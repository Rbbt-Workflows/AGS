require 'set'

module AGS

  #{{{ CONSTANTS

  DAG_COMPARISONS = {
    'pdpi_vs_singles'  => ['INT_PD_PI', 'PD', 'PI'],
    'fivezpi_vs_singles' => ['INT_FiveZ_PI', 'FiveZ', 'PI'],
    'pdpi_vs_fivezpi'  => ['INT_PD_PI', 'INT_FiveZ_PI'],
  }

  DAG_COMPARISON_LABELS = {
    'pdpi_vs_singles'  => 'PD+PI vs singles',
    'fivezpi_vs_singles' => '5Z+PI vs singles',
    'pdpi_vs_fivezpi'  => 'PD+PI vs 5Z+PI',
  }

  DAG_TIMEPOINT_ORDER = ['T1', 'T2', 'T4', 'T8', 'T24']

  DAG_TP_X = {
    'T1'  => 0,
    'T2'  => 1,
    'T4'  => 2,
    'T8'  => 3,
    'T24' => 4,
  }

  DAG_TP_RANK = {
    'T1'  => 0,
    'T2'  => 1,
    'T4'  => 2,
    'T8'  => 3,
    'T24' => 4,
  }

  DAG_ELIGIBLE_PAIRS = Set.new([
    ["T1", "T2"], ["T1", "T4"],
    ["T2", "T2"], ["T2", "T4"],
    ["T4", "T4"], ["T4", "T8"],
    ["T8", "T8"], ["T8", "T24"],
  ])

  DAG_SAME_TIMEPOINT_PAIRS = Set.new([
    ["T2", "T2"], ["T4", "T4"], ["T8", "T8"],
  ])

  DAG_TREATMENTS_FOR_COMPARISONS = DAG_COMPARISONS.values.flatten.uniq

  #{{{ HELPERS: edge extraction and pruning

  helper :dag_extract_filtered_edges do |seq, activity_cutoff, include_same_timepoint|
    seq = seq.to_double unless seq.type == :double
    edges = []
    seq.through do |key, values|
      source      = values["Source"].first
      source_tp   = values["Source timepoint"].first
      source_act  = values["Source activity"].first.to_f
      target      = values["Target"].first
      target_tp   = values["Target timepoint"].first
      target_act  = values["Target activity"].first.to_f
      source_sc   = values["Source self-consistent"].first.to_s == "true" || values["Source self-consistent"].first.to_s == "1"
      target_sc   = values["Target self-consistent"].first.to_s == "true" || values["Target self-consistent"].first.to_s == "1"
      effect      = values["Effect"].first.to_f
      offset      = values["Offset"].first.to_i
      type_str    = values["Type"].first

      # Step A: eligible pair check
      next unless DAG_ELIGIBLE_PAIRS.include?([source_tp, target_tp])

      # Optional: exclude same-timepoint sequences
      if !include_same_timepoint && DAG_SAME_TIMEPOINT_PAIRS.include?([source_tp, target_tp])
        next
      end

      # Step B: self-consistency rules
      next unless target_sc
      if source_tp == "T4" || source_tp == "T8"
        next unless source_sc
      end

      # Step B: activity cutoff
      next unless source_act.abs >= activity_cutoff && target_act.abs >= activity_cutoff

      edges << {
        key: key,
        source: source, source_tp: source_tp, source_act: source_act,
        target: target, target_tp: target_tp, target_act: target_act,
        source_sc: source_sc, target_sc: target_sc,
        effect: effect, offset: offset, type_str: type_str,
      }
    end
    edges
  end

  # Prune edges to unbroken chains from T1/T2 to T8/T24.
  # Forward pass: ensure every edge source is reachable from a T1/T2 start node.
  # Backward pass (optional): ensure every edge target can reach a T8/T24 end node.
  # Returns a Set of surviving edge indices.
  helper :dag_prune do |edges, prune_dead_ends|
    # Build adjacency lists keyed on [tf, timepoint]
    forward_adj = {}
    reverse_adj = {}

    edges.each_with_index do |e, i|
      s_node = [e[:source], e[:source_tp]]
      t_node = [e[:target], e[:target_tp]]
      (forward_adj[s_node] ||= []) << i
      (reverse_adj[t_node] ||= []) << i
    end

    # --- Forward reachability from T1/T2 source nodes ---
    start_nodes = Set.new
    edges.each do |e|
      if ['T1', 'T2'].include?(e[:source_tp])
        start_nodes << [e[:source], e[:source_tp]]
      end
    end

    fwd_reach = Set.new(start_nodes)
    queue = start_nodes.to_a
    until queue.empty?
      node = queue.shift
      (forward_adj[node] || []).each do |ei|
        t_node = [edges[ei][:target], edges[ei][:target_tp]]
        unless fwd_reach.include?(t_node)
          fwd_reach << t_node
          queue << t_node
        end
      end
    end

    if prune_dead_ends
      # --- Backward reachability from T8/T24 target nodes ---
      end_nodes = Set.new
      edges.each do |e|
        if ['T8', 'T24'].include?(e[:target_tp])
          end_nodes << [e[:target], e[:target_tp]]
        end
      end

      bwd_reach = Set.new(end_nodes)
      queue = end_nodes.to_a
      until queue.empty?
        node = queue.shift
        (reverse_adj[node] || []).each do |ei|
          s_node = [edges[ei][:source], edges[ei][:source_tp]]
          unless bwd_reach.include?(s_node)
            bwd_reach << s_node
            queue << s_node
          end
        end
      end

      # Edge survives if source is forward-reachable AND target is backward-reachable
      surviving = Set.new
      edges.each_with_index do |e, i|
        s_node = [e[:source], e[:source_tp]]
        t_node = [e[:target], e[:target_tp]]
        surviving << i if fwd_reach.include?(s_node) && bwd_reach.include?(t_node)
      end
    else
      # Forward-only: keep edge if source is forward-reachable
      surviving = Set.new
      edges.each_with_index do |e, i|
        s_node = [e[:source], e[:source_tp]]
        surviving << i if fwd_reach.include?(s_node)
      end
    end

    surviving
  end

  #{{{ TASK: dag_sequences (per-treatment unbroken DAG)

  desc "Generate unbroken DAG chains from T1/T2 to T8/T24 targets with filtering and pruning"
  dep :sequence, jobname: 'Default'
  input :activity_cutoff, :float, "Minimum absolute activity for both source and target TFs", 3.0
  input :prune_dead_ends, :boolean, "Remove edges that do not lead to T8/T24 endpoints", true
  input :include_same_timepoint, :boolean, "Include same-timepoint sequences (T2-T2, T4-T4, T8-T8)", true
  task :dag_sequences => :tsv do |activity_cutoff, prune_dead_ends, include_same_timepoint|
    treatment = recursive_inputs[:treatment]
    seq = step(:sequence).load
    seq = seq.to_double unless seq.type == :double

    edges = dag_extract_filtered_edges(seq, activity_cutoff, include_same_timepoint)
    surviving = dag_prune(edges, prune_dead_ends)

    result = TSV.setup({}, key_field: "ID", fields: [
      'Source', 'Source timepoint', 'Source activity',
      'Target', 'Target timepoint', 'Target activity',
      'Effect', 'Offset', 'Type',
      'Source self-consistent', 'Target self-consistent'
    ], type: :list)

    surviving.sort.each do |i|
      e = edges[i]
      key = "dag-#{treatment}-#{e[:source]}_#{e[:source_tp]}-#{e[:target]}_#{e[:target_tp]}"
      result[key] = [
        e[:source], e[:source_tp], e[:source_act],
        e[:target], e[:target_tp], e[:target_act],
        e[:effect], e[:offset], e[:type_str],
        e[:source_sc], e[:target_sc]
      ]
    end

    result
  end

  dep :dag_sequences do |jobname,options|
    TREATMENTS.collect{|treatment| {treatment: treatment} }
  end
  task :dag_suite => :array do
    dependencies.each do |dep|
      treatment = dep.recursive_inputs[:treatment]
      activity_cutoff = dep.recursive_inputs[:activity_cutoff]

      target = file("#{treatment}-#{activity_cutoff}.tsv")
      Open.cp dep.path, target
    end
    files
  end


  #{{{ TASK: dag_common_edges (comparison common edges + D3 filtering)

  desc "Find common DAG edges across treatments for each comparison, with D3 time-compatibility filtering"
  dep :dag_sequences, treatment: :placeholder do |jobname, options|
    DAG_TREATMENTS_FOR_COMPARISONS.collect do |treatment|
      options.merge(treatment: treatment)
    end
  end
  input :activity_cutoff, :float, "Minimum absolute activity for both source and target TFs", 3.0
  input :prune_dead_ends, :boolean, "Remove edges that do not lead to T8/T24 endpoints", true
  input :include_same_timepoint, :boolean, "Include same-timepoint sequences (T2-T2, T4-T4, T8-T8)", true
  task :dag_common_edges => :tsv do |activity_cutoff, prune_dead_ends, include_same_timepoint|
    # Collect DAG outputs per treatment
    treatment_dags = {}
    dag_deps = dependencies.select { |d| d.task_name == :dag_sequences }
    dag_deps.each do |dep|
      treatment = dep.recursive_inputs[:treatment]
      tsv = dep.load
      tsv = tsv.to_double unless tsv.type == :double
      treatment_dags[treatment] = tsv
    end

    comparison_fields = ['Source', 'Source timepoint', 'Target', 'Target timepoint', 'Effect', 'Comparison']
    # Add per-treatment activity columns
    all_comparison_treatments = DAG_COMPARISONS.values.flatten.uniq.sort
    all_comparison_treatments.each do |t|
      comparison_fields << "#{t} source activity"
      comparison_fields << "#{t} target activity"
    end
    comparison_fields << 'Chain step'

    result = TSV.setup({}, key_field: "ID", fields: comparison_fields, type: :list)

    DAG_COMPARISONS.each do |comp_name, comp_treatments|
      # For each treatment, build a set of edge tuples and a lookup of activities
      treatment_edge_maps = {}
      comp_treatments.each do |t|
        tsv = treatment_dags[t]
        next unless tsv
        edge_map = {}
        tsv.through do |key, values|
          tuple = [
            values["Source"].first,
            values["Source timepoint"].first,
            values["Target"].first,
            values["Target timepoint"].first
          ]
          edge_map[tuple] = {
            source_act: values["Source activity"].first.to_f,
            target_act: values["Target activity"].first.to_f,
            effect: values["Effect"].first.to_f,
          }
        end
        treatment_edge_maps[t] = edge_map
      end

      # Find intersection of edge tuples across all treatments
      if treatment_edge_maps.empty?
        next
      end

      common_tuples = treatment_edge_maps.values.first.keys.to_set
      treatment_edge_maps.values[1..-1].each do |em|
        common_tuples &= em.keys.to_set
      end

      # Convert common tuples to edge hashes for pruning
      common_edges = common_tuples.map do |tuple|
        source, source_tp, target, target_tp = tuple
        first_map = treatment_edge_maps[comp_treatments.first]
        info = first_map[tuple]
        {
          source: source, source_tp: source_tp,
          target: target, target_tp: target_tp,
          source_act: info[:source_act],
          target_act: info[:target_act],
          effect: info[:effect],
        }
      end

      # D3: apply same pruning to common edges
      surviving_indices = dag_prune(common_edges, true)

      surviving_indices.sort.each do |i|
        e = common_edges[i]
        tuple = [e[:source], e[:source_tp], e[:target], e[:target_tp]]

        # Determine chain step label
        step_label = "#{e[:target_tp]}<-#{e[:source_tp]}"

        # Collect per-treatment activities
        row = [e[:source], e[:source_tp], e[:target], e[:target_tp], e[:effect], comp_name]
        all_comparison_treatments.each do |t|
          if treatment_edge_maps[t] && treatment_edge_maps[t][tuple]
            row << treatment_edge_maps[t][tuple][:source_act]
            row << treatment_edge_maps[t][tuple][:target_act]
          else
            row << ""
            row << ""
          end
        end
        row << step_label

        key = "#{comp_name}-#{e[:source]}_#{e[:source_tp]}-#{e[:target]}_#{e[:target_tp]}"
        result[key] = row
      end
    end

    result
  end

  #{{{ TASK: dag_common_edges_by_tf (time-agnostic common edges)

  desc "Find common DAG edges across treatments matched by TF-pair identity (ignoring timepoint), recording per-treatment timing"
  dep :dag_sequences, treatment: :placeholder do |jobname, options|
    DAG_TREATMENTS_FOR_COMPARISONS.collect do |treatment|
      options.merge(treatment: treatment)
    end
  end
  input :activity_cutoff, :float, "Minimum absolute activity for both source and target TFs", 3.0
  input :prune_dead_ends, :boolean, "Remove edges that do not lead to T8/T24 endpoints", true
  input :include_same_timepoint, :boolean, "Include same-timepoint sequences (T2-T2, T4-T4, T8-T8)", true
  task :dag_common_edges_by_tf => :tsv do |activity_cutoff, prune_dead_ends, include_same_timepoint|
    # Collect DAG outputs per treatment
    treatment_dags = {}
    dag_deps = dependencies.select { |d| d.task_name == :dag_sequences }
    dag_deps.each do |dep|
      treatment = dep.recursive_inputs[:treatment]
      tsv = dep.load
      tsv = tsv.to_double unless tsv.type == :double
      treatment_dags[treatment] = tsv
    end

    all_comparison_treatments = DAG_COMPARISONS.values.flatten.uniq.sort

    out_fields = ['Source', 'Target', 'Effect', 'Comparison']
    all_comparison_treatments.each do |t|
      out_fields << "#{t} timepoints"
      out_fields << "#{t} source activities"
      out_fields << "#{t} target activities"
    end
    out_fields << 'Earliest single timepoint'
    out_fields << 'Earliest combo timepoint'
    out_fields << 'Classification'

    result = TSV.setup({}, key_field: "Edge ID", fields: out_fields, type: :list)

    DAG_COMPARISONS.each do |comp_name, comp_treatments|
      # For each treatment, build a map keyed on (Source TF, Target TF) => array of timing entries
      treatment_tf_edges = {}
      comp_treatments.each do |t|
        tsv = treatment_dags[t]
        next unless tsv
        tf_map = Hash.new { |h, k| h[k] = [] }
        tsv.through do |key, values|
          source    = values["Source"].first
          source_tp = values["Source timepoint"].first
          target    = values["Target"].first
          target_tp = values["Target timepoint"].first
          source_act = values["Source activity"].first.to_f
          target_act = values["Target activity"].first.to_f
          effect     = values["Effect"].first.to_f

          tf_pair = [source, target]
          tf_map[tf_pair] << {
            source_tp: source_tp,
            target_tp: target_tp,
            source_act: source_act,
            target_act: target_act,
            effect: effect,
          }
        end
        treatment_tf_edges[t] = tf_map
      end

      next if treatment_tf_edges.empty?

      # Find intersection of TF-pair keys across all treatments
      common_tf_pairs = treatment_tf_edges.values.first.keys.to_set
      treatment_tf_edges.values[1..-1].each do |tf_map|
        common_tf_pairs &= tf_map.keys.to_set
      end

      # Identify which treatments are "singles" vs "combination"
      combo_treatments = comp_treatments.select { |t| t.start_with?('INT_') }
      single_treatments = comp_treatments.reject { |t| t.start_with?('INT_') }

      common_tf_pairs.sort.each do |tf_pair|
        source, target = tf_pair

        # Gather per-treatment data
        per_treatment_timepoints = {}
        per_treatment_source_acts = {}
        per_treatment_target_acts = {}
        effect_val = nil

        comp_treatments.each do |t|
          entries = treatment_tf_edges[t][tf_pair] || []
          if entries.any?
            effect_val ||= entries.first[:effect]
            per_treatment_timepoints[t] = entries.map { |e| "#{e[:source_tp]}->#{e[:target_tp]}" }.sort.uniq
            per_treatment_source_acts[t] = entries.map { |e| "%.2f" % e[:source_act] }
            per_treatment_target_acts[t] = entries.map { |e| "%.2f" % e[:target_act] }
          else
            per_treatment_timepoints[t] = []
            per_treatment_source_acts[t] = []
            per_treatment_target_acts[t] = []
          end
        end

        # Compute earliest timepoint for classification
        tp_rank = lambda do |tp_str|
          parts = tp_str.split('->')
          src_rank = DAG_TP_RANK[parts[0]] || 99
          tgt_rank = DAG_TP_RANK[parts[1]] || 99
          src_rank + tgt_rank / 10.0
        end

        earliest_single = nil
        single_treatments.each do |t|
          per_treatment_timepoints[t].each do |tp_str|
            r = tp_rank.call(tp_str)
            earliest_single = r if earliest_single.nil? || r < earliest_single
          end
        end

        earliest_combo = nil
        combo_treatments.each do |t|
          per_treatment_timepoints[t].each do |tp_str|
            r = tp_rank.call(tp_str)
            earliest_combo = r if earliest_combo.nil? || r < earliest_combo
          end
        end

        # Classification
        classification = if earliest_single.nil? && earliest_combo.nil?
                           'absent'
                         elsif earliest_single.nil?
                           'combination_specific'
                         elsif earliest_combo.nil?
                           'single_only'
                         elsif earliest_combo < earliest_single
                           'combo_before_singles'
                         elsif earliest_combo > earliest_single
                           'single_before_combo'
                         else
                           'same_time'
                         end

        # Build row
        row = [source, target, effect_val, comp_name]
        all_comparison_treatments.each do |t|
          row << (per_treatment_timepoints[t] || []).join(';')
          row << (per_treatment_source_acts[t] || []).join(';')
          row << (per_treatment_target_acts[t] || []).join(';')
        end

        # Earliest labels
        earliest_single_str = single_treatments.map do |t|
          next if per_treatment_timepoints[t].empty?
          "#{t}:#{per_treatment_timepoints[t].min_by { |tp_str| tp_rank.call(tp_str) }}"
        end.compact.join(';')
        earliest_combo_str = combo_treatments.map do |t|
          next if per_treatment_timepoints[t].empty?
          "#{t}:#{per_treatment_timepoints[t].min_by { |tp_str| tp_rank.call(tp_str) }}"
        end.compact.join(';')

        row << (earliest_single_str || '')
        row << (earliest_combo_str || '')
        row << classification

        key = "#{comp_name}-#{source}-#{target}"
        result[key] = row
      end
    end

    result
  end

  #{{{ TASK: dag_combination_specific_edges (edges in combo but not in any single)

  desc "Find TF-pair edges present in the combination treatment but absent from ALL single treatments"
  dep :dag_sequences, treatment: :placeholder do |jobname, options|
    DAG_TREATMENTS_FOR_COMPARISONS.collect do |treatment|
      options.merge(treatment: treatment)
    end
  end
  input :activity_cutoff, :float, "Minimum absolute activity for both source and target TFs", 3.0
  input :prune_dead_ends, :boolean, "Remove edges that do not lead to T8/T24 endpoints", true
  input :include_same_timepoint, :boolean, "Include same-timepoint sequences (T2-T2, T4-T4, T8-T8)", true
  task :dag_combination_specific_edges => :tsv do |activity_cutoff, prune_dead_ends, include_same_timepoint|
    treatment_dags = {}
    dag_deps = dependencies.select { |d| d.task_name == :dag_sequences }
    dag_deps.each do |dep|
      treatment = dep.recursive_inputs[:treatment]
      tsv = dep.load
      tsv = tsv.to_double unless tsv.type == :double
      treatment_dags[treatment] = tsv
    end

    out_fields = ['Source', 'Target', 'Effect', 'Combination', 'Combination timepoints',
                  'Combination source activities', 'Combination target activities']
    result = TSV.setup({}, key_field: "Edge ID", fields: out_fields, type: :list)

    DAG_COMPARISONS.each do |comp_name, comp_treatments|
      combo_treatments = comp_treatments.select { |t| t.start_with?('INT_') }
      single_treatments = comp_treatments.reject { |t| t.start_with?('INT_') }

      # Build TF-pair sets for singles
      single_tf_pairs = Set.new
      single_treatments.each do |t|
        tsv = treatment_dags[t]
        next unless tsv
        tsv.through do |key, values|
          single_tf_pairs << [values["Source"].first, values["Target"].first]
        end
      end

      # Find combo edges not in singles
      combo_treatments.each do |t|
        tsv = treatment_dags[t]
        next unless tsv
        tf_map = Hash.new { |h, k| h[k] = [] }
        tsv.through do |key, values|
          source    = values["Source"].first
          source_tp = values["Source timepoint"].first
          target    = values["Target"].first
          target_tp = values["Target timepoint"].first
          source_act = values["Source activity"].first.to_f
          target_act = values["Target activity"].first.to_f
          effect     = values["Effect"].first.to_f

          tf_pair = [source, target]
          next if single_tf_pairs.include?(tf_pair)

          tf_map[tf_pair] << {
            source_tp: source_tp,
            target_tp: target_tp,
            source_act: source_act,
            target_act: target_act,
            effect: effect,
          }
        end

        tf_map.each do |tf_pair, entries|
          source, target = tf_pair
          effect_val = entries.first[:effect]
          timepoints = entries.map { |e| "#{e[:source_tp]}->#{e[:target_tp]}" }.sort.uniq
          source_acts = entries.map { |e| "%.2f" % e[:source_act] }
          target_acts = entries.map { |e| "%.2f" % e[:target_act] }

          key = "#{comp_name}-#{t}-#{source}-#{target}"
          result[key] = [source, target, effect_val, t, timepoints.join(';'), source_acts.join(';'), target_acts.join(';')]
        end
      end
    end

    result
  end

  #{{{ TASK: dag_sif (SIF export for individual + comparison DAGs)

  desc "Export DAG networks (per-treatment and comparison) in SIF format"
  dep :dag_sequences, treatment: :placeholder do |jobname, options|
    DAG_TREATMENTS_FOR_COMPARISONS.collect do |treatment|
      options.merge(treatment: treatment)
    end
  end
  dep :dag_common_edges do |jobname, options|
    options
  end
  input :activity_cutoff, :float, "Minimum absolute activity for both source and target TFs", 3.0
  input :prune_dead_ends, :boolean, "Remove edges that do not lead to T8/T24 endpoints", true
  input :include_same_timepoint, :boolean, "Include same-timepoint sequences (T2-T2, T4-T4, T8-T8)", true
  task :dag_sif => :text do |activity_cutoff, prune_dead_ends, include_same_timepoint|
    lines = []

    # Per-treatment DAGs
    dag_deps = dependencies.select { |d| d.task_name == :dag_sequences }
    dag_deps.each do |dep|
      treatment = dep.recursive_inputs[:treatment]
      tsv = dep.load
      tsv = tsv.to_double unless tsv.type == :double
      lines << "# Treatment: #{treatment}"
      tsv.through do |key, values|
        source = values["Source"].first
        target = values["Target"].first
        effect = values["Effect"].first.to_f
        interaction = effect > 0 ? "activates" : "represses"
        lines << [source, interaction, target] * "\t"
      end
      lines << ""
    end

    # Comparison DAGs
    comp_tsv = step(:dag_common_edges).load
    comp_tsv = comp_tsv.to_double unless comp_tsv.type == :double

    current_comp = nil
    comp_tsv.through do |key, values|
      comp = values["Comparison"].first
      if comp != current_comp
        lines << "" if current_comp
        lines << "# Comparison: #{comp}"
        current_comp = comp
      end
      source = values["Source"].first
      target = values["Target"].first
      effect = values["Effect"].first.to_f
      interaction = effect > 0 ? "activates" : "represses"
      lines << [source, interaction, target] * "\t"
    end

    sif_text = lines * "\n" + "\n"
    Open.write(file("dag_network.sif"), sif_text)
    sif_text
  end

  #{{{ TASK: dag_layout (layered layout)

  desc "Compute a fixed layered node layout for DAGs (left-to-right by timepoint)"
  dep :dag_sequences, treatment: :placeholder do |jobname, options|
    DAG_TREATMENTS_FOR_COMPARISONS.collect do |treatment|
      options.merge(treatment: treatment)
    end
  end
  dep :dag_common_edges do |jobname, options|
    options
  end
  input :activity_cutoff, :float, "Minimum absolute activity for both source and target TFs", 3.0
  input :prune_dead_ends, :boolean, "Remove edges that do not lead to T8/T24 endpoints", true
  input :include_same_timepoint, :boolean, "Include same-timepoint sequences (T2-T2, T4-T4, T8-T8)", true
  task :dag_layout => :tsv do |activity_cutoff, prune_dead_ends, include_same_timepoint|
    # Collect all unique nodes (TF@timepoint) from all DAGs and comparisons
    all_nodes = Set.new

    dag_deps = dependencies.select { |d| d.task_name == :dag_sequences }
    dag_deps.each do |dep|
      tsv = dep.load
      tsv = tsv.to_double unless tsv.type == :double
      tsv.through do |key, values|
        source = values["Source"].first
        source_tp = values["Source timepoint"].first
        target = values["Target"].first
        target_tp = values["Target timepoint"].first
        all_nodes << "#{source}@#{source_tp}"
        all_nodes << "#{target}@#{target_tp}"
      end
    end

    comp_tsv = step(:dag_common_edges).load
    comp_tsv = comp_tsv.to_double unless comp_tsv.type == :double
    comp_tsv.through do |key, values|
      source = values["Source"].first
      source_tp = values["Source timepoint"].first
      target = values["Target"].first
      target_tp = values["Target timepoint"].first
      all_nodes << "#{source}@#{source_tp}"
      all_nodes << "#{target}@#{target_tp}"
    end

    # Build layered layout: group nodes by timepoint, sort within layer
    layers = {}
    all_nodes.each do |node_id|
      tf, tp = node_id.split('@')
      (layers[tp] ||= []) << tf
    end

    layout = TSV.setup({}, key_field: "Node", fields: ["X", "Y"], type: :list)

    # Assign coordinates: X by timepoint column, Y by alphabetical order within column
    col_width = 200.0
    row_height = 40.0
    max_layer_size = layers.values.map(&:length).max || 1

    DAG_TIMEPOINT_ORDER.each_with_index do |tp, col_idx|
      tfs = (layers[tp] || []).sort
      tfs.each_with_index do |tf, row_idx|
        node_id = "#{tf}@#{tp}"
        x = col_idx * col_width
        # Center the layer vertically
        y = (max_layer_size - tfs.length) / 2.0 * row_height + row_idx * row_height
        layout[node_id] = [x, y]
      end
    end

    Open.write(file("dag_layout.tsv"), layout.to_s)
    layout
  end

  #{{{ SVG HELPERS

  helper :dag_activity_to_color do |activity|
    return "#F5F5F5" unless activity && activity != ""

    activity = activity.to_f
    threshold = 2.5

    if activity.abs >= threshold
      if activity > 0
        intensity = [[activity.abs / 8.0, 1.0].min, 0.3].max
        r = (215 + (127 - 215) * intensity).round
        g = (48 + (0 - 48) * intensity).round
        b = (39 + (0 - 39) * intensity).round
        sprintf("#%02X%02X%02X", r, g, b)
      else
        intensity = [[activity.abs / 8.0, 1.0].min, 0.3].max
        r = (69 + (8 - 69) * intensity).round
        g = (117 + (30 - 117) * intensity).round
        b = (180 + (107 - 180) * intensity).round
        sprintf("#%02X%02X%02X", r, g, b)
      end
    else
      "#F5F5F5"
    end
  end

  helper :build_dag_svg do |layout, edges, node_activities, title, svg_w, svg_h|
    # Collect all nodes from edges and layout
    all_nodes = Set.new
    edges.each do |e|
      all_nodes << "#{e[:source]}@#{e[:source_tp]}"
      all_nodes << "#{e[:target]}@#{e[:target_tp]}"
    end

    # Determine bounding box from layout
    xs = all_nodes.map { |n| layout[n] ? layout[n][0].first.to_f : 0.0 }
    ys = all_nodes.map { |n| layout[n] ? layout[n][1].first.to_f : 0.0 }
    min_x, max_x = xs.minmax
    min_y, max_y = ys.minmax

    margin = 80
    node_r = 16.0
    range_x = (max_x - min_x); range_x = 1.0 if range_x == 0
    range_y = (max_y - min_y); range_y = 1.0 if range_y == 0

    inner_w = svg_w - 2 * margin
    inner_h = svg_h - 2 * margin - 30  # extra top margin for title

    # Compute screen positions
    pos = {}
    all_nodes.each do |n|
      lx = layout[n] ? layout[n][0].first.to_f : 0.0
      ly = layout[n] ? layout[n][1].first.to_f : 0.0
      sx = margin + ((lx - min_x) / range_x) * inner_w
      sy = margin + 30 + ((max_y - ly) / range_y) * inner_h
      pos[n] = [sx, sy]
    end

    # Helper: clip line endpoints to node boundaries
    clip_edge = lambda do |sx, sy, tx, ty, r|
      dx = tx - sx
      dy = ty - sy
      dist = Math.sqrt(dx * dx + dy * dy)
      return [sx, sy, tx, ty] if dist == 0
      ux = dx / dist
      uy = dy / dist
      [sx + ux * r, sy + uy * r, tx - ux * r, ty - uy * r]
    end

    svg = []
    svg << %(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)
    svg << %(<svg xmlns="http://www.w3.org/2000/svg" width="#{svg_w}" height="#{svg_h}" viewBox="0 0 #{svg_w} #{svg_h}">)
    svg << %(<rect width="#{svg_w}" height="#{svg_h}" fill="white"/>)

    # Title
    svg << %(<text x="#{svg_w / 2}" y="22" text-anchor="middle" font-family="Helvetica" font-size="16" font-weight="bold">#{title}</text>)

    # Timepoint column headers
    DAG_TIMEPOINT_ORDER.each do |tp|
      tp_nodes = all_nodes.select { |n| n.end_with?("@#{tp}") }
      next if tp_nodes.empty?
      sample_node = tp_nodes.first
      next unless pos[sample_node]
      x = pos[sample_node][0]
      svg << %(<text x="#{x.round(1)}" y="#{(margin + 18).round(1)}" text-anchor="middle" font-family="Helvetica" font-size="12" fill="#999">#{tp}</text>)
    end

    # Arrow marker definitions
    svg << %(<defs>)
    svg << %(<marker id="arrow-act-dag" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto"><polygon points="0,0 10,4 0,8" fill="#D73027"/></marker>)
    svg << %(<marker id="arrow-rep-dag" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto"><polygon points="0,0 10,4 0,8" fill="#4575B4"/></marker>)
    svg << %(</defs>)

    # Draw edges
    edges.each do |e|
      s_node = "#{e[:source]}@#{e[:source_tp]}"
      t_node = "#{e[:target]}@#{e[:target_tp]}"
      next unless pos[s_node] && pos[t_node]

      sx, sy = pos[s_node]
      tx, ty = pos[t_node]
      csx, csy, ctx, cty = clip_edge.call(sx, sy, tx, ty, node_r)

      effect = e[:effect].to_f
      color = effect > 0 ? "#D73027" : "#4575B4"
      marker_id = effect > 0 ? "arrow-act-dag" : "arrow-rep-dag"

      svg << %(<line x1="#{csx.round(1)}" y1="#{csy.round(1)}" x2="#{ctx.round(1)}" y2="#{cty.round(1)}" stroke="#{color}" stroke-width="1.5" marker-end="url(##{marker_id})"/>)
    end

    # Draw nodes
    all_nodes.sort.each do |n|
      next unless pos[n]
      sx, sy = pos[n]
      tf, tp = n.split('@')

      activity = node_activities[n]
      fill = dag_activity_to_color(activity)

      svg << %(<circle cx="#{sx.round(1)}" cy="#{sy.round(1)}" r="#{node_r}" fill="#{fill}" stroke="#333" stroke-width="0.8"/>)
      svg << %(<text x="#{sx.round(1)}" y="#{(sy + node_r + 11).round(1)}" text-anchor="middle" font-family="Helvetica" font-size="8">#{tf}</text>)

      # Show activity value if present
      if activity && activity != ""
        act_str = "%.1f" % activity.to_f
        svg << %(<text x="#{sx.round(1)}" y="#{(sy + 3).round(1)}" text-anchor="middle" font-family="Helvetica" font-size="7" fill="white" font-weight="bold">#{act_str}</text>)
      end
    end

    svg << %(</svg>)
    svg * "\n"
  end

  #{{{ TASK: dag_svg_panels (individual treatment + comparison SVGs)

  desc "Generate SVG panels for individual treatment DAGs and comparison DAGs using a fixed layout"
  dep :dag_sequences, treatment: :placeholder do |jobname, options|
    DAG_TREATMENTS_FOR_COMPARISONS.collect do |treatment|
      options.merge(treatment: treatment)
    end
  end
  dep :dag_common_edges do |jobname, options|
    options
  end
  dep :dag_layout do |jobname, options|
    options
  end
  dep :tf_predictions, jobname: 'Default'
  input :activity_cutoff, :float, "Minimum absolute activity for both source and target TFs", 3.0
  input :prune_dead_ends, :boolean, "Remove edges that do not lead to T8/T24 endpoints", true
  input :include_same_timepoint, :boolean, "Include same-timepoint sequences (T2-T2, T4-T4, T8-T8)", true
  task :dag_svg_panels => :array do |activity_cutoff, prune_dead_ends, include_same_timepoint|
    layout = step(:dag_layout).load
    layout = layout.to_double unless layout.type == :double
    preds  = step(:tf_predictions).load
    preds  = preds.to_double unless preds.type == :double

    svg_files = []

    # --- Per-treatment DAG panels ---
    dag_deps = dependencies.select { |d| d.task_name == :dag_sequences }
    dag_deps.each do |dep|
      treatment = dep.recursive_inputs[:treatment]
      tsv = dep.load
      tsv = tsv.to_double unless tsv.type == :double

      edges = []
      tsv.through do |key, values|
        edges << {
          source: values["Source"].first,
          source_tp: values["Source timepoint"].first,
          target: values["Target"].first,
          target_tp: values["Target timepoint"].first,
          effect: values["Effect"].first.to_f,
        }
      end

      # Get activity values for this treatment
      node_activities = {}
      DAG_TIMEPOINT_ORDER.each do |tp|
        col = "#{treatment}-T#{tp.sub('T','')}"
        col_idx = preds.fields.index(col)
        next unless col_idx
        preds.each do |tf, values|
          val = values[col_idx].first rescue values[col_idx]
          val = val.to_s if val
          node_activities["#{tf}@#{tp}"] = val if val && val != ""
        end
      end

      label = AGS.figure_treatment_label(treatment.to_s)
      svg = build_dag_svg(layout, edges, node_activities, label, 700, 500)
      filename = "dag_#{treatment}.svg"
      Open.write(file(filename), svg)
      svg_files << file(filename).to_s
    end

    # --- Comparison DAG panels ---
    comp_tsv = step(:dag_common_edges).load
    comp_tsv = comp_tsv.to_double unless comp_tsv.type == :double

    DAG_COMPARISONS.each do |comp_name, comp_treatments|
      edges = []
      node_activities = {}

      comp_tsv.through do |key, values|
        next unless values["Comparison"].first == comp_name

        source = values["Source"].first
        source_tp = values["Source timepoint"].first
        target = values["Target"].first
        target_tp = values["Target timepoint"].first
        effect = values["Effect"].first.to_f

        edges << {
          source: source, source_tp: source_tp,
          target: target, target_tp: target_tp,
          effect: effect,
        }

        # Average activity across treatments for coloring
        comp_treatments.each do |t|
          sa_field = "#{t} source activity"
          ta_field = "#{t} target activity"
          sa = values[sa_field] ? values[sa_field].first : nil
          ta = values[ta_field] ? values[ta_field].first : nil

          s_node = "#{source}@#{source_tp}"
          t_node = "#{target}@#{target_tp}"

          if sa && sa != ""
            node_activities[s_node] ||= []
            node_activities[s_node] << sa.to_f
          end
          if ta && ta != ""
            node_activities[t_node] ||= []
            node_activities[t_node] << ta.to_f
          end
        end
      end

      next if edges.empty?

      # Compute average activity for each node
      avg_activities = {}
      node_activities.each do |node, vals|
        avg_activities[node] = vals.sum / vals.length if vals.any?
      end

      label = DAG_COMPARISON_LABELS[comp_name] || comp_name
      svg = build_dag_svg(layout, edges, avg_activities, label, 700, 500)
      filename = "dag_comp_#{comp_name}.svg"
      Open.write(file(filename), svg)
      svg_files << file(filename).to_s
    end

    Open.write(file("dag_index.txt"), svg_files * "\n" + "\n")
    svg_files
  end

  #{{{ TASK: dag_grid_svg (combined grid visualization)

  desc "Combine all DAG SVG panels (individual treatments + comparisons) into a single grid SVG"
  dep :dag_svg_panels do |jobname, options|
    options
  end
  input :activity_cutoff, :float, "Minimum absolute activity for both source and target TFs", 3.0
  input :prune_dead_ends, :boolean, "Remove edges that do not lead to T8/T24 endpoints", true
  input :include_same_timepoint, :boolean, "Include same-timepoint sequences (T2-T2, T4-T4, T8-T8)", true
  task :dag_grid_svg => :binary do |activity_cutoff, prune_dead_ends, include_same_timepoint|
    panel_files = step(:dag_svg_panels).load

    # Separate into treatment and comparison panels
    treatment_panels = {}
    comparison_panels = {}

    panel_files.each do |path|
      basename = File.basename(path, ".svg")
      if basename =~ /^dag_comp_(.+)$/
        comparison_panels[$1] = path
      elsif basename =~ /^dag_(.+)$/
        treatment_panels[$1] = path
      end
    end

    # Layout: treatment panels in rows, comparison panels below
    treatment_order = DAG_COMPARISONS.values.flatten.uniq
    comparison_order = DAG_COMPARISONS.keys

    panels = []
    treatment_order.each { |t| panels << treatment_panels[t] }
    comparison_order.each { |c| panels << comparison_panels[c] }
    panels.compact!

    panel_w = 700
    panel_h = 500
    margin = 40
    title_h = 25
    cell_w = panel_w + margin
    cell_h = panel_h + margin + title_h

    # Arrange in a grid: up to 3 columns
    cols = [panels.length, 3].min
    rows = (panels.length / cols.to_f).ceil

    total_w = (cell_w * cols).round
    total_h = (cell_h * rows).round

    grid_parts = []
    grid_parts << %(<?xml version="1.0" encoding="UTF-8" standalone="no"?>)
    grid_parts << %(<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="#{total_w}" height="#{total_h}" viewBox="0 0 #{total_w} #{total_h}">)
    grid_parts << %(<rect width="#{total_w}" height="#{total_h}" fill="white"/>)

    panels.each_with_index do |path, i|
      row = i / cols
      col = i % cols
      x = col * cell_w
      y = row * cell_h

      basename = File.basename(path, ".svg")
      if basename =~ /^dag_comp_(.+)$/
        label = DAG_COMPARISON_LABELS[$1] || $1
      elsif basename =~ /^dag_(.+)$/
        label = AGS.figure_treatment_label($1)
      else
        label = basename
      end

      grid_parts << %(<text x="#{x + cell_w / 2}" y="#{y + 18}" text-anchor="middle" font-family="Helvetica" font-size="14" font-weight="bold">#{label}</text>)

      if path && File.exist?(path)
        inner_svg = Open.read(path)
        inner_content = inner_svg
          .sub(/<\?xml[^>]*\?>/, "")
          .sub(/^<svg[^>]*>/, "")
          .sub(/<\/svg>\s*$/, "")

        offset_x = x + margin / 2
        offset_y = y + title_h

        grid_parts << %(<g transform="translate(#{offset_x}, #{offset_y})">)
        grid_parts << inner_content
        grid_parts << %(</g>)
      end
    end

    grid_parts << %(</svg>)

    grid_svg = grid_parts * "\n"
    Open.write(file("dag_grid.svg"), grid_svg)
    grid_svg
  end

end
