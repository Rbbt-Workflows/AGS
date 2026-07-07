require 'json'
require 'set'

module AGS

  #{{{ GRN CONSTRUCTION

  dep :change_offsets, jobname: 'Default'
  dep :tf_predictions, jobname: 'Default'
  task :grns => :array do
    change_offsets = step(:change_offsets).load
    decoupler_jobs = rec_dependencies.select do |dep|
      dep.task_name == :timepoint_decoupler
    end

    self.progress_bar max: decoupler_jobs.length do |bar|
      decoupler_jobs.each do |job|
        treatment, time_point = job.recursive_inputs.values_at :treatment, :time_point

        tfs = {}
        tgs = {}
        edges = {}

        regulome = job.step(:filtered_regulome).path.tsv type: :double
        regulome = regulome.change_key 'source', merge: true
        regulome.unnamed = false

        job.load.each do |tf,activity|
          tfs[tf] = activity.first
          NamedArray.zip_fields(regulome[tf]).each do |tg,weight|
            changes = change_offsets[tg][treatment]
            direction = nil

            # Iterate in reverse and allow matching in previous sampled
            # timepoint if no more recent change is found. Sampled timepoints
            # are not equally spaced: previous to 4 h is 2 h, not 3 h.
            time_index = TIME_POINTS.index(time_point)
            previous_time_point = time_index && time_index > 0 ? TIME_POINTS[time_index - 1] : nil

            changes.reverse.each do |change|
              match = false
              dir,_sep, time = change.partition ' '
              time = time.sub 'h', ''
              if time.include? '-'
                start, eend = time.split '-'
                match = true if (start.to_i..eend.to_i).include? time_point
              else
                match = true if time.to_i == time_point
                match = true if previous_time_point && time.to_i == previous_time_point
              end

              if match
                direction = dir
                break
              end
            end
            tgs[tg] = direction
            edges[[tf, tg]*'~'] = weight.to_i
          end
        end


        grn = {tfs: tfs, tgs: tgs, edges: edges}

        target = [treatment, time_point] * '-'
        Open.write file(target + '.json'), grn.to_json
        bar.tick
      end
    end
    files
  end

  #{{{ GRN ANALYSIS HELPERS

  helper :grn_load_json do |treatment, time_point|
    filename = [treatment, time_point] * '-'
    json_file = step(:grns).file(filename + '.json')
    content = Open.read(json_file)
    data = JSON.parse(content)
    {
      :tfs   => data['tfs'] || {},
      :tgs   => data['tgs'] || {},
      :edges => data['edges'] || {}
    }
  end

  helper :grn_load_modules do
    Rbbt.data["grn_modules.tsv"].tsv type: :flat, fields: ["Gene"]
  end

  helper :grn_gene_modules do |modules|
    gene_modules = Hash.new { |h, k| h[k] = [] }
    modules.each do |mod, genes|
      genes.each do |gene|
        gene_modules[gene] << mod
      end
    end
    gene_modules
  end

  helper :grn_edge_pressure do |tf_activity, edge_sign|
    return nil if tf_activity.nil? || edge_sign.nil?
    tf_activity.to_f * edge_sign.to_i
  end

  helper :grn_coherent? do |pressure, tg_direction|
    return false if pressure.nil? || tg_direction.nil?
    if pressure > 0 && tg_direction == "increase"
      true
    elsif pressure < 0 && tg_direction == "decrease"
      true
    else
      false
    end
  end


  helper :grn_tg_state do |tgs, gene|
    return "absent" unless tgs.key?(gene)
    direction = tgs[gene]
    return "absent" if direction.nil?
    direction == "increase" ? "up" : "down"
  end

  helper :grn_transition_counts_for_genes do |tgs1, tgs2, genes|
    counts = Hash.new(0)
    genes.each do |gene|
      s1 = grn_tg_state(tgs1, gene)
      s2 = grn_tg_state(tgs2, gene)
      counts["#{s1}-to-#{s2}"] += 1
    end
    counts
  end

  helper :grn_transition_counts_string do |counts|
    counts.sort.collect { |cat, cnt| "#{cat}:#{cnt}" } * ","
  end

  helper :grn_bounded_pressure do |normalized_net_pressure|
    # Normalize pressure into approximately [-1, 1]. The raw pressure values
    # depend on edge counts and decoupleR score scale; tanh keeps module scores
    # comparable across modules while preserving sign and ordering.
    Math.tanh(normalized_net_pressure.to_f / 2.0)
  end

  helper :grn_module_state_score do |score, module_size|
    module_size = module_size.to_f
    module_size = 1.0 if module_size == 0
    up = score[:tg_increase].to_f
    down = score[:tg_decrease].to_f
    coverage = score[:coverage].to_f
    coherence_ratio = score[:coherence_ratio].to_f
    pressure = grn_bounded_pressure(score[:normalized_net_pressure].to_f)

    target_activation = (up - down) / module_size
    target_suppression = (down - up) / module_size

    activation_score = [target_activation, 0.0].max + [pressure, 0.0].max + 0.15 * coverage + 0.10 * coherence_ratio
    suppression_score = [target_suppression, 0.0].max + [-pressure, 0.0].max + 0.15 * coverage + 0.10 * coherence_ratio

    {
      :activation_score => activation_score.round(4),
      :suppression_score => suppression_score.round(4),
      :target_activation => target_activation.round(4),
      :target_suppression => target_suppression.round(4),
      :bounded_pressure => pressure.round(4)
    }
  end

  helper :grn_module_detail do |grn, mod_genes|
    tfs = grn[:tfs]
    tgs = grn[:tgs]
    edges = grn[:edges]

    target_edges = Hash.new { |h, k| h[k] = [] }
    edges.each do |edge_key, sign|
      tf, tg = edge_key.split("~")
      next unless tf && tg
      target_edges[tg] << [tf, sign.to_i, edge_key]
    end

    positive_targets = []
    negative_targets = []
    target_net_pressure = Hash.new(0.0)
    target_pos_pressure = Hash.new(0.0)
    target_neg_pressure = Hash.new(0.0)
    target_top_drivers = Hash.new { |h, k| h[k] = [] }
    tf_coherent_positive = Hash.new(0.0)
    tf_coherent_negative = Hash.new(0.0)
    tf_all_abs = Hash.new(0.0)
    coherent_positive_edges = 0
    coherent_negative_edges = 0
    conflict_edges = 0

    mod_genes.each do |gene|
      direction = tgs[gene]
      positive_targets << gene if direction == "increase"
      negative_targets << gene if direction == "decrease"
      next unless target_edges.key?(gene)

      edge_driver_values = []
      target_edges[gene].each do |tf, sign, _edge_key|
        tf_act = tfs[tf]
        next if tf_act.nil?
        pressure = tf_act.to_f * sign
        target_net_pressure[gene] += pressure
        if pressure > 0
          target_pos_pressure[gene] += pressure
        elsif pressure < 0
          target_neg_pressure[gene] += pressure.abs
        end
        tf_all_abs[tf] += pressure.abs
        edge_driver_values << [tf, pressure]

        if grn_coherent?(pressure, direction)
          if pressure > 0
            coherent_positive_edges += 1
            tf_coherent_positive[tf] += pressure
          elsif pressure < 0
            coherent_negative_edges += 1
            tf_coherent_negative[tf] += pressure.abs
          end
        elsif !direction.nil?
          conflict_edges += 1
        end
      end
      target_top_drivers[gene] = edge_driver_values.sort_by { |_, p| -p.abs }.first(5)
    end

    top_coherent_positive_tfs = tf_coherent_positive.sort_by { |_, v| -v }.first(10).collect { |tf, v| "#{tf}(#{v.round(2)})" } * ","
    top_coherent_negative_tfs = tf_coherent_negative.sort_by { |_, v| -v }.first(10).collect { |tf, v| "#{tf}(#{v.round(2)})" } * ","
    top_all_tfs = tf_all_abs.sort_by { |_, v| -v }.first(10).collect { |tf, v| "#{tf}(#{v.round(2)})" } * ","

    {
      :positive_targets => positive_targets,
      :negative_targets => negative_targets,
      :target_net_pressure => target_net_pressure,
      :target_pos_pressure => target_pos_pressure,
      :target_neg_pressure => target_neg_pressure,
      :target_top_drivers => target_top_drivers,
      :coherent_positive_edges => coherent_positive_edges,
      :coherent_negative_edges => coherent_negative_edges,
      :conflict_edges => conflict_edges,
      :top_coherent_positive_tfs => top_coherent_positive_tfs,
      :top_coherent_negative_tfs => top_coherent_negative_tfs,
      :top_all_tfs => top_all_tfs
    }
  end

  helper :grn_format_gene_list do |genes, max = 25|
    genes = genes.compact.uniq.sort
    suffix = genes.length > max ? "...(+#{genes.length - max})" : ""
    (genes.first(max) * ",") + suffix
  end

  helper :grn_format_target_drivers do |genes, details, max_genes = 12|
    genes.compact.uniq.sort.first(max_genes).collect do |gene|
      drivers = (details[:target_top_drivers][gene] || []).first(3).collect do |tf, pressure|
        "#{tf}:#{pressure.round(2)}"
      end * ";"
      drivers.empty? ? gene : "#{gene}[#{drivers}]"
    end * ","
  end

  # Compute module-level scores from a single GRN hash.
  # Returns a hash keyed by module name, each value is a hash of computed fields.
  helper :grn_compute_module_scores do |grn, modules|
    tfs = grn[:tfs]
    tgs = grn[:tgs]
    edges = grn[:edges]

    # Build target -> edges index
    target_edges = Hash.new { |h, k| h[k] = [] }
    edges.each do |edge_key, sign|
      tf, tg = edge_key.split("~")
      next unless tf && tg
      target_edges[tg] << [tf, sign.to_i, edge_key]
    end

    results = {}
    modules.each do |mod_name, mod_genes|
      # TG counts
      present_genes = mod_genes.select { |g| tgs.key?(g) }
      tg_increase = present_genes.count { |g| tgs[g] == "increase" }
      tg_decrease = present_genes.count { |g| tgs[g] == "decrease" }
      tg_absent = mod_genes.length - tg_increase - tg_decrease
      coverage = mod_genes.length > 0 ? present_genes.length.to_f / mod_genes.length : 0.0

      # Edge-based module scores. Pressure > 0 predicts target increase;
      # pressure < 0 predicts target decrease.
      positive_pressure = 0.0
      negative_pressure = 0.0
      coherence = 0
      conflict = 0
      edge_count = 0
      tf_abs_pressure = Hash.new(0.0)
      tf_positive_pressure = Hash.new(0.0)
      tf_negative_pressure = Hash.new(0.0)

      mod_genes.each do |gene|
        next unless target_edges.key?(gene)
        tg_dir = tgs[gene]
        target_edges[gene].each do |tf, sign, edge_key|
          tf_act = tfs[tf]
          next if tf_act.nil?
          pressure = tf_act.to_f * sign
          edge_count += 1
          if pressure > 0
            positive_pressure += pressure
            tf_positive_pressure[tf] += pressure
          elsif pressure < 0
            negative_pressure += pressure
            tf_negative_pressure[tf] += pressure.abs
          end
          tf_abs_pressure[tf] += pressure.abs
          if grn_coherent?(pressure, tg_dir)
            coherence += 1
          elsif !tg_dir.nil?
            conflict += 1
          end
        end
      end

      net_pressure = positive_pressure + negative_pressure
      normalized_net_pressure = edge_count > 0 ? net_pressure / edge_count : 0.0
      coherence_ratio = (coherence + conflict) > 0 ? coherence.to_f / (coherence + conflict) : 0.0

      top_tfs = tf_abs_pressure.sort_by { |_, v| -v }.first(10).collect { |tf, v| "#{tf}(#{v.round(2)})" } * ","
      top_positive_tfs = tf_positive_pressure.sort_by { |_, v| -v }.first(10).collect { |tf, v| "#{tf}(#{v.round(2)})" } * ","
      top_negative_tfs = tf_negative_pressure.sort_by { |_, v| -v }.first(10).collect { |tf, v| "#{tf}(#{v.round(2)})" } * ","

      results[mod_name] = {
        :tg_increase       => tg_increase,
        :tg_decrease       => tg_decrease,
        :tg_absent         => tg_absent,
        :coverage          => coverage.round(4),
        :edge_count        => edge_count,
        :positive_pressure => positive_pressure.round(4),
        :negative_pressure => negative_pressure.round(4),
        :net_pressure      => net_pressure.round(4),
        :normalized_net_pressure => normalized_net_pressure.round(4),
        :coherence         => coherence,
        :conflict          => conflict,
        :coherence_ratio   => coherence_ratio.round(4),
        :top_tfs           => top_tfs,
        :top_positive_tfs  => top_positive_tfs,
        :top_negative_tfs  => top_negative_tfs
      }
    end

    results
  end

  #{{{ GRN MODULE SCORES

  dep :grns, jobname: 'Default'
  input :treatment, :select, "Treatment", nil, :select_options => TREATMENTS, :required => true
  input :time_point, :select, "Timepoint", nil, :select_options => TIME_POINTS, :required => true
  task :grn_module_scores => :tsv do |treatment, time_point|
    grn = grn_load_json(treatment, time_point)
    modules = grn_load_modules
    scores = grn_compute_module_scores(grn, modules)

    fields = %w(TG_increase TG_decrease TG_absent Coverage Edge_count Positive_pressure Negative_pressure Net_pressure Normalized_net_pressure Coherence Conflict Coherence_ratio Top_TFs Top_positive_TFs Top_negative_TFs)
    result = TSV.setup({}, :key_field => "Module", :fields => fields, :type => :list)

    modules.keys.sort.each do |mod|
      s = scores[mod]
      result[mod] = [
        s[:tg_increase].to_s,
        s[:tg_decrease].to_s,
        s[:tg_absent].to_s,
        s[:coverage].to_s,
        s[:edge_count].to_s,
        s[:positive_pressure].to_s,
        s[:negative_pressure].to_s,
        s[:net_pressure].to_s,
        s[:normalized_net_pressure].to_s,
        s[:coherence].to_s,
        s[:conflict].to_s,
        s[:coherence_ratio].to_s,
        s[:top_tfs],
        s[:top_positive_tfs],
        s[:top_negative_tfs]
      ]
    end

    result
  end

  #{{{ GRN TRANSITION SCORES

  dep :grns, jobname: 'Default'
  input :treatment, :select, "Treatment", nil, :select_options => TREATMENTS, :required => true
  input :time_point, :select, "Earlier timepoint (T1)", 4, :select_options => TIME_POINTS, :required => true
  input :time_point2, :select, "Later timepoint (T2)", 8, :select_options => TIME_POINTS, :required => true
  task :grn_transition_scores => :tsv do |treatment, time_point, time_point2|
    grn1 = grn_load_json(treatment, time_point)
    grn2 = grn_load_json(treatment, time_point2)
    modules = grn_load_modules

    tgs1 = grn1[:tgs]
    tgs2 = grn2[:tgs]
    tfs1 = grn1[:tfs]
    tfs2 = grn2[:tfs]
    edges1 = grn1[:edges]
    edges2 = grn2[:edges]

    # Helper to classify a TG's state
    tg_state = lambda do |tgs, gene|
      return "absent" unless tgs.key?(gene)
      d = tgs[gene]
      return "absent" if d.nil?
      d == "increase" ? "up" : "down"
    end

    # Compute module-level scores at each timepoint
    scores1 = grn_compute_module_scores(grn1, modules)
    scores2 = grn_compute_module_scores(grn2, modules)

    # Build result TSV
    fields = %w(Transition_counts TG_total Net_pressure_T1 Net_pressure_T2 Pressure_delta Coherence_T1 Coherence_T2 Coherence_delta Top_TFs_T1 Top_TFs_T2)
    result = TSV.setup({}, :key_field => "Module", :fields => fields, :type => :list)

    all_tgs = Set.new(tgs1.keys) | Set.new(tgs2.keys)

    modules.keys.sort.each do |mod|
      mod_genes = modules[mod]
      mod_gene_set = Set.new(mod_genes)

      transition_counts = Hash.new(0)
      mod_genes.each do |gene|
        s1 = tg_state.call(tgs1, gene)
        s2 = tg_state.call(tgs2, gene)
        transition = "#{s1}-to-#{s2}"
        transition_counts[transition] += 1
      end

      # Format transition counts as "category:count,..."
      trans_str = transition_counts.sort.collect { |cat, cnt| "#{cat}:#{cnt}" } * ","

      s1 = scores1[mod]
      s2 = scores2[mod]

      net1 = s1[:net_pressure]
      net2 = s2[:net_pressure]
      pressure_delta = (net2 - net1).round(4)

      coh1 = s1[:coherence]
      coh2 = s2[:coherence]
      coh_delta = coh2 - coh1

      result[mod] = [
        trans_str,
        mod_genes.length.to_s,
        net1.to_s,
        net2.to_s,
        pressure_delta.to_s,
        coh1.to_s,
        coh2.to_s,
        coh_delta.to_s,
        s1[:top_tfs],
        s2[:top_tfs]
      ]
    end

    result
  end

  #{{{ GRN COMBINATION CONTRAST

  dep :grns, jobname: 'Default'
  input :combo_treatment, :select, "Combination treatment", "INT_PD_PI", :select_options => TREATMENTS
  input :component1, :select, "Component 1 treatment", "PI", :select_options => TREATMENTS
  input :component2, :select, "Component 2 treatment", "PD", :select_options => TREATMENTS
  input :time_point, :select, "Timepoint", 8, :select_options => TIME_POINTS
  task :grn_combination_contrast => :tsv do |combo_treatment, component1, component2, time_point|
    grn_combo = grn_load_json(combo_treatment, time_point)
    grn_c1 = grn_load_json(component1, time_point)
    grn_c2 = grn_load_json(component2, time_point)

    combo_tgs = grn_combo[:tgs]
    c1_tgs = grn_c1[:tgs]
    c2_tgs = grn_c2[:tgs]

    combo_tfs = grn_combo[:tfs]
    c1_tfs = grn_c1[:tfs]
    c2_tfs = grn_c2[:tfs]

    combo_edges = grn_combo[:edges]
    c1_edges = grn_c1[:edges]
    c2_edges = grn_c2[:edges]

    modules = grn_load_modules
    gene_modules = grn_gene_modules(modules)

    # Helper: classify a TG across combo vs two components
    classify_tg = lambda do |combo_dir, c1_dir, c2_dir|
      combo_present = !combo_dir.nil?
      c1_present = !c1_dir.nil?
      c2_present = !c2_dir.nil?

      # cooperative-latch: absent in both components but present in combo
      if combo_present && !c1_present && !c2_present
        return "cooperative-latch"
      end

      # combo-specific: present in combo, absent in one component, present in other but different direction
      if combo_present
        if c1_present && combo_dir == c1_dir && !c2_present
          return "component1-inherited"
        end
        if c2_present && combo_dir == c2_dir && !c1_present
          return "component2-inherited"
        end
        if c1_present && combo_dir == c1_dir && c2_present && combo_dir == c2_dir
          return "shared"
        end
        # sign-inverted relative to a component that was present
        if c1_present && combo_dir != c1_dir
          return "sign-inverted"
        end
        if c2_present && combo_dir != c2_dir
          return "sign-inverted"
        end
        # combo-specific: present in combo but neither component matches
        return "combo-specific"
      else
        # combo absent
        if (c1_present || c2_present)
          return "component-suppressed"
        else
          return "absent-in-all"
        end
      end
    end

    # Helper: classify a TF based on activity comparison
    classify_tf = lambda do |combo_act, c1_act, c2_act|
      combo_present = !combo_act.nil? && combo_act.to_f != 0
      c1_present = !c1_act.nil? && c1_act.to_f != 0
      c2_present = !c2_act.nil? && c2_act.to_f != 0

      combo_v = combo_act.nil? ? 0.0 : combo_act.to_f
      c1_v = c1_act.nil? ? 0.0 : c1_act.to_f
      c2_v = c2_act.nil? ? 0.0 : c2_act.to_f

      max_comp = [c1_v.abs, c2_v.abs].max

      if combo_present
        if c1_present && (combo_v <=> 0) == (c1_v <=> 0) && !c2_present
          return "component1-inherited"
        end
        if c2_present && (combo_v <=> 0) == (c2_v <=> 0) && !c1_present
          return "component2-inherited"
        end
        if c1_present && c2_present
          same_sign_c1 = (combo_v <=> 0) == (c1_v <=> 0)
          same_sign_c2 = (combo_v <=> 0) == (c2_v <=> 0)
          if same_sign_c1 && same_sign_c2
            if combo_v.abs > max_comp * 1.2
              return "stronger"
            elsif combo_v.abs < max_comp * 0.8
              return "weaker"
            else
              return "shared"
            end
          elsif !same_sign_c1 || !same_sign_c2
            return "sign-inverted"
          end
        end
        if !c1_present && !c2_present
          return "combo-specific"
        end
        return "combo-specific"
      else
        if c1_present || c2_present
          return "suppressed"
        else
          return "absent-in-all"
        end
      end
    end

    # Helper: classify an edge
    classify_edge = lambda do |combo_key, combo_sign, c1_sign, c2_sign|
      combo_present = !combo_sign.nil?
      c1_present = !c1_sign.nil?
      c2_present = !c2_sign.nil?

      if combo_present
        if c1_present && combo_sign == c1_sign && !c2_present
          return "#{component1}-inherited"
        end
        if c2_present && combo_sign == c2_sign && !c1_present
          return "#{component2}-inherited"
        end
        if c1_present && c2_present
          if combo_sign == c1_sign && combo_sign == c2_sign
            return "shared"
          end
          return "sign-inverted"
        end
        if c1_present && combo_sign != c1_sign
          return "sign-inverted"
        end
        if c2_present && combo_sign != c2_sign
          return "sign-inverted"
        end
        return "combo-specific"
      else
        if c1_present || c2_present
          return "component-suppressed"
        else
          return "absent-in-all"
        end
      end
    end

    # Helper: compute edge pressure for a single GRN
    edge_pressure = lambda do |edges, tfs, edge_key, sign|
      tf, _tg = edge_key.split("~")
      tf_act = tfs[tf]
      return nil if tf_act.nil?
      tf_act.to_f * sign.to_i
    end

    fields = %w(Element_type Modules Combo_state Component1_state Component2_state Classification Combo_pressure Component1_pressure Component2_pressure)
    result = TSV.setup({}, :key_field => "Element", :fields => fields, :type => :list)

    # --- Target genes ---
    all_tgs = Set.new(combo_tgs.keys) | Set.new(c1_tgs.keys) | Set.new(c2_tgs.keys)
    all_tgs.sort.each do |gene|
      combo_dir = combo_tgs[gene]
      c1_dir = c1_tgs[gene]
      c2_dir = c2_tgs[gene]
      classification = classify_tg.call(combo_dir, c1_dir, c2_dir)
      result["TG:#{gene}"] = [
        "TG",
        (gene_modules[gene] || []) * ",",
        combo_dir.nil? ? "absent" : combo_dir,
        c1_dir.nil? ? "absent" : c1_dir,
        c2_dir.nil? ? "absent" : c2_dir,
        classification,
        "NA",
        "NA",
        "NA"
      ]
    end

    # --- TFs ---
    all_tfs = Set.new(combo_tfs.keys) | Set.new(c1_tfs.keys) | Set.new(c2_tfs.keys)
    all_tfs.sort.each do |tf|
      combo_act = combo_tfs[tf]
      c1_act = c1_tfs[tf]
      c2_act = c2_tfs[tf]
      classification = classify_tf.call(combo_act, c1_act, c2_act)
      result["TF:#{tf}"] = [
        "TF",
        (gene_modules[tf] || []) * ",",
        combo_act.nil? ? "absent" : combo_act.to_s,
        c1_act.nil? ? "absent" : c1_act.to_s,
        c2_act.nil? ? "absent" : c2_act.to_s,
        classification,
        combo_act.nil? ? "NA" : combo_act.to_s,
        c1_act.nil? ? "NA" : c1_act.to_s,
        c2_act.nil? ? "NA" : c2_act.to_s
      ]
    end

    # --- Edges ---
    all_edges = Set.new(combo_edges.keys) | Set.new(c1_edges.keys) | Set.new(c2_edges.keys)
    all_edges.sort.each do |edge_key|
      combo_sign = combo_edges[edge_key]
      c1_sign = c1_edges[edge_key]
      c2_sign = c2_edges[edge_key]
      classification = classify_edge.call(edge_key, combo_sign, c1_sign, c2_sign)

      combo_p = combo_sign.nil? ? "NA" : edge_pressure.call(combo_edges, combo_tfs, edge_key, combo_sign).to_s
      c1_p = c1_sign.nil? ? "NA" : edge_pressure.call(c1_edges, c1_tfs, edge_key, c1_sign).to_s
      c2_p = c2_sign.nil? ? "NA" : edge_pressure.call(c2_edges, c2_tfs, edge_key, c2_sign).to_s

      _tf, edge_tg = edge_key.split("~")
      result["EDGE:#{edge_key}"] = [
        "Edge",
        (gene_modules[edge_tg] || []) * ",",
        combo_sign.nil? ? "absent" : combo_sign.to_s,
        c1_sign.nil? ? "absent" : c1_sign.to_s,
        c2_sign.nil? ? "absent" : c2_sign.to_s,
        classification,
        combo_p,
        c1_p,
        c2_p
      ]
    end

    result
  end


  #{{{ GRN SYNERGY WINDOW

  dep :grns, jobname: 'Default'
  input :combo_treatment, :select, "Combination treatment", "INT_PD_PI", :select_options => TREATMENTS, :required => true
  input :component1, :select, "Component 1 treatment", "PI", :select_options => TREATMENTS, :required => true
  input :component2, :select, "Component 2 treatment", "PD", :select_options => TREATMENTS, :required => true
  input :time_point, :select, "Earlier timepoint (default: T4, interval 2-4 h)", 4, :select_options => TIME_POINTS, :required => true
  input :time_point2, :select, "Later timepoint (default: T8, interval 4-8 h)", 8, :select_options => TIME_POINTS, :required => true
  task :grn_synergy_window => :tsv do |combo_treatment, component1, component2, time_point, time_point2|
    modules = grn_load_modules

    grn_combo_1 = grn_load_json(combo_treatment, time_point)
    grn_combo_2 = grn_load_json(combo_treatment, time_point2)
    grn_c1_1 = grn_load_json(component1, time_point)
    grn_c1_2 = grn_load_json(component1, time_point2)
    grn_c2_1 = grn_load_json(component2, time_point)
    grn_c2_2 = grn_load_json(component2, time_point2)

    combo_scores_1 = grn_compute_module_scores(grn_combo_1, modules)
    combo_scores_2 = grn_compute_module_scores(grn_combo_2, modules)
    c1_scores_1 = grn_compute_module_scores(grn_c1_1, modules)
    c1_scores_2 = grn_compute_module_scores(grn_c1_2, modules)
    c2_scores_1 = grn_compute_module_scores(grn_c2_1, modules)
    c2_scores_2 = grn_compute_module_scores(grn_c2_2, modules)

    arrest_modules = Set.new %w(checkpoint_arrest replication_licensing mitotic_execution myc_e2f_foxm1)
    escape_modules = Set.new %w(ap1_ets_srf plasticity_adhesion inflammatory_stat_nfkb isr_upr_proteostasis metabolism)

    fields = %w(
      Module_class Rank_score Event_call Interpretation_hint
      Combo_transition Component1_transition Component2_transition
      Combo_T1_activation Combo_T2_activation Component1_T2_activation Component2_T2_activation
      Combo_T1_suppression Combo_T2_suppression Component1_T2_suppression Component2_T2_suppression
      Combo_pressure_T1 Combo_pressure_T2 Component1_pressure_T2 Component2_pressure_T2
      Latch_advantage Escape_block Component_escape_mean Combo_escape
      Combo_coherent_negative_T2 Combo_coherent_positive_T2 Component1_coherent_positive_T2 Component2_coherent_positive_T2
      Combo_latch_targets Component_escape_targets Combo_specific_down_targets Component_released_targets
      Combo_top_negative_TFs_T2 Combo_top_positive_TFs_T2 Component1_top_positive_TFs_T2 Component2_top_positive_TFs_T2
      Combo_latch_target_drivers Component1_escape_target_drivers Component2_escape_target_drivers
    )
    result = TSV.setup({}, :key_field => "Module", :fields => fields, :type => :list)

    modules.keys.sort.each do |mod|
      genes = modules[mod]
      module_size = genes.length

      cs1 = combo_scores_1[mod]
      cs2 = combo_scores_2[mod]
      c1s1 = c1_scores_1[mod]
      c1s2 = c1_scores_2[mod]
      c2s1 = c2_scores_1[mod]
      c2s2 = c2_scores_2[mod]

      combo_state_1 = grn_module_state_score(cs1, module_size)
      combo_state_2 = grn_module_state_score(cs2, module_size)
      c1_state_1 = grn_module_state_score(c1s1, module_size)
      c1_state_2 = grn_module_state_score(c1s2, module_size)
      c2_state_1 = grn_module_state_score(c2s1, module_size)
      c2_state_2 = grn_module_state_score(c2s2, module_size)

      combo_detail_2 = grn_module_detail(grn_combo_2, genes)
      c1_detail_2 = grn_module_detail(grn_c1_2, genes)
      c2_detail_2 = grn_module_detail(grn_c2_2, genes)

      combo_trans = grn_transition_counts_for_genes(grn_combo_1[:tgs], grn_combo_2[:tgs], genes)
      c1_trans = grn_transition_counts_for_genes(grn_c1_1[:tgs], grn_c1_2[:tgs], genes)
      c2_trans = grn_transition_counts_for_genes(grn_c2_1[:tgs], grn_c2_2[:tgs], genes)

      comp_t2_activation_mean = (c1_state_2[:activation_score] + c2_state_2[:activation_score]) / 2.0
      comp_t2_suppression_mean = (c1_state_2[:suppression_score] + c2_state_2[:suppression_score]) / 2.0
      component_escape_mean = comp_t2_activation_mean
      combo_escape = combo_state_2[:activation_score]

      latch_advantage = combo_state_2[:suppression_score] - comp_t2_suppression_mean
      escape_block = component_escape_mean - combo_escape

      module_class = if arrest_modules.include?(mod)
                       "arrest_latch"
                     elsif escape_modules.include?(mod)
                       "escape_priming"
                     else
                       "other"
                     end

      # Rank score prioritizes different evidence depending on module class.
      # Arrest modules score highly if the combination maintains stronger T8
      # suppression than its components. Escape modules score highly if the
      # components show T8 activation that is blocked in the combination.
      rank_score = if module_class == "arrest_latch"
                     combo_state_2[:suppression_score] + [latch_advantage, 0.0].max + [escape_block, 0.0].max
                   elsif module_class == "escape_priming"
                     [escape_block, 0.0].max + component_escape_mean + [combo_state_2[:suppression_score] - comp_t2_suppression_mean, 0.0].max
                   else
                     [latch_advantage, 0.0].max + [escape_block, 0.0].max
                   end
      rank_score = rank_score.round(4)

      combo_down = Set.new(combo_detail_2[:negative_targets])
      combo_up = Set.new(combo_detail_2[:positive_targets])
      c1_down = Set.new(c1_detail_2[:negative_targets])
      c1_up = Set.new(c1_detail_2[:positive_targets])
      c2_down = Set.new(c2_detail_2[:negative_targets])
      c2_up = Set.new(c2_detail_2[:positive_targets])

      component_up = c1_up | c2_up
      component_down = c1_down | c2_down

      combo_latch_targets = combo_down.to_a
      component_escape_targets = (component_up - combo_up).to_a
      combo_specific_down_targets = (combo_down - component_down).to_a
      component_released_targets = genes.select do |gene|
        grn_tg_state(grn_combo_1[:tgs], gene) == "down" &&
          grn_tg_state(grn_combo_2[:tgs], gene) != "down" &&
          (grn_tg_state(grn_c1_2[:tgs], gene) == "up" || grn_tg_state(grn_c2_2[:tgs], gene) == "up")
      end

      event_call = if module_class == "arrest_latch" && latch_advantage > 0.25
                     "combination_latch"
                   elsif module_class == "escape_priming" && escape_block > 0.25
                     "blocked_escape_priming"
                   elsif combo_state_2[:activation_score] > component_escape_mean + 0.25
                     "combo_reopening"
                   elsif combo_state_2[:suppression_score] > 1.0
                     "maintained_suppression"
                   elsif combo_state_2[:activation_score] > 1.0
                     "reopening"
                   else
                     "mixed_or_weak"
                   end

      interpretation_hint = case event_call
                            when "combination_latch"
                              "#{combo_treatment} maintains stronger #{mod} suppression at #{time_point2}h than #{component1}/#{component2}; candidate arrest latch."
                            when "blocked_escape_priming"
                              "#{component1}/#{component2} show more #{mod} activation at #{time_point2}h than #{combo_treatment}; candidate blocked escape route."
                            when "combo_reopening"
                              "#{combo_treatment} reopens #{mod} more than its components; inspect for combination-specific stress or rebound."
                            when "maintained_suppression"
                              "#{combo_treatment} keeps #{mod} suppressed at #{time_point2}h, but not clearly more than components."
                            when "reopening"
                              "#{combo_treatment} shows #{mod} reopening at #{time_point2}h."
                            else
                              "No dominant differential window behavior detected for #{mod}."
                            end

      result[mod] = [
        module_class,
        rank_score.to_s,
        event_call,
        interpretation_hint,
        grn_transition_counts_string(combo_trans),
        grn_transition_counts_string(c1_trans),
        grn_transition_counts_string(c2_trans),
        combo_state_1[:activation_score].to_s,
        combo_state_2[:activation_score].to_s,
        c1_state_2[:activation_score].to_s,
        c2_state_2[:activation_score].to_s,
        combo_state_1[:suppression_score].to_s,
        combo_state_2[:suppression_score].to_s,
        c1_state_2[:suppression_score].to_s,
        c2_state_2[:suppression_score].to_s,
        cs1[:normalized_net_pressure].to_s,
        cs2[:normalized_net_pressure].to_s,
        c1s2[:normalized_net_pressure].to_s,
        c2s2[:normalized_net_pressure].to_s,
        latch_advantage.round(4).to_s,
        escape_block.round(4).to_s,
        component_escape_mean.round(4).to_s,
        combo_escape.round(4).to_s,
        combo_detail_2[:coherent_negative_edges].to_s,
        combo_detail_2[:coherent_positive_edges].to_s,
        c1_detail_2[:coherent_positive_edges].to_s,
        c2_detail_2[:coherent_positive_edges].to_s,
        grn_format_gene_list(combo_latch_targets),
        grn_format_gene_list(component_escape_targets),
        grn_format_gene_list(combo_specific_down_targets),
        grn_format_gene_list(component_released_targets),
        combo_detail_2[:top_coherent_negative_tfs],
        combo_detail_2[:top_coherent_positive_tfs],
        c1_detail_2[:top_coherent_positive_tfs],
        c2_detail_2[:top_coherent_positive_tfs],
        grn_format_target_drivers(combo_latch_targets, combo_detail_2),
        grn_format_target_drivers(component_escape_targets, c1_detail_2),
        grn_format_target_drivers(component_escape_targets, c2_detail_2)
      ]
    end

    # Return modules in descending rank-score order while preserving TSV shape.
    ordered = TSV.setup({}, :key_field => "Module", :fields => fields, :type => :list)
    result.sort_by { |_mod, values| -values[1].to_f }.each do |mod, values|
      ordered[mod] = values
    end
    ordered
  end


  #{{{ GRN DRIVER SET

  dep :grns, jobname: 'Default'
  input :treatment, :select, "Treatment", nil, :select_options => TREATMENTS, :required => true
  input :time_point, :select, "Timepoint", nil, :select_options => TIME_POINTS, :required => true
  input :module_filter, :string, "Optional module name, or 'all'", "all"
  input :max_drivers, :integer, "Maximum drivers to report per module and direction", 8
  input :coverage_threshold, :float, "Stop once this fraction of dynamic targets is covered", 0.7
  task :grn_driver_set => :tsv do |treatment, time_point, module_filter, max_drivers, coverage_threshold|
    grn = grn_load_json(treatment, time_point)
    modules = grn_load_modules
    module_filter = module_filter.to_s
    max_drivers = max_drivers.to_i
    max_drivers = 8 if max_drivers <= 0
    coverage_threshold = coverage_threshold.to_f
    coverage_threshold = 0.7 if coverage_threshold <= 0

    selected_modules = if module_filter.nil? || module_filter.empty? || module_filter == "all"
                         modules.keys.sort
                       else
                         module_filter.split(/[,;\s]+/).select { |m| modules.key?(m) }
                       end

    tfs = grn[:tfs]
    tgs = grn[:tgs]
    edges = grn[:edges]

    # Build TF -> edges index to avoid repeatedly scanning all edges for each module.
    tf_edges = Hash.new { |h, k| h[k] = [] }
    edges.each do |edge_key, sign|
      tf, tg = edge_key.split("~")
      next unless tf && tg
      next unless tfs.key?(tf)
      next unless tgs.key?(tg)
      pressure = tfs[tf].to_f * sign.to_i
      next if pressure == 0
      coherent = grn_coherent?(pressure, tgs[tg])
      tf_edges[tf] << [tg, sign.to_i, pressure, coherent]
    end

    fields = %w(
      Direction Step Driver_TF TF_activity Marginal_targets Cumulative_targets Total_dynamic_targets
      Marginal_coverage Cumulative_coverage Marginal_abs_pressure Marginal_net_pressure
      Coherent_positive_edges Coherent_negative_edges Targets Top_target_pressures
    )
    result = TSV.setup({}, :key_field => "Module", :fields => fields, :type => :list)

    selected_modules.each do |mod|
      genes = modules[mod]
      gene_set = Set.new(genes)
      dynamic_by_direction = {
        "increase" => genes.select { |g| tgs[g] == "increase" },
        "decrease" => genes.select { |g| tgs[g] == "decrease" }
      }

      %w(increase decrease).each do |direction|
        dynamic_targets = dynamic_by_direction[direction]
        next if dynamic_targets.empty?
        target_universe = Set.new(dynamic_targets)

        candidates = {}
        tf_edges.each do |tf, edge_list|
          coherent_targets = {}
          positive_edges = 0
          negative_edges = 0
          edge_list.each do |tg, _sign, pressure, coherent|
            next unless coherent
            next unless target_universe.include?(tg)
            next unless (direction == "increase" && pressure > 0) || (direction == "decrease" && pressure < 0)
            coherent_targets[tg] ||= 0.0
            coherent_targets[tg] += pressure
            positive_edges += 1 if pressure > 0
            negative_edges += 1 if pressure < 0
          end
          next if coherent_targets.empty?
          candidates[tf] = {
            :targets => coherent_targets,
            :positive_edges => positive_edges,
            :negative_edges => negative_edges,
            :abs_pressure => coherent_targets.values.collect { |v| v.abs }.inject(0.0, :+),
            :net_pressure => coherent_targets.values.inject(0.0, :+)
          }
        end

        covered = Set.new
        step = 0
        while step < max_drivers && covered.length.to_f / target_universe.length < coverage_threshold
          best_tf = nil
          best_data = nil
          best_uncovered = []
          candidates.each do |tf, data|
            uncovered_targets = data[:targets].keys.reject { |tg| covered.include?(tg) }
            next if uncovered_targets.empty?
            uncovered_pressure = uncovered_targets.collect { |tg| data[:targets][tg].abs }.inject(0.0, :+)
            if best_tf.nil? || uncovered_targets.length > best_uncovered.length || (uncovered_targets.length == best_uncovered.length && uncovered_pressure > best_data[:_uncovered_pressure].to_f)
              best_tf = tf
              best_data = data.merge(:_uncovered_pressure => uncovered_pressure)
              best_uncovered = uncovered_targets
            end
          end
          break if best_tf.nil?

          step += 1
          best_uncovered.each { |tg| covered << tg }
          marginal_targets = best_uncovered.sort
          marginal_abs_pressure = marginal_targets.collect { |tg| best_data[:targets][tg].abs }.inject(0.0, :+)
          marginal_net_pressure = marginal_targets.collect { |tg| best_data[:targets][tg] }.inject(0.0, :+)
          cumulative_coverage = covered.length.to_f / target_universe.length
          marginal_coverage = marginal_targets.length.to_f / target_universe.length
          target_pressures = marginal_targets.collect { |tg| "#{tg}:#{best_data[:targets][tg].round(2)}" } * ","

          key = [mod, direction, step, best_tf] * "~"
          result[key] = [
            direction,
            step.to_s,
            best_tf,
            tfs[best_tf].to_s,
            marginal_targets.length.to_s,
            covered.length.to_s,
            target_universe.length.to_s,
            marginal_coverage.round(4).to_s,
            cumulative_coverage.round(4).to_s,
            marginal_abs_pressure.round(4).to_s,
            marginal_net_pressure.round(4).to_s,
            best_data[:positive_edges].to_s,
            best_data[:negative_edges].to_s,
            grn_format_gene_list(marginal_targets, 40),
            target_pressures
          ]
        end
      end
    end

    result
  end


  #{{{ GRN MOTIFS

  dep :grns, jobname: 'Default'
  input :treatment, :select, "Treatment", nil, :select_options => TREATMENTS, :required => true
  input :time_point, :select, "Timepoint", nil, :select_options => TIME_POINTS, :required => true
  input :module_filter, :string, "Optional target module name, or 'all'", "all"
  input :max_motifs, :integer, "Maximum motifs to return", 500
  task :grn_motifs => :tsv do |treatment, time_point, module_filter, max_motifs|
    grn = grn_load_json(treatment, time_point)
    modules = grn_load_modules
    gene_modules = grn_gene_modules(modules)
    module_filter = module_filter.to_s
    max_motifs = max_motifs.to_i
    max_motifs = 500 if max_motifs <= 0

    selected_modules = if module_filter.nil? || module_filter.empty? || module_filter == "all"
                         modules.keys.sort
                       else
                         module_filter.split(/[,;\s]+/).select { |m| modules.key?(m) }
                       end
    selected_target_genes = Set.new
    selected_modules.each { |mod| modules[mod].each { |g| selected_target_genes << g } }

    tfs = grn[:tfs]
    tgs = grn[:tgs]
    edges = grn[:edges]

    # Index outgoing edges for active TFs.
    outgoing = Hash.new { |h, k| h[k] = {} }
    edges.each do |edge_key, sign|
      tf, tg = edge_key.split("~")
      next unless tf && tg
      next unless tfs.key?(tf)
      outgoing[tf][tg] = sign.to_i
    end

    fields = %w(
      Motif_class Target_modules TF_A TF_B Target A_activity B_activity
      A_to_B_sign A_to_Target_sign B_to_Target_sign Target_direction B_direction
      A_to_B_pressure A_to_Target_pressure B_to_Target_pressure
      Target_coherence B_transcript_coherence Motif_score Interpretation_hint
    )
    result = TSV.setup({}, :key_field => "Motif", :fields => fields, :type => :list)

    motifs = []
    outgoing.each do |tf_a, targets_a|
      a_activity = tfs[tf_a].to_f
      targets_a.each do |tf_b, sign_ab|
        next if tf_b == tf_a
        next unless tfs.key?(tf_b)
        next unless outgoing.key?(tf_b)
        shared_targets = targets_a.keys & outgoing[tf_b].keys
        shared_targets.each do |target|
          next if target == tf_a || target == tf_b
          next unless selected_target_genes.include?(target)
          next unless tgs.key?(target)
          sign_at = targets_a[target]
          sign_bt = outgoing[tf_b][target]
          b_activity = tfs[tf_b].to_f
          p_ab = a_activity * sign_ab
          p_at = a_activity * sign_at
          p_bt = b_activity * sign_bt
          target_direction = tgs[target]
          b_direction = tgs[tf_b]

          target_coherence_count = 0
          target_coherence_count += 1 if grn_coherent?(p_at, target_direction)
          target_coherence_count += 1 if grn_coherent?(p_bt, target_direction)
          target_coherence = case target_coherence_count
                             when 2 then "both_coherent"
                             when 1 then "one_coherent"
                             else "conflict"
                             end
          b_transcript_coherence = if b_direction.nil?
                                     "b_not_dynamic_tg"
                                   elsif grn_coherent?(p_ab, b_direction)
                                     "coherent"
                                   else
                                     "conflict"
                                   end

          motif_class = if p_at > 0 && p_bt > 0
                          "coherent_activation_ffl"
                        elsif p_at < 0 && p_bt < 0
                          "coherent_repression_ffl"
                        else
                          "incoherent_ffl"
                        end

          motif_score = p_at.abs + p_bt.abs + (target_coherence_count * 2.0)
          motif_score += p_ab.abs * 0.25
          motif_score += 1.0 if b_transcript_coherence == "coherent"
          hint = if motif_class == "coherent_repression_ffl" && target_direction == "decrease"
                   "Feed-forward repression coherently explains target decrease."
                 elsif motif_class == "coherent_activation_ffl" && target_direction == "increase"
                   "Feed-forward activation coherently explains target increase."
                 elsif motif_class == "incoherent_ffl"
                   "Incoherent feed-forward motif may buffer or create delayed reversal."
                 else
                   "Feed-forward motif has mixed relation to observed target direction."
                 end

          motifs << [motif_score, tf_a, tf_b, target, [
            motif_class,
            (gene_modules[target] || []) * ",",
            tf_a,
            tf_b,
            target,
            a_activity.round(4).to_s,
            b_activity.round(4).to_s,
            sign_ab.to_s,
            sign_at.to_s,
            sign_bt.to_s,
            target_direction.to_s,
            b_direction.nil? ? "absent" : b_direction.to_s,
            p_ab.round(4).to_s,
            p_at.round(4).to_s,
            p_bt.round(4).to_s,
            target_coherence,
            b_transcript_coherence,
            motif_score.round(4).to_s,
            hint
          ]]
        end
      end
    end

    motifs.sort_by { |score, _a, _b, _t, _values| -score }.first(max_motifs).each_with_index do |(_score, tf_a, tf_b, target, values), i|
      result[[i + 1, tf_a, tf_b, target] * "~"] = values
    end

    result
  end

  helper :grn_treatment_targets do |treatment|
    target_map = {
      "FiveZ" => %w(MAP3K7),
      "PD" => %w(MAP2K1 MAP2K2),
      "PI" => %w(PIK3CA PIK3CB PIK3CG PIK3R1),
      "INT_PD_PI" => %w(PIK3CA PIK3CB PIK3CG PIK3R1 MAP2K1 MAP2K2),
      "INT_FiveZ_PI" => %w(PIK3CA PIK3CB PIK3CG PIK3R1 MAP3K7),
      "DMSO" => []
    }
    target_map[treatment.to_s] || []
  end

  helper :grn_signor_effect_sign do |effect|
    effect = effect.to_s.downcase
    return -1 if effect.include?("down-regulates") || effect.include?("decreases") || effect.include?("inhibits")
    return 1 if effect.include?("up-regulates") || effect.include?("increases") || effect.include?("activates")
    nil
  end

  helper :grn_signor_edge_info do |index, pair|
    raw = index[pair]
    info = {}
    if raw.respond_to?(:fields) && raw.respond_to?(:[])
      raw.fields.each do |field|
        value = raw[field]
        value = value.first if Array === value && value.length == 1
        info[field] = value
      end
    elsif Hash === raw
      info = raw
    end
    info
  end

  helper :grn_signor_out_edges do |index, source|
    matches = index.match(source)
    out = []
    matches.each do |pair|
      src, tgt = pair.split("~")
      next unless src == source && tgt && !tgt.empty?
      info = grn_signor_edge_info(index, pair)
      effect = info["Effect"] || info[:Effect]
      sign = grn_signor_effect_sign(effect)
      next if sign.nil?
      mechanism = info["Mechanism"] || info[:Mechanism]
      residue = info["Residue"] || info[:Residue]
      out << {
        :pair => pair,
        :source => src,
        :target => tgt,
        :sign => sign,
        :effect => effect.to_s,
        :mechanism => mechanism.to_s,
        :residue => residue.to_s
      }
    end
    out
  end

  helper :grn_selected_module_genes do |modules, module_filter|
    module_filter = module_filter.to_s
    selected_modules = if module_filter.nil? || module_filter.empty? || module_filter == "all"
                         modules.keys.sort
                       else
                         module_filter.split(/[,;\s]+/).select { |m| modules.key?(m) }
                       end
    selected_genes = Set.new
    selected_modules.each { |mod| modules[mod].each { |g| selected_genes << g } }
    [selected_modules, selected_genes]
  end

  helper :grn_module_driver_tfs do |grn, selected_genes|
    tfs = grn[:tfs]
    driver_tfs = Set.new
    grn[:edges].each do |edge_key, _sign|
      tf, tg = edge_key.split("~")
      next unless tf && tg
      next unless selected_genes.include?(tg)
      next unless tfs.key?(tf)
      driver_tfs << tf
    end
    driver_tfs
  end


  #{{{ GRN SIGNOR BRIDGES

  dep :grns, jobname: 'Default'
  input :treatment, :select, "Treatment", nil, :select_options => TREATMENTS, :required => true
  input :time_point, :select, "Timepoint", nil, :select_options => TIME_POINTS, :required => true
  input :module_filter, :string, "Optional module name, or 'all'", "all"
  input :max_len, :integer, "Maximum SIGNOR path length", 3
  input :max_paths, :integer, "Maximum bridge paths to return", 1000
  task :grn_signor_bridges => :tsv do |treatment, time_point, module_filter, max_len, max_paths|
    grn = grn_load_json(treatment, time_point)
    modules = grn_load_modules
    gene_modules = grn_gene_modules(modules)
    selected_modules, selected_genes = grn_selected_module_genes(modules, module_filter)
    active_tfs = grn[:tfs]
    target_driver_tfs = module_filter.to_s == "all" ? Set.new(active_tfs.keys) : grn_module_driver_tfs(grn, selected_genes)

    max_len = max_len.to_i
    max_len = 3 if max_len <= 0
    max_paths = max_paths.to_i
    max_paths = 1000 if max_paths <= 0

    drug_targets = grn_treatment_targets(treatment)
    kb = gene_kb
    signor = kb.get_index 'Signor'

    fields = %w(
      Drug_target TF TF_activity TF_modules Selected_modules Path_length
      Target_activation_sign Inhibitor_expected_sign Observed_TF_sign Agreement
      Path Effects Mechanisms Residues Score Interpretation_hint
    )
    result = TSV.setup({}, :key_field => "Bridge", :fields => fields, :type => :list)

    bridges = []
    drug_targets.each do |drug_target|
      # Each path stores nodes, cumulative target-activation sign, and edge annotations.
      queue = [[drug_target, [drug_target], 1, [], [], []]]
      until queue.empty?
        current, nodes, cumulative_sign, effects, mechanisms, residues = queue.shift
        depth = nodes.length - 1
        next if depth >= max_len

        grn_signor_out_edges(signor, current).each do |edge|
          nxt = edge[:target]
          next if nodes.include?(nxt)
          new_nodes = nodes + [nxt]
          new_sign = cumulative_sign * edge[:sign]
          new_effects = effects + [edge[:effect]]
          new_mechanisms = mechanisms + [edge[:mechanism]]
          new_residues = residues + [edge[:residue]]
          new_depth = new_nodes.length - 1

          if target_driver_tfs.include?(nxt)
            tf_activity = active_tfs[nxt].to_f
            observed_sign = tf_activity > 0 ? 1 : -1
            inhibitor_expected = -1 * new_sign
            agreement = observed_sign == inhibitor_expected ? "match" : "mismatch"
            # Shorter paths, matching signs and larger TF activity rank higher.
            score = tf_activity.abs + (agreement == "match" ? 5.0 : 0.0) + (max_len - new_depth + 1)
            hint = if agreement == "match"
                     "Inhibition of #{drug_target} predicts #{nxt} #{inhibitor_expected > 0 ? 'activation' : 'repression'}, matching GRN TF activity."
                   else
                     "SIGNOR path predicts #{nxt} #{inhibitor_expected > 0 ? 'activation' : 'repression'}, but GRN TF activity has opposite sign."
                   end
            bridges << [score, drug_target, nxt, [
              drug_target,
              nxt,
              tf_activity.round(4).to_s,
              (gene_modules[nxt] || []) * ",",
              selected_modules * ",",
              new_depth.to_s,
              new_sign.to_s,
              inhibitor_expected.to_s,
              observed_sign.to_s,
              agreement,
              new_nodes * "->",
              new_effects * "|",
              new_mechanisms * "|",
              new_residues * "|",
              score.round(4).to_s,
              hint
            ]]
          end

          queue << [nxt, new_nodes, new_sign, new_effects, new_mechanisms, new_residues] if new_depth < max_len
        end
      end
    end

    bridges.sort_by { |score, _drug_target, _tf, _values| -score }.first(max_paths).each_with_index do |(_score, drug_target, tf, values), i|
      result[[i + 1, drug_target, tf] * "~"] = values
    end

    result
  end


  #{{{ GRN CORE STORY

  dep :grns, jobname: 'Default'
  input :latch_treatment, :select, "Latch / durable-suppression treatment", "INT_PD_PI", :select_options => TREATMENTS, :required => true
  input :escape_treatment, :select, "Escape / rebound-prone treatment", "INT_FiveZ_PI", :select_options => TREATMENTS, :required => true
  input :time_point, :select, "Timepoint to summarize", 8, :select_options => TIME_POINTS, :required => true
  input :max_edges_per_module, :integer, "Maximum coherent GRN edges per module/context", 12
  input :max_signor_paths, :integer, "Maximum SIGNOR bridge paths to include", 30
  task :grn_core_story => :string do |latch_treatment, escape_treatment, time_point, max_edges_per_module, max_signor_paths|
    max_edges_per_module = max_edges_per_module.to_i
    max_edges_per_module = 12 if max_edges_per_module <= 0
    max_signor_paths = max_signor_paths.to_i
    max_signor_paths = 30 if max_signor_paths <= 0

    modules = grn_load_modules
    gene_modules = grn_gene_modules(modules)
    latch_modules = %w(replication_licensing mitotic_execution myc_e2f_foxm1 checkpoint_arrest)
    escape_modules = %w(ap1_ets_srf plasticity_adhesion inflammatory_stat_nfkb isr_upr_proteostasis)

    contexts = {
      "latch" => {
        :treatment => latch_treatment,
        :grn => grn_load_json(latch_treatment, time_point),
        :modules => latch_modules
      },
      "escape" => {
        :treatment => escape_treatment,
        :grn => grn_load_json(escape_treatment, time_point),
        :modules => escape_modules
      }
    }

    nodes = {}
    grn_edges = []

    add_node = lambda do |gene, role, context_name, grn|
      nodes[gene] ||= {
        :id => gene,
        :roles => [],
        :modules => gene_modules[gene] || [],
        :tf_activity => {},
        :tg_direction => {},
        :contexts => []
      }
      nodes[gene][:roles] << role unless nodes[gene][:roles].include?(role)
      nodes[gene][:contexts] << context_name unless nodes[gene][:contexts].include?(context_name)
      nodes[gene][:tf_activity][context_name] = grn[:tfs][gene].to_f.round(4) if grn[:tfs].key?(gene)
      nodes[gene][:tg_direction][context_name] = grn[:tgs][gene] if grn[:tgs].key?(gene)
    end

    contexts.each do |context_name, info|
      grn = info[:grn]
      info[:modules].each do |mod|
        next unless modules.key?(mod)
        mod_genes = Set.new(modules[mod])
        candidates = []
        grn[:edges].each do |edge_key, sign|
          tf, tg = edge_key.split("~")
          next unless tf && tg
          next unless mod_genes.include?(tg)
          next unless grn[:tfs].key?(tf)
          next unless grn[:tgs].key?(tg)
          pressure = grn[:tfs][tf].to_f * sign.to_i
          next unless grn_coherent?(pressure, grn[:tgs][tg])
          score = pressure.abs
          # Emphasize the two intended story axes.
          score += 3.0 if context_name == "latch" && %w(replication_licensing mitotic_execution).include?(mod)
          score += 3.0 if context_name == "escape" && %w(ap1_ets_srf plasticity_adhesion).include?(mod)
          score += 1.0 if gene_modules[tf] && !(gene_modules[tf] & info[:modules]).empty?
          candidates << [score, tf, tg, sign.to_i, pressure, grn[:tgs][tg], mod]
        end

        candidates.sort_by { |score, _tf, _tg, _sign, _pressure, _dir, _mod| -score }.first(max_edges_per_module).each do |score, tf, tg, sign, pressure, direction, edge_module|
          add_node.call(tf, "TF", context_name, grn)
          add_node.call(tg, "TG", context_name, grn)
          grn_edges << {
            :type => "GRN",
            :context => context_name,
            :treatment => info[:treatment],
            :time_point => time_point,
            :source => tf,
            :target => tg,
            :module => edge_module,
            :sign => sign,
            :tf_activity => grn[:tfs][tf].to_f.round(4),
            :target_direction => direction,
            :pressure => pressure.round(4),
            :score => score.round(4)
          }
        end
      end
    end

    # Add drug target nodes.
    contexts.each do |context_name, info|
      grn_treatment_targets(info[:treatment]).each do |drug_target|
        nodes[drug_target] ||= {
          :id => drug_target,
          :roles => [],
          :modules => gene_modules[drug_target] || [],
          :tf_activity => {},
          :tg_direction => {},
          :contexts => []
        }
        nodes[drug_target][:roles] << "drug_target" unless nodes[drug_target][:roles].include?("drug_target")
        nodes[drug_target][:contexts] << context_name unless nodes[drug_target][:contexts].include?(context_name)
      end
    end

    # SIGNOR bridge paths from drug targets to selected TF nodes.
    selected_tfs = Set.new(nodes.values.select { |n| n[:roles].include?("TF") }.collect { |n| n[:id] })
    signor_edges = []
    begin
      signor = gene_kb.get_index 'Signor'
      contexts.each do |context_name, info|
        drug_targets = grn_treatment_targets(info[:treatment])
        grn = info[:grn]
        bridges = []
        drug_targets.each do |drug_target|
          queue = [[drug_target, [drug_target], 1, [], [], []]]
          until queue.empty?
            current, path_nodes, cumulative_sign, effects, mechanisms, residues = queue.shift
            depth = path_nodes.length - 1
            next if depth >= 3
            grn_signor_out_edges(signor, current).each do |edge|
              nxt = edge[:target]
              next if path_nodes.include?(nxt)
              new_nodes = path_nodes + [nxt]
              new_sign = cumulative_sign * edge[:sign]
              new_effects = effects + [edge[:effect]]
              new_mechanisms = mechanisms + [edge[:mechanism]]
              new_residues = residues + [edge[:residue]]
              new_depth = new_nodes.length - 1
              if selected_tfs.include?(nxt) && grn[:tfs].key?(nxt)
                tf_activity = grn[:tfs][nxt].to_f
                observed_sign = tf_activity > 0 ? 1 : -1
                inhibitor_expected = -1 * new_sign
                agreement = observed_sign == inhibitor_expected ? "match" : "mismatch"
                score = tf_activity.abs + (agreement == "match" ? 5.0 : 0.0) + (4 - new_depth)
                bridges << [score, {
                  :type => "SIGNOR_bridge",
                  :context => context_name,
                  :treatment => info[:treatment],
                  :drug_target => drug_target,
                  :target_tf => nxt,
                  :path => new_nodes,
                  :path_length => new_depth,
                  :target_activation_sign => new_sign,
                  :inhibitor_expected_sign => inhibitor_expected,
                  :observed_tf_sign => observed_sign,
                  :agreement => agreement,
                  :tf_activity => tf_activity.round(4),
                  :effects => new_effects,
                  :mechanisms => new_mechanisms,
                  :residues => new_residues,
                  :score => score.round(4)
                }]
              end
              queue << [nxt, new_nodes, new_sign, new_effects, new_mechanisms, new_residues] if new_depth < 3
            end
          end
        end
        signor_edges.concat bridges.sort_by { |score, _edge| -score }.first(max_signor_paths).collect { |_score, edge| edge }
      end
    rescue Exception => e
      signor_edges << {:type => "SIGNOR_bridge_error", :message => e.message}
    end

    # A compact summary for immediate inspection.
    summary = {
      :time_point => time_point,
      :latch_treatment => latch_treatment,
      :escape_treatment => escape_treatment,
      :latch_modules => latch_modules,
      :escape_modules => escape_modules,
      :node_count => nodes.length,
      :grn_edge_count => grn_edges.length,
      :signor_bridge_count => signor_edges.length,
      :interpretation => "Core story graph contrasts a durable INT_PD_PI-like replication/mitotic latch with an INT_FiveZ_PI-like AP1/plasticity escape-priming state."
    }

    JSON.pretty_generate({
      :summary => summary,
      :nodes => nodes.values,
      :edges => grn_edges + signor_edges
    })
  end


  helper :gene_kb do
    KnowledgeBase.new Scout.var.Agent.Gene.knowledge_base
  end
end
