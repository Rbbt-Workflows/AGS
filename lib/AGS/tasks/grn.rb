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


  helper :gene_kb do
    KnowledgeBase.new Scout.var.Agent.Gene.knowledge_base
  end

end

require_relative 'grn/prior'
require_relative 'grn/agnostic'
