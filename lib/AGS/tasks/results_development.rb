module AGS

  helper :scalar_value do |value|
    if Array === value
      value = value.compact.first
    end
    value = nil if value.respond_to?(:empty?) && value.empty?
    value
  end

  helper :numeric_value do |value, default=nil|
    value = scalar_value(value)
    value.nil? ? default : value.to_f
  end

  helper :context_field do |prefix, treatment, time_point|
    "#{prefix}_#{treatment}.T#{time_point}"
  end

  helper :parse_treatment_time_field do |field|
    if field =~ /^(.+)-T(\d+)$/
      [$1, $2.to_i]
    elsif field =~ /^(?:PvalueFC1|FC1|Pvalue|FC)_(.+)\.T(\d+)$/
      [$1, $2.to_i]
    else
      nil
    end
  end

  helper :onset_events do |label|
    values = [label].flatten.compact.collect(&:to_s).reject{|v| v.empty? }
    values = values.collect{|v| v.split('|') }.flatten
    values.reject!{|v| v.nil? || v.empty? || v == 'unclassified' }
    values.collect do |entry|
      if entry =~ /^(increase|decrease)\s+(\d+)(?:-(\d+))?h$/
        direction = $1 == 'increase' ? 'up' : 'down'
        start_time = $2.to_i
        end_time = ($3 || $2).to_i
        {
          :label => entry,
          :direction => direction,
          :start_time => start_time,
          :end_time => end_time
        }
      end
    end.compact
  end

  helper :onset_direction_at_time do |label, time_point|
    events = onset_events(label).select{|event| event[:start_time] == time_point.to_i }
    dirs = events.collect{|event| event[:direction] }.uniq
    if dirs.length == 1
      dirs.first
    elsif dirs.length > 1
      'both'
    else
      nil
    end
  end

  helper :fc1_value_from_fc0_values do |values, time_point|
    time_point = time_point.to_i
    index = AGS::TIME_POINTS.index(time_point)
    current = values[index]
    current = current.to_f unless current.nil?
    if index == 0
      current
    else
      previous = values[index - 1]
      previous = previous.to_f unless previous.nil?
      current.nil? || previous.nil? ? nil : current - previous
    end
  end

  helper :fc1_pvalue_surrogate_from_values do |values, time_point|
    time_point = time_point.to_i
    index = AGS::TIME_POINTS.index(time_point)
    current = values[index]
    current = current.nil? ? 1.0 : current.to_f
    if index == 0
      current
    else
      previous = values[index - 1]
      previous = previous.nil? ? 1.0 : previous.to_f
      [current, previous].max
    end
  end

  helper :activity_sign_label do |value|
    value = value.to_f
    value > 0 ? 'positive' : (value < 0 ? 'negative' : 'zero')
  end

  helper :activity_direction_label do |value|
    value = value.to_f
    value > 0 ? 'up' : (value < 0 ? 'down' : nil)
  end

  helper :boolean_value do |value|
    value = scalar_value(value)
    value == true || value.to_s == 'true'
  end

  helper :parse_go_obo_terms do |text|
    terms = {}
    current = nil
    in_term = false

    finalize = lambda do
      if current && current[:id] && ! current[:obsolete]
        terms[current[:id]] = current
      end
    end

    text.each_line do |line|
      line = line.chomp
      case line
      when '[Term]'
        finalize.call
        current = {:parents => []}
        in_term = true
      when /^\[/
        finalize.call if in_term
        current = nil
        in_term = false
      else
        next unless in_term && current
        case line
        when /^id:\s+(GO:\d+)/
          current[:id] = $1
        when /^name:\s+(.+)/
          current[:name] = $1.strip
        when /^namespace:\s+(.+)/
          current[:namespace] = $1.strip
        when /^is_a:\s+(GO:\d+)/
          current[:parents] << $1
        when /^relationship:\s+part_of\s+(GO:\d+)/
          current[:parents] << $1
        when /^is_obsolete:\s+true/
          current[:obsolete] = true
        end
      end
    end
    finalize.call
    terms
  end

  helper :go_ancestors_for_term do |term_id, terms, cache, visiting = nil|
    visiting ||= {}
    return cache[term_id] if cache.include?(term_id)
    return [] if visiting[term_id]
    term = terms[term_id]
    return cache[term_id] = [] if term.nil?

    visiting[term_id] = true
    ancestors = [term_id]
    term[:parents].each do |parent|
      ancestors.concat(go_ancestors_for_term(parent, terms, cache, visiting))
    end
    visiting.delete(term_id)
    cache[term_id] = ancestors.uniq
  end

  helper :hypergeom_upper_tail do |k, m, n, population|
    k = k.to_i
    m = m.to_i
    n = n.to_i
    population = population.to_i
    return 1.0 if k <= 0 || m <= 0 || n <= 0 || population <= 0
    max_k = [m, n].min
    return 1.0 if k > max_k

    log_choose = lambda do |a,b|
      return -Float::INFINITY if b < 0 || b > a
      Math.lgamma(a + 1)[0] - Math.lgamma(b + 1)[0] - Math.lgamma(a - b + 1)[0]
    end

    denom = log_choose.call(population, n)
    logs = (k..max_k).collect do |i|
      log_choose.call(m, i) + log_choose.call(population - m, n - i) - denom
    end
    max_log = logs.max
    return 0.0 if max_log == -Float::INFINITY
    pvalue = Math.exp(max_log) * logs.collect{|l| Math.exp(l - max_log) }.inject(0.0, &:+)
    [[pvalue, 1.0].min, 0.0].max
  end

  helper :bh_adjust_values do |pvalues|
    indexed = pvalues.each_with_index.select{|pvalue,i| ! pvalue.nil? }.sort_by{|pvalue,i| pvalue.to_f }
    adjusted = Array.new(pvalues.length)
    m = indexed.length
    min_q = 1.0
    indexed.reverse.each_with_index do |(pvalue, i), rank_from_end|
      rank = m - rank_from_end
      qvalue = pvalue.to_f * m / rank
      min_q = [min_q, qvalue].min
      adjusted[i] = [min_q, 1.0].min
    end
    adjusted
  end

  helper :result_treatment_order do
    %w(DMSO FiveZ INT_FiveZ_PI INT_PD_PI PD PI)
  end

  helper :standard_treatment_sort_index do |treatment|
    order = result_treatment_order
    order.index(treatment) || AGS::TREATMENTS.index(treatment) || 999
  end

  
  # ScoutCoder: this function does nothing useful, just step(task_name).load
  # suffices
  #helper :load_dependency_tsv do |task_name|
  #  step = step(task_name)
  #  step.load
  #  step.path.tsv
  #end

  dep :fold_changes, :fc_source => 'NTNU'
  task :interval_fold_changes_fc1 => :tsv do
    fc0 = step(:fold_changes).load.transpose('Associated Gene Name')

    treatments = fc0.fields.collect{|field| field.sub(/^FC_/, '').split('.T').first }.uniq
    fields = treatments.collect{|treatment| AGS::TIME_POINTS.collect{|time_point| context_field('FC1', treatment, time_point) } }.flatten

    dumper = TSV::Dumper.new :key_field => 'Associated Gene Name', :fields => fields, :type => :list, :namespace => AGS.organism
    dumper.init

    TSV.traverse fc0, :into => dumper, :bar => self.progress_bar('Computing interval fold changes') do |gene, values|
      values = NamedArray.setup(values, fc0.fields)
      row = fields.collect do |field|
        treatment, time_point = parse_treatment_time_field(field)
        fc_fields = AGS::TIME_POINTS.collect{|t| context_field('FC', treatment, t) }
        fc_values = values.values_at(*fc_fields).collect{|v| scalar_value(v) }
        fc1_value_from_fc0_values(fc_values, time_point)
      end
      [gene, row]
    end
  end

  dep :pvalues, :fc_source => 'NTNU'
  task :interval_pvalue_surrogate_fc1 => :tsv do
    p0 = step(:pvalues).load.transpose('Associated Gene Name')

    p_prefix = p0.fields.first.to_s.start_with?('Pvalue_') ? 'Pvalue' : 'FC'
    treatments = p0.fields.collect{|field| field.sub(/^(?:Pvalue|FC)_/, '').split('.T').first }.uniq
    fields = treatments.collect{|treatment| AGS::TIME_POINTS.collect{|time_point| context_field('PvalueFC1', treatment, time_point) } }.flatten

    dumper = TSV::Dumper.new :key_field => 'Associated Gene Name', :fields => fields, :type => :list, :namespace => AGS.organism
    dumper.init

    TSV.traverse p0, :into => dumper, :bar => self.progress_bar('Computing interval p-value surrogates') do |gene, values|
      values = NamedArray.setup(values, p0.fields)
      row = fields.collect do |field|
        treatment, time_point = parse_treatment_time_field(field)
        p_fields = AGS::TIME_POINTS.collect{|t| context_field(p_prefix, treatment, t) }
        p_values = values.values_at(*p_fields).collect{|v| scalar_value(v) }
        fc1_pvalue_surrogate_from_values(p_values, time_point)
      end
      [gene, row]
    end
  end

  dep :interval_fold_changes_fc1
  dep :interval_pvalue_surrogate_fc1
  input :fc1_threshold, :float, 'Absolute interval log2 fold-change threshold', 0.3
  input :pvalue_threshold, :float, 'Surrogate interval adjusted p-value threshold, used as a strict upper bound', 0.05
  task :interval_de_gene_calls_fc1 => :tsv do |fc1_threshold, pvalue_threshold|
    fc1 = step(:interval_fold_changes_fc1).load
    p1 = step(:interval_pvalue_surrogate_fc1).load

    fields = %w(Gene Treatment Time FC1 PvalueFC1 Direction AbsFC1Threshold PvalueThreshold)
    dumper = TSV::Dumper.new :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism
    dumper.init
    id = 0
    TSV.traverse fc1, :bar => self.progress_bar('Finding interval DE genes'), into: dumper do |gene, values|
      values = NamedArray.setup(values, fc1.fields)
      pvalues = NamedArray.setup(p1[gene], p1.fields)
      res = fc1.fields.collect do |field|
        treatment, time_point = parse_treatment_time_field(field)
        fc_value = numeric_value(values[field])
        next if fc_value.nil?
        p_field = context_field('PvalueFC1', treatment, time_point)
        p_value = numeric_value(pvalues[p_field], 1.0)
        next unless fc_value.abs >= fc1_threshold && p_value < pvalue_threshold
        direction = fc_value > 0 ? 'up' : 'down'
        id += 1
        [id, [gene, treatment, time_point, fc_value, p_value, direction, fc1_threshold, pvalue_threshold]]
      end.compact
      
      res.extend MultipleResult

      res
    end
  end

  dep :interval_de_gene_calls_fc1
  task :interval_de_gene_counts_fc1 => :tsv do
    # ScoutCoder: if dependency returns a tsv, calling 'load' already loads the
    # TSV, no need to use the load_dependency_tsv
    calls = step(:interval_de_gene_calls_fc1).load
    counts = Hash.new(0)
    calls.through do |id, values|
      treatment = values['Treatment']
      time_point = values['Time'].to_i
      direction = values['Direction']
      counts[[treatment, time_point, direction]] += 1
      counts[[treatment, time_point, 'both']] += 1
    end

    tsv = TSV.setup({}, 'ID~Treatment,Time,Direction,Genes')
    id = 0
    treatments = AGS::TREATMENTS & counts.keys.collect{|k| k[0] }.uniq
    treatments.each do |treatment|
      AGS::TIME_POINTS.each do |time_point|
        %w(up down both).each do |direction|
          id += 1
          tsv[id] = [treatment, time_point, direction, counts[[treatment, time_point, direction]] || 0]
        end
      end
    end
    tsv
  end

  dep :interval_de_gene_calls_fc1
  task :interval_de_gene_lists_fc1 => :array do
    calls = step(:interval_de_gene_calls_fc1).load
    lists = Hash.new{|h,k| h[k] = [] }
    calls.through do |id, values|
      treatment = values['Treatment']
      time_point = values['Time'].to_i
      direction = values['Direction']
      gene = values['Gene']
      lists[[treatment, time_point, direction]] << gene
      lists[[treatment, time_point, 'both']] << gene
    end

    paths = []
    lists.each do |(treatment, time_point, direction), genes|
      path = file([treatment, "T#{time_point}", direction] * '-' + '.list')
      Open.write(path, genes.uniq.sort * "\n" + "\n")
      paths << path
    end
    paths
  end

  dep :full_gene_info
  dep :interval_fold_changes_fc1
  dep :interval_pvalue_surrogate_fc1
  input :fc1_threshold, :float, 'Absolute interval log2 fold-change threshold', 0.3
  input :pvalue_threshold, :float, 'Surrogate interval adjusted p-value threshold, used as a strict upper bound', 0.05
  task :fc1_onset_relationship_long => :tsv do |fc1_threshold, pvalue_threshold|
    info = step(:full_gene_info).load
    fc1 = step(:interval_fold_changes_fc1).load
    p1 = step(:interval_pvalue_surrogate_fc1).load

    fields = %w(Gene Treatment Time FC1 PvalueFC1 FC1Direction FC1IntervalDE OnsetLabel OnsetDirectionAtTime FirstSameDirectionOnsetTime Relationship ProteinCoding INSPEcTGene DynamicGene)
    dumper = TSV::Dumper.new :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism
    dumper.init

    id = 0
    TSV.traverse info, :bar => self.progress_bar('Comparing interval DE and onset calls'), into: dumper do |gene, values|
      next unless fc1.include?(gene)
      info_values = NamedArray.setup(values, info.fields)
      fc_values = NamedArray.setup(fc1[gene], fc1.fields)
      p_values = NamedArray.setup(p1[gene], p1.fields)

      res = []
      AGS::TREATMENTS.each do |treatment|
        onset_label = info_values["#{treatment}: FC clusters"]
        events = onset_events(onset_label)
        AGS::TIME_POINTS.each do |time_point|
          fc_field = context_field('FC1', treatment, time_point)
          p_field = context_field('PvalueFC1', treatment, time_point)
          fc_value = numeric_value(fc_values[fc_field])
          p_value = numeric_value(p_values[p_field], 1.0)
          next if fc_value.nil?

          fc1_de = fc_value.abs >= fc1_threshold && p_value < pvalue_threshold
          fc1_direction = fc_value > 0 ? 'up' : (fc_value < 0 ? 'down' : nil)
          onset_dir_time = onset_direction_at_time(onset_label, time_point)
          onset_starts_here = events.select{|event| event[:start_time] == time_point }

          same_dir_events = fc1_direction ? events.select{|event| event[:direction] == fc1_direction } : []
          first_same = same_dir_events.collect{|event| event[:start_time] }.min
          relationship = nil

          if fc1_de
            if same_dir_events.any?{|event| event[:start_time] == time_point }
              relationship = 'fc1_DE_and_matching_onset'
            elsif same_dir_events.any?{|event| event[:start_time] < time_point }
              relationship = 'fc1_DE_but_onset_earlier'
            elsif same_dir_events.any?{|event| event[:start_time] > time_point }
              relationship = 'fc1_DE_but_onset_later'
            elsif onset_starts_here.any?{|event| event[:direction] != fc1_direction }
              relationship = 'fc1_DE_but_opposite_onset'
            elsif events.empty?
              relationship = 'fc1_DE_but_unclassified_onset'
            else
              relationship = 'fc1_DE_but_other_onset_pattern'
            end
          elsif onset_starts_here.any?
            relationship = 'onset_without_fc1_DE'
          else
            next
          end

          id += 1
          res << [id, [gene, treatment, time_point, fc_value, p_value, fc1_direction, fc1_de, [onset_label].flatten.compact * '|', onset_dir_time, first_same, relationship, scalar_value(info_values['Protein coding']), scalar_value(info_values['INSPEcT gene']), scalar_value(info_values['Dynamic gene'])]]
        end
      end

      res.extend MultipleResult
      res
    end
  end

  dep :fc1_onset_relationship_long
  task :fc1_onset_relationship_summary => :tsv do
    rel = step(:fc1_onset_relationship_long).load
    counts = Hash.new(0)
    rel.through do |id, values|
      treatment = values['Treatment']
      time_point = values['Time'].to_i
      direction = values['FC1Direction']
      direction = values['OnsetDirectionAtTime'] if direction.nil? || direction.empty?
      relationship = values['Relationship']
      counts[[treatment, time_point, direction, relationship]] += 1
    end

    tsv = TSV.setup({}, 'ID~Treatment,Time,Direction,Relationship,Genes')
    id = 0
    counts.keys.sort_by{|treatment,time,direction,relationship| [AGS::TREATMENTS.index(treatment) || 999, time, direction.to_s, relationship.to_s] }.each do |key|
      id += 1
      tsv[id] = key + [counts[key]]
    end
    tsv
  end

  dep :fc1_onset_relationship_long
  task :fc1_onset_relationship_8h_focus => :tsv do
    rel = step(:fc1_onset_relationship_long).load
    counts = Hash.new(0)
    rel.through do |id, values|
      next unless values['Time'].to_i == 8
      next unless values['FC1IntervalDE'].to_s == 'true'
      treatment = values['Treatment']
      direction = values['FC1Direction']
      relationship = values['Relationship']
      counts[[treatment, direction, relationship]] += 1
    end

    tsv = TSV.setup({}, 'ID~Treatment,Direction,Relationship,Genes')
    id = 0
    counts.keys.sort_by{|treatment,direction,relationship| [AGS::TREATMENTS.index(treatment) || 999, direction.to_s, relationship.to_s] }.each do |key|
      id += 1
      tsv[id] = key + [counts[key]]
    end
    tsv
  end

  dep :tf_predictions
  dep :full_gene_info
  dep :regulome
  dep :interval_fold_changes_fc1
  dep :interval_pvalue_surrogate_fc1
  dep :treatment_tf_consistency, treatment: :placeholder do |jobname, options|
    AGS::TREATMENTS.collect{|treatment| options.merge(:treatment => treatment) }
  end
  input :fc1_threshold, :float, 'Absolute interval log2 fold-change threshold for target support', 0.3
  input :pvalue_threshold, :float, 'Surrogate interval adjusted p-value threshold for target support, used as a strict upper bound', 0.05
  task :tf_timepoint_report_card => :tsv do |fc1_threshold, pvalue_threshold|
    predictions = step(:tf_predictions).load
    info = step(:full_gene_info).load
    regulome = step(:regulome).load
    fc1 = step(:interval_fold_changes_fc1).load
    p1 = step(:interval_pvalue_surrogate_fc1).load

    consistency_by_treatment = {}
    dependencies.select{|dep| dep.task_name == :treatment_tf_consistency }.each do |dep|
      consistency_by_treatment[dep.recursive_inputs[:treatment]] = dep.load
    end

    targets_by_tf = Hash.new{|h,k| h[k] = [] }
    regulome.through do |id, values|
      tf, target, weight = values.values_at('source', 'target', 'weight')
      tf ||= values[0]
      target ||= values[1]
      weight ||= values[2]
      targets_by_tf[tf] << [target, weight.to_f]
    end

    fields = %w(TF Treatment Time TFActivityScore TFActivitySign TFActivityCalled TFGeneFC0 TFGeneFC1 TFGenePvalueFC0 TFGenePvalueFC1Surrogate TFGeneOnsetLabel TFGeneOnsetDirectionAtTime SelfConsistency TotalRegulomeTargets TargetsInGeneInfo DynamicTargetsAnyTime DynamicTargetsThisTime PositiveEdgesTotal NegativeEdgesTotal PositiveEdgesDynamicThisTime NegativeEdgesDynamicThisTime TargetsExpectedUp TargetsExpectedDown TargetsObservedUpFC1 TargetsObservedDownFC1 TargetsObservedUpOnset TargetsObservedDownOnset TargetSignConcordantCount TargetSignDiscordantCount TargetSignUnknownCount TargetConcordanceFraction DynamicTargetConcordanceFraction TopConcordantTargets TopDiscordantTargets)
    dumper = TSV::Dumper.new :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism
    dumper.init

    id = 0
    TSV.traverse predictions, :bar => self.progress_bar('Building TF timepoint report card'), :into => dumper do |tf, pred_values|
      res = []
      pred_values = NamedArray.setup(pred_values, predictions.fields)
      predictions.fields.each do |field|
        parsed = parse_treatment_time_field(field)
        next if parsed.nil?
        treatment, time_point = parsed
        activity = numeric_value(pred_values[field], 0.0)
        activity_sign = activity > 0 ? 'positive' : (activity < 0 ? 'negative' : 'zero')
        activity_called = activity != 0

        tf_info = info[tf]
        tf_info = NamedArray.setup(tf_info, info.fields) if tf_info
        tf_fc0 = tf_info ? numeric_value(tf_info[context_field('FC', treatment, time_point)]) : nil
        tf_p0 = tf_info ? numeric_value(tf_info[context_field('Pvalue', treatment, time_point)]) : nil
        tf_fc1 = fc1.include?(tf) ? numeric_value(fc1[tf][context_field('FC1', treatment, time_point)]) : nil
        tf_p1 = p1.include?(tf) ? numeric_value(p1[tf][context_field('PvalueFC1', treatment, time_point)], 1.0) : nil
        tf_onset = tf_info ? [tf_info["#{treatment}: FC clusters"]].flatten.compact * '|' : nil
        tf_onset_dir = tf_info ? onset_direction_at_time(tf_info["#{treatment}: FC clusters"], time_point) : nil

        consistency = nil
        if consistency_by_treatment[treatment] && consistency_by_treatment[treatment].include?(tf)
          consistency = scalar_value(consistency_by_treatment[treatment][tf]["Consistent at #{time_point}h"])
        end

        target_edges = targets_by_tf[tf]
        total_targets = target_edges.length
        positive_edges_total = target_edges.select{|target, weight| weight > 0 }.length
        negative_edges_total = target_edges.select{|target, weight| weight < 0 }.length

        targets_in_info = 0
        dynamic_any = 0
        dynamic_this = 0
        pos_edges_dynamic = 0
        neg_edges_dynamic = 0
        expected_up = 0
        expected_down = 0
        obs_up_fc1 = 0
        obs_down_fc1 = 0
        obs_up_onset = 0
        obs_down_onset = 0
        concordant = 0
        discordant = 0
        unknown = 0
        dyn_concordant = 0
        dyn_discordant = 0
        concordant_targets = []
        discordant_targets = []

        target_edges.each do |target, weight|
          target_info = info[target]
          next if target_info.nil?
          targets_in_info += 1
          target_info = NamedArray.setup(target_info, info.fields)
          events = onset_events(target_info["#{treatment}: FC clusters"])
          target_dynamic_any = events.any?
          target_dynamic_this = events.any?{|event| event[:start_time] == time_point }
          dynamic_any += 1 if target_dynamic_any
          dynamic_this += 1 if target_dynamic_this
          if target_dynamic_this
            pos_edges_dynamic += 1 if weight > 0
            neg_edges_dynamic += 1 if weight < 0
          end

          target_fc1 = fc1.include?(target) ? numeric_value(fc1[target][context_field('FC1', treatment, time_point)]) : nil
          target_p1 = p1.include?(target) ? numeric_value(p1[target][context_field('PvalueFC1', treatment, time_point)], 1.0) : nil
          target_fc1_de = target_fc1 && target_fc1.abs >= fc1_threshold && target_p1 < pvalue_threshold
          observed_fc1_dir = target_fc1_de ? (target_fc1 > 0 ? 'up' : 'down') : nil
          obs_up_fc1 += 1 if observed_fc1_dir == 'up'
          obs_down_fc1 += 1 if observed_fc1_dir == 'down'

          observed_onset_dir = onset_direction_at_time(target_info["#{treatment}: FC clusters"], time_point)
          obs_up_onset += 1 if observed_onset_dir == 'up'
          obs_down_onset += 1 if observed_onset_dir == 'down'

          expected_dir = nil
          if activity > 0
            expected_dir = weight > 0 ? 'up' : 'down'
          elsif activity < 0
            expected_dir = weight > 0 ? 'down' : 'up'
          end
          expected_up += 1 if expected_dir == 'up'
          expected_down += 1 if expected_dir == 'down'

          observed_dir = observed_onset_dir || observed_fc1_dir
          if expected_dir.nil? || observed_dir.nil? || observed_dir == 'both'
            unknown += 1
          elsif expected_dir == observed_dir
            concordant += 1
            dyn_concordant += 1 if target_dynamic_this
            concordant_targets << target if concordant_targets.length < 20
          else
            discordant += 1
            dyn_discordant += 1 if target_dynamic_this
            discordant_targets << target if discordant_targets.length < 20
          end
        end

        informative = concordant + discordant
        dyn_informative = dyn_concordant + dyn_discordant
        target_fraction = informative == 0 ? nil : concordant.to_f / informative
        dynamic_fraction = dyn_informative == 0 ? nil : dyn_concordant.to_f / dyn_informative

        id += 1
        res << [id, [tf, treatment, time_point, activity, activity_sign, activity_called, tf_fc0, tf_fc1, tf_p0, tf_p1, tf_onset, tf_onset_dir, consistency, total_targets, targets_in_info, dynamic_any, dynamic_this, positive_edges_total, negative_edges_total, pos_edges_dynamic, neg_edges_dynamic, expected_up, expected_down, obs_up_fc1, obs_down_fc1, obs_up_onset, obs_down_onset, concordant, discordant, unknown, target_fraction, dynamic_fraction, concordant_targets * '|', discordant_targets * '|']]
      end

      res.extend MultipleResult
      res
    end
  end

  dep :tf_predictions
  dep :full_gene_info
  dep :regulome
  dep :interval_fold_changes_fc1
  dep :interval_pvalue_surrogate_fc1
  input :called_only, :boolean, 'Only report targets for non-zero TF activity calls', true
  input :dynamic_targets_only, :boolean, 'Only report targets with an onset at the treatment-timepoint', false
  input :fc1_threshold, :float, 'Absolute interval log2 fold-change threshold for target support', 0.3
  input :pvalue_threshold, :float, 'Surrogate interval adjusted p-value threshold for target support, used as a strict upper bound', 0.05
  task :tf_target_report_card => :tsv do |called_only, dynamic_targets_only, fc1_threshold, pvalue_threshold|
    predictions = step(:tf_predictions).load
    info = step(:full_gene_info).load
    regulome = step(:regulome).load
    fc1 = step(:interval_fold_changes_fc1).load
    p1 = step(:interval_pvalue_surrogate_fc1).load

    target_edges_by_tf = Hash.new{|h,k| h[k] = [] }
    regulome.through do |id, values|
      tf, target, weight = values.values_at('source', 'target', 'weight')
      tf ||= values[0]
      target ||= values[1]
      weight ||= values[2]
      target_edges_by_tf[tf] << [target, weight.to_f]
    end

    fields = %w(TF Target Treatment Time TFActivityScore TFActivitySign RegulomeWeight RegulomeEffectSign ExpectedTargetDirection TargetFC0 TargetFC1 TargetPvalueFC0 TargetPvalueFC1Surrogate TargetOnsetLabel TargetOnsetDirectionAtTime TargetDynamicThisTime TargetDynamicAnyTime TargetFC1IntervalDE ConcordantWithTFActivity ProteinCoding INSPEcTGene DynamicGene)
    dumper = TSV::Dumper.new :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism
    dumper.init

    id = 0
    TSV.traverse predictions, :bar => self.progress_bar('Building TF target report card'), :into => dumper do |tf, pred_values|
      res = []
      pred_values = NamedArray.setup(pred_values, predictions.fields)
      target_edges = target_edges_by_tf[tf]
      next if target_edges.empty?

      predictions.fields.each do |field|
        parsed = parse_treatment_time_field(field)
        next if parsed.nil?
        treatment, time_point = parsed
        activity = numeric_value(pred_values[field], 0.0)
        next if called_only && activity == 0
        activity_sign = activity > 0 ? 'positive' : (activity < 0 ? 'negative' : 'zero')

        target_edges.each do |target, weight|
          target_info = info[target]
          next if target_info.nil?
          target_info = NamedArray.setup(target_info, info.fields)
          onset_label = target_info["#{treatment}: FC clusters"]
          events = onset_events(onset_label)
          target_dynamic_this = events.any?{|event| event[:start_time] == time_point }
          target_dynamic_any = events.any?
          next if dynamic_targets_only && ! target_dynamic_this

          expected_dir = nil
          if activity > 0
            expected_dir = weight > 0 ? 'up' : 'down'
          elsif activity < 0
            expected_dir = weight > 0 ? 'down' : 'up'
          end

          target_fc0 = numeric_value(target_info[context_field('FC', treatment, time_point)])
          target_p0 = numeric_value(target_info[context_field('Pvalue', treatment, time_point)])
          target_fc1 = fc1.include?(target) ? numeric_value(fc1[target][context_field('FC1', treatment, time_point)]) : nil
          target_p1 = p1.include?(target) ? numeric_value(p1[target][context_field('PvalueFC1', treatment, time_point)], 1.0) : nil
          target_fc1_de = target_fc1 && target_fc1.abs >= fc1_threshold && target_p1 < pvalue_threshold
          observed_onset_dir = onset_direction_at_time(onset_label, time_point)
          observed_fc1_dir = target_fc1_de ? (target_fc1 > 0 ? 'up' : 'down') : nil
          observed_dir = observed_onset_dir || observed_fc1_dir

          concordant = nil
          unless expected_dir.nil? || observed_dir.nil? || observed_dir == 'both'
            concordant = expected_dir == observed_dir
          end

          id += 1
          res << [id, [tf, target, treatment, time_point, activity, activity_sign, weight, weight > 0 ? 'positive' : 'negative', expected_dir, target_fc0, target_fc1, target_p0, target_p1, [onset_label].flatten.compact * '|', observed_onset_dir, target_dynamic_this, target_dynamic_any, target_fc1_de, concordant, scalar_value(target_info['Protein coding']), scalar_value(target_info['INSPEcT gene']), scalar_value(target_info['Dynamic gene'])]]
        end
      end

      res.extend MultipleResult
      res
    end
  end

  dep :tf_target_report_card
  task :tf_target_edge_consistency_summary => :tsv do

    # ScoutCoder: the tf_target_report_card is a very large file
    # se we stream it using TSV.traverse instead of loading
    # the tsv to memory and traversing the TSV.
    #
    # Also, when using TSV.traverse values are not extended with 
    # NamedArray, which also improves performance, but we can't use
    # that inteface to extract the correct field from the result.
    # we get the field index before hand from the file headers.
    #
    # We use the TSV::Parser object to be able to read the header
    # lines at the start of the stream and still be able to traverse
    # the reminder stream, which avoids reopening or rewinding the
    # stream, which may not always be posible.
    
    parser = TSV::Parser.new step(:tf_target_report_card)

    treatment_field_num = parser.fields.index 'Treatment'
    time_field_num = parser.fields.index 'Time'
    dynamic_target_field_num = parser.fields.index 'TargetDynamicThisTime'
    concordant_field_num = parser.fields.index 'ConcordantWithTFActivity'

    counts = Hash.new{|h,k| h[k] = Hash.new(0) }
    TSV.traverse parser, bar: self.progress_bar("Processing tf_target_report_card") do |id, values| 
      treatment = values[treatment_field_num]
      time_point = values[time_field_num].to_i
      dynamic = values[dynamic_target_field_num].to_s == 'true' ? 'dynamic_target' : 'other_target'
      concordant = values[concordant_field_num]
      category = if concordant.to_s == 'true'
                   'concordant'
                 elsif concordant.to_s == 'false'
                   'discordant'
                 else
                   'unknown'
                 end
      counts[[treatment, time_point, dynamic]][category] += 1
    end

    tsv = TSV.setup({}, 'ID~Treatment,Time,TargetClass,Concordant,Discordant,Unknown,Informative,ConcordanceFraction')
    id = 0
    counts.keys.sort_by{|treatment,time,target_class| [AGS::TREATMENTS.index(treatment) || 999, time, target_class] }.each do |key|
      c = counts[key]
      informative = c['concordant'] + c['discordant']
      fraction = informative == 0 ? nil : c['concordant'].to_f / informative
      id += 1
      tsv[id] = key + [c['concordant'], c['discordant'], c['unknown'], informative, fraction]
    end
    tsv
  end

  dep :fold_changes, :fc_source => 'NTNU'
  dep :pvalues, :fc_source => 'NTNU'
  input :fc0_threshold, :float, 'Absolute cumulative log2 fold-change threshold', 0.0
  input :pvalue_threshold, :float, 'Adjusted p-value threshold, used as a strict upper bound', 0.05
  task :de_gene_counts_fc0 => :tsv do |fc0_threshold, pvalue_threshold|
    fc0 = step(:fold_changes).load.transpose('Associated Gene Name')
    p0 = step(:pvalues).load.transpose('Associated Gene Name')
    p_prefix = p0.fields.first.to_s.start_with?('Pvalue_') ? 'Pvalue' : 'FC'

    counts = Hash.new(0)
    fc0.fields.each do |field|
      parsed = parse_treatment_time_field(field)
      next if parsed.nil?
      treatment, time_point = parsed
      p_field = context_field(p_prefix, treatment, time_point)
      fc0.through do |gene, values|
        fc_value = numeric_value(values[field])
        next if fc_value.nil?
        p_value = p0.include?(gene) ? numeric_value(p0[gene][p_field], 1.0) : 1.0
        next unless fc_value.abs >= fc0_threshold && p_value < pvalue_threshold
        direction = fc_value > 0 ? 'up' : 'down'
        counts[[treatment, time_point, direction]] += 1
        counts[[treatment, time_point, 'both']] += 1
      end
    end

    tsv = TSV.setup({}, 'ID~Treatment,Time,Direction,Genes,AbsFC0Threshold,PvalueThreshold')
    id = 0
    result_treatment_order.each do |treatment|
      AGS::TIME_POINTS.each do |time_point|
        %w(up down both).each do |direction|
          id += 1
          tsv[id] = [treatment, time_point, direction, counts[[treatment, time_point, direction]] || 0, fc0_threshold, pvalue_threshold]
        end
      end
    end
    tsv
  end

  dep :full_gene_info
  task :onset_first_counts => :tsv do
    info = step(:full_gene_info).load
    counts = Hash.new(0)
    info.through do |gene, values|
      values = NamedArray.setup(values, info.fields)
      result_treatment_order.each do |treatment|
        events = onset_events(values["#{treatment}: FC clusters"])
        next if events.empty?
        first = events.sort_by{|event| [AGS::TIME_POINTS.index(event[:start_time]) || 999, event[:label]] }.first
        counts[[treatment, first[:start_time], first[:direction]]] += 1
        counts[[treatment, first[:start_time], 'both']] += 1
      end
    end

    tsv = TSV.setup({}, 'ID~Treatment,Time,Direction,Genes')
    id = 0
    result_treatment_order.each do |treatment|
      AGS::TIME_POINTS.each do |time_point|
        %w(up down both).each do |direction|
          id += 1
          tsv[id] = [treatment, time_point, direction, counts[[treatment, time_point, direction]] || 0]
        end
      end
    end
    tsv
  end

  dep :full_gene_info
  task :onset_episode_counts => :tsv do
    info = step(:full_gene_info).load
    counts = Hash.new(0)
    info.through do |gene, values|
      values = NamedArray.setup(values, info.fields)
      result_treatment_order.each do |treatment|
        onset_events(values["#{treatment}: FC clusters"]).each do |event|
          counts[[treatment, event[:start_time], event[:direction]]] += 1
          counts[[treatment, event[:start_time], 'both']] += 1
        end
      end
    end

    tsv = TSV.setup({}, 'ID~Treatment,Time,Direction,Episodes')
    id = 0
    result_treatment_order.each do |treatment|
      AGS::TIME_POINTS.each do |time_point|
        %w(up down both).each do |direction|
          id += 1
          tsv[id] = [treatment, time_point, direction, counts[[treatment, time_point, direction]] || 0]
        end
      end
    end
    tsv
  end

  dep :full_gene_info
  task :onset_direction_switch_summary => :tsv do
    info = step(:full_gene_info).load
    counts = Hash.new(0)
    info.through do |gene, values|
      values = NamedArray.setup(values, info.fields)
      result_treatment_order.each do |treatment|
        events = onset_events(values["#{treatment}: FC clusters"])
        episode_count = events.length
        directions = events.collect{|event| event[:direction] }.uniq.length
        category = if episode_count == 0
                     'unclassified'
                   elsif episode_count == 1
                     'single_episode'
                   elsif directions == 1
                     'multiple_same_direction'
                   else
                     'direction_switch'
                   end
        counts[[treatment, category]] += 1
      end
    end

    tsv = TSV.setup({}, 'ID~Treatment,Category,Genes')
    id = 0
    result_treatment_order.each do |treatment|
      %w(unclassified single_episode multiple_same_direction direction_switch).each do |category|
        id += 1
        tsv[id] = [treatment, category, counts[[treatment, category]] || 0]
      end
    end
    tsv
  end

  dep :tf_predictions
  task :tf_activity_call_counts_dynamic => :tsv do
    predictions = step(:tf_predictions).load
    counts = Hash.new(0)
    predictions.fields.each do |field|
      parsed = parse_treatment_time_field(field)
      next if parsed.nil?
      treatment, time_point = parsed
      predictions.through do |tf, values|
        activity = numeric_value(values[field], 0.0)
        next if activity == 0
        sign = activity > 0 ? 'positive' : 'negative'
        counts[[treatment, time_point, sign]] += 1
        counts[[treatment, time_point, 'both']] += 1
      end
    end

    tsv = TSV.setup({}, 'ID~Treatment,Time,Sign,TFActivityCalls')
    id = 0
    result_treatment_order.each do |treatment|
      AGS::TIME_POINTS.each do |time_point|
        %w(positive negative both).each do |sign|
          id += 1
          tsv[id] = [treatment, time_point, sign, counts[[treatment, time_point, sign]] || 0]
        end
      end
    end
    tsv
  end

  dep :tf_predictions
  task :combination_tf_categories => :tsv do
    predictions = step(:tf_predictions).load
    combinations = {
      'INT_PD_PI' => ['PI', 'PD'],
      'INT_FiveZ_PI' => ['PI', 'FiveZ']
    }

    first_call = Hash.new{|h,k| h[k] = {} }
    predictions.through do |tf, values|
      values = NamedArray.setup(values, predictions.fields)
      result_treatment_order.each do |treatment|
        AGS::TIME_POINTS.each do |time_point|
          field = "#{treatment}-T#{time_point}"
          next unless predictions.fields.include?(field)
          activity = numeric_value(values[field], 0.0)
          next if activity == 0
          first_call[[tf, treatment]][activity_sign_label(activity)] ||= time_point
          first_call[[tf, treatment]]['any'] ||= time_point
        end
      end
    end

    fields = %w(TF Combination Component1 Component2 Time CombinationActivity CombinationSign Component1Activity Component1Sign Component2Activity Component2Sign Category CombinationEarlierThanBoth StrongerThanBoth SignReversal)
    dumper = TSV::Dumper.new :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism
    dumper.init
    id = 0

    TSV.traverse predictions, :bar => self.progress_bar('Classifying combination TF activity calls'), :into => dumper do |tf, values|
      values = NamedArray.setup(values, predictions.fields)
      res = []
      combinations.each do |combination, components|
        comp1, comp2 = components
        AGS::TIME_POINTS.each do |time_point|
          combo_field = "#{combination}-T#{time_point}"
          next unless predictions.fields.include?(combo_field)
          combo_activity = numeric_value(values[combo_field], 0.0)
          next if combo_activity == 0
          combo_sign = activity_sign_label(combo_activity)
          comp1_field = "#{comp1}-T#{time_point}"
          comp2_field = "#{comp2}-T#{time_point}"
          comp1_activity = predictions.fields.include?(comp1_field) ? numeric_value(values[comp1_field], 0.0) : 0.0
          comp2_activity = predictions.fields.include?(comp2_field) ? numeric_value(values[comp2_field], 0.0) : 0.0
          comp1_sign = activity_sign_label(comp1_activity)
          comp2_sign = activity_sign_label(comp2_activity)

          comp1_same = comp1_activity != 0 && comp1_sign == combo_sign
          comp2_same = comp2_activity != 0 && comp2_sign == combo_sign
          sign_reversal = (comp1_activity != 0 && comp1_sign != combo_sign) || (comp2_activity != 0 && comp2_sign != combo_sign)
          stronger = combo_activity.abs > comp1_activity.abs && combo_activity.abs > comp2_activity.abs
          comp1_first = first_call[[tf, comp1]][combo_sign] || first_call[[tf, comp1]]['any']
          comp2_first = first_call[[tf, comp2]][combo_sign] || first_call[[tf, comp2]]['any']
          earlier = (comp1_first.nil? || AGS::TIME_POINTS.index(time_point) < AGS::TIME_POINTS.index(comp1_first)) && (comp2_first.nil? || AGS::TIME_POINTS.index(time_point) < AGS::TIME_POINTS.index(comp2_first))

          category = if earlier
                       'combination_earlier_than_both'
                     elsif comp1_same && comp2_same
                       'shared_with_both_same_sign'
                     elsif comp1_same
                       'shared_with_component1_same_sign'
                     elsif comp2_same
                       'shared_with_component2_same_sign'
                     elsif sign_reversal
                       'sign_reversed_relative_to_component'
                     else
                       'combination_specific_at_time'
                     end

          id += 1
          res << [id, [tf, combination, comp1, comp2, time_point, combo_activity, combo_sign, comp1_activity, comp1_sign, comp2_activity, comp2_sign, category, earlier, stronger, sign_reversal]]
        end
      end
      res.extend MultipleResult
      res
    end
  end

  dep :combination_tf_categories
  task :combination_tf_category_counts => :tsv do
    categories = step(:combination_tf_categories).load
    counts = Hash.new(0)
    categories.through do |id, values|
      combination = values['Combination']
      time_point = values['Time'].to_i
      category = values['Category']
      counts[[combination, time_point, category]] += 1
    end

    tsv = TSV.setup({}, 'ID~Combination,Time,Category,TFActivityCalls')
    id = 0
    %w(INT_FiveZ_PI INT_PD_PI).each do |combination|
      AGS::TIME_POINTS.each do |time_point|
        counts.keys.collect{|k| k[2] }.uniq.sort.each do |category|
          id += 1
          tsv[id] = [combination, time_point, category, counts[[combination, time_point, category]] || 0]
        end
      end
    end
    tsv
  end


  input :scheme, :select, 'TF prediction scheme to summarize', 'dynamic', :select_options => %w(dynamic non-dynamic)
  dep :tf_predictions, :scheme => 'dynamic'
  dep :tf_predictions, :scheme => 'non-dynamic'
  task :tf_activity_call_counts_by_scheme => :tsv do |scheme|
    dep = dependencies.find{|d| d.recursive_inputs[:scheme].to_s == scheme.to_s } || dependencies.first
    predictions = dep.load
    counts = Hash.new(0)
    predictions.fields.each do |field|
      parsed = parse_treatment_time_field(field)
      next if parsed.nil?
      treatment, time_point = parsed
      predictions.through do |tf, values|
        activity = numeric_value(values[field], 0.0)
        next if activity == 0
        sign = activity > 0 ? 'positive' : 'negative'
        counts[[treatment, time_point, sign]] += 1
        counts[[treatment, time_point, 'both']] += 1
      end
    end

    tsv = TSV.setup({}, 'ID~Scheme,Treatment,Time,Sign,TFActivityCalls')
    id = 0
    result_treatment_order.each do |treatment|
      AGS::TIME_POINTS.each do |time_point|
        %w(positive negative both).each do |sign|
          id += 1
          tsv[id] = [scheme, treatment, time_point, sign, counts[[treatment, time_point, sign]] || 0]
        end
      end
    end
    tsv
  end

  input :scheme, :select, 'TF prediction scheme to export', 'dynamic', :select_options => %w(dynamic non-dynamic)
  input :normalization, :select, 'Matrix normalization', 'raw', :select_options => %w(raw row_zscore column_zscore)
  dep :tf_predictions, :scheme => 'dynamic'
  dep :tf_predictions, :scheme => 'non-dynamic'
  task :tf_activity_heatmap_matrix => :tsv do |scheme, normalization|
    dep = dependencies.find{|d| d.recursive_inputs[:scheme].to_s == scheme.to_s } || dependencies.first
    predictions = dep.load
    fields = predictions.fields.select{|field| parse_treatment_time_field(field) }

    matrix = {}
    predictions.through do |tf, values|
      matrix[tf] = fields.collect{|field| numeric_value(values[field], 0.0) }
    end

    case normalization
    when 'row_zscore'
      matrix.each do |tf, row|
        mean = row.inject(0.0, &:+) / row.length
        variance = row.collect{|v| (v - mean) ** 2 }.inject(0.0, &:+) / row.length
        sd = Math.sqrt(variance)
        matrix[tf] = sd == 0 ? row.collect{0.0} : row.collect{|v| (v - mean) / sd }
      end
    when 'column_zscore'
      means = []
      sds = []
      fields.each_index do |i|
        col = matrix.values.collect{|row| row[i] }
        mean = col.inject(0.0, &:+) / col.length
        variance = col.collect{|v| (v - mean) ** 2 }.inject(0.0, &:+) / col.length
        means[i] = mean
        sds[i] = Math.sqrt(variance)
      end
      matrix.each do |tf, row|
        matrix[tf] = row.each_with_index.collect{|v,i| sds[i] == 0 ? 0.0 : (v - means[i]) / sds[i] }
      end
    end

    tsv = TSV.setup({}, :key_field => 'Associated Gene Name', :fields => fields, :type => :list, :namespace => AGS.organism)
    matrix.each{|tf, row| tsv[tf] = row }
    tsv
  end


  dep :neko_bootstrap_sweep
  task :neko_dynamic_non_dynamic_summary => :tsv do
    sweep = step(:neko_bootstrap_sweep).load
    idx = Hash[sweep.fields.each_with_index.to_a]
    tsv = TSV.setup({}, 'ID~Treatment,Scheme,Vetting,Target,Match,Miss,Total,MatchFraction,MatchOdds')
    id = 0
    sweep.through do |row_id, values|
      scheme = scalar_value(scalar_value(values[idx['Scheme']]))
      next unless %w(dynamic non-dynamic).include?(scheme)
      treatment = scalar_value(scalar_value(values[idx['Treatment']]))
      vetting = scalar_value(scalar_value(values[idx['Vetting']]))
      target = scalar_value(scalar_value(values[idx['Target']]))
      match = scalar_value(scalar_value(values[idx['Match']])).to_f
      miss = scalar_value(scalar_value(values[idx['Miss']])).to_f
      total = scalar_value(scalar_value(values[idx['Total']])).to_f
      next if total == 0
      match_fraction = match / total
      match_odds = (total - match) == 0 ? nil : match / (total - match)
      id += 1
      tsv[id] = [treatment, scheme, vetting, target, match.to_i, miss.to_i, total.to_i, match_fraction, match_odds]
    end
    tsv
  end

  dep :neko_dynamic_non_dynamic_summary
  task :neko_dynamic_vs_non_dynamic_odds => :tsv do
    data = step(:neko_dynamic_non_dynamic_summary).load
    idx = Hash[data.fields.each_with_index.to_a]
    by_key = Hash.new{|h,k| h[k] = {} }
    data.through do |id, values|
      key = [scalar_value(scalar_value(values[idx['Treatment']])), scalar_value(scalar_value(values[idx['Vetting']])), scalar_value(scalar_value(values[idx['Target']]))]
      by_key[key][scalar_value(scalar_value(values[idx['Scheme']]))] = values
    end

    tsv = TSV.setup({}, 'ID~Treatment,Vetting,Target,DynamicMatch,DynamicMiss,DynamicTotal,NonDynamicMatch,NonDynamicMiss,NonDynamicTotal,DynamicMatchFraction,NonDynamicMatchFraction,DynamicVsNonDynamicOddsRatio')
    id = 0
    by_key.keys.sort.each do |key|
      dyn = by_key[key]['dynamic']
      nd = by_key[key]['non-dynamic']
      next if dyn.nil? || nd.nil?
      dm = scalar_value(scalar_value(dyn[idx['Match']])).to_f
      dt = scalar_value(scalar_value(dyn[idx['Total']])).to_f
      nm = scalar_value(scalar_value(nd[idx['Match']])).to_f
      nt = scalar_value(scalar_value(nd[idx['Total']])).to_f
      dyn_odds = dt == dm ? nil : dm / (dt - dm)
      nd_odds = nt == nm ? nil : nm / (nt - nm)
      oratio = dyn_odds.nil? || nd_odds.nil? || nd_odds == 0 ? nil : dyn_odds / nd_odds
      id += 1
      tsv[id] = key + [scalar_value(scalar_value(dyn[idx['Match']])), scalar_value(scalar_value(dyn[idx['Miss']])), scalar_value(scalar_value(dyn[idx['Total']])), scalar_value(scalar_value(nd[idx['Match']])), scalar_value(scalar_value(nd[idx['Miss']])), scalar_value(scalar_value(nd[idx['Total']])), scalar_value(scalar_value(dyn[idx['MatchFraction']])), scalar_value(scalar_value(nd[idx['MatchFraction']])), oratio]
    end
    tsv
  end

  dep :consistency_counts
  task :self_consistency_dynamic_non_dynamic_summary => :tsv do
    counts = step(:consistency_counts).load
    idx = Hash[counts.fields.each_with_index.to_a]
    tsv = TSV.setup({}, 'ID~Treatment,Time,Scheme,Vetting,Matches,Miss,Total,MatchFraction,MatchOdds')
    id = 0
    counts.through do |row_id, values|
      scheme = scalar_value(scalar_value(values[idx['Scheme']]))
      next unless %w(dynamic non-dynamic).include?(scheme)
      treatment = scalar_value(scalar_value(values[idx['Treatment']]))
      time = scalar_value(scalar_value(values[idx['Time']]))
      vetting = scalar_value(scalar_value(values[idx['Vetting']]))
      matches = scalar_value(scalar_value(values[idx['Matches']])).to_f
      miss = scalar_value(scalar_value(values[idx['Miss']])).to_f
      total = scalar_value(scalar_value(values[idx['Total']])).to_f
      next if total == 0
      match_fraction = matches / total
      match_odds = (total - matches) == 0 ? nil : matches / (total - matches)
      id += 1
      tsv[id] = [treatment, time, scheme, vetting, matches.to_i, miss.to_i, total.to_i, match_fraction, match_odds]
    end
    tsv
  end

  dep :self_consistency_dynamic_non_dynamic_summary
  task :self_consistency_dynamic_vs_non_dynamic_odds => :tsv do
    data = step(:self_consistency_dynamic_non_dynamic_summary).load
    idx = Hash[data.fields.each_with_index.to_a]
    by_key = Hash.new{|h,k| h[k] = {} }
    data.through do |id, values|
      key = [scalar_value(scalar_value(values[idx['Treatment']])), scalar_value(scalar_value(values[idx['Time']])), scalar_value(scalar_value(values[idx['Vetting']]))]
      by_key[key][scalar_value(scalar_value(values[idx['Scheme']]))] = values
    end

    tsv = TSV.setup({}, 'ID~Treatment,Time,Vetting,DynamicMatches,DynamicMiss,DynamicTotal,NonDynamicMatches,NonDynamicMiss,NonDynamicTotal,DynamicMatchFraction,NonDynamicMatchFraction,DynamicVsNonDynamicOddsRatio')
    id = 0
    by_key.keys.sort_by{|treatment,time,vetting| [standard_treatment_sort_index(treatment), time.to_i, vetting.to_s] }.each do |key|
      dyn = by_key[key]['dynamic']
      nd = by_key[key]['non-dynamic']
      next if dyn.nil? || nd.nil?
      dm = scalar_value(scalar_value(dyn[idx['Matches']])).to_f
      dt = scalar_value(scalar_value(dyn[idx['Total']])).to_f
      nm = scalar_value(scalar_value(nd[idx['Matches']])).to_f
      nt = scalar_value(scalar_value(nd[idx['Total']])).to_f
      dyn_odds = dt == dm ? nil : dm / (dt - dm)
      nd_odds = nt == nm ? nil : nm / (nt - nm)
      oratio = dyn_odds.nil? || nd_odds.nil? || nd_odds == 0 ? nil : dyn_odds / nd_odds
      id += 1
      tsv[id] = key + [scalar_value(scalar_value(dyn[idx['Matches']])), scalar_value(scalar_value(dyn[idx['Miss']])), scalar_value(scalar_value(dyn[idx['Total']])), scalar_value(scalar_value(nd[idx['Matches']])), scalar_value(scalar_value(nd[idx['Miss']])), scalar_value(scalar_value(nd[idx['Total']])), scalar_value(scalar_value(dyn[idx['MatchFraction']])), scalar_value(scalar_value(nd[idx['MatchFraction']])), oratio]
    end
    tsv
  end


  input :go_slim_url, :string, 'URL for the generic GO slim OBO subset', 'https://current.geneontology.org/ontology/subsets/goslim_generic.obo'
  task :go_slim_bp_terms => :tsv do |go_slim_url|
    require 'open-uri'
    slim_text = URI.open(go_slim_url, 'User-Agent' => 'Mozilla/5.0').read
    slim_terms = parse_go_obo_terms(slim_text)

    tsv = TSV.setup({}, 'SlimGOID~SlimName,Namespace')
    slim_terms.keys.sort.each do |go_id|
      term = slim_terms[go_id]
      next unless term[:namespace] == 'biological_process'
      tsv[go_id] = [term[:name], term[:namespace]]
    end
    tsv
  end

  input :go_basic_url, :string, 'URL for go-basic.obo', 'https://current.geneontology.org/ontology/go-basic.obo'
  input :go_slim_url, :string, 'URL for the generic GO slim OBO subset', 'https://current.geneontology.org/ontology/subsets/goslim_generic.obo'
  task :go_bp_to_goslim_map => :tsv do |go_basic_url, go_slim_url|
    require 'open-uri'
    require 'set'

    go_text = URI.open(go_basic_url, 'User-Agent' => 'Mozilla/5.0').read
    slim_text = URI.open(go_slim_url, 'User-Agent' => 'Mozilla/5.0').read

    terms = parse_go_obo_terms(go_text)
    slim_terms = parse_go_obo_terms(slim_text)
    slim_bp_ids = slim_terms.keys.select{|go_id| slim_terms[go_id][:namespace] == 'biological_process' }
    slim_bp_set = slim_bp_ids.to_set
    ancestor_cache = {}

    tsv = TSV.setup({}, :key_field => 'GOID', :fields => %w(GOName Namespace SlimGOIDs SlimNames SlimGOIDsAll SlimNamesAll), :type => :list)
    terms.keys.sort.each do |go_id|
      term = terms[go_id]
      next unless term[:namespace] == 'biological_process'
      all_slim = (go_ancestors_for_term(go_id, terms, ancestor_cache) & slim_bp_ids)
      next if all_slim.empty?

      reduced = all_slim.reject do |candidate|
        all_slim.any? do |other|
          other != candidate && go_ancestors_for_term(other, terms, ancestor_cache).include?(candidate)
        end
      end
      reduced = all_slim if reduced.empty?
      reduced = reduced.sort
      all_slim = all_slim.sort

      tsv[go_id] = [
        term[:name],
        term[:namespace],
        reduced * '|',
        reduced.collect{|slim_id| slim_terms[slim_id] ? slim_terms[slim_id][:name] : terms[slim_id][:name] } * '|',
        all_slim * '|',
        all_slim.collect{|slim_id| slim_terms[slim_id] ? slim_terms[slim_id][:name] : terms[slim_id][:name] } * '|'
      ]
    end
    tsv
  end

  dep :go_bp_to_goslim_map
  task :gene_goslim_bp_annotations => :tsv do
    go_to_slim = step(:go_bp_to_goslim_map).load
    gene_go = Organism.gene_go_bp(AGS.organism).tsv(type: :flat)
    identifiers = Organism.identifiers(AGS.organism).tsv(:key_field => 'Ensembl Gene ID', :fields => ['Associated Gene Name'], :type => :flat)

    slim_name_by_id = {}
    go_to_slim.through do |go_id, values|
      values = NamedArray.setup(values, go_to_slim.fields)
      ids = scalar_value(values['SlimGOIDs']).to_s.split('|').reject{|v| v.empty? }
      names = scalar_value(values['SlimNames']).to_s.split('|')
      ids.each_with_index{|slim_id, i| slim_name_by_id[slim_id] ||= names[i] || slim_id }
    end

    annotations_by_gene = Hash.new{|h,k| h[k] = [] }
    gene_go.through do |ensembl_gene, go_ids|
      next unless identifiers.include?(ensembl_gene)
      gene_names = [identifiers[ensembl_gene]].flatten.compact.collect(&:to_s).reject{|gene| gene.empty? }
      next if gene_names.empty?
      slim_ids = []
      [go_ids].flatten.compact.each do |go_id|
        next unless go_to_slim.include?(go_id)
        slim_ids.concat scalar_value(go_to_slim[go_id]['SlimGOIDs']).to_s.split('|').reject{|v| v.empty? }
      end
      slim_ids = slim_ids.uniq
      next if slim_ids.empty?
      gene_names.each{|gene| annotations_by_gene[gene].concat(slim_ids) }
    end

    tsv = TSV.setup({}, :key_field => 'Associated Gene Name', :fields => %w(SlimGOIDs SlimNames), :type => :list, :namespace => AGS.organism)
    annotations_by_gene.keys.sort.each do |gene|
      slim_ids = annotations_by_gene[gene].uniq.sort
      tsv[gene] = [slim_ids * '|', slim_ids.collect{|slim_id| slim_name_by_id[slim_id] || slim_id } * '|']
    end
    tsv
  end

  dep :full_gene_info
  input :source_type, :select, 'Gene sets to summarize', 'cluster', :select_options => %w(cluster)
  task :goslim_bp_gene_sets => :tsv do |source_type|
    info = step(:full_gene_info).load
    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(Treatment Time Direction SourceType Genes QuerySize), :type => :list, :namespace => AGS.organism)
    id = 0

    result_treatment_order.each do |treatment|
      AGS::TIME_POINTS.each do |time_point|
        %w(up down both).each do |direction|
          genes = []
          info.through do |gene, values|
            values = NamedArray.setup(values, info.fields)
            events = onset_events(values["#{treatment}: FC clusters"])
            selected = events.select{|event| event[:start_time] == time_point }
            selected = selected.select{|event| event[:direction] == direction } unless direction == 'both'
            genes << gene if selected.any?
          end
          id += 1
          genes = genes.uniq.sort
          tsv[id] = [treatment, time_point, direction, source_type, genes * '|', genes.length]
        end
      end
    end
    tsv
  end

  dep :goslim_bp_gene_sets
  dep :gene_goslim_bp_annotations
  dep :expressed_coding_genes
  input :min_query_size, :integer, 'Minimum genes in query set', 10
  input :min_intersection, :integer, 'Minimum genes in GO slim term intersection', 3
  task :goslim_bp_enrichment => :tsv do |min_query_size, min_intersection|
    require 'set'
    gene_sets = step(:goslim_bp_gene_sets).load
    annotations = step(:gene_goslim_bp_annotations).load
    background = step(:expressed_coding_genes).load.collect(&:to_s).to_set
    annotated_background = background & annotations.keys.collect(&:to_s).to_set

    genes_by_slim = Hash.new{|h,k| h[k] = [] }
    slim_name = {}
    annotations.through do |gene, values|
      values = NamedArray.setup(values, annotations.fields)
      next unless annotated_background.include?(gene.to_s)
      ids = scalar_value(values['SlimGOIDs']).to_s.split('|').reject{|v| v.empty? }
      names = scalar_value(values['SlimNames']).to_s.split('|')
      ids.each_with_index do |slim_id, i|
        genes_by_slim[slim_id] << gene.to_s
        slim_name[slim_id] ||= names[i] || slim_id
      end
    end
    genes_by_slim.each{|slim_id, genes| genes.uniq! }

    raw_rows = []
    gene_sets.through do |set_id, values|
      values = NamedArray.setup(values, gene_sets.fields)
      treatment = scalar_value(values['Treatment'])
      time_point = scalar_value(values['Time']).to_i
      direction = scalar_value(values['Direction'])
      source_type = scalar_value(values['SourceType'])
      query_genes = scalar_value(values['Genes']).to_s.split('|').reject{|v| v.empty? }.to_set & annotated_background
      next if query_genes.length < min_query_size

      pvalues = []
      rows = []
      genes_by_slim.keys.sort.each do |slim_id|
        term_genes = genes_by_slim[slim_id].to_set
        intersection = (query_genes & term_genes).to_a.sort
        next if intersection.length < min_intersection
        pvalue = hypergeom_upper_tail(intersection.length, term_genes.length, query_genes.length, annotated_background.length)
        pvalues << pvalue
        rows << [treatment, time_point, direction, source_type, slim_id, slim_name[slim_id], query_genes.length, term_genes.length, intersection.length, pvalue, nil, intersection.length.to_f / query_genes.length, intersection.length.to_f / term_genes.length, intersection * '|']
      end
      qvalues = bh_adjust_values(pvalues)
      rows.each_with_index do |row, row_i|
        row[10] = qvalues[row_i]
        raw_rows << row
      end
    end

    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(Treatment Time Direction SourceType SlimGOID SlimName QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue Precision Recall Genes), :type => :list, :namespace => AGS.organism)
    raw_rows.sort_by{|row| [standard_treatment_sort_index(row[0]), row[1].to_i, row[2].to_s, row[10].to_f, row[5].to_s] }.each_with_index do |row, i|
      tsv[i + 1] = row
    end
    tsv
  end

  dep :goslim_bp_enrichment
  input :adjusted_pvalue_threshold, :float, 'Adjusted p-value threshold', 0.05
  input :top_n_per_context, :integer, 'Top GO slim terms to retain per treatment, time, and direction', 5
  task :goslim_bp_top_terms => :tsv do |adjusted_pvalue_threshold, top_n_per_context|
    enrichment = step(:goslim_bp_enrichment).load
    idx = Hash[enrichment.fields.each_with_index.to_a]
    grouped = Hash.new{|h,k| h[k] = [] }
    enrichment.through do |row_id, values|
      qvalue = numeric_value(values[idx['AdjustedPValue']], 1.0)
      next if qvalue > adjusted_pvalue_threshold
      key = [scalar_value(values[idx['Treatment']]), scalar_value(values[idx['Time']]).to_i, scalar_value(values[idx['Direction']])]
      grouped[key] << values
    end

    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(Treatment Time Direction SlimGOID SlimName QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue Precision Recall Genes), :type => :list, :namespace => AGS.organism)
    id = 0
    grouped.keys.sort_by{|treatment,time,direction| [standard_treatment_sort_index(treatment), time, direction.to_s] }.each do |key|
      grouped[key].sort_by{|values| [numeric_value(values[idx['AdjustedPValue']], 1.0), -numeric_value(values[idx['IntersectionSize']], 0), scalar_value(values[idx['SlimName']]).to_s] }.first(top_n_per_context).each do |values|
        id += 1
        tsv[id] = [
          scalar_value(values[idx['Treatment']]),
          scalar_value(values[idx['Time']]),
          scalar_value(values[idx['Direction']]),
          scalar_value(values[idx['SlimGOID']]),
          scalar_value(values[idx['SlimName']]),
          scalar_value(values[idx['QuerySize']]),
          scalar_value(values[idx['TermBackgroundSize']]),
          scalar_value(values[idx['IntersectionSize']]),
          scalar_value(values[idx['PValue']]),
          scalar_value(values[idx['AdjustedPValue']]),
          scalar_value(values[idx['Precision']]),
          scalar_value(values[idx['Recall']]),
          values[idx['Genes']]
        ]
      end
    end
    tsv
  end

end
