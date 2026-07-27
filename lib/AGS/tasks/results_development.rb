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
    when 'column_min_max'
      max = []
      min = []
      fields.each_index do |i|
        col = matrix.values.collect{|row| row[i] }
        max[i] = col.max
        min[i] = col.min.abs
      end
      matrix.each do |tf, row|
        matrix[tf] = row.each_with_index.collect{|v,i| 
          v == 0 ? 0 : (v < 0 ? v / min[i] : v / max[i])
        }
      end
    when 'column_min_max_r'
      max = []
      min = []
      fields.each_index do |i|
        col = matrix.values.collect{|row| row[i] }
        col = col.sort[3..-4]
        max[i] = col.max
        min[i] = col.min.abs
      end
      matrix.each do |tf, row|
        matrix[tf] = row.each_with_index.collect{|v,i| 
          v == 0 ? 0 : (v < 0 ? v / min[i] : v / max[i])
        }
      end
    else
      raise "Unknown normalization #{normalization}"
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


  # Manuscript-ready result tables
  #
  # These tasks reshape the long-form workflow summaries into wide tables that
  # can be pasted directly into the working manuscript. They deliberately read
  # from the primary workflow outputs rather than from hand-copied markdown
  # tables, so they remain synchronized with changes in onset thresholds such as
  # the reverted default 24 h cutoff.

  dep :de_gene_counts_fc0, :fc0_threshold => 0.0, :pvalue_threshold => 0.05
  task :results_table_de_gene_counts_fc0 => :tsv do
    counts = step(:de_gene_counts_fc0).load
    by_key = {}
    counts.through do |id, values|
      values = NamedArray.setup(values, counts.fields)
      next unless scalar_value(values['Direction']) == 'both'
      by_key[[scalar_value(values['Treatment']), scalar_value(values['Time']).to_i]] = scalar_value(values['Genes']).to_i
    end

    fields = AGS::TIME_POINTS.collect{|time| "#{time} h" }
    tsv = TSV.setup({}, :key_field => 'Treatment', :fields => fields, :type => :list)
    result_treatment_order.each do |treatment|
      tsv[treatment] = AGS::TIME_POINTS.collect{|time| by_key[[treatment, time]] || 0 }
    end
    tsv
  end

  dep :onset_first_counts
  task :results_table_onset_first_counts => :tsv do
    counts = step(:onset_first_counts).load
    by_key = {}
    counts.through do |id, values|
      values = NamedArray.setup(values, counts.fields)
      by_key[[scalar_value(values['Treatment']), scalar_value(values['Time']).to_i, scalar_value(values['Direction'])]] = scalar_value(values['Genes']).to_i
    end

    fields = AGS::TIME_POINTS.collect{|time| ["inc #{time} h", "dec #{time} h"] }.flatten + ['Total dynamic']
    tsv = TSV.setup({}, :key_field => 'Treatment', :fields => fields, :type => :list)
    result_treatment_order.each do |treatment|
      row = []
      total = 0
      AGS::TIME_POINTS.each do |time|
        up = by_key[[treatment, time, 'up']] || 0
        down = by_key[[treatment, time, 'down']] || 0
        total += up + down
        row.concat [up, down]
      end
      row << total
      tsv[treatment] = row
    end
    tsv
  end

  dep :onset_episode_counts
  task :results_table_onset_episode_counts => :tsv do
    counts = step(:onset_episode_counts).load
    by_key = {}
    counts.through do |id, values|
      values = NamedArray.setup(values, counts.fields)
      by_key[[scalar_value(values['Treatment']), scalar_value(values['Time']).to_i, scalar_value(values['Direction'])]] = scalar_value(values['Episodes']).to_i
    end

    fields = AGS::TIME_POINTS.collect{|time| ["inc #{time} h", "dec #{time} h"] }.flatten + ['Total episodes']
    tsv = TSV.setup({}, :key_field => 'Treatment', :fields => fields, :type => :list)
    result_treatment_order.each do |treatment|
      row = []
      total = 0
      AGS::TIME_POINTS.each do |time|
        up = by_key[[treatment, time, 'up']] || 0
        down = by_key[[treatment, time, 'down']] || 0
        total += up + down
        row.concat [up, down]
      end
      row << total
      tsv[treatment] = row
    end
    tsv
  end

  dep :tf_activity_call_counts_dynamic
  task :results_table_tf_activity_counts_dynamic => :tsv do
    counts = step(:tf_activity_call_counts_dynamic).load
    by_key = {}
    counts.through do |id, values|
      values = NamedArray.setup(values, counts.fields)
      by_key[[scalar_value(values['Treatment']), scalar_value(values['Time']).to_i, scalar_value(values['Sign'])]] = scalar_value(values['TFActivityCalls']).to_i
    end

    fields = AGS::TIME_POINTS.collect{|time| ["T#{time} positive", "T#{time} negative", "T#{time} total"] }.flatten
    tsv = TSV.setup({}, :key_field => 'Treatment', :fields => fields, :type => :list)
    result_treatment_order.each do |treatment|
      row = []
      AGS::TIME_POINTS.each do |time|
        row << (by_key[[treatment, time, 'positive']] || 0)
        row << (by_key[[treatment, time, 'negative']] || 0)
        row << (by_key[[treatment, time, 'both']] || 0)
      end
      tsv[treatment] = row
    end
    tsv
  end

  dep :sequence_with_changes, :treatment => :placeholder do |jobname, options|
    %w(DMSO FiveZ INT_FiveZ_PI INT_PD_PI PD PI).collect{|treatment| options.merge(:treatment => treatment) }
  end
  task :results_table_sequence_edge_counts => :tsv do
    tsv = TSV.setup({}, :key_field => 'Treatment', :fields => ['Total sequence edges', 'Edges with both TFs self-consistent'], :type => :list)

    dependencies.each do |dep|
      treatment = dep.recursive_inputs[:treatment].to_s
      data = dep.load
      total = 0
      both_self_consistent = 0
      data.through do |id, values|
        values = NamedArray.setup(values, data.fields)
        total += 1
        source_sc = scalar_value(values['Source self-consistent']).to_s == 'true'
        target_sc = scalar_value(values['Target self-consistent']).to_s == 'true'
        both_self_consistent += 1 if source_sc && target_sc
      end
      tsv[treatment] = [total, both_self_consistent]
    end
    tsv
  end

  dep :results_table_de_gene_counts_fc0
  dep :results_table_onset_first_counts
  dep :results_table_onset_episode_counts
  dep :results_table_tf_activity_counts_dynamic
  dep :results_table_sequence_edge_counts
  task :results_manuscript_tables => :array do
    paths = []
    {
      'table_de_gene_counts_fc0.tsv' => step(:results_table_de_gene_counts_fc0),
      'table_onset_first_counts.tsv' => step(:results_table_onset_first_counts),
      'table_onset_episode_counts.tsv' => step(:results_table_onset_episode_counts),
      'table_tf_activity_counts_dynamic.tsv' => step(:results_table_tf_activity_counts_dynamic),
      'table_sequence_edge_counts.tsv' => step(:results_table_sequence_edge_counts)
    }.each do |filename, dep|
      path = file(filename)
      Open.write(path, dep.load.to_s)
      paths << path
    end
    paths
  end


  # Dynamic versus non-dynamic TF activity timing comparisons
  #
  # These tasks support the manuscript claim that the dynamic-onset strategy
  # improves regulatory chronology. The unit of comparison is a
  # TF-treatment-sign event, where sign is positive or negative inferred TF
  # activity. For each scheme we record the first sampled window where the event
  # appears and the number of sampled windows where it is present.

  dep :tf_predictions, :scheme => 'dynamic'
  dep :tf_predictions, :scheme => 'non-dynamic'
  task :dynamic_vs_nondynamic_tf_timing_events => :tsv do
    scheme_data = {}
    dependencies.each do |dep|
      scheme = dep.recursive_inputs[:scheme].to_s
      predictions = dep.load
      first = {}
      persistence = Hash.new(0)

      predictions.through do |tf, values|
        values = NamedArray.setup(values, predictions.fields)
        result_treatment_order.each do |treatment|
          AGS::TIME_POINTS.each do |time_point|
            field = "#{treatment}-T#{time_point}"
            next unless predictions.fields.include?(field)
            activity = numeric_value(values[field], 0.0)
            next if activity == 0
            sign = activity > 0 ? 'positive' : 'negative'
            key = [tf, treatment, sign]
            first[key] ||= time_point
            persistence[key] += 1
          end
        end
      end
      scheme_data[scheme] = {:first => first, :persistence => persistence}
    end

    dynamic = scheme_data['dynamic'] || {:first => {}, :persistence => {}}
    nondynamic = scheme_data['non-dynamic'] || {:first => {}, :persistence => {}}
    keys = (dynamic[:first].keys + nondynamic[:first].keys).uniq
    time_index = Hash[AGS::TIME_POINTS.each_with_index.to_a]

    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(TF Treatment Sign DynamicFirst NonDynamicFirst Category DynamicPersistence NonDynamicPersistence), :type => :list, :namespace => AGS.organism)
    keys.sort_by{|tf,treatment,sign| [standard_treatment_sort_index(treatment), sign, tf] }.each_with_index do |key, i|
      tf, treatment, sign = key
      d_first = dynamic[:first][key]
      n_first = nondynamic[:first][key]
      category = if d_first && n_first
                   if time_index[d_first] < time_index[n_first]
                     'dynamic earlier'
                   elsif time_index[d_first] > time_index[n_first]
                     'non-dynamic earlier'
                   else
                     'same time'
                   end
                 elsif d_first
                   'dynamic only'
                 else
                   'non-dynamic only'
                 end
      tsv[i + 1] = [tf, treatment, sign, d_first, n_first, category, dynamic[:persistence][key] || 0, nondynamic[:persistence][key] || 0]
    end
    tsv
  end

  dep :dynamic_vs_nondynamic_tf_timing_events
  task :dynamic_vs_nondynamic_tf_timing_summary => :tsv do
    events = step(:dynamic_vs_nondynamic_tf_timing_events).load
    counts = Hash.new(0)
    category_order = ['dynamic earlier', 'same time', 'non-dynamic earlier', 'dynamic only', 'non-dynamic only']
    events.through do |id, values|
      values = NamedArray.setup(values, events.fields)
      treatment = scalar_value(values['Treatment'])
      sign = scalar_value(values['Sign'])
      category = scalar_value(values['Category'])
      counts[[treatment, sign, category]] += 1
      counts[[treatment, 'both', category]] += 1
    end

    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(Treatment Sign Category TFEvents SharedCategory), :type => :list, :namespace => AGS.organism)
    id = 0
    result_treatment_order.each do |treatment|
      %w(positive negative both).each do |sign|
        category_order.each do |category|
          id += 1
          shared = ['dynamic earlier', 'same time', 'non-dynamic earlier'].include?(category)
          tsv[id] = [treatment, sign, category, counts[[treatment, sign, category]] || 0, shared]
        end
      end
    end
    tsv
  end

  dep :dynamic_vs_nondynamic_tf_timing_events
  task :dynamic_vs_nondynamic_tf_persistence_distribution => :tsv do
    events = step(:dynamic_vs_nondynamic_tf_timing_events).load
    counts = Hash.new(0)
    events.through do |id, values|
      values = NamedArray.setup(values, events.fields)
      treatment = scalar_value(values['Treatment'])
      sign = scalar_value(values['Sign'])
      d_persistence = scalar_value(values['DynamicPersistence']).to_i
      n_persistence = scalar_value(values['NonDynamicPersistence']).to_i
      if d_persistence > 0
        counts[['dynamic', treatment, sign, d_persistence]] += 1
        counts[['dynamic', treatment, 'both', d_persistence]] += 1
      end
      if n_persistence > 0
        counts[['non-dynamic', treatment, sign, n_persistence]] += 1
        counts[['non-dynamic', treatment, 'both', n_persistence]] += 1
      end
    end

    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(Scheme Treatment Sign ActiveTimepoints TFEvents), :type => :list, :namespace => AGS.organism)
    id = 0
    %w(dynamic non-dynamic).each do |scheme|
      result_treatment_order.each do |treatment|
        %w(positive negative both).each do |sign|
          (1..AGS::TIME_POINTS.length).each do |active_timepoints|
            id += 1
            tsv[id] = [scheme, treatment, sign, active_timepoints, counts[[scheme, treatment, sign, active_timepoints]] || 0]
          end
        end
      end
    end
    tsv
  end

  dep :dynamic_vs_nondynamic_tf_timing_events
  task :dynamic_vs_nondynamic_tf_persistence_comparison => :tsv do
    events = step(:dynamic_vs_nondynamic_tf_timing_events).load
    counts = Hash.new(0)
    events.through do |id, values|
      values = NamedArray.setup(values, events.fields)
      d_first = scalar_value(values['DynamicFirst'])
      n_first = scalar_value(values['NonDynamicFirst'])
      next if d_first.nil? || d_first.to_s.empty? || n_first.nil? || n_first.to_s.empty?
      treatment = scalar_value(values['Treatment'])
      sign = scalar_value(values['Sign'])
      d_persistence = scalar_value(values['DynamicPersistence']).to_i
      n_persistence = scalar_value(values['NonDynamicPersistence']).to_i
      category = if d_persistence < n_persistence
                   'dynamic shorter'
                 elsif d_persistence > n_persistence
                   'non-dynamic shorter'
                 else
                   'same persistence'
                 end
      counts[[treatment, sign, category]] += 1
      counts[[treatment, 'both', category]] += 1
    end

    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(Treatment Sign Category TFEvents), :type => :list, :namespace => AGS.organism)
    id = 0
    result_treatment_order.each do |treatment|
      %w(positive negative both).each do |sign|
        ['dynamic shorter', 'same persistence', 'non-dynamic shorter'].each do |category|
          id += 1
          tsv[id] = [treatment, sign, category, counts[[treatment, sign, category]] || 0]
        end
      end
    end
    tsv
  end


  dep :tf_predictions
  task :dmso_reference_tf_activity_similarity => :tsv do
    predictions = step(:tf_predictions).load
    refs = ['DMSO-T8', 'DMSO-T24']
    inhibitor_treatments = result_treatment_order.reject{|treatment| treatment == 'DMSO' }

    sign = lambda do |value|
      value = value.to_f
      value > 0 ? 1 : (value < 0 ? -1 : 0)
    end

    pearson = lambda do |a, b|
      n = a.length
      return nil if n == 0
      ma = a.inject(0.0, &:+) / n
      mb = b.inject(0.0, &:+) / n
      va = a.collect{|v| (v - ma) ** 2 }.inject(0.0, &:+)
      vb = b.collect{|v| (v - mb) ** 2 }.inject(0.0, &:+)
      return nil if va == 0 || vb == 0
      cov = a.zip(b).collect{|x,y| (x - ma) * (y - mb) }.inject(0.0, &:+)
      cov / Math.sqrt(va * vb)
    end

    cosine = lambda do |a, b|
      na = Math.sqrt(a.collect{|v| v ** 2 }.inject(0.0, &:+))
      nb = Math.sqrt(b.collect{|v| v ** 2 }.inject(0.0, &:+))
      return nil if na == 0 || nb == 0
      a.zip(b).collect{|x,y| x * y }.inject(0.0, &:+) / (na * nb)
    end

    fields = %w(Reference InhibitorT24 ReferenceActive InhibitorActive SharedActive UnionActive Jaccard SharedSameSign SharedOppositeSign SharedSignConcordance InhibitorActiveAlsoReferenceFraction ReferenceActiveReappearsFraction PearsonAllTFs CosineAllTFs PearsonUnionActive CosineUnionActive)
    tsv = TSV.setup({}, :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism)
    id = 0

    refs.each do |ref_field|
      next unless predictions.fields.include?(ref_field)
      inhibitor_treatments.each do |treatment|
        inhib_field = "#{treatment}-T24"
        next unless predictions.fields.include?(inhib_field)

        ref_values = {}
        inhib_values = {}
        predictions.through do |tf, values|
          values = NamedArray.setup(values, predictions.fields)
          ref_values[tf] = numeric_value(values[ref_field], 0.0)
          inhib_values[tf] = numeric_value(values[inhib_field], 0.0)
        end

        ref_active = ref_values.keys.select{|tf| ref_values[tf] != 0 }
        inhib_active = inhib_values.keys.select{|tf| inhib_values[tf] != 0 }
        shared = ref_active & inhib_active
        union = ref_active | inhib_active
        same = shared.count{|tf| sign.call(ref_values[tf]) == sign.call(inhib_values[tf]) }
        opposite = shared.count{|tf| sign.call(ref_values[tf]) == - sign.call(inhib_values[tf]) }

        all_tfs = predictions.keys
        ref_all = all_tfs.collect{|tf| ref_values[tf] || 0.0 }
        inhib_all = all_tfs.collect{|tf| inhib_values[tf] || 0.0 }
        ref_union = union.collect{|tf| ref_values[tf] || 0.0 }
        inhib_union = union.collect{|tf| inhib_values[tf] || 0.0 }

        id += 1
        tsv[id] = [ref_field, inhib_field, ref_active.length, inhib_active.length, shared.length, union.length,
                   union.empty? ? nil : shared.length.to_f / union.length,
                   same, opposite,
                   shared.empty? ? nil : same.to_f / shared.length,
                   inhib_active.empty? ? nil : shared.length.to_f / inhib_active.length,
                   ref_active.empty? ? nil : shared.length.to_f / ref_active.length,
                   pearson.call(ref_all, inhib_all), cosine.call(ref_all, inhib_all),
                   union.length > 1 ? pearson.call(ref_union, inhib_union) : nil,
                   union.length > 1 ? cosine.call(ref_union, inhib_union) : nil]
      end
    end
    tsv
  end

  dep :tf_predictions
  task :dmso_reference_tf_activity_overlap_details => :tsv do
    predictions = step(:tf_predictions).load
    refs = ['DMSO-T8', 'DMSO-T24']
    inhibitor_treatments = result_treatment_order.reject{|treatment| treatment == 'DMSO' }

    sign = lambda do |value|
      value = value.to_f
      value > 0 ? 'positive' : (value < 0 ? 'negative' : 'inactive')
    end

    fields = %w(Reference InhibitorT24 TF ReferenceActivity InhibitorActivity ReferenceSign InhibitorSign Category)
    dumper = TSV::Dumper.new :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism
    dumper.init
    id = 0

    TSV.traverse predictions, :into => dumper, :bar => self.progress_bar('Comparing DMSO reference and inhibitor T24 TF scoreboards') do |tf, values|
      values = NamedArray.setup(values, predictions.fields)
      res = []
      refs.each do |ref_field|
        next unless predictions.fields.include?(ref_field)
        ref_activity = numeric_value(values[ref_field], 0.0)
        ref_sign = sign.call(ref_activity)
        inhibitor_treatments.each do |treatment|
          inhib_field = "#{treatment}-T24"
          next unless predictions.fields.include?(inhib_field)
          inhib_activity = numeric_value(values[inhib_field], 0.0)
          inhib_sign = sign.call(inhib_activity)
          next if ref_sign == 'inactive' && inhib_sign == 'inactive'
          category = if ref_sign != 'inactive' && inhib_sign != 'inactive'
                       ref_sign == inhib_sign ? 'shared_same_sign' : 'shared_opposite_sign'
                     elsif ref_sign != 'inactive'
                       'reference_only'
                     else
                       'inhibitor_only'
                     end
          id += 1
          res << [id, [ref_field, inhib_field, tf, ref_activity, inhib_activity, ref_sign, inhib_sign, category]]
        end
      end
      res.extend MultipleResult
      res
    end
  end


end
