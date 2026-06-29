module AGS

  PROTEOMICS_TIMES = %w(0.5 2 8) unless const_defined?(:PROTEOMICS_TIMES)
  PROTEOMICS_TREATMENTS = %w(DMSO PI PD FiveZ INT_PD_PI INT_FiveZ_PI) unless const_defined?(:PROTEOMICS_TREATMENTS)
  PROTEOMICS_DRUG_TREATMENTS = %w(PI PD FiveZ INT_PD_PI INT_FiveZ_PI) unless const_defined?(:PROTEOMICS_DRUG_TREATMENTS)

  helper :proteomics_data_dir do
    Path.setup('data/proteomics')
  end

  helper :proteomics_float do |value|
    value = value.first if Array === value
    return nil if value.nil?
    value = value.to_s
    return nil if value.empty? || value == 'NA' || value == 'NaN'
    value.to_f
  end

  helper :proteomics_scalar do |value|
    value = value.first if Array === value
    value.nil? ? nil : value.to_s
  end

  helper :proteomics_mean do |values|
    values = values.collect{|v| proteomics_float(v) }.compact
    values.empty? ? nil : values.inject(0.0){|acc,v| acc + v } / values.length
  end

  helper :proteomics_normalize_treatment do |treatment|
    case treatment.to_s
    when '5Z', 'FiveZ', 'X5Z'
      'FiveZ'
    when 'PIPD', 'INT_PD_PI', 'XPIPD'
      'INT_PD_PI'
    when 'PI5Z', 'INT_FiveZ_PI', 'XPI5Z'
      'INT_FiveZ_PI'
    when 'CTRL', 'CONTROL', 'DMSO', 'XDMSO'
      'DMSO'
    else
      treatment.to_s.sub(/^X/, '')
    end
  end

  helper :proteomics_raw_treatment_code do |treatment|
    case proteomics_normalize_treatment(treatment)
    when 'FiveZ'
      '5Z'
    when 'INT_PD_PI'
      'PIPD'
    when 'INT_FiveZ_PI'
      'PI5Z'
    when 'DMSO'
      'DMSO'
    else
      treatment.to_s
    end
  end

  helper :proteomics_normalize_time do |time|
    time.to_s.sub('0-5', '0.5').sub('0_5', '0.5').sub(/h$/, '')
  end

  helper :proteomics_time_label do |time|
    proteomics_normalize_time(time) + 'h'
  end

  helper :proteomics_prediction_sign do |prediction|
    p = proteomics_float(prediction)
    return 0 if p.nil? || p == 0
    p > 0 ? 1 : -1
  end

  helper :proteomics_observed_sign do |value|
    v = proteomics_float(value)
    return 0 if v.nil? || v == 0
    v > 0 ? 1 : -1
  end

  helper :proteomics_direction_from_sign do |sign|
    case sign.to_i
    when 1
      'up'
    when -1
      'down'
    else
      'no_change'
    end
  end

  helper :proteomics_prediction_support do |prediction, observed, effect_threshold=0.25|
    observed = proteomics_float(observed)
    return 'not_measured' if observed.nil?
    pred_sign = proteomics_prediction_sign(prediction)
    obs_sign = proteomics_observed_sign(observed)

    if pred_sign == 0
      observed.abs < effect_threshold ? 'no_change_match' : 'unexpected_change'
    else
      if observed.abs < effect_threshold
        'weak_or_no_effect'
      elsif pred_sign == obs_sign
        observed.abs >= 0.5 ? 'strong_directional_match' : 'directional_match'
      else
        observed.abs >= 0.5 ? 'strong_directional_miss' : 'directional_miss'
      end
    end
  end

  helper :proteomics_parse_abundance_filename do |path|
    basename = path.basename.to_s.sub(/\.tsv$/, '')
    if basename =~ /^(.+)-T(.+)$/
      treatment = proteomics_normalize_treatment($1)
      time = proteomics_normalize_time($2)
      [treatment, time]
    else
      nil
    end
  end

  helper :proteomics_abundance_paths do
    proteomics_data_dir['abundance'].glob('*.tsv').sort
  end

  helper :proteomics_select_field do |fields, pattern|
    fields.find{|field| field.to_s =~ pattern }
  end

  helper :proteomics_abundance_replicate_fields do |fields, raw_code, dmso=false|
    if dmso
      fields.select{|field| field.to_s =~ /^XDMSO[\-_]/ && field.to_s =~ /Rep\d+/ }
    else
      fields.select{|field| field.to_s =~ /^X#{Regexp.escape(raw_code)}[\-_]/ && field.to_s =~ /Rep\d+/ }
    end
  end

  task :proteomics_abundance_file_index => :tsv do
    tsv = TSV.setup({}, 'ID~Treatment,Time,File,Rows,TreatmentReplicates,DMSOReplicates,DifferenceField,QValueField')
    id = 0
    proteomics_abundance_paths.each do |path|
      treatment, time = proteomics_parse_abundance_filename(path)
      raw_code = proteomics_raw_treatment_code(treatment)
      table = path.tsv :type => :list
      treatment_fields = proteomics_abundance_replicate_fields(table.fields, raw_code, false)
      dmso_fields = proteomics_abundance_replicate_fields(table.fields, raw_code, true)
      diff_field = proteomics_select_field(table.fields, /Student's T-test Difference/)
      q_field = proteomics_select_field(table.fields, /Student's T-test q-value/)
      id += 1
      tsv[id] = [treatment, time, path.to_s, table.keys.length, treatment_fields.length, dmso_fields.length, diff_field, q_field]
    end
    tsv
  end

  task :proteomics_abundance_long => :tsv do
    fields = %w(Gene Treatment Time TreatmentMean DMSOMean Difference QValue NegLogP TestStatistic Significant TreatmentValidValues DMSOValidValues)
    tsv = TSV.setup({}, 'ID~' + fields * ',')
    id = 0

    proteomics_abundance_paths.each do |path|
      treatment, time = proteomics_parse_abundance_filename(path)
      raw_code = proteomics_raw_treatment_code(treatment)
      table = path.tsv :type => :list
      diff_field = proteomics_select_field(table.fields, /Student's T-test Difference/)
      q_field = proteomics_select_field(table.fields, /Student's T-test q-value/)
      neglogp_field = proteomics_select_field(table.fields, /-Log Student's T-test p-value/)
      statistic_field = proteomics_select_field(table.fields, /Student's T-test Test statistic/)
      significant_field = proteomics_select_field(table.fields, /Student's T-test Significant/)
      treatment_valid_field = proteomics_select_field(table.fields, /Valid values X#{Regexp.escape(raw_code)}/)
      dmso_valid_field = proteomics_select_field(table.fields, /Valid values XDMSO/)
      treatment_fields = proteomics_abundance_replicate_fields(table.fields, raw_code, false)
      dmso_fields = proteomics_abundance_replicate_fields(table.fields, raw_code, true)

      table.through bar: self.progress_bar(path) do |gene, values|
        values = NamedArray.setup(values, table.fields)
        id += 1
        treatment_mean = proteomics_mean(treatment_fields.collect{|f| values[f] })
        dmso_mean = proteomics_mean(dmso_fields.collect{|f| values[f] })
        difference = proteomics_float(values[diff_field])
        qvalue = proteomics_float(values[q_field])
        neglogp = proteomics_float(values[neglogp_field])
        statistic = proteomics_float(values[statistic_field])
        significant_value = significant_field ? values[significant_field].to_s : ''
        significant = significant_value == '+' ? 'true' : 'false'
        treatment_valid = treatment_valid_field ? proteomics_float(values[treatment_valid_field]) : nil
        dmso_valid = dmso_valid_field ? proteomics_float(values[dmso_valid_field]) : nil
        tsv[id] = [gene, treatment, time, treatment_mean, dmso_mean, difference, qvalue, neglogp, statistic, significant, treatment_valid, dmso_valid]
      end
    end

    tsv
  end

  dep :proteomics_abundance_long
  task :proteomics_abundance_matrix => :tsv do
    long = step(:proteomics_abundance_long).load
    context_fields = PROTEOMICS_DRUG_TREATMENTS.collect{|treatment| PROTEOMICS_TIMES.collect{|time| "#{treatment}-T#{time}" } }.flatten
    data = {}
    long.through do |id, values|
      gene = proteomics_scalar(values['Gene'])
      treatment = proteomics_scalar(values['Treatment'])
      time = proteomics_scalar(values['Time'])
      field = "#{treatment}-T#{time}"
      next unless context_fields.include?(field)
      data[gene] ||= Array.new(context_fields.length)
      data[gene][context_fields.index(field)] = proteomics_float(values['Difference'])
    end
    TSV.setup(data, :key_field => 'Associated Gene Name', :fields => context_fields, :type => :list, :namespace => AGS.organism)
  end

  dep :proteomics_abundance_long
  input :qvalue_threshold, :float, 'q-value threshold', 0.05
  input :abs_difference_threshold, :float, 'Absolute log2 protein abundance difference threshold', 0.0
  task :proteomics_abundance_calls => :tsv do |qvalue_threshold, abs_difference_threshold|
    long = step(:proteomics_abundance_long).load
    tsv = TSV.setup({}, 'ID~Gene,Treatment,Time,Difference,QValue,Direction,Significant')
    new_id = 0
    long.through do |id, values|
      difference = proteomics_float(values['Difference'])
      qvalue = proteomics_float(values['QValue'])
      next if difference.nil? || qvalue.nil?
      next unless qvalue < qvalue_threshold && difference.abs >= abs_difference_threshold
      direction = difference > 0 ? 'up' : 'down'
      new_id += 1
      tsv[new_id] = [proteomics_scalar(values['Gene']), proteomics_scalar(values['Treatment']), proteomics_scalar(values['Time']), difference, qvalue, direction, proteomics_scalar(values['Significant'])]
    end
    tsv
  end

  dep :proteomics_abundance_calls
  task :proteomics_abundance_call_counts => :tsv do
    calls = step(:proteomics_abundance_calls).load
    counts = Hash.new(0)
    calls.through do |id, values|
      treatment = proteomics_scalar(values['Treatment'])
      time = proteomics_scalar(values['Time'])
      direction = proteomics_scalar(values['Direction'])
      counts[[treatment, time, direction]] += 1
      counts[[treatment, time, 'both']] += 1
    end

    tsv = TSV.setup({}, 'ID~Treatment,Time,Direction,Proteins')
    id = 0
    PROTEOMICS_DRUG_TREATMENTS.each do |treatment|
      PROTEOMICS_TIMES.each do |time|
        %w(up down both).each do |direction|
          id += 1
          tsv[id] = [treatment, time, direction, counts[[treatment, time, direction]]]
        end
      end
    end
    tsv
  end
end
