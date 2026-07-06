require 'distribution'
module AGS

  helper :proteomics_prediction_file do |prediction_set|
    case prediction_set.to_s
    when 'full', 'all'
      proteomics_data_dir['ptm_predictions.tsv']
    when 'early'
      proteomics_data_dir['ptm_predictions_early.tsv']
    else
      raise ParameterException, "Unknown prediction set: #{prediction_set}"
    end
  end

  helper :proteomics_prediction_tsv do |prediction_set|
    tsv = proteomics_prediction_file(prediction_set).tsv header_hash: '', type: :double, merge: true
    time_field = proteomics_prediction_time_field tsv.fields
    tsv = tsv.unzip "Treatment", merge: true, delete: false
    tsv = tsv.unzip time_field, delete: false
    tsv
  end

  helper :proteomics_prediction_time_field do |prediction_fields|
    prediction_fields.include?('Time') ? 'Time' : 'MS time'
  end

  helper :proteomics_prediction_module_field do |prediction_fields|
    prediction_fields.include?('Mechanism theme') ? 'Mechanism theme' : 'Mechanism'
  end

  helper :proteomics_prediction_extra_fields do |prediction_fields|
    prediction_fields - ['Protein', 'Residue', 'Treatment', 'Time', 'MS time', 'Prediction', 'Direction']
  end

  helper :proteomics_match_support? do |support|
    %w(directional_match strong_directional_match no_change_match).include?(support.to_s)
  end

  helper :proteomics_miss_support? do |support|
    %w(directional_miss strong_directional_miss unexpected_change).include?(support.to_s)
  end

  helper :proteomics_directional_match_support? do |support|
    %w(directional_match strong_directional_match).include?(support.to_s)
  end

  helper :proteomics_directional_miss_support? do |support|
    %w(directional_miss strong_directional_miss).include?(support.to_s)
  end



  helper :proteomics_binomial_p_greater_half do |matches, trials|
	matches = matches.to_i
	trials  = trials.to_i

	return nil if trials <= 0
	return 1.0 if matches <= 0
	return 0.0 if matches > trials

	1.0 - Distribution::Binomial.cdf(matches - 1, trials, 0.5)
  end

  helper :proteomics_build_ptm_lookup do |ptm_long|
	lookup = {}
	ptm_long.through do |id, values|
	  site = proteomics_scalar(values['Site'])
	  treatment = proteomics_scalar(values['Treatment'])
	  time = proteomics_scalar(values['Time'])
	  time = proteomics_normalize_time(time)
	  lookup[[site, treatment, time]] = NamedArray.setup(values, ptm_long.fields)
	end
	lookup
  end

  helper :proteomics_build_abundance_lookup do |abundance_long|
	lookup = {}
	abundance_long.through do |id, values|
	  gene = proteomics_scalar(values['Gene'])
	  treatment = proteomics_scalar(values['Treatment'])
	  time = proteomics_scalar(values['Time'])
	  time = proteomics_normalize_time(time)
	  lookup[[gene, treatment, time]] = NamedArray.setup(values, abundance_long.fields)
	end
	lookup
  end

  helper :proteomics_prediction_summary_table do |report, support_field='Support'|
	groups = Hash.new{|h,k| h[k] = Hash.new(0) }
	report.through do |id, values|
	  values = NamedArray.setup(values, report.fields)
	  support = proteomics_scalar(values[support_field])
	  prediction = proteomics_float(values['Prediction'])
	  direction_class = prediction.nil? || prediction == 0 ? 'no_change_prediction' : 'directional_prediction'
	  [
		['overall', 'all'],
		['mechanism_theme', proteomics_scalar(values['MechanismTheme'])],
		['confidence', proteomics_scalar(values['Confidence'])],
		['treatment', proteomics_scalar(values['Treatment'])],
		['time', proteomics_scalar(values['Time'])],
		['direction_class', direction_class]
	  ].each do |section, item|
		key = [section, item]
		groups[key]['Total'] += 1
		groups[key][support] += 1
		groups[key]['Matches'] += 1 if proteomics_match_support?(support)
		groups[key]['Misses'] += 1 if proteomics_miss_support?(support)
		groups[key]['DirectionalMatches'] += 1 if proteomics_directional_match_support?(support)
		groups[key]['DirectionalMisses'] += 1 if proteomics_directional_miss_support?(support)
		groups[key]['DirectionalTested'] += 1 if proteomics_directional_match_support?(support) || proteomics_directional_miss_support?(support)
	  end
	end

	fields = %w(Section Item Total Matches Misses MatchFraction DirectionalTested DirectionalMatches DirectionalMisses DirectionalMatchFraction BinomialPAboveHalf NoChangeMatch UnexpectedChange WeakOrNoEffect NotMeasured)
	tsv = TSV.setup({}, 'ID~' + fields * ',')
	id = 0
	groups.keys.sort.each do |section, item|
	  counts = groups[[section, item]]
	  total = counts['Total']
	  matches = counts['Matches']
	  misses = counts['Misses']
	  directional_tested = counts['DirectionalTested']
	  directional_matches = counts['DirectionalMatches']
	  directional_misses = counts['DirectionalMisses']
	  match_fraction = (matches + misses) == 0 ? nil : matches.to_f / (matches + misses)
	  directional_fraction = directional_tested == 0 ? nil : directional_matches.to_f / directional_tested
	  pvalue = proteomics_binomial_p_greater_half(directional_matches, directional_tested)
	  id += 1
	  tsv[id] = [section, item, total, matches, misses, match_fraction, directional_tested, directional_matches, directional_misses, directional_fraction, pvalue, counts['no_change_match'], counts['unexpected_change'], counts['weak_or_no_effect'], counts['not_measured']]
	end
	tsv
  end

  dep :proteomics_ptm_long
  input :prediction_set, :select, 'Prediction file to use', 'full', :select_options => %w(full early)
  input :ptm_source, :select, 'PTM source file', 'collapsed', :select_options => %w(collapsed filtered)
  input :comparison_basis, :select, 'Observed PTM difference used for scoring', 'dmso', :select_options => %w(baseline dmso)
  input :effect_threshold, :float, 'Minimum absolute observed log2 difference to call a measured direction', 0.25
  task :proteomics_prediction_comparison_ptm => :tsv do |prediction_set, ptm_source, comparison_basis, effect_threshold|
	predictions = proteomics_prediction_tsv(prediction_set)
	time_field = proteomics_prediction_time_field(predictions.fields)
	module_field = proteomics_prediction_module_field(predictions.fields)
	ptm_lookup = proteomics_build_ptm_lookup(step(:proteomics_ptm_long).load)

	fields = %w(PredictionSet PTMSource ComparisonBasis Site Protein Residue Treatment Time Prediction Direction Confidence MechanismTheme ObservedDifference ObservedDirection Support DifferenceVsBaseline DifferenceVsDMSO TreatmentMean BaselineMean DMSOMean TreatmentValidValues BaselineValidValues DMSOValidValues PhosphositeMeaning PriorKnowledge PredictionRationale)
	tsv = TSV.setup({}, 'ID~' + fields * ',')
	id = 0


    predictions.through do |key, values|
      site, *_rest = key.split(':')
      values = NamedArray.setup(values, predictions.fields)
      protein = proteomics_scalar(values['Protein']) || proteomics_parse_ptm_key(site).first
      residue = proteomics_scalar(values['Residue']) || proteomics_parse_ptm_key(site)[1]
      treatment = proteomics_normalize_treatment(proteomics_scalar(values['Treatment']))
      time = proteomics_normalize_time(proteomics_scalar(values[time_field]))
      prediction = proteomics_float(values['Prediction'])
      direction = proteomics_scalar(values['Direction'])
      confidence = proteomics_scalar(values['Confidence'])
      mechanism = proteomics_scalar(values[module_field])

      ptm = ptm_lookup[[site, treatment, time]]

      difference_baseline = ptm ? proteomics_float(ptm['DifferenceVsBaseline']) : nil
      difference_dmso = ptm ? proteomics_float(ptm['DifferenceVsDMSO']) : nil
      observed = if comparison_basis.to_s == 'dmso' && treatment != 'DMSO'
                   difference_dmso
                 else
                   difference_baseline
                 end
      support = proteomics_prediction_support(prediction, observed, effect_threshold)
      observed_direction = proteomics_direction_from_sign(proteomics_observed_sign(observed))

      id += 1
      tsv[id] = [prediction_set, ptm_source, comparison_basis, site, protein, residue, treatment, time, prediction, direction, confidence, mechanism, observed, observed_direction, support, difference_baseline, difference_dmso, ptm ? proteomics_float(ptm['TreatmentMean']) : nil, ptm ? proteomics_float(ptm['BaselineMean']) : nil, ptm ? proteomics_float(ptm['DMSOMean']) : nil, ptm ? proteomics_float(ptm['TreatmentValidValues']) : nil, ptm ? proteomics_float(ptm['BaselineValidValues']) : nil, ptm ? proteomics_float(ptm['DMSOValidValues']) : nil, proteomics_scalar(values['Phosphosite meaning']), proteomics_scalar(values['Prior knowledge used']) || proteomics_scalar(values['SIGNOR upstream']), proteomics_scalar(values['Prediction rationale'])]
    end
    tsv
  end

  dep :proteomics_prediction_comparison_ptm
  task :proteomics_prediction_summary_ptm => :tsv do
    proteomics_prediction_summary_table(step(:proteomics_prediction_comparison_ptm).load, 'Support')
  end

  dep :proteomics_ptm_long
  dep :proteomics_abundance_long
  input :prediction_set, :select, 'Prediction file to use', 'full', :select_options => %w(full early)
  input :ptm_source, :select, 'PTM source file', 'collapsed', :select_options => %w(collapsed filtered)
  input :effect_threshold, :float, 'Minimum absolute abundance-adjusted log2 difference to call a measured direction', 0.25
  task :proteomics_prediction_comparison_ptm_abundance => :tsv do |prediction_set, ptm_source, effect_threshold|
    predictions = proteomics_prediction_tsv(prediction_set)
    time_field = proteomics_prediction_time_field(predictions.fields)
    module_field = proteomics_prediction_module_field(predictions.fields)
    ptm_lookup = proteomics_build_ptm_lookup(step(:proteomics_ptm_long).load)
    abundance_lookup = proteomics_build_abundance_lookup(step(:proteomics_abundance_long).load)

    fields = %w(PredictionSet PTMSource Site Protein Residue Treatment Time Prediction Direction Confidence MechanismTheme PTMDifferenceVsDMSO ProteinAbundanceDifference AbundanceAdjustedPTMDifference ObservedDirection Support PTMDifferenceVsBaseline AbundanceQValue AbundanceSignificant PhosphositeMeaning PriorKnowledge PredictionRationale)
    tsv = TSV.setup({}, 'ID~' + fields * ',')
    id = 0

    predictions.through do |key, values|
      site, *_rest = key.split(':')
      values = NamedArray.setup(values, predictions.fields)
      protein = proteomics_scalar(values['Protein']) || proteomics_parse_ptm_key(site).first
      residue = proteomics_scalar(values['Residue']) || proteomics_parse_ptm_key(site)[1]
      treatment = proteomics_normalize_treatment(proteomics_scalar(values['Treatment']))
      time = proteomics_normalize_time(proteomics_scalar(values[time_field]))
      prediction = proteomics_float(values['Prediction'])
      direction = proteomics_scalar(values['Direction'])
      confidence = proteomics_scalar(values['Confidence'])
      mechanism = proteomics_scalar(values[module_field])

      ptm = ptm_lookup[[site, treatment, time]]
      abundance = abundance_lookup[[protein, treatment, time]]

      ptm_baseline = ptm ? proteomics_float(ptm['DifferenceVsBaseline']) : nil
      ptm_dmso = ptm ? proteomics_float(ptm['DifferenceVsDMSO']) : nil
      abundance_difference = abundance ? proteomics_float(abundance['Difference']) : nil
      adjusted = if treatment == 'DMSO'
                   ptm_baseline
                 elsif ptm_dmso.nil? || abundance_difference.nil?
                   nil
                 else
                   ptm_dmso - abundance_difference
                 end
      support = proteomics_prediction_support(prediction, adjusted, effect_threshold)
      observed_direction = proteomics_direction_from_sign(proteomics_observed_sign(adjusted))

      id += 1
      tsv[id] = [prediction_set, ptm_source, site, protein, residue, treatment, time, prediction, direction, confidence, mechanism, ptm_dmso, abundance_difference, adjusted, observed_direction, support, ptm_baseline, abundance ? proteomics_float(abundance['QValue']) : nil, abundance ? proteomics_scalar(abundance['Significant']) : nil, proteomics_scalar(values['Phosphosite meaning']), proteomics_scalar(values['Prior knowledge used']) || proteomics_scalar(values['SIGNOR upstream']), proteomics_scalar(values['Prediction rationale'])]
    end
    tsv
  end

  dep :proteomics_prediction_comparison_ptm_abundance
  task :proteomics_prediction_summary_ptm_abundance => :tsv do
    proteomics_prediction_summary_table(step(:proteomics_prediction_comparison_ptm_abundance).load, 'Support')
  end

  dep :proteomics_prediction_summary_ptm do |jobname,options|
    jobs = []
    ['full', 'early'].each do |prediction_set|
      ['collapsed', 'filtered'].each do |ptm_source|
        ['baseline', 'dmso'].each do |comparison_basis|
          jobs << {prediction_set: prediction_set, ptm_source: ptm_source, comparison_basis: comparison_basis}
        end
      end
    end
    jobs
  end
  dep :proteomics_prediction_summary_ptm_abundance do |jobname,options|
    jobs = []
    ['full', 'early'].each do |prediction_set|
      ['collapsed', 'filtered'].each do |ptm_source|
        ['baseline', 'dmso'].each do |comparison_basis|
          jobs << {prediction_set: prediction_set, ptm_source: ptm_source, comparison_basis: comparison_basis}
        end
      end
    end
    jobs
  end
  task :proteomics_prediction_suite => :array do
    dependencies.each do |dep|
      prediction_set, ptm_source, comparison_basis = dep.recursive_inputs.values_at \
        :prediction_set, :ptm_source, :comparison_basis

      abundance = dep.task_name.to_s.include? 'abundance'

      target = [prediction_set, ptm_source, comparison_basis].compact
      target << 'abundance' if abundance

      Open.cp dep.path, file("ptm_summary.#{target*'.'}.tsv")
      Open.cp dep.dependencies.first.path, file("ptm_comparison.#{target*'.'}.tsv")
    end
    files
  end

end
