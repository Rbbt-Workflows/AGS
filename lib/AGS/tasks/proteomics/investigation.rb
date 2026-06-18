module AGS

  PROTEOMICS_MECHANISM_EXPECTATIONS = [
    ['FOXO_CDk_inhibitor_checkpoint', 'PI', '8', 'CDKN1B', 'up', 'FOXO release by PI3K/AKT inhibition should increase CDKN1B/p27'],
    ['FOXO_CDk_inhibitor_checkpoint', 'INT_PD_PI', '8', 'CDKN1B', 'up', 'INT_PD_PI should combine FOXO release with MEK/AP1 repression'],
    ['FOXO_CDk_inhibitor_checkpoint', 'INT_FiveZ_PI', '8', 'CDKN1B', 'up', 'PI-containing combination should increase CDKN1B/p27'],
    ['FOXO_CDk_inhibitor_checkpoint', 'PD', '8', 'CDKN1B', 'up', 'Delayed cell-cycle arrest may also increase CDKN1B'],
    ['FOXO_CDk_inhibitor_checkpoint', 'FiveZ', '8', 'CDKN1B', 'up', 'Delayed cell-cycle arrest may also increase CDKN1B'],

    ['Cyclin_D_suppression', 'PI', '8', 'CCND1', 'down', 'FOXO and ERK-output repression should reduce cyclin D1'],
    ['Cyclin_D_suppression', 'PD', '8', 'CCND1', 'down', 'MEK/ERK/AP1 repression should reduce cyclin D1'],
    ['Cyclin_D_suppression', 'FiveZ', '8', 'CCND1', 'down', 'FiveZ early PD-like AP1 repression predicts cyclin D1 loss'],
    ['Cyclin_D_suppression', 'INT_PD_PI', '8', 'CCND1', 'down', 'Combination should strongly reduce cyclin D1'],
    ['Cyclin_D_suppression', 'INT_FiveZ_PI', '8', 'CCND1', 'down', 'Combination should reduce cyclin D1 before possible late rebound'],

    ['AP1_output_suppression', 'PD', '8', 'JUN', 'down', 'MEK/ERK inhibition should reduce AP1 output'],
    ['AP1_output_suppression', 'FiveZ', '8', 'JUN', 'down', 'FiveZ PD-like early state predicts reduced JUN protein at 8 h'],
    ['AP1_output_suppression', 'INT_PD_PI', '8', 'JUN', 'down', 'INT_PD_PI should suppress AP1 rescue'],
    ['AP1_output_suppression', 'INT_FiveZ_PI', '8', 'JUN', 'down', 'INT_FiveZ_PI should suppress AP1 before late rebound'],

    ['MYC_E2F_cell_cycle_collapse', 'PI', '8', 'E2F3', 'down', 'E2F/cell-cycle repression should lower E2F3 protein'],
    ['MYC_E2F_cell_cycle_collapse', 'PD', '8', 'E2F3', 'down', 'Delayed E2F collapse expected at 8 h'],
    ['MYC_E2F_cell_cycle_collapse', 'FiveZ', '8', 'E2F3', 'down', 'Delayed E2F collapse expected at 8 h'],
    ['MYC_E2F_cell_cycle_collapse', 'INT_PD_PI', '8', 'E2F3', 'down', 'INT_PD_PI should strongly collapse E2F/cell-cycle module'],
    ['MYC_E2F_cell_cycle_collapse', 'INT_FiveZ_PI', '8', 'E2F3', 'down', 'INT_FiveZ_PI should suppress E2F before late rebound'],

    ['CDK_inhibitor_checkpoint_extension', 'PD', '8', 'CDKN1C', 'up', 'PD late growth arrest may increase CDK inhibitor CDKN1C'],
    ['CDK_inhibitor_checkpoint_extension', 'INT_PD_PI', '8', 'CDKN1C', 'up', 'INT_PD_PI strong arrest should increase CDKN1C'],
    ['CDK_inhibitor_checkpoint_extension', 'INT_FiveZ_PI', '8', 'CDKN2B', 'up', 'INT_FiveZ_PI should show CDK-inhibitor/checkpoint support'],
    ['CDK_inhibitor_checkpoint_extension', 'PD', '8', 'CDKN2B', 'up', 'PD should show CDK-inhibitor/checkpoint support'],

    ['Late_stress_inflammatory_axis', 'INT_PD_PI', '8', 'NFKB2', 'up', 'Late INT_PD_PI stress/remodeling may engage NF-kB2'],
    ['Late_stress_inflammatory_axis', 'INT_PD_PI', '8', 'IRF3', 'up', 'IRF/interferon module predicted transcriptionally; protein abundance may not follow'],
    ['Late_stress_inflammatory_axis', 'INT_PD_PI', '8', 'CEBPB', 'up', 'CEBPB transcriptional activity/onset predicts stress response, but abundance may not follow'],
    ['Late_stress_inflammatory_axis', 'INT_FiveZ_PI', '8', 'NFKB2', 'up', 'NF-kB2 may participate in FiveZ_PI late rebound, but 8 h may be pre-rebound'],
  ] unless const_defined?(:PROTEOMICS_MECHANISM_EXPECTATIONS)

  helper :proteomics_support_label do |expected_direction, difference, qvalue|
    return 'not_measured' if difference.nil?
    observed_direction = difference > 0 ? 'up' : (difference < 0 ? 'down' : 'zero')
    match = observed_direction == expected_direction.to_s
    if match && qvalue && qvalue < 0.05
      'strong_support'
    elsif match && qvalue && qvalue < 0.2
      'moderate_support'
    elsif match
      'directional_support_not_significant'
    elsif qvalue && qvalue < 0.05
      'significant_opposite'
    else
      'no_support_or_opposite'
    end
  end

  helper :proteomics_observed_direction do |difference|
    return nil if difference.nil?
    difference > 0 ? 'up' : (difference < 0 ? 'down' : 'zero')
  end

  dep :proteomics_abundance_long
  dep :full_gene_info
  task :proteomics_candidate_marker_report => :tsv do
    long = step(:proteomics_abundance_long).load
    info = step(:full_gene_info).load

    lookup = {}
    long.through do |id, values|
      lookup[[proteomics_scalar(values['Gene']), proteomics_scalar(values['Treatment']), proteomics_scalar(values['Time'])]] = values
    end

    tsv = TSV.setup({}, 'ID~Mechanism,Treatment,Time,Gene,ExpectedDirection,ObservedDifference,QValue,ObservedDirection,Support,MRNAOnset,Rationale')
    id = 0
    PROTEOMICS_MECHANISM_EXPECTATIONS.each do |mechanism, treatment, time, gene, expected, rationale|
      values = lookup[[gene, treatment, time]]
      difference = values ? proteomics_float(values['Difference']) : nil
      qvalue = values ? proteomics_float(values['QValue']) : nil
      observed = proteomics_observed_direction(difference)
      support = proteomics_support_label(expected, difference, qvalue)
      onset = nil
      if info.include?(gene)
        row = NamedArray.setup(info[gene], info.fields)
        onset_value = row["#{treatment}: FC clusters"]
        onset = [onset_value].flatten.compact.collect(&:to_s).reject{|v| v.empty? } * '|'
      end
      id += 1
      tsv[id] = [mechanism, treatment, time, gene, expected, difference, qvalue, observed, support, onset, rationale]
    end
    tsv
  end

  dep :proteomics_abundance_long
  task :proteomics_targeted_marker_matrix => :tsv do
    markers = PROTEOMICS_MECHANISM_EXPECTATIONS.collect{|m| m[3] }.uniq
    long = step(:proteomics_abundance_long).load
    fields = PROTEOMICS_TREATMENTS.collect{|treatment| PROTEOMICS_TIMES.collect{|time| "#{treatment}-T#{time}" } }.flatten
    data = {}
    long.through do |id, values|
      gene = proteomics_scalar(values['Gene'])
      next unless markers.include?(gene)
      treatment = proteomics_scalar(values['Treatment'])
      time = proteomics_scalar(values['Time'])
      field = "#{treatment}-T#{time}"
      next unless fields.include?(field)
      data[gene] ||= Array.new(fields.length)
      data[gene][fields.index(field)] = proteomics_float(values['Difference'])
    end
    TSV.setup(data, :key_field => 'Associated Gene Name', :fields => fields, :type => :list, :namespace => AGS.organism)
  end

  dep :proteomics_candidate_marker_report
  dep :proteomics_abundance_call_counts
  dep :proteomics_ptm_status
  task :proteomics_mechanism_summary => :tsv do
    report = step(:proteomics_candidate_marker_report).load
    counts = step(:proteomics_abundance_call_counts).load
    ptm = step(:proteomics_ptm_status).load

    support_counts = Hash.new(0)
    report.through do |id, values|
      support_counts[[proteomics_scalar(values['Mechanism']), proteomics_scalar(values['Support'])]] += 1
    end

    tsv = TSV.setup({}, 'ID~Section,Item,Value,Comment')
    id = 0
    support_counts.keys.collect{|k| k[0] }.uniq.sort.each do |mechanism|
      %w(strong_support moderate_support directional_support_not_significant no_support_or_opposite significant_opposite not_measured).each do |support|
        value = support_counts[[mechanism, support]]
        next if value == 0
        id += 1
        tsv[id] = ['mechanism_support', mechanism + ':' + support, value, 'Counts over targeted protein markers']
      end
    end

    counts.through do |cid, values|
      next unless proteomics_scalar(values['Direction']) == 'both'
      id += 1
      tsv[id] = ['global_abundance_calls', "#{proteomics_scalar(values['Treatment'])}-T#{proteomics_scalar(values['Time'])}", proteomics_scalar(values['Proteins']), 'Proteins with q < 0.05']
    end

    ptm.through do |pid, values|
      id += 1
      tsv[id] = ['ptm_status', File.basename(proteomics_scalar(values['File'])), proteomics_scalar(values['Status']), "Data rows: #{proteomics_scalar(values['DataRows'])}"]
    end

    tsv
  end
end
