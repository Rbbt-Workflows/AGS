module AGS

  PROTEOMICS_PTM_EXPECTATIONS = [
    ['PI3K_AKT_FOXO_release', 'PI', '0.5', 'FOXO3_T32_M1', 'down', 'AKT inhibitory phosphorylation site; PI should reduce FOXO3 T32 phosphorylation'],
    ['PI3K_AKT_FOXO_release', 'PI', '2', 'FOXO3_T32_M1', 'down', 'AKT inhibitory phosphorylation site; PI should reduce FOXO3 T32 phosphorylation'],
    ['PI3K_AKT_FOXO_release', 'INT_PD_PI', '0.5', 'FOXO3_T32_M1', 'down', 'PI-containing combination should reduce inhibitory FOXO3 phosphorylation'],
    ['PI3K_AKT_FOXO_release', 'INT_PD_PI', '2', 'FOXO3_T32_M1', 'down', 'PI-containing combination should reduce inhibitory FOXO3 phosphorylation'],
    ['PI3K_AKT_FOXO_release', 'INT_FiveZ_PI', '0.5', 'FOXO3_T32_M1', 'down', 'PI-containing combination should reduce inhibitory FOXO3 phosphorylation'],

    ['GSK3B_release_from_AKT', 'PI', '0.5', 'GSK3B_S9_M1', 'down', 'AKT inhibitory GSK3B S9 phosphorylation should decrease if PI reduces AKT activity'],
    ['GSK3B_release_from_AKT', 'INT_PD_PI', '0.5', 'GSK3B_S9_M1', 'down', 'PI-containing combination should decrease inhibitory GSK3B S9 phosphorylation'],
    ['GSK3B_release_from_AKT', 'INT_FiveZ_PI', '0.5', 'GSK3B_S9_M1', 'down', 'PI-containing combination should decrease inhibitory GSK3B S9 phosphorylation'],

    ['MEK_ERK_output_suppression', 'PD', '0.5', 'MAPK1_Y187_M1', 'down', 'PD should reduce ERK2/MAPK1 activating phosphorylation'],
    ['MEK_ERK_output_suppression', 'INT_PD_PI', '0.5', 'MAPK1_Y187_M1', 'down', 'PD-containing combination should reduce ERK2/MAPK1 activating phosphorylation'],
    ['MEK_ERK_output_suppression', 'PD', '0.5', 'MAPK3_Y204_M1', 'down', 'PD should reduce ERK1/MAPK3 activating phosphorylation'],
    ['MEK_ERK_output_suppression', 'INT_PD_PI', '0.5', 'MAPK3_Y204_M1', 'down', 'PD-containing combination should reduce ERK1/MAPK3 activating phosphorylation'],
    ['MEK_ERK_output_suppression', 'INT_FiveZ_PI', '0.5', 'MAPK3_Y204_M1', 'down', 'FiveZ_PI may suppress ERK/AP1-like output early'],
    ['MEK_ERK_output_suppression', 'FiveZ', '0.5', 'MAPK3_Y204_M1', 'down', 'FiveZ early PD-like TF state may involve reduced ERK/AP1-like signaling'],
    ['MEK_ERK_output_suppression', 'PD', '0.5', 'MAPK3_T202_M2', 'down', 'PD should reduce ERK1 activating phosphorylation'],
    ['MEK_ERK_output_suppression', 'INT_PD_PI', '0.5', 'MAPK3_T202_M2', 'down', 'PD-containing combination should reduce ERK1 activating phosphorylation'],
    ['MEK_ERK_output_suppression', 'INT_FiveZ_PI', '0.5', 'MAPK3_T202_M2', 'down', 'FiveZ_PI may suppress ERK/AP1-like output early'],

    ['AP1_MYC_phosphorylation_suppression', 'PD', '2', 'MYC_S62_M2', 'down', 'ERK stabilizing/activating MYC S62 phosphorylation should decrease with MEK inhibition'],
    ['AP1_MYC_phosphorylation_suppression', 'PD', '8', 'MYC_S62_M2', 'down', 'ERK stabilizing/activating MYC S62 phosphorylation should remain reduced'],
    ['AP1_MYC_phosphorylation_suppression', 'INT_PD_PI', '2', 'MYC_S62_M2', 'down', 'INT_PD_PI should suppress ERK/MYC phosphorylation'],
    ['AP1_MYC_phosphorylation_suppression', 'INT_PD_PI', '8', 'MYC_S62_M2', 'down', 'INT_PD_PI should suppress ERK/MYC phosphorylation'],
    ['AP1_MYC_phosphorylation_suppression', 'INT_FiveZ_PI', '2', 'MYC_S62_M2', 'down', 'INT_FiveZ_PI should suppress MYC phosphorylation before late rebound'],
    ['AP1_MYC_phosphorylation_suppression', 'FiveZ', '2', 'MYC_S62_M2', 'down', 'FiveZ PD-like AP1 output suggests reduced MYC S62 phosphorylation'],
    ['AP1_MYC_phosphorylation_suppression', 'FiveZ', '0.5', 'JUN_S63_M1', 'down', 'FiveZ should reduce JNK/ERK-associated JUN S63 phosphorylation early'],
    ['AP1_MYC_phosphorylation_suppression', 'INT_FiveZ_PI', '0.5', 'JUN_S63_M1', 'down', 'FiveZ_PI should reduce JUN S63 phosphorylation early'],

    ['MTOR_translation_axis', 'PI', '0.5', 'EIF4EBP1_T46_M2', 'down', 'PI should reduce mTOR-dependent EIF4EBP1 phosphorylation'],
    ['MTOR_translation_axis', 'INT_PD_PI', '0.5', 'EIF4EBP1_T46_M2', 'down', 'PI-containing combination should reduce EIF4EBP1 phosphorylation'],
    ['MTOR_translation_axis', 'INT_FiveZ_PI', '0.5', 'EIF4EBP1_T46_M2', 'down', 'PI-containing combination should reduce EIF4EBP1 phosphorylation'],
    ['MTOR_translation_axis', 'PI', '8', 'EIF4EBP1_T46_M2', 'down', 'PI should reduce mTOR-dependent EIF4EBP1 phosphorylation at 8 h'],
    ['MTOR_translation_axis', 'INT_PD_PI', '8', 'EIF4EBP1_T46_M2', 'down', 'INT_PD_PI should reduce EIF4EBP1 phosphorylation at 8 h'],
    ['MTOR_translation_axis', 'INT_FiveZ_PI', '8', 'EIF4EBP1_T46_M2', 'down', 'INT_FiveZ_PI should reduce EIF4EBP1 phosphorylation at 8 h'],

    ['YAP_CTNNB1_plasticity_axis', 'PI', '0.5', 'CTNNB1_S552_M1', 'down', 'AKT-dependent beta-catenin S552 phosphorylation might decrease after PI'],
    ['YAP_CTNNB1_plasticity_axis', 'INT_PD_PI', '0.5', 'CTNNB1_S552_M1', 'down', 'PI-containing combination might decrease AKT-dependent beta-catenin S552 phosphorylation'],
    ['YAP_CTNNB1_plasticity_axis', 'INT_FiveZ_PI', '8', 'YAP1_S127_M2', 'down', 'Reduced YAP S127 inhibitory phosphorylation could support late YAP/TEAD rebound'],
    ['YAP_CTNNB1_plasticity_axis', 'INT_PD_PI', '8', 'YAP1_S127_M2', 'down', 'Reduced YAP S127 inhibitory phosphorylation could indicate YAP release'],

    ['STAT_IRF_stress_axis', 'PD', '8', 'STAT3_Y705_M1', 'up', 'Late stress/inflammatory state may involve STAT3 activating phosphorylation'],
    ['STAT_IRF_stress_axis', 'FiveZ', '8', 'STAT3_Y705_M1', 'up', 'Late stress/inflammatory state may involve STAT3 activating phosphorylation'],
    ['STAT_IRF_stress_axis', 'INT_PD_PI', '8', 'STAT3_Y705_M1', 'up', 'Late stress/inflammatory state may involve STAT3 activating phosphorylation'],
    ['STAT_IRF_stress_axis', 'INT_FiveZ_PI', '8', 'STAT3_Y705_M1', 'up', 'Late stress/inflammatory state may involve STAT3 activating phosphorylation'],
    ['STAT_IRF_stress_axis', 'INT_PD_PI', '8', 'IRF3_S175_M1', 'up', 'IRF/stress axis may involve increased IRF3 phosphorylation'],
  ] unless const_defined?(:PROTEOMICS_PTM_EXPECTATIONS)

  helper :proteomics_ptm_file do |source|
    case source.to_s
    when 'collapsed'
      proteomics_data_dir['ptm/PTM_collapsed_log2_norm_imp_18517.txt']
    when 'filtered'
      proteomics_data_dir['ptm/dDIA_PHOS_log2_filter75perc_18517psites.txt']
    else
      raise ParameterException, "Unknown PTM source: #{source}"
    end
  end

  helper :proteomics_time_pattern do |time|
    case time.to_s
    when '0.5'
      '0[-_]5h'
    when '2', '8'
      "#{time}h"
    else
      "#{Regexp.escape(time.to_s)}h"
    end
  end

  helper :proteomics_ptm_replicate_fields do |fields, raw_code, time, dmso=false|
    code = dmso ? 'DMSO' : raw_code
    pattern = proteomics_time_pattern(time)
    fields.select do |field|
      field.to_s =~ /^X#{Regexp.escape(code)}[_-]#{pattern}_Rep\d+/ 
    end
  end

  helper :proteomics_parse_ptm_key do |key|
    if key.to_s =~ /^(.+)_([STY])(\d+)_([^_]+)$/
      [$1, $2, $3.to_i, $2 + $3, $4]
    else
      [key.to_s.split('_').first, nil, nil, nil, nil]
    end
  end

  helper :proteomics_ptm_source_fields do |source|
    source.to_s == 'filtered' ? ['Localization'] : []
  end

  input :ptm_source, :select, 'PTM source file', 'collapsed', :select_options => %w(collapsed filtered)
  task :proteomics_ptm_file_index => :tsv do |ptm_source|
    path = proteomics_ptm_file(ptm_source)
    table = path.tsv :type => :list
    tsv = TSV.setup({}, 'ID~Source,File,Rows,Fields,ExtraFields')
    tsv[1] = [ptm_source, path.to_s, table.keys.length, table.fields.length, (table.fields - table.fields.select{|f| f =~ /^X/ }).length]
    tsv
  end

  input :ptm_source, :select, 'PTM source file', 'collapsed', :select_options => %w(collapsed filtered)
  task :proteomics_ptm_long => :tsv do |ptm_source|
    path = proteomics_ptm_file(ptm_source)
    table = path.tsv :type => :list
    extra_localization_field = table.fields.find{|f| f.to_s.include?('PTM_localization') }
    fields = %w(Site Gene Residue Position Modification Treatment Time TreatmentMean DMSOMean Difference TreatmentValidValues DMSOValidValues Localization Source)
    tsv = TSV.setup({}, 'ID~' + fields * ',')
    id = 0

    table.through do |site_key, values|
      values = NamedArray.setup(values, table.fields)
      gene, residue, position, site, modification = proteomics_parse_ptm_key(site_key)
      localization = extra_localization_field ? proteomics_scalar(values[extra_localization_field]) : nil
      PROTEOMICS_TREATMENTS.each do |treatment|
        raw = proteomics_raw_treatment_code(treatment)
        PROTEOMICS_TIMES.each do |time|
          treatment_fields = proteomics_ptm_replicate_fields(table.fields, raw, time, false)
          dmso_fields = proteomics_ptm_replicate_fields(table.fields, raw, time, true)
          treatment_values = treatment_fields.collect{|f| proteomics_float(values[f]) }.compact
          dmso_values = dmso_fields.collect{|f| proteomics_float(values[f]) }.compact
          treatment_mean = treatment_values.empty? ? nil : treatment_values.inject(0.0){|acc,v| acc + v } / treatment_values.length
          dmso_mean = dmso_values.empty? ? nil : dmso_values.inject(0.0){|acc,v| acc + v } / dmso_values.length
          difference = treatment_mean.nil? || dmso_mean.nil? ? nil : treatment_mean - dmso_mean
          id += 1
          tsv[id] = [site_key, gene, residue, position, modification, treatment, time, treatment_mean, dmso_mean, difference, treatment_values.length, dmso_values.length, localization, ptm_source]
        end
      end
    end
    tsv
  end

  dep :proteomics_ptm_long
  task :proteomics_ptm_matrix => :tsv do
    long = step(:proteomics_ptm_long).load
    fields = PROTEOMICS_TREATMENTS.collect{|treatment| PROTEOMICS_TIMES.collect{|time| "#{treatment}-T#{time}" } }.flatten
    data = {}
    long.through do |id, values|
      site = proteomics_scalar(values['Site'])
      field = "#{proteomics_scalar(values['Treatment'])}-T#{proteomics_scalar(values['Time'])}"
      next unless fields.include?(field)
      data[site] ||= Array.new(fields.length)
      data[site][fields.index(field)] = proteomics_float(values['Difference'])
    end
    TSV.setup(data, :key_field => 'PTM collapse key', :fields => fields, :type => :list, :namespace => AGS.organism)
  end

  helper :proteomics_ptm_support_label do |expected_direction, difference|
    return 'not_measured' if difference.nil?
    observed = difference > 0 ? 'up' : (difference < 0 ? 'down' : 'zero')
    match = observed == expected_direction.to_s
    if match && difference.abs >= 0.5
      'strong_support'
    elsif match && difference.abs >= 0.25
      'moderate_support'
    elsif match
      'directional_support_small_effect'
    elsif difference.abs >= 0.5
      'strong_opposite'
    elsif difference.abs >= 0.25
      'moderate_opposite'
    else
      'no_support_or_small_opposite'
    end
  end

  dep :proteomics_ptm_long
  input :ptm_source, :select, 'PTM source file', 'collapsed', :select_options => %w(collapsed filtered)
  task :proteomics_ptm_candidate_site_report => :tsv do |ptm_source|
    long = step(:proteomics_ptm_long).load
    lookup = {}
    long.through do |id, values|
      lookup[[proteomics_scalar(values['Site']), proteomics_scalar(values['Treatment']), proteomics_scalar(values['Time'])]] = values
    end

    fields = %w(Mechanism Treatment Time Site Gene ExpectedDirection ObservedDifference ObservedDirection Support TreatmentValidValues DMSOValidValues Localization Rationale)
    tsv = TSV.setup({}, 'ID~' + fields * ',')
    id = 0
    PROTEOMICS_PTM_EXPECTATIONS.each do |mechanism, treatment, time, site, expected, rationale|
      values = lookup[[site, treatment, time]]
      difference = values ? proteomics_float(values['Difference']) : nil
      observed = difference.nil? ? nil : (difference > 0 ? 'up' : (difference < 0 ? 'down' : 'zero'))
      support = proteomics_ptm_support_label(expected, difference)
      gene = values ? proteomics_scalar(values['Gene']) : proteomics_parse_ptm_key(site).first
      id += 1
      tsv[id] = [mechanism, treatment, time, site, gene, expected, difference, observed, support, values ? proteomics_scalar(values['TreatmentValidValues']) : nil, values ? proteomics_scalar(values['DMSOValidValues']) : nil, values ? proteomics_scalar(values['Localization']) : nil, rationale]
    end
    tsv
  end

  dep :proteomics_ptm_candidate_site_report
  dep :proteomics_ptm_file_index
  task :proteomics_ptm_mechanism_summary => :tsv do
    report = step(:proteomics_ptm_candidate_site_report).load
    index = step(:proteomics_ptm_file_index).load
    counts = Hash.new(0)
    report.through do |id, values|
      counts[[proteomics_scalar(values['Mechanism']), proteomics_scalar(values['Support'])]] += 1
    end
    tsv = TSV.setup({}, 'ID~Section,Item,Value,Comment')
    id = 0
    counts.keys.collect{|k| k[0] }.uniq.sort.each do |mechanism|
      %w(strong_support moderate_support directional_support_small_effect no_support_or_small_opposite moderate_opposite strong_opposite not_measured).each do |support|
        value = counts[[mechanism, support]]
        next if value == 0
        id += 1
        tsv[id] = ['ptm_mechanism_support', mechanism + ':' + support, value, 'Counts over targeted PTM expectations']
      end
    end
    index.through do |idx, values|
      id += 1
      tsv[id] = ['ptm_file_index', proteomics_scalar(values['Source']), proteomics_scalar(values['Rows']), "Fields: #{proteomics_scalar(values['Fields'])}; extra fields: #{proteomics_scalar(values['ExtraFields'])}"]
    end
    tsv
  end
end
