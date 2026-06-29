module AGS

  helper :proteomics_ptm_file do |source|
    case source.to_s
    when 'collapsed', 'imputed'
      proteomics_data_dir['ptm/PTM_collapsed_log2_norm_imp_18517.txt']
    when 'filtered', 'raw'
      proteomics_data_dir['ptm/dDIA_PHOS_log2_filter75perc_18517psites.txt']
    when 'sites'
      proteomics_data_dir['ptm/sites']
    else
      raise ParameterException, "Unknown PTM source: #{source}"
    end
  end

  helper :proteomics_ptm_control_fields do |fields|
    fields.select do |field|
      field.to_s =~ /^XCONTROL_0h_Rep\d+/ || field.to_s =~ /^XCTRL_0h_Rep\d+/
    end
  end

  helper :proteomics_ptm_time_pattern do |time, source='collapsed'|
    normalized = proteomics_normalize_time(time)
    if normalized == '0.5'
      source.to_s == 'filtered' || source.to_s == 'raw' ? '0-5h' : '0_5h'
    else
      "#{normalized}h"
    end
  end

  helper :proteomics_ptm_replicate_fields do |fields, treatment, time, source='collapsed'|
    raw_code = proteomics_raw_treatment_code(treatment)
    pattern = proteomics_ptm_time_pattern(time, source)
    fields.select do |field|
      field.to_s =~ /^X#{Regexp.escape(raw_code)}_#{Regexp.escape(pattern)}_Rep\d+/
    end
  end

  helper :proteomics_parse_ptm_key do |key|
    if key.to_s =~ /^(.+)_([STY])(\d+)_([^_]+)$/
      [$1, $2 + $3, $2, $3.to_i, $4]
    else
      [key.to_s.split('_').first, nil, nil, nil, nil]
    end
  end

  input :ptm_source, :select, 'PTM source file', 'collapsed', :select_options => %w(collapsed filtered)
  task :proteomics_ptm_file_index => :tsv do |ptm_source|
    path = proteomics_ptm_file(ptm_source)
    table = path.tsv :type => :list
    measurement_fields = table.fields.select{|f| f.to_s =~ /^X/ }
    extra_fields = table.fields - measurement_fields
    tsv = TSV.setup({}, 'ID~Source,File,Rows,Fields,MeasurementFields,ExtraFields')
    tsv[1] = [ptm_source, path.to_s, table.keys.length, table.fields.length, measurement_fields.length, extra_fields * '|']
    tsv
  end

  input :ptm_source, :select, 'PTM source file', 'collapsed', :select_options => %w(collapsed filtered)
  task :proteomics_ptm_long => :tsv do |ptm_source|
    path = proteomics_ptm_file(ptm_source)
    table = path.tsv :type => :list
    localization_field = table.fields.find{|f| f.to_s.include?('PTM_localization') }
    control_fields = proteomics_ptm_control_fields(table.fields)

    fields = %w(Site Protein Residue ResidueAA Position Modification Treatment Time TreatmentMean BaselineMean DMSOMean DifferenceVsBaseline DifferenceVsDMSO TreatmentValidValues BaselineValidValues DMSOValidValues Localization Source)
    tsv = TSV.setup({}, 'ID~' + fields * ',')
    id = 0

    table.through bar: self.progress_bar do |site_key, values|
      values = NamedArray.setup(values, table.fields)
      protein, residue, residue_aa, position, modification = proteomics_parse_ptm_key(site_key)
      baseline_values = control_fields.collect{|f| proteomics_float(values[f]) }.compact
      baseline_mean = proteomics_mean(baseline_values)
      localization = localization_field ? proteomics_scalar(values[localization_field]) : nil

      PROTEOMICS_TREATMENTS.each do |treatment|
        PROTEOMICS_TIMES.each do |time|
          treatment_fields = proteomics_ptm_replicate_fields(table.fields, treatment, time, ptm_source)
          dmso_fields = proteomics_ptm_replicate_fields(table.fields, 'DMSO', time, ptm_source)
          treatment_values = treatment_fields.collect{|f| proteomics_float(values[f]) }.compact
          dmso_values = dmso_fields.collect{|f| proteomics_float(values[f]) }.compact
          treatment_mean = proteomics_mean(treatment_values)
          dmso_mean = proteomics_mean(dmso_values)
          difference_baseline = treatment_mean.nil? || baseline_mean.nil? ? nil : treatment_mean - baseline_mean
          difference_dmso = if treatment == 'DMSO'
                              0.0
                            else
                              treatment_mean.nil? || dmso_mean.nil? ? nil : treatment_mean - dmso_mean
                            end
          id += 1
          tsv[id] = [site_key, protein, residue, residue_aa, position, modification, treatment, time, treatment_mean, baseline_mean, dmso_mean, difference_baseline, difference_dmso, treatment_values.length, baseline_values.length, dmso_values.length, localization, ptm_source]
        end
      end
    end

    tsv
  end

#  dep :proteomics_ptm_long
#  task :proteomics_ptm_matrix => :tsv do
#    long = step(:proteomics_ptm_long).load
#    fields = PROTEOMICS_TREATMENTS.collect{|treatment| PROTEOMICS_TIMES.collect{|time| "#{treatment}-T#{time}" } }.flatten
#    data = {}
#    long.through do |id, values|
#      site = proteomics_scalar(values['Site'])
#      field = "#{proteomics_scalar(values['Treatment'])}-T#{proteomics_scalar(values['Time'])}"
#      next unless fields.include?(field)
#      data[site] ||= Array.new(fields.length)
#      data[site][fields.index(field)] = proteomics_float(values['DifferenceVsBaseline'])
#    end
#    TSV.setup(data, :key_field => 'PTM collapse key', :fields => fields, :type => :list, :namespace => AGS.organism)
#  end
#
#  dep :proteomics_ptm_long
#  task :proteomics_ptm_matrix_vs_dmso => :tsv do
#    long = step(:proteomics_ptm_long).load
#    fields = PROTEOMICS_DRUG_TREATMENTS.collect{|treatment| PROTEOMICS_TIMES.collect{|time| "#{treatment}-T#{time}" } }.flatten
#    data = {}
#    long.through do |id, values|
#      treatment = proteomics_scalar(values['Treatment'])
#      next if treatment == 'DMSO'
#      site = proteomics_scalar(values['Site'])
#      field = "#{treatment}-T#{proteomics_scalar(values['Time'])}"
#      next unless fields.include?(field)
#      data[site] ||= Array.new(fields.length)
#      data[site][fields.index(field)] = proteomics_float(values['DifferenceVsDMSO'])
#    end
#    TSV.setup(data, :key_field => 'PTM collapse key', :fields => fields, :type => :list, :namespace => AGS.organism)
#  end
end
