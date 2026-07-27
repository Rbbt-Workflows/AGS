
module AGS
  dep :regulome
  dep :full_gene_info
  dep :full_gene_info, low_min_24h:  1, mid_min_24h: 1, high_min_24h: 1, next_min_24h: 1 
  dep :full_gene_info, low_min_24h:  0.5, mid_min_24h: 0.5, high_min_24h: 0.5, next_min_24h: 0.5
  dep :treatment_tf_consistency do
    TREATMENTS.collect do |treatment|
      {treatment: treatment}
    end
  end
  dep :sequence_with_changes do
    TREATMENTS.collect do |treatment|
      {treatment: treatment}
    end
  end
  dep :consistency_counts
  dep :consistency_counts, low_min_24h:  1, mid_min_24h: 1, high_min_24h: 1, next_min_24h: 1 
  dep :consistency_counts, low_min_24h:  0.5, mid_min_24h: 0.5, high_min_24h: 0.5, next_min_24h: 0.5
  dep :neko_bootstrap_sweep
  dep :neko_bootstrap_consistency do
    %w(PD PI FiveZ).collect do |treatment|
      {treatment: treatment}
    end
  end
  dep :tf_predictions, scheme: :placeholder do |jobname,options|
    IndiferentHash.setup(options)
    jobs = %w(dynamic non-dynamic fc0 diff).collect do |scheme|
      if scheme == 'dynamic'
        [
          options.merge({scheme: scheme}),
          options.merge({scheme: scheme, low_min_24h:  1, mid_min_24h: 1, high_min_24h: 1, next_min_24h: 1 }), 
          options.merge({scheme: scheme, low_min_24h:  0.5, mid_min_24h: 0.5, high_min_24h: 0.5, next_min_24h: 0.5 }) 
        ]
      else
        options.merge({scheme: scheme})
      end
    end.flatten
    jobs
  end
  dep :valid_TFs
  dep :interval_fold_changes_fc1
  dep :interval_pvalue_surrogate_fc1
  dep :interval_de_gene_counts_fc1
  dep :fc1_onset_relationship_summary
  dep :fc1_onset_relationship_8h_focus
  dep :de_gene_counts_fc0
  dep :onset_first_counts
  dep :onset_episode_counts
  dep :onset_direction_switch_summary
  dep :tf_activity_call_counts_dynamic
  dep :tf_activity_call_counts_by_scheme, scheme: :placeholder do |jobname,options|
    %w(dynamic non-dynamic).collect{|scheme| options.merge(scheme: scheme) }
  end
  dep :tf_activity_heatmap_matrix, scheme: :placeholder, normalization: :placeholder do |jobname,options|
    %w(dynamic non-dynamic).collect do |scheme|
      %w(raw row_zscore column_zscore).collect do |normalization|
        options.merge(scheme: scheme, normalization: normalization)
      end
    end.flatten
  end
  dep :neko_dynamic_non_dynamic_summary
  dep :neko_dynamic_vs_non_dynamic_odds
  dep :self_consistency_dynamic_non_dynamic_summary
  dep :self_consistency_dynamic_vs_non_dynamic_odds
  dep :tf_timepoint_report_card
  dep :tf_target_report_card
  dep :tf_target_edge_consistency_summary
  dep :combination_tf_categories
  dep :combination_tf_category_counts
  dep :functional_enrichment_suite
  dep :proteomics_prediction_suite
  dep :grns
  task :freeze => :array do
    dependencies.each do |dep|
      other = dependencies.select{|d| d.task_name == dep.task_name }.length > 1
      filename = case dep.task_name
                 when :functional_enrichment_suite
                   Open.cp dep.files_dir, file("functional_enrichment_suite")
                 when :grns
                   Open.cp dep.files_dir, file("grns")
                 when :proteomics_prediction_suite
                   Open.cp dep.files_dir, file("proteomics_prediction_suite")
                 when :valid_TFs
                   'valid_TFs.list'
                 when :list_tfs
                   treatment, time_point, direction, threshold, scheme = dep.recursive_inputs.values_at :treatment, :time_point, :direction, :high_min_24h, :scheme
                   scheme = nil if scheme.to_s == 'dynamic'
                   threshold = nil if threshold == 0.25
                   threshold = "T24_fc_cutoff_#{threshold}" if threshold
                   time_point = "#{time_point}h" if time_point
                   ["list_tfs", treatment, time_point, direction, threshold, scheme].compact * "-" + ".list"
                 when :tf_predictions
                   treatment, time_point, direction, threshold, scheme = dep.recursive_inputs.values_at :treatment, :time_point, :direction, :high_min_24h, :scheme
                   scheme = nil if scheme.to_s == 'dynamic'
                   treatment = nil
                   threshold = nil if threshold && threshold == 0.25
                   threshold = "T24_fc_cutoff_#{threshold}" if threshold
                   time_point = "#{time_point}h" if time_point
                   time_point = nil
                   [dep.task_name.to_s, treatment, time_point, direction, threshold, scheme].compact * "-" + ".tsv"
                 when :consistency_counts
                   treatment, time_point, direction, threshold, scheme = dep.recursive_inputs.values_at :treatment, :time_point, :direction, :high_min_24h, :scheme
                   scheme = nil if scheme.to_s == 'dynamic'
                   treatment = nil
                   threshold = nil if threshold && threshold == 0.25
                   threshold = "T24_fc_cutoff_#{threshold}" if threshold
                   time_point = "#{time_point}h" if time_point
                   time_point = nil
                   [dep.task_name.to_s, treatment, time_point, direction, threshold, scheme].compact * "-" + ".tsv"
                 else
                   treatment, time_point, direction, threshold, scheme, normalization = dep.recursive_inputs.values_at :treatment, :time_point, :direction, :high_min_24h, :scheme, :normalization
                   treatment = nil unless other
                   scheme = nil if scheme.to_s == 'dynamic'
                   threshold = nil if threshold && threshold.to_f == 0.25
                   threshold = "T24_fc_cutoff_#{threshold}" if threshold
                   time_point = "#{time_point}h" if time_point
                   time_point = nil
                   normalization = nil if normalization.nil? || normalization.to_s.empty?
                   [dep.task_name.to_s, treatment, time_point, direction, threshold, scheme, normalization].compact * "-" + ".tsv"
                 end
      
      next if filename.nil?

      Open.cp dep.path, file(filename)
      file(filename)
    end

    # Full gene info

    info = file('full_gene_info.tsv').tsv

    cluster_fields = info.fields.select{|f| f.include?('FC clusters') }

    info1 = file('full_gene_info-T24_fc_cutoff_1.tsv').tsv fields: cluster_fields
    info05 = file('full_gene_info-T24_fc_cutoff_0.5.tsv').tsv fields: cluster_fields

    info1.fields = info1.fields.collect{|f| f + ' T24_fc_cutoff_1' }
    info05.fields = info05.fields.collect{|f| f + ' T24_fc_cutoff_0.5' }

    info.attach info1
    info.attach info05

    file('full_gene_info_extended.tsv').write info.to_s

    # TF Predictions

    preds = file('tf_predictions.tsv').tsv
    %w( T24_fc_cutoff_0.5 T24_fc_cutoff_1 diff fc0 non-dynamic).each do |tag|
      new = file("tf_predictions-#{tag}.tsv").tsv
      new.fields = new.fields.collect{|f| f + " #{tag}" }
      preds.attach new
    end
    file('tf_predictions_extended.tsv').write preds.to_s

    files
  end
end
