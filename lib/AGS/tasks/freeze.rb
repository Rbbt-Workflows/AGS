
module AGS
  #dep :treatment_tfs do 
  #  TREATMENTS.collect{|t| {treatment: t} }
  #end
  #dep :treatment_tfs, vetting: :synthesis do 
  #  TREATMENTS.collect{|t| {treatment: t} }
  #end
  #dep :treatment_tfs, vetting: :degradation do 
  #  TREATMENTS.collect{|t| {treatment: t} }
  #end
  #dep :treatment_tfs, vetting: :relaxed_degradation do 
  #  TREATMENTS.collect{|t| {treatment: t} }
  #end
  #dep :treatment_tfs, low_min_24h:  1, mid_min_24h: 1, high_min_24h: 1, next_min_24h: 1 do 
  #  TREATMENTS.collect{|t| {treatment: t} }
  #end
  #dep :treatment_tfs, low_min_24h:  0.5, mid_min_24h: 0.5, high_min_24h: 0.5, next_min_24h: 0.5 do 
  #  TREATMENTS.collect{|t| {treatment: t} }
  #end
  #dep :treatment_tfs_diff do 
  #  TREATMENTS.collect{|t| {treatment: t} }
  #end
  #dep :treatment_tfs_non_dynamic do 
  #  TREATMENTS.collect{|t| {treatment: t} }
  #end
  #dep :treatment_tfs_fc0 do 
  #  TREATMENTS.collect{|t| {treatment: t} }
  #end
  #dep :list_tfs do
  #  jobs = []
  #  TREATMENTS.each do |treatment|
  #    TIME_POINTS.each do |time_point|
  #      %w(up down).each do |direction|
  #        jobs << {treatment: treatment, time_point: time_point, direction: direction}
  #      end
  #    end
  #  end
  #  jobs
  #end
  #dep :list_tfs do
  #  jobs = []
  #  TREATMENTS.each do |treatment|
  #    time_point = 24
  #    %w(up down).each do |direction|
  #      [0.5, 1].each do |threshold|
  #        jobs << {
  #          treatment: treatment, time_point: time_point, direction: direction,
  #          low_min_24h:  threshold, mid_min_24h: threshold, high_min_24h: threshold, next_min_24h: threshold
  #        }
  #      end
  #    end
  #  end
  #  jobs
  #end
  dep :regulome
  dep :full_gene_info
  dep :full_gene_info, low_min_24h:  1, mid_min_24h: 1, high_min_24h: 1, next_min_24h: 1 
  dep :full_gene_info, low_min_24h:  0.5, mid_min_24h: 0.5, high_min_24h: 0.5, next_min_24h: 0.5
  dep :treatment_tf_consistency do
    TREATMENTS.collect do |treatment|
      {treatment: treatment}
    end
  end
  dep :consistency_counts
  dep :neko_bootstrap_sweep
  dep :neko_bootstrap_consistency do
    %w(PD PI FiveZ).collect do |treatment|
      {treatment: treatment}
    end
  end
  dep :tf_predictions, scheme: :placeholder do |jobname,options|
    %w(dynamic non-dynamic fc0 diff).collect do |scheme|
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
  end
  dep :valid_TFs
  task :freeze => :array do
    dependencies.each do |dep|
      other = dependencies.select{|d| d.task_name == dep.task_name }.length > 1
      filename = case dep.task_name
                 when :valid_TFs
                   'valid_TFs.list'
                 when :list_tfs
                   treatment, time_point, direction, threshold, scheme = dep.recursive_inputs.values_at :treatment, :time_point, :direction, :high_min_24h, :scheme
                   scheme = nil if scheme.to_s == 'dynamic'
                   threshold = nil if threshold <= 0.25
                   threshold = "T24_fc_cutoff_#{threshold}" if threshold
                   time_point = "#{time_point}h" if time_point
                   ["list_tfs", treatment, time_point, direction, threshold, scheme].compact * "-" + ".list"
                 when :tf_predictions
                   treatment, time_point, direction, threshold, scheme = dep.recursive_inputs.values_at :treatment, :time_point, :direction, :high_min_24h, :scheme
                   scheme = nil if scheme.to_s == 'dynamic'
                   treatment = nil
                   threshold = nil if threshold && threshold <= 0.25
                   threshold = "T24_fc_cutoff_#{threshold}" if threshold
                   time_point = "#{time_point}h" if time_point
                   time_point = nil
                   [dep.task_name.to_s, treatment, time_point, direction, threshold, scheme].compact * "-" + ".tsv"
                 else
                   treatment, time_point, direction, threshold, scheme = dep.recursive_inputs.values_at :treatment, :time_point, :direction, :high_min_24h, :scheme
                   treatment = nil unless other
                   scheme = nil if scheme.to_s == 'dynamic'
                   threshold = nil if threshold && threshold <= 0.25
                   threshold = "T24_fc_cutoff_#{threshold}" if threshold
                   time_point = "#{time_point}h" if time_point
                   time_point = nil
                   [dep.task_name.to_s, treatment, time_point, direction, threshold, scheme].compact * "-" + ".tsv"
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
