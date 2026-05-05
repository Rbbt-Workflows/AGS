
module AGS
  dep :regulome
  dep :full_gene_info
  dep :treatment_tfs do 
    TREATMENTS.collect{|t| {treatment: t} }
  end
  dep :treatment_tfs, vetting: :synthesis do 
    TREATMENTS.collect{|t| {treatment: t} }
  end
  dep :treatment_tfs, vetting: :degradation do 
    TREATMENTS.collect{|t| {treatment: t} }
  end
  dep :treatment_tfs, vetting: :relaxed_degradation do 
    TREATMENTS.collect{|t| {treatment: t} }
  end
  dep :treatment_tfs, low_min_24h:  1, mid_min_2eh: 1, high_min_24h: 1, next_min_24h: 1 do 
    TREATMENTS.collect{|t| {treatment: t} }
  end
  dep :treatment_tfs, low_min_24h:  0.5, mid_min_2eh: 0.5, high_min_24h: 0.5, next_min_24h: 0.5 do 
    TREATMENTS.collect{|t| {treatment: t} }
  end
  dep :treatment_tfs_diff do 
    TREATMENTS.collect{|t| {treatment: t} }
  end
  dep :treatment_tfs_non_dynamic do 
    TREATMENTS.collect{|t| {treatment: t} }
  end
  dep :treatment_tfs_fc0 do 
    TREATMENTS.collect{|t| {treatment: t} }
  end
  dep :list_tfs do
    jobs = []
    TREATMENTS.each do |treatment|
      TIME_POINTS.each do |time_point|
        %w(up down).each do |direction|
          jobs << {treatment: treatment, time_point: time_point, direction: direction}
        end
      end
    end
    jobs
  end
  dep :list_tfs do
    jobs = []
    TREATMENTS.each do |treatment|
      time_point = 24
      %w(up down).each do |direction|
        [0.5, 1].each do |threshold|
          jobs << {
            treatment: treatment, time_point: time_point, direction: direction,
            low_min_24h:  threshold, mid_min_2eh: threshold, high_min_24h: threshold, next_min_24h: threshold
          }
        end
      end
    end
    jobs
  end
  dep :treatment_tf_consistency do
    TREATMENTS.collect do |treatment|
      {treatment: treatment}
    end
  end
  dep :neko_bootstrap_consistency do
    %w(PD PI FiveZ).collect do |treatment|
      {treatment: treatment}
    end
  end
  task :freeze => :array do
    dependencies.collect do |dep|
      filename = case dep.task_name
                 when :list_tfs
                   treatment, time_point, direction, threshold = dep.recursive_inputs.values_at :treatment, :time_point, :direction, :high_min_24h 
                   threshold = nil if threshold <= 0.25
                   threshold = "fc_#{threshold}" if threshold
                   time_point = "#{time_point}h" if time_point
                   ["list_tfs", treatment, time_point, direction, threshold].compact * "-" + ".list"
                 else
                   treatment, time_point, direction, threshold = dep.recursive_inputs.values_at :treatment, :time_point, :direction, :high_min_24h 
                   threshold = nil if threshold && threshold <= 0.25
                   threshold = "fc_#{threshold}" if threshold
                   time_point = "#{time_point}h" if time_point
                   time_point = nil
                   [dep.task_name.to_s, treatment, time_point, direction, threshold].compact * "-" + ".tsv"
                 end
      
      next if filename.nil?

      Open.cp dep.path, file(filename)
      file(filename)
    end.compact
  end
end
