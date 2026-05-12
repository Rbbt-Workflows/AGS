module AGS

  dep :timepoint_decoupler, time_point: :placeholder, data_type: :fc0, dynamic: false do |jobname,options|
    AGS::TIME_POINTS.collect do |time_point,index|
      options.merge({time_point: time_point})
    end
  end
  task :treatment_tfs_diff => :tsv do
    tsvs = []
    AGS::TIME_POINTS.each_with_index do |time_point,index|
      job = dependencies.select{|d| d.recursive_inputs[:time_point] == time_point }.first
      if index > 0
        prev_time_point = AGS::TIME_POINTS[index-1]
        previous_job = dependencies.select{|d| d.recursive_inputs[:time_point] == prev_time_point }.first
      end

      tsv = job.load

      if previous_job
        previous_tsv = previous_job.load
        tsv.fields do |field|
          tsv.process field do |k,v|
            v.to_f - previous_tsv[k][field].to_f
          end
        end
      end

      tsvs << tsv
    end

    tsvs.inject do |acc,tsv| acc.attach tsv, complete: true end
  end

  dep :timepoint_decoupler, time_point: :placeholder, dynamic: false, data_type: :fc do |jobname,options|
    AGS::TIME_POINTS.collect do |time_point,index|
      options.merge({time_point: time_point})
    end
  end
  task :treatment_tfs_non_dynamic => :tsv do
    tsvs = []
    AGS::TIME_POINTS.each_with_index do |time_point,index|
      job = dependencies.select{|d| d.recursive_inputs[:time_point] == time_point }.first

      tsv = job.load

      tsvs << tsv
    end

    tsvs.inject do |acc,tsv| acc.attach tsv, complete: true end
  end

  dep :timepoint_decoupler, time_point: :placeholder, dynamic: false, data_type: :fc0 do |jobname,options|
    AGS::TIME_POINTS.collect do |time_point,index|
      options.merge({time_point: time_point})
    end
  end
  task :treatment_tfs_fc0 => :tsv do
    tsvs = []
    AGS::TIME_POINTS.each_with_index do |time_point,index|
      job = dependencies.select{|d| d.recursive_inputs[:time_point] == time_point }.first

      tsv = job.load

      tsvs << tsv
    end

    tsvs.inject do |acc,tsv| acc.attach tsv, complete: true end
  end

  input :data_type, :select, "Values to use", :range, :select_options => %w(fc fc0 ext_fc binary range)
  dep :timepoint_decoupler, time_point: :placeholder, dynamic: true do |jobname,options|
    AGS::TIME_POINTS.collect do |time_point,index|
      options.merge({time_point: time_point})
    end
  end
  task :treatment_tfs_dynamic => :tsv do
    tsvs = []
    AGS::TIME_POINTS.each_with_index do |time_point,index|
      job = dependencies.select{|d| d.recursive_inputs[:time_point] == time_point }.first

      tsv = job.load

      tsvs << tsv
    end

    tsvs.inject do |acc,tsv| acc.attach tsv, complete: true end
  end


  dep :treatment_tfs_dynamic, threshold: 0.05, synthesis: false, vetting: 'degradation', data_type: :range

  #Exclude 
  dep :treatment_tfs_dynamic, threshold: 0.05, vetting: 'none', data_type: :range
  dep :treatment_tfs_dynamic, threshold: 0.05, vetting: 'synthesis', data_type: :range
  dep :treatment_tfs_dynamic, threshold: 0.05, vetting: 'degradation', data_type: :range
  dep :treatment_tfs_dynamic, threshold: 0.05, vetting: 'relaxed_degradation', data_type: :range

  #Score 
  dep :treatment_tfs_dynamic, threshold: 0.01, vetting: 'degradation', data_type: :range
  dep :treatment_tfs_dynamic, threshold: 0.05, vetting: 'extended_degradation', data_type: :range
  dep :treatment_tfs_non_dynamic, threshold: 0.05, vetting: 'degradation'
  dep :treatment_tfs_diff, threshold: 0.05, vetting: 'degradation'
  dep :treatment_tfs_dynamic, threshold: 0.05, vetting: 'none', data_type: :range
  input :score_points, :array, "Score points for each scoring job", [1, 2, 0.5, 0.5, 0]
  task :treatment_tfs_priority => :tsv do |score_points|
    orig, exclude1, exclude2, exclude3, exclude4, *score = dependencies

    orig = orig.load

    [exclude1, exclude2, exclude3, exclude4].each do |exclude|
      exclude_tsv = exclude.load
      exclude_tsv.through do |gene,values|
        current = orig[gene]
        next if current.nil?
        values.to_hash.each do |field,value|
          value = value.to_f
          next if value == 0
          if (value > 0) != (current[field].to_f > 0)
            current[field] = 0
          end
        end
        orig[gene] = current
      end
    end

    orig.fields.each do |field|
      orig.add_field "Score #{field}" do
        0
      end
    end

    score.zip(score_points).each do |score_job,score_points|
      score_points = score_points.to_f
      exclude_tsv = score_job.load
      exclude_tsv.through do |gene,values|
        current = orig[gene]
        next if current.nil?
        values.to_hash.each do |field,value|
          value = value.to_f
          score = current["Score #{field}"]
          score = score.to_f
          next if value == 0
          next if current[field].nil?
          next if current[field] == 0
          if (value > 0) != (current[field].to_f > 0)
            score -= score_points
          else
            score += score_points
          end
          current["Score #{field}"] = score
        end
        orig[gene] = current
      end
    end

    orig
  end

  input :scheme, :select, "Scheme to use: diff, dynamic, non-dynamic", :dynamic, select_options: %w(priority dynamic non-dynamic diff fc0)
  task_alias :treatment_tfs, AGS, :treatment_tfs_dynamic do |jobname,options|
    case options[:scheme].to_s
    when "priority"
      {task: :treatment_tfs_priority, jobname: jobname, options: options}
    when "dynamic"
      {task: :treatment_tfs_dynamic, jobname: jobname, options: options}
    when "non-dynamic"
      {task: :treatment_tfs_non_dynamic, jobname: jobname, options: options}
    when "diff"
      {task: :treatment_tfs_diff, jobname: jobname, options: options}
    when "fc0"
      {task: :treatment_tfs_fc0, jobname: jobname, options: options}
    end
  end


  dep :treatment_tfs, vetting: :none, data_type: :range, treatment: :placeholder do |jobname,options|
    jobs = []
    AGS::TREATMENTS.each do |treatment|
      %w(ulm).each do |method|
        jobs << options.merge(treatment: treatment)
      end
    end
    jobs
  end
  task :tf_predictions => :tsv do
    dependencies.inject(nil) do |acc,job|
      tsv = job.load
      if acc.nil?
        acc = tsv
      else
        acc.attach tsv, complete: true
      end
    end
  end
end
