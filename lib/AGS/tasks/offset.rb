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

  dep :timepoint_decoupler, time_point: :placeholder, data_type: :fc, dynamic: false do |jobname,options|
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
          next if value == 0
          if (value > 0) != (current[field] > 0)
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
      exclude_tsv = score_job.load
      exclude_tsv.through do |gene,values|
        current = orig[gene]
        next if current.nil?
        values.to_hash.each do |field,value|
          score = current["Score #{field}"]
          next if value == 0
          next if current[field].nil?
          next if current[field] == 0
          if (value > 0) != (current[field] > 0)
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


  dep :treatment_tfs_priority
  dep :change_offsets
  task :treatment_tfs_priority_consistency => :tsv do
    treatment =recursive_inputs[:treatment]
    changes = step(:change_offsets).load
    fields = AGS::TIME_POINTS.collect{|t| "Consistent at #{t}h" }
    consistency = TSV.setup({}, key_field: "Associated Gene Name", fields: fields, type: :list)
    traverse step(:treatment_tfs_priority), into: consistency do |gene,values|
      next unless changes.include?(gene)

      gene_changes = changes[gene][treatment]
      negative_activities = gene_changes.select{|c| c.include?("decrease")}.
        collect{|c| c.split(" ").last.to_i}.collect{|t| AGS::TIME_POINTS.index(t)}.
        collect{|i| [i, i+1]}.flatten.reject{|i| i > 5}.uniq.
        collect{|i| AGS::TIME_POINTS[i] }

      positive_activities = gene_changes.select{|c| c.include?("increase")}.
        collect{|c| c.split(" ").last.to_i}.collect{|t| AGS::TIME_POINTS.index(t)}.
        collect{|i| [i, i+1]}.flatten.reject{|i| i > 5}.uniq.
        collect{|i| AGS::TIME_POINTS[i] }

      
      res = []
      AGS::TIME_POINTS.each_with_index do |time_point,i|
        res << (values[i] > 0 && positive_activities.include?(time_point)) || (values[i] < 0 && negative_activities.include?(time_point))
      end
      res = res.collect{|c| c ? 1 : 0 }
      [gene, res]
    end
    tsv = step(:treatment_tfs_priority).load.attach consistency
    tsv.cast = nil
    tsv
  end

end
