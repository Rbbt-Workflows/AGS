module AGS

  #input :timepoint_suite, :tsv
  #task :consistent => :tsv do |suite|
  #  [true, false].each do |synthesis|
  #    [true, false].each do |dynamic|
  #      [:fc, :binary, :ext_fc].each do |data_type|
  #        [:mayority, :one, :two, :three].each do |synthesis_criteria|
  #          TREATMENTS.each do |treatment|
  #            dir_name = [data_type, synthesis_criteria, "synth=#{synthesis}", "dynamic=#{dynamic}"] * "-"
  #            dir = file(dir_name)
  #            Open.mkdir dir

  #            target = dir[time_point.to_s + "h.png"]
  #            fields = TIME_POINTS.collect do |time_point|
  #              field_name = [treatment, time_point, data_type, synthesis_criteria, "synth=#{synthesis}", "dynamic=#{dynamic}"] * "-"
  #              field_name
  #            end
  #            data = suite.slice(fields)
  #            Log.tsv data
  #          end
  #        end
  #      end
  #    end
  #  end
  #  Dir.glob(files_dir + "/**")
  #end

  helper :prev_time_point do |time_point|
    index = TIME_POINTS.index time_point
    if index > 0
      TIME_POINTS[index-1]
    else
      "NO TIMEPOINT"
    end
  end

  dep :timepoint_decoupler_suite
  dep :change_offsets, :jobname => "Default"
  task :consistent => :tsv do 
    suite = step(:timepoint_decoupler_suite).load
    clusters = step(:change_offsets).load
    consistent = TSV.setup({}, :key_field => "Associated Gene Name", :fields => [], :type => :list)

    fields = suite.fields
    TSV.traverse fields, :bar => true do |field|
      treatment, time_point = field.split("-") 
      time_point = time_point.to_i
      prev_time_point = prev_time_point(time_point)
      treatment_clusters = clusters.column(treatment).to_flat

      treatment_consistency = TSV.setup({}, :key_field => consistent.key_field, :fields => [field], :type => :single)

      suite.through nil, [field] do |gene,vs|
        v = vs.first.to_f
        gene_clusters = treatment_clusters[gene] || []
        consistency = if v.to_f == 0.0
                        next
                      elsif v.to_f > 0
                        gene_clusters.include?("increase #{time_point}h") || 
                          gene_clusters.include?("increase #{prev_time_point}h") ||
                          gene_clusters.start_with?("increase #{time_point}-") ||
                          gene_clusters.start_with?("increase #{prev_time_point}-")
                      else
                        gene_clusters.include?("decrease #{time_point}h") || 
                          gene_clusters.include?("decrease #{prev_time_point}h") ||
                          gene_clusters.start_with?("decrease #{time_point}-") || 
                          gene_clusters.start_with?("decrease #{prev_time_point}-")
                      end
        treatment_consistency[gene] = consistency
      end
      consistent = consistent.attach treatment_consistency, :complete => true
    end
    consistent
  end

  dep :consistent
  task :consistency_summary => :tsv do
    tsv = step(:consistent).load

    stats = TSV.setup({}, :key_field => "Experimen", :fields => %w(Consistent Inconsistent Proportion), :type => :list, :cast => :to_f)

    tsv.fields.each do |field|
      counts = Misc.counts(tsv.column(field).values.flatten.compact.reject{|v| v.empty?})
      counts["true"] ||= 0
      counts["false"] ||= 0

      stats[field] = counts["true"], counts["false"], counts["true"].to_f / (counts["true"] + counts["false"])
    end

    stats
  end
end
