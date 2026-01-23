module AGS

  helper :prev_time_point do |time_point|
    index = TIME_POINTS.index time_point
    if index > 0
      TIME_POINTS[index-1]
    else
      "NO TIMEPOINT"
    end
  end

  dep :treatment_tfs
  dep :change_offsets
  dep :dbTFs
  input :bona_fide, :boolean, "Use only bona fide TFs", false
  input :remove_consecutive, :boolean, "Remove consecutive changes", true
  task :treatment_tf_consistency => :tsv do |bona_fide,remove_consecutive|
    dbTFs = step(:dbTFs).load
    activities = dependencies.first
    treatment = recursive_inputs[:treatment]
    changes = step(:change_offsets).load
    fields = AGS::TIME_POINTS.collect{|t| "Consistent at #{t}h" }
    consistency = TSV.setup({}, key_field: "Associated Gene Name", fields: fields, type: :list)
    traverse activities, into: consistency do |gene,values|
      next unless changes.include?(gene)
      next if bona_fide and not dbTFs.include?(gene)

      gene_changes = changes[gene][treatment]
      #negative_activities = gene_changes.select{|c| c.include?("decrease")}.
      #  collect{|c| c.split(" ").last.to_i}.collect{|t| AGS::TIME_POINTS.index(t)}.
      #  collect{|i| [i, i+1]}.flatten.reject{|i| i > 5}.uniq.
      #  collect{|i| AGS::TIME_POINTS[i] }

      #positive_activities = gene_changes.select{|c| c.include?("increase")}.
      #  collect{|c| c.split(" ").last.to_i}.collect{|t| AGS::TIME_POINTS.index(t)}.
      #  collect{|i| [i, i+1]}.flatten.reject{|i| i > 5}.uniq.
      #  collect{|i| AGS::TIME_POINTS[i] }

      negative_changes = gene_changes.select{|c| c.include?("decrease")}.
        collect{|c| c.split(" ").last.to_i}.collect{|t| AGS::TIME_POINTS.index(t)}

      positive_changes = gene_changes.select{|c| c.include?("increase")}.
        collect{|c| c.split(" ").last.to_i}.collect{|t| AGS::TIME_POINTS.index(t)}

      if remove_consecutive
        negative_changes = negative_changes.uniq.sort.each_with_object([]) {|n, res| 
          res << n unless negative_changes.include?(n-1) 
        }
        positive_changes = positive_changes.uniq.sort.each_with_object([]) {|n, res| 
          res << n unless positive_changes.include?(n-1)
        }
      end

      #res = []
      #AGS::TIME_POINTS.each_with_index do |time_point,i|
      #  res << (values[i].to_f > 0 && positive_activities.include?(time_point)) || (values[i].to_f < 0 && negative_activities.include?(time_point))
      #end
      #res = res.collect{|c| c ? 1 : 0 }

      res = []
      AGS::TIME_POINTS.each_with_index do |time_point,i|
        res << begin
                 if values[i].to_f > 0
                   neg = negative_changes.reject{|v| v >= i }.max || -1
                   pos = positive_changes.reject{|v| v > i }.max || -1

                   if pos > neg
                     if i >= 4 && pos <= 1
                       0
                     elsif i >= 3 && pos == 0
                       0
                     else
                       1
                     end
                   elsif neg == pos
                     0
                   else
                     -1
                   end
                 elsif values[i].to_f < 0
                   neg = negative_changes.reject{|v| v > i }.max || -1
                   pos = positive_changes.reject{|v| v >= i }.max || -1

                   if neg > pos
                     if i >= 4 && neg <= 1
                       0
                     elsif i >= 3 && neg == 0
                       0
                     else
                       1
                     end
                   elsif neg == pos
                     0
                   else
                     -1
                   end
                 else
                   nil
                 end
        end
      end

      [gene, res]
    end
    tsv = activities.load.attach consistency
    tsv.cast = nil
    tsv
  end

  dep :dbTFs
  dep :treatment_tf_consistency, treatment: :placeholder do |jobname,options|
    AGS::TREATMENTS.collect do |treatment|
      next if treatment == "DMSO"
      options.merge(treatment: treatment)
    end.compact
  end
  task :consistency_summary => :float do
    consistency_tfs = step(:dbTFs).load
    dep_consistency = dependencies[1..-1].collect do |dep| 
      tsv = dep.load
      tsv = tsv.subset(consistency_tfs)
      fields = tsv.fields.select{|f| f.start_with?("Consistent") } 
      Misc.mean(fields.collect{|field| Misc.sum(tsv.column(field).values.flatten.compact.reject{|v| v.empty? || v.to_f <= 0 }.map(&:to_f)) }.reject{|v| v.nan? })
    end
    Misc.mean(dep_consistency)
  end

  dep :treatment_tf_consistency, treatment: :placeholder, scheme: :placeholder, vetting: :placeholder do |jobname,options|
    %w(dynamic non-dynamic fc0).collect do |scheme|
      %w(none relaxed_degradation).collect do |vetting|
        AGS::TREATMENTS.collect do |treatment|
          next if treatment == "DMSO"
          options.merge(treatment: treatment, scheme: scheme, vetting: vetting)
        end.compact
      end
    end.flatten
  end
  task :consistency_sweep => :tsv do
    dep_consistencies = dependencies.collect do |dep| 
      tsv = dep.load
      treatment = dep.recursive_inputs[:treatment]
      scheme = dep.recursive_inputs[:scheme]
      vetting = dep.recursive_inputs[:vetting]
      tsv.fields = tsv.fields.collect{|f|
        f = treatment + "." + f unless f.include? treatment
        f = [f, scheme, vetting] * "."
        f
      }
      tsv
    end
    dep_consistencies.inject(nil){|acc,tsv| acc = acc.nil? ? tsv : acc.attach(tsv, complete: true) }
  end

  dep :consistency_sweep
  task :consistency_counts => :tsv do
    data = step(:consistency_sweep).load
    consistency_fields = data.fields.select{|field| field.include? 'Consistent' }

    tsv = TSV.setup({}, "ID~Treatment,Time,Scheme,Vetting,Matches,Miss,Total,Match/Miss Odds,Match/Total Odds")
    id = 1
    consistency_fields.each do |field|
      treatment, title, scheme, vetting = field.split('.')
      time = title.split(' ').last.to_i
      counts = Misc.counts(data.column(field).values)
      matches = counts['1'] || 0
      miss = counts['-1'] || 0
      zero = counts['0'] || 0
      total = matches+miss+zero
      tsv[id] = [treatment, time, scheme, vetting, matches, miss, total, matches.to_f/miss, matches.to_f/total]
      id += 1
    end

    tsv
  end



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

  #dep :timepoint_decoupler_suite
  #dep :change_offsets, :jobname => "Default"
  #task :consistent => :tsv do 
  #  suite = step(:timepoint_decoupler_suite).load
  #  clusters = step(:change_offsets).load
  #  consistent = TSV.setup({}, :key_field => "Associated Gene Name", :fields => [], :type => :list)

  #  fields = suite.fields
  #  TSV.traverse fields, :bar => true do |field|
  #    treatment, time_point = field.split("-") 
  #    time_point = time_point.to_i
  #    prev_time_point = prev_time_point(time_point)
  #    treatment_clusters = clusters.column(treatment).to_flat

  #    treatment_consistency = TSV.setup({}, :key_field => consistent.key_field, :fields => [field], :type => :single)

  #    # ToDo: check the start_with
  #    suite.through nil, [field] do |gene,vs|
  #      v = vs.first.to_f
  #      gene_clusters = treatment_clusters[gene] || []
  #      consistency = if v.to_f == 0.0
  #                      next
  #                    elsif v.to_f > 0
  #                      gene_clusters.include?("increase #{time_point}h") || 
  #                        gene_clusters.include?("increase #{prev_time_point}h") ||
  #                        gene_clusters.start_with?("increase #{time_point}-") ||
  #                        gene_clusters.start_with?("increase #{prev_time_point}-")
  #                    else
  #                      gene_clusters.include?("decrease #{time_point}h") || 
  #                        gene_clusters.include?("decrease #{prev_time_point}h") ||
  #                        gene_clusters.start_with?("decrease #{time_point}-") || 
  #                        gene_clusters.start_with?("decrease #{prev_time_point}-")
  #                    end
  #      treatment_consistency[gene] = consistency
  #    end
  #    consistent = consistent.attach treatment_consistency, :complete => true
  #  end
  #  consistent
  #end

  #dep :consistent
  #task :consistency_summary => :tsv do
  #  tsv = step(:consistent).load

  #  stats = TSV.setup({}, :key_field => "Experimen", :fields => %w(Consistent Inconsistent Proportion), :type => :list, :cast => :to_f)

  #  tsv.fields.each do |field|
  #    counts = Misc.counts(tsv.column(field).values.flatten.compact.reject{|v| v.empty?})
  #    counts["true"] ||= 0
  #    counts["false"] ||= 0

  #    stats[field] = counts["true"], counts["false"], counts["true"].to_f / (counts["true"] + counts["false"])
  #  end

  #  stats
  #end
end
