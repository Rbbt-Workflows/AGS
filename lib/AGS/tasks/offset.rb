module AGS

  input :pvalue_threshold, :float, "Threshold for decoupler p-value", 0.01
  input :reported_method, :select, "Method to report", :prev_clust, :select_options => %w(best prev_clust prev_all diff)
  dep :timepoint_decoupler,
    :treatment => :placeholder, :time_point => :placeholder,
    :threshold => :placeholder, :data_type => :placeholder, 
    :dynamic => :placeholder  do |jobname,options|
      jobs = []
      pvalue_threshold = options[:pvalue_threshold]
      AGS::TREATMENTS.each do |treatment|
        AGS::TIME_POINTS.each_with_index do |time_point,index|
          jobs << AGS.job(:timepoint_decoupler, treatment: treatment, time_point: time_point, data_type: :fc0, :threshold => pvalue_threshold)
          jobs << AGS.job(:timepoint_decoupler, treatment: treatment, time_point: AGS::TIME_POINTS[index - 1], data_type: :fc0, :threshold => pvalue_threshold)
          jobs << AGS.job(:timepoint_decoupler, treatment: treatment, time_point: time_point, data_type: :fc, :threshold => pvalue_threshold)
          jobs << AGS.job(:timepoint_decoupler, treatment: treatment, time_point: time_point, data_type: :range, dynamic: true, :threshold => pvalue_threshold)
        end
      end
      jobs.uniq
    end
  task :tf_offsets => :array do |pvalue_threshold, reported_method|
    fc0_all = nil
    fc_all = nil
    fcf_all = nil

    AGS::TREATMENTS.each do |treatment|
      AGS::TIME_POINTS.each_with_index do |time_point,index|
        fc0 = AGS.job(:timepoint_decoupler, treatment: treatment, time_point: time_point, data_type: :fc0, :threshold => pvalue_threshold).run
        if index > 0
          fc0_prev = AGS.job(:timepoint_decoupler, treatment: treatment, time_point: AGS::TIME_POINTS[index - 1], data_type: :fc0, :threshold => pvalue_threshold).run
          fc0.fields do |field|
            fc0.process field do |k,v|
              v.to_f - fc0_prev[k][field].to_f
            end
          end
        end

        fc = AGS.job(:timepoint_decoupler, treatment: treatment, time_point: time_point, data_type: :fc, :threshold => pvalue_threshold).run

        fcf = AGS.job(:timepoint_decoupler, treatment: treatment, time_point: time_point, data_type: :range, dynamic: true, :threshold => pvalue_threshold).run

        fc_all = fc_all.nil? ? fc : fc_all.attach(fc, :complete => true)
        fcf_all = fcf_all.nil? ? fcf : fcf_all.attach(fcf, :complete => true)
        fc0_all = fc0_all.nil? ? fc0 : fc0_all.attach(fc0, :complete => true)
      end
    end

    AGS::TREATMENTS.each do |treatment|
      tpfc = fc_all.slice(fc_all.fields.select{|f| f.split("-").include?(treatment) })
      tpfcf = fcf_all.slice(fcf_all.fields.select{|f| f.split("-").include?(treatment) })
      tpfc0 = fc0_all.slice(fc0_all.fields.select{|f| f.split("-").include?(treatment) })

      tpfc = tpfc.select{|k,l| l.select{|v| v.to_f != 0.0 }.any? }
      tpfcf = tpfcf.select{|k,l| l.select{|v| v.to_f != 0.0 }.any? }
      tpfc0 = tpfc0.select{|k,l| l.select{|v| v.to_f != 0.0 }.any? }

      tpfc.fields = tpfc.fields.collect{|f| f + ' prev all' } 
      tpfcf.fields = tpfcf.fields.collect{|f| f + ' prev clust' }
      tpfc0.fields = tpfc0.fields.collect{|f| f + ' diff' }

      data = tpfcf.attach tpfc, complete: true
      data = data.attach tpfc0, complete: true

      gene_pos_offsets = {}
      gene_neg_offsets = {}
      data.through do |gene,values|
        prev = 0
        values.each_with_index do |v,i|
          if ! v.nil?
            if v > 0 && prev <= 0
              gene_pos_offsets[gene] ||= []
              gene_pos_offsets[gene] << i
            end
            if v < 0 && prev >= 0
              gene_neg_offsets[gene] ||= []
              gene_neg_offsets[gene] << i
            end
          end
          prev = (((i + 1) % 5 == 0) || i < 5) ? 0 : v.to_f
        end
      end

      genes = (gene_neg_offsets.keys + gene_pos_offsets.keys).uniq

      gene_timepoint_info = {}
      genes.each do |gene|
        offsets = (gene_pos_offsets[gene] || []) + (gene_neg_offsets[gene] || [])
        offset_values = data[gene].values_at(*offsets)
        offset_names = data.fields.values_at(*offsets)
        offset_time_points = offset_names.collect{|n| n.match(/T(\d+)/)[1] }
        offset_names.zip(offset_time_points, offset_values).each do |n, t, v|
          method = n.partition(" ").last
          gene_timepoint_info[gene] ||= {}
          gene_timepoint_info[gene][t] ||= []
          gene_timepoint_info[gene][t] << {method: method, value: v}
        end
      end

      AGS::TIME_POINTS.each do |time_point|
        tsv = TSV.setup({}, :key_field =>"Transcription factor (Associated Gene Name)", :fields => %w(Value Best All Conflict?), :type => :list)
        gene_timepoint_info.each do |gene,info|
          list = info[time_point.to_s]
          next unless list && list.any?
          all_values = list.collect{|e| e[:value]}
          all = list.collect{|e| e[:method] }
          case reported_method.to_sym
          when :best
            max = all_values.max
            best = list.select{|e| e[:value] == max}.first[:method]
          when :prev_clust
            next unless all.include?("prev clust")
            max = list.select{|e| e[:method] == "prev clust"}.first[:value]
            best = reported_method
          when :prev_all
            next unless all.include?("prev all")
            max = list.select{|e| e[:method] == "prev all"}.first[:value]
            best = reported_method
          when :diff
            next unless all.include?("diff")
            max = list.select{|e| e[:method] == "diff"}.first[:value]
            best = reported_method
          end
          confict = all_values.select{|v| v > 0 }.any? && all_values.select{|v| v < 0 }.any?
          tsv[gene] = [max, best, all * "|", confict]
        end
        file([treatment, time_point] * "-T" + ".tsv").write(tsv.to_s)
      end
    end
    files_dir.glob("*.tsv")
  end

  dep :tf_offsets
  input :treatment, :select, "Treatment", nil, :select_options => TREATMENTS, :required => true
  task :tf_treatment_offsets => :tsv do |treatment|
    tf_offset_job = step(:tf_offsets)

    AGS::TIME_POINTS.inject(nil) do |acc,timepoint|
      tsv = tf_offset_job.file([treatment, timepoint] * "-T" + ".tsv").tsv :fields => %w(Value)
      tsv.fields = ["T" + timepoint.to_s]
      acc = acc.nil? ? tsv : acc.attach(tsv)
    end
  end
end
