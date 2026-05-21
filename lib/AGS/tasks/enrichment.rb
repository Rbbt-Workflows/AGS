require 'rbbt/statistics/hypergeometric'

module AGS

  dep :expressed_coding_genes, jobname: "Default"
  input :list, :array, "Gene set"
  input :background, :array, "Background genes", nil
  input :database, :select, "Annotation to find enrichment", :go_bp, select_options: %w(go_bp)
  task :gs_hyper => :tsv do |list,background,database|
    background = step(:expressed_coding_genes).load if background.nil?

    tsv = case database.to_sym
          when :go_bp
            Organism.gene_go_bp(AGS.organism).tsv type: :flat
          else
            raise ParameterException, "Unkown database parameter #{database}"
          end

    tsv = tsv.change_key "Associated Gene Name"

    tsv.enrichment list, nil, background: background
  end

  dep :expressed_coding_genes, jobname: "Default"
  input :list, :array, "Gene set"
  input :background, :array, "Background genes", nil
  input :database, :select, "Annotation to find enrichment", :go_bp, select_options: %w(go_bp)
  task :gprofiler => :tsv do |list,background,database|
    require 'rbbt/util/python'
    background = step(:expressed_coding_genes).load if background.nil?

    list = list.load if Step === list
    list.shift if list[1] && list[1].start_with?(">")
    RbbtPython.run 'gprofiler' do
      gp = gprofiler.GProfiler.new(return_dataframe:true)
      if background && background.any?
        res = gp.profile(organism:'hsapiens', query: list, sources: %w(GO:BP), background: background)
      else
        res = gp.profile(organism:'hsapiens', query: list, sources: %w(GO:BP))
      end
      tsv = RbbtPython.df2tsv res
      tsv.type = :list
      tsv.key_field = "Position"
      tsv = tsv.reorder "native"
      tsv.fields = tsv.fields.collect{|f| f == "p_value" ? "p-value" : f }
      tsv.each do |k,v|
        v[13] = v[13..-1] * "|"
        v.replace v.slice(0, 14)
      end
      tsv
    end
  end

  dep :expressed_coding_genes
  input :queries, :yaml, "Gene set"
  input :background, :array, "Background genes", nil
  input :database, :select, "Annotation to find enrichment", :go_bp, select_options: %w(go_bp)
  task :gprofiler_multiple => :tsv do |queries,background,database|
    require 'rbbt/util/python'
    background = step(:expressed_coding_genes).load if background.nil?

    queries = queries.select{|k,v|
      k.include? "_both"
    }
    RbbtPython.run 'gprofiler' do
      gp = gprofiler.GProfiler.new(return_dataframe:true)
      if background && background.any?
        #res = gp.profile(organism:'hsapiens', query: queries, sources: %w(GO:BP), background: background)
        res = gp.profile(organism:'hsapiens', query: queries, sources: %w(GO:BP), background: background)
      else
        res = gp.profile(organism:'hsapiens', query: queries, sources: %w(GO:BP))
      end
      tsv = RbbtPython.df2tsv res
      tsv = tsv.to_double
      tsv.key_field = "Position"
      tsv.fields = tsv.fields.collect{|f| f == "p_value" ? "p-value" : f }
      tsv.each do |k,v|
        v[13] = v[13..-1].flatten
        v.replace v.slice(0, 14)
      end
      tsv
    end
  end

  dep :change_offsets_simplified
  dep :tf_predictions
  task :gprofiler_queries => :array do
    change_offsets = step(:change_offsets_simplified).load
    AGS::TIME_POINTS.each do |timepoint|
      AGS::TREATMENTS.each do |treatment|
        up_cluster = change_offsets.select(treatment => "increase #{timepoint}h").keys
        down_cluster = change_offsets.select(treatment => "decrease #{timepoint}h").keys

        tsv = AGS.job(:timepoint_matrix, treatment: treatment, time_point: timepoint, data_type: :fc).run.transpose
        fc_up_03 = tsv.select do |k,values|
          values.flatten.first.to_f > 0.3
        end.keys
        fc_down_03 = tsv.select do |k,values|
          values.flatten.first.to_f < -0.3
        end.keys

        fc_up_07 = tsv.select do |k,values|
          values.flatten.first.to_f > 0.7
        end.keys
        fc_down_07 = tsv.select do |k,values|
          values.flatten.first.to_f < -0.7
        end.keys


        tsv = AGS.job(:timepoint_matrix, treatment: treatment, time_point: timepoint, data_type: :fc0).run.transpose
        fc0_up_03 = tsv.select do |k,values|
          values.flatten.first.to_f > 0.3
        end.keys
        fc0_down_03 = tsv.select do |k,values|
          values.flatten.first.to_f < -0.3
        end.keys

        fc0_up_07 = tsv.select do |k,values|
          values.flatten.first.to_f > 0.7
        end.keys
        fc0_down_07 = tsv.select do |k,values|
          values.flatten.first.to_f < -0.7
        end.keys

        [
          [up_cluster, 'cluster', 'up'],
          [down_cluster, 'cluster', 'down'],
          [fc_up_03, 'fc_03', 'up'],
          [fc_down_03, 'fc_03', 'down'],
          [fc_up_07, 'fc_07', 'up'],
          [fc_down_07, 'fc_07', 'down'],
          [fc0_up_03, 'fc0_03', 'up'],
          [fc0_down_03, 'fc0_03', 'down'],
          [fc0_up_07, 'fc0_07', 'up'],
          [fc0_down_07, 'fc0_07', 'down'],
        ].each do |list,type,direction|

          name = [treatment, timepoint.to_s+'h', type, direction] * "-"
          file(name + '.txt').write <<-EOF
>#{name}
#{list*"\n"}
          EOF
        end
      end
    end

    predictions = step(:tf_predictions).load
    
    predictions.fields.each do |treatment_tp|
      treatment, timepoint = treatment_tp.split('-')
      timepoint = timepoint.sub('T','').to_s + 'h'
      up = predictions.select(treatment_tp){|v| v.to_f > 0 }.keys
      down = predictions.select(treatment_tp){|v| v.to_f < 0 }.keys

      up_name = [treatment, timepoint, 'TF', 'up'] * '-'
      down_name = [treatment, timepoint, 'TF', 'down'] * '-'
      file(up_name + '.txt').write <<-EOF
>#{up_name}
#{up*"\n"}
      EOF

      file(down_name + '.txt').write <<-EOF
>#{down_name}
#{down*"\n"}
      EOF
    end

    files
  end

  dep :gprofiler_queries, compute: :produce
  dep :gprofiler do |jobname,options,dependencies|
    queries = dependencies.flatten.first
    queries.files.collect do |file|
      options.merge(list: queries.file(file).list, jobname: file.sub('.txt', ''))
    end
  end
  task :gprofiler_suite => :tsv do
    tsv = TSV.setup({}, key_field: "Treatment:Time", fields: ["Up", "Down"], type: :double)
    dependencies[1..-1].collect do |dep|
      name = dep.clean_name
      Open.cp dep.path, file(name + '.tsv')
      next unless name.include?("cluster")
      treatment, hour, type, direction = name.split("-")
      tp = [treatment, hour + "h"] * ":"
      dep.load.each do |id, values|
        name = values[2]
        pvalue = values[3]
        next unless pvalue.to_f < 0.05
        tsv[tp] ||= [[], []]
        case direction
        when "up"
          tsv[tp][0] << name 
        when "down"
          tsv[tp][1] << name 
        else
          next
        end
      end
    end
    tsv
  end

end
