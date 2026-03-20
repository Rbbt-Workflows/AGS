module AGS
  dep ExTRI2, :ExTRI2_clean, :only_authoritative_tfs => false
  task :coTF_sign => :tsv do
    coTFs = Rbbt.data.coTFs.list
    coTF_sign = TSV.setup({}, "CoTF~Sign#:type=:single")
    traverse step(:ExTRI2_clean) do |k,values|
      tf, tg,_, sign = values
      #next unless coTFs.include?(tf)
      next unless %w(UP DOWN).include? sign
      coTF_sign[tf] ||= 0
      coTF_sign[tf] += 1
    end
    coTF_sign
  end

  dep :coTF_sign
  extension :png
  task :coTF_plot => :binary do
    data = step(:coTF_sign).load
    data.R <<-EOF
    EOF
  end

  dep :treatment_tfs_dynamic, :ExTRI2_regulome => true, vetting: 'none' do |jobname,options|
    TREATMENTS.collect do |treatment|
      {:inputs => options.merge(treatment: treatment)}
    end
  end
  task :activities_for_signaling => :tsv do
    tsv = dependencies.inject(nil){|acc,dep| acc = acc.nil? ? dep.load : acc.attach(dep.path, complete: true) }
    tsv.subset(Rbbt.data["DNA_binding_and_some_likely_coTFs.13012025"].list) 
  end

  dep :timepoint_matrix, data_type: :range do |jobname,options|
    jobs = []
    AGS::TREATMENTS.each do |treatment|
      AGS::TIME_POINTS.each do |time_point|
        jobs << options.merge(treatment: treatment, time_point: time_point)
      end
    end
    jobs
  end
  task :range_matrix => :tsv do
    dependencies.inject(nil) do |acc,dep| 
      tsv = dep.load.transpose "Associated Gene Name"
      acc = acc.nil? ? tsv : acc.attach(tsv)
    end
  end


  dep :dbTFs
  dep :tf_predictions
  input :bona_fide, :boolean, "Use only bona fide TFs", false
  task :faster_in_combinations => :array do |bona_fide|
    dbTFs = step(:dbTFs).load
    tsv = step(:tf_predictions).load
    tsv.select do |gene,values|
      next if bona_fide and not dbTFs.include?(genes)
      up_in_combo = values.values_at("INT_PD_PI-T2",  "INT_PD_PI-T4",  "INT_FiveZ_PI-T2",  "INT_FiveZ_PI-T4").select{|v| v.to_f != 0 }.any?
      up_in_single = values.values_at(
        "PD-T2", 
        "PD-T4", 
        "FiveZ-T2", 
        "FiveZ-T4", 
      ).select{|v| v.to_f != 0 }.any?
      up_in_combo && ! up_in_single
    end.keys
  end
  
  dep :dbTFs
  dep :faster_in_combinations
  task_alias :faster_in_combination_enrichment, AGS, :gprofiler, list: :faster_in_combinations, background: :dbTFs, bona_fide: true

  dep :dbTFs
  dep :tf_predictions
  task :not_recovered_in_PD_PI => :array do |bona_fide|
    dbTFs = step(:dbTFs).load
    tsv = step(:tf_predictions).load
    other_fields = tsv.fields - ["INT_PD_PI-T24"]
    tsv.select do |gene,values|
      next if bona_fide and not dbTFs.include?(genes)
      up_in_combo = values.values_at("INT_PD_PI-T24").select{|v| v.to_f != 0 }.any?
      up_in_others = values.values_at(*other_fields).select{|v| v.to_f != 0 }.any?
      ! up_in_combo && up_in_others
    end.keys
  end

  dep :dbTFs
  dep :not_recovered_in_PD_PI
  task_alias :not_recovered_in_PD_PI_enrichment, AGS, :gprofiler, list: :not_recovered_in_PD_PI, background: :dbTFs, bona_fide: true

  dep :tf_predictions
  task :statistics => :yaml do
    pred = step(:tf_predictions).load
    
    res = {}
    res[:total] = pred.keys.length
    res[:bona_fide] = (AGS::BONAFIDE_TFS & pred.keys).length
    
    pred.fields.each do |field|
      res["#{field}_total"] = pred.select(field){|v| v  }.keys.length
      res["#{field}_bona_fide"] = (AGS::BONAFIDE_TFS & pred.select(field){|v| v }.keys).length
    end

    res
  end

  desc "Return a list of the top derregulated transcription factor activities inferred from the gene expression levels of a particular treatment at each of the different timepoints"
  dep :treatment_tfs
  input :time_point, :integer, "Time point", nil, required: true
  input :direction, :select, "Activity regulation direction, up or down", nil, select_options: %w(up down), required: true
  input :max, :integer, "Maximum number of top derregulated transcription factors returned", 100
  task :list_tfs => :array do |time_point,direction,max|
    tsv = step(:treatment_tfs).load
    field = tsv.fields.select{|field| field.end_with?("-T#{time_point}") }.first
    case direction.to_s
    when 'up'
      tsv = tsv.select(field){|v| v.to_f > 0 }
    when 'down'
      tsv = tsv.select(field){|v| v.to_f < 0 }
    end

    keys = tsv.column(field).to_single.sort_by{|k,v| v.to_f.abs }.reverse

    keys[0..max.to_i].collect{|k,v| k }
  end

  dep :list_tfs
  dep :change_offsets_simplified
  input :tfs, :array, "Restrict genes to targets of these transcription factors", nil
  desc "Return a list of the target genes that are starting to be derragulated in expression, up or down, in a particular treatment at each of the different timepoints"
  task :list_tgs => :array do |tfs|
    time_point, direction, max = self.recursive_inputs.values_at :time_point, :direction, :max

    tfs  ||= step(:list_tfs).load
    treatment = recursive_inputs[:treatment]
    change_offsets = step(:change_offsets_simplified).load
    regulome = step(:regulome).path.tsv key_field: 'source', fields: ['target'], type: :flat

    target_genes = regulome.values_at(*tfs).flatten.uniq

    case direction.to_s
    when 'up'
      tsv = change_offsets.select(treatment){|v| v.include? "increase #{time_point}h" }
    when 'down'
      tsv = change_offsets.select(treatment){|v| v.include? "decrease #{time_point}h" }
    else
      raise
    end

    tsv.subset(target_genes).keys.uniq
  end

  dep :regulome
  input :gene, :string, "Transcription factor"
  task :tf_targets => :tsv do |gene|
    regulome = step(:regulome).path.tsv key_field: "source", fields: ['target', 'weight'], type: :double
    tsv = regulome.subset([gene]).reorder 'target', ['weight'], one2one: true
    tsv = tsv.add_field "Mode of regulation" do |k,values|
      values.last.last.to_f < 0 ? "inhibit" : "activate"
    end

    tsv.reorder('target', ['Mode of regulation']).to_single
  end

  dep :tf_targets
  dep :list_tgs, direction: "up"
  dep :list_tgs, direction: "down"
  task :target_analysis => :tsv do |tfs|
    targets, up, down = dependencies.collect{|dep| dep.load }
    targets = targets.to_list

    targets.add_field "Regulated" do |tg,values|
      if up.include? tg
        "up"
      elsif up.include? tg
        "down"
      else
        ""
      end
    end

    targets = targets.select("Regulated"){|v| not (v.nil? or v.empty?) }

    targets.add_field "Consistent" do |tg,values|
      mor, direction = values
      if mor == 'activate' 
        if direction == 'up'
          "consistent"
        else
          "inconsistent"
        end
      else
        if direction == 'up'
          "inconsistent"
        else
          "consistent"
        end
      end
    end
  end

  dep :change_offsets_simplified
  task :gprofiler_queries => :string do
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

          name = [treatment, timepoint, type, direction] * "_"
          file(name + '.txt').write <<-EOF
>#{name}
#{up_cluster*"\n"}
          EOF
        end
      end
    end
    'DONE'
  end

end

