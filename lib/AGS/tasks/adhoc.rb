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
rbbt.png.plot("de
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

  dep :treatment_tfs
  input :time_point, :string, nil, 1
  input :direction, :select, nil, :up, select_options: %w(up down)
  input :max, :integer, nil, 20, select_options: %w(up down)
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
end

