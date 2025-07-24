module AGS
  

  dep :treatment_tfs
  dep :dbTFs
  input :dbTFs, :boolean, "Use only  dbTFs", true
  task :neko_consistency => :tsv do |use_dbtf|
    dbTFs = step(:dbTFs).load


    treatment = recursive_inputs[:treatment]
    predicted = case treatment
                when "PI"
                  Rbbt.data.neko_benchmarks['PI_benchmark_020425'].tsv type: :list 
                when "PD"
                  Rbbt.data.neko_benchmarks['PD_benchmark_020425'].tsv  type: :list
                when "FiveZ"
                  Rbbt.data.neko_benchmarks['5Z_benchmark_020425'].tsv  type: :list
                else
                  raise ParameterException, "Treatment not supported #{treatment}"
                end

    tfs = step(:treatment_tfs).load

    tsv = TSV.setup({}, "TF~Type,Effect,Predicted,T1 Activity,T2 Activity,T1 Match,T2 Match#:type=:list")
    tfs.each do |tf,timepoints|
      next if use_dbtf and not dbTFs.include?(tf)
      t1, t2 = timepoints
      p = predicted[tf]
      next if p.nil?
      type, target, effect, value = p

      value = - value.to_i

      t1_match = (value <=> 0) == (t1.to_f <=> 0) if t1
      t2_match = (value <=> 0) == (t2.to_f <=> 0) if t2

      next unless t1_match || t2_match
      tsv[tf] = [type, effect, value, t1, t2, t1_match, t2_match]
    end

    tsv
  end
  
  dep :neko_consistency
  input :target, :select, 'What to compare', :strict, select_options: %w(strict relaxed T1 T2)
  task :neko_summary => :array do |target|
    tsv = step(:neko_consistency).load
    
    match = 0.0
    miss_match = 0.0

    tsv.each do |tf,values|
      type, effect, predicted, t1, t2, m1, m2 = values
      m1 = nil if m1 && m1.empty?
      m2 = nil if m2 && m2.empty?
      
      m1 = m1.to_s == "true" if m1
      m2 = m2.to_s == "true" if m2

      case target.to_sym
      when :relaxed
        if m1 && m2.nil?
          match += 1
        elsif m2 && m1.nil?
          match += 1
        elsif m1 && m2
          match += 1
          match += 1
        else
          match += 1 if m1
          match += 1 if m2

          miss_match += 1 if !m1
          miss_match += 1 if !m2
        end
      when :strict
        match += 1 if m1
        match += 1 if m2

        miss_match += 1 if !m1
        miss_match += 1 if !m2
      when :T1
        match += 1 if m1
        miss_match += 1 if !m1
      when :T2
        match += 1 if m2
        miss_match += 1 if !m2
      end
    end

    [match, miss_match]
  end

  dep :neko_summary, treatment: :placeholder, scheme: :placeholder, vetting: :placeholder do |jobname,options|
    %w(strict relaxed T1 T2).collect do |target|
      %w(dynamic non-dynamic).collect do |scheme|
        %w(none relaxed_degradation).collect do |vetting|
          %w(PD PI FiveZ).collect do |treatment|
            options.merge(treatment: treatment, scheme: scheme, vetting: vetting, target: target)
          end.compact
        end
      end
    end.flatten
  end
  task :neko_sweep => :tsv do
    tsv = TSV.setup({}, "ID~Treatment,Scheme,Vetting,Target,Match,Miss,Total,Odds")
    id = 0
    dependencies.each do |dep| 
      match, miss_match = dep.load
      total = match.to_f + miss_match.to_f
      odds = match.to_f / total
      treatment = dep.recursive_inputs[:treatment]
      scheme = dep.recursive_inputs[:scheme]
      vetting = dep.recursive_inputs[:vetting]
      target = dep.recursive_inputs[:target]
      tsv[id]=[treatment, scheme, vetting, target, match, miss_match, total, odds]
      id += 1 
    end
    tsv
  end

  input :treatment, :select, "Treatment", nil, :select_options => %w(PI PD FiveZ PD-MAP2K2 PD-MAP2K1), :required => true
  task :neko_bootstrap => :tsv do |treatment|
    treatment = treatment == "FiveZ" ? "5Z" : treatment
    files = Rbbt.modules['TFe-nekofinder']["bootstrap_#{treatment.split('-').first}_test"].glob("*.csv")

    main_targets = {
      "PD-MAP2K2" => ['MAP2K2'],
      "PD-MAP2K1" => ['MAP2K1'],
      "PD" => ['MAP2K2', 'MAP2K1'],
      "PI" => ['PIK3CA'],
      "5Z" => ['MAP3K7']
    }
    
    main = main_targets[treatment]

    res = files.inject({}) do |acc,file|
      tsv = TSV.csv file, merge: true, type: :double

      tsv.each do |target,values|
        Misc.zip_fields(values).each do |tf,sign|
          acc[tf] ||= []
          acc[tf] << [target, sign]
        end
      end
      acc
    end

    new = {}
    res.each do |tf,pairs|
      main_matches = pairs.select{|target,sign| main.include?(target)}
      if true # We are only looking at main targets
        new[tf] = main_matches.collect{|target,sign| sign }
      else
        new[tf] = pairs.collect{|target,sign| sign }
      end
    end

    TSV.setup(new, key_field: "Target", fields: ["Signs"], type: :flat)
  end

  dep :neko_consistency
  input :target, :select, 'What to compare', :strict, select_options: %w(strict relaxed T1 T2)
  task :neko_summary => :array do |target|
    tsv = step(:neko_consistency).load
    
    match = 0.0
    miss_match = 0.0

    tsv.each do |tf,values|
      type, effect, predicted, t1, t2, m1, m2 = values
      m1 = nil if m1 && m1.empty?
      m2 = nil if m2 && m2.empty?
      
      m1 = m1.to_s == "true" if m1
      m2 = m2.to_s == "true" if m2

      case target.to_sym
      when :relaxed
        if m1 && m2.nil?
          match += 1
        elsif m2 && m1.nil?
          match += 1
        elsif m1 && m2
          match += 1
          match += 1
        else
          match += 1 if m1
          match += 1 if m2

          miss_match += 1 if !m1
          miss_match += 1 if !m2
        end
      when :strict
        match += 1 if m1
        match += 1 if m2

        miss_match += 1 if !m1
        miss_match += 1 if !m2
      when :T1
        match += 1 if m1
        miss_match += 1 if !m1
      when :T2
        match += 1 if m2
        miss_match += 1 if !m2
      end
    end

    [match, miss_match]
  end

  dep :neko_bootstrap
  task :neko_bootstrap_overview => :tsv do
    tsv = step(:neko_bootstrap).load
    tsv = tsv.to_double
    
    size = tsv.column("Signs").values.collect{|l| l.length }.max

    min = size * 0.8

    tsv = tsv.select do |k,values|
      values["Signs"].length > min
    end

    tsv.add_field "Positive" do |k,values|
      Misc.counts(values["Signs"])['1'].to_f / values["Signs"].length
    end

    tsv.add_field "Negative" do |k,values|
      Misc.counts(values["Signs"])['-1'].to_f / values["Signs"].length
    end

    tsv
  end

  dep :neko_bootstrap_overview
  task :neko_bootstrap_final => :tsv do
    tsv = step(:neko_bootstrap_overview).load

    final = TSV.setup({}, key_field: "TF", fields: ["Sign"], type: :single)

    tsv.through do |tf,values|
      pos = values["Positive"].first.to_f
      if pos >= 0.8
        final[tf] = 1
      elsif pos <= 0.2
        final[tf] = -1
      else
        next
      end
    end

    final
  end

  input :treatment, :select, "Treatment to compare, PD is separated in two for convenience", nil, select_options: %w(PI PD 5Z PD-MAP2K1 PD-MAP2K2) 
  dep :treatment_tfs do |jobname,options,dependencies|
    options.merge(treatment: options[:treatment].split("-").first)
  end
  dep :dbTFs
  dep :neko_bootstrap_final
  input :dbTFs, :boolean, "Use only  dbTFs", true
  task :neko_bootstrap_consistency => :tsv do |use_dbtf|
    dbTFs = step(:dbTFs).load
    bootstrap = step(:neko_bootstrap_final).load

    tfs = step(:treatment_tfs).load

    tsv = TSV.setup({}, "TF~Predicted,T1 Activity,T2 Activity,T1 Match,T2 Match#:type=:list")
    tfs.each do |tf,timepoints|
      next if use_dbtf and not dbTFs.include?(tf)
      t1, t2 = timepoints
      value = bootstrap[tf]
      next if value.nil?
      next if t1.nil? && t2.nil?

      value = - value.to_i

      t1_match = (value <=> 0) == (t1.to_f <=> 0) if t1
      t2_match = (value <=> 0) == (t2.to_f <=> 0) if t2

      #next unless t1_match || t2_match
      tsv[tf] = [value, t1, t2, t1_match, t2_match]
    end

    tsv
  end

  dep :neko_bootstrap_consistency
  input :target, :select, 'What to compare', :relaxed, select_options: %w(strict relaxed T1 T2)
  task :neko_bootstrap_summary => :array do |target|
    tsv = step(:neko_bootstrap_consistency).load
    
    match = 0
    miss_match = 0

    tsv.each do |tf,values|
      predicted, t1, t2, m1, m2 = values
      m1 = nil if m1 && m1.empty?
      m2 = nil if m2 && m2.empty?
      
      m1 = m1.to_s == "true" if m1
      m2 = m2.to_s == "true" if m2

      case target.to_sym
      when :relaxed
        if m1 || m2
          match += 1
        else
          miss_match += 1
        end
      when :strict
        match += 1 if m1
        match += 1 if m2

        miss_match += 1 if !m1
        miss_match += 1 if !m2
      when :T1
        next if m1.nil?
        match += 1 if m1
        miss_match += 1 if !m1
      when :T2
        next if m2.nil?
        match += 1 if m2
        miss_match += 1 if !m2
      when :T1_old
        match += 1 if m1
        miss_match += 1 if !m1
      when :T2_old
        match += 1 if m2
        miss_match += 1 if !m2
      when :relaxed_old
        if m1 && m2.nil?
          match += 1
        elsif m2 && m1.nil?
          match += 1
        elsif m1 && m2
          match += 1
          match += 1
        else
          match += 1 if m1
          match += 1 if m2

          miss_match += 1 if !m1
          miss_match += 1 if !m2
        end
      end
    end

    [match, miss_match]
  end

  dep :neko_bootstrap_summary, treatment: :placeholder, scheme: :placeholder, vetting: :placeholder do |jobname,options|
    %w(relaxed T1 T2).collect do |target|
      %w(dynamic non-dynamic).collect do |scheme|
        %w(none relaxed_degradation).collect do |vetting|
          %w(PD PI FiveZ PD-MAP2K1 PD-MAP2K2).collect do |treatment|
            options.merge(treatment: treatment, scheme: scheme, vetting: vetting, target: target)
          end.compact
        end
      end
    end.flatten
  end
  task :neko_bootstrap_sweep => :tsv do
    tsv = TSV.setup({}, "ID~Treatment,Scheme,Vetting,Target,Match,Miss,Total,Odds")
    id = 0
    dependencies.each do |dep| 
      match, miss_match = dep.load
      total = match.to_f + miss_match.to_f
      odds = match.to_f / total
      treatment = dep.recursive_inputs[:treatment]
      scheme = dep.recursive_inputs[:scheme]
      vetting = dep.recursive_inputs[:vetting]
      target = dep.recursive_inputs[:target]
      tsv[id]=[treatment, scheme, vetting, target, match, miss_match, total, odds]
      id += 1 
    end
    tsv
  end

end
