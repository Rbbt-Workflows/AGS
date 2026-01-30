module AGS

  dep :neko_bootstrap_summary, treatment: :placeholder, scheme: :placeholder, vetting: :placeholder do |jobname,options|
    [0.07, 0.05, 0.1].collect do |low_min_1h|
      [0.6, 10].collect do |mid_multiplier_1h|
        %w(relaxed T1 T2).collect do |target|
          %w(PD PI FiveZ PD-MAP2K1 PD-MAP2K2).collect do |treatment|
            options.merge(treatment: treatment, scheme: 'dynamic', vetting: 'none', target: target, low_min_1h: low_min_1h, mid_multiplier_1h: mid_multiplier_1h)
          end
        end
      end
    end.flatten
  end
  task :neko_bootstrap_rule_variants_t1 => :tsv do
    tsv = TSV.setup({}, "ID~Treatment,low_min_1h,mid_multiplier_1h,Target,Match,Miss,Total,Odds")
    id = 0
    dependencies.each do |dep| 
      match, miss_match = dep.load
      total = match.to_f + miss_match.to_f
      odds = match.to_f / total
      treatment = dep.recursive_inputs[:treatment]
      low_min_1h = dep.recursive_inputs[:low_min_1h]
      mid_multiplier_1h = dep.recursive_inputs[:mid_multiplier_1h]
      target = dep.recursive_inputs[:target]
      tsv[id]=[treatment, low_min_1h, mid_multiplier_1h, target, match, miss_match, total, odds]
      id += 1 
    end
    tsv
  end

  dep :treatment_tf_consistency, treatment: :placeholder, scheme: :placeholder, vetting: :placeholder do |jobname,options|
    #[0.25, 0.3,0.35,0.4,0.45].collect do |low_min_24h|
      #[0.25, 0.3,0.35,0.4,0.45].collect do |mid_min_24h|
    [0.25, 0.3,0.35,0.4,0.45].collect do |high_min_24h|
      %w(relaxed T1 T2).collect do |target|
        AGS::TREATMENTS.collect do |treatment|
          options.merge(treatment: treatment, scheme: 'dynamic', vetting: 'none', target: target, high_min_24h: high_min_24h)
        end
      end
    end.flatten.compact
  end
  task :consistency_summary_rule_variants_t24 => :tsv do
    tsv = TSV.setup({}, "ID~Treatment,low_min_24h,mid_min_24h,high_min_24h,Total,Consistent,Ratio")
    id = 0
    dependencies.each do |dep| 
      match, miss_match = dep.load
      values = Misc.counts(dep.column('Consistent at 24h').values.compact)
      total = Misc.sum(values.values)
      consistent = total['1']
      ratio = consistent.to_f / total
      treatment = dep.recursive_inputs[:treatment]
      low_min_24h = dep.recursive_inputs[:low_min_24h]
      mid_min_24h = dep.recursive_inputs[:mid_min_24h]
      high_min_24h = dep.recursive_inputs[:high_min_24h]
      tsv[id]=[treatment, low_min_24h, mid_min_24h, high_min_24h, total, consistent, ratio]
      id += 1 
    end
    tsv
  end
end
