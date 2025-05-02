module AGS

  dep :treatment_tfs_dynamic, ExTRI2_regulome: false, vetting: 'none'
  dep :treatment_tfs_dynamic, ExTRI2_regulome: false, vetting: 'synthesis'
  dep :treatment_tfs_dynamic, ExTRI2_regulome: false, vetting: 'degradation'
  dep :treatment_tfs_dynamic, ExTRI2_regulome: false, vetting: 'relaxed_degradation'
  dep :treatment_tfs_dynamic, ExTRI2_regulome: true, vetting: 'none'
  dep :treatment_tfs_dynamic, ExTRI2_regulome: true, vetting: 'synthesis'
  dep :treatment_tfs_dynamic, ExTRI2_regulome: true, vetting: 'degradation'
  dep :treatment_tfs_dynamic, ExTRI2_regulome: true, vetting: 'relaxed_degradation'
  dep :treatment_tfs_diff, ExTRI2_regulome: false, vetting: 'none'
  dep :treatment_tfs_diff, ExTRI2_regulome: false, vetting: 'synthesis'
  dep :treatment_tfs_diff, ExTRI2_regulome: false, vetting: 'degradation'
  dep :treatment_tfs_diff, ExTRI2_regulome: false, vetting: 'relaxed_degradation'
  dep :treatment_tfs_diff, ExTRI2_regulome: true, vetting: 'none'
  dep :treatment_tfs_diff, ExTRI2_regulome: true, vetting: 'synthesis'
  dep :treatment_tfs_diff, ExTRI2_regulome: true, vetting: 'degradation'
  dep :treatment_tfs_diff, ExTRI2_regulome: true, vetting: 'relaxed_degradation'
  dep :treatment_tfs_non_dynamic, ExTRI2_regulome: false, vetting: 'none'
  dep :treatment_tfs_non_dynamic, ExTRI2_regulome: false, vetting: 'synthesis'
  dep :treatment_tfs_non_dynamic, ExTRI2_regulome: false, vetting: 'degradation'
  dep :treatment_tfs_non_dynamic, ExTRI2_regulome: false, vetting: 'relaxed_degradation'
  dep :treatment_tfs_non_dynamic, ExTRI2_regulome: true, vetting: 'none'
  dep :treatment_tfs_non_dynamic, ExTRI2_regulome: true, vetting: 'synthesis'
  dep :treatment_tfs_non_dynamic, ExTRI2_regulome: true, vetting: 'degradation'
  dep :treatment_tfs_non_dynamic, ExTRI2_regulome: true, vetting: 'relaxed_degradation'
  task :treatment_overview => :tsv do
    dependencies.inject(nil) do |acc,dep|
      tsv = dep.load
      vetting = dep.recursive_inputs[:vetting]
      extri2 = dep.recursive_inputs[:ExTRI2_regulome]
      name = extri2 ? "ExTRI2" : "ExTRI1"
      name += " #{vetting}"
      name += " #{dep.task_name}"
      tsv.fields = tsv.fields.collect{|f| "#{f} #{name}" }
      acc = acc.nil? ? tsv : acc.attach(tsv, complete: true)
    end
  end

  dep :treatment_overview, treatment: :placeholder, method: :placeholder do |jobname,options|
    jobs = []
    AGS::TREATMENTS.each do |treatment|
      %w(ulm).each do |method|
        jobs << options.merge(treatment: treatment, method: method)
      end
    end
    jobs
  end
  task :treatment_overview_sweep => :tsv do
    dependencies.collect do |job|
      treatment, use_ExTRI2_regulome, method = job.recursive_inputs.values_at :treatment, :ExTRI2_regulome, :method
      name = ["treatment_overview", treatment, method, use_ExTRI2_regulome ? "ExTRI2" : "ExTRI1", 'tsv'] * "."
      Open.ln(job.path, file(name))
      file(name)
    end
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
