module AGS
  #task_alias :regulome, ExTRI2, :regulome, jobname: "Default", 
  #  filter_valid_tfs: true,
  #  only_authoritative_tfs: false,
  #  remove_auto_regulation: false,
  #  no_MoR: false do |jobname,options|
  #    if options[:ExTRI2_regulome] 
  #      Workflow.require_workflow "ExTRI2"
  #      {workflow: ExTRI2, inputs: options}
  #    else
  #      {workflow: SaezLab, task: :regulome, inputs: options}
  #    end
  #  end

  input :ExTRI2_regulome, :boolean, "Use ExTRI2 regulome", true
  task_alias :regulome, ExTRI2, :regulome, jobname: "Default" do |jobname,options|
      if options[:ExTRI2_regulome] 
        Workflow.require_workflow "ExTRI2"
        {workflow: ExTRI2, inputs: options}
      else
        {workflow: SaezLab, task: :regulome, inputs: options}
      end
    end
end
