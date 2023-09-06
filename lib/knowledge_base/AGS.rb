require 'rbbt/knowledge_base'

module AGS
  KB = KnowledgeBase.new Rbbt.var.knowledge_base.AGS.find, AGS.organism
  KB.format = {"Gene" => "Associated Gene Name"}

  KB.register :regulome, nil, :source => 'source (Associated Gene Name)', :target => 'target (Associated Gene Name)', :source_format => "Associated Gene Name", :target_format => "Associated Gene Name" do
    regulome = SaezLab.job(:regulome).run
    regulome.fields = ["source (Associated Gene Name)", "target (Associated Gene Name)", "weight"]
    regulome
  end
end
