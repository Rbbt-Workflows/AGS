require 'rbbt/sources/organism'

module Entity::REST
  USE_ENSEMBL = false

end

module Object::Gene
  extend Entity
  include Entity::REST
  add_identifiers Organism.identifiers(AGS.organism), "Associated Gene Name", "Associated Gene Name"
end

module Object::Cluster
  extend Entity
  include Entity::REST

  property :genes do
    genes = AGS.gene_clusters.select("Cluster" => self).keys
    genes = AGS.gene_clusters.select("MetaCluster" => self).keys if genes.empty?
    genes = AGS.gene_clusters.select("SuperCluster" => self).keys if genes.empty?
    genes
  end

  property :name do
    if name = AGS.super_cluster_names[self]
      "(#{ self }) #{name}"
    else
      self
    end
  end

  self.format = ["MetaCluster", "SuperCluster"]
end

Workflow.require_workflow "Genomics"
Workflow.require_workflow "ExTRI"
Genomics.knowledge_base.register "ExTRI", ExTRI.job(:pairs).produce.path, :source => "Transcription Factor (Associated Gene Name)", :target => "Target Gene (Associated Gene Name)"
Workflow.require_workflow "Enrichment"
