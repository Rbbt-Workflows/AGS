require 'safe_ruby'
require 'AGS/tasks/gene_clusters/rules'
require 'AGS/tasks/gene_clusters'
module AGS

  # ToDo: Check fcs_one[0], should it be equal to fcs[0]?
  def self.apply_cluster_rules(fcs, pvalues, cluster_rules)
    fcs_one = [fcs[0]] 
    fcs_two = [0,fcs[0]] 
    fcs_ratio = [1,1] 

    # make pvalues consider the latest value
    
    (1..fcs.length-1).each do |i|
      fcs_one << fcs[i] - fcs[i-1]
    end

    (2..fcs.length-1).each do |i|
      fcs_two << fcs[i] - fcs[i-2]
    end

    (2..fcs.length-1).each do |i|
      fcs_ratio << fcs_one[i] / fcs_one[i-1]
    end

    #{{{ UP RULES

    threshold = 0.1

    clusters = []
    cluster_rules.each do |name,rules|
      rules.each do |rule|
        if res = eval(rule)
          break if res == 'break'
          clusters << name 
          break
        end
      end
    end

    clusters = ["unclassified"] if clusters.empty?
    clusters.sort_by{|c| c.scan(/\d+/)[0].to_i}
  end

  def self.cluster(fcs, pvalues, cluster_rules = nil)
    cluster_rules ||= self.cluster_rules
    apply_cluster_rules(fcs, pvalues, cluster_rules)
  end

  dep :fold_changes, :fc_source => "NTNU"
  dep :pvalues, :fc_source => "NTNU"
  input :custom_rule, :text, "Custom rule", nil
  input :protein_coding, :boolean, "Only keep protein coding genes", true
  task :gene_clusters => :tsv do |custom_rule,protein_coding|
    fc_tsv = step(:fold_changes).load.transpose "Associated Gene Name"
    pvalue_tsv = step(:pvalues_NTNU).load.transpose "Associated Gene Name"
    experiments = fc_tsv.fields.collect{|f| f.split(".").first.sub("FC_", "") }.uniq 

    cluster_rules = AGS.cluster_rules
    cluster_rules["custom"] = custom_rule.split("\n") if custom_rule

    cpus = config :cpus, :gene_clusters, :default => 1

    if protein_coding
      gene_info = Organism.identifiers(AGS.organism).tsv(:key_field => "Associated Gene Name", :fields => ["Ensembl Gene ID"], :type => :list)
      gene_info = gene_info.attach Organism.transcripts(AGS.organism), :fields => ["Ensembl Transcript ID"]
      gene_info = gene_info.attach Organism.transcript_biotype(AGS.organism), :fields => ["Ensembl Transcript Biotype"]

      protein_coding_genes = Set.new(gene_info.select("Ensembl Transcript Biotype" => 'protein_coding').keys)
    end

    clusters = TSV::Dumper.new :key_field => "Associated Gene Name", :fields => experiments, :type => :double
    clusters.init
    TSV.traverse fc_tsv, :into => clusters, :bar => self.progress_bar("Processing genes"), :cpus => cpus do |gene,values|
      if protein_coding && ! protein_coding_genes.include?(gene)
        next
      end
      clusters = experiments.collect do |experiment|
        experiment_fields = [1, 2, 4, 8, 24].collect{|t| "FC_#{experiment}.T#{t}" }
        fcs = values.values_at *experiment_fields
        pvalues = pvalue_tsv[gene].values_at *experiment_fields

        AGS.cluster(fcs, pvalues, cluster_rules)
      end
      [gene, clusters]
    end
  end

  dep :gene_clusters
  task :gene_clusters_unfolded => :tsv do
    tsv = step(:gene_clusters).load
    clusters = AGS.cluster_rules.keys + ["unclassified"]
    experiments = tsv.fields
    experiments.each do |experiment|
      clusters.each do |cluster|
        name = [experiment, cluster] * "-"
        tsv.add_field name do |k,v|
          v[experiment].include? cluster
        end
      end
    end
    tsv
  end

  dep :fold_changes, :fc_source => "NTNU"
  dep :gene_clusters_unfolded
  task :gene_clusters_fc => :tsv do
    fcs = step(:fold_changes).load.transpose "Associated Gene Name"
    clusters = step(:gene_clusters_unfolded).load
    Log.tsv fcs
    Log.tsv clusters
    clusters.attach fcs, :complete => true
  end

  dep :fold_changes, :fc_source => "NTNU"
  dep :gene_clusters_unfolded
  task :gene_clusters_fc_2 => :tsv do
    fcs = step(:fold_changes).load.transpose "Associated Gene Name"
    clusters = step(:gene_clusters_unfolded).load
    clusters_2 = TSV.setup(clusters.keys, :key_field => "Associated Gene Name", :fields => [], :type => :list)
    types = clusters.values.flatten.uniq.sort.sort_by{|s| s.scan(/\d+/)[0].to_i } 
    experiments = clusters.fields

    types.each do |type|
      experiments.each do |experiment|
        field = [experiment, type] * " - "

        clusters_2.add_field field do |g,v|
          clusters[g][experiment].include? type
        end
      end
    end


    clusters.attach clusters_2, :complete => true
    clusters.attach fcs, :complete => true
  end
end
