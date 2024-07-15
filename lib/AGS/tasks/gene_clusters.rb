require 'safe_ruby'
require 'AGS/tasks/gene_clusters/rules'
require 'AGS/tasks/gene_clusters'
module AGS

  # ToDo: Check fcs_one[0], should it be equal to fcs[0]?
  #def self.apply_cluster_rules(fcs, pvalues, cluster_rules)
  helper :apply_cluster_rules do |fcs, pvalues, cluster_rules,binding|
    fcs_one = [fcs[0]] 
    fcs_two = [0,fcs[0]] 
    fcs_three = [0, 0, fcs[0]] 
    fcs_ratio = [1,1] 

    (1..fcs.length-1).each do |i|
      one_step_change = fcs[i] - fcs[i-1]
      #one_step_change = 0.01 if one_step_change == 0
      fcs_one << one_step_change
    end

    (2..fcs.length-1).each do |i|
      fcs_two << fcs[i] - fcs[i-2]
    end

    (2..fcs.length-1).each do |i|
      fcs_three << fcs[i] - fcs[i-3]
    end

    (2..fcs.length-1).each do |i|
      fcs_ratio << fcs_one[i] / fcs_one[i-1]
    end

    # make pvalues consider the latest value
    p_one = [pvalues[0]]
    p_two = [1, pvalues[1]]
    p_three = [1, 1, pvalues[2]]
    
    (1..fcs.length-1).each do |i|
      p_one << [pvalues[i], pvalues[i-1]].max
    end

    (2..fcs.length-1).each do |i|
      p_two << [pvalues[i], pvalues[i-2]].max
    end

    (2..fcs.length-1).each do |i|
      p_three << [pvalues[i], pvalues[i-3]].max
    end


    threshold = 0.05

    s_one = p_one.collect{|v| v.to_f <= threshold }
    s_two = p_two.collect{|v| v.to_f <= threshold }
    s_three = p_three.collect{|v| v.to_f <= threshold }

    #{{{ UP RULES

    self.binding.local_variables.each do |v|
      binding.local_variable_set(v, self.binding.local_variable_get(v))
    end

    clusters = []
    cluster_rules.each do |name,rules|
      rules.each do |rule|
        binding.local_variable_set(:clusters, clusters)
        if res = eval(rule,binding)
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
  extension :png
  task :vulcano_plot => :binary do
    fc_tsv = step(:fold_changes).load.transpose "Associated Gene Name"
    pvalue_tsv = step(:pvalues).load.transpose "Associated Gene Name"

    tsv = TSV.setup({}, "Key~FC,pvalue#:type=:list")

    timepoints = fc_tsv.fields.reject{|f| f.include? "DMSO" }
    timepoints = fc_tsv.fields.select{|f| f.include? "PD_PI" }
    TSV.traverse fc_tsv.keys, bar: self.progress_bar("Processing genes") do |gene|
      timepoints.each do |tp|
        fc = fc_tsv[gene][tp]
        p = pvalue_tsv[gene][tp]
        key = Misc.digest([gene, tp] * ":")
        tsv[key] = [fc.to_f, p.to_f]
      end
    end

    R::PNG.ggplot self.tmp_path, tsv, <<-EOF
ggplot(data = data, aes(x = FC, y = -log10(pvalue))) +
   geom_vline(xintercept = c(-0.6, 0.6), col = 'gray', linetype = 'dashed') + 
   geom_hline(yintercept = -log10(0.05), col = 'red', linetype = 'dashed') + 
   geom_hline(yintercept = -log10(0.1), col = 'yellow', linetype = 'dashed') + 
   geom_point(alpha=0.01) 
    EOF

    nil
  end

  #dep :fold_changes, :fc_source => "NTNU"
  #dep :pvalues, :fc_source => "NTNU"
  #input :custom_rule, :text, "Custom rule", nil
  #input :protein_coding, :boolean, "Only keep protein coding genes", true
  #input :low_min_1h, :float, "", 0.1
  #input :mid_min_1h, :float, "", 0.15
  #input :high_min_1h, :float, "", 0.25
  #input :next_min_1h, :float, "", 0.03
  #input :low_multiplier_1h, :float, "", 8
  #input :mid_multiplier_1h, :float, "", 0.6
  #input :low_min_2h, :float, "", 0.1
  #input :mid_min_2h, :float, "", 0.2
  #input :high_min_2h, :float, "", 0.3
  #input :next_min_2h, :float, "", 0.1
  #input :low_multiplier_2h, :float, "", 8
  #input :mid_multiplier_2h, :float, "", 0.7
  #input :low_min_4h, :float, "", 0.1
  #input :mid_min_4h, :float, "", 0.2
  #input :high_min_4h, :float, "", 0.3
  #input :next_min_4h, :float, "", 0.1
  #input :low_multiplier_4h, :float, "", 10_000
  #input :mid_multiplier_4h, :float, "", 0.9
  #input :low_min_8h, :float, "", 0.1
  #input :mid_min_8h, :float, "", 0.2
  #input :high_min_8h, :float, "", 0.3
  #input :next_min_8h, :float, "", 0.15
  #input :low_multiplier_8h, :float, "", 8
  #input :mid_multiplier_8h, :float, "", 0.9
  #input :low_min_24h, :float, "", 0.3
  #input :mid_min_24h, :float, "", 0.3
  #input :high_min_24h, :float, "", 0.3
  #input :next_min_24h, :float, "", 0.15
  #input :low_multiplier_24h, :float, "", 8
  #input :mid_multiplier_24h, :float, "", 0.9
  #task :gene_clusters => :tsv do |custom_rule,protein_coding, 
  #  low_min_1h, mid_min_1h, high_min_1h, next_min_1h, low_multiplier_1h, mid_multiplier_1h,
  #  low_min_2h, mid_min_2h, high_min_2h, next_min_2h, low_multiplier_2h, mid_multiplier_2h, 
  #  low_min_4h, mid_min_4h, high_min_4h, next_min_4h, low_multiplier_4h, mid_multiplier_4h,
  #  low_min_8h, mid_min_8h, high_min_8h, next_min_8h, low_multiplier_8h, mid_multiplier_8h,
  #  low_min_24h, mid_min_24h, high_min_24h, next_min_24h, low_multiplier_24h, mid_multiplier_24h|

  #  if protein_coding
  #    gene_info = Organism.identifiers(AGS.organism).tsv(:key_field => "Associated Gene Name", :fields => ["Ensembl Gene ID"], :type => :double)
  #    gene_info = gene_info.attach Organism.transcripts(AGS.organism), :fields => ["Ensembl Transcript ID"]
  #    gene_info = gene_info.attach Organism.transcript_biotype(AGS.organism), :fields => ["Ensembl Transcript Biotype"], one2one: true

  #    protein_coding_genes = Set.new(gene_info.select("Ensembl Transcript Biotype" => 'protein_coding').keys)
  #  end

  #  fc_tsv = step(:fold_changes).load.transpose "Associated Gene Name"
  #  pvalue_tsv = step(:pvalues).load.transpose "Associated Gene Name"
  #  experiments = fc_tsv.fields.collect{|f| f.split(".").first.sub("FC_", "") }.uniq 

  #  cluster_rules = AGS.cluster_rules
  #  cluster_rules["custom"] = custom_rule.split("\n") if custom_rule

  #  cpus = config :cpus, :gene_clusters, :default => 1

  #  clusters = TSV::Dumper.new :key_field => "Associated Gene Name", :fields => experiments, :type => :double, :namespace => AGS.organism
  #  clusters.init
  #  TSV.traverse fc_tsv, :into => clusters, :bar => self.progress_bar("Processing genes"), :cpus => cpus do |gene,values|
  #    if protein_coding && ! protein_coding_genes.include?(gene)
  #      next
  #    end
  #    clusters = experiments.collect do |experiment|
  #      experiment_fields = [1, 2, 4, 8, 24].collect{|t| fc_tsv.fields.index("FC_#{experiment}.T#{t}") }
  #      fcs = values.values_at *experiment_fields
  #      pvalues = pvalue_tsv[gene].values_at *experiment_fields

  #      apply_cluster_rules(fcs, pvalues, AGS.cluster_rules, self.binding)
  #    end
  #    [gene, clusters]
  #  end
  #end

  #dep :gene_clusters
  #task :gene_clusters_unfolded => :tsv do
  #  tsv = step(:gene_clusters).load
  #  clusters = AGS.cluster_rules.keys + ["unclassified"]
  #  experiments = tsv.fields
  #  experiments.each do |experiment|
  #    clusters.each do |cluster|
  #      name = [experiment, cluster] * "-"
  #      tsv.add_field name do |k,v|
  #        v[experiment].include? cluster
  #      end
  #    end
  #  end
  #  tsv
  #end

  #dep :fold_changes, :fc_source => "NTNU"
  #dep :gene_clusters_unfolded
  #task :gene_clusters_fc => :tsv do
  #  fcs = step(:fold_changes).load.transpose "Associated Gene Name"
  #  clusters = step(:gene_clusters_unfolded).load
  #  clusters.attach fcs, :complete => true
  #end

  #dep :fold_changes, :fc_source => "NTNU"
  #dep :gene_clusters
  #task :gene_clusters_fc_simple => :tsv do
  #  fcs = step(:fold_changes).load.transpose "Associated Gene Name"
  #  clusters = step(:gene_clusters).load
  #  clusters.attach fcs
  #end


  #dep :fold_changes, :fc_source => "NTNU"
  #dep :gene_clusters_unfolded
  #task :gene_clusters_fc_2 => :tsv do
  #  fcs = step(:fold_changes).load.transpose "Associated Gene Name"
  #  clusters = step(:gene_clusters_unfolded).load
  #  clusters_2 = TSV.setup(clusters.keys, :key_field => "Associated Gene Name", :fields => [], :type => :list)
  #  types = clusters.values.flatten.uniq.sort.sort_by{|s| s.scan(/\d+/)[0].to_i } 
  #  experiments = clusters.fields

  #  types.each do |type|
  #    experiments.each do |experiment|
  #      field = [experiment, type] * " - "

  #      clusters_2.add_field field do |g,v|
  #        clusters[g][experiment].include? type
  #      end
  #    end
  #  end


  #  clusters.attach clusters_2, :complete => true
  #  clusters.attach fcs, :complete => true
  #end

  #task :cluster_rules => :text do
  #  Path.setup('lib/AGS/tasks/gene_clusters/rules.rb').find(:lib).read
  #end
end
