module AGS

  def self.low_rule(fc_current, fc_next, current_min, next_min, next_multiplier)
    fc_current > current_min && fc_next > next_min && fc_current * next_multiplier > fc_next
  end

  def self.mid_rule(fc_current, fc_next, current_min, next_multiplier)
    fc_current > current_min && fc_next > - fc_current * next_multiplier * (fc_current/current_min)**2
  end

  def self.high_rule(fc_current, current_min)
    fc_current > current_min
  end

  def self.offset_rule(fc_current, fc_next, low_min, mid_min, high_min, next_min, low_multiplier, mid_multiplier)
    low_rule(fc_current, fc_next, low_min, next_min, low_multiplier) ||
      mid_rule(fc_current, fc_next, mid_min, mid_multiplier) ||
      high_rule(fc_current, high_min)
  end

  def self.collapse_offsets(offsets)
    last = nil
    start = nil
    segments = []
    current = []
    while i = offsets.shift
      if current.any? && current.last + 1 == i
        current << i
      elsif current.any?
        segments << [current.first, current.last]
        current = []
      else
        current = []
      end
      current << i
    end
    segments << [current.first, current.last] if current.any?

    segments.collect{|s,e| ([s,e].uniq.collect{|i| AGS::TIME_POINTS[i] } * "-") + "h" }
  end

  def self.clean_offsets(offsets)
    clean_offsets = []
    last = nil
    while offset = offsets.shift
      type = offset.split(" ").first

      clean_offsets << offset unless last == type

      last = type
    end
    clean_offsets
  end

  helper :offsets do |fcs, pvalues, inputs|
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

    increase_offsets = []
    decrease_offsets = []
    AGS::TIME_POINTS.each_with_index do |time_point,i|
      fc_current = fcs_one[i]
      fc_next = fcs_one[i+1] || 0
      low_min, mid_min, high_min, next_min, low_multiplier, mid_multiplier = 
        inputs.values_at(*%w(low_min mid_min high_min next_min low_multiplier mid_multiplier).collect{|v| v + "_#{time_point}h" })
      increase_offsets << i if AGS.offset_rule(fc_current, fc_next, low_min, mid_min, high_min, next_min, low_multiplier, mid_multiplier)
      decrease_offsets << i if AGS.offset_rule(-fc_current, -fc_next, low_min, mid_min, high_min, next_min, low_multiplier, mid_multiplier)
    end

    clusters = []
    clusters += AGS.collapse_offsets(increase_offsets.dup).collect{|t| "increase #{t}" }
    clusters += AGS.collapse_offsets(decrease_offsets.dup).collect{|t| "decrease #{t}" }

    clusters = clusters.sort_by{|c| c.scan(/\d+/)[0].to_i}

    clusters = AGS.clean_offsets(clusters)

    clusters = ["unclassified"] if clusters.empty?

    clusters
  end

  dep :fold_changes, :fc_source => "NTNU"
  dep :pvalues, :fc_source => "NTNU"
  input :protein_coding, :boolean, "Only keep protein coding genes", true
  input :low_min_1h, :float, "", 0.1
  input :mid_min_1h, :float, "", 0.15
  input :high_min_1h, :float, "", 0.25
  input :next_min_1h, :float, "", 0.02
  input :low_multiplier_1h, :float, "", 8
  input :mid_multiplier_1h, :float, "", 0.6
  input :low_min_2h, :float, "", 0.1
  input :mid_min_2h, :float, "", 0.15
  input :high_min_2h, :float, "", 0.3
  input :next_min_2h, :float, "", 0.05
  input :low_multiplier_2h, :float, "", 8
  input :mid_multiplier_2h, :float, "", 0.5
  input :low_min_4h, :float, "", 0.1
  input :mid_min_4h, :float, "", 0.15
  input :high_min_4h, :float, "", 0.3
  input :next_min_4h, :float, "", 0.05
  input :low_multiplier_4h, :float, "", 10
  input :mid_multiplier_4h, :float, "", 0.50
  input :low_min_8h, :float, "", 0.1
  input :mid_min_8h, :float, "", 0.15
  input :high_min_8h, :float, "", 0.3
  input :next_min_8h, :float, "", 0.10
  input :low_multiplier_8h, :float, "", 12
  input :mid_multiplier_8h, :float, "", 0.5
  input :low_min_24h, :float, "", 0.3
  input :mid_min_24h, :float, "", 0.3
  input :high_min_24h, :float, "", 0.3
  input :next_min_24h, :float, "", 0.15
  input :low_multiplier_24h, :float, "", 8
  input :mid_multiplier_24h, :float, "", 0.9
  task :change_offsets => :tsv do |protein_coding, 
    low_min_1h, mid_min_1h, high_min_1h, next_min_1h, low_multiplier_1h, mid_multiplier_1h,
    low_min_2h, mid_min_2h, high_min_2h, next_min_2h, low_multiplier_2h, mid_multiplier_2h, 
    low_min_4h, mid_min_4h, high_min_4h, next_min_4h, low_multiplier_4h, mid_multiplier_4h,
    low_min_8h, mid_min_8h, high_min_8h, next_min_8h, low_multiplier_8h, mid_multiplier_8h,
    low_min_24h, mid_min_24h, high_min_24h, next_min_24h, low_multiplier_24h, mid_multiplier_24h|

    if protein_coding
      log :gene_info, "Identifiers"
      gene_info = Organism.identifiers(AGS.organism).tsv(:key_field => "Associated Gene Name", :fields => ["Ensembl Gene ID"], :type => :double)
      log :gene_info, "Transcripts"
      gene_info = gene_info.attach Organism.transcripts(AGS.organism), :fields => ["Ensembl Transcript ID"]
      log :gene_info, "Biotype"
      gene_info = gene_info.attach Organism.transcript_biotype(AGS.organism), :fields => ["Ensembl Transcript Biotype"], one2one: true

      protein_coding_genes = Set.new(gene_info.select("Ensembl Transcript Biotype" => 'protein_coding').keys)
    end

    log :transpose, "Foldchanges"
    fc_tsv = step(:fold_changes).load.transpose "Associated Gene Name"
    log :transpose, "P-values"
    pvalue_tsv = step(:pvalues).load.transpose "Associated Gene Name"
    experiments = fc_tsv.fields.collect{|f| f.split(".").first.sub("FC_", "") }.uniq 

    cpus = config :cpus, :gene_clusters, :default => 1

    clusters = TSV::Dumper.new :key_field => "Associated Gene Name", :fields => experiments, :type => :double, :namespace => AGS.organism
    clusters.init
    TSV.traverse fc_tsv, :into => clusters, :bar => self.progress_bar("Processing genes"), :cpus => cpus do |gene,values|
      if protein_coding && ! protein_coding_genes.include?(gene)
        next
      end
      clusters = experiments.collect do |experiment|
        experiment_fields = AGS::TIME_POINTS.collect{|t| fc_tsv.fields.index("FC_#{experiment}.T#{t}") }
        fcs = values.values_at *experiment_fields
        pvalues = pvalue_tsv[gene].values_at *experiment_fields

        offsets(fcs, pvalues, self.inputs)
      end
      [gene, clusters]
    end
  end

end
