require 'rbbt/statistics/hypergeometric'

module AGS

  #dep :expressed_coding_genes, jobname: "Default"
  #input :list, :array, "Gene set"
  #input :background, :array, "Background genes", nil
  #input :database, :select, "Annotation to find enrichment", :go_bp, select_options: %w(go_bp)
  #task :gs_hyper => :tsv do |list,background,database|
  #  background = step(:expressed_coding_genes).load if background.nil?

  #  tsv = case database.to_sym
  #        when :go_bp
  #          Organism.gene_go_bp(AGS.organism).tsv type: :flat
  #        else
  #          raise ParameterException, "Unkown database parameter #{database}"
  #        end

  #  tsv = tsv.change_key "Associated Gene Name"

  #  tsv.enrichment list, nil, background: background
  #end

  #{{{ G:PROFILER

  dep :expressed_coding_genes, jobname: "Default"
  input :list, :array, "Gene set"
  input :background, :array, "Background genes", nil
  input :database, :select, "Annotation to find enrichment", :go_bp, select_options: %w(go_bp)
  task :gprofiler => :tsv do |list,background,database|
    require 'rbbt/util/python'
    background = step(:expressed_coding_genes).load if background.nil?

    list = list.load if Step === list
    list.shift if list[1] && list[1].start_with?(">")
    RbbtPython.run 'gprofiler' do
      gp = gprofiler.GProfiler.new(return_dataframe:true)
      if background && background.any?
        res = gp.profile(organism:'hsapiens', query: list, sources: %w(GO:BP), background: background)
      else
        res = gp.profile(organism:'hsapiens', query: list, sources: %w(GO:BP))
      end
      tsv = RbbtPython.df2tsv res
      tsv.type = :list
      tsv.key_field = "Position"
      tsv = tsv.reorder "native"
      tsv.fields = tsv.fields.collect{|f| f == "p_value" ? "p-value" : f }
      tsv.each do |k,v|
        v[13] = v[13..-1] * "|"
        v.replace v.slice(0, 14)
      end
      tsv
    end
  end

  dep :expressed_coding_genes
  input :queries, :yaml, "Gene set"
  input :background, :array, "Background genes", nil
  input :database, :select, "Annotation to find enrichment", :go_bp, select_options: %w(go_bp)
  task :gprofiler_multiple => :tsv do |queries,background,database|
    require 'rbbt/util/python'
    background = step(:expressed_coding_genes).load if background.nil?

    queries = queries.select{|k,v|
      k.include? "_both"
    }
    RbbtPython.run 'gprofiler' do
      gp = gprofiler.GProfiler.new(return_dataframe:true)
      if background && background.any?
        #res = gp.profile(organism:'hsapiens', query: queries, sources: %w(GO:BP), background: background)
        res = gp.profile(organism:'hsapiens', query: queries, sources: %w(GO:BP), background: background)
      else
        res = gp.profile(organism:'hsapiens', query: queries, sources: %w(GO:BP))
      end
      tsv = RbbtPython.df2tsv res
      tsv = tsv.to_double
      tsv.key_field = "Position"
      tsv.fields = tsv.fields.collect{|f| f == "p_value" ? "p-value" : f }
      tsv.each do |k,v|
        v[13] = v[13..-1].flatten
        v.replace v.slice(0, 14)
      end
      tsv
    end
  end

  dep :change_offsets_simplified
  dep :tf_predictions
  task :gprofiler_queries => :array do
    change_offsets = step(:change_offsets_simplified).load
    AGS::TIME_POINTS.each do |timepoint|
      AGS::TREATMENTS.each do |treatment|
        up_cluster = change_offsets.select(treatment => "increase #{timepoint}h").keys
        down_cluster = change_offsets.select(treatment => "decrease #{timepoint}h").keys

        tsv = AGS.job(:timepoint_matrix, treatment: treatment, time_point: timepoint, data_type: :fc).run.transpose
        fc_up_03 = tsv.select do |k,values|
          values.flatten.first.to_f > 0.3
        end.keys
        fc_down_03 = tsv.select do |k,values|
          values.flatten.first.to_f < -0.3
        end.keys

        fc_up_07 = tsv.select do |k,values|
          values.flatten.first.to_f > 0.7
        end.keys
        fc_down_07 = tsv.select do |k,values|
          values.flatten.first.to_f < -0.7
        end.keys


        tsv = AGS.job(:timepoint_matrix, treatment: treatment, time_point: timepoint, data_type: :fc0).run.transpose
        fc0_up_03 = tsv.select do |k,values|
          values.flatten.first.to_f > 0.3
        end.keys
        fc0_down_03 = tsv.select do |k,values|
          values.flatten.first.to_f < -0.3
        end.keys

        fc0_up_07 = tsv.select do |k,values|
          values.flatten.first.to_f > 0.7
        end.keys
        fc0_down_07 = tsv.select do |k,values|
          values.flatten.first.to_f < -0.7
        end.keys

        [
          [up_cluster, 'cluster', 'up'],
          [down_cluster, 'cluster', 'down'],
          [fc_up_03, 'fc_03', 'up'],
          [fc_down_03, 'fc_03', 'down'],
          [fc_up_07, 'fc_07', 'up'],
          [fc_down_07, 'fc_07', 'down'],
          [fc0_up_03, 'fc0_03', 'up'],
          [fc0_down_03, 'fc0_03', 'down'],
          [fc0_up_07, 'fc0_07', 'up'],
          [fc0_down_07, 'fc0_07', 'down'],
        ].each do |list,type,direction|

          name = [treatment, timepoint.to_s+'h', type, direction] * "-"
          file(name + '.txt').write <<-EOF
>#{name}
#{list*"\n"}
          EOF
        end
      end
    end

    predictions = step(:tf_predictions).load
    
    predictions.fields.each do |treatment_tp|
      treatment, timepoint = treatment_tp.split('-')
      timepoint = timepoint.sub('T','').to_s + 'h'
      up = predictions.select(treatment_tp){|v| v.to_f > 0 }.keys
      down = predictions.select(treatment_tp){|v| v.to_f < 0 }.keys

      up_name = [treatment, timepoint, 'TF', 'up'] * '-'
      down_name = [treatment, timepoint, 'TF', 'down'] * '-'
      file(up_name + '.txt').write <<-EOF
>#{up_name}
#{up*"\n"}
      EOF

      file(down_name + '.txt').write <<-EOF
>#{down_name}
#{down*"\n"}
      EOF
    end

    files
  end

  dep :gprofiler_queries, compute: :produce
  dep :gprofiler do |jobname,options,dependencies|
    queries = dependencies.flatten.first
    queries.files.collect do |file|
      options.merge(list: queries.file(file).list, jobname: file.sub('.txt', ''))
    end
  end
  task :gprofiler_suite => :tsv do
    tsv = TSV.setup({}, key_field: "Treatment:Time", fields: ["Up", "Down"], type: :double)
    dependencies[1..-1].collect do |dep|
      name = dep.clean_name
      Open.cp dep.path, file(name + '.tsv')
      next unless name.include?("cluster")
      treatment, hour, type, direction = name.split("-")
      tp = [treatment, hour + "h"] * ":"
      dep.load.each do |id, values|
        name = values[2]
        pvalue = values[3]
        next unless pvalue.to_f < 0.05
        tsv[tp] ||= [[], []]
        case direction
        when "up"
          tsv[tp][0] << name 
        when "down"
          tsv[tp][1] << name 
        else
          next
        end
      end
    end
    tsv
  end

  #}}} G:PROFILER

  #{{{ HELPERS

  # Paper-oriented functional enrichment utilities
  #
  # These tasks keep the broad g:Profiler GO-BP analysis available while adding
  # two more manuscript-facing layers: GO-Slim enrichment and MSigDB Hallmark
  # enrichment. Both use onset-defined dynamic gene sets as the query layer and
  # an expressed-coding gene background.

  helper :enrichment_scalar_value do |value|
    value = value.compact.first if Array === value
    value = nil if value.respond_to?(:empty?) && value.empty?
    value
  end

  helper :enrichment_treatment_order do
    %w(DMSO FiveZ INT_FiveZ_PI INT_PD_PI PD PI)
  end

  helper :enrichment_treatment_sort_index do |treatment|
    enrichment_treatment_order.index(treatment.to_s) || AGS::TREATMENTS.index(treatment.to_s) || 999
  end

  helper :enrichment_onset_events do |label|
    values = [label].flatten.compact.collect(&:to_s).reject{|v| v.empty? }
    values = values.collect{|v| v.split('|') }.flatten
    values.reject!{|v| v.nil? || v.empty? || v == 'unclassified' }
    values.collect do |entry|
      if entry =~ /^(increase|decrease)\s+(\d+)(?:-(\d+))?h$/
        direction = $1 == 'increase' ? 'up' : 'down'
        start_time = $2.to_i
        end_time = ($3 || $2).to_i
        {:label => entry, :direction => direction, :start_time => start_time, :end_time => end_time}
      end
    end.compact
  end

  helper :enrichment_parse_go_obo_terms do |text|
    terms = {}
    current = nil
    in_term = false
    finalize = lambda do
      if current && current[:id] && ! current[:obsolete]
        terms[current[:id]] = current
      end
    end
    text.each_line do |line|
      line = line.chomp
      case line
      when '[Term]'
        finalize.call
        current = {:parents => []}
        in_term = true
      when /^\[/
        finalize.call if in_term
        current = nil
        in_term = false
      else
        next unless in_term && current
        case line
        when /^id:\s+(GO:\d+)/
          current[:id] = $1
        when /^name:\s+(.+)/
          current[:name] = $1.strip
        when /^namespace:\s+(.+)/
          current[:namespace] = $1.strip
        when /^is_a:\s+(GO:\d+)/
          current[:parents] << $1
        when /^relationship:\s+part_of\s+(GO:\d+)/
          current[:parents] << $1
        when /^is_obsolete:\s+true/
          current[:obsolete] = true
        end
      end
    end
    finalize.call
    terms
  end

  helper :enrichment_go_ancestors_for_term do |term_id, terms, cache, visiting = nil|
    visiting ||= {}
    return cache[term_id] if cache.include?(term_id)
    return [] if visiting[term_id]
    term = terms[term_id]
    return cache[term_id] = [] if term.nil?
    visiting[term_id] = true
    ancestors = [term_id]
    term[:parents].each do |parent|
      ancestors.concat(enrichment_go_ancestors_for_term(parent, terms, cache, visiting))
    end
    visiting.delete(term_id)
    cache[term_id] = ancestors.uniq
  end

  helper :enrichment_hypergeom_upper_tail do |k, m, n, population|
    k = k.to_i
    m = m.to_i
    n = n.to_i
    population = population.to_i
    return 1.0 if k <= 0 || m <= 0 || n <= 0 || population <= 0
    max_k = [m, n].min
    return 1.0 if k > max_k
    log_choose = lambda do |a,b|
      return -Float::INFINITY if b < 0 || b > a
      Math.lgamma(a + 1)[0] - Math.lgamma(b + 1)[0] - Math.lgamma(a - b + 1)[0]
    end
    denom = log_choose.call(population, n)
    logs = (k..max_k).collect do |i|
      log_choose.call(m, i) + log_choose.call(population - m, n - i) - denom
    end
    max_log = logs.max
    return 0.0 if max_log == -Float::INFINITY
    pvalue = Math.exp(max_log) * logs.collect{|l| Math.exp(l - max_log) }.inject(0.0, &:+)
    [[pvalue, 1.0].min, 0.0].max
  end

  helper :enrichment_bh_adjust_values do |pvalues|
    indexed = pvalues.each_with_index.select{|pvalue,i| ! pvalue.nil? }.sort_by{|pvalue,i| pvalue.to_f }
    adjusted = Array.new(pvalues.length)
    m = indexed.length
    min_q = 1.0
    indexed.reverse.each_with_index do |(pvalue, i), rank_from_end|
      rank = m - rank_from_end
      qvalue = pvalue.to_f * m / rank
      min_q = [min_q, qvalue].min
      adjusted[i] = [min_q, 1.0].min
    end
    adjusted
  end

  helper :enrichment_overlap_score do |genes_a, genes_b, method|
    a = genes_a.to_set
    b = genes_b.to_set
    return 0.0 if a.empty? || b.empty?
    intersection = (a & b).length.to_f
    case method.to_s
    when 'jaccard'
      intersection / (a | b).length.to_f
    else
      intersection / [a.length, b.length].min.to_f
    end
  end

  #}}} HELPERS

  #{{{ ONSET GENE SETS

  dep :full_gene_info
  input :source_type, :select, 'Gene sets to summarize', 'cluster', :select_options => %w(cluster)
  task :functional_onset_gene_sets => :tsv do |source_type|
    info = step(:full_gene_info).load
    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(Treatment Time Direction SourceType Genes QuerySize), :type => :list, :namespace => AGS.organism)
    id = 0
    enrichment_treatment_order.each do |treatment|
      AGS::TIME_POINTS.each do |time_point|
        %w(up down both).each do |direction|
          genes = []
          info.through do |gene, values|
            values = NamedArray.setup(values, info.fields)
            events = enrichment_onset_events(values["#{treatment}: FC clusters"])
            selected = events.select{|event| event[:start_time] == time_point }
            selected = selected.select{|event| event[:direction] == direction } unless direction == 'both'
            genes << gene if selected.any?
          end
          id += 1
          genes = genes.uniq.sort
          tsv[id] = [treatment, time_point, direction, source_type, genes * '|', genes.length]
        end
      end
    end
    tsv
  end

  #}}} ONSET GENE SETS

  #{{{ GOSLIM
  dep :functional_onset_gene_sets
  task :goslim_bp_gene_sets => :tsv do
    step(:functional_onset_gene_sets).load
  end

  input :go_slim_url, :string, 'URL for the generic GO slim OBO subset', 'https://current.geneontology.org/ontology/subsets/goslim_generic.obo'
  task :go_slim_bp_terms => :tsv do |go_slim_url|
    require 'open-uri'
    slim_text = URI.open(go_slim_url, 'User-Agent' => 'Mozilla/5.0').read
    slim_terms = enrichment_parse_go_obo_terms(slim_text)
    tsv = TSV.setup({}, 'SlimGOID~SlimName,Namespace')
    slim_terms.keys.sort.each do |go_id|
      term = slim_terms[go_id]
      next unless term[:namespace] == 'biological_process'
      tsv[go_id] = [term[:name], term[:namespace]]
    end
    tsv
  end

  input :go_basic_url, :string, 'URL for go-basic.obo', 'https://current.geneontology.org/ontology/go-basic.obo'
  input :go_slim_url, :string, 'URL for the generic GO slim OBO subset', 'https://current.geneontology.org/ontology/subsets/goslim_generic.obo'
  task :go_bp_to_goslim_map => :tsv do |go_basic_url, go_slim_url|
    require 'open-uri'
    require 'set'
    go_text = URI.open(go_basic_url, 'User-Agent' => 'Mozilla/5.0').read
    slim_text = URI.open(go_slim_url, 'User-Agent' => 'Mozilla/5.0').read
    terms = enrichment_parse_go_obo_terms(go_text)
    slim_terms = enrichment_parse_go_obo_terms(slim_text)
    slim_bp_ids = slim_terms.keys.select{|go_id| slim_terms[go_id][:namespace] == 'biological_process' }
    ancestor_cache = {}
    tsv = TSV.setup({}, :key_field => 'GOID', :fields => %w(GOName Namespace SlimGOIDs SlimNames SlimGOIDsAll SlimNamesAll), :type => :list)
    terms.keys.sort.each do |go_id|
      term = terms[go_id]
      next unless term[:namespace] == 'biological_process'
      all_slim = (enrichment_go_ancestors_for_term(go_id, terms, ancestor_cache) & slim_bp_ids)
      next if all_slim.empty?
      reduced = all_slim.reject do |candidate|
        all_slim.any? do |other|
          other != candidate && enrichment_go_ancestors_for_term(other, terms, ancestor_cache).include?(candidate)
        end
      end
      reduced = all_slim if reduced.empty?
      reduced = reduced.sort
      all_slim = all_slim.sort
      tsv[go_id] = [term[:name], term[:namespace], reduced * '|', reduced.collect{|slim_id| slim_terms[slim_id] ? slim_terms[slim_id][:name] : terms[slim_id][:name] } * '|', all_slim * '|', all_slim.collect{|slim_id| slim_terms[slim_id] ? slim_terms[slim_id][:name] : terms[slim_id][:name] } * '|']
    end
    tsv
  end

  dep :go_bp_to_goslim_map
  task :gene_goslim_bp_annotations => :tsv do
    go_to_slim = step(:go_bp_to_goslim_map).load
    gene_go = Organism.gene_go_bp(AGS.organism).tsv(type: :flat)
    identifiers = Organism.identifiers(AGS.organism).tsv(:key_field => 'Ensembl Gene ID', :fields => ['Associated Gene Name'], :type => :flat)
    slim_name_by_id = {}
    go_to_slim.through do |go_id, values|
      values = NamedArray.setup(values, go_to_slim.fields)
      ids = enrichment_scalar_value(values['SlimGOIDs']).to_s.split('|').reject{|v| v.empty? }
      names = enrichment_scalar_value(values['SlimNames']).to_s.split('|')
      ids.each_with_index{|slim_id, i| slim_name_by_id[slim_id] ||= names[i] || slim_id }
    end
    annotations_by_gene = Hash.new{|h,k| h[k] = [] }
    gene_go.through do |ensembl_gene, go_ids|
      next unless identifiers.include?(ensembl_gene)
      gene_names = [identifiers[ensembl_gene]].flatten.compact.collect(&:to_s).reject{|gene| gene.empty? }
      next if gene_names.empty?
      slim_ids = []
      [go_ids].flatten.compact.each do |go_id|
        next unless go_to_slim.include?(go_id)
        slim_ids.concat enrichment_scalar_value(go_to_slim[go_id]['SlimGOIDs']).to_s.split('|').reject{|v| v.empty? }
      end
      slim_ids = slim_ids.uniq
      next if slim_ids.empty?
      gene_names.each{|gene| annotations_by_gene[gene].concat(slim_ids) }
    end
    tsv = TSV.setup({}, :key_field => 'Associated Gene Name', :fields => %w(SlimGOIDs SlimNames), :type => :list, :namespace => AGS.organism)
    annotations_by_gene.keys.sort.each do |gene|
      slim_ids = annotations_by_gene[gene].uniq.sort
      tsv[gene] = [slim_ids * '|', slim_ids.collect{|slim_id| slim_name_by_id[slim_id] || slim_id } * '|']
    end
    tsv
  end

  dep :goslim_bp_gene_sets
  dep :gene_goslim_bp_annotations
  dep :expressed_coding_genes
  input :min_query_size, :integer, 'Minimum genes in query set', 10
  input :min_intersection, :integer, 'Minimum genes in GO slim term intersection', 3
  task :goslim_bp_enrichment => :tsv do |min_query_size, min_intersection|
    require 'set'
    gene_sets = step(:goslim_bp_gene_sets).load
    annotations = step(:gene_goslim_bp_annotations).load
    background = step(:expressed_coding_genes).load.collect(&:to_s).to_set
    annotated_background = background & annotations.keys.collect(&:to_s).to_set
    genes_by_slim = Hash.new{|h,k| h[k] = [] }
    slim_name = {}
    annotations.through do |gene, values|
      values = NamedArray.setup(values, annotations.fields)
      next unless annotated_background.include?(gene.to_s)
      ids = enrichment_scalar_value(values['SlimGOIDs']).to_s.split('|').reject{|v| v.empty? }
      names = enrichment_scalar_value(values['SlimNames']).to_s.split('|')
      ids.each_with_index do |slim_id, i|
        genes_by_slim[slim_id] << gene.to_s
        slim_name[slim_id] ||= names[i] || slim_id
      end
    end
    genes_by_slim.each{|slim_id, genes| genes.uniq! }
    raw_rows = []
    gene_sets.through do |set_id, values|
      values = NamedArray.setup(values, gene_sets.fields)
      treatment = enrichment_scalar_value(values['Treatment'])
      time_point = enrichment_scalar_value(values['Time']).to_i
      direction = enrichment_scalar_value(values['Direction'])
      source_type = enrichment_scalar_value(values['SourceType'])
      query_genes = enrichment_scalar_value(values['Genes']).to_s.split('|').reject{|v| v.empty? }.to_set & annotated_background
      next if query_genes.length < min_query_size
      pvalues = []
      rows = []
      genes_by_slim.keys.sort.each do |slim_id|
        term_genes = genes_by_slim[slim_id].to_set
        intersection = (query_genes & term_genes).to_a.sort
        next if intersection.length < min_intersection
        pvalue = enrichment_hypergeom_upper_tail(intersection.length, term_genes.length, query_genes.length, annotated_background.length)
        pvalues << pvalue
        rows << [treatment, time_point, direction, source_type, slim_id, slim_name[slim_id], query_genes.length, term_genes.length, intersection.length, pvalue, nil, intersection.length.to_f / query_genes.length, intersection.length.to_f / term_genes.length, intersection * '|']
      end
      qvalues = enrichment_bh_adjust_values(pvalues)
      rows.each_with_index do |row, row_i|
        row[10] = qvalues[row_i]
        raw_rows << row
      end
    end
    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(Treatment Time Direction SourceType SlimGOID SlimName QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue Precision Recall Genes), :type => :list, :namespace => AGS.organism)
    raw_rows.sort_by{|row| [enrichment_treatment_sort_index(row[0]), row[1].to_i, row[2].to_s, row[10].to_f, row[5].to_s] }.each_with_index{|row, i| tsv[i + 1] = row }
    tsv
  end

  dep :goslim_bp_enrichment
  input :adjusted_pvalue_threshold, :float, 'Adjusted p-value threshold', 0.05
  input :top_n_per_context, :integer, 'Top GO slim terms to retain per treatment, time, and direction', 5
  task :goslim_bp_top_terms => :tsv do |adjusted_pvalue_threshold, top_n_per_context|
    enrichment = step(:goslim_bp_enrichment).load
    idx = Hash[enrichment.fields.each_with_index.to_a]
    grouped = Hash.new{|h,k| h[k] = [] }
    enrichment.through do |row_id, values|
      qvalue = enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f
      next if qvalue > adjusted_pvalue_threshold
      key = [enrichment_scalar_value(values[idx['Treatment']]), enrichment_scalar_value(values[idx['Time']]).to_i, enrichment_scalar_value(values[idx['Direction']])]
      grouped[key] << values
    end
    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(Treatment Time Direction SlimGOID SlimName QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue Precision Recall Genes), :type => :list, :namespace => AGS.organism)
    id = 0
    grouped.keys.sort_by{|treatment,time,direction| [enrichment_treatment_sort_index(treatment), time, direction.to_s] }.each do |key|
      grouped[key].sort_by{|values| [enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f, -enrichment_scalar_value(values[idx['IntersectionSize']]).to_i, enrichment_scalar_value(values[idx['SlimName']]).to_s] }.first(top_n_per_context).each do |values|
        id += 1
        tsv[id] = [enrichment_scalar_value(values[idx['Treatment']]), enrichment_scalar_value(values[idx['Time']]), enrichment_scalar_value(values[idx['Direction']]), enrichment_scalar_value(values[idx['SlimGOID']]), enrichment_scalar_value(values[idx['SlimName']]), enrichment_scalar_value(values[idx['QuerySize']]), enrichment_scalar_value(values[idx['TermBackgroundSize']]), enrichment_scalar_value(values[idx['IntersectionSize']]), enrichment_scalar_value(values[idx['PValue']]), enrichment_scalar_value(values[idx['AdjustedPValue']]), enrichment_scalar_value(values[idx['Precision']]), enrichment_scalar_value(values[idx['Recall']]), values[idx['Genes']]]
      end
    end
    tsv
  end

  #}}} GOSLIM

  #{{{ HALMARKS

  input :library_url, :string, 'Enrichr text export for MSigDB Hallmark gene sets', 'https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020'
  task :hallmark_gene_sets => :tsv do |library_url|
    require 'open-uri'
    text = URI.open(library_url, 'User-Agent' => 'Mozilla/5.0').read
    tsv = TSV.setup({}, :key_field => 'Hallmark', :fields => %w(Description Genes Size), :type => :list, :namespace => AGS.organism)
    text.each_line do |line|
      parts = line.chomp.split("\t")
      next if parts.length < 3
      hallmark = parts.shift
      description = parts.shift.to_s
      genes = parts.collect(&:to_s).reject{|gene| gene.empty? }.uniq.sort
      tsv[hallmark] = [description, genes * '|', genes.length]
    end
    tsv
  end

  dep :functional_onset_gene_sets
  dep :hallmark_gene_sets
  dep :expressed_coding_genes
  input :min_query_size, :integer, 'Minimum genes in query set', 10
  input :min_intersection, :integer, 'Minimum genes in Hallmark intersection', 3
  task :hallmark_enrichment => :tsv do |min_query_size, min_intersection|
    require 'set'
    gene_sets = step(:functional_onset_gene_sets).load
    hallmark = step(:hallmark_gene_sets).load
    background = step(:expressed_coding_genes).load.collect(&:to_s).to_set
    hallmark_genes = {}
    hallmark.through do |term, values|
      values = NamedArray.setup(values, hallmark.fields)
      genes = enrichment_scalar_value(values['Genes']).to_s.split('|').reject{|v| v.empty? }.to_set & background
      hallmark_genes[term] = genes.to_a.sort
    end
    effective_background = background & hallmark_genes.values.flatten.to_set
    raw_rows = []
    gene_sets.through do |set_id, values|
      values = NamedArray.setup(values, gene_sets.fields)
      treatment = enrichment_scalar_value(values['Treatment'])
      time_point = enrichment_scalar_value(values['Time']).to_i
      direction = enrichment_scalar_value(values['Direction'])
      source_type = enrichment_scalar_value(values['SourceType'])
      query_genes = enrichment_scalar_value(values['Genes']).to_s.split('|').reject{|v| v.empty? }.to_set & effective_background
      next if query_genes.length < min_query_size
      pvalues = []
      rows = []
      hallmark_genes.keys.sort.each do |term|
        term_genes = hallmark_genes[term].to_set & effective_background
        intersection = (query_genes & term_genes).to_a.sort
        next if intersection.length < min_intersection
        pvalue = enrichment_hypergeom_upper_tail(intersection.length, term_genes.length, query_genes.length, effective_background.length)
        pvalues << pvalue
        rows << [treatment, time_point, direction, source_type, term, term, query_genes.length, term_genes.length, intersection.length, pvalue, nil, intersection.length.to_f / query_genes.length, intersection.length.to_f / term_genes.length, intersection * '|']
      end
      qvalues = enrichment_bh_adjust_values(pvalues)
      rows.each_with_index do |row, row_i|
        row[10] = qvalues[row_i]
        raw_rows << row
      end
    end
    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(Treatment Time Direction SourceType TermID TermName QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue Precision Recall Genes), :type => :list, :namespace => AGS.organism)
    raw_rows.sort_by{|row| [enrichment_treatment_sort_index(row[0]), row[1].to_i, row[2].to_s, row[10].to_f, row[5].to_s] }.each_with_index{|row, i| tsv[i + 1] = row }
    tsv
  end

  dep :hallmark_enrichment
  input :adjusted_pvalue_threshold, :float, 'Adjusted p-value threshold', 0.05
  input :top_n_per_context, :integer, 'Top Hallmark terms to retain per treatment, time, and direction', 5
  task :hallmark_top_terms => :tsv do |adjusted_pvalue_threshold, top_n_per_context|
    enrichment = step(:hallmark_enrichment).load
    idx = Hash[enrichment.fields.each_with_index.to_a]
    grouped = Hash.new{|h,k| h[k] = [] }
    enrichment.through do |row_id, values|
      qvalue = enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f
      next if qvalue > adjusted_pvalue_threshold
      key = [enrichment_scalar_value(values[idx['Treatment']]), enrichment_scalar_value(values[idx['Time']]).to_i, enrichment_scalar_value(values[idx['Direction']])]
      grouped[key] << values
    end
    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(Treatment Time Direction TermID TermName QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue Precision Recall Genes), :type => :list, :namespace => AGS.organism)
    id = 0
    grouped.keys.sort_by{|treatment,time,direction| [enrichment_treatment_sort_index(treatment), time, direction.to_s] }.each do |key|
      grouped[key].sort_by{|values| [enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f, -enrichment_scalar_value(values[idx['IntersectionSize']]).to_i, enrichment_scalar_value(values[idx['TermName']]).to_s] }.first(top_n_per_context).each do |values|
        id += 1
        tsv[id] = [enrichment_scalar_value(values[idx['Treatment']]), enrichment_scalar_value(values[idx['Time']]), enrichment_scalar_value(values[idx['Direction']]), enrichment_scalar_value(values[idx['TermID']]), enrichment_scalar_value(values[idx['TermName']]), enrichment_scalar_value(values[idx['QuerySize']]), enrichment_scalar_value(values[idx['TermBackgroundSize']]), enrichment_scalar_value(values[idx['IntersectionSize']]), enrichment_scalar_value(values[idx['PValue']]), enrichment_scalar_value(values[idx['AdjustedPValue']]), enrichment_scalar_value(values[idx['Precision']]), enrichment_scalar_value(values[idx['Recall']]), values[idx['Genes']]]
      end
    end
    tsv
  end

  #}}} HALLMARKS

  #{{{ REDUCTION

  input :source, :select, 'Enrichment source to reduce', 'hallmark', :select_options => %w(hallmark goslim)
  input :adjusted_pvalue_threshold, :float, 'Adjusted p-value threshold', 0.05
  input :overlap_method, :select, 'Overlap score used for redundancy filtering', 'overlap_coefficient', :select_options => %w(overlap_coefficient jaccard)
  input :overlap_threshold, :float, 'Remove lower-ranked terms with overlap at or above this value', 0.5
  dep :hallmark_enrichment
  dep :goslim_bp_enrichment
  task :reduced_functional_enrichment => :tsv do |source, adjusted_pvalue_threshold, overlap_method, overlap_threshold|
    table = source.to_s == 'goslim' ? step(:goslim_bp_enrichment).load : step(:hallmark_enrichment).load
    idx = Hash[table.fields.each_with_index.to_a]
    term_field = idx['TermName'] ? 'TermName' : 'SlimName'
    id_field = idx['TermID'] ? 'TermID' : 'SlimGOID'
    grouped = Hash.new{|h,k| h[k] = [] }
    table.through do |row_id, values|
      qvalue = enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f
      next if qvalue > adjusted_pvalue_threshold
      key = [enrichment_scalar_value(values[idx['Treatment']]), enrichment_scalar_value(values[idx['Time']]).to_i, enrichment_scalar_value(values[idx['Direction']])]
      genes = enrichment_scalar_value(values[idx['Genes']]).to_s.split('|').reject{|v| v.empty? }
      grouped[key] << {
        :id => enrichment_scalar_value(values[idx[id_field]]),
        :name => enrichment_scalar_value(values[idx[term_field]]),
        :qvalue => qvalue,
        :pvalue => enrichment_scalar_value(values[idx['PValue']]).to_f,
        :intersection => enrichment_scalar_value(values[idx['IntersectionSize']]).to_i,
        :query_size => enrichment_scalar_value(values[idx['QuerySize']]).to_i,
        :term_size => enrichment_scalar_value(values[idx['TermBackgroundSize']]).to_i,
        :genes => genes
      }
    end
    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(Source Treatment Time Direction TermID TermName QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue MostSimilarKept MaxOverlap Genes), :type => :list, :namespace => AGS.organism)
    out_id = 0
    grouped.keys.sort_by{|treatment,time,direction| [enrichment_treatment_sort_index(treatment), time, direction.to_s] }.each do |key|
      selected = []
      grouped[key].sort_by{|row| [row[:qvalue], -row[:intersection], row[:name].to_s] }.each do |row|
        redundant_with = nil
        max_overlap = 0.0
        selected.each do |kept|
          overlap = enrichment_overlap_score(row[:genes], kept[:genes], overlap_method)
          if overlap > max_overlap
            max_overlap = overlap
            redundant_with = kept
          end
        end
        next if redundant_with && max_overlap >= overlap_threshold
        selected << row
        out_id += 1
        tsv[out_id] = [source, key[0], key[1], key[2], row[:id], row[:name], row[:query_size], row[:term_size], row[:intersection], row[:pvalue], row[:qvalue], redundant_with ? redundant_with[:name] : nil, max_overlap, row[:genes] * '|']
      end
    end
    tsv
  end

  #}}} REDUCTION

  #{{{ TF and TARGETS

  # TF-regulatory functional enrichment
  #
  # These tasks annotate sets of TF activity calls not only by the TF genes
  # themselves, but by the union of the TFs and their CollecTRI2 target genes.
  # This avoids the common TF-only enrichment problem where most terms collapse
  # to transcriptional regulation.

  helper :enrichment_parse_tf_context do |field|
    if field.to_s =~ /^(.+)-T(\d+)$/
      [$1, $2.to_i]
    else
      nil
    end
  end

  dep :tf_predictions
  dep :regulome
  dep :full_gene_info
  input :gene_mode, :select, 'Genes used to represent each TF activity set', 'tf_and_targets', :select_options => %w(tf_only tf_targets tf_and_targets)
  input :onset_targets_only, :boolean, 'Restrict CollecTRI2 targets to genes with an onset in the same treatment-time context', false
  task :tf_activity_regulatory_gene_sets => :tsv do |gene_mode, onset_targets_only|
    require 'set'
    predictions = step(:tf_predictions).load
    regulome = step(:regulome).load

    onset_targets_by_context = Hash.new{|h,k| h[k] = Set.new }
    if onset_targets_only
      info = step(:full_gene_info).load
      info.through do |gene, values|
        values = NamedArray.setup(values, info.fields)
        enrichment_treatment_order.each do |treatment|
          enrichment_onset_events(values["#{treatment}: FC clusters"]).each do |event|
            onset_targets_by_context[[treatment, event[:start_time]]] << gene.to_s
          end
        end
      end
    end

    targets_by_tf = Hash.new{|h,k| h[k] = [] }
    regulome.through do |edge_id, values|
      values = NamedArray.setup(values, regulome.fields)
      tf = enrichment_scalar_value(values['source']) || enrichment_scalar_value(values[0])
      target = enrichment_scalar_value(values['target']) || enrichment_scalar_value(values[1])
      next if tf.nil? || tf.to_s.empty? || target.nil? || target.to_s.empty?
      targets_by_tf[tf.to_s] << target.to_s
    end
    targets_by_tf.each{|tf, targets| targets.uniq! }

    rows = []
    predictions.fields.each do |field|
      parsed = enrichment_parse_tf_context(field)
      next if parsed.nil?
      treatment, time_point = parsed
      sign_tfs = Hash.new{|h,k| h[k] = [] }
      predictions.through do |tf, values|
        values = NamedArray.setup(values, predictions.fields)
        activity = enrichment_scalar_value(values[field]).to_f
        next if activity == 0
        sign = activity > 0 ? 'positive' : 'negative'
        sign_tfs[sign] << tf.to_s
        sign_tfs['both'] << tf.to_s
      end

      %w(positive negative both).each do |sign|
        tfs = sign_tfs[sign].uniq.sort
        targets = tfs.collect{|tf| targets_by_tf[tf] }.flatten.compact.uniq
        targets = targets.select{|target| onset_targets_by_context[[treatment, time_point]].include?(target) } if onset_targets_only
        targets = targets.sort
        genes = case gene_mode.to_s
                when 'tf_only'
                  tfs
                when 'tf_targets'
                  targets
                else
                  (tfs + targets).uniq.sort
                end
        rows << [treatment, time_point, sign, gene_mode, tfs, targets, genes]
      end
    end

    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(Treatment Time Sign GeneMode TargetScope TFs Targets Genes TFCount TargetCount GeneCount), :type => :list, :namespace => AGS.organism)
    rows.sort_by{|treatment,time,sign,mode,tfs,targets,genes| [enrichment_treatment_sort_index(treatment), time.to_i, sign] }.each_with_index do |row, i|
      treatment, time_point, sign, mode, tfs, targets, genes = row
      target_scope = onset_targets_only ? 'onset_targets_this_time' : 'all_collectri_targets'
      tsv[i + 1] = [treatment, time_point, sign, mode, target_scope, tfs * '|', targets * '|', genes * '|', tfs.length, targets.length, genes.length]
    end
    tsv
  end

  dep :tf_activity_regulatory_gene_sets
  dep :gene_goslim_bp_annotations
  dep :expressed_coding_genes
  input :min_query_size, :integer, 'Minimum genes in TF regulatory query set', 10
  input :min_intersection, :integer, 'Minimum genes in GO slim term intersection', 3
  task :tf_regulatory_goslim_bp_enrichment => :tsv do |min_query_size, min_intersection|
    require 'set'
    gene_sets = step(:tf_activity_regulatory_gene_sets).load
    annotations = step(:gene_goslim_bp_annotations).load
    background = step(:expressed_coding_genes).load.collect(&:to_s).to_set
    annotated_background = background & annotations.keys.collect(&:to_s).to_set

    genes_by_slim = Hash.new{|h,k| h[k] = [] }
    slim_name = {}
    annotations.through do |gene, values|
      values = NamedArray.setup(values, annotations.fields)
      next unless annotated_background.include?(gene.to_s)
      ids = enrichment_scalar_value(values['SlimGOIDs']).to_s.split('|').reject{|v| v.empty? }
      names = enrichment_scalar_value(values['SlimNames']).to_s.split('|')
      ids.each_with_index do |slim_id, i|
        genes_by_slim[slim_id] << gene.to_s
        slim_name[slim_id] ||= names[i] || slim_id
      end
    end
    genes_by_slim.each{|slim_id, genes| genes.uniq! }

    raw_rows = []
    gene_sets.through do |set_id, values|
      values = NamedArray.setup(values, gene_sets.fields)
      treatment = enrichment_scalar_value(values['Treatment'])
      time_point = enrichment_scalar_value(values['Time']).to_i
      sign = enrichment_scalar_value(values['Sign'])
      gene_mode = enrichment_scalar_value(values['GeneMode'])
      tf_count = enrichment_scalar_value(values['TFCount']).to_i
      target_count = enrichment_scalar_value(values['TargetCount']).to_i
      query_genes = enrichment_scalar_value(values['Genes']).to_s.split('|').reject{|v| v.empty? }.to_set & annotated_background
      next if query_genes.length < min_query_size

      pvalues = []
      rows = []
      genes_by_slim.keys.sort.each do |slim_id|
        term_genes = genes_by_slim[slim_id].to_set
        intersection = (query_genes & term_genes).to_a.sort
        next if intersection.length < min_intersection
        pvalue = enrichment_hypergeom_upper_tail(intersection.length, term_genes.length, query_genes.length, annotated_background.length)
        pvalues << pvalue
        rows << [treatment, time_point, sign, gene_mode, slim_id, slim_name[slim_id], tf_count, target_count, query_genes.length, term_genes.length, intersection.length, pvalue, nil, intersection.length.to_f / query_genes.length, intersection.length.to_f / term_genes.length, intersection * '|']
      end
      qvalues = enrichment_bh_adjust_values(pvalues)
      rows.each_with_index do |row, row_i|
        row[12] = qvalues[row_i]
        raw_rows << row
      end
    end

    fields = %w(Treatment Time Sign GeneMode SlimGOID SlimName TFCount TargetCount QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue Precision Recall Genes)
    tsv = TSV.setup({}, :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism)
    raw_rows.sort_by{|row| [enrichment_treatment_sort_index(row[0]), row[1].to_i, row[2].to_s, row[12].to_f, row[5].to_s] }.each_with_index{|row, i| tsv[i + 1] = row }
    tsv
  end

  dep :tf_regulatory_goslim_bp_enrichment
  input :adjusted_pvalue_threshold, :float, 'Adjusted p-value threshold', 0.05
  input :top_n_per_context, :integer, 'Top GO-Slim terms to retain per treatment, time, and sign', 5
  task :tf_regulatory_goslim_bp_top_terms => :tsv do |adjusted_pvalue_threshold, top_n_per_context|
    enrichment = step(:tf_regulatory_goslim_bp_enrichment).load
    idx = Hash[enrichment.fields.each_with_index.to_a]
    grouped = Hash.new{|h,k| h[k] = [] }
    enrichment.through do |row_id, values|
      qvalue = enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f
      next if qvalue > adjusted_pvalue_threshold
      key = [enrichment_scalar_value(values[idx['Treatment']]), enrichment_scalar_value(values[idx['Time']]).to_i, enrichment_scalar_value(values[idx['Sign']])]
      grouped[key] << values
    end
    fields = %w(Treatment Time Sign GeneMode SlimGOID SlimName TFCount TargetCount QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue Precision Recall Genes)
    tsv = TSV.setup({}, :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism)
    id = 0
    grouped.keys.sort_by{|treatment,time,sign| [enrichment_treatment_sort_index(treatment), time, sign.to_s] }.each do |key|
      grouped[key].sort_by{|values| [enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f, -enrichment_scalar_value(values[idx['IntersectionSize']]).to_i, enrichment_scalar_value(values[idx['SlimName']]).to_s] }.first(top_n_per_context).each do |values|
        id += 1
        tsv[id] = fields.collect{|field| enrichment_scalar_value(values[idx[field]]) }
      end
    end
    tsv
  end

  dep :tf_activity_regulatory_gene_sets
  dep :hallmark_gene_sets
  dep :expressed_coding_genes
  input :min_query_size, :integer, 'Minimum genes in TF regulatory query set', 10
  input :min_intersection, :integer, 'Minimum genes in Hallmark intersection', 3
  task :tf_regulatory_hallmark_enrichment => :tsv do |min_query_size, min_intersection|
    require 'set'
    gene_sets = step(:tf_activity_regulatory_gene_sets).load
    hallmark = step(:hallmark_gene_sets).load
    background = step(:expressed_coding_genes).load.collect(&:to_s).to_set
    hallmark_genes = {}
    hallmark.through do |term, values|
      values = NamedArray.setup(values, hallmark.fields)
      genes = enrichment_scalar_value(values['Genes']).to_s.split('|').reject{|v| v.empty? }.to_set & background
      hallmark_genes[term] = genes.to_a.sort
    end
    effective_background = background & hallmark_genes.values.flatten.to_set

    raw_rows = []
    gene_sets.through do |set_id, values|
      values = NamedArray.setup(values, gene_sets.fields)
      treatment = enrichment_scalar_value(values['Treatment'])
      time_point = enrichment_scalar_value(values['Time']).to_i
      sign = enrichment_scalar_value(values['Sign'])
      gene_mode = enrichment_scalar_value(values['GeneMode'])
      tf_count = enrichment_scalar_value(values['TFCount']).to_i
      target_count = enrichment_scalar_value(values['TargetCount']).to_i
      query_genes = enrichment_scalar_value(values['Genes']).to_s.split('|').reject{|v| v.empty? }.to_set & effective_background
      next if query_genes.length < min_query_size

      pvalues = []
      rows = []
      hallmark_genes.keys.sort.each do |term|
        term_genes = hallmark_genes[term].to_set & effective_background
        intersection = (query_genes & term_genes).to_a.sort
        next if intersection.length < min_intersection
        pvalue = enrichment_hypergeom_upper_tail(intersection.length, term_genes.length, query_genes.length, effective_background.length)
        pvalues << pvalue
        rows << [treatment, time_point, sign, gene_mode, term, term, tf_count, target_count, query_genes.length, term_genes.length, intersection.length, pvalue, nil, intersection.length.to_f / query_genes.length, intersection.length.to_f / term_genes.length, intersection * '|']
      end
      qvalues = enrichment_bh_adjust_values(pvalues)
      rows.each_with_index do |row, row_i|
        row[12] = qvalues[row_i]
        raw_rows << row
      end
    end

    fields = %w(Treatment Time Sign GeneMode TermID TermName TFCount TargetCount QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue Precision Recall Genes)
    tsv = TSV.setup({}, :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism)
    raw_rows.sort_by{|row| [enrichment_treatment_sort_index(row[0]), row[1].to_i, row[2].to_s, row[12].to_f, row[5].to_s] }.each_with_index{|row, i| tsv[i + 1] = row }
    tsv
  end

  dep :tf_regulatory_hallmark_enrichment
  input :adjusted_pvalue_threshold, :float, 'Adjusted p-value threshold', 0.05
  input :top_n_per_context, :integer, 'Top Hallmark terms to retain per treatment, time, and sign', 5
  task :tf_regulatory_hallmark_top_terms => :tsv do |adjusted_pvalue_threshold, top_n_per_context|
    enrichment = step(:tf_regulatory_hallmark_enrichment).load
    idx = Hash[enrichment.fields.each_with_index.to_a]
    grouped = Hash.new{|h,k| h[k] = [] }
    enrichment.through do |row_id, values|
      qvalue = enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f
      next if qvalue > adjusted_pvalue_threshold
      key = [enrichment_scalar_value(values[idx['Treatment']]), enrichment_scalar_value(values[idx['Time']]).to_i, enrichment_scalar_value(values[idx['Sign']])]
      grouped[key] << values
    end
    fields = %w(Treatment Time Sign GeneMode TermID TermName TFCount TargetCount QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue Precision Recall Genes)
    tsv = TSV.setup({}, :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism)
    id = 0
    grouped.keys.sort_by{|treatment,time,sign| [enrichment_treatment_sort_index(treatment), time, sign.to_s] }.each do |key|
      grouped[key].sort_by{|values| [enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f, -enrichment_scalar_value(values[idx['IntersectionSize']]).to_i, enrichment_scalar_value(values[idx['TermName']]).to_s] }.first(top_n_per_context).each do |values|
        id += 1
        tsv[id] = fields.collect{|field| enrichment_scalar_value(values[idx[field]]) }
      end
    end
    tsv
  end

  dep :tf_regulatory_hallmark_enrichment
  input :adjusted_pvalue_threshold, :float, 'Adjusted p-value threshold', 0.05
  input :overlap_method, :select, 'Overlap score used for redundancy filtering', 'overlap_coefficient', :select_options => %w(overlap_coefficient jaccard)
  input :overlap_threshold, :float, 'Remove lower-ranked terms with overlap at or above this value', 0.5
  task :reduced_tf_regulatory_hallmark_enrichment => :tsv do |adjusted_pvalue_threshold, overlap_method, overlap_threshold|
    table = step(:tf_regulatory_hallmark_enrichment).load
    idx = Hash[table.fields.each_with_index.to_a]
    grouped = Hash.new{|h,k| h[k] = [] }
    table.through do |row_id, values|
      qvalue = enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f
      next if qvalue > adjusted_pvalue_threshold
      key = [enrichment_scalar_value(values[idx['Treatment']]), enrichment_scalar_value(values[idx['Time']]).to_i, enrichment_scalar_value(values[idx['Sign']])]
      genes = enrichment_scalar_value(values[idx['Genes']]).to_s.split('|').reject{|v| v.empty? }
      grouped[key] << {
        :id => enrichment_scalar_value(values[idx['TermID']]),
        :name => enrichment_scalar_value(values[idx['TermName']]),
        :qvalue => qvalue,
        :pvalue => enrichment_scalar_value(values[idx['PValue']]).to_f,
        :intersection => enrichment_scalar_value(values[idx['IntersectionSize']]).to_i,
        :query_size => enrichment_scalar_value(values[idx['QuerySize']]).to_i,
        :term_size => enrichment_scalar_value(values[idx['TermBackgroundSize']]).to_i,
        :tf_count => enrichment_scalar_value(values[idx['TFCount']]).to_i,
        :target_count => enrichment_scalar_value(values[idx['TargetCount']]).to_i,
        :gene_mode => enrichment_scalar_value(values[idx['GeneMode']]),
        :genes => genes
      }
    end
    fields = %w(Source Treatment Time Sign GeneMode TermID TermName TFCount TargetCount QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue MostSimilarKept MaxOverlap Genes)
    tsv = TSV.setup({}, :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism)
    out_id = 0
    grouped.keys.sort_by{|treatment,time,sign| [enrichment_treatment_sort_index(treatment), time, sign.to_s] }.each do |key|
      selected = []
      grouped[key].sort_by{|row| [row[:qvalue], -row[:intersection], row[:name].to_s] }.each do |row|
        most_similar = nil
        max_overlap = 0.0
        selected.each do |kept|
          overlap = enrichment_overlap_score(row[:genes], kept[:genes], overlap_method)
          if overlap > max_overlap
            max_overlap = overlap
            most_similar = kept
          end
        end
        next if most_similar && max_overlap >= overlap_threshold
        selected << row
        out_id += 1
        tsv[out_id] = ['hallmark', key[0], key[1], key[2], row[:gene_mode], row[:id], row[:name], row[:tf_count], row[:target_count], row[:query_size], row[:term_size], row[:intersection], row[:pvalue], row[:qvalue], most_similar ? most_similar[:name] : nil, max_overlap, row[:genes] * '|']
      end
    end
    tsv
  end

  dep :tf_regulatory_goslim_bp_enrichment
  input :adjusted_pvalue_threshold, :float, 'Adjusted p-value threshold', 0.05
  input :overlap_method, :select, 'Overlap score used for redundancy filtering', 'overlap_coefficient', :select_options => %w(overlap_coefficient jaccard)
  input :overlap_threshold, :float, 'Remove lower-ranked terms with overlap at or above this value', 0.5
  task :reduced_tf_regulatory_goslim_bp_enrichment => :tsv do |adjusted_pvalue_threshold, overlap_method, overlap_threshold|
    table = step(:tf_regulatory_goslim_bp_enrichment).load
    idx = Hash[table.fields.each_with_index.to_a]
    grouped = Hash.new{|h,k| h[k] = [] }
    table.through do |row_id, values|
      qvalue = enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f
      next if qvalue > adjusted_pvalue_threshold
      key = [enrichment_scalar_value(values[idx['Treatment']]), enrichment_scalar_value(values[idx['Time']]).to_i, enrichment_scalar_value(values[idx['Sign']])]
      genes = enrichment_scalar_value(values[idx['Genes']]).to_s.split('|').reject{|v| v.empty? }
      grouped[key] << {
        :id => enrichment_scalar_value(values[idx['SlimGOID']]),
        :name => enrichment_scalar_value(values[idx['SlimName']]),
        :qvalue => qvalue,
        :pvalue => enrichment_scalar_value(values[idx['PValue']]).to_f,
        :intersection => enrichment_scalar_value(values[idx['IntersectionSize']]).to_i,
        :query_size => enrichment_scalar_value(values[idx['QuerySize']]).to_i,
        :term_size => enrichment_scalar_value(values[idx['TermBackgroundSize']]).to_i,
        :tf_count => enrichment_scalar_value(values[idx['TFCount']]).to_i,
        :target_count => enrichment_scalar_value(values[idx['TargetCount']]).to_i,
        :gene_mode => enrichment_scalar_value(values[idx['GeneMode']]),
        :genes => genes
      }
    end
    fields = %w(Source Treatment Time Sign GeneMode TermID TermName TFCount TargetCount QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue MostSimilarKept MaxOverlap Genes)
    tsv = TSV.setup({}, :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism)
    out_id = 0
    grouped.keys.sort_by{|treatment,time,sign| [enrichment_treatment_sort_index(treatment), time, sign.to_s] }.each do |key|
      selected = []
      grouped[key].sort_by{|row| [row[:qvalue], -row[:intersection], row[:name].to_s] }.each do |row|
        most_similar = nil
        max_overlap = 0.0
        selected.each do |kept|
          overlap = enrichment_overlap_score(row[:genes], kept[:genes], overlap_method)
          if overlap > max_overlap
            max_overlap = overlap
            most_similar = kept
          end
        end
        next if most_similar && max_overlap >= overlap_threshold
        selected << row
        out_id += 1
        tsv[out_id] = ['goslim', key[0], key[1], key[2], row[:gene_mode], row[:id], row[:name], row[:tf_count], row[:target_count], row[:query_size], row[:term_size], row[:intersection], row[:pvalue], row[:qvalue], most_similar ? most_similar[:name] : nil, max_overlap, row[:genes] * '|']
      end
    end
    tsv
  end


  dep :tf_activity_regulatory_gene_sets
  dep :regulome
  dep :gene_goslim_bp_annotations
  task :tf_regulatory_goslim_bp_annotation_counts => :tsv do
    require 'set'
    gene_sets = step(:tf_activity_regulatory_gene_sets).load
    regulome = step(:regulome).load
    annotations = step(:gene_goslim_bp_annotations).load

    slim_name_by_id = {}
    annotations.through do |gene, values|
      values = NamedArray.setup(values, annotations.fields)
      ids = enrichment_scalar_value(values['SlimGOIDs']).to_s.split('|').reject{|v| v.empty? }
      names = enrichment_scalar_value(values['SlimNames']).to_s.split('|')
      ids.each_with_index{|slim_id, i| slim_name_by_id[slim_id] ||= names[i] || slim_id }
    end

    targets_by_tf = Hash.new{|h,k| h[k] = [] }
    regulome.through do |edge_id, values|
      values = NamedArray.setup(values, regulome.fields)
      tf = enrichment_scalar_value(values['source']) || enrichment_scalar_value(values[0])
      target = enrichment_scalar_value(values['target']) || enrichment_scalar_value(values[1])
      next if tf.nil? || tf.to_s.empty? || target.nil? || target.to_s.empty?
      targets_by_tf[tf.to_s] << target.to_s
    end
    targets_by_tf.each{|tf, targets| targets.uniq! }

    fields = %w(Treatment Time Sign GeneMode SlimGOID SlimName TFCount TargetCount UniqueGenesWithTerm RepeatedGeneOccurrences QueryUniqueGenes QueryRepeatedOccurrences)
    tsv = TSV.setup({}, :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism)
    out_id = 0

    gene_sets.through do |set_id, values|
      values = NamedArray.setup(values, gene_sets.fields)
      treatment = enrichment_scalar_value(values['Treatment'])
      time_point = enrichment_scalar_value(values['Time']).to_i
      sign = enrichment_scalar_value(values['Sign'])
      gene_mode = enrichment_scalar_value(values['GeneMode'])
      tfs = enrichment_scalar_value(values['TFs']).to_s.split('|').reject{|v| v.empty? }
      tf_count = enrichment_scalar_value(values['TFCount']).to_i
      target_count = enrichment_scalar_value(values['TargetCount']).to_i

      occurrences = []
      occurrences.concat tfs if gene_mode != 'tf_targets'
      if gene_mode != 'tf_only'
        tfs.each do |tf|
          occurrences.concat targets_by_tf[tf]
        end
      end
      occurrences = occurrences.collect(&:to_s).reject{|gene| gene.empty? }
      unique_genes = occurrences.uniq

      repeated_counts = Hash.new(0)
      unique_counts = Hash.new{|h,k| h[k] = Set.new }
      occurrences.each do |gene|
        next unless annotations.include?(gene)
        row = NamedArray.setup(annotations[gene], annotations.fields)
        ids = enrichment_scalar_value(row['SlimGOIDs']).to_s.split('|').reject{|v| v.empty? }
        names = enrichment_scalar_value(row['SlimNames']).to_s.split('|')
        ids.each_with_index do |slim_id, i|
          slim_name_by_id[slim_id] ||= names[i] || slim_id
          repeated_counts[slim_id] += 1
          unique_counts[slim_id] << gene
        end
      end

      repeated_counts.keys.sort_by{|slim_id| [-repeated_counts[slim_id], slim_name_by_id[slim_id].to_s] }.each do |slim_id|
        out_id += 1
        tsv[out_id] = [treatment, time_point, sign, gene_mode, slim_id, slim_name_by_id[slim_id], tf_count, target_count, unique_counts[slim_id].length, repeated_counts[slim_id], unique_genes.length, occurrences.length]
      end
    end
    tsv
  end

  #}}} TF and TARGETS

  # Generic functional enrichment entry points
  #
  # These are the preferred paper-development tasks. They separate the query
  # layer from the annotation layer so the same enrichment code can be used for
  # onset genes, TF activity calls, TF targets, or TF plus target regulatory
  # sets, and for GO-Slim, MSigDB Hallmark, and Menyhart cancer hallmark
  # annotations.

  helper :functional_cancer_hallmark_url do |annotation|
    case annotation.to_s
    when 'cancerhallmarks_core'
      'https://cancerhallmarks.com/download_file/Menyhart_JPA_CancerHallmarks_core.txt'
    when 'cancerhallmarks_integrated'
      'https://cancerhallmarks.com/download_file/Menyhart_JPA_CancerHallmarks_integrated.txt'
    else
      nil
    end
  end

  helper :functional_parse_gene_set_text do |text|
    sets = {}
    text.each_line do |line|
      parts = line.chomp.split("\t").collect{|part| part.to_s.strip }
      next if parts.length < 2
      term = parts.shift
      next if term.empty?
      genes = parts.collect do |entry|
        entry.to_s.split(/[\/;, ]+/)
      end.flatten.collect(&:strip).reject{|gene| gene.empty? }.uniq.sort
      next if genes.empty?
      sets[term] = genes
    end
    sets
  end

  helper :functional_term_field_names do |annotation|
    case annotation.to_s
    when 'goslim'
      ['SlimGOID', 'SlimName']
    else
      ['TermID', 'TermName']
    end
  end

  input :annotation, :select, 'Annotation source', 'cancerhallmarks_core', :select_options => %w(goslim msigdb_hallmark cancerhallmarks_core cancerhallmarks_integrated)
  dep :gene_goslim_bp_annotations
  dep :hallmark_gene_sets
  dep :expressed_coding_genes
  task :functional_annotation_gene_sets => :tsv do |annotation|
    require 'open-uri'
    require 'set'
    background = step(:expressed_coding_genes).load.collect(&:to_s).to_set
    sets = {}
    names = {}

    case annotation.to_s
    when 'goslim'
      ann = step(:gene_goslim_bp_annotations).load
      ann.through do |gene, values|
        values = NamedArray.setup(values, ann.fields)
        next unless background.include?(gene.to_s)
        ids = enrichment_scalar_value(values['SlimGOIDs']).to_s.split('|').reject{|v| v.empty? }
        term_names = enrichment_scalar_value(values['SlimNames']).to_s.split('|')
        ids.each_with_index do |term_id, i|
          sets[term_id] ||= []
          sets[term_id] << gene.to_s
          names[term_id] ||= term_names[i] || term_id
        end
      end
    when 'msigdb_hallmark'
      hallmark = step(:hallmark_gene_sets).load
      hallmark.through do |term, values|
        values = NamedArray.setup(values, hallmark.fields)
        genes = enrichment_scalar_value(values['Genes']).to_s.split('|').reject{|v| v.empty? } & background.to_a
        next if genes.empty?
        sets[term.to_s] = genes.uniq.sort
        names[term.to_s] = term.to_s
      end
    when 'cancerhallmarks_core', 'cancerhallmarks_integrated'
      url = functional_cancer_hallmark_url(annotation)
      parsed = functional_parse_gene_set_text(URI.open(url, 'User-Agent' => 'Mozilla/5.0').read)
      parsed.each do |term, genes|
        filtered = genes & background.to_a
        next if filtered.empty?
        sets[term] = filtered.uniq.sort
        names[term] = term
      end
    else
      raise ParameterException, "Unknown annotation source #{annotation}"
    end

    tsv = TSV.setup({}, :key_field => 'TermID', :fields => %w(TermName Annotation Genes Size), :type => :list, :namespace => AGS.organism)
    sets.keys.sort.each do |term|
      genes = sets[term].uniq.sort
      tsv[term] = [names[term] || term, annotation, genes * '|', genes.length]
    end
    tsv
  end

  input :query_type, :select, 'Query gene-set source', 'onset_genes', :select_options => %w(onset_genes tfs targets tf_and_targets)
  dep :full_gene_info
  dep :tf_predictions
  dep :regulome
  task :functional_query_gene_sets => :tsv do |query_type|
    require 'set'
    tsv = TSV.setup({}, :key_field => 'ID', :fields => %w(QueryType Treatment Time Direction Genes TFs Targets QuerySize TFCount TargetCount), :type => :list, :namespace => AGS.organism)
    id = 0

    if query_type.to_s == 'onset_genes'
      info = step(:full_gene_info).load
      enrichment_treatment_order.each do |treatment|
        AGS::TIME_POINTS.each do |time_point|
          %w(up down both).each do |direction|
            genes = []
            info.through do |gene, values|
              values = NamedArray.setup(values, info.fields)
              events = enrichment_onset_events(values["#{treatment}: FC clusters"])
              selected = events.select{|event| event[:start_time] == time_point }
              selected = selected.select{|event| event[:direction] == direction } unless direction == 'both'
              genes << gene.to_s if selected.any?
            end
            genes = genes.uniq.sort
            id += 1
            tsv[id] = [query_type, treatment, time_point, direction, genes * '|', '', '', genes.length, 0, 0]
          end
        end
      end
      next tsv
    end

    predictions = step(:tf_predictions).load
    regulome = step(:regulome).load
    targets_by_tf = Hash.new{|h,k| h[k] = [] }
    regulome.through do |edge_id, values|
      values = NamedArray.setup(values, regulome.fields)
      tf = enrichment_scalar_value(values['source']) || enrichment_scalar_value(values[0])
      target = enrichment_scalar_value(values['target']) || enrichment_scalar_value(values[1])
      next if tf.nil? || tf.to_s.empty? || target.nil? || target.to_s.empty?
      targets_by_tf[tf.to_s] << target.to_s
    end
    targets_by_tf.each{|tf, targets| targets.uniq! }

    predictions.fields.each do |field|
      parsed = enrichment_parse_tf_context(field)
      next if parsed.nil?
      treatment, time_point = parsed
      tfs_by_sign = Hash.new{|h,k| h[k] = [] }
      predictions.through do |tf, values|
        values = NamedArray.setup(values, predictions.fields)
        activity = enrichment_scalar_value(values[field]).to_f
        next if activity == 0
        sign = activity > 0 ? 'positive' : 'negative'
        tfs_by_sign[sign] << tf.to_s
        tfs_by_sign['both'] << tf.to_s
      end

      %w(positive negative both).each do |sign|
        tfs = tfs_by_sign[sign].uniq.sort
        targets = tfs.collect{|tf| targets_by_tf[tf] }.flatten.compact.uniq.sort
        genes = case query_type.to_s
                when 'tfs'
                  tfs
                when 'targets'
                  targets
                else
                  (tfs + targets).uniq.sort
                end
        id += 1
        tsv[id] = [query_type, treatment, time_point, sign, genes * '|', tfs * '|', targets * '|', genes.length, tfs.length, targets.length]
      end
    end
    tsv
  end

  input :query_type, :select, 'Query gene-set source', 'onset_genes', :select_options => %w(onset_genes tfs targets tf_and_targets)
  input :annotation, :select, 'Annotation source', 'cancerhallmarks_core', :select_options => %w(goslim msigdb_hallmark cancerhallmarks_core cancerhallmarks_integrated)
  input :min_query_size, :integer, 'Minimum query size', 10
  input :min_intersection, :integer, 'Minimum term intersection', 3
  dep :functional_query_gene_sets do |jobname, options|
    options.merge(:query_type => options[:query_type])
  end
  dep :functional_annotation_gene_sets do |jobname, options|
    options.merge(:annotation => options[:annotation])
  end
  dep :expressed_coding_genes
  task :functional_enrichment => :tsv do |query_type, annotation, min_query_size, min_intersection|
    require 'set'
    queries = dependencies.find{|dep| dep.task_name == :functional_query_gene_sets }.load
    annotations = dependencies.find{|dep| dep.task_name == :functional_annotation_gene_sets }.load
    background = step(:expressed_coding_genes).load.collect(&:to_s).to_set

    term_genes = {}
    term_names = {}
    annotations.through do |term_id, values|
      values = NamedArray.setup(values, annotations.fields)
      genes = enrichment_scalar_value(values['Genes']).to_s.split('|').reject{|v| v.empty? }.to_set & background
      next if genes.empty?
      term_genes[term_id.to_s] = genes.to_a.sort
      term_names[term_id.to_s] = enrichment_scalar_value(values['TermName']).to_s
    end
    effective_background = background & term_genes.values.flatten.to_set

    raw_rows = []
    queries.through do |query_id, values|
      values = NamedArray.setup(values, queries.fields)
      qtype = enrichment_scalar_value(values['QueryType'])
      next unless qtype.to_s == query_type.to_s
      treatment = enrichment_scalar_value(values['Treatment'])
      time_point = enrichment_scalar_value(values['Time']).to_i
      direction = enrichment_scalar_value(values['Direction'])
      query_genes = enrichment_scalar_value(values['Genes']).to_s.split('|').reject{|v| v.empty? }.to_set & effective_background
      next if query_genes.length < min_query_size
      tf_count = enrichment_scalar_value(values['TFCount']).to_i
      target_count = enrichment_scalar_value(values['TargetCount']).to_i

      rows = []
      pvalues = []
      term_genes.keys.sort.each do |term_id|
        genes = term_genes[term_id].to_set
        intersection = (query_genes & genes).to_a.sort
        next if intersection.length < min_intersection
        pvalue = enrichment_hypergeom_upper_tail(intersection.length, genes.length, query_genes.length, effective_background.length)
        pvalues << pvalue
        rows << [qtype, annotation, treatment, time_point, direction, term_id, term_names[term_id], tf_count, target_count, query_genes.length, genes.length, intersection.length, pvalue, nil, intersection.length.to_f / query_genes.length, intersection.length.to_f / genes.length, intersection * '|']
      end
      qvalues = enrichment_bh_adjust_values(pvalues)
      rows.each_with_index do |row, row_i|
        row[13] = qvalues[row_i]
        raw_rows << row
      end
    end

    fields = %w(QueryType Annotation Treatment Time Direction TermID TermName TFCount TargetCount QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue Precision Recall Genes)
    tsv = TSV.setup({}, :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism)
    raw_rows.sort_by{|row| [row[0].to_s, row[1].to_s, enrichment_treatment_sort_index(row[2]), row[3].to_i, row[4].to_s, row[13].to_f, row[6].to_s] }.each_with_index do |row, i|
      tsv[i + 1] = row
    end
    tsv
  end

  input :query_type, :select, 'Query gene-set source', 'onset_genes', :select_options => %w(onset_genes tfs targets tf_and_targets)
  input :annotation, :select, 'Annotation source', 'cancerhallmarks_core', :select_options => %w(goslim msigdb_hallmark cancerhallmarks_core cancerhallmarks_integrated)
  input :adjusted_pvalue_threshold, :float, 'Adjusted p-value threshold', 0.05
  input :top_n_per_context, :integer, 'Top terms per context', 5
  dep :functional_enrichment do |jobname, options|
    options.merge(:query_type => options[:query_type], :annotation => options[:annotation])
  end
  task :functional_top_terms => :tsv do |query_type, annotation, adjusted_pvalue_threshold, top_n_per_context|
    enrichment = step(:functional_enrichment).load
    idx = Hash[enrichment.fields.each_with_index.to_a]
    grouped = Hash.new{|h,k| h[k] = [] }
    enrichment.through do |row_id, values|
      qvalue = enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f
      next if qvalue > adjusted_pvalue_threshold
      key = [enrichment_scalar_value(values[idx['QueryType']]), enrichment_scalar_value(values[idx['Annotation']]), enrichment_scalar_value(values[idx['Treatment']]), enrichment_scalar_value(values[idx['Time']]).to_i, enrichment_scalar_value(values[idx['Direction']])]
      grouped[key] << values
    end
    fields = %w(QueryType Annotation Treatment Time Direction TermID TermName TFCount TargetCount QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue Precision Recall Genes)
    tsv = TSV.setup({}, :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism)
    id = 0
    grouped.keys.sort_by{|qtype,ann,treatment,time,direction| [qtype.to_s, ann.to_s, enrichment_treatment_sort_index(treatment), time, direction.to_s] }.each do |key|
      grouped[key].sort_by{|values| [enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f, -enrichment_scalar_value(values[idx['IntersectionSize']]).to_i, enrichment_scalar_value(values[idx['TermName']]).to_s] }.first(top_n_per_context).each do |values|
        id += 1
        tsv[id] = fields.collect{|field| enrichment_scalar_value(values[idx[field]]) }
      end
    end
    tsv
  end

  input :query_type, :select, 'Query gene-set source', 'onset_genes', :select_options => %w(onset_genes tfs targets tf_and_targets)
  input :annotation, :select, 'Annotation source', 'cancerhallmarks_core', :select_options => %w(goslim msigdb_hallmark cancerhallmarks_core cancerhallmarks_integrated)
  input :adjusted_pvalue_threshold, :float, 'Adjusted p-value threshold', 0.05
  input :overlap_method, :select, 'Overlap score used for redundancy filtering', 'overlap_coefficient', :select_options => %w(overlap_coefficient jaccard)
  input :overlap_threshold, :float, 'Remove lower-ranked terms with overlap at or above this value', 0.5
  dep :functional_enrichment do |jobname, options|
    options.merge(:query_type => options[:query_type], :annotation => options[:annotation])
  end
  task :functional_enrichment_reduced => :tsv do |query_type, annotation, adjusted_pvalue_threshold, overlap_method, overlap_threshold|
    enrichment = step(:functional_enrichment).load
    idx = Hash[enrichment.fields.each_with_index.to_a]
    grouped = Hash.new{|h,k| h[k] = [] }
    enrichment.through do |row_id, values|
      qvalue = enrichment_scalar_value(values[idx['AdjustedPValue']]).to_f
      next if qvalue > adjusted_pvalue_threshold
      key = [enrichment_scalar_value(values[idx['QueryType']]), enrichment_scalar_value(values[idx['Annotation']]), enrichment_scalar_value(values[idx['Treatment']]), enrichment_scalar_value(values[idx['Time']]).to_i, enrichment_scalar_value(values[idx['Direction']])]
      genes = enrichment_scalar_value(values[idx['Genes']]).to_s.split('|').reject{|v| v.empty? }
      grouped[key] << {
        :term_id => enrichment_scalar_value(values[idx['TermID']]),
        :term_name => enrichment_scalar_value(values[idx['TermName']]),
        :tf_count => enrichment_scalar_value(values[idx['TFCount']]).to_i,
        :target_count => enrichment_scalar_value(values[idx['TargetCount']]).to_i,
        :query_size => enrichment_scalar_value(values[idx['QuerySize']]).to_i,
        :term_size => enrichment_scalar_value(values[idx['TermBackgroundSize']]).to_i,
        :intersection => enrichment_scalar_value(values[idx['IntersectionSize']]).to_i,
        :pvalue => enrichment_scalar_value(values[idx['PValue']]).to_f,
        :qvalue => qvalue,
        :precision => enrichment_scalar_value(values[idx['Precision']]),
        :recall => enrichment_scalar_value(values[idx['Recall']]),
        :genes => genes
      }
    end
    fields = %w(QueryType Annotation Treatment Time Direction TermID TermName TFCount TargetCount QuerySize TermBackgroundSize IntersectionSize PValue AdjustedPValue Precision Recall MostSimilarKept MaxOverlap Genes)
    tsv = TSV.setup({}, :key_field => 'ID', :fields => fields, :type => :list, :namespace => AGS.organism)
    id = 0
    grouped.keys.sort_by{|qtype,ann,treatment,time,direction| [qtype.to_s, ann.to_s, enrichment_treatment_sort_index(treatment), time, direction.to_s] }.each do |key|
      kept = []
      grouped[key].sort_by{|row| [row[:qvalue], -row[:intersection], row[:term_name].to_s] }.each do |row|
        most_similar = nil
        max_overlap = 0.0
        kept.each do |prev|
          overlap = enrichment_overlap_score(row[:genes], prev[:genes], overlap_method)
          if overlap > max_overlap
            max_overlap = overlap
            most_similar = prev
          end
        end
        next if most_similar && max_overlap >= overlap_threshold
        kept << row
        id += 1
        tsv[id] = key + [row[:term_id], row[:term_name], row[:tf_count], row[:target_count], row[:query_size], row[:term_size], row[:intersection], row[:pvalue], row[:qvalue], row[:precision], row[:recall], most_similar ? most_similar[:term_name] : nil, max_overlap, row[:genes] * '|']
      end
    end
    tsv
  end

  dep :functional_top_terms do |jobname,options|
    %w(onset_genes tfs targets tf_and_targets).collect do |query_type|
      %w(goslim msigdb_hallmark cancerhallmarks_core cancerhallmarks_integrated).collect do |annotation|
        {annotation: annotation, query_type: query_type}
      end
    end.flatten
  end
  dep :functional_enrichment_reduced do |jobname,options|
    %w(onset_genes tfs targets tf_and_targets).collect do |query_type|
      %w(goslim msigdb_hallmark cancerhallmarks_core cancerhallmarks_integrated).collect do |annotation|
        {annotation: annotation, query_type: query_type}
      end
    end.flatten
  end
  task :functional_enrichment_suite => :array do
    dependencies.collect do |dep|
      task_name = dep.task_name
      annotation = dep.recursive_inputs[:annotation]
      query_type = dep.recursive_inputs[:query_type]

      target = file("#{task_name}.#{annotation}.#{query_type}.tsv")
      Open.cp dep.path, target
      target
    end
  end


end
