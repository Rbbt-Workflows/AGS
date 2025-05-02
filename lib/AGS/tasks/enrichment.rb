require 'rbbt/statistics/hypergeometric'

module AGS

  dep :expressed_coding_genes, jobname: "Default"
  input :list, :array, "Gene set"
  input :background, :array, "Background genes", nil
  input :database, :select, "Annotation to find enrichment", :go_bp, select_options: %w(go_bp)
  task :gs_hyper => :tsv do |list,background,database|
    background = step(:expressed_coding_genes).load if background.nil?

    tsv = case database.to_sym
          when :go_bp
            Organism.gene_go_bp(AGS.organism).tsv type: :flat
          else
            raise ParameterException, "Unkown database parameter #{database}"
          end

    tsv = tsv.change_key "Associated Gene Name"

    tsv.enrichment list, nil, background: background
  end

  dep :expressed_coding_genes, jobname: "Default"
  input :list, :array, "Gene set"
  input :background, :array, "Background genes", nil
  input :database, :select, "Annotation to find enrichment", :go_bp, select_options: %w(go_bp)
  task :gprofiler => :tsv do |list,background,database|
    require 'rbbt/util/python'
    background = step(:expressed_coding_genes).load if background.nil?

    list = list.load if Step === list
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

    RbbtPython.run 'gprofiler' do
      gp = gprofiler.GProfiler.new(return_dataframe:true)
      if background && background.any?
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
end
