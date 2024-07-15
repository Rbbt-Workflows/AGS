module AGS
  input :subset, :array, "Subset to a list of genes"
  task :fold_changes_NTNU => :tsv do |subset|
    tsv = Rbbt.data.FC["deseq2_log2foldchange_ashrshrink_outliersremoved_untreatedT0asreference_independentfiltering.txt"].tsv
    tsv.identifiers = Organism.identifiers(AGS.organism)

    log :translate, "Translate to gene name"
    tsv = tsv.change_key("Associated Gene Name")

    tsv.delete("")

    if subset && subset.any?
      log :subset, "Subset"
      tsv = tsv.select(subset)
    end

    tsv.transpose
  end

  input :subset, :array, "Subset to a list of genes"
  task :pvalues_NTNU => :tsv do |subset|
    tsv = TSV.open Rbbt.data.FC["deseq2_pvalue_fdradjusted_outliersremoved_untreatedT0asreference_independentfiltering.txt"], :header_hash => '', :type => :list
    genes = Rbbt.data.FC["deseq2_log2foldchange_ashrshrink_outliersremoved_untreatedT0asreference_independentfiltering.txt"].tsv.keys

    ids = {}
    i=0
    genes.each do |gene|
      i += 1
      ids[i] = gene
    end

    fixed = tsv.annotate({})

    TSV.traverse tsv do |id,values|
      gene = ids[id.to_i]
      values = values.collect{|v| v == "NA" ? nil : v.to_f}
      fixed[gene] = values
    end

    tsv = fixed

    tsv.key_field = "Ensembl Gene ID"

    tsv = tsv.select(subset) if subset && subset.any?

    tsv = tsv.change_key("Associated Gene Name", :identifiers => Organism.identifiers(AGS.organism))

    tsv.delete("")

    tsv.cast = :to_f

    tsv.transpose
  end
end
