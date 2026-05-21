module AGS
  input :ExTRI2_regulome, :boolean, "Use ExTRI2 regulome", true
  task_alias :regulome_pre, ExTRI2, :regulome, jobname: "Default" do |jobname,options|
      if options[:ExTRI2_regulome] 
        Workflow.require_workflow "ExTRI2"
        {workflow: ExTRI2, inputs: options}
      else
        {workflow: SaezLab, task: :regulome, inputs: options}
      end
    end

  dep :full_gene_info
  task :translations => :tsv do
    old = Organism.identifiers(AGS.organism).tsv key_field: 'Associated Gene Name', field: 'Ensembl Gene ID', type: :single, persist: true
    new = Organism.identifiers("Hsa").tsv key_field: 'Ensembl Gene ID', field: 'Associated Gene Name', type: :single, persist: true

    translations = {}
    genes = step(:full_gene_info).load.keys
    genes.each do |orig|
      ens = old[orig]
      name = new[ens]
      next if ens.nil?
      next if orig == name
      next if name.nil?
      translations[name] = orig
    end
    TSV.setup(translations, key_field: 'Modern name', fields: ['Old name'], type: :single)
  end

  dep :regulome_pre
  dep :translations
  task :regulome => :tsv do
    translations = step(:translations).load
    tsv = step(:regulome).load
    tsv.process 'source' do |gene|
      if translations[gene]
        translations[gene]
      else
        gene
      end
    end
    tsv.process 'target' do |gene|
      if translations[gene]
        translations[gene]
      else
        gene
      end
    end

    clean = tsv.annotate({})
    pairs = Set.new
    tsv.each do |k,v|
      pair = v.values_at(0, 1) * ':'
      next if pairs.include? pair
      pairs << pair
      clean[k] = v
    end
    clean
  end
end
