module AGS
  task :INSPEcT => :tsv do
    files = Rbbt.data.INSPEcT.glob("*")
    bar = self.progress_bar("Loading INSPEcT files", max: files.length)
    treatment_tsvs = files.collect do |dirname|
      Path.setup(dirname)
      treatment = File.basename(dirname)
      data_file = dirname.glob("*_data.csv").first
      gene_info_file = dirname.glob("*_gene_info.csv").first

      data = TSV.open(data_file, fix: proc{|l| l.gsub('"','') }, :header_hash => '')
      data.key_field = "Ensembl Gene ID"
      gene_info = TSV.open(gene_info_file, fix: proc{|l| l.gsub('"','') }, :header_hash => '')
      gene_info.key_field = "Ensembl Gene ID"

      data.fields = data.fields.collect{|f| [treatment, f.gsub('"','')] * ": " }
      gene_info.fields = gene_info.fields.collect{|f| [treatment, f.gsub('"','')] * ": " }

      gene_info.attach data, :complete => true, insitu: true

      bar.tick
      gene_info
    end
    bar.remove

    bar = self.progress_bar("Merging INSPEcT files", max: files.length)
    res = treatment_tsvs.inject(nil) do |acc,tsv|
      if acc.nil?
        acc = tsv
      else
        acc.attach tsv, :complete => true, insitu: true
      end
      bar.tick
      acc
    end
    bar.remove
    res
  end

  dep :INSPEcT
  task :synthesis_changes => :tsv do
    tsv = step(:INSPEcT).load

    time_points = %w(0 1 2 4 8 24)

    treatments = tsv.fields.collect{|f| f.split(": ").first }.uniq
    treatments.each do |treatment|
      field = [treatment, 'INSPEcT synthesis clusters'] * ": "
      synthesis_fields = time_points.collect{|t| [treatment, "synthesis_#{t}"] * ": " }
      synthesis_fields_index = synthesis_fields.collect{|f| tsv.fields.index f }
      tsv.add_field field do |ens,values|
        synthesis_values = values.values_at(*synthesis_fields_index).collect{|v| v.first }
        elbows = (synthesis_values.length - 1).times.collect do |i| 
          current = synthesis_values[i+1]
          next if current.nil?
          prev = synthesis_values[i]

          current = current.to_f
          prev = prev.to_f

          if current > prev * 1.05
            :up
          elsif current < prev * 0.95
            :down
          end
        end

        clusters = []
        last = nil
        elbows.each do |e|
          if e != last
            clusters << e
          else
            clusters << nil
          end
          last = e
        end

        time_points[1..-1].zip(elbows).collect do |t,e|
          next if e.nil?
          [e, t] * " "
        end.compact
      end
    end
    tsv
  end

  dep :INSPEcT
  task :degradation_changes => :tsv do
    tsv = step(:INSPEcT).load

    time_points = %w(0 1 2 4 8 24)

    treatments = tsv.fields.collect{|f| f.split(": ").first }.uniq
    treatments.each do |treatment|
      field = [treatment, 'INSPEcT degradation clusters'] * ": "
      synthesis_fields = time_points.collect{|t| [treatment, "degradation_#{t}"] * ": " }
      synthesis_fields_index = synthesis_fields.collect{|f| tsv.fields.index f }
      tsv.add_field field do |ens,values|
        synthesis_values = values.values_at(*synthesis_fields_index).collect{|v| v.first }
        elbows = (synthesis_values.length - 1).times.collect do |i| 
          current = synthesis_values[i+1]
          next if current.nil?
          prev = synthesis_values[i]

          current = current.to_f
          prev = prev.to_f

          if current > prev * 1.05
            :up
          elsif current < prev * 0.95
            :down
          end
        end

        clusters = []
        last = nil
        elbows.each do |e|
          if e != last
            clusters << e
          else
            clusters << nil
          end
          last = e
        end

        time_points[1..-1].zip(elbows).collect do |t,e|
          next if e.nil?
          [e, t] * " "
        end.compact
      end
    end
    tsv
  end

  dep :synthesis_changes
  dep :degradation_changes
  dep :fold_changes_NTNU
  dep :change_offsets
  task :full_gene_info => :tsv do
    s = step(:synthesis_changes).load
    d = step(:degradation_changes).load
    f = step(:fold_changes_NTNU).load
    f = f.transpose("Associated Gene Name")

    s = s.change_key "Associated Gene Name", identifiers: Organism.identifiers(AGS.organism)
    d = d.change_key "Associated Gene Name", identifiers: Organism.identifiers(AGS.organism)

    s.attach d, :identifiers => Organism.identifiers(AGS.organism), complete: true
    s.attach f, :identifiers => Organism.identifiers(AGS.organism), complete: true

    index = Organism.identifiers(AGS.organism).index :target => "Ensembl Gene ID", 
      :fields => ["Associated Gene Name"], 
      :order => true, 
      :persist => true
    s.add_field "Ensembl Gene ID" do |n|
      index[n]
    end

    s = s.reorder :key, ["Ensembl Gene ID"] + s.fields[0..-2] 
    clusters = step(:change_offsets).load
    clusters.fields = clusters.fields.collect{|f| [f, " FC clusters"] * ":" }
    s.attach clusters
  end

  desc "Synthesis genes are regulated by syntesis in more treatments than processing and degradation"
  dep :full_gene_info
  input :synthesis_criteria, :select, "Criteria to define synthesis", :mayority, :select_options => %w(mayority one two three)
  task :synthesis_genes => :tsv do |criteria|
    criteria = criteria.to_sym
    tsv = step(:full_gene_info).load

    synthesis_pos = tsv.fields.select{|f| f =~ /synthesis$/ }.collect{|f| tsv.fields.index f }
    processing_pos = tsv.fields.select{|f| f =~ /processing$/ }.collect{|f| tsv.fields.index f }
    degradation_pos = tsv.fields.select{|f| f =~ /degradation$/ }.collect{|f| tsv.fields.index f }
    threshold = 0.05
    tsv.add_field "Synthesis gene" do |k,values|
      synthesis_values = values.values_at(*synthesis_pos).collect{|f| f.first.to_f  < threshold }
      processing_values = values.values_at(*processing_pos).collect{|f| f.first.to_f  < threshold }
      degradation_values = values.values_at(*degradation_pos).collect{|f| f.first.to_f  < threshold }

      synthesis_count = synthesis_values[1..-1].select{|v| v }.length
      processing_count = processing_values[1..-1].select{|v| v }.length
      degradation_count = degradation_values[1..-1].select{|v| v }.length

      degradation = case criteria
                  when :mayority
                    synthesis_count >= 1 && (synthesis_count == [synthesis_count, processing_count, degradation_count].max)
                  when :one
                    synthesis_count >= 1
                  when :two
                    synthesis_count >= 2
                  when :three
                    synthesis_count >= 3
                  else
                    raise "Criteria not understood #{criteria}"
                  end
      [degradation]
    end
  end

  desc "Degradation genes are regulated by degradation in more treatments than processing and synthesis"
  dep :full_gene_info
  input :degradation_criteria, :select, "Criteria to define synthesis", :mayority, :select_options => %w(mayority one two three)
  task :degradation_genes => :tsv do |criteria|
    criteria = criteria.to_sym
    tsv = step(:full_gene_info).load

    synthesis_pos = tsv.fields.select{|f| f =~ /synthesis$/ }.collect{|f| tsv.fields.index f }
    processing_pos = tsv.fields.select{|f| f =~ /processing$/ }.collect{|f| tsv.fields.index f }
    degradation_pos = tsv.fields.select{|f| f =~ /degradation$/ }.collect{|f| tsv.fields.index f }
    threshold = 0.05
    tsv.add_field "Degradation gene" do |k,values|
      synthesis_values = values.values_at(*synthesis_pos).collect{|f| f.first.to_f  < threshold }
      processing_values = values.values_at(*processing_pos).collect{|f| f.first.to_f  < threshold }
      degradation_values = values.values_at(*degradation_pos).collect{|f| f.first.to_f  < threshold }

      synthesis_count = synthesis_values[1..-1].select{|v| v }.length
      processing_count = processing_values[1..-1].select{|v| v }.length
      degradation_count = degradation_values[1..-1].select{|v| v }.length

      degradation = case criteria
                  when :mayority
                    degradation_count >= 1 && (degradation_count == [synthesis_count, processing_count, degradation_count].max)
                  when :one
                    degradation_count >= 1
                  when :two
                    degradation_count >= 2
                  when :three
                    degradation_count >= 3
                  else
                    raise "Criteria not understood #{criteria}"
                  end
      [degradation]
    end
  end

  desc "Vetted genes are systhesis genes where the up or down synthesis trends coincide with heuristic clustering in at least one timepoint in at least on treatment"
  dep :synthesis_genes
  task :vetted_genes => :tsv do
    tsv = step(:synthesis_genes).load
    treatments = %w(FiveZ INT_FiveZ_PI INT_PD_PI PD PI)
    time_points = [0,1,2,4,8,24]

    tsv.add_field "Degradation timepoints" do |k,values|
      treatments.collect do |treatment|
        clusters = values[treatment + ": FC clusters"]
        inspect = values[treatment + ": INSPEcT degradation clusters"]

        inspect_directions = {}
        inspect.each do |ic|
          direction, time = ic.split(" ")
          inspect_directions[direction] ||= []
          inspect_directions[direction] << time_points.index(time.to_i)
        end

        inspect_directions.each{|direction, values| values.sort! }

        clean_inspect_directions = {}
        inspect_directions.each do |direction,values|
          clean_inspect_directions[direction] = values.reject{|v| values.include?(v-1) }
        end

        matches = []
        (clean_inspect_directions["down"] || []).each do |time_index|
          matches << time_index if clusters.select{|value| value.include?("increase") }.collect{|value| time_points.index(value.split(" ").last.to_i) }.include?(time_index)
          matches << time_index if clusters.select{|value| value.include?("increase") }.collect{|value| time_points.index(value.split(" ").last.to_i) }.include?(time_index - 1)
          matches << time_index if clusters.select{|value| value.include?("increase") }.collect{|value| time_points.index(value.split(" ").last.to_i) }.include?(time_index + 1)
        end

        (clean_inspect_directions["up"] || []).each do |time_index|
          matches << time_index if clusters.select{|value| value.include?("decrease") }.collect{|value| time_points.index(value.split(" ").last.to_i) }.include?(time_index)
          matches << time_index if clusters.select{|value| value.include?("decrease") }.collect{|value| time_points.index(value.split(" ").last.to_i) }.include?(time_index - 1)
          matches << time_index if clusters.select{|value| value.include?("decrease") }.collect{|value| time_points.index(value.split(" ").last.to_i) }.include?(time_index + 1)
        end

        matches.collect{|m| [treatment, time_points[m]] * ":" }
      end.flatten
    end

    tsv.add_field "Treatment synthesis profile match" do |k,values|
      treatments.select do |treatment|
        clusters = values[treatment + ": FC clusters"]
        inspect = values[treatment + ": INSPEcT synthesis clusters"]

        inspect_directions = {}
        inspect.each do |ic|
          direction, time = ic.split(" ")
          inspect_directions[direction] ||= []
          inspect_directions[direction] << time_points.index(time.to_i)
        end

        inspect_directions.each{|direction, values| values.sort! }

        clean_inspect_directions = {}
        inspect_directions.each do |direction,values|
          clean_inspect_directions[direction] = values.reject{|v| values.include?(v-1)}
        end

        match = false
        (clean_inspect_directions["up"] || []).each do |time_index|
          match = true if clusters.select{|value| value.include?("increase") }.collect{|value| time_points.index(value.split(" ").last.to_i) }.include?(time_index)
          match = true if clusters.select{|value| value.include?("increase") }.collect{|value| time_points.index(value.split(" ").last.to_i) }.include?(time_index - 1)
          match = true if clusters.select{|value| value.include?("increase") }.collect{|value| time_points.index(value.split(" ").last.to_i) }.include?(time_index + 1)
        end

        (clean_inspect_directions["down"] || []).each do |time_index|
          match = true if clusters.select{|value| value.include?("decrease") }.collect{|value| time_points.index(value.split(" ").last.to_i) }.include?(time_index)
          match = true if clusters.select{|value| value.include?("decrease") }.collect{|value| time_points.index(value.split(" ").last.to_i) }.include?(time_index - 1)
          match = true if clusters.select{|value| value.include?("decrease") }.collect{|value| time_points.index(value.split(" ").last.to_i) }.include?(time_index + 1)
        end


        match
      end
    end

    tsv.add_field "Vetted synthesis gene" do |k,values|
      values["Treatment synthesis profile match"].any? && values["Synthesis gene"].first.to_s == "true"
    end

  end

  dep :vetted_genes
  dep ExTRI, :CollecTRI
  task :vetted_genes_in_CollecTRI => :tsv do
    collectri = step(:CollecTRI).load
    tgs = collectri.column("Target Gene (Associated Gene Name)").values.flatten.compact.uniq

    tsv = step(:vetted_genes).load
    tsv.add_field "In CollecTRI" do |ens,values|
      gene = values.flatten.first
      tgs.include?(gene)
    end
  end

  dep :vetted_genes
  task :vetted_gene_list => :array do
    step(:vetted_genes).load.select("Vetted synthesis gene" => "true").column("Associated Gene Name").values.flatten.uniq
  end

  dep :vetted_gene_list
  dep ExTRI, :CollecTRI
  task :vetted_gene_enrichment => :tsv do
    require 'rbbt/statistics/hypergeometric'
    tsv = step(:CollecTRI).load.reorder "Transcription Factor (Associated Gene Name)", ["Target Gene (Associated Gene Name)"]
    vetted_genes = step(:vetted_gene_list).load
    all_genes = step(:vetted_genes).load.column("Associated Gene Name").values.flatten.uniq
    tsv.enrichment(vetted_genes, nil, :fdr => false, :background => all_genes)
  end

  dep :vetted_gene_list
  dep ExTRI, :CollecTRI
  task :vetted_gene_overlap => :tsv do
    tsv = step(:CollecTRI).load.reorder "Transcription Factor (Associated Gene Name)", ["Target Gene (Associated Gene Name)"]
    vetted_genes = step(:vetted_gene_list).load
    all_genes = step(:vetted_genes).load.column("Associated Gene Name").values.flatten.uniq
    res = TSV.setup({}, :key_field => "Statistic", :fields => ["Value"], :type => :single, :cast => :to_f)

    res["Total genes"] = all_genes.length
    res["Total vetted"] = vetted_genes.length
    res["CollecTRI targets"] = tsv.values.flatten.uniq.length
    res["Overlap"] = (tsv.values.flatten.uniq & vetted_genes).length
    res["Overlap %"] = ((100.0) * (tsv.values.flatten.uniq & vetted_genes).length) / vetted_genes.length
    res["All CollecTRI targets"] = (tsv.values.flatten.uniq & all_genes).length
    res["All Overlap %"] = ((100.0) * (tsv.values.flatten.uniq & all_genes).length) / all_genes.length

    res
  end

  dep :degradation_genes
  task :degradation_gene_list => :array do
    step(:degradation_genes).load.select("Degradation gene" => "true").column("Associated Gene Name").values.flatten.uniq
  end

  dep :degradation_gene_list
  dep ExTRI, :CollecTRI
  task :degradation_gene_overlap => :tsv do
    tsv = step(:CollecTRI).load.reorder "Transcription Factor (Associated Gene Name)", ["Target Gene (Associated Gene Name)"]
    vetted_genes = step(:degradation_gene_list).load
    all_genes = step(:degradation_genes).load.column("Associated Gene Name").values.flatten.uniq
    res = TSV.setup({}, :key_field => "Statistic", :fields => ["Value"], :type => :single, :cast => :to_f)

    res["Total genes"] = all_genes.length
    res["Total vetted"] = vetted_genes.length
    res["CollecTRI targets"] = tsv.values.flatten.uniq.length
    res["Overlap"] = (tsv.values.flatten.uniq & vetted_genes).length
    res["Overlap %"] = ((100.0) * (tsv.values.flatten.uniq & vetted_genes).length) / vetted_genes.length
    res["All CollecTRI targets"] = (tsv.values.flatten.uniq & all_genes).length
    res["All Overlap %"] = ((100.0) * (tsv.values.flatten.uniq & all_genes).length) / all_genes.length

    res
  end

end
