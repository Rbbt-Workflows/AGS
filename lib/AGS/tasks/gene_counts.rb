require 'rbbt/matrix'
require 'rbbt/matrix/differential'

module AGS
  task :gene_counts_pre => :tsv do
    header = nil
    tsv = {}
    TSV.traverse Rbbt.data.gene_counts["gene_counts.tsv"], :type => :array do |line|
      parts = line.split("\t")
      if header.nil?
        header = parts
      else
        key, *values = parts
        tsv[key] = values
      end
    end

    tsv = TSV.setup(tsv, :key_field => "Ensembl Gene ID", :fields => header)
    tsv.slice(tsv.fields - ["63", "101"])
  end

  task :samples => :tsv do
    header = nil
    tsv = {}
    TSV.traverse Rbbt.data.gene_counts["sample_info.tsv"], :type => :array do |line|
      parts = line.split("\t")
      if header.nil?
        header = parts[1..-1]
      else
        key, *values = parts
        tsv[key] = values
      end
    end

    tsv = TSV.setup(tsv, :key_field => "code", :fields => header, :type => :list)
  end

  dep :samples
  dep :gene_counts_pre
  task :gene_counts => :tsv do
    samples = step(:samples).load.slice %w(treatment Timepoint)
    samples.add_field "Treatment code" do |k,values|
      treatment, timepoint = values
      [treatment, "T" + timepoint.to_i.to_s] * "-"
    end
    sample_treatments = samples.slice("Treatment code").to_single
    counts = step(:gene_counts_pre).load

    treatments = sample_treatments.values.uniq.sort

    tsv = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => treatments, :type => :double, :cast => :to_f)

    samples = counts.fields
    counts.each do |gene,values|
      res = [[]] * treatments.length
      samples.zip(values).each do |sample,value|

        treatment = sample_treatments[sample]
        index = treatments.index(treatment)
        res[index] += [value]
      end
      tsv[gene] = res
    end
    tsv
  end


  helper :unzip_repeats do |tsv,sample|
    sample_column = tsv.column(sample)
    sample_column.unnamed = true
    sample_tsv = sample_column.to_flat
    sample_tsv.type = :list
    sample_tsv.fields = (1..sample_tsv.values.first.length).to_a.collect{|rep| sample + '_' + rep.to_s}
    sample_tsv
  end

  dep :gene_counts, :jobname => "Default"
  input :case_sample, :string, "Case sample name"
  input :control_sample, :string, "Control sample name"
  task :differential => :tsv do |case_sample,control_sample|

    gene_counts = step(:gene_counts).load

    case_tsv = unzip_repeats gene_counts, case_sample
    control_tsv = unzip_repeats gene_counts, control_sample

    case_sample_reps = case_tsv.fields.dup
    control_sample_reps = control_tsv.fields.dup

    gene_counts_by_rep = case_tsv.attach control_tsv

    matrix_file = file('matrix.tsv')
    Open.write(matrix_file, gene_counts_by_rep.to_s)

    matrix = RbbtMatrix.new matrix_file, nil, 'fpkm'

    matrix.differential(case_sample_reps, control_sample_reps, self.tmp_path)

    nil
  end

  dep :differential, :case_sample => :placeholder, :control_sample => :placeholder do |jobname,options|
    jobs = []
    %w(DMSO 5Z PD_PI 5Z_PI PD PI ).each do |treatment|
      %w(1 2 4 8 24).each do |case_time|
        %w(0 1 2 4 8 24).each do |control_time|
          next if case_time.to_i <= control_time.to_i

          case_sample = [treatment, "T" + case_time] * "-"

          if control_time == "0"
            control_sample = "untreated-T0"
          else
            control_sample = [treatment, "T" + control_time] * "-"
          end

          name = [case_sample, control_sample] * " <-> "
          jobs << {:inputs => options.merge(:case_sample => case_sample, :control_sample => control_sample), :jobname => name}
        end
      end
    end

    jobs
  end
  task :all_differential => :tsv do 
    tsv = Log::ProgressBar.with_bar dependencies.length do |bar|
      dependencies.inject(nil) do |acc,dep|
        tsv = dep.load
        ratio = tsv.column("log2FoldChange")
        pvalues = tsv.column("pvalue")

        name = dep.clean_name
        ratio.fields = [name]
        pvalues.fields = [name + " (pvalues)"]

        if acc.nil?
          acc = ratio
        else
          acc = acc.attach ratio, :complete => true
        end

        acc = acc.attach pvalues, :complete => true

        bar.tick

        acc
      end
    end

    tsv = tsv.change_key "Associated Gene Name"

    tsv
  end

  dep :all_differential
  task :log2foldchanges => :tsv do
    tsv_file = step(:all_differential).join.path
    fields = TSV.parse_header(tsv_file).fields
    fields = fields.select{|f| f !~ /pvalues/ }
    TSV.open(tsv_file, :fields => fields)
  end

  dep :all_differential
  task :pvalues_BSC => :tsv do
    tsv_file = step(:all_differential).join.path
    fields = TSV.parse_header(tsv_file).fields
    fields = fields.select{|f| f =~ /pvalues/ }
    tsv = TSV.open(tsv_file, :fields => fields)
    tsv.fields = tsv.fields.collect{|f| f.sub(' (pvalues)', '') }
    tsv
  end

  dep :all_differential
  input :threshold, :float, "P-value threshold", 0.1
  task :significant_log2foldchanges => :tsv do |threshold|
    tsv_file = step(:all_differential).join.path

    parser = TSV::Parser.new tsv_file
    fc_fields = parser.fields.select{|f| f !~ /pvalues/ }

    dumper = TSV::Dumper.new parser.options.merge(:fields => fc_fields)
    dumper.init
    TSV.traverse parser, :into => dumper do |key,values|
      experiments = values.length / 2
      res = []
      experiments.times do |i|
        fc = values[i]
        pv = values[i + 1]
        res << ((pv && pv > 0 && pv < threshold) ? fc : nil)
      end
      [key, res]
    end
  end

  dep :significant_log2foldchanges
  task :significant_foldchanges => :tsv do 
    tsv_file = step(:significant_log2foldchanges).join.path

    parser = TSV::Parser.new tsv_file

    dumper = TSV::Dumper.new parser.options
    dumper.init
    TSV.traverse parser, :into => dumper do |key,values|
      res = values.collect{|v| v.nil? ? nil : 2**v }
      [key, res]
    end

  end

end
