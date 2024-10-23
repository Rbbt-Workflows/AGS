
#require 'AGS/tasks/decoupler_old'
module AGS

  def self.setup_name(treatment, time_point, data_type, synthesis_criteria, synthesis, dynamic, finegrained_degradation)
    [treatment, time_point, data_type, synthesis_criteria, "synth=#{synthesis}", "dynamic=#{dynamic}", "finegrained_degradation=#{finegrained_degradation}"] * "-"
  end

  def self.range_value(fc_values, cluster, time_point)
    range = cluster.split(" ").last
    range_start, range_end = range.split("-")
    range_start = range_start.to_i
    range_end = range_end.nil? ? range_start : range_end.to_i
    return nil unless time_point == range_start
    i_start = TIME_POINTS.index(range_start)
    i_end = TIME_POINTS.index(range_end)

    fc = fc_values[i_start].to_f
    fc -= fc_values[i_start-1].to_f unless i_start == 0

    (i_start+1..i_end).each do |ext_i|
      fc += (fc_values[ext_i].to_f - fc_values[ext_i-1].to_f) / (TIME_POINTS[ext_i] - TIME_POINTS[ext_i-1])
    end
    return fc
  end

  dep :vetted_genes
  input :treatment, :select, "Treatment", nil, :select_options => TREATMENTS, :required => true
  input :time_point, :select, "Timepoint", nil, :select_options => TIME_POINTS, :required => true
  input :data_type, :select, "Values to use", :ext_fc, :select_options => %w(fc fc0 ext_fc binary range)
  task :timepoint_matrix => :tsv do |treatment,time_point,data_type|
    time_point = time_point.to_i
    sample_name = [treatment, time_point] * "-T"
    matrix = TSV.setup({}, :key_field => "Associated Gene Name", :fields => [sample_name], :type => :single, :cast => :to_f)

    fc_fields = TIME_POINTS.collect{|t| ["FC_#{treatment}", "T#{t}"] * "." }
    TSV.traverse step(:vetted_genes), :into => matrix do |name,values,fields|
      name = name.first if Array === name
      NamedArray.setup(values, fields)
      fc_values = values.values_at(*fc_fields).collect{|v| v.first}
      clusters = values[treatment + ": FC clusters"]

      time_point_index =  TIME_POINTS.index time_point

      next if fc_values[time_point_index].nil?

      value = case data_type.to_sym
              when :fc0
                fc = fc_values[time_point_index].to_f
              when :fc
                fc = fc_values[time_point_index].to_f
                fc -= fc_values[time_point_index-1].to_f unless time_point_index == 0
                fc
              when :ext_fc
                fc = (time_point_index == fc_values.length - 1) ? fc_values[time_point_index].to_f : fc_values[time_point_index+1].to_f
                fc -= fc_values[time_point_index-1].to_f unless time_point_index == 0
                fc
              when :binary
                if increase
                  1
                elsif decrease
                  -1
                else
                  0
                end
              when :range
                fc = nil
                clusters.each{|cluster| 
                  fc = AGS.range_value(fc_values, cluster, time_point)
                  break if ! fc.nil?
                  #range = cluster.split(" ").last
                  #range_start, range_end = range.split("-")
                  #range_start = range_start.to_i
                  #range_end = range_end.nil? ? range_start : range_end.to_i
                  #next unless time_point == range_start
                  #i_start = TIME_POINTS.index(range_start)
                  #i_end = TIME_POINTS.index(range_end)

                  #fc = fc_values[i_start].to_f
                  #fc -= fc_values[i_start-1].to_f unless i_start == 0

                  #(i_start+1..i_end).each do |ext_i|
                  #  fc += (fc_values[ext_i].to_f - fc_values[ext_i-1].to_f) / (TIME_POINTS[ext_i] - TIME_POINTS[ext_i-1])
                  #end
                }
                fc.nil? ? 0 : fc
              else
                raise
              end

      [name, value]
    end

    matrix.transpose("Sample")
  end

  dep :vetted_genes, :jobname => "Default"
  input :treatment, :select, "Treatment", nil, :select_options => TREATMENTS, :required => true
  input :time_point, :select, "Timepoint", nil, :select_options => TIME_POINTS, :required => true
  input :synthesis, :boolean, "Filter for synthesis genes", false
  input :dynamic, :boolean, "Filter for dynamic cluster profiles"
  input :dynamic_subset, :select, "Increase, decrease or both", :both, :select_options => %w(increase decrease both)
  input :finegrained_degradation, :boolean, "Filter for finegrained dynamic degradation", true
  input :target_subset, :array, "Subset of genes to consider" 
  task :decoupler_targets => :array do |treatment,time_point,synthesis,dynamic,dynamic_subset,finegrained_degradation,target_subset|

    noisy_genes = Rbbt.data.noisy_genes.list
    TSV.traverse step(:vetted_genes), :into => :stream do |name,values,fields|
      name = name.first if Array === name
      NamedArray.setup(values, fields)
      next if synthesis && ! values["Vetted synthesis gene"].include?("true")
      next if noisy_genes.include?(name)
      next if target_subset && ! target_subset.include?(name)
      cluters = values[treatment + ": FC clusters"]
      #increase = cluters.include? "increase #{time_point}"
      #decrease = cluters.include? "decrease #{time_point}"
      increase = cluters.select{|c| c == "increase #{time_point}" || c.include?("increase #{time_point}-") }.any?
      decrease = cluters.select{|c| c == "decrease #{time_point}" || c.include?("decrease #{time_point}-") }.any?

      if dynamic
        case dynamic_subset.to_sym
        when :both
          next unless increase || decrease
        when :increase
          next unless increase
        when :decrease
          next unless decrease
        else
          raise ParameterException, "Unkown dynamic_subset #{dynamic_subset}"
        end
      end

      if finegrained_degradation and values["Degradation timepoints"].include?([treatment, time_point]*":")
        next
      end

      name
    end
  end

  dep SaezLab, :regulome
  dep :decoupler_targets
  task :filtered_regulome => :tsv do
    targets = step(:decoupler_targets).load
    interesting_tfs = Rbbt.data.interesting_tfs.list
    dumper = TSV::Dumper.new step(:regulome).load.options
    dumper.init
    TSV.traverse step(:regulome), :into => dumper do |id,values|
      tf, tg, weight = values
      next unless targets.include?(tg)
      next unless interesting_tfs.include?(tf)
      [id, [tf, tg, weight]]
    end
  end

  dep :timepoint_matrix
  dep :filtered_regulome
  dep_task :timepoint_decoupler_pre, SaezLab, :decoupler, :matrix => :timepoint_matrix, :network => :filtered_regulome, :min_n => 5

  dep :timepoint_decoupler_pre
  input :threshold, :float, "P-value threshold", 0.05
  task :timepoint_decoupler => :tsv do |threshold|
    tsv = step(:timepoint_decoupler_pre).load.transpose("Associated Gene Name")

    genes = tsv.keys.collect{|f| f.split(" ").first}.uniq

    new = tsv.annotate({})
    genes.each do |gene|
      if tsv.include? gene
        values = tsv[gene]
        pvalues = tsv[gene + " (pvalue)"]
      else
        values = tsv[gene + " (consensus_estimate)"]
        pvalues = tsv[gene + " (consensus_pvals)"]
      end
      new_values = values.zip(pvalues).collect do |v,p|
        p.to_f < threshold ? v : 0
      end
      new[gene] = new_values
    end
    new
  end

  dep :timepoint_decoupler, :canfail => true do |jobname,options|
    jobs = []
    TREATMENTS.each do |treatment|
      TIME_POINTS.each do |time_point|
        [true, false].each do |synthesis|
          next if synthesis
          [true, false].each do |dynamic|
            next if dynamic
            [true, false].each do |finegrained_degradation|
              next if finegrained_degradation
              [:fc, :fc0, :binary, :ext_fc].each do |data_type|
                [:mayority, :one, :two, :three].each do |synthesis_criteria|
                  next if synthesis_criteria != :mayority
                  setup_name = AGS.setup_name(treatment, time_point, data_type, synthesis_criteria, synthesis, dynamic, finegrained_degradation)
                  setup_options = {:treatment => treatment, :time_point => time_point, :synthesis => synthesis, :dynamic => dynamic, :data_type => data_type, :synthesis_criteria => synthesis_criteria, :finegrained_degradation => finegrained_degradation}
                  jobs << {:inputs => options.merge(setup_options), :jobname => setup_name}
                end
              end
            end
          end
        end
      end
    end
    jobs
  end
  task :timepoint_decoupler_suite => :tsv do
    not_expressed = Rbbt.data.not_expressed.list
    noisy_genes = Rbbt.data.noisy_genes.list
    dependencies.inject(nil) do |acc,dep|
      next if dep.error?
      tsv = dep.load
      tsv = tsv.subset(tsv.keys - not_expressed - noisy_genes)
      tsv.fields = [dep.clean_name]
      if acc.nil?
        acc = tsv
      else
        acc = acc.attach tsv, :complete => true
      end
      acc
    end
  end

  dep :timepoint_decoupler_suite
  task :timepoint_decoupler_excel => :array do
    require 'rbbt/tsv/excel'
    tsv = step(:timepoint_decoupler_suite).load
    output = file('excel')
    Open.mkdir output
    TREATMENTS.each do |treatment|
      treatment_fields = tsv.fields.select{|f| f.start_with?(treatment) }
      target = output[treatment + ".xlsx"]
      tsv.slice(treatment_fields).xlsx(target)
    end

    output.glob("*")
  end

  input :condition1, :string
  input :condition2, :string
  input :direction, :select, "Direction of regulation", :up, :select_options => %w(up down both)
  extension :png
  task :timepoint_decoupler_venn => :binary do |condition1,condition2,direction|
    tsv = TSV.open("/home/mvazque2/.rbbt/var/jobs/AGS/timepoint_decoupler_suite/Default.tsv")
    up1 = tsv.select(condition1){|v| v.to_f > 0}.keys
    down1 = tsv.select(condition1){|v| v.to_f < 0}.keys
    up2 = tsv.select(condition2){|v| v.to_f > 0}.keys
    down2 = tsv.select(condition2){|v| v.to_f < 0}.keys

    set_info :up1, up1
    set_info :down1, down1
    set_info :up2, up2
    set_info :down2, down2

    data = TSV.setup({}, :key_field => "Gene", :fields => [condition1, condition2], :type => :list)

    case direction.to_sym
    when :up
      up1.each{|g| v = data[g] ||= [false, false]; v[0] = true; data[g] = v }
      up2.each{|g| v = data[g] ||= [false, false]; v[1] = true; data[g] = v }
    when :down
      down1.each{|g| v = data[g] ||= [false, false]; v[0] = true; data[g] = v }
      down2.each{|g| v = data[g] ||= [false, false]; v[1] = true; data[g] = v }
    when :both
      up1.each{|g| v = data[g] ||= [false, false]; v[0] = true; data[g] = v }
      up2.each{|g| v = data[g] ||= [false, false]; v[1] = true; data[g] = v }
      down1.each{|g| v = data[g] ||= [false, false]; v[0] = true; data[g] = v }
      down2.each{|g| v = data[g] ||= [false, false]; v[1] = true; data[g] = v }
    end

    R::PNG.plot(self.tmp_path, data, <<-EOF)
rbbt.plot.venn(data)
    EOF
    nil
  end


end
