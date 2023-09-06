
module AGS

  TREATMENTS = %w(DMSO FiveZ INT_FiveZ_PI INT_PD_PI PD PI)
  TIME_POINTS = [1,2,4,8,24]
  

  dep :vetted_genes, :jobname => "Default"
  input :treatment, :select, "Treatment", nil, :select_options => TREATMENTS, :required => true
  input :time_point, :select, "Timepoint", nil, :select_options => TIME_POINTS, :required => true
  input :synthesis, :boolean, "Filter for synthesis genes"
  input :dynamic, :boolean, "Filter for dynamic cluster profiles"
  input :data_type, :select, "Values to use", :ext_fc, :select_options => %w(fc ext_fc binary)
  input :dynamic_subset, :select, "Increase, decrease or both", :both, :select_options => %w(increase decrease both)
  task :decoupler_matrix => :tsv do |treatment,time_point,synthesis,dynamic,data_type,dynamic_subset|
    time_point = time_point.to_i
    sample_name = [treatment, time_point] * "-T"
    matrix = TSV.setup({}, :key_field => "Associated Gene Name", :fields => [sample_name], :type => :single, :cast => :to_f)

    fc_fields = TIME_POINTS.collect{|t| ["FC_#{treatment}", "T#{t}"] * "." }
    TSV.traverse step(:vetted_genes), :into => matrix do |ens,values,fields|
      NamedArray.setup(values, fields)
      next if synthesis && ! values["Vetted gene"].include?("true")
      name = values["Associated Gene Name"].first
      fc_values = values.values_at(*fc_fields).collect{|v| v.first.to_f}
      cluters = values[treatment + ": FC clusters"]
      increase = cluters.include? "start increase #{time_point}h"
      decrease = cluters.include? "start decrease #{time_point}h"

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

      time_point_index =  TIME_POINTS.index time_point

      value = case data_type.to_sym
              when :fc
                fc = fc_values[time_point_index]
                fc -= fc_values[time_point_index-1] unless time_point_index == 0
                fc
              when :ext_fc
                fc = (time_point_index == fc_values.length - 1) ? fc_values[time_point_index] : fc_values[time_point_index+1]
                fc -= fc_values[time_point_index-1] unless time_point_index == 0
                fc
              when :binary
                if increase
                  1
                elsif decrease
                  -1
                else
                  0
                end
              else
                raise
              end

      [name, value]
    end

    matrix.transpose("Sample")
  end

  dep :decoupler_matrix
  dep SaezLab, :regulome, :jobname => "Default"
  dep_task :timepoint_decoupler_pre, SaezLab, :decoupler, :matrix => :decoupler_matrix, :network => :regulome, :min_n => 3

  dep :timepoint_decoupler_pre
  input :threshold, :float, "P-value threshold", 0.5
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
          [true, false].each do |dynamic|
            [:fc, :binary, :ext_fc].each do |data_type|
              [:mayority, :one, :two, :three].each do |synthesis_criteria|
                setup_name = [treatment, time_point, data_type, synthesis_criteria, "synth=#{synthesis}", "dynamic=#{dynamic}"] * "-"
                setup_options = {:treatment => treatment, :time_point => time_point, :synthesis => synthesis, :dynamic => dynamic, :data_type => data_type, :synthesis_criteria => synthesis_criteria}
                jobs << {:inputs => options.merge(setup_options), :jobname => setup_name}
              end
            end
          end
        end
      end
    end
    jobs
  end
  task :timepoint_decoupler_suite => :tsv do
    dependencies.inject(nil) do |acc,dep|
      next if dep.error?
      tsv = dep.load
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
