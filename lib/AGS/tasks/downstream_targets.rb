module AGS
  TREATMENT_TARGETS ||= {
    'FiveZ' => %w(MAP3K7),
    'PD' => %w(MAP2K1 MAP2K2 ),
    'PI' => %w(PIK3CG PIK3CA PIK3R1 PIK3CB)
  }

  TREATMENT_TARGETS["INT_FiveZ_PI"] = TREATMENT_TARGETS["FiveZ"] + TREATMENT_TARGETS["PI"]
  TREATMENT_TARGETS["INT_PD_PI"] = TREATMENT_TARGETS["PD"] + TREATMENT_TARGETS["PI"]

  dep :regulome
  input :steps_max, :integer, "Max steps", 3
  task :downstream_targets => :yaml do |steps|

    require 'rbbt/network/paths'

    sources = TREATMENT_TARGETS.values.flatten.uniq
    reg = step(:regulome).load
    targets = reg.column("source").values.flatten.uniq

    signor_all = Signor.protein_protein.tsv 
    signor_all.namespace = "Hsa/may2017"
    signor_all.identifiers = Organism.identifiers(signor_all.namespace)
    signor_all = signor_all.change_key "Associated Gene Name"
    signor_all = signor_all.swap_id 'Target (UniProt/SwissProt Accession)', "Associated Gene Name"

    signor = signor_all.slice([0])
    signor = signor.to_flat

    signor_pairs = signor_all.unzip.slice("Effect").to_single

    found_targets = {}
    sources.each do |source|
      step = 0
      found = []
      current = [source]
      type = [1]
      seen = []
      while current.any?
        iif [current, found, seen]
        ncurrent = []
        ntype = []
        current.zip(type).each do |c,t|
          next unless signor.include?(c)
          ns = signor[c] - seen
          nt = ns.collect{|n| 
            dir = signor_pairs[[c,n] * ":"]
            dir.split(/\W/).first == "down" ? -1 : 1
          }
          ncurrent += ns
          if t == 1
            ntype += nt
          else
            ntype += nt.collect{|v| - v }
          end
        end
        break if ncurrent.empty?
        seen += current
        current = ncurrent
        type = ntype
        current, type = Misc.zip_fields current.zip(type).uniq
        step += 1
        current.zip(type).each do |gene,t|
          found << [gene, t] if targets.include?(gene)
        end
        found = found.uniq
        break if step > steps
      end

      found_targets[source] = found
    end

    found_targets
  end

  dep :downstream_targets
  input :timepoint_suite, :tsv
  task :downstream_consistency => :array do |suite|
    downstream_targets = step(:downstream_targets).load
    TREATMENT_TARGETS.collect do |treatment,sources|
      target_signs = {} 
      downstream_targets.values_at(*sources).uniq.each do |list|
        list.each do |target,sign|
          target_signs[target] ||= []
          target_signs[target] << sign
        end
      end
      targets = target_signs.keys
      [true, false].each do |synthesis|
        [true, false].each do |dynamic|
          [:fc, :binary, :ext_fc].each do |data_type|
            [:mayority, :one, :two, :three].each do |synthesis_criteria|
              dir_name = [data_type, synthesis_criteria, "synth=#{synthesis}", "dynamic=#{dynamic}"] * "-"
              dir = file(dir_name)
              Open.mkdir dir
              output = dir[treatment + '.tsv']

              time_point = "1"
              field_name = [treatment, time_point, data_type, synthesis_criteria, "synth=#{synthesis}", "dynamic=#{dynamic}"] * "-"

              data = suite.slice([field_name])
              data = data.subset(targets)
              data.add_field "Sign" do |target|
                target_signs[target] * "|"
              end
              Open.write(output, data.to_s)
            end
          end
        end
      end
    end
    Dir.glob(files_dir + "/**/*")
  end
end
