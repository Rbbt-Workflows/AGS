# .source
# .source/abundance
# .source/abundance/5Z_0-5h.txt
# .source/abundance/PI5Z_2h.txt
# .source/abundance/PI5Z_8h.txt
# .source/abundance/PI_8h.txt
# .source/abundance/PI_2h.txt
# .source/abundance/5Z_8h.txt
# .source/abundance/PD_2h.txt
# .source/abundance/PI5Z_0.5h.txt
# .source/abundance/PIPD_0.5h.txt
# .source/abundance/PIPD_8h.txt
# .source/abundance/PI_0.5h.txt
# .source/abundance/5Z_2h.txt
# .source/abundance/PD_0-5h.txt
# .source/abundance/PIPD_2h.txt
# .source/abundance/PD_8h.txt
# .source/RE__Proteomics_data.zip
# .source/ptm
# .source/ptm/PTM_collapsed_log2_norm_imp_18517.txt
# .source/ptm/dDIA_PHOS_log2_filter75perc_18517psites.txt

require "scout"

Workflow.require_workflow 'AGS'

source = Path.setup '.source'
sss 0

uni2name = Organism.identifiers(AGS.organism).index target: 'Associated Gene Name', fields:  ['UniProt/SwissProt Accession'], persist: true

file = Path.setup ".source/ptm/PTM_collapsed_log2_norm_imp_18517.txt"
tsv = file.tsv header_hash: '', key_field: 'PTM_collapse_key', type: :list
new = tsv.annotate({})
tsv.each do |uni_multi, values|
  uni_multi.split(";").each do |uni|
    name = uni2name[uni]
    next if name.nil?
    new[name] = values
  end
end
new.key_field = "Associated Gene Name"

file = Path.setup ".source/ptm/dDIA_PHOS_log2_filter75perc_18517psites.txt"
tsv = file.tsv header_hash: '', key_field: 'T: T: PTM_collapse_key', type: :list
Log.tsv tsv

source.ptm.glob.each do |file|
  tsv = begin
          file.tsv header_hash: '', key_field: 'PTM_collapse_key', type: :list
        rescue
          file.tsv header_hash: '', key_field: 'T: T: PTM_collapse_key', type: :list
        end

  new = tsv.annotate({})
  tsv.each do |key, values|
    names, ptm, mod  = key.split('_')

    names.split(";").each do |name|
      new[[name, ptm, mod]*"_"] = values.collect{|v| v == 'NaN' ? nil : v }
    end
  end
  new.key_field = "PTM collapse key"

  Open.write "ptm/#{file.basename}", new.to_s
end

raise

source.abundance.glob.each do |file|
  treatment, time = file.basename.remove_extension.split('_')
  tsv = file.tsv header_hash: '', key_field: 'T: UniprotID', type: :list
  new = tsv.annotate({})
  tsv.each do |uni_multi, values|
    uni_multi.split(";").each do |uni|
      name = uni2name[uni]
      next if name.nil?
      new[name] = values
    end
  end
  new.key_field = "Associated Gene Name"

  treatment = 'INT_PD_PI' if treatment == 'PIPD'
  treatment = 'INT_FiveZ_PI' if treatment == 'PI5Z'

  target = [treatment, time[0..-2]]  * '-T' + '.tsv'
  Open.write 'abundance/'+target, new.to_s
end


