module AGS

  desc <<-EOF
Returns, for each treatment, the time timepoints in  which a 
gene starts increasing or decreasing
  EOF
  dep :change_offsets
  input :gene, :string, "Gene symbol to get offset information from"
  task :gene_offsets => :json do |gene|
    tsv = step(:change_offsets).path.tsv persist: true
    raise ParameterException, "Gene not found" unless tsv.include? gene
    tsv[gene].to_hash
  end
end
