module AGS
  task :dbTFs => :array do
    tsv = Rbbt.data["EXTRI2-regulome_dbTFs_coTFs_230425"].tsv
    tsv.select("dbTF" => "dbTF").keys
  end  
end
