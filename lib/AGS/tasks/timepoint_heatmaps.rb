module AGS


  dep :timepoint_decoupler_suite
  task :timepoint_heatmaps => :array do 
    suite = step(:timepoint_decoupler_suite).path.tsv :cast => :to_f
    [true, false].each do |synthesis|
      [true, false].each do |dynamic|
        [:fc, :binary, :ext_fc].each do |data_type|
          [:mayority, :one, :two, :three].each do |synthesis_criteria|
            dir_name = [data_type, synthesis_criteria, "synth=#{synthesis}", "dynamic=#{dynamic}"] * "-"
            setup_name = [treatment, time_point, data_type, synthesis_criteria, "synth=#{synthesis}", "dynamic=#{dynamic}", "finegrained_degradation=#{finegrained_degradation}"] * "-"

            TIME_POINTS.each do |time_point|
              dir = file(dir_name)
              Open.mkdir dir
              target = dir[time_point.to_s + "h.png"]
              fields = TREATMENTS.collect do |treatment|
                field_name = [treatment, time_point, data_type, synthesis_criteria, "synth=#{synthesis}", "dynamic=#{dynamic}"] * "-"
                field_name
              end
              data = suite.slice(fields)
              data = data.select{|k,vs| vs.reject{|v| v.to_f == 0 }.any? }
              next if data.size < 2
              data.fields = data.fields.collect{|f| f.split("-").first }
              R::PNG.plot(target, data, <<-EOF, 1000, 500 + data.size * 12)
rbbt.require('gplots')
data[is.na(data)] = 0
m = as.matrix(data)
heatmap.2(m, dendrogram='row', Colv=FALSE, trace='none', col = terrain.colors(256), density.info="none", keysize = 0.75)
              EOF
            end
          end
        end
      end
    end
    Dir.glob(files_dir + "/**")
  end
end
