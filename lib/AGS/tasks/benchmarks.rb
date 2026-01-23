#module AGS
#
#
#  dep :decoupler_matrix, :jobname => "Default"
#  input :ExTRI2_regulome, :boolean, "Use ExTRI2 regulome", true
#  input :test_regime, :array, "Test regime", [90,80,65,50,30]
#  input :experiments, :integer, "Number of experiments per regime", 10
#  task :suite => :tsv do |test_regime,experiments|
#    matrix = step(:decoupler_matrix).load
#    genes = matrix.fields
#    size = genes.length
#    input = file(:input)
#    test_regime.each do |percent|
#      experiments.times do |i|
#        filename = input[[percent, i] * "-"]
#        subset = genes.sample((size * percent) / 100)
#        Open.write(filename, matrix.slice(subset).to_s)
#      end
#    end
#  end
#
#  dep SaezLab, :regulome, :jobname => "Default"
#  input :test_regime, :array, "Test regime", [90,80,65,50,30]
#  input :experiments, :integer, "Number of experiments per regime", 10
#  dep :suite, :test_regime => :test_regime, :experiments => :experiments
#  dep SaezLab, :decoupler_activities, :matrix => :placeholder, :network => :regulome, :min_n => 3 do |jobname,options,dependencies|
#    jobs = []
#    suite = dependencies.flatten.last
#    test_regime = options[:test_regime]
#    experiments = options[:experiments]
#    test_regime.each do |percent|
#      experiments.times do |i|
#        name = [percent, i] * "-"
#        matrix_file = suite.file(:input)[name]
#        jobs << {inputs: options.merge(:matrix => matrix_file), jobname: jobname + "-large-#{name}"}
#      end
#    end
#    jobs
#  end
#  task :run_suite => :tsv do
#    dependencies[2..-1].inject(nil) do |acc,dep|
#      tsv = dep.load
#      tsv.fields = [dep.clean_name]
#      if acc.nil?
#        acc = tsv
#      else
#        acc = acc.attach tsv, complete: true
#      end
#    end
#  end
#
#  dep :run_suite, :syntesis => false, :dynamic => false, :jobname => "All"
#  dep :run_suite, :syntesis => false, :dynamic => true, :jobname => "Dynamic"
#  task :side_by_side => :tsv do
#    large, small = dependencies.collect{|d| d.load}
#    tsv = large.attach small
#  end
#
#  SIDE_BY_SIDE_OPTIONS={:treatment => "INT_PD_PI", :time_point => 1, :data_type => :ext_fc}
#  dep :side_by_side, SIDE_BY_SIDE_OPTIONS
#  extension :png
#  task :side_by_side_plot => :binary do
#    tsv = step(:side_by_side).path.tsv :cast => :to_f
#    R::PNG.plot(self.tmp_path, tsv, <<-EOF, 1000, tsv.size * 12)
#rbbt.require('gplots')
#m = as.matrix(data)
#heatmap.2(m,dendrogram='none', Colv=FALSE, trace='none', col = terrain.colors(256), density.info="none", keysize = 0.75)
#    EOF
#    nil
#  end
#
#  dep :side_by_side, SIDE_BY_SIDE_OPTIONS
#  extension :svg
#  task :side_by_side_svg => :text do
#    tsv = step(:side_by_side).path.tsv :cast => :to_f
#    R::SVG.plot(self.tmp_path, tsv, <<-EOF, 8, tsv.size / 5 )
#rbbt.require('gplots')
#m = as.matrix(data)
#heatmap.2(m,dendrogram='none', Colv=FALSE, trace='none', col = terrain.colors(256), density.info="none", keysize = 0.75)
#    EOF
#    nil
#  end
#
#  dep :side_by_side, SIDE_BY_SIDE_OPTIONS
#  task :side_by_side_binary => :tsv do
#    tsv = step(:side_by_side).load
#    fields = tsv.fields
#    fields.each do |f|
#      tsv.process f do |v|
#        if v.nil?
#          nil
#        elsif v > 0
#          1
#        else
#          -1
#        end
#      end
#    end
#    tsv.add_field "Positive" do |k,v|
#      values = fields.collect{|f| v[f] }
#      values.select{|v| v == 1 }.length
#    end
#    tsv.add_field "Negative" do |k,v|
#      values = fields.collect{|f| v[f] }
#      values.select{|v| v == -1 }.length
#    end
#
#    tsv.add_field "Discrepancies" do |k,v|
#      pos, neg = v.values_at "Positive", "Negative"
#
#      [pos, neg].min.to_f / (pos + neg)
#    end
#  end
#end
