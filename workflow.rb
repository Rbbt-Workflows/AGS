require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/sources/organism'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/AGS'

Workflow.require_workflow "SaezLab"
Workflow.require_workflow "ExTRI2"

module AGS
  extend Workflow

  #TREATMENTS = %w(DMSO FiveZ INT_FiveZ_PI INT_PD_PI PD PI)
  TREATMENTS = %w(INT_PD_PI INT_FiveZ_PI PI FiveZ PD DMSO)
  TIME_POINTS = [1,2,4,8,24]
  
  def self.gprofiler
    @gprofiler ||= begin
                     RbbtPython.pyimport 'gprofiler'
                     RbbtPython.gprofiler.GProfiler.new(return_dataframe: true)
                   end
  end

  #self.gprofiler

  def self.organism
    "Hsa/feb2014"
  end

  def self.regulome
    @@regulome ||= SaezLab.job(:regulome).produce.path.tsv :key_field => "source", :fields => ["target"], :type => :flat
  end

  def self.meta_clusters
    @@meta_clusters ||= begin
                          meta = Rbbt.data["cluster_groups.tsv"].tsv :fields => %w(5Z_PI PD_PI PI 5Z PD DMSO)
                          meta.add_field "MetaCluster" do |cluster,values|
                            values.zip(meta.fields).collect do |meta,experiment|
                              next if meta.nil? || meta.empty?
                              [meta, experiment] * "-"
                            end.compact
                          end
                          meta
                        end
  end

  def self.super_clusters
    @@super_clusters ||= begin
                          meta = Rbbt.data["cluster_groups.tsv"].tsv 
                          meta.add_field "SuperCluster" do |cluster,values|
                            values.flatten.select{|v| v =~ /^[A-Z]\d?$/}
                          end
                          meta
                        end
  end

  def self.super_cluster_names
    @@super_cluster_names ||= Rbbt.data["super_cluster_names.tsv"].tsv :type => :single
  end


  def self.gene_clusters
    @@gene_clusters ||= begin
                          tsv = Rbbt.data["gene_modules.tsv"].tsv :key_field => "Associated Gene Name", :fields => ["Cluster"], :type => :double
                          tsv.attach super_clusters, :fields => ["SuperCluster"]
                          tsv.attach meta_clusters, :fields => ["MetaCluster"]
                        end
  end

  def self.gene_counts
    @@gene_counts ||= AGS.job(:gene_counts).run #.change_key "Associated Gene Name", :identifiers => Organism.identifiers(AGS.organism)
  end


  def self.fold_changes(fc_source = "NTNU")
    @@fold_changes ||= AGS.job(:fold_changes, nil, :fc_source => fc_source).run.transpose "Associated Gene Name"
  end

  def self.decoupler(fc_source = "NTNU")
    @@decoupler ||= AGS.job(:decoupler, nil, :fc_source => fc_source).run
  end


  input :fc_source, :select, "Source of fold_changes", "NTNU", :select_options => %w(BSC NTNU)
  dep :fold_changes_NTNU do |jobname,options|
    if options[:fc_source] == "BSC"
      nil
    else
      {:inputs => options, :jobname => jobname}
    end
  end
  dep :log2foldchanges do |jobname, options|
    if options[:fc_source] == "BSC"
      {:inputs => options, :jobname => jobname}
    else
      nil
    end
  end
  input :subset_effect, :select, "Effect considered", :both, :select_options => %w(both up down)
  input :conditions, :array, "Subset conditions" 
  task :fold_changes => :tsv do  |fc_source,subset_effect,conditions|
    tsv = case fc_source.to_s
          when "BSC"
            tsv = step(:log2foldchanges).load
            tsv.fields.each do |field|
              tsv.process field do |v|
                v.nil? ? 0 : v
              end
            end
            tsv.transpose("ID")
          when "NTNU"
            step(:fold_changes_NTNU).load
          else
            raise ParameterException
          end

    tsv = tsv.subset(conditions) if conditions && conditions.any?

    case subset_effect
    when 'up'
      tsv.fields.each do |field|
        tsv.process field do |v|
          v < 0 ? 0 : v
        end
      end
      tsv
    when 'down'
      tsv.fields.each do |field|
        tsv.process field do |v|
          v > 0 ? 0 : v
        end
      end
      tsv
    else
      tsv
    end
  end

  input :fc_source, :select, "Source of fold_changes", "NTNU", :select_options => %w(BSC NTNU)
  dep :pvalues_NTNU do |jobname,options|
    if options[:fc_source] == "BSC"
      nil
    else
      {:inputs => options, :jobname => jobname}
    end
  end
  dep :pvalues_BSC do |jobname,options|
    if options[:fc_source] == "BSC"
      {:inputs => options, :jobname => jobname}
    else
      nil
    end
  end
  input :subset_effect, :select, "Effect considered", :both, :select_options => %w(both up down)
  input :conditions, :array, "Subset conditions" 
  task :pvalues => :tsv do  |fc_source,subset_effect,conditions|
    tsv = case fc_source.to_s
          when "BSC"
            tsv = step(:pvalues).load
            tsv.fields.each do |field|
              tsv.process field do |v|
                v.nil? ? 0 : v
              end
            end
            tsv.transpose("ID")
          when "NTNU"
            step(:pvalues_NTNU).load
          else
            raise ParameterException
          end

    case subset_effect
    when 'up'
      tsv.fields.each do |field|
        tsv.process field do |v|
          v < 0 ? 0 : v
        end
      end
      tsv
    when 'down'
      tsv.fields.each do |field|
        tsv.process field do |v|
          v > 0 ? 0 : v
        end
      end
      tsv
    else
      tsv
    end

    tsv = tsv.slice(conditions) if conditions && conditions.any?

    tsv
  end


  dep :fold_changes
  dep SaezLab, :regulome
  dep_task :decoupler_pre, SaezLab, :decoupler, :matrix => :fold_changes, :network => :regulome

  dep :decoupler_pre
  input :threshold, :float, "P-value threshold", 0.5
  task :decoupler => :tsv do |threshold|
    tsv = step(:decoupler_pre).load.transpose("Associated Gene Name")

    genes = tsv.keys.reject{|g| g.include? 'pvalue'}

    new = tsv.annotate({})
    genes.each do |gene|
      values = tsv[gene]
      pvalues = tsv[gene + " (pvalue)"]
      new_values = values.zip(pvalues).collect do |v,p|
        p.to_f < threshold ? v : 0
      end
      new[gene] = new_values
    end
    new
  end

  dep :fold_changes, :fc_source => "NTNU"
  task :fold_changes_fc_one => :tsv do
    tsv = step(:fold_changes).load
    new = tsv.annotate({})
    tsv.keys.each do |k|
      experiment, time = k.split(".T")
      prev_time = case time
                  when "1"
                    next
                  when "2"
                    "1"
                  when "4"
                    "2"
                  when "8"
                    "4"
                  when "24"
                    "8"
                  end
      prev = [experiment, prev_time] * ".T"
      new_values = tsv[k].zip(tsv[prev]).collect{|c,p| c - p }
      new[k] = new_values
    end
    new
  end

  dep :fold_changes_fc_one
  dep SaezLab, :regulome
  dep_task :decoupler_pre_fc_one, SaezLab, :decoupler, :matrix => :fold_changes_fc_one, :network => :regulome

  dep :decoupler_pre_fc_one
  input :threshold, :float, "P-value threshold", 0.5
  task :decoupler_fc_one => :tsv do |threshold|
    tsv = step(:decoupler_pre_fc_one).load.transpose("Associated Gene Name")

    genes = tsv.keys.reject{|g| g.include? 'pvalue'}

    new = tsv.annotate({})
    genes.each do |gene|
      values = tsv[gene]
      pvalues = tsv[gene + " (pvalue)"]
      new_values = values.zip(pvalues).collect do |v,p|
        p.to_f < threshold ? v : 0
      end
      new[gene] = new_values
    end
    new
  end
end

require 'AGS/tasks/NTNU'
require 'AGS/tasks/gene_counts'
require 'AGS/tasks/change_starts'
require 'AGS/tasks/gene_clusters'
require 'AGS/tasks/INSPEcT'
require 'AGS/tasks/decoupler'
require 'AGS/tasks/offset'
require 'AGS/tasks/timepoint_heatmaps'
require 'AGS/tasks/downstream_targets'
require 'AGS/tasks/consistency'

require 'AGS/tasks/enrichment'

require 'AGS/tasks/excel'

require 'AGS/tasks/adhoc'

#require 'AGS/tasks/benchmarks'
require 'knowledge_base/AGS'
#require 'rbbt/entity/AGS'

Workflow.main = AGS
