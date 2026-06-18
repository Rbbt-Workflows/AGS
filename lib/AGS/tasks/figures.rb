module AGS

  FIGURE_TREATMENT_ORDER = %w(DMSO FiveZ INT_FiveZ_PI INT_PD_PI PD PI) unless const_defined?(:FIGURE_TREATMENT_ORDER)
  FIGURE_DRUG_TREATMENT_ORDER = %w(FiveZ INT_FiveZ_PI INT_PD_PI PD PI) unless const_defined?(:FIGURE_DRUG_TREATMENT_ORDER)
  FIGURE_TIME_POINTS = [1, 2, 4, 8, 24] unless const_defined?(:FIGURE_TIME_POINTS)

  FIGURE_TREATMENT_LABELS = {
    'DMSO' => 'DMSO',
    'FiveZ' => '5Z',
    'INT_FiveZ_PI' => '5Z+PI',
    'INT_PD_PI' => 'PD+PI',
    'PD' => 'PD',
    'PI' => 'PI'
  } unless const_defined?(:FIGURE_TREATMENT_LABELS)

  FIGURE_TREATMENT_COLORS = {
    'DMSO' => '#7F7F7F',
    'FiveZ' => '#4DAF4A',
    'INT_FiveZ_PI' => '#FF7F00',
    'INT_PD_PI' => '#E41A1C',
    'PD' => '#984EA3',
    'PI' => '#377EB8'
  } unless const_defined?(:FIGURE_TREATMENT_COLORS)

  FIGURE_DIRECTION_COLORS = {
    'up' => '#D73027',
    'down' => '#4575B4',
    'positive' => '#D73027',
    'negative' => '#4575B4',
    'both' => '#4D4D4D'
  } unless const_defined?(:FIGURE_DIRECTION_COLORS)

  FIGURE_ENRICHMENT_FILL_COLORS = {
    'none' => '#FFFFFF',
    'up_q05' => '#F6B6B6',
    'up_q1e4' => '#D73027',
    'up_q1e6' => '#7F0000',
    'down_q05' => '#B8D4F0',
    'down_q1e4' => '#4575B4',
    'down_q1e6' => '#08306B',
    'positive_q05' => '#F6B6B6',
    'positive_q1e4' => '#D73027',
    'positive_q1e6' => '#7F0000',
    'negative_q05' => '#B8D4F0',
    'negative_q1e4' => '#4575B4',
    'negative_q1e6' => '#08306B'
  } unless const_defined?(:FIGURE_ENRICHMENT_FILL_COLORS)

  FIGURE_COMBINATION_CATEGORY_LABELS = {
    'combination_earlier_than_both' => 'combination earlier than both',
    'combination_specific_at_time' => 'combination specific at time',
    'shared_with_both_same_sign' => 'shared with both, same sign',
    'shared_with_component1_same_sign' => 'shared with PI, same sign',
    'shared_with_component2_same_sign' => 'shared with second component, same sign',
    'sign_reversed_relative_to_component' => 'sign reversed relative to component'
  } unless const_defined?(:FIGURE_COMBINATION_CATEGORY_LABELS)

  FIGURE_SELECTED_GOSLIM_TERMS = [
    'DNA-templated transcription',
    'regulation of DNA-templated transcription',
    'ribosome biogenesis',
    'tRNA metabolic process',
    'mRNA metabolic process',
    'protein folding',
    'mitochondrion organization',
    'generation of precursor metabolites and energy',
    'mitotic cell cycle',
    'DNA replication',
    'DNA repair',
    'DNA recombination',
    'chromosome segregation',
    'lipid metabolic process',
    'transmembrane transport',
    'amino acid metabolic process',
    'cell adhesion',
    'signaling',
    'cell differentiation',
    'autophagy',
    'extracellular matrix organization',
    'cell motility'
  ] unless const_defined?(:FIGURE_SELECTED_GOSLIM_TERMS)

  def self.figure_treatment_label(treatment)
    FIGURE_TREATMENT_LABELS[treatment.to_s] || treatment.to_s
  end

  def self.figure_treatment_color(treatment)
    FIGURE_TREATMENT_COLORS[treatment.to_s] || '#000000'
  end

  def self.figure_treatment_levels(treatments = FIGURE_TREATMENT_ORDER)
    treatments.collect{|t| "'#{AGS.figure_treatment_label(t)}'" } * ','
  end

  def self.figure_treatment_color_scale(aesthetic = 'color', treatments = FIGURE_TREATMENT_ORDER)
    values = treatments.collect{|t| "'#{AGS.figure_treatment_label(t)}'='#{AGS.figure_treatment_color(t)}'" } * ','
    "scale_#{aesthetic}_manual(values=c(#{values}))"
  end

  def self.figure_direction_color_scale(aesthetic = 'fill')
    values = FIGURE_DIRECTION_COLORS.collect{|k,v| "'#{k}'='#{v}'" } * ','
    "scale_#{aesthetic}_manual(values=c(#{values}))"
  end

  def self.figure_enrichment_fill_scale
    values = FIGURE_ENRICHMENT_FILL_COLORS.collect{|k,v| "'#{k}'='#{v}'" } * ','
    "scale_fill_manual(values=c(#{values}), drop=FALSE, name='Direction/FDR')"
  end

  def self.figure_build_tsv(rows, key_field, fields)
    tsv = TSV.setup({}, :key_field => key_field, :fields => fields, :type => :list)
    rows.each do |row|
      row = row.dup
      key = row.shift
      tsv[key] = row
    end
    tsv
  end

  def self.figure_safe_float(value, default = 0.0)
    value = value.compact.first if Array === value
    return default if value.nil? || value.to_s.empty?
    value.to_f
  end

  def self.figure_safe_int(value, default = 0)
    value = value.compact.first if Array === value
    return default if value.nil? || value.to_s.empty?
    value.to_i
  end

  def self.figure_context_fields(prefix = 'FC')
    FIGURE_TREATMENT_ORDER.collect do |treatment|
      FIGURE_TIME_POINTS.collect{|time_point| "#{prefix}_#{treatment}.T#{time_point}" }
    end.flatten
  end

  def self.figure_tf_context_fields(tsv)
    FIGURE_TREATMENT_ORDER.collect do |treatment|
      FIGURE_TIME_POINTS.collect{|time_point| "#{treatment}-T#{time_point}" }
    end.flatten.select{|field| tsv.fields.include?(field) }
  end

  #def self.figure_ggplot_svg(tsv, r_code, width = 8, height = 5)
  #  R::SVG.ggplot tsv, r_code, width, height
  #end

  #def self.figure_base_svg(file, tsv, r_code, width = 11, height = 8)
  #  R::SVG.plot file, tsv, r_code, width, height
  #  nil
  #end

  helper :figure_ggplot do |tsv,r_code,width=8,height=5|
    Open.mkdir files_dir
    R::PNG.ggplot file('plot.png'), tsv, r_code, width, height
    R::SVG.ggplot tsv, r_code, width, height
  end

  helper :figure_base do |tsv,r_code,width=11,height=8|
    Open.mkdir files_dir
    pw = width
    ph = height
    begin
      Misc.insist do
        R::PNG.plot file('plot.png'), tsv, r_code, pw, ph
        pw *= 2 
        ph *= 2
      end
    rescue
    end
    R::SVG.plot self.tmp_path, tsv, r_code, width, height
    nil
  end

  def self.figure_functional_term_group(term, annotation = nil)
    term_up = term.to_s.upcase
    term_down = term.to_s.downcase

    if annotation.to_s =~ /^cancerhallmarks/
      return 'Proliferation and growth' if term_up =~ /PROLIFERATIVE|GROWTH|IMMORTALITY/
      return 'Cell death and suppression' if term_up =~ /CELL DEATH|GROWTH SUPPRESSORS/
      return 'Genome maintenance' if term_up =~ /GENOME/
      return 'Metabolism' if term_up =~ /METABOLISM/
      return 'Inflammation and immunity' if term_up =~ /IMMUNE|INFLAMMATION/
      return 'Invasion, metastasis and angiogenesis' if term_up =~ /INVASION|METASTASIS|ANGIOGENESIS/
      return 'Other cancer hallmarks'
    end

    return 'Cell cycle and genome' if term_down =~ /cell cycle|mitotic|chromosome|spindle|cytokinesis|replication|dna repair|dna recombination|checkpoint/
    return 'RNA, ribosome and protein homeostasis' if term_down =~ /ribosome|rrna|trna|mrna|rna |transcription|translation|splicing|protein folding|proteasome|unfolded/
    return 'Metabolism and mitochondria' if term_down =~ /metabolic|metabolism|glycolysis|lipid|amino acid|mitochond|respiration|energy|oxidative phosphorylation|mtorc/
    return 'Stress, death and autophagy' if term_down =~ /stress|apopt|death|autophagy|hypoxia|p53|uv response|unfolded/
    return 'Inflammation and immunity' if term_down =~ /immune|inflammatory|cytokine|interferon|tnf|nf.kb|complement|il6|jak|stat/
    return 'Signaling' if term_down =~ /signaling|signal|phosphorylation|kinase|mapk|pi3k|wnt|tgf|notch|hedgehog|kras/
    return 'Adhesion, ECM and plasticity' if term_down =~ /adhesion|migration|motility|extracellular matrix|emt|epithelial|mesenchymal|angiogenesis|coagulation|apical/
    return 'Development and differentiation' if term_down =~ /differentiation|development|morphogenesis|anatomical|nervous system/
    'Other processes'
  end

  def self.figure_functional_contexts(column_grouping)
    contexts = []
    gap = 0.8
    if column_grouping.to_s == 'time_treatment'
      FIGURE_TIME_POINTS.each_with_index do |time, outer_i|
        group_start = nil
        group_end = nil
        FIGURE_TREATMENT_ORDER.each_with_index do |treatment, inner_i|
          x = outer_i * (FIGURE_TREATMENT_ORDER.length + gap) + inner_i + 1
          group_start ||= x
          group_end = x
          contexts << {
            :treatment => treatment,
            :time => time,
            :x => x,
            :inner_label => figure_treatment_label(treatment),
            :inner_color => figure_treatment_label(treatment),
            :group_label => "#{time}h",
            :group_color => 'black',
            :group_key => time.to_s
          }
        end
      end
    else
      FIGURE_TREATMENT_ORDER.each_with_index do |treatment, outer_i|
        FIGURE_TIME_POINTS.each_with_index do |time, inner_i|
          x = outer_i * (FIGURE_TIME_POINTS.length + gap) + inner_i + 1
          contexts << {
            :treatment => treatment,
            :time => time,
            :x => x,
            :inner_label => "#{time}h",
            :inner_color => 'black',
            :group_label => figure_treatment_label(treatment),
            :group_color => figure_treatment_label(treatment),
            :group_key => treatment
          }
        end
      end
    end
    contexts
  end

  def self.figure_functional_fill_class(direction, qvalue, qvalue_strong, qvalue_very_strong)
    direction = direction.to_s
    return 'none' if direction.empty? || direction == 'none'
    qvalue = qvalue.to_f
    suffix = if qvalue <= qvalue_very_strong.to_f
               'q1e6'
             elsif qvalue <= qvalue_strong.to_f
               'q1e4'
             else
               'q05'
             end
    direction + '_' + suffix
  end

  def self.figure_cap(value, cap)
    [[value.to_f, cap.to_f].min, -cap.to_f].max
  end

  helper :figure do |tsv, code, w, h|
  end

  extension :svg
  task :figure_overall_intermediate_phenotypes => :text do
    tsv = AGS.figure_build_tsv([['one', 1]], 'ID', ['Value'])
    figure_ggplot(tsv, <<-RCODE, 9, 4)
ggplot(data.frame(x=1,y=1), aes(x,y)) +
  xlim(0, 10) + ylim(0, 5) +
  annotate('rect', xmin=0.4, xmax=2.4, ymin=2.1, ymax=3.5, fill='#F0F0F0', color='grey30') +
  annotate('rect', xmin=3.0, xmax=5.1, ymin=2.1, ymax=3.5, fill='#E6F2FF', color='grey30') +
  annotate('rect', xmin=5.7, xmax=7.8, ymin=2.1, ymax=3.5, fill='#FFF2CC', color='grey30') +
  annotate('rect', xmin=8.4, xmax=9.8, ymin=2.1, ymax=3.5, fill='#FCE4D6', color='grey30') +
  annotate('text', x=1.4, y=2.8, label='Kinase\ninhibition', size=4) +
  annotate('text', x=4.05, y=2.8, label='TF activity\nintermediate\nphenotype', size=4) +
  annotate('text', x=6.75, y=2.8, label='Gene expression\nintermediate\nphenotype', size=4) +
  annotate('text', x=9.1, y=2.8, label='Cellular\nphenotype', size=4) +
  annotate('segment', x=2.45, xend=2.95, y=2.8, yend=2.8, arrow=arrow(length=unit(0.15,'inches'))) +
  annotate('segment', x=5.15, xend=5.65, y=2.8, yend=2.8, arrow=arrow(length=unit(0.15,'inches'))) +
  annotate('segment', x=7.85, xend=8.35, y=2.8, yend=2.8, arrow=arrow(length=unit(0.15,'inches'))) +
  annotate('text', x=5.1, y=1.3, label='Molecular causalities are inferred from time-resolved intermediate phenotypes', size=4) +
  theme_void()
RCODE
  end

  extension :svg
  task :figure_overall_dynamic_workflow => :text do
    tsv = AGS.figure_build_tsv([['one', 1]], 'ID', ['Value'])
    figure_ggplot(tsv, <<-RCODE, 10, 4.5)
ggplot(data.frame(x=1,y=1), aes(x,y)) +
  xlim(0, 10) + ylim(0, 5) +
  annotate('rect', xmin=0.3, xmax=1.8, ymin=2.4, ymax=3.5, fill='#F0F0F0', color='grey30') +
  annotate('rect', xmin=2.3, xmax=4.1, ymin=2.4, ymax=3.5, fill='#E2F0D9', color='grey30') +
  annotate('rect', xmin=4.6, xmax=6.3, ymin=2.4, ymax=3.5, fill='#DDEBF7', color='grey30') +
  annotate('rect', xmin=6.8, xmax=8.2, ymin=2.4, ymax=3.5, fill='#FFF2CC', color='grey30') +
  annotate('rect', xmin=8.6, xmax=9.8, ymin=2.4, ymax=3.5, fill='#FCE4D6', color='grey30') +
  annotate('text', x=1.05, y=2.95, label='RNA-seq\ntime series', size=3.7) +
  annotate('text', x=3.2, y=2.95, label='Onset-defined\ndynamic genes', size=3.7) +
  annotate('text', x=5.45, y=2.95, label='CollecTRI2\nregulome', size=3.7) +
  annotate('text', x=7.5, y=2.95, label='TF activity\nscoreboards', size=3.7) +
  annotate('text', x=9.2, y=2.95, label='Process\nscoreboards', size=3.7) +
  annotate('segment', x=1.85, xend=2.25, y=2.95, yend=2.95, arrow=arrow(length=unit(0.13,'inches'))) +
  annotate('segment', x=4.15, xend=4.55, y=2.95, yend=2.95, arrow=arrow(length=unit(0.13,'inches'))) +
  annotate('segment', x=6.35, xend=6.75, y=2.95, yend=2.95, arrow=arrow(length=unit(0.13,'inches'))) +
  annotate('segment', x=8.25, xend=8.55, y=2.95, yend=2.95, arrow=arrow(length=unit(0.13,'inches'))) +
  annotate('text', x=5.0, y=1.4, label='Dynamic inference asks which regulators explain genes entering a new expression regime', size=3.6) +
  theme_void()
RCODE
  end

  dep :fold_changes, :fc_source => 'NTNU'
  extension :svg
  task :figure_mrna_pca_by_treatment => :text do
    fc0 = step(:fold_changes).load.transpose('Associated Gene Name')
    fields = AGS.figure_context_fields('FC').select{|field| fc0.fields.include?(field) }
    data_tsv = fc0.reorder('Associated Gene Name', fields)
    figure_ggplot(data_tsv, <<-RCODE, 7, 5.5)
mat <- data[, c(#{fields.collect{|f| "'#{f}'"} * ', '}), drop=FALSE]
mat[] <- lapply(mat, as.numeric)
mat[is.na(mat)] <- 0
mat <- t(as.matrix(mat))
keep <- apply(mat, 2, sd) > 0
mat <- mat[, keep, drop=FALSE]
pca <- prcomp(mat, center=TRUE, scale.=TRUE)
contexts <- rownames(mat)
plot_df <- data.frame(Context=contexts, PC1=pca$x[,1], PC2=pca$x[,2], stringsAsFactors=FALSE)
tmp_context <- sub('^FC_', '', plot_df$Context)
parts <- strsplit(tmp_context, '.T', fixed=TRUE)
plot_df$Treatment <- vapply(parts, function(x) x[1], character(1))
plot_df$Time <- as.numeric(vapply(parts, function(x) x[2], character(1)))
plot_df$TreatmentLabel <- factor(plot_df$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
var_exp <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
#{AGS.figure_treatment_color_scale('color')}
ggplot(plot_df, aes(PC1, PC2, color=TreatmentLabel)) + geom_point(aes(size=Time), alpha=0.9) +
  scale_size_continuous(breaks=c(1,2,4,8,24), name='Time') +
  theme_bw() + labs(x=paste0('PC1 (', var_exp[1], '%)'), y=paste0('PC2 (', var_exp[2], '%)'), color='Treatment') +
  theme(legend.position='bottom')
RCODE
  end

  dep :fold_changes, :fc_source => 'NTNU'
  extension :svg
  task :figure_mrna_pca_by_time => :text do
    fc0 = step(:fold_changes).load.transpose('Associated Gene Name')
    fields = AGS.figure_context_fields('FC').select{|field| fc0.fields.include?(field) }
    data_tsv = fc0.reorder('Associated Gene Name', fields)
    figure_ggplot(data_tsv, <<-RCODE, 7, 5.5)
mat <- data[, c(#{fields.collect{|f| "'#{f}'"} * ', '}), drop=FALSE]
mat[] <- lapply(mat, as.numeric)
mat[is.na(mat)] <- 0
mat <- t(as.matrix(mat))
keep <- apply(mat, 2, sd) > 0
mat <- mat[, keep, drop=FALSE]
pca <- prcomp(mat, center=TRUE, scale.=TRUE)
contexts <- rownames(mat)
plot_df <- data.frame(Context=contexts, PC1=pca$x[,1], PC2=pca$x[,2], stringsAsFactors=FALSE)
tmp_context <- sub('^FC_', '', plot_df$Context)
parts <- strsplit(tmp_context, '.T', fixed=TRUE)
plot_df$Treatment <- vapply(parts, function(x) x[1], character(1))
plot_df$Time <- factor(as.numeric(vapply(parts, function(x) x[2], character(1))), levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
plot_df$TreatmentLabel <- factor(plot_df$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
var_exp <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
ggplot(plot_df, aes(PC1, PC2, color=Time, shape=TreatmentLabel)) + geom_point(size=3, alpha=0.9) +
  scale_color_brewer(palette='YlOrRd') + theme_bw() +
  labs(x=paste0('PC1 (', var_exp[1], '%)'), y=paste0('PC2 (', var_exp[2], '%)'), shape='Treatment') +
  theme(legend.position='bottom')
RCODE
  end

  dep :fold_changes, :fc_source => 'NTNU'
  extension :svg
  task :figure_mrna_context_distance_clustered => :text do
    fc0 = step(:fold_changes).load.transpose('Associated Gene Name')
    fields = AGS.figure_context_fields('FC').select{|field| fc0.fields.include?(field) }
    data_tsv = fc0.reorder('Associated Gene Name', fields)
    figure_base(data_tsv, <<-RCODE, 11, 10)
rbbt.require('gplots')
mat <- data[, c(#{fields.collect{|f| "'#{f}'"} * ', '}), drop=FALSE]
mat[] <- lapply(mat, as.numeric)
mat[is.na(mat)] <- 0
mat <- as.matrix(mat)
cor_mat <- cor(mat, use='pairwise.complete.obs')
dist_mat <- 1 - cor_mat
labs <- colnames(dist_mat)
labs <- sub('^FC_', '', labs)
labs <- gsub('INT_FiveZ_PI', '5Z+PI', labs)
labs <- gsub('INT_PD_PI', 'PD+PI', labs)
labs <- gsub('FiveZ', '5Z', labs)
#labs <- gsub('\\.T', ' T', labs)
rownames(dist_mat) <- labs
colnames(dist_mat) <- labs
heatmap.2(dist_mat, trace='none', col=colorRampPalette(c('white','#2166AC'))(101),
          dendrogram='both', density.info='none', key.title='1-r', margins=c(9,9), symm=TRUE)
RCODE
  end

  dep :de_gene_counts_fc0
  extension :svg
  task :figure_de_counts_fc0_bar => :text do
    counts = step(:de_gene_counts_fc0).load
    figure_ggplot(counts, <<-RCODE, 10, 5)
data <- subset(data, Direction == 'both')
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Genes <- as.numeric(data$Genes)
#{AGS.figure_treatment_color_scale('fill')}
ggplot(data, aes(Time, Genes, fill=TreatmentLabel)) + geom_col(position='dodge') +
  theme_bw() + labs(x='Time point', y='DE genes relative to baseline', fill='Treatment') + theme(legend.position='bottom')
RCODE
  end

  dep :interval_de_gene_counts_fc1
  extension :svg
  task :figure_de_counts_fc1_bar => :text do
    counts = step(:interval_de_gene_counts_fc1).load
    figure_ggplot(counts, <<-RCODE, 10, 5)
data <- subset(data, Direction == 'both')
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Genes <- as.numeric(data$Genes)
#{AGS.figure_treatment_color_scale('fill')}
ggplot(data, aes(Time, Genes, fill=TreatmentLabel)) + geom_col(position='dodge') +
  theme_bw() + labs(x='Interval ending at time point', y='Interval DE genes', fill='Treatment') + theme(legend.position='bottom')
RCODE
  end

  dep :onset_first_counts
  extension :svg
  task :figure_dynamic_onset_first_bar => :text do
    tsv = step(:onset_first_counts).load
    figure_ggplot(tsv, <<-RCODE, 10, 5)
data <- subset(data, Direction != 'both')
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Direction <- factor(data$Direction, levels=c('up','down'))
data$Genes <- as.numeric(data$Genes)
#{AGS.figure_direction_color_scale('fill')}
ggplot(data, aes(Time, Genes, fill=Direction)) + geom_col() + facet_wrap(~TreatmentLabel, nrow=2) +
  theme_bw() + labs(x='First onset time', y='Genes', fill='Direction') + theme(legend.position='bottom')
RCODE
  end

  dep :onset_episode_counts
  extension :svg
  task :figure_dynamic_onset_episode_bar => :text do
    tsv = step(:onset_episode_counts).load
    figure_ggplot(tsv, <<-RCODE, 10, 5)
data <- subset(data, Direction != 'both')
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Direction <- factor(data$Direction, levels=c('up','down'))
data$Episodes <- as.numeric(data$Episodes)
#{AGS.figure_direction_color_scale('fill')}
ggplot(data, aes(Time, Episodes, fill=Direction)) + geom_col() + facet_wrap(~TreatmentLabel, nrow=2) +
  theme_bw() + labs(x='Onset episode time', y='Episodes', fill='Direction') + theme(legend.position='bottom')
RCODE
  end

  dep :onset_direction_switch_summary
  extension :svg
  task :figure_dynamic_onset_switch_summary => :text do
    tsv = step(:onset_direction_switch_summary).load
    figure_ggplot(tsv, <<-RCODE, 8, 5)
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Category <- factor(data$Category, levels=c('unclassified','single_episode','multiple_same_direction','direction_switch'))
data$Genes <- as.numeric(data$Genes)
ggplot(data, aes(TreatmentLabel, Genes, fill=Category)) + geom_col() +
  scale_fill_brewer(palette='Set2') + theme_bw() + labs(x='', y='Genes', fill='Expression-regime class') +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
RCODE
  end

  dep :fc1_onset_relationship_summary
  extension :svg
  task :figure_dynamic_vs_interval_onset_relationship => :text do
    tsv = step(:fc1_onset_relationship_summary).load
    figure_ggplot(tsv, <<-RCODE, 12, 6)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Genes <- as.numeric(data$Genes)
data$Relationship <- factor(data$Relationship)
ggplot(data, aes(Time, Genes, fill=Relationship)) + geom_col() +
  facet_grid(Direction ~ TreatmentLabel, scales='free_y') + scale_fill_brewer(palette='Paired') +
  theme_bw() + labs(x='Time point', y='Genes', fill='FC1 versus onset') +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom', legend.text=element_text(size=7))
RCODE
  end

  input :query_type, :select, 'Query gene-set source', 'onset_genes', :select_options => %w(onset_genes tfs targets tf_and_targets)
  input :annotation, :select, 'Annotation source', 'cancerhallmarks_core', :select_options => %w(goslim msigdb_hallmark cancerhallmarks_core cancerhallmarks_integrated)
  input :column_grouping, :select, 'Column grouping', 'time_treatment', :select_options => %w(time_treatment treatment_time)
  input :max_terms, :integer, 'Maximum terms to show', 30
  input :adjusted_pvalue_threshold, :float, 'FDR cutoff for display', 0.05
  input :qvalue_strong, :float, 'Darker shade cutoff', 1e-4
  input :qvalue_very_strong, :float, 'Darkest shade cutoff', 1e-6
  input :reduced_terms, :boolean, 'Use redundancy-reduced enrichment table', false
  dep :functional_enrichment_reduced do |jobname, options|
    options.merge(:query_type => options[:query_type], :annotation => options[:annotation], :adjusted_pvalue_threshold => options[:adjusted_pvalue_threshold])
  end
  dep :functional_enrichment do |jobname, options|
    options.merge(:query_type => options[:query_type], :annotation => options[:annotation])
  end
  extension :svg
  task :figure_functional_enrichment_scoreboard => :text do |query_type, annotation, column_grouping, max_terms, adjusted_pvalue_threshold, qvalue_strong, qvalue_very_strong, reduced_terms|
    table = if reduced_terms
              dependencies.find{|dep| dep.task_name == :functional_enrichment_reduced }.load
            else
              dependencies.find{|dep| dep.task_name == :functional_enrichment }.load
            end
    idx = Hash[table.fields.each_with_index.to_a]
    allowed_directions = query_type.to_s == 'onset_genes' ? %w(up down) : %w(positive negative)
    best = {}
    term_best = Hash.new(1.0)
    table.through do |row_id, values|
      values = NamedArray.setup(values, table.fields)
      direction = values['Direction'].to_s
      next unless allowed_directions.include?(direction)
      qvalue = AGS.figure_safe_float(values['AdjustedPValue'], 1.0)
      next if qvalue > adjusted_pvalue_threshold.to_f
      treatment = values['Treatment'].to_s
      time = AGS.figure_safe_int(values['Time'])
      term = values['TermName'].to_s
      next if term.empty?
      key = [term, treatment, time]
      previous = best[key]
      if previous.nil? || qvalue < previous[:qvalue]
        best[key] = {
          :direction => direction,
          :qvalue => qvalue,
          :pvalue => AGS.figure_safe_float(values['PValue'], 1.0),
          :intersection => AGS.figure_safe_int(values['IntersectionSize'], 0),
          :query_size => AGS.figure_safe_int(values['QuerySize'], 0),
          :term_size => AGS.figure_safe_int(values['TermBackgroundSize'], 0),
          :term_id => values['TermID'].to_s,
          :group => AGS.figure_functional_term_group(term, annotation)
        }
      end
      term_best[term] = [term_best[term], qvalue].min
    end

    selected_terms = term_best.keys.sort_by{|term| [term_best[term], AGS.figure_functional_term_group(term, annotation), term] }.first(max_terms.to_i)
    selected_terms = selected_terms.sort_by{|term| [AGS.figure_functional_term_group(term, annotation), term_best[term], term] }
    contexts = AGS.figure_functional_contexts(column_grouping)
    group_order = selected_terms.collect{|term| AGS.figure_functional_term_group(term, annotation) }.uniq
    term_index = {}
    selected_terms.each_with_index{|term, i| term_index[term] = i + 1 }

    rows = []
    id = 0
    selected_terms.each do |term|
      contexts.each do |context|
        entry = best[[term, context[:treatment], context[:time]]]
        direction = entry ? entry[:direction] : 'none'
        qvalue = entry ? entry[:qvalue] : 1.0
        fill_class = entry ? AGS.figure_functional_fill_class(direction, qvalue, qvalue_strong, qvalue_very_strong) : 'none'
        id += 1
        rows << [id, term, term_index[term], AGS.figure_functional_term_group(term, annotation), context[:treatment], AGS.figure_treatment_label(context[:treatment]), context[:time], "#{context[:time]}h", context[:x], context[:inner_label], context[:inner_color], context[:group_label], context[:group_color], context[:group_key], direction, fill_class, qvalue, entry ? entry[:intersection] : 0, entry ? entry[:query_size] : 0, entry ? entry[:term_size] : 0]
      end
    end
    plot_tsv = AGS.figure_build_tsv(rows, 'ID', %w(Term TermIndex TermGroup Treatment TreatmentLabel Time TimeLabel X InnerLabel InnerColor GroupLabel GroupColor GroupKey Direction FillClass QValue Intersection QuerySize TermSize))
    fill_scale = AGS.figure_enrichment_fill_scale
    figure_ggplot(plot_tsv, <<-RCODE, 15, [5, selected_terms.length * 0.28 + 2.8].max)
data$TermIndex <- as.numeric(data$TermIndex)
data$X <- as.numeric(data$X)
data$QValue <- as.numeric(data$QValue)
data$Intersection <- as.numeric(data$Intersection)
data$Term <- factor(data$Term, levels=rev(unique(data$Term[order(data$TermIndex)])))
label_df <- unique(data[, c('X','InnerLabel','InnerColor','GroupLabel','GroupColor','GroupKey')])
group_df <- aggregate(X ~ GroupKey + GroupLabel + GroupColor, data=label_df, FUN=function(x) mean(range(x)))
sep_df <- aggregate(X ~ GroupKey, data=label_df, FUN=max)
term_group_df <- unique(data[, c('TermGroup','TermIndex')])
term_group_mid <- aggregate(TermIndex ~ TermGroup, data=term_group_df, FUN=function(x) mean(range(x)))
term_group_sep <- aggregate(TermIndex ~ TermGroup, data=term_group_df, FUN=max)
max_y <- max(data$TermIndex)
max_x <- max(data$X)
ggplot(data, aes(X, TermIndex)) +
  geom_tile(aes(fill=FillClass), color='grey88', width=0.92, height=0.92) +
  geom_vline(data=sep_df, aes(xintercept=X + 0.55), inherit.aes=FALSE, color='grey70', linewidth=0.35) +
  geom_hline(data=term_group_sep, aes(yintercept=TermIndex + 0.5), inherit.aes=FALSE, color='grey80', linewidth=0.3) +
  geom_text(data=label_df, aes(x=X, y=-0.20, label=InnerLabel, color=InnerColor), inherit.aes=FALSE, angle=90, hjust=1, size=2.6) +
  geom_text(data=group_df, aes(x=X, y=-1.05, label=GroupLabel, color=GroupColor), inherit.aes=FALSE, fontface='bold', size=3.0) +
  geom_text(data=term_group_mid, aes(x=max_x + 0.8, y=TermIndex, label=TermGroup), inherit.aes=FALSE, hjust=0, size=2.7, color='grey30') +
  #{fill_scale} + scale_color_manual(values=c('black'='black','DMSO'='#7F7F7F','5Z'='#4DAF4A','5Z+PI'='#FF7F00','PD+PI'='#E41A1C','PD'='#984EA3','PI'='#377EB8'), guide='none') +
  scale_y_continuous(breaks=data$TermIndex[match(levels(data$Term), data$Term)], labels=levels(data$Term), limits=c(-1.45, max_y + 0.5), expand=c(0,0)) +
  scale_x_continuous(limits=c(min(data$X)-0.6, max_x + 4.8), expand=c(0,0)) +
  coord_cartesian(clip='off') + theme_bw() + labs(x='', y='', title='#{query_type} / #{annotation}') +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=7), panel.grid=element_blank(), plot.margin=margin(8, 130, 45, 8), legend.position='bottom')
RCODE
  end

  input :query_type, :select, 'Query gene-set source', 'onset_genes', :select_options => %w(onset_genes tfs targets tf_and_targets)
  input :annotation, :select, 'Annotation source', 'cancerhallmarks_core', :select_options => %w(goslim msigdb_hallmark cancerhallmarks_core cancerhallmarks_integrated)
  dep :functional_enrichment_reduced do |jobname, options|
    options.merge(:query_type => options[:query_type], :annotation => options[:annotation])
  end
  extension :svg
  task :figure_functional_enrichment_term_frequency => :text do |query_type, annotation|
    table = step(:functional_enrichment_reduced).load
    idx = Hash[table.fields.each_with_index.to_a]
    counts = Hash.new{|h,k| h[k] = [0, 1.0, ''] }
    table.through do |row_id, values|
      values = NamedArray.setup(values, table.fields)
      direction = values['Direction'].to_s
      next if direction == 'both'
      term = values['TermName'].to_s
      q = AGS.figure_safe_float(values['AdjustedPValue'], 1.0)
      counts[term][0] += 1
      counts[term][1] = [counts[term][1], q].min
      counts[term][2] = AGS.figure_functional_term_group(term, annotation)
    end
    rows = counts.keys.sort_by{|term| [-counts[term][0], counts[term][1], term] }.first(30).each_with_index.collect do |term, i|
      [i + 1, term, counts[term][2], counts[term][0], -Math.log10([counts[term][1], 1e-300].max)]
    end
    plot_tsv = AGS.figure_build_tsv(rows, 'ID', %w(Term TermGroup Contexts BestNegLogQ))
    figure_ggplot(plot_tsv, <<-RCODE, 8, 8)
data$Contexts <- as.numeric(data$Contexts)
data$BestNegLogQ <- as.numeric(data$BestNegLogQ)
data$Term <- factor(data$Term, levels=rev(data$Term[order(data$Contexts, data$BestNegLogQ)]))
ggplot(data, aes(Term, Contexts, fill=BestNegLogQ)) + geom_col() + coord_flip() +
  facet_grid(TermGroup ~ ., scales='free_y', space='free_y') + scale_fill_viridis_c(name='best -log10(FDR)') +
  theme_bw() + labs(x='', y='Enriched treatment-time-direction contexts')
RCODE
  end

  dep :tf_activity_call_counts_by_scheme, :scheme => 'dynamic'
  extension :svg
  task :figure_tf_activity_counts_bar => :text do
    tsv = step(:tf_activity_call_counts_by_scheme).load
    figure_ggplot(tsv, <<-RCODE, 10, 5)
data <- subset(data, Sign != 'both')
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$TFActivityCalls <- as.numeric(data$TFActivityCalls)
data$SignedCalls <- ifelse(data$Sign == 'negative', -data$TFActivityCalls, data$TFActivityCalls)
ggplot(data, aes(TreatmentLabel, SignedCalls, fill=Sign)) + geom_col() + facet_wrap(~Time, nrow=1, scales='free_y') +
  scale_fill_manual(values=c('positive'='#D73027','negative'='#4575B4')) + coord_flip() +
  theme_bw() + labs(x='', y='TF activity-change calls', fill='Activity sign') + theme(legend.position='bottom')
RCODE
  end

  dep :tf_activity_heatmap_matrix, :scheme => 'dynamic', :normalization => 'raw'
  extension :svg
  task :figure_tf_activity_pca_by_treatment => :text do
    matrix = step(:tf_activity_heatmap_matrix).load
    fields = AGS.figure_tf_context_fields(matrix)
    data_tsv = matrix.reorder('Associated Gene Name', fields)
    figure_ggplot(data_tsv, <<-RCODE, 7, 5.5)
mat <- data[, c(#{fields.collect{|f| "'#{f}'"} * ', '}), drop=FALSE]
mat[] <- lapply(mat, as.numeric)
mat[is.na(mat)] <- 0
mat <- t(as.matrix(mat))
keep <- apply(mat, 2, sd) > 0
mat <- mat[, keep, drop=FALSE]
pca <- prcomp(mat, center=TRUE, scale.=TRUE)
plot_df <- data.frame(Context=rownames(mat), PC1=pca$x[,1], PC2=pca$x[,2], stringsAsFactors=FALSE)
plot_df$Treatment <- sub('-T[0-9]+$', '', plot_df$Context)
plot_df$Time <- as.numeric(sub('.*-T', '', plot_df$Context))
plot_df$TreatmentLabel <- factor(plot_df$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
var_exp <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
#{AGS.figure_treatment_color_scale('color')}
ggplot(plot_df, aes(PC1, PC2, color=TreatmentLabel)) + geom_point(aes(size=Time), alpha=0.9) +
  scale_size_continuous(breaks=c(1,2,4,8,24), name='Time') + theme_bw() +
  labs(x=paste0('PC1 (', var_exp[1], '%)'), y=paste0('PC2 (', var_exp[2], '%)'), color='Treatment') + theme(legend.position='bottom')
RCODE
  end

  dep :tf_activity_call_counts_by_scheme, :scheme => 'dynamic'
  extension :svg
  task :figure_tf_activity_sign_balance_bar => :text do
    tsv = step(:tf_activity_call_counts_by_scheme).load
    figure_ggplot(tsv, <<-RCODE, 9, 5)
data <- subset(data, Sign != 'both')
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$TFActivityCalls <- as.numeric(data$TFActivityCalls)
ggplot(data, aes(Time, TFActivityCalls, fill=Sign)) + geom_col(position='fill') + facet_wrap(~TreatmentLabel, nrow=2) +
  scale_fill_manual(values=c('positive'='#D73027','negative'='#4575B4')) + theme_bw() +
  labs(x='Time point', y='Fraction of TF activity calls', fill='Sign') + theme(legend.position='bottom')
RCODE
  end

  input :scheme, :select, 'TF prediction scheme', 'dynamic', :select_options => %w(dynamic non-dynamic)
  input :normalization, :select, 'Normalization', 'row_zscore', :select_options => %w(raw row_zscore column_zscore)
  input :max_abs, :float, 'Maximum absolute value for color scale', 3.0
  input :cluster_columns, :boolean, 'Cluster columns and show column dendrogram', false
  input :show_tf_labels, :boolean, 'Show TF labels', true
  dep :tf_activity_heatmap_matrix, :scheme => :placeholder, :normalization => :placeholder do |jobname, options|
    { :scheme => options[:scheme], :normalization => options[:normalization] }
  end
  extension :svg
  task :figure_tf_activity_heatmap_dendrogram => :text do |scheme, normalization, max_abs, cluster_columns, show_tf_labels|
    matrix = step(:tf_activity_heatmap_matrix).load
    fields = AGS.figure_tf_context_fields(matrix)
    data_tsv = matrix.reorder('Associated Gene Name', fields)
    dendrogram = cluster_columns ? 'both' : 'row'
    colv = cluster_columns ? 'TRUE' : 'FALSE'
    labrow = show_tf_labels ? 'rownames(mat)' : 'rep("", nrow(mat))'
    figure_base(data_tsv, <<-RCODE, 12, 10)
rbbt.require('gplots')
mat <- data[, c(#{fields.collect{|f| "'#{f}'"} * ', '}), drop=FALSE]
mat[] <- lapply(mat, as.numeric)
mat[is.na(mat)] <- 0
mat <- as.matrix(mat)
rownames(mat) <- rownames(data)
mat[mat > #{max_abs}] <- #{max_abs}
mat[mat < -#{max_abs}] <- -#{max_abs}
colnames(mat) <- gsub('INT_FiveZ_PI', '5Z+PI', colnames(mat))
colnames(mat) <- gsub('INT_PD_PI', 'PD+PI', colnames(mat))
colnames(mat) <- gsub('FiveZ', '5Z', colnames(mat))
heatmap.2(mat, trace='none', col=colorRampPalette(c('#4575B4','white','#D73027'))(101),
          dendrogram='#{dendrogram}', Colv=#{colv}, density.info='none', key.title='activity',
          symbreaks=FALSE, labRow=#{labrow}, margins=c(9,4), cexCol=0.65, cexRow=0.35)
RCODE
  end

  dep :tf_predictions
  input :combination, :select, 'Combination treatment', 'INT_FiveZ_PI', :select_options => %w(INT_FiveZ_PI INT_PD_PI)
  input :time_point, :integer, 'Time point', 2
  input :sign_mode, :select, 'Activity sign subset', 'all', :select_options => %w(all positive negative)
  extension :svg
  task :figure_combination_tf_activity_upset_bar => :text do |combination, time_point, sign_mode|
    predictions = step(:tf_predictions).load
    components = combination == 'INT_PD_PI' ? ['PI', 'PD'] : ['PI', 'FiveZ']
    labels = components + [combination]
    sets = Hash.new{|h,k| h[k] = [] }
    predictions.through do |tf, values|
      values = NamedArray.setup(values, predictions.fields)
      labels.each do |treatment|
        field = "#{treatment}-T#{time_point}"
        next unless predictions.fields.include?(field)
        value = AGS.figure_safe_float(values[field], 0.0)
        active = case sign_mode
                 when 'positive' then value > 0
                 when 'negative' then value < 0
                 else value != 0
                 end
        sets[treatment] << tf if active
      end
    end
    counts = Hash.new(0)
    all_tfs = labels.collect{|l| sets[l] }.flatten.uniq
    all_tfs.each do |tf|
      pattern = labels.collect{|l| sets[l].include?(tf) ? '1' : '0' } * ''
      counts[pattern] += 1
    end
    rows = counts.keys.sort_by{|pattern| [-counts[pattern], pattern] }.each_with_index.collect do |pattern, i|
      display = labels.each_with_index.collect{|l,j| pattern[j] == '1' ? AGS.figure_treatment_label(l) : nil }.compact * ' & '
      display = 'none' if display.empty?
      [i + 1, pattern, display, counts[pattern]]
    end
    plot_tsv = AGS.figure_build_tsv(rows, 'ID', %w(Pattern Display Count))
    title = "#{AGS.figure_treatment_label(combination)} T#{time_point} #{sign_mode} TF activity exact intersections"
    figure_ggplot(plot_tsv, <<-RCODE, 8, 5)
data$Count <- as.numeric(data$Count)
data <- data[order(data$Count, decreasing=TRUE),]
data$Display <- factor(data$Display, levels=rev(data$Display))
ggplot(data, aes(Display, Count)) + geom_col(fill='grey35') + coord_flip() +
  theme_bw() + labs(x='Exact active-TF intersection', y='TFs', title='#{title}')
RCODE
  end

  dep :combination_tf_category_counts
  extension :svg
  task :figure_combination_category_counts_bar => :text do
    tsv = step(:combination_tf_category_counts).load
    figure_ggplot(tsv, <<-RCODE, 10, 5)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$CombinationLabel <- factor(data$Combination, levels=c('INT_FiveZ_PI','INT_PD_PI'), labels=c('5Z+PI','PD+PI'))
data$TFActivityCalls <- as.numeric(data$TFActivityCalls)
data$CategoryLabel <- factor(data$Category, levels=c(#{FIGURE_COMBINATION_CATEGORY_LABELS.keys.collect{|c| "'#{c}'"} * ','}), labels=c(#{FIGURE_COMBINATION_CATEGORY_LABELS.values.collect{|c| "'#{c}'"} * ','}))
ggplot(data, aes(Time, TFActivityCalls, fill=CategoryLabel)) + geom_col() + facet_wrap(~CombinationLabel, nrow=1) +
  scale_fill_brewer(palette='Set3') + theme_bw() + labs(x='Time point', y='TF activity calls in combination', fill='Category') +
  theme(legend.position='bottom', legend.text=element_text(size=7))
RCODE
  end

  dep :tf_predictions
  input :time_point, :integer, 'Early time point', 1
  extension :svg
  task :figure_opposite_early_single_agent_tf_signs => :text do |time_point|
    predictions = step(:tf_predictions).load
    pairs = [['PI','PD'], ['PI','FiveZ'], ['PD','FiveZ']]
    counts = Hash.new(0)
    predictions.through do |tf, values|
      values = NamedArray.setup(values, predictions.fields)
      pairs.each do |a,b|
        fa = "#{a}-T#{time_point}"; fb = "#{b}-T#{time_point}"
        next unless predictions.fields.include?(fa) && predictions.fields.include?(fb)
        va = AGS.figure_safe_float(values[fa], 0.0); vb = AGS.figure_safe_float(values[fb], 0.0)
        next if va == 0 && vb == 0
        category = if va != 0 && vb != 0
                     va * vb > 0 ? 'shared same sign' : 'shared opposite sign'
                   elsif va != 0
                     "#{a} only"
                   else
                     "#{b} only"
                   end
        counts[["#{AGS.figure_treatment_label(a)} vs #{AGS.figure_treatment_label(b)}", category]] += 1
      end
    end
    rows = counts.each_with_index.collect{|((pair, category), count), i| [i + 1, pair, category, count] }
    plot_tsv = AGS.figure_build_tsv(rows, 'ID', %w(Pair Category Count))
    figure_ggplot(plot_tsv, <<-RCODE, 8, 5)
data$Count <- as.numeric(data$Count)
ggplot(data, aes(Pair, Count, fill=Category)) + geom_col(position='stack') + coord_flip() +
  scale_fill_brewer(palette='Set2') + theme_bw() + labs(x='', y='TFs', fill='Sign-aware overlap', title='Early single-agent TF signs at T#{time_point}') +
  theme(legend.position='bottom')
RCODE
  end

  dep :tf_predictions
  input :time_point, :integer, 'Time point', 1
  input :treatment_x, :select, 'Treatment on x axis', 'PI', :select_options => FIGURE_TREATMENT_ORDER
  input :treatment_y, :select, 'Treatment on y axis', 'FiveZ', :select_options => FIGURE_TREATMENT_ORDER
  input :label_tfs, :integer, 'Number of extreme TFs to label', 12
  extension :svg
  task :figure_tf_activity_pair_scatter => :text do |time_point, treatment_x, treatment_y, label_tfs|
    predictions = step(:tf_predictions).load
    fx = "#{treatment_x}-T#{time_point}"
    fy = "#{treatment_y}-T#{time_point}"
    rows = []
    predictions.through do |tf, values|
      values = NamedArray.setup(values, predictions.fields)
      vx = predictions.fields.include?(fx) ? AGS.figure_safe_float(values[fx], 0.0) : 0.0
      vy = predictions.fields.include?(fy) ? AGS.figure_safe_float(values[fy], 0.0) : 0.0
      next if vx == 0 && vy == 0
      category = if vx != 0 && vy != 0
                   vx * vy > 0 ? 'same sign' : 'opposite sign'
                 elsif vx != 0
                   "#{treatment_x} only"
                 else
                   "#{treatment_y} only"
                 end
      rows << [tf, vx, vy, vx.abs + vy.abs, category]
    end
    top = rows.sort_by{|row| -row[3] }.first(label_tfs.to_i).collect{|row| row[0] }
    rows = rows.each_with_index.collect{|row,i| [i + 1] + row + [top.include?(row[0]) ? row[0] : ''] }
    plot_tsv = AGS.figure_build_tsv(rows, 'ID', %w(TF ActivityX ActivityY Magnitude Category Label))
    figure_ggplot(plot_tsv, <<-RCODE, 6, 6)
data$ActivityX <- as.numeric(data$ActivityX)
data$ActivityY <- as.numeric(data$ActivityY)
ggplot(data, aes(ActivityX, ActivityY, color=Category)) + geom_hline(yintercept=0, color='grey75') + geom_vline(xintercept=0, color='grey75') +
  geom_point(alpha=0.75, size=1.8) + geom_text(data=subset(data, Label != ''), aes(label=Label), size=2.5, vjust=-0.5, check_overlap=TRUE) +
  scale_color_brewer(palette='Set2') + theme_bw() + labs(x='#{AGS.figure_treatment_label(treatment_x)} T#{time_point}', y='#{AGS.figure_treatment_label(treatment_y)} T#{time_point}', color='Relationship')
RCODE
  end

  dep :neko_dynamic_non_dynamic_summary
  extension :svg
  task :figure_neko_match_fraction_bar => :text do
    tsv = step(:neko_dynamic_non_dynamic_summary).load
    figure_ggplot(tsv, <<-RCODE, 9, 5)
data <- subset(data, Vetting == 'none')
data$MatchFraction <- as.numeric(data$MatchFraction)
data$Target <- factor(data$Target, levels=c('T1','T2','relaxed'))
ggplot(data, aes(Target, MatchFraction, fill=Scheme)) + geom_col(position='dodge') + facet_wrap(~Treatment) +
  scale_fill_manual(values=c('dynamic'='#D73027','non-dynamic'='#4575B4')) + ylim(0,1) +
  theme_bw() + labs(x='Benchmark evaluation', y='Match fraction', fill='Scheme') + theme(legend.position='bottom')
RCODE
  end

  dep :neko_dynamic_non_dynamic_summary
  extension :svg
  task :figure_neko_counts_bar => :text do
    tsv = step(:neko_dynamic_non_dynamic_summary).load
    figure_ggplot(tsv, <<-RCODE, 9, 5)
data <- subset(data, Vetting == 'none')
data$Match <- as.numeric(data$Match)
data$Miss <- as.numeric(data$Miss)
data$Target <- factor(data$Target, levels=c('T1','T2','relaxed'))
long <- rbind(data.frame(data[, c('Treatment','Scheme','Target')], Class='Match', Count=data$Match), data.frame(data[, c('Treatment','Scheme','Target')], Class='Miss', Count=data$Miss))
ggplot(long, aes(Target, Count, fill=Class)) + geom_col() + facet_grid(Scheme ~ Treatment) +
  scale_fill_manual(values=c('Match'='#1A9850','Miss'='#D73027')) + theme_bw() + labs(x='Benchmark evaluation', y='TFs', fill='') +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
RCODE
  end

  dep :self_consistency_dynamic_non_dynamic_summary
  extension :svg
  task :figure_self_consistency_match_fraction_bar => :text do
    tsv = step(:self_consistency_dynamic_non_dynamic_summary).load
    figure_ggplot(tsv, <<-RCODE, 10, 5)
data <- subset(data, Vetting == 'none')
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$MatchFraction <- as.numeric(data$MatchFraction)
ggplot(data, aes(Time, MatchFraction, fill=Scheme)) + geom_col(position='dodge') + facet_wrap(~TreatmentLabel, nrow=2) +
  scale_fill_manual(values=c('dynamic'='#D73027','non-dynamic'='#4575B4')) + ylim(0,1) +
  theme_bw() + labs(x='Time point', y='Self-consistency match fraction', fill='Scheme') + theme(legend.position='bottom')
RCODE
  end

  dep :self_consistency_dynamic_non_dynamic_summary
  extension :svg
  task :figure_self_consistency_counts_bar => :text do
    tsv = step(:self_consistency_dynamic_non_dynamic_summary).load
    figure_ggplot(tsv, <<-RCODE, 10, 5)
data <- subset(data, Vetting == 'none')
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Matches <- as.numeric(data$Matches)
data$Miss <- as.numeric(data$Miss)
long <- rbind(data.frame(data[, c('TreatmentLabel','Time','Scheme')], Class='Match', Count=data$Matches), data.frame(data[, c('TreatmentLabel','Time','Scheme')], Class='Miss', Count=data$Miss))
ggplot(long, aes(Time, Count, fill=Class)) + geom_col() + facet_grid(Scheme ~ TreatmentLabel, scales='free_y') +
  scale_fill_manual(values=c('Match'='#1A9850','Miss'='#D73027')) + theme_bw() + labs(x='Time point', y='TF activity calls', fill='') +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
RCODE
  end

  dep :tf_predictions, :scheme => 'dynamic'
  dep :tf_predictions, :scheme => 'non-dynamic'
  extension :svg
  task :figure_dynamic_vs_nondynamic_earliest_detection => :text do
    dyn = dependencies.find{|d| d.recursive_inputs[:scheme].to_s == 'dynamic' }.load
    nd = dependencies.find{|d| d.recursive_inputs[:scheme].to_s == 'non-dynamic' }.load
    maps = {}
    {'dynamic' => dyn, 'non-dynamic' => nd}.each do |scheme_name, tsv|
      map = {}
      tsv.through do |tf, values|
        values = NamedArray.setup(values, tsv.fields)
        FIGURE_TREATMENT_ORDER.each do |treatment|
          %w(positive negative).each do |sign|
            FIGURE_TIME_POINTS.each do |time|
              field = "#{treatment}-T#{time}"
              next unless tsv.fields.include?(field)
              value = AGS.figure_safe_float(values[field], 0.0)
              next unless (sign == 'positive' && value > 0) || (sign == 'negative' && value < 0)
              map[[tf, treatment, sign]] ||= time
            end
          end
        end
      end
      maps[scheme_name] = map
    end
    counts = Hash.new(0)
    (maps['dynamic'].keys | maps['non-dynamic'].keys).each do |tf,treatment,sign|
      d = maps['dynamic'][[tf,treatment,sign]]
      n = maps['non-dynamic'][[tf,treatment,sign]]
      category = if d && n
                   FIGURE_TIME_POINTS.index(d) < FIGURE_TIME_POINTS.index(n) ? 'dynamic earlier' : (FIGURE_TIME_POINTS.index(d) > FIGURE_TIME_POINTS.index(n) ? 'non-dynamic earlier' : 'same time')
                 elsif d
                   'dynamic only'
                 else
                   'non-dynamic only'
                 end
      counts[[treatment, sign, category]] += 1
    end
    rows = counts.each_with_index.collect{|((treatment, sign, category), count), i| [i+1, treatment, sign, category, count] }
    plot_tsv = AGS.figure_build_tsv(rows, 'ID', %w(Treatment Sign Category Count))
    figure_ggplot(plot_tsv, <<-RCODE, 10, 5)
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Count <- as.numeric(data$Count)
ggplot(data, aes(TreatmentLabel, Count, fill=Category)) + geom_col() + facet_wrap(~Sign, nrow=1) +
  scale_fill_brewer(palette='Set2') + coord_flip() + theme_bw() + labs(x='', y='TF-treatment-sign calls', fill='Earliest detection') + theme(legend.position='bottom')
RCODE
  end

  dep :tf_predictions, :scheme => 'dynamic'
  dep :tf_predictions, :scheme => 'non-dynamic'
  extension :svg
  task :figure_dynamic_vs_nondynamic_persistence => :text do
    rows = []
    id = 0
    dependencies.each do |dep|
      scheme_name = dep.recursive_inputs[:scheme].to_s
      tsv = dep.load
      tsv.through do |tf, values|
        values = NamedArray.setup(values, tsv.fields)
        FIGURE_TREATMENT_ORDER.each do |treatment|
          %w(positive negative).each do |sign|
            count = 0
            FIGURE_TIME_POINTS.each do |time|
              field = "#{treatment}-T#{time}"
              next unless tsv.fields.include?(field)
              value = AGS.figure_safe_float(values[field], 0.0)
              count += 1 if (sign == 'positive' && value > 0) || (sign == 'negative' && value < 0)
            end
            next if count == 0
            id += 1
            rows << [id, scheme_name, treatment, sign, count]
          end
        end
      end
    end
    plot_tsv = AGS.figure_build_tsv(rows, 'ID', %w(Scheme Treatment Sign ActiveTimepoints))
    figure_ggplot(plot_tsv, <<-RCODE, 9, 5)
data$ActiveTimepoints <- as.numeric(data$ActiveTimepoints)
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
ggplot(data, aes(factor(ActiveTimepoints), fill=Scheme)) + geom_bar(position='dodge') + facet_grid(Sign ~ TreatmentLabel) +
  scale_fill_manual(values=c('dynamic'='#D73027','non-dynamic'='#4575B4')) + theme_bw() +
  labs(x='Active sampled windows per TF-treatment-sign', y='Count', fill='Scheme') + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
RCODE
  end

  dep :results_table_sequence_edge_counts
  extension :svg
  task :figure_sequence_edge_counts_bar => :text do
    tsv = step(:results_table_sequence_edge_counts).load
    rows = []
    id = 0
    tsv.through do |treatment, values|
      values = NamedArray.setup(values, tsv.fields)
      id += 1
      rows << [id, treatment, 'all sequence edges', AGS.figure_safe_int(values['Total sequence edges'], 0)]
      id += 1
      rows << [id, treatment, 'both TFs self-consistent', AGS.figure_safe_int(values['Edges with both TFs self-consistent'], 0)]
    end
    plot_tsv = AGS.figure_build_tsv(rows, 'ID', %w(Treatment Class Count))
    figure_ggplot(plot_tsv, <<-RCODE, 8, 5)
data$Count <- as.numeric(data$Count)
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
ggplot(data, aes(TreatmentLabel, Count, fill=Class)) + geom_col(position='dodge') + coord_flip() +
  scale_fill_manual(values=c('all sequence edges'='grey60','both TFs self-consistent'='#1B7837')) +
  theme_bw() + labs(x='', y='Candidate TF-to-TF edges', fill='Sequence class') + theme(legend.position='bottom')
RCODE
  end

  extension :svg
  task :figure_three_layer_interpretation_scaffold => :text do
    tsv = AGS.figure_build_tsv([['one', 1]], 'ID', ['Value'])
    figure_ggplot(tsv, <<-RCODE, 10, 5)
ggplot(data.frame(x=1,y=1), aes(x,y)) + xlim(0,10) + ylim(0,5) +
  annotate('rect', xmin=0.6, xmax=3.0, ymin=3.2, ymax=4.4, fill='#DDEBF7', color='grey30') +
  annotate('rect', xmin=3.8, xmax=6.2, ymin=3.2, ymax=4.4, fill='#FFF2CC', color='grey30') +
  annotate('rect', xmin=7.0, xmax=9.4, ymin=3.2, ymax=4.4, fill='#E2F0D9', color='grey30') +
  annotate('text', x=1.8, y=3.8, label='Signaling\ncascade states', size=4) +
  annotate('text', x=5.0, y=3.8, label='TF activity-change\nscoreboards', size=4) +
  annotate('text', x=8.2, y=3.8, label='Functional process\nscoreboards', size=4) +
  annotate('segment', x=3.05, xend=3.75, y=3.8, yend=3.8, arrow=arrow(length=unit(0.15,'inches'))) +
  annotate('segment', x=6.25, xend=6.95, y=3.8, yend=3.8, arrow=arrow(length=unit(0.15,'inches'))) +
  annotate('rect', xmin=1.2, xmax=8.8, ymin=1.0, ymax=2.2, fill='#FCE4D6', color='grey30') +
  annotate('text', x=5.0, y=1.6, label='Structured mechanistic narrative: time-resolved, evidence-linked, and testable', size=4) +
  annotate('segment', x=5.0, xend=5.0, y=3.15, yend=2.25, arrow=arrow(length=unit(0.15,'inches'))) +
  theme_void()
RCODE
  end

  dep :figure_overall_intermediate_phenotypes
  dep :figure_overall_dynamic_workflow
  dep :figure_mrna_pca_by_treatment
  dep :figure_mrna_pca_by_time
  dep :figure_mrna_context_distance_clustered
  dep :figure_de_counts_fc0_bar
  dep :figure_de_counts_fc1_bar
  dep :figure_dynamic_onset_first_bar
  dep :figure_dynamic_onset_episode_bar
  dep :figure_dynamic_onset_switch_summary
  dep :figure_dynamic_vs_interval_onset_relationship
  dep :figure_functional_enrichment_scoreboard, :query_type => 'onset_genes', :annotation => 'cancerhallmarks_core', :column_grouping => 'time_treatment', :max_terms => 20
  dep :figure_functional_enrichment_scoreboard, :query_type => 'onset_genes', :annotation => 'cancerhallmarks_core', :column_grouping => 'treatment_time', :max_terms => 20
  dep :figure_functional_enrichment_scoreboard, :query_type => 'onset_genes', :annotation => 'goslim', :column_grouping => 'time_treatment', :max_terms => 28
  dep :figure_functional_enrichment_scoreboard, :query_type => 'tf_and_targets', :annotation => 'cancerhallmarks_core', :column_grouping => 'treatment_time', :max_terms => 20
  dep :figure_functional_enrichment_term_frequency, :query_type => 'onset_genes', :annotation => 'cancerhallmarks_core'
  dep :figure_tf_activity_counts_bar
  dep :figure_tf_activity_pca_by_treatment
  dep :figure_tf_activity_sign_balance_bar
  dep :figure_tf_activity_heatmap_dendrogram, :scheme => 'dynamic', :normalization => 'raw', :max_abs => 50.0, :cluster_columns => false
  dep :figure_tf_activity_heatmap_dendrogram, :scheme => 'dynamic', :normalization => 'row_zscore', :max_abs => 3.0, :cluster_columns => false
  dep :figure_tf_activity_heatmap_dendrogram, :scheme => 'dynamic', :normalization => 'row_zscore', :max_abs => 3.0, :cluster_columns => true
  dep :figure_combination_tf_activity_upset_bar, :combination => 'INT_FiveZ_PI', :time_point => 2, :sign_mode => 'all'
  dep :figure_combination_tf_activity_upset_bar, :combination => 'INT_FiveZ_PI', :time_point => 24, :sign_mode => 'all'
  dep :figure_combination_tf_activity_upset_bar, :combination => 'INT_PD_PI', :time_point => 1, :sign_mode => 'all'
  dep :figure_combination_tf_activity_upset_bar, :combination => 'INT_PD_PI', :time_point => 24, :sign_mode => 'all'
  dep :figure_combination_category_counts_bar
  dep :figure_opposite_early_single_agent_tf_signs, :time_point => 1
  dep :figure_opposite_early_single_agent_tf_signs, :time_point => 2
  dep :figure_tf_activity_pair_scatter, :time_point => 1, :treatment_x => 'PI', :treatment_y => 'FiveZ'
  dep :figure_tf_activity_pair_scatter, :time_point => 1, :treatment_x => 'PI', :treatment_y => 'PD'
  dep :figure_tf_activity_pair_scatter, :time_point => 1, :treatment_x => 'FiveZ', :treatment_y => 'PD'
  dep :figure_tf_activity_pair_scatter, :time_point => 2, :treatment_x => 'PI', :treatment_y => 'FiveZ'
  dep :figure_tf_activity_pair_scatter, :time_point => 2, :treatment_x => 'PI', :treatment_y => 'PD'
  dep :figure_tf_activity_pair_scatter, :time_point => 2, :treatment_x => 'FiveZ', :treatment_y => 'PD'
  dep :figure_neko_match_fraction_bar
  dep :figure_neko_counts_bar
  dep :figure_self_consistency_match_fraction_bar
  dep :figure_self_consistency_counts_bar
  dep :figure_dynamic_vs_nondynamic_earliest_detection
  dep :figure_dynamic_vs_nondynamic_persistence
  dep :figure_sequence_edge_counts_bar
  dep :figure_three_layer_interpretation_scaffold
  task :figures_manuscript => :array do
    require 'fileutils'
    FileUtils.rm_rf(files_dir) if File.directory?(files_dir)
    FileUtils.mkdir_p(files_dir)
    dependencies.each do |dep|
      next unless dep.done?
      inputs = dep.recursive_inputs
      task_name = dep.task_name.to_s
      tags = []
      if task_name == 'figure_tf_activity_heatmap_dendrogram'
        tags << inputs[:scheme]
        tags << inputs[:normalization]
        tags << ('columns_clustered' if inputs[:cluster_columns].to_s == 'true')
      elsif task_name == 'figure_combination_tf_activity_upset_bar'
        tags << inputs[:combination]
        tags << "T#{inputs[:time_point]}"
        tags << inputs[:sign_mode]
      elsif task_name == 'figure_opposite_early_single_agent_tf_signs'
        tags << "T#{inputs[:time_point]}"
      elsif task_name == 'figure_tf_activity_pair_scatter'
        tags << inputs[:treatment_x]
        tags << 'vs'
        tags << inputs[:treatment_y]
        tags << "T#{inputs[:time_point]}"
      elsif task_name == 'figure_functional_enrichment_scoreboard'
        tags << inputs[:query_type]
        tags << inputs[:annotation]
        tags << inputs[:column_grouping]
      elsif task_name == 'figure_functional_enrichment_term_frequency'
        tags << inputs[:query_type]
        tags << inputs[:annotation]
      end
      name = ([task_name.sub(/^figure_/, '')] + tags.compact.collect(&:to_s)) * '-'
      target = file(name + '.svg')
      png = file(name + '.png')
      Open.cp(dep.path, target)
      Open.cp(dep.file('plot.png'), png) if dep.file('plot.png').exists?
    end
    files
  end

  dep_task :figures_current, AGS, :figures_manuscript
end
