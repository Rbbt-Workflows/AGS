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

  FIGURE_COMBINATION_CATEGORY_LABELS = {
    'combination_earlier_than_both' => 'combination earlier than both',
    'combination_specific_at_time' => 'combination specific at time',
    'shared_with_both_same_sign' => 'shared with both, same sign',
    'shared_with_component1_same_sign' => 'shared with PI, same sign',
    'shared_with_component2_same_sign' => 'shared with second component, same sign',
    'sign_reversed_relative_to_component' => 'sign reversed relative to component'
  } unless const_defined?(:FIGURE_COMBINATION_CATEGORY_LABELS)

  FIGURE_GO_THEME_KEYWORDS = {
    'cell cycle' => %w(cell cycle mitotic chromosome spindle cytokinesis replication),
    'DNA damage and repair' => %w(repair checkpoint damage recombination),
    'RNA and ribosome' => %w(ribosome ribosomal rRNA RNA transcription splicing),
    'translation and protein homeostasis' => %w(translation proteasome ubiquitin folding unfolded endoplasmic),
    'stress response' => %w(stress oxidative hypoxia heat shock),
    'inflammatory signaling' => %w(inflammatory cytokine interferon NF-kappa immune),
    'metabolism' => %w(metabolic metabolism glycolysis lipid amino mitochondrial respiration),
    'cell death' => %w(apoptotic apoptosis death autophagy),
    'adhesion and migration' => %w(adhesion migration motility extracellular matrix morphogenesis),
    'signaling' => %w(signaling signal phosphorylation kinase MAPK PI3K)
  } unless const_defined?(:FIGURE_GO_THEME_KEYWORDS)

  FIGURE_GENERAL_GO_TERMS = [
    'biological process',
    'cellular process',
    'metabolic process',
    'biological regulation',
    'regulation of biological process',
    'regulation of cellular process',
    'response to stimulus',
    'multicellular organismal process',
    'developmental process',
    'localization',
    'cellular component organization'
  ] unless const_defined?(:FIGURE_GENERAL_GO_TERMS)

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

  def self.figure_parse_context(field)
    case field.to_s
    when /^(.+)-T(\d+)$/
      [$1, $2.to_i]
    when /^(?:FC|Pvalue|FC1|PvalueFC1)_(.+)\.T(\d+)$/
      [$1, $2.to_i]
    else
      nil
    end
  end

  def self.figure_ggplot(file, tsv, r_code, width = 8, height = 5)
    R::PNG.ggplot file, tsv, r_code, width, height
    nil
  end

  def self.figure_cap(value, cap)
    [[value.to_f, cap.to_f].min, -cap.to_f].max
  end

  def self.figure_gprofiler_rows(files_dir, source_type, pvalue_cutoff, max_terms, selected_terms = nil)
    all_rows = []
    term_counts = Hash.new(0)
    term_best = Hash.new(1.0)

    FIGURE_TREATMENT_ORDER.each do |treatment|
      FIGURE_TIME_POINTS.each do |time_point|
        %w(up down).each do |direction|
          path = File.join(files_dir, [treatment, "#{time_point}h", source_type, direction] * '-' + '.tsv')
          next unless File.exist?(path)
          tsv = TSV.open(path, :type => :list)
          name_idx = tsv.fields.index('name') || 2
          p_idx = tsv.fields.index('p-value') || tsv.fields.index('p_value') || 3
          sig_idx = tsv.fields.index('significant')
          inter_idx = tsv.fields.index('intersection_size')
          query_idx = tsv.fields.index('query_size')
          tsv.through do |id, values|
            name = values[name_idx].to_s
            next if name.empty?
            next if FIGURE_GENERAL_GO_TERMS.include?(name)
            pvalue = AGS.figure_safe_float(values[p_idx], 1.0)
            next if pvalue <= 0 || pvalue > pvalue_cutoff
            if sig_idx
              sig = values[sig_idx].to_s
              next if sig == 'false'
            end
            term_counts[name] += 1
            term_best[name] = [term_best[name], pvalue].min
            all_rows << {
              :term => name,
              :treatment => treatment,
              :time => time_point,
              :direction => direction,
              :pvalue => pvalue,
              :neglogp => -Math.log10(pvalue),
              :intersection => inter_idx ? AGS.figure_safe_int(values[inter_idx], 0) : 0,
              :query_size => query_idx ? AGS.figure_safe_int(values[query_idx], 0) : 0
            }
          end
        end
      end
    end

    selected = if selected_terms && selected_terms.any?
                 selected_terms
               else
                 term_counts.keys.sort_by{|term| [-term_counts[term], term_best[term], term] }.first(max_terms)
               end

    all_rows.select{|row| selected.include?(row[:term]) }
  end

  def self.figure_go_theme(term)
    down = term.to_s.downcase
    FIGURE_GO_THEME_KEYWORDS.each do |theme, words|
      return theme if words.any?{|word| down.include?(word.downcase) }
    end
    'other'
  end

  def self.figure_format_number(value)
    return '' if value.nil?
    if value.to_f.nan? || value.to_f.infinite?
      ''
    else
      value.to_s
    end
  end

  #############################################################################
  # Figure 1 concept panels
  #############################################################################

  extension :png
  task :figure_01b_concept_intermediate_phenotypes => :binary do
    tsv = AGS.figure_build_tsv([['one', 1]], 'ID', ['Value'])
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 9, 4)
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

  extension :png
  task :figure_01c_concept_dynamic_workflow => :binary do
    tsv = AGS.figure_build_tsv([['one', 1]], 'ID', ['Value'])
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 10, 4.5)
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
  annotate('text', x=7.5, y=2.95, label='TF activity\ncalls', size=3.7) +
  annotate('text', x=9.2, y=2.95, label='Process\nchronology', size=3.7) +
  annotate('segment', x=1.85, xend=2.25, y=2.95, yend=2.95, arrow=arrow(length=unit(0.13,'inches'))) +
  annotate('segment', x=4.15, xend=4.55, y=2.95, yend=2.95, arrow=arrow(length=unit(0.13,'inches'))) +
  annotate('segment', x=6.35, xend=6.75, y=2.95, yend=2.95, arrow=arrow(length=unit(0.13,'inches'))) +
  annotate('segment', x=8.25, xend=8.55, y=2.95, yend=2.95, arrow=arrow(length=unit(0.13,'inches'))) +
  annotate('text', x=5.0, y=1.4, label='Only the dynamic TF activity scheme is shown here; fc0 and diff variants are left out of the main figure.', size=3.6) +
  theme_void()
RCODE
  end

  #############################################################################
  # Figure 2: transcriptome PCA and treatment distance panels C, D, E
  #############################################################################

  dep :fold_changes, :fc_source => 'NTNU'
  extension :png
  task :figure_02c_mrna_pca_treatment => :binary do
    fc0 = step(:fold_changes).load.transpose('Associated Gene Name')
    fields = AGS.figure_context_fields('FC').select{|field| fc0.fields.include?(field) }
    data_tsv = fc0.reorder('Associated Gene Name', fields)
    labels = AGS.figure_treatment_levels
    colors = AGS.figure_treatment_color_scale('color')
    AGS.figure_ggplot(self.tmp_path, data_tsv, <<-RCODE, 7, 5.5)
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
plot_df$TreatmentLabel <- factor(plot_df$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{labels}))
var_exp <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
#{colors}
ggplot(plot_df, aes(PC1, PC2, color=TreatmentLabel, group=TreatmentLabel)) +
  geom_path(alpha=0.6) + geom_point(aes(size=Time), alpha=0.9) +
  scale_size_continuous(breaks=c(1,2,4,8,24), name='Time') +
  theme_bw() + labs(x=paste0('PC1 (', var_exp[1], '%)'), y=paste0('PC2 (', var_exp[2], '%)'), color='Treatment') +
  theme(legend.position='bottom')
RCODE
  end

  dep :fold_changes, :fc_source => 'NTNU'
  extension :png
  task :figure_02d_mrna_pca_time => :binary do
    fc0 = step(:fold_changes).load.transpose('Associated Gene Name')
    fields = AGS.figure_context_fields('FC').select{|field| fc0.fields.include?(field) }
    data_tsv = fc0.reorder('Associated Gene Name', fields)
    labels = AGS.figure_treatment_levels
    AGS.figure_ggplot(self.tmp_path, data_tsv, <<-RCODE, 7, 5.5)
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
plot_df$TreatmentLabel <- factor(plot_df$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{labels}))
var_exp <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
ggplot(plot_df, aes(PC1, PC2, color=Time, shape=TreatmentLabel)) +
  geom_point(size=3, alpha=0.9) +
  scale_color_brewer(palette='YlOrRd') +
  theme_bw() + labs(x=paste0('PC1 (', var_exp[1], '%)'), y=paste0('PC2 (', var_exp[2], '%)'), shape='Treatment') +
  theme(legend.position='bottom')
RCODE
  end

  dep :fold_changes, :fc_source => 'NTNU'
  extension :png
  task :figure_02e_mrna_treatment_distance => :binary do
    fc0 = step(:fold_changes).load.transpose('Associated Gene Name')
    fields = AGS.figure_context_fields('FC').select{|field| fc0.fields.include?(field) }
    data_tsv = fc0.reorder('Associated Gene Name', fields)
    labels = FIGURE_TREATMENT_ORDER.collect{|t| AGS.figure_treatment_label(t) }
    AGS.figure_ggplot(self.tmp_path, data_tsv, <<-RCODE, 6, 5.5)
mat0 <- data[, c(#{fields.collect{|f| "'#{f}'"} * ', '}), drop=FALSE]
mat0[] <- lapply(mat0, as.numeric)
mat0[is.na(mat0)] <- 0
mat <- t(as.matrix(mat0))
treatments <- c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','})
labels <- c(#{labels.collect{|l| "'#{l}'"} * ','})
traj <- lapply(treatments, function(tr) as.vector(mat[startsWith(rownames(mat), paste0('FC_', tr, '.T')), , drop=FALSE]))
names(traj) <- labels
traj <- do.call(rbind, traj)
d <- as.matrix(as.dist(1 - cor(t(traj), use='pairwise.complete.obs')))
df <- as.data.frame(as.table(d))
names(df) <- c('Treatment1','Treatment2','Distance')
df$Treatment1 <- factor(df$Treatment1, levels=labels)
df$Treatment2 <- factor(df$Treatment2, levels=rev(labels))
ggplot(df, aes(Treatment1, Treatment2, fill=Distance)) + geom_tile(color='white') +
  scale_fill_gradient(low='white', high='#1F78B4', name='1 - r') +
  coord_fixed() + theme_bw() + labs(x='', y='') +
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank())
RCODE
  end

  #############################################################################
  # Figure 3: DE count panels B and C
  #############################################################################

  dep :de_gene_counts_fc0
  extension :png
  task :figure_03b_fc0_de_counts => :binary do
    counts = step(:de_gene_counts_fc0).load
    AGS.figure_ggplot(self.tmp_path, counts, <<-RCODE, 8, 5)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Genes <- as.numeric(data$Genes)
data <- subset(data, Direction == 'both')
#{AGS.figure_treatment_color_scale('color')}
ggplot(data, aes(Time, Genes, color=TreatmentLabel, group=TreatmentLabel)) + geom_line(linewidth=1) + geom_point(size=2) +
  theme_bw() + labs(x='Time point', y='DE genes relative to baseline', color='Treatment') +
  theme(legend.position='bottom')
RCODE
  end

  dep :interval_de_gene_counts_fc1
  extension :png
  task :figure_03c_fc1_de_counts => :binary do
    counts = step(:interval_de_gene_counts_fc1).load
    AGS.figure_ggplot(self.tmp_path, counts, <<-RCODE, 8, 5)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Genes <- as.numeric(data$Genes)
data <- subset(data, Direction == 'both')
#{AGS.figure_treatment_color_scale('color')}
ggplot(data, aes(Time, Genes, color=TreatmentLabel, group=TreatmentLabel)) + geom_line(linewidth=1) + geom_point(size=2) +
  theme_bw() + labs(x='Interval ending at time point', y='Interval DE genes using FC1 p-value surrogate', color='Treatment') +
  theme(legend.position='bottom')
RCODE
  end

  #############################################################################
  # Figure 4: dynamic onset panels
  #############################################################################

  dep :onset_first_counts
  extension :png
  task :figure_04a_onset_first_counts => :binary do
    tsv = step(:onset_first_counts).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 10, 5)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Direction <- factor(data$Direction, levels=c('up','down'))
data$Genes <- as.numeric(data$Genes)
data <- subset(data, Direction != 'both')
#{AGS.figure_direction_color_scale('fill')}
ggplot(data, aes(Time, Genes, fill=Direction)) + geom_col(position='stack') +
  facet_wrap(~TreatmentLabel, nrow=2) + theme_bw() + labs(x='First onset time', y='Genes', fill='Direction') +
  theme(legend.position='bottom')
RCODE
  end

  dep :onset_episode_counts
  extension :png
  task :figure_04b_onset_episode_counts => :binary do
    tsv = step(:onset_episode_counts).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 10, 5)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Direction <- factor(data$Direction, levels=c('up','down'))
data$Episodes <- as.numeric(data$Episodes)
data <- subset(data, Direction != 'both')
#{AGS.figure_direction_color_scale('fill')}
ggplot(data, aes(Time, Episodes, fill=Direction)) + geom_col(position='stack') +
  facet_wrap(~TreatmentLabel, nrow=2) + theme_bw() + labs(x='Onset episode time', y='Onset episodes', fill='Direction') +
  theme(legend.position='bottom')
RCODE
  end

  dep :onset_direction_switch_summary
  extension :png
  task :figure_04c_onset_switch_summary => :binary do
    tsv = step(:onset_direction_switch_summary).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 8, 5)
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Category <- factor(data$Category, levels=c('unclassified','single_episode','multiple_same_direction','direction_switch'))
data$Genes <- as.numeric(data$Genes)
ggplot(data, aes(TreatmentLabel, Genes, fill=Category)) + geom_col() +
  scale_fill_brewer(palette='Set2') + theme_bw() + labs(x='', y='Genes', fill='Trajectory class') +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom')
RCODE
  end

  dep :full_gene_info
  input :example_treatment, :select, 'Treatment for onset example trajectories', 'PI', :select_options => FIGURE_TREATMENT_ORDER
  input :example_onset, :string, 'Onset label to show', 'increase 2h'
  input :max_genes, :integer, 'Maximum genes to show', 40
  extension :png
  task :figure_04d_onset_example_profiles => :binary do |example_treatment, example_onset, max_genes|
    fields = FIGURE_TIME_POINTS.collect{|time| "FC_#{example_treatment}.T#{time}" } + ["#{example_treatment}: FC clusters"]
    info = step(:full_gene_info).load.reorder('Associated Gene Name', fields)
    cluster_field = "#{example_treatment}: FC clusters"
    cluster_idx = info.fields.index(cluster_field)
    genes = []
    info.through do |gene, values|
      label = values[cluster_idx].to_s
      genes << gene if label.split('|').include?(example_onset)
      break if genes.length >= max_genes
    end
    rows = []
    genes.each do |gene|
      values = info[gene]
      FIGURE_TIME_POINTS.each do |time|
        field = "FC_#{example_treatment}.T#{time}"
        idx = info.fields.index(field)
        rows << ["#{gene}-#{time}", gene, time, AGS.figure_safe_float(values[idx]), 'gene']
      end
    end
    FIGURE_TIME_POINTS.each do |time|
      field = "FC_#{example_treatment}.T#{time}"
      idx = info.fields.index(field)
      vals = genes.collect{|gene| AGS.figure_safe_float(info[gene][idx]) }
      mean = vals.empty? ? 0.0 : vals.inject(0.0, &:+) / vals.length
      rows << ["mean-#{time}", 'mean', time, mean, 'mean']
    end
    plot_tsv = AGS.figure_build_tsv(rows, 'ID', %w(Gene Time FC Type))
    AGS.figure_ggplot(self.tmp_path, plot_tsv, <<-RCODE, 6.5, 5)
data$Time <- as.numeric(data$Time)
data$FC <- as.numeric(data$FC)
ggplot(data, aes(Time, FC, group=Gene)) +
  geom_line(data=subset(data, Type == 'gene'), color='grey70', alpha=0.5) +
  geom_line(data=subset(data, Type == 'mean'), color='#D73027', linewidth=1.2) +
  geom_point(data=subset(data, Type == 'mean'), color='#D73027', size=2) +
  geom_hline(yintercept=0, linetype='dashed', color='grey50') +
  scale_x_continuous(breaks=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h')) +
  theme_bw() + labs(x='Time point', y='log2 fold change', title='#{example_treatment}: #{example_onset}', subtitle='#{genes.length} example genes')
RCODE
  end

  #############################################################################
  # Figure 5: FC1 and onset relationship panels
  #############################################################################

  dep :interval_de_gene_counts_fc1
  extension :png
  task :figure_05a_fc1_de_direction_counts => :binary do
    tsv = step(:interval_de_gene_counts_fc1).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 10, 5)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Direction <- factor(data$Direction, levels=c('up','down'))
data$Genes <- as.numeric(data$Genes)
data <- subset(data, Direction != 'both')
#{AGS.figure_direction_color_scale('fill')}
ggplot(data, aes(Time, Genes, fill=Direction)) + geom_col(position='stack') +
  facet_wrap(~TreatmentLabel, nrow=2) + theme_bw() +
  labs(x='Interval ending at time point', y='Interval DE genes', fill='FC1 direction') +
  theme(legend.position='bottom')
RCODE
  end

  dep :fc1_onset_relationship_summary
  extension :png
  task :figure_05b_fc1_onset_relationship => :binary do
    tsv = step(:fc1_onset_relationship_summary).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 12, 6)
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

  dep :fc1_onset_relationship_8h_focus
  extension :png
  task :figure_05e_fc1_onset_8h_focus => :binary do
    tsv = step(:fc1_onset_relationship_8h_focus).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 8, 5)
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$Genes <- as.numeric(data$Genes)
ggplot(data, aes(TreatmentLabel, Genes, fill=Relationship)) + geom_col() +
  facet_wrap(~Direction, nrow=1) + scale_fill_brewer(palette='Paired') +
  theme_bw() + labs(x='', y='8 h FC1 DE genes', fill='Onset relationship') +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom', legend.text=element_text(size=7))
RCODE
  end

  #############################################################################
  # Figure 6: GO enrichment investigation panels
  #############################################################################

  dep :gprofiler_suite
  input :source_type, :select, 'gProfiler query type', 'cluster', :select_options => %w(cluster TF fc0_03 fc0_07 fc_03 fc_07)
  input :pvalue_cutoff, :float, 'Maximum enrichment p-value', 1e-6
  input :max_terms, :integer, 'Maximum GO terms to show', 25
  extension :png
  task :figure_06a_go_top_terms_dotplot => :binary do |source_type, pvalue_cutoff, max_terms|
    suite = step(:gprofiler_suite)
    rows = AGS.figure_gprofiler_rows(suite.files_dir, source_type, pvalue_cutoff, max_terms)
    plot_rows = rows.each_with_index.collect do |row, i|
      id = i + 1
      [id, row[:term], AGS.figure_treatment_label(row[:treatment]), row[:time], row[:direction], row[:neglogp], row[:intersection], row[:query_size]]
    end
    plot_tsv = AGS.figure_build_tsv(plot_rows, 'ID', %w(Term Treatment Time Direction NegLogP Intersection QuerySize))
    AGS.figure_ggplot(self.tmp_path, plot_tsv, <<-RCODE, 14, 8)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$Treatment <- factor(data$Treatment, levels=c(#{AGS.figure_treatment_levels}))
data$Direction <- factor(data$Direction, levels=c('up','down'))
data$NegLogP <- as.numeric(data$NegLogP)
data$Intersection <- as.numeric(data$Intersection)
term_order <- names(sort(tapply(data$NegLogP, data$Term, max), decreasing=FALSE))
data$Term <- factor(data$Term, levels=term_order)
data$Context <- factor(paste(data$Treatment, data$Time, data$Direction, sep=' '), levels=unique(paste(data$Treatment, data$Time, data$Direction, sep=' ')))
ggplot(data, aes(Context, Term, size=Intersection, color=NegLogP)) + geom_point(alpha=0.85) +
  scale_color_viridis_c(name='-log10(p)') + scale_size_continuous(name='Genes', range=c(1,6)) +
  theme_bw() + labs(x='', y='') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6), axis.text.y=element_text(size=7), legend.position='right')
RCODE
  end

  dep :gprofiler_suite
  input :source_type, :select, 'gProfiler query type', 'cluster', :select_options => %w(cluster TF fc0_03 fc0_07 fc_03 fc_07)
  input :pvalue_cutoff, :float, 'Maximum enrichment p-value', 1e-6
  input :max_terms, :integer, 'Maximum GO terms to show', 40
  extension :png
  task :figure_06b_go_term_frequency => :binary do |source_type, pvalue_cutoff, max_terms|
    rows = AGS.figure_gprofiler_rows(step(:gprofiler_suite).files_dir, source_type, pvalue_cutoff, max_terms * 4)
    counts = Hash.new{|h,k| h[k] = [0, 1.0] }
    rows.each do |row|
      counts[row[:term]][0] += 1
      counts[row[:term]][1] = [counts[row[:term]][1], row[:pvalue]].min
    end
    plot_rows = counts.keys.sort_by{|term| [-counts[term][0], counts[term][1], term] }.first(max_terms).each_with_index.collect do |term, i|
      [i + 1, term, counts[term][0], -Math.log10(counts[term][1])]
    end
    plot_tsv = AGS.figure_build_tsv(plot_rows, 'ID', %w(Term Occurrences BestNegLogP))
    AGS.figure_ggplot(self.tmp_path, plot_tsv, <<-RCODE, 8, 8)
data$Occurrences <- as.numeric(data$Occurrences)
data$BestNegLogP <- as.numeric(data$BestNegLogP)
data$Term <- factor(data$Term, levels=rev(data$Term[order(data$Occurrences, data$BestNegLogP)]))
ggplot(data, aes(Term, Occurrences, fill=BestNegLogP)) + geom_col() + coord_flip() +
  scale_fill_viridis_c(name='best -log10(p)') + theme_bw() + labs(x='', y='Number of enriched treatment-time-direction contexts')
RCODE
  end

  dep :gprofiler_suite
  input :source_type, :select, 'gProfiler query type', 'cluster', :select_options => %w(cluster TF fc0_03 fc0_07 fc_03 fc_07)
  input :pvalue_cutoff, :float, 'Maximum enrichment p-value', 1e-6
  extension :png
  task :figure_06c_go_theme_heatmap => :binary do |source_type, pvalue_cutoff|
    rows = AGS.figure_gprofiler_rows(step(:gprofiler_suite).files_dir, source_type, pvalue_cutoff, 1000)
    theme_scores = Hash.new(0.0)
    rows.each do |row|
      theme = AGS.figure_go_theme(row[:term])
      next if theme == 'other'
      key = [theme, AGS.figure_treatment_label(row[:treatment]), row[:time], row[:direction]]
      theme_scores[key] = [theme_scores[key], row[:neglogp]].max
    end
    plot_rows = []
    id = 0
    theme_scores.keys.sort.each do |theme, treatment, time, direction|
      id += 1
      plot_rows << [id, theme, treatment, time, direction, theme_scores[[theme, treatment, time, direction]]]
    end
    plot_tsv = AGS.figure_build_tsv(plot_rows, 'ID', %w(Theme Treatment Time Direction Score))
    AGS.figure_ggplot(self.tmp_path, plot_tsv, <<-RCODE, 12, 6)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$Treatment <- factor(data$Treatment, levels=c(#{AGS.figure_treatment_levels}))
data$Direction <- factor(data$Direction, levels=c('up','down'))
data$Score <- as.numeric(data$Score)
data$Context <- factor(paste(data$Treatment, data$Time, data$Direction, sep=' '), levels=unique(paste(data$Treatment, data$Time, data$Direction, sep=' ')))
ggplot(data, aes(Context, Theme, fill=Score)) + geom_tile(color='white') +
  scale_fill_viridis_c(name='max -log10(p)') + theme_bw() + labs(x='', y='GO theme') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6), panel.grid=element_blank())
RCODE
  end

  dep :gprofiler_suite
  input :source_type, :select, 'gProfiler query type', 'cluster', :select_options => %w(cluster TF fc0_03 fc0_07 fc_03 fc_07)
  input :pvalue_cutoff, :float, 'Maximum enrichment p-value', 1e-6
  input :selected_terms, :array, 'GO terms to show', ['cell cycle process', 'DNA replication', 'DNA repair', 'RNA processing', 'ribosome biogenesis', 'apoptotic process', 'response to stress', 'cell adhesion']
  extension :png
  task :figure_06d_go_selected_terms => :binary do |source_type, pvalue_cutoff, selected_terms|
    rows = AGS.figure_gprofiler_rows(step(:gprofiler_suite).files_dir, source_type, pvalue_cutoff, selected_terms.length, selected_terms)
    plot_rows = rows.each_with_index.collect do |row, i|
      [i + 1, row[:term], AGS.figure_treatment_label(row[:treatment]), row[:time], row[:direction], row[:neglogp], row[:intersection]]
    end
    plot_tsv = AGS.figure_build_tsv(plot_rows, 'ID', %w(Term Treatment Time Direction NegLogP Intersection))
    AGS.figure_ggplot(self.tmp_path, plot_tsv, <<-RCODE, 12, 5)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$Treatment <- factor(data$Treatment, levels=c(#{AGS.figure_treatment_levels}))
data$Term <- factor(data$Term, levels=rev(c(#{selected_terms.collect{|t| "'#{t}'"} * ','})))
data$NegLogP <- as.numeric(data$NegLogP)
data$Intersection <- as.numeric(data$Intersection)
data$Context <- factor(paste(data$Treatment, data$Time, data$Direction, sep=' '), levels=unique(paste(data$Treatment, data$Time, data$Direction, sep=' ')))
ggplot(data, aes(Context, Term, size=Intersection, color=NegLogP)) + geom_point() +
  scale_color_viridis_c(name='-log10(p)') + scale_size_continuous(name='Genes', range=c(1,6)) +
  theme_bw() + labs(x='', y='') + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6))
RCODE
  end

  #############################################################################
  # Figure 7: dynamic TF activity call panels B, C, D
  #############################################################################

  dep :tf_activity_call_counts_by_scheme, :scheme => 'dynamic'
  extension :png
  task :figure_07b_tf_activity_call_counts => :binary do
    tsv = step(:tf_activity_call_counts_by_scheme).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 10, 5)
data <- subset(data, Sign != 'both')
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$TFActivityCalls <- as.numeric(data$TFActivityCalls)
data$SignedCalls <- ifelse(data$Sign == 'negative', -data$TFActivityCalls, data$TFActivityCalls)
ggplot(data, aes(TreatmentLabel, SignedCalls, fill=Sign)) + geom_col() +
  facet_wrap(~Time, nrow=1, scales='free_y') +
  scale_fill_manual(values=c('positive'='#D73027','negative'='#4575B4')) + coord_flip() +
  theme_bw() + labs(x='', y='TF activity calls, positive and negative', fill='Activity sign') +
  theme(legend.position='bottom')
RCODE
  end

  dep :tf_activity_heatmap_matrix, :scheme => 'dynamic', :normalization => 'raw'
  extension :png
  task :figure_07c_tf_activity_pca => :binary do
    matrix = step(:tf_activity_heatmap_matrix).load
    fields = AGS.figure_tf_context_fields(matrix)
    data_tsv = matrix.reorder('Associated Gene Name', fields)
    AGS.figure_ggplot(self.tmp_path, data_tsv, <<-RCODE, 7, 5.5)
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
ggplot(plot_df, aes(PC1, PC2, color=TreatmentLabel, group=TreatmentLabel)) + geom_path(alpha=0.6) + geom_point(aes(size=Time)) +
  scale_size_continuous(breaks=c(1,2,4,8,24), name='Time') + theme_bw() +
  labs(x=paste0('PC1 (', var_exp[1], '%)'), y=paste0('PC2 (', var_exp[2], '%)'), color='Treatment') +
  theme(legend.position='bottom')
RCODE
  end

  dep :tf_activity_call_counts_by_scheme, :scheme => 'dynamic'
  extension :png
  task :figure_07d_tf_activity_sign_balance => :binary do
    tsv = step(:tf_activity_call_counts_by_scheme).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 8, 5)
wide <- reshape(data[, c('Treatment','Time','Sign','TFActivityCalls')], idvar=c('Treatment','Time'), timevar='Sign', direction='wide')
wide$positive <- as.numeric(wide$TFActivityCalls.positive)
wide$negative <- as.numeric(wide$TFActivityCalls.negative)
wide$both <- as.numeric(wide$TFActivityCalls.both)
wide$PositiveFraction <- ifelse(wide$both > 0, wide$positive / wide$both, NA)
wide$Time <- factor(wide$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
wide$TreatmentLabel <- factor(wide$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
#{AGS.figure_treatment_color_scale('color')}
ggplot(wide, aes(Time, PositiveFraction, color=TreatmentLabel, group=TreatmentLabel)) + geom_hline(yintercept=0.5, linetype='dashed', color='grey50') + geom_line(linewidth=1) + geom_point(size=2) +
  ylim(0,1) + theme_bw() + labs(x='Time point', y='Fraction positive among TF activity calls', color='Treatment') +
  theme(legend.position='bottom')
RCODE
  end

  #############################################################################
  # Figure 8: TF report-card and target consistency panels
  #############################################################################

  dep :tf_timepoint_report_card
  input :report_treatment, :select, 'Treatment', 'PI', :select_options => FIGURE_TREATMENT_ORDER
  input :report_time_point, :integer, 'Time point', 2
  input :top_n_tfs, :integer, 'Number of TFs to show', 30
  extension :png
  task :figure_08a_tf_report_card_example => :binary do |report_treatment, report_time_point, top_n_tfs|
    report = step(:tf_timepoint_report_card).load
    idx = Hash[report.fields.each_with_index.to_a]
    selected = []
    report.through do |id, values|
      next unless values[idx['Treatment']].to_s == report_treatment
      next unless values[idx['Time']].to_i == report_time_point.to_i
      next unless values[idx['TFActivityCalled']].to_s == 'true'
      activity = AGS.figure_safe_float(values[idx['TFActivityScore']])
      selected << [values[idx['TF']], activity, values]
    end
    selected = selected.sort_by{|tf, activity, values| -activity.abs }.first(top_n_tfs)
    metrics = ['TFActivityScore', 'TFGeneFC1', 'DynamicTargetsThisTime', 'TargetConcordanceFraction', 'DynamicTargetConcordanceFraction']
    rows = []
    row_id = 0
    selected.each do |tf, activity, values|
      metrics.each do |metric|
        row_id += 1
        rows << [row_id, tf, metric, AGS.figure_safe_float(values[idx[metric]], 0.0)]
      end
      row_id += 1
      consistency = values[idx['SelfConsistency']].to_s
      consistency_value = consistency == '1' ? 1 : (consistency == '-1' ? -1 : 0)
      rows << [row_id, tf, 'SelfConsistency', consistency_value]
    end
    plot_tsv = AGS.figure_build_tsv(rows, 'ID', %w(TF Metric Value))
    AGS.figure_ggplot(self.tmp_path, plot_tsv, <<-RCODE, 8, 8)
data$Value <- as.numeric(data$Value)
data$TF <- factor(data$TF, levels=rev(unique(data$TF)))
data$Metric <- factor(data$Metric, levels=c('TFActivityScore','TFGeneFC1','SelfConsistency','DynamicTargetsThisTime','TargetConcordanceFraction','DynamicTargetConcordanceFraction'))
ggplot(data, aes(Metric, TF, fill=Value)) + geom_tile(color='white') +
  scale_fill_gradient2(low='#4575B4', mid='white', high='#D73027', midpoint=0, name='Value') +
  theme_bw() + labs(x='', y='', title='#{report_treatment} T#{report_time_point} TF report card') +
  theme(axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank())
RCODE
  end

  dep :tf_target_report_card, :called_only => true, :dynamic_targets_only => false
  input :selected_tf, :string, 'TF to show', 'FOXO3'
  input :report_treatment, :select, 'Treatment', 'PI', :select_options => FIGURE_TREATMENT_ORDER
  input :report_time_point, :integer, 'Time point', 2
  input :max_targets, :integer, 'Maximum targets to show', 80
  extension :png
  task :figure_08b_tf_targets_dynamic_highlight => :binary do |selected_tf, report_treatment, report_time_point, max_targets|
    report = step(:tf_target_report_card).load
    idx = Hash[report.fields.each_with_index.to_a]
    rows = []
    report.through do |id, values|
      next unless values[idx['TF']].to_s == selected_tf
      next unless values[idx['Treatment']].to_s == report_treatment
      next unless values[idx['Time']].to_i == report_time_point.to_i
      target = values[idx['Target']].to_s
      fc1 = AGS.figure_safe_float(values[idx['TargetFC1']], 0.0)
      dynamic = values[idx['TargetDynamicThisTime']].to_s == 'true'
      concordant = values[idx['ConcordantWithTFActivity']].to_s
      rows << [target, target, fc1, dynamic ? 'dynamic target' : 'other target', concordant]
      break if rows.length >= max_targets
    end
    plot_tsv = AGS.figure_build_tsv(rows, 'ID', %w(Target FC1 TargetClass Concordant))
    AGS.figure_ggplot(self.tmp_path, plot_tsv, <<-RCODE, 8, 8)
data$FC1 <- as.numeric(data$FC1)
data$Target <- factor(data$Target, levels=rev(data$Target[order(data$FC1)]))
ggplot(data, aes(FC1, Target, color=TargetClass, shape=Concordant)) + geom_point(size=2) +
  geom_vline(xintercept=0, linetype='dashed', color='grey50') +
  scale_color_manual(values=c('dynamic target'='#D73027','other target'='grey60')) +
  theme_bw() + labs(x='Target interval log2 fold change', y='', color='', shape='Concordant', title='#{selected_tf} targets in #{report_treatment} T#{report_time_point}')
RCODE
  end

  dep :tf_target_edge_consistency_summary
  extension :png
  task :figure_08c_tf_target_edge_consistency => :binary do
    tsv = step(:tf_target_edge_consistency_summary).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 9, 5)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$ConcordanceFraction <- as.numeric(data$ConcordanceFraction)
data$TargetClass <- factor(data$TargetClass, levels=c('dynamic_target','other_target'))
ggplot(data, aes(Time, TreatmentLabel, fill=ConcordanceFraction)) + geom_tile(color='white') +
  facet_wrap(~TargetClass, nrow=1) + scale_fill_gradient(low='white', high='#1B7837', limits=c(0,1), name='Concordance') +
  theme_bw() + labs(x='', y='') + theme(panel.grid=element_blank())
RCODE
  end

  #############################################################################
  # Figure 9: TF heatmaps with raw and normalized scores
  #############################################################################

  input :scheme, :select, 'TF prediction scheme', 'dynamic', :select_options => %w(dynamic non-dynamic)
  input :normalization, :select, 'Normalization', 'row_zscore', :select_options => %w(raw row_zscore column_zscore)
  input :max_abs, :float, 'Maximum absolute value for color scale', 3.0
  input :cluster_columns, :boolean, 'Cluster columns', false
  input :show_tf_labels, :boolean, 'Show TF labels', false
  dep :tf_activity_heatmap_matrix, :scheme => :placeholder, :normalization => :placeholder do |jobname, options|
    { :scheme => options[:scheme], :normalization => options[:normalization] }
  end
  extension :png
  task :figure_09_tf_activity_heatmap => :binary do |scheme, normalization, max_abs, cluster_columns, show_tf_labels|
    matrix = step(:tf_activity_heatmap_matrix).load
    fields = AGS.figure_tf_context_fields(matrix)
    data_tsv = matrix.reorder('Associated Gene Name', fields)
    label_code = show_tf_labels ? "element_text(size=4)" : "element_blank()"
    column_order_code = cluster_columns ? "if (ncol(mat) > 2) mat <- mat[, hclust(dist(t(mat)))$order, drop=FALSE]" : "mat <- mat[, c(#{fields.collect{|f| "'#{f}'"} * ','}), drop=FALSE]"
    AGS.figure_ggplot(self.tmp_path, data_tsv, <<-RCODE, 12, 10)
mat <- data[, c(#{fields.collect{|f| "'#{f}'"} * ', '}), drop=FALSE]
mat[] <- lapply(mat, as.numeric)
mat[is.na(mat)] <- 0
mat <- as.matrix(mat)
rownames(mat) <- rownames(data)
if (nrow(mat) > 2) mat <- mat[hclust(dist(mat))$order, , drop=FALSE]
#{column_order_code}
mat[mat > #{max_abs}] <- #{max_abs}
mat[mat < -#{max_abs}] <- -#{max_abs}
df <- as.data.frame(as.table(mat))
names(df) <- c('TF','Context','Activity')
df$TF <- factor(df$TF, levels=rev(unique(df$TF)))
df$Context <- factor(df$Context, levels=unique(df$Context))
ggplot(df, aes(Context, TF, fill=Activity)) + geom_raster() +
  scale_fill_gradient2(low='#4575B4', mid='white', high='#D73027', midpoint=0, limits=c(-#{max_abs}, #{max_abs}), name='Activity') +
  theme_minimal() + labs(x='', y='', title='#{scheme}, #{normalization}') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6), axis.text.y=#{label_code}, panel.grid=element_blank())
RCODE
  end

  dep :tf_activity_heatmap_matrix, :scheme => 'dynamic', :normalization => 'raw'
  extension :png
  task :figure_09d_tf_context_correlation => :binary do
    matrix = step(:tf_activity_heatmap_matrix).load
    fields = AGS.figure_tf_context_fields(matrix)
    data_tsv = matrix.reorder('Associated Gene Name', fields)
    AGS.figure_ggplot(self.tmp_path, data_tsv, <<-RCODE, 9, 8)
mat <- data[, c(#{fields.collect{|f| "'#{f}'"} * ', '}), drop=FALSE]
mat[] <- lapply(mat, as.numeric)
mat[is.na(mat)] <- 0
cor_mat <- cor(as.matrix(mat), use='pairwise.complete.obs')
ord <- hclust(as.dist(1 - cor_mat))$order
cor_mat <- cor_mat[ord, ord]
df <- as.data.frame(as.table(cor_mat))
names(df) <- c('Context1','Context2','Correlation')
df$Context1 <- factor(df$Context1, levels=colnames(cor_mat))
df$Context2 <- factor(df$Context2, levels=rev(colnames(cor_mat)))
ggplot(df, aes(Context1, Context2, fill=Correlation)) + geom_tile() + coord_fixed() +
  scale_fill_gradient2(low='#4575B4', mid='white', high='#D73027', midpoint=0, limits=c(-1,1), name='Pearson r') +
  theme_minimal() + labs(x='', y='') + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6), axis.text.y=element_text(size=6), panel.grid=element_blank())
RCODE
  end

  #############################################################################
  # Figure 10: column-clustered TF heatmap and context similarities
  #############################################################################

  input :scheme, :select, 'TF prediction scheme', 'dynamic', :select_options => %w(dynamic non-dynamic)
  dep :tf_activity_heatmap_matrix, :scheme => :placeholder, :normalization => 'row_zscore' do |jobname, options|
    { :scheme => options[:scheme], :normalization => 'row_zscore' }
  end
  extension :png
  task :figure_10c_tf_heatmap_columns_clustered => :binary do |scheme|
    matrix = step(:tf_activity_heatmap_matrix).load
    fields = AGS.figure_tf_context_fields(matrix)
    data_tsv = matrix.reorder('Associated Gene Name', fields)
    AGS.figure_ggplot(self.tmp_path, data_tsv, <<-RCODE, 12, 10)
mat <- data[, c(#{fields.collect{|f| "'#{f}'"} * ', '}), drop=FALSE]
mat[] <- lapply(mat, as.numeric)
mat[is.na(mat)] <- 0
mat <- as.matrix(mat)
rownames(mat) <- rownames(data)
if (nrow(mat) > 2) mat <- mat[hclust(dist(mat))$order, , drop=FALSE]
if (ncol(mat) > 2) mat <- mat[, hclust(dist(t(mat)))$order, drop=FALSE]
mat[mat > 3] <- 3
mat[mat < -3] <- -3
df <- as.data.frame(as.table(mat))
names(df) <- c('TF','Context','Activity')
df$TF <- factor(df$TF, levels=rev(unique(df$TF)))
df$Context <- factor(df$Context, levels=unique(df$Context))
ggplot(df, aes(Context, TF, fill=Activity)) + geom_raster() +
  scale_fill_gradient2(low='#4575B4', mid='white', high='#D73027', midpoint=0, limits=c(-3,3), name='row z-score') +
  theme_minimal() + labs(x='', y='', title='#{scheme}: columns clustered') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6), axis.text.y=element_blank(), panel.grid=element_blank())
RCODE
  end

  dep :tf_activity_heatmap_matrix, :scheme => 'dynamic', :normalization => 'row_zscore'
  extension :png
  task :figure_10d_tf_context_similarity_heatmap => :binary do
    matrix = step(:tf_activity_heatmap_matrix).load
    fields = AGS.figure_tf_context_fields(matrix)
    data_tsv = matrix.reorder('Associated Gene Name', fields)
    AGS.figure_ggplot(self.tmp_path, data_tsv, <<-RCODE, 9, 8)
mat <- data[, c(#{fields.collect{|f| "'#{f}'"} * ', '}), drop=FALSE]
mat[] <- lapply(mat, as.numeric)
mat[is.na(mat)] <- 0
cor_mat <- cor(as.matrix(mat), use='pairwise.complete.obs')
ord <- hclust(as.dist(1 - cor_mat))$order
cor_mat <- cor_mat[ord, ord]
df <- as.data.frame(as.table(cor_mat))
names(df) <- c('Context1','Context2','Correlation')
df$Context1 <- factor(df$Context1, levels=colnames(cor_mat))
df$Context2 <- factor(df$Context2, levels=rev(colnames(cor_mat)))
ggplot(df, aes(Context1, Context2, fill=Correlation)) + geom_tile() + coord_fixed() +
  scale_fill_gradient2(low='#4575B4', mid='white', high='#D73027', midpoint=0, limits=c(-1,1), name='Pearson r') +
  theme_minimal() + labs(x='', y='', title='TF activity context similarity') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6), axis.text.y=element_text(size=6), panel.grid=element_blank())
RCODE
  end

  dep :tf_activity_call_counts_by_scheme, :scheme => 'dynamic'
  extension :png
  task :figure_10e_tf_activity_call_burden => :binary do
    tsv = step(:tf_activity_call_counts_by_scheme).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 8, 5)
data <- subset(data, Sign == 'both')
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$TFActivityCalls <- as.numeric(data$TFActivityCalls)
ggplot(data, aes(Time, TreatmentLabel, fill=TFActivityCalls)) + geom_tile(color='white') +
  scale_fill_viridis_c(name='TF activity calls') + theme_bw() + labs(x='', y='') + theme(panel.grid=element_blank())
RCODE
  end

  #############################################################################
  # Figure 11: combination TF activity categories
  #############################################################################

  dep :combination_tf_category_counts
  extension :png
  task :figure_11b_combination_category_counts => :binary do
    tsv = step(:combination_tf_category_counts).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 10, 5)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$CombinationLabel <- factor(data$Combination, levels=c('INT_FiveZ_PI','INT_PD_PI'), labels=c('5Z+PI','PD+PI'))
data$TFActivityCalls <- as.numeric(data$TFActivityCalls)
data$CategoryLabel <- factor(data$Category, levels=c(#{FIGURE_COMBINATION_CATEGORY_LABELS.keys.collect{|c| "'#{c}'"} * ','}), labels=c(#{FIGURE_COMBINATION_CATEGORY_LABELS.values.collect{|c| "'#{c}'"} * ','}))
ggplot(data, aes(Time, TFActivityCalls, fill=CategoryLabel)) + geom_col() + facet_wrap(~CombinationLabel, nrow=1) +
  scale_fill_brewer(palette='Set3') + theme_bw() + labs(x='Time point', y='TF activity calls in combination', fill='Category') +
  theme(legend.position='bottom', legend.text=element_text(size=7))
RCODE
  end

  dep :combination_tf_category_counts
  extension :png
  task :figure_11c_combination_category_heatmap => :binary do
    tsv = step(:combination_tf_category_counts).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 9, 6)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$CombinationLabel <- factor(data$Combination, levels=c('INT_FiveZ_PI','INT_PD_PI'), labels=c('5Z+PI','PD+PI'))
data$TFActivityCalls <- as.numeric(data$TFActivityCalls)
data$CategoryLabel <- factor(data$Category, levels=c(#{FIGURE_COMBINATION_CATEGORY_LABELS.keys.collect{|c| "'#{c}'"} * ','}), labels=c(#{FIGURE_COMBINATION_CATEGORY_LABELS.values.collect{|c| "'#{c}'"} * ','}))
ggplot(data, aes(Time, CategoryLabel, fill=TFActivityCalls)) + geom_tile(color='white') + facet_wrap(~CombinationLabel, nrow=1) +
  scale_fill_viridis_c(name='TF calls') + theme_bw() + labs(x='', y='') + theme(panel.grid=element_blank())
RCODE
  end

  dep :combination_tf_category_counts
  extension :png
  task :figure_11d_combination_category_fraction => :binary do
    tsv = step(:combination_tf_category_counts).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 10, 5)
data$TFActivityCalls <- as.numeric(data$TFActivityCalls)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$CombinationLabel <- factor(data$Combination, levels=c('INT_FiveZ_PI','INT_PD_PI'), labels=c('5Z+PI','PD+PI'))
data$CategoryLabel <- factor(data$Category, levels=c(#{FIGURE_COMBINATION_CATEGORY_LABELS.keys.collect{|c| "'#{c}'"} * ','}), labels=c(#{FIGURE_COMBINATION_CATEGORY_LABELS.values.collect{|c| "'#{c}'"} * ','}))
ggplot(data, aes(Time, TFActivityCalls, fill=CategoryLabel)) + geom_col(position='fill') + facet_wrap(~CombinationLabel, nrow=1) +
  scale_fill_brewer(palette='Set3') + theme_bw() + labs(x='Time point', y='Fraction of combination TF activity calls', fill='Category') +
  theme(legend.position='bottom', legend.text=element_text(size=7))
RCODE
  end

  dep :combination_tf_category_counts
  extension :png
  task :figure_11f_combination_earlier_counts => :binary do
    tsv = step(:combination_tf_category_counts).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 7, 4.5)
data <- subset(data, Category == 'combination_earlier_than_both')
data$TFActivityCalls <- as.numeric(data$TFActivityCalls)
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$CombinationLabel <- factor(data$Combination, levels=c('INT_FiveZ_PI','INT_PD_PI'), labels=c('5Z+PI','PD+PI'))
ggplot(data, aes(Time, TFActivityCalls, fill=CombinationLabel)) + geom_col(position='dodge') +
  scale_fill_manual(values=c('5Z+PI'='#{FIGURE_TREATMENT_COLORS['INT_FiveZ_PI']}','PD+PI'='#{FIGURE_TREATMENT_COLORS['INT_PD_PI']}')) +
  theme_bw() + labs(x='Time point', y='Combination earlier than both components', fill='Combination') + theme(legend.position='bottom')
RCODE
  end

  #############################################################################
  # Figure 12: Neko validation, dynamic versus non-dynamic
  #############################################################################

  dep :neko_dynamic_non_dynamic_summary
  extension :png
  task :figure_12a_neko_match_odds => :binary do
    tsv = step(:neko_dynamic_non_dynamic_summary).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 9, 5)
data <- subset(data, Vetting == 'none')
data$MatchOdds <- as.numeric(data$MatchOdds)
data$Treatment <- factor(data$Treatment)
data$Target <- factor(data$Target, levels=c('T1','T2','relaxed'))
ggplot(data, aes(Target, MatchOdds, color=Scheme, group=Scheme)) + geom_hline(yintercept=1, linetype='dashed', color='grey50') + geom_point(size=2) + geom_line() +
  facet_wrap(~Treatment, scales='free_y') + scale_color_manual(values=c('dynamic'='#D73027','non-dynamic'='#4575B4')) +
  theme_bw() + labs(x='Benchmark evaluation', y='Match odds: matches / non-matches', color='Scheme') + theme(legend.position='bottom')
RCODE
  end

  dep :neko_dynamic_vs_non_dynamic_odds
  extension :png
  task :figure_12b_neko_dynamic_vs_nondynamic_odds => :binary do
    tsv = step(:neko_dynamic_vs_non_dynamic_odds).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 9, 5)
data <- subset(data, Vetting == 'none')
data$DynamicVsNonDynamicOddsRatio <- as.numeric(data$DynamicVsNonDynamicOddsRatio)
data$Target <- factor(data$Target, levels=c('T1','T2','relaxed'))
ggplot(data, aes(Target, DynamicVsNonDynamicOddsRatio, fill=Treatment)) + geom_hline(yintercept=1, linetype='dashed', color='grey50') + geom_col(position='dodge') +
  theme_bw() + labs(x='Benchmark evaluation', y='Dynamic / non-dynamic match odds ratio', fill='Benchmark') + theme(legend.position='bottom')
RCODE
  end

  dep :neko_dynamic_non_dynamic_summary
  extension :png
  task :figure_12c_neko_match_fraction => :binary do
    tsv = step(:neko_dynamic_non_dynamic_summary).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 9, 5)
data <- subset(data, Vetting == 'none')
data$MatchFraction <- as.numeric(data$MatchFraction)
data$Target <- factor(data$Target, levels=c('T1','T2','relaxed'))
ggplot(data, aes(Target, MatchFraction, color=Scheme, group=Scheme)) + geom_point(size=2) + geom_line() +
  facet_wrap(~Treatment) + scale_color_manual(values=c('dynamic'='#D73027','non-dynamic'='#4575B4')) +
  ylim(0,1) + theme_bw() + labs(x='Benchmark evaluation', y='Match fraction', color='Scheme') + theme(legend.position='bottom')
RCODE
  end

  dep :neko_dynamic_non_dynamic_summary
  extension :png
  task :figure_12d_neko_counts => :binary do
    tsv = step(:neko_dynamic_non_dynamic_summary).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 9, 5)
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

  #############################################################################
  # Figure 13: self-consistency validation, dynamic versus non-dynamic
  #############################################################################

  dep :self_consistency_dynamic_non_dynamic_summary
  extension :png
  task :figure_13b_self_consistency_match_odds => :binary do
    tsv = step(:self_consistency_dynamic_non_dynamic_summary).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 10, 5)
data <- subset(data, Vetting == 'none')
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$MatchOdds <- as.numeric(data$MatchOdds)
ggplot(data, aes(Time, TreatmentLabel, fill=MatchOdds)) + geom_tile(color='white') + facet_wrap(~Scheme, nrow=1) +
  scale_fill_viridis_c(name='Match odds', trans='sqrt', na.value='grey90') + theme_bw() + labs(x='', y='') + theme(panel.grid=element_blank())
RCODE
  end

  dep :self_consistency_dynamic_non_dynamic_summary
  extension :png
  task :figure_13c_self_consistency_counts => :binary do
    tsv = step(:self_consistency_dynamic_non_dynamic_summary).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 10, 5)
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

  dep :self_consistency_dynamic_vs_non_dynamic_odds
  extension :png
  task :figure_13d_self_consistency_dynamic_vs_nondynamic_odds => :binary do
    tsv = step(:self_consistency_dynamic_vs_non_dynamic_odds).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 8, 5)
data <- subset(data, Vetting == 'none')
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$DynamicVsNonDynamicOddsRatio <- as.numeric(data$DynamicVsNonDynamicOddsRatio)
ggplot(data, aes(Time, TreatmentLabel, fill=DynamicVsNonDynamicOddsRatio)) + geom_tile(color='white') +
  scale_fill_gradient2(low='#4575B4', mid='white', high='#D73027', midpoint=1, name='Odds ratio', na.value='grey90') +
  theme_bw() + labs(x='', y='') + theme(panel.grid=element_blank())
RCODE
  end

  dep :self_consistency_dynamic_non_dynamic_summary
  extension :png
  task :figure_13e_self_consistency_match_fraction => :binary do
    tsv = step(:self_consistency_dynamic_non_dynamic_summary).load
    AGS.figure_ggplot(self.tmp_path, tsv, <<-RCODE, 10, 5)
data <- subset(data, Vetting == 'none')
data$Time <- factor(data$Time, levels=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h'))
data$TreatmentLabel <- factor(data$Treatment, levels=c(#{FIGURE_TREATMENT_ORDER.collect{|t| "'#{t}'"} * ','}), labels=c(#{AGS.figure_treatment_levels}))
data$MatchFraction <- as.numeric(data$MatchFraction)
ggplot(data, aes(Time, MatchFraction, color=Scheme, group=Scheme)) + geom_line() + geom_point(size=2) +
  facet_wrap(~TreatmentLabel, nrow=2) + scale_color_manual(values=c('dynamic'='#D73027','non-dynamic'='#4575B4')) +
  ylim(0,1) + theme_bw() + labs(x='Time point', y='Match fraction among all TF activity calls', color='Scheme') + theme(legend.position='bottom')
RCODE
  end

  #############################################################################
  # Convenience collectors
  #############################################################################

  dep :figure_01b_concept_intermediate_phenotypes
  dep :figure_01c_concept_dynamic_workflow
  dep :figure_02c_mrna_pca_treatment
  dep :figure_02d_mrna_pca_time
  dep :figure_02e_mrna_treatment_distance
  dep :figure_03b_fc0_de_counts
  dep :figure_03c_fc1_de_counts
  dep :figure_04a_onset_first_counts
  dep :figure_04b_onset_episode_counts
  dep :figure_04c_onset_switch_summary
  dep :figure_04d_onset_example_profiles
  dep :figure_05a_fc1_de_direction_counts
  dep :figure_05b_fc1_onset_relationship
  dep :figure_05e_fc1_onset_8h_focus
  dep :figure_06a_go_top_terms_dotplot
  dep :figure_06b_go_term_frequency
  dep :figure_06c_go_theme_heatmap
  dep :figure_06d_go_selected_terms
  dep :figure_07b_tf_activity_call_counts
  dep :figure_07c_tf_activity_pca
  dep :figure_07d_tf_activity_sign_balance
  dep :figure_08a_tf_report_card_example
  dep :figure_08b_tf_targets_dynamic_highlight
  dep :figure_08c_tf_target_edge_consistency
  dep :figure_09_tf_activity_heatmap, :scheme => 'dynamic', :normalization => 'raw', :max_abs => 50.0
  dep :figure_09_tf_activity_heatmap, :scheme => 'dynamic', :normalization => 'row_zscore', :max_abs => 3.0
  dep :figure_09_tf_activity_heatmap, :scheme => 'dynamic', :normalization => 'column_zscore', :max_abs => 3.0
  dep :figure_09d_tf_context_correlation
  dep :figure_10c_tf_heatmap_columns_clustered, :scheme => 'dynamic'
  dep :figure_10d_tf_context_similarity_heatmap
  dep :figure_10e_tf_activity_call_burden
  dep :figure_11b_combination_category_counts
  dep :figure_11c_combination_category_heatmap
  dep :figure_11d_combination_category_fraction
  dep :figure_11f_combination_earlier_counts
  dep :figure_12a_neko_match_odds
  dep :figure_12b_neko_dynamic_vs_nondynamic_odds
  dep :figure_12c_neko_match_fraction
  dep :figure_12d_neko_counts
  dep :figure_13b_self_consistency_match_odds
  dep :figure_13c_self_consistency_counts
  dep :figure_13d_self_consistency_dynamic_vs_nondynamic_odds
  dep :figure_13e_self_consistency_match_fraction
  task :figures_current => :array do
    dependencies.each do |dep|
      inputs = dep.recursive_inputs
      tags = []
      tags << inputs[:scheme] if inputs[:scheme] && ! inputs[:scheme].to_s.empty?
      tags << inputs[:normalization] if inputs[:normalization] && ! inputs[:normalization].to_s.empty?
      name = ([dep.task_name.to_s] + tags).compact * '-'
      target = file(name + '.png')
      Open.cp dep.path, target
    end
    files
  end

end
