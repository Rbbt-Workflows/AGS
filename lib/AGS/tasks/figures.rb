module AGS

  DISPLAY_TREATMENTS = %w(INT_PD_PI INT_FiveZ_PI PI FiveZ PD DMSO)
  DISPLAY_TIME_POINTS = [1, 2, 4, 8, 24]

  TREATMENT_LABELS = {
    'INT_PD_PI'    => 'PD+PI',
    'INT_FiveZ_PI' => '5Z+PI',
    'PI'           => 'PI',
    'FiveZ'        => '5Z',
    'PD'           => 'PD',
    'DMSO'         => 'DMSO',
  }

  TREATMENT_COLORS = {
    'INT_PD_PI'    => '#E41A1C',
    'INT_FiveZ_PI' => '#FF7F00',
    'PI'           => '#377EB8',
    'FiveZ'        => '#4DAF4A',
    'PD'           => '#984EA3',
    'DMSO'         => '#999999',
  }

  def self.results_dir
    @results_dir ||= File.expand_path(File.join(File.dirname(__FILE__), '..', '..', '..', 'results'))
  end

  def self.load_result(name, opts = {})
    TSV.open(File.join(results_dir, name), opts)
  end

  # Returns treatment-timepoint fields in canonical order, filtered to those present in tsv
  def self.ordered_fields(tsv)
    DISPLAY_TREATMENTS.flat_map do |trt|
      DISPLAY_TIME_POINTS.map { |tp| "#{trt}-T#{tp}" }
    end.select { |f| tsv.fields.include?(f) }
  end

  # Build a :list TSV from rows. Each row is [key, val1, val2, ...]
  def self.build_tsv(rows, key_field, fields)
    tsv = TSV.setup({}, :key_field => key_field, :fields => fields, :type => :list)
    rows.each { |r| k = r.shift; tsv[k] = r }
    tsv
  end

  # R code fragments for ggplot color/fill scales and factor levels
  def self.r_color_scale(treatments = DISPLAY_TREATMENTS)
    vals = treatments.map { |t| "'#{TREATMENT_LABELS[t]}'='#{TREATMENT_COLORS[t]}'" }.join(',')
    "scale_color_manual(values=c(#{vals}))"
  end

  def self.r_fill_scale(treatments = DISPLAY_TREATMENTS)
    vals = treatments.map { |t| "'#{TREATMENT_LABELS[t]}'='#{TREATMENT_COLORS[t]}'" }.join(',')
    "scale_fill_manual(values=c(#{vals}))"
  end

  def self.r_label_levels(treatments = DISPLAY_TREATMENTS)
    treatments.map { |t| "'#{TREATMENT_LABELS[t]}'" }.join(',')
  end

  # Human-readable column name for a treatment-timepoint field (use _ not space for R safety)
  def self.nice_colname(field)
    if field =~ /^(.+)-T(\d+)$/
      trt = $1; tp = $2
      TREATMENT_LABELS.key?(trt) ? "#{TREATMENT_LABELS[trt]}_T#{tp}" : field
    else
      field
    end
  end

  # Save a ggplot figure using the rbbt R integration.
  # +tsv+ is the data TSV; +r_code+ is R code that uses the data frame +data+
  # and evaluates to a ggplot object.
  def self.ggplot_png(path, tsv, r_code, width = 8, height = 4)
    FileUtils.mkdir_p File.dirname(path)
    R::PNG.ggplot path, tsv, r_code, width, height
  end

  # Same but returns binary content for workflow :binary tasks
  def self.ggplot_figure(tsv, file, r_code, width = 8, height = 4)
    ggplot_png(file, tsv, r_code, width, height)
    nil
  end

  #############################################################################
  # Figure R2: DE gene count trajectories across treatments and time points
  #############################################################################

  input :pvalue_threshold, :float, "Adjusted p-value threshold for DE calling", 0.05
  input :fc_threshold, :float, "Minimum |log2FC| threshold (0 = none)", 0
  extension :png
  task :figure_R2_de_gene_counts => :binary do |pvalue_threshold, fc_threshold|
    # Select only the FC and Pvalue columns we need
    fc_fields = DISPLAY_TREATMENTS.flat_map { |trt| DISPLAY_TIME_POINTS.map { |tp| "FC_#{trt}.T#{tp}" } }
    pv_fields = DISPLAY_TREATMENTS.flat_map { |trt| DISPLAY_TIME_POINTS.map { |tp| "Pvalue_#{trt}.T#{tp}" } }
    keep_fields = fc_fields + pv_fields
    info = AGS.load_result('full_gene_info.tsv', :type => :list, :fields => keep_fields)

    rows = []
    DISPLAY_TREATMENTS.each do |trt|
      DISPLAY_TIME_POINTS.each do |tp|
        fc_f   = "FC_#{trt}.T#{tp}"
        pv_f   = "Pvalue_#{trt}.T#{tp}"
        next unless info.fields.include?(fc_f) && info.fields.include?(pv_f)

        fc_col = info.identify_field fc_f
        pv_col = info.identify_field pv_f
        count = 0
        info.each do |gene, values|
          pv_s = values[pv_col]; pv = pv_s.nil? || pv_s.empty? ? 1.0 : pv_s.to_f
          fc_s = values[fc_col]; fc = fc_s.nil? || fc_s.empty? ? 0.0 : fc_s.to_f
          next if pv > pvalue_threshold
          next if fc_threshold > 0 && fc.abs < fc_threshold
          count += 1
        end
        rows << ["#{trt}_#{tp}", trt, tp.to_s, count.to_s, TREATMENT_LABELS[trt]]
      end
    end

    tsv = AGS.build_tsv(rows, "ID", %w(Treatment TimePoint Count Label))
    threshold_note = fc_threshold > 0 ? ", |log2FC| > #{fc_threshold}" : ""
    color_scale = AGS.r_color_scale
    label_levels = AGS.r_label_levels

    AGS.ggplot_figure(tsv, self.tmp_path, <<-RCODE, 8, 4)
data$TimePoint <- as.numeric(data$TimePoint)
data$Count     <- as.numeric(data$Count)
data$Label     <- factor(data$Label, levels=c(#{label_levels}))
#{color_scale}
ggplot(data, aes(x=TimePoint, y=Count, color=Label, group=Label)) +
  geom_line(linewidth=1) + geom_point(size=2) +
  scale_x_continuous(breaks=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h')) +
  theme_bw() +
  labs(x='Time point', y='DE genes (adj. p < #{pvalue_threshold}#{threshold_note})',
       color='Treatment') +
  theme(legend.position='bottom', text=element_text(size=10))
RCODE
    nil
  end

  #############################################################################
  # Figure R3: Dynamic gene sets bar plot
  #############################################################################

  input :min_genes, :integer, "Minimum gene set size to show", 0
  extension :png
  task :figure_R3_dynamic_gene_sets => :binary do |min_genes|
    cluster_fields = DISPLAY_TREATMENTS.map { |trt| "#{trt}: FC clusters" }
    info = AGS.load_result('full_gene_info.tsv', :type => :list, :fields => cluster_fields)

    rows = []
    DISPLAY_TREATMENTS.each do |trt|
      cf = "#{trt}: FC clusters"
      next unless info.fields.include?(cf)
      col = info.identify_field cf
      onset_counts = Hash.new(0)
      info.each do |gene, values|
        label = values[col].to_s
        next if label.nil? || label.empty? || label == 'unclassified'
        label.split('|').each do |part|
          part = part.strip
          if part =~ /^(increase|decrease)\s+(\d+)(?:-(\d+))?h$/
            direction = $1; start_t = $2.to_i
            onset_counts["#{direction} #{start_t}h"] += 1
          end
        end
      end
      onset_counts.each do |onset, count|
        next if count < min_genes
        direction = onset.start_with?('increase') ? 'Increase' : 'Decrease'
        time = onset.sub(/^(increase|decrease)\s+/, '').sub(/h$/, '')
        rows << ["#{trt}_#{onset}", trt, time, direction, count.to_s, TREATMENT_LABELS[trt]]
      end
    end

    tsv = AGS.build_tsv(rows, "ID", %w(Treatment TimePoint Direction Count Label))
    fill_scale = AGS.r_fill_scale
    label_levels = AGS.r_label_levels

    AGS.ggplot_figure(tsv, self.tmp_path, <<-RCODE, 12, 5)
data$TimePoint <- factor(data$TimePoint, levels=c('1','2','4','8','24'))
data$Count     <- as.numeric(data$Count)
data$Label     <- factor(data$Label, levels=c(#{label_levels}))
#{fill_scale}
ggplot(data, aes(x=TimePoint, y=Count, fill=Label)) +
  geom_bar(stat='identity', position=position_dodge2(preserve='single')) +
  facet_wrap(~ Direction, nrow=1) +
  theme_bw() +
  labs(x='Onset interval start', y='Number of genes', fill='Treatment') +
  theme(legend.position='bottom', text=element_text(size=10),
        axis.text.x=element_text(angle=45, hjust=1))
RCODE
  end

  #############################################################################
  # Figure R4: Example gene trajectories for a selected dynamic set
  #############################################################################

  input :example_treatment, :string, "Treatment to show examples for", "PI"
  input :example_onset, :string, "Onset label to highlight", "increase 2h"
  input :max_genes, :integer, "Maximum number of example gene trajectories", 30
  extension :png
  task :figure_R4_example_trajectories => :binary do |example_treatment, example_onset, max_genes|
    trt = example_treatment
    fc_fields = DISPLAY_TIME_POINTS.map { |tp| "FC_#{trt}.T#{tp}" }
    cf = "#{trt}: FC clusters"
    info = AGS.load_result('full_gene_info.tsv', :type => :list, :fields => fc_fields + [cf])

    cf_col = info.identify_field cf
    matching_genes = []
    info.each do |gene, values|
      label = values[cf_col].to_s
      next unless label && label.include?(example_onset)
      matching_genes << gene
    end
    matching_genes = matching_genes.first(max_genes)

    fc_cols = fc_fields.map { |f| info.identify_field f }

    rows = []
    matching_genes.each do |gene|
      values = info[gene]; next unless values
      DISPLAY_TIME_POINTS.each_with_index do |tp, i|
        fc = values[fc_cols[i]]; fc = fc.nil? || fc.empty? ? 0.0 : fc.to_f
        rows << ["#{gene}_#{tp}", gene, tp.to_s, fc.to_s, 'gene']
      end
    end

    # Add mean trajectory
    unless matching_genes.empty?
      DISPLAY_TIME_POINTS.each_with_index do |tp, i|
        vals = matching_genes.map do |gene|
          values = info[gene]; next unless values
          fc = values[fc_cols[i]]; fc.nil? || fc.empty? ? 0.0 : fc.to_f
        end.compact
        mean_fc = vals.empty? ? 0 : vals.sum / vals.length
        rows << ["_MEAN_#{tp}", "__MEAN__", tp.to_s, mean_fc.to_s, 'mean']
      end
    end

    tsv = AGS.build_tsv(rows, "ID", %w(Gene TimePoint FC Type))
    n_genes = matching_genes.length

    AGS.ggplot_figure(tsv, self.tmp_path, <<-RCODE, 7, 5)
data$TimePoint <- as.numeric(data$TimePoint)
data$FC <- as.numeric(data$FC)
ggplot(data, aes(x=TimePoint, y=FC, group=Gene)) +
  geom_line(data=subset(data, Type=='gene'), color='grey70', alpha=0.5) +
  geom_line(data=subset(data, Type=='mean'), color='#E41A1C', linewidth=1.5) +
  geom_point(data=subset(data, Type=='mean'), color='#E41A1C', size=2) +
  geom_hline(yintercept=0, linetype='dashed', color='grey40') +
  scale_x_continuous(breaks=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h')) +
  theme_bw() +
  labs(x='Time point', y='log2 fold change',
       title='#{trt}: #{example_onset} (#{n_genes} genes)',
       subtitle='Red = mean, grey = individual genes') +
  theme(text=element_text(size=10))
RCODE
  end

  #############################################################################
  # Figure R6: TF activity counts by treatment, timepoint, direction
  #############################################################################

  extension :png
  task :figure_R6_tf_activity_counts => :binary do
    preds = AGS.load_result('tf_predictions.tsv', :type => :list)

    rows = []
    DISPLAY_TREATMENTS.each do |trt|
      DISPLAY_TIME_POINTS.each do |tp|
        field = "#{trt}-T#{tp}"
        next unless preds.fields.include?(field)
        col = preds.identify_field field
        up_count = 0; dn_count = 0
        preds.each do |tf, values|
          v = values[col]; v = v.nil? ? 0.0 : v.to_f
          up_count += 1 if v > 0
          dn_count += 1 if v < 0
        end
        rows << ["#{trt}_#{tp}", trt, tp.to_s, up_count.to_s, dn_count.to_s, TREATMENT_LABELS[trt]]
      end
    end

    tsv = AGS.build_tsv(rows, "ID", %w(Treatment TimePoint Up Down Label))
    label_levels = AGS.r_label_levels

    AGS.ggplot_figure(tsv, self.tmp_path, <<-RCODE, 10, 5)
data$TimePoint <- factor(data$TimePoint, levels=c('1','2','4','8','24'))
data$Up   <- as.numeric(data$Up)
data$Down <- -as.numeric(data$Down)
data$Label <- factor(data$Label, levels=c(#{label_levels}))
ggplot(data) +
  geom_bar(aes(x=Label, y=Up,   fill='Activated'), stat='identity') +
  geom_bar(aes(x=Label, y=Down, fill='Repressed'), stat='identity') +
  facet_wrap(~ TimePoint, nrow=1, scales='free_y') +
  scale_fill_manual(values=c('Activated'='#E41A1C', 'Repressed'='#377EB8')) +
  coord_flip() +
  theme_bw() +
  labs(x='', y='Number of TFs (+activated / -repressed)', fill='') +
  theme(legend.position='bottom', text=element_text(size=9))
RCODE
  end

  #############################################################################
  # Figure R7: Opposing early TF activities (scatter plots at a given time point)
  #############################################################################

  input :time_point, :integer, "Time point for comparison", 1
  input :treatment_x, :string, "Treatment on x-axis", "PI"
  input :treatment_y, :string, "Treatment on y-axis", "FiveZ"
  input :label_tfs, :array, "TFs to label in the scatter",
    %w(FOXO3 FOXO1 JUN ELK1 MYC TP53 PPARG CTNNB1 E2F1 STAT3)
  extension :png
  task :figure_R7_tf_scatter => :binary do |time_point, treatment_x, treatment_y, label_tfs|
    preds = AGS.load_result('tf_predictions.tsv', :type => :list)
    field_x = "#{treatment_x}-T#{time_point}"
    field_y = "#{treatment_y}-T#{time_point}"
    raise "Field not found: #{field_x}" unless preds.fields.include?(field_x)
    raise "Field not found: #{field_y}" unless preds.fields.include?(field_y)
    col_x = preds.identify_field field_x
    col_y = preds.identify_field field_y

    rows = []
    preds.each do |tf, values|
      vx = values[col_x]; vx = vx.nil? ? 0.0 : vx.to_f
      vy = values[col_y]; vy = vy.nil? ? 0.0 : vy.to_f
      next if vx.abs < 0.01 && vy.abs < 0.01
      lbl = label_tfs.include?(tf) ? tf : ""
      rows << [tf, vx.to_s, vy.to_s, lbl]
    end

    tsv = AGS.build_tsv(rows, "TF", %w(ActivityX ActivityY LabelTF))
    label_x = TREATMENT_LABELS[treatment_x]
    label_y = TREATMENT_LABELS[treatment_y]

    AGS.ggplot_figure(tsv, self.tmp_path, <<-RCODE, 6, 6)
library(ggrepel)
data$ActivityX <- as.numeric(data$ActivityX)
data$ActivityY <- as.numeric(data$ActivityY)
r_val <- tryCatch(round(cor(data$ActivityX, data$ActivityY), 2), error=function(e) NA)
ggplot(data, aes(x=ActivityX, y=ActivityY)) +
  geom_point(alpha=0.4, size=1.5) +
  geom_hline(yintercept=0, linetype='dashed', color='grey50') +
  geom_vline(xintercept=0, linetype='dashed', color='grey50') +
  geom_text_repel(aes(label=LabelTF), size=3, max.overlaps=20,
    box.padding=0.5, segment.color='grey60') +
  annotate('text', x=Inf, y=Inf, hjust=1.1, vjust=1.5,
    label=paste('r =', r_val), size=4) +
  theme_bw() +
  labs(x='#{label_x} T#{time_point} activity',
       y='#{label_y} T#{time_point} activity') +
  theme(text=element_text(size=10))
RCODE
  end

  #############################################################################
  # Figure R8: Correlation heatmap of TF activities across all treatment-timepoints
  #############################################################################

  extension :png
  task :figure_R8_correlation_heatmap => :binary do
    preds = AGS.load_result('tf_predictions.tsv', :type => :list)
    ordered = AGS.ordered_fields(preds)
    # Slice to ordered columns only
    data_tsv = preds.reorder("Associated Gene Name", ordered)

    # Rename columns to human-readable names (safe for R: use underscores)
    nice_names = ordered.map { |f| AGS.nice_colname(f) }
    data_tsv.fields = nice_names

    # Build R code that creates the correlation matrix and plots it
    col_r_vector = nice_names.map { |n| "'#{n}'" }.join(',')

    AGS.ggplot_figure(data_tsv, self.tmp_path, <<-RCODE, 14, 12)
mat_data <- data[, c(#{col_r_vector})]
mat_data[] <- lapply(mat_data, as.numeric)
mat_data[is.na(mat_data)] <- 0
mat <- as.matrix(mat_data)
cor_mat <- cor(mat, use='pairwise.complete.obs')
cor_df <- as.data.frame(as.table(cor_mat))
names(cor_df) <- c('Var1', 'Var2', 'Correlation')
cor_df$Var1 <- factor(cor_df$Var1, levels=unique(cor_df$Var1))
cor_df$Var2 <- factor(cor_df$Var2, levels=rev(unique(cor_df$Var1)))
ggplot(cor_df, aes(x=Var1, y=Var2, fill=Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low='#377EB8', mid='white', high='#E41A1C', midpoint=0,
    limits=c(-1,1), name='Pearson r') +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6),
        axis.text.y=element_text(size=6),
        axis.title=element_blank(), panel.grid=element_blank()) +
  coord_fixed()
RCODE
  end

  #############################################################################
  # Figure R9: Dynamic vs fc0 scheme comparison for selected TFs
  #############################################################################

  input :compare_tfs, :array, "TFs to compare", %w(FOXO3 MYC ELK1 TP53)
  input :selected_treatment, :string, "Treatment to show", "PI"
  extension :png
  task :figure_R9_scheme_comparison => :binary do |compare_tfs, selected_treatment|
    trt = selected_treatment
    preds_dyn = AGS.load_result('tf_predictions.tsv', :type => :list)
    preds_fc0 = AGS.load_result('tf_predictions-fc0.tsv', :type => :list)

    rows = []
    compare_tfs.each do |tf|
      next unless preds_dyn[tf] && preds_fc0[tf]
      DISPLAY_TIME_POINTS.each do |tp|
        fd = "#{trt}-T#{tp}"
        next unless preds_dyn.fields.include?(fd) && preds_fc0.fields.include?(fd)
        vd = preds_dyn[tf][preds_dyn.identify_field(fd)]
        vf = preds_fc0[tf][preds_fc0.identify_field(fd)]
        vd = vd.nil? ? 0.0 : vd.to_f
        vf = vf.nil? ? 0.0 : vf.to_f
        rows << ["#{tf}_dyn_#{tp}", tf, 'dynamic', tp.to_s, vd.to_s]
        rows << ["#{tf}_fc0_#{tp}", tf, 'fc0', tp.to_s, vf.to_s]
      end
    end

    tsv = AGS.build_tsv(rows, "ID", %w(TF Scheme TimePoint Activity))
    ncol = [compare_tfs.length, 4].min
    tf_levels = compare_tfs.join("','")
    trt_label = TREATMENT_LABELS[trt]
    width = ncol * 3.5
    height = (compare_tfs.length.to_f / ncol).ceil * 2.5 + 1

    AGS.ggplot_figure(tsv, self.tmp_path, <<-RCODE, width, height)
data$TimePoint <- as.numeric(data$TimePoint)
data$Activity  <- as.numeric(data$Activity)
data$TF <- factor(data$TF, levels=c('#{tf_levels}'))
ggplot(data, aes(x=TimePoint, y=Activity, color=Scheme, group=Scheme)) +
  geom_line(linewidth=0.8) + geom_point(size=1.5) +
  geom_hline(yintercept=0, linetype='dashed', color='grey50') +
  scale_x_continuous(breaks=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h')) +
  scale_color_manual(values=c('dynamic'='#E41A1C', 'fc0'='#377EB8')) +
  facet_wrap(~ TF, scales='free_y', ncol=#{ncol}) +
  theme_bw() +
  labs(x='', y='Inferred activity', color='Scheme',
       subtitle='Treatment: #{trt_label}') +
  theme(legend.position='bottom', text=element_text(size=9),
        strip.text=element_text(face='bold'))
RCODE
  end

  #############################################################################
  # Figure R10: Full TF activity heatmap with hierarchical clustering
  #############################################################################

  input :scheme, :select, "TF prediction scheme", "dynamic",
    :select_options => %w(dynamic non-dynamic fc0 diff)
  input :show_rownames, :boolean, "Show TF names on rows", false
  input :max_abs, :float, "Cap absolute activity for coloring", 50
  extension :png
  task :figure_R10_tf_heatmap => :binary do |scheme, show_rownames, max_abs|
    fname = scheme == 'dynamic' ? 'tf_predictions.tsv' : "tf_predictions-#{scheme}.tsv"
    preds = AGS.load_result(fname, :type => :list)
    ordered = AGS.ordered_fields(preds)
    data_tsv = preds.reorder("Associated Gene Name", ordered)

    # Filter to TFs with at least one non-zero value
    data_tsv = data_tsv.select do |tf, values|
      values.any? { |v| !v.nil? && v.to_f.abs > 0.01 }
    end

    nice_names = ordered.map { |f| AGS.nice_colname(f) }
    data_tsv.fields = nice_names

    cap = max_abs.to_f
    show_labels = show_rownames ? 'element_text(size=4)' : 'element_blank()'
    n_tfs = data_tsv.keys.length
    height = [8 + n_tfs * 0.06, 30].min

    col_r_vector = nice_names.map { |n| "'#{n}'" }.join(',')

    AGS.ggplot_figure(data_tsv, self.tmp_path, <<-RCODE, 14, height)
mat_data <- data[, c(#{col_r_vector})]
mat_data[] <- lapply(mat_data, as.numeric)
mat_data[is.na(mat_data)] <- 0
mat <- as.matrix(mat_data)
rownames(mat) <- rownames(data)
mat[mat >  #{cap}] <-  #{cap}
mat[mat < -#{cap}] <- -#{cap}
if (nrow(mat) > 2) mat <- mat[hclust(dist(mat))$order, ]
mat_df <- as.data.frame(as.table(mat))
names(mat_df) <- c('TF', 'Condition', 'Activity')
mat_df$Condition <- factor(mat_df$Condition, levels=unique(mat_df$Condition))
mat_df$TF <- factor(mat_df$TF, levels=rev(unique(mat_df$TF)))
ggplot(mat_df, aes(x=Condition, y=TF, fill=Activity)) +
  geom_raster() +
  scale_fill_gradient2(low='#377EB8', mid='white', high='#E41A1C', midpoint=0,
    limits=c(-#{cap}, #{cap}), name='Activity') +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7),
        axis.text.y=#{show_labels}, axis.ticks.y=element_blank(),
        axis.title=element_blank(), panel.grid=element_blank())
RCODE
  end

  #############################################################################
  # Figure R12: TF activity trajectory plots for key TFs across treatments
  #############################################################################

  input :highlight_tfs, :array, "TFs to highlight",
    %w(FOXO3 FOXO1 ELK1 ETS1 MYC TP53 E2F1 FOXM1 STAT3 RELA)
  input :scheme, :select, "TF prediction scheme", "dynamic",
    :select_options => %w(dynamic non-dynamic fc0 diff)
  extension :png
  task :figure_R12_tf_trajectories => :binary do |highlight_tfs, scheme|
    fname = scheme == 'dynamic' ? 'tf_predictions.tsv' : "tf_predictions-#{scheme}.tsv"
    preds = AGS.load_result(fname, :type => :list)
    selected = preds.keys.select { |tf| highlight_tfs.include?(tf) }
    sub = TSV.setup({}, :key_field => preds.key_field, :fields => preds.fields, :type => :list)
    selected.each { |tf| sub[tf] = preds[tf] }
    preds = sub
    ordered = AGS.ordered_fields(preds)




    rows = []
    preds.each do |tf, all_values|
      ordered.each do |field|
        col = preds.identify_field field
        v = all_values[col]; v = v.nil? ? 0.0 : v.to_f
        trt = field.split("-T")[0]; tp = field.split("-T")[1]
        rows << ["#{tf}_#{field}", tf, trt, tp, v.to_s, TREATMENT_LABELS[trt]]
      end
    end

    tsv = AGS.build_tsv(rows, "ID", %w(TF Treatment TimePoint Activity Label))
    ncol = [highlight_tfs.length, 5].min
    width = ncol * 3
    height = (highlight_tfs.length.to_f / ncol).ceil * 2.5 + 1.5
    tf_levels = highlight_tfs.join("','")
    color_scale = AGS.r_color_scale
    label_levels = AGS.r_label_levels

    AGS.ggplot_figure(tsv, self.tmp_path, <<-RCODE, width, height)
data$TimePoint <- as.numeric(data$TimePoint)
data$Activity  <- as.numeric(data$Activity)
data$Label <- factor(data$Label, levels=c(#{label_levels}))
data$TF    <- factor(data$TF, levels=c('#{tf_levels}'))
#{color_scale}
ggplot(data, aes(x=TimePoint, y=Activity, color=Label, group=Label)) +
  geom_line(linewidth=0.8) + geom_point(size=1.5) +
  geom_hline(yintercept=0, linetype='dashed', color='grey50') +
  scale_x_continuous(breaks=c(1,2,4,8,24), labels=c('1h','2h','4h','8h','24h')) +
  facet_wrap(~ TF, scales='free_y', ncol=#{ncol}) +
  theme_bw() +
  labs(x='', y='Inferred activity', color='Treatment') +
  theme(legend.position='bottom', text=element_text(size=9),
        strip.text=element_text(face='bold'))
RCODE
  end

  #############################################################################
  # Figure R13: Self-consistency of TF activities over time
  #############################################################################

  extension :png
  task :figure_R13_self_consistency => :binary do
    counts = AGS.load_result('consistency_counts.tsv', :type => :list)

    rows = []
    counts.each do |id, values|
      treatment = values[0].to_s
      time      = values[1].to_s.to_i
      scheme    = values[2].to_s
      vetting   = values[3].to_s
      matches   = values[4].to_s.to_i
      miss      = values[5].to_s.to_i

      next unless scheme == 'dynamic' && vetting == 'none'
      next unless DISPLAY_TREATMENTS.include?(treatment)
      rows << ["#{treatment}_#{time}", treatment, time.to_s, matches.to_s, miss.to_s,
               TREATMENT_LABELS[treatment]]
    end

    tsv = AGS.build_tsv(rows, "ID", %w(Treatment Time Matches Miss Label))
    label_levels = AGS.r_label_levels

    AGS.ggplot_figure(tsv, self.tmp_path, <<-RCODE, 9, 4)
data$Time    <- factor(data$Time, levels=c('1','2','4','8','24'))
data$Matches <- as.numeric(data$Matches)
data$Miss    <- as.numeric(data$Miss)
data$Label   <- factor(data$Label, levels=c(#{label_levels}))
ggplot(data, aes(x=Time)) +
  geom_bar(aes(y=Matches, fill='Match'), stat='identity') +
  geom_bar(aes(y=-Miss,    fill='Miss'),  stat='identity') +
  facet_wrap(~ Label, nrow=2) +
  scale_fill_manual(values=c('Match'='#4DAF4A', 'Miss'='#E41A1C')) +
  theme_bw() +
  labs(x='Time point', y='Count (+matches / -misses)', fill='Consistency') +
  theme(legend.position='bottom', text=element_text(size=9))
RCODE
  end

  #############################################################################
  # Figure R14: External benchmark agreement (Neko)
  #############################################################################

  extension :png
  task :figure_R14_neko_benchmark => :binary do
    neko_treatments = %w(PI FiveZ PD)
    rows = []
    neko_treatments.each do |trt|
      neko = AGS.load_result("neko_bootstrap_consistency-#{trt}.tsv", :type => :list)
      next if neko.nil?
      t1_match = 0; t2_match = 0; t1_total = 0; t2_total = 0
      neko.each do |tf, values|
        t1_a = values[1].to_s.to_f rescue 0.0
        t2_a = values[2].to_s.to_f rescue 0.0
        t1_m = values[3].to_s
        t2_m = values[4].to_s
        t1_total += 1 if t1_a.abs > 0.01
        t1_match += 1 if t1_m == 'true' && t1_a.abs > 0.01
        t2_total += 1 if t2_a.abs > 0.01
        t2_match += 1 if t2_m == 'true' && t2_a.abs > 0.01
      end
      t1_rate = t1_total > 0 ? (t1_match.to_f / t1_total * 100).round(1) : 0
      t2_rate = t2_total > 0 ? (t2_match.to_f / t2_total * 100).round(1) : 0
      rows << ["#{trt}_T1", trt, 'T1', t1_match.to_s, (t1_total - t1_match).to_s,
               t1_total.to_s, t1_rate.to_s, TREATMENT_LABELS[trt]]
      rows << ["#{trt}_T2", trt, 'T2', t2_match.to_s, (t2_total - t2_match).to_s,
               t2_total.to_s, t2_rate.to_s, TREATMENT_LABELS[trt]]
    end

    tsv = AGS.build_tsv(rows, "ID", %w(Treatment TP Match Miss Total Rate Label))
    fill_scale = AGS.r_fill_scale(neko_treatments)
    neko_levels = AGS.r_label_levels(neko_treatments)

    AGS.ggplot_figure(tsv, self.tmp_path, <<-RCODE, 6, 4)
data$Match <- as.numeric(data$Match)
data$Miss  <- as.numeric(data$Miss)
data$Rate  <- as.numeric(data$Rate)
data$Label <- factor(data$Label, levels=c(#{neko_levels}))
#{fill_scale}
ggplot(data, aes(x=TP, y=Rate, fill=Label)) +
  geom_bar(stat='identity', position=position_dodge()) +
  geom_text(aes(label=paste0(Rate, '%')), vjust=-0.5, size=3,
    position=position_dodge(width=0.9)) +
  ylim(0, 100) +
  theme_bw() +
  labs(x='Time point', y='Agreement with Neko (%)', fill='Treatment') +
  theme(legend.position='bottom', text=element_text(size=10))
RCODE
  end

  #############################################################################
  # Figure R15: TF-TF sequence edge summary by temporal offset
  #############################################################################

  extension :png
  task :figure_R15_sequence_summary => :binary do
    hc_counts = Hash.new { |h, k| h[k] = Hash.new(0) }
    total_counts = Hash.new { |h, k| h[k] = Hash.new(0) }

    DISPLAY_TREATMENTS.each do |trt|
      seq = AGS.load_result("sequence_with_changes-#{trt}.tsv", :type => :list)
      next if seq.nil?
      offset_col = seq.identify_field 'Offset'
      src_sc_col = seq.identify_field 'Source self-consistent'
      tgt_sc_col = seq.identify_field 'Target self-consistent'
      next unless offset_col

      seq.each do |id, values|
        offset = values[offset_col].to_s
        next if offset.nil? || offset.empty?
        total_counts[trt][offset] += 1
        if src_sc_col && tgt_sc_col
          s = values[src_sc_col].to_s
          t = values[tgt_sc_col].to_s
          hc_counts[trt][offset] += 1 if s == 'true' && t == 'true'
        end
      end
    end

    rows = []
    DISPLAY_TREATMENTS.each do |trt|
      %w(0 1 2).each do |offset|
        total = total_counts[trt][offset] || 0
        hc    = hc_counts[trt][offset] || 0
        next if total == 0
        ol = case offset; when '0' then 'Coincident'; when '1' then 'Delay 1'; when '2' then 'Delay 2'; end
        rows << ["#{trt}_#{offset}_all", trt, ol, 'All', total.to_s, TREATMENT_LABELS[trt]]
        rows << ["#{trt}_#{offset}_hc",  trt, ol, 'High conf.', hc.to_s, TREATMENT_LABELS[trt]]
      end
    end

    tsv = AGS.build_tsv(rows, "ID", %w(Treatment Offset Type Count Label))
    fill_scale = AGS.r_fill_scale
    label_levels = AGS.r_label_levels

    AGS.ggplot_figure(tsv, self.tmp_path, <<-RCODE, 8, 4)
data$Count <- as.numeric(data$Count)
data$Label  <- factor(data$Label, levels=c(#{label_levels}))
data$Offset <- factor(data$Offset, levels=c('Coincident', 'Delay 1', 'Delay 2'))
#{fill_scale}
ggplot(subset(data, Type=='High conf.'), aes(x=Offset, y=Count, fill=Label)) +
  geom_bar(stat='identity', position=position_dodge()) +
  theme_bw() +
  labs(x='Temporal offset', y='High-confidence TF-TF edges', fill='Treatment') +
  theme(legend.position='bottom', text=element_text(size=10))
RCODE
  end

  #############################################################################
  # Convenience: run all figures and collect outputs into files_dir
  #############################################################################

  dep :figure_R2_de_gene_counts
  dep :figure_R3_dynamic_gene_sets
  dep :figure_R4_example_trajectories
  dep :figure_R6_tf_activity_counts
  dep :figure_R7_tf_scatter
  dep :figure_R7_tf_scatter, treatment_x: "PI", treatment_y: "PD"
  dep :figure_R7_tf_scatter, treatment_x: "FiveZ", treatment_y: "PD"
  dep :figure_R8_correlation_heatmap
  dep :figure_R9_scheme_comparison
  dep :figure_R10_tf_heatmap
  dep :figure_R12_tf_trajectories
  dep :figure_R13_self_consistency
  dep :figure_R14_neko_benchmark
  dep :figure_R15_sequence_summary
  task :all_figures => :array do
    dependencies.each do |dep|
      target = file(dep.task_name.to_s + ".png")
      Open.ln_h dep.path, target
    end
    dependencies.collect { |d| d.task_name.to_s }.uniq
  end

end
