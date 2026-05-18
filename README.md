# AGS workflow

Time-resolved transcriptional and signaling analysis of AGS gastric cancer cells treated with kinase inhibitors.

## Overview

This workflow implements the full analysis used to interpret a time course experiment in the gastric cancer cell line **AGS** exposed to kinase inhibitors. Cells were treated with:

- PI: PIK3CA inhibitor
- PD: MAP2K2/MAP2K1 inhibitor
- FiveZ: MAP3K7 inhibitor
- INT_PD_PI: PI + PD at half doses
- INT_FiveZ_PI: PI + FiveZ at half doses
- DMSO: vehicle control

Transcriptome measurements were collected at **1, 2, 4, 8, and 24 hours** after treatment. The workflow is designed to explain how cells initially respond and then either adapt or succumb to treatment, in terms of:

- Differential gene expression and its temporal structure
- Transcription factor (TF) activities inferred with `decoupler`
- RNA synthesis/degradation dynamics inferred with INSPEcT
- Regulatory networks (regulome) and downstream TF–TF causal chains
- Enrichment of biological processes (GO:BP) and benchmarking against external TF predictions (NeKo/CollecTRI)

Conceptually, the pipeline:

1. **Prepares expression data**: loads count data and differential expression results from two sources (BSC and NTNU), and normalizes them into fold-changes and p-values across all treatments and time points.
2. **Segments gene expression dynamics**: applies heuristic rules to fold-change trajectories to identify intervals of increasing or decreasing expression for each gene and treatment.
3. **Integrates INSPEcT rates**: classifies genes into synthesis-driven vs degradation-driven regulation and assembles a vetted gene set where expression clusters agree with INSPEcT kinetics.
4. **Builds a regulome and filters it**: obtains TF → target regulatory edges from ExTRI2 or SaezLab, optionally restricting to DNA-binding TFs and expressed coding genes.
5. **Infers TF activities**: runs decoupler on different data matrices (per-timepoint and per-treatment time series) to infer TF activity profiles, optionally restricted to vetted, dynamic genes.
6. **Checks temporal consistency and causal chains**: compares TF activity profiles with gene-level change offsets, evaluates consistency, and derives sequences of TF–TF activation events.
7. **Benchmarks and enrichment**: compares inferred TF activities to external NeKo predictions, explores rule variants, and performs GO/TF enrichments.
8. **Exports results**: generates heatmaps, volcano plots, Excel files and helper lists for downstream manual inspection.

All tasks live in the `AGS` workflow module defined in `workflow.rb` and in `lib/AGS/tasks/*.rb`.

## Usage

You can run this workflow either via the Rbbt command line or programmatically from Ruby. Below are minimal examples; adjust paths and options to your environment.

### Command line (Rbbt)

```bash
# Change to the workflow directory
cd /bulk/mvazque2/git/workflows/AGS

# List tasks
rbbt workflow AGS

# Compute fold-changes from NTNU differential expression
rbbt workflow AGS fold_changes --fc_source=NTNU > fold_changes.tsv

# Infer TF activities for PI at 4 h using dynamic scheme
rbbt workflow AGS treatment_tfs \
  --treatment=PI \
  --scheme=dynamic \
  > PI_dynamic_tfs.tsv

# Run INSPEcT integration and build vetted gene table
rbbt workflow AGS vetted_genes > vetted_genes.tsv

# Build per-timepoint decoupler suite and heatmaps
rbbt workflow AGS timepoint_decoupler_suite > decoupler_suite.tsv
rbbt workflow AGS timepoint_heatmaps
```

### Ruby API

```ruby
require 'rbbt/workflow'
require './workflow'   # loads AGS

# 1. Fold-changes for NTNU data
fc_job = AGS.job(:fold_changes, 'NTNU_default', fc_source: 'NTNU')
fc_job.run
fc = fc_job.path.tsv   # key: gene, fields: treatment.Ttime

# 2. Vetted dynamic genes (INSPEcT + fold-change clusters)
vg_job = AGS.job(:vetted_genes, 'vetted_default')
vg_job.run
vetted = vg_job.path.tsv

# 3. TF activities for PI at 4 h, dynamic scheme
pi_tfs_job = AGS.job(:treatment_tfs, 'PI_dynamic_TFs',
  treatment: 'PI',
  scheme: 'dynamic',
  vetting: 'none'
)
pi_tfs_job.run
pi_tfs = pi_tfs_job.path.tsv

# 4. Time-resolved TF–TF sequence for a treatment
seq_job = AGS.job(:sequence_with_changes, 'PI_sequence', treatment: 'PI')
seq_job.run
seq = seq_job.path.tsv
```

The following section documents each task available in the AGS workflow.

# Tasks

## fold_changes
Compute fold-changes for all genes across treatments and time points, combining BSC or NTNU differential expression sources.

This is the central task for obtaining per-gene log2 fold-changes under different experimental conditions. It wraps either the `log2foldchanges` task from the BSC pipeline or the `fold_changes_NTNU` task and applies optional filtering on direction and conditions.

Inputs include:
- `fc_source` (`BSC` or `NTNU`, default `NTNU`): which upstream DE pipeline to use.
- `subset_effect` (`both`, `up`, `down`): restrict fold-changes by sign (values outside the chosen direction are set to 0).
- `conditions` (array of field names): optional subset of experimental columns to keep.

The output is a TSV with key field as gene identifier (Associated Gene Name) and one column per `TREATMENT.Ttime` combination.

Example (NTNU fold-changes restricted to up-regulation):

```bash
rbbt workflow AGS fold_changes \
  --fc_source=NTNU \
  --subset_effect=up \
  > fold_changes_up.tsv
```

## pvalues
Compute p-values for all genes across treatments and time points, combining BSC or NTNU differential expression sources.

This task parallels `fold_changes` but returns p-values instead of fold-changes. For BSC, it uses the `pvalues_BSC` task; for NTNU, it uses `pvalues_NTNU`.

Inputs:
- `fc_source` (`BSC` or `NTNU`, default `NTNU`).
- `subset_effect` (`both`, `up`, `down`): if set, p-values not matching the chosen direction are zeroed.
- `conditions` (array): optional subset of experimental columns.

The output is a TSV keyed by gene, with one p-value column per experimental condition.

## decoupler_pre
Run decoupler on overall fold-change matrix using the SaezLab regulome (internal helper).

This `dep_task` delegates to `SaezLab.decoupler` using the `fold_changes` matrix and the `regulome` network as inputs. It produces an intermediate TSV with decoupler scores and p-values for each TF.

You rarely need to call this task directly; instead, use `decoupler`, which thresholds the p-values and zeros out non-significant activities.

## decoupler
Apply p-value filtering to the global decoupler results over the fold-change matrix.

This task takes the output of `decoupler_pre`, transposes it to a gene-centric view, and sets activity values to zero wherever the associated p-value is above a user-defined threshold.

Inputs:
- `threshold` (float, default `0.5`): p-value cutoff.

The output is a TSV of filtered decoupler scores (one column per experimental condition, one row per TF).

## fold_changes_fc_one
Compute one-step differences in fold-change between consecutive time points (ΔFC per interval).

Given the NTNU-derived `fold_changes`, this task computes, for each treatment and time point, the difference between the current fold-change and the previous one, effectively capturing the incremental change between time points.

The output is a TSV where keys are `TREATMENT.Ttime` columns and values are ΔFC vs the previous time point. The 1 h time point has no previous point and is skipped.

## fold_changes_fc_one_transpose
Transpose the one-step fold-change matrix to be keyed by gene.

This task simply transposes the result from `fold_changes_fc_one` so that the key field is `Associated Gene Name` and columns correspond to `TREATMENT.Ttime` intervals.

Useful as an input matrix for decoupler analyses that focus on incremental changes rather than absolute fold-changes.

## decoupler_pre_fc_one
Run decoupler on the one-step fold-change matrix (internal helper).

Analogous to `decoupler_pre`, but uses `fold_changes_fc_one` as the input matrix. It calls `SaezLab.decoupler` with the one-step FC matrix and the regulome.

Use `decoupler_fc_one` for a thresholded, user-facing result.

## decoupler_fc_one
Apply p-value filtering to decoupler results over one-step fold-changes.

This task filters the output of `decoupler_pre_fc_one` by p-value, zeroing non-significant entries.

Inputs:
- `threshold` (float, default `0.5`): p-value cutoff.

The output is a TSV of TF activities corresponding to incremental FC changes rather than absolute levels.

## fold_changes_NTNU
Load NTNU differential expression fold-changes and map them to gene symbols.

This task reads the NTNU fold-change file, converts Ensembl gene IDs to `Associated Gene Name`, optionally subsets to a given list of genes, removes empty entries, and transposes to have experiments as keys.

Inputs:
- `subset` (array of gene names): optional subset of genes to keep.

The output is a TSV with key field as experiment identifiers (e.g. `FC_PD.T1`) and fields as genes.

## fold_changes_NTNU_tranposed
Transpose NTNU fold-changes to be gene-centric.

This is a convenience wrapper around `fold_changes_NTNU` that returns a TSV keyed by `Associated Gene Name`, with one column per experimental condition.

## pvalues_NTNU
Load NTNU p-values, convert IDs and align them with the fold-change experiments.

This task parses the NTNU p-value file, maps internal numeric IDs back to genes using the fold-change file, converts `NA` values to `nil`, remaps Ensembl IDs to `Associated Gene Name`, and transposes the result.

Inputs:
- `subset` (array of gene names): optional subset of genes to keep.

The output is a TSV of p-values per experiment, keyed by experiment.

## pvalues_NTNU_transposed
Transpose NTNU p-values to be gene-centric.

This is the gene-centric counterpart of `pvalues_NTNU`, returning a TSV keyed by `Associated Gene Name`.

## gene_counts_pre
Load raw gene count matrix and drop problematic samples.

This task parses `gene_counts.tsv` (Ensembl gene IDs by sample), sets `"Ensembl Gene ID"` as key, and removes hard-coded sample columns `63` and `101`.

The result is a TSV used by downstream tasks to aggregate counts by treatment and timepoint.

## samples
Load the sample metadata table.

This task parses `sample_info.tsv`, using the sample code as key and all other columns (e.g. treatment, timepoint) as fields.

The result is used to map each sample to its treatment/timepoint combination.

## gene_counts
Aggregate per-sample counts into treatment–timepoint groups and average replicates.

Using `gene_counts_pre` and `samples`, this task constructs a matrix keyed by `Ensembl Gene ID`, with fields representing `Treatment code` (e.g. `PI-T1`, `PD-T8`). For each group of replicates, it aggregates counts into a list per gene.

The output is a TSV suitable for downstream differential analysis (`differential`).

## differential
Run differential expression between two aggregated samples using RbbtMatrix.

Given two named treatment–timepoint sample groups (`case_sample` and `control_sample`), this task unpacks replicate columns, merges them, writes a temporary matrix file, and calls `RbbtMatrix#differential`.

Inputs:
- `case_sample` (string): e.g. `PI-T2`.
- `control_sample` (string): e.g. `PI-T0` or `untreated-T0`.

The task writes out DE results to the job’s temp directory via `RbbtMatrix` and returns `nil` as the TSV output is created by the downstream matrix methods.

## all_differential
Aggregate all pairwise differential analyses across treatments and timepoints.

This task depends on multiple `differential` jobs covering all allowed case–control pairs (respecting temporal order). It collects log2 fold-changes and p-values, attaches them column-wise, and converts to `Associated Gene Name` keys.

The resulting TSV has one or more columns per (case,control) comparison, with paired log2FC and p-value columns.

## log2foldchanges
Extract only log2 fold-change columns from `all_differential`.

This task filters out all p-value columns from `all_differential`, returning only log2 fold-change columns.

Useful when you want a dense fold-change matrix without p-values.

## pvalues_BSC
Extract only p-value columns from `all_differential`.

This task filters log2 fold-change columns out of `all_differential`, keeps only p-value columns, and renames them to drop the ` (pvalues)` suffix.

The output aligns with the BSC differential experiments.

## significant_log2foldchanges
Zero out log2 fold-changes that do not pass a p-value threshold.

This task re-parses the `all_differential` TSV and, for each experiment, sets log2 fold-change values to `nil` unless the corresponding p-value is below a user-provided threshold.

Input:
- `threshold` (float, default `0.1`): p-value cutoff for significance.

The output is a TSV of filtered log2 fold-changes.

## significant_foldchanges
Convert significant log2 fold-changes to fold-changes.

This task exponentiates base-2 the values from `significant_log2foldchanges` to obtain fold-changes, leaving `nil` where no significant change is present.

## dbTFs
Return the list of DNA-binding transcription factors (dbTFs) from the ExTRI2 regulome.

This task loads the `EXTRI2-regulome_dbTFs_coTFs_230425` table and selects entries annotated as dbTFs.

The output is an array of TF gene symbols usable as a filter in other tasks.

## vulcano_plot
Draw a volcano plot of fold-change vs p-value for selected timepoints.

This task uses NTNU `fold_changes` and `pvalues`, focusing on selected timepoints (by default PD_PI conditions in the code) and generates a volcano plot using R’s `ggplot2`.

The output is a PNG image (binary) of the volcano plot.

## change_offsets
Classify genes into temporal intervals of increasing or decreasing expression based on fold-change dynamics.

This is a core task for temporal segmentation. It:

- Transposes NTNU `fold_changes` and `pvalues` to be gene-centric.
- Optionally restricts to protein-coding genes using Ensembl biotype annotations.
- For each gene and treatment, computes first-step changes between consecutive timepoints and uses rule parameters (`low_min_Xh`, `mid_min_Xh`, `high_min_Xh`, `next_min_Xh`, `low_multiplier_Xh`, `mid_multiplier_Xh`) to decide whether a given time point initiates an **increase** or **decrease** interval.
- Collapses consecutive indices into readable intervals like `increase 1-4h` or `decrease 2-8h`.
- Optionally extends intervals based on a fold-change p-value threshold (`fc_one_threshold`).

Inputs include many per-timepoint tuning parameters; defaults in the code capture a calibrated rule set.

The output is a TSV keyed by `Associated Gene Name`, with one field per experiment indicating a list of intervals (e.g. `["increase 1-4h", "decrease 8-24h"]`).

## change_offsets_simplified
Simplify change intervals to single timepoints.

This task post-processes `change_offsets` by truncating the end of intervals (e.g. `increase 1-4h` → `increase 1h`), yielding a simpler representation of "change onset" time.

The output is structurally identical to `change_offsets`, but each interval ends at its start time.

## INSPEcT
Merge INSPEcT RNA rate files across treatments into a single table.

This task scans the INSPEcT data directory, loading for each treatment the `*_data.csv` and `*_gene_info.csv` files, aligning them by `Ensembl Gene ID`, labeling fields as `treatment: field`, and merging all treatments into one TSV.

The result is a multi-field per-gene table containing synthesis, degradation and processing rates across time.

## synthesis_changes
Call up/down changes in RNA synthesis rates over time, per treatment.

Using the `INSPEcT` table, this task computes discrete changes (up, down, or none) in synthesis between consecutive timepoints and encodes them as values such as `"up 2"` or `"down 4"` for each treatment.

One new field per treatment is added, named `treatment: INSPEcT synthesis clusters`, containing a list of change events.

## degradation_changes
Call up/down changes in RNA degradation rates over time, per treatment.

This task mirrors `synthesis_changes` but uses `degradation_t` fields instead of synthesis.

One new field per treatment is added, named `treatment: INSPEcT degradation clusters`.

## full_gene_info
Integrate INSPEcT clusters, fold-changes, change_offsets and biotype into a single per-gene table.

This task:

- Loads `synthesis_changes`, `degradation_changes`, `fold_changes_NTNU` (transposed), and `change_offsets`.
- Harmonizes keys to `Associated Gene Name` and then adds an `Ensembl Gene ID` column.
- Attaches gene identifiers and Ensembl transcript annotations to mark protein-coding genes.
- Adds boolean flags:
  - `Protein coding`
  - `FC gene` (has any fold-change in any treatment/timepoint)
  - `INSPEcT gene` (present in INSPEcT with a geneClass)
  - `Dynamic gene` (has at least one non-unclassified FC cluster)

The output is a rich per-gene profile used for downstream classification.

## synthesis_genes
Identify genes whose regulation is dominated by synthesis in INSPEcT.

Based on `full_gene_info`, this task counts, for each gene, how many treatments show significant synthesis, processing, or degradation regulation (p-value < 0.05), and labels a gene as a "Synthesis gene" if it meets a chosen criterion.

Inputs:
- `synthesis_criteria` (`mayority`, `one`, `two`, `three`; default `mayority`): decides how many synthesis-regulated treatments are required and whether synthesis must dominate over processing/degradation.

A new field `Synthesis gene` is added to the TSV from `full_gene_info`.

## vetted_genes
Mark synthesis genes whose expression clusters agree with INSPEcT-derived synthesis/degradation trends.

Using `synthesis_genes`, this task:

- For each gene and treatment, compares INSPEcT degradation clusters with FC change offsets to detect consistent pairs of "up" vs "decrease" or "down" vs "increase".
- Computes strict and relaxed extended degradation timepoints by spreading matches across treatments under direction constraints.
- Determines, per treatment, whether the overall synthesis profile agrees with FC clusters (`Treatment synthesis profile match`).
- Adds a final `Vetted synthesis gene` flag when a gene is both a synthesis gene and has at least one treatment where INSPEcT and FC clusters match.

The output is the main vetted gene table used in TF activity inference and filtering.

## inspect_genes
Return the list of genes present in INSPEcT.

This task selects `INSPEcT gene == true` from `vetted_genes` and returns their `Associated Gene Name` list.

## dynamic_genes
Return the list of genes with dynamic FC clusters.

This task selects `Dynamic gene == true` from `vetted_genes` and returns their `Associated Gene Name` list.

## expressed_coding_genes
Return the list of protein-coding genes with non-zero fold-changes.

Based on `vetted_genes`, this task selects genes where `FC gene == true` and `Protein coding == true`.

This list is used as a filter for TFs and targets in `filtered_regulome`.

## inspect_dynamic_genes
Return the intersection of INSPEcT genes and dynamic genes.

This task returns a list of genes that are both INSPEcT-measured and dynamically clustered (`inspect_genes ∩ dynamic_genes`).

## vetted_genes_in_CollecTRI
Annotate vetted genes with whether they appear as targets in CollecTRI.

This task loads CollecTRI from the `ExTRI` workflow, retrieves all CollecTRI target genes and adds an `In CollecTRI` boolean flag to the `vetted_genes` table.

## vetted_gene_list
Return the list of vetted synthesis genes.

This task selects genes with `Vetted synthesis gene == true` from `vetted_genes` and returns their `Associated Gene Name` list.

## vetted_gene_enrichment
Compute TF enrichment of vetted genes using CollecTRI as TF–target mapping.

This task uses CollecTRI’s TF→target mapping as a gene set collection and performs a hypergeometric enrichment of vetted genes against all genes found in `vetted_genes`.

The output TSV summarizes enrichment statistics per TF.

## vetted_gene_overlap
Summarize overlap between vetted genes and all CollecTRI targets.

This task computes simple counts and percentages of vetted genes overlapping CollecTRI targets and writes them as a single-row TSV keyed by statistic name.

## regulome
Provide a regulome (TF→target mapping) from ExTRI2 or SaezLab.

This task alias resolves to:

- `ExTRI2.regulome` by default (`ExTRI2_regulome = true`), or
- `SaezLab.regulome` otherwise.

Inputs:
- `ExTRI2_regulome` (boolean, default `true`).
- Additional options are passed through to the underlying regulome task (e.g. `only_authoritative_tfs`, `remove_auto_regulation`, `no_MoR`).

The output is a TSV with key field `source` and target and weight columns.

## timepoint_matrix
Build a single-sample expression matrix for a given treatment and timepoint.

Using `vetted_genes`, this task constructs a one-column matrix containing gene-level values for a specific `(treatment, time_point)` and `data_type`. It only includes:

- Protein-coding genes.
- Genes marked as dynamic and INSPEcT genes.

Inputs:
- `treatment` (select from `TREATMENTS`, required).
- `time_point` (select from `TIME_POINTS`, required).
- `data_type` (`fc`, `fc0`, `ext_fc`, `binary`, `range`; default `ext_fc`).

The output is a TSV with one column named `treatment-Ttime` and key field `Associated Gene Name`. It is mainly used as input to `timepoint_decoupler_pre`.

## decoupler_targets
Define the subset of genes to use as decoupler targets for a specific treatment and timepoint.

This task filters vetted genes to those that are:

- Not in the noisy gene list.
- Optionally restricted to a specified target subset.
- Optionally restricted to dynamic genes at this timepoint and direction.
- Optionally filtered by vetting scheme (`synthesis`, `degradation`, `relaxed_degradation`, `none`).

Inputs include:
- `treatment` (select).
- `time_point` (select).
- `vetting` (`none`, `synthesis`, `degradation`, `relaxed_degradation`).
- `dynamic` (boolean): whether to require dynamic changes at this timepoint.
- `dynamic_subset` (`increase`, `decrease`, `both`).
- `target_subset` (array): optional list of genes to restrict to.

The output is an array of target gene names.

## filtered_regulome
Filter the regulome to keep only relevant edges for decoupler.

This task intersects the regulome with `decoupler_targets` and `expressed_coding_genes`, and optionally restricts TFs to dbTFs.

Inputs:
- `only_dbTF` (boolean, default `false`).

The output is a TF→target mapping suitable as a network input for decoupler.

## timepoint_decoupler_pre
Run decoupler on a single timepoint matrix and filtered regulome (internal helper).

This `dep_task` delegates to `SaezLab.decoupler`, passing:

- `matrix`: the `timepoint_matrix` for a given treatment and timepoint.
- `network`: the `filtered_regulome`.

The raw decoupler output is then processed by `timepoint_decoupler`.

## timepoint_decoupler
Apply p-value filtering to decoupler results at a single timepoint.

This task thresholds the output of `timepoint_decoupler_pre`: for each TF, it merges effect/p-value columns (or consensus estimates), zeroing non-significant entries.

Inputs:
- `threshold` (float, default `0.05`).

The output is a TSV keyed by TF, with the column representing the selected `(treatment, timepoint)` sample.

## timepoint_decoupler_suite
Run decoupler across all treatments and timepoints and merge into a single TF activity matrix.

This task depends on a grid of `timepoint_decoupler` jobs that vary treatment, timepoint, data type, vetting and other options following a hard-coded scheme. It:

- Loads each `timepoint_decoupler` output.
- Removes not-expressed and noisy genes.
- Attaches all columns into a single TSV keyed by TF.

This suite is a central input for heatmaps and Venn diagrams.

## timepoint_decoupler_excel
Export per-treatment TF activity tables to Excel.

Using `timepoint_decoupler_suite`, this task creates one Excel workbook per treatment, each containing TF activities across time and parameter settings.

The output is the list of produced `.xlsx` files.

## timepoint_decoupler_venn
Plot a Venn diagram of TF activity overlap between two conditions.

This task loads a precomputed `timepoint_decoupler_suite` TSV from disk and computes TFs that are active (positive or negative) in each of two specified condition columns. It then calls an R plotting function to draw a Venn diagram.

Inputs:
- `condition1`, `condition2` (strings): column names in the suite.
- `direction` (`up`, `down`, `both`): whether to consider up-regulation, down-regulation or both.

The output is a PNG binary with the Venn diagram.

## treatment_tfs_diff
Compute differences in TF activities between consecutive timepoints for a single treatment.

This task chains multiple `timepoint_decoupler` calls per timepoint (with `data_type = fc0`, dynamic disabled), and for each timepoint subtracts the previous one, yielding per-TF activity deltas.

The output is a TSV keyed by TF, with one column per `Ttime` difference for the chosen treatment.

## treatment_tfs_non_dynamic
Collect TF activities across timepoints without dynamic filtering for a single treatment.

This task collects `timepoint_decoupler` outputs (with `data_type = fc`, dynamic disabled) across timepoints and concatenates them into a TF × timepoint activity matrix.

## treatment_tfs_fc0
Collect TF activities on absolute fold-changes (fc0) across time for a single treatment.

Similar to `treatment_tfs_non_dynamic` but based on `fc0` instead of one-step FC changes.

## treatment_tfs_dynamic
Collect TF activities across time for dynamic genes only for a single treatment.

This task runs `timepoint_decoupler` with `dynamic = true` across timepoints and concatenates results into TF × timepoint matrices. It is the main input to higher-level `treatment_tfs` and TF predictions.

Inputs include `data_type` for interaction with upstream dep options.

## treatment_tfs_priority
Score TF activities across multiple vetting and thresholding schemes to build a priority score.

This task:

- Starts from a baseline dynamic decoupler result.
- Applies multiple exclusion schemes (e.g. different vetting/thresholds) to zero out inconsistent TF activities.
- Adds `Score <field>` columns and updates them according to how consistent each TF’s activity direction is across different scoring jobs.

Inputs:
- `score_points` (array of floats, default `[1, 2, 0.5, 0.5, 0]`): weights for each scoring job.

The output is a TF × timepoint activity matrix with additional score columns.

## treatment_tfs
High-level entry point for per-treatment TF time series under different schemes.

This is a `task_alias` over `treatment_tfs_dynamic`, selecting the underlying implementation based on the `scheme` input:

- `priority` → `treatment_tfs_priority`.
- `dynamic` → `treatment_tfs_dynamic`.
- `non-dynamic` → `treatment_tfs_non_dynamic`.
- `diff` → `treatment_tfs_diff`.
- `fc0` → `treatment_tfs_fc0`.

Inputs:
- `scheme` (`priority`, `dynamic`, `non-dynamic`, `diff`, `fc0`; default `dynamic`).
- `treatment` (select from `TREATMENTS`).

The output is a TF × timepoint activity matrix.

Example:

```bash
rbbt workflow AGS treatment_tfs \
  --treatment=INT_PD_PI \
  --scheme=dynamic \
  > INT_PD_PI_tfs.tsv
```

## tf_predictions
Merge dynamic TF activity matrices across all treatments.

This task runs `treatment_tfs` (dynamic scheme, `vetting = none`, `data_type = range`) for all treatments and attaches the resulting matrices into a single TSV keyed by TF.

The output can be used to query, for any TF, its activity across all treatments and timepoints.

## timepoint_heatmaps
Generate heatmaps of TF activities across treatments at each timepoint and parameter setting.

This task uses `timepoint_decoupler_suite` and, for each combination of synthesis/dynamic/data_type/synthesis_criteria settings, slices the suite to build treatment × TF matrices at each timepoint and plots them via R’s `heatmap.2`.

The output is a collection of PNG files organized in directories named after the parameter combinations.

## downstream_targets
Propagate treatment targets downstream through the SIGNOR network and find regulome targets.

This task starts from a predefined map of direct drug targets per treatment (`TREATMENT_TARGETS`), then:

- Loads the SIGNOR protein–protein interaction network and maps it to gene symbols.
- For each source target gene, performs a breadth-first search up to `steps_max` to find reachable genes.
- Restricts to genes that appear as targets in the regulome.
- Records the sign (+1 or -1) of regulation based on SIGNOR effect annotations.

Input:
- `steps_max` (integer, default `3`): maximum path length.

The output is a YAML mapping each drug target gene to a list of (target, sign) pairs.

## downstream_consistency
Export downstream target consistency tables for each treatment.

Using `downstream_targets` and an externally provided `timepoint_suite` TSV (e.g. `timepoint_decoupler_suite`), this task:

- Derives, for each treatment, a target set with associated sign information.
- For each parameter combination, extracts decoupler activities at 1 h and attaches the expected sign per target gene.
- Writes treatment-specific TSV files summarizing sign consistency.

The output is a list of generated files under the job directory.

## treatment_tf_consistency
Assess consistency between TF activities and gene-level change offsets for a given treatment.

This task compares per-timepoint TF activities (from `treatment_tfs`) with the per-gene change offsets (from `change_offsets`) for a selected treatment. For each TF and timepoint, it returns:

- `1` for consistent evidence.
- `-1` for inconsistent evidence.
- `0` or `nil` otherwise.

Inputs:
- `bona_fide` (boolean, default `false`): restrict to dbTFs only.
- `remove_consecutive` (boolean, default `true`): collapse consecutive change points to avoid over-counting.

The output is a TSV where some columns are consistency scores named `Consistent at Xh`.

## consistency_summary
Compute a global summary consistency score across treatments for dbTFs.

This task filters `treatment_tf_consistency` to dbTFs and computes the average number of consistent occurrences across treatments and timepoints.

The output is a single float value.

## consistency_sweep
Sweep consistency across multiple `scheme` and `vetting` combinations for all treatments.

This task runs `treatment_tf_consistency` for several `(scheme, vetting, treatment)` combinations (dynamic/non-dynamic/fc0 and none/relaxed_degradation). It then attaches all resulting TSVs, renaming fields to encode treatment, scheme and vetting.

The output is a large TSV of consistency scores for all tested setups.

## consistency_counts
Summarize counts of consistent vs inconsistent TFs for each experiment.

Using `consistency_sweep`, this task iterates over all `Consistent` columns and, for each, counts matches (`1`), mismatches (`-1`) and totals, computing odds and proportions.

The output is a TSV with one row per `(treatment, time, scheme, vetting)` combination.

## neko_consistency
Compare inferred TF activities with NeKo benchmark predictions for a treatment.

This task loads NeKo benchmark prediction tables specific to `PI`, `PD`, or `FiveZ`, and compares sign predictions with TF activities from `treatment_tfs`. It records, per TF, whether 1 h and 2 h activities match the predicted sign.

Input:
- `dbTFs` (boolean, default `true`): restrict to dbTFs.

The output is a TSV with per-TF match flags and activity values.

## neko_summary
Summarize match vs mismatch counts against NeKo predictions.

Given `neko_consistency`, this task counts matches and mismatches according to a chosen comparison target:

- `strict`: both timepoints counted independently.
- `relaxed`: allow partial matches.
- `T1` / `T2`: focus on a single timepoint.

Input:
- `target` (`strict`, `relaxed`, `T1`, `T2`).

It returns an array `[match_count, mismatch_count]`.

## neko_sweep
Sweep NeKo consistency across treatments, schemes, vettings and targets.

This task runs `neko_summary` across combinations of `target`, `scheme` (dynamic vs non-dynamic), `vetting` (none vs relaxed_degradation) and treatment (`PI`, `PD`, `FiveZ`). It then assembles a TSV of match/mismatch counts and odds.

## neko_bootstrap
Load and aggregate NeKo bootstrap cumulative effect tables for a treatment.

This task reads multiple bootstrap `.csv` files from the `TFe-nekofinder` module, aggregates predicted signs per TF and target, and keeps either only the main targets or all targets depending on configuration.

Input:
- `treatment` (select: `PI`, `PD`, `FiveZ`, `PD-MAP2K2`, `PD-MAP2K1`).

The output is a TSV keyed by TF, with a `Signs` list of predicted sign values over bootstraps.

## neko_bootstrap_overview
Summarize NeKo bootstrap sign distributions per TF.

This task filters TFs with enough bootstrap samples, then adds `Positive` and `Negative` fields representing the fraction of positive vs negative effects across bootstraps.

## neko_bootstrap_final
Select robust TF predictions from NeKo bootstrap results.

Using `neko_bootstrap_overview`, this task selects TFs with strong positive or negative consensus (positive fraction ≥ 0.8 or ≤ 0.2) and encodes them as +1 or −1 in a single-column TSV keyed by TF.

## neko_bootstrap_consistency
Compare inferred TF activities to robust NeKo bootstrap predictions.

This task aligns predicted signs from `neko_bootstrap_final` with TF activities from `treatment_tfs`, computing match indicators at early timepoints.

Inputs:
- `treatment` (select; PD may be split into MAP2K1/2 subcases).
- `dbTFs` (boolean, default `true`).

The output is a TSV of per-TF matches.

## neko_bootstrap_summary
Summarize overall consistency with NeKo bootstrap predictions.

Given `neko_bootstrap_consistency`, this task counts matches and mismatches under different `target` strategies (`relaxed`, `T1`, `T2`, etc.) and returns `[match, mismatch]`.

## neko_bootstrap_sweep
Sweep NeKo bootstrap consistency across treatments, schemes and vettings.

This task runs `neko_bootstrap_summary` for multiple combinations of `target`, `scheme`, `vetting` and treatment, then assembles a TSV of match/miss/total/odds per combination.

## gs_hyper
Perform simple gene-set enrichment using a hypergeometric test.

This task uses `Organism.gene_go_bp` annotations (mapped to `Associated Gene Name`) to compute GO:BP enrichment for an input gene list against a chosen background (default: `expressed_coding_genes`).

Inputs:
- `list` (array): query gene set.
- `background` (array, default `expressed_coding_genes`).
- `database` (select, currently only `go_bp`).

The output is a TSV of GO terms with enrichment statistics.

## gprofiler
Run g:Profiler enrichment for a single gene set.

This task calls the Python g:Profiler client via `RbbtPython`, performing GO:BP enrichment for a gene list, optionally with a given background.

Inputs:
- `list` (array or Step): query genes.
- `background` (array, default `expressed_coding_genes`).
- `database` (currently `go_bp`, mapped to GO:BP in g:Profiler).

The output is a TSV similar to g:Profiler’s tabular output, with p-values, term names, and compacted gene lists.

## gprofiler_multiple
Run g:Profiler enrichment for multiple named queries.

This task is similar to `gprofiler` but takes a YAML mapping of query names to gene lists and runs g:Profiler in multiple-query mode.

Inputs:
- `queries` (YAML): mapping name → [genes].
- `background` and `database` as in `gprofiler`.

The output TSV contains results for all queries, with query identifiers in the fields.

## treatment_overview
Merge multiple TF activity matrices (dynamic, diff, non-dynamic) under different regulome/vetting options.

This task depends on `treatment_tfs_dynamic`, `treatment_tfs_diff`, and `treatment_tfs_non_dynamic` computed under both ExTRI1 (SaezLab) and ExTRI2 regulomes and multiple vettings. It renames columns to encode `ExTRI1/ExTRI2`, vetting scheme, and task name, then attaches all matrices into a single overview TSV.

## treatment_overview_sweep
Create per-treatment overview files and a global merged TSV.

For each treatment and method, this task:

- Links the `treatment_overview` result into a named file under the job directory.
- Attaches all treatment-specific overviews into a global TSV.

The result is both side-effect files and a merged TF activity overview.

## coTF_sign
Count sign occurrences for co-regulator TFs in the ExTRI2 regulome.

This task loads `ExTRI2_clean`, reads a list of coTFs, and counts how many signed interactions (`UP` or `DOWN`) each coTF participates in.

The output is a TSV keyed by coTF name with a `Sign` count.

## coTF_plot
Plot coTF sign statistics (incomplete helper).

Nominally, this task should plot `coTF_sign` as a PNG using R. The current implementation skeleton writes an R call but is likely incomplete; treat it as experimental.

## activities_for_signaling
Extract TF activities for a curated set of signaling-related TFs.

This task runs `treatment_tfs_dynamic` with `ExTRI2_regulome = true`, `vetting = 'none'` for all treatments, merges them into a single matrix, and then subsets to TFs listed in `DNA_binding_and_some_likely_coTFs.13012025`.

The output is a TSV of activities for signaling-relevant TFs.

## range_matrix
Build a matrix of range-based FC values across all treatments and timepoints.

This task calls `timepoint_matrix` with `data_type = :range` for all treatments and timepoints, transposes each result, and attaches them into a large gene × (treatment,timepoint) matrix.

Useful for custom downstream processing or visualization.

## faster_in_combinations
List TFs that respond faster in combination treatments than in single-agent treatments.

This task uses `dbTFs` and `tf_predictions` and returns TFs whose activity is present in early timepoints (e.g. 2–4 h) for combinations (INT_PD_PI, INT_FiveZ_PI) but absent in the corresponding single treatments.

Input:
- `bona_fide` (boolean): restrict to dbTFs.

The output is an array of TF gene names.

## faster_in_combination_enrichment
Perform g:Profiler enrichment for TFs that respond faster in combinations.

This `task_alias` wraps `gprofiler` with:

- `list` = `faster_in_combinations`.
- `background` = `dbTFs`.

It returns a GO:BP enrichment TSV for these candidate TFs.

## not_recovered_in_PD_PI
List TFs active in other conditions but not recovered in INT_PD_PI.

This task scans `tf_predictions`, identifying TFs active in any field except the `INT_PD_PI-T24` column but inactive in `INT_PD_PI-T24` itself.

Input:
- `bona_fide` (boolean): restrict to dbTFs.

The output is a list of TF gene names.

## not_recovered_in_PD_PI_enrichment
Perform g:Profiler enrichment for TFs not recovered in INT_PD_PI.

This `task_alias` wraps `gprofiler` with:

- `list` = `not_recovered_in_PD_PI`.
- `background` = `dbTFs`.

The output is a GO:BP enrichment TSV.

## statistics
Compute summary statistics on TF predictions.

This task uses `tf_predictions` to compute basic counts: number of TFs overall, number of bona fide TFs, and per-field counts of TFs inferred as active, split by dbTF membership.

The output is a YAML hash of statistics.

## list_tfs
List the top deregulated TFs for a treatment and timepoint.

This user-facing helper uses `treatment_tfs` to select TFs that show positive or negative activity at a given timepoint and ranks them by absolute activity.

Inputs:
- `time_point` (integer, default `1`).
- `direction` (`up` or `down`, default `up`).
- `max` (integer, default `20`): maximum number of TFs to return.

The output is an array of TF names.

Example: top 20 up-regulated TFs for PI at 4 h:

```bash
rbbt workflow AGS list_tfs \
  --treatment=PI \
  --time_point=4 \
  --direction=up \
  --max=20
```

## list_tgs
List target genes whose expression starts to change at a given timepoint, for TFs of interest.

This task combines `list_tfs`, `change_offsets_simplified`, and `regulome` to derive genes which both:

- Are targets of the selected TFs.
- Belong to change-offset clusters (`increase` or `decrease`) at a specified timepoint and treatment.

Inputs:
- Indirectly uses `time_point`, `direction`, `max` and `treatment` from recursive inputs.
- `tfs` (array, optional): explicit TFs to use instead of those from `list_tfs`.

The output is a list of target genes.

## tf_targets
Return all targets and modes of regulation for a TF.

This task slices the regulome to keep only edges for a given TF, then annotates each target with its mode of regulation (`activate` or `inhibit`) based on the sign of the weight.

Input:
- `gene` (string): TF gene symbol.

The output is a TSV keyed by target gene, with mode of regulation as the main field.

## target_analysis
Compare target-level fold-change starts with inferred mode of regulation for a TF.

This higher-level helper combines:

- `tf_targets` for modes of regulation.
- `list_tgs` for up and down regulated targets at a given timepoint.

It marks each target as `Regulated` (`up` or `down`) and `Consistent` (whether regulation matches the expected effect given the TF’s mode of regulation).

## gprofiler_queries
Generate gene list files for offline GO enrichment with g:Profiler.

This task iterates over all `TREATMENTS` and `TIME_POINTS`, builds multiple gene sets per combination (change-offset clusters and FC-based thresholds, for `fc` and `fc0` data), and writes them as text files with a FASTA-like header format suitable for g:Profiler’s web interface.

The output is the string `DONE` after all files are written.

## sequence_old
Infer simple TF–TF activation events from TF activity time series and regulome.

This early implementation uses `treatment_tfs` and `regulome` to:

- Identify TFs with strong activation/repression at each timepoint.
- For each TF, scan its targets for coincident or next-timepoint activity with consistent sign given the edge effect.

The output is an array of hashes, each describing a TF–TF event with source, target, timepoints, activities, and effect.

## sequence
Infer a time-resolved sequence of TF activations with self-consistency constraints.

This is the main sequence inference task. It uses:

- `treatment_tfs` (with fields normalized to `T1`, `T2`, etc.).
- `treatment_tf_consistency` (without collapsing consecutive changes) to identify self-consistent TFs at each timepoint.
- `regulome` edges between TFs.

For each treatment, it:

- Builds TF activity series per timepoint.
- Identifies sustained TF activities that are strong and consistent over at least two timepoints.
- For each TF pair (A→B) connected in the regulome, scans offsets 0–2 (same time, +1, +2 timepoints) and records events where source and target activities match the sign of the edge.
- Requires self-consistency of both source and target at the relevant timepoints.

The output is a TSV keyed by an event ID, with fields describing source, target, timepoints, activities, effect, offset, type (`sustained/coincident` vs `delay_1/2`), and self-consistency flags.

## sequence_with_changes
Augment TF–TF sequence events with gene-level change offsets.

This task attaches `change_offsets` to the source and target TFs in the `sequence` output, adding `Source Changes` and `Target Changes` fields containing the per-treatment change intervals.

This makes it easier to relate TF–TF causal events to the onset of transcriptional changes.

## neko_bootstrap_rule_variants_t1
Explore rule parameter variants and their impact on NeKo bootstrap consistency at early timepoints.

This task sweeps selected `change_offsets` rule parameters (`low_min_1h`, `mid_multiplier_1h`) and evaluates `neko_bootstrap_summary` for multiple targets and treatments.

The output is a TSV summarizing match/miss/total/odds per parameter set.

## consistency_summary_rule_variants_t24
Explore rule parameter variants and their impact on consistency at 24 h.

This task sweeps `high_min_24h` and uses `treatment_tf_consistency` to compute, for each `(treatment, high_min_24h)` combination, how many TFs are consistent at 24 h.

The output is a TSV with counts and proportions per configuration.
