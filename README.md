# 🧬 Methylation Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15150516.svg)](https://doi.org/10.5281/zenodo.15150516)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A Python pipeline for **longitudinal analysis of cell-free DNA (cfDNA) methylation** from enzymatic methyl-seq (EM-seq) data. It takes CpG island (CGI) fragment counts from one or more sequencing runs, converts them into normalized methylation signals, summarizes them at the gene level, and tracks how methylation changes across treatment timepoints within each patient and across a cohort, with comparisons to healthy controls.

The pipeline works at three levels:

| Level | Question it answers | Example outputs |
|:--|:--|:--|
| **Global** (sample) | How does overall methylation signal change over the course of treatment? | Condition dotplots, per-patient trajectories, box/violin plots |
| **Locus** (CpG island / chromosome) | Which regions and chromosomes change most between timepoints? | Bubble plots, per-chromosome change, top differentially methylated CGIs |
| **Gene** | Which genes carry the strongest signal, and how does each gene's rank shift over treatment? | Gene rank trajectories, top-gene heatmaps, paired t-tests |

> **Early detection.** The [`early_detection/`](early_detection/) folder applies the same normalization to a different design: plasma drawn 15-9 and 9-3 months before a cancer diagnosis and near diagnosis, compared with healthy controls. It reconstructs the pipeline's first (2024) use case and runs end to end on a synthetic demo dataset. See [`early_detection/README.md`](early_detection/README.md).

---

## 📑 Contents

- [Workflow overview](#-workflow-overview)
- [Early detection (separate folder)](early_detection/README.md)
- [Quick start](#-quick-start)
- [Input requirements](#-input-requirements)
- [Preprocessing (Steps 1–6)](#-preprocessing-steps-16)
- [Global-level analysis](#-global-level-analysis)
- [Locus-level analysis](#-locus-level-analysis)
- [Gene-level analysis and ranking](#-gene-level-analysis-and-ranking)
- [Script versions](#-script-versions)
- [Project structure](#-project-structure)
- [Known limitations](#-known-limitations)
- [CGI naming format](#-cgi-naming-format)
- [Citation](#-citation)
- [License](#-license)

---

## 🗺 Workflow overview

```
data/ (EM-seq CGI count files + patient ID list)
   │
   ├─ Step 1  Filter samples by patient ID
   ├─ Step 2  Merge runs from multiple EM-seq batches        (optional)
   ├─ Step 3  Glob20 / GlobMin80 ratio → scaled methylation signal
   ├─ Step 4  Parse CGI names into gene annotations           (needed for gene-level)
   ├─ Step 5  Build CGI → gene map                            (needed for gene-level)
   └─ Step 6  Sum CGI signals into a gene × sample matrix     (needed for gene ranking)
   │
   ├─▶ Global analysis   (scripts/global/)
   ├─▶ Locus analysis    (scripts/locus/)
   └─▶ Gene analysis     (scripts/locus/ — gene rank and top-gene scripts)
```

---

## ▶️ Quick start

```bash
git clone https://github.com/codebyesther/methylation-pipeline.git
cd methylation-pipeline
pip install -r requirements.txt

mkdir -p data output plots
# place your patient ID list and EM-seq Excel files in data/

python scripts/step_1_filter_patients_local.py
python scripts/step_2_merge_filtered_files.py        # if you have multiple runs
python scripts/step_3_convert_to_aberrant_signals.py
python scripts/step_4_generate_gene_annotation.py
python scripts/step_5_build_gene_cgi_map.py
python scripts/step_6_gene-methylation-matrix.py
```

Then run any of the analysis scripts below. All scripts are run from the repository root and read from `data/` and `output/`, writing results to `output/` and `plots/`.

**Requirements:** Python 3 with `pandas`, `numpy`, `openpyxl`, `matplotlib`, `seaborn`, `scipy`, `tqdm` (see `requirements.txt`).

---

## 📁 Input requirements

Place all inputs in `data/`.

**1. Patient ID list** (`.xlsx`)
- Patient IDs in the first column.
- The filename must contain `patient` (most scripts) or `patient id` (the `top10dm-plots*` scripts). A name like `Patient ID list.xlsx` satisfies both.
- The gene-rank scripts currently expect a specific filename; see [Known limitations](#-known-limitations).

| Patient ID |
|:--|
| PT001 |
| PT002 |
| PT003 |

**2. EM-seq methylation files** (`.xlsx`)
- One row per CpG island, with CGI names in the first column (e.g. `CGI_chr1_778604_779167_LOC100288069_0`).
- One column per sample. Each column name must contain the patient ID and a timepoint label.
- A summary row labeled `Total CpG island fragments counts for this particular spreadsheet` holds the total CGI fragment count per sample. The global-level scripts use this row.
- Each run is expected as a pair of files: a **Glob20** file and a **GlobMin80** file, used together in Step 3.

**Sample naming conventions.** Timepoints are read from the sample (column) names:

| Label in sample name | Interpreted as |
|:--|:--|
| `Baseline` | Baseline (pre-treatment) |
| `C<n>` (e.g. `C2D1`, `C1D15`) | Treatment cycle *n* (On-Treatment) |
| `Off-tx` | Off treatment (Post-Treatment) |
| `INNOV` | Healthy control |

Example sample names: `PT001_Baseline`, `PT001_C2D1`, `PT001_C7D1-Off-tx`, `INNOV_control_01`.

---

## 🔧 Preprocessing (Steps 1–6)

### Step 1 — Filter by patient ID
`scripts/step_1_filter_patients_local.py`

- Reads the patient ID list and every methylation Excel file in `data/`.
- Keeps only the sample columns that belong to the listed patients.
- Writes one filtered file per input to `output/` (e.g. `EMSeq-NN-Run-1_..._Glob20_cpgi_counts.xlsx_samples-of-interest`).

### Step 2 — Merge multiple EM-seq runs *(optional)*
`scripts/step_2_merge_filtered_files.py`

- Needed only if your samples span more than one sequencing run.
- Finds the filtered files in `output/` with `Glob20` or `GlobMin80` in their names and merges each set by row label.
- Outputs `output/merged_output_glob20.xlsx` and `output/merged_output_globmin80.xlsx`.

### Step 3 — Convert to normalized methylation signal
`scripts/step_3_convert_to_aberrant_signals.py`

- Reads `output_glob20*.xlsx` and `output_globmin80*.xlsx` from `output/`.
- For each CGI and sample, divides the Glob20 fragment count by the matching GlobMin80 count and multiplies by a scaling factor, producing a normalized methylation signal.
- Keeps the total-count summary row so global-level scripts can use it.
- Outputs `output/scaled_fragment_ratios_matrix.xlsx`.

Two alternative versions are included (see [Script versions](#-script-versions)). The `v3` version aligns rows and columns by name rather than position, marks division-by-zero cells as `INF` (highlighted red in Excel), and logs them to `output/inf_locations_log.csv` and `output/inf_summary.csv`.

### Step 4 — Annotate genes
`scripts/step_4_generate_gene_annotation.py`

- Reads an Excel file with `matrix` in its name from `output/`.
- Parses each `CGI_...` name into chromosome, start, end, gene name(s) and index.
- Outputs `output/structured_gene_annotation.csv`.

### Step 5 — Build the CGI → gene map
`scripts/step_5_build_gene_cgi_map.py`

- Reads the `gene_annotation` file from Step 4.
- Reshapes it into one row per CGI–gene pair, with CGI IDs in `chr:start-end` format.
- Outputs `output/gene_cgi_map.csv`.

| cgi_id | gene_name |
|:--|:--|
| chr1:778604-779167 | LOC100288069 |
| chr1:827036-828065 | LINC00115 |
| chr1:940225-942983 | SAMD11 |

### Step 6 — Build the gene methylation matrix
`scripts/step_6_gene-methylation-matrix.py`

- Reads `output/merged_output_glob20.xlsx` and `output/gene_cgi_map.csv`.
- Sums the signal of all CpG islands assigned to each gene, giving one value per gene per sample.
- Outputs `output/gene_methylation_matrix.csv` (genes × samples). This is the input for the gene-rank scripts.

---

## 🌐 Global-level analysis

These scripts use each sample's total CGI signal (the summary row from Step 3) to look at overall trends.

### Dotplots by condition
`scripts/global/dotplots-by-condition.py`

- Input: `output/*fragment_ratios_matrix*.xlsx` (from Step 3).
- Groups samples into Healthy, Baseline, On-Treatment and Post-Treatment, and computes the n, mean, median and SD for each group.
- Outputs to `plots/dotplots/`: `*_summary_stats.csv`, `*_median_dotplot.png`, `*_mean_sd_dotplot.png`.

<img src="sampleimages/scaled_fragment_ratios_matrix_median_dotplot.png" width="450"/> <img src="sampleimages/scaled_fragment_ratios_matrix_mean_sd_dotplot.png" width="450"/>

### Per-patient trajectories
`scripts/global/lineplots-perpatient_v3.py` (grouped timepoints)
`scripts/global/lineplots-perpatient_v4.py` (individual treatment cycles)

- Input: `output/scaled_fragment_ratios_matrix.xlsx`.
- `v3` groups samples into Baseline / On-Treatment / Post-Treatment. `v4` keeps each treatment cycle as its own point on the x-axis.
- Plots: longitudinal trajectories for all patients, the average trajectory, box plots, violin + swarm plots, and one plot per patient.
- Tables: `sample_metadata.csv`, `summary_statistics.csv`, `boxplot_summary_by_timepoint.csv`, and per-patient replicate tables in `output/per_patient_tables/`.
- `v3` writes plots to `plots/lineplots/`. `v4` writes them to `plots/global-lineplots/`.

<img src="sampleimages/methylation_longitudinal_plot.png" width="600"/>
<img src="sampleimages/average_trajectory.png" width="450"/> <img src="sampleimages/violinplot_by_timepoint.png" width="450"/>

### Plot by cycle *(utility)*
`scripts/plot-by-cycle.py`

A short helper that draws line plots from a manually assembled table (patients × cycle groups). It reads a fixed file path, so edit `data_file` before using it.

---

## 🔍 Locus-level analysis

### Average methylation change per chromosome
`scripts/locus/avg-methylation-change-per-chromosome.py`

- Inputs: a `patient` file in `data/` and a `scaled` `.xlsx` file in `output/`.
- Averages CGI signal per chromosome for each patient and timepoint (healthy controls excluded), then computes the change between timepoints.
- Outputs to `plots/avg-methylation-change-per-chromosome/`: `chr_avg_overlay_<name>_aligned.png` and `chr_avg_summary_<name>.xlsx`. Chromosomes are shown in karyotype order.

<img src="sampleimages/chr_avg_overlay_scaled_fragment_ratios_matrix_aligned.png" width="900"/>

### Bubble plots of CGI methylation along each chromosome
`scripts/locus/bubbleplot_generator_v9_gridsoff.py` (latest)

- Inputs: a `*ratios_matrix*` file in `output/` and a `patient` file in `data/`.
- For each patient and chromosome, plots every CpG island at its genomic midpoint, sized and colored by signal, with one row per timepoint. Also makes cohort-averaged plots for each chromosome.
- Saves PNG and SVG files, bundled into `plots/bubbleplots.zip`.
- `merged-glob20-bubbleplot.py` produces the same plots directly from `output/merged_output_glob20.xlsx`.

<img src="sampleimages/bubbleplot_chr1.png" width="900"/>

### Top 10 differentially methylated CGIs and multi-CGI genes
`scripts/locus/top10dm-plots_using-map_v4.py` (latest)

- Inputs: a `patient id` file in `data/`, a `matrix` file in `output/`, and `output/gene_cgi_map.csv`.
- Averages each CGI across patients at each timepoint, then finds the 10 CGIs with the largest change for Baseline → On-Treatment, Baseline → Post-Treatment, and On-Treatment → Post-Treatment.
- Also identifies genes with more than one affected CGI.
- Outputs bar charts plus the intermediate tables (`v3`/`v4`), bundled into a zip file in `plots/`.

<img src="sampleimages/top10_diff_CGIsubregions_baseline_on.png" width="300"/> <img src="sampleimages/top10_diff_CGIsubregions_on_post.png" width="300"/> <img src="sampleimages/top10_diff_CGIsubregions_baseline_post.png" width="300"/>

---

## 🧪 Gene-level analysis and ranking

### Gene rank trajectories
`scripts/locus/gene_rank_slopeplot_top10_colored_v2.py` (recommended)

This is the gene-ranking part of the pipeline. Instead of comparing raw signal values, which can vary in scale between samples and runs, it compares each gene's **rank** within its sample.

1. **Rank genes within each sample.** Genes are ranked by methylation signal (rank 1 = most methylated). Genes ranked 501 or lower are grouped together at rank 501 to reduce noise from genes with little signal.
2. **Assign patients and timepoints.** Each sample is mapped to its patient and timepoint (Baseline, each treatment cycle, Off-Tx, or Healthy).
3. **Average replicates.** Ranks from replicate samples (same patient, same timepoint) are averaged.
4. **Track each gene over time.** For each patient, a slope plot follows every gene's rank from baseline through each treatment cycle to the end of treatment. Genes that rank in the top 10 at any of that patient's timepoints are highlighted.

- Input: `output/gene_methylation_matrix.csv` (from Step 6) and the patient ID list.
- Outputs to `plots/rank-slopeplot/`:
  - `rank_slopeplot_patient_<ID>.png`: rank trajectories for each patient (log scale)
  - `gene_methylation_ranks.csv`: the full gene × sample rank matrix
  - `melted_gene_methylation_ranks_avg.csv`: ranks in long format, averaged over replicates
  - `top_genes_per_patient.csv`: the highlighted top genes for each patient

Variants (see [Script versions](#-script-versions)) differ in the rank cap (none, 21, or 501) and in whether replicates are averaged.

### Top genes: deltas, paired t-tests and heatmaps
`scripts/locus/top10genes-heatmap-barplot.py`

- Inputs: a `matrix` file in `output/`, a `patient` file in `data/`, and `output/gene_cgi_map.csv`.
- Builds a gene-level matrix (genes with multiple CGIs), computes each gene's change between timepoints for every patient, and runs **paired t-tests** for Baseline vs On-Treatment, Baseline vs Post-Treatment, and On-Treatment vs Post-Treatment.
- Selects the top 10 genes by size of the Baseline → Post-Treatment change and plots heatmaps, line plots and bar plots.
- Outputs to `plots/heatmaps-lineplots/`: `gene_deltas_all_comparisons.csv`, `*_ttest.csv`, `patient_deltas_*.csv`, `top10_gene_deltas_*.csv`, `gene_methylation_matrix.csv`, and figures.

<img src="sampleimages/avg_methylation_heatmap.png" width="500"/> <img src="sampleimages/avg_methylation_lineplot.png" width="500"/>

### Top genes across individual treatment cycles
`scripts/locus/top10genes-barplot-heatmap-lineplot.py`

- Inputs: `output/merged_output_glob20.xlsx`, a `patient` file, and `output/gene_cgi_map.csv`.
- Selects the 10 genes with the highest baseline signal and follows them across each treatment cycle (not grouped timepoints).
- Outputs to `plots/heatmaps-lineplots/`: `heatmap_top10_genes_detailed_timepoints.png`, plus cohort and per-patient line plots in `gene-lineplots/`.

---

## 🗂 Script versions

Several scripts have more than one version. The table shows which one to use.

| Purpose | Recommended | Other versions |
|:--|:--|:--|
| Step 3 signal conversion | `step_3_convert_to_aberrant_signals_v3-matchingRowColumnNames.py` when sample rows or columns may not line up between the Glob20 and GlobMin80 files | `step_3_..._signals.py` (×1,000 scaling, position-based); `_v2` (×100,000 scaling, row alignment) |
| Global trajectories | `lineplots-perpatient_v3.py` (grouped) or `_v4.py` (per cycle) | — |
| Bubble plots | `bubbleplot_generator_v9_gridsoff.py` | `_v8.py`, `_v8_gridsoff.py`, `merged-glob20-bubbleplot.py` |
| Top differentially methylated CGIs | `top10dm-plots_using-map_v4.py` | `top10dm-plots.py`, `_using-map_v2.py`, `_v3.py` |
| Gene heatmaps / t-tests | `top10genes-heatmap-barplot.py` | `heatmap-lineplot-barplot_v2.py`, `_v3.py` |
| Gene rank trajectories | `gene_rank_slopeplot_top10_colored_v2.py` | `generank-slopeplot.py` (no cap), `-21cap.py`, `-501cap.py`, `gene_rank_slopeplot_top10_colored.py` (no replicate averaging) |

> **Note on Step 3 scaling:** the original Step 3 script multiplies ratios by 1,000; `v2` and `v3` multiply by 100,000. Use the same version for every sample in an analysis so values stay comparable.

---

## 📂 Project structure

```
methylation-pipeline/
├── data/                         # Inputs (not tracked): patient list, EM-seq Excel files
├── output/                       # Intermediate matrices and tables
├── plots/                        # Figures and summary tables
├── sampleimages/                 # Example figures used in this README
├── early_detection/              # Pre-diagnosis vs healthy-control version (own README)
├── scripts/
│   ├── step_1_filter_patients_local.py
│   ├── step_2_merge_filtered_files.py
│   ├── step_3_convert_to_aberrant_signals.py
│   ├── step_3_convert_to_aberrant_signals_v2.py
│   ├── step_3_convert_to_aberrant_signals_v3-matchingRowColumnNames.py
│   ├── step_4_generate_gene_annotation.py
│   ├── step_5_build_gene_cgi_map.py
│   ├── step_6_gene-methylation-matrix.py
│   ├── plot-by-cycle.py
│   ├── global/
│   │   ├── dotplots-by-condition.py
│   │   ├── lineplots-perpatient_v3.py
│   │   └── lineplots-perpatient_v4.py
│   └── locus/
│       ├── avg-methylation-change-per-chromosome.py
│       ├── bubbleplot_generator_v8.py
│       ├── bubbleplot_generator_v8_gridsoff.py
│       ├── bubbleplot_generator_v9_gridsoff.py
│       ├── merged-glob20-bubbleplot.py
│       ├── top10dm-plots.py
│       ├── top10dm-plots_using-map_v2.py
│       ├── top10dm-plots_using-map_v3.py
│       ├── top10dm-plots_using-map_v4.py
│       ├── heatmap-lineplot-barplot_v2.py
│       ├── heatmap-lineplot-barplot_v3.py
│       ├── top10genes-heatmap-barplot.py
│       ├── top10genes-barplot-heatmap-lineplot.py
│       ├── generank-slopeplot.py
│       ├── generank-slopeplot-21cap.py
│       ├── generank-slopeplot-501cap.py
│       ├── gene_rank_slopeplot_top10_colored.py
│       └── gene_rank_slopeplot_top10_colored_v2.py
├── CITATION.cff
├── LICENSE
├── README.md
└── requirements.txt
```

---

## ⚠️ Known limitations

- **Fixed patient-list filename in the gene-rank scripts.** The `generank-*` and `gene_rank_slopeplot_*` scripts read `data/Patient ID list fot EMseq16-18-20.xlsx`. Rename your file to match, or edit `patient_list_path` at the top of the script.
- **`generank-slopeplot.py`** reads its gene matrix from `plots/heatmaps-lineplots/gene_methylation_matrix.csv` (written by the heatmap scripts). The other rank scripts read `output/gene_methylation_matrix.csv` (Step 6).
- **File auto-detection uses the first match.** If more than one file in a folder matches a keyword (e.g. `matrix`), scripts use the first one found. Keep only one candidate file per keyword in `data/` and `output/`.
- **Timepoint parsing depends on sample names** following the conventions in [Input requirements](#-input-requirements).
- **Within-patient ranking.** The gene-rank scripts track rank changes within each patient. They don't run a statistical test that combines ranks across patients.

---

## 🧠 CGI naming format

`CGI_chr1_778604_779167_LOC100288069_0`

| Part | Meaning |
|:--|:--|
| `CGI_` | CpG island prefix |
| `chr1` | Chromosome |
| `778604_779167` | Start and end coordinates |
| `LOC100288069` | Associated gene |
| `_0` | CGI index / probe ID |

---

## 📖 Citation

If you use this code, please cite:

> Choi, E. (2025). *Methylation Pipeline* (Version 1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.15150516

Citation metadata is also available in [`CITATION.cff`](CITATION.cff) (use **"Cite this repository"** on GitHub).

---

## 📜 License

Released under the [MIT License](LICENSE).
