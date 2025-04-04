# 🧬 Methylation Pipeline

This repository contains preprocessing and visualization tools for analyzing DNA methylation changes across treatment timepoints using filtered EMseq data. It is useful for following longitudinal studies.

## 🔧 Preprocessing Steps

These scripts help prepare your input methylation data.

### Step 1: Filter Methylation Files by Patient ID
- Script: `scripts/step_1_filter_patients_local.py`
- Applies to all methylation Excel files and a patient ID list.
- Filters large Excel files to retain only samples from selected patients.
- Outputs one filtered Excel file per original input in the output directory (e.g. `EMSeq-NN-Run-1_CG80_CH20_Tot10_Glob20_cpgi_counts.xlsx_samples-of-interest`)

#### 📁 Input File Requirements
- One Excel file with patient IDs (e.g. `patient_list.xlsx`)
  - IDs must be in the first column.
- One or more methylation Excel files
  - Columns should contain patient IDs in their names.

These input files are used at the start of Step 1. Place them in a `data/` folder. Filtered outputs will be saved to an `output/` folder and will be needed for downstream analysis.

### Step 2 (Optional): Merge Filtered Files from Multiple EMseq Runs
- Script: `scripts/step_2_merge_filtered_files.py`
- Required only if you have `.xlsx` outputs from multiple different EMseq batches.
- Produces one merged file for unified analysis in the output directory (e.g. `merged_output_glob20.xlsx`, `merged_output_globmin80.xlsx`)

### Step 3: Convert to Aberrant Signals
- Script: `scripts/step_4_generate_gene_annotation.py`
- Takes the cpgi methylation fragment counts in Glob20 Excel files and divides them by corresponding values in GlobMin80 Excel files to create methylation fragment ratios. Then, scales those numbers up to >1 by multiplying 1000 each.
- Produces `scaled_fragment_ratios_matrix.xlsx` in the output directory.

### Step 4: Annotate Genes
- Script: `scripts/step_3_convert_to_aberrant_signals.py`
- Reads an Excel file containing "matrix" in its name and processes its content to generate a CSV file with structured gene annotation data:
  - Extracts CGI names from the first column, filtering those starting with "CGI_".
  - Processes each CGI name by splitting it into parts and extracting chromosome, start and end genomic coordinates, gene names, and probe IDs.
- Produces `structured_gene_annotation.csv` in the output directory.

## 🔍 Locus-Level (CpG Island / Gene) Analysis

This workflow zooms into locus-specific (CpG island or gene-level) methylation dynamics.

### Step 3 (Update Coming): Visualize Methylation by Region
- Script: `scripts/locus/step_3_cpg_subregion_gene_plots_darkblue_local.py`
- Generate plots per CpG island, gene, or genomic region.
- Includes bar plots, heatmaps, line plots, and longitudinal summaries.

### Step 3 (Update Coming): Visualize Chromosome-Level Methylation Change
- Visualizes average methylation change across chromosomes using bar plots (Baseline → On-Tx, On-Tx → Post-Tx) and a line plot (Baseline → Post-Tx).
- Plots are ordered by chromosome number (karyotype order) and designed for intuitive interpretation.

## 🌐 Global-Level (Patient / Timepoint) Analysis

This workflow focuses on overall methylation trends per patient or treatment condition.

### Generate Dotplots by Condition/Timepoint
- Script: `scripts/global/dotplots-by-condition.py`
- Processes Excel file(s) located in the output directory that contains "fragment_ratios_matrix" in its name.
  - Classifies samples into conditions (Healthy, Baseline, On-Treatment, Post-Treatment).
  - Calculates summary statistics (Mean, Median, Standard Deviation), and generates two types of scatter plots: Median Scatter Plot and Mean ± SD Scatter Plot.
- A CSV file containing summary statistics (Mean, Median, Standard Deviation) for each condition is saved as `*_summary_stats.csv` in the plots/dotplots directory.
- Resulting plots are saved as `*_median_dotplot.png` and `*_mean_sd_dotplot.png` in the plots/dotplots directory.

### Generate Trajectory Lineplots, Boxplots and Violin + Swarm Overlay Plots by Condition/Timepoint
- Script: `scripts/global/lineplots-perpatient_v3.py`
- Input file(s): Processes .xlsx or .xls file(s) located in the output directory that contains "scaled_fragment_ratios_matrix" in its name.
- Generated file(s):
  - Main Plot Directory: plots/lineplots
    - methylation_longitudinal_plot.png
    - average_trajectory.png
    - boxplot_by_timepoint.png
    - violinplot_by_timepoint.png
  - Per Patient Plot Directory: plots/lineplots/per_patient
    - Individual patient plots in the format {Patient_ID}.png
  - Replicate Table Directory: output/per_patient_tables
    - Replicate tables for each patient in the format {Patient_ID}.csv
  - Metadata and Summary Files Directory: output
    - sample_metadata.csv
    - summary_statistics.csv
    - per_patient_summary.csv
    - boxplot_summary_by_timepoint.csv

### Generate Heatmaps and Lineplots of Top 10 Genes by Condition/Timepoint
- Script: `scripts/locus/deltagene-heatmaps-lineplots-bothfonts-labels.py`
- Input file(s):
  - CpG Methylation Matrix
    - Directory: output
    - Identifier: Contains the keyword "matrix"
    - Format: .xlsx or .txt
  - Patient List
    - Directory: data
    - Identifier: Contains the keyword "patient"
    - Format: .xlsx
  - Gene Annotation Map
    - Directory: output
    - Identifier: Contains the keyword "cgi_map"
    - Format: .xlsx or .csv
- Data Processing Steps:
  - Loads CpG methylation matrix, patient list, and gene annotation map.
  - Matches CpGs to genes and filters genes with multiple CpGs.
  - Calculates delta methylation values between specified timepoints and performs paired t-tests.
  - Generates delta values and statistical comparison results for:
    - Baseline vs Post-Treatment
    - Baseline vs On-Treatment
    - On-Treatment vs Post-Treatment
  - Saves delta values and t-test results to CSV files in the plots/heatmaps-lineplots directory.
  - Identifies top 10 genes with highest delta values (Baseline to Post-Treatment).
  - Computes average gene methylation values and constructs a gene matrix.
- Data Visualization Steps:
  - Generates heatmaps and line plots for average gene methylation across timepoints.
  - Saves heatmaps and line plots to the plots/heatmaps-lineplots directory.
- Generated file(s): plots/heatmaps-lineplots directory
  - Delta Values and T-Test Results
    - gene_deltas_all_comparisons.csv
    - baseline_vs_post_ttest.csv
    - baseline_vs_on_ttest.csv
    - on_vs_post_ttest.csv
  - Visualizations
    - avg_methylation_heatmap.png
    - avg_methylation_lineplot.png
  - Per-Patient Visualizations and Matrices
    - heatmap_{patient}.png
    - lineplot_{patient}.png
    - methylation_matrix_{patient}.csv

---

## ▶️ How to Use

### Option 1: Run Locally
- Clone the repository:
  ```bash
  git clone https://github.com/codebyesther/methylation-pipeline.git
  cd methylation-pipeline
  ```
- Install dependencies:
  ```bash
  pip install pandas numpy openpyxl matplotlib seaborn scipy tqdm
  ```
- Place your input `.xlsx` files in a `data/` folder.
- Run the script (e.g. `step_1_filter_patients.py`):
  ```bash
  python scripts/step_1_filter_patients.py
  ```

### Option 2: Google Colab (if using the Colab version)
- Open the `.ipynb` file in Colab
- Upload your files when prompted
- Download the filtered results

---

## 📦 Dependencies

Install all requirements at once:

```bash
pip install -r requirements.txt
```

- `pandas`
- `numpy`
- `openpyxl`
- `matplotlib`
- `seaborn`
- `scipy`
- `tqdm`

## 📂 Project Structure

```
methylation-pipeline/
├── data/                      # Input Excel files
├── output/                    # Filtered and merged outputs
├── plots/                     # Generated plots and Excel summaries
├── scripts/
│   ├── step_1_filter_patients_local.py
│   ├── step_2_merge_filtered_files.py
│   └── global/
│       └── step_3_chr_avg_overlay_global_local.py
│   └── locus/
│       └── step_3_cpg_subregion_gene_plots_darkblue_local.py
├── README.md
├── .gitignore
├── requirements.txt
├── CITATION.cff
```

## 🧠 Additional Context: Data Annotation Format
A breakdown of CGI names (i.e. `CGI_chr1_778604_779167_LOC100288069_0`) can be done to interpret CGI identifiers:
- `CGI_` → prefix 
- `chr1` → chromosome
- `778604_779167` → genomic coordinates
- `LOC100288069` → gene name
- `_0` → sometimes denotes CpG island index or probe ID

Hopefully, this helps you understand the CpG region labels.

---

## 📜 License

This project is licensed under the MIT License.

## 🧾 Citation

Please cite this work if you use any scripts or plots from this repository.

**Suggested citation:**

```
Choi, E. (2025). Methylation Pipeline for Visualizing Longitudinal Methylation Changes. GitHub repository. https://github.com/codebyesther/methylation-pipeline
```
You can also find formal citation formats in the `CITATION.cff` file or by clicking **"Cite this repository"** on GitHub.

If you use or adapt specific components, please mention the relevant script in your methods or supplementary material.
