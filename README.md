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
- Auto-detect file(s):
- Required only if you have `.xlsx` outputs from multiple different EMseq batches.
- Produces one merged file for unified analysis in the output directory (e.g. `merged_output_glob20.xlsx`, `merged_output_globmin80.xlsx`)

### Step 3: Convert to Aberrant Signals
- Script: `scripts/step_3_convert_to_aberrant_signals.py` 
- Auto-detect file(s):
- Takes the cpgi methylation fragment counts in Glob20 Excel files and divides them by corresponding values in GlobMin80 Excel files to create methylation fragment ratios. Then, scales those numbers up to >1 by multiplying 1000 each.
- Produces `scaled_fragment_ratios_matrix.xlsx` in the output directory.

### Step 4 (Optional): Annotate Genes
- Script: `scripts/step_4_generate_gene_annotation.py`
- Required only if you are going to run `scripts/locus/deltagene-heatmaps-lineplots-bothfonts-labels.py` later.
- Auto-detect file(s): an Excel file containing "matrix" in its name
- Reads an Excel file containing "matrix" in its name and processes its content to generate a CSV file with structured gene annotation data:
  - Extracts CGI names from the first column, filtering those starting with "CGI_".
  - Processes each CGI name by splitting it into parts and extracting chromosome, start and end genomic coordinates, gene names, and probe IDs.
- Produces `structured_gene_annotation.csv` in the output directory.

### Step 5 (Optional): Build a CGI Map
- Script: `scripts/step_5_build_gene_cgi_map.py`
- Required only if you are going to run `scripts/locus/deltagene-heatmaps-lineplots-bothfonts-labels.py` later.
- Auto-detect file(s): .xlsx and .csv files in the output directory that contain the term "gene_annotation" in their filename
  - You had to have run Step 4 of the Preprocessing Steps (`scripts/step_4_generate_gene_annotation.py`) for this file to have been generated in your output directory.
- Processes these files to create a mapping of genes to CGI identifiers:
  - The script selects the first file from the list of matching files and prints its name.
  - The script identifies all columns in the DataFrame that start with "Gene".
  - It reshapes the DataFrame into a long format where each row represents a gene and its associated CGI.
  - The script creates a cgi_id column in the format chr:start-end
  - The final gene_cgi_map DataFrame is created by selecting the cgi_id and gene_name columns and removing duplicates.
  - The script saves the gene_cgi_map DataFrame as a CSV file in the output directory.
- Generated file(s): `gene_cgi_map.csv` in the output directory

## 🌐 Global-Level (Patient / Timepoint) Analysis

This workflow focuses on overall methylation trends per patient or treatment condition.

### Generate Dotplots by Condition/Timepoint
- Script: `scripts/global/dotplots-by-condition.py`
- Auto-detect file(s): Excel file(s) located in the output directory that contains "fragment_ratios_matrix" in its name
  - You had to have run Step 3 of the Preprocessing Steps (scripts/step_3_convert_to_aberrant_signals.py) for this file to have been generated in your output directory.
- Data Processing Steps:
  - Classifies samples into conditions (Healthy, Baseline, On-Treatment, Post-Treatment).
  - Calculates summary statistics (Mean, Median, Standard Deviation), and generates two types of scatter plots: Median Scatter Plot and Mean ± SD Scatter Plot.
- Generated file(s):
  - A CSV file containing summary statistics (Mean, Median, Standard Deviation) for each condition is saved as `*_summary_stats.csv` in the plots/dotplots directory.
  - Resulting plots are saved as `*_median_dotplot.png` and `*_mean_sd_dotplot.png` in the plots/dotplots directory.

### Generate Trajectory Lineplots, Boxplots and Violin + Swarm Overlay Plots by Condition/Timepoint
- Script: `scripts/global/lineplots-perpatient_v3.py`
- Auto-detect file(s): .xlsx or .xls file(s) located in the output directory that contains "scaled_fragment_ratios_matrix" in its name
  - You had to have run Step 3 of the Preprocessing Steps (scripts/step_3_convert_to_aberrant_signals.py) for this file to have been generated in your output directory.
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


## 🔍 Locus-Level (CpG Island / Gene) Analysis
This workflow zooms into locus-specific (CpG island or gene-level) methylation dynamics.

### Generate Average Methylation Change per Chromosome
Looking across chromosomes, we can identify global trends—such as whether certain chromosomes are more epigenetically active or suppressed in response to treatment.​
- Script: `scripts/locus/intuitive_storytelling_avg_methylation_local.py`
- Auto-detect file(s):
  - Patient file: a "patient" file in the data directory
  - Methylation file: a "scaled" methylation file in the output directory
    - You had to have run Step 3 of the Preprocessing Steps (`scripts/step_3_convert_to_aberrant_signals.py`) for this file to have been generated in your output directory.
- Data Processing Steps:
  - Loads the patient file and extracts patient IDs.
  - Loads the methylation data file and prepares it for analysis by extracting CpG island data.
  - Identifies and normalizes timepoints in the sample data (e.g., Baseline, On-Treatment, Post-Treatment, Healthy).
  - Filters and structures the data for valid samples, excluding "Healthy" samples.
  - Calculates the average methylation levels for each chromosome.
  - Computes the changes (deltas) in methylation levels between different timepoints for each patient and chromosome.
  - Summarizes the mean changes in methylation levels for each chromosome and comparison. Saves the summary data as an Excel file in the plots directory.
  - Creates bar and line plots to visualize the average methylation changes per chromosome for different comparisons. Plots are ordered by chromosome number (karyotype order) and designed for intuitive interpretation. Saves the plots as PNG files in the plots directory.
- Generated file(s):
  - Plot: A PNG file named chr_avg_overlay_<base_fname>_aligned.png saved in the plots directory.
  - Excel Summary: An Excel file named chr_avg_summary_<base_fname>.xlsx saved in the plots directory.

### Generate Bubble Plots for Visualization of Longitudinal DNA Hypermethylation Profiles per Chromosome
- Script: `scripts/locus/bubbleplot_generator_v8_gridsoff.py`
- Auto-detect file(s):
  - Patient file: a "patient" file in the data directory (both .xlsx and .csv formats)
  - Methylation file: a "ratios_matrix" methylation file in the output directory (both .xlsx and .csv formats)
    - You had to have run Step 3 of the Preprocessing Steps (`scripts/step_3_convert_to_aberrant_signals.py`) for this file to have been generated in your output directory.
- Data Processing Steps:
  - Normalizes sample names to standard timepoints such as "Baseline", "On-Treatment", and "Post-Treatment".
  - Processes each methylation data file to extract CpG Island coordinates and associated values, and merges this data with sample metadata.
- Data Visualization Steps:
- Per Patient Per Chromosome: Generates bubble plots for each patient and chromosome combination, showing the DNA hypermethylation profiles across different timepoints.
- Per Chromosome (Averaged Across Patients): Generates bubble plots for each chromosome, averaged across all patients, to visualize overall methylation patterns.
- Saves the generated plots as PNG and SVG files in the plots directory.
- Creates a ZIP file (bubbleplots.zip) containing all the plot files.
- Deletes the individual plot files after zipping to save space.
- Generated file(s):
  - Bubble plots saved as PNG and SVG files in the plots directory.
  - A ZIP file (bubbleplots.zip) containing all the plot files.

### Generate Heatmaps and Lineplots of Top 10 Genes by Condition/Timepoint
- Script: `scripts/locus/deltagene-heatmaps-lineplots-bothfonts-labels.py`
- Auto-detect file(s):
  - CpG Methylation Matrix
    - Directory: output
    - Identifier: Contains the keyword "matrix"
    - Format: .xlsx or .txt
    - You had to have run Step 3 of the Preprocessing Steps (`scripts/step_3_convert_to_aberrant_signals.py`) for this file to have been generated in your output directory.
  - Patient List
    - Directory: data
    - Identifier: Contains the keyword "patient"
    - Format: .xlsx
  - Gene Annotation Map
    - Directory: output
    - Identifier: Contains the keyword "cgi_map"
    - Format: .xlsx or .csv
    - You have to have run Step 5 of the Preprocessing Steps (`scripts/step_5_build_gene_cgi_map.py`) for this file to have been generated in your output directory.
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
  - Average Visualizations (Across Patients)
    - avg_methylation_heatmap.png
    - avg_methylation_lineplot.png
  - Per-Patient Visualizations and Matrices
    - heatmap_{patient}.png
    - lineplot_{patient}.png
    - methylation_matrix_{patient}.csv

### Identify Top 10 Differentially Methylated CGI Subregions and Top 10 Differentially Methylated Genes
- Script: `scripts/locus/alltimepoint-comparisons_top10dm-cgi-subregions_alldm-genes-with-multiple-affected-cgi.py`
- Auto-detect file(s):
  - Patient Files: It searches for Excel files containing the term "patient id" in the data directory.
  - Methylation Matrix Files: It searches for Excel files containing the term "matrix" in the output directory.
    - You had to have run Step 3 of the Preprocessing Steps (scripts/step_3_convert_to_aberrant_signals.py) for this file to have been generated in your output directory.
- Data Processing Steps:
  - Extracts CpG island data from the matrix file.
  - Generates metadata for each sample, including patient ID and treatment timepoint (Baseline, On-Treatment, Post-Treatment, Healthy).
  - Filters samples to exclude those labeled as "Healthy".
  - Constructs a matrix with methylation data, grouped by patient and timepoint.
  - Collapses the matrix to average methylation levels for each CpG island across patients and timepoints.
  - Calculates average changes in methylation levels between different treatment timepoints (Baseline vs Post-Treatment, Baseline vs On-Treatment, On-Treatment vs Post-Treatment).
- Data Visualization Steps:
  - Top 10 Differentially Methylated CpG Subregions:
    - Plots bar charts for the top 10 CpG islands with the highest average changes in methylation levels between treatment timepoints. Saves these plots as PNG files in the plots directory.
  - Genes with Multiple Affected CpG Islands:
    - Identifies genes with more than one affected CpG island and plots bar charts for these genes. Saves these plots as PNG files in the plots directory.
- Generated file(s): plots directory
  - top10_diff_CGIsubregions_baseline_post.png
  - multi_CpG_genes_baseline_post.png
  - top10_diff_CGIsubregions_baseline_on.png
  - multi_CpG_genes_baseline_on.png
  - top10_diff_CGIsubregions_on_post.png
  - multi_CpG_genes_on_post.png

---

## ▶️ How to Use

### Run Locally
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


---

## 📂 Project Structure

```
methylation-pipeline/
├── data/                      # Input Excel files
├── output/                    # Filtered and merged outputs
├── plots/                     # Generated plots and Excel summaries
├── scripts/
│   ├── step_1_filter_patients_local.py
│   ├── step_2_merge_filtered_files.py
│   ├── step_3_convert_to_aberrant_signals.py
│   ├── step_4_generate_gene_annotation.py
│   ├── step_5_build_gene_cgi_map.py
│   ├── global/
│   │   ├── dotplots-by-condition.py
│   │   └── lineplots-perpatient_v3.py
│   └── locus/
│       ├── intuitive_storytelling_avg_methylation_local.py
│       ├── bubbleplot_generator_v8_gridsoff.py
│       ├── bubbleplot_generator_v8.py
│       ├── deltagene-heatmaps-lineplots-bothfonts-labels.py
│       └── alltimepoint-comparisons_top10dm-cgi-subregions_alldm-genes-with-multiple-affected-cgi.py
├── .gitignore
├── CITATION.cff
├── LICENSE
├── README.md
└── requirements.txt

```

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
