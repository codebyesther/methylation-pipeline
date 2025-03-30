# 🧬 Methylation Pipeline

This repository contains preprocessing and visualization tools for analyzing DNA methylation changes across treatment timepoints using filtered EMseq data. It is useful for following longitudinal studies.

## 🔧 Preprocessing Steps

These scripts help prepare your input methylation data.

### Step 1: Filter Methylation Files by Patient ID
- Script: `scripts/step_1_filter_patients_local.py`
- Applies to all methylation Excel files and a patient ID list.
- Filters large Excel files to retain only samples from selected patients.
- Outputs one filtered Excel file per original input.

#### 📁 Input File Requirements
- One Excel file with patient IDs (e.g. `patient_list.xlsx`)
  - IDs must be in the first column.
- One or more methylation Excel files
  - Columns should contain patient IDs in their names.

These input files are used at the start of Step 1. Place them in a `data/` folder. Filtered outputs will be saved to an `output/` folder and will be needed for downstream analysis.

### Step 2 (Optional): Merge Filtered Files from Multiple EMseq Runs
- Script: `scripts/step_2_merge_filtered_files.py`
- Required only if you have `.xlsx` outputs from multiple different EMseq batches.
- Produces one merged file for unified analysis.

## 🔍 Locus-Level (CpG Island / Gene) Analysis

This workflow zooms into locus-specific (CpG island or gene-level) methylation dynamics.

### Step 3 (Update Coming): Visualize Methylation by Region
- Script: `scripts/locus/step_3_cpg_subregion_gene_plots_darkblue_local.py`
- Generate plots per CpG island, gene, or genomic region.
- Includes bar plots, heatmaps, line plots, and longitudinal summaries.

## 🌐 Global-Level (Patient / Timepoint) Analysis

This workflow focuses on overall methylation trends per patient or treatment condition.

### Step 3 (Update Coming): Visualize Chromosome-Level Methylation Change
- Script: `scripts/global/step_3_chr_avg_overlay_global_local.py`
- Visualizes average methylation change across chromosomes using bar plots (Baseline → On-Tx, On-Tx → Post-Tx) and a line plot (Baseline → Post-Tx).
- Plots are ordered by chromosome number (karyotype order) and designed for intuitive interpretation.

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
  pip install pandas numpy openpyxl matplotlib seaborn
  ```
- Place your input `.xlsx` files in a `data/` folder.
- Run the script:
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
