# 🧬 Methylation Pipeline

This repository contains preprocessing and visualization tools for analyzing DNA methylation changes across treatment timepoints using filtered EMseq data.

## 🔧 Preprocessing Steps

These scripts help prepare your input methylation data.

### Step 1: Filter Samples of Interest
- Script: `scripts/step_1_filter_patients_local.py`
- Filters large Excel files to retain only samples from selected patients.
- Outputs one filtered Excel file per original input.

### Step 2 (Optional): Merge Filtered Files
- Script: `scripts/step_2_merge_filtered_files.py`
- Merges filtered Excel files from separate EMseq runs into a single dataset.
- Only needed if your filtered files are split by run.

## 📊 Locus-Level (CpG Island / Gene) Analysis

This workflow zooms into locus-specific (CpG island or gene-level) methylation dynamics.

> (Coming soon)

## 🌐 Global-Level (Patient / Timepoint) Analysis

This workflow focuses on overall methylation trends per patient or treatment condition.

### Step 3: Visualize Chromosome-Level Methylation Change
- Script: `scripts/global/step_3_chr_avg_overlay_global_local.py`
- Visualizes average methylation change across chromosomes using bar plots (Baseline → On-Tx, On-Tx → Post-Tx) and a line plot (Baseline → Post-Tx).
- Plots are ordered by chromosome number (karyotype order) and designed for intuitive interpretation.

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
├── README.md
├── .gitignore
├── requirements.txt
├── CITATION.cff
```

## 📦 Dependencies

Install all requirements at once:

```bash
pip install -r requirements.txt
```

- `pandas`
- `openpyxl`
- `matplotlib`
- `seaborn`

## 📜 License

This project is licensed under the MIT License.

## 🧾 Citation

Please cite this work if you use any scripts or plots from this repository.
