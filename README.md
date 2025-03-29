# Methylation Pipeline

This repository contains Python-based tools for filtering and visualizing DNA methylation data from patient samples. It is designed for cfDNA or cancer methylation studies where preprocessing and longitudinal analysis are needed.

---

## 🔧 What It Does

### ✅ Step 1: Filter Methylation Files by Patient ID
- Loads a patient ID Excel file and multiple methylation Excel files.
- Filters each methylation dataset to include only the columns matching specified patient IDs.
- Saves the filtered outputs to new Excel files for downstream analysis.

### 🔁 Step 2 (Optional): Merge Filtered Files from Multiple EMseq Runs
This step is optional for both global and CGI-level methylation analysis.

Use it when you have separate filtered Excel files from **multiple EMseq runs** and want to combine them into a single file for unified downstream analysis.

- Merges filtered `.xlsx` files from the `output/` directory
- Produces a single combined file: `merged_output.xlsx`
- Useful when samples were processed across different EMseq runs

### 📊 Step 3: Visualize Methylation Data *(coming soon)*
- Plot longitudinal methylation trends by treatment timepoint.
- Generate chromosome-level summaries and CpG island heatmaps.
- Identify genes with significant methylation changes.

---

## 📁 Input File Requirements

- One Excel file with patient IDs (e.g. `patient_list.xlsx`)
  - IDs must be in the first column.
- One or more methylation Excel files
  - Columns should contain patient IDs in their names.

These input files are used at the start of Step 1. Place them in a `data/` folder. Filtered outputs will be saved to an `output/` folder.

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
  pip install pandas openpyxl
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

Install required Python packages with:

```bash
pip install -r requirements.txt
```

## 📂 Project Structure

```
methylation-pipeline/
├── data/                  # Input Excel files
├── output/                # Filtered and merged outputs
├── scripts/
│   ├── step_1_filter_patients_local.py
│   └── step_2_merge_filtered_files.py
├── README.md
├── .gitignore
├── requirements.txt
├── CITATION.cff
methylation-pipeline/
├── data/                  # Input Excel files
├── output/                # Filtered and merged outputs
├── scripts/
│   ├── step_1_filter_patients_local.py
│   └── step_2_merge_filtered_files.py
├── README.md
├── .gitignore
├── requirements.txt
├── CITATION.cff
```

---

## ✍️ Author

Developed by **Esther Choi**  
If you use or adapt this code in a publication or project, **please include appropriate attribution or authorship**.  

Suggested citation:  
Esther Choi (2025), *Methylation Pipeline*, GitHub: [https://github.com/codebyesther/methylation-pipeline](https://github.com/codebyesther/methylation-pipeline)

You can also find formal citation formats in the `CITATION.cff` file or by clicking **"Cite this repository"** on GitHub.

---

## 📜 License

MIT License
