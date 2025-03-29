# Methylation Pipeline

This repository contains Python-based tools for filtering and visualizing DNA methylation data from patient samples. It is designed for cfDNA or cancer methylation studies where preprocessing and longitudinal analysis are needed.

---

## 🔧 What It Does

### ✅ Step 1: Filter Methylation Files by Patient ID
- Loads a patient ID Excel file and multiple methylation Excel files.
- Filters each methylation dataset to include only the columns matching specified patient IDs.
- Saves the filtered outputs to new Excel files for downstream analysis.

### 🔁 Step 2 (Optional): Merge Filtered Files for Global-Scale Analysis
This step is only required for global/sample-level methylation analysis, not CGI-level analysis.

- Merges the filtered Excel files generated in Step 1
- Produces a single combined file for downstream global metrics
- Looks for all `.xlsx` files in the `output/` directory

### 📊 Step 3: Visualize Methylation Data *(coming soon)*
- Plot longitudinal methylation trends by treatment timepoint.
- Generate chromosome-level summaries and CpG island heatmaps.
- Identify genes with significant methylation changes.

---

## 📁 File Requirements

- One Excel file with patient IDs (e.g. `patient_list.xlsx`)
  - IDs must be in the first column.
- One or more methylation Excel files
  - Columns should contain patient IDs in their names.

Place all input files in a `data/` folder and results will be saved to an `output/` folder.

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
├── output/                # Filtered Excel outputs
├── scripts/
│   ├── step_1_filter_patients.py
│   └── step_2_visualization.py (coming soon)
├── README.md
├── .gitignore
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
