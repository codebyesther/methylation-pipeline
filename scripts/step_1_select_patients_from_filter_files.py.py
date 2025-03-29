# Step 1: Filter Methylation Files by Patient ID (Local Version)

import pandas as pd
import os

# Set input and output directories
input_dir = "data"
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

# Initialize containers
methylation_dfs = {}
patient_ids = []

# Load all Excel files from input_dir
for fname in os.listdir(input_dir):
    if fname.endswith(".xlsx") or fname.endswith(".xls"):
        fpath = os.path.join(input_dir, fname)
        if "patient" in fname.lower():
            # Load patient ID file
            patient_df = pd.read_excel(fpath)
            patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()
        else:
            # Load methylation data file
            methylation_dfs[fname] = pd.read_excel(fpath)

# Filter methylation files by patient IDs
filtered_methylation_dfs = {}
for fname, df in methylation_dfs.items():
    filtered = df[
        [df.columns[0]] + 
        [col for col in df.columns[1:] if any(pid in str(col) for pid in patient_ids)]
    ]
    filtered_methylation_dfs[fname] = filtered

# Preview filtered output
for name, df in filtered_methylation_dfs.items():
    print(f"Preview of {name} (filtered):")
    print(df.head())

# Save filtered DataFrames to Excel
for name, df in filtered_methylation_dfs.items():
    output_path = os.path.join(output_dir, f"{name}_samples-of-interest.xlsx")
    df.to_excel(output_path, index=False)

print("Filtering complete. Files saved in the 'output/' directory.")
