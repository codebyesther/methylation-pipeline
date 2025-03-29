# Step 2: Merge Filtered Files (Local Version)

import pandas as pd
import os

# Set input and output paths
input_dir = "output"
output_file = "output/merged_output.xlsx"

# Collect all Excel files from the output folder
dfs = []

for filename in os.listdir(input_dir):
    if filename.endswith(".xlsx"):
        fpath = os.path.join(input_dir, filename)
        print(f"Processing {filename}...")
        df = pd.read_excel(fpath, sheet_name="Sheet1")
        df_t = df.set_index("Header").T
        dfs.append(df_t)

# Combine and reset index
merged_df = pd.concat(dfs, axis=0).T.reset_index()
print("Merge complete. Preview:")
print(merged_df.head())

# Save merged file
merged_df.to_excel(output_file, index=False)
print(f"Merged file saved as: {output_file}")
