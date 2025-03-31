# Step 2: Merge Filtered Files (Local Version)

import pandas as pd
import os

# Set input and output paths
input_dir = "output"
output_file_glob20 = "output/merged_output_glob20.xlsx"
output_file_globmin80 = "output/merged_output_globmin80.xlsx"

# Collect all Excel files from the output folder
dfs_glob20 = []
dfs_globmin80 = []

for filename in os.listdir(input_dir):
    if filename.endswith(".xlsx"):
        fpath = os.path.join(input_dir, filename)
        if "Glob20" in filename:
            print(f"Processing {filename} for Glob20...")
            df = pd.read_excel(fpath, sheet_name="Sheet1")
            df_t = df.set_index("Header").T
            dfs_glob20.append(df_t)
        elif "GlobMin80" in filename:
            print(f"Processing {filename} for GlobMin80...")
            df = pd.read_excel(fpath, sheet_name="Sheet1")
            df_t = df.set_index("Header").T
            dfs_globmin80.append(df_t)

# Combine and reset index for Glob20 files
if dfs_glob20:
    merged_df_glob20 = pd.concat(dfs_glob20, axis=0).T.reset_index()
    print("Merge complete for Glob20. Preview:")
    print(merged_df_glob20.head())
    # Save merged file for Glob20
    merged_df_glob20.to_excel(output_file_glob20, index=False)
    print(f"Merged file saved as: {output_file_glob20}")

# Combine and reset index for GlobMin80 files
if dfs_globmin80:
    merged_df_globmin80 = pd.concat(dfs_globmin80, axis=0).T.reset_index()
    print("Merge complete for GlobMin80. Preview:")
    print(merged_df_globmin80.head())
    # Save merged file for GlobMin80
    merged_df_globmin80.to_excel(output_file_globmin80, index=False)
    print(f"Merged file saved as: {output_file_globmin80}")
