import os
import pandas as pd

# STEP1: Auto-detect input files in the "output" directory
output_dir = 'output'
files_in_dir = os.listdir(output_dir)
glob20_file = next((os.path.join(output_dir, f) for f in files_in_dir if "output_glob20" in f and f.endswith(".xlsx")), None)
globmin80_file = next((os.path.join(output_dir, f) for f in files_in_dir if "output_globmin80" in f and f.endswith(".xlsx")), None)

if not glob20_file or not globmin80_file:
    raise FileNotFoundError("One or both input files ('output_glob20_*.xlsx', 'output_globmin80_*.xlsx') not found in the 'output' directory.")

# STEP2: Load the Excel files
glob20_df = pd.read_excel(glob20_file)
globmin80_df = pd.read_excel(globmin80_file)

# STEP3: Filter rows independently for each DataFrame
filter_condition_glob20 = glob20_df.iloc[:, 0].str.contains("CGI_chr") | (glob20_df.iloc[:, 0] == "Total CpG island fragments counts for this particular spreadsheet")
filtered_glob20 = glob20_df[filter_condition_glob20]

filter_condition_globmin80 = globmin80_df.iloc[:, 0].str.contains("CGI_chr") | (globmin80_df.iloc[:, 0] == "Total CpG island fragments counts for this particular spreadsheet")
filtered_globmin80 = globmin80_df[filter_condition_globmin80]

# STEP4: Align rows by their first column (ensure they have the same labels)
aligned_glob20 = filtered_glob20[filtered_glob20.iloc[:, 0].isin(filtered_globmin80.iloc[:, 0])]
aligned_globmin80 = filtered_globmin80[filtered_globmin80.iloc[:, 0].isin(filtered_glob20.iloc[:, 0])]

# Debugging: Print shapes to ensure alignment
print("Shape of aligned_glob20:", aligned_glob20.shape)
print("Shape of aligned_globmin80:", aligned_globmin80.shape)

# Ensure the aligned DataFrames have the same shape
assert aligned_glob20.shape == aligned_globmin80.shape, "Aligned DataFrames do not have the same shape."

# STEP5: Divide the values in the output_glob20 file by the corresponding values in the output_globmin80 file and multiply by 1000
ratio_df = (aligned_glob20.iloc[:, 1:].astype(float).reset_index(drop=True) /
            aligned_globmin80.iloc[:, 1:].astype(float).reset_index(drop=True)) * 1000

# STEP6: Reassemble the DataFrame with the original row labels
result_df = pd.concat([aligned_glob20.iloc[:, 0].reset_index(drop=True), ratio_df], axis=1)

# STEP7: Add specific rows from glob20_df and globmin80_df
# Extract the "Total CpG island fragments counts for this particular spreadsheet" row from both DataFrames
total_glob20_row = glob20_df[glob20_df.iloc[:, 0] == "Total CpG island fragments counts for this particular spreadsheet"]
total_globmin80_row = globmin80_df[globmin80_df.iloc[:, 0] == "Total CpG island fragments counts for this particular spreadsheet"]

# Rename the rows
if not total_glob20_row.empty:
    total_glob20_row.iloc[0, 0] = "Total CpG island fragments counts for Glob20"
if not total_globmin80_row.empty:
    total_globmin80_row.iloc[0, 0] = "Total CpG island fragments counts for GlobMin80"

# Combine the renamed rows and append them above the original summary row
summary_rows = pd.concat([total_glob20_row, total_globmin80_row], ignore_index=True)
result_df = pd.concat([summary_rows, result_df], ignore_index=True)

# STEP8: Export results
output_file = os.path.join(output_dir, "scaled_fragment_ratios_matrix.xlsx")
result_df.to_excel(output_file, index=False)

print(f"Saved output to {output_file}")
