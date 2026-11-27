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

# STEP3: Filter rows to keep only those containing "Total CpG island fragments counts for this particular spreadsheet" or "CGI_chr" in the first column
filter_condition = glob20_df.iloc[:, 0].str.contains("CGI_chr") | (glob20_df.iloc[:, 0] == "Total CpG island fragments counts for this particular spreadsheet")
glob20_filtered = glob20_df[filter_condition]
globmin80_filtered = globmin80_df[filter_condition]

# Ensure the filtered DataFrames have the same shape
assert glob20_filtered.shape == globmin80_filtered.shape, "Filtered DataFrames do not have the same shape."

# STEP4: Divide the values in the output_glob20 file by the corresponding values in the output_globmin80 file and multiply by 1000
ratio_df = (glob20_filtered.iloc[:, 1:].astype(float) / globmin80_filtered.iloc[:, 1:].astype(float)) * 1000

# STEP5: Reassemble the DataFrame with the original row labels
result_df = pd.concat([glob20_filtered.iloc[:, 0].reset_index(drop=True), ratio_df.reset_index(drop=True)], axis=1)

# STEP6: Export results
output_file = os.path.join(output_dir, "scaled_fragment_ratios_matrix.xlsx")
result_df.to_excel(output_file, index=False)

print(f"Saved output to {output_file}")
