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

# STEP3: Filter for CGI rows only
cpg_label = "Total CpG island fragments counts for this particular spreadsheet"
filtered_glob20 = glob20_df[glob20_df.iloc[:, 0].str.startswith("CGI_")].copy()
filtered_globmin80 = globmin80_df[globmin80_df.iloc[:, 0].str.startswith("CGI_")].copy()

# STEP4: Align rows and columns
# Set row labels as index
filtered_glob20.set_index(filtered_glob20.columns[0], inplace=True)
filtered_globmin80.set_index(filtered_globmin80.columns[0], inplace=True)

# Get intersection of row and column labels
shared_rows = filtered_glob20.index.intersection(filtered_globmin80.index)
shared_columns = filtered_glob20.columns.intersection(filtered_globmin80.columns)

# Subset and sort for exact alignment
aligned_glob20 = filtered_glob20.loc[shared_rows, shared_columns].sort_index().sort_index(axis=1)
aligned_globmin80 = filtered_globmin80.loc[shared_rows, shared_columns].sort_index().sort_index(axis=1)

# STEP5: Safe division — avoid divide-by-zero and preserve alignment
with pd.option_context('mode.use_inf_as_na', True):
    ratio_df = (aligned_glob20 / aligned_globmin80.replace(0, pd.NA)) * 100000

# Reset index to turn row labels back into the first column
ratio_df.reset_index(inplace=True)

# STEP6: Add total rows from input files
total_glob20_row = glob20_df[glob20_df.iloc[:, 0] == cpg_label].copy()
total_globmin80_row = globmin80_df[globmin80_df.iloc[:, 0] == cpg_label].copy()

if not total_glob20_row.empty:
    total_glob20_row.iloc[0, 0] = "Total CpG island fragments counts for Glob20"
if not total_globmin80_row.empty:
    total_globmin80_row.iloc[0, 0] = "Total CpG island fragments counts for GlobMin80"

# Combine summary and ratio rows
result_df = pd.concat([total_glob20_row, total_globmin80_row, ratio_df], ignore_index=True)

# STEP7: Export results
output_file = os.path.join(output_dir, "scaled_fragment_ratios_matrix.xlsx")
result_df.to_excel(output_file, index=False)

print(f"✅ Saved output to: {output_file}")
print(f"📊 Output includes {len(shared_rows)} aligned CGI rows and {len(shared_columns)} shared sample columns.")
