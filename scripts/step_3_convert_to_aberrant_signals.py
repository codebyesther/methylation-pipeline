import os
import pandas as pd

# STEP1: Auto-detect input files in the "output" directory
output_dir = 'output'
files_in_dir = os.listdir(output_dir)
glob20_file = next((os.path.join(output_dir, f) for f in files_in_dir if "output_glob20" in f and f.endswith(".xlsx")), None)
globmin80_file = next((os.path.join(output_dir, f) for f in files_in_dir if "output_globmin80" in f and f.endswith(".xlsx")), None)

if not glob20_file or not globmin80_file:
    raise FileNotFoundError("One or both input files ('Glob20_*.xlsx', 'GlobMin80_*.xlsx') not found in the 'output' directory.")

# STEP2: Load the Excel files
glob20_df = pd.read_excel(glob20_file)
globmin80_df = pd.read_excel(globmin80_file)

# STEP3: Extract the target rows
target_rows_glob20 = glob20_df[glob20_df['Header'].str.contains("CGI_chr") | (glob20_df['Header'] == "Total CpG island fragments counts for this particular spreadsheet")]
glob20_rows = target_rows_glob20.set_index("Header")

target_rows_globmin80 = globmin80_df[globmin80_df['Header'].str.contains("CGI_chr") | (globmin80_df['Header'] == "Total CpG island fragments counts for this particular spreadsheet")]
globmin80_rows = target_rows_globmin80.set_index("Header")

# Ensure we only work with common columns
common_columns = glob20_rows.columns.intersection(globmin80_rows.columns)
glob20_rows = glob20_rows[common_columns]
globmin80_rows = globmin80_rows[common_columns]

# STEP4: Clean and align sample names
glob20_cleaned = glob20_rows.rename(lambda x: str(x).strip(), axis='columns')
globmin80_cleaned = globmin80_rows.rename(lambda x: str(x).strip(), axis='columns')
glob20_numeric, globmin80_numeric = glob20_cleaned.astype(float).align(globmin80_cleaned.astype(float), join="inner")

# STEP5: Compute scaled ratio
scaled_ratio = (glob20_numeric / globmin80_numeric) * 100000

# STEP6: Assemble final DataFrame
result_df = pd.concat([glob20_numeric, globmin80_numeric, scaled_ratio], keys=['Glob20', 'GlobMin80', 'Scaled Ratio (x100K)']).unstack(level=1)

# STEP7: Export results
output_file = os.path.join(output_dir, "scaled_fragment_ratios_matrix.xlsx")

# Split DataFrame into chunks if too large for a single sheet
max_columns = 16384
num_chunks = (result_df.shape[1] // max_columns) + 1

with pd.ExcelWriter(output_file) as writer:
    for i in range(num_chunks):
        start_col = i * max_columns
        end_col = min((i + 1) * max_columns, result_df.shape[1])
        result_df.iloc[:, start_col:end_col].to_excel(writer, sheet_name=f'Sheet_{i+1}')

print(f"Saved output to {output_file}")
