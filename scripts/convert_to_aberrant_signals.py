import os
import pandas as pd

# STEP1: Auto-detect input files in the "output" directory
output_dir = 'output'
files_in_dir = os.listdir(output_dir)
glob20_file = next((os.path.join(output_dir, f) for f in files_in_dir if "Glob20" in f and f.endswith(".xlsx")), None)
globmin80_file = next((os.path.join(output_dir, f) for f in files_in_dir if "GlobMin80" in f and f.endswith(".xlsx")), None)

if not glob20_file or not globmin80_file:
    raise FileNotFoundError("One or both input files ('Glob20_*.xlsx', 'GlobMin80_*.xlsx') not found in the 'output' directory.")

# STEP2: Load the Excel files
glob20_df = pd.read_excel(glob20_file)
globmin80_df = pd.read_excel(globmin80_file)

# STEP3: Extract the target rows
target_rows = glob20_df[glob20_df['Header'].str.contains("CGI_chr") | (glob20_df['Header'] == "Total CpG island fragments counts for this particular spreadsheet")]
glob20_rows = target_rows.set_index("Header").loc[:, target_rows.columns != "Header"]
target_rows = globmin80_df[globmin80_df['Header'].str.contains("CGI_chr") | (globmin80_df['Header'] == "Total CpG island fragments counts for this particular spreadsheet")]
globmin80_rows = target_rows.set_index("Header").loc[:, target_rows.columns != "Header"]

# STEP4: Clean and align sample names
glob20_cleaned = glob20_rows.rename(lambda x: str(x).strip(), axis='columns')
globmin80_cleaned = globmin80_rows.rename(lambda x: str(x).strip(), axis='columns')
glob20_numeric, globmin80_numeric = glob20_cleaned.astype(float).align(globmin80_cleaned.astype(float), join="inner")

# STEP5: Compute scaled ratio
scaled_ratio = (glob20_numeric / globmin80_numeric) * 100000

# STEP6: Assemble final DataFrame
result_df = pd.DataFrame({
    "Glob20": glob20_numeric.stack(),
    "GlobMin80": globmin80_numeric.stack(),
    "Scaled Ratio (x100K)": scaled_ratio.stack()
}).reset_index()

# STEP7: Export results
output_file = os.path.join(output_dir, "scaled_fragment_ratios.xlsx")
result_df.to_excel(output_file, index=False)

print(f"Saved output to {output_file}")
