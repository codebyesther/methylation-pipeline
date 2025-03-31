import os
import pandas as pd

# STEP1: Auto-detect input files in current directory
files_in_dir = os.listdir()
glob20_file = next((f for f in files_in_dir if "Glob20" in f and f.endswith(".xlsx")), None)
globmin80_file = next((f for f in files_in_dir if "GlobMin80" in f and f.endswith(".xlsx")), None)

if not glob20_file or not globmin80_file:
    raise FileNotFoundError("One or both input files ('Glob20_*.xlsx', 'GlobMin80_*.xlsx') not found in current directory.")

# STEP2: Load the Excel files
glob20_df = pd.read_excel(glob20_file)
globmin80_df = pd.read_excel(globmin80_file)

# STEP3: Extract the target fragment count row
target_row = "Total CpG island fragments counts for this particular spreadsheet"
glob20_row = glob20_df.set_index("Header").loc[target_row]
globmin80_row = globmin80_df.set_index("Header").loc[target_row]

# STEP4: Clean and align sample names
glob20_cleaned = glob20_row.rename(lambda x: str(x).strip())
globmin80_cleaned = globmin80_row.rename(lambda x: str(x).strip())
glob20_numeric, globmin80_numeric = glob20_cleaned.astype(float).align(globmin80_cleaned.astype(float), join="inner")

# STEP5: Compute scaled ratio
scaled_ratio = (glob20_numeric / globmin80_numeric) * 100000

# STEP6: Assemble final DataFrame
result_df = pd.DataFrame({
    "Glob20": glob20_numeric,
    "GlobMin80": globmin80_numeric,
    "Scaled Ratio (x100K)": scaled_ratio
})

# STEP7: Export results
output_file = "scaled_fragment_ratios.xlsx"
result_df.to_excel(output_file)

print(f"Saved output to {output_file}")
