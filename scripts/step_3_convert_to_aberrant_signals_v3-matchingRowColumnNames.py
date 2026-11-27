"""Step 3: Glob20 / GlobMin80 ratio for every CpG island and sample.

Inputs (output/):
  merged_output_glob20.xlsx and merged_output_globmin80.xlsx from Step 2, or,
  with a single sequencing run, the one pair of Step 1 files
  (*Glob20*_samples-of-interest.xlsx and *GlobMin80*_samples-of-interest.xlsx).

Output: output/scaled_fragment_ratios_matrix.xlsx, laid out as
  row 1  Scaled ratio of totals (Glob20 / GlobMin80 x 100000)  <- sample-level signal
  row 2  Total CpG island fragments counts for Glob20           (raw)
  row 3  Total CpG island fragments counts for GlobMin80        (raw)
  rows   CGI_...  Glob20 count / GlobMin80 count of the same CGI x 100000

Cells are matched by CGI name and sample name, not by position.
  * GlobMin80 = 0 and Glob20 > 0  ->  "INF" (highlighted red, logged)
  * GlobMin80 = 0 and Glob20 = 0  ->  0 (no fragments at all is not infinite)
"""

import os

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import PatternFill

from common import GLOB20_TOTAL_ROW, GLOBMIN80_TOTAL_ROW, SCALED_TOTAL_ROW, TOTAL_ROW

output_dir = "output"
SCALE = 100000


def pick(keyword_merged, keyword_run):
    files = sorted(os.listdir(output_dir))
    merged = [f for f in files if keyword_merged in f.lower() and f.endswith(".xlsx") and not f.startswith("~$")]
    if merged:
        return os.path.join(output_dir, merged[0])
    runs = [f for f in files if keyword_run in f and "samples-of-interest" in f and f.endswith(".xlsx")]
    if len(runs) == 1:
        print(f"[note] no merged file for {keyword_run}; using the single run {runs[0]}")
        return os.path.join(output_dir, runs[0])
    if len(runs) > 1:
        raise FileNotFoundError(f"{len(runs)} {keyword_run} runs found; run Step 2 to merge them first.")
    raise FileNotFoundError(f"No {keyword_run} input in output/; run Steps 1-2 first.")


# STEP 1: find inputs ("Glob20" never matches GlobMin80 file names)
glob20_file = pick("output_glob20", "Glob20")
globmin80_file = pick("output_globmin80", "GlobMin80")
print(f"Glob20:    {glob20_file}\nGlobMin80: {globmin80_file}")

# STEP 2: load
glob20_df = pd.read_excel(glob20_file)
globmin80_df = pd.read_excel(globmin80_file)
label_col = glob20_df.columns[0]
glob20_df[label_col] = glob20_df[label_col].astype(str).str.strip()
globmin80_df[globmin80_df.columns[0]] = globmin80_df[globmin80_df.columns[0]].astype(str).str.strip()
glob20_df = glob20_df.set_index(label_col)
globmin80_df = globmin80_df.set_index(globmin80_df.columns[0])

# STEP 3: CGI rows only, numeric
num = glob20_df[glob20_df.index.str.startswith("CGI_")].apply(pd.to_numeric, errors="coerce").fillna(0)
ref = globmin80_df[globmin80_df.index.str.startswith("CGI_")].apply(pd.to_numeric, errors="coerce").fillna(0)

# STEP 4: align by row and column labels
rows = num.index.intersection(ref.index)
cols = num.columns.intersection(ref.columns)
dropped = sorted(set(num.columns) ^ set(ref.columns))
if dropped:
    print(f"[warn] samples in only one of the two files, left out: {dropped}")
num = num.loc[rows, cols].sort_index().sort_index(axis=1)
ref = ref.loc[rows, cols].sort_index().sort_index(axis=1)
print(f"{len(rows)} CpG islands x {len(cols)} samples")

# STEP 5: ratio; x/0 -> INF, 0/0 -> 0
ratio = num / ref.replace(0, np.nan) * SCALE
inf_mask = (ref == 0) & (num > 0)
ratio = ratio.mask((ref == 0) & (num == 0), 0.0)
ratio_out = ratio.astype(object).where(~inf_mask, "INF")
print(f"'INF' cells (GlobMin80 = 0, Glob20 > 0): {int(inf_mask.values.sum())}")

inf_locations = inf_mask.stack()
inf_locations = inf_locations[inf_locations].reset_index().iloc[:, :2]
inf_locations.columns = ["CGI_Region", "Sample"]
inf_locations.to_csv(os.path.join(output_dir, "inf_locations_log.csv"), index=False)
inf_summary = inf_mask.sum().reset_index()
inf_summary.columns = ["Sample", "INF_Count"]
inf_summary.to_csv(os.path.join(output_dir, "inf_summary.csv"), index=False)


# STEP 6: totals. Use the lab's total row when present, else the CGI sum.
def total(df, cgi):
    if TOTAL_ROW in df.index:
        return pd.to_numeric(df.loc[TOTAL_ROW], errors="coerce").reindex(cols)
    print("[note] no total row in an input file; using the sum of its CGI rows")
    return cgi.sum()


tot20, tot80 = total(glob20_df, num), total(globmin80_df, ref)
scaled_total = tot20 / tot80.replace(0, np.nan) * SCALE

summary = pd.DataFrame([scaled_total, tot20, tot80],
                       index=[SCALED_TOTAL_ROW, GLOB20_TOTAL_ROW, GLOBMIN80_TOTAL_ROW])
result = pd.concat([summary.astype(object), ratio_out])
result.index.name = "Header"

# STEP 7: export, INF cells in red
output_file = os.path.join(output_dir, "scaled_fragment_ratios_matrix.xlsx")
result.reset_index().to_excel(output_file, index=False)
wb = load_workbook(output_file)
ws = wb.active
ws.title = "Scaled Fragment Ratios"
red = PatternFill(start_color="FF0000", end_color="FF0000", fill_type="solid")
for row in ws.iter_rows(min_row=5, min_col=2, max_col=ws.max_column):
    for cell in row:
        if cell.value == "INF":
            cell.fill = red
wb.save(output_file)
print(f"Saved {output_file}")
