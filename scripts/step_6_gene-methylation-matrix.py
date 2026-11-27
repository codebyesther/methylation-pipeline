"""Step 6: gene x sample matrix of Glob20 fragment counts.

Input:  output/merged_output_glob20.xlsx (Step 2), or the single Step 1 Glob20
        file when there is only one run.
Output: output/gene_methylation_matrix.csv  (index "Gene", one column per sample)

Each gene's value is the sum of the Glob20 fragment counts of the CpG islands
whose name lists that gene. Genes are matched by exact symbol, so GATA2 does
not pick up GATA2-AS1's islands. A CpG island listing two genes counts toward
both. These are raw counts; the gene-rank scripts rank genes within each
sample, which does not depend on sequencing depth.
"""

import os

import pandas as pd

from common import gene_matrix

output_folder = "output"
os.makedirs(output_folder, exist_ok=True)

files = sorted(os.listdir(output_folder))
merged = [f for f in files if "merged_output_glob20" in f and f.endswith(".xlsx")]
single = [f for f in files if "Glob20" in f and "samples-of-interest" in f and f.endswith(".xlsx")]
if merged:
    src = merged[0]
elif len(single) == 1:
    src = single[0]
    print(f"[note] no merged file; using the single run {src}")
else:
    raise FileNotFoundError("No merged_output_glob20.xlsx (Step 2) or single Glob20 run (Step 1) in output/")

cpg_matrix = pd.read_excel(os.path.join(output_folder, src), index_col=0)
cpg_matrix.index = cpg_matrix.index.astype(str).str.strip()

gm = gene_matrix(cpg_matrix, how="sum")
out_path = os.path.join(output_folder, "gene_methylation_matrix.csv")
gm.to_csv(out_path)
print(f"{gm.shape[0]} genes x {gm.shape[1]} samples from {src}")
print(f"Saved gene methylation matrix to: {out_path}")
