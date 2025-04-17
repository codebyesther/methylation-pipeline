# This script generates the gene_methylation_matrix from the merged_output_glob20.xlsx and gene_cgi_map.csv files

import os
import glob
import pandas as pd
from tqdm import tqdm

# === Settings ===
data_folder = "data"
output_folder = "output"
os.makedirs(output_folder, exist_ok=True)

# === Helper Function ===
def find_file(directory, keyword):
    files = glob.glob(os.path.join(directory, f"*{keyword}*"))
    if not files:
        raise FileNotFoundError(f"No file containing '{keyword}' found in '{directory}'")
    return files[0]

# === Locate Files ===
cpg_matrix_file = find_file(output_folder, "merged_output_glob20")
gene_annotation_file = find_file(output_folder, "cgi_map")

# === Load Files ===
cpg_matrix = pd.read_excel(cpg_matrix_file, index_col=0) if cpg_matrix_file.endswith('.xlsx') else pd.read_csv(cpg_matrix_file, sep="\t", index_col=0)
gene_annot_raw = pd.read_excel(gene_annotation_file) if gene_annotation_file.endswith('.xlsx') else pd.read_csv(gene_annotation_file)

# === Prepare Gene Annotations ===
gene_annot = gene_annot_raw[gene_annot_raw['gene_name'].notna()].copy()
gene_annot['gene_name'] = gene_annot['gene_name'].astype(str)
cpg_headers = cpg_matrix.index.astype(str).tolist()

matched = []
for _, row in tqdm(gene_annot.iterrows(), total=gene_annot.shape[0], desc="Matching CpGs"):
    gene = row['gene_name']
    matched_cpgs = [h for h in cpg_headers if gene in h]
    for cpg in matched_cpgs:
        matched.append({'cgi_id': cpg, 'gene_name': gene})

gene_annot = pd.DataFrame(matched)
all_genes = gene_annot['gene_name'].unique().tolist()

# === Build Methylation Matrix ===
gene_rows = []
for gene in tqdm(all_genes, desc="Building gene methylation matrix"):
    cpgs = gene_annot[gene_annot['gene_name'] == gene]['cgi_id']
    gene_data = cpg_matrix.loc[cpg_matrix.index.intersection(cpgs)]
    if gene_data.empty:
        continue
    summed_row = gene_data.sum(axis=0)
    summed_row.name = gene
    gene_rows.append(summed_row)

gene_methylation_matrix = pd.DataFrame(gene_rows)
gene_methylation_matrix.index.name = "Gene"
gene_methylation_matrix.columns.name = "Sample"

# === Save Output ===
out_path = os.path.join("plots", "heatmaps-lineplots")
os.makedirs(out_path, exist_ok=True)
gene_methylation_matrix.to_csv(os.path.join(out_path, "gene_methylation_matrix.csv"))
print(f"Saved gene methylation matrix to: {os.path.join(out_path, 'gene_methylation_matrix.csv')}")
