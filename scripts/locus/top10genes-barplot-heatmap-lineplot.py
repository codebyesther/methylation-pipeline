import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
from tqdm import tqdm
import zipfile
import re

# === Argument Parser ===
parser = argparse.ArgumentParser(description='Generate gene-level methylation barplots, heatmaps, and line plots based on delta values.')
parser.add_argument('--output_dir', type=str, default='plots/heatmaps-lineplots', help='Directory to save plots')
args = parser.parse_args()
os.makedirs(args.output_dir, exist_ok=True)

# === Helper Functions ===
def find_file(directory, keyword):
    files = glob.glob(os.path.join(directory, f"*{keyword}*"))
    if not files:
        raise FileNotFoundError(f"No file containing '{keyword}' found in directory '{directory}'")
    return files[0]

def classify_detailed_timepoint(sample):
    if match := re.search(r'(Baseline|Off-tx|INNOV|C\d+)', sample):
        token = match.group(1)
        if token == "INNOV":
            return "Healthy"
        elif token == "Off-tx":
            return "Off-Tx"
        else:
            return token  # Baseline or C1-C9
    return "On-Treatment"

def get_patient(sample):
    return next((pid for pid in patient_ids if pid in sample), None)

def sort_timepoints(tp):
    if tp == "Healthy": return -2
    if tp == "Baseline": return -1
    if tp == "Off-Tx": return 99
    if tp.startswith("C") and tp[1:].isdigit(): return int(tp[1:])
    return 98

# === Load Data ===
output_folder = 'output'
data_folder = 'data'
cpg_matrix_file = find_file(output_folder, "merged_output_glob20")
patient_list_file = find_file(data_folder, "patient")
gene_annotation_file = find_file(output_folder, "cgi_map")

cpg_matrix = pd.read_excel(cpg_matrix_file, index_col=0) if cpg_matrix_file.endswith('.xlsx') else pd.read_csv(cpg_matrix_file, sep="\t", index_col=0)
patient_df = pd.read_excel(patient_list_file)
gene_annot_raw = pd.read_excel(gene_annotation_file) if gene_annotation_file.endswith('.xlsx') else pd.read_csv(gene_annotation_file)

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
cpg_gene_counts = gene_annot['gene_name'].value_counts()
all_genes = cpg_gene_counts.index.tolist()

patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()

# === Generate Gene Methylation Matrix (Raw Fragment Counts) ===
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
gene_methylation_matrix.to_csv(os.path.join(args.output_dir, "gene_methylation_matrix.csv"))
print(f"Gene methylation matrix saved to {os.path.join(args.output_dir, 'gene_methylation_matrix.csv')}")

# === Plot Top Genes Heatmap Across Detailed Timepoints ===
baseline_means = gene_methylation_matrix.filter(like="Baseline").mean(axis=1)
top_genes = baseline_means.abs().sort_values(ascending=False).head(10).index.tolist()
ordered_top_genes = top_genes

heatmap_df = gene_methylation_matrix.loc[ordered_top_genes]
column_meta = pd.DataFrame({
    'Sample': heatmap_df.columns,
    'Timepoint': [classify_detailed_timepoint(col) for col in heatmap_df.columns]
})
avg_methylation = heatmap_df.T.join(column_meta.set_index("Sample")).groupby("Timepoint").mean().T
avg_methylation = avg_methylation[sorted(avg_methylation.columns, key=sort_timepoints)]

plt.figure(figsize=(12, 6))
sns.heatmap(avg_methylation, cmap="coolwarm", linewidths=0.5,
            cbar_kws={"label": "Scaled Methylated Fragment Count Ratio"})
plt.title("Average Methylation per Gene Across Treatment Cycles (Top 10 Genes)", fontsize=12)
plt.xlabel("Timepoint")
plt.ylabel("Gene")
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "heatmap_top10_genes_detailed_timepoints.png"))
plt.close()
print("Heatmap saved.")

# === Generate Multi-Gene Line Plots per Detailed Timepoint ===
lineplot_dir = os.path.join(args.output_dir, "gene-lineplots")
os.makedirs(lineplot_dir, exist_ok=True)

melted = gene_methylation_matrix.loc[ordered_top_genes].reset_index().melt(id_vars='Gene', var_name='Sample', value_name='Methylation')
melted['Timepoint'] = melted['Sample'].map(classify_detailed_timepoint)
melted['Patient'] = melted['Sample'].map(get_patient)
melted.dropna(subset=['Patient'], inplace=True)
melted['Timepoint'] = pd.Categorical(melted['Timepoint'], categories=sorted(melted['Timepoint'].unique(), key=sort_timepoints), ordered=True)

# Average Line Plot Across Patients
avg = melted.groupby(['Timepoint', 'Gene'])['Methylation'].mean().reset_index()
plt.figure(figsize=(12, 6))
sns.lineplot(data=avg, x='Timepoint', y='Methylation', hue='Gene', marker='o')
plt.title("Avg Methylation per Gene Across Detailed Timepoints")
plt.tight_layout()
plt.savefig(os.path.join(lineplot_dir, "multi_gene_avg_lineplot_detailed.png"))
plt.close()

# Line Plot Per Patient
for patient_id, subdf in melted.groupby('Patient'):
    plt.figure(figsize=(12, 6))
    sns.lineplot(data=subdf, x='Timepoint', y='Methylation', hue='Gene', marker='o')
    plt.title(f"Patient {patient_id} - Methylation per Gene (Detailed Timepoints)")
    plt.tight_layout()
    plt.savefig(os.path.join(lineplot_dir, f"lineplot_patient_{patient_id}_detailed.png"))
    plt.close()

print(f"Line plots saved to {lineplot_dir}")
