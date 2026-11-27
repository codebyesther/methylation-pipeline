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
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (detailed_timepoint, find_file, find_patient_file, gene_matrix,  # noqa: E402
                    get_patient as _get_patient, read_patient_ids, timepoint_sort_key)

# === Argument Parser ===
parser = argparse.ArgumentParser(description='Generate gene-level methylation barplots, heatmaps, and line plots based on delta values.')
parser.add_argument('--output_dir', type=str, default='plots/heatmaps-lineplots', help='Directory to save plots')
args = parser.parse_args()
os.makedirs(args.output_dir, exist_ok=True)

# === Helper Functions ===
def classify_detailed_timepoint(sample):
    return detailed_timepoint(sample)

def get_patient(sample):
    return _get_patient(sample, patient_ids)

def sort_timepoints(tp):
    return timepoint_sort_key(tp)

# === Load Data ===
output_folder = 'output'
data_folder = 'data'
cpg_matrix_file = find_file(output_folder, "merged_output_glob20", exts=(".xlsx",))
patient_list_file = find_patient_file(data_folder)

cpg_matrix = pd.read_excel(cpg_matrix_file, index_col=0)
cpg_matrix.index = cpg_matrix.index.astype(str).str.strip()
patient_ids = read_patient_ids(patient_list_file)

# === Gene x sample matrix of raw Glob20 fragment counts (exact gene matching) ===
gene_methylation_matrix = gene_matrix(cpg_matrix, how="sum")
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
            cbar_kws={"label": "Glob20 Fragment Count (summed per gene)"})
plt.title("Average Methylation per Gene Across Treatment Cycles (Top 10 Genes by Baseline Signal)", fontsize=12)
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
