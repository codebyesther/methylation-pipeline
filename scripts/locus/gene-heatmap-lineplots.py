# gene_heatmap_lineplots_local.py

import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

# Argument parser for local execution
parser = argparse.ArgumentParser(description='Generate gene-level methylation heatmaps and line plots.')
parser.add_argument('--input_dir', type=str, required=True, help='Directory containing input files')
parser.add_argument('--output_dir', type=str, default='outputs', help='Directory to save plots')
parser.add_argument('--data_dir', type=str, default='data', help='Directory containing data files')
args = parser.parse_args()

# Create output directory if it doesn't exist
os.makedirs(args.output_dir, exist_ok=True)

# Helper to find file containing a keyword
def find_file(directory, keyword):
    return glob.glob(os.path.join(directory, f"*{keyword}*"))[0]

# Load inputs based on known filename patterns
cpg_matrix_file = find_file(args.output_dir, "matrix")
patient_list_file = find_file(args.data_dir, "patient")
gene_annotation_file = find_file(args.output_dir, "gene_annotation")

cpg_matrix = pd.read_csv(cpg_matrix_file, sep="\t", index_col=0)
patient_df = pd.read_excel(patient_list_file)
gene_annot = pd.read_csv(gene_annotation_file)

# Function to classify sample timepoints
def classify_timepoint(sample_name):
    if "Baseline" in sample_name:
        return "Baseline"
    elif "Off-tx" in sample_name:
        return "Post-Treatment"
    elif "INNOV" in sample_name:
        return "Healthy"
    else:
        return "On-Treatment"

# Annotate sample timepoints
sample_timepoints = {s: classify_timepoint(s) for s in cpg_matrix.columns}
timepoint_df = pd.DataFrame.from_dict(sample_timepoints, orient='index', columns=['Timepoint'])

# Filter for genes with multiple CpG islands
cpg_gene_counts = gene_annot['gene_name'].value_counts()
multicpg_genes = cpg_gene_counts[cpg_gene_counts > 1].index.tolist()

# Average methylation for each gene
gene_means = {}
for gene in multicpg_genes:
    cpgs = gene_annot[gene_annot['gene_name'] == gene]['cgi_id']
    gene_data = cpg_matrix.loc[cpg_matrix.index.isin(cpgs)]
    gene_means[gene] = gene_data.mean()

gene_matrix = pd.DataFrame(gene_means).T

# Plot heatmap (aggregated by timepoint)
gene_matrix_T = gene_matrix.T
merged = gene_matrix_T.merge(timepoint_df, left_index=True, right_index=True)
avg_by_tp = merged.groupby("Timepoint").mean().T

plt.figure(figsize=(15, len(avg_by_tp)))
sns.heatmap(avg_by_tp, cmap="coolwarm")
plt.title("Average Methylation per Gene Across Timepoints")
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "avg_methylation_heatmap.png"))
plt.close()

# Plot line plot (aggregated)
avg_by_tp_ordered = avg_by_tp[['Baseline', 'On-Treatment', 'Post-Treatment'] if 'Healthy' not in avg_by_tp.columns else ['Healthy', 'Baseline', 'On-Treatment', 'Post-Treatment']]

plt.figure(figsize=(15, 6))
for gene in avg_by_tp_ordered.index:
    plt.plot(avg_by_tp_ordered.columns, avg_by_tp_ordered.loc[gene], label=gene)
plt.title("Methylation Trends Across Timepoints")
plt.ylabel("Average Methylation")
plt.xticks(rotation=45)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "avg_methylation_lineplot.png"))
plt.close()

print("Plots saved to", args.output_dir)
