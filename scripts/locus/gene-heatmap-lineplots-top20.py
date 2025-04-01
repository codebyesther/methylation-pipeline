
import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import argparse

# Argument parser for local execution
parser = argparse.ArgumentParser(description='Generate gene-level methylation heatmaps and line plots.')
parser.add_argument('--output_dir', type=str, default='plots/heatmaps-lineplots', help='Directory to save plots')
args = parser.parse_args()

# Create output directory if it doesn't exist
os.makedirs(args.output_dir, exist_ok=True)

# Helper to find file containing a keyword
def find_file(directory, keyword):
    files = glob.glob(os.path.join(directory, f"*{keyword}*"))
    if not files:
        raise FileNotFoundError(f"No file containing '{keyword}' found in directory '{directory}'")
    return files[0]

# Define the directories
output_folder = 'output'
data_folder = 'data'

# Load inputs based on known filename patterns
try:
    cpg_matrix_file = find_file(output_folder, "matrix")
    patient_list_file = find_file(data_folder, "patient")
    gene_annotation_file = find_file(output_folder, "cgi_map")
except FileNotFoundError as e:
    print(e)
    exit(1)

# Load cpg_matrix
if cpg_matrix_file.endswith('.csv'):
    cpg_matrix = pd.read_csv(cpg_matrix_file, sep="\t", index_col=0, encoding='ISO-8859-1')
elif cpg_matrix_file.endswith('.xlsx'):
    cpg_matrix = pd.read_excel(cpg_matrix_file, index_col=0)
else:
    raise ValueError("Unsupported methylation matrix format")

# Load patient list
patient_df = pd.read_excel(patient_list_file)

# Load gene annotation
if gene_annotation_file.endswith('.csv'):
    gene_annot_raw = pd.read_csv(gene_annotation_file, encoding='ISO-8859-1')
elif gene_annotation_file.endswith('.xlsx'):
    gene_annot_raw = pd.read_excel(gene_annotation_file)
else:
    raise ValueError("Unsupported gene annotation format")

# Match gene names from gene_annotation_file to CpG matrix headers
gene_annot = gene_annot_raw.copy()
gene_annot = gene_annot[gene_annot['gene_name'].notna()]
gene_annot['gene_name'] = gene_annot['gene_name'].astype(str)

cpg_headers = cpg_matrix.index.astype(str).tolist()
matched = []
for _, row in gene_annot.iterrows():
    gene = row['gene_name']
    matched_cpgs = [h for h in cpg_headers if gene in h]
    for cpg in matched_cpgs:
        matched.append({'cgi_id': cpg, 'gene_name': gene})

gene_annot = pd.DataFrame(matched)

# Classify timepoints
def classify_timepoint(sample):
    if "Baseline" in sample:
        return "Baseline"
    elif "Off-tx" in sample:
        return "Post-Treatment"
    elif "INNOV" in sample:
        return "Healthy"
    else:
        return "On-Treatment"

sample_timepoints = {sample: classify_timepoint(sample) for sample in cpg_matrix.columns}
timepoint_df = pd.DataFrame.from_dict(sample_timepoints, orient='index', columns=['Timepoint'])

# Filter for genes with multiple CpG islands
cpg_gene_counts = gene_annot['gene_name'].value_counts()
multicpg_genes = cpg_gene_counts[cpg_gene_counts > 1].index.tolist()

# Weighted average methylation per gene
gene_means = {}
for gene in multicpg_genes:
    cpgs = gene_annot[gene_annot['gene_name'] == gene]['cgi_id']
    gene_data = cpg_matrix.loc[cpg_matrix.index.intersection(cpgs)]
    if gene_data.empty:
        continue
    weighted_avg = gene_data.mean(axis=0)
    gene_means[gene] = weighted_avg

if not gene_means:
    print("No valid gene methylation data found. Check gene names and CpG headers.")
    exit(1)

gene_matrix_all = pd.DataFrame(gene_means).T

gene_matrix = gene_matrix_all.copy()

# Select top 20 genes by variance
top_genes = gene_matrix.var(axis=1).sort_values(ascending=False).head(20).index.tolist()
gene_matrix = gene_matrix.loc[top_genes]

# Aggregate by timepoint
gene_matrix_T = gene_matrix.T
merged = gene_matrix_T.merge(timepoint_df, left_index=True, right_index=True)
avg_by_tp = merged.groupby("Timepoint").mean().T

# Plot heatmap (average)
max_height = 20  # cap the figure height
fig_height = min(max_height, len(avg_by_tp))
plt.figure(figsize=(15, fig_height))
sns.heatmap(avg_by_tp, cmap="coolwarm")
plt.title("Average Methylation per Gene Across Timepoints")
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "avg_methylation_heatmap.png"))
plt.close()

# Line plot (average)
tp_order = ['Healthy', 'Baseline', 'On-Treatment', 'Post-Treatment']
avg_by_tp = avg_by_tp[[tp for tp in tp_order if tp in avg_by_tp.columns]]
plt.figure(figsize=(15, 6))
for gene in avg_by_tp.index:
    plt.plot(avg_by_tp.columns, avg_by_tp.loc[gene], label=gene)
plt.title("Methylation Trends Across Timepoints")
plt.ylabel("Average Methylation")
plt.xticks(rotation=45)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "avg_methylation_lineplot.png"))
plt.close()


# Per-patient heatmap and line plot
patients = patient_df.iloc[:, 0].dropna().unique().tolist()
for patient in patients:
    sample_cols = [s for s in cpg_matrix.columns if patient in s]
    if not sample_cols:
        continue
    patient_gene_matrix_full = gene_matrix_all[sample_cols].T
    timepoints = [classify_timepoint(s) for s in patient_gene_matrix_full.index]
    patient_gene_matrix_full['Timepoint'] = timepoints
    grouped_full = patient_gene_matrix_full.groupby("Timepoint").mean().T

    if grouped_full.empty:
        continue

    # Select top 20 genes by variance for this patient
    top_genes = grouped_full.var(axis=1).sort_values(ascending=False).head(20).index.tolist()
    grouped = grouped_full.loc[top_genes]

    # Track color scale
    vmin = grouped.min().min()
    vmax = grouped.max().max()

    # Heatmap
    fig_height = min(max_height, len(grouped))
    plt.figure(figsize=(12, fig_height))
    ax = sns.heatmap(grouped, cmap="coolwarm", vmin=vmin, vmax=vmax, cbar_kws={"label": "Methylation"})
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([vmin, vmax])
    colorbar.set_ticklabels([f"Min: {vmin:.2f}", f"Max: {vmax:.2f}"])
    plt.title(f"Gene Methylation - {patient}")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, f"heatmap_{patient}.png"))
    plt.close()

    # Line plot
    grouped = grouped[[tp for tp in tp_order if tp in grouped.columns]]
    plt.figure(figsize=(12, 5))
    for gene in grouped.index:
        plt.plot(grouped.columns, grouped.loc[gene], label=gene)
    plt.title(f"Methylation Trends - {patient}")
    plt.ylabel("Methylation")
    plt.xticks(rotation=45)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, f"lineplot_{patient}.png"))
    plt.close()

print("✅ All plots saved to:", args.output_dir)


# Save final gene_matrix as CSV
gene_matrix.to_csv(os.path.join("output", "top20_gene_matrix.csv"))
