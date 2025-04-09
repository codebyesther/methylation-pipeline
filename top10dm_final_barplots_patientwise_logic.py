import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
from scipy.stats import ttest_rel
from tqdm import tqdm

parser = argparse.ArgumentParser(description='Generate gene-level methylation heatmaps and line plots based on delta values.')
parser.add_argument('--output_dir', type=str, default='plots/heatmaps-lineplots', help='Directory to save plots')
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

# Helper functions
def find_file(directory, keyword):
    files = glob.glob(os.path.join(directory, f"*{keyword}*"))
    if not files:
        raise FileNotFoundError(f"No file containing '{keyword}' found in directory '{directory}'")
    return files[0]

def classify_timepoint(sample):
    if "Baseline" in sample:
        return "Baseline"
    elif "Off-tx" in sample:
        return "Post-Treatment"
    elif "INNOV" in sample:
        return "Healthy"
    else:
        return "On-Treatment"

def get_patient(sample):
    return next((pid for pid in patient_ids if pid in sample), None)

# Load data
output_folder = 'output'
data_folder = 'data'
cpg_matrix_file = find_file(output_folder, "matrix")
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
multicpg_genes = cpg_gene_counts[cpg_gene_counts > 1].index.tolist()

patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()

# Delta calculations
def calculate_deltas(tp1, tp2):
    deltas, stats = {}, []
    for gene in tqdm(multicpg_genes, desc=f"Calculating deltas for {tp1} vs {tp2}"):
        cpgs = gene_annot[gene_annot['gene_name'] == gene]['cgi_id']
        gene_data = cpg_matrix.loc[cpg_matrix.index.intersection(cpgs)]
        if gene_data.empty:
            continue
        avg_gene_methylation = gene_data.mean(axis=0)
        df = avg_gene_methylation.reset_index()
        df.columns = ['Sample', 'Methylation']
        df['Timepoint'] = df['Sample'].map(classify_timepoint)
        df['Patient'] = df['Sample'].map(get_patient)
        df.dropna(subset=['Patient'], inplace=True)
        grouped = df.groupby(['Patient', 'Timepoint'])['Methylation'].mean().unstack()
        if tp1 in grouped.columns and tp2 in grouped.columns:
            grouped = grouped.dropna(subset=[tp1, tp2])
            if len(grouped) >= 2:
                avg_delta = (grouped[tp2] - grouped[tp1]).mean()
                deltas[gene] = avg_delta
                t_stat, p_val = ttest_rel(grouped[tp2], grouped[tp1])
                stats.append({'Gene': gene, 'Delta': avg_delta, 'T-stat': t_stat, 'P-value': p_val})
    return deltas, pd.DataFrame(stats)

baseline_post, stats_bp = calculate_deltas("Baseline", "Post-Treatment")
baseline_on, stats_bo = calculate_deltas("Baseline", "On-Treatment")
on_post, stats_op = calculate_deltas("On-Treatment", "Post-Treatment")

# Save stats
pd.DataFrame({
    'Baseline → Post-Treatment': pd.Series(baseline_post),
    'Baseline → On-Treatment': pd.Series(baseline_on),
    'On-Treatment → Post-Treatment': pd.Series(on_post),
}).to_csv(os.path.join(args.output_dir, "gene_deltas_all_comparisons.csv"))
stats_bp.to_csv(os.path.join(args.output_dir, "baseline_vs_post_ttest.csv"), index=False)
stats_bo.to_csv(os.path.join(args.output_dir, "baseline_vs_on_ttest.csv"), index=False)
stats_op.to_csv(os.path.join(args.output_dir, "on_vs_post_ttest.csv"), index=False)

# Gene matrix for all genes
gene_means = {}
for gene in tqdm(multicpg_genes, desc="Computing gene matrix for all genes"):
    cpgs = gene_annot[gene_annot['gene_name'] == gene]['cgi_id']
    gene_data = cpg_matrix.loc[cpg_matrix.index.intersection(cpgs)]
    if not gene_data.empty:
        gene_means[gene] = gene_data.mean(axis=0)

gene_matrix = pd.DataFrame(gene_means).T
gene_matrix.to_csv(os.path.join(args.output_dir, "full_gene_matrix.csv"))

# Per-patient gene matrix export
for patient in tqdm(patient_ids, desc="Saving per-patient matrices for all genes"):
    sample_cols = [s for s in cpg_matrix.columns if patient in s]
    if not sample_cols:
        continue
    patient_gene_matrix = gene_matrix[sample_cols].T
    patient_gene_matrix.dropna(axis=1, how='all').T.to_csv(
        os.path.join(args.output_dir, f"methylation_matrix_{patient}_all_genes.csv"))

# Combine per-patient matrices
patient_matrix_files = [os.path.join(args.output_dir, f"methylation_matrix_{patient}_all_genes.csv") for patient in patient_ids]
combined_matrix = pd.concat([pd.read_csv(file, index_col=0) for file in patient_matrix_files], axis=1)
combined_matrix.to_csv(os.path.join(args.output_dir, "combined_methylation_matrix_all_genes.csv"))

# Filter top 10 genes by delta
top_genes = pd.Series(baseline_post).abs().sort_values(ascending=False).head(10).index.tolist()
filtered_gene_matrix = gene_matrix.loc[top_genes]
filtered_gene_matrix.to_csv(os.path.join(args.output_dir, "top10_gene_matrix.csv"))

# Visualization: heatmap + lineplot
sample_timepoints = {sample: classify_timepoint(sample) for sample in cpg_matrix.columns}
timepoint_df = pd.DataFrame.from_dict(sample_timepoints, orient='index', columns=['Timepoint'])
gene_matrix_T = filtered_gene_matrix.T
merged = gene_matrix_T.merge(timepoint_df, left_index=True, right_index=True)
avg_by_tp = merged.groupby("Timepoint").mean().T

# Reorder and fix timepoint x-axis
tp_order = ['Healthy', 'Baseline', 'On-Treatment', 'Post-Treatment']
avg_by_tp = avg_by_tp[[tp for tp in tp_order if tp in avg_by_tp.columns]]
avg_by_tp.columns = pd.CategoricalIndex(avg_by_tp.columns, categories=tp_order, ordered=True)
avg_by_tp = avg_by_tp.sort_index(axis=1)

# Heatmap
fig_height = min(20, len(avg_by_tp))
plt.figure(figsize=(15, fig_height))
ax = sns.heatmap(avg_by_tp, cmap="coolwarm", cbar_kws={'label': 'Methylation', 'shrink': 1})
ax.set_xlabel("Timepoint", fontsize=14)
ax.set_ylabel("Gene", fontsize=14)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=14)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=14)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=14)
cbar.set_label("Scaled Methylation Fragment Count Ratio", fontsize=14)
plt.title("Average Methylation per Gene Across Timepoints (Top 10 by ∆ Baseline → Post-Tx)", fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "avg_methylation_heatmap_top10_genes.png"))
plt.close()

# Line plot
plt.figure(figsize=(15, 6))
for gene in avg_by_tp.index:
    plt.plot(avg_by_tp.columns, avg_by_tp.loc[gene], label=gene)
plt.title("Methylation Trends Across Timepoints (Top 10 by ∆ Baseline → Post-Tx)", fontsize=14)
plt.ylabel("Average Scaled Methylation Fragment Count Ratio", fontsize=14)
plt.xlabel("Timepoint", fontsize=14)
plt.xticks(rotation=0, fontsize=14)
plt.yticks(fontsize=14)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "avg_methylation_lineplot_top10_genes.png"))
plt.close()
