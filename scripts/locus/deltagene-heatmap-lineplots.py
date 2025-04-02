# Updated script: gene-heatmap-lineplots-top10.py (delta-based logic + per-patient exports + stats)

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

parser = argparse.ArgumentParser(description='Generate gene-level methylation heatmaps and line plots based on delta values.')
parser.add_argument('--output_dir', type=str, default='plots/heatmaps-lineplots', help='Directory to save plots')
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

def find_file(directory, keyword):
    files = glob.glob(os.path.join(directory, f"*{keyword}*"))
    if not files:
        raise FileNotFoundError(f"No file containing '{keyword}' found in directory '{directory}'")
    return files[0]

output_folder = 'output'
data_folder = 'data'

try:
    cpg_matrix_file = find_file(output_folder, "matrix")
    patient_list_file = find_file(data_folder, "patient")
    gene_annotation_file = find_file(output_folder, "cgi_map")
except FileNotFoundError as e:
    print(e)
    exit(1)

# Load files
cpg_matrix = pd.read_excel(cpg_matrix_file, index_col=0) if cpg_matrix_file.endswith('.xlsx') else pd.read_csv(cpg_matrix_file, sep="\t", index_col=0)
patient_df = pd.read_excel(patient_list_file)
gene_annot_raw = pd.read_excel(gene_annotation_file) if gene_annotation_file.endswith('.xlsx') else pd.read_csv(gene_annotation_file)

# Gene annotation filtering
gene_annot = gene_annot_raw[gene_annot_raw['gene_name'].notna()].copy()
gene_annot['gene_name'] = gene_annot['gene_name'].astype(str)
cpg_headers = cpg_matrix.index.astype(str).tolist()

matched = []
for _, row in gene_annot.iterrows():
    gene = row['gene_name']
    matched_cpgs = [h for h in cpg_headers if gene in h]
    for cpg in matched_cpgs:
        matched.append({'cgi_id': cpg, 'gene_name': gene})

gene_annot = pd.DataFrame(matched)
cpg_gene_counts = gene_annot['gene_name'].value_counts()
multicpg_genes = cpg_gene_counts[cpg_gene_counts > 1].index.tolist()

# Timepoint classification
def classify_timepoint(sample):
    if "Baseline" in sample:
        return "Baseline"
    elif "Off-tx" in sample:
        return "Post-Treatment"
    elif "INNOV" in sample:
        return "Healthy"
    else:
        return "On-Treatment"

sample_meta = pd.DataFrame({
    'Sample': cpg_matrix.columns,
    'Timepoint': [classify_timepoint(s) for s in cpg_matrix.columns]
})

# Helper for patient ID
patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()
def get_patient(sample):
    return next((pid for pid in patient_ids if pid in sample), None)

# Calculate deltas and t-tests for a pair of timepoints
def calculate_deltas(tp1, tp2):
    deltas = {}
    stats = []
    for gene in multicpg_genes:
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
                baseline_vals = grouped[tp1]
                post_vals = grouped[tp2]
                avg_delta = (post_vals - baseline_vals).mean()
                deltas[gene] = avg_delta
                t_stat, p_val = ttest_rel(post_vals, baseline_vals)
                stats.append({'Gene': gene, 'Delta': avg_delta, 'T-stat': t_stat, 'P-value': p_val})
    return deltas, pd.DataFrame(stats)

# Compute deltas and stats
baseline_post, stats_bp = calculate_deltas("Baseline", "Post-Treatment")
baseline_on, stats_bo = calculate_deltas("Baseline", "On-Treatment")
on_post, stats_op = calculate_deltas("On-Treatment", "Post-Treatment")

# Save all deltas and stats to CSVs
pd.DataFrame({
    'Baseline→Post-Treatment': pd.Series(baseline_post),
    'Baseline→On-Treatment': pd.Series(baseline_on),
    'On-Treatment→Post-Treatment': pd.Series(on_post),
}).to_csv(os.path.join(args.output_dir, "gene_deltas_all_comparisons.csv"))
stats_bp.to_csv(os.path.join(args.output_dir, "baseline_vs_post_ttest.csv"), index=False)
stats_bo.to_csv(os.path.join(args.output_dir, "baseline_vs_on_ttest.csv"), index=False)
stats_op.to_csv(os.path.join(args.output_dir, "on_vs_post_ttest.csv"), index=False)

# Use Baseline→Post-Treatment for top 10 selection
top_genes = pd.Series(baseline_post).abs().sort_values(ascending=False).head(10).index.tolist()

# Compute gene matrix
gene_means = {}
for gene in top_genes:
    cpgs = gene_annot[gene_annot['gene_name'] == gene]['cgi_id']
    gene_data = cpg_matrix.loc[cpg_matrix.index.intersection(cpgs)]
    gene_means[gene] = gene_data.mean(axis=0)

gene_matrix = pd.DataFrame(gene_means).T

# Save full methylation matrix of top genes
gene_matrix.to_csv(os.path.join(args.output_dir, "top10_gene_methylation_matrix.csv"))

sample_timepoints = {sample: classify_timepoint(sample) for sample in cpg_matrix.columns}
timepoint_df = pd.DataFrame.from_dict(sample_timepoints, orient='index', columns=['Timepoint'])
gene_matrix_T = gene_matrix.T
merged = gene_matrix_T.merge(timepoint_df, left_index=True, right_index=True)
avg_by_tp = merged.groupby("Timepoint").mean().T

# Heatmap
max_height = 20
fig_height = min(max_height, len(avg_by_tp))
plt.figure(figsize=(15, fig_height))
sns.heatmap(avg_by_tp, cmap="coolwarm")
plt.title("Average Methylation per Gene Across Timepoints (Top 10 by Δ Baseline→Post-Tx)")
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "avg_methylation_heatmap.png"))
plt.close()

# Line plot
tp_order = ['Healthy', 'Baseline', 'On-Treatment', 'Post-Treatment']
avg_by_tp = avg_by_tp[[tp for tp in tp_order if tp in avg_by_tp.columns]]
plt.figure(figsize=(15, 6))
for gene in avg_by_tp.index:
    plt.plot(avg_by_tp.columns, avg_by_tp.loc[gene], label=gene)
plt.title("Methylation Trends Across Timepoints (Top 10 by Δ Baseline→Post-Tx)")
plt.ylabel("Average Methylation")
plt.xticks(rotation=45)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "avg_methylation_lineplot.png"))
plt.close()

# Per-patient heatmap, line plot, and matrix export
for patient in patient_ids:
    sample_cols = [s for s in cpg_matrix.columns if patient in s]
    if not sample_cols:
        continue
    patient_gene_matrix = gene_matrix[sample_cols].T
    timepoints = [classify_timepoint(s) for s in patient_gene_matrix.index]
    patient_gene_matrix['Timepoint'] = timepoints
    grouped = patient_gene_matrix.groupby("Timepoint").mean().T
    if grouped.empty:
        continue
    fig_height = min(max_height, len(grouped))
    plt.figure(figsize=(12, fig_height))
    sns.heatmap(grouped, cmap="coolwarm")
    plt.title(f"Gene Methylation - {patient} (Top 10 by Δ Baseline→Post-Tx)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, f"heatmap_{patient}.png"))
    plt.close()

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

    # Save patient-specific methylation matrix
    patient_gene_matrix.drop(columns='Timepoint').T.to_csv(os.path.join(args.output_dir, f"methylation_matrix_{patient}.csv"))

print("✅ All top gene delta-based plots, per-patient matrices, and stats saved to:", args.output_dir)
