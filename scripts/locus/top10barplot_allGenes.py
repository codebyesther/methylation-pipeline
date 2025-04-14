# -*- coding: utf-8 -*-
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

# === Helper Functions ===
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

# === Load Data ===
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
all_genes = cpg_gene_counts.index.tolist()

patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()

# === Delta Calculation Function ===
def calculate_deltas(tp1, tp2):
    deltas = {}
    stats = []
    gene_patient_deltas = {}

    for gene in tqdm(all_genes, desc=f"Calculating deltas for {tp1} vs {tp2}"):
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
                delta_per_patient = grouped[tp2] - grouped[tp1]
                avg_delta = delta_per_patient.mean()
                deltas[gene] = avg_delta
                gene_patient_deltas[gene] = delta_per_patient
                t_stat, p_val = ttest_rel(grouped[tp2], grouped[tp1])
                stats.append({'Gene': gene, 'Delta': avg_delta, 'T-stat': t_stat, 'P-value': p_val})
    return deltas, pd.DataFrame(stats), pd.DataFrame(gene_patient_deltas).T

# === Run Delta Comparisons ===
baseline_post, stats_bp, bp_patient_deltas = calculate_deltas("Baseline", "Post-Treatment")
baseline_on, stats_bo, bo_patient_deltas = calculate_deltas("Baseline", "On-Treatment")
on_post, stats_op, op_patient_deltas = calculate_deltas("On-Treatment", "Post-Treatment")

# === Save Delta and T-Test Results ===
pd.DataFrame({
    'Baseline → Post-Treatment': pd.Series(baseline_post),
    'Baseline → On-Treatment': pd.Series(baseline_on),
    'On-Treatment → Post-Treatment': pd.Series(on_post),
}).to_csv(os.path.join(args.output_dir, "gene_deltas_all_comparisons.csv"))

stats_bp.to_csv(os.path.join(args.output_dir, "baseline_vs_post_ttest.csv"), index=False)
stats_bo.to_csv(os.path.join(args.output_dir, "baseline_vs_on_ttest.csv"), index=False)
stats_op.to_csv(os.path.join(args.output_dir, "on_vs_post_ttest.csv"), index=False)

# === Save Per-Gene, Per-Patient Delta Tables ===
bp_patient_deltas.to_csv(os.path.join(args.output_dir, "patient_deltas_baseline_to_post.csv"))
bo_patient_deltas.to_csv(os.path.join(args.output_dir, "patient_deltas_baseline_to_on.csv"))
op_patient_deltas.to_csv(os.path.join(args.output_dir, "patient_deltas_on_to_post.csv"))

# === Generate Gene Methylation Matrix (Raw Fragment Counts) ===
gene_methylation_matrix = pd.DataFrame()

for gene in all_genes:
    cpgs = gene_annot[gene_annot['gene_name'] == gene]['cgi_id']
    gene_data = cpg_matrix.loc[cpg_matrix.index.intersection(cpgs)]
    if gene_data.empty:
        continue
    gene_methylation_matrix[gene] = gene_data.sum(axis=0)

gene_methylation_matrix = gene_methylation_matrix.T
gene_methylation_matrix.index.name = "Gene"
gene_methylation_matrix.columns.name = "Sample"

gene_methylation_matrix.to_csv(os.path.join(args.output_dir, "gene_methylation_matrix.csv"))
print(f"Gene methylation matrix saved to {os.path.join(args.output_dir, 'gene_methylation_matrix.csv')}")

# === Plot Barplots and Save Delta Tables for Top 10 Genes ===
top_genes = pd.Series(baseline_post).abs().sort_values(ascending=False).head(10).index.tolist()
comparisons = {
    "Baseline → Post-Treatment": (baseline_post, bp_patient_deltas),
    "Baseline → On-Treatment": (baseline_on, bo_patient_deltas),
    "On-Treatment → Post-Treatment": (on_post, op_patient_deltas),
}

for comparison, (delta_values, patient_deltas_df) in comparisons.items():
    delta_series = pd.Series(delta_values)
    top_gene_deltas = delta_series.reindex(top_genes).fillna(0).sort_values()
    top_patient_deltas = patient_deltas_df.loc[top_gene_deltas.index]

    # Save top 10 delta values and patient-level breakdown
    comparison_name = comparison.replace(' ', '_').replace('→', 'to')
    top_gene_deltas.to_csv(os.path.join(args.output_dir, f"top10_gene_deltas_{comparison_name}.csv"))
    top_patient_deltas.to_csv(os.path.join(args.output_dir, f"top10_patient_deltas_{comparison_name}.csv"))

    # n = number of patients used (in any of the top 10 genes for this comparison)
    patients_used = top_patient_deltas.dropna(how='all', axis=1).shape[1]

    # Plot
    plt.figure(figsize=(10, 6))
    sns.barplot(x=top_gene_deltas.values, y=top_gene_deltas.index, color="darkblue")
    plt.axvline(0, color="gray", linestyle="--")
    plt.xlabel(f"Avg Change in Methylation ({comparison.split(' → ')[1]} − {comparison.split(' → ')[0]})", fontsize=14)
    plt.ylabel("Gene", fontsize=14)
    plt.title(f"Top 10 Genes by Avg Methylation Change ({comparison})\n(n = {patients_used} patients)", fontsize=14)
    plt.tight_layout()
    filename = f"barplot_top10_gene_deltas_{comparison_name}.png"
    plt.savefig(os.path.join(args.output_dir, filename))
    plt.close()
