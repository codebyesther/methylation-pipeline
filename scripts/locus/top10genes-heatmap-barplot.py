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
import zipfile

# === Argument Parser ===
parser = argparse.ArgumentParser(description='Generate gene-level methylation barplots and heatmaps based on delta values.')
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
    available_genes = [gene for gene in top_gene_deltas.index if gene in patient_deltas_df.index]
    top_patient_deltas = patient_deltas_df.loc[available_genes]
    top_gene_deltas = top_gene_deltas.loc[available_genes]

    comparison_name = comparison.replace(' ', '_').replace('→', 'to')
    top_gene_deltas.to_csv(os.path.join(args.output_dir, f"top10_gene_deltas_{comparison_name}.csv"))
    top_patient_deltas.to_csv(os.path.join(args.output_dir, f"top10_patient_deltas_{comparison_name}.csv"))

    patients_used = top_patient_deltas.dropna(how='all', axis=1).shape[1]

    plt.figure(figsize=(10, 6))
    sns.barplot(x=top_gene_deltas.values, y=top_gene_deltas.index, color="darkblue")
    plt.axvline(0, color="gray", linestyle="--")
    plt.xlabel(f"Avg Change in Methylation ({comparison.split(' → ')[1]} − {comparison.split(' → ')[0]})", fontsize=14)
    plt.ylabel("Gene", fontsize=14)
    plt.title(f"Top 10 Genes by Avg Methylation Change ({comparison})\n(n = {patients_used} patients)", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, f"barplot_top10_gene_deltas_{comparison_name}.png"))
    plt.close()

# === Generate Heatmap of Raw Methylation for Top Genes ===
ordered_top_genes = pd.Series(baseline_post).reindex(top_genes).fillna(0).sort_values().index.tolist()
heatmap_df = gene_methylation_matrix.loc[ordered_top_genes]

column_meta = pd.DataFrame({
    'Sample': heatmap_df.columns,
    'Timepoint': [classify_timepoint(col) for col in heatmap_df.columns]
})
avg_methylation = heatmap_df.T.join(column_meta.set_index("Sample")).groupby("Timepoint").mean().T

timepoint_order = ["Healthy", "Baseline", "On-Treatment", "Post-Treatment"]
avg_methylation = avg_methylation[[tp for tp in timepoint_order if tp in avg_methylation.columns]]

plt.figure(figsize=(10, 6))
sns.heatmap(avg_methylation, cmap="coolwarm", linewidths=0.5,
            cbar_kws={"label": "Scaled Methylated Fragment Count Ratio"})
plt.title("Average Methylation per Gene Across Timepoints (Top 10 by Δ Baseline → Post-Tx)", fontsize=12)
plt.xlabel("Timepoint")
plt.ylabel("Gene")
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "heatmap_top10_genes_avg_methylation.png"))
plt.close()
print("Heatmap saved.")

# === Generate Per-Patient Heatmaps for Top 10 Genes ===
per_patient_output_dir = os.path.join(args.output_dir, "per_patient_heatmaps")
os.makedirs(per_patient_output_dir, exist_ok=True)

sample_metadata = pd.DataFrame({
    'Sample': gene_methylation_matrix.columns,
    'Timepoint': [classify_timepoint(col) for col in gene_methylation_matrix.columns],
    'Patient': [get_patient(col) for col in gene_methylation_matrix.columns]
}).dropna(subset=["Patient"])

for patient_id in sample_metadata["Patient"].unique():
    patient_samples = sample_metadata[sample_metadata["Patient"] == patient_id]
    patient_cols = patient_samples["Sample"].tolist()

    if not patient_cols:
        continue

    patient_data = gene_methylation_matrix.loc[
        gene_methylation_matrix.index.intersection(ordered_top_genes),
        gene_methylation_matrix.columns.intersection(patient_cols)
    ]
    if patient_data.empty:
        continue

    patient_data = patient_data[
        sorted(patient_data.columns, key=lambda x: timepoint_order.index(classify_timepoint(x)) if classify_timepoint(x) in timepoint_order else 99)
    ]
    patient_data = patient_data.loc[pd.Index(ordered_top_genes).intersection(patient_data.index)]

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(patient_data, cmap="coolwarm", linewidths=0.5,
                cbar_kws={"label": "Methylated Fragment Count"}, ax=ax)
    ax.set_title(f"Gene Methylation for Patient {patient_id}", fontsize=12)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Gene")
    plt.tight_layout()
    fig_path = os.path.join(per_patient_output_dir, f"heatmap_top10_genes_patient_{patient_id}.png")
    plt.savefig(fig_path)
    plt.close()

zip_path = os.path.join(args.output_dir, "per_patient_heatmaps.zip")
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
    for root, _, files in os.walk(per_patient_output_dir):
        for file in files:
            if file.endswith(".png"):
                zipf.write(os.path.join(root, file), arcname=file)

print(f"Per-patient heatmaps saved to: {per_patient_output_dir}")
print(f"All heatmaps zipped at: {zip_path}")
