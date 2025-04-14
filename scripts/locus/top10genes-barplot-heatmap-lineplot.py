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

# === (Remaining content stays the same until the end of heatmap section) ===

# === Generate Multi-Gene Line Plots ===
lineplot_dir = os.path.join(args.output_dir, "lineplots")
os.makedirs(lineplot_dir, exist_ok=True)

# Melt matrix for top genes only

# Ensure timepoints are in correct order
melted_timepoint_order = ["Healthy", "Baseline", "On-Treatment", "Post-Treatment"]
melted = gene_methylation_matrix.loc[ordered_top_genes].reset_index().melt(id_vars='Gene', var_name='Sample', value_name='Methylation')
melted['Timepoint'] = melted['Sample'].map(classify_timepoint)
melted['Patient'] = melted['Sample'].map(get_patient)
melted.dropna(subset=['Patient'], inplace=True)

# Apply categorical order to Timepoint
melted['Timepoint'] = pd.Categorical(melted['Timepoint'], categories=melted_timepoint_order, ordered=True)

# Average Line Plot Across Patients
avg = melted.groupby(['Timepoint', 'Gene'])['Methylation'].mean().reset_index()
plt.figure(figsize=(10, 6))
sns.lineplot(data=avg, x='Timepoint', y='Methylation', hue='Gene', marker='o')
plt.title("Avg Methylation per Gene Across Timepoints")
plt.tight_layout()
plt.savefig(os.path.join(lineplot_dir, "multi_gene_avg_lineplot.png"))
plt.close()

# Line Plot Per Patient
for patient_id, subdf in melted.groupby('Patient'):
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=subdf, x='Timepoint', y='Methylation', hue='Gene', marker='o')
    plt.title(f"Patient {patient_id} - Methylation per Gene")
    plt.tight_layout()
    plt.savefig(os.path.join(lineplot_dir, f"lineplot_patient_{patient_id}.png"))
    plt.close()

print(f"Line plots saved to {lineplot_dir}")
