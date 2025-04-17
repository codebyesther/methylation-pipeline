# This version ensures replicates (same patient, same timepoint) are collapsed via averaging before plotting.

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re

# === Settings ===
input_path = os.path.join("output", "gene_methylation_matrix.csv")
patient_list_path = os.path.join("data", "Patient ID list fot EMseq16-18-20.xlsx")
output_dir = os.path.join("plots", "rank-slopeplot")
os.makedirs(output_dir, exist_ok=True)

# === Helper Functions ===
def classify_detailed_timepoint(sample_name):
    if "INNOV" in sample_name:
        return "Healthy"
    if "Baseline" in sample_name:
        return "Baseline"
    if re.search(r'Off[-_]?tx', sample_name, re.IGNORECASE):
        return "Off-Tx"
    match = re.search(r'C(\d{1,2})(?!\d)', sample_name)
    if match:
        return f"C{int(match.group(1))}"
    return "On-Treatment"

def sort_timepoints(tp):
    if tp == "Healthy": return -2
    if tp == "Baseline": return -1
    if tp == "Off-Tx": return 99
    if tp.startswith("C") and tp[1:].isdigit(): return int(tp[1:])
    return 98

def get_patient(sample_name, patient_ids):
    return next((pid for pid in patient_ids if pid in sample_name), None)

# === Load Data ===
matrix = pd.read_csv(input_path, index_col=0)
patients = pd.read_excel(patient_list_path).iloc[:, 0].dropna().astype(str).tolist()

# === Create Rankings Per Sample ===
ranks = matrix.rank(axis=0, method='min', ascending=False)
ranks[ranks >= 501] = 501
ranks.to_csv(os.path.join(output_dir, "gene_methylation_ranks.csv"))

# === Metadata Mapping ===
timepoint_map = {col: classify_detailed_timepoint(col) for col in ranks.columns}
patient_map = {col: get_patient(col, patients) for col in ranks.columns}

# === Melt and Annotate ===
melted = ranks.reset_index().melt(id_vars='Gene', var_name='Sample', value_name='Rank')
melted['Timepoint'] = melted['Sample'].map(timepoint_map)
melted['Patient'] = melted['Sample'].map(patient_map)
melted.dropna(subset=['Patient', 'Timepoint'], inplace=True)
all_timepoints_ordered = sorted(set(melted['Timepoint']), key=sort_timepoints)

# === Average over replicates ===
grouped_avg = melted.groupby(['Gene', 'Patient', 'Timepoint'], as_index=False).agg({'Rank': 'mean'})

# === Store top genes summary ===
top_genes_summary = []

# === Plot Per-Patient Slope Charts ===
for patient_id, subdf in grouped_avg.groupby("Patient"):
    if subdf['Timepoint'].nunique() < 2:
        continue

    # Determine highlight genes from all raw (non-averaged) ranks
    patient_samples = melted[melted['Patient'] == patient_id]['Sample'].unique()
    highlight_genes = set()
    for sample in patient_samples:
        top10_genes = ranks[sample].nsmallest(10).index.tolist()
        highlight_genes.update(top10_genes)

    top_genes_summary.append({
        "Patient": patient_id,
        "Top_Genes": ", ".join(sorted(highlight_genes))
    })

    # Assign colors
    highlight_genes = sorted(highlight_genes)
    color_palette = sns.color_palette("husl", len(highlight_genes))
    gene_color_dict = dict(zip(highlight_genes, color_palette))

    plt.figure(figsize=(14, 8))
    subdf = subdf.copy()
    subdf['Timepoint'] = pd.Categorical(subdf['Timepoint'], categories=all_timepoints_ordered, ordered=True)

    for gene, gene_df in subdf.groupby("Gene"):
        gene_df = gene_df.sort_values('Timepoint')
        is_highlighted = gene in gene_color_dict
        plt.plot(
            gene_df['Timepoint'], gene_df['Rank'],
            alpha=0.6,
            linewidth=2 if is_highlighted else 0.5,
            color=gene_color_dict[gene] if is_highlighted else 'black',
            linestyle='-' if is_highlighted else ':',
            label=gene if is_highlighted else None
        )

    plt.yscale("log")
    plt.gca().invert_yaxis()
    collapsed_count = (subdf['Rank'] == 501).sum()
    plt.text(0.99, 0.02, f"{collapsed_count} gene-timepoints at rank 501", ha='right', va='bottom',
             transform=plt.gca().transAxes, fontsize=10, color='gray')
    plt.title(f"Gene Ranking Trajectories for Patient {patient_id}")
    plt.ylabel("Gene Rank (smaller = more methylated, log scale)")
    plt.xlabel("Timepoint")

    handles, labels = plt.gca().get_legend_handles_labels()
    if handles:
        unique = dict(zip(labels, handles))
        plt.legend(
            unique.values(), unique.keys(),
            title="Top 10 Genes (per sample)",
            fontsize=8, title_fontsize=9,
            bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.
        )

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"rank_slopeplot_patient_{patient_id}.png"), bbox_inches='tight')
    plt.close()

# Save full melted table and top gene list
grouped_avg.to_csv(os.path.join(output_dir, "melted_gene_methylation_ranks_avg.csv"), index=False)
pd.DataFrame(top_genes_summary).to_csv(os.path.join(output_dir, "top_genes_per_patient.csv"), index=False)

print(f"Finished! Slope plots and summaries saved in: {output_dir}")

