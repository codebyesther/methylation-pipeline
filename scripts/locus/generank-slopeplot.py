import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re

# === Settings ===
input_path = os.path.join("plots", "heatmaps-lineplots", "gene_methylation_matrix.csv")
patient_list_path = os.path.join("data", "Patient ID list fot EMseq16-18-20.xlsx")
output_dir = os.path.join("plots", "rank-slopeplot")
os.makedirs(output_dir, exist_ok=True)

# === Helper Functions ===
def classify_detailed_timepoint(sample_name):
    if "INNOV" in sample_name:
        return "Healthy"
    if "Baseline" in sample_name:
        return "Baseline"
    if "Off-tx" in sample_name or "Off-Tx" in sample_name:
        return "Off-Tx"
    match = re.search(r'C(\\d+)[^\\d]?', sample_name)
    if match:
        return f\"C{match.group(1)}\"
    return \"On-Treatment\"

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
ranks.to_csv(os.path.join(output_dir, "gene_methylation_ranks.csv"))
timepoint_map = {col: classify_detailed_timepoint(col) for col in ranks.columns}
patient_map = {col: get_patient(col, patients) for col in ranks.columns}

# === Melt and Annotate ===
melted = ranks.reset_index().melt(id_vars='Gene', var_name='Sample', value_name='Rank')
melted['Timepoint'] = melted['Sample'].map(timepoint_map)
melted['Patient'] = melted['Sample'].map(patient_map)
melted.dropna(subset=['Patient'], inplace=True)
melted['Timepoint'] = pd.Categorical(melted['Timepoint'], categories=sorted(set(timepoint_map.values()), key=sort_timepoints), ordered=True)
melted.to_csv(os.path.join(output_dir, "melted_gene_methylation_ranks.csv"))

# === Plot Per-Patient Slope Charts ===
for patient_id, subdf in melted.groupby("Patient"):
    if subdf['Timepoint'].nunique() < 2:
        continue  # skip patients with < 2 timepoints
    plt.figure(figsize=(14, 8))
    for gene, gene_df in subdf.groupby("Gene"):
        plt.plot(gene_df['Timepoint'], gene_df['Rank'], alpha=0.5, linewidth=0.7)
    plt.yscale("log")
    plt.gca().invert_yaxis()
    plt.title(f"Gene Ranking Trajectories for Patient {patient_id}")
    plt.ylabel("Gene Rank (lower = more methylated, log scale)")
    plt.xlabel("Timepoint")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"rank_slopeplot_patient_{patient_id}.png"))
    plt.close()

print(f"All per-patient rank slope plots saved to: {output_dir}")
