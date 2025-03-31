#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob

class Args:
    patients = ""
    methylation = ""
    outdir = "plots"
args = Args()

# Locate the files based on the specified criteria
data_folder = "data"
output_folder = "output"

patient_files = glob.glob(os.path.join(data_folder, '*[pP][aA][tT][iI][eE][nN][tT] [iI][dD]*.xlsx'))
matrix_files = glob.glob(os.path.join(output_folder, '*[mM][aA][tT][rR][iI][xX]*.xlsx'))

if patient_files:
    args.patients = patient_files[0]
else:
    raise FileNotFoundError("No patient ID file found in the data folder.")

if matrix_files:
    args.methylation = matrix_files[0]
else:
    raise FileNotFoundError("No matrix file found in the output folder.")

methylation_dfs = {}
patient_ids = []

os.makedirs(args.outdir, exist_ok=True)

# Read files
if "patient" in args.patients.lower():
    patient_df = pd.read_excel(args.patients)
    patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()

df = pd.read_excel(args.methylation)
methylation_dfs[os.path.basename(args.methylation)] = df

# Metadata functions
def get_patient(sample):
    for pid in patient_ids:
        if pid in sample:
            return pid
    return None

def normalize_timepoint(sample):
    if "Baseline" in sample:
        return "Baseline"
    elif "Off-tx" in sample:
        return "Post-Treatment"
    elif "INNOV" in sample:
        return "Healthy"
    else:
        return "On-Treatment"

# Process the methylation data
for fname, df in methylation_dfs.items():
    print(f"\n=== Processing file: {fname} ===")
    base_fname = os.path.splitext(fname)[0]  # For cleaner filenames

    # Locate locus-level matrix
    start_idx = df[df.iloc[:, 0].astype(str).str.contains("CGI_chr", na=False)].index[0]
    cpg_island_df = df.iloc[start_idx:].reset_index(drop=True)
    cpg_island_df.rename(columns={cpg_island_df.columns[0]: "CpG_Island"}, inplace=True)
    cpg_island_df.dropna(how="all", subset=cpg_island_df.columns[1:], inplace=True)

    samples = cpg_island_df.columns[1:]
    sample_meta = pd.DataFrame({
        "Sample": samples,
        "Patient": [get_patient(s) for s in samples],
        "Timepoint": [normalize_timepoint(s) for s in samples]
    })
    valid_samples = sample_meta.dropna()
    valid_samples = valid_samples[valid_samples["Timepoint"] != "Healthy"]

    matrix = cpg_island_df.set_index("CpG_Island")[valid_samples["Sample"].tolist()]
    matrix.columns = pd.MultiIndex.from_frame(valid_samples[["Patient", "Timepoint"]])
    matrix = matrix.apply(pd.to_numeric, errors='coerce')
    collapsed = matrix.T.groupby(level=[0, 1]).mean().T

    # === Top10 Differentially Methylated Subregions of CpG Islands ===
    top_changes = []
    for cpg in collapsed.index:
        deltas = []
        for patient in collapsed.columns.levels[0]:
            if (patient, "Baseline") in collapsed.columns and (patient, "Post-Treatment") in collapsed.columns:
                b = collapsed.loc[cpg, (patient, "Baseline")]
                p = collapsed.loc[cpg, (patient, "Post-Treatment")]
                if not pd.isna(b) and not pd.isna(p):
                    deltas.append(p - b)
        if len(deltas) >= 2:
            top_changes.append((cpg, np.mean(deltas), len(deltas)))

    top_df = pd.DataFrame(top_changes, columns=["CpG_Island", "Avg_Delta", "n"]).sort_values("Avg_Delta")

    plt.figure(figsize=(10, 6))
    sns.barplot(data=top_df.head(10), x="Avg_Delta", y="CpG_Island", color='darkblue')
    plt.xlabel("Avg Change in Methylated Fragment Count")
    plt.axvline(0, color="gray", linestyle="--")
    plt.title(f"Top10 Differentially Methylated Subregions of CpG Islands")
    plt.tight_layout()
    plot1_path = os.path.join(args.outdir, "top10_diff_CGIsubregions.png")
    plt.savefig(plot1_path)
    plt.show()

    # === Genes with More than One Affected CpG Island ===
    top_df["Gene"] = top_df["CpG_Island"].str.extract(r"chr\w+_\d+_\d+_(.+?)_")
    gene_df = top_df.groupby("Gene").agg(count=("CpG_Island", "count"), avg_delta=("Avg_Delta", "mean")).reset_index()
    
    # Debugging statements
    print("Gene DataFrame:")
    print(gene_df)
    print("Multi CpG Genes DataFrame:")
    multi_cpg_genes = gene_df[gene_df["count"] > 1].sort_values("avg_delta")
    print(multi_cpg_genes)

    # Modify this part to plot only the top 10 genes
    plt.figure(figsize=(10, 6))
    sns.barplot(data=multi_cpg_genes.head(10).sort_values("avg_delta"), x="avg_delta", y="Gene", color='darkblue')
    plt.xlabel("Avg Change in Methylated Fragment Count")
    plt.axvline(0, color="gray", linestyle="--")
    plt.title(f"Top 10 Differentially Methylated Genes with Affected CpG Islands")
    plt.tight_layout()
    plot2_path = os.path.join(args.outdir, "top10_CpG_genes.png")
    plt.savefig(plot2_path)
    plt.show()
