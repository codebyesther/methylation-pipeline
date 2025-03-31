#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from google.colab import files
uploaded = files.upload()


# In[ ]:


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# In[13]:

class Args:
    methylation = "matrix.xlsx"
    patients = "Patient ID list.xlsx"
    outdir = "plots"
args = Args()

methylation_dfs = {}
patient_ids = []

import os
os.makedirs(args.outdir, exist_ok=True)

# Read files
if "patient" in args.patients.lower():
    patient_df = pd.read_excel(args.patients)
    patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()

df = pd.read_excel(args.methylation)
methylation_dfs[os.path.basename(args.methylation)] = df


# This script takes your Patient ID list.xlsx (has to have "patient" in filename) and Methylation data.xlsx (has to have "matrix" in filename) to plot "Top Differentially Methylated CpG Islands" and "Genes with More than One Affected CpG Island".

"""STEP1: Install Dependencies.
"""


"""STEP2: Upload Excel files (Methylation matrix data + Patient list)."""


"""STEP3: Import packages."""

import io
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os
import time
import matplotlib.pyplot as plt

"""STEP4: Extract methylation dataframes and patient IDs."""

methylation_dfs = {}
patient_ids = []

for fname, fdata in uploaded.items():
    if "patient" in fname.lower():
        patient_df = pd.read_excel(io.BytesIO(fdata))
        patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()
    else:
        df = pd.read_excel(io.BytesIO(fdata))
        methylation_dfs[fname] = df

"""STEP5: Loop through each methylation file.
1. Locate locus-level matrix starting with "CGI_chr"
2. Extract metadata
3. Filter out healthy samples
4. Format matrix and group by patient + timepoint
5. Show the top differentially methylated CpG islands
6. Show genes with more than one affected CpG island
"""

os.makedirs("plots", exist_ok=True)

for fname, df in methylation_dfs.items():
    print(f"\n=== Processing file: {fname} ===")
    base_fname = os.path.splitext(fname)[0]  # For cleaner filenames

    # Locate locus-level matrix
    start_idx = df[df.iloc[:, 0].astype(str).str.contains("CGI_chr", na=False)].index[0]
    cpg_island_df = df.iloc[start_idx:].reset_index(drop=True)
    cpg_island_df.rename(columns={cpg_island_df.columns[0]: "CpG_Island"}, inplace=True)
    cpg_island_df.dropna(how="all", subset=cpg_island_df.columns[1:], inplace=True)

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
    collapsed = matrix.groupby(axis=1, level=[0, 1]).mean()

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
    multi_cpg_genes = gene_df[gene_df["count"] > 1].sort_values("avg_delta")

    plt.figure(figsize=(10, 6))
    sns.barplot(data=multi_cpg_genes, x="avg_delta", y="Gene", color='darkblue')
    plt.xlabel("Avg Change in Methylated Fragment Count")
    plt.axvline(0, color="gray", linestyle="--")
    plt.title(f"Genes with More than One Affected CpG Island")
    plt.tight_layout()
    plot2_path = os.path.join(args.outdir, "multi_CpG_genes.png")
    plt.savefig(plot2_path)
    plt.show()

    # Download all saved plots
    for plot_path in [plot1_path, plot2_path]:
        time.sleep(2)  # slight delay to ensure Colab has time to process download
        files.download(plot_path)

