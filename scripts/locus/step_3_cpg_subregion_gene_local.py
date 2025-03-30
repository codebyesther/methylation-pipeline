import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description="CGI and Gene-level methylation visualization.")
parser.add_argument("--methylation", required=True, help="Path to the methylation Excel file")
parser.add_argument("--patients", required=True, help="Path to the patient IDs Excel file")
parser.add_argument("--outdir", default="plots", help="Output directory for plots")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

patient_df = pd.read_excel(args.patients)
patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()

df = pd.read_excel(args.methylation)

start_idx = df[df.iloc[:, 0].astype(str).str.contains("CGI_chr", na=False)].index[0]
cpg_island_df = df.iloc[start_idx:].reset_index(drop=True)
cpg_island_df.rename(columns={cpg_island_df.columns[0]: "CpG_Island"}, inplace=True)
cpg_island_df.dropna(how="all", subset=cpg_island_df.columns[1:], inplace=True)

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
plt.close()

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
plt.close()