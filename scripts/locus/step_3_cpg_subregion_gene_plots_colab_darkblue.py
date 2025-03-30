import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description="CpG subregion and gene-level methylation analysis with auto file detection.")
parser.add_argument("--datadir", default="data", help="Directory containing input Excel files")
parser.add_argument("--outdir", default="plots", help="Output directory for plots")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

patient_file, methylation_file = None, None
for fname in os.listdir(args.datadir):
    if fname.endswith(".xlsx"):
        if "patient" in fname.lower():
            patient_file = os.path.join(args.datadir, fname)
        else:
            methylation_file = os.path.join(args.datadir, fname)

if not patient_file or not methylation_file:
    raise FileNotFoundError("Could not find required input files in the specified data directory.")

patient_df = pd.read_excel(patient_file)
patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()

df = pd.read_excel(methylation_file)

start_idx = df[df.iloc[:, 0].astype(str).str.contains("CGI_chr", na=False)].index[0]
cpg_island_df = df.iloc[start_idx:].reset_index(drop=True)
cpg_island_df.rename(columns={cpg_island_df.columns[0]: "CpG_Island"}, inplace=True)
cpg_island_df.dropna(how="all", subset=cpg_island_df.columns[1:], inplace=True)

def get_patient(sample):
    return next((pid for pid in patient_ids if pid in sample), None)

def normalize_timepoint(sample):
    if "Baseline" in sample:
        return "Baseline"
    elif "Off-tx" in sample:
        return "Post-Treatment"
    elif "INNOV" in sample:
        return "Healthy"
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

# Top Differentially Methylated CpG Islands
for comparison in [("Baseline", "On-Treatment"), ("On-Treatment", "Post-Treatment"), ("Baseline", "Post-Treatment")]:
    top_changes = []
    for cpg in collapsed.index:
        deltas = []
        for patient in collapsed.columns.levels[0]:
            if (patient, comparison[0]) in collapsed.columns and (patient, comparison[1]) in collapsed.columns:
                t0, t1 = collapsed.loc[cpg, (patient, comparison[0])], collapsed.loc[cpg, (patient, comparison[1])]
                if not pd.isna(t0) and not pd.isna(t1):
                    deltas.append(t1 - t0)
        if len(deltas) >= 2:
            top_changes.append((cpg, np.mean(deltas), len(deltas)))

    top_df = pd.DataFrame(top_changes, columns=["CpG_Island", "Avg_Delta", "n"]).sort_values("Avg_Delta").head(10)

    plt.figure(figsize=(10, 6))
    sns.barplot(data=top_df, x="Avg_Delta", y="CpG_Island", color='#000080')
    plt.xlabel("Avg Change in Methylated Fragment Count")
    plt.axvline(0, color="gray", linestyle="--")
    plt.title(f"Top10 Differentially Methylated CpG Subregions ({comparison[0]} → {comparison[1]})")
    plt.tight_layout()
    plot_path = os.path.join(args.outdir, f"top10_diff_CGI_{comparison[0]}_{comparison[1]}.png")
    plt.savefig(plot_path)
    plt.close()

# Genes with More than One Affected CpG Island
top_df_all = pd.DataFrame(top_changes, columns=["CpG_Island", "Avg_Delta", "n"])
top_df_all["Gene"] = top_df_all["CpG_Island"].str.extract(r"chr\w+_\d+_\d+_(.+?)_")
gene_df = top_df_all.groupby("Gene").agg(count=("CpG_Island", "count"), avg_delta=("Avg_Delta", "mean")).reset_index()
multi_cpg_genes = gene_df[gene_df["count"] > 1].sort_values("avg_delta")

plt.figure(figsize=(10, 6))
sns.barplot(data=multi_cpg_genes, x="avg_delta", y="Gene", color='#000080')
plt.xlabel("Avg Change in Methylated Fragment Count")
plt.axvline(0, color="gray", linestyle="--")
plt.title("Genes with More than One Affected CpG Island")
plt.tight_layout()
plot_path = os.path.join(args.outdir, "multi_CpG_genes.png")
plt.savefig(plot_path)
plt.close()
