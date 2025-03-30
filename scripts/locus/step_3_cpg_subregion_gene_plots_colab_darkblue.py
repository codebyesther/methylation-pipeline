import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--datadir", type=str, required=True, help="Directory containing Excel files")
    parser.add_argument("--outdir", type=str, required=False, default="plots")
    return parser.parse_args()

args = parse_args()
os.makedirs(args.outdir, exist_ok=True)

files = os.listdir(args.datadir)
patient_file = [f for f in files if "patient" in f.lower()][0]
methylation_file = [f for f in files if "patient" not in f.lower()][0]

patient_df = pd.read_excel(os.path.join(args.datadir, patient_file))
patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()

df = pd.read_excel(os.path.join(args.datadir, methylation_file))
start_idx = df[df.iloc[:, 0].astype(str).str.contains("CGI_chr", na=False)].index[0]
cpg_island_df = df.iloc[start_idx:].reset_index(drop=True)
cpg_island_df.rename(columns={cpg_island_df.columns[0]: "CpG_Island"}, inplace=True)

# Metadata functions
def get_patient(sample):
    return next((pid for pid in patient_ids if pid in sample), None)

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
matrix = matrix.apply(pd.to_numeric, errors='coerce').fillna(0)
collapsed = matrix.groupby(axis=1, level=[0, 1]).mean()

comparisons = [("Baseline", "On-Treatment"), ("On-Treatment", "Post-Treatment"), ("Baseline", "Post-Treatment")]
all_gene_dfs = []

# Loop over comparisons
for comp in comparisons:
    gene_changes = []
    for cpg in collapsed.index:
        deltas = []
        for patient in collapsed.columns.levels[0]:
            if (patient, comp[0]) in collapsed.columns and (patient, comp[1]) in collapsed.columns:
                val1, val2 = collapsed.loc[cpg, (patient, comp[0])], collapsed.loc[cpg, (patient, comp[1])]
                deltas.append(val2 - val1)
        if len(deltas) >= 2:
            avg_delta = np.mean(deltas)
            gene = cpg.split("_")[3]  # extract gene name
            gene_changes.append((gene, avg_delta))

    gene_df = pd.DataFrame(gene_changes, columns=["Gene", "Avg_Delta"])
    multi_gene_df = gene_df.groupby("Gene").agg(count=("Avg_Delta", "count"), avg_delta=("Avg_Delta", "mean"))
    multi_gene_df = multi_gene_df[multi_gene_df["count"] > 1].sort_values("avg_delta")

    if not multi_gene_df.empty:
        plt.figure(figsize=(10, 6))
        sns.barplot(x="avg_delta", y=multi_gene_df.index, data=multi_gene_df, color="#000080")
        plt.axvline(0, color="gray", linestyle="--")
        plt.xlabel("Avg Change in Methylated Fragment Count")
        plt.ylabel("Gene")
        plt.title(f"Genes with Multiple Affected CpG Islands ({comp[0]} → {comp[1]})")
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, f"multi_cpg_genes_{comp[0]}_{comp[1]}.png"))
        plt.close()

    multi_gene_df['Comparison'] = f"{comp[0]} → {comp[1]}"
    all_gene_dfs.append(multi_gene_df)

# Aggregate data across all comparisons
agg_gene_df = pd.concat(all_gene_dfs).groupby('Gene').agg(count=("count", "sum"), avg_delta=("avg_delta", "mean"))
agg_gene_df = agg_gene_df.sort_values("avg_delta")

if not agg_gene_df.empty:
    plt.figure(figsize=(10, 6))
    sns.barplot(x="avg_delta", y=agg_gene_df.index, data=agg_gene_df, color="#000080")
    plt.axvline(0, color="gray", linestyle="--")
    plt.xlabel("Avg Change in Methylated Fragment Count")
    plt.ylabel("Gene")
    plt.title("Aggregated Genes with Multiple Affected CpG Islands (All Comparisons)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "multi_cpg_genes_all_comparisons.png"))
    plt.close()
