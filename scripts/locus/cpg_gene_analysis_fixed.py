
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
import argparse

# Argument parsing
parser = argparse.ArgumentParser(description="CpG subregion and gene-level methylation analysis")
parser.add_argument("--datadir", type=str, required=True, help="Directory with methylation data and patient files")
parser.add_argument("--outdir", type=str, default="plots", help="Output directory for plots")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# Automated file detection
methylation_file = next(f for f in os.listdir(args.datadir) if "patient" not in f.lower() and f.lower().endswith('.xlsx'))
patient_file = next(f for f in os.listdir(args.datadir) if "patient" in f.lower())

print("Detected files:", methylation_file, patient_file)

methylation_df = pd.read_excel(os.path.join(args.datadir, methylation_file))
patient_df = pd.read_excel(os.path.join(args.datadir, patient_file))
patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()

# Extract locus-level matrix
start_idx = methylation_df[methylation_df.iloc[:, 0].astype(str).str.contains("CGI_chr", na=False)].index[0]
cpg_island_df = methylation_df.iloc[start_idx:].reset_index(drop=True)
cpg_island_df.rename(columns={cpg_island_df.columns[0]: "CpG_Island"}, inplace=True)
cpg_island_df.dropna(how="all", subset=cpg_island_df.columns[1:], inplace=True)

# Metadata extraction
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
matrix = matrix.apply(pd.to_numeric, errors='coerce')
collapsed = matrix.groupby(level=[0, 1], axis=1).mean()

# Genes with multiple CpG islands
multi_cpg_genes_list = []
comparisons = [("Baseline", "On-Treatment"), 
               ("On-Treatment", "Post-Treatment"), 
               ("Baseline", "Post-Treatment")]

for comp in comparisons:
    comp_label = f"{comp[0]} → {comp[1]}"
    gene_changes = []

    for cpg in collapsed.index:
        deltas = []
        for patient in collapsed.columns.levels[0]:
            if (patient, comp[0]) in collapsed.columns and (patient, comp[1]) in collapsed.columns:
                before = collapsed.loc[cpg, (patient, comp[0])]
                after = collapsed.loc[cpg, (patient, comp[1])]
                if pd.notna(before) and pd.notna(after):
                    deltas.append(after - before)
        
        if len(deltas) >= 1:
            avg_delta = np.mean(deltas)
            gene_name = re.search(r"CGI_chr[\w]+_[\d]+_[\d]+_(.+?)_", cpg)
            gene = gene_name.group(1) if gene_name else "Unknown"
            gene_changes.append((gene, avg_delta))

    gene_df = pd.DataFrame(gene_changes, columns=["Gene", "Avg_Delta"])
    gene_counts = gene_df.groupby("Gene").size()
    multi_gene_df = gene_df[gene_df["Gene"].isin(gene_counts[gene_counts > 1].index)]
    multi_gene_agg = multi_gene_df.groupby("Gene")["Avg_Delta"].mean().reset_index().sort_values("Avg_Delta")

    multi_cpg_genes_list.append(multi_gene_agg.assign(Comparison=comp_label))

    plt.figure(figsize=(10, 6))
    sns.barplot(data=multi_gene_agg, x="Avg_Delta", y="Gene", color="#000080")
    plt.axvline(0, color="gray", linestyle="--")
    plt.title(f"Genes with Multiple Affected CpG Islands ({comp_label})")
    plt.xlabel("Avg Change in Methylated Fragment Count")
    plt.ylabel("Gene")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, f"multi_CpG_genes_{comp[0]}_{comp[1]}.png"))
    plt.close()

# Aggregate all comparisons
aggregate_genes = pd.concat(multi_cpg_genes_list).groupby("Gene")["Avg_Delta"].mean().reset_index().sort_values("Avg_Delta")

plt.figure(figsize=(10, 6))
sns.barplot(data=aggregate_genes, x="Avg_Delta", y="Gene", color="#000080")
plt.axvline(0, color="gray", linestyle="--")
plt.title("Genes with Multiple Affected CpG Islands (All Comparisons)")
plt.xlabel("Avg Change in Methylated Fragment Count")
plt.ylabel("Gene")
plt.tight_layout()
plt.savefig(os.path.join(args.outdir, "multi_CpG_genes_all_comparisons.png"))
plt.close()

print("Plots generated successfully.")
