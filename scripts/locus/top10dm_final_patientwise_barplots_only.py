import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import ttest_rel

output_dir = 'plots/top10dm_barplots'
os.makedirs(output_dir, exist_ok=True)

# Helper functions
def classify_timepoint(sample):
    if "Baseline" in sample:
        return "Baseline"
    elif "Off-tx" in sample:
        return "Post-Treatment"
    elif "INNOV" in sample:
        return "Healthy"
    else:
        return "On-Treatment"

def get_patient(sample, patient_ids):
    return next((pid for pid in patient_ids if pid in sample), None)

# Load data
output_folder = 'output'
data_folder = 'data'

cpg_matrix_file = glob.glob(os.path.join(output_folder, "*matrix*"))[0]
patient_list_file = glob.glob(os.path.join(data_folder, "*patient*"))[0]
cgi_map_file = os.path.join(output_folder, "gene_cgi_map.csv")

cpg_matrix = pd.read_excel(cpg_matrix_file, index_col=0) if cpg_matrix_file.endswith('.xlsx') else pd.read_csv(cpg_matrix_file, sep="\t", index_col=0)
patient_df = pd.read_excel(patient_list_file)
gene_map_df = pd.read_csv(cgi_map_file)

gene_map_df.columns = gene_map_df.columns.str.strip()
gene_map_df["cgi_id"] = gene_map_df["cgi_id"].astype(str).str.strip()
cgi_to_gene = dict(zip(gene_map_df["cgi_id"], gene_map_df["gene_name"]))

patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()

# Delta calculation logic
def calculate_deltas(cpg_matrix, timepoint1, timepoint2, patient_ids):
    deltas = []
    for cpg in cpg_matrix.index:
        values = []
        for patient in patient_ids:
            t1 = cpg_matrix.loc[cpg, cpg_matrix.columns.str.contains(f"{timepoint1}_{patient}")]
            t2 = cpg_matrix.loc[cpg, cpg_matrix.columns.str.contains(f"{timepoint2}_{patient}")]
            if not t1.empty and not t2.empty:
                values.append(t2.mean() - t1.mean())
        if len(values) > 1:
            deltas.append((cpg, np.mean(values)))
    return pd.DataFrame(deltas, columns=['CpG_Island', 'Avg_Delta']).sort_values('Avg_Delta', ascending=False)

# Visualization functions
def plot_top10_diff_cgi_subregions(df, title, filename):
    plt.figure(figsize=(10, 6))
    sns.barplot(data=df.head(10), x="Avg_Delta", y="CpG_Island", color='darkblue')
    plt.axvline(0, color="gray", linestyle="--")
    plt.xlabel("Avg Change in Scaled Methylation Fragment Count Ratio")
    plt.title(title, loc='center')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, filename))
    plt.close()

def plot_multi_cpg_genes(df, title, filename):
    df["Gene"] = df["CpG_Island"].map(cgi_to_gene)
    gene_df = df.groupby("Gene").agg(
        count=("CpG_Island", "count"),
        avg_delta=("Avg_Delta", "mean")
    ).reset_index()

    multi_cpg_genes = gene_df[gene_df["count"] > 1].sort_values("avg_delta", ascending=False)

    plt.figure(figsize=(10, 6))
    sns.barplot(data=multi_cpg_genes.head(10), x="avg_delta", y="Gene", color='darkblue')
    plt.axvline(0, color="gray", linestyle="--")
    plt.xlabel("Avg Change in Scaled Methylation Fragment Count Ratio")
    plt.title(title, loc='center')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, filename))
    plt.close()

# Generate plots
comparisons = [
    ("Baseline", "Post-Treatment", "baseline_post"),
    ("Baseline", "On-Treatment", "baseline_on"),
    ("On-Treatment", "Post-Treatment", "on_post"),
]

for t1, t2, suffix in comparisons:
    deltas_df = calculate_deltas(cpg_matrix, t1, t2, patient_ids)
    plot_top10_diff_cgi_subregions(
        deltas_df,
        f"Top 10 Differentially Methylated Subregions of CpG Islands ({t1} vs {t2})",
        f"top10_diff_CGIsubregions_{suffix}.png"
    )
    plot_multi_cpg_genes(
        deltas_df,
        f"Genes with More than One Affected CpG Island ({t1} vs {t2})",
        f"multi_CpG_genes_{suffix}.png"
    )
