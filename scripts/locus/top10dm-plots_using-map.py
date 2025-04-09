import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import zipfile

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

# Load CpG → Gene mapping
map_file = os.path.join(output_folder, "gene_cgi_map.csv")
if not os.path.exists(map_file):
    raise FileNotFoundError("The gene_cgi_map.csv file is missing from the 'output/' folder.")
gene_map_df = pd.read_csv(map_file)
gene_map_df.columns = gene_map_df.columns.str.strip()
cgi_to_gene = dict(zip(gene_map_df["cgi_id"], gene_map_df["gene_name"]))

methylation_dfs = {}
patient_ids = []
os.makedirs(args.outdir, exist_ok=True)

# Read patient IDs
if "patient" in args.patients.lower():
    patient_df = pd.read_excel(args.patients)
    patient_ids = patient_df.iloc[:, 0].dropna().astype(str).tolist()

# Read methylation matrix
df = pd.read_excel(args.methylation)
methylation_dfs[os.path.basename(args.methylation)] = df

# Metadata helpers
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

# Output plot tracking
top10dmplot_filenames = []

# Process methylation files
for fname, df in methylation_dfs.items():
    print(f"\n=== Processing file: {fname} ===")
    base_fname = os.path.splitext(fname)[0]

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

    def calculate_deltas(collapsed, timepoint1, timepoint2):
        changes = []
        for cpg in collapsed.index:
            deltas = []
            for patient in collapsed.columns.levels[0]:
                if (patient, timepoint1) in collapsed.columns and (patient, timepoint2) in collapsed.columns:
                    t1 = collapsed.loc[cpg, (patient, timepoint1)]
                    t2 = collapsed.loc[cpg, (patient, timepoint2)]
                    if not pd.isna(t1) and not pd.isna(t2):
                        deltas.append(t2 - t1)
            if len(deltas) >= 2:
                changes.append((cpg, np.mean(deltas), len(deltas)))
        return pd.DataFrame(changes, columns=["CpG_Island", "Avg_Delta", "n"]).sort_values("Avg_Delta")

    def plot_top10_diff_cgi_subregions(df, title, filename):
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.barplot(data=df.head(10), x="Avg_Delta", y="CpG_Island", color='darkblue', ax=ax)
        ax.set_xlabel("Avg Change in Scaled Methylated Fragment Count Ratio")
        ax.axvline(0, color="gray", linestyle="--")
        ax.set_title("")
        fig.text(0.17, 0.98, title, ha='left', va='top', fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        plot_path = os.path.join(args.outdir, filename)
        fig.savefig(plot_path)
        plt.close(fig)
        top10dmplot_filenames.append(plot_path)

    def plot_multi_cpg_genes(df, title, filename):
        df["Gene"] = df["CpG_Island"].map(cgi_to_gene)
        df = df.dropna(subset=["Gene"])

        gene_df = df.groupby("Gene").agg(
            count=("CpG_Island", "count"),
            avg_delta=("Avg_Delta", "mean")
        ).reset_index()
        multi_cpg_genes = gene_df[gene_df["count"] > 1].sort_values("avg_delta")

        if multi_cpg_genes.empty:
            print(f"No multi-CpG genes found for: {title}")
            return

        plt.figure(figsize=(10, 6))
        sns.barplot(data=multi_cpg_genes, x="avg_delta", y="Gene", color='darkblue')
        plt.xlabel("Avg Change in Scaled Methylated Fragment Count Ratio")
        plt.axvline(0, color="gray", linestyle="--")
        plt.title(title, loc='center')
        plt.tight_layout()
        plot_path = os.path.join(args.outdir, filename)
        plt.savefig(plot_path)
        plt.close()
        top10dmplot_filenames.append(plot_path)

    # Comparisons
    comparisons = [
        ("Baseline", "Post-Treatment", "baseline_post"),
        ("Baseline", "On-Treatment", "baseline_on"),
        ("On-Treatment", "Post-Treatment", "on_post"),
    ]

    for t1, t2, suffix in comparisons:
        top_df = calculate_deltas(collapsed, t1, t2)
        plot_top10_diff_cgi_subregions(
            top_df,
            f"Top 10 Differentially Methylated Subregions of CpG Islands ({t1} vs {t2})",
            f"top10_diff_CGIsubregions_{suffix}.png"
        )
        plot_multi_cpg_genes(
            top_df,
            f"Genes with More than One Affected CpG Island ({t1} vs {t2})",
            f"multi_CpG_genes_{suffix}.png"
        )

# Zip plots
zip_filename = os.path.join(args.outdir, 'top-10-differential-methylation-plots.zip')
with zipfile.ZipFile(zip_filename, 'w') as zipf:
    for plot_filename in top10dmplot_filenames:
        zipf.write(plot_filename, arcname=os.path.basename(plot_filename))

# Clean up individual plots
for plot_filename in top10dmplot_filenames:
    os.remove(plot_filename)

print(f'Saved plots and zipped them in {zip_filename}')
