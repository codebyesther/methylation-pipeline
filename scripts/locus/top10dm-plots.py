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
top10dmplot_filenames = []

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
        plt.figure(figsize=(10, 6))
        sns.barplot(data=df.head(10), x="Avg_Delta", y="CpG_Island", color='darkblue')
        plt.xlabel("Avg Change in Scaled Methylated Fragment Count Ratio")
        plt.axvline(0, color="gray", linestyle="--")
        plt.title(title)
        plt.tight_layout(pad=10.0)  # Add padding to avoid clipping
        plot_path = os.path.join(args.outdir, filename)
        plt.savefig(plot_path)
        plt.close()
        top10dmplot_filenames.append(plot_path)

    def plot_multi_cpg_genes(df, title, filename):
        df["Gene"] = df["CpG_Island"].str.extract(r"chr\w+_\d+_\d+_(.+?)_")
        gene_df = df.groupby("Gene").agg(count=("CpG_Island", "count"), avg_delta=("Avg_Delta", "mean")).reset_index()
        multi_cpg_genes = gene_df[gene_df["count"] > 1].sort_values("avg_delta")
        
        plt.figure(figsize=(10, 6))
        sns.barplot(data=multi_cpg_genes, x="avg_delta", y="Gene", color='darkblue')
        plt.xlabel("Avg Change in Scaled Methylated Fragment Count Ratio")
        plt.axvline(0, color="gray", linestyle="--")
        plt.title(title)
        plt.tight_layout()
        plot_path = os.path.join(args.outdir, filename)
        plt.savefig(plot_path)
        plt.close()
        top10dmplot_filenames.append(plot_path)

    # Generate and plot for baseline vs post-treatment
    top_df_baseline_post = calculate_deltas(collapsed, "Baseline", "Post-Treatment")
    plot_top10_diff_cgi_subregions(top_df_baseline_post, "Top 10 Differentially Methylated Subregions of CpG Islands (Baseline vs Post-Treatment)", "top10_diff_CGIsubregions_baseline_post.png")
    plot_multi_cpg_genes(top_df_baseline_post, "Genes with More than One Affected CpG Island (Baseline vs Post-Treatment)", "multi_CpG_genes_baseline_post.png")

    # Generate and plot for baseline vs on-treatment
    top_df_baseline_on = calculate_deltas(collapsed, "Baseline", "On-Treatment")
    plot_top10_diff_cgi_subregions(top_df_baseline_on, "Top 10 Differentially Methylated Subregions of CpG Islands (Baseline vs On-Treatment)", "top10_diff_CGIsubregions_baseline_on.png")
    plot_multi_cpg_genes(top_df_baseline_on, "Genes with More than One Affected CpG Island (Baseline vs On-Treatment)", "multi_CpG_genes_baseline_on.png")

    # Generate and plot for on-treatment vs post-treatment
    top_df_on_post = calculate_deltas(collapsed, "On-Treatment", "Post-Treatment")
    plot_top10_diff_cgi_subregions(top_df_on_post, "Top 10 Differentially Methylated Subregions of CpG Islands (On-Treatment vs Post-Treatment)", "top10_diff_CGIsubregions_on_post.png")
    plot_multi_cpg_genes(top_df_on_post, "Genes with More than One Affected CpG Island (On-Treatment vs Post-Treatment)", "multi_CpG_genes_on_post.png")

# Create a zip file containing all plot files
zip_filename = os.path.join(args.outdir, 'top-10-differential-methylation-plots.zip')
with zipfile.ZipFile(zip_filename, 'w') as zipf:
    for plot_filename in top10dmplot_filenames:
        zipf.write(plot_filename, arcname=os.path.basename(plot_filename))

# Delete the individual plot files after zipping
for plot_filename in top10dmplot_filenames:
    os.remove(plot_filename)

print(f'Saved plots and zipped them in {zip_filename}')
