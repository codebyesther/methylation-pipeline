import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Global-level methylation visualization with auto file detection.")
parser.add_argument("--datadir", default="data", help="Directory containing input Excel files")
parser.add_argument("--outdir", default="plots", help="Output directory for plots")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# Auto-detect patient and methylation files
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

# Load methylation data
df = pd.read_excel(methylation_file)

# Locate methylation data
start_idx = df[df.iloc[:, 0].astype(str).str.contains("CGI_chr", na=False)].index[0]
cpg_island_df = df.iloc[start_idx:].reset_index(drop=True)
cpg_island_df.rename(columns={cpg_island_df.columns[0]: "CpG_Island"}, inplace=True)
cpg_island_df.dropna(how="all", subset=cpg_island_df.columns[1:], inplace=True)

# Define helper functions
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
valid_samples = sample_meta.dropna().query('Timepoint != "Healthy"')

matrix = cpg_island_df.set_index("CpG_Island")[valid_samples["Sample"].tolist()]
matrix.columns = pd.MultiIndex.from_frame(valid_samples[["Patient", "Timepoint"]])
matrix = matrix.apply(pd.to_numeric, errors='coerce').fillna(0)
collapsed = matrix.T.groupby(level=[0, 1]).mean().T

chromosome_lookup = pd.Series(collapsed.index.to_series().str.extract(r"chr([^_]+)")[0], index=collapsed.index)
chr_order = [str(i) for i in range(1, 23)] + ["X", "Y"]

comparisons = [
    ("Baseline", "On-Treatment", "Baseline → On-Tx"),
    ("On-Treatment", "Post-Treatment", "On-Tx → Post-Tx"),
    ("Baseline", "Post-Treatment", "Baseline → Post-Tx")
]

summaries = []
for t0, t1, label in comparisons:
    deltas = []
    for patient in collapsed.columns.levels[0]:
        if (patient, t0) in collapsed and (patient, t1) in collapsed:
            delta = collapsed[(patient, t1)] - collapsed[(patient, t0)]
            for cpg, val in delta.items():
                chrom = chromosome_lookup.at[cpg]
                if pd.notna(val) and pd.notna(chrom):
                    deltas.append((chrom, val))
    df_delta = pd.DataFrame(deltas, columns=["Chromosome", "Delta"])
    summary = df_delta.groupby("Chromosome").agg(Mean_Delta=("Delta", "mean"))
    summary = summary.reindex(chr_order).fillna(0).reset_index()
    summary["Comparison"] = label
    summaries.append(summary)

merged = pd.concat(summaries)
merged["Chromosome"] = pd.Categorical(merged["Chromosome"], categories=chr_order, ordered=True)

bar_data = merged.query('Comparison != "Baseline → Post-Tx"')
line_data = merged.query('Comparison == "Baseline → Post-Tx"')

x = np.arange(len(chr_order))
width = 0.35

bar_pivot = bar_data.pivot(index="Chromosome", columns="Comparison", values="Mean_Delta").fillna(0)
line_vals = line_data.set_index("Chromosome").reindex(chr_order).fillna(0)["Mean_Delta"]

fig, ax = plt.subplots(figsize=(16, 6))
ax.bar(x - width/2, bar_pivot["Baseline → On-Tx"], width, label="Baseline → On-Tx", color="#1f77b4", alpha=0.8)
ax.bar(x + width/2, bar_pivot["On-Tx → Post-Tx"], width, label="On-Tx → Post-Tx", color="#8B0000", alpha=0.6)
ax.plot(x, line_vals, marker="o", color="#000000", label="Baseline → Post-Tx", linewidth=2)

ax.set_xticks(x)
ax.set_xticklabels(chr_order)
ax.set_xlabel("Chromosome")
ax.set_ylabel("Mean Methylation Delta")
ax.set_title("Avg Methylation Change per Chromosome")
ax.axhline(0, color="gray", linestyle="--")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)
plt.tight_layout()

plot_path = os.path.join(args.outdir, "chr_avg_overlay.png")
plt.savefig(plot_path, bbox_inches='tight')
plt.close()

excel_path = os.path.join(args.outdir, "chr_avg_summary.xlsx")
merged.to_excel(excel_path, index=False, sheet_name="Chromosome_Deltas")
