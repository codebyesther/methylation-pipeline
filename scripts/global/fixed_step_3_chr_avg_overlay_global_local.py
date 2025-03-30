import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description="Global-level methylation visualization.")
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
matrix = matrix.apply(pd.to_numeric, errors='coerce').fillna(0)
collapsed = matrix.T.groupby(level=[0, 1]).mean().T

chromosome_lookup = pd.Series(collapsed.index.to_series().str.extract(r"chr([^_]+)")[0], index=collapsed.index)

chr_order = [str(i) for i in range(1, 23)] + ["X", "Y"]
all_comparisons = [
    ("Baseline", "On-Treatment", "Baseline → On-Tx"),
    ("On-Treatment", "Post-Treatment", "On-Tx → Post-Tx"),
    ("Baseline", "Post-Treatment", "Baseline → Post-Tx")
]

summaries = []
for t0, t1, label in all_comparisons:
    deltas = []
    for patient in collapsed.columns.levels[0]:
        if (patient, t0) in collapsed.columns and (patient, t1) in collapsed.columns:
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

bar_data = merged[merged["Comparison"].isin(["Baseline → On-Tx", "On-Tx → Post-Tx"])]
line_data = merged[merged["Comparison"] == "Baseline → Post-Tx"]

colors = {
    "Baseline → On-Tx": "#1f77b4",
    "On-Tx → Post-Tx": "#8B0000",
    "Baseline → Post-Tx": "#000000"
}

x = np.arange(len(chr_order))
width = 0.35

bar_pivot = bar_data.pivot(index="Chromosome", columns="Comparison", values="Mean_Delta").reindex(chr_order).fillna(0)
line_vals = line_data.set_index("Chromosome").reindex(chr_order).fillna(0)["Mean_Delta"].values

fig, ax = plt.subplots(figsize=(16, 6))
ax.bar(x - width/2, bar_pivot["Baseline → On-Tx"], width, label="Baseline → On-Tx", color=colors["Baseline → On-Tx"], alpha=0.8)
ax.bar(x + width/2, bar_pivot["On-Tx → Post-Tx"], width, label="On-Tx → Post-Tx", color=colors["On-Tx → Post-Tx"], alpha=0.6)
ax.plot(x, line_vals, marker="o", color=colors["Baseline → Post-Tx"], label="Baseline → Post-Tx", linewidth=2)

ax.set_xticks(x)
ax.set_xticklabels(chr_order)
ax.set_xlabel("Chromosome")
ax.set_ylabel("Mean Methylation Delta")
ax.set_title("Avg Methylation Change per Chromosome")
ax.axhline(0, color="gray", linestyle="--")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)
plt.tight_layout()

overlay_path = os.path.join(args.outdir, "chr_avg_overlay_aligned.png")
plt.savefig(overlay_path, bbox_inches='tight')
plt.close()

excel_path = os.path.join(args.outdir, "chr_avg_summary.xlsx")
merged.to_excel(excel_path, index=False, sheet_name="Chromosome_Deltas")
