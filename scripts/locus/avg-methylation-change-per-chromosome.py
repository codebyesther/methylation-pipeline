import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import sys

import matplotlib
matplotlib.use("Agg")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common import (classify_timepoint, find_patient_file, find_scaled_matrix,  # noqa: E402
                    get_patient as _get_patient, read_patient_ids, read_scaled_matrix)

# Auto-detect files
patient_file = find_patient_file()
methylation_path = find_scaled_matrix()
methylation_file = os.path.basename(methylation_path)
print(f"Using patient file: {os.path.basename(patient_file)}")
print(f"Using methylation file: {methylation_file}")

patient_ids = read_patient_ids(patient_file)
base_fname = os.path.splitext(methylation_file)[0]

# CGI x sample ratios; "INF" cells (GlobMin80 = 0) become NaN and are left out
cgi, _ = read_scaled_matrix(methylation_path)
cpg_island_df = cgi.reset_index().rename(columns={cgi.index.name or "index": "CpG_Island"})
cpg_island_df.rename(columns={cpg_island_df.columns[0]: "CpG_Island"}, inplace=True)

def get_patient(sample):
    return _get_patient(sample, patient_ids)

def normalize_timepoint(sample):
    return classify_timepoint(sample)

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
matrix = matrix.apply(pd.to_numeric, errors='coerce').replace([np.inf, -np.inf], np.nan)  # INF cells skipped by mean()
collapsed = matrix.T.groupby(level=[0, 1]).mean().T

chromosome_lookup = pd.Series(
    collapsed.index.to_series().str.extract(r"chr([^_]+)")[0],
    index=collapsed.index
)

chr_order = [str(i) for i in range(1, 23)] + ["X", "Y"]
all_comparisons = [
    ("Baseline", "On-Treatment", "Baseline → On-Tx"),
    ("On-Treatment", "Post-Treatment", "On-Tx → Post-Tx"),
    ("Baseline", "Post-Treatment", "Baseline → Post-Tx"),
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

os.makedirs("plots/avg-methylation-change-per-chromosome", exist_ok=True)

fig, ax = plt.subplots(figsize=(16, 6))
ax.bar(x - width / 2, bar_pivot["Baseline → On-Tx"], width, label="Baseline → On-Tx", color=colors["Baseline → On-Tx"], alpha=0.8)
ax.bar(x + width / 2, bar_pivot["On-Tx → Post-Tx"], width, label="On-Tx → Post-Tx", color=colors["On-Tx → Post-Tx"], alpha=0.8)
ax.plot(x, line_vals, marker="o", color=colors["Baseline → Post-Tx"], label="Baseline → Post-Tx", linewidth=2)

ax.set_xticks(x)
ax.set_xticklabels(chr_order)
ax.set_xlabel("Chromosome")
ax.set_ylabel("Mean Scaled Methylation Fragment Count Ratio Delta")
ax.set_title("Avg Methylation Change per Chromosome")
ax.axhline(0, color="gray", linestyle="--")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)
plt.tight_layout()

# Save results
plot_path = f"plots/avg-methylation-change-per-chromosome/chr_avg_overlay_{base_fname}_aligned.png"
excel_path = f"plots/avg-methylation-change-per-chromosome/chr_avg_summary_{base_fname}.xlsx"
plt.savefig(plot_path, bbox_inches="tight")
plt.close(fig)

merged.to_excel(excel_path, index=False)

print(f"Saved plot to {plot_path}")
print(f"Saved Excel summary to {excel_path}")
