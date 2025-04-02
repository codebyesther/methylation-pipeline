# -*- coding: utf-8 -*-
"""Local-scale figures_Part2_BubblePlots_03292025_InterRowSpacingAdjustment.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1-fKb2KGMaxSqcah5i0NgBFbH85ynNxiG
"""

import io, os, zipfile, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from matplotlib.gridspec import GridSpec
from tqdm import tqdm  # Import tqdm for progress bars

sns.set(style="whitegrid")
os.makedirs("plots", exist_ok=True)

# === Load Files ===
methylation_dfs = {}
patient_ids = []

# Automatically detect files
output_dir = "output/"
data_dir = "data/"

# Find files with "ratios_matrix" in their name from the output directory
ratio_files = glob.glob(os.path.join(output_dir, "*ratios_matrix*.xlsx")) + glob.glob(os.path.join(output_dir, "*ratios_matrix*.csv"))

# Find files with "patient" in their name from the data directory
patient_files = glob.glob(os.path.join(data_dir, "*patient*.xlsx")) + glob.glob(os.path.join(data_dir, "*patient*.csv"))

# Process patient files
for file_path in tqdm(patient_files, desc="Processing patient files"):
    df = pd.read_excel(file_path, header=None) if file_path.endswith('.xlsx') else pd.read_csv(file_path, header=None)
    values = df[0].dropna().astype(str).tolist()
    patient_ids.extend([v for v in values if not v.lower().startsWith("unnamed")])

# Process ratio files
for file_path in tqdm(ratio_files, desc="Processing ratio files"):
    file_name = os.path.basename(file_path)
    df = pd.read_excel(file_path) if file_path.endswith('.xlsx') else pd.read_csv(file_path)
    methylation_dfs[file_name] = df

# === Normalize sample names ===
def normalize_timepoint(sample):
    if "Baseline" in sample:
        return "Baseline"
    elif "Off-tx" in sample:
        return "Post-Treatment"
    elif "INNOV" in sample:
        return "Healthy"
    else:
        return "On-Treatment"

# Define timepoints
timepoints_patient = ["Baseline", "On-Treatment", "Post-Treatment"]  # Per-patient plots
timepoints_chromosome = ["Baseline", "On-Treatment", "Post-Treatment", "Healthy"]  # Per-chromosome plots

# Define timepoint positions for plotting
timepoint_positions_patient = {
    "Baseline": 1.0,
    "On-Treatment": 0.9,
    "Post-Treatment": 1.1
}

timepoint_positions_chromosome = {
    "Baseline": 1.0,
    "On-Treatment": 0.9,
    "Post-Treatment": 1.1,
    "Healthy": 1.2
}

# === Process Files ===
for fname, df in tqdm(methylation_dfs.items(), desc="Processing methylation files"):
    print(f"\n=== Processing file: {fname} ===")
    start_idx = df[df.iloc[:, 0].astype(str).str.contains("CGI_chr", na=False)].index[0]
    cpg_df = df.iloc[start_idx:].reset_index(drop=True)
    cpg_df.rename(columns={cpg_df.columns[0]: "CpG_Island"}, inplace=True)
    cpg_df.dropna(how="all", subset=cpg_df.columns[1:], inplace=True)

    samples = cpg_df.columns[1:]
    def get_patient(s): return next((pid for pid in patient_ids if pid in s), None)
    sample_meta = pd.DataFrame({
        "Sample": samples,
        "Patient": [get_patient(s) for s in samples],
        "Timepoint": [normalize_timepoint(s) for s in samples]
    }).dropna()

    # Remove Healthy samples from matrix
    non_healthy = sample_meta[sample_meta["Timepoint"] != "Healthy"]
    matrix = cpg_df.set_index("CpG_Island")[non_healthy["Sample"]]
    matrix.columns = pd.MultiIndex.from_frame(non_healthy[["Patient", "Timepoint"]])
    matrix = matrix.apply(pd.to_numeric, errors='coerce')
    collapsed = matrix.groupby(axis=1, level=[0, 1]).mean()

    # === Extract CpG Coordinates ===
    cpg_coords = cpg_df["CpG_Island"].str.extract(r"CGI_(chr\w+)_(\d+)_(\d+)", expand=True)
    cpg_coords.columns = ["Chr", "Start", "End"]
    cpg_coords = cpg_coords.dropna()
    cpg_coords["Start"] = cpg_coords["Start"].astype(int)
    cpg_coords["End"] = cpg_coords["End"].astype(int)
    cpg_coords["Midpoint"] = (cpg_coords["Start"] + cpg_coords["End"]) // 2
    cpg_coords["CpG_Island"] = cpg_df["CpG_Island"]
    coords_df = cpg_coords.set_index("CpG_Island")

    # === Merge data with coordinates ===
    collapsed_flat = collapsed.copy()
    collapsed_flat.columns = ['{}_{}'.format(p, tp) for p, tp in collapsed.columns]
    collapsed_flat = collapsed_flat.reset_index()

    bubble_data = pd.merge(collapsed_flat, coords_df.reset_index(), on="CpG_Island").set_index("CpG_Island")

# === Bubble plots per patient per chromosome ===
for patient in tqdm(collapsed.columns.levels[0], desc="Generating bubble plots per patient"):
    for chrom in coords_df["Chr"].unique():
        # Build a subset for each timepoint, merging them
        chr_data = bubble_data[bubble_data["Chr"] == chrom]

        # Collect rows for all timepoints in a single DataFrame
        all_rows = []
        for tp in timepoints_patient:
            col = f"{patient}_{tp}"
            if col not in chr_data.columns:
                print(f"[WARNING] Column {col} not found for {patient}, {chrom}, skipping.")
                continue
            sub = chr_data[[col, "Midpoint"]].dropna().rename(columns={col: "value"})
            if sub.empty:
                continue
            sub["Timepoint"] = tp
            all_rows.append(sub)

        if not all_rows:
            # No data for any timepoint
            continue

        subset_df = pd.concat(all_rows, ignore_index=True)

        # Create figure with 2 columns:
        # - Left col = main bubble plot
        # - Right col = sub-gridspec for colorbar (top) + bubble legend (bottom)
        fig = plt.figure(figsize=(14, 6))
        gs = GridSpec(nrows=1, ncols=2, width_ratios=[5, 1], figure=fig)

        # Main axis on the left
        ax_main = fig.add_subplot(gs[0, 0])

        # Sub-gridspec on the right: 2 rows (colorbar top, legend bottom)
        gs_right = gs[0, 1].subgridspec(nrows=2, ncols=1, height_ratios=[0.5, 0.5])
        ax_cbar = fig.add_subplot(gs_right[0, 0])
        ax_legend = fig.add_subplot(gs_right[1, 0])

        # Plot each row in subset_df
        sc = None
        for _, row in subset_df.iterrows():
            tp = row["Timepoint"]
            y = timepoint_positions_patient[tp]
            # Increase bubble scale factor for bigger bubbles
            bubble_size = row["value"]**0.5 * 150
            sc = ax_main.scatter(
                row["Midpoint"],
                y,
                s=bubble_size,
                c=row["value"],
                cmap="viridis",
                alpha=0.6,
                vmin=0,  # lower bound of color scale
                vmax=170    # upper bound of color scale
            )

        # Format main axis
        ax_main.set_yticks(list(timepoint_positions_patient.values()))  # Ensure the number of ticks matches the number of labels
        ax_main.set_yticklabels(timepoints_patient)
        ax_main.set_ylim(0.75, 1.25) # set y-limits so large bubbles have padding above & below
        ax_main.set_xlabel("CpG Island Genomic Coordinate Midpoint (bp)")
        ax_main.set_ylabel("Timepoint")
        ax_main.set_title(f"DNA Hypermethylation Profiles Throughout Treatment\nPatient: {patient}, Chromosome: {chrom}")

        # If we got at least one scatter point, make the colorbar in ax_cbar
        if sc is not None:
            fig.colorbar(sc, cax=ax_cbar, label="Fragment Count")

        # Create bubble-size legend in ax_legend
        ax_legend.axis("off")  # hide ticks and background
        handles, labels = [], []
        for s in [1, 15, 150]:
            h = ax_legend.scatter([], [], s=s**0.5 * 150, color="gray", alpha=0.5)
            handles.append(h)
            labels.append(str(s))

        ax_legend.legend(
            handles,
            labels,
            title="Bubble Size\n(Fragment Count)",
            labelspacing=1.5,
            handletextpad=2.0,
            borderpad=1.3,
            loc="center",
            bbox_to_anchor=(0.5, 0.5)
        )

        plt.tight_layout()
        filename_base = os.path.join("plots", f"bubble_{patient}_{chrom}")
        plt.savefig(f"{filename_base}.png")
        plt.savefig(f"{filename_base}.svg")
        plt.close()

# === Bubble plots per chromosome (averaged across patients) ===
for chrom in tqdm(coords_df["Chr"].unique(), desc="Generating bubble plots per chromosome"):
    chr_data = bubble_data[bubble_data["Chr"] == chrom]

    all_rows = []
    for tp in timepoints_chromosome:
        # Find columns that end with e.g. "_Baseline"
        cols = [c for c in chr_data.columns if c.endsWith(f"_{tp}")]
        if not cols:
            continue
        avg = chr_data[cols].mean(axis=1)
        sub = chr_data.loc[avg.notna(), ["Midpoint"]].copy()
        sub["value"] = avg[avg.notna()]
        sub["Timepoint"] = tp
        all_rows.append(sub)

    if not all_rows:
        continue

    subset_df = pd.concat(all_rows, ignore_index=True)

    fig = plt.figure(figsize=(14, 6))
    gs = GridSpec(nrows=1, ncols=2, width_ratios=[5, 1], figure=fig)

    ax_main = fig.add_subplot(gs[0, 0])
    gs_right = gs[0, 1].subgridspec(nrows=2, ncols=1, height_ratios=[0.5, 0.5])
    ax_cbar = fig.add_subplot(gs_right[0, 0])
    ax_legend = fig.add_subplot(gs_right[1, 0])

    sc = None
    for _, row in subset_df.iterrows():
        tp = row["Timepoint"]
        y = timepoint_positions_chromosome[tp]
        bubble_size = row["value"]**0.5 * 50
        sc = ax_main.scatter(
            row["Midpoint"],
            y,
            s=bubble_size,
            c=row["value"],
            cmap="viridis",
            alpha=0.6,
            vmin=0,  # lower bound of color scale
            vmax=170    # upper bound of color scale
        )

    ax_main.set_yticks(list(timepoint_positions_chromosome.values()))  # Ensure the number of ticks matches the number of labels
    ax_main.set_yticklabels(timepoints_chromosome)
    ax_main.set_ylim(0.75, 1.25) # set y-limits so large bubbles have padding above & below
    ax_main.set_xlabel("CpG Island Genomic Coordinate Midpoint (bp)")
    ax_main.set_ylabel("Timepoint")
    ax_main.set_title("DNA Hypermethylation Profiles Throughout Treatment (Averaged Across Patients)")

    if sc is not None:
        fig.colorbar(sc, cax=ax_cbar, label="Fragment Count")

    ax_legend.axis("off")
    handles, labels = [], []
    for s in [1, 15, 150]:
        h = ax_legend.scatter([], [], s=s**0.5 * 50, color="gray", alpha=0.5)
        handles.append(h)
        labels.append(str(s))

    ax_legend.legend(
        handles,
        labels,
        title="Bubble Size\n(Fragment Count)",
        labelspacing=1.5,
        handletextpad=2.0,
        borderpad=1.3,
        loc="center",
        bbox_to_anchor=(0.5, 0.5)
    )

    plt.tight_layout()
    filename_base = os.path.join("plots", f"bubble_{chrom}")
    plt.savefig(f"{filename_base}.png")
    plt.savefig(f"{filename_base}.svg")
    plt.close()

# === Zip All Plots ===
zipf = zipfile.ZipFile("bubble_plots.zip", "w", zipfile.ZIP_DEFLATED)
for root, dirs, files in os.walk("plots"):
    for f in files:
        zipf.write(os.path.join(root, f), arcname=f)
zipf.close()

# === Create ZIP of all plots ===
zip_path = os.path.join("plots", "bubbleplots.zip")
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
    for root, _, files in os.walk("plots"):
        for file in files:
            file_path = os.path.join(root, file)
            arcname = os.path.relpath(file_path, start="plots")
            zipf.write(file_path, arcname=arcname)

print(f"All plots zipped and saved to: {zip_path}")
