# -*- coding: utf-8 -*-

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
    patient_ids.extend([v for v in values if not v.lower().startswith("unnamed")])
    
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
timepoints_chromosome = ["Baseline", "On-Treatment", "Post-Treatment"]  # Per-chromosome plots (removed "Healthy")

# Define timepoint positions for plotting with increased spacing
timepoint_positions_patient = {
    "Baseline": 0.6,
    "On-Treatment": 1.0,
    "Post-Treatment": 1.4
}

timepoint_positions_chromosome = {
    "Baseline": 0.6,
    "On-Treatment": 1.0,
    "Post-Treatment": 1.4
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
        fig = plt.figure(figsize=(18, 10))  # Increased figure size to prevent cropping
        gs = GridSpec(nrows=1, ncols=2, width_ratios=[5, 1], figure=fig)

        # Main axis on the left
        ax_main = fig.add_subplot(gs[0, 0])
        ax_main.grid(False)  # Hide grid

        # Sub-gridspec on the right: 2 rows (colorbar top, legend bottom)
        gs_right = gs[0, 1].subgridspec(nrows=2, ncols=1, height_ratios=[0.5, 0.5])
        ax_cbar = fig.add_subplot(gs_right[0, 0])
        ax_cbar.grid(False)  # Hide grid
        ax_legend = fig.add_subplot(gs_right[1, 0])

        # Plot each row in subset_df
        sc = None
        for _, row in subset_df.iterrows():
            tp = row["Timepoint"]
            y = timepoint_positions_patient[tp]
            # Increase bubble scale factor for bigger bubbles
            bubble_size = row["value"]**0.5 * 50
            sc = ax_main.scatter(
                row["Midpoint"],
                y,
                s=bubble_size,
                c=row["value"],
                cmap="viridis",
                alpha=0.6,
                vmin=0,  # lower bound of color scale
                vmax=2000    # increased upper bound of color scale
            )

        # Format main axis
        ax_main.set_yticks(list(timepoint_positions_patient.values()))  # Ensure the number of ticks matches the number of labels
        ax_main.set_yticklabels(timepoints_patient, fontsize=14)
        ax_main.set_ylim(0.3, 1.7)  # set y-limits so large bubbles have padding above & below
        # Add x-axis padding to avoid bubble clipping
        x_min, x_max = subset_df['Midpoint'].min(), subset_df['Midpoint'].max()
        x_range = x_max - x_min
        ax_main.set_xlim(x_min - 0.1 * x_range, x_max + 0.1 * x_range)  # Increase x-axis padding

        ax_main.set_xlabel("CpG Island Genomic Coordinate Midpoint (bp)", fontsize=16)
        ax_main.set_ylabel("Timepoint", fontsize=16)
        ax_main.set_title(f"DNA Hypermethylation Profiles Throughout Treatment\nPatient: {patient}, Chromosome: {chrom}", fontsize=18)

        # If we got at least one scatter point, make the colorbar in ax_cbar
        if sc is not None:
            cb = fig.colorbar(sc, cax=ax_cbar, label="Scaled Fragment Count Ratio")
            cb.set_label("Scaled Fragment Count Ratio", fontsize=14)

        # Define the x-coordinate for the legend title and scatter plot positions
        legend_x_coord = 0.25

        # Create bubble-size legend in ax_legend
        ax_legend.axis("off")  # hide ticks and background

        # Define the sizes and calculate bubble sizes
        sizes = [1, 80, 800, 8000]
        bubble_sizes = [size**0.5 * 50 for size in sizes]

        # Calculate proportional vertical positions based on bubble radii
        cumulative_height = np.cumsum([size**0.5 for size in sizes])
        total_height = cumulative_height[-1]
        positions = np.array([0.1, 1, 2, 3.4]) * 9/ 1000 * total_height / len(sizes)    # vertical spacing between gray bubble markers

        # Set the x-axis limits explicitly for the legend axis
        ax_legend.set_xlim(0, 1)

        # Print the axis limits to check if legend_x_coord is within range
        x_min, x_max = ax_legend.get_xlim()
        print(f"Legend x-axis limits: min={x_min}, max={x_max}")
        print(f"legend_x_coord: {legend_x_coord}")

        if legend_x_coord < x_min or legend_x_coord > x_max:
            print(f"Warning: legend_x_coord ({legend_x_coord}) is out of bounds!")

        # Manually draw the legend using scatter and text
        for size, pos in zip(sizes, positions):
            ax_legend.scatter(legend_x_coord, pos, s=size**0.5 * 50, color="gray", alpha=0.5)    # Use legend_x_coord for gray bubble x-coordinate adjustment
            ax_legend.text(legend_x_coord + 0.4, pos, str(size), verticalalignment='center', horizontalalignment='center', fontsize=14)    # Adjust text position based on legend_x_coord

        # Legend title slightly above top bubble, adjust the position as needed
        ax_legend.text(legend_x_coord + 0.2, positions[-1] + 0.3, "Bubble Size\n(Scaled Fragment Count Ratio)",    # title is 0.3 above the top bubble
                    horizontalalignment='center', verticalalignment='center', fontweight='bold', fontsize=14)

        # Adjust y-limits to ensure no clipping
        ax_legend.set_ylim(0, positions[-1] + 0.3)    # legend y-axis limit should be larger than the distance between top bubble and title

        fig.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.3, hspace=0.5)
        filename_base = os.path.join("plots", f"bubbleplot_{patient}_{chrom}")
        plt.savefig(f"{filename_base}.png")
        plt.savefig(f"{filename_base}.svg")
        plt.close()

# === Bubble plots per chromosome (averaged across patients) ===
for chrom in tqdm(coords_df["Chr"].unique(), desc="Generating bubble plots per chromosome"):
    chr_data = bubble_data[bubble_data["Chr"] == chrom]

    all_rows = []
    for tp in timepoints_chromosome:
        # Find columns that end with e.g. "_Baseline"
        cols = [c for c in chr_data.columns if c.endswith(f"_{tp}")]
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

    fig = plt.figure(figsize=(18, 10))  # Increased figure size to prevent cropping
    gs = GridSpec(nrows=1, ncols=2, width_ratios=[5, 1], figure=fig)

    ax_main = fig.add_subplot(gs[0, 0])
    ax_main.grid(False)  # Hide grid
    gs_right = gs[0, 1].subgridspec(nrows=2, ncols=1, height_ratios=[0.5, 0.5])
    ax_cbar = fig.add_subplot(gs_right[0, 0])
    ax_cbar.grid(False)  # Hide grid
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
            vmax=2000    # increased upper bound of color scale
        )

    ax_main.set_yticks(list(timepoint_positions_chromosome.values()))  # Ensure the number of ticks matches the number of labels
    ax_main.set_yticklabels(timepoints_chromosome, fontsize=14)
    ax_main.set_ylim(0.3, 1.7)  # set y-limits so large bubbles have padding above & below
    # Add x-axis padding to avoid bubble clipping
    x_min, x_max = subset_df['Midpoint'].min(), subset_df['Midpoint'].max()
    x_range = x_max - x_min
    ax_main.set_xlim(x_min - 0.1 * x_range, x_max + 0.1 * x_range)  # Increase x-axis padding

    ax_main.set_xlabel("CpG Island Genomic Coordinate Midpoint (bp)", fontsize=16)
    ax_main.set_ylabel("Timepoint", fontsize=16)
    ax_main.set_title(f"DNA Hypermethylation Profiles Throughout Treatment (Averaged Across Patients)\nChromosome: {chrom}", fontsize=18)

    if sc is not None:
        cbar = fig.colorbar(sc, cax=ax_cbar, label="Scaled Fragment Count Ratio")
        cbar.set_label("Scaled Fragment Count Ratio", fontsize=14)

    ax_legend.axis("off")

    # Use the same logic as the per-patient plots for legend spacing
    sizes = [1, 80, 800, 8000]
    bubble_sizes = [size**0.5 * 50 for size in sizes]

    # Calculate proportional vertical positions based on bubble radii
    cumulative_height = np.cumsum([size**0.5 for size in sizes])
    total_height = cumulative_height[-1]
    positions = np.array([0.1, 1, 2, 3.4]) * 9/ 1000 * total_height / len(sizes)    # vertical spacing between gray bubble markers

    # Set the x-axis limits explicitly for the legend axis
    ax_legend.set_xlim(0, 1)

    # Print the axis limits to check if legend_x_coord is within range
    x_min, x_max = ax_legend.get_xlim()
    print(f"Legend x-axis limits: min={x_min}, max={x_max}")
    print(f"legend_x_coord: {legend_x_coord}")

    if legend_x_coord < x_min or legend_x_coord > x_max:
        print(f"Warning: legend_x_coord ({legend_x_coord}) is out of bounds!")

    # Manually draw the legend using scatter and text
    for size, pos in zip(sizes, positions):
        ax_legend.scatter(legend_x_coord, pos, s=size**0.5 * 50, color="gray", alpha=0.5)    # Use legend_x_coord for gray bubble x-coordinate adjustment
        ax_legend.text(legend_x_coord + 0.4, pos, str(size), verticalalignment='center', horizontalalignment='center', fontsize=14)    # Adjust text position based on legend_x_coord

    # Legend title slightly above top bubble, adjust the position as needed
    ax_legend.text(legend_x_coord + 0.2, positions[-1] + 0.3, "Bubble Size\n(Scaled Fragment Count Ratio)",    # title is 0.3 above the top bubble
                horizontalalignment='center', verticalalignment='center', fontweight='bold', fontsize=14)

    # Adjust y-limits to ensure no clipping
    ax_legend.set_ylim(0, positions[-1] + 0.3)    # legend y-axis limit should be larger than the distance between top bubble and title

    fig.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.3, hspace=0.5)
    filename_base = os.path.join("plots", f"bubbleplot_{chrom}")
    plt.savefig(f"{filename_base}.png")
    plt.savefig(f"{filename_base}.svg")
    plt.close()

# === Create ZIP of all plots ===
zip_path = os.path.join("plots", "bubbleplots.zip")
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
    for root, _, files in os.walk("plots"):
        for file in files:
            if file.startswith("bubbleplot_") and (file.endswith(".png") or file.endswith(".svg")):
                file
