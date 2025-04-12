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
    })

    # === Plotting Section ===
    fig, ax = plt.subplots(figsize=(10, 6))  # Adjust figure size if needed

    # Example Plot (replace with your actual plotting code)
    some_plot = ax.plot(range(len(samples)), np.random.rand(len(samples)), label="Example Data")

    # Adjust layout to reduce blank space
    plt.subplots_adjust(left=0.1, right=0.85, top=0.9, bottom=0.1)  # Adjust these values as needed

    # Add colorbar (if applicable)
    cbar = fig.colorbar(some_plot[0], ax=ax)  # Replace `some_plot[0]` with your actual mappable object
    cbar.ax.set_position([0.88, 0.15, 0.03, 0.7])  # [left, bottom, width, height]

    # Adjust legend position
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)

    # Use tight_layout as a fallback
    plt.tight_layout()

    # Save plot
    plot_path = os.path.join("plots", f"{os.path.splitext(fname)[0]}_adjusted_plot.png")
    fig.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

# === Create ZIP of all plots ===
zip_path = os.path.join("plots", "bubbleplots.zip")
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
    for root, _, files in os.walk("plots"):
        for file in files:
            if file.startswith("bubbleplot_") and (file.endswith(".png") or file.endswith(".svg")):
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, start="plots")
                zipf.write(file_path, arcname=arcname)

print(f"All bubble plot files zipped and saved to: {zip_path}")

# === Remove individual plot files after zipping ===
for root, _, files in os.walk("plots"):
    for file in files:
        if file.startswith("bubbleplot_") and (file.endswith(".png") or file.endswith(".svg")):
            file_path = os.path.join(root, file)
            if file_path != zip_path:
                os.remove(file_path)

print(f"Individual bubble plot files have been removed after zipping.")
