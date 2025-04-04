import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === Directories ===
input_dir = "output"
base_plot_dir = os.path.join("plots", "lineplots")
metadata_path = os.path.join("output", "sample_metadata.csv")

# === Helper Functions ===
def assign_patient_id(name):
    if "INNOV" in name:
        return name
    return name.split('_')[0]

def simplify_timepoint(name):
    if "INNOV" in name:
        return "Healthy"
    elif name.endswith("Baseline"):
        return "Baseline"
    elif name.endswith("Off-tx"):
        return "Post-Treatment"
    return "On-Treatment"

def make_plot(df, innov_mean, save_path, overlay_innov):
    unique_patients = df["Patient_ID"].unique()
    palette_cb = sns.color_palette("colorblind", len(unique_patients))
    markers = ['o', 's', 'D', '^', 'v', 'P', '*', 'X', 'H', '8', '<', '>']
    linestyles = ['-', '--', '-.', ':']

    fig, ax = plt.subplots(figsize=(10, 6))
    for i, (pid, group) in enumerate(df.groupby("Patient_ID")):
        group = group.sort_values("Timepoint")
        ax.plot(
            group["Timepoint"],
            group["Scaled_Ratio"],
            marker=markers[i % len(markers)],
            linestyle=linestyles[i % len(linestyles)],
            linewidth=2,
            markersize=8,
            color=palette_cb[i % len(palette_cb)],
            label=pid
        )
    if overlay_innov and innov_mean is not None:
        ax.axhline(innov_mean, color="gray", linestyle="--", linewidth=2, label="INNOV Avg")

    ax.set_ylabel("Scaled Ratio (x100K)", fontsize=12)
    ax.set_xlabel("Treatment Timepoint", fontsize=12)
    ax.set_title("Longitudinal CpG Methylation (Baseline → On-Tx → Post-Tx)", fontsize=14)
    ax.tick_params(axis='x', labelrotation=0, labelsize=11)
    fig.tight_layout()
    ax.legend(title="Patient ID", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11, title_fontsize=12)
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

def make_per_patient_plots(df, innov_mean, save_dir, overlay_innov):
    os.makedirs(save_dir, exist_ok=True)
    for pid, group in df.groupby("Patient_ID"):
        fig, ax = plt.subplots(figsize=(7, 4))
        group = group.sort_values("Timepoint")
        ax.plot(
            group["Timepoint"],
            group["Scaled_Ratio"],
            marker='o',
            linestyle='-',
            linewidth=2,
            markersize=8,
            color='navy',
            label=pid
        )
        if overlay_innov and innov_mean is not None:
            ax.axhline(innov_mean, color="gray", linestyle="--", linewidth=2, label="INNOV Avg")
        ax.set_title(f"Patient: {pid}")
        ax.set_ylabel("Scaled Ratio (x100K)")
        ax.set_xlabel("Treatment Timepoint")
        ax.set_ylim(0, max(df["Scaled_Ratio"].max(), innov_mean or 0) * 1.2)
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(save_dir, f"{pid}.png"), dpi=300)
        plt.close(fig)

# === Step 1: Locate File ===
file_to_use = None
for fname in os.listdir(input_dir):
    if "scaled_fragment_ratios_matrix" in fname.lower() and fname.endswith((".xlsx", ".xls")):
        file_to_use = os.path.join(input_dir, fname)
        break
if file_to_use is None:
    raise FileNotFoundError("No file with 'scaled_fragment_ratios_matrix' found in the output/ directory.")

# === Step 2: Load first two rows ===
df_raw = pd.read_excel(file_to_use, header=None, nrows=2)
sample_names = df_raw.iloc[0, 1:].tolist()
scaled_ratios = df_raw.iloc[1, 1:].tolist()

# === Step 3: Create tidy DataFrame ===
df_long = pd.DataFrame({
    "Sample": sample_names,
    "Scaled_Ratio": scaled_ratios
})
df_long["Patient_ID"] = df_long["Sample"].apply(assign_patient_id)
df_long["Timepoint"] = df_long["Sample"].apply(simplify_timepoint)

# === Step 4: Save sample metadata ===
df_long[["Sample", "Patient_ID", "Timepoint"]].to_csv(metadata_path, index=False)
print(f"✅ Saved sample metadata to: {metadata_path}")

# === Step 5: Filter non-INNOV with ≥2 timepoints ===
df_non_innov = df_long[df_long["Timepoint"] != "Healthy"]
valid = df_non_innov["Patient_ID"].value_counts()
df_non_innov = df_non_innov[df_non_innov["Patient_ID"].isin(valid[valid > 1].index)]
df_non_innov["Timepoint"] = pd.Categorical(df_non_innov["Timepoint"], categories=["Baseline", "On-Treatment", "Post-Treatment"], ordered=True)

# === Step 6: Get INNOV average ===
df_innov = df_long[df_long["Timepoint"] == "Healthy"]
innov_mean = df_innov["Scaled_Ratio"].mean() if not df_innov.empty else None

# === Step 7: Generate plots with & without INNOV overlay ===
for overlay_innov in [True, False]:
    label = "with_innov" if overlay_innov else "without_innov"
    plot_dir = os.path.join(base_plot_dir, label)
    os.makedirs(plot_dir, exist_ok=True)
    per_patient_dir = os.path.join(plot_dir, "per_patient")
    os.makedirs(per_patient_dir, exist_ok=True)

    # Main plot
    main_path = os.path.join(plot_dir, "methylation_longitudinal_plot.png")
    make_plot(df_non_innov, innov_mean, main_path, overlay_innov)

    # Per-patient plots
    make_per_patient_plots(df_non_innov, innov_mean, per_patient_dir, overlay_innov)

    print(f"✅ Saved {'WITH' if overlay_innov else 'WITHOUT'} INNOV overlay plots to: {plot_dir}")
