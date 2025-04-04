import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === Directory Setup ===
input_dir = "output"
plot_dir = os.path.join("plots", "lineplots")
per_patient_dir = os.path.join(plot_dir, "per_patient")
os.makedirs(per_patient_dir, exist_ok=True)

metadata_path = os.path.join("output", "sample_metadata.csv")
summary_stats_path = os.path.join("output", "summary_statistics.csv")

time_order = ["Baseline", "On-Treatment", "Post-Treatment"]

# === Functions ===
def assign_patient_id(name):
    if "INNOV" in name:
        return None
    elif name.startswith("LOI") and name.count('_') >= 2:
        return name.split('_')[1]
    elif "_" in name:
        return name.split('_')[0]
    else:
        return None

def simplify_timepoint(name):
    if "INNOV" in name:
        return None
    elif "Baseline" in name:
        return "Baseline"
    elif "Off-tx" in name:
        return "Post-Treatment"
    else:
        return "On-Treatment"

def make_main_plot(df, save_path):
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
    ax.set_ylabel("Scaled Ratio (x100K)", fontsize=12)
    ax.set_xlabel("Treatment Timepoint", fontsize=12)
    ax.set_title("Longitudinal CpG Methylation (Baseline → On-Tx → Post-Tx)", fontsize=14)
    ax.tick_params(axis='x', labelrotation=0, labelsize=11)
    fig.tight_layout()
    ax.legend(title="Patient ID", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11, title_fontsize=12)
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

def make_per_patient_plots(df, save_dir):
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
        ax.set_title(f"Patient: {pid}")
        ax.set_ylabel("Scaled Ratio (x100K)")
        ax.set_xlabel("Treatment Timepoint")
        ax.set_ylim(0, df["Scaled_Ratio"].max() * 1.2)
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(save_dir, f"{pid}.png"), dpi=300)
        plt.close(fig)

def make_average_trajectory_plot(df, save_path):
    agg = df.groupby("Timepoint")["Scaled_Ratio"].agg(["mean", "std"]).reindex(time_order)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(
        agg.index,
        agg["mean"],
        yerr=agg["std"],
        marker='o',
        linestyle='-',
        color='darkgreen',
        linewidth=2,
        capsize=5,
        label="Average"
    )
    ax.set_title("Average Methylation Trajectory Across Patients")
    ax.set_ylabel("Scaled Ratio (x100K)")
    ax.set_xlabel("Treatment Timepoint")
    ax.set_ylim(0, max(agg["mean"].max(), df["Scaled_Ratio"].max()) * 1.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_path, dpi=300)
    plt.close(fig)

# === Load Excel File ===
file_to_use = next(
    (os.path.join(input_dir, f) for f in os.listdir(input_dir)
     if "scaled_fragment_ratios_matrix" in f.lower() and f.endswith((".xlsx", ".xls"))),
    None
)
if file_to_use is None:
    raise FileNotFoundError("No file with 'scaled_fragment_ratios_matrix' found in the output/ directory.")

df_raw = pd.read_excel(file_to_use, header=None, nrows=2)
sample_names = df_raw.iloc[0, 1:].tolist()
scaled_ratios = df_raw.iloc[1, 1:].tolist()

# === Build and Clean Data ===
df_long = pd.DataFrame({
    "Sample": sample_names,
    "Scaled_Ratio": scaled_ratios
})
df_long["Patient_ID"] = df_long["Sample"].apply(assign_patient_id)
df_long["Timepoint"] = df_long["Sample"].apply(simplify_timepoint)
df_long = df_long.dropna(subset=["Patient_ID", "Timepoint"])
df_long["Timepoint"] = pd.Categorical(df_long["Timepoint"], categories=time_order, ordered=True)

# === Save Sample Metadata ===
df_long[["Sample", "Patient_ID", "Timepoint"]].to_csv(metadata_path, index=False)
print(f"✅ Saved sample metadata to: {metadata_path}")

# === Average over duplicate Patient_ID + Timepoint ===
df_long = df_long.groupby(["Patient_ID", "Timepoint"], as_index=False)["Scaled_Ratio"].mean()

# === Filter for patients with ≥2 timepoints ===
valid = df_long["Patient_ID"].value_counts()
df_long = df_long[df_long["Patient_ID"].isin(valid[valid > 1].index)]

# === Save Summary Stats ===
summary_stats = df_long.groupby("Timepoint")["Scaled_Ratio"].agg(["count", "mean", "median", "std"]).reindex(time_order)
summary_stats.to_csv(summary_stats_path)
print(f"✅ Saved summary stats to: {summary_stats_path}")

# === Generate Plots ===
make_main_plot(df_long, os.path.join(plot_dir, "methylation_longitudinal_plot.png"))
make_per_patient_plots(df_long, per_patient_dir)
make_average_trajectory_plot(df_long, os.path.join(plot_dir, "average_trajectory.png"))

print(f"✅ All plots saved in: {plot_dir}")
