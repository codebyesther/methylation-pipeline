import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === Directory Setup ===
input_dir = "output"
plot_dir = os.path.join("plots", "lineplots")
per_patient_dir = os.path.join(plot_dir, "per_patient")
replicate_table_dir = os.path.join("output", "per_patient_tables")
os.makedirs(per_patient_dir, exist_ok=True)
os.makedirs(replicate_table_dir, exist_ok=True)

metadata_path = os.path.join("output", "sample_metadata.csv")
summary_stats_path = os.path.join("output", "summary_statistics.csv")
per_patient_summary_path = os.path.join("output", "per_patient_summary.csv")

time_order = ["Baseline", "On-Treatment", "Post-Treatment"]

# === Functions ===
def assign_patient_id(name):
    if "INNOV" in name:
        return None
    elif name.startswith("LOI") and name.count('_') >= 2:
        base = name.split('_')[1]
    elif "_" in name:
        base = name.split('_')[0]
    else:
        return None
    return base.split("-MVS")[0] if "-MVS" in base else base

def simplify_timepoint(name):
    if "INNOV" in name:
        return None
    elif "Baseline" in name:
        return "Baseline"
    elif "Off-tx" in name:
        return "Post-Treatment"
    else:
        return "On-Treatment"

def extract_replicate(name):
    parts = name.split('_')
    if len(parts) > 2:
        return "_".join(parts[2:])
    elif len(parts) == 2:
        return parts[1] if parts[1] not in ["Baseline", "Off-tx"] else None
    return None

def make_main_plot_with_box(df, save_path):
    fig, ax = plt.subplots(figsize=(10, 6))
    palette_cb = sns.color_palette("colorblind", len(df["Patient_ID"].unique()))
    for i, (pid, group) in enumerate(df.groupby("Patient_ID")):
        group = group.sort_values("Timepoint")
        ax.plot(
            group["Timepoint"],
            group["Scaled_Ratio"],
            marker='o',
            linestyle='-',
            linewidth=2,
            markersize=8,
            color=palette_cb[i % len(palette_cb)],
            label=pid
        )
    ax.set_ylabel("Scaled Methylation Fragment Count Ratio", fontsize=12)
    ax.set_xlabel("Treatment Timepoint", fontsize=12)
    ax.set_title("CpG Methylation Trajectories by Patient", fontsize=14)
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
        ax.set_ylabel("Scaled Methylation Fragment Count Ratio")
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
    ax.set_ylabel("Scaled Methylation Fragment Count Ratio")
    ax.set_xlabel("Treatment Timepoint")
    ax.set_ylim(0, max(agg["mean"].max(), df["Scaled_Ratio"].max()) * 1.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(save_path, dpi=300)
    plt.close(fig)

def make_standalone_boxplot(df, save_path):
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.boxplot(data=df, x="Timepoint", y="Scaled_Ratio", palette="pastel", ax=ax)
    ax.set_title("Distribution of Methylation by Timepoint (Boxplot)")
    ax.set_ylabel("Scaled Methylation Fragment Count Ratio")
    ax.set_xlabel("Treatment Timepoint")
    fig.tight_layout()
    fig.savefig(save_path, dpi=300)
    plt.close(fig)

def make_violin_plot(df, save_path):
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.violinplot(data=df, x="Timepoint", y="Scaled_Ratio", palette="pastel", ax=ax, cut=0, inner=None)
    sns.swarmplot(data=df, x="Timepoint", y="Scaled_Ratio", color="k", size=4, ax=ax)
    ax.set_title("Distribution of Methylation by Timepoint (Violin + Swarm)")
    ax.set_ylabel("Scaled Methylation Fragment Count Ratio")
    ax.set_xlabel("Treatment Timepoint")
    fig.tight_layout()
    fig.savefig(save_path, dpi=300)
    plt.close(fig)

def export_boxplot_summary(df, save_path):
    summary = df.groupby("Timepoint")["Scaled_Ratio"].agg(["count", "mean", "median", "std"]).reindex(time_order)
    summary.to_csv(save_path)

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

# === Build Full Raw Table
df_full = pd.DataFrame({
    "Sample": sample_names,
    "Scaled_Ratio": scaled_ratios
})
df_full["Patient_ID"] = df_full["Sample"].apply(assign_patient_id)
df_full["Timepoint"] = df_full["Sample"].apply(simplify_timepoint)
df_full["Replicate_ID"] = df_full["Sample"].apply(extract_replicate)
df_full = df_full.dropna(subset=["Patient_ID", "Timepoint"])
df_full["Timepoint"] = pd.Categorical(df_full["Timepoint"], categories=time_order, ordered=True)

# === Save sample metadata
df_full[["Sample", "Patient_ID", "Timepoint", "Replicate_ID"]].to_csv(metadata_path, index=False)
print(f"✅ Saved sample metadata to: {metadata_path}")

# === Per-patient timepoint summary
summary = df_full.groupby(["Patient_ID", "Timepoint"]).size().unstack(fill_value=0)
summary.to_csv(per_patient_summary_path)
print(f"✅ Saved per-patient timepoint summary to: {per_patient_summary_path}")

# === Replicate tables per patient
for pid, group in df_full.groupby("Patient_ID"):
    group[["Sample", "Timepoint", "Replicate_ID", "Scaled_Ratio"]].to_csv(
        os.path.join(replicate_table_dir, f"{pid}.csv"), index=False
    )
print(f"✅ Saved per-patient replicate tables to: {replicate_table_dir}")

# === Average across replicates for plotting
df_plot = df_full.groupby(["Patient_ID", "Timepoint"], as_index=False)["Scaled_Ratio"].mean()
valid = df_plot["Patient_ID"].value_counts()
df_plot = df_plot[df_plot["Patient_ID"].isin(valid[valid > 1].index)]

# === Summary stats (used for longitudinal trajectory)
summary_stats = df_plot.groupby("Timepoint")["Scaled_Ratio"].agg(["count", "mean", "median", "std"]).reindex(time_order)
summary_stats.to_csv(summary_stats_path)
print(f"✅ Saved summary stats to: {summary_stats_path}")

# === Plot Generation
make_main_plot_with_box(df_plot, os.path.join(plot_dir, "methylation_longitudinal_plot.png"))
make_per_patient_plots(df_plot, per_patient_dir)
make_average_trajectory_plot(df_plot, os.path.join(plot_dir, "average_trajectory.png"))
make_standalone_boxplot(df_full, os.path.join(plot_dir, "boxplot_by_timepoint.png"))
make_violin_plot(df_full, os.path.join(plot_dir, "violinplot_by_timepoint.png"))
export_boxplot_summary(df_full, os.path.join("output", "boxplot_summary_by_timepoint.csv"))

print("✅ All plots and tables saved successfully.")
