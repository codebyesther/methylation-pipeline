import pandas as pd
import matplotlib.pyplot as plt
import os
import glob

# === Setup ===
input_dir = "output"
plot_dir = os.path.join("plots", "dotplots")
os.makedirs(plot_dir, exist_ok=True)

# === Sample classification function ===
def classify_condition(name):
    if "INNOV" in name:
        return "Healthy"
    elif "Baseline" in name:
        return "Baseline"
    elif "Off-tx" in name:
        return "Post-Treatment"
    else:
        return "On-Treatment"

# === Plot settings ===
order = ["Healthy", "Baseline", "On-Treatment", "Post-Treatment"]
palette = {
    "Healthy": "#3B0A45",
    "Baseline": "#0A42A3",
    "On-Treatment": "#1A7771",
    "Post-Treatment": "#D7542D",
}

# === Find all valid Excel files ===
excel_files = [
    f for f in glob.glob(os.path.join(input_dir, "*fragment_ratios_matrix*.xlsx"))
    if not os.path.basename(f).startswith("~$")
]

if not excel_files:
    print("⚠️ No matching Excel files found in the 'output/' directory.")
else:
    for filepath in excel_files:
        filename = os.path.basename(filepath)
        base_name = filename.replace('.xlsx', '')
        print(f"\n🔍 Processing: {filename}")

        # === Load and reformat the scaled matrix ===
        raw_df = pd.read_excel(filepath, header=None)

        samples = raw_df.iloc[0, 1:]  # skip "Header" column
        ratios = raw_df.iloc[1, 1:]

        df = pd.DataFrame({
            "Sample": samples.values,
            "Scaled_Ratio": ratios.values
        })

        # === Annotate and prepare for plotting ===
        df["Condition"] = df["Sample"].apply(classify_condition)
        x_labels = [f"{cond}\n(n={len(df[df['Condition'] == cond])})" for cond in order]

        # === Export Summary Stats ===
        stats = []
        for cond in order:
            group = df[df["Condition"] == cond]["Scaled_Ratio"]
            stats.append({
                "Condition": cond,
                "n": len(group),
                "Mean": round(group.mean(), 2),
                "Median": round(group.median(), 2),
                "SD": round(group.std(), 2)
            })
        stats_df = pd.DataFrame(stats)
        stats_path = os.path.join(plot_dir, f"{base_name}_summary_stats.csv")
        stats_df.to_csv(stats_path, index=False)

        # === Plot 1: Median Scatter Plot ===
        plt.figure(figsize=(8, 6))
        for i, cond in enumerate(order):
            group = df[df['Condition'] == cond]['Scaled_Ratio']
            plt.scatter([i]*len(group), group, color=palette[cond], s=60,
                        alpha=0.8, edgecolors='k', linewidth=0.5)
            plt.plot([i - 0.2, i + 0.2], [group.median()] * 2,
                     color=palette[cond], lw=3)
        plt.xticks(range(len(order)), x_labels)
        plt.ylabel("CpG Methylation\n(Scaled Ratio)", fontsize=12)
        plt.xlabel("Sample Condition", fontsize=12)
        plt.ylim(bottom=0)
        plt.title(f"​CpG Methylation by Condition (Median)", fontsize=14)
        plt.tight_layout()
        median_path = os.path.join(plot_dir, f"{base_name}_median_dotplot.png")
        plt.savefig(median_path, dpi=300)
        plt.close()

        # === Plot 2: Mean ± SD Scatter Plot ===
        plt.figure(figsize=(8, 6))
        for i, cond in enumerate(order):
            group = df[df['Condition'] == cond]['Scaled_Ratio']
            mean = group.mean()
            sd = group.std()
            plt.scatter([i]*len(group), group, color=palette[cond], s=60,
                        alpha=0.8, edgecolors='k', linewidth=0.5)
            plt.plot([i - 0.2, i + 0.2], [mean] * 2,
                     color=palette[cond], lw=3)
            plt.errorbar(i, mean, yerr=sd, fmt='none', ecolor='gray', capsize=5, lw=1.5)
        plt.xticks(range(len(order)), x_labels)
        plt.ylabel("CpG Methylation\n(Scaled Ratio)", fontsize=12)
        plt.xlabel("Sample Condition", fontsize=12)
        plt.ylim(bottom=0)
        plt.title(f"​CpG Methylation by Condition (Mean ± SD)", fontsize=14)
        plt.tight_layout()
        mean_sd_path = os.path.join(plot_dir, f"{base_name}_mean_sd_dotplot.png")
        plt.savefig(mean_sd_path, dpi=300)
        plt.close()

        print(f"✅ Saved outputs:")
        print(f"  - {median_path}")
        print(f"  - {mean_sd_path}")
        print(f"  - {stats_path}")

print("\n🎉 All done!")
