import pandas as pd
import matplotlib.pyplot as plt
import os
import glob

# Set up directories
input_dir = "output"
plot_dir = os.path.join("plots", "dotplots")
os.makedirs(plot_dir, exist_ok=True)

# Define condition classifier
def classify_condition(name):
    if "INNOV" in name:
        return "Healthy"
    elif "Baseline" in name:
        return "Baseline"
    elif "Off-tx" in name:
        return "Post-Treatment"
    else:
        return "On-Treatment"

# Define plotting order and color palette
order = ["Healthy", "Baseline", "On-Treatment", "Post-Treatment"]
palette = {
    "Healthy": "#3B0A45",
    "Baseline": "#0A42A3",
    "On-Treatment": "#1A7771",
    "Post-Treatment": "#D7542D",
}

# Find all matching Excel files
excel_files = glob.glob(os.path.join(input_dir, "*scaled_fragment_ratios_matrix*.xlsx"))

if not excel_files:
    print("⚠️ No matching Excel files found in the 'output/' directory.")
else:
    for filepath in excel_files:
        filename = os.path.basename(filepath)
        print(f"Processing: {filename}")

        # Load and clean data
        df = pd.read_excel(filepath)
        df.columns = ['Sample', 'Glob20', 'GlobMin80', 'Scaled_Ratio']
        df['Condition'] = df['Sample'].apply(classify_condition)
        x_labels = [f"{cond}\n(n={len(df[df['Condition'] == cond])})" for cond in order]

        # Plot 1: Median scatter plot
        plt.figure(figsize=(8, 6))
        for i, cond in enumerate(order):
            group = df[df['Condition'] == cond]['Scaled_Ratio']
            plt.scatter([i]*len(group), group, color=palette[cond], s=60, alpha=0.8, edgecolors='k', linewidth=0.5)
            plt.plot([i - 0.2, i + 0.2], [group.median()] * 2, color=palette[cond], lw=3)
        plt.xticks(range(len(order)), x_labels)
        plt.ylabel("CpG Methylation\n(Scaled Ratio x100K)", fontsize=12)
        plt.xlabel("Sample Condition", fontsize=12)
        plt.ylim(bottom=0)
        plt.title(f"{filename} - Median Lines", fontsize=14)
        plt.tight_layout()
        median_path = os.path.join(plot_dir, f"{filename.replace('.xlsx', '')}_median_dotplot.png")
        plt.savefig(median_path, dpi=300)
        plt.close()

        # Plot 2: Mean ± SD scatter plot
        plt.figure(figsize=(8, 6))
        for i, cond in enumerate(order):
            group = df[df['Condition'] == cond]['Scaled_Ratio']
            mean = group.mean()
            sd = group.std()
            plt.scatter([i]*len(group), group, color=palette[cond], s=60, alpha=0.8, edgecolors='k', linewidth=0.5)
            plt.plot([i - 0.2, i + 0.2], [mean] * 2, color=palette[cond], lw=3)
            plt.errorbar(i, mean, yerr=sd, fmt='none', ecolor='gray', capsize=5, lw=1.5)
        plt.xticks(range(len(order)), x_labels)
        plt.ylabel("CpG Methylation\n(Scaled Ratio x100K)", fontsize=12)
        plt.xlabel("Sample Condition", fontsize=12)
        plt.ylim(bottom=0)
        plt.title(f"{filename} - Mean ± SD", fontsize=14)
        plt.tight_layout()
        mean_sd_path = os.path.join(plot_dir, f"{filename.replace('.xlsx', '')}_mean_sd_dotplot.png")
        plt.savefig(mean_sd_path, dpi=300)
        plt.close()

        print(f"✔️ Saved plots to:\n  - {median_path}\n  - {mean_sd_path}")

print("✅ Done.")
