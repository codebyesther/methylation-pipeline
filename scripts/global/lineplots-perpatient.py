import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set directories
input_dir = "output"
plot_dir = os.path.join("plots", "lineplots")
os.makedirs(plot_dir, exist_ok=True)

# Detect the correct Excel file
file_to_use = None
for fname in os.listdir(input_dir):
    if "scaled_fragment_ratios_matrix" in fname.lower() and fname.endswith((".xlsx", ".xls")):
        file_to_use = os.path.join(input_dir, fname)
        break

if file_to_use is None:
    raise FileNotFoundError("No file with 'scaled_fragment_ratios_matrix' found in the output/ directory.")

# Load data
df = pd.read_excel(file_to_use)

# Clean column names
df.columns = ['Sample', 'Glob20', 'GlobMin80', 'Scaled_Ratio']

# Assign patient ID
def assign_patient_id(name):
    return name if "INNOV" in name else name.split('_')[0]

# Assign simplified timepoint
def simplify_timepoint(name):
    if "INNOV" in name:
        return None
    elif "Baseline" in name:
        return "Baseline"
    elif "Off-tx" in name:
        return "Post-Treatment"
    else:
        return "On-Treatment"

# Apply transformations
df['Patient_ID'] = df['Sample'].apply(assign_patient_id)
df['Timepoint'] = df['Sample'].apply(simplify_timepoint)

# Drop healthy and invalid entries
df = df.dropna(subset=['Timepoint'])

# Average duplicate On-Tx entries
df_avg = df.groupby(['Patient_ID', 'Timepoint'], as_index=False)['Scaled_Ratio'].mean()

# Keep only patients with 2+ timepoints
valid = df_avg['Patient_ID'].value_counts()
df_avg = df_avg[df_avg['Patient_ID'].isin(valid[valid > 1].index)]

# Sort timepoints
time_order = ['Baseline', 'On-Treatment', 'Post-Treatment']
df_avg['Timepoint'] = pd.Categorical(df_avg['Timepoint'], categories=time_order, ordered=True)

# Plot
unique_patients = df_avg['Patient_ID'].unique()
palette_cb = sns.color_palette("colorblind", len(unique_patients))
markers = ['o', 's', 'D', '^', 'v', 'P', '*', 'X', 'H', '8', '<', '>']
linestyles = ['-', '--', '-.', ':']

fig, ax = plt.subplots(figsize=(10, 6))
for i, (pid, group) in enumerate(df_avg.groupby('Patient_ID')):
    group = group.sort_values('Timepoint')
    ax.plot(
        group['Timepoint'],
        group['Scaled_Ratio'],
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
ax.legend(
    title="Patient ID",
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    fontsize=11,
    title_fontsize=12
)

# Save figure
plot_path = os.path.join(plot_dir, "methylation_longitudinal_plot.png")
fig.savefig(plot_path, dpi=300, bbox_inches='tight')
print(f"Saved plot to {plot_path}")
