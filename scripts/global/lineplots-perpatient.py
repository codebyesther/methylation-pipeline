import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set directories
input_dir = "output"
plot_dir = os.path.join("plots", "lineplots")
os.makedirs(plot_dir, exist_ok=True)

# Step 1: Find the file
file_to_use = None
for fname in os.listdir(input_dir):
    if "scaled_fragment_ratios_matrix" in fname.lower() and fname.endswith((".xlsx", ".xls")):
        file_to_use = os.path.join(input_dir, fname)
        break

if file_to_use is None:
    raise FileNotFoundError("No file with 'scaled_fragment_ratios_matrix' found in the output/ directory.")

# Step 2: Load first two rows only
df_raw = pd.read_excel(file_to_use, header=None, nrows=2)

# Extract sample names (row 0, skip first column if index)
sample_names = df_raw.iloc[0, 1:].tolist()

# Extract scaled methylation values (row 1, skip first column if index)
scaled_ratios = df_raw.iloc[1, 1:].tolist()

# Build DataFrame
df_long = pd.DataFrame({
    "Sample": sample_names,
    "Scaled_Ratio": scaled_ratios
})

# Step 3: Extract patient and timepoint
def assign_patient_id(name):
    if "INNOV" in name:
        return name
    else:
        return "_".join(name.split('_')[:2])  # ← fixed logic here

def simplify_timepoint(name):
    if "INNOV" in name:
        return None
    elif "Baseline" in name:
        return "Baseline"
    elif "Off-tx" in name:
        return "Post-Treatment"
    else:
        return "On-Treatment"

df_long['Patient_ID'] = df_long['Sample'].apply(assign_patient_id)
df_long['Timepoint'] = df_long['Sample'].apply(simplify_timepoint)

# Drop rows with undefined timepoints (e.g., INNOV)
df_long = df_long.dropna(subset=['Timepoint'])

# Keep only patients with at least 2 timepoints
valid = df_long['Patient_ID'].value_counts()
df_long = df_long[df_long['Patient_ID'].isin(valid[valid > 1].index)]

# Order timepoints
time_order = ['Baseline', 'On-Treatment', 'Post-Treatment']
df_long['Timepoint'] = pd.Categorical(df_long['Timepoint'], categories=time_order, ordered=True)

# Step 4: Plot
unique_patients = df_long['Patient_ID'].unique()
palette_cb = sns.color_palette("colorblind", len(unique_patients))
markers = ['o', 's', 'D', '^', 'v', 'P', '*', 'X', 'H', '8', '<', '>']
linestyles = ['-', '--', '-.', ':']

fig, ax = plt.subplots(figsize=(10, 6))
for i, (pid, group) in enumerate(df_long.groupby('Patient_ID')):
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

# Step 5: Save
plot_path = os.path.join(plot_dir, "methylation_longitudinal_plot.png")
fig.savefig(plot_path, dpi=300, bbox_inches='tight')
print(f"✅ Saved plot to: {plot_path}")
