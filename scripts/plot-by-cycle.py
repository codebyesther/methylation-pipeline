# this is a temporary script to create line plots of manually assembled data (Scaled_LOI-in-EMseq-16-18-20_by-Cycle.xlsx)
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your data (replace 'your_data.csv' with the actual data file)
# Make sure the data is structured with rows as patients and columns containing timepoints
data_file = "Scaled_LOI-in-EMseq-16-18-20_by-Cycle.xlsx"  # Update with your data file path
data = pd.read_excel(data_file)

# Example structure of the input data (data.head()):
# Patient ID | Baseline | C1-2 | C3-4 | C5+
# ------------------------------------------
# Patient 1  | 0.8      | 1.2  | 1.0  | 1.4
# Patient 2  | 0.6      | 1.1  | 0.9  | 1.3
# ...

# Melt the data for easier plotting
melted_data = data.melt(id_vars=["Patient ID"], var_name="Timepoint", value_name="Scaled Fragment Count Ratio")

# Create a line plot using seaborn
plt.figure(figsize=(12, 6))
sns.lineplot(data=melted_data, x="Timepoint", y="Scaled Fragment Count Ratio", hue="Patient ID", marker="o")

# Customize the plot
plt.title("Scaled Fragment Count Ratio Over Timepoints", fontsize=16)
plt.xlabel("Timepoint", fontsize=12)
plt.ylabel("Scaled Fragment Count Ratio", fontsize=12)
plt.legend(title="Patient ID", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()

# Save the plot or display it
output_file = "line_plot_by_patient.png"  # Update with your desired output file name
plt.savefig(output_file)
plt.show()
