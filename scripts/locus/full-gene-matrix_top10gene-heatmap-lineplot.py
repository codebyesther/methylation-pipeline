# Gene matrix for all genes
gene_means = {}
for gene in tqdm(multicpg_genes, desc="Computing gene matrix for all genes"):
    cpgs = gene_annot[gene_annot['gene_name'] == gene]['cgi_id']
    gene_data = cpg_matrix.loc[cpg_matrix.index.intersection(cpgs)]
    if not gene_data.empty:
        gene_means[gene] = gene_data.mean(axis=0)

# Create a gene matrix DataFrame with all genes
gene_matrix = pd.DataFrame(gene_means).T

# Save the full gene matrix to a CSV file
gene_matrix.to_csv(os.path.join(args.output_dir, "full_gene_matrix.csv"))

# Per-patient gene matrix export for all genes
for patient in tqdm(patient_ids, desc="Saving per-patient matrices for all genes"):
    sample_cols = [s for s in cpg_matrix.columns if patient in s]
    if not sample_cols:
        continue
    patient_gene_matrix = gene_matrix[sample_cols].T
    patient_gene_matrix.dropna(axis=1, how='all').T.to_csv(
        os.path.join(args.output_dir, f"methylation_matrix_{patient}_all_genes.csv"))

# Combine all patient methylation matrices into one CSV
patient_matrix_files = [os.path.join(args.output_dir, f"methylation_matrix_{patient}_all_genes.csv") for patient in patient_ids]
combined_matrix = pd.concat([pd.read_csv(file, index_col=0) for file in patient_matrix_files], axis=1)
combined_matrix.to_csv(os.path.join(args.output_dir, "combined_methylation_matrix_all_genes.csv"))

# **Filter for the top 10 genes AFTER creating the full gene matrix**
top_genes = pd.Series(baseline_post).abs().sort_values(ascending=False).head(10).index.tolist()
filtered_gene_matrix = gene_matrix.loc[top_genes]

# Save the top 10 gene matrix to a CSV file
filtered_gene_matrix.to_csv(os.path.join(args.output_dir, "top10_gene_matrix.csv"))

# Visualizations and analysis for top 10 genes
sample_timepoints = {sample: classify_timepoint(sample) for sample in cpg_matrix.columns}
timepoint_df = pd.DataFrame.from_dict(sample_timepoints, orient='index', columns=['Timepoint'])
gene_matrix_T = filtered_gene_matrix.T
merged = gene_matrix_T.merge(timepoint_df, left_index=True, right_index=True)
avg_by_tp = merged.groupby("Timepoint").mean().T

# Reorder timepoints
tp_order = ['Healthy', 'Baseline', 'On-Treatment', 'Post-Treatment']
avg_by_tp = avg_by_tp[[tp for tp in tp_order if tp in avg_by_tp.columns]]

# Heatmap for top 10 genes
fig_height = min(20, len(avg_by_tp))
plt.figure(figsize=(15, fig_height))
ax = sns.heatmap(avg_by_tp, cmap="coolwarm", cbar_kws={'label': 'Methylation', 'shrink': 1})
ax.set_xlabel("Timepoint", fontsize=14)
ax.set_ylabel("Gene", fontsize=14)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=14)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=14)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=14)
cbar.set_label("Scaled Methylation Fragment Count Ratio", fontsize=14)
plt.title("Average Methylation per Gene Across Timepoints (Top 10 by ∆ Baseline → Post-Tx)", fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "avg_methylation_heatmap_top10_genes.png"))
plt.close()

# Average Line plot for top 10 genes
plt.figure(figsize=(15, 6))
for gene in avg_by_tp.index:
    plt.plot(avg_by_tp.columns, avg_by_tp.loc[gene], label=gene)
plt.title("Methylation Trends Across Timepoints (Top 10 by ∆ Baseline → Post-Tx)", fontsize=14)
plt.ylabel("Average Scaled Methylation Fragment Count Ratio", fontsize=14)
plt.xlabel("Timepoint", fontsize=14)
plt.xticks(rotation=0, fontsize=14)
plt.yticks(fontsize=14)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(args.output_dir, "avg_methylation_lineplot_top10_genes.png"))
plt.close()
