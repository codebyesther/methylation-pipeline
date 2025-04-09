def plot_multi_cpg_genes(df, title, filename):
    df["CpG_Island"] = df["CpG_Island"].astype(str).str.strip()

    # Convert to chr:start-end format to match gene mapping keys
    matches = df["CpG_Island"].str.extract(r"CGI_(chr\d+)_(\d+)_(\d+)_")
    df["Map_ID"] = matches[0] + ':' + matches[1] + '-' + matches[2]

    df["Gene"] = df["Map_ID"].map(cgi_to_gene)
    df = df.dropna(subset=["Gene"])

    print(f"[DEBUG] {df['Gene'].nunique()} unique genes mapped in: {title}")
    print(df[["CpG_Island", "Map_ID", "Gene"]].dropna().head(5))

    # Aggregate gene-level stats
    gene_df = df.groupby("Gene").agg(
        count=("CpG_Island", "count"),
        avg_delta=("Avg_Delta", "mean")
    ).reset_index()

    print("\n[DEBUG] Full gene count summary:")
    print(gene_df.sort_values("count", ascending=False).head(10))

    # Primary: genes with >1 CpG island
    multi_cpg_genes = gene_df[gene_df["count"] > 1].sort_values("avg_delta")

    if multi_cpg_genes.empty:
        print(f"⚠️ No genes with >1 CpG island — plotting top 10 genes with highest CpG count instead.")
        multi_cpg_genes = gene_df.sort_values("count", ascending=False).head(10)

    print("[DEBUG] Genes selected for plotting:")
    print(multi_cpg_genes)

    # Plot
    plt.figure(figsize=(10, 6))
    sns.barplot(data=multi_cpg_genes, x="avg_delta", y="Gene", color='darkblue')
    plt.xlabel("Avg Change in Scaled Methylated Fragment Count Ratio")
    plt.axvline(0, color="gray", linestyle="--")
    plt.title(title, loc='center')
    plt.tight_layout()
    plot_path = os.path.join(args.outdir, filename)
    plt.savefig(plot_path)
    plt.close()
    top10dmplot_filenames.append(plot_path)
