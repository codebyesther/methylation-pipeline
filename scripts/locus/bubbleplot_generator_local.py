# === Bubble plots per chromosome (averaged across patients) ===
for chrom in tqdm(coords_df["Chr"].unique(), desc="Generating bubble plots per chromosome"):
    chr_data = bubble_data[bubble_data["Chr"] == chrom]

    all_rows = []
    for tp in timepoints_chromosome:
        # Find columns that end with e.g. "_Baseline"
        cols = [c for c in chr_data.columns if c.endswith(f"_{tp}")]
        if not cols:
            continue
        avg = chr_data[cols].mean(axis=1)
        sub = chr_data.loc[avg.notna(), ["Midpoint"]].copy()
        sub["value"] = avg[avg.notna()]
        sub["Timepoint"] = tp
        all_rows.append(sub)

    if not all_rows:
        continue

    subset_df = pd.concat(all_rows, ignore_index=True)

    fig = plt.figure(figsize=(14, 6))
    gs = GridSpec(nrows=1, ncols=2, width_ratios=[5, 1], figure=fig)

    ax_main = fig.add_subplot(gs[0, 0])
    gs_right = gs[0, 1].subgridspec(nrows=2, ncols=1, height_ratios=[0.5, 0.5])
    ax_cbar = fig.add_subplot(gs_right[0, 0])
    ax_legend = fig.add_subplot(gs_right[1, 0])

    sc = None
    for _, row in subset_df.iterrows():
        tp = row["Timepoint"]
        y = timepoint_positions_chromosome[tp]
        bubble_size = row["value"]**0.5 * 50
        sc = ax_main.scatter(
            row["Midpoint"],
            y,
            s=bubble_size,
            c=row["value"],
            cmap="viridis",
            alpha=0.6,
            vmin=0,  # lower bound of color scale
            vmax=170    # upper bound of color scale
        )

    ax_main.set_yticks(list(timepoint_positions_chromosome.values()))  # Ensure the number of ticks matches the number of labels
    ax_main.set_yticklabels(timepoints_chromosome)
    ax_main.set_ylim(0.75, 1.25) # set y-limits so large bubbles have padding above & below
    ax_main.set_xlabel("CpG Island Genomic Coordinate Midpoint (bp)")
    ax_main.set_ylabel("Timepoint")
    ax_main.set_title(f"DNA Hypermethylation Profiles Throughout Treatment (Averaged Across Patients)\nChromosome: {chrom}")

    if sc is not None:
        fig.colorbar(sc, cax=ax_cbar, label="Fragment Count")

    ax_legend.axis("off")
    handles, labels = [], []
    for s in [1, 15, 150]:
        h = ax_legend.scatter([], [], s=s**0.5 * 50, color="gray", alpha=0.5)
        handles.append(h)
        labels.append(str(s))

    ax_legend.legend(
        handles,
        labels,
        title="Bubble Size\n(Fragment Count)",
        labelspacing=1.5,
        handletextpad=2.0,
        borderpad=1.3,
        loc="center",
        bbox_to_anchor=(0.5, 0.5)
    )

    plt.tight_layout()
    filename_base = os.path.join("plots", f"bubble_{chrom}")
    plt.savefig(f"{filename_base}.png")
    plt.savefig(f"{filename_base}.svg")
    plt.close()
