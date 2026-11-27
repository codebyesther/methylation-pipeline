"""Gene-level analysis: which genes carry the most signal earliest, and how
do they change toward diagnosis compared with healthy controls?

Cohort mode (default):
  Takes the top N genes by mean signal in the ranking window (default: the
  earliest pre-diagnosis window, pre12) and follows the same genes through
  the later windows, with the healthy-control mean alongside.

Subject mode (--subject ID):
  Same, but ranks and plots one subject's own samples. Useful when windows
  are collected from a single person.

Outputs:
  output/top<N>_genes_<rank-window>[_<subject>].csv
  plots/top<N>_heatmap_<...>.png      genes x windows, annotated
  plots/top<N>_lines_<...>.png        one line per gene across windows

Ranking by mean signal favors genes covered by many CGIs (signals are summed
per gene). Use --per-cgi to rank by mean signal per CGI instead.
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from edcommon import OUTPUT_DIR, WINDOW_ORDER, save_fig, window_label  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=OUTPUT_DIR)
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--rank-window", default="pre12", choices=["pre12", "pre6", "near"])
    ap.add_argument("--subject", default=None, help="rank and plot one subject's samples only")
    ap.add_argument("--per-cgi", action="store_true", help="divide each gene's signal by its CGI count")
    args = ap.parse_args()

    long = pd.read_csv(os.path.join(args.out, "gene_signal_long.csv"))
    if args.per_cgi:
        cmap = pd.read_csv(os.path.join(args.out, "cgi_gene_map.csv"))
        n = cmap.groupby("gene")["CGI"].nunique()
        long["signal"] = long["signal"] / long["gene"].map(n)

    controls = long[long["window"] == "control"]
    cases = long[long["window"] != "control"]
    if args.subject:
        cases = cases[cases["subject"].astype(str) == args.subject]
        if cases.empty:
            raise SystemExit(f"No pre-diagnosis samples for subject {args.subject}")

    means = cases.groupby(["gene", "window"])["signal"].mean().unstack("window")
    if args.rank_window not in means:
        raise SystemExit(f"No samples in window '{args.rank_window}'")
    top = means[args.rank_window].dropna().sort_values(ascending=False).head(args.top).index

    cols = [w for w in WINDOW_ORDER if w in means.columns and w != "control"]
    table = means.loc[top, cols]
    if not controls.empty:
        table.insert(0, "control", controls.groupby("gene")["signal"].mean().reindex(top))
    table.columns = [window_label(c) for c in table.columns]

    tag = f"top{args.top}_genes_{args.rank_window}" + (f"_{args.subject}" if args.subject else "")
    table.to_csv(os.path.join(args.out, f"{tag}.csv"))
    print(table.round(1).head(10).to_string())

    who = f"subject {args.subject}" if args.subject else "cohort mean"
    unit = "signal per CGI" if args.per_cgi else "summed CGI signal"

    # heatmap
    fig, ax = plt.subplots(figsize=(1.6 + 1.3 * table.shape[1], 0.32 * len(table) + 1.5))
    sns.heatmap(table, cmap="Blues", linewidths=0.5, annot=True, fmt=".0f",
                annot_kws={"size": 7}, cbar_kws={"label": unit}, ax=ax)
    ax.set_title(f"Top {len(table)} genes in '{window_label(args.rank_window)}' ({who})")
    ax.set_ylabel("")
    plt.xticks(rotation=15)
    save_fig(fig, f"{tag.replace('_genes', '_heatmap')}.png")
    plt.close(fig)

    # line chart: top third dotted, middle dashed, bottom solid (as in 2024)
    pre = [c for c in table.columns if c != window_label("control")]
    fig, ax = plt.subplots(figsize=(8, 5.5))
    third = max(1, len(table) // 3)
    for i, (gname, row) in enumerate(table[pre].iterrows()):
        style = ":" if i < third else ("--" if i < 2 * third else "-")
        ax.plot(range(len(pre)), row.values, linestyle=style, marker="o", markersize=3, label=gname)
    ax.set_xticks(range(len(pre)), pre)
    ax.set_xlabel("Blood collection window")
    ax.set_ylabel(unit)
    ax.set_title(f"Top {len(table)} genes toward diagnosis ({who})")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=6, ncol=2, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    save_fig(fig, f"{tag.replace('_genes', '_lines')}.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
