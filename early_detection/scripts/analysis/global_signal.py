"""Sample-level analysis: does the overall scaled score rise toward diagnosis?

Reads output/sample_scores.csv and writes:

  plots/global_score_by_window.png      every sample's scaled score by window
  plots/global_score_trajectories.png   one line per subject across the
                                        pre-diagnosis windows, with the
                                        healthy-control range shaded
  output/global_window_stats.csv        each pre-diagnosis window vs controls:
                                        n, median, Mann-Whitney U, two-sided
                                        p, and AUC = U / (n_case x n_control)

Read the statistics as descriptive when groups are small: with fewer than
about 5 samples per group, a p-value cannot reach conventional thresholds
and the AUC has a very wide uncertainty. The script prints a warning then.
"""

import argparse
import os
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from edcommon import (OUTPUT_DIR, PRE_DX_KEYS, WINDOW_COLORS, WINDOW_ORDER,  # noqa: E402
                      save_fig, window_label)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=OUTPUT_DIR)
    ap.add_argument("--log", action="store_true", help="log10 y-axis")
    args = ap.parse_args()

    s = pd.read_csv(os.path.join(args.out, "sample_scores.csv")).dropna(subset=["window", "scaled_score"])
    present = [w for w in WINDOW_ORDER if w in set(s["window"])]
    ctrl = s.loc[s["window"] == "control", "scaled_score"]

    # --- statistics -------------------------------------------------------
    rows = []
    for w in [k for k in PRE_DX_KEYS if k in present]:
        case = s.loc[s["window"] == w, "scaled_score"]
        row = {"window": w, "label": window_label(w), "n_case": len(case), "n_control": len(ctrl),
               "median_case": case.median(), "median_control": ctrl.median() if len(ctrl) else np.nan}
        if len(case) and len(ctrl):
            u, p = mannwhitneyu(case, ctrl, alternative="two-sided")
            row.update(U=u, p_two_sided=p, AUC=u / (len(case) * len(ctrl)))
        rows.append(row)
    stats = pd.DataFrame(rows)
    stats.to_csv(os.path.join(args.out, "global_window_stats.csv"), index=False)
    print(stats.round(4).to_string(index=False))
    if len(stats) and (stats["n_case"].min() < 5 or len(ctrl) < 5):
        print("[warn] fewer than 5 samples in a group: treat p-values and AUCs as descriptive only")

    # --- strip + box by window ---------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 4.5))
    data = [s.loc[s["window"] == w, "scaled_score"].values for w in present]
    ax.boxplot(data, widths=0.5, showfliers=False, medianprops={"color": "black"})
    rng = np.random.default_rng(0)
    for i, (w, vals) in enumerate(zip(present, data), start=1):
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(vals)), vals, s=28,
                   color=WINDOW_COLORS[w], edgecolor="black", linewidth=0.4, zorder=3)
    ax.set_xticks(range(1, len(present) + 1), [textwrap.fill(window_label(w), 12) + f"\n(n={len(d)})" for w, d in zip(present, data)])
    ax.set_ylabel("Scaled score (Glob20 / GlobMin80 x 1e5)")
    ax.set_title("Sample-level methylation score by collection window")
    if args.log:
        ax.set_yscale("log")
    ax.spines[["top", "right"]].set_visible(False)
    save_fig(fig, "global_score_by_window.png")
    plt.close(fig)

    # --- per-subject trajectories ------------------------------------------
    pre = [w for w in PRE_DX_KEYS if w in present]
    if pre:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        x = {w: i for i, w in enumerate(pre)}
        if len(ctrl):
            ax.axhspan(ctrl.min(), ctrl.max(), color=WINDOW_COLORS["control"], alpha=0.18,
                       label=f"Healthy control range (n={len(ctrl)})")
        cases = s[s["window"].isin(pre)]
        for subj, g in cases.groupby("subject"):
            g = g.groupby("window")["scaled_score"].mean().reindex(pre).dropna()
            ax.plot([x[w] for w in g.index], g.values, marker="o", linewidth=1.4, label=subj)
        ax.set_xticks(range(len(pre)), [window_label(w) for w in pre])
        ax.set_ylabel("Scaled score")
        ax.set_title("Per-subject trajectory before diagnosis")
        if args.log:
            ax.set_yscale("log")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
        save_fig(fig, "global_score_trajectories.png")
        plt.close(fig)


if __name__ == "__main__":
    main()
