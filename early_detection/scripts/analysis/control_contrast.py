"""Gene-level contrast with healthy controls.

For early detection the question is not only "which genes are high before
diagnosis" but "which genes are high before diagnosis AND low in people
without cancer". For every gene this script reports, per pre-diagnosis
window:

  mean signal in cases and in controls
  log2 fold change vs controls (pseudocount 1)
  fraction of case samples above the highest control value

Outputs:
  output/control_contrast.csv
  output/control_contrast_<window>.rnk    two-column ranked list (gene,
                                          log2FC) for GSEA preranked
                                          (GSEA desktop or gseapy.prerank)
  plots/control_contrast_<window>.png     log2FC vs fraction above control

Run step 1 with healthy control IDs in the subject list; without controls
this script stops.
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from edcommon import OUTPUT_DIR, PRE_DX_KEYS, WINDOW_COLORS, save_fig, window_label  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=OUTPUT_DIR)
    ap.add_argument("--min-mean", type=float, default=1.0,
                    help="skip genes whose mean signal is below this in every group")
    ap.add_argument("--label-top", type=int, default=10)
    args = ap.parse_args()

    long = pd.read_csv(os.path.join(args.out, "gene_signal_long.csv"))
    ctrl = long[long["window"] == "control"]
    if ctrl.empty:
        raise SystemExit("No healthy-control samples found; add control IDs to the subject list.")

    c_mean = ctrl.groupby("gene")["signal"].mean()
    c_max = ctrl.groupby("gene")["signal"].max()
    out = pd.DataFrame({"control_mean": c_mean, "control_max": c_max, "n_control": ctrl["sample"].nunique()})

    windows = [w for w in PRE_DX_KEYS if w in set(long["window"])]
    if not windows:
        raise SystemExit("No pre-diagnosis samples found (windows pre12, pre6, near).")
    for w in windows:
        g = long[long["window"] == w]
        m = g.groupby("gene")["signal"].mean()
        above = (g.assign(above=g["signal"] > g["gene"].map(c_max))
                 .groupby("gene")["above"].mean())
        out[f"{w}_mean"] = m
        out[f"{w}_log2FC"] = np.log2((m + 1) / (out["control_mean"] + 1))
        out[f"{w}_frac_above_control_max"] = above
        out[f"n_{w}"] = g["sample"].nunique()

    keep = out[[c for c in out if c.endswith("_mean")]].max(axis=1) >= args.min_mean
    out = out[keep].sort_values(f"{windows[0]}_log2FC", ascending=False)
    out.index.name = "gene"
    out.to_csv(os.path.join(args.out, "control_contrast.csv"))
    print(f"{len(out)} genes with mean signal >= {args.min_mean} in some group")

    for w in windows:
        rnk = out[f"{w}_log2FC"].dropna().sort_values(ascending=False)
        rnk.to_csv(os.path.join(args.out, f"control_contrast_{w}.rnk"), sep="\t", header=False)

        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        x, y = out[f"{w}_log2FC"], out[f"{w}_frac_above_control_max"]
        ax.scatter(x, y, s=10, alpha=0.5, color=WINDOW_COLORS[w], edgecolor="none")
        best = (x.clip(lower=0) * y).sort_values(ascending=False).head(args.label_top).index
        placed = []  # skip labels that would sit on top of one already drawn
        xr = (x.max() - x.min()) or 1.0
        for gname in best:
            px, py = x[gname], y[gname]
            if any(abs(px - qx) < 0.08 * xr and abs(py - qy) < 0.04 for qx, qy in placed):
                continue
            placed.append((px, py))
            ax.annotate(gname, (px, py), fontsize=7, xytext=(3, 2), textcoords="offset points")
        ax.axvline(0, color="grey", linewidth=0.8)
        ax.set_xlabel(f"log2 fold change vs healthy controls ({window_label(w)})")
        ax.set_ylabel("Fraction of case samples above control max")
        ax.set_title(f"Genes separating '{window_label(w)}' from controls "
                     f"(n={int(out[f'n_{w}'].iloc[0])} vs {int(out['n_control'].iloc[0])})")
        ax.spines[["top", "right"]].set_visible(False)
        save_fig(fig, f"control_contrast_{w}.png")
        plt.close(fig)

    print(out.head(args.label_top).round(2).to_string())


if __name__ == "__main__":
    main()
