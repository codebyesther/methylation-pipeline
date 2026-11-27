"""Step 5: group samples by collection window and average per gene.

Reads output/gene_signal_matrix.csv and output/sample_scores.csv and writes:

  gene_signal_long.csv     one row per gene x sample with subject and window
  gene_window_means.csv    gene x window: mean signal across the samples in
                           each window, with the sample count per window in
                           the file header comment
  gene_window_ranks.csv    gene x window: rank of the gene within the window
                           (1 = highest mean signal)

Windows (edit WINDOWS in edcommon.py to change labels or patterns):
  pre12    15-9 months before diagnosis
  pre6     9-3 months before diagnosis
  near     near diagnosis
  control  healthy controls
"""

import argparse
import os

import pandas as pd

from edcommon import OUTPUT_DIR, WINDOW_ORDER


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=OUTPUT_DIR)
    args = ap.parse_args()

    gene = pd.read_csv(os.path.join(args.out, "gene_signal_matrix.csv"), index_col=0)
    scores = pd.read_csv(os.path.join(args.out, "sample_scores.csv"))
    meta = scores.set_index("sample")[["subject", "window"]]

    long = (gene.rename_axis("gene").reset_index()
            .melt(id_vars="gene", var_name="sample", value_name="signal")
            .join(meta, on="sample"))
    long.to_csv(os.path.join(args.out, "gene_signal_long.csv"), index=False)

    present = [w for w in WINDOW_ORDER if w in set(long["window"])]
    means = (long.groupby(["gene", "window"])["signal"].mean()
             .unstack("window").reindex(columns=present))
    means.to_csv(os.path.join(args.out, "gene_window_means.csv"))
    ranks = means.rank(ascending=False, method="min")
    ranks.to_csv(os.path.join(args.out, "gene_window_ranks.csv"))

    counts = meta["window"].value_counts().reindex(present).fillna(0).astype(int)
    print("samples per window:", ", ".join(f"{w}={n}" for w, n in counts.items()))
    print(f"{means.shape[0]} genes summarized")


if __name__ == "__main__":
    main()
