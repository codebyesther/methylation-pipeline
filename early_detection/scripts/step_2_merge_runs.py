"""Step 2 (optional): merge filtered files from several sequencing runs.

Merges every output/*_filtered.xlsx whose name contains the keyword into one
table, aligned by CGI name. CGIs missing from a run count as 0 for that run's
samples. If the same sample appears in two runs, the later file is kept and a
warning is printed.

Outputs: output/merged_<keyword>.xlsx for the numerator and reference keywords.
"""

import argparse
import glob
import os

import pandas as pd

from edcommon import OUTPUT_DIR, keyword_in, load_count_file, write_count_file


def merge(keyword: str, out_dir: str) -> None:
    paths = [p for p in sorted(glob.glob(os.path.join(out_dir, "*_filtered.xlsx")))
             if keyword_in(os.path.basename(p), keyword)
             and not os.path.basename(p).lower().startswith("merged_")]
    if not paths:
        raise FileNotFoundError(f"No *_filtered.xlsx containing '{keyword}' in {out_dir}/")

    cgis, totals = [], []
    for p in paths:
        c, t = load_count_file(p)
        dup = set(c.columns) & set().union(*[set(x.columns) for x in cgis]) if cgis else set()
        if dup:
            print(f"  [warn] samples in more than one run, keeping {os.path.basename(p)}: {sorted(dup)}")
            cgis = [x.drop(columns=list(dup & set(x.columns))) for x in cgis]
            totals = [None if x is None else x.drop(list(dup & set(x.index))) for x in totals]
        cgis.append(c)
        totals.append(t)

    merged = pd.concat(cgis, axis=1).fillna(0)
    tot = None if any(t is None for t in totals) else pd.concat(totals)
    out = os.path.join(out_dir, f"merged_{keyword}.xlsx")
    write_count_file(merged, tot, out)
    print(f"{keyword}: merged {len(paths)} run(s) -> {merged.shape[1]} samples, {merged.shape[0]} CGIs")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=OUTPUT_DIR)
    ap.add_argument("--numerator", default="Glob20")
    ap.add_argument("--reference", default="GlobMin80")
    args = ap.parse_args()
    merge(args.numerator, args.out)
    merge(args.reference, args.out)


if __name__ == "__main__":
    main()
