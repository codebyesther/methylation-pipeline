"""Step 3: turn fragment counts into depth-scaled methylation signals.

Numerator = fragments from CGIs that are lowly methylated in healthy cells
            (default file keyword "Glob20"): unexpected methylation.
Reference = fragments from CGIs that are normally methylated in healthy cells
            (default keyword "GlobMin80"): grows with input and depth, not
            with tumor signal.

Two scaling modes:

  --mode total   (default; the 2024 HIV-lymphoma method)
      CGI signal   = numerator count of that CGI / sample's total reference
                     count x 1e6
      sample score = sample's total numerator / total reference x 1e5
      Every CGI in a sample is divided by the same number, so no division by
      zero occurs unless the whole sample has no reference fragments.

  --mode per-cgi (the 2025 treatment-monitoring method)
      CGI signal = numerator count / reference count of the same CGI x 1e5.
      Cells whose reference count is 0 are left empty and logged.

Outputs (output/):
  scaled_signal_matrix.csv   CGI x sample
  sample_scores.csv          one row per sample: subject, window, totals,
                             scaled score
  zero_reference_log.csv     (per-cgi mode) CGI/sample pairs left empty
"""

import argparse
import glob
import os

import numpy as np
import pandas as pd

from edcommon import (DATA_DIR, OUTPUT_DIR, find_file, keyword_in, load_count_file,
                      read_id_list, subject_of, window_of)


def pick_input(out_dir: str, keyword: str) -> str:
    merged = os.path.join(out_dir, f"merged_{keyword}.xlsx")
    if os.path.exists(merged):
        return merged
    hits = [p for p in sorted(glob.glob(os.path.join(out_dir, "*_filtered.xlsx")))
            if keyword_in(os.path.basename(p), keyword)]
    if len(hits) == 1:
        return hits[0]
    if not hits:
        raise FileNotFoundError(f"No filtered or merged '{keyword}' file in {out_dir}/ (run steps 1-2)")
    raise RuntimeError(f"Several '{keyword}' runs found; run step 2 to merge them first")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=DATA_DIR)
    ap.add_argument("--out", default=OUTPUT_DIR)
    ap.add_argument("--numerator", default="Glob20")
    ap.add_argument("--reference", default="GlobMin80")
    ap.add_argument("--mode", choices=["total", "per-cgi"], default="total")
    args = ap.parse_args()

    num_cgi, num_tot = load_count_file(pick_input(args.out, args.numerator))
    ref_cgi, ref_tot = load_count_file(pick_input(args.out, args.reference))

    samples = num_cgi.columns.intersection(ref_cgi.columns)
    dropped = set(num_cgi.columns) ^ set(ref_cgi.columns)
    if dropped:
        print(f"[warn] samples present in only one of the two files, skipped: {sorted(dropped)}")
    num_cgi, ref_cgi = num_cgi[samples], ref_cgi[samples]

    # Use the lab's summary row when present; otherwise sum the CGI rows.
    if num_tot is None:
        print("[note] numerator file has no total row; using the sum of CGI rows")
        num_tot = num_cgi.sum()
    if ref_tot is None:
        print("[note] reference file has no total row; using the sum of CGI rows")
        ref_tot = ref_cgi.sum()
    num_tot, ref_tot = num_tot[samples].astype(float), ref_tot[samples].astype(float)

    if args.mode == "total":
        denom = ref_tot.replace(0, np.nan)
        signal = num_cgi.div(denom, axis=1) * 1e6
        log = None
    else:
        rows = num_cgi.index.intersection(ref_cgi.index)
        print(f"per-cgi mode: {len(rows)} CGIs appear in both files")
        if len(rows) == 0:
            raise SystemExit("No CGI is in both files; per-cgi mode needs shared CGI rows. Use --mode total.")
        n, r = num_cgi.loc[rows], ref_cgi.loc[rows]
        signal = n / r.replace(0, np.nan) * 1e5
        zero = signal.isna() & (n > 0)
        log = zero.stack().loc[lambda s: s].reset_index().iloc[:, :2]
        log.columns = ["CGI", "Sample"]
        signal = signal.where(~(signal.isna() & (n == 0)), 0.0)  # 0/0 -> 0

    bad = ref_tot[ref_tot == 0].index.tolist()
    if bad:
        print(f"[warn] samples with zero reference fragments (left empty): {bad}")

    ids_path = None
    try:
        ids_path = find_file(args.data, "subject", exts=(".xlsx", ".csv"))
    except FileNotFoundError:
        try:
            ids_path = find_file(args.data, "patient", exts=(".xlsx", ".csv"))
        except FileNotFoundError:
            pass
    ids = read_id_list(ids_path) if ids_path else None

    scores = pd.DataFrame({
        "sample": samples,
        "subject": [subject_of(s, ids) for s in samples],
        "window": [window_of(s) for s in samples],
        "numerator_total": num_tot.values,
        "reference_total": ref_tot.values,
    })
    scores["scaled_score"] = scores["numerator_total"] / scores["reference_total"].replace(0, np.nan) * 1e5

    signal.to_csv(os.path.join(args.out, "scaled_signal_matrix.csv"))
    scores.to_csv(os.path.join(args.out, "sample_scores.csv"), index=False)
    print(f"mode={args.mode}: {signal.shape[0]} CGIs x {signal.shape[1]} samples")
    print(scores.groupby("window", dropna=False)["sample"].count().rename("samples").to_string())
    if log is not None:
        log.to_csv(os.path.join(args.out, "zero_reference_log.csv"), index=False)
        print(f"{len(log)} CGI/sample cells had numerator > 0 and reference = 0 (left empty, logged)")


if __name__ == "__main__":
    main()
