"""Step 1: keep only the samples that belong to the study.

Inputs (data/):
  * one ID list whose file name contains "subject" (or "patient"):
    first column = subject IDs (cases and healthy controls)
  * the lab's CGI count files, as pairs: a numerator file (default keyword
    "Glob20") and a reference file (default keyword "GlobMin80")

Output (output/):
  * <original name>_filtered.xlsx for every count file, same layout

A sample (column) is kept if its name contains one of the listed IDs and it
can be assigned to a collection window (see WINDOWS in edcommon.py). Samples
that match an ID but no window are reported, not silently dropped.
"""

import argparse
import glob
import os

from edcommon import (DATA_DIR, OUTPUT_DIR, contains_id, ensure_dirs, load_count_file,
                      read_id_list, window_of, write_count_file)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=DATA_DIR)
    ap.add_argument("--out", default=OUTPUT_DIR)
    args = ap.parse_args()
    ensure_dirs(args.out)

    files = sorted(glob.glob(os.path.join(args.data, "*.xlsx")))
    id_files = [f for f in files if any(k in os.path.basename(f).lower() for k in ("subject", "patient"))]
    if not id_files:
        raise FileNotFoundError("Put an ID list named like 'subject_ids.xlsx' in data/")
    ids = read_id_list(id_files[0])
    print(f"{len(ids)} subject IDs from {os.path.basename(id_files[0])}")

    count_files = [f for f in files if f not in id_files and not os.path.basename(f).startswith("~$")]
    if not count_files:
        raise FileNotFoundError("No count files (.xlsx) found in data/")

    for path in count_files:
        cgi, totals = load_count_file(path)
        keep, no_window = [], []
        for col in cgi.columns:
            if any(contains_id(col, sid) for sid in ids):
                (keep if window_of(col) else no_window).append(col)
        name = os.path.splitext(os.path.basename(path))[0]
        out_path = os.path.join(args.out, f"{name}_filtered.xlsx")
        write_count_file(cgi[keep], None if totals is None else totals[keep], out_path)
        print(f"{os.path.basename(path)}: kept {len(keep)} of {cgi.shape[1]} samples, {cgi.shape[0]} CGIs")
        if no_window:
            print(f"  [warn] matched an ID but no collection window, skipped: {no_window}")


if __name__ == "__main__":
    main()
