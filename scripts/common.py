"""Shared helpers for the treatment-monitoring scripts.

Scripts in scripts/, scripts/global/ and scripts/locus/ import this module so
that file discovery, patient/timepoint parsing, gene matching and reading the
Step 3 matrix behave the same way everywhere. Run every script from the
repository root, as before.
"""

from __future__ import annotations

import glob
import os
import re

import numpy as np
import pandas as pd

DATA_DIR = "data"
OUTPUT_DIR = "output"

TOTAL_ROW = "Total CpG island fragments counts for this particular spreadsheet"
GLOB20_TOTAL_ROW = "Total CpG island fragments counts for Glob20"
GLOBMIN80_TOTAL_ROW = "Total CpG island fragments counts for GlobMin80"
SCALED_TOTAL_ROW = "Scaled ratio of totals (Glob20 / GlobMin80 x 100000)"
LEGACY_PATIENT_FILE = "Patient ID list fot EMseq16-18-20.xlsx"


# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------
def find_file(directory: str, keyword: str, exts=(".xlsx", ".xls", ".csv"), prefer: str | None = None) -> str:
    """Find a file whose name contains `keyword`, ignoring case.

    Excel lock files (~$...) are skipped. If several files match, the one
    whose name also contains `prefer` wins; otherwise the first in
    alphabetical order, with a note so the choice is visible.
    """
    hits = sorted(
        f for f in glob.glob(os.path.join(directory, "*"))
        if keyword.lower() in os.path.basename(f).lower()
        and f.lower().endswith(exts)
        and not os.path.basename(f).startswith("~$")
    )
    if not hits:
        raise FileNotFoundError(f"No file containing '{keyword}' found in '{directory}/'")
    if prefer:
        preferred = [f for f in hits if prefer.lower() in os.path.basename(f).lower()]
        if preferred:
            hits = preferred
    if len(hits) > 1:
        print(f"[note] {len(hits)} files contain '{keyword}'; using {os.path.basename(hits[0])}")
    return hits[0]


def find_scaled_matrix(directory: str = OUTPUT_DIR) -> str:
    """The Step 3 output (scaled_fragment_ratios_matrix.xlsx)."""
    return find_file(directory, "fragment_ratios_matrix", exts=(".xlsx",))


def find_patient_file(directory: str = DATA_DIR) -> str:
    """The patient ID list: any .xlsx/.csv with 'patient' in its name."""
    return find_file(directory, "patient", exts=(".xlsx", ".xls", ".csv"), prefer="patient id")


def read_table(path: str, **kw) -> pd.DataFrame:
    return pd.read_excel(path, **kw) if path.lower().endswith((".xlsx", ".xls")) else pd.read_csv(path, **kw)


def read_patient_ids(path: str | None = None) -> list[str]:
    """Patient IDs from the first column of the patient list (header row skipped)."""
    path = path or find_patient_file()
    ids = read_table(path).iloc[:, 0].dropna().astype(str).str.strip()
    return [i for i in ids if i]


# ---------------------------------------------------------------------------
# Samples: patient and timepoint
# ---------------------------------------------------------------------------
def contains_id(sample: str, pid: str) -> bool:
    """True if the patient ID appears in the sample name as a whole token,
    so "PT001" matches "PT001_Baseline" but "PT01" does not match "PT010_C2D1"."""
    if not pid:
        return False
    return re.search(r"(?<![A-Za-z0-9])" + re.escape(pid) + r"(?![A-Za-z0-9])", str(sample)) is not None


def get_patient(sample: str, patient_ids: list[str]) -> str | None:
    """Longest patient ID found in the sample name as a whole token."""
    for pid in sorted(patient_ids, key=len, reverse=True):
        if contains_id(sample, pid):
            return pid
    return None


def is_healthy(sample: str) -> bool:
    return "INNOV" in str(sample).upper()


def is_off_tx(sample: str) -> bool:
    return re.search(r"off[-_ ]?tx", str(sample), re.IGNORECASE) is not None


def classify_timepoint(sample: str) -> str:
    """Healthy / Baseline / On-Treatment / Post-Treatment."""
    s = str(sample)
    if is_healthy(s):
        return "Healthy"
    if "baseline" in s.lower():
        return "Baseline"
    if is_off_tx(s):
        return "Post-Treatment"
    return "On-Treatment"


def detailed_timepoint(sample: str) -> str:
    """Healthy / Baseline / C<n> / Off-Tx (On-Treatment if no cycle is found)."""
    grouped = classify_timepoint(sample)
    if grouped == "Post-Treatment":
        return "Off-Tx"
    if grouped != "On-Treatment":
        return grouped
    m = re.search(r"(?<![A-Za-z])C(\d{1,2})(?!\d)", str(sample))
    return f"C{int(m.group(1))}" if m else "On-Treatment"


def timepoint_sort_key(tp: str) -> int:
    if tp == "Healthy":
        return -2
    if tp == "Baseline":
        return -1
    if tp in ("Off-Tx", "Post-Treatment"):
        return 99
    if tp.startswith("C") and tp[1:].isdigit():
        return int(tp[1:])
    return 98


# ---------------------------------------------------------------------------
# CGI names and genes
# ---------------------------------------------------------------------------
def cgi_genes(cgi_name: str) -> list[str]:
    """Gene symbols in a CGI name: CGI_chr1_778604_779167_GENE1_GENE2_0."""
    parts = str(cgi_name).split("_")
    if parts and parts[0] == "CGI":
        parts = parts[1:]
    if parts and parts[-1].isdigit():
        parts = parts[:-1]
    return [p for p in parts[3:] if p]


def cgi_gene_pairs(cgi_names) -> pd.DataFrame:
    """One row per (CGI, gene) pair, matched by exact gene symbol.

    Replaces substring matching, which assigned e.g. GATA2-AS1's CpG islands
    to GATA2 and APCDD1's to APC.
    """
    rows = [(c, g) for c in cgi_names for g in cgi_genes(c)]
    return pd.DataFrame(rows, columns=["cgi_id", "gene_name"])


def gene_matrix(cgi_matrix: pd.DataFrame, how: str = "sum") -> pd.DataFrame:
    """Gene x sample matrix from a CGI x sample matrix (exact gene matching)."""
    cgi_rows = [c for c in cgi_matrix.index.astype(str) if c.startswith("CGI_")]
    pairs = cgi_gene_pairs(cgi_rows)
    data = cgi_matrix.loc[pairs["cgi_id"]].apply(pd.to_numeric, errors="coerce")
    data.index = pairs["gene_name"].values
    agg = data.groupby(level=0).sum(min_count=1) if how == "sum" else data.groupby(level=0).mean()
    agg.index.name = "Gene"
    agg.columns.name = "Sample"
    return agg


# ---------------------------------------------------------------------------
# The Step 3 matrix
# ---------------------------------------------------------------------------
def read_scaled_matrix(path: str | None = None):
    """Read scaled_fragment_ratios_matrix.xlsx.

    Returns (cgi, scaled_total):
      cgi           CGI x sample, numeric; "INF" cells (reference count 0,
                    numerator > 0) become NaN so they are left out of means
      scaled_total  per-sample ratio of totals (Glob20 / GlobMin80 x 1e5)

    The scaled total is read from its labeled row if present, otherwise
    computed from the two raw total rows written by Step 3.
    """
    path = path or find_scaled_matrix()
    raw = pd.read_excel(path)
    label = raw.columns[0]
    raw[label] = raw[label].astype(str).str.strip()
    raw = raw.set_index(label)

    # pd.to_numeric reads the text "INF" as infinity; treat it as missing instead
    cgi = (raw[raw.index.str.startswith("CGI_")]
           .apply(pd.to_numeric, errors="coerce")
           .replace([np.inf, -np.inf], np.nan))

    def row(name):
        return pd.to_numeric(raw.loc[name], errors="coerce") if name in raw.index else None

    scaled = row(SCALED_TOTAL_ROW)
    if scaled is None:
        g20, gm80 = row(GLOB20_TOTAL_ROW), row(GLOBMIN80_TOTAL_ROW)
        if g20 is None or gm80 is None:
            raise ValueError(f"{os.path.basename(path)} has no scaled-total row and no raw total rows; rerun Step 3.")
        scaled = g20 / gm80.replace(0, np.nan) * 100000
    scaled.name = "Scaled_Ratio"
    return cgi, scaled
