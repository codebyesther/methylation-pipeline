"""Shared helpers for the early-detection methylation pipeline.

Every step script imports from here so that file discovery, CGI-name parsing,
sample-to-window assignment and the normalization rule are defined once.
"""

from __future__ import annotations

import glob
import os
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Folders (all scripts are run from the repository root)
# ---------------------------------------------------------------------------
DATA_DIR = "data"
OUTPUT_DIR = "output"
PLOTS_DIR = "plots"

TOTAL_ROW_KEY = "total cpg island fragments counts"

# ---------------------------------------------------------------------------
# Collection windows relative to diagnosis
# ---------------------------------------------------------------------------
# Each window is matched against the sample (column) name with a regex,
# case-insensitively, in the order listed. The first match wins, so the more
# specific patterns come first. Edit this table if your lab labels samples
# differently.


@dataclass(frozen=True)
class Window:
    key: str      # short machine name, used in file names
    label: str    # axis / legend label
    pattern: str  # regex matched against the sample name
    is_control: bool = False


WINDOWS: list[Window] = [
    Window("control", "Healthy control", r"control|healthy|innov|ctrl", True),
    Window("pre12", "15-9 months before Dx", r"12[-_ ]?m(onths?)?\b|12m"),
    Window("pre6", "9-3 months before Dx", r"(?<!1)6[-_ ]?m(onths?)?\b|(?<!1)6m"),
    Window("near", "Near Dx", r"at[-_ ]?cancer|near[-_ ]?dx|at[-_ ]?dx|\bdx\b"),
]

# Order used on every x-axis: earliest pre-diagnosis window first, control last
WINDOW_ORDER = ["pre12", "pre6", "near", "control"]
PRE_DX_KEYS = ["pre12", "pre6", "near"]


def window_of(sample_name: str) -> str | None:
    """Return the window key for a sample name, or None if no pattern matches."""
    for w in WINDOWS:
        if re.search(w.pattern, str(sample_name), flags=re.IGNORECASE):
            return w.key
    return None


def window_label(key: str) -> str:
    return next(w.label for w in WINDOWS if w.key == key)


def subject_of(sample_name: str, subject_ids: list[str] | None = None) -> str:
    """Return the subject ID a sample belongs to.

    If a subject ID list is given, the first ID contained in the sample name is
    returned. Otherwise the sample name up to the first separator is used.
    """
    name = str(sample_name)
    if subject_ids:
        for sid in subject_ids:
            if sid and sid in name:
                return sid
    return re.split(r"[_\- ]", name, maxsplit=1)[0]


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------
def ensure_dirs(*dirs: str) -> None:
    for d in dirs:
        os.makedirs(d, exist_ok=True)


def find_file(directory: str, keyword: str, exts=(".xlsx", ".csv")) -> str:
    """Find exactly one file whose name contains `keyword` (case-insensitive)."""
    hits = [
        f for f in sorted(glob.glob(os.path.join(directory, "*")))
        if keyword.lower() in os.path.basename(f).lower() and f.endswith(exts)
        and not os.path.basename(f).startswith("~$")
    ]
    if not hits:
        raise FileNotFoundError(f"No file containing '{keyword}' in '{directory}/'")
    if len(hits) > 1:
        print(f"[note] {len(hits)} files match '{keyword}'; using {os.path.basename(hits[0])}")
    return hits[0]


def read_table(path: str, **kw) -> pd.DataFrame:
    return pd.read_excel(path, **kw) if path.endswith((".xlsx", ".xls")) else pd.read_csv(path, **kw)


def read_id_list(path: str) -> list[str]:
    df = read_table(path)
    return df.iloc[:, 0].dropna().astype(str).str.strip().tolist()


# ---------------------------------------------------------------------------
# Count files from the lab preprocessing
# ---------------------------------------------------------------------------
def load_count_file(path: str) -> tuple[pd.DataFrame, pd.Series | None]:
    """Load a CGI fragment-count spreadsheet.

    Expected layout (the lab preprocessing output): first column holds row
    labels; rows whose label starts with ``CGI_`` are CpG islands; one summary
    row whose label contains "Total CpG island fragments counts" holds the
    per-sample total. Any other rows (pool, cancer type, notes) are ignored.

    Returns (cgi_counts, totals). ``totals`` is None if no summary row exists.
    """
    raw = read_table(path)
    label_col = raw.columns[0]
    labels = raw[label_col].astype(str)

    cgi = raw[labels.str.startswith("CGI_")].set_index(label_col)
    cgi.index.name = "CGI"
    cgi = cgi.apply(pd.to_numeric, errors="coerce").fillna(0)

    total_rows = raw[labels.str.lower().str.contains(TOTAL_ROW_KEY)]
    totals = None
    if not total_rows.empty:
        totals = pd.to_numeric(total_rows.iloc[0, 1:], errors="coerce")
        totals.index = raw.columns[1:]
        totals = totals.reindex(cgi.columns)
    return cgi, totals


def write_count_file(cgi: pd.DataFrame, totals: pd.Series | None, path: str) -> None:
    """Write counts back in the same layout load_count_file reads."""
    out = cgi.copy()
    if totals is not None:
        top = pd.DataFrame([totals.reindex(out.columns).values], columns=out.columns,
                           index=["Total CpG island fragments counts for this particular spreadsheet"])
        out = pd.concat([top, out])
    out.index.name = "Header"
    out.reset_index().to_excel(path, index=False)


# ---------------------------------------------------------------------------
# CGI names
# ---------------------------------------------------------------------------
def parse_cgi_name(name: str) -> dict:
    """Parse ``CGI_chr1_778604_779167_GENE1_GENE2_0`` into its parts.

    The trailing integer is the CGI index; anything between the coordinates
    and the index is a gene symbol (zero or more).
    """
    parts = str(name).split("_")
    if parts and parts[0] == "CGI":
        parts = parts[1:]
    index = parts.pop() if parts and parts[-1].isdigit() else ""
    chrom, start, end = (parts + ["", "", ""])[:3]
    genes = [g for g in parts[3:] if g]
    return {
        "CGI": name,
        "chr": chrom,
        "start": int(start) if str(start).isdigit() else np.nan,
        "end": int(end) if str(end).isdigit() else np.nan,
        "genes": genes,
        "cgi_index": index,
    }


# ---------------------------------------------------------------------------
# Plot style shared by the analysis scripts
# ---------------------------------------------------------------------------
WINDOW_COLORS = {
    "pre12": "#9ecae1",
    "pre6": "#4292c6",
    "near": "#08306b",
    "control": "#969696",
}


FIG_FORMAT = os.environ.get("FIG_FORMAT", "png")  # e.g. FIG_FORMAT=svg


def save_fig(fig, name: str) -> str:
    import matplotlib
    matplotlib.rcParams["svg.fonttype"] = "none"  # keep SVG text as text, not outlines
    matplotlib.rcParams["svg.hashsalt"] = "ed"     # stable element ids
    ensure_dirs(PLOTS_DIR)
    path = os.path.join(PLOTS_DIR, os.path.splitext(name)[0] + "." + FIG_FORMAT)
    fig.savefig(path, dpi=200, bbox_inches="tight", metadata={"Date": None} if FIG_FORMAT == "svg" else None)
    print(f"  saved {path}")
    return path
