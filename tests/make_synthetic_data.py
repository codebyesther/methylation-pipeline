"""Make a small SYNTHETIC treatment-monitoring dataset in the lab file layout.

No real samples: 6 patients (Baseline, C2D1, C4D1, C7D1-Off-tx; one replicate),
4 healthy controls (INNOV), two sequencing runs, Glob20 and GlobMin80 files,
some zero GlobMin80 cells, chrX/chrY islands, and gene pairs such as GATA2 /
GATA2-AS1 that a substring match would confuse.

Usage (from the repository root):  python tests/make_synthetic_data.py data
"""
import os
import sys

import numpy as np
import pandas as pd

out = sys.argv[1] if len(sys.argv) > 1 else "data"
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(7)

GENES = ["GATA2", "GATA2-AS1", "APC", "APCDD1", "SEPT9", "SHOX2", "RASSF1", "SOX17", "LYPD1",
         "NKX2-2", "PLD6", "ZBTB16", "TWIST1", "VIM", "SFRP1", "HOXA9", "PAX6", "OTX2", "TBX5",
         "EN2", "HAND2", "ALX4", "GDNF", "RUNX3", "IKZF1", "BCAT1", "SDC2", "NDRG4", "TAC1", "NPY"]
chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
cgis = []
for i in range(300):
    c = chroms[i % len(chroms)]
    s = int(rng.integers(1_000_000, 100_000_000))
    k = rng.choice([0, 1, 1, 2])
    g = list(rng.choice(GENES, size=k, replace=False)) if k else []
    cgis.append("_".join(["CGI", c, str(s), str(s + int(rng.integers(300, 1500))), *g, str(i % 3)]))
cgis = list(dict.fromkeys(cgis))

patients = ["PT001", "PT002", "PT003", "PT010", "PT011", "PT012"]   # PT001 vs PT010-style prefixes
tps = ["Baseline", "C2D1", "C4D1", "C7D1-Off-tx"]
response = {"Baseline": 1.0, "C2D1": 0.5, "C4D1": 0.7, "C7D1-Off-tx": 1.3}


def sample(name, factor):
    depth = rng.uniform(20, 60)
    ref = rng.poisson(depth * rng.gamma(2, 1, len(cgis)))
    ref[rng.random(len(cgis)) < 0.05] = 0            # some zero reference cells
    num = rng.poisson(depth * 0.05 * factor * rng.gamma(1.5, 1, len(cgis)))
    return name, num, ref


runs = {"EMSeq-16-Run-1": [], "EMSeq-18-Run-1": []}
for j, p in enumerate(patients):
    run = "EMSeq-16-Run-1" if j < 3 else "EMSeq-18-Run-1"
    for t in tps:
        runs[run].append(sample(f"{p}_{t}", response[t] * rng.uniform(0.5, 2)))
    if p == "PT002":
        runs[run].append(sample(f"{p}_C2D1_rep2", response["C2D1"]))
runs["EMSeq-16-Run-1"] += [sample(f"INNOV_control_0{i}", 0.6) for i in (1, 2)]
runs["EMSeq-18-Run-1"] += [sample(f"INNOV_control_0{i}", 0.6) for i in (3, 4)]
runs["EMSeq-18-Run-1"].append(sample("OTHER99_Baseline", 1.0))   # not in the patient list

for run, samples in runs.items():
    for kind, idx in (("Glob20", 1), ("GlobMin80", 2)):
        cols = {name: v[idx - 1] for name, *v in [(s[0], s[1], s[2]) for s in samples]}
        df = pd.DataFrame(cols, index=cgis)
        tot = pd.DataFrame([df.sum().values], columns=df.columns,
                           index=["Total CpG island fragments counts for this particular spreadsheet"])
        info = pd.DataFrame([["x"] * df.shape[1]], columns=df.columns, index=["Info_Pool"])
        full = pd.concat([info, tot, df])
        full.index.name = "Header"
        full.reset_index().to_excel(os.path.join(out, f"{run}_CG80_CH20_Tot10_{kind}_cpgi_counts.xlsx"),
                                    index=False, sheet_name="Sheet1")

ids = pd.DataFrame({"Patient ID": patients + ["INNOV"]})  # INNOV keeps the healthy controls
ids.to_excel(os.path.join(out, "Patient ID list.xlsx"), index=False)
print("wrote", out)
