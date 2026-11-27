"""Make a small SYNTHETIC dataset in the lab's spreadsheet layout.

No real samples are used. Gene symbols are real (public) so the CGI names
look realistic; counts, subjects and effects are simulated.

Design:
  * 6 cases, each sampled 15-9 months before diagnosis, 9-3 months before,
    and near diagnosis; 6 healthy controls sampled once
  * two sequencing runs (cases 1-3 + controls 1-3 in run A, the rest in
    run B) to exercise the merge step
  * 25 "informative" genes whose Glob20 CGIs gain fragments toward
    diagnosis in cases only; everything else is background
  * one extra sample not in the ID list, and one with no window label, to
    show what step 1 filters out

Writes into demo/data/:
  subject_ids.xlsx
  RunA_..._Glob20_cpgi_counts.xlsx, RunA_..._GlobMin80_cpgi_counts.xlsx
  RunB_..._Glob20_cpgi_counts.xlsx, RunB_..._GlobMin80_cpgi_counts.xlsx
"""

import os

import numpy as np
import pandas as pd

rng = np.random.default_rng(2024)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data")
os.makedirs(OUT, exist_ok=True)

GENES = """PRR25 LBX2 GLI4 GATA2 NDST3 FAM53A NECAB3 PDE2A PAX2 TMEM151A PLD6 RASIP1
NKX1-1 FENDRR TNPO2 OLIG2 SOX17 ZNF503 NTNG1 SIM2 ZBTB16 LYPD1 CBX8 NKX2-2 NAPSA
SEPT9 SHOX2 RASSF1 APC CDKN2A MGMT GSTP1 RARB TWIST1 VIM SFRP1 SFRP2 WIF1 DKK1
HOXA9 HOXD13 PCDH10 TFPI2 BCAT1 IKZF1 NDRG4 BMP3 SDC2 ZNF331 HIC1 PENK NPY TAC1
SOX1 PAX6 OTX2 LHX2 FOXA1 FOXE1 TBX5 IRX1 EN2 HAND2 ALX4 NKX6-1 GDNF RUNX3""".split()
INFORMATIVE = GENES[:25]

CASES = [f"HL{i:02d}" for i in range(1, 7)]
CONTROLS = [f"CTRL{i:02d}" for i in range(1, 7)]
WINDOWS = {"12m": 0.15, "6m": 0.5, "near-Dx": 1.0}  # effect size multiplier


def make_cgis(n, genes_pool, prefix_seed):
    r = np.random.default_rng(prefix_seed)
    names, gene_of = [], []
    for i in range(n):
        chrom = f"chr{r.integers(1, 23)}"
        start = int(r.integers(1_000_000, 200_000_000))
        end = start + int(r.integers(200, 2000))
        k = r.choice([0, 1, 1, 1, 2], p=[0.1, 0.3, 0.3, 0.2, 0.1])
        genes = list(r.choice(genes_pool, size=k, replace=False)) if k else []
        names.append("_".join(["CGI", chrom, str(start), str(end), *genes, str(i % 4)]))
        gene_of.append(genes)
    return names, gene_of


glob20_names, glob20_genes = make_cgis(400, GENES, 1)
globmin80_names, _ = make_cgis(500, GENES, 2)
informative_mask = np.array([any(g in INFORMATIVE for g in gs) for gs in glob20_genes])


def sample_counts(subject, window, depth):
    """Return (glob20 counts, globmin80 counts) for one sample."""
    ref = rng.poisson(depth * rng.gamma(2.0, 1.0, len(globmin80_names)))
    person = rng.gamma(6.0, 1 / 6.0)  # person-to-person background, sd ~0.4
    lam = depth * 0.02 * person * rng.gamma(1.5, 1.0, len(glob20_names))
    if subject in CASES:
        effect = WINDOWS[window] * rng.uniform(0.6, 1.4)
        lam = lam + informative_mask * depth * 0.06 * effect
    return rng.poisson(lam), ref


def build_run(cases, controls, extra=()):
    g20, gm80 = {}, {}
    for s in cases:
        for w in WINDOWS:
            name = f"{s}_{w}"
            a, b = sample_counts(s, w, depth=rng.uniform(15, 60))
            g20[name], gm80[name] = a, b
    for s in controls:
        name = f"{s}_control"
        a, b = sample_counts(s, None, depth=rng.uniform(15, 60))
        g20[name], gm80[name] = a, b
    for name in extra:
        a, b = sample_counts("x", None, depth=30)
        g20[name], gm80[name] = a, b
    return pd.DataFrame(g20, index=glob20_names), pd.DataFrame(gm80, index=globmin80_names)


def write(df, path, pool):
    info = pd.DataFrame([[pool] * df.shape[1]], columns=df.columns, index=["Info_Pool"])
    total = pd.DataFrame([df.sum().values], columns=df.columns,
                         index=["Total CpG island fragments counts for this particular spreadsheet"])
    out = pd.concat([info, total, df])
    out.index.name = "Header"
    out.reset_index().to_excel(path, index=False)


runA = build_run(CASES[:3], CONTROLS[:3], extra=["OTHERPROJ07_baseline"])
runB = build_run(CASES[3:], CONTROLS[3:], extra=["HL06_unlabeled"])
for tag, (g20, gm80) in {"RunA": runA, "RunB": runB}.items():
    write(g20, os.path.join(OUT, f"{tag}_CG80_CH20_Tot10_Glob20_cpgi_counts.xlsx"), tag)
    write(gm80, os.path.join(OUT, f"{tag}_CG80_CH20_Tot10_GlobMin80_cpgi_counts.xlsx"), tag)

pd.DataFrame({"Subject ID": CASES + CONTROLS}).to_excel(os.path.join(OUT, "subject_ids.xlsx"), index=False)
print(f"Synthetic data written to {OUT}")
