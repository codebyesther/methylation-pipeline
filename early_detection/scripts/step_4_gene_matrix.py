"""Step 4: annotate CGIs with genes and sum CGI signals per gene.

CGI names carry their genes: CGI_chr1_778604_779167_LOC100288069_0.
This step parses every name in output/scaled_signal_matrix.csv and writes:

  cgi_annotation.csv       one row per CGI: chr, start, end, genes, index
  cgi_gene_map.csv         one row per CGI-gene pair
  gene_signal_matrix.csv   gene x sample, sum of the CGI signals

--gene-mode first  (default, as in the 2024 notebooks)
    each CGI counts toward its first listed gene only, so every CGI is
    counted once.
--gene-mode all
    each CGI counts toward every gene it lists. Signal of a CGI shared by two
    genes is counted in both.

Genes are matched by exact symbol, never by substring, so GATA2 does not
collect GATA2-AS1's CGIs. CGIs without a gene symbol are left out.
"""

import argparse
import os

import pandas as pd

from edcommon import OUTPUT_DIR, parse_cgi_name


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=OUTPUT_DIR)
    ap.add_argument("--gene-mode", choices=["first", "all"], default="first")
    args = ap.parse_args()

    signal = pd.read_csv(os.path.join(args.out, "scaled_signal_matrix.csv"), index_col=0)

    ann = pd.DataFrame([parse_cgi_name(c) for c in signal.index])
    ann["n_genes"] = ann["genes"].str.len()
    ann_out = ann.assign(genes=ann["genes"].str.join(";"))
    ann_out.to_csv(os.path.join(args.out, "cgi_annotation.csv"), index=False)

    pairs = ann[["CGI", "genes"]].explode("genes").dropna(subset=["genes"])
    if args.gene_mode == "first":
        pairs = pairs.groupby("CGI", sort=False).head(1)
    pairs = pairs.rename(columns={"genes": "gene"})
    pairs.to_csv(os.path.join(args.out, "cgi_gene_map.csv"), index=False)

    gene = (signal.loc[pairs["CGI"]]
            .assign(gene=pairs["gene"].values)
            .groupby("gene").sum(min_count=1))
    gene.index.name = "gene"
    gene.to_csv(os.path.join(args.out, "gene_signal_matrix.csv"))

    no_gene = int((ann["n_genes"] == 0).sum())
    print(f"{len(ann)} CGIs parsed; {no_gene} without a gene symbol (left out)")
    print(f"gene-mode={args.gene_mode}: {gene.shape[0]} genes x {gene.shape[1]} samples")


if __name__ == "__main__":
    main()
