# Legacy scripts

Earlier versions of pipeline scripts, kept unchanged so figures made with them can be reproduced. For new analyses use the scripts in `scripts/`, `scripts/global/` and `scripts/locus/`, which include the fixes listed in the main README's changelog.

Run them from the repository root, as before, e.g. `python scripts/legacy/locus/generank-slopeplot-21cap.py`.

| Script | Superseded by | Difference |
|:--|:--|:--|
| `step_3_convert_to_aberrant_signals.py` | `step_3_..._v3-matchingRowColumnNames.py` | x1,000 scaling, rows matched by position |
| `step_3_convert_to_aberrant_signals_v2.py` | `step_3_..._v3-matchingRowColumnNames.py` | x100,000 scaling, rows matched by label, columns by position |
| `locus/bubbleplot_generator_v8.py`, `_v8_gridsoff.py` | `locus/bubbleplot_generator_v9_gridsoff.py` | earlier layout |
| `locus/merged-glob20-bubbleplot.py` | `locus/bubbleplot_generator_v9_gridsoff.py` | plots raw Glob20 counts from `merged_output_glob20.xlsx` |
| `locus/top10dm-plots.py`, `_using-map_v2.py`, `_v3.py` | `locus/top10dm-plots_using-map_v4.py` | fewer intermediate tables |
| `locus/heatmap-lineplot-barplot_v2.py`, `_v3.py` | `locus/top10genes-heatmap-barplot.py` | earlier versions |
| `locus/generank-slopeplot.py` | `locus/gene_rank_slopeplot_top10_colored_v2.py` | no rank cap; reads `plots/heatmaps-lineplots/gene_methylation_matrix.csv` |
| `locus/generank-slopeplot-21cap.py`, `-501cap.py` | `locus/gene_rank_slopeplot_top10_colored_v2.py` | rank cap 21 or 501, no top-10 coloring |
| `locus/gene_rank_slopeplot_top10_colored.py` | `locus/gene_rank_slopeplot_top10_colored_v2.py` | replicates not averaged |

Known issues in these versions (fixed in the current scripts): patient IDs and gene symbols are matched as substrings (e.g. `GATA2` also collects `GATA2-AS1`), the patient-list file name is fixed in the gene-rank scripts, and the pandas `groupby(axis=1)` used by the bubble-plot scripts no longer exists in pandas 3.
