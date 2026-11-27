#!/usr/bin/env bash
# Run every recommended script, in order, on the synthetic dataset.
# Works in a temporary copy, so your data/, output/ and plots/ are untouched.
#
# Usage (from the repository root):  bash tests/run_pipeline.sh [workdir]
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
W="${1:-$(mktemp -d)}"
mkdir -p "$W" && cd "$W"
rm -rf scripts data output plots
cp -r "$REPO/scripts" .
python "$REPO/tests/make_synthetic_data.py" data > /dev/null
export MPLBACKEND=Agg

STEPS=(
  scripts/step_1_filter_patients_local.py
  scripts/step_2_merge_filtered_files.py
  scripts/step_3_convert_to_aberrant_signals_v3-matchingRowColumnNames.py
  scripts/step_4_generate_gene_annotation.py
  scripts/step_5_build_gene_cgi_map.py
  scripts/step_6_gene-methylation-matrix.py
  scripts/global/dotplots-by-condition.py
  scripts/global/lineplots-perpatient_v3.py
  scripts/global/lineplots-perpatient_v4.py
  scripts/locus/avg-methylation-change-per-chromosome.py
  scripts/locus/bubbleplot_generator_v9_gridsoff.py
  scripts/locus/top10dm-plots_using-map_v4.py
  scripts/locus/gene_rank_slopeplot_top10_colored_v2.py
  scripts/locus/top10genes-heatmap-barplot.py
  scripts/locus/top10genes-barplot-heatmap-lineplot.py
)
failed=0
for s in "${STEPS[@]}"; do
  log="log_$(basename "$s" .py).txt"
  if python "$s" > "$log" 2>&1; then
    echo "ok    $s"
  else
    echo "FAIL  $s  (see $W/$log)"; tail -3 "$log"; failed=1
  fi
done
echo "Outputs in $W/output and $W/plots"
exit $failed
