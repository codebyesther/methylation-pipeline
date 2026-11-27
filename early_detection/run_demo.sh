#!/usr/bin/env bash
# Run the whole early-detection pipeline on the synthetic demo data.
# Usage (from this folder):  bash run_demo.sh
set -euo pipefail
cd "$(dirname "$0")"

python demo/make_synthetic_data.py
rm -rf data output plots
mkdir -p data
cp demo/data/*.xlsx data/

python scripts/step_1_filter_samples.py
python scripts/step_2_merge_runs.py
python scripts/step_3_scaled_signal.py --mode total
python scripts/step_4_gene_matrix.py --gene-mode first
python scripts/step_5_window_summary.py

python scripts/analysis/global_signal.py
python scripts/analysis/top_genes.py --top 30 --rank-window pre12
python scripts/analysis/top_genes.py --top 15 --rank-window pre12 --subject HL02
python scripts/analysis/control_contrast.py

echo "Done. Tables in output/, figures in plots/."
