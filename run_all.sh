#!/usr/bin/env bash
set -euo pipefail

# run_all.sh
# Execute the full scRNA-seq workflow in the intended order.
# Assumes you are running from the repository root.

echo "[INFO] Running scRNA-seq pipeline..."

python scripts/01_build_and_label_tulane.py
python scripts/02_label_egfp_adi_tulane.py
python scripts/03_leiden_umap_visualization.py
python scripts/04_degs_and_similarity.py

echo "[INFO] Done."

