#!/usr/bin/env bash
set -euo pipefail

# run_all.sh
# Execute the full scRNA-seq workflow in the intended order.
# Must be run from the repository root.

echo "[INFO] Starting scRNA-seq EGFP / Krt8 ADI analysis pipeline"

echo "[STEP 1] QC, normalization, and initial labeling"
python scripts/01_build_and_label_tulane.py

echo "[STEP 2] EGFP+ and Krt8+ ADI-like cell annotation"
python scripts/02_label_egfp_adi_tulane.py

echo "[STEP 3] Leiden clustering (resolution = 0.17) and UMAP visualization"
python scripts/03_leiden_umap_visualization.py

echo "[STEP 4] Differential expression and transcriptomic similarity analysis"
python scripts/04_degs_and_similarity.py

echo "[INFO] Pipeline completed successfully"

