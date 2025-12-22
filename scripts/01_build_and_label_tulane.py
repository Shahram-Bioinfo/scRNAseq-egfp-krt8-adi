from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
INPUT_DIR = DATA_DIR / "input"
OUTPUT_DIR = DATA_DIR / "processed"

INPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MTX_FILE = INPUT_DIR / "matrix.mtx"
FEATURES_FILE = INPUT_DIR / "features.tsv"
BARCODES_FILE = INPUT_DIR / "barcodes.tsv"

# Output is named to be consumed by scripts/02_label_egfp_adi_tulane.py
OUT_H5AD = OUTPUT_DIR / "adata_tulane_clean.h5ad"

ADI_MARKERS = [
    "Krt8", "Krt18", "Krt7", "Krt19",
    "Sprr1a", "Areg", "Hbegf", "Tnfrsf12a",
    "Itgb6", "Cldn18", "Cldn4", "Cryab", "Cdkn1a"
]


def main() -> None:
    print("[INFO] Loading scRNA-seq matrix input files")

    adata = sc.read_mtx(MTX_FILE).T
    features = pd.read_csv(FEATURES_FILE, header=None, sep="\t")
    barcodes = pd.read_csv(BARCODES_FILE, header=None, sep="\t")

    adata.var_names = features[1].values
    adata.obs_names = barcodes[0].values
    adata.var["gene_ids"] = features[0].values
    adata.var_names_make_unique()

    print("[INFO] Running basic QC filters")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    print("[INFO] Normalizing to 1e4 counts per cell and log-transforming")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print("[INFO] Labeling EGFP status (EGFP > 0)")
    if "EGFP" not in adata.var_names:
        raise RuntimeError("EGFP gene not found in dataset.")
    egfp_expr = adata[:, "EGFP"].X.toarray().ravel()
    adata.obs["EGFP_Status"] = np.where(egfp_expr > 0, "EGFP+", "EGFP-").astype(str)
    adata.obs["EGFP_Status"] = adata.obs["EGFP_Status"].astype("category")

    print("[INFO] Labeling Krt8+ ADI-like cells (marker expression > 1)")
    available_markers = [g for g in ADI_MARKERS if g in adata.var_names]
    if len(available_markers) == 0:
        raise RuntimeError("None of the ADI marker genes were found in the dataset.")
    X_markers = adata[:, available_markers].X.toarray()
    adata.obs["Krt8_ADI_like"] = np.where(
        np.any(X_markers > 1, axis=1),
        "Krt8+ ADI-like",
        "Other"
    ).astype(str)
    adata.obs["Krt8_ADI_like"] = adata.obs["Krt8_ADI_like"].astype("category")

    adata.write(OUT_H5AD)
    print(f"[INFO] Saved processed AnnData: {OUT_H5AD}")


if __name__ == "__main__":
    main()
