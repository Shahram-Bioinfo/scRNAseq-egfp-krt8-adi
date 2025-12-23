#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = PROJECT_ROOT / "data" / "input"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

OUT_H5AD = PROCESSED_DIR / "adata_tulane_normalized.h5ad"

MIN_GENES = 200
MIN_CELLS = 3
TARGET_SUM = 1e4
MT_PCT_MAX = 20.0
EGFP_CANDIDATES = ["EGFP", "Egfp", "eGFP", "GFP"]

def _to_dense_1d(x):
    if hasattr(x, "toarray"):
        return x.toarray().ravel()
    if hasattr(x, "A"):
        return x.A.ravel()
    return np.asarray(x).ravel()

def _find_first(patterns, root: Path):
    for pat in patterns:
        hits = list(root.rglob(pat))
        if hits:
            return hits[0]
    return None

def _read_10x_h5(h5_path: Path):
    adata = sc.read_10x_h5(str(h5_path))
    adata.var_names_make_unique()
    return adata

def _read_mtx_triplet(mtx_path: Path, features_path: Path, barcodes_path: Path):
    adata = sc.read_mtx(str(mtx_path)).T
    features = pd.read_csv(str(features_path), header=None, sep="\t", compression="infer")
    barcodes = pd.read_csv(str(barcodes_path), header=None, sep="\t", compression="infer")
    adata.var_names = features[1].astype(str).values
    adata.obs_names = barcodes[0].astype(str).values
    adata.var["gene_ids"] = features[0].astype(str).values
    adata.var_names_make_unique()
    return adata

def main():
    if not INPUT_DIR.exists():
        raise FileNotFoundError(str(INPUT_DIR))

    h5 = _find_first(["filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5"], INPUT_DIR)
    if h5 is not None:
        adata = _read_10x_h5(h5)
    else:
        mtx = _find_first(["matrix.mtx", "matrix.mtx.gz"], INPUT_DIR)
        features = _find_first(["features.tsv", "features.tsv.gz"], INPUT_DIR)
        barcodes = _find_first(["barcodes.tsv", "barcodes.tsv.gz"], INPUT_DIR)
        if mtx is None or features is None or barcodes is None:
            raise FileNotFoundError(
                "Expected either a 10x .h5 file or matrix.mtx(.gz)+features.tsv(.gz)+barcodes.tsv(.gz) under data/input."
            )
        adata = _read_mtx_triplet(mtx, features, barcodes)

    adata.obs["CellType"] = "Unknown"

    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)

    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs["pct_counts_mt"] < MT_PCT_MAX].copy()

    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)

    egfp_gene = next((g for g in EGFP_CANDIDATES if g in adata.var_names), None)
    if egfp_gene is None:
        raise RuntimeError("EGFP gene not found in var_names; cannot derive GFP_Status.")
    xg = _to_dense_1d(adata[:, egfp_gene].X)
    adata.obs["GFP_Status"] = pd.Categorical(np.where(xg > 0, "GFP+", "GFP-"))

    adata.write(str(OUT_H5AD), compression="gzip")
    print(str(OUT_H5AD))

if __name__ == "__main__":
    main()


