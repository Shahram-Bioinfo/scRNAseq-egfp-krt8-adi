#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

IN_H5AD = PROCESSED_DIR / "adata_tulane_normalized.h5ad"
OUT_H5AD = PROCESSED_DIR / "adata_tulane_normalized_labeled.h5ad"
OUT_CROSSTAB = PROCESSED_DIR / "Tulane_Krt8ADI_vs_GFP_crosstab.csv"

KRT8_MARKERS = [
    "Krt8","Krt18","Krt7","Krt19",
    "Sprr1a","Areg","Hbegf","Tnfrsf12a","Itgb6",
    "Cldn18","Cldn4","Cryab","Cdkn1a"
]

THRESH = 1.0

def to_dense(x):
    if hasattr(x, "toarray"):
        return x.toarray()
    if hasattr(x, "A"):
        return x.A
    return np.asarray(x)

def main():
    adata = sc.read_h5ad(str(IN_H5AD))

    present = [g for g in KRT8_MARKERS if g in adata.var_names]
    if len(present) < 5:
        raise RuntimeError(f"Too few ADI markers present ({len(present)}).")

    X = to_dense(adata[:, present].X)

    is_adi = np.any(X > THRESH, axis=1)

    adata.obs["Krt8+ADI-like"] = pd.Categorical(
        np.where(is_adi, "Krt8+ADI-like", "Other")
    )

    adata.write(str(OUT_H5AD), compression="gzip")

    pd.crosstab(
        adata.obs["Krt8+ADI-like"],
        adata.obs["GFP_Status"]
    ).to_csv(str(OUT_CROSSTAB))

    print(str(OUT_H5AD))

if __name__ == "__main__":
    main()

