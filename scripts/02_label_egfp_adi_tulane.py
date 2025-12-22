from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
PROCESSED_DIR = DATA_DIR / "processed"

TULANE_H5AD = PROCESSED_DIR / "adata_tulane_clean.h5ad"
OUT_PATH = PROCESSED_DIR / "adata_tulane_labeled.h5ad"

GFP_GENE_CANDIDATES = ["EGFP", "Egfp", "eGFP", "GFP"]
EPITHELIAL_MARKERS = ["Epcam", "Krt8", "Krt18", "Krt19", "Krt7"]

ADI_SIGNATURE = [
    "Krt8", "Krt18", "Krt7", "Krt19",
    "Cldn4", "Areg", "Hbegf", "Edn1", "Cyr61", "Plaur", "Tnip3",
    "Cdkn1a", "S100a6", "S100a10", "Lgals3", "Tnfrsf12a", "Fn1", "Ctgf", "Itgb6"
]


def present(varnames, genes):
    s = set(varnames)
    return [g for g in genes if g in s]


def to_dense(X):
    if hasattr(X, "A"):
        return X.A
    if hasattr(X, "toarray"):
        return X.toarray()
    return np.asarray(X)


def main() -> None:
    print(f"[INFO] Reading AnnData: {TULANE_H5AD}")
    ad = sc.read_h5ad(TULANE_H5AD)
    ad.var_names_make_unique()

    egfp = next((g for g in GFP_GENE_CANDIDATES if g in ad.var_names), None)
    if egfp is None:
        raise RuntimeError(f"EGFP/GFP gene not found. Checked: {GFP_GENE_CANDIDATES}")

    print("[INFO] Labeling GFP status (EGFP > 0)")
    eg = to_dense(ad[:, egfp].X).ravel().astype(float)
    ad.obs["GFP_Status"] = pd.Categorical(
        np.where(eg > 0, "GFP+", "GFP-"),
        categories=["GFP-", "GFP+"],
        ordered=True,
    )

    epi = present(ad.var_names, EPITHELIAL_MARKERS)
    if len(epi) >= 2:
        Xe = to_dense(ad[:, epi].X)
        epi_mask = (Xe > 0).sum(axis=1) >= 2
    else:
        epi_mask = np.ones(ad.n_obs, dtype=bool)

    print("[INFO] Scoring ADI signature (score_genes)")
    sig = present(ad.var_names, ADI_SIGNATURE)
    if len(sig) < 5:
        print(f"[WARN] Only {len(sig)} ADI signature genes present; scoring anyway.")

    ad_epi = ad[epi_mask].copy()
    sc.tl.score_genes(ad_epi, gene_list=sig, score_name="ADI_score")

    ad.obs["ADI_score"] = np.nan
    ad.obs.loc[ad_epi.obs_names, "ADI_score"] = ad_epi.obs["ADI_score"].values

    print("[INFO] Defining Krt8+ ADI-like cells (ADI_score >= 90th percentile)")
    thr = np.nanpercentile(ad.obs["ADI_score"].values, 90)
    ad.obs["Krt8_ADI_like"] = np.where(
        ad.obs["ADI_score"].values >= thr,
        "Krt8+ ADI-like",
        "Other",
    )
    ad.obs["Krt8_ADI_like"] = pd.Categorical(
        ad.obs["Krt8_ADI_like"],
        categories=["Other", "Krt8+ ADI-like"],
        ordered=True,
    )

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    ad.write(OUT_PATH)
    print(f"[INFO] Saved labeled AnnData: {OUT_PATH}")

    crosstab_path = OUT_PATH.with_name("Tulane_Krt8ADI_vs_GFP_crosstab.csv")
    pd.crosstab(ad.obs["Krt8_ADI_like"], ad.obs["GFP_Status"]).to_csv(crosstab_path)
    print(f"[INFO] Saved crosstab: {crosstab_path}")


if __name__ == "__main__":
    main()
