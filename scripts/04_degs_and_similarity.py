#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import json
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import cosine
from scipy.stats import pearsonr, spearmanr

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
OUT_DIR = PROJECT_ROOT / "outputs" / "degs_similarity"
OUT_DIR.mkdir(parents=True, exist_ok=True)

IN_H5AD = PROCESSED_DIR / "adata_tulane_normalized_labeled.h5ad"

TOP_N = 50
ALPHA = 0.05

def _to_dense_2d(x):
    if hasattr(x, "toarray"):
        return x.toarray()
    if hasattr(x, "A"):
        return x.A
    return np.asarray(x)

def top_degs_to_df(adata, group_key, group_name, top_n):
    df = sc.get.rank_genes_groups_df(adata, group=group_name)
    df = df[df["pvals_adj"] < ALPHA].copy()
    df = df.sort_values(["scores"], ascending=False).head(top_n)
    df["group"] = group_name
    df["group_key"] = group_key
    return df

def main():
    adata = sc.read_h5ad(str(IN_H5AD))

    if "GFP_Status" not in adata.obs.columns:
        raise RuntimeError("GFP_Status missing.")
    if "Krt8+ADI-like" not in adata.obs.columns:
        raise RuntimeError("Krt8+ADI-like missing.")

    adata_egfp = adata.copy()
    adata_egfp.obs["deg_group"] = pd.Categorical(
        np.where(adata_egfp.obs["GFP_Status"].astype(str) == "GFP+", "EGFP+", "Other")
    )
    sc.tl.rank_genes_groups(adata_egfp, groupby="deg_group", method="wilcoxon")
    deg_egfp = top_degs_to_df(adata_egfp, "deg_group", "EGFP+", TOP_N)
    deg_egfp.to_csv(OUT_DIR / "deg_EGFPplus_top50.csv", index=False)

    adata_krt8 = adata.copy()
    adata_krt8.obs["deg_group"] = pd.Categorical(
        np.where(adata_krt8.obs["Krt8+ADI-like"].astype(str) == "Krt8+ADI-like", "Krt8+ADI-like", "Other")
    )
    sc.tl.rank_genes_groups(adata_krt8, groupby="deg_group", method="wilcoxon")
    deg_krt8 = top_degs_to_df(adata_krt8, "deg_group", "Krt8+ADI-like", TOP_N)
    deg_krt8.to_csv(OUT_DIR / "deg_Krt8ADI_top50.csv", index=False)

    genes_egfp = [g for g in deg_egfp["names"].tolist() if g in adata.var_names]
    genes_krt8 = [g for g in deg_krt8["names"].tolist() if g in adata.var_names]
    genes_union = sorted(set(genes_egfp) | set(genes_krt8))

    egfp_cells = adata.obs["GFP_Status"].astype(str) == "GFP+"
    krt8_cells = adata.obs["Krt8+ADI-like"].astype(str) == "Krt8+ADI-like"

    X_egfp = _to_dense_2d(adata[egfp_cells, genes_union].X)
    X_krt8 = _to_dense_2d(adata[krt8_cells, genes_union].X)

    v1 = np.asarray(X_egfp.mean(axis=0)).ravel()
    v2 = np.asarray(X_krt8.mean(axis=0)).ravel()

    pear = float(pearsonr(v1, v2)[0]) if len(v1) > 1 else float("nan")
    spear = float(spearmanr(v1, v2)[0]) if len(v1) > 1 else float("nan")
    cos_sim = float(1.0 - cosine(v1, v2)) if (np.any(v1) or np.any(v2)) else float("nan")

    metrics = {"n_union_genes": int(len(genes_union)), "pearson": pear, "spearman": spear, "cosine_similarity": cos_sim}

    with open(OUT_DIR / "similarity_metrics.json", "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)

    pd.DataFrame({"gene": genes_union, "mean_EGFPplus": v1, "mean_Krt8ADI": v2}).to_csv(
        OUT_DIR / "shared_genes_mean_expression.csv", index=False
    )

    print(str(OUT_DIR))

if __name__ == "__main__":
    main()
