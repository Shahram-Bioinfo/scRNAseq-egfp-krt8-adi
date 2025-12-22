from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import cosine as cosine_distance
from scipy.stats import pearsonr, spearmanr

import matplotlib.pyplot as plt
import seaborn as sns


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
PROCESSED_DIR = DATA_DIR / "processed"
OUT_DIR = PROJECT_ROOT / "outputs" / "degs_similarity"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TULANE_H5AD = PROCESSED_DIR / "adata_tulane_labeled.h5ad"

GFP_COL_CANDIDATES = ["GFP_Status", "EGFP_Status"]
ADI_COL_CANDIDATES = ["Krt8_ADI_like", "Krt8ADI_highlight", "Krt8ADI_label", "Krt8+ADI-like"]

TOPK = 50


def _to_dense(x):
    if hasattr(x, "A"):
        return x.A
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x)


def _find_obs_col(adata: sc.AnnData, candidates: list[str]) -> str:
    for c in candidates:
        if c in adata.obs.columns:
            return c
    raise KeyError(f"None of the expected columns found in adata.obs: {candidates}")


def _ensure_binary_labels(series: pd.Series, positive_values: set[str], out_pos="True", out_neg="False") -> pd.Categorical:
    s = series.astype(str)
    is_pos = s.isin(positive_values)
    return pd.Categorical(np.where(is_pos, out_pos, out_neg), categories=[out_neg, out_pos], ordered=True)


def run_deg_topk(adata: sc.AnnData, groupby: str, group: str, reference: str, topk: int) -> pd.DataFrame:
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=[group],
        reference=reference,
        method="wilcoxon",
    )
    df = sc.get.rank_genes_groups_df(adata, group=group)

    if "pvals_adj" in df.columns:
        df = df[df["pvals_adj"] < 0.05].copy()

    sort_col = "pvals_adj" if "pvals_adj" in df.columns else "pvals"
    return df.sort_values(sort_col, ascending=True).head(topk).copy()


def mean_expression_vector(adata: sc.AnnData, genes: list[str], mask: np.ndarray) -> np.ndarray:
    X = _to_dense(adata[:, genes].X)
    return X[mask].mean(axis=0)


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    return float(1.0 - cosine_distance(a, b))


def main() -> None:
    if not TULANE_H5AD.exists():
        raise FileNotFoundError(f"Missing input: {TULANE_H5AD}")

    print(f"[INFO] Reading AnnData: {TULANE_H5AD}")
    ad = sc.read_h5ad(TULANE_H5AD)
    ad.var_names_make_unique()

    gfp_col = _find_obs_col(ad, GFP_COL_CANDIDATES)
    adi_col = _find_obs_col(ad, ADI_COL_CANDIDATES)

    ad.obs["__GFP_BOOL__"] = _ensure_binary_labels(
        ad.obs[gfp_col],
        positive_values={"GFP+", "EGFP+", "True", "1"},
    )

    ad.obs["__ADI_BOOL__"] = _ensure_binary_labels(
        ad.obs[adi_col],
        positive_values={"Krt8+ ADI-like", "Krt8+ ADI", "Krt8+ADI-like", "True", "1"},
    )

    print("[INFO] Differential expression: EGFP+ vs EGFP-")
    deg_gfp = run_deg_topk(ad, groupby="__GFP_BOOL__", group="True", reference="False", topk=TOPK)
    deg_gfp.to_csv(OUT_DIR / "deg_EGFPplus_top50.csv", index=False)

    print("[INFO] Differential expression: Krt8+ ADI-like vs Other")
    deg_adi = run_deg_topk(ad, groupby="__ADI_BOOL__", group="True", reference="False", topk=TOPK)
    deg_adi.to_csv(OUT_DIR / "deg_Krt8ADI_top50.csv", index=False)

    genes_gfp = deg_gfp["names"].dropna().astype(str).tolist() if "names" in deg_gfp.columns else []
    genes_adi = deg_adi["names"].dropna().astype(str).tolist() if "names" in deg_adi.columns else []

    gene_union = sorted(set(genes_gfp) | set(genes_adi))
    gene_union = [g for g in gene_union if g in ad.var_names]
    if len(gene_union) == 0:
        raise RuntimeError("No DE genes available for similarity computation after filtering.")

    gfp_mask = (ad.obs["__GFP_BOOL__"].astype(str).values == "True")
    adi_mask = (ad.obs["__ADI_BOOL__"].astype(str).values == "True")

    gfp_avg = mean_expression_vector(ad, gene_union, gfp_mask)
    adi_avg = mean_expression_vector(ad, gene_union, adi_mask)

    metrics = {
        "n_union_genes_used": int(len(gene_union)),
        "pearson": float(pearsonr(gfp_avg, adi_avg)[0]),
        "spearman": float(spearmanr(gfp_avg, adi_avg)[0]),
        "cosine": float(cosine_similarity(gfp_avg, adi_avg)),
    }
    with open(OUT_DIR / "similarity_metrics.json", "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)

    df_means = pd.DataFrame(
        {"gene": gene_union, "EGFP_plus_mean": gfp_avg, "Krt8_ADI_like_mean": adi_avg}
    )
    df_means.to_csv(OUT_DIR / "shared_genes_mean_expression.csv", index=False)

    hm = pd.DataFrame(
        np.vstack([gfp_avg, adi_avg]),
        index=["EGFP+", "Krt8+ ADI-like"],
        columns=gene_union,
    )

    plt.figure(figsize=(max(10, 0.25 * len(gene_union)), 3.2))
    sns.heatmap(hm, cmap="viridis", cbar=True)
    plt.title("Mean expression (Top50 DEG union)")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "heatmap_mean_expression.png", dpi=300)
    plt.close()

    print(f"[INFO] Outputs saved to: {OUT_DIR}")


if __name__ == "__main__":
    main()
