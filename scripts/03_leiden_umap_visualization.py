from __future__ import annotations

from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


PROJECT_ROOT = Path(__file__).resolve().parents[1]
IN_H5AD = PROJECT_ROOT / "data" / "processed" / "adata_tulane_labeled.h5ad"

STAMP = datetime.now().strftime("%Y%m%d_%H%M")
OUT_DIR = PROJECT_ROOT / "outputs" / "umap_leiden_res017" / STAMP
OUT_DIR.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = str(OUT_DIR)
sc.settings.autoshow = False
sc.settings.set_figure_params(dpi=220, dpi_save=220, facecolor="white")

LEIDEN_RES = 0.17
RANDOM_STATE = 0

CURATED_MARKERS = {
    "AT1": ["Pdpn", "Ager", "Cav1"],
    "AT2": ["Sftpc", "Sftpb", "Abca3", "Lamp3"],
    "Injury_AT2": ["S100a10", "Lcn2", "Lgals3"],
    "Krt8_ADI": [
        "Krt8", "Krt18", "Krt7", "Krt19", "Sprr1a", "Areg", "Hbegf", "Tnfrsf12a", "Itgb6",
        "Cldn18", "Cldn4", "Cryab", "Cdkn1a", "S100a6", "S100a11", "Emp2", "Icam1", "Gprc5a", "Nupr1", "Sfn",
    ],
    "Basal": ["Krt5", "Krt14", "Trp63"],
    "Club": ["Scgb1a1", "Scgb3a2"],
    "Ciliated": ["Foxj1", "Pifo", "Tekt1"],
    "Endothelial": ["Pecam1", "Kdr"],
    "Fibroblast": ["Col1a1", "Col1a2", "Dcn"],
    "Macrophage": ["Lyz2", "Adgre1"],
    "T_cell": ["Cd3d", "Cd3e", "Trac"],
    "B_cell": ["Ms4a1", "Cd79a"],
    "NK": ["Nkg7", "Klrb1c"],
    "DC": ["Itgax", "H2-Aa"],
}


def to_numpy(X):
    if hasattr(X, "A"):
        return X.A
    if hasattr(X, "toarray"):
        return X.toarray()
    return np.asarray(X)


def present_genes(var_names, genes):
    s = set(var_names)
    return [g for g in genes if g in s]


def safe_savefig(fig, path, **kwargs):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(p), **kwargs)


def leiden_with_fallback(adata, resolution, key_added="leiden"):
    try:
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=key_added,
            random_state=RANDOM_STATE,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
    except TypeError:
        sc.tl.leiden(adata, resolution=resolution, key_added=key_added, random_state=RANDOM_STATE)


def main() -> None:
    if not IN_H5AD.exists():
        raise FileNotFoundError(f"Input file not found: {IN_H5AD}")

    adata = sc.read_h5ad(IN_H5AD)

    if "GFP_Status" not in adata.obs.columns:
        egfp_gene = next((g for g in ["EGFP", "Egfp", "eGFP", "GFP"] if g in adata.var_names), None)
        if egfp_gene is None:
            raise RuntimeError("Missing obs['GFP_Status'] and no EGFP gene found to derive it.")
        eg = to_numpy(adata[:, egfp_gene].X).ravel().astype(float)
        adata.obs["GFP_Status"] = pd.Categorical(np.where(eg > 0, "GFP+", "GFP-"))

    krt8_label_candidates = ["Krt8_ADI_like", "Krt8ADI_highlight", "Krt8+ADI-like", "Krt8+ ADI", "Krt8_ADI_bool", "Krt8_ADI"]
    krt8_col = next((c for c in krt8_label_candidates if c in adata.obs.columns), None)
    if krt8_col is None:
        raise RuntimeError("No Krt8+ ADI label found in adata.obs.")

    sc.pp.pca(adata, n_comps=30, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(30, adata.obsm["X_pca"].shape[1]))
    sc.tl.umap(adata, random_state=RANDOM_STATE)
    leiden_with_fallback(adata, resolution=LEIDEN_RES, key_added="leiden")

    gfp_mask = (adata.obs["GFP_Status"] == "GFP+").to_numpy()
    krt8_mask = adata.obs[krt8_col].astype(str).isin(["Krt8+ ADI-like", "Krt8+ ADI", "True", "1"]).to_numpy()

    adata.obs["EGFP_highlight"] = pd.Categorical(
        np.where(gfp_mask, "GFP+", "Other"),
        categories=["Other", "GFP+"],
        ordered=True,
    )
    adata.obs["Krt8ADI_highlight"] = pd.Categorical(
        np.where(krt8_mask, "Krt8+ ADI", "Other"),
        categories=["Other", "Krt8+ ADI"],
        ordered=True,
    )
    adata.uns["EGFP_highlight_colors"] = ["lightgrey", "#1f77b4"]
    adata.uns["Krt8ADI_highlight_colors"] = ["lightgrey", "#1f77b4"]

    fig, axs = plt.subplots(1, 3, figsize=(18, 5.8))
    sc.pl.umap(adata, color="leiden", legend_loc="on data", frameon=False, ax=axs[0], show=False, title=f"(A) Leiden (res={LEIDEN_RES:.2f})")
    sc.pl.umap(adata, color="EGFP_highlight", frameon=False, ax=axs[1], show=False, title="(B) GFP+ highlighted")
    sc.pl.umap(adata, color="Krt8ADI_highlight", frameon=False, ax=axs[2], show=False, title="(C) Krt8+ ADI highlighted")
    for ax in axs:
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
    plt.tight_layout()
    safe_savefig(fig, OUT_DIR / "Figure_Composite_UMAPs_res017.png", dpi=300)
    plt.close(fig)

    markers_present = {k: present_genes(adata.var_names, v) for k, v in CURATED_MARKERS.items()}
    missing = []
    for k, v in CURATED_MARKERS.items():
        mis = sorted(set(v) - set(markers_present[k]))
        if mis:
            missing.append({"group": k, "missing_genes": ",".join(mis)})
    if missing:
        pd.DataFrame(missing).to_csv(OUT_DIR / "Missing_Markers_Report.csv", index=False)

    if any(len(v) for v in markers_present.values()):
        sc.pl.dotplot(
            adata,
            var_names=markers_present,
            groupby="leiden",
            standard_scale="var",
            dendrogram=False,
            show=False,
            save="_DotPlot_CuratedMarkers_by_Leiden_res017.png",
        )

    stats = (
        adata.obs.assign(
            is_gfp=(adata.obs["GFP_Status"] == "GFP+").astype(int),
            is_krt8=adata.obs["Krt8ADI_highlight"].astype(str).eq("Krt8+ ADI").astype(int),
        )
        .groupby("leiden", observed=True)
        .agg(n_cells=("leiden", "size"), n_gfp_plus=("is_gfp", "sum"), n_krt8_pos=("is_krt8", "sum"))
        .reset_index()
    )
    stats["frac_GFP_plus"] = stats["n_gfp_plus"] / stats["n_cells"]
    stats["frac_Krt8_pos"] = stats["n_krt8_pos"] / stats["n_cells"]
    stats.to_csv(OUT_DIR / "Table_Fraction_GFPplus_Krt8pos_by_Leiden_res017.csv", index=False)

    plt.figure(figsize=(9.8, 4.2))
    ax = stats.set_index("leiden")["frac_GFP_plus"].sort_index().plot(kind="bar", rot=0)
    ax.set_ylabel("Fraction EGFP+")
    ax.set_xlabel("Leiden subcluster")
    ax.set_ylim(0, 1)
    plt.tight_layout()
    safe_savefig(plt.gcf(), OUT_DIR / "Figure_Bar_Fraction_GFPplus_by_Leiden_res017.png", dpi=300)
    plt.close()

    egfp_gene = next((g for g in ["EGFP", "Egfp", "eGFP", "GFP"] if g in adata.var_names), None)
    if egfp_gene:
        sc.pl.violin(
            adata,
            keys=egfp_gene,
            groupby="leiden",
            stripplot=False,
            jitter=False,
            log=True,
            show=False,
            save="_Violin_EGFP_by_Leiden_res017.png",
        )

    pd.crosstab(adata.obs["Krt8ADI_highlight"].astype(str), adata.obs["GFP_Status"]).to_csv(
        OUT_DIR / "Table_Crosstab_Krt8ADI_vs_GFPStatus_res017.csv"
    )

    print(f"[INFO] Outputs saved to: {OUT_DIR}")


if __name__ == "__main__":
    main()

