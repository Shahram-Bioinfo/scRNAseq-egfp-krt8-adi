#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
OUT_BASE = PROJECT_ROOT / "outputs" / "umap_leiden_res017"
OUT_BASE.mkdir(parents=True, exist_ok=True)

IN_H5AD = PROCESSED_DIR / "adata_tulane_normalized_labeled.h5ad"

LEIDEN_RES = 0.17
RANDOM_STATE = 0
N_PCS = 30
N_NEIGHBORS = 15

GRAY = "#D3D3D3"
BLUE = "#1f77b4"

KRT8_PROG: List[str] = [
    "Krt8","Krt18","Krt7","Krt19",
    "Sprr1a","Areg","Hbegf","Tnfrsf12a","Itgb6",
    "Cldn18","Cldn4","Cryab","Cdkn1a",
    "S100a6","S100a11","Emp2","Icam1","Gprc5a","Nupr1","Sfn"
]

CURATED_MARKERS: Dict[str, List[str]] = {
    "AT1": ["Pdpn","Ager"],
    "AT2": ["Sftpc","Sftpa1","Sftpb","Lamp3"],
    "Krt8_ADI": KRT8_PROG,
    "Basal": ["Krt5","Krt14","Trp63"],
    "Club": ["Scgb1a1","Scgb3a2"],
    "Ciliated": ["Foxj1","Tekt1"],
    "Prolif": ["Mki67","Top2a","Ccnb1"],
    "Fibro": ["Col1a1","Col3a1","Dcn","Pdgfra"],
    "Endoth": ["Pecam1","Kdr"],
    "Macro": ["Adgre1","Lyz2"],
    "Neutro": ["S100a8","S100a9"],
    "T_cells": ["Cd3d","Cd3e","Trac"],
    "B_cells": ["Ms4a1","Cd79a"],
    "NK": ["Nkg7","Klrb1c"],
    "DC": ["Itgax","H2-Aa"]
}

def to_numpy(x):
    if hasattr(x, "toarray"):
        return x.toarray()
    if hasattr(x, "A"):
        return x.A
    return np.asarray(x)

def present_genes(var_names, genes):
    s = set(map(str, var_names))
    return [g for g in genes if g in s]

def ensure_umap(adata):
    if "X_umap" in adata.obsm:
        return
    sc.pp.pca(adata, n_comps=N_PCS, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=RANDOM_STATE)
    sc.tl.umap(adata, random_state=RANDOM_STATE)

def main():
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = OUT_BASE / stamp
    out_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(out_dir)
    sc.settings.autoshow = False

    adata = sc.read_h5ad(str(IN_H5AD))

    if "GFP_Status" not in adata.obs:
        egfp_gene = next((g for g in ["EGFP","Egfp","eGFP","GFP"] if g in adata.var_names), None)
        if egfp_gene is None:
            raise RuntimeError("GFP_Status missing and EGFP gene not found in var_names.")
        xg = to_numpy(adata[:, egfp_gene].X).ravel()
        adata.obs["GFP_Status"] = np.where(xg > 0, "GFP+", "GFP-").astype("category")

    raw_krt8_col = None
    krt8_bool = None
    if "Krt8+ADI-like" in adata.obs:
        raw_krt8_col = "Krt8+ADI-like"
        ss = adata.obs[raw_krt8_col].astype(str)
        krt8_bool = (ss == "Krt8+ADI-like") | (ss.str.lower().isin(["true","1","yes","y","t"]))
        if int(np.sum(krt8_bool)) == 0:
            krt8_bool = None

    if krt8_bool is None:
        genes_present = present_genes(adata.var_names, KRT8_PROG)
        if len(genes_present) < 5:
            raise RuntimeError(f"Too few Krt8 program genes present ({len(genes_present)}).")
        sc.tl.score_genes(adata, gene_list=genes_present, score_name="Krt8ADI_score", use_raw=False)
        thr = float(np.quantile(adata.obs["Krt8ADI_score"].to_numpy(), 0.90))
        krt8_bool = adata.obs["Krt8ADI_score"].to_numpy() > thr
        ensure_umap(adata)
        sc.pl.umap(adata, color="Krt8ADI_score", show=False, save="_Krt8ADI_score.png")

    adata.obs["Krt8ADI_highlight"] = pd.Categorical(
        np.where(krt8_bool, "Krt8+ ADI", "Other"),
        categories=["Other","Krt8+ ADI"], ordered=True
    )

    ensure_umap(adata)
    sc.tl.leiden(adata, resolution=LEIDEN_RES, random_state=RANDOM_STATE)

    adata.obs["EGFP_highlight"] = pd.Categorical(
        np.where(adata.obs["GFP_Status"].astype(str) == "GFP+", "GFP+", "Other"),
        categories=["Other","GFP+"], ordered=True
    )

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    sc.pl.umap(adata, color="leiden", legend_loc="on data", title="Leiden", ax=axs[0], show=False)
    sc.pl.umap(adata, color="EGFP_highlight", palette=[GRAY, BLUE], title="EGFP+", ax=axs[1], show=False)
    sc.pl.umap(adata, color="Krt8ADI_highlight", palette=[GRAY, BLUE], title="Krt8+ ADI", ax=axs[2], show=False)
    for ax in axs:
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
    plt.tight_layout()
    plt.savefig(out_dir / "PanelC_res017_Highlighted.png", dpi=300)
    plt.close()

    missing_rows = []
    markers_present = {}
    for group, genes in CURATED_MARKERS.items():
        pres = present_genes(adata.var_names, genes)
        markers_present[group] = pres
        mis = sorted(set(genes) - set(pres))
        if mis:
            missing_rows.append({"group": group, "missing_genes": ",".join(mis)})
    if missing_rows:
        pd.DataFrame(missing_rows).to_csv(out_dir / "Missing_Markers_Report.csv", index=False)

    sc.pl.dotplot(
        adata,
        var_names=markers_present,
        groupby="leiden",
        standard_scale="var",
        dendrogram=False,
        show=False,
        save="__DotPlot_CuratedMarkers_by_Leiden.png",
    )

    egfp_gene = next((g for g in ["EGFP","Egfp","eGFP","GFP"] if g in adata.var_names), None)
    if egfp_gene is not None:
        sc.pl.violin(
            adata,
            keys=[egfp_gene],
            groupby="leiden",
            stripplot=False,
            jitter=False,
            log=True,
            show=False,
            save="__Violin_EGFP_by_Leiden.png",
        )

    stats = (
        adata.obs[["leiden","EGFP_highlight"]]
        .assign(is_gfp=(adata.obs["EGFP_highlight"] == "GFP+").astype(int))
        .groupby("leiden", observed=True)
        .agg(n_cells=("leiden","size"), n_gfp_plus=("is_gfp","sum"))
        .reset_index()
    )
    stats["fraction_GFP_plus"] = stats["n_gfp_plus"] / stats["n_cells"]
    stats["cluster_weight"] = stats["n_cells"] / stats["n_cells"].sum()
    stats["weighted_contrib"] = stats["fraction_GFP_plus"] * stats["cluster_weight"]

    overall = pd.DataFrame({
        "leiden": ["ALL"],
        "n_cells": [int(stats["n_cells"].sum())],
        "n_gfp_plus": [int(stats["n_gfp_plus"].sum())],
        "fraction_GFP_plus": [stats["n_gfp_plus"].sum() / stats["n_cells"].sum()],
        "cluster_weight": [1.0],
        "weighted_contrib": [stats["n_gfp_plus"].sum() / stats["n_cells"].sum()],
    })
    pd.concat([stats, overall], ignore_index=True).to_csv(out_dir / "Leiden_GFP_fraction.csv", index=False)

    plt.figure(figsize=(9, 4))
    stats.set_index("leiden")["fraction_GFP_plus"].sort_index().plot(kind="bar")
    plt.ylabel("Fraction EGFP+")
    plt.xlabel("Leiden cluster")
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(out_dir / "Bar_Fraction_GFPplus_by_Leiden.png", dpi=300)
    plt.close()

    ov = (
        adata.obs[["leiden","EGFP_highlight","Krt8ADI_highlight"]]
        .assign(is_gfp=(adata.obs["EGFP_highlight"] == "GFP+").astype(int))
        .assign(is_krt8=(adata.obs["Krt8ADI_highlight"] == "Krt8+ ADI").astype(int))
    )
    ov["is_both"] = (ov["is_gfp"] & ov["is_krt8"]).astype(int)
    ov = (
        ov.groupby("leiden", observed=True)
          .agg(
              n_cells=("leiden","size"),
              n_gfp_plus=("is_gfp","sum"),
              n_krt8_plus=("is_krt8","sum"),
              n_both=("is_both","sum"),
          )
          .reset_index()
    )
    ov.to_csv(out_dir / "Leiden_EGFP_Krt8_overlap.csv", index=False)

    score_cols: List[str] = []
    for group, genes in CURATED_MARKERS.items():
        genes_present = present_genes(adata.var_names, genes)
        if len(genes_present) == 0:
            continue
        score_name = f"score__{group}"
        sc.tl.score_genes(adata, gene_list=genes_present, score_name=score_name, use_raw=False)
        score_cols.append(score_name)

    if len(score_cols) > 0:
        mean_scores = (
            adata.obs.groupby("leiden", observed=True)[score_cols]
            .mean()
            .fillna(-np.inf)
        )
        best_types = mean_scores.idxmax(axis=1)
        cluster_celltype_map = {str(k): v.replace("score__","") for k, v in best_types.to_dict().items()}
        mean_scores.to_csv(out_dir / "Leiden_celltype_mean_scores.csv")
    else:
        cluster_celltype_map = {str(k): "unknown" for k in adata.obs["leiden"].astype(str).unique()}

    adata.obs["celltype"] = pd.Categorical(adata.obs["leiden"].astype(str).map(cluster_celltype_map))
    adata.uns["cluster_celltype_map"] = cluster_celltype_map
    adata.uns["marker_scores_columns"] = score_cols

    pd.DataFrame({"leiden": list(cluster_celltype_map.keys()), "celltype": list(cluster_celltype_map.values())}).to_csv(
        out_dir / "Leiden_to_CellType_map.csv", index=False
    )

    adata.write(out_dir / f"adata_with_leiden_res017_{stamp}.h5ad", compression="gzip")
    print(str(out_dir))

if __name__ == "__main__":
    main()


