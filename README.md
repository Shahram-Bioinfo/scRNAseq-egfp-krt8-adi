**scRNA-seq computational analysis (EGFP / Krt8 ADI)**

This repository contains the computational scripts used to perform the single-cell RNA sequencing (scRNA-seq) analyses described in the Methods section of the study entitled:

Transcriptomic Signature-Guided Depletion of Intermediate Alveolar Epithelial Cells Ameliorates Pulmonary Fibrosis in Mice.

All analyses were conducted in Python using established open-source libraries (including Scanpy, SciPy, Seaborn, and Matplotlib). The scripts are provided to reproduce the analysis workflow by integrating existing software packages; no novel algorithms were developed.

Scope of the analysis

The workflow implemented in this repository includes the following computational steps:

Quality control and preprocessing of scRNA-seq data, including filtering of low-quality cells, normalization to a total of 10,000 reads per cell, and log-transformation.

Identification of EGFP-positive (EGFP⁺) cells based on non-zero expression of the EGFP gene.

Dimensionality reduction using principal component analysis (PCA), construction of a k-nearest neighbor graph, unsupervised clustering using the Leiden algorithm (resolution = 0.17), and two-dimensional visualization using Uniform Manifold Approximation and Projection (UMAP).

Annotation of Krt8⁺ alveolar differentiation intermediate (ADI)-like cells using a curated panel of structural cytokeratin and injury-response marker genes based on previously published criteria.

Differential gene expression analysis using the Wilcoxon rank-sum test with multiple-testing correction.

Assessment of transcriptomic similarity between EGFP⁺ and Krt8⁺ ADI-like cell populations using Pearson correlation, Spearman correlation, and cosine similarity metrics.

Generation of visualization outputs, including UMAP embeddings, violin plots, dot plots, and summary tables.

Execution order

The scripts are intended to be executed sequentially, with each step producing outputs that serve as inputs for subsequent analyses:

scripts/01_build_tulane_normalized.py
Loads the scRNA-seq dataset, performs quality control (including mitochondrial-content filtering), normalization, and log-transformation, and writes a processed AnnData object.

scripts/02_label_krt8_adi_like.py
Identifies EGFP⁺ cells and annotates Krt8⁺ ADI-like cells using the predefined marker gene panel.

scripts/03_leiden_umap_and_plots.py
Performs PCA, constructs the neighborhood graph, applies Leiden clustering (resolution = 0.17), computes UMAP embeddings, and generates visualization outputs (UMAP panels, violin plots, dot plots, and per-cluster summaries).

scripts/04_degs_and_similarity.py
Conducts differential gene expression analysis, computes transcriptomic similarity metrics, and produces summary tables and heatmaps.

Reproducibility and usage

All file paths within the scripts are defined relative to the project root to ensure portability across computing environments. Users can reproduce the analysis by providing the required input datasets in the designated directories and executing the scripts in the specified order. The repository is intended to support transparency and reproducibility of the computational analyses.
