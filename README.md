[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18022192.svg)](https://doi.org/10.5281/zenodo.18022192)

# scRNA-seq computational analysis (EGFP / Krt8 ADI)

This repository contains the computational scripts used to perform the single-cell RNA sequencing (scRNA-seq) analyses described in the Methods section of the study entitled:

**Transcriptomic Signature-Guided Depletion of Intermediate Alveolar Epithelial Cells Ameliorates Pulmonary Fibrosis in Mice**.

All analyses were conducted in Python using established open-source libraries (including Scanpy, SciPy, Seaborn, and Matplotlib). The scripts are provided to reproduce the analysis workflow by integrating existing software packages; no novel algorithms were developed.


## Scope of the analysis

The workflow implemented in this repository includes the following computational steps:

1- Quality control and preprocessing of scRNA-seq data, including filtering of low-quality cells, normalization to a total of 10,000 reads per cell, and log-transformation.

2- Identification of EGFP-positive (EGFP⁺) cells based on non-zero expression of the EGFP gene.

3- Dimensionality reduction using principal component analysis (PCA), construction of a k-nearest neighbor graph, unsupervised clustering using the Leiden algorithm (resolution = 0.17), and two-dimensional visualization using Uniform Manifold Approximation and Projection (UMAP).

4- Annotation of Krt8⁺ alveolar differentiation intermediate (ADI)-like cells using a curated panel of structural cytokeratin and injury-response marker genes based on previously published criteria.

5- Differential gene expression analysis using the Wilcoxon rank-sum test with multiple-testing correction.

6- Assessment of transcriptomic similarity between EGFP⁺ and Krt8⁺ ADI-like cell populations using Pearson correlation, Spearman correlation, and cosine similarity metrics.

7- Generation of visualization outputs, including UMAP embeddings, violin plots, and heatmaps.

## Execution order

The scripts are intended to be executed sequentially, with each step producing outputs that serve as inputs for subsequent analyses:

1. **scripts/01_build_and_label_tulane.py**  
   Loads the scRNA-seq dataset, performs quality control, normalization, and log-transformation, and generates a processed AnnData object.

2. **scripts/02_label_egfp_adi_tulane.py**  
   Identifies EGFP⁺ cells and annotates Krt8⁺ ADI-like cells using the predefined marker gene panel.

3. **scripts/03_leiden_umap_visualization.py**  
   Performs PCA, constructs the neighborhood graph, applies Leiden clustering (resolution = 0.17), computes UMAP embeddings, and generates visualization outputs.

4. **scripts/04_degs_and_similarity.py**  
   Conducts differential gene expression analysis, computes transcriptomic similarity metrics, and produces summary tables and heatmaps.


## Reproducibility and usage

All file paths within the scripts are defined relative to the project root to ensure portability across computing environments. Users can reproduce the analysis by providing the required input datasets in the designated directories and executing the scripts in the specified order. The repository is intended to support transparency and reproducibility of the computational analyses.