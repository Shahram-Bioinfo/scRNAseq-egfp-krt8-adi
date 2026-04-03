# scRNA-seq Computational Analysis (EGFP / Krt8 ADI)

This repository contains the computational workflows used for single-cell RNA sequencing (scRNA-seq) analysis described in the study:

Transcriptomic Signature-Guided Depletion of Intermediate Alveolar Epithelial Cells Ameliorates Pulmonary Fibrosis in Mice

All analyses were implemented in Python using Scanpy, SciPy, Seaborn, and Matplotlib. No novel algorithms were developed.

## Scope of Analysis

### 1. Quality Control and Preprocessing
- Filtering of low-quality cells  
- Mitochondrial content filtering  
- Normalization to 10,000 reads per cell  
- Log-transformation  

### 2. EGFP⁺ Cell Identification
- Based on non-zero expression of EGFP  

### 3. Dimensionality Reduction and Clustering
- PCA  
- kNN graph  
- Leiden clustering (resolution = 0.17)  
- UMAP visualization  

### 4. Annotation of Krt8⁺ ADI-like Cells
- Based on cytokeratin and injury-response markers  

### 5. Differential Expression Analysis
- Wilcoxon rank-sum test with correction  

### 6. Transcriptomic Similarity
- Pearson, Spearman, cosine similarity  

### 7. Visualization Outputs
- UMAP, violin plots, dot plots, heatmaps  

## Execution

Run in order:

```
python scripts/01_build_tulane_normalized.py
python scripts/02_label_krt8_adi_like.py
python scripts/03_leiden_umap_and_plots.py
python scripts/04_degs_and_similarity.py
```

## Requirements

```
pip install scanpy scipy seaborn matplotlib
```
