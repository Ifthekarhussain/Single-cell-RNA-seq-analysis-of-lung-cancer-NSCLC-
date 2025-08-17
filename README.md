# Single-cell RNA-seq Analysis of Lung Cancer (NSCLC)

This project showcases my skills in **single-cell RNA-seq data analysis** using R and popular bioinformatics tools.  
It is based on a dataset from **non-small cell lung cancer (NSCLC)** and demonstrates a complete end-to-end workflow including:

---

## 🧰 Tools & Skills Demonstrated

- **Seurat** – Data QC, normalization, clustering, UMAP, marker gene analysis  
- **SingleR** – Automated cell type annotation using Human Primary Cell Atlas (HPCA)  
- **Monocle3** – Pseudotime and trajectory inference  
- **clusterProfiler** – GO and KEGG enrichment analysis  
- **R / Bioconductor / Tidyverse / HPC workflows**

---

## 📌 Project Overview

This analysis includes:
- Preprocessing of 10X Genomics data  
- Quality control (`nFeature_RNA`, `percent.mt`)
- Identification of highly variable genes  
- Clustering (SNN-based) and UMAP visualization  
- Cell type annotation with **SingleR + celldex**  
- Marker gene detection and heatmaps  
- Functional enrichment (GO/KEGG)  
- T-cell subclustering  
- Trajectory inference with Monocle3

---

## 📂 Structure

NSCLC_scRNAseq_project/
├── scripts/ # R scripts (pipeline, installer)
├── data/ # Input files (10X .h5)
├── results/ # Output plots, marker tables
├── docs/ # Notes, slides, etc.
├── README.md # This file
├── LICENSE # MIT license
└── .gitignore

## ▶️ How to Run the Pipeline

1. Install dependencies:
```r
source("scripts/install_packages.R")

Run the main pipeline:
source("scripts/NSCLC_scRNAseq_pipeline.R")

h5_path <- "data/your_file.h5"
source("scripts/NSCLC_scRNAseq_pipeline.R")

```

Presentation

Slides explaining the workflow and key results: 📊 [Click here to view the project presentation (PPTX)](https://github.com/Ifthekarhussain/Single-cell-RNA-seq-analysis-of-lung-cancer-NSCLC-/blob/main/Introduction%20to%20Single-Cell%20RNA-Seq%20with%20Seurat.pptx?raw=true)






