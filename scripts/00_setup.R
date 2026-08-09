# Install CRAN Packages
install.packages("pak")

CRAN_REQUIRED <- c(
  "tidyverse", "rio", "ggplot2", "ggrepel", "ggpubr", "pheatmap", "patchwork",
  "igraph", "ggraph", "scales", "httr", "jsonlite", "Matrix", "metafor",
  "survival", "survminer", "glmnet", "randomForest", "pROC", "caret",
  "VennDiagram", "ggVennDiagram", "msigdbr", "gridExtra",
  "circlize", "future"
)

pak::pkg_install(CRAN_REQUIRED)


# Bioconductor Packages
BIOC_REQUIRED <- c(
  "tximport",
  "DESeq2",
  "sva",
  "edgeR",
  "limma",
  "RUVSeq",
  "clusterProfiler",
  "enrichplot",
  "org.Hs.eg.db",
  "ReactomePA",
  "GSVA",
  "recount3",
  "SummarizedExperiment",
  "MultiAssayExperiment",
  "curatedTCGAData",
  "ComplexHeatmap",
  "genekitr",
  "TCGAbiolinks",
  "biomaRt",
  "GEOquery",
  "ExperimentHub",
  "AnnotationDbi",
  "ensembldb",
  "EnsDb.Hsapiens.v86"
)

pak::pkg_install(BIOC_REQUIRED)


# GitHub Packages
yes


# Optional Packages
BIOC_OPTIONAL <- c(
  "decoupleR",
  "multiMiR",
  "depmap",
  "cBioPortalData",
  "SingleR",
  "celldex",
  "slingshot",
  "AUCell"
)

pak::pkg_install(BIOC_OPTIONAL)