#!/usr/bin/env Rscript

# ---------------------------------------------------------------------
# install_packages.R
# Unified installer for CRAN, Bioconductor, and GitHub R packages
# ---------------------------------------------------------------------

options(error = function(e) { message("FATAL ERROR: ", e$message); quit(status = 1) })
options(Ncpus = parallel::detectCores())

# ---------------------------------------------------------------------
# 1. Set repository options (stable CRAN snapshot)
# ---------------------------------------------------------------------
options(
  repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/2024-01-01")
)

# ---------------------------------------------------------------------
# 2. Ensure pak is available
# ---------------------------------------------------------------------
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = "https://cran.rstudio.com/")
}

# Add Bioconductor repo to pak so it can see both CRAN and Bioc
pak::repo_add(bioc = "https://bioconductor.org/packages/3.19/bioc")

# ---------------------------------------------------------------------
# 3. Define package groups
# ---------------------------------------------------------------------
cran_and_github_pkgs <- c(
  "pacman", "dplyr", "readr", "stringr", "tidyr", "tibble",
  "ggplot2", "gridExtra", "png", "magick", "colorspace", "pracma",
  "kableExtra", "Rtsne", "argparse", "docstring", "R.devices",
  "doParallel", "readxl", "ggpubr", "itertools", "ComplexHeatmap",
  "ggtext", "ranger", "fairadapt", "rliger", "huge", "umap",
  "mclust", "future", "bigmemory", "gplots", "RColorBrewer",
  "corrplot", "DT", "plotly", "heatmaply", "moments",
  "matrixStats", "MCMCpack", "nnls", "glmnet", "e1071",
  "MLmetrics", "tidyverse", "data.table", "caret", "ROCR",
  "reshape2", "plyr", "rpart", "genefilter", "nnet", "RcppArmadillo",
  "foreach", "parallel", "doParallel", "scales", "Seurat",
  "devtools", "MatrixModels", "bmbolstad/preprocessCore"
)

bioc_pkgs <- c(
  "BiocManager", "sva", "SCAN.UPC", "SummarizedExperiment",
  "ExperimentHub", "limma", "vsn", "batchelor", "BatchQC",
  "GEOquery", "pd.hg.u133a", "pd.hg.u133.plus.2", "AnnotationDbi",
  "biomaRt", "DESeq2", "illuminaHumanv4.db",
  "hugene11sttranscriptcluster.db", "annotate", "DelayedMatrixStats",
  "RUVSeq", "EDASeq", "polyester"
)

# ---------------------------------------------------------------------
# 4. Install CRAN + GitHub packages via pak
# ---------------------------------------------------------------------
cat("--- Installing CRAN and GitHub packages with pak ---\n")
pak::pkg_install(unique(cran_and_github_pkgs))

# ---------------------------------------------------------------------
# 5. Install Bioconductor packages via BiocManager
# ---------------------------------------------------------------------
cat("--- Installing Bioconductor packages ---\n")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cran.rstudio.com/")

BiocManager::install(unique(bioc_pkgs), ask = FALSE, update = TRUE)

# ---------------------------------------------------------------------
# 6. Install custom annotation packages from MBNI
# ---------------------------------------------------------------------
cat("--- Installing custom annotation packages from mbni.org ---\n")
tmpDir <- tempdir()

pkgUrl1 <- "http://mbni.org/customcdf/25.0.0/entrezg.download/hgu133ahsentrezgprobe_25.0.0.tar.gz"
pkgUrl2 <- "http://mbni.org/customcdf/25.0.0/entrezg.download/hgu133plus2hsentrezgprobe_25.0.0.tar.gz"

pkgFiles <- file.path(tmpDir, basename(c(pkgUrl1, pkgUrl2)))

for (i in seq_along(pkgFiles)) {
  download.file(c(pkgUrl1, pkgUrl2)[i], pkgFiles[i], quiet = TRUE)
  pak::pkg_install(pkgFiles[i])
}

cat("✅ All R packages installed successfully.\n")
