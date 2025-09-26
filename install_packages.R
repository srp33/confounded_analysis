# install_packages.R
# This script consolidates all R package installations into a single step.

# Install pak if it's not already available.
if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak", repos = "https://cran.rstudio.com/")
}

# Use all available CPU cores for parallel installation.
options(Ncpus = parallel::detectCores())

# Define all packages that can be installed together.
# Seurat and SeuratData are excluded here and will be installed sequentially
# to resolve a dependency issue during the Docker build process.
all_other_packages <- c(
    # From install_main_packages.R
    "pacman", "dplyr", "readr", "stringr", "tidyr", "tibble", "ggplot2", 
    "gridExtra", "png", "magick", "colorspace", "pracma", "kableExtra", 
    "Rtsne", "argparse", "docstring", "R.devices", "doParallel", "readxl", 
    "sva", "SCAN.UPC", "SummarizedExperiment", "ExperimentHub", "limma", 
    "vsn", "ggpubr", "itertools", "ComplexHeatmap", "ggtext",

    # From install_adjuster_specific_packages.R
    "BatchQC", "batchelor", "ranger", "fairadapt", "rliger", "huge",
    "MASS", "umap", "mclust", "future", "bigmemory",
    
    # From install_annotation_packages.R
    "GEOquery", "pd.hg.u133a", "pd.hg.u133.plus.2", "AnnotationDbi",

    # Other GitHub packages
    "bmbolstad/preprocessCore",
    "Seurat"
)

# Install the bulk of packages first.
cat("--- Installing all CRAN, Bioconductor, and GitHub packages ---\n")
pak::pkg_install(unique(all_other_packages))

# Handle the custom source packages from mbni.org.
# Download first, then use pak to install from the local file.
cat("--- Installing custom annotation packages from source ---\n")
tmpDir <- tempdir()

# Package 1: hgu133ahsentrezgprobe
pkgUrl1 <- "http://mbni.org/customcdf/25.0.0/entrezg.download/hgu133ahsentrezgprobe_25.0.0.tar.gz"
pkgFilePath1 <- file.path(tmpDir, basename(pkgUrl1))
download.file(pkgUrl1, pkgFilePath1)
pak::pkg_install(pkgFilePath1)

# Package 2: hgu133plus2hsentrezgprobe
pkgUrl2 <- "http://mbni.org/customcdf/25.0.0/entrezg.download/hgu133plus2hsentrezgprobe_25.0.0.tar.gz"
pkgFilePath2 <- file.path(tmpDir, basename(pkgUrl2))
download.file(pkgUrl2, pkgFilePath2)
pak::pkg_install(pkgFilePath2)

cat("--- All R packages installed successfully. ---\n")
