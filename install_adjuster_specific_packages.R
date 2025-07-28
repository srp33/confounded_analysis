# Install pak if it's not already available
if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak")
}

# Install packages in dependency order to avoid conflicts
# First install base packages
base_packages <- c(
    # Bioconductor Packages
    "limma", "vsn", "doParallel", "ggplot2", "ggpubr", 
    "BatchQC", "batchelor",
    
    # CRAN Packages
    "ranger", "fairadapt", "rliger", "huge",
    "MASS", "Rtsne", "umap"
)

pak::pkg_install(base_packages)

# Then install Seurat
pak::pkg_install("satijalab/seurat@seurat5")

# Finally install SeuratData (depends on Seurat)
pak::pkg_install("satijalab/seurat-data@seurat5")