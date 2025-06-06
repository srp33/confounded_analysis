if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

tampor_packages = c("limma", "vsn", "doParallel", "ggplot2", "ggpubr")

BiocManager::install(tampor_packages)

BiocManager::install("BatchQC")

devtools::install_github("satijalab/seurat", "seurat5")
devtools::install_github("satijalab/seurat-data", "seurat5")

