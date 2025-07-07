if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

tampor_packages = c("limma", "vsn", "doParallel", "ggplot2", "ggpubr")

BiocManager::install(tampor_packages)

BiocManager::install("BatchQC")

BiocManager::install(c("batchelor", "rliger"))

devtools::install_github("satijalab/seurat", "seurat5")
devtools::install_github("satijalab/seurat-data", "seurat5")


# Needed for FairAdapt
install.packages("ranger")
install.packages("fairadapt")