if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

tampor_packages = c("limma", "vsn", "doParallel", "ggplot2", "ggpubr")

BiocManager::install(tampor_packages)
