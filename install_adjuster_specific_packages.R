# Install BiocManager if it's not already present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# --- Bioconductor Packages ---
# A single, combined list of all Bioconductor packages to install
bioc_packages <- c(
    "limma", "vsn", "doParallel", "ggplot2", "ggpubr", 
    "BatchQC", "batchelor" # rliger was removed from here
)
BiocManager::install(bioc_packages, update = FALSE, ask = FALSE)

# --- CRAN Packages ---
# Install standard R packages from CRAN
install.packages(c("ranger", "fairadapt"))
install.packages("rliger") 

# --- GitHub Packages ---
# Install packages from GitHub using devtools
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("satijalab/seurat", ref = "seurat5")
devtools::install_github("satijalab/seurat-data", ref = "seurat5")