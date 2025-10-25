# install_combatseq_packages.R

# Set number of parallel workers
options(Ncpus = parallel::detectCores())

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Set Bioconductor version explicitly
suppressMessages(BiocManager::install(version = "3.11", ask = FALSE, update = FALSE))

# Optional: Set CRAN snapshot mirror for reproducibility (June 2020)
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2020-06-01"))

# ----- Install CRAN packages -----
cran_pkgs <- c(
  "pacman", "dplyr", "readr", "stringr", "tidyr", "tibble", "ggplot2",
  "gridExtra", "png", "magick", "colorspace", "pracma", "kableExtra",
  "Rtsne", "argparse", "docstring", "R.devices", "doParallel", "readxl",
  "itertools", "ggpubr", "MCMCpack", "nnls", "glmnet", "e1071", "MLmetrics",
  "tidyverse", "data.table", "caret", "annotate", "ROCR", "reshape2",
  "plyr", "rpart", "nnet", "foreach", "parallel", "scales", "vroom", 
  "future", "huge", "preprocessCore"
)

install.packages(cran_pkgs)

# ----- Install Bioconductor packages -----
bioc_pkgs <- c(
  "BatchQC", "sva", "SCAN.UPC", "SummarizedExperiment", "ExperimentHub",
  "limma", "vsn", "ComplexHeatmap", "DelayedMatrixStats", "RUVSeq",
  "GEOquery", "pd.hg.u133a", "pd.hg.u133.plus.2", "AnnotationDbi",
  "biomaRt", "DESeq2", "illuminaHumanv4.db", "hugene11sttranscriptcluster.db", 
  "batchelor"
)

BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

# ----- GitHub packages -----
# Install remotes (for GitHub installs)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install preprocessCore (if not covered above)
if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  remotes::install_github("bmbolstad/preprocessCore")
}

# Install Seurat — specific version from ~2020 if possible
# Use CRAN version from 2020 snapshot; no need for GitHub
install.packages("Seurat")

# Try to install fairadapt (no CRAN version)
tryCatch({
  remotes::install_github("YosefLab/fairadapt")
}, error = function(e) {
  message("WARNING: fairadapt failed to install. You may need to install it manually or skip it.")
})

# Try to install rliger (may fail in R 4.0)
tryCatch({
  remotes::install_github("MacoskoLab/liger")
}, error = function(e) {
  message("WARNING: rliger failed to install. You may need to install it manually or exclude methods using it.")
})

message("--- All R packages installation attempted ---")