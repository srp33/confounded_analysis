# install_combatseq_packages.R

# Set number of parallel workers
options(Ncpus = parallel::detectCores())

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Set Bioconductor version explicitly
BiocManager::install(version = "3.11", ask = FALSE, update = FALSE)

# Optional: Set CRAN snapshot mirror for reproducibility (June 2020)
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2020-06-01"))

# ----- Install CRAN packages -----
cran_pkgs <- c(
  "pacman", "dplyr", "readr", "stringr", "tidyr", "tibble", "ggplot2",
  "gridExtra", "png", "magick", "colorspace", "pracma", "kableExtra",
  "Rtsne", "argparse", "docstring", "R.devices", "doParallel", "readxl",
  "itertools", "ggpubr", "MCMCpack", "nnls", "glmnet", "e1071", "MLmetrics",
  "tidyverse", "data.table", "caret", "annotate", "ROCR", "reshape2",
  "plyr", "rpart", "nnet", "foreach", "parallel", "scales"
)

install.packages(cran_pkgs)

# ----- Install Bioconductor packages -----
bioc_pkgs <- c(
  "BatchQC", "sva", "SCAN.UPC", "SummarizedExperiment", "ExperimentHub",
  "limma", "vsn", "ComplexHeatmap", "DelayedMatrixStats", "RUVSeq",
  "GEOquery", "pd.hg.u133a", "pd.hg.u133.plus.2", "AnnotationDbi",
  "biomaRt", "DESeq2", "illuminaHumanv4.db", "hugene11sttranscriptcluster.db"
)

BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

# ----- GitHub packages -----
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# These may or may not be available for Bioc 3.11 – try with caution
# Check if still needed
remotes::install_github("bmbolstad/preprocessCore")

# Note: rliger is not installable in R 4.0/Bioc 3.11 without custom handling.
# You may want to skip or install manually if you need it.


message("--- All R packages installed successfully ---")
