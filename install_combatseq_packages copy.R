# install_combatseq_packages.R

# Install pak from CRAN if it's not already installed.
if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak", repos = "https://cloud.r-project.org")
}

# Load the pak library
library(pak)

# --- Set Version Constraints for Reproducibility (Crucial Step RE-ADDED) ---

# CRAN Snapshot: Re-add the specific date to find package versions that 
# were compatible with Bioconductor 3.11 when it was current.
pak::repo_set(CRAN = "https://packagemanager.posit.co/cran/2020-06-01")

# --- Define Package Lists (Keep as is) ---
all_pkgs <- c(
  "pacman", "tidyverse", "vroom", "data.table", "readxl", "reshape2",
  "gridExtra", "png", "magick", "colorspace", "scales", "ggpubr", 
  "pracma", "kableExtra", "Rtsne", "argparse", "docstring", 
  "R.devices", "doParallel", "itertools", "future", "parallel", "foreach", 

  "MCMCpack", "nnls", "glmnet", "e1071", "MLmetrics", "caret", "ROCR", 
  "plyr", "rpart", "nnet", "huge", 

  "BatchQC", "sva", "SCAN.UPC", "SummarizedExperiment", "ExperimentHub",
  "limma", "vsn", "ComplexHeatmap", "DelayedMatrixStats", "RUVSeq",
  "GEOquery", "pd.hg.u133a", "pd.hg.u133.plus.2", "AnnotationDbi",
  "biomaRt", "DESeq2", "illuminaHumanv4.db", "hugene11sttranscriptcluster.db", 
  "batchelor", "clusterProfiler", "DOSE", "fgsea", "annotate",

  # GitHub Package
  "bmbolstad/preprocessCore"
)

# --- Execute Single, Unified Installation ---

message("--- Starting unified package installation with pak ---")
message(paste("--- Using CRAN Snapshot:", pak::pak_opt("repos")["CRAN"], "---"))
message("--- Assuming R 4.0 for Bioconductor 3.11 compatibility ---")

# pkg_install handles CRAN, Bioconductor, and GitHub in parallel
pak::pkg_install(all_pkgs)

message("--- All R packages installation attempted and resolved by pak ---")