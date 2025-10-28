# install_base_packages.R
# This script installs the core R packages that take the longest to compile
# These will be baked into the base image for faster subsequent builds

cat("=== STARTING R PACKAGE INSTALLATION ===\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("Installation started at:", as.character(Sys.time()), "\n")

# Configure repositories to prefer binary packages
cat("Configuring repositories...\n")
options(repos = c(
    CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest",
    CRAN_source = "https://cran.rstudio.com/"
))
cat("Repositories configured:", getOption("repos"), "\n")

# Install pak first
cat("Installing pak package manager...\n")
if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak", repos = getOption("repos"))
}
cat("pak installation completed\n")

# Use all available CPU cores for parallel installation
num_cores <- parallel::detectCores()
cat("Detected", num_cores, "CPU cores for base package installation\n")
options(Ncpus = num_cores)

# Configure pak for maximum parallelization
options(pak.no_extra_messages = TRUE)
Sys.setenv(PKG_BUILD_EXTRA_FLAGS = paste0("-j", num_cores))
Sys.setenv(MAKEFLAGS = paste0("-j", num_cores))

# OPTIMIZED: Install packages in strategic batches for maximum parallelization

# Batch 1: Core infrastructure (install first as many packages depend on these)
cat("--- Batch 1: Core infrastructure packages ---\n")
cat("Batch 1 started at:", as.character(Sys.time()), "\n")
core_infrastructure <- c(
    "Rcpp", "RcppArmadillo", "BH", "data.table", "Matrix", "foreach", "parallel", 
    "doParallel", "future", "matrixStats", "remotes"  # remotes for GitHub installations
)
cat("Installing packages:", paste(core_infrastructure, collapse = ", "), "\n")
pak::pkg_install(core_infrastructure)
cat("Batch 1 completed at:", as.character(Sys.time()), "\n")

# Batch 2: Tidyverse + heavy CRAN packages (can be installed together efficiently)
cat("--- Batch 2: Tidyverse + heavy CRAN packages ---\n")
cat("Batch 2 started at:", as.character(Sys.time()), "\n")
heavy_cran_batch <- c(
    "tidyverse",  # This installs the whole ecosystem
    # Heavy statistical packages
    "glmnet", "e1071", "ranger", "caret", "magick", "stringi", 
    "RcppEigen", "lme4", "bigmemory", "mclust", "MCMCpack", "huge",
    # Heavy visualization packages
    "gridExtra", "gplots", "RColorBrewer", "plotly", "heatmaply", 
    "png", "colorspace", "scales", "R.devices", "readxl"
)
cat("Installing packages:", paste(heavy_cran_batch, collapse = ", "), "\n")
pak::pkg_install(heavy_cran_batch)
cat("Batch 2 completed at:", as.character(Sys.time()), "\n")

# Batch 3: Bioconductor packages
cat("--- Batch 3: Bioconductor packages ---\n")
cat("Batch 3 started at:", as.character(Sys.time()), "\n")
all_bioc_base <- c(
    # Core Bioconductor packages
    "bioc::SummarizedExperiment", "bioc::limma", "bioc::vsn", 
    "bioc::AnnotationDbi", "bioc::biomaRt", "bioc::ExperimentHub", 
    "bioc::ComplexHeatmap", "bmbolstad/preprocessCore",

    "bioc::Biostrings", "bioc::IRanges", "bioc::S4Vectors", "bioc::zlibbioc",
    "logspline"
)
cat("Installing packages:", paste(all_bioc_base, collapse = ", "), "\n")
pak::pkg_install(all_bioc_base)
cat("Batch 3 completed at:", as.character(Sys.time()), "\n")


cat("=== BASE R PACKAGES INSTALLED SUCCESSFULLY ===\n")
cat("Installation completed at:", as.character(Sys.time()), "\n")
cat("Total packages installed:", length(installed.packages()[,1]), "\n")
cat("Installed packages can be seen with: installed.packages()[,1]\n")