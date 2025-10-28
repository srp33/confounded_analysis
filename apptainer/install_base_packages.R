# install_base_packages.R
# This script installs the core R packages that take the longest to compile
# These will be baked into the base image for faster subsequent builds

# Configure repositories to prefer binary packages
options(repos = c(
    CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest",
    CRAN_source = "https://cran.rstudio.com/"
))

# Install pak first
if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak", repos = getOption("repos"))
}

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
core_infrastructure <- c(
    "Rcpp", "RcppArmadillo", "BH", "data.table", "Matrix", "foreach", "parallel", 
    "doParallel", "future", "matrixStats", "remotes"  # remotes for GitHub installations
)
pak::pkg_install(core_infrastructure)

# Batch 2: Tidyverse + heavy CRAN packages (can be installed together efficiently)
cat("--- Batch 2: Tidyverse + heavy CRAN packages ---\n")
heavy_cran_batch <- c(
    "tidyverse",  # This installs the whole ecosystem
    # Heavy statistical packages
    "glmnet", "e1071", "ranger", "caret", "magick", "stringi", 
    "RcppEigen", "lme4", "bigmemory", "mclust", "MCMCpack", "huge",
    # Heavy visualization packages
    "gridExtra", "gplots", "RColorBrewer", "plotly", "heatmaply", 
    "png", "colorspace", "scales", "R.devices", "readxl"
)
pak::pkg_install(heavy_cran_batch)

# Batch 3: All Bioconductor packages together (including annotation packages)
cat("--- Batch 3: All Bioconductor packages + annotation packages ---\n")
all_bioc_base <- c(
    # Core Bioconductor packages
    "bioc::SummarizedExperiment", "bioc::limma", "bioc::vsn", 
    "bioc::AnnotationDbi", "bioc::biomaRt", "bioc::ExperimentHub", 
    "bioc::ComplexHeatmap", "bmbolstad/preprocessCore",
    # Polyester dependencies (moved to base for better performance)
    "bioc::Biostrings", "bioc::IRanges", "bioc::S4Vectors", "bioc::zlibbioc",
    "logspline"
)
pak::pkg_install(all_bioc_base)

# Note: Custom annotation packages from mbni.org moved to separate build stage
# due to server reliability issues. They will be installed in the fast build stage.

cat("--- Base R packages installed successfully ---\n")
cat("Installed packages can be seen with: installed.packages()[,1]\n")