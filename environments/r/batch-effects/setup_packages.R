#!/usr/bin/env Rscript
# setup_packages.R
# Install additional R packages not included in the minimal Nix environment
# Run this script after activating the Nix environment

cat("=== Setting up additional R packages ===\n\n")

# Determine library path
lib <- Sys.getenv("R_LIBS_USER")
if (lib == "" || !dir.exists(dirname(lib))) {
  lib <- "/grphome/grp_batch_effects/environments/r/r-libs"
  cat(sprintf("Using shared library: %s\n", lib))
}

# Ensure library directory exists with group permissions
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
Sys.chmod(lib, mode = "0775")

# Set library path
.libPaths(c(lib, .libPaths()))

cat(sprintf("Library path: %s\n", lib))
cat(sprintf("R version: %s\n\n", R.version.string))

# Function to install if not present
install_if_missing <- function(pkg, source = "CRAN", repo_spec = NULL) {
  pkg_name <- if (is.null(repo_spec)) pkg else basename(repo_spec)
  
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    cat(sprintf("Installing %s from %s...\n", pkg_name, source))
    
    tryCatch({
      if (source == "CRAN") {
        install.packages(pkg, lib = lib, repos = "https://cloud.r-project.org", 
                        dependencies = TRUE, quiet = FALSE)
      } else if (source == "Bioconductor") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", lib = lib, repos = "https://cloud.r-project.org")
        }
        BiocManager::install(pkg, lib = lib, update = FALSE, ask = FALSE)
      } else if (source == "GitHub") {
        if (!requireNamespace("remotes", quietly = TRUE)) {
          install.packages("remotes", lib = lib, repos = "https://cloud.r-project.org")
        }
        remotes::install_github(repo_spec, lib = lib, upgrade = "never")
      }
      cat(sprintf("  ✓ %s installed successfully\n", pkg_name))
    }, error = function(e) {
      cat(sprintf("  ✗ Failed to install %s: %s\n", pkg_name, e$message))
    })
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg_name))
  }
}

cat("=== Installing additional packages ===\n\n")

# === Tidyverse Extensions ===
cat("--- Tidyverse Extensions ---\n")
install_if_missing("purrr", "CRAN")
install_if_missing("readxl", "CRAN")
install_if_missing("stringi", "CRAN")

# === Statistical Packages ===
cat("\n--- Statistical Packages ---\n")
install_if_missing("ranger", "CRAN")
install_if_missing("lme4", "CRAN")
install_if_missing("mclust", "CRAN")
install_if_missing("MCMCpack", "CRAN")
install_if_missing("huge", "CRAN")
install_if_missing("nnls", "CRAN")
install_if_missing("moments", "CRAN")
install_if_missing("pracma", "CRAN")
install_if_missing("logspline", "CRAN")

# === Visualization ===
cat("\n--- Visualization Packages ---\n")
install_if_missing("gplots", "CRAN")
install_if_missing("plotly", "CRAN")
install_if_missing("heatmaply", "CRAN")
install_if_missing("ggpubr", "CRAN")
install_if_missing("ggtext", "CRAN")
install_if_missing("corrplot", "CRAN")
install_if_missing("png", "CRAN")
install_if_missing("colorspace", "CRAN")
install_if_missing("scales", "CRAN")
install_if_missing("DT", "CRAN")
install_if_missing("kableExtra", "CRAN")

# === Bioconductor Packages ===
cat("\n--- Bioconductor Packages ---\n")
install_if_missing("vsn", "Bioconductor")
install_if_missing("AnnotationDbi", "Bioconductor")
install_if_missing("biomaRt", "Bioconductor")
install_if_missing("ExperimentHub", "Bioconductor")
install_if_missing("Biostrings", "Bioconductor")
install_if_missing("IRanges", "Bioconductor")
install_if_missing("S4Vectors", "Bioconductor")
install_if_missing("zlibbioc", "Bioconductor")

# === Machine Learning ===
cat("\n--- Machine Learning Packages ---\n")
install_if_missing("MLmetrics", "CRAN")
install_if_missing("lightgbm", "CRAN")
install_if_missing("parsnip", "CRAN")
install_if_missing("tidymodels", "CRAN")

# === Dimensionality Reduction ===
cat("\n--- Dimensionality Reduction ---\n")
install_if_missing("Rtsne", "CRAN")
install_if_missing("umap", "CRAN")

# === Additional Utilities ===
cat("\n--- Additional Utilities ---\n")
install_if_missing("RcppArmadillo", "CRAN")
install_if_missing("RcppEigen", "CRAN")
install_if_missing("BH", "CRAN")
install_if_missing("matrixStats", "CRAN")
install_if_missing("future", "CRAN")
install_if_missing("pak", "CRAN")

# === Optional Specialized Packages ===
cat("\n--- Optional Specialized Packages ---\n")
cat("(These may fail if dependencies are missing - that's OK)\n")
install_if_missing("Seurat", "CRAN")
install_if_missing("seqgendiff", "CRAN")

# === GitHub Packages (if not already in Nix) ===
cat("\n--- GitHub Packages ---\n")
# Note: preprocessCore and polyester should be in Nix via git_pkgs
# But we can install them here as backup if Nix installation failed
if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  cat("preprocessCore not found in Nix environment, installing from GitHub...\n")
  install_if_missing("preprocessCore", "GitHub", "bmbolstad/preprocessCore")
}

if (!requireNamespace("polyester", quietly = TRUE)) {
  cat("polyester not found in Nix environment, installing from GitHub...\n")
  install_if_missing("polyester", "GitHub", "alyssafrazee/polyester")
}

cat("\n=== Package setup complete! ===\n")
cat("\nInstalled packages are in: ", lib, "\n")
cat("\nTo verify installation, run:\n")
cat("  library(tidyverse)\n")
cat("  library(sva)\n")
cat("  library(limma)\n")
cat("  library(xgboost)\n")
cat("\nNote: Some optional packages may have failed - that's OK if you don't need them.\n")
