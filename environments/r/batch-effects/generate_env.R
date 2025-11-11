# generate_env.R
# Phase 1 (Authoring): Generate default.nix and .Rprofile for batch-effects environment
# This script uses rix to create a reproducible R environment specification

library(rix)

cat("=== Generating batch-effects R environment ===\n")
cat("This will create default.nix and .Rprofile files\n")
cat("Date pinning: 2023-04-01 (Latest available in rstats-on-nix)\n\n")

# Generate the Nix environment specification
rix(
  # Use date pinning for deep reproducibility with rstats-on-nix fork
  # This pins R version, CRAN packages, and Bioconductor to a specific snapshot
  # Using 2023-04-01 - the latest snapshot available in rstats-on-nix/nixpkgs
  # This ensures binary cache coverage (packages won't need to build from source)
  r_ver = "4.5.1",
  
  # R packages from CRAN and Bioconductor
  # Organized by category for maintainability
  r_pkgs = c(
    # === Core Infrastructure ===
    "Rcpp", "RcppArmadillo", "RcppEigen", "BH", 
    "data.table", "Matrix", "matrixStats",
    "foreach", "doParallel", "future",
    "remotes", "pak",
    
    # === Tidyverse Ecosystem ===
    "tidyverse",  # Includes dplyr, ggplot2, tidyr, purrr, readr, stringr, etc.
    "readxl",
    
    # === Statistical Packages ===
    "glmnet", "e1071", "ranger", "caret", "lme4",
    "mclust", "MCMCpack", "huge", "nnls", "moments", "pracma",
    
    # === Visualization ===
    "gridExtra", "gplots", "RColorBrewer", "plotly", "heatmaply",
    "ggpubr", "ggtext", "corrplot", "png", "colorspace", "scales",
    "R.devices", "DT", "kableExtra",
    
    # === Bioconductor Core ===
    "BiocManager",
    "SummarizedExperiment", "limma", "vsn",
    "AnnotationDbi", "biomaRt", "ExperimentHub",
    "ComplexHeatmap",
    "Biostrings", "IRanges", "S4Vectors", "zlibbioc",
    
    # === Batch Correction ===
    "sva",           # ComBat and other batch correction methods
    "batchelor",     # Bioconductor batch correction methods
    "RUVSeq",        # Remove Unwanted Variation
    # Note: BatchQC excluded due to NCmisc dependency issues
    # Note: preprocessCore will be added via git_pkgs if needed
    
    # === Machine Learning ===
    "MLmetrics", "ROCR", "xgboost", "lightgbm",
    "parsnip", "tidymodels",
    
    # === Dimensionality Reduction ===
    "Rtsne", "umap",
    
    # === Utilities ===
    "argparse", "docstring", "itertools", "pacman",
    "stringi", "magick", "bigmemory",
    
    # === Specialized Packages ===
    "Seurat",        # Single-cell analysis
    "fairadapt",     # Fairness-aware data preprocessing
    "seqgendiff",    # RNA-seq simulation (alternative to polyester)
    "logspline",     # For polyester dependency
    
    # === VS Code Integration ===
    "languageserver"  # REQUIRED for VS Code R extension
  ),
  
  # System dependencies (handled by Nix)
  # These are the underlying system libraries needed by R packages
  system_pkgs = c(
    # SSL and networking
    "openssl", "curl", "libxml2",
    
    # Graphics and fonts
    "fontconfig", "harfbuzz", "fribidi", "freetype",
    "libpng", "libtiff", "libjpeg", "libwebp",
    
    # Image processing
    "imagemagick",
    
    # Math and optimization
    "glpk", "gsl",
    
    # Geospatial (if needed)
    "udunits"
  ),
  
  # GitHub packages (for packages not in CRAN/Bioconductor)
  # Format: list(package_name = "user/repo@commit_or_tag")
  # Note: polyester is not in Bioconductor 3.21, may need GitHub installation
  # Note: preprocessCore may need GitHub installation
  # We'll add these in Task 8.5 if needed after testing
  git_pkgs = NULL,
  
  # HPC configuration
  # Use "none" for HPC/SLURM compatibility (not "code" or "rstudio")
  ide = "none",
  
  # Generate files in current directory
  project_path = ".",
  
  # Overwrite existing files
  overwrite = TRUE,
  
  # Print the generated default.nix for verification
  print = TRUE
)

cat("\n=== Generation complete ===\n")
cat("Files created:\n")
cat("  - default.nix: Nix expression defining the R environment\n")
cat("  - .Rprofile: R profile for library path isolation\n\n")
cat("Next steps:\n")
cat("  1. Review default.nix to verify package list\n")
cat("  2. Run Phase 2 to build the environment: ./run_generator.sh\n")
cat("  3. Test the environment: nix-shell --run 'R --version'\n")
