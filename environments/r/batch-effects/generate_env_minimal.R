# generate_env_minimal.R
# Phase 1 (Authoring): Generate default.nix and .Rprofile for batch-effects environment
# MINIMAL VERSION: Core packages only, with GitHub sources for essential missing packages
# Additional packages can be installed via R's package manager after environment activation

library(rix)

cat("=== Generating MINIMAL batch-effects R environment ===\n")
cat("This will create default.nix and .Rprofile files\n")
cat("Date pinning: 2024-12-14 (R 4.4.x, Bioconductor 3.21)\n")
cat("Strategy: Core packages in Nix + additional packages via R\n\n")

# Generate the Nix environment specification
rix(
  # Use date pinning for deep reproducibility with rstats-on-nix fork
  # Using 2024-12-14 which has better binary cache coverage
  date = "2024-12-14",
  
  # MINIMAL R packages - only core essentials
  # Additional packages will be installed via R's package manager
  r_pkgs = c(
    # === Core Infrastructure (Essential) ===
    "Rcpp",          # C++ integration (many packages depend on this)
    "data.table",    # Fast data manipulation
    "Matrix",        # Sparse matrices
    "foreach",       # Parallel loops
    "doParallel",    # Parallel backend
    "remotes",       # Install from GitHub
    
    # === Tidyverse (Essential subset) ===
    "dplyr",         # Data manipulation
    "ggplot2",       # Plotting
    "tidyr",         # Data tidying
    "readr",         # Reading data
    "stringr",       # String manipulation
    
    # === Statistics (Core) ===
    "glmnet",        # Regularized regression
    "e1071",         # SVM and other ML
    "caret",         # ML framework
    
    # === Bioconductor (Batch Correction Core) ===
    "BiocManager",   # Bioconductor package manager
    "sva",           # ComBat and other batch correction
    "limma",         # Linear models for genomics
    "batchelor",     # Batch correction methods
    "RUVSeq",        # Remove Unwanted Variation
    "SummarizedExperiment",  # Bioconductor data structure
    "ComplexHeatmap", # Heatmap visualization
    
    # === Machine Learning (Core) ===
    "ROCR",          # ROC curves
    "xgboost",       # Gradient boosting
    
    # === Visualization (Core) ===
    "gridExtra",     # Multiple plots
    "RColorBrewer",  # Color palettes
    
    # === Utilities (Essential) ===
    "argparse",      # Command-line arguments
    
    # === VS Code Integration (Required) ===
    "languageserver" # REQUIRED for VS Code R extension
  ),
  
  # System dependencies (minimal)
  system_pkgs = c(
    "openssl",       # SSL/TLS
    "curl",          # HTTP client
    "libxml2",       # XML parsing
    "fontconfig",    # Font configuration
    "harfbuzz",      # Text shaping
    "fribidi"        # Bidirectional text
  ),
  
  # GitHub packages for essential missing packages
  # These are packages not available in rstats-on-nix that are critical for the pipeline
  git_pkgs = list(
    list(
      package_name = "preprocessCore",
      repo_url = "https://github.com/bmbolstad/preprocessCore",
      branch_name = "master",
      commit = "HEAD"  # Use latest; can pin to specific commit for reproducibility
    ),
    list(
      package_name = "polyester",
      repo_url = "https://github.com/alyssafrazee/polyester",
      branch_name = "master",
      commit = "HEAD"
    )
  ),
  
  # HPC configuration
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
cat("  - default.nix: Nix expression with MINIMAL core packages\n")
cat("  - .Rprofile: R profile for library path isolation\n\n")
cat("IMPORTANT: Additional packages\n")
cat("  This environment includes only core packages (~30 packages)\n")
cat("  Additional packages can be installed via R after activation:\n")
cat("    - Run: Rscript setup_packages.R\n")
cat("    - Or install manually: install.packages('package_name')\n\n")
cat("Next steps:\n")
cat("  1. Review default.nix to verify package list\n")
cat("  2. Run Phase 2 to build the environment: ./build_with_cache.sh\n")
cat("  3. Install additional packages: Rscript setup_packages.R\n")
cat("  4. Test the environment: nix-shell --run 'R --version'\n")
