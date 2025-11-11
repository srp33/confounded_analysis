# generate_env_minimal.R
# Phase 1 (Authoring): Generate default.nix and .Rprofile for combatseq environment
# MINIMAL VERSION: Core ComBat-seq packages only
# This is a focused subset for ComBat-seq workflows

library(rix)

cat("=== Generating MINIMAL combatseq R environment ===\n")
cat("This will create default.nix and .Rprofile files\n")
cat("Date pinning: 2024-10-15 (R 4.4.x, Bioconductor 3.21)\n")
cat("Strategy: Core ComBat-seq packages only\n\n")

# Generate the Nix environment specification
rix(
  # Use same date as batch-effects for consistency
  date = "2024-10-15",
  
  # MINIMAL R packages for ComBat-seq workflows
  r_pkgs = c(
    # === Bioconductor Core ===
    "BiocManager",
    "SummarizedExperiment",
    
    # === Batch Correction (ComBat-seq specific) ===
    "sva",           # ComBat and ComBat-seq
    "limma",         # Linear models for microarray/RNA-seq
    "vsn",           # Variance stabilization
    "ComplexHeatmap", # Visualization
    "RUVSeq",        # Remove Unwanted Variation
    "DESeq2",        # Differential expression
    "batchelor",     # Batch correction methods
    
    # === Core Infrastructure ===
    "data.table",
    "Matrix",
    
    # === Tidyverse (minimal) ===
    "dplyr",
    "ggplot2",
    
    # === Statistical Packages ===
    "caret",
    "glmnet",
    "e1071",
    
    # === Machine Learning Metrics ===
    "ROCR",
    
    # === Utilities ===
    "argparse",      # Command-line argument parsing
    
    # === VS Code Integration ===
    "languageserver"  # REQUIRED for VS Code R extension
  ),
  
  # System dependencies (minimal for ComBat-seq)
  system_pkgs = c(
    "openssl",
    "curl",
    "libxml2"
  ),
  
  # GitHub packages (if needed)
  # preprocessCore may be needed for some ComBat-seq workflows
  git_pkgs = list(
    list(
      package_name = "preprocessCore",
      repo_url = "https://github.com/bmbolstad/preprocessCore",
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
cat("  - default.nix: Nix expression with MINIMAL ComBat-seq packages\n")
cat("  - .Rprofile: R profile for library path isolation\n\n")
cat("IMPORTANT: Additional packages\n")
cat("  This environment includes only core ComBat-seq packages (~20 packages)\n")
cat("  Additional packages can be installed via R after activation if needed\n\n")
cat("Next steps:\n")
cat("  1. Review default.nix to verify package list\n")
cat("  2. Run Phase 2 to build the environment: ./build_with_cache.sh\n")
cat("  3. Test the environment: nix-shell --run 'R --version'\n")
