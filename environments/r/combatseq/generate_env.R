# generate_env.R
# Phase 1 (Authoring): Generate default.nix and .Rprofile for combatseq environment
# This script uses rix to create a reproducible R environment specification
# for ComBat-seq workflows (Bioconductor 3.11 compatibility)

library(rix)

cat("=== Generating combatseq R environment ===\n")
cat("This will create default.nix and .Rprofile files\n")
cat("Date pinning: 2024-12-14 (R 4.4.x, Bioconductor 3.21)\n\n")
cat("Note: Using same R/Bioc version as batch-effects environment\n")
cat("      ComBat-seq methods are stable across Bioconductor versions\n")
cat("      This is a focused subset for ComBat-seq workflows\n\n")

# Generate the Nix environment specification
rix(
  # Use same date as batch-effects for consistency
  # ComBat-seq methods are stable across Bioconductor versions
  date = "2024-12-14",
  
  # R packages for ComBat-seq workflows
  # Focused on batch correction and differential expression
  r_pkgs = c(
    # === Bioconductor Core (Bioc 3.11) ===
    "BiocManager",
    "SummarizedExperiment",
    "ExperimentHub",
    
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
    
    # === Tidyverse (for data manipulation) ===
    "tidyverse",     # Includes dplyr, ggplot2, tidyr, etc.
    "dplyr",
    "ggplot2",
    
    # === Statistical Packages ===
    "caret",
    "glmnet",
    "e1071",
    
    # === Machine Learning Metrics ===
    "MLmetrics",
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
