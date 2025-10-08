#!/usr/bin/env Rscript

# Install BatchQC and dependencies

cat("Installing BatchQC and dependencies...\n")

# Check if BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cran.r-project.org")
}

library(BiocManager)

# Install BatchQC and its dependencies
packages_to_install <- c(
  "BatchQC",
  "sva",
  "limma", 
  "preprocessCore",
  "gplots",
  "RColorBrewer",
  "corrplot",
  "DT",
  "ggplot2",
  "plotly",
  "heatmaply",
  "moments",
  "matrixStats"
)

cat("Installing packages:", paste(packages_to_install, collapse = ", "), "\n")

for (pkg in packages_to_install) {
  cat("Installing", pkg, "...\n")
  tryCatch({
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    cat("✓", pkg, "installed successfully\n")
  }, error = function(e) {
    cat("✗ Failed to install", pkg, ":", e$message, "\n")
  })
}

# Test if BatchQC loads
cat("\nTesting BatchQC installation...\n")
tryCatch({
  library(BatchQC)
  cat("✓ BatchQC loaded successfully!\n")
  
  # Check if main function exists
  if (exists("batchQC")) {
    cat("✓ batchQC function is available\n")
  } else {
    cat("✗ batchQC function not found\n")
  }
  
}, error = function(e) {
  cat("✗ Failed to load BatchQC:", e$message, "\n")
  cat("Trying alternative approach...\n")
  
  # Try installing from GitHub as fallback
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cran.r-project.org")
  }
  
  tryCatch({
    devtools::install_github("mani2012/BatchQC")
    library(BatchQC)
    cat("✓ BatchQC installed from GitHub successfully!\n")
  }, error = function(e2) {
    cat("✗ GitHub installation also failed:", e2$message, "\n")
  })
})

cat("\nInstallation complete!\n")