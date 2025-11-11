#!/usr/bin/env Rscript
# Test script for ComBat-seq environment

cat("=== Testing ComBat-seq R Environment ===\n")
cat("R version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n\n")

# Test core packages
test_packages <- c(
  "tidyverse", "dplyr", "ggplot2",
  "BiocManager", "sva", "BatchQC", "limma",
  "SummarizedExperiment", "batchelor", "DESeq2",
  "caret", "glmnet", "e1071"
)

cat("Testing package availability:\n")
success_count <- 0
fail_count <- 0

for (pkg in test_packages) {
  result <- tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    cat("  ✓", pkg, "\n")
    success_count <- success_count + 1
    TRUE
  }, error = function(e) {
    cat("  ✗", pkg, "- FAILED\n")
    fail_count <- fail_count + 1
    FALSE
  })
}

cat("\n=== Test Summary ===\n")
cat("Successful:", success_count, "/", length(test_packages), "\n")
cat("Failed:", fail_count, "/", length(test_packages), "\n")

if (fail_count == 0) {
  cat("\n✓ All core packages loaded successfully!\n")
  quit(status = 0)
} else {
  cat("\n✗ Some packages failed to load\n")
  quit(status = 1)
}
