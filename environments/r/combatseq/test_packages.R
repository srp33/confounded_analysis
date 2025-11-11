#!/usr/bin/env Rscript
# test_packages.R
# Test script to verify all ComBat-seq packages load correctly

cat("=== Testing ComBat-seq R Environment ===\n\n")

# Track success/failure
failed_packages <- c()
success_count <- 0
total_count <- 0

# Function to test package loading
test_package <- function(pkg_name) {
  total_count <<- total_count + 1
  cat(sprintf("Testing %s... ", pkg_name))
  
  result <- tryCatch({
    suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
    cat("✓ OK\n")
    success_count <<- success_count + 1
    TRUE
  }, error = function(e) {
    cat(sprintf("✗ FAILED: %s\n", e$message))
    failed_packages <<- c(failed_packages, pkg_name)
    FALSE
  })
  
  return(result)
}

cat("--- Bioconductor Core Packages ---\n")
test_package("BiocManager")
test_package("SummarizedExperiment")
test_package("ExperimentHub")

cat("\n--- Batch Correction Packages (ComBat-seq specific) ---\n")
test_package("sva")
test_package("limma")
test_package("vsn")
test_package("ComplexHeatmap")
test_package("RUVSeq")
test_package("DESeq2")
test_package("batchelor")

cat("\n--- Core Infrastructure ---\n")
test_package("data.table")
test_package("Matrix")

cat("\n--- Tidyverse ---\n")
test_package("tidyverse")
test_package("dplyr")
test_package("ggplot2")

cat("\n--- Statistical Packages ---\n")
test_package("caret")
test_package("glmnet")
test_package("e1071")

cat("\n--- Machine Learning Metrics ---\n")
test_package("MLmetrics")
test_package("ROCR")

cat("\n--- Utilities ---\n")
test_package("argparse")

cat("\n--- VS Code Integration ---\n")
test_package("languageserver")

# Summary
cat("\n=== Test Summary ===\n")
cat(sprintf("Total packages tested: %d\n", total_count))
cat(sprintf("Successful: %d\n", success_count))
cat(sprintf("Failed: %d\n", length(failed_packages)))

if (length(failed_packages) > 0) {
  cat("\nFailed packages:\n")
  for (pkg in failed_packages) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\nEnvironment test FAILED\n")
  quit(status = 1)
} else {
  cat("\n✓ All packages loaded successfully!\n")
  cat("Environment test PASSED\n")
  quit(status = 0)
}
