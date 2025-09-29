#!/usr/bin/env Rscript

# Test script to verify the GMM streamlined fix
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

source("scripts/adjust/gmm_adjust_streamlined.R")

# Create test data that might cause the original error
set.seed(42)

# Test case 1: Normal case
cat("Testing normal case...\n")
normal_data <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
colnames(normal_data) <- paste0("gene_", 1:10)
batch <- rep(c("batch1", "batch2"), each = 50)

result1 <- gmm_adjust_streamlined(normal_data, batch, debug = TRUE)
cat("Normal case: SUCCESS\n\n")

# Test case 2: Problematic case with very low variance genes
cat("Testing low variance case...\n")
low_var_data <- matrix(rnorm(100 * 10, sd = 0.001), nrow = 100, ncol = 10)
colnames(low_var_data) <- paste0("gene_", 1:10)

result2 <- gmm_adjust_streamlined(low_var_data, batch, debug = TRUE)
cat("Low variance case: SUCCESS\n\n")

# Test case 3: Constant genes
cat("Testing constant genes case...\n")
constant_data <- matrix(1, nrow = 100, ncol = 10)
colnames(constant_data) <- paste0("gene_", 1:10)

result3 <- gmm_adjust_streamlined(constant_data, batch, debug = TRUE)
cat("Constant genes case: SUCCESS\n\n")

# Test case 4: Mixed case with some problematic genes
cat("Testing mixed case...\n")
mixed_data <- cbind(
  matrix(rnorm(100 * 5), nrow = 100, ncol = 5),  # Normal genes
  matrix(rep(1, 100 * 2), nrow = 100, ncol = 2),  # Constant genes
  matrix(rnorm(100 * 3, sd = 0.001), nrow = 100, ncol = 3)  # Low variance genes
)
colnames(mixed_data) <- paste0("gene_", 1:10)

result4 <- gmm_adjust_streamlined(mixed_data, batch, debug = TRUE)
cat("Mixed case: SUCCESS\n\n")

cat("All tests completed successfully!\n")