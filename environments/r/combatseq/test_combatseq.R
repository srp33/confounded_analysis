#!/usr/bin/env Rscript
# test_combatseq.R
# Test ComBat-seq functionality with a simple example

cat("=== Testing ComBat-seq Functionality ===\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(sva)
  library(limma)
  library(DESeq2)
})

cat("--- Environment Information ---\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("sva version: %s\n", packageVersion("sva")))
cat(sprintf("limma version: %s\n", packageVersion("limma")))
cat(sprintf("DESeq2 version: %s\n", packageVersion("DESeq2")))

# Check if BiocManager is available and get Bioconductor version
if (requireNamespace("BiocManager", quietly = TRUE)) {
  bioc_version <- BiocManager::version()
  cat(sprintf("Bioconductor version: %s\n", bioc_version))
} else {
  cat("BiocManager not available\n")
}

cat("\n--- Testing ComBat-seq Function ---\n")

# Create simple test data
set.seed(42)
n_genes <- 100
n_samples <- 20

# Simulate count data
counts <- matrix(
  rpois(n_genes * n_samples, lambda = 100),
  nrow = n_genes,
  ncol = n_samples
)
rownames(counts) <- paste0("Gene", 1:n_genes)
colnames(counts) <- paste0("Sample", 1:n_samples)

# Create batch variable (2 batches)
batch <- rep(c(1, 2), each = n_samples / 2)

# Create biological covariate
bio_group <- rep(c("A", "B"), times = n_samples / 2)

cat(sprintf("Test data: %d genes x %d samples\n", n_genes, n_samples))
cat(sprintf("Batches: %s\n", paste(unique(batch), collapse = ", ")))
cat(sprintf("Biological groups: %s\n", paste(unique(bio_group), collapse = ", ")))

# Test ComBat_seq function
cat("\nRunning ComBat_seq...\n")
result <- tryCatch({
  adjusted_counts <- ComBat_seq(
    counts = counts,
    batch = batch,
    group = bio_group
  )
  cat("✓ ComBat_seq completed successfully\n")
  cat(sprintf("Output dimensions: %d genes x %d samples\n", 
              nrow(adjusted_counts), ncol(adjusted_counts)))
  
  # Verify output
  if (all(dim(adjusted_counts) == dim(counts))) {
    cat("✓ Output dimensions match input\n")
  } else {
    cat("✗ Output dimensions do not match input\n")
    stop("Dimension mismatch")
  }
  
  if (all(rownames(adjusted_counts) == rownames(counts))) {
    cat("✓ Gene names preserved\n")
  } else {
    cat("✗ Gene names not preserved\n")
    stop("Gene name mismatch")
  }
  
  TRUE
}, error = function(e) {
  cat(sprintf("✗ ComBat_seq failed: %s\n", e$message))
  FALSE
})

# Summary
cat("\n=== Test Summary ===\n")
if (result) {
  cat("✓ ComBat-seq functionality test PASSED\n")
  cat("Environment is ready for ComBat-seq workflows\n")
  quit(status = 0)
} else {
  cat("✗ ComBat-seq functionality test FAILED\n")
  quit(status = 1)
}
