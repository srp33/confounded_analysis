#!/usr/bin/env Rscript

# Test BatchQC setup with a small dataset

suppressPackageStartupMessages({
  library(BatchQC)
  library(dplyr)
  library(readr)
})

# Test with one dataset
test_file <- "grp_batch_effects/data/paired_datasets/gse115577_gse58644/unadjusted.csv"
output_dir <- "grp_batch_effects/outputs/batchqc_test"

if (!file.exists(test_file)) {
  stop("Test file not found: ", test_file)
}

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Testing BatchQC with:", test_file, "\n")

# Read a subset of the data for testing (first 1000 genes)
data <- read_csv(test_file, show_col_types = FALSE)
cat("Original data dimensions:", nrow(data), "samples x", ncol(data), "columns\n")

# Take first 1000 genes for testing
meta_cols <- grep("^meta_", colnames(data), value = TRUE)
expr_cols <- setdiff(colnames(data), meta_cols)
test_genes <- head(expr_cols, 1000)

data_subset <- data[, c(test_genes, meta_cols)]
cat("Test subset dimensions:", nrow(data_subset), "samples x", ncol(data_subset), "columns\n")

# Prepare data for BatchQC
expr_matrix <- as.matrix(data_subset[, test_genes])
rownames(expr_matrix) <- paste0("Sample_", 1:nrow(expr_matrix))
expr_matrix <- t(expr_matrix)

# Extract metadata
metadata <- data_subset[, meta_cols, drop = FALSE]
rownames(metadata) <- colnames(expr_matrix)
colnames(metadata) <- gsub("^meta_", "", colnames(metadata))

# Set up batch and condition
batch <- as.factor(metadata$source)
condition <- if ("er_status" %in% colnames(metadata)) as.factor(metadata$er_status) else NULL

cat("Batch levels:", paste(levels(batch), collapse = ", "), "\n")
if (!is.null(condition)) {
  cat("Condition levels:", paste(levels(condition), collapse = ", "), "\n")
}

# Run BatchQC
cat("Running BatchQC test...\n")
output_file <- file.path(output_dir, "test_batchqc_report.html")

tryCatch({
  if (!is.null(condition)) {
    batchQC(expr_matrix, 
            batch = batch, 
            condition = condition,
            report_file = output_file,
            report_dir = output_dir)
  } else {
    batchQC(expr_matrix, 
            batch = batch,
            report_file = output_file,
            report_dir = output_dir)
  }
  cat("✓ BatchQC test successful!\n")
  cat("Report saved to:", output_file, "\n")
}, error = function(e) {
  cat("✗ BatchQC test failed:\n")
  cat("Error:", e$message, "\n")
})