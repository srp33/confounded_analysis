#!/usr/bin/env Rscript

# BatchQC Analysis Script
# Usage: Rscript batchqc_analysis.R <input_csv> <output_dir>

suppressPackageStartupMessages({
  library(BatchQC)
  library(dplyr)
  library(readr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript batchqc_analysis.R <input_csv> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read data
cat("Reading data from:", input_file, "\n")
data <- read_csv(input_file, show_col_types = FALSE)

# Separate expression data from metadata
meta_cols <- grep("^meta_", colnames(data), value = TRUE)
expr_cols <- setdiff(colnames(data), meta_cols)

# Extract expression matrix (samples x genes) and transpose to (genes x samples)
expr_matrix <- as.matrix(data[, expr_cols])
rownames(expr_matrix) <- paste0("Sample_", 1:nrow(expr_matrix))
expr_matrix <- t(expr_matrix)

# Extract metadata
metadata <- data[, meta_cols, drop = FALSE]
rownames(metadata) <- colnames(expr_matrix)

# Clean up column names (remove meta_ prefix)
colnames(metadata) <- gsub("^meta_", "", colnames(metadata))

# Ensure batch column exists (use source as batch)
if ("source" %in% colnames(metadata)) {
  batch <- as.factor(metadata$source)
} else {
  stop("No 'meta_source' column found for batch information")
}

# Create condition variable if ER status is available
condition <- NULL
if ("er_status" %in% colnames(metadata)) {
  condition <- as.factor(metadata$er_status)
}

cat("Data dimensions:\n")
cat("  Expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
cat("  Batch levels:", length(levels(batch)), "-", paste(levels(batch), collapse = ", "), "\n")
if (!is.null(condition)) {
  cat("  Condition levels:", length(levels(condition)), "-", paste(levels(condition), collapse = ", "), "\n")
}

# Run BatchQC analysis
cat("Running BatchQC analysis...\n")

# Set output file path
output_file <- file.path(output_dir, "batchqc_report.html")

# Run BatchQC with condition if available, otherwise just batch
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

cat("BatchQC analysis complete!\n")
cat("Report saved to:", output_file, "\n")