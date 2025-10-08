#!/usr/bin/env Rscript

# Test BatchQC setup with a small dataset using individual analysis functions

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(SummarizedExperiment)
  library(S4Vectors)
})

library(BatchQC)

# Test with one dataset
test_file <- "/data/paired_datasets/gse115577_gse58644/unadjusted.csv"
output_dir <- "/outputs/batchqc_test"

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

# Create clean metadata DataFrame for SummarizedExperiment
cat("Creating SummarizedExperiment object...\n")

# Create a clean metadata data frame with only the columns we need
clean_metadata <- data.frame(
  batch = batch,
  stringsAsFactors = FALSE
)

# Add condition if available
if (!is.null(condition)) {
  clean_metadata$condition <- condition
}

# Add sample names
rownames(clean_metadata) <- colnames(expr_matrix)

# Create the SummarizedExperiment with clean metadata
se_data <- SummarizedExperiment(
  assays = list(counts = expr_matrix),
  colData = DataFrame(clean_metadata)
)

cat("SummarizedExperiment created successfully!\n")
cat("Dimensions:", nrow(se_data), "genes x", ncol(se_data), "samples\n")

cat("Batch levels:", paste(levels(batch), collapse = ", "), "\n")
if (!is.null(condition)) {
  cat("Condition levels:", paste(levels(condition), collapse = ", "), "\n")
}

# Run BatchQC analysis using individual functions
cat("Running BatchQC analysis...\n")

# Run BatchQC analysis
cat("Running BatchQC analysis...\n")

tryCatch({
  # 1. PCA Analysis
  cat("Creating PCA plot...\n")
  pca_plot <- PCA_plotter(
    se = se_data, 
    nfeature = 500,  # Number of features to use for PCA
    batch = "batch", 
    color = "batch", 
    assays = "counts",
    xaxisPC = 1,
    yaxisPC = 2
  )
  ggsave(file.path(output_dir, "pca_plot.png"), pca_plot, width = 10, height = 8)
  
  # 2. Heatmap
  cat("Creating heatmap...\n")
  heatmap_plot <- heatmap_plotter(se_data, batch = "batch", assays = "counts")
  png(file.path(output_dir, "heatmap.png"), width = 1200, height = 1000)
  print(heatmap_plot)
  dev.off()
  
  # 3. Explained Variation Analysis
  cat("Calculating explained variation...\n")
  ev_result <- batchqc_explained_variation(se_data, batch = "batch", condition = if(!is.null(condition)) "condition" else NULL)
  ev_plot <- EV_plotter(ev_result)
  ggsave(file.path(output_dir, "explained_variation.png"), ev_plot, width = 10, height = 6)
  
  # 4. Summary statistics
  cat("Generating summary...\n")
  summary_text <- paste0(
    "BatchQC Analysis Summary\n",
    "========================\n",
    "Dataset: ", basename(test_file), "\n",
    "Samples: ", ncol(se_data), "\n",
    "Genes: ", nrow(se_data), "\n",
    "Batches: ", length(levels(batch)), " (", paste(levels(batch), collapse = ", "), ")\n",
    "Conditions: ", if(!is.null(condition)) length(levels(condition)) else "None", "\n",
    "Output directory: ", output_dir, "\n"
  )
  
  writeLines(summary_text, file.path(output_dir, "summary.txt"))
  cat(summary_text)
  
  cat("✓ BatchQC analysis successful!\n")
  cat("Output files saved to:", output_dir, "\n")
  
}, error = function(e) {
  cat("✗ BatchQC analysis failed:\n")
  cat("Error:", e$message, "\n")
  cat("Traceback:\n")
  print(traceback())
})