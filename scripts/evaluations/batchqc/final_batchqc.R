#!/usr/bin/env Rscript

# Final working BatchQC analysis script

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
output_dir <- "/outputs/batchqc_final"

if (!file.exists(test_file)) {
  stop("Test file not found: ", test_file)
}

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Final BatchQC analysis with:", test_file, "\n")

# Read data and ensure we get samples from both batches
data <- read_csv(test_file, show_col_types = FALSE)
cat("Original data dimensions:", nrow(data), "samples x", ncol(data), "columns\n")

# Get metadata columns
meta_cols <- grep("^meta_", colnames(data), value = TRUE)
expr_cols <- setdiff(colnames(data), meta_cols)

# Sample data to ensure we have both batches
batch_col <- data$meta_source
batch_counts <- table(batch_col)
cat("Available batches:", paste(names(batch_counts), "=", batch_counts, collapse = ", "), "\n")

# Take 100 samples from each batch for a more robust analysis
samples_per_batch <- 100
selected_indices <- c()

for (batch_name in names(batch_counts)) {
  batch_indices <- which(batch_col == batch_name)
  selected_batch_indices <- head(batch_indices, samples_per_batch)
  selected_indices <- c(selected_indices, selected_batch_indices)
}

# Take first 500 genes for testing
test_genes <- head(expr_cols, 500)

# Subset data
data_subset <- data[selected_indices, c(test_genes, meta_cols)]
cat("Balanced subset dimensions:", nrow(data_subset), "samples x", ncol(data_subset), "columns\n")

# Check batch distribution in subset
subset_batch_counts <- table(data_subset$meta_source)
cat("Subset batch distribution:", paste(names(subset_batch_counts), "=", subset_batch_counts, collapse = ", "), "\n")

# Prepare expression matrix
expr_matrix <- as.matrix(data_subset[, test_genes])
rownames(expr_matrix) <- paste0("Sample_", 1:nrow(expr_matrix))
expr_matrix <- t(expr_matrix)  # Genes as rows, samples as columns

# Extract and clean metadata
metadata <- data_subset[, meta_cols, drop = FALSE]
rownames(metadata) <- colnames(expr_matrix)
colnames(metadata) <- gsub("^meta_", "", colnames(metadata))

# Set up batch and condition
batch <- as.factor(metadata$source)
condition <- if ("er_status" %in% colnames(metadata)) as.factor(metadata$er_status) else NULL

cat("Final batch distribution:", paste(names(table(batch)), "=", table(batch), collapse = ", "), "\n")
if (!is.null(condition)) {
  cat("Condition distribution:", paste(names(table(condition)), "=", table(condition), collapse = ", "), "\n")
}

# Create SummarizedExperiment
cat("Creating SummarizedExperiment object...\n")

clean_metadata <- DataFrame(
  batch = batch,
  row.names = colnames(expr_matrix)
)

if (!is.null(condition)) {
  clean_metadata$condition <- condition
}

se_data <- SummarizedExperiment(
  assays = list(counts = expr_matrix),
  colData = clean_metadata
)

cat("SummarizedExperiment created successfully!\n")
cat("Dimensions:", nrow(se_data), "genes x", ncol(se_data), "samples\n")

# Run BatchQC analysis
cat("\n=== Running BatchQC Analysis ===\n")

# 1. PCA Analysis
cat("1. Creating PCA plot...\n")
tryCatch({
  pca_data <- t(assay(se_data, "counts"))
  pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
  
  # Create manual PCA plot
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    batch = colData(se_data)$batch,
    condition = if(!is.null(condition)) colData(se_data)$condition else "None"
  )
  
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = batch, shape = condition)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = "PCA Plot - Batch Effects Visualization", 
         x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "% variance)"),
         y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "% variance)")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(output_dir, "pca_plot.png"), pca_plot, width = 10, height = 8)
  cat("✓ PCA plot created successfully\n")
  
}, error = function(e) {
  cat("✗ PCA analysis failed:", e$message, "\n")
})

# 2. Explained Variation Analysis
cat("2. Calculating explained variation...\n")
tryCatch({
  ev_result <- batchqc_explained_variation(
    se = se_data, 
    batch = "batch", 
    condition = if(!is.null(condition)) "condition" else NULL,
    assay_name = "counts"
  )
  
  ev_plot <- EV_plotter(ev_result)
  ggsave(file.path(output_dir, "explained_variation.png"), ev_plot, width = 10, height = 6)
  cat("✓ Explained variation analysis successful\n")
  
}, error = function(e) {
  cat("✗ Explained variation failed:", e$message, "\n")
})

# 3. Alternative heatmap approaches
cat("3. Creating correlation heatmap...\n")
tryCatch({
  # Create sample correlation heatmap manually
  library(pheatmap)
  
  # Calculate sample correlations
  expr_data <- assay(se_data, "counts")
  sample_cor <- cor(expr_data, method = "pearson")
  
  # Create annotation for heatmap
  annotation_df <- data.frame(
    Batch = colData(se_data)$batch,
    row.names = colnames(expr_data)
  )
  
  if (!is.null(condition)) {
    annotation_df$Condition <- colData(se_data)$condition
  }
  
  # Create heatmap
  png(file.path(output_dir, "sample_correlation_heatmap.png"), width = 1200, height = 1000)
  pheatmap(sample_cor, 
           annotation_col = annotation_df,
           annotation_row = annotation_df,
           main = "Sample Correlation Heatmap",
           show_rownames = FALSE,
           show_colnames = FALSE)
  dev.off()
  
  cat("✓ Sample correlation heatmap created successfully\n")
  
}, error = function(e) {
  cat("✗ Correlation heatmap failed:", e$message, "\n")
})

# 4. Boxplot for batch effects
cat("4. Creating batch effect boxplots...\n")
tryCatch({
  # Create boxplot of first PC by batch
  pca_data <- t(assay(se_data, "counts"))
  pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
  
  boxplot_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    batch = colData(se_data)$batch
  )
  
  pc1_boxplot <- ggplot(boxplot_df, aes(x = batch, y = PC1, fill = batch)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = "PC1 Distribution by Batch", x = "Batch", y = "PC1 Score") +
    theme_minimal()
  
  pc2_boxplot <- ggplot(boxplot_df, aes(x = batch, y = PC2, fill = batch)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = "PC2 Distribution by Batch", x = "Batch", y = "PC2 Score") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "pc1_boxplot.png"), pc1_boxplot, width = 8, height = 6)
  ggsave(file.path(output_dir, "pc2_boxplot.png"), pc2_boxplot, width = 8, height = 6)
  
  cat("✓ Batch effect boxplots created successfully\n")
  
}, error = function(e) {
  cat("✗ Boxplot creation failed:", e$message, "\n")
})

# 5. Summary statistics
cat("5. Generating comprehensive summary...\n")

# Calculate some basic statistics
batch_table <- table(colData(se_data)$batch)
condition_table <- if(!is.null(condition)) table(colData(se_data)$condition) else NULL

summary_text <- paste0(
  "BatchQC Analysis Summary\n",
  "========================\n",
  "Dataset: ", basename(test_file), "\n",
  "Analysis Date: ", Sys.Date(), "\n",
  "Samples: ", ncol(se_data), "\n",
  "Genes: ", nrow(se_data), "\n",
  "Batches: ", length(levels(batch)), " (", paste(names(batch_table), "=", batch_table, collapse = ", "), ")\n",
  if(!is.null(condition_table)) paste0("Conditions: ", length(levels(condition)), " (", paste(names(condition_table), "=", condition_table, collapse = ", "), ")\n") else "",
  "Output directory: ", output_dir, "\n",
  "\nGenerated Files:\n",
  "- pca_plot.png: PCA visualization showing batch separation\n",
  "- explained_variation.png: Variance explained by batch vs condition\n",
  "- sample_correlation_heatmap.png: Sample-to-sample correlation matrix\n",
  "- pc1_boxplot.png: PC1 distribution by batch\n",
  "- pc2_boxplot.png: PC2 distribution by batch\n",
  "- summary.txt: This summary file\n"
)

writeLines(summary_text, file.path(output_dir, "summary.txt"))
cat(summary_text)

cat("\n✓ Complete BatchQC analysis finished successfully!\n")
cat("All output files saved to:", output_dir, "\n")