# Usage:
#   Rscript run_scaling_experiment.R <adjuster> <subset_index>
# Example:
#   Rscript run_scaling_experiment.R log_combat 2

# ---- Load Libraries ----
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(argparse)
})

# ---- Parse Arguments ----
parser <- ArgumentParser(description = "Adjust each subset with a given adjuster and subset index k.")

parser$add_argument('--adjuster', required=TRUE, help='Desired adjuster (gmm, min_mean, combat, mnn, or log_transformed)')
parser$add_argument('--subset-path', required=TRUE, help='Subset .csv file to be adjusted.')
parser$add_argument('--k', required=TRUE, help='Number of datasets in the subset, k.')
parser$add_argument('--test', required=TRUE, help='Test source name/ID.')
parser$add_argument('--output-dir', required=TRUE, help='Output directory for adjusted datasets.')
parser$add_argument('--adjust-script', required=TRUE, help='Path to adjust.R script.')
parser$add_argument('--metadata-file', required=TRUE, help='GEO metadata CSV to create MNN order list.')

args <- parser$parse_args()

adjuster <- args$adjuster
subset_path <- args$subset_path
subset_index <- args$k
test_source <- args$test
output_dir <- args$output_dir
adjust_script <- args$adjust_script
metadata_file <- args$metadata_file

cat("Running scaling experiment with adjuster:", adjuster, "on file: ", subset_path, "\n")

# ---- Source adjustment functions ----
source(adjust_script)

# ---- Set random seed ----
set.seed(1)

# ---- Adjustment wrapper ----
apply_adjustment <- function(df, method, test_source, metadata_file) {
  meta_cols <- df %>% select(starts_with("meta_"))
  num_cols <- df %>% select(where(is.numeric), -starts_with("meta_"))
  if (ncol(num_cols) == 0) stop("No numeric columns found in dataset.")

  # Convert to matrix
  num_cols_matrix <- as.matrix(num_cols)
  num_cols_matrix <- t(num_cols_matrix) # Should be genes x samples

  # Check that there are more genes than samples
  stopifnot(nrow(num_cols_matrix) > ncol(num_cols_matrix))

  batch_vec <- df$meta_source
  design <- model.matrix(~1, data=df)

  # ---- Gene-wise centering ----
  # Only center genes if adjuster is NOT the baseline
  if (adjuster != "log_transformed") { 
      gene_means <- rowMeans(num_cols_matrix, na.rm = TRUE)
      num_cols_matrix <- sweep(num_cols_matrix, 1, gene_means, FUN = "-")
      cat(" [preprocess] Applied gene-wise centering\n")
  } else {
      cat(" [preprocess] Skipping gene-wise centering for log_transformed baseline\n")
  }

  cat("[debug] dim(num_cols_matrix):", dim(num_cols_matrix), "\n")

  # HVG Selection for MNN
  if (method == "mnn") {

    cat(" [mnn] Selecting top 3000 HVGs from training data\n")

    # num_cols_matrix is genes x samples (already log-scaled and centered)
    train_idx <- which(df$meta_source != test_source)
    train_mat <- num_cols_matrix[, train_idx, drop = FALSE]

    cat("[debug] length(train_idx):", length(train_idx), "\n")
    cat("[debug] ncol(train_mat):", ncol(train_mat), "\n")

    gene_vars <- rowMeans(train_mat^2)

    # Remove genes with NA / zero variance
    valid <- is.finite(gene_vars) & gene_vars > 0
    gene_vars <- gene_vars[valid]

    top_n <- min(3000, length(gene_vars))
    top_idx <- order(gene_vars, decreasing = TRUE)[seq_len(top_n)]
    hvg_genes <- names(gene_vars)[top_idx]

    # Subset full matrix (train + test) to HVGs
    num_cols_matrix <- num_cols_matrix[hvg_genes, , drop = FALSE]

    cat(" [mnn] Retained", nrow(num_cols_matrix), "genes for MNN\n")
  }

  # Create list of order from largest to smallest for MNN
  if (method == "mnn") {
    train_datasets <- df %>%
      filter(meta_source != test_source) %>%
      pull(meta_source) %>% 
      unique()

    geo_meta <- read_csv(metadata_file, col_types = cols()) %>%
      filter(gse_id %in% train_datasets) %>%
      arrange(desc(sample_size))
    
    cat("  [apply_adjustment] MNN order (train datasets by sample size:\n")

    batch_levels <- c(geo_meta$gse_id, test_source)
  }

  # Ensure numeric matrix has proper row and column names *after* log transform
  # genes × samples matrix
  if (is.null(rownames(num_cols_matrix))) {
    rownames(num_cols_matrix) <- colnames(num_cols)  # gene names
  }
  if (is.null(colnames(num_cols_matrix))) {
    colnames(num_cols_matrix) <- df$meta_Sample_ID   # or rownames(df)
  }


  cat("Applying adjustment method:", method, "\n")
  adjusted <- switch(method,
    min_mean = adjust_min_mean(num_cols_matrix, batch = df$meta_source),
    log_combat = adjust_log_combat(num_cols_matrix, batch = batch_vec, design = design),
    mnn = adjust_mnn(df_=num_cols_matrix, batch=batch_vec, test_source=test_source, data_are_counts=FALSE, batch_levels=batch_levels, debug = FALSE),
    gmm = adjust_gmm(matrix_ = num_cols_matrix, batch = batch_vec, debug = FALSE),
    log_transformed = num_cols_matrix,
    stop("Unknown adjuster: ", method)
  )

  cat("  [apply_adjustment] Adjustment done.\n")
  cat("  [apply_adjustment] Adjusted matrix dimensions:", dim(adjusted), "\n")
  cat("  [apply_adjustment] Adjusted value summary:\n")

  # Random sampling from the full matrix for a value summary
  vals <- as.vector(adjusted)
  vals <- vals[is.finite(vals)]

  n_show <- min(10000, length(vals))
  sample_vals <- sample(vals, n_show)

  print(summary(sample_vals))


  adjusted <- t(adjusted)

  # Recombine metadata safely
  colnames(adjusted) <- colnames(num_cols)
  adjusted_df <- bind_cols(meta_cols, adjusted)
  return(adjusted_df)
}

# ---- Process single subset ----
if (!file.exists(subset_path)) {
  stop("Missing subset file:", subset_path)
}

df <- read_csv(subset_path, show_col_types = FALSE)
cat("  Loaded subset with", nrow(df), "rows and", ncol(df), "columns.\n")
cat("  Example columns:", paste(head(colnames(df), 5), collapse = ", "), "\n")
cat("  Numeric columns:", sum(sapply(df, is.numeric)), "\n\n")


cat("\n--- Processing test source:", test_source, "---\n")
  tryCatch({
    adjusted_df <- apply_adjustment(df, method = adjuster, test_source = test_source, metadata_file=metadata_file)
    
    out_dir <- file.path(output_dir, adjuster)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    cat("  Preparing to write CSV: ", nrow(adjusted_df), "rows x", ncol(adjusted_df), "columns\n")
    cat("  Column names (first few):", paste(head(colnames(adjusted_df), 5), collapse = ", "), "\n")

    out_path <- file.path(out_dir, sprintf("%s_%sstudies_test_%s.csv", adjuster, subset_index, test_source))
    write_csv(adjusted_df, out_path)
    cat("Saved adjusted dataset to:", out_path, "\n")
  }, error = function(e) {
    cat("⚠️  Error while processing subset", subset_index, "test source", test_source, ":", conditionMessage(e), "\n")
  })

cat("\n=== Finished subset", subset_index, "for adjuster:", adjuster, "===\n")
