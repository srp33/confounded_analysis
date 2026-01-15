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


  batch_vec <- df$meta_source
  design <- model.matrix(~1, data=df)

  # FIX: Add back later with verbose flag
  # cat("  [apply_adjustment] Input matrix:", nrow(num_cols), "rows x", ncol(num_cols), "cols\n")
  # cat("  [apply_adjustment] Numeric colnames example:", paste(head(colnames(num_cols), 5), collapse = ", "), "\n")
  # cat("  [apply_adjustment] Value summary (pre-transform):\n")
  # print(summary(as.vector(num_cols_matrix)[1:min(1000, length(num_cols_matrix))]))
  
  # ---- Log-scale harmonization ----
  # Heuristic: if data looksk liek counts (large max, strong skew), apply log2(x + 1)

  # Check counts
  idx <- sample(length(num_cols_matrix),
                min(100000, length(num_cols_matrix)))
  num_vals <- num_cols_matrix[idx]
  num_vals <- num_vals[is.finite(num_vals)]


  looks_like_counts <- (
    max(num_vals, na.rm = TRUE) > 100 ||
    quantile(num_vals, 0.99, na.rm = TRUE) > 50
  )

  if (looks_like_counts) {
    cat(" [preprocess] Detected count-like data: applying log2(x+1\n)")
    num_cols_matrix <- log2(num_cols + 1)
  } else {
    cat(" [preprocess] Data appear already log-scaled; no log transform applied\n")
  }

  # ---- Gene-wise centering ----
  # Center each gene (row) across samples

  if (method != "log_transformed") { 
    gene_means <- rowMeans(num_cols_matrix, na.rm = TRUE)
    num_cols_matrix <- sweep(num_cols_matrix, 1, gene_means, FUN = "-")

    cat(" [preprocess] Applied gene-wise centering\n")
    cat("  [preprocess] Value summary after log+centering:\n")
    print(summary(as.vector(num_cols_matrix)[
      1:min(1000, length(as.vector(num_cols_matrix)))
    ]))
  }

  # HVG Selection for MNN
  if (method == "mnn") {

    cat(" [mnn] Selecting top 3000 HVGs from training data\n")

    # num_cols_matrix is genes x samples (already log-scaled and centered)
    train_idx <- which(df$meta_source != test_source)
    train_mat <- num_cols_matrix[, train_idx, drop = FALSE]

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
  if (is.null(rownames(num_cols_matrix))) {
    rownames(num_cols_matrix) <- paste0("gene_", seq_len(nrow(num_cols_matrix)))
  }
  if (is.null(colnames(num_cols_matrix))) {
    colnames(num_cols_matrix) <- colnames(num_cols)
  }

  num_cols_matrix_t <- t(num_cols_matrix)

  cat("Applying adjustment method:", method, "\n")
  adjusted <- switch(method,
    min_mean = adjust_min_mean(num_cols_matrix_t, batch = df$meta_source),
    log_combat = adjust_log_combat(num_cols_matrix_t, batch = batch_vec, design = design),
    mnn = adjust_mnn(df_=num_cols_matrix_t, batch=batch_vec, test_source=test_source, data_are_counts=FALSE, batch_levels=batch_levels, debug = FALSE),
    gmm = adjust_gmm(matrix_ = num_cols_matrix_t, batch = batch_vec, debug = FALSE),
    log_transformed = num_cols_matrix_t,
    stop("Unknown adjuster: ", method)
  )

  cat("  [apply_adjustment] Adjustment done.\n")
  cat("  [apply_adjustment] Adjusted matrix dimensions:", dim(adjusted), "\n")
  cat("  [apply_adjustment] Adjusted value summary:\n")

  # Random sampling from the full matrix for a value summary
  vals <- as.vector(num_cols_matrix)
  vals <- vals[is.finite(vals)]

  n_show <- min(10000, length(vals))
  sample_vals <- sample(vals, n_show)

  summary(sample_vals)


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
