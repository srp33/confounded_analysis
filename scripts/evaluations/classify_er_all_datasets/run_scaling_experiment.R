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

parser$add_argument('--adjuster', required=TRUE, help='Desired adjuster (gmm, min_mean, combat, mnn, or log_transform)')
parser$add_argument('--subset-path', required=TRUE, help='Subset .csv file to be adjusted.')
parser$add_argument('--k', required=TRUE, help='Number of datasets in the subset, k.')
parser$add_argument('--test', required=TRUE, help='Test source name/ID.')
parser$add_argument('--output-dir', required=TRUE, help='Output directory for adjusted datasets.')
parser$add_argument('--adjust-script', required=TRUE, help='Path to adjust.R script.')

args <- parser$parse_args()

adjuster <- args$adjuster
subset_path <- args$subset_path
subset_index <- args$k
test_source <- args$test
output_dir <- args$output_dir
adjust_script <- args$adjust_script

cat("Running scaling experiment with adjuster:", adjuster, "on file: ", subset_path, "\n")

# ---- Source adjustment functions ----
source(adjust_script)

# ---- Adjustment wrapper ----
apply_adjustment <- function(df, method, test_source) {
  meta_cols <- df %>% select(starts_with("meta_"))
  num_cols <- df %>% select(where(is.numeric), -starts_with("meta_"))
  if (ncol(num_cols) == 0) stop("No numeric columns found in dataset.")

  batch_vec <- df$meta_source
  design <- model.matrix(~1, data=df)

  cat("  [apply_adjustment] Input matrix:", nrow(num_cols), "rows x", ncol(num_cols), "cols\n")
  cat("  [apply_adjustment] Numeric colnames example:", paste(head(colnames(num_cols), 5), collapse = ", "), "\n")
  cat("  [apply_adjustment] Value summary (pre-transform):\n")
  print(summary(as.vector(as.matrix(num_cols))[1:min(1000, length(as.matrix(num_cols)))]))
  
  if (method != "gmm" ) {
    if (any(num_cols < 0, na.rm = TRUE)) {
      warning("Negative values found in numeric columns; consider shifting data.")
    }
    # Safe log transform (shift by min if negatives exist)
    shift <- ifelse(any(num_cols < 0, na.rm = TRUE), abs(min(num_cols, na.rm = TRUE)) + 1, 0)
    num_cols <- log1p(num_cols + shift)
    cat("Applying log transform before ", method, "\n")
  }

  # Ensure numeric matrix has proper row and column names *after* log transform
  num_cols_matrix <- as.matrix(num_cols)
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
    mnn = adjust_mnn(df_ = num_cols_matrix_t, batch = batch_vec, test_source=test_source, data_are_counts = FALSE, debug = FALSE),
    gmm = adjust_gmm(matrix_ = num_cols_matrix_t, batch = batch_vec, debug = FALSE),
    log_transformed = num_cols_matrix_t,
    stop("Unknown adjuster: ", method)
  )

  cat("  [apply_adjustment] Adjustment done.\n")
  cat("  [apply_adjustment] Adjusted matrix dimensions:", dim(adjusted), "\n")
  cat("  [apply_adjustment] Adjusted value summary:\n")
  print(summary(as.vector(adjusted)[1:min(1000, length(as.vector(adjusted)))]))

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
    adjusted_df <- apply_adjustment(df, method = adjuster, test_source = test_source)
    
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
