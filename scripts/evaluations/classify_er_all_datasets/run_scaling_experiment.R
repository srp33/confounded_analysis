# Usage:
#   Rscript run_scaling_experiment.R <adjuster> <subset_index>
# Example:
#   Rscript run_scaling_experiment.R log_combat 2

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ---- Parse Arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_scaling_experiment.R <adjuster> <subset_index>")
}

adjuster <- args[1]
subset_index <- as.integer(args[2])

cat("Running scaling experiment with adjuster:", adjuster, "on subset:", subset_index, "\n")

# ---- Source adjustment functions ----
source("/scripts/adjust/adjust.R")

# ---- Adjustment wrapper ----
apply_adjustment <- function(df, method, test_source) {
  meta_cols <- df %>% select(starts_with("meta_"))
  num_cols <- df %>% select(where(is.numeric))

  if (ncol(num_cols) == 0) stop("No numeric columns found in dataset.")

  if (method != "gmm") {
    if (any(num_cols < 0, na.rm = TRUE)) {
      warning("Negative values found in numeric columns; consider shifting data.")
    }
    # Safe log transform (shift by min if negatives exist)
    shift <- ifelse(any(num_cols < 0, na.rm = TRUE), abs(min(num_cols, na.rm = TRUE)) + 1, 0)
    num_cols <- log1p(num_cols + shift)
    cat("Applying log transform before ", method, "\n")
  }

  cat("Applying adjustment method:", method, "\n")
  adjusted <- switch(method,
                     min_mean = adjust_min_mean(num_cols),
                     log_combat = adjust_log_combat(num_cols, ref_batch = unique(df$meta_source)[1]),
                     mnn = adjust_mnn(num_cols, merge_last = test_source),
                     gmm = adjust_gmm(num_cols),
                     stop("Unknown adjuster: ", method)
  )

  # Recombine metadata safely
  adjusted_df <- bind_cols(meta_cols, adjusted)
  return(adjusted_df)
}

# ---- Process single subset ----
subset_path <- sprintf("data/all_combined_subsets/subset_%dstudies.csv", subset_index)
cat("\nProcessing:", subset_path, "\n")

if (!file.exists(subset_path)) {
  stop("Missing subset file:", subset_path)
}

df <- read_csv(subset_path, show_col_types = FALSE)
sources <- unique(df$meta_source)

for (test_source in sources) {
  cat("  -> Test source:", test_source, "\n")
  tryCatch({
    adjusted_df <- apply_adjustment(df, method = adjuster, test_source = test_source)

    out_dir <- file.path("data/adjusted_datasets", adjuster)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    out_path <- file.path(out_dir, sprintf("subset%dstudies_test_%s.csv", subset_index, test_source))
    write_csv(adjusted_df, out_path)
    cat("Saved adjusted dataset to:", out_path, "\n")
  }, error = function(e) {
    cat("⚠️  Error while processing subset", subset_index, "test source", test_source, ":", conditionMessage(e), "\n")
  })
}

cat("\n=== Finished subset", subset_index, "for adjuster:", adjuster, "===\n")
