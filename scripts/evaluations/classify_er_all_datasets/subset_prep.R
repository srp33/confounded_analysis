library(readr)
library(dplyr)
library(argparse)

# -------------------------
# 1. Load combined data
# -------------------------
load_combined <- function(input_path) {
  df <- read.csv(input_path, stringsAsFactors = FALSE)
  if (!"meta_source" %in% colnames(df)) stop("Missing meta_source column")
  return(df)
}

# -------------------------
# 2. Create study subset
# -------------------------
create_subset <- function(df, test_source, order_vector, k) {
  if (k < 1 || k > length(order_vector)) stop("Invalid k")
  selected_studies <- unique(c(test_source, order_vector[1:k]))
  subset_df <- df %>% filter(meta_source %in% selected_studies)
  
  missing_studies <- setdiff(selected_studies, unique(subset_df$meta_source))
  if (length(missing_studies) > 0) warning("Missing studies: ", paste(missing_studies, collapse=", "))
  
  return(subset_df)
}

# -------------------------
# 3. Per-dataset log transform
# -------------------------
log_transform_per_dataset <- function(df) {
  num_cols <- df %>% select(where(is.numeric), -starts_with("meta_"))
  meta_cols <- df %>% select(starts_with("meta_"))
  num_mat <- as.matrix(num_cols)
  
  for(ds in unique(df$meta_source)) {
    idx <- which(df$meta_source == ds)
    mat_ds <- num_mat[idx, , drop = FALSE]
    
    if (all(mat_ds >= 0) && (max(mat_ds, na.rm=TRUE) > 100 || quantile(mat_ds, 0.99, na.rm=TRUE) > 50)) {
      message(">>> Applying log1p to dataset: ", ds)
      mat_ds <- log1p(mat_ds)
    } else {
      message(">>> Skipping log transform for dataset: ", ds)
    }
    
    num_mat[idx, ] <- mat_ds
  }
  
  return(cbind(meta_cols, as.data.frame(num_mat)))
}

# -------------------------
# 4. Write output
# -------------------------
write_subset <- function(df, output_path) {
  dir.create(dirname(output_path), recursive=TRUE, showWarnings=FALSE)
  write_csv(df, output_path)
  message(">>> Subset written to: ", output_path)
}

# -------------------------
# Main function
# -------------------------
main <- function() {
  parser <- ArgumentParser(description = "Create study subset with preprocessing")
  parser$add_argument('--input', required=TRUE)
  parser$add_argument('--test', required=TRUE)
  parser$add_argument('--order', required=TRUE)
  parser$add_argument('--k', required=TRUE, type='integer')
  parser$add_argument('--output', required=TRUE)
  args <- parser$parse_args()
  
  df <- load_combined(args$input)
  order_vector <- read_csv(args$order, col_types = cols())$train_source
  subset <- create_subset(df, args$test, order_vector, args$k)
  processed <- log_transform_per_dataset(subset)
  write_subset(processed, args$output)
}

if (!interactive()) main()
