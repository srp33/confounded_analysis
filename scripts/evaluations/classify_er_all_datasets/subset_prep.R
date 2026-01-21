# Prepping breast cancer study data
# Create new files that subset all_combined.csv based on the meta_source column

# ---- Load Libraries ----
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(argparse)
})

# ---- Parse Arguments ----
parser <- ArgumentParser(description = "Create a subset of all_combined.csv using the top-K studies by sample count")

parser$add_argument('--input', required=TRUE, help='Path to all_combined.csv')
parser$add_argument('--test',required=TRUE, help='Test source name')
parser$add_argument('--order', required=TRUE, help='Randomized, fixed vector of the order to add datasets')
parser$add_argument('--k', required=TRUE, type='integer', help='Number of studies to include (k)')
parser$add_argument('--output', required=TRUE, help='Output CSV path for the subset')
args <- parser$parse_args()

input_path <- args$input
test_source <- args$test
order_file <- args$order
k <- args$k 
output_path <- args$output

message(">>> Input combined file: ", input_path)
message(">>> Test source: ", test_source)
message(">>> Order file: ", order_file)
message(">>> Requested K = ", k)
message(">>> Output subset file: ", output_path)

# ---- Load combined data ----
# If the first column is gene names, set col_names = TRUE and then set rownames
combined <- read.csv(input_path, stringsAsFactors=FALSE) 

# Confirm that the dataset has a column called 'meta_source'
if (!"meta_source" %in% colnames(combined)) {
        stop("The file must have a 'meta_source' column with study identifiers.")
}

# Create order vector from the order file
order_df <- read_csv(order_file, col_types = cols())
order_vector <- order_df$train_source

# Check that 
if (k < 1 || k > length(order_vector)) {
        stop("Requested k=", k, " is outside valid range: 2-", length(order_vector))
}

# ---- Create subset ----
selected_studies <- unique(c(test_source, order_vector[1:k]))

message(">>> Selected studies (k=", k, "): ", paste(selected_studies, collapse=", "))

subset_data <- combined %>%
        filter(meta_source %in% selected_studies)

# Check for studies with no rows
missing_studies <- setdiff(selected_studies, unique(subset_data$meta_source))
if (length(missing_studies) > 0) {
        warning("The following studies had no rows in combined data: ", paste(missing_studies, collapse=", "))
}

# ---- Preprocessing: per-dataset log transform ----

# Identify numeric columns (gene expression)
num_cols <- subset_data %>% select(where(is.numeric), -starts_with("meta_"))
meta_cols <- subset_data %>% select(starts_with("meta_"))

if(ncol(num_cols) == 0) stop("No numeric columns found in subset.")

# Convert to matrix
num_mat <- as.matrix(num_cols)

# Apply per-dataset log transform
for(ds in unique(subset_data$meta_source)) {
    idx <- which(subset_data$meta_source == ds)
    mat_ds <- num_mat[idx, , drop = FALSE]

    # Decide if log-transform is needed
    # Simple heuristic: RNA-seq counts are non-negative and have high max/quantile
    if (all(mat_ds >= 0) && (max(mat_ds, na.rm=TRUE) > 100 || quantile(mat_ds, 0.99, na.rm=TRUE) > 50)) {
        message(">>> Applying log1p to dataset: ", ds)
        mat_ds <- log1p(mat_ds)
    } else {
        message(">>> Skipping log transform for dataset: ", ds)
    }

    num_mat[idx, ] <- mat_ds
}

# Recombine metadata with numeric matrix
subset_data_processed <- cbind(meta_cols, as.data.frame(num_mat))


# ---- Write output ----
out_dir <- dirname(output_path)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

write_csv(subset_data_processed, output_path)

message(">>> Subset processed and created successfully: ", output_path)