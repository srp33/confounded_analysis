# Prepping breast cancer study data
# Create new files that subset all_combined.csv based on the meta_source column

# Libraries
library(readr)
library(dplyr)
library(argparse)

# Parse command-line arguments
parser <- ArgumentParser(description = "Create a subset of all_combined.csv using the top-K studies by sample count")

parser$add_argument('--input', required=True, help='Path to all_combined.csv')
parser$add_argument('--test',required=True, help='Test source name')
parser$add_argumet('--order', required=True, help='Randomized, fixed vector of the order to add datasets')
parser$add_argument('--k', required=True, type="integer", required=True, help='Number of studies to include (k)')
parser$add_argument('--output', type=Path, required=True, help='Output CSV path for the subset')
args <- parser$parse_args()

input_path <- args$input
test_source <- args$test
order_file <- args$order
k <- args$k 
output_path <- args$output

message(">>> Input combined file: ", input_path)
message(">>> Test source: ", test_source)
message(">>> Training order: ", dataset_order)
message(">>> Requested K = ", k)
message(">>> Output subset file: ", output_path)

# ---- Load combined data ----
# If the first column is gene names, set col_names = TRUE and then set rownames
combined <- read.csv(input_path, stringsAsFactors=FALSE) 

# Confirm that the dataset has a column called 'meta_source'
if (!"meta_source" %in% colnames(combined)) {
        stop("The file must have a 'meta_source' column with study identifiers.")
}

# Load dataset metadata
metadata_path <- "/scripts/evaluations/geo_metadata.csv"

if (!file.exists(metadata_path)) {
    stop("Metadata file not found at: ", metadata_path)
}

metadata <- read_csv(metadata_path, show_col_types = FALSE)

if (!all(c("gse_id", "sample_size") %in% colnames(metadata))) {
        stop("Metadata file must have columns: ", paste(required_cols, collapse=", "))
}

# Drop specific rows
metadata <- metadata %>%
        filter(!gse_id %in% c('gse115577', 'gse123845', 'gse163882'))

# Create order vector from the order file
order_vector <- readLines(order_file)

# Check that 
if (k < 1 || k > length(order_vector)) {
        stop("Requested k=", k, " is outside valid range: 2-", length(order_vector))
}

# ---- Create subset ----
selected_studies <- order_vector[0:k]

message(">>> Selected studies (k=", k, "): ", paste(selected_studies, collapse=", "))

subset_data <- combined %>%
        filter(meta_source %in% selected_studies)

# ---- Write output ----
out_dir <- dirname(output_path)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

write_csv(subset_data, output_path)

message(">>> Subset created successfully: ", output_path)