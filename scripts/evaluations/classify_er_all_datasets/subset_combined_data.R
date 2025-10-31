# Prepping breast cancer study data
# Create new files that subset all_combined.tsv based on the meta_source column

# Libraries
library(readr)
library(dplyr)

# ---- Load combined data ----
# If the first column is gene names, set col_names = TRUE and then set rownames
combined <- read.csv("/data/all_combined_data/all_combined.csv", stringsAsFactors=FALSE) 

# Confirm that the dataset has a column called 'meta_source'
if (!"meta_source" %in% colnames(combined)) {
        stop("The file must have a 'meta_source' column with study identifiers.")
}

# Load dataset metadata and order from smallest sample size to largest
metadata <- read_csv("/scripts/evaluations/geo_metadata.csv")
if (!all(c("gse_id", "sample_size") %in% colnames(metadata))) {
        stop("Metadata file must have columns: gse_id and sample_size")
}

# Drop specific rows
metadata <- metadata %>%
        filter(!gse_id %in% c('gse115577', 'gse123845', 'gse163882', 'metabric'))
      
# Order by sample size (ascending)
metadata_ordered <- metadata %>% 
        arrange(sample_size)

# ---- Define study names ----
study_names <- metadata_ordered$gse_id
test_name <- study_names[1]
train_names <- study_names[-1]

# ---- Output directory ----
# Output folder
out_dir <- "/data/all_combined_subsets"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# ---- Generate subsets ----
set.seed(123) # ensuring reproducibility

for (k in 2:(length(study_names))) {
        selected_studies <- study_names[1:k]
        message("Creating subset with ", k, " studies: ", paste(selected_studies, collapse = ", "))

        subset_data <- combined %>% 
                filter(meta_source %in% selected_studies)

        out_file <- file.path(out_dir, paste0("subset_", k, "studies.csv"))

        write_csv(subset_data, out_file)
}

message("Done! Created one subset for each number of studies (2–", length(study_names), ") in:\n", out_dir)
