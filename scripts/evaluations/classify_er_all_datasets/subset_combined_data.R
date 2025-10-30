# Prepping breast cancer study data
# Create new files that subset all_combined.tsv based on the meta_source column

# Libraries
library(readr)
library(dplyr)

# ---- Load combined data ----
# If the first column is gene names, set col_names = TRUE and then set rownames
combined <- read.csv("~/confounded_analysis/grp_batch_effects/data/all_combined_data/all_combined.csv", stringsAsFactors=FALSE) 

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

# ---- Output directory ----
out_dir <- "/data/all_combined_subsets"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# ---- Training set size tiers ----
training_tiers <- list(
        small = 2:4,
        medium = 5:10,
        large = 11:13
)

set.seed(123) # ensuring reproducibility

# ---- Helper functions to sample combinations ---- 
sample_combinations <- function(studies, size, n_samples) {
    all_combos <- combn(studies, size, simplify = FALSE)
    if(length(all_combos) <= n_samples) return(all_combos)
    sample(all_combos, n_samples)
}

# ---- Generate subsets ----
for (tier_name in names(training_tiers)) {
        tier_dir <- file.path(out_dir, tier_name)
        dir.create(tier_dir, recursive=TRUE, showWarnings=FALSE)
        sizes <- training_tiers[[tier_name]]
        for (size in sizes) {
                # Decide number of combinations based on tier
                n_combos <- switch(
                        tier_name, 
                        small = Inf,
                        medium = 30,
                        large = 10
                )

                # Generate all possible combinations
                if(size <= length(study_names)) {
                        all_combos <- sample_combinations(study_names, size, n_combos)

                        # LOSO: each study not in training is test
                        for (train_set in all_combos) {
                                test_candidates <- setdiff(study_names, train_set)

                                for (test_name in test_candidates) {
                                        subset_train <- combined %>% filter(meta_source %in% train_set)
                                        subset_test <- combined %>% filter(meta_source == test_name)

                                        # Save files
                                        train_idx <- paste(match(train_set, study_names), collapse="-")
                                        train_file <- file.path(tier_dir, paste0("train_", train_idx, "_size", length(train_set), ".csv"))
                                        test_file <- file.path(tier_dir, paste0("test_", match(test_name, study_names), "_train", length(train_set), ".csv"))

                                        write_csv(subset_train, train_file)
                                        write_csv(subset_test, test_file)

                                        message("Tier: ", tier_name,
                                                "; Train size: ", size,
                                                "; Train: ", paste(train_set, collapse=", "),
                                                "; Test: ", test_name)
                                }
                        }
                }
        }
}

message("Done! Created one subset for each number of studies (2–", length(study_names), ") in:\n", out_dir)
