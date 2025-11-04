# This script will take arguments
# Usage: Rscript run_scaling_experiment.R --adjuster log_combat

# Loop through subsets 2-14 and load the csv, 
#  - identify meta_source values,
#  - split train/test, 
#  - apply chosen adjuster, 
#  - call modeling function, 
#  - and write to csv

suppressPackageStartupMessages({
    library(optparse)
    library(readr)
    library(dplyr)
    library(reticulate)
})

# ---- Parse Arguments ----
option_list <- list(
    make_option("--adjuster", type = "character", help = "Adjustment method name")
)
opt <- parse_args(OptionParser(option_list = option_list))
adjuster <- opt$adjuster

cat("Running scaling experiment with adjuster:", adjuster, "\n")

# ---- Source R scripts ----
source("scripts/adjust/adjust.R")
source("scripts/evaluations/classify_er_all_datasets/sklearn_helpers.R", local = TRUE)

apply_adjustment <- function(df, method, test_source) {
    if (method != "gmm") {
        df[, sapply(df, is.numeric)] <- log1p(df[, sapply(df, is.numeric)])
        cat("Applying log transform before ", method, "\n")
    }
    
    cat("Applying adjustment method: ", method, "\n")
    if (method == "min_mean") {
        adjusted <- adjust_min_mean(df)
    } else if (method == "log_combat") {
        adjusted <- adjust_log_combat(df, ref_batch = unique(df$meta_source)[1])
    } else if (method == "mnn") {
        adjusted <- adjust_mnn(df, merge_last = test_source)
    } else if (method == "gmm") {
        adjusted <- adjust_gmm(df)
    } else {
        stop("Unknown adjuster: ", method)
    }
    return(adjusted)
}

# Leave one out loop
for (n_studies in 2:14) {
    subset_path <- sprintf("grp_batch_effects/data/all_combined_subsets/subset_%dstudies.csv", n_studies)
    cat("\nProcessing:", subset_path, "\n")
    df <- read_csv(subset_path)

    sources <- unique(df$meta_source)
    for (test_source in sources) {
        adjusted_df <- apply_adjustment(df, method = adjuster, test_source = test_source)

        df_train <- adjusted_df %>% filter(meta_source != test_source)
        df_test <- adjusted_df %>% filter(meta_source == test_source)

        # Run classification (reusing function)
        metrics_df <- run_sklearn_model(
            X_train = df_train %>% select(where(is.numeric)),
            y_train = df_train$meta_er_status,
            X_test = df_test %>% select(where(is.numeric)),
            y_test = df_test$meta_er_status
        )

        # Add metadata to metrics_df
        metrics_df$adjuster <- adjuster
        metrics_df$subset <- n_studies
        metrics_df$test <- test_source

        # Save result
        out_path <- sprintf("grp_batch_effects/outputs/er_all_datasets/%dstudies_test%s.csv",
                            n_studies, test_source)
        if (!dir.exists(dirname(out_path))) dir.create(dirname(out_path), recursive = TRUE)
        write_csv(metrics_df, out_path)
    }
}