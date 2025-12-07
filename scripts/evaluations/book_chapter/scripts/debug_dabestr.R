#!/usr/bin/env Rscript

library(dplyr)
library(dabestr)

# Load the actual data
data <- read.csv("/home/phr23/confounded_analysis/grp_batch_effects/outputs/book_chapter/adjusters_on_classifiers.csv", stringsAsFactors = FALSE)

# Filter to MCC and one classifier
mxe_data <- data[data$metric == "mcc" & !is.na(data$n_datasets), ]
classifier_data <- mxe_data[mxe_data$classifier == "logistic", ]

cat("Total rows:", nrow(classifier_data), "\n")
cat("Unique adjusters:", length(unique(classifier_data$adjuster)), "\n")
cat("Unique n_datasets:", unique(classifier_data$n_datasets), "\n")

# Prepare data like the script does
plot_data <- classifier_data %>%
  select(adjuster, classifier, n_datasets, test_study, value) %>%
  mutate(condition_id = paste(classifier, n_datasets, test_study, sep = "_"))

# Filter to 3 studies
subset_data <- plot_data %>% filter(n_datasets == 3)

cat("\nSubset data (3 studies):\n")
cat("Rows:", nrow(subset_data), "\n")
cat("Adjusters:", paste(unique(subset_data$adjuster), collapse = ", "), "\n")
cat("Condition IDs:", length(unique(subset_data$condition_id)), "\n")

# Try to create dabest object
adjusters_to_compare <- unique(subset_data$adjuster)
cat("\nAttempting dabestr::load with adjusters:", paste(adjusters_to_compare, collapse = ", "), "\n")

tryCatch({
  dabest_obj <- dabestr::load(
    data = subset_data,
    x = adjuster,
    y = value,
    idx = adjusters_to_compare,
    paired = "baseline",
    id_col = condition_id
  )
  
  cat("SUCCESS: dabest_obj created\n")
  cat("Class:", class(dabest_obj), "\n")
  
  # Try mean_diff
  dabest_diff <- dabestr::mean_diff(dabest_obj, perm_count = 100)
  cat("SUCCESS: mean_diff computed\n")
  cat("Class:", class(dabest_diff), "\n")
  cat("Has boot_result:", !is.null(dabest_diff$boot_result), "\n")
  cat("Has permtest_pvals:", !is.null(dabest_diff$permtest_pvals), "\n")
  
  if (!is.null(dabest_diff$boot_result)) {
    cat("boot_result rows:", nrow(dabest_diff$boot_result), "\n")
    print(dabest_diff$boot_result)
  }
  
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  cat("Traceback:\n")
  print(traceback())
})
