#!/usr/bin/env Rscript

# classify_class_imbalanced.R - Classification on systematically class-imbalanced TB data
# Tests Combat supervised vs unsupervised vs unadjusted across multiple imbalance levels

suppressMessages(suppressWarnings({
  required_packages <- c("glmnet", "SummarizedExperiment", "sva", "DESeq2", 
                        "ROCR", "ggplot2", "gridExtra", "reshape2", 
                        "dplyr", "purrr", "nnls", "batchelor",
                        "argparse", "class", "xgboost", "sda")
  sapply(required_packages, require, character.only=TRUE, quietly=TRUE)
}))

# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Execute class-imbalanced TB classification")

parser$add_argument("--adjuster", type = "character", required = TRUE,
                   help = "Batch correction method: unadjusted, combat, or combat_sup")
parser$add_argument("--classifier", type = "character", required = TRUE,
                   help = "Classifier type: logistic, elnet, elasticnet, svm, rf, nnet, knn, xgboost, or shrinkageLDA")
parser$add_argument("--scenario-name", type = "character", required = TRUE,
                   help = "Scenario name (e.g., GSE37250_SA_USA_India_imbal20_testGSE37250_SA)")
parser$add_argument("--data-dir", type = "character", required = TRUE,
                   help = "Directory containing imbalanced datasets")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                   help = "Output CSV file path")

args <- parser$parse_args()

# Parameter validation
valid_adjusters <- c("unadjusted", "combat", "combat_sup")
valid_classifiers <- c("rda", "logistic", "elnet", "elasticnet", "svm", "rf", "nnet", "knn", "xgboost", "shrinkageLDA")

if (!args$adjuster %in% valid_adjusters) {
  cat(sprintf("Error: Invalid adjuster '%s'. Must be one of: %s\n", 
              args$adjuster, paste(valid_adjusters, collapse=", ")))
  quit(status=1)
}

if (!args$classifier %in% valid_classifiers) {
  cat(sprintf("Error: Invalid classifier '%s'. Must be one of: %s\n", 
              args$classifier, paste(valid_classifiers, collapse=", ")))
  quit(status=1)
}

adjuster <- args$adjuster
classifier <- args$classifier
scenario_name <- args$scenario_name
data_dir <- args$data_dir
output_file <- args$output

# ====================================================================
# HELPER FUNCTIONS
# ====================================================================

# Parse scenario name to extract components
parse_scenario_name <- function(scenario_name, dat_lst_names) {
  # Expected format: TRAININGPAIR-imbalXX-testDATASET-repN
  # Split by the known markers: -imbal, -test, -rep
  
  parts <- strsplit(scenario_name, "-imbal|-test|-rep")[[1]]
  
  if (length(parts) < 3) {
    stop(sprintf("Cannot parse scenario name (expected at least 3 parts): %s", scenario_name))
  }
  
  training_pair <- parts[1]
  imbalance_pct <- as.numeric(parts[2]) / 100
  test_dataset <- parts[3]
  replicate <- if (length(parts) >= 4) as.numeric(parts[4]) else 1
  
  # Infer train_datasets from the actual dat_lst keys
  # The training pair should match exactly 2 datasets from dat_lst (excluding test)
  available_train_datasets <- setdiff(dat_lst_names, test_dataset)
  
  if (length(available_train_datasets) != 2) {
    stop(sprintf("Expected exactly 2 training datasets in loaded data (excluding test=%s), found %d: %s",
                test_dataset, length(available_train_datasets), paste(available_train_datasets, collapse=", ")))
  }
  
  train_datasets <- available_train_datasets
  
  return(list(
    train_datasets = train_datasets,
    test_dataset = test_dataset,
    imbalance_pct = imbalance_pct,
    training_pair = training_pair,
    replicate = replicate
  ))
}

# Define the train_and_evaluate_classifier function
train_and_evaluate_classifier <- function(classifier_type, train_data, train_labels, test_data, test_labels) {
  
  # Source helper functions
  source("scripts/helper.R")
  
  # Initialize variables
  trained_model <- NULL
  test_predictions <- NULL
  
  if (classifier_type == "shrinkageLDA") {
    cat("Training shrinkage LDA (sda)...\n")
    
    # Transpose data: R (features x samples) -> sda expects (samples x features)
    X_train <- t(train_data)
    X_test <- t(test_data)
    y_train <- as.factor(train_labels)
    
    # Ensure data is in matrix format
    if (!is.matrix(X_train)) X_train <- as.matrix(X_train)
    if (!is.matrix(X_test)) X_test <- as.matrix(X_test)
    
    # Train shrinkage LDA model
    lda_fit <- sda(Xtrain = X_train, L = y_train, diagonal = FALSE)
    
    # Generate predictions on test set
    pred <- predict(lda_fit, X_test)
    
    # For binary classification, extract probability of positive class
    if (nlevels(y_train) == 2) {
      test_predictions <- pred$posterior[, 2]
    } else {
      test_predictions <- apply(pred$posterior, 1, max)
    }
    
    trained_model <- list(mod = lda_fit)
    
  } else {
    # Get prediction function for classifier type
    learner_fit <- getPredFunctions(classifier_type)
    
    # Validate data before training
    if (!is.matrix(train_data)) train_data <- as.matrix(train_data)
    if (!is.matrix(test_data)) test_data <- as.matrix(test_data)
    
    # Train model
    cat(sprintf("Training %s classifier...\n", classifier_type))
    trained_model <- trainPipe(train_set = train_data, train_label = train_labels, lfit = learner_fit)
    
    # Generate predictions on test set
    cat(sprintf("Generating predictions...\n"))
    test_predictions <- predWrapper(trained_model$mod, test_data, classifier_type)
  }
  
  # Calculate performance metrics
  perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc")
  
  # Create predictions list in format expected by perf_wrapper
  predictions_list <- list()
  predictions_list[[adjuster]] <- test_predictions
  
  # Calculate performance using original perf_wrapper function
  perf_results <- perf_wrapper(perf_measures, predictions_list, test_labels)
  
  # Calculate confusion matrix elements and derived metrics
  confusion_results <- confusion_matrix_wrapper(predictions_list, test_labels)
  
  # Combine performance metrics and confusion matrix metrics
  combined_results <- rbind(perf_results, confusion_results)
  
  # Extract performance values for this method
  perf_values <- combined_results[, adjuster]
  names(perf_values) <- rownames(combined_results)
  
  return(list(
    model = trained_model,
    predictions = test_predictions,
    performance = perf_values
  ))
}

# ====================================================================
# MAIN ANALYSIS FUNCTION
# ====================================================================

main_class_imbalanced_analysis <- function() {
  
  cat("=== CLASS-IMBALANCED TB CLASSIFICATION ===\n")
  cat(sprintf("Adjuster: %s\n", adjuster))
  cat(sprintf("Classifier: %s\n", classifier))
  cat(sprintf("Scenario: %s\n", scenario_name))
  cat("==========================================\n\n")
  
  # Load imbalanced data
  data_file <- file.path(data_dir, sprintf("%s.RData", scenario_name))
  if (!file.exists(data_file)) {
    stop(sprintf("Data file not found: %s", data_file))
  }
  
  load(data_file)  # Loads dat_lst and label_lst
  
  # Parse scenario name using actual dataset names from loaded data
  scenario_info <- parse_scenario_name(scenario_name, names(dat_lst))
  
  cat(sprintf("Parsed scenario:\n"))
  cat(sprintf("  Training pair: %s\n", scenario_info$training_pair))
  cat(sprintf("  Train datasets: %s\n", paste(scenario_info$train_datasets, collapse = ", ")))
  cat(sprintf("  Test dataset: %s\n", scenario_info$test_dataset))
  cat(sprintf("  Imbalance level: %.0f%%\n", scenario_info$imbalance_pct * 100))
  cat(sprintf("  Replicate: %d\n", scenario_info$replicate))
  
  # Validate that all required datasets are present
  required_datasets <- c(scenario_info$train_datasets, scenario_info$test_dataset)
  missing_datasets <- setdiff(required_datasets, names(dat_lst))
  if (length(missing_datasets) > 0) {
    stop(sprintf("Missing datasets in loaded data: %s", paste(missing_datasets, collapse = ", ")))
  }
  
  # Prepare training and test data
  train_datasets <- scenario_info$train_datasets
  test_dataset <- scenario_info$test_dataset
  
  # Ensure all datasets have the same genes by taking intersection
  all_datasets <- c(train_datasets, test_dataset)
  common_genes <- Reduce(intersect, lapply(dat_lst[all_datasets], rownames))
  
  cat(sprintf("Gene intersection: %d common genes across all datasets\n", length(common_genes)))
  
  # Prepare training data (combine training datasets)
  dat_lst_train_subset <- lapply(dat_lst[train_datasets], function(x) x[common_genes, , drop=FALSE])
  dat <- do.call(cbind, dat_lst_train_subset)
  batch <- rep(1:length(train_datasets), times=sapply(dat_lst_train_subset, ncol))
  group <- do.call(c, label_lst[train_datasets])
  
  # Prepare test data
  dat_test <- dat_lst[[test_dataset]][common_genes, , drop=FALSE]
  group_test <- label_lst[[test_dataset]]
  
  # Feature reduction (top 1000 most variable genes)
  n_highvar_genes <- min(1000, nrow(dat))
  genes_sel_names <- order(rowVars(dat), decreasing=TRUE)[1:n_highvar_genes]
  dat <- dat[genes_sel_names, ]
  dat_test <- dat_test[genes_sel_names, ]
  
  # Log transformation
  needs_log_transform_train <- max(dat) > 100 || mean(dat) > 20 || (max(dat) / median(dat)) > 50
  needs_log_transform_test <- max(dat_test) > 100 || mean(dat_test) > 20 || (max(dat_test) / median(dat_test)) > 50
  
  if (needs_log_transform_train) {
    train_min <- min(dat)
    if (train_min < 0) dat <- dat - train_min
    dat <- log2(dat + 1)
  }
  
  if (needs_log_transform_test) {
    test_min <- min(dat_test)
    if (test_min < 0) dat_test <- dat_test - test_min
    dat_test <- log2(dat_test + 1)
  }
  
  cat(sprintf("Data preparation completed:\n"))
  cat(sprintf("  Training samples: %d\n", ncol(dat)))
  cat(sprintf("  Test samples: %d\n", ncol(dat_test)))
  cat(sprintf("  Features (genes): %d\n", nrow(dat)))
  cat(sprintf("  Training batches: %d\n", length(unique(batch))))
  
  # Print class distribution to show the imbalance effect
  train_active <- sum(group == 1)
  train_latent <- sum(group == 0)
  test_active <- sum(group_test == 1)
  test_latent <- sum(group_test == 0)
  
  cat(sprintf("TRAINING class distribution: LTBI=%d (%.1f%%), Active TB=%d (%.1f%%)\n", 
              train_latent, 100 * train_latent / (train_active + train_latent),
              train_active, 100 * train_active / (train_active + train_latent)))
  cat(sprintf("TEST class distribution: LTBI=%d (%.1f%%), Active TB=%d (%.1f%%)\n", 
              test_latent, 100 * test_latent / (test_active + test_latent),
              test_active, 100 * test_active / (test_active + test_latent)))
  
  # Print batch-specific class distributions
  cat("Training batch-specific class distributions:\n")
  for (i in 1:length(train_datasets)) {
    batch_mask <- batch == i
    batch_group <- group[batch_mask]
    batch_active <- sum(batch_group == 1)
    batch_latent <- sum(batch_group == 0)
    batch_total <- batch_active + batch_latent
    
    cat(sprintf("  %s (batch %d): LTBI=%d (%.1f%%), Active TB=%d (%.1f%%), Total=%d\n",
                train_datasets[i], i, batch_latent, 100 * batch_latent / batch_total,
                batch_active, 100 * batch_active / batch_total, batch_total))
  }
  
  # ====================================================================
  # BATCH CORRECTION
  # ====================================================================
  
  if (adjuster == "unadjusted") {
    cat("Using unadjusted data (no batch correction)\n")
    dat_corrected <- dat
    dat_test_corrected <- dat_test
    
  } else if (adjuster == "combat") {
    cat("Applying Combat (unsupervised) batch correction\n")
    
    # Step 1: Correct training data without using labels
    dat_corrected <- ComBat(dat, batch=batch, mod=NULL)
    
    # Step 2: Adjust test data to match corrected training distribution
    ref_batch_id <- 1
    test_batch_id <- max(batch) + 1
    combined_dat <- cbind(dat_corrected, dat_test)
    combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), 
                       rep(test_batch_id, ncol(dat_test)))
    
    combat_combined <- ComBat(combined_dat, batch=combined_batch, 
                             mod=NULL, ref.batch=ref_batch_id)
    
    dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
    
  } else if (adjuster == "combat_sup") {
    cat("Applying Combat (supervised) batch correction\n")
    
    # Step 1: Correct training data while preserving biological signal
    dat_corrected <- ComBat(dat, batch=batch, mod=model.matrix(~group))
    
    # Step 2: Adjust test data to match corrected training distribution
    ref_batch_id <- 1
    test_batch_id <- max(batch) + 1
    combined_dat <- cbind(dat_corrected, dat_test)
    combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), 
                       rep(test_batch_id, ncol(dat_test)))
    
    # Apply Combat without using test labels
    combat_combined <- ComBat(combined_dat, batch=combined_batch, 
                             mod=NULL, ref.batch=ref_batch_id)
    
    dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
  }
  
  # ====================================================================
  # CLASSIFICATION
  # ====================================================================
  
  cat(sprintf("Training %s classifier...\n", classifier))
  
  # Train and evaluate classifier
  result <- train_and_evaluate_classifier(
    classifier_type = classifier,
    train_data = dat_corrected,
    train_labels = group,
    test_data = dat_test_corrected,
    test_labels = group_test
  )
  
  # ====================================================================
  # SAVE RESULTS
  # ====================================================================
  
  # Extract performance metrics from the result
  performance <- result$performance
  
  # Calculate batch-class confounding metrics
  batch_tb_ratios <- sapply(1:length(train_datasets), function(i) {
    batch_mask <- batch == i
    batch_group <- group[batch_mask]
    sum(batch_group == 1) / length(batch_group)
  })
  confounding_strength <- sd(batch_tb_ratios)
  
  # Create results dataframe
  results_df <- data.frame(
    scenario_name = scenario_name,
    training_pair = scenario_info$training_pair,
    test_dataset = scenario_info$test_dataset,
    train_dataset_1 = train_datasets[1],
    train_dataset_2 = train_datasets[2],
    imbalance_pct = scenario_info$imbalance_pct,
    replicate = scenario_info$replicate,
    adjuster = adjuster,
    classifier = classifier,
    mcc = performance["mcc"],
    accuracy = performance["acc"],
    sensitivity = performance["sensitivity"],
    specificity = performance["specificity"],
    auc = performance["auc"],
    f1_score = performance["f"],
    train_samples = ncol(dat),
    test_samples = ncol(dat_test),
    train_active = train_active,
    train_latent = train_latent,
    test_active = test_active,
    test_latent = test_latent,
    train_active_pct = train_active / (train_active + train_latent),
    test_active_pct = test_active / (test_active + test_latent),
    batch_confounding_strength = confounding_strength,
    analysis_type = "class_imbalanced",
    stringsAsFactors = FALSE
  )
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write.csv(results_df, output_file, row.names = FALSE)
  
  cat(sprintf("Results saved to: %s\n", output_file))
  cat(sprintf("MCC: %.4f, Accuracy: %.4f, AUC: %.4f\n", 
              performance["mcc"], performance["acc"], performance["auc"]))
  
  return(results_df)
}

# ====================================================================
# EXECUTE WITH ERROR HANDLING
# ====================================================================

tryCatch({
  result <- main_class_imbalanced_analysis()
}, error = function(e) {
  cat(sprintf("[ERROR] %s\n", e$message), file = stderr())
  quit(status = 1)
})