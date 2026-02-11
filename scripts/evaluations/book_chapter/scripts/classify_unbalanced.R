#!/usr/bin/env Rscript

# classify_unbalanced.R - Classification on unbalanced TB data
# Tests Combat supervised vs unsupervised vs unadjusted on deliberately imbalanced data

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

parser <- ArgumentParser(description = "Execute unbalanced TB classification comparison")

parser$add_argument("--adjuster", type = "character", required = TRUE,
                   help = "Batch correction method: unadjusted, combat, or combat_sup")
parser$add_argument("--classifier", type = "character", required = TRUE,
                   help = "Classifier type: logistic, elnet, elasticnet, svm, rf, nnet, knn, xgboost, or shrinkageLDA")
parser$add_argument("--num-datasets", type = "integer", required = TRUE,
                   help = "Number of datasets to include: 3, 4, 5, or 6")
parser$add_argument("--test-study", type = "character", required = TRUE,
                   help = "Test study name (e.g., GSE37250_SA, USA, India, etc.)")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                   help = "Output CSV file path")

args <- parser$parse_args()

# Parameter validation
valid_adjusters <- c("unadjusted", "combat", "combat_sup")
valid_classifiers <- c("rda", "logistic", "elnet", "elasticnet", "svm", "rf", "nnet", "knn", "xgboost", "shrinkageLDA")
valid_num_datasets <- c(3, 4, 5, 6)

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

if (!args$num_datasets %in% valid_num_datasets) {
  cat(sprintf("Error: Invalid num-datasets '%d'. Must be one of: %s\n", 
              args$num_datasets, paste(valid_num_datasets, collapse=", ")))
  quit(status=1)
}

adjuster <- args$adjuster
classifier <- args$classifier
num_datasets <- args$num_datasets
test_study <- args$test_study
output_file <- args$output

# ====================================================================
# MAIN ANALYSIS FUNCTION
# ====================================================================

main_unbalanced_analysis <- function() {
  
  cat("=== UNBALANCED TB CLASSIFICATION ===\n")
  cat(sprintf("Adjuster: %s\n", adjuster))
  cat(sprintf("Classifier: %s\n", classifier))
  cat(sprintf("Num datasets: %d\n", num_datasets))
  cat(sprintf("Test study: %s\n", test_study))
  cat("====================================\n\n")
  
  # Load unbalanced data for this num_datasets + test_study combination
  data_path <- sprintf("data/TB_unbalanced_data_%d_%s.RData", num_datasets, test_study)
  if (!file.exists(data_path)) {
    stop(sprintf("Unbalanced data file not found: %s", data_path))
  }
  
  load(data_path)  # Loads dat_lst and label_lst
  source("scripts/helper.R")
  
  # Define the train_and_evaluate_classifier function locally
  # (copied from classify_adjusters.R to avoid sourcing the entire file)
  train_and_evaluate_classifier <- function(classifier_type, train_data, train_labels, test_data, test_labels) {
    
    # Initialize variables
    trained_model <- NULL
    test_predictions <- NULL
    
    if (classifier_type == "shrinkageLDA") {
      cat("Training shrinkage LDA (sda)...\n")
      
      # Transpose data: R (features x samples) -> sda expects (samples x features)
      X_train <- t(train_data)
      X_test <- t(test_data)
      y_train <- as.factor(train_labels)
      
      # Ensure data is in matrix format (sda requires matrix, not data.frame)
      if (!is.matrix(X_train)) {
        X_train <- as.matrix(X_train)
      }
      if (!is.matrix(X_test)) {
        X_test <- as.matrix(X_test)
      }
      
      cat(sprintf("  X_train dimensions: %d x %d (class: %s)\n", nrow(X_train), ncol(X_train), class(X_train)[1]))
      cat(sprintf("  X_test dimensions: %d x %d (class: %s)\n", nrow(X_test), ncol(X_test), class(X_test)[1]))
      
      # Train shrinkage LDA model
      lda_fit <- sda(
        Xtrain = X_train,
        L = y_train,
        diagonal = FALSE      # Allow correlations between features
      )
      
      # Generate predictions on test set
      pred <- predict(lda_fit, X_test)
      
      # For binary classification, extract probability of positive class
      if (nlevels(y_train) == 2) {
        # Probability of positive class (second column)
        test_predictions <- pred$posterior[, 2]
      } else {
        # Multiclass: pick max probability
        test_predictions <- apply(pred$posterior, 1, max)
      }
      
      trained_model <- list(mod = lda_fit)
      
    } else {
      # Get prediction function for classifier type
      learner_fit <- getPredFunctions(classifier_type)
      
      # Validate data before training
      if (!is.matrix(train_data)) {
        cat(sprintf("Converting train_data to matrix (was %s)\n", class(train_data)))
        train_data <- as.matrix(train_data)
      }
      if (!is.matrix(test_data)) {
        cat(sprintf("Converting test_data to matrix (was %s)\n", class(test_data)))
        test_data <- as.matrix(test_data)
      }
      
      # Train model
      cat(sprintf("Training %s classifier...\n", classifier_type))
      trained_model <- trainPipe(train_set = train_data, train_label = train_labels, 
                                lfit = learner_fit)
      
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
  
  # Use the same data preparation logic as classify_adjusters.R
  # but with unbalanced data
  
  # Filter studies based on num_datasets parameter
  all_studies <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
  selected_studies <- all_studies[1:num_datasets]
  
  dat_lst_filtered <- dat_lst[selected_studies]
  label_lst_filtered <- label_lst[selected_studies]
  study_names <- names(dat_lst_filtered)
  
  cat(sprintf("Running %d-study analysis with studies: %s\n", 
              num_datasets, paste(study_names, collapse=", ")))
  
  # Validate test study is in the filtered list
  if (!test_study %in% study_names) {
    stop(sprintf("Test study '%s' not found in selected studies: %s", 
                 test_study, paste(study_names, collapse=", ")))
  }
  
  # Prepare datasets with CRITICAL MODIFICATION:
  # - Use IMBALANCED training data (from filtered dataset)
  # - Use ORIGINAL test data (from original dataset for consistent evaluation)
  
  # Load ORIGINAL data for test set
  original_data_path <- "data/TB_real_data.RData"
  if (!file.exists(original_data_path)) {
    stop(sprintf("Original data file not found: %s", original_data_path))
  }
  
  # Load original data into separate variables
  original_env <- new.env()
  load(original_data_path, envir = original_env)
  dat_lst_original <- original_env$dat_lst
  label_lst_original <- original_env$label_lst
  
  cat("CRITICAL: Using IMBALANCED training data + ORIGINAL test data\n")
  cat("This tests how methods handle training imbalance while evaluating on consistent test sets\n\n")
  
  test_name <- test_study
  train_name <- setdiff(study_names, test_name)
  
  # Ensure all datasets have the same genes by taking intersection
  all_datasets <- c(train_name, test_name)
  
  # Get common genes from BOTH filtered and original data
  common_genes_filtered <- Reduce(intersect, lapply(dat_lst_filtered[train_name], rownames))
  common_genes_original <- Reduce(intersect, lapply(dat_lst_original[all_datasets], rownames))
  common_genes <- intersect(common_genes_filtered, common_genes_original)
  
  cat(sprintf("Gene intersection: %d common genes across all datasets\n", length(common_genes)))
  
  # Prepare TRAINING data from FILTERED (imbalanced) dataset
  dat_lst_train_subset <- lapply(dat_lst_filtered[train_name], function(x) x[common_genes, , drop=FALSE])
  dat <- do.call(cbind, dat_lst_train_subset)
  batch <- rep(1:length(train_name), times=sapply(dat_lst_train_subset, ncol))
  group <- do.call(c, label_lst_filtered[train_name])
  
  # Prepare TEST data from ORIGINAL (unfiltered) dataset  
  dat_test <- dat_lst_original[[test_name]][common_genes, , drop=FALSE]
  group_test <- label_lst_original[[test_name]]
  
  # Feature reduction (top 1000 most variable genes)
  n_highvar_genes <- 1000
  genes_sel_names <- order(rowVars(dat), decreasing=TRUE)[1:n_highvar_genes]
  dat <- dat[genes_sel_names, ]
  dat_test <- dat_test[genes_sel_names, ]
  
  # Log transformation (same logic as classify_adjusters.R)
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
  cat(sprintf("TRAINING class distribution (imbalanced): LTBI=%d, Active TB=%d\n", 
              sum(group == 0), sum(group == 1)))
  cat(sprintf("TEST class distribution (original): LTBI=%d, Active TB=%d\n", 
              sum(group_test == 0), sum(group_test == 1)))
  
  train_imbalance_ratio <- sum(group == 1) / sum(group == 0)
  test_imbalance_ratio <- sum(group_test == 1) / sum(group_test == 0)
  cat(sprintf("Training TB/LTBI ratio: %.3f, Test TB/LTBI ratio: %.3f\n", 
              train_imbalance_ratio, test_imbalance_ratio))
  
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
  
  # Use the same classification logic as classify_adjusters.R
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
  
  # Create results dataframe
  results_df <- data.frame(
    adjuster = adjuster,
    classifier = classifier,
    num_datasets = num_datasets,
    test_study = test_study,
    mcc = performance["mcc"],
    accuracy = performance["acc"],
    sensitivity = performance["sensitivity"],
    specificity = performance["specificity"],
    auc = performance["auc"],
    train_imbalance_ratio = sum(group == 1) / sum(group == 0),
    test_imbalance_ratio = sum(group_test == 1) / sum(group_test == 0),
    analysis_type = "unbalanced_tb",
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
  result <- main_unbalanced_analysis()
}, error = function(e) {
  cat(sprintf("[ERROR] %s\n", e$message), file = stderr())
  quit(status = 1)
})