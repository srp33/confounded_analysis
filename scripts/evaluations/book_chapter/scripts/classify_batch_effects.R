#!/usr/bin/env Rscript

# classify_batch_effects.R
# Single job script for batch effects analysis extracted from 1_simpipe.R
# Executes one classifier × one parameter set × one seed combination

# Suppress all output and warnings for cleaner job execution
options(warn = -1)
suppressMessages(suppressWarnings({
  
  # Load required packages
  required_packages <- c("SummarizedExperiment", "plyr", "sva", "MCMCpack", "ROCR", "ggplot2", 
                        "limma", "nnls", "glmnet", "rpart", "genefilter", "nnet", "e1071", 
                        "RcppArmadillo", "foreach", "parallel", "doParallel", "ranger", "scales",
                        "purrr", "dplyr", "batchelor", "reticulate")
  
  package_results <- sapply(required_packages, require, character.only=TRUE, quietly=TRUE)
  failed_packages <- names(package_results)[!package_results]
  if(length(failed_packages) > 0) {
    cat("Warning: Failed to load packages:", paste(failed_packages, collapse=", "), "\n", file=stderr())
  }
}))

# ====================================================================
# [RVC MODIFICATION] RETICULATE / RVC SETUP
# ====================================================================
# This section attempts to load the necessary Python modules for RVC.
# We define them as NULL globally and use <<- to assign them if found.
RVC_py <- NULL
np_py <- NULL

tryCatch({
  cat("Attempting to import Python modules for RVC...\n")
  rvm_module <- import("scikit_rvm")
  RVC_py <<- rvm_module$RVC
  np_py <<- import("numpy")
  cat("Successfully imported scikit_rvm and numpy.\n")
}, error = function(e) {
  cat("[WARNING] Could not import Python modules 'scikit_rvm' or 'numpy'.\n")
  cat("[WARNING] The 'rvc' classifier will be unavailable.\n")
  cat(sprintf("[WARNING] Python Error: %s\n", e$message))
})

####  Command Line Interface  ####

#' Parse command line arguments with validation and help text
parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Help text
  help_text <- "
Usage: 
  ./run_in_apptainer.sh scripts/evaluations/book_chapter/scripts/classify_batch_effects.R --classifier <type> --mean <value> --var <value> --seed <value> -o <output>
  
  Or directly: Rscript classify_batch_effects.R --classifier <type> --mean <value> --var <value> --seed <value> -o <output>

Arguments:
  --classifier    Classifier type (logistic, elnet, elasticnet, svm, rf, nnet, knn, xgboost, rvc)
  --mean         Batch effect mean parameter (typically 5)
  --var          Batch effect variance parameter (1, 3, or 5)  
  --seed         Random seed for reproducibility (integer)
  -o, --output   Output CSV file path
  --usage        Show this help message (use instead of -h to avoid conflict with apptainer)

Examples:
  ./run_in_apptainer.sh scripts/evaluations/book_chapter/scripts/classify_batch_effects.R --classifier logistic --mean 5 --var 1 --seed 42 -o results/batch_effects_logistic_5_1_42.csv
  ./run_in_apptainer.sh scripts/evaluations/book_chapter/scripts/classify_batch_effects.R --usage

This script extracts single job functionality from 1_simpipe.R for parallel execution.
"
  
  # Check for help request (avoid conflict with apptainer script's -h/--help)
  if(length(args) == 0 || any(args %in% c("--usage", "--script-help"))) {
    cat(help_text)
    quit(status = 0)
  }
  
  # Parse arguments
  if(length(args) != 10) {
    cat("Error: Expected 10 arguments, got", length(args), "\n", file=stderr())
    cat("Use --help for usage information\n", file=stderr())
    quit(status = 1)
  }
  
  # Extract parameter values
  classifier_idx <- which(args == "--classifier")
  mean_idx <- which(args == "--mean")
  var_idx <- which(args == "--var")
  seed_idx <- which(args == "--seed")
  output_idx <- which(args %in% c("-o", "--output"))
  
  if(length(classifier_idx) != 1 || length(mean_idx) != 1 || length(var_idx) != 1 || 
     length(seed_idx) != 1 || length(output_idx) != 1) {
    cat("Error: Missing or duplicate required arguments\n", file=stderr())
    cat("Use --help for usage information\n", file=stderr())
    quit(status = 1)
  }
  
  classifier <- args[classifier_idx + 1]
  mean_val <- as.numeric(args[mean_idx + 1])
  var_val <- as.numeric(args[var_idx + 1])
  seed_val <- as.integer(args[seed_idx + 1])
  output_path <- args[output_idx + 1]
  
  # Validate parameters
  valid_classifiers <- c("logistic", "elnet", "elasticnet", "svm", "rf", "nnet", "knn", "xgboost", "rvc")
  if(!classifier %in% valid_classifiers) {
    cat("Error: Invalid classifier. Must be one of:", paste(valid_classifiers, collapse=", "), "\n", file=stderr())
    quit(status = 1)
  }
  
  if(is.na(mean_val) || mean_val < 0) {
    cat("Error: Mean parameter must be a non-negative number\n", file=stderr())
    quit(status = 1)
  }
  
  if(is.na(var_val) || !var_val %in% c(1, 3, 5)) {
    cat("Error: Variance parameter must be 1, 3, or 5\n", file=stderr())
    quit(status = 1)
  }
  
  if(is.na(seed_val) || seed_val < 0) {
    cat("Error: Seed must be a non-negative integer\n", file=stderr())
    quit(status = 1)
  }
  
  if(is.null(output_path) || output_path == "") {
    cat("Error: Output path cannot be empty\n", file=stderr())
    quit(status = 1)
  }
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_path)
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  list(
    classifier = classifier,
    mean = mean_val,
    var = var_val,
    seed = seed_val,
    output = output_path
  )
}

# Parse command line arguments
params <- parse_arguments()

# Load dependencies
source("scripts/helper.R")

# Validate RVC dependencies if RVC classifier is requested
if (params$classifier == "rvc" && (is.null(RVC_py) || is.null(np_py))) {
  cat("Error: Classifier 'rvc' was requested, but Python dependencies 'scikit-rvm' or 'numpy' could not be imported.\n", file=stderr())
  cat("Please install them using: pip install scikit-rvm numpy\n", file=stderr())
  quit(status = 1)
}

#### Data Preparation Logic (Extracted from 1_simpipe.R) ####

#' Load and prepare training/test data with gene selection
prepare_data <- function() {
  # Load the combined dataset
  load("data/combined_sub.RData")
  
  # Select 1000 genes with largest variance in training set (Africa)
  var_trn <- rowVars(train_expr)
  genes_sel <- rownames(train_expr)[order(var_trn, decreasing=TRUE)[1:1000]]
  train_expr <- train_expr[genes_sel, ]
  test_expr <- test_expr[genes_sel, ]
  
  list(
    train_expr = train_expr,
    test_expr = test_expr,
    y_train = y_train,
    y_test = y_test
  )
}

#' Subset training data into batches and filter genes
#' @param train_expr Training expression matrix
#' @param y_train Training labels
#' @param N_sample_size Sample size per batch
#' @param N_batch Number of batches (fixed at 3)
subset_and_filter_data <- function(train_expr, y_train, N_sample_size, N_batch = 3) {
  
  # Subset training set in batches using subsetBatch function from helper.R
  batches_ind <- subsetBatch(condition = y_train, N_sample_size = N_sample_size, N_batch = N_batch)
  
  # Create batch assignments
  batch <- rep(0, ncol(train_expr))
  iwalk(batches_ind, ~{batch[.x] <<- .y})
  
  # Extract samples for selected batches
  curr_train_expr <- train_expr[, do.call(c, batches_ind)]
  curr_y_train <- y_train[do.call(c, batches_ind)]
  batch <- batch[do.call(c, batches_ind)]  
  batches_ind <- map(1:N_batch, ~which(batch == .x))
  
  # Remove genes with only 0 values in any batch
  # This is critical for preventing numerical issues in batch correction
  g_keep <- map(1:N_batch, ~{
    which(apply(curr_train_expr[, batch == .x], 1, function(x){!all(x==0)}))
  })
  g_keep <- Reduce(intersect, g_keep)  
  
  curr_train_expr <- curr_train_expr[g_keep, ]
  
  list(
    train_expr = curr_train_expr,
    y_train = curr_y_train,
    batch = batch,
    batches_ind = batches_ind,
    gene_filter = g_keep
  )
}

#### Batch Effect Simulation Logic (Extracted from 1_simpipe.R) ####

#' Generate hyperparameters for batch effect simulation
#' @param max_batch_mean Maximum batch effect mean
#' @param max_batch_var Maximum batch effect variance  
#' @param N_batch Number of batches
generate_hyperparameters <- function(max_batch_mean, max_batch_var, N_batch) {
  
  # Generate hyperparameters using identical logic from 1_simpipe.R
  hyper_pars <- list(
    hyper_mu = seq(from = -max_batch_mean, to = max_batch_mean, length.out = N_batch),  
    hyper_sd = sqrt(rep(0.01, N_batch)),
    hyper_alpha = mv2ab(m = seq(from = 1/max_batch_var, to = max_batch_var, length.out = N_batch), 
                       v = rep(0.01, N_batch))$alpha,
    hyper_beta = mv2ab(m = seq(from = 1/max_batch_var, to = max_batch_var, length.out = N_batch), 
                      v = rep(0.01, N_batch))$beta
  )
  
  return(hyper_pars)
}

#' Simulate batch effects on training data
#' @param train_expr Training expression matrix
#' @param y_train Training labels
#' @param batch Batch assignments
#' @param batches_ind List of batch indices
#' @param hyper_pars Hyperparameters for simulation
simulate_batch_effects <- function(train_expr, y_train, batch, batches_ind, hyper_pars) {
  
  # Set random seed for reproducibility
  set.seed(params$seed)
  
  # Simulate batch effect using simBatch function from helper.R
  # This preserves the exact same methodology as the original implementation
  sim_batch_res <- simBatch(dat = train_expr, 
                           condition = y_train, 
                           batches_ind = batches_ind, 
                           batch = batch, 
                           hyper_pars = hyper_pars)
  
  return(sim_batch_res$new_dat)
}

#### Single Classifier Training and Evaluation (Extracted from 1_simpipe.R) ####

#' Normalize datasets before training
#' @param train_expr Training expression matrix
#' @param test_expr Test expression matrix
#' @param train_expr_batch Training expression matrix with batch effects
#' @param batches_ind List of batch indices
#' @param norm_data Whether to normalize data
normalize_datasets <- function(train_expr, test_expr, train_expr_batch, batches_ind, norm_data = TRUE) {
  
  if(norm_data) {
    train_expr_norm <- normalizeData(train_expr)
    test_expr_norm <- normalizeData(test_expr)
    train_expr_batch_whole_norm <- normalizeData(train_expr_batch)
    
    # Normalize within each batch
    train_expr_batch_norm <- matrix(NA, nrow = nrow(train_expr_batch), ncol = ncol(train_expr_batch), 
                                   dimnames = dimnames(train_expr_batch))
    iwalk(batches_ind, ~{
      train_expr_batch_norm[, .x] <<- normalizeData(train_expr_batch[, .x])
    })
  } else {
    train_expr_norm <- train_expr
    test_expr_norm <- test_expr
    train_expr_batch_whole_norm <- train_expr_batch_norm <- train_expr_batch
  }
  
  list(
    train_expr_norm = train_expr_norm,
    test_expr_norm = test_expr_norm,
    train_expr_batch_whole_norm = train_expr_batch_whole_norm,
    train_expr_batch_norm = train_expr_batch_norm
  )
}

#' Apply batch correction methods
#' @param train_expr_batch Training data with batch effects
#' @param test_expr Test expression data
#' @param batch Batch assignments
#' @param y_train Training labels
apply_batch_corrections <- function(train_expr_batch, test_expr, batch, y_train) {
  
  # ComBat with reference batch (use first batch as reference)
  ref_batch <- min(batch)
  combined_dat <- cbind(train_expr_batch, test_expr)
  combined_batch <- c(batch, rep(ref_batch, ncol(test_expr)))  # Test set uses reference batch
  combined_labels <- c(y_train, rep(0, ncol(test_expr)))  # Use dummy labels for test
  train_expr_combat <- ComBat(combined_dat, batch = combined_batch, 
                             mod = model.matrix(~combined_labels), ref.batch = ref_batch)
  
  # Split back into training and test
  train_expr_combat_adj <- train_expr_combat[, 1:ncol(train_expr_batch)]
  test_expr_combat_adj <- train_expr_combat[, (ncol(train_expr_batch) + 1):ncol(train_expr_combat)]
  
  # MNN correction with test set last in merge order
  library(batchelor, quietly = TRUE)
  combined_dat <- cbind(train_expr_batch, test_expr)
  combined_batch <- c(batch, rep(max(batch) + 1, ncol(test_expr)))  # Test set gets new batch ID
  
  # Create batch list for MNN (each batch as separate matrix)
  unique_batches <- sort(unique(combined_batch))
  batch_list <- map(unique_batches, ~{
    combined_dat[, combined_batch == .x]
  })
  
  # Apply MNN correction with test set last in merge order
  mnn_result <- do.call(fastMNN, c(batch_list, list(merge.order = seq_along(unique_batches))))
  
  # Get corrected data (handle different batchelor versions)
  if("corrected" %in% assayNames(mnn_result)) {
    corrected_combined <- assay(mnn_result, "corrected")
  } else if("reconstructed" %in% assayNames(mnn_result)) {
    corrected_combined <- assay(mnn_result, "reconstructed")
  } else {
    # Fallback to first assay
    corrected_combined <- assay(mnn_result, 1)
  }
  
  # Split back into training and test
  train_expr_mnn_adj <- corrected_combined[, 1:ncol(train_expr_batch)]
  test_expr_mnn_adj <- corrected_combined[, (ncol(train_expr_batch) + 1):ncol(corrected_combined)]
  
  list(
    train_combat = train_expr_combat_adj,
    test_combat = test_expr_combat_adj,
    train_mnn = train_expr_mnn_adj,
    test_mnn = test_expr_mnn_adj
  )
}

#' Train and evaluate single classifier
#' @param classifier_type Type of classifier to use
#' @param normalized_data List of normalized datasets
#' @param corrected_data List of batch-corrected datasets
#' @param y_train Training labels
#' @param y_test Test labels
train_and_evaluate_classifier <- function(classifier_type, normalized_data, corrected_data, y_train, y_test) {
  
  # Define training configurations
  training_configs <- list(
    NoBatch = list(
      train = normalized_data$train_expr_norm,
      test = normalized_data$test_expr_norm
    ),
    Batch = list(
      train = normalized_data$train_expr_batch_whole_norm,
      test = normalized_data$test_expr_norm
    ),
    ComBat = list(
      train = normalizeData(corrected_data$train_combat),
      test = normalizeData(corrected_data$test_combat)
    ),
    MNNcorrect = list(
      train = normalizeData(corrected_data$train_mnn),
      test = normalizeData(corrected_data$test_mnn)
    )
  )
  
  # Train and predict for each configuration
  if (classifier_type == "rvc") {
    # RVC needs special handling due to Python interop
    cat("Training rvc classifier (using reticulate)...\n")
    
    predictions <- map(training_configs, ~{
      # Transpose data: R (features x samples) -> Python (samples x features)
      X_train_t <- t(.x$train)
      X_test_t <- t(.x$test)
      
      # Convert to Python
      X_train_py <- r_to_py(X_train_t)
      y_train_r <- as.numeric(as.factor(y_train)) - 1
      y_train_py <- r_to_py(y_train_r)
      X_test_py <- r_to_py(X_test_t)
      
      # Train RVC model with rbf kernel
      model_py <- RVC_py(kernel = "rbf")
      model_py$fit(X_train_py, y_train_py)
      
      # Get predictions (returns n_samples x n_classes)
      preds_py <- model_py$predict_proba(X_test_py)
      
      # Convert back to R and extract positive class probability
      preds_r <- py_to_r(preds_py)
      preds_r[, 2]  # Column 2 is positive class "1"
    })
    
  } else {
    # Standard R classifiers
    learner_fit <- getPredFunctions(classifier_type)
    
    predictions <- map(training_configs, ~{
      pred_res <- trainPipe(train_set = .x$train, train_label = y_train, 
                           test_set = .x$test, lfit = learner_fit)
      pred_res$pred_tst_prob
    })
  }
  
  # Calculate performance metrics
  perf_measures <- c("mxe", "auc", "acc", "f", "err")
  perf_df <- perf_wrapper(perf_measures, predictions, y_test)
  
  # Calculate confusion matrix elements (TP, FP, TN, FN) and derived metrics
  confusion_df <- confusion_matrix_wrapper(predictions, y_test)
  
  # Combine performance metrics and confusion matrix metrics
  combined_perf <- rbind(perf_df, confusion_df)
  
  list(
    predictions = predictions,
    performance = combined_perf
  )
}

#### CSV Output Format (Extracted from 1_simpipe.R) ####

#' Convert performance results to CSV format
#' @param performance_df Performance data frame from perf_wrapper
#' @param classifier Classifier type
#' @param mean_val Mean parameter value
#' @param var_val Variance parameter value
#' @param seed_val Seed value
create_output_dataframe <- function(performance_df, classifier, mean_val, var_val, seed_val) {
  
  # Extract performance metrics and methods
  methods <- colnames(performance_df)
  metrics <- rownames(performance_df)
  
  # Create long format data frame
  output_rows <- list()
  row_idx <- 1
  
  for(method in methods) {
    for(metric in metrics) {
      output_rows[[row_idx]] <- data.frame(
        classifier = classifier,
        mean = mean_val,
        variance = var_val,
        seed = seed_val,
        method = method,
        metric = metric,
        value = performance_df[metric, method],
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1
    }
  }
  
  # Combine all rows
  output_df <- do.call(rbind, output_rows)
  
  return(output_df)
}

#' Write results to CSV file
#' @param output_df Output data frame
#' @param output_path Path to output CSV file
write_results_csv <- function(output_df, output_path) {
  
  # Ensure output directory exists
  output_dir <- dirname(output_path)
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write CSV file
  write.csv(output_df, output_path, row.names = FALSE)
  
  cat(sprintf("Results written to: %s\n", output_path))
  cat(sprintf("Output contains %d rows with %d methods and %d metrics\n", 
             nrow(output_df), 
             length(unique(output_df$method)), 
             length(unique(output_df$metric))))
}

#### Error Handling and Logging ####

#' Main execution wrapper with comprehensive error handling
main_execution <- function() {
  
  # Create job identifier for logging
  job_id <- sprintf("batch_effects_%s_%s_%s_%s", params$classifier, params$mean, params$var, params$seed)
  
  tryCatch({
    
    # Log job start
    start_time <- Sys.time()
    cat(sprintf("[INFO] Job %s started at %s\n", job_id, start_time))
    cat(sprintf("Starting batch effects analysis: classifier=%s, mean=%s, var=%s, seed=%s\n", 
               params$classifier, params$mean, params$var, params$seed))
    
    # Load and prepare data
    cat("Loading and preparing data...\n")
    data_prep <- prepare_data()
    train_expr <- data_prep$train_expr
    test_expr <- data_prep$test_expr
    y_train <- data_prep$y_train
    y_test <- data_prep$y_test
    
    # Subset data into batches and filter genes
    cat("Subsetting data into batches and filtering genes...\n")
    N_batch <- 3
    N_sample_size <- 20  # Fixed sample size per batch (not tied to batch effect mean!)
    
    # Validate we have enough samples
    n_cases <- sum(y_train == 1)
    n_controls <- sum(y_train == 0)
    required_cases <- N_batch * (N_sample_size / 2)
    required_controls <- N_batch * (N_sample_size / 2)
    
    if(n_cases < required_cases || n_controls < required_controls) {
      cat(sprintf("Warning: Insufficient samples. Need %d cases + %d controls, have %d cases + %d controls\n",
                 required_cases, required_controls, n_cases, n_controls))
      # Reduce sample size if needed
      max_possible <- min(floor(n_cases / N_batch) * 2, floor(n_controls / N_batch) * 2)
      if(max_possible >= 4) {  # Minimum 2 cases + 2 controls per batch
        N_sample_size <- max_possible
        cat(sprintf("Reducing sample size to %d per batch\n", N_sample_size))
      } else {
        stop("Insufficient samples for meaningful analysis")
      }
    }
    
    subset_data <- subset_and_filter_data(train_expr, y_train, N_sample_size, N_batch)
    curr_train_expr <- subset_data$train_expr
    curr_y_train <- subset_data$y_train
    batch <- subset_data$batch
    batches_ind <- subset_data$batches_ind
    
    # Apply gene filter to test data
    curr_test_expr <- test_expr[subset_data$gene_filter, ]
    
    cat(sprintf("Data preparation complete: %d genes, %d training samples (%d per batch), %d test samples\n",
               nrow(curr_train_expr), ncol(curr_train_expr), N_sample_size, ncol(curr_test_expr)))
    
    # Generate hyperparameters for batch effect simulation
    cat("Generating hyperparameters for batch effect simulation...\n")
    hyper_pars <- generate_hyperparameters(max_batch_mean = params$mean, 
                                          max_batch_var = params$var, 
                                          N_batch = N_batch)
    
    cat(sprintf("Hyperparameters: mu range [%.2f, %.2f], var range [%.2f, %.2f]\n",
               min(hyper_pars$hyper_mu), max(hyper_pars$hyper_mu),
               min(hyper_pars$hyper_alpha/hyper_pars$hyper_beta), max(hyper_pars$hyper_alpha/hyper_pars$hyper_beta)))
    
    # Simulate batch effects on training data
    cat("Simulating batch effects...\n")
    train_expr_batch <- simulate_batch_effects(curr_train_expr, curr_y_train, batch, batches_ind, hyper_pars)
    
    cat(sprintf("Batch effect simulation complete: seed=%d\n", params$seed))
    
    # Normalize datasets
    cat("Normalizing datasets...\n")
    norm_data <- TRUE
    normalized_data <- normalize_datasets(curr_train_expr, curr_test_expr, train_expr_batch, 
                                         batches_ind, norm_data)
    
    # Apply batch corrections
    cat("Applying batch corrections (ComBat and MNN)...\n")
    corrected_data <- apply_batch_corrections(train_expr_batch, curr_test_expr, batch, curr_y_train)
    
    # Train and evaluate the specified classifier
    cat(sprintf("Training and evaluating %s classifier...\n", params$classifier))
    results <- train_and_evaluate_classifier(params$classifier, normalized_data, corrected_data, 
                                            curr_y_train, y_test)
    
    cat("Classifier training and evaluation complete\n")
    
    # Create output data frame
    cat("Creating output data frame...\n")
    output_df <- create_output_dataframe(results$performance, 
                                        params$classifier, 
                                        params$mean, 
                                        params$var, 
                                        params$seed)
    
    # Write results to CSV
    cat("Writing results to CSV...\n")
    write_results_csv(output_df, params$output)
    
    # Display summary of results
    cat("\nResults Summary:\n")
    print(output_df)
    
    # Success logging
    end_time <- Sys.time()
    duration <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
    cat(sprintf("[SUCCESS] Job %s completed successfully at %s\n", job_id, end_time))
    cat(sprintf("[SUCCESS] Total execution time: %.2f minutes\n", duration))
    cat(sprintf("[SUCCESS] Output file: %s\n", params$output))
    
    return(0)  # Success exit code
    
  }, error = function(e) {
    
    # Detailed error logging to stderr
    error_time <- Sys.time()
    
    cat(sprintf("[ERROR] Job %s failed at %s\n", job_id, error_time), file = stderr())
    cat(sprintf("[ERROR] Error message: %s\n", e$message), file = stderr())
    cat(sprintf("[ERROR] Error class: %s\n", class(e)[1]), file = stderr())
    
    # Log job parameters for debugging
    cat(sprintf("[ERROR] Job parameters:\n"), file = stderr())
    cat(sprintf("[ERROR]   Classifier: %s\n", params$classifier), file = stderr())
    cat(sprintf("[ERROR]   Mean: %s\n", params$mean), file = stderr())
    cat(sprintf("[ERROR]   Variance: %s\n", params$var), file = stderr())
    cat(sprintf("[ERROR]   Seed: %s\n", params$seed), file = stderr())
    cat(sprintf("[ERROR]   Output: %s\n", params$output), file = stderr())
    
    # Log system information for debugging
    cat(sprintf("[ERROR] Working directory: %s\n", getwd()), file = stderr())
    cat(sprintf("[ERROR] R version: %s\n", R.version.string), file = stderr())
    
    # Check input file existence
    input_files <- c(
      "data/combined_sub.RData",
      "scripts/helper.R"
    )
    
    cat(sprintf("[ERROR] Input file status:\n"), file = stderr())
    for(file in input_files) {
      exists <- file.exists(file)
      cat(sprintf("[ERROR]   %s: %s\n", file, ifelse(exists, "EXISTS", "MISSING")), file = stderr())
    }
    
    # Memory usage information
    gc_info <- gc()
    cat(sprintf("[ERROR] Memory usage: %.1f MB used, %.1f MB max\n", 
               sum(gc_info[, "used"]), sum(gc_info[, "max used"])), file = stderr())
    
    # Print traceback to stderr
    cat("[ERROR] Traceback:\n", file = stderr())
    traceback_lines <- capture.output(traceback())
    for(line in traceback_lines) {
      cat(sprintf("[ERROR] %s\n", line), file = stderr())
    }
    
    return(1)  # Error exit code
    
  })
}

# Execute main function and exit with appropriate code
exit_code <- main_execution()
quit(status = exit_code)