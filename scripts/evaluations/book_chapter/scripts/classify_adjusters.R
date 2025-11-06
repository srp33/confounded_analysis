#!/usr/bin/env Rscript

# classify_adjusters.R - Single job adjuster comparison script
# Executes single adjuster × classifier × dataset × seed combination

# Suppress warnings and messages for cleaner output
options(warn = -1)
suppressMessages(suppressWarnings({
  rm(list=ls())
}))

# Load required libraries
suppressMessages(suppressWarnings({
  required_packages <- c("glmnet", "SummarizedExperiment", "sva", "DESeq2", "ROCR", "ggplot2", 
                        "gridExtra", "reshape2", "dplyr", "purrr", "nnls", "lightgbm", "batchelor",
                        "argparse")
  sapply(required_packages, require, character.only=TRUE, quietly=TRUE)
}))

# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Execute single adjuster comparison job for batch correction analysis")

parser$add_argument("--adjuster", type = "character", required = TRUE,
                   help = "Batch correction method: unadjusted, combat, or mnn")
parser$add_argument("--classifier", type = "character", required = TRUE,
                   help = "Classifier type: logistic, elnet, svm, rf, lightgbm, or nnet")
parser$add_argument("--num-datasets", type = "integer", required = TRUE,
                   help = "Number of datasets to include: 3, 4, 5, or 6")
parser$add_argument("--seed", type = "integer", required = TRUE,
                   help = "Random seed for reproducibility")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                   help = "Output CSV file path")

# Parse arguments
args <- parser$parse_args()

# Arguments are automatically validated as required by argparse

# Parameter validation
valid_adjusters <- c("unadjusted", "combat", "mnn")
valid_classifiers <- c("logistic", "elnet", "svm", "rf", "lightgbm", "nnet")
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

if (args$seed < 1 || args$seed > 1000) {
  cat("Error: Seed must be between 1 and 1000\n")
  quit(status=1)
}

# Extract validated parameters
adjuster <- args$adjuster
classifier <- args$classifier
num_datasets <- args$num_datasets
seed <- args$seed
output_file <- args$output

# Validate output directory exists
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  cat(sprintf("Error: Output directory does not exist: %s\n", output_dir))
  quit(status=1)
}

# ====================================================================
# ERROR HANDLING AND LOGGING WRAPPER
# ====================================================================

# Create job ID for logging
job_id <- sprintf("adjuster_%s_%s_%d_%d", adjuster, classifier, num_datasets, seed)

# Main job wrapper with comprehensive error handling
main_job_wrapper <- function() {
  tryCatch({
    # Print job parameters for logging
    cat("=== ADJUSTER COMPARISON JOB ===\n")
    cat(sprintf("Job ID: %s\n", job_id))
    cat(sprintf("Adjuster: %s\n", adjuster))
    cat(sprintf("Classifier: %s\n", classifier))
    cat(sprintf("Num datasets: %d\n", num_datasets))
    cat(sprintf("Seed: %d\n", seed))
    cat(sprintf("Output: %s\n", output_file))
    cat(sprintf("Start time: %s\n", Sys.time()))
    cat("===============================\n\n")
    
    # Execute main analysis
    result <- main_analysis_function()
    
    # Success logging
    cat(sprintf("\n[SUCCESS] Job %s completed at %s\n", job_id, Sys.time()))
    
    return(result)
    
  }, error = function(e) {
    # Detailed error logging
    cat(sprintf("[ERROR] Job %s failed at %s\n", job_id, Sys.time()), file = stderr())
    cat(sprintf("[ERROR] Error: %s\n", e$message), file = stderr())
    cat(sprintf("[ERROR] Parameters: adjuster=%s, classifier=%s, num_datasets=%d, seed=%d\n", 
                adjuster, classifier, num_datasets, seed), file = stderr())
    
    # Additional debugging information
    cat(sprintf("[ERROR] Working directory: %s\n", getwd()), file = stderr())
    cat(sprintf("[ERROR] Output file: %s\n", output_file), file = stderr())
    
    # Check input files
    data_path <- "/scripts/evaluations/book_chapter/data/TB_real_data.RData"
    helper_path <- "/scripts/evaluations/book_chapter/scripts/helper.R"
    common_path <- "/scripts/evaluations/book_chapter/scripts/common_functions.R"
    
    cat(sprintf("[ERROR] Data file exists: %s\n", file.exists(data_path)), file = stderr())
    cat(sprintf("[ERROR] Helper file exists: %s\n", file.exists(helper_path)), file = stderr())
    cat(sprintf("[ERROR] Common functions file exists: %s\n", file.exists(common_path)), file = stderr())
    
    # Memory usage
    gc_info <- capture.output(gc())
    cat(sprintf("[ERROR] Memory usage: %s\n", paste(gc_info, collapse="; ")), file = stderr())
    
    # Session info
    session_info <- capture.output(sessionInfo())
    cat("[ERROR] Session info:\n", file = stderr())
    cat(paste(session_info[1:5], collapse="\n"), file = stderr())
    cat("\n", file = stderr())
    
    # Exit with error code for Snakemake
    quit(status = 1)
  })
}

# ====================================================================
# MAIN ANALYSIS FUNCTION
# ====================================================================

main_analysis_function <- function() {
  # Load data and dependencies
  data_path <- "/scripts/evaluations/book_chapter/data/TB_real_data.RData"
  if (!file.exists(data_path)) {
    stop(sprintf("Data file not found: %s", data_path))
  }
  
  load(data_path)
  source("/scripts/evaluations/book_chapter/scripts/helper.R")
  source("/scripts/evaluations/book_chapter/scripts/common_functions.R")
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # ====================================================================
  # REAL DATA PREPARATION LOGIC (extracted from common_functions.R)
  # ====================================================================
  
  #' Filter studies based on analysis type
  filter_studies <- function(dat_lst, label_lst, n_studies) {
    all_studies <- c("GSE37250_SA", "US", "India", "GSE37250_M", "Africa", "GSE39941_M")
    selected_studies <- all_studies[1:n_studies]
    
    # Filter data and labels to keep only selected studies
    dat_lst <- dat_lst[selected_studies]
    label_lst <- label_lst[selected_studies]
    study_names <- names(dat_lst)
    cat(sprintf("Running %d-study analysis with studies: %s\n", 
                n_studies, paste(study_names, collapse=", ")))
    
    list(dat_lst = dat_lst, label_lst = label_lst, study_names = study_names)
  }
  
  #' Prepare training and test data
  prepare_datasets <- function(dat_lst, label_lst, test_name, study_names) {
    train_name <- setdiff(study_names, test_name)
    
    dat <- do.call(cbind, dat_lst[train_name])
    batch <- rep(1:length(train_name), times=sapply(dat_lst[train_name], ncol))
    batches_ind <- lapply(1:length(train_name), function(x) which(batch == x))
    batch_names <- levels(factor(batch))
    group <- do.call(c, label_lst[train_name])
    y_sgbatch_train <- lapply(batch_names, function(x) group[batch == x])
    
    dat_test <- dat_lst[[test_name]]
    group_test <- label_lst[[test_name]]
    
    list(dat=dat, batch=batch, batches_ind=batches_ind, batch_names=batch_names, 
         group=group, y_sgbatch_train=y_sgbatch_train, 
         dat_test=dat_test, group_test=group_test)
  }
  
  #' Select highly variable genes and reduce feature space
  reduce_features <- function(dat, dat_test, n_genes=1000) {
    genes_sel_names <- order(rowVars(dat), decreasing=TRUE)[1:n_genes]
    list(dat=dat[genes_sel_names, ], 
         dat_test=dat_test[genes_sel_names, ])
  }
  
  # ====================================================================
  # EXECUTE DATA PREPARATION
  # ====================================================================
  
  cat("Starting data preparation...\n")
  
  # Filter studies based on num_datasets parameter
  filtered_data <- filter_studies(dat_lst, label_lst, num_datasets)
  dat_lst_filtered <- filtered_data$dat_lst
  label_lst_filtered <- filtered_data$label_lst
  study_names <- filtered_data$study_names
  
  # For single job execution, we need to select one test study
  # Use seed to deterministically select test study
  set.seed(seed)
  test_study_index <- ((seed - 1) %% length(study_names)) + 1
  test_name <- study_names[test_study_index]
  
  cat(sprintf("Selected test study: %s (index %d based on seed %d)\n", 
              test_name, test_study_index, seed))
  
  # Prepare datasets
  datasets <- prepare_datasets(dat_lst_filtered, label_lst_filtered, test_name, study_names)
  
  # Validate datasets
  if(is.null(datasets$dat_test)) {
    stop(sprintf("Test dataset '%s' is NULL or missing from data", test_name))
  }
  if(ncol(datasets$dat_test) == 0) {
    stop(sprintf("Test dataset '%s' has no samples", test_name))
  }
  
  # Feature reduction (top 1000 most variable genes)
  n_highvar_genes <- 1000
  feat_reduced <- reduce_features(datasets$dat, datasets$dat_test, n_highvar_genes)
  dat <- feat_reduced$dat
  dat_test <- feat_reduced$dat_test
  
  # Validate feature-reduced data
  if(is.null(dat_test)) {
    stop("Test data became NULL after feature reduction")
  }
  
  cat(sprintf("Data preparation completed:\n"))
  cat(sprintf("  Training samples: %d\n", ncol(dat)))
  cat(sprintf("  Test samples: %d\n", ncol(dat_test)))
  cat(sprintf("  Features (genes): %d\n", nrow(dat)))
  cat(sprintf("  Training batches: %d\n", length(unique(datasets$batch))))
  
  # ====================================================================
  # BATCH CORRECTION APPLICATION LOGIC
  # ====================================================================
  
  #' Apply batch correction methods with sophisticated test set handling
  apply_batch_corrections <- function(dat, batch, group, dat_test, method) {
    if (method == "unadjusted") {
      # No correction - return original data
      return(list(
        dat_corrected = dat,
        dat_test_corrected = dat_test
      ))
    } else if (method == "combat") {
      # ComBat with reference batch (use first batch as reference)
      ref_batch <- min(batch)
      
      # Apply ComBat to combined training and test data
      combined_dat <- cbind(dat, dat_test)
      combined_batch <- c(batch, rep(ref_batch, ncol(dat_test)))  # Test set uses reference batch
      combined_labels <- c(group, rep(0, ncol(dat_test)))  # Use dummy labels for test
      
      # Apply ComBat correction
      combat_combined <- ComBat(combined_dat, batch=combined_batch, 
                               mod=model.matrix(~combined_labels), ref.batch=ref_batch)
      
      # Split back into training and test
      dat_corrected <- combat_combined[, 1:ncol(dat)]
      dat_test_corrected <- combat_combined[, (ncol(dat) + 1):ncol(combat_combined)]
      
      return(list(
        dat_corrected = dat_corrected,
        dat_test_corrected = dat_test_corrected
      ))
    } else if (method == "mnn") {
      # MNNcorrect - merge order with test set last
      library(batchelor, quietly = TRUE)
      
      # Combine training and test data for MNN correction
      combined_dat <- cbind(dat, dat_test)
      combined_batch <- c(batch, rep(max(batch) + 1, ncol(dat_test)))  # Test set gets new batch ID
      
      # Create batch list for MNN (each batch as separate matrix)
      unique_batches <- sort(unique(combined_batch))
      batch_list <- lapply(unique_batches, function(b) {
        combined_dat[, combined_batch == b]
      })
      
      # Apply MNN correction with test set last in merge order
      mnn_result <- do.call(fastMNN, c(batch_list, list(merge.order = seq_along(unique_batches))))
      
      # Extract corrected data (handle different batchelor versions)
      if("corrected" %in% assayNames(mnn_result)) {
        corrected_combined <- assay(mnn_result, "corrected")
      } else if("reconstructed" %in% assayNames(mnn_result)) {
        corrected_combined <- assay(mnn_result, "reconstructed")
      } else {
        # Fallback to first assay
        corrected_combined <- assay(mnn_result, 1)
      }
      
      # Split back into training and test
      n_train <- ncol(dat)
      dat_corrected <- corrected_combined[, 1:n_train]
      dat_test_corrected <- corrected_combined[, (n_train + 1):ncol(corrected_combined)]
      
      return(list(
        dat_corrected = dat_corrected,
        dat_test_corrected = dat_test_corrected
      ))
    } else {
      stop(sprintf("Unknown batch correction method: %s", method))
    }
  }
  
  #' Normalize data matrices within batches
  normalize_within_batches <- function(dat, batch, batch_names) {
    dat_whole_norm <- normalizeData(dat)
    
    # Normalize each batch separately
    batch_normalized_data <- lapply(batch_names, function(b) {
      batch_indices <- batch == b
      normalizeData(dat[, batch_indices])
    })
    
    # Reconstruct the full matrix
    dat_batch_norm <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat), dimnames=dimnames(dat))
    for (i in seq_along(batch_names)) {
      batch_indices <- batch == batch_names[i]
      dat_batch_norm[, batch_indices] <- batch_normalized_data[[i]]
    }
    
    list(whole_norm=dat_whole_norm, batch_norm=dat_batch_norm)
  }
  
  # ====================================================================
  # EXECUTE BATCH CORRECTION
  # ====================================================================
  
  cat(sprintf("Applying batch correction method: %s\n", adjuster))
  
  # Apply batch correction
  batch_corr_result <- apply_batch_corrections(dat, datasets$batch, datasets$group, dat_test, adjuster)
  dat_corrected <- batch_corr_result$dat_corrected
  dat_test_corrected <- batch_corr_result$dat_test_corrected
  
  # Normalize data
  norm_data <- TRUE  # Always normalize as in original implementation
  
  if (norm_data) {
    # Normalize training data
    norm_result <- normalize_within_batches(dat_corrected, datasets$batch, datasets$batch_names)
    dat_train_norm <- norm_result$whole_norm
    
    # Normalize test data with error handling
    cat(sprintf("Debug - About to normalize test data. dat_test_corrected is NULL: %s\n", is.null(dat_test_corrected)))
    
    if(is.null(dat_test_corrected)) {
      stop("Test data is NULL before normalization - batch correction failed")
    }
    
    dat_test_norm <- normalizeData(dat_test_corrected)
    
    cat(sprintf("Debug - After normalization. dat_test_norm is NULL: %s\n", is.null(dat_test_norm)))
    if(!is.null(dat_test_norm)) {
      cat(sprintf("Debug - dat_test_norm dimensions after normalization: %d x %d\n", nrow(dat_test_norm), ncol(dat_test_norm)))
    }
    
    # Check for normalization issues
    if (any(is.na(dat_train_norm)) || any(is.na(dat_test_norm))) {
      stop("Normalization produced NA values")
    }
    if (is.null(dat_test_norm)) {
      stop("Test data became NULL during normalization")
    }
  } else {
    dat_train_norm <- dat_corrected
    dat_test_norm <- dat_test_corrected
  }
  
  cat(sprintf("Batch correction completed successfully\n"))
  cat(sprintf("  Method: %s\n", adjuster))
  cat(sprintf("  Training data shape: %d x %d\n", nrow(dat_train_norm), ncol(dat_train_norm)))
  cat(sprintf("  Test data shape: %d x %d\n", nrow(dat_test_norm), ncol(dat_test_norm)))
  
  # ====================================================================
  # SINGLE CLASSIFIER TRAINING AND EVALUATION
  # ====================================================================
  
  #' Train single classifier and evaluate performance
  train_and_evaluate_classifier <- function(classifier_type, train_data, train_labels, test_data, test_labels) {
    # Get prediction function for classifier type
    learner_fit <- getPredFunctions(classifier_type)
    
    # Train model
    cat(sprintf("Training %s classifier...\n", classifier_type))
    trained_model <- trainPipe(train_set = train_data, train_label = train_labels, 
                              test_set = NULL, lfit = learner_fit)
    
    # Generate predictions on test set
    cat(sprintf("Generating predictions...\n"))
    test_predictions <- predWrapper(trained_model$mod, test_data, classifier_type)
    
    # Calculate performance metrics
    perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc")
    
    # Create predictions list in format expected by perf_wrapper
    predictions_list <- list()
    predictions_list[[adjuster]] <- test_predictions
    
    # Calculate performance using original perf_wrapper function
    perf_results <- perf_wrapper(perf_measures, predictions_list, test_labels)
    
    # Extract performance values for this method
    perf_values <- perf_results[, adjuster]
    names(perf_values) <- rownames(perf_results)
    
    return(list(
      model = trained_model,
      predictions = test_predictions,
      performance = perf_values
    ))
  }
  
  # ====================================================================
  # EXECUTE CLASSIFIER TRAINING AND EVALUATION
  # ====================================================================
  
  cat(sprintf("Training and evaluating classifier: %s\n", classifier))
  
  # Debug: Check data before passing to classifier
  cat(sprintf("Debug - dat_train_norm is NULL: %s\n", is.null(dat_train_norm)))
  cat(sprintf("Debug - dat_test_norm is NULL: %s\n", is.null(dat_test_norm)))
  cat(sprintf("Debug - dat_test_corrected is NULL: %s\n", is.null(dat_test_corrected)))
  
  if(!is.null(dat_test_norm)) {
    cat(sprintf("Debug - dat_test_norm dimensions: %d x %d\n", nrow(dat_test_norm), ncol(dat_test_norm)))
  }
  if(!is.null(dat_test_corrected)) {
    cat(sprintf("Debug - dat_test_corrected dimensions: %d x %d\n", nrow(dat_test_corrected), ncol(dat_test_corrected)))
  }
  
  # Additional validation before classifier training
  if(is.null(dat_test_norm)) {
    stop("Test data is NULL after normalization - this should not happen")
  }
  if(ncol(dat_test_norm) == 0) {
    stop("Test data has zero columns after normalization")
  }
  if(nrow(dat_test_norm) == 0) {
    stop("Test data has zero rows after normalization")
  }
  
  # Classifier-specific early validation
  n_train_samples <- ncol(dat_train_norm)
  n_test_samples <- ncol(dat_test_norm)
  n_features <- nrow(dat_train_norm)
  
  cat(sprintf("Dataset summary before classifier training:\n"))
  cat(sprintf("  Training samples: %d\n", n_train_samples))
  cat(sprintf("  Test samples: %d\n", n_test_samples))
  cat(sprintf("  Features: %d\n", n_features))
  cat(sprintf("  Training labels: %d unique values\n", length(unique(datasets$group))))
  
  # Early validation for problematic cases
  if(n_train_samples < 10) {
    warning(sprintf("Very small training set (%d samples) - results may be unreliable", n_train_samples))
  }
  
  if(classifier == "lightgbm" && n_train_samples < 20) {
    warning(sprintf("LightGBM with small dataset (%d samples) - using simplified configuration", n_train_samples))
  }
  
  # Check for class imbalance
  class_counts <- table(datasets$group)
  if(min(class_counts) < 3) {
    warning(sprintf("Severe class imbalance detected: %s", paste(names(class_counts), class_counts, sep="=", collapse=", ")))
  }
  
  # Train classifier and evaluate performance
  result <- train_and_evaluate_classifier(
    classifier_type = classifier,
    train_data = dat_train_norm,
    train_labels = datasets$group,
    test_data = dat_test_norm,
    test_labels = datasets$group_test
  )
  
  cat("Classification completed successfully\n")
  cat(sprintf("Performance metrics for %s + %s:\n", adjuster, classifier))
  for (metric in names(result$performance)) {
    cat(sprintf("  %s: %.6f\n", metric, result$performance[metric]))
  }
  
  # ====================================================================
  # CSV OUTPUT FORMAT IMPLEMENTATION
  # ====================================================================
  
  #' Create output data frame with required columns
  create_output_dataframe <- function(adjuster, classifier, n_datasets, seed, performance_metrics) {
    # Create one row per metric
    output_rows <- lapply(names(performance_metrics), function(metric) {
      data.frame(
        adjuster = adjuster,
        classifier = classifier,
        n_datasets = n_datasets,
        seed = seed,
        metric = metric,
        value = performance_metrics[metric],
        stringsAsFactors = FALSE
      )
    })
    
    # Combine all rows
    output_df <- do.call(rbind, output_rows)
    return(output_df)
  }
  
  # ====================================================================
  # GENERATE AND WRITE OUTPUT
  # ====================================================================
  
  cat("Generating output CSV...\n")
  
  # Create output data frame
  output_df <- create_output_dataframe(
    adjuster = adjuster,
    classifier = classifier,
    n_datasets = num_datasets,
    seed = seed,
    performance_metrics = result$performance
  )
  
  # Write to CSV file
  write.csv(output_df, file = output_file, row.names = FALSE)
  
  cat(sprintf("Results written to: %s\n", output_file))
  cat(sprintf("Output contains %d rows (one per metric)\n", nrow(output_df)))
  
  # Display output for verification
  cat("\nOutput preview:\n")
  print(output_df)
  
  cat(sprintf("\n=== JOB COMPLETED SUCCESSFULLY ===\n"))
  cat(sprintf("Adjuster: %s\n", adjuster))
  cat(sprintf("Classifier: %s\n", classifier))
  cat(sprintf("Datasets: %d\n", num_datasets))
  cat(sprintf("Seed: %d\n", seed))
  cat(sprintf("Output: %s\n", output_file))
  cat("==================================\n")
  
  return(output_df)
}

# ====================================================================
# EXECUTE MAIN JOB
# ====================================================================

# Run the main job with error handling
main_job_wrapper()