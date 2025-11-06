# Common Functions for Real Data Analysis Pipeline
# Enhanced with functional programming and improved error handling
# Preserves all scientific logic while improving code maintainability

options(warn = -1)
suppressMessages(suppressWarnings({
  # Enhanced package loading with functional programming support
  required_packages <- c("glmnet", "SummarizedExperiment", "sva", "DESeq2", "ROCR", "ggplot2", 
                        "gridExtra", "reshape2", "dplyr", "purrr", "nnls", "lightgbm", "batchelor")
  sapply(required_packages, require, character.only=TRUE, quietly=TRUE)
}))

####  Common Parameters and Setup  ####

#' Initialize analysis parameters and directories
#' @param n_studies Number of studies to include (3, 4, 5, or 6)
#' @param debug_mode Whether to run in debug mode with fewer iterations
#' @return List of parameters and setup information
initialize_analysis <- function(n_studies, debug_mode = FALSE) {
  # Create results directory
  results_dir <- sprintf("~/confounded_analysis/scripts/evaluations/robustifying/results_real_%dstudies", n_studies)
  if(!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  # Set parameters
  norm_data <- TRUE

  n_highvar_genes <- 1000
  B <- if(debug_mode) 3 else 100
  
  learner_types <- c("logistic", "elnet", "svm", "rf", "lightgbm", "nnet", "knn", "xgboost")
  perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc")
  perf_measures_names <- c("Mean cross-entropy loss", "AUC", "Root-mean-squared error", 
                           "F1 score", "Error rate", "Accuracy")
  names(perf_measures_names) <- perf_measures
  
  if(debug_mode) {
    cat("=== DEBUG MODE ENABLED ===\n")
    cat("Running with", B, "iterations instead of 100\n")
    cat("==========================\n\n")
  }
  
  list(
    results_dir = results_dir,
    norm_data = norm_data,

    n_highvar_genes = n_highvar_genes,
    B = B,
    learner_types = learner_types,
    perf_measures = perf_measures,
    perf_measures_names = perf_measures_names
  )
}

#' Filter studies based on analysis type
#' @param dat_lst List of data matrices
#' @param label_lst List of labels
#' @param n_studies Number of studies to include
#' @return List with filtered data and labels
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

####  Data Processing Functions  ####

#' Prepare training and test data
prepare_datasets <- function(dat_lst, label_lst, test_name, study_names) {
  train_name <- setdiff(study_names, test_name)
  
  dat <- do.call(cbind, dat_lst[train_name])
  batch <- rep(1:length(train_name), times=map_int(dat_lst[train_name], ncol))
  batches_ind <- map(1:length(train_name), ~which(batch == .x))
  batch_names <- levels(factor(batch))
  group <- do.call(c, label_lst[train_name])
  y_sgbatch_train <- map(batch_names, ~group[batch == .x])
  
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

#' Apply batch correction methods
apply_batch_corrections <- function(dat, batch, group, dat_test) {
  # ComBat with reference batch (use first batch as reference)
  ref_batch <- min(batch)
  
  # Apply ComBat to combined training and test data
  combined_dat <- cbind(dat, dat_test)
  combined_batch <- c(batch, rep(ref_batch, ncol(dat_test)))  # Test set uses reference batch
  combined_labels <- c(group, rep(0, ncol(dat_test)))  # Use dummy labels for test
  combat_combined <- ComBat(combined_dat, batch=combined_batch, mod=model.matrix(~combined_labels), ref.batch=ref_batch)
  
  # Split back into training and test
  dat_combat <- combat_combined[, 1:ncol(dat)]
  dat_test_combat <- combat_combined[, (ncol(dat) + 1):ncol(combat_combined)]
  
  # MNNcorrect - merge order with test set last
  # Convert to SingleCellExperiment format for batchelor
  library(batchelor, quietly = TRUE)
  
  # Combine training and test data for MNN correction
  combined_dat <- cbind(dat, dat_test)
  combined_batch <- c(batch, rep(max(batch) + 1, ncol(dat_test)))  # Test set gets new batch ID
  
  # Create batch list for MNN (each batch as separate matrix)
  unique_batches <- sort(unique(combined_batch))
  batch_list <- map(unique_batches, ~{
    combined_dat[, combined_batch == .x]
  })
  
  # Apply MNN correction with test set last in merge order
  mnn_result <- do.call(fastMNN, c(batch_list, list(merge.order = seq_along(unique_batches))))
  
  # Extract corrected data
  corrected_combined <- assay(mnn_result, "corrected")
  
  # Split back into training and test
  n_train <- ncol(dat)
  dat_mnn <- corrected_combined[, 1:n_train]
  dat_test_mnn <- corrected_combined[, (n_train + 1):ncol(corrected_combined)]
  
  list(dat_combat=dat_combat, dat_test_combat=dat_test_combat, dat_mnn=dat_mnn, dat_test_mnn=dat_test_mnn)
}

#' Normalize data matrices within batches (enhanced with functional programming)
normalize_within_batches <- function(dat, batch, batch_names) {
  dat_whole_norm <- normalizeData(dat)
  
  # Enhanced: Use map instead of for loop
  batch_normalized_data <- map(batch_names, ~{
    batch_indices <- batch == .x
    normalizeData(dat[, batch_indices])
  })
  
  # Reconstruct the full matrix
  dat_batch_norm <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat), dimnames=dimnames(dat))
  iwalk(batch_names, ~{
    batch_indices <- batch == .x
    dat_batch_norm[, batch_indices] <<- batch_normalized_data[[.y]]
  })
  
  list(whole_norm=dat_whole_norm, batch_norm=dat_batch_norm)
}

####  Model Training Functions  ####

#' Train model for single learner type
train_single_learner <- function(l_type, dat_batch_whole_norm, dat_combat_whole_norm, 
                                 dat_mnn_whole_norm, dat_batch_norm, group, batch, batch_names) {
  learner_fit <- getPredFunctions(l_type)
  
  # Define training configurations
  training_configs <- list(
    unadj = list(data = dat_batch_whole_norm, name = "unadj"),
    combat = list(data = dat_combat_whole_norm, name = "combat"),
    mnn = list(data = dat_mnn_whole_norm, name = "mnn")
  )
  
  # Train models using map instead of repetitive code
  trained_models <- map(training_configs, ~{
    trainPipe(train_set = .x$data, train_label = group, 
              test_set = NULL, lfit = learner_fit)
  })
  
  # Train single-batch models
  sgbatch_models <- map(batch_names, ~{
    batch_indices <- batch == .x
    trainPipe(train_set = dat_batch_norm[, batch_indices], 
              train_label = group[batch_indices], 
              test_set = NULL, lfit = learner_fit)
  })
  names(sgbatch_models) <- paste0("Batch", batch_names)
  
  result <- list(
    unadj_mod = trained_models$unadj$mod,
    combat_mod = trained_models$combat$mod, 
    mnn_mod = trained_models$mnn$mod,
    sgbatch_mod = map(sgbatch_models, "mod")
  )
  
  result
}

#' Train models for all learner types (enhanced with functional programming)
train_all_learners <- function(learner_types, dat_batch_whole_norm, dat_combat_whole_norm, 
                               dat_mnn_whole_norm, dat_batch_norm, group, batch, batch_names) {
  
  # Use map instead of for loop
  all_trained <- map(learner_types, ~{
    train_single_learner(.x, dat_batch_whole_norm, dat_combat_whole_norm,
                        dat_mnn_whole_norm, dat_batch_norm, group, batch, batch_names)
  })
  names(all_trained) <- learner_types
  
  # Reorganize results by model type instead of learner type
  list(
    unadj_mod_lst = map(all_trained, "unadj_mod"),
    combat_mod_lst = map(all_trained, "combat_mod"),
    mnn_mod_lst = map(all_trained, "mnn_mod"),
    sgbatch_mod_lst = map(all_trained, "sgbatch_mod")
  )
}

####  Prediction Functions  ####

#' Generate predictions for single learner type in bootstrap iteration
generate_predictions_single_learner <- function(l_type, unadj_mod_lst, combat_mod_lst, 
                                                mnn_mod_lst, sgbatch_mod_lst, dat_test_norm,
                                                perf_measures, group_test) {
  
  # Define prediction configurations
  pred_configs <- list(
    Unadjusted = list(model = unadj_mod_lst[[l_type]]),
    ComBat = list(model = combat_mod_lst[[l_type]]),
    MNNcorrect = list(model = mnn_mod_lst[[l_type]])
  )
  
  # Generate predictions using map
  predictions <- map(pred_configs, ~{
    predWrapper(.x$model, dat_test_norm, l_type)
  })
  
  # Calculate performance metrics
  perf_df <- perf_wrapper(perf_measures, predictions, group_test)
  
  list(
    perf_df = perf_df, 
    tst_scores = predictions,
    unadj_tst_prob = predictions$Unadjusted,
    combat_tst_prob = predictions$ComBat, 
    mnn_tst_prob = predictions$MNNcorrect
  )
}

#' Generate predictions for all learner types in bootstrap iteration
generate_predictions <- function(learner_types, unadj_mod_lst, combat_mod_lst, 
                                mnn_mod_lst, sgbatch_mod_lst, dat_test_norm,
                                perf_measures, group_test) {
  
  # Use map instead of for loop
  all_predictions <- map(learner_types, ~{
    generate_predictions_single_learner(.x, unadj_mod_lst, combat_mod_lst,
                                       mnn_mod_lst, sgbatch_mod_lst, dat_test_norm,
                                       perf_measures, group_test)
  })
  names(all_predictions) <- learner_types
  
  # Extract components
  list(
    perf_df_lst = map(all_predictions, "perf_df"),
    tst_scores_modlst = map(all_predictions, "tst_scores"),
    unadj_tst_prob = all_predictions[[1]]$unadj_tst_prob,  # Take from first learner
    combat_tst_prob = all_predictions[[1]]$combat_tst_prob,
    mnn_tst_prob = all_predictions[[1]]$mnn_tst_prob
  )
}

#' Generate predictions for single learner type with method-specific test sets
generate_predictions_single_learner_specific_tests <- function(l_type, unadj_mod_lst, combat_mod_lst, 
                                                              mnn_mod_lst, sgbatch_mod_lst, 
                                                              dat_test_unadj, dat_test_combat, dat_test_mnn,
                                                              perf_measures, group_test) {
  
  # Define prediction configurations with specific test sets
  pred_configs <- list(
    Unadjusted = list(model = unadj_mod_lst[[l_type]], test_data = dat_test_unadj),
    ComBat = list(model = combat_mod_lst[[l_type]], test_data = dat_test_combat),
    MNNcorrect = list(model = mnn_mod_lst[[l_type]], test_data = dat_test_mnn)
  )
  
  # Generate predictions using map with specific test sets
  predictions <- map(pred_configs, ~{
    predWrapper(.x$model, .x$test_data, l_type)
  })
  
  # Calculate performance metrics
  perf_df <- map_dfr(names(predictions), ~{
    pred_prob <- predictions[[.x]]
    perf_vals <- map_dbl(perf_measures, ~perfWrapper(pred_prob, group_test, .x))
    names(perf_vals) <- perf_measures
    data.frame(Method = .x, t(perf_vals), stringsAsFactors = FALSE)
  })
  
  list(
    perf_df = perf_df,
    tst_scores = predictions,
    unadj_tst_prob = predictions$Unadjusted,
    combat_tst_prob = predictions$ComBat, 
    mnn_tst_prob = predictions$MNNcorrect
  )
}

#' Generate predictions with method-specific test sets
generate_predictions_with_specific_test_sets <- function(learner_types, unadj_mod_lst, combat_mod_lst, 
                                                        mnn_mod_lst, sgbatch_mod_lst, 
                                                        dat_test_unadj, dat_test_combat, dat_test_mnn,
                                                        perf_measures, group_test) {
  
  # Use map instead of for loop
  all_predictions <- map(learner_types, ~{
    generate_predictions_single_learner_specific_tests(.x, unadj_mod_lst, combat_mod_lst,
                                                      mnn_mod_lst, sgbatch_mod_lst, 
                                                      dat_test_unadj, dat_test_combat, dat_test_mnn,
                                                      perf_measures, group_test)
  })
  names(all_predictions) <- learner_types
  
  # Extract components
  list(
    perf_df_lst = map(all_predictions, "perf_df"),
    tst_scores_modlst = map(all_predictions, "tst_scores"),
    unadj_tst_prob = all_predictions[[1]]$unadj_tst_prob,  # Take from first learner
    combat_tst_prob = all_predictions[[1]]$combat_tst_prob,
    mnn_tst_prob = all_predictions[[1]]$mnn_tst_prob
  )
}

####  Results Processing Functions  ####

#' Extract performance metrics by model and measure
extract_performance_data <- function(perf_df_lst, perf_measure, n_batches) {
  
  # Use map instead of for loop
  perf_data <- map(perf_df_lst, ~{
    perf_res <- .x
    if(perf_measure %in% rownames(perf_res)) {
      # Handle different model types
      if(names(perf_df_lst)[which(map_lgl(perf_df_lst, identical, .x))] == "crossmod") {
        perf_res[perf_measure, , drop=FALSE]
      } else {
        start_col <- 2 + n_batches
        if(ncol(perf_res) >= start_col) {
          selected_cols <- c(1, start_col:ncol(perf_res))
          perf_res[perf_measure, selected_cols]
        } else {
          perf_res[perf_measure, 1, drop=FALSE]
        }
      }
    } else {
      NULL
    }
  })
  
  # Remove NULL entries
  compact(perf_data)
}

#' Write performance results to CSV
write_performance_results <- function(perf_data, test_name, perf_measure, iteration, results_dir) {
  if(length(perf_data) == 0) return(NULL)
  
  summary_df <- reshape2::melt(perf_data) %>%
    mutate(iteration = iteration) %>%
    rename(Method = Var1, value = value, Model = L1, Iteration = iteration) %>%
    select(Method, value, Model, Iteration)
  
  # Ensure results directory exists
  if(!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  filepath <- file.path(results_dir, sprintf('test%s_%s.csv', test_name, perf_measure))
  first_file <- !file.exists(filepath)
  write.table(summary_df, filepath, append=!first_file, col.names=first_file,
              row.names=FALSE, sep=",")
  
  summary_df
}

####  Bootstrap Evaluation  ####

#' PRESERVED: Bootstrap evaluation with original retry logic
#' This maintains the exact statistical behavior of the original while loop
run_bootstrap_evaluation <- function(params, trained_models, dat_testOri, group_testOri, 
                                   test_name, n_batches) {
  
  # PRESERVED: Original while loop structure with retry logic
  # This is critical for maintaining exact bootstrap behavior
  b <- 1
  successful_iterations <- 0
  
  while(b <= params$B) {
    boot_ind <- sample(1:ncol(dat_testOri), ncol(dat_testOri), replace=TRUE)
    dat_test <- dat_testOri[, boot_ind]
    group_test <- group_testOri[boot_ind]
    
    # Create different test sets for different methods
    dat_test_combat_boot <- dat_test_combat[, boot_ind]
    dat_test_mnn_boot <- dat_test_mnn[, boot_ind]
    
    if(params$norm_data) {
      dat_test_unadj_norm <- normalizeData(dat_test)
      dat_test_combat_norm <- normalizeData(dat_test_combat_boot)
      dat_test_mnn_norm <- normalizeData(dat_test_mnn_boot)
    } else {
      dat_test_unadj_norm <- dat_test
      dat_test_combat_norm <- dat_test_combat_boot
      dat_test_mnn_norm <- dat_test_mnn_boot
    }
    
    # PRESERVED: Original retry logic - decrement b if normalization fails
    if(any(is.na(dat_test_unadj_norm)) || any(is.na(dat_test_combat_norm)) || any(is.na(dat_test_mnn_norm))) {
      b <- b - 1  # This is the critical retry logic from original
    } else {
      # Generate predictions with method-specific test sets
      single_preds <- generate_predictions_with_specific_test_sets(
        params$learner_types, trained_models$unadj_mod_lst,
        trained_models$combat_mod_lst, trained_models$mnn_mod_lst,
        trained_models$sgbatch_mod_lst, 
        dat_test_unadj_norm, dat_test_combat_norm, dat_test_mnn_norm,
        params$perf_measures, group_test
      )
      
      # Write results for each performance measure
      walk(params$perf_measures, ~{
        perf_data <- extract_performance_data(single_preds$perf_df_lst, .x, n_batches)
        write_performance_results(perf_data, test_name, .x, b, params$results_dir)
      })
      
      successful_iterations <- successful_iterations + 1
      
      # PRESERVED: Explicit cleanup from original
      rm(boot_ind, dat_test, group_test, single_preds, envir=environment())
    }
    
    b <- b + 1
  }
  
  cat(sprintf("Completed %d/%d bootstrap iterations for %s\n", 
              successful_iterations, params$B, test_name))
  
  successful_iterations
}

####  Main Pipeline Function  ####

#' Run complete analysis pipeline for specified number of studies
#' @param n_studies Number of studies (3, 4, 5, or 6)
#' @param dat_lst List of data matrices
#' @param label_lst List of labels  
#' @param debug_mode Whether to run in debug mode
run_analysis_pipeline <- function(n_studies, dat_lst, label_lst, debug_mode = FALSE) {
  # Initialize parameters and setup
  params <- initialize_analysis(n_studies, debug_mode)
  
  # Filter studies based on analysis type
  filtered_data <- filter_studies(dat_lst, label_lst, n_studies)
  dat_lst <- filtered_data$dat_lst
  label_lst <- filtered_data$label_lst
  study_names <- filtered_data$study_names
  
  # Main analysis using map instead of for loop
  study_results <- map(study_names, ~{
    test_name <- .x
    cat(sprintf("Processing study: %s\n", test_name))
    
    ## Prepare datasets
    datasets <- prepare_datasets(dat_lst, label_lst, test_name, study_names)
    
    ## Feature reduction
    feat_reduced <- reduce_features(datasets$dat, datasets$dat_test, params$n_highvar_genes)
    dat <- feat_reduced$dat
    dat_test <- feat_reduced$dat_test
    
    ## Batch correction
    batch_corr <- apply_batch_corrections(dat, datasets$batch, datasets$group, dat_test)
    dat_combat <- batch_corr$dat_combat
    dat_test_combat <- batch_corr$dat_test_combat
    dat_mnn <- batch_corr$dat_mnn
    dat_test_mnn <- batch_corr$dat_test_mnn
    
    ## Normalize data
    if(params$norm_data) {
      norm_raw <- normalize_within_batches(dat, datasets$batch, datasets$batch_names)
      norm_combat <- normalize_within_batches(dat_combat, datasets$batch, datasets$batch_names)
      norm_mnn <- normalize_within_batches(dat_mnn, datasets$batch, datasets$batch_names)
    } else {
      norm_raw <- list(whole_norm=dat, batch_norm=dat)
      norm_combat <- list(whole_norm=dat_combat, batch_norm=dat_combat)
      norm_mnn <- list(whole_norm=dat_mnn, batch_norm=dat_mnn)
    }
    
    ## Train all models
    trained_models <- train_all_learners(
      params$learner_types, norm_raw$whole_norm, norm_combat$whole_norm, 
      norm_mnn$whole_norm, norm_raw$batch_norm, datasets$group, datasets$batch,
      datasets$batch_names
    )
    
    ## Run bootstrap evaluation
    n_batches <- length(unique(datasets$batch))
    successful_iterations <- run_bootstrap_evaluation(
      params, trained_models, datasets$dat_test, datasets$group_test, 
      test_name, n_batches
    )
    
    list(study = test_name, successful_iterations = successful_iterations)
  })
  
  names(study_results) <- study_names
  
  # Summary
  total_successful <- sum(map_int(study_results, "successful_iterations"))
  total_expected <- length(study_names) * params$B
  
  cat(sprintf("\n%d-study analysis completed successfully\n", n_studies))
  cat(sprintf("Total successful bootstrap iterations: %d/%d\n", total_successful, total_expected))
  
  study_results
}