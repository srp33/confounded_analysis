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
  # REMOVE "reticulate" from this list
  required_packages <- c("glmnet", "SummarizedExperiment", "sva", "DESeq2", 
                        "ROCR", "ggplot2", "gridExtra", "reshape2", 
                        "dplyr", "purrr", "nnls", "batchelor",
                        "argparse", "class", "xgboost")
  sapply(required_packages, require, character.only=TRUE, quietly=TRUE)
}))

# ====================================================================
# [RVC MODIFICATION] RETICULATE / RVC SETUP
# ====================================================================
# This section attempts to load the necessary Python modules for RVC.
# We define them as NULL globally and use <<- to assign them if found.
RVC_py <- NULL
np_py <- NULL

import_reticulate <- function() {
  tryCatch({
    # LOAD THE LIBRARY HERE
    library(reticulate) 
    
    cat("Attempting to import Python modules for RVC...\n")
    rvm_module <- import("sklearn_rvm")
    RVC_py <<- rvm_module$em_rvm$EMRVC
    np_py <<- import("numpy")
    cat("Successfully imported sklearn_rvm and numpy.\n")
  }, error = function(e) {
    cat("[WARNING] Could not import Python modules 'sklearn_rvm' or 'numpy'.\n")
    cat("[WARNING] The 'rvc' classifier will be unavailable.\n")
    cat(sprintf("[WARNING] Python Error: %s\n", e$message))
  })
}

# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Execute single adjuster comparison job for batch correction analysis")

parser$add_argument("--adjuster", type = "character", required = TRUE,
                   help = "Batch correction method: unadjusted, combat, combat_mean, combat_sup, mnn, mnn_centered, ruvr, ruvg, or gmm")
parser$add_argument("--classifier", type = "character", required = TRUE,
                   help = "Classifier type: logistic, elnet, elasticnet, svm, rf, nnet, knn, xgboost, or rvc")
parser$add_argument("--num-datasets", type = "integer", required = TRUE,
                   help = "Number of datasets to include: 3, 4, 5, or 6")
parser$add_argument("--test-study", type = "character", required = TRUE,
                   help = "Test study name (e.g., GSE37250_SA, USA, India, etc.)")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                   help = "Output CSV file path")

# Parse arguments
args <- parser$parse_args()

# Arguments are automatically validated as required by argparse

# Parameter validation
valid_adjusters <- c("unadjusted", "combat", "combat_mean", "combat_sup", "mnn", "mnn_centered", "ruvr", "ruvg", "gmm")
valid_classifiers <- c("logistic", "elnet", "elasticnet", "svm", "rf", "nnet", "knn", "xgboost", "rvc")
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

# Extract validated parameters
adjuster <- args$adjuster
classifier <- args$classifier
num_datasets <- args$num_datasets
test_study <- args$test_study
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
job_id <- sprintf("adjuster_%s_%s_%d_%s", adjuster, classifier, num_datasets, test_study)

# Main job wrapper with comprehensive error handling
main_job_wrapper <- function() {
  tryCatch({
    # Print job parameters for logging
    cat("=== ADJUSTER COMPARISON JOB ===\n")
    cat(sprintf("Job ID: %s\n", job_id))
    cat(sprintf("Adjuster: %s\n", adjuster))
    cat(sprintf("Classifier: %s\n", classifier))
    cat(sprintf("Num datasets: %d\n", num_datasets))
    cat(sprintf("Test study: %s\n", test_study))
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
    cat(sprintf("[ERROR] Parameters: adjuster=%s, classifier=%s, num_datasets=%d, test_study=%s\n", 
                adjuster, classifier, num_datasets, test_study), file = stderr())
    
    # Additional debugging information
    cat(sprintf("[ERROR] Working directory: %s\n", getwd()), file = stderr())
    cat(sprintf("[ERROR] Output file: %s\n", output_file), file = stderr())
    
    # Check input files
    data_path <- "data/TB_real_data.RData"
    helper_path <- "scripts/helper.R"
    
    cat(sprintf("[ERROR] Data file exists: %s\n", file.exists(data_path)), file = stderr())
    cat(sprintf("[ERROR] Helper file exists: %s\n", file.exists(helper_path)), file = stderr())
    
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
  data_path <- "data/TB_real_data.RData"
  if (!file.exists(data_path)) {
    stop(sprintf("Data file not found: %s", data_path))
  }
  
  load(data_path)
  source("scripts/helper.R")
  
  if (classifier == "rvc"){
    import_reticulate()
    if (is.null(RVC_py) || is.null(np_py)) {
      stop("Classifier 'rvc' was requested, but Python dependencies 'sklearn_rvm' or 'numpy' could not be imported. Please install them.")
    }
  }
  
  # ====================================================================
  # REAL DATA PREPARATION LOGIC
  # ====================================================================
  
  filter_studies <- function(dat_lst, label_lst, n_studies) {
    all_studies <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
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
  
  # Validate test study is in the filtered list
  if (!test_study %in% study_names) {
    stop(sprintf("Test study '%s' not found in selected studies: %s", 
                 test_study, paste(study_names, collapse=", ")))
  }
  
  test_name <- test_study
  cat(sprintf("Using test study: %s\n", test_name))
  
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
  # BATCH CORRECTION 
  # ====================================================================
  
  apply_batch_corrections <- function(dat, batch, group, dat_test, method) {
    if (method == "unadjusted") {
      # No correction - return original data
      return(list(
        dat_corrected = dat,
        dat_test_corrected = dat_test
      ))
    } else if (method == "combat") {
      # ComBat correction without labels (unsupervised)
      # Step 1: Correct batch effects within training data without using labels
      dat_corrected <- ComBat(dat, batch=batch, mod=NULL)
      
      # Step 2: Adjust test data to match corrected training distribution
      # Use entire corrected training set as reference batch
      combined_dat <- cbind(dat_corrected, dat_test)
      ref_batch_id <- 1  # Training data batch ID
      test_batch_id <- 2  # Test data batch ID
      combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), 
                         rep(test_batch_id, ncol(dat_test)))
      
      # Apply ComBat with training as reference (no mod matrix)
      combat_combined <- ComBat(combined_dat, batch=combined_batch, 
                               mod=NULL, ref.batch=ref_batch_id)
      
      # Extract corrected test data (training data unchanged as it's the reference)
      dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
      
      return(list(
        dat_corrected = dat_corrected,
        dat_test_corrected = dat_test_corrected
      ))
    } else if (method == "combat_mean") {
      # ComBat correction with mean adjustment only (no variance adjustment)
      # Step 1: Correct batch effects within training data (mean only)
      dat_corrected <- ComBat(dat, batch=batch, mod=NULL, mean.only=TRUE)
      
      # Step 2: Adjust test data to match corrected training distribution
      combined_dat <- cbind(dat_corrected, dat_test)
      ref_batch_id <- 1
      test_batch_id <- 2
      combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), 
                         rep(test_batch_id, ncol(dat_test)))
      
      # Apply ComBat with training as reference (mean only)
      combat_combined <- ComBat(combined_dat, batch=combined_batch, 
                               mod=NULL, ref.batch=ref_batch_id, mean.only=TRUE)
      
      dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
      
      return(list(
        dat_corrected = dat_corrected,
        dat_test_corrected = dat_test_corrected
      ))
    } else if (method == "combat_sup") {
      # ComBat correction with labels (supervised)
      # Step 1: Correct batch effects within training data while preserving biological signal
      dat_corrected <- ComBat(dat, batch=batch, mod=model.matrix(~group))
      
      # Step 2: Adjust test data to match corrected training distribution
      # Use entire corrected training set as reference batch
      combined_dat <- cbind(dat_corrected, dat_test)
      ref_batch_id <- 1  # Training data batch ID
      test_batch_id <- 2  # Test data batch ID
      combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), 
                         rep(test_batch_id, ncol(dat_test)))
      
      # Apply ComBat with training as reference (no mod matrix to avoid using test labels)
      combat_combined <- ComBat(combined_dat, batch=combined_batch, 
                               mod=NULL, ref.batch=ref_batch_id)
      
      # Extract corrected test data (training data unchanged as it's the reference)
      dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
      
      return(list(
        dat_corrected = dat_corrected,
        dat_test_corrected = dat_test_corrected
      ))
    } else if (method == "mnn") {
      library(batchelor, quietly = TRUE)
      library(SummarizedExperiment, quietly = TRUE)

      # MNN without pre-centering
      combined_dat <- cbind(dat, dat_test)
      
      # Create batch vector (Test set gets a new unique ID)
      test_id <- max(batch) + 1
      combined_batch <- c(batch, rep(test_id, ncol(dat_test)))
      
      # Determine merge order: Training batches by size (descending), Test set last.
      u_batches <- sort(unique(batch))
      train_sizes <- table(batch)[as.character(u_batches)]
      train_ord <- order(train_sizes, decreasing = TRUE)
      # merge.order needs actual batch IDs, not indices
      merge_ord <- c(u_batches[train_ord], test_id)
      
      cat(sprintf("DEBUG MNN: batch range [%s], test_id=%s\n", paste(range(batch), collapse="-"), test_id))
      cat(sprintf("DEBUG MNN: unique batches: %s\n", paste(u_batches, collapse=",")))
      cat(sprintf("DEBUG MNN: merge order: %s\n", paste(merge_ord, collapse=",")))
      cat(sprintf("DEBUG MNN: combined_batch length=%d, unique values: %s\n", 
                  length(combined_batch), paste(unique(combined_batch), collapse=",")))
      cat(sprintf("DEBUG MNN: combined_dat dims: %d x %d\n", nrow(combined_dat), ncol(combined_dat)))
      
      tryCatch({
        mnn_object <- batchelor::mnnCorrect(combined_dat, batch = combined_batch, merge.order = merge_ord)
        mnn_matrix <- SummarizedExperiment::assay(mnn_object)
      }, error = function(e) {
        cat(sprintf("[ERROR] mnnCorrect failed: %s\n", e$message), file = stderr())
        cat(sprintf("[ERROR] Error class: %s\n", paste(class(e), collapse=", ")), file = stderr())
        cat(sprintf("[ERROR] Batch info: unique=%s, length=%d\n", 
                    paste(unique(combined_batch), collapse=","), length(combined_batch)), file = stderr())
        cat(sprintf("[ERROR] Merge order: %s (class: %s)\n", 
                    paste(merge_ord, collapse=","), class(merge_ord)), file = stderr())
        stop(sprintf("MNN correction failed: %s", e$message))
      })
      
      # Split result
      dat_corrected <- mnn_matrix[, 1:ncol(dat)]
      dat_test_corrected <- mnn_matrix[, (ncol(dat) + 1):ncol(mnn_matrix)]
      
      return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))
    } else if (method == "mnn_centered") {
      library(batchelor, quietly = TRUE)
      library(SummarizedExperiment, quietly = TRUE)

      # Pre-center each gene within each batch before MNN
      cat("Pre-centering data for MNN...\n")
      
      # Center training data within each batch
      dat_centered <- dat
      for (b in unique(batch)) {
        batch_idx <- batch == b
        gene_means <- rowMeans(dat[, batch_idx, drop = FALSE])
        dat_centered[, batch_idx] <- dat[, batch_idx] - gene_means
      }
      
      # Center test data (as its own batch)
      test_gene_means <- rowMeans(dat_test)
      dat_test_centered <- dat_test - test_gene_means
      
      # Combine centered data
      combined_dat <- cbind(dat_centered, dat_test_centered)
      
      # Create batch vector (Test set gets a new unique ID)
      test_id <- max(batch) + 1
      combined_batch <- c(batch, rep(test_id, ncol(dat_test)))
      
      # Determine merge order: Training batches by size (descending), Test set last.
      u_batches <- sort(unique(batch))
      train_sizes <- table(batch)[as.character(u_batches)]
      train_ord <- order(train_sizes, decreasing = TRUE)
      merge_ord <- c(u_batches[train_ord], test_id)
      
      cat(sprintf("DEBUG MNN_CENTERED: batch range [%s], test_id=%s\n", paste(range(batch), collapse="-"), test_id))
      cat(sprintf("DEBUG MNN_CENTERED: merge order: %s\n", paste(merge_ord, collapse=",")))
      
      tryCatch({
        mnn_object <- batchelor::mnnCorrect(combined_dat, batch = combined_batch, merge.order = merge_ord)
        mnn_matrix <- SummarizedExperiment::assay(mnn_object)
      }, error = function(e) {
        cat(sprintf("[ERROR] mnnCorrect (centered) failed: %s\n", e$message), file = stderr())
        stop(sprintf("MNN centered correction failed: %s", e$message))
      })
      
      # Split result
      dat_corrected <- mnn_matrix[, 1:ncol(dat)]
      dat_test_corrected <- mnn_matrix[, (ncol(dat) + 1):ncol(mnn_matrix)]
      
      return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))
    } else if (method == "ruvr") {
      # RUVr: Remove Unwanted Variation using Residuals
      # Custom implementation without ruv package dependency
      cat("Applying RUVr correction...\n")
      
      # Step 1: Fit initial GLM on training data to get residuals
      # Create design matrix with TB status and batch
      design <- model.matrix(~ group + batch)
      
      # Fit gene-wise linear models
      cat("Fitting initial GLM to estimate residuals...\n")
      residuals <- matrix(NA, nrow = nrow(dat), ncol = ncol(dat))
      for (i in 1:nrow(dat)) {
        fit <- lm(dat[i, ] ~ group + batch)
        residuals[i, ] <- residuals(fit)
      }
      
      # Step 2: Estimate unwanted variation factors from residuals using SVD
      k <- 3  # Number of unwanted variation factors (can be tuned)
      cat(sprintf("Estimating %d unwanted variation factors...\n", k))
      
      svd_res <- svd(residuals)
      W <- svd_res$u[, 1:k, drop = FALSE]  # Factor loadings (genes x k)
      alpha <- svd_res$v[, 1:k, drop = FALSE] %*% diag(svd_res$d[1:k])  # Factor scores (samples x k)
      
      # Step 3: Correct training data by regressing out the factors
      dat_corrected <- dat
      for (i in 1:nrow(dat)) {
        fit <- lm(dat[i, ] ~ alpha)
        dat_corrected[i, ] <- residuals(fit) + mean(dat[i, ])
      }
      
      # Step 4: Project test data onto the learned factors and correct
      cat("Projecting test data onto learned factors...\n")
      # Project test data onto factor space: alpha_test = dat_test^T %*% W
      alpha_test <- t(dat_test) %*% W
      
      dat_test_corrected <- dat_test
      for (i in 1:nrow(dat_test)) {
        fit <- lm(dat_test[i, ] ~ alpha_test)
        dat_test_corrected[i, ] <- residuals(fit) + mean(dat_test[i, ])
      }
      
      cat("RUVr correction complete\n")
      
      return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))
    } else if (method == "ruvg") {
      # RUVg: Remove Unwanted Variation using housekeeping genes
      # Uses negative control genes to estimate unwanted variation
      cat("Applying RUVg correction...\n")
      
      # Define housekeeping genes
      housekeeping_genes <- c("GAPDH", "ACTG1", "RPS18", "POM121C", "MRPL18", 
                             "TOMM5", "YTHDF1", "TPT1", "RPS27")
      
      # Find which housekeeping genes are present in the data
      available_hk <- intersect(housekeeping_genes, rownames(dat))
      
      if (length(available_hk) == 0) {
        stop("None of the housekeeping genes found in data. Cannot apply RUVg.")
      }
      
      cat(sprintf("Using %d housekeeping genes: %s\n", 
                  length(available_hk), paste(available_hk, collapse=", ")))
      
      # Extract housekeeping gene expression
      hk_dat <- dat[available_hk, , drop = FALSE]
      
      # Step 1: Estimate unwanted variation factors from housekeeping genes using SVD
      k <- min(3, length(available_hk) - 1)  # Number of factors (limited by number of HK genes)
      cat(sprintf("Estimating %d unwanted variation factors from housekeeping genes...\n", k))
      
      # Center housekeeping gene data
      hk_centered <- hk_dat - rowMeans(hk_dat)
      
      svd_res <- svd(hk_centered)
      W <- svd_res$u[, 1:k, drop = FALSE]  # Factor loadings (HK genes x k)
      alpha <- svd_res$v[, 1:k, drop = FALSE] %*% diag(svd_res$d[1:k])  # Factor scores (samples x k)
      
      # Step 2: Correct training data by regressing out the factors
      dat_corrected <- dat
      for (i in 1:nrow(dat)) {
        fit <- lm(dat[i, ] ~ alpha)
        dat_corrected[i, ] <- residuals(fit) + mean(dat[i, ])
      }
      
      # Step 3: Estimate factors for test data using housekeeping genes
      cat("Estimating unwanted variation in test data...\n")
      hk_test <- dat_test[available_hk, , drop = FALSE]
      hk_test_centered <- hk_test - rowMeans(hk_dat)  # Use training means for centering
      
      # Project test housekeeping genes onto factor space
      alpha_test <- t(hk_test_centered) %*% W
      
      # Step 4: Correct test data
      dat_test_corrected <- dat_test
      for (i in 1:nrow(dat_test)) {
        fit <- lm(dat_test[i, ] ~ alpha_test)
        dat_test_corrected[i, ] <- residuals(fit) + mean(dat_test[i, ])
      }
      
      cat("RUVg correction complete\n")
      
      return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))
    } else if (method == "gmm") {
      # GMM adjustment: fits 2-component GMM to each gene within each batch
      cat("Applying GMM adjustment...\n")
      
      # Source the GMM adjustment function (use absolute path from workspace root)
      gmm_script <- file.path(getwd(), "..", "..", "adjust", "gmm_adjust.R")
      if (!file.exists(gmm_script)) {
        gmm_script <- "scripts/adjust/gmm_adjust.R"  # Fallback to relative from workspace root
      }
      source(gmm_script)
      
      # Apply GMM to training data (genes are rows, samples are columns)
      dat_corrected <- gmm_adjust(
        data = dat,
        batch = batch,
        genes_are_columns = FALSE,  # Our data has genes as rows
        mean_mean_zero = TRUE,
        unit_var = TRUE,
        log_transform = FALSE,  # Data is already log-transformed
        debug = FALSE,
        num_workers = 1
      )
      
      # Apply GMM to test data (treat as single batch)
      dat_test_corrected <- gmm_adjust(
        data = dat_test,
        batch = rep(1, ncol(dat_test)),  # Single batch
        genes_are_columns = FALSE,
        mean_mean_zero = TRUE,
        unit_var = TRUE,
        log_transform = FALSE,
        debug = FALSE,
        num_workers = 1
      )
      
      cat("GMM adjustment complete\n")
      
      return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))
    } else {
      stop(sprintf("Unknown batch correction method: %s", method))
    }
  }
  
  #' Global scaling: scale entire dataset to have overall variance = 1
  #' Preserves relative gene importance while putting data on consistent scale
  #' @param dat_train Training data matrix
  #' @param dat_test Test data matrix
  #' @return List with scaled training and test data
  global_scale <- function(dat_train, dat_test) {
    # Compute global mean and SD from training data
    train_mean <- mean(dat_train)
    train_sd <- sd(as.vector(dat_train))
    
    # Apply same transformation to both train and test
    dat_train_scaled <- (dat_train - train_mean) / train_sd
    dat_test_scaled <- (dat_test - train_mean) / train_sd
    
    cat(sprintf("Global scaling: mean=%.4f, sd=%.4f\n", train_mean, train_sd))
    
    list(
      dat_train = dat_train_scaled,
      dat_test = dat_test_scaled
    )
  }
  
  # ====================================================================
  # EXECUTE BATCH CORRECTION
  # ====================================================================
  
  cat(sprintf("Applying batch correction method: %s\n", adjuster))
  
  # Apply batch correction
  batch_corr_result <- apply_batch_corrections(dat, datasets$batch, datasets$group, dat_test, adjuster)
  dat_corrected <- batch_corr_result$dat_corrected
  dat_test_corrected <- batch_corr_result$dat_test_corrected
  
  # Global scaling (not per-gene normalization)
  cat("Applying global scaling to training and test data...\n")
  
  if(is.null(dat_test_corrected)) {
    stop("Test data is NULL after batch correction")
  }
  
  scaled_data <- global_scale(dat_corrected, dat_test_corrected)
  dat_train_norm <- scaled_data$dat_train
  dat_test_norm <- scaled_data$dat_test
  
  # Validation
  if (any(is.na(dat_train_norm)) || any(is.na(dat_test_norm))) {
    stop("Scaling produced NA values")
  }
  if (is.null(dat_test_norm)) {
    stop("Test data became NULL during scaling")
  }
  
  cat(sprintf("Scaled data dimensions - Train: %d x %d, Test: %d x %d\n", 
              nrow(dat_train_norm), ncol(dat_train_norm),
              nrow(dat_test_norm), ncol(dat_test_norm)))
  
  cat(sprintf("Batch correction completed successfully\n"))
  cat(sprintf("  Method: %s\n", adjuster))
  cat(sprintf("  Training data shape: %d x %d\n", nrow(dat_train_norm), ncol(dat_train_norm)))
  cat(sprintf("  Test data shape: %d x %d\n", nrow(dat_test_norm), ncol(dat_test_norm)))
  
  # ====================================================================
  # SINGLE CLASSIFIER TRAINING AND EVALUATION
  # ====================================================================
  
  #' Train single classifier and evaluate performance
  # RVC needs it's own branch since it's from python
  train_and_evaluate_classifier <- function(classifier_type, train_data, train_labels, test_data, test_labels) {
    
    # Initialize variables
    trained_model <- NULL
    test_predictions <- NULL
    
    if (classifier_type == "rvc") {
      cat("Training rvc classifier (using reticulate)...\n")
      
      # Transpose data: R (features x samples) -> Python (samples x features)
      X_train_t <- t(train_data)
      X_test_t <- t(test_data)
      
      X_train_py <- r_to_py(X_train_t)
      y_train_r <- as.numeric(as.factor(train_labels)) - 1
      y_train_py <- r_to_py(y_train_r)
      X_test_py <- r_to_py(X_test_t)
      
      # Using rbf kernel
      model_py <- RVC_py(kernel = "rbf")
      model_py$fit(X_train_py, y_train_py)
      
      # This returns (n_samples, n_classes)
      preds_py <- model_py$predict_proba(X_test_py)
      
      # Convert predictions back to R
      # We need the probability of the positive class (class "1"), which is the 2nd column
      preds_r <- py_to_r(preds_py)
      test_predictions <- preds_r[, 2] 
      
      trained_model <- list(mod = model_py)
      
    } else {
      # Get prediction function for classifier type
      learner_fit <- getPredFunctions(classifier_type)
      
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
  create_output_dataframe <- function(adjuster, classifier, n_datasets, test_study, performance_metrics) {
    # Create one row per metric
    output_rows <- lapply(names(performance_metrics), function(metric) {
      data.frame(
        adjuster = adjuster,
        classifier = classifier,
        n_datasets = n_datasets,
        test_study = test_study,
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
    test_study = test_study,
    performance_metrics = result$performance
  )
  
  # Write to CSV file with error handling
  tryCatch({
    write.csv(output_df, file = output_file, row.names = FALSE)
    
    # Verify file was created
    if (!file.exists(output_file)) {
      stop(sprintf("File was not created: %s", output_file))
    }
    
    # Verify file has content
    file_size <- file.info(output_file)$size
    if (is.na(file_size) || file_size == 0) {
      stop(sprintf("File was created but is empty: %s", output_file))
    }
    
    cat(sprintf("Results written to: %s\n", output_file))
    cat(sprintf("Output contains %d rows (one per metric)\n", nrow(output_df)))
    cat(sprintf("File size: %d bytes\n", file_size))
    
  }, error = function(e) {
    cat(sprintf("[ERROR] Failed to write output file: %s\n", e$message), file = stderr())
    cat(sprintf("[ERROR] Output file path: %s\n", output_file), file = stderr())
    cat(sprintf("[ERROR] Output directory exists: %s\n", dir.exists(dirname(output_file))), file = stderr())
    cat(sprintf("[ERROR] Output directory writable: %s\n", file.access(dirname(output_file), 2) == 0), file = stderr())
    stop(sprintf("Failed to write output file: %s", e$message))
  })
  
  # Display output for verification
  cat("\nOutput preview:\n")
  print(output_df)
  
  cat(sprintf("\n=== JOB COMPLETED SUCCESSFULLY ===\n"))
  cat(sprintf("Adjuster: %s\n", adjuster))
  cat(sprintf("Classifier: %s\n", classifier))
  cat(sprintf("Datasets: %d\n", num_datasets))
  cat(sprintf("Test study: %s\n", test_study))
  cat(sprintf("Output: %s\n", output_file))
  cat("==================================\n")
  
  return(output_df)
}

# ====================================================================
# EXECUTE MAIN JOB
# ====================================================================

# Store result (or NULL if it crashed inside wrapper, though wrapper handles that)
res <- main_job_wrapper()

# Force a clean exit to signal to Snakemake that we are happy
quit(save = "no", status = 0)
