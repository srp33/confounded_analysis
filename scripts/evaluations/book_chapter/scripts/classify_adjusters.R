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
  required_packages <- c("glmnet", "SummarizedExperiment", "sva", "DESeq2", 
                        "ROCR", "ggplot2", "gridExtra", "reshape2", 
                        "dplyr", "purrr", "nnls", "batchelor",
                        "argparse", "class", "xgboost", "sda", "klaR")
  sapply(required_packages, require, character.only=TRUE, quietly=TRUE)
}))

# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Execute single adjuster comparison job for batch correction analysis")

parser$add_argument("--adjuster", type = "character", required = TRUE,
                   help = "Batch correction method: unadjusted, naive, rank_samples, rank_twice, npn, combat, combat_mean, combat_sup, mnn, fast_mnn, ruvr, ruvg, gmm, or posse")
parser$add_argument("--classifier", type = "character", required = TRUE,
                   help = "Classifier type: rda, elnet, elasticnet, svm, rf, nnet, knn, xgboost, or shrinkageLDA")
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
valid_adjusters <- c("unadjusted", "naive", "rank_samples", "rank_twice", "npn", "combat", "combat_mean", "combat_sup", "mnn", "fast_mnn", "ruvr", "ruvg", "gmm", "posse_default", "posse_aggressive", "posse_focused", "posse_conservative", "posse_housekeeping", "posse_ultra_aggressive", "posse_two_phase", "posse_sniper")
valid_classifiers <- c("rda", "elnet", "elasticnet", "svm", "rf", "nnet", "knn", "xgboost", "shrinkageLDA")
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
# RANK ADJUSTMENT HELPER FUNCTIONS
# ====================================================================

rank_normalized <- function(matrix_, dim) {
  if (dim < 1 || dim > 2) {
    stop("Invalid dimension. Must be 1 for rows or 2 for columns.")
  }
  ranked = apply(matrix_, dim, rank, ties.method = "average")
  
  # apply() transposes the result when dim=1, so we need to transpose it back
  # When dim=1: apply ranks across columns (samples) for each row (feature)
  # When dim=2: apply ranks across rows (features) for each column (sample)
  if (dim == 1 && is.matrix(ranked)) {
    ranked = t(ranked)
  }
  
  return(ranked / max(ranked, na.rm = TRUE))
}

adjust_ranked_with_batch_info <- function(matrix_, batch, debug = FALSE) {
  #' Normalize sample-wise by ranking the genes within the sample, and then by batch.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param debug Logical flag for debug output.
  #' @return The adjusted matrix (features x samples).
  
  cat("Adjusting with ranked with batch info.\n")
  ranked = rank_normalized(matrix_, 1)
  
  if (debug) {
    cat("DEBUG: matrix_ dimensions: ", nrow(matrix_), " x ", ncol(matrix_), "\n")
    cat("DEBUG: ranked dimensions: ", nrow(ranked), " x ", ncol(ranked), "\n")
  }
  
  batch_levels <- unique(batch)
  ranked2 <- matrix(NA, nrow = nrow(ranked), ncol = ncol(ranked))
  
  for (b in batch_levels) {
    # For each batch, we rank by sample.
    batch_indices <- which(batch == b)
    batch_data <- ranked[, batch_indices, drop = FALSE]
    
    if (debug) {
      cat("DEBUG: Processing batch '", b, "' with ", length(batch_indices), " samples\n")
      cat("DEBUG: batch_data dimensions: ", nrow(batch_data), " x ", ncol(batch_data), "\n")
    }
    
    # Only apply ranking if there's more than one sample in the batch
    if (ncol(batch_data) > 1) {
      batch_ranked <- rank_normalized(batch_data, 2)
      if (debug) {
        cat("DEBUG: batch_ranked dimensions: ", nrow(batch_ranked), " x ", ncol(batch_ranked), "\n")
      }
      ranked2[, batch_indices] <- batch_ranked
    } else {
      # For single-sample batches, just use the original ranked values
      ranked2[, batch_indices] <- batch_data
    }
  }
  
  # Handle any remaining NA values
  if (any(is.na(ranked2))) {
    cat("WARNING: Found NA values in ranked2 matrix. Replacing with original ranked values.\n")
    ranked2[is.na(ranked2)] <- ranked[is.na(ranked2)]
  }
  
  max_val <- max(ranked2, na.rm = TRUE)
  if (max_val == 0) {
    cat("WARNING: Maximum value in ranked2 is 0. Using 1 as denominator.\n")
    max_val <- 1
  }
  
  return(ranked2 / max_val)
}

adjust_npn <- function(matrix_, batch, debug = FALSE) {
  #' Adjust matrix using Nonparanormal (NPN) transformation.
  #' If batch is NULL, the entire matrix is adjusted at once.
  #' Assumes the batch vector contains no NA values.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @return The adjusted matrix.
  
  if (is.null(batch)) {
    # If no batch is provided, adjust the whole matrix.
    cat("Batch is NULL. Adjusting entire matrix with NPN transformation.\n")
    
    # Transpose to (samples x features) for huge.npn.
    matrix_t <- t(matrix_)
    
    npn_transformed_t <- huge::huge.npn(matrix_t, verbose = FALSE)
    
    return(t(npn_transformed_t))
    
  } else {
    cat("Adjusting using Nonparanormal (NPN) transformation by batch.\n")
    
    # Split the matrix by batch.
    batch_levels <- unique(batch)
    matrix_by_batch <- list()
    
    for (b in batch_levels) {
      batch_indices <- which(batch == b)
      if (length(batch_indices) > 0) {
        matrix_by_batch[[as.character(b)]] <- matrix_[, batch_indices, drop = FALSE]
      }
    }
    
    # Apply NPN transformation to each batch.
    for (b in names(matrix_by_batch)) {
      matrix_t <- t(matrix_by_batch[[b]])
      npn_transformed_t <- huge::huge.npn(matrix_t, verbose = FALSE)
      matrix_by_batch[[b]] <- t(npn_transformed_t)
    }
    
    # Reassemble the matrix from the adjusted batches.
    result_matrix <- matrix_
    for (b in names(matrix_by_batch)) {
      batch_indices <- which(batch == as.character(b))
      result_matrix[, batch_indices] <- matrix_by_batch[[b]]
    }
    
    return(result_matrix)
  }
}

adjust_ranked_samples_with_batch_info <- function(matrix_, batch, debug = FALSE) {
  #' Rank samples within each gene (across samples), then merge batches with normalized ranks.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param debug Logical flag for debug output.
  #' @return The adjusted matrix (features x samples).
  
  cat("Adjusting with rank_samples (rank samples within genes, batch-aware).\n")
  
  batch_levels <- unique(batch)
  result_matrix <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  
  for (b in batch_levels) {
    # For each batch, rank samples within each gene
    batch_indices <- which(batch == b)
    batch_data <- matrix_[, batch_indices, drop = FALSE]
    
    if (debug) {
      cat("DEBUG: Processing batch '", b, "' with ", length(batch_indices), " samples\n")
      cat("DEBUG: batch_data dimensions: ", nrow(batch_data), " x ", ncol(batch_data), "\n")
    }
    
    # Only apply ranking if there's more than one sample in the batch
    if (ncol(batch_data) > 1) {
      # Rank samples within each gene (dim=2 means rank across columns/samples for each row/gene)
      batch_ranked <- rank_normalized(batch_data, 2)
      if (debug) {
        cat("DEBUG: batch_ranked dimensions: ", nrow(batch_ranked), " x ", ncol(batch_ranked), "\n")
      }
      result_matrix[, batch_indices] <- batch_ranked
    } else {
      # For single-sample batches, just normalize to [0,1]
      result_matrix[, batch_indices] <- batch_data / max(batch_data, na.rm = TRUE)
    }
  }
  
  # Handle any remaining NA values
  if (any(is.na(result_matrix))) {
    cat("WARNING: Found NA values in result_matrix. Replacing with original values.\n")
    result_matrix[is.na(result_matrix)] <- matrix_[is.na(result_matrix)]
  }
  
  max_val <- max(result_matrix, na.rm = TRUE)
  if (max_val == 0) {
    cat("WARNING: Maximum value in result_matrix is 0. Using 1 as denominator.\n")
    max_val <- 1
  }
  
  return(result_matrix / max_val)
}

adjust_ranked_twice_with_batch_info <- function(matrix_, batch, debug = FALSE) {
  #' Double ranking: first rank genes within samples, then rank samples within genes, batch-aware.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param debug Logical flag for debug output.
  #' @return The adjusted matrix (features x samples).
  
  cat("Adjusting with rank_twice (genes within samples, then samples within genes, batch-aware).\n")
  
  batch_levels <- unique(batch)
  result_matrix <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  
  for (b in batch_levels) {
    # For each batch, apply double ranking
    batch_indices <- which(batch == b)
    batch_data <- matrix_[, batch_indices, drop = FALSE]
    
    if (debug) {
      cat("DEBUG: Processing batch '", b, "' with ", length(batch_indices), " samples\n")
      cat("DEBUG: batch_data dimensions: ", nrow(batch_data), " x ", ncol(batch_data), "\n")
    }
    
    # Only apply ranking if there's more than one sample in the batch
    if (ncol(batch_data) > 1) {
      # First: rank genes within each sample (dim=1)
      # Then: rank samples within each gene (dim=2)
      batch_ranked <- rank_normalized(rank_normalized(batch_data, 1), 2)
      if (debug) {
        cat("DEBUG: batch_ranked dimensions: ", nrow(batch_ranked), " x ", ncol(batch_ranked), "\n")
      }
      result_matrix[, batch_indices] <- batch_ranked
    } else {
      # For single-sample batches, just apply single ranking (genes within sample)
      batch_ranked <- rank_normalized(batch_data, 1)
      result_matrix[, batch_indices] <- batch_ranked
    }
  }
  
  # Handle any remaining NA values
  if (any(is.na(result_matrix))) {
    cat("WARNING: Found NA values in result_matrix. Replacing with original values.\n")
    result_matrix[is.na(result_matrix)] <- matrix_[is.na(result_matrix)]
  }
  
  max_val <- max(result_matrix, na.rm = TRUE)
  if (max_val == 0) {
    cat("WARNING: Maximum value in result_matrix is 0. Using 1 as denominator.\n")
    max_val <- 1
  }
  
  return(result_matrix / max_val)
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
  
  # DIAGNOSTIC: Print loaded data characteristics
  cat("[DATA DIAGNOSTIC] Loaded TB data characteristics:\n")
  if (exists("dat_lst")) {
    cat(sprintf("  dat_lst contains %d studies: %s\n", length(dat_lst), paste(names(dat_lst), collapse=", ")))
    
    # Check first study's data characteristics
    first_study <- names(dat_lst)[1]
    first_data <- dat_lst[[first_study]]
    cat(sprintf("  First study (%s) data shape: %d x %d\n", first_study, nrow(first_data), ncol(first_data)))
    cat(sprintf("  First study data range: [%.3f, %.3f]\n", min(first_data), max(first_data)))
    cat(sprintf("  First study data mean: %.3f, sd: %.3f\n", mean(first_data), sd(first_data)))
    cat(sprintf("  First study has negative values: %s\n", any(first_data < 0)))
    cat(sprintf("  First study appears to be integers: %s\n", all(first_data == round(first_data))))
    
    # Sample values from first study
    sample_values <- first_data[1:min(3, nrow(first_data)), 1:min(3, ncol(first_data))]
    cat("  Sample values from first study:\n")
    print(sample_values)
  } else {
    cat("  ERROR: dat_lst not found in loaded data\n")
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
    
    # Ensure all datasets have the same genes by taking intersection
    all_datasets <- c(train_name, test_name)
    common_genes <- Reduce(intersect, lapply(dat_lst[all_datasets], rownames))
    
    cat(sprintf("  Gene intersection: %d common genes across all datasets\n", length(common_genes)))
    
    # Subset all datasets to common genes
    dat_lst_subset <- lapply(dat_lst[all_datasets], function(x) x[common_genes, , drop=FALSE])
    
    # Combine training datasets
    dat <- do.call(cbind, dat_lst_subset[train_name])
    batch <- rep(1:length(train_name), times=sapply(dat_lst_subset[train_name], ncol))
    batches_ind <- lapply(1:length(train_name), function(x) which(batch == x))
    batch_names <- levels(factor(batch))
    group <- do.call(c, label_lst[train_name])
    y_sgbatch_train <- lapply(batch_names, function(x) group[batch == x])
    
    dat_test <- dat_lst_subset[[test_name]]
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
  
  # DIAGNOSTIC: Check data after preparation
  cat(sprintf("[POST-PREPARE DIAGNOSTIC] Data after prepare_datasets:\n"))
  cat(sprintf("  dat shape: %d x %d\n", nrow(datasets$dat), ncol(datasets$dat)))
  cat(sprintf("  dat range: [%.3f, %.3f]\n", min(datasets$dat), max(datasets$dat)))
  cat(sprintf("  dat mean: %.3f, sd: %.3f\n", mean(datasets$dat), sd(datasets$dat)))
  
  # Validate datasets
  if(is.null(datasets$dat_test)) {
    stop(sprintf("Test dataset '%s' is NULL or missing from data", test_name))
  }
  if(ncol(datasets$dat_test) == 0) {
    stop(sprintf("Test dataset '%s' has no samples", test_name))
  }
  
  # Feature reduction (top 1000 most variable genes)
  # EXCEPTION: POSSE methods need access to all genes for pathway analysis
  if (grepl("^posse", adjuster, ignore.case = TRUE)) {
    cat(sprintf("[POSSE EXCEPTION] Skipping feature reduction - POSSE needs all genes for pathway analysis\n"))
    dat <- datasets$dat
    dat_test <- datasets$dat_test
    cat(sprintf("  POSSE data shape: %d genes x %d samples (training)\n", nrow(dat), ncol(dat)))
    cat(sprintf("  POSSE test shape: %d genes x %d samples (test)\n", nrow(dat_test), ncol(dat_test)))
  } else {
    n_highvar_genes <- 1000
    feat_reduced <- reduce_features(datasets$dat, datasets$dat_test, n_highvar_genes)
    dat <- feat_reduced$dat
    dat_test <- feat_reduced$dat_test
    
    # DIAGNOSTIC: Check data after feature reduction
    cat(sprintf("[POST-FEATURE DIAGNOSTIC] Data after feature reduction:\n"))
    cat(sprintf("  dat shape: %d x %d\n", nrow(dat), ncol(dat)))
    cat(sprintf("  dat range: [%.3f, %.3f]\n", min(dat), max(dat)))
    cat(sprintf("  dat mean: %.3f, sd: %.3f\n", mean(dat), sd(dat)))
  }
  
  # ====================================================================
  # DATA PREPROCESSING: LOG TRANSFORMATION FOR RAW INTENSITY DATA
  # EXCEPTION: POSSE methods handle their own transformations
  # ====================================================================
  
  if (grepl("^posse", adjuster, ignore.case = TRUE)) {
    cat(sprintf("[POSSE EXCEPTION] Skipping log transformation - POSSE handles its own arcsinh transformation\n"))
    cat(sprintf("  Raw data preserved for POSSE: range=[%.3f, %.3f], mean=%.3f\n", 
                min(dat), max(dat), mean(dat)))
  } else if (adjuster == "gmm") {
    cat(sprintf("[GMM EXCEPTION] Skipping log transformation - GMM handles its own log transformation\n"))
    cat(sprintf("  Raw data preserved for GMM: range=[%.3f, %.3f], mean=%.3f\n", 
                min(dat), max(dat), mean(dat)))
  } else {
    # Detect if training data appears to be raw intensity values that need log transformation
    # Heuristics: 
    # 1. Max value > 100 (log-transformed data typically < 20)
    # 2. Mean > 20 (log-transformed data typically < 10)
    # 3. Large dynamic range (max/median > 50)
    
    needs_log_transform_train <- max(dat) > 100 || mean(dat) > 20 || (max(dat) / median(dat)) > 50
    needs_log_transform_test <- max(dat_test) > 100 || mean(dat_test) > 20 || (max(dat_test) / median(dat_test)) > 50
    
    cat(sprintf("[PREPROCESSING] Log transformation assessment:\n"))
    cat(sprintf("  Training data needs log transform: %s\n", needs_log_transform_train))
    cat(sprintf("  Test data needs log transform: %s\n", needs_log_transform_test))
    
    if (needs_log_transform_train) {
      cat(sprintf("[PREPROCESSING] Applying log transformation to TRAINING data\n"))
      cat(sprintf("  Before: range=[%.3f, %.3f], mean=%.3f, median=%.3f\n", 
                  min(dat), max(dat), mean(dat), median(dat)))
      
      # Handle negative values in training data
      train_min <- min(dat)
      if (train_min < 0) {
        dat <- dat - train_min
        cat(sprintf("  After subtracting min (%.3f): range=[%.3f, %.3f]\n", 
                    train_min, min(dat), max(dat)))
      }
      
      # Add 1 before log transform
      dat <- dat + 1
      cat(sprintf("  After adding 1: range=[%.3f, %.3f]\n", min(dat), max(dat)))
      
      # Log2 transform
      dat <- log2(dat)
      cat(sprintf("  After log2: range=[%.3f, %.3f], mean=%.3f\n", 
                  min(dat), max(dat), mean(dat)))
      
      # Check for problems
      n_na <- sum(is.na(dat))
      n_inf <- sum(is.infinite(dat))
      if (n_na > 0 || n_inf > 0) {
        stop(sprintf("Log transformation of training data produced %d NA and %d Inf values", n_na, n_inf))
      }
    }
    
    if (needs_log_transform_test) {
      cat(sprintf("[PREPROCESSING] Applying log transformation to TEST data\n"))
      cat(sprintf("  Before: range=[%.3f, %.3f], mean=%.3f, median=%.3f\n", 
                  min(dat_test), max(dat_test), mean(dat_test), median(dat_test)))
      
      # Handle negative values in test data
      test_min <- min(dat_test)
      if (test_min < 0) {
        dat_test <- dat_test - test_min
        cat(sprintf("  After subtracting min (%.3f): range=[%.3f, %.3f]\n", 
                    test_min, min(dat_test), max(dat_test)))
      }
      
      # Add 1 before log transform
      dat_test <- dat_test + 1
      cat(sprintf("  After adding 1: range=[%.3f, %.3f]\n", min(dat_test), max(dat_test)))
      
      # Log2 transform
      dat_test <- log2(dat_test)
      cat(sprintf("  After log2: range=[%.3f, %.3f], mean=%.3f\n", 
                  min(dat_test), max(dat_test), mean(dat_test)))
      
      # Check for problems
      n_na <- sum(is.na(dat_test))
      n_inf <- sum(is.infinite(dat_test))
      if (n_na > 0 || n_inf > 0) {
        stop(sprintf("Log transformation of test data produced %d NA and %d Inf values", n_na, n_inf))
      }
    }
    
    if (!needs_log_transform_train && !needs_log_transform_test) {
      cat(sprintf("[PREPROCESSING] Both datasets appear already log-transformed - no transformation applied\n"))
    }
  }
  
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
  
  #' Save diagnostic parameters for batch correction methods
  #' @param method Correction method name
  #' @param raw_data Original data before correction
  #' @param corrected_data Data after correction
  #' @param gene_names Gene identifiers
  #' @param output_file Path to save diagnostic CSV
  save_correction_diagnostics <- function(method, raw_data, corrected_data, gene_names = NULL, output_file) {
    tryCatch({
      # DIAGNOSTIC: Print data characteristics
      cat(sprintf("[%s DIAGNOSTIC] Input data characteristics:\n", toupper(method)))
      cat(sprintf("  Raw data shape: %d genes x %d samples\n", nrow(raw_data), ncol(raw_data)))
      cat(sprintf("  Raw data range: [%.3f, %.3f]\n", min(raw_data), max(raw_data)))
      cat(sprintf("  Raw data mean: %.3f, sd: %.3f\n", mean(raw_data), sd(raw_data)))
      cat(sprintf("  Raw data has negative values: %s\n", any(raw_data < 0)))
      cat(sprintf("  Raw data appears to be integers: %s\n", all(raw_data == round(raw_data))))
      
      cat(sprintf("  Corrected data range: [%.3f, %.3f]\n", min(corrected_data), max(corrected_data)))
      cat(sprintf("  Corrected data mean: %.3f, sd: %.3f\n", mean(corrected_data), sd(corrected_data)))
      
      # Sample a few genes to show the transformation
      sample_genes <- min(3, nrow(raw_data))
      for (i in 1:sample_genes) {
        raw_samples <- raw_data[i, 1:min(3, ncol(raw_data))]
        corr_samples <- corrected_data[i, 1:min(3, ncol(corrected_data))]
        cat(sprintf("  Gene %d samples: raw=[%.3f, %.3f, %.3f], corr=[%.3f, %.3f, %.3f]\n", 
                   i, raw_samples[1], raw_samples[2], raw_samples[3],
                   corr_samples[1], corr_samples[2], corr_samples[3]))
      }
      
      # Create diagnostic directory
      diagnostic_dir <- dirname(output_file)
      if (!dir.exists(diagnostic_dir)) {
        dir.create(diagnostic_dir, recursive = TRUE)
      }
      
      # Calculate per-gene correction parameters
      n_genes <- nrow(raw_data)
      
      if (is.null(gene_names)) {
        gene_names <- paste0("gene_", 1:n_genes)
      }
      
      # Calculate per-gene alpha (scale) and beta (shift) parameters
      # These represent what correction was actually applied: corrected = alpha * raw + beta
      diagnostic_data <- data.frame(
        gene_id = gene_names,
        gene_index = 1:n_genes,
        stringsAsFactors = FALSE
      )
      
      # Calculate correction parameters for each gene
      for (i in 1:n_genes) {
        raw_gene <- raw_data[i, ]
        corrected_gene <- corrected_data[i, ]
        
        # Fit linear model: corrected = alpha * raw + beta
        if (var(raw_gene) > 1e-10) {  # Avoid division by zero
          alpha <- cov(corrected_gene, raw_gene) / var(raw_gene)
          beta <- mean(corrected_gene) - alpha * mean(raw_gene)
        } else {
          # If raw gene has no variance, use identity transformation
          alpha <- 1.0
          beta <- mean(corrected_gene) - mean(raw_gene)
        }
        
        diagnostic_data$alpha_final[i] <- alpha
        diagnostic_data$beta_final[i] <- beta
        
        # DIAGNOSTIC: Print first few genes' parameters
        if (i <= 3) {
          cat(sprintf("  Gene %d: alpha=%.6f, beta=%.3f, raw_mean=%.3f, raw_sd=%.3f\n", 
                     i, alpha, beta, mean(raw_gene), sd(raw_gene)))
        }
        
        # Calculate relative corrections compared to naive correction
        # Naive correction would set alpha=1 (no scale change) and beta=0 (no shift)
        # Measure deviation from naive as relative measures
        eps <- 1e-9
        
        # Alpha deviation: how much does the scale correction deviate from identity (1.0)?
        alpha_deviation_from_naive <- abs(alpha - 1.0)
        
        # Beta deviation: how much shift relative to the typical data scale?
        # Use the standard deviation of raw data as the scale reference
        data_scale <- sd(raw_gene) + eps
        beta_deviation_from_naive <- abs(beta) / data_scale
        
        # Combined relative correction magnitude
        relative_correction_magnitude <- alpha_deviation_from_naive + beta_deviation_from_naive
        
        diagnostic_data$alpha_deviation_from_naive[i] <- alpha_deviation_from_naive
        diagnostic_data$beta_deviation_from_naive[i] <- beta_deviation_from_naive
        diagnostic_data$relative_correction_magnitude[i] <- relative_correction_magnitude
        diagnostic_data$data_scale[i] <- data_scale
      }
      
      # DIAGNOSTIC: Print summary of calculated parameters
      cat(sprintf("[%s DIAGNOSTIC] Calculated parameters summary:\n", toupper(method)))
      cat(sprintf("  Alpha range: [%.6f, %.6f]\n", min(diagnostic_data$alpha_final), max(diagnostic_data$alpha_final)))
      cat(sprintf("  Beta range: [%.3f, %.3f]\n", min(diagnostic_data$beta_final), max(diagnostic_data$beta_final)))
      cat(sprintf("  Alpha mean: %.6f, sd: %.6f\n", mean(diagnostic_data$alpha_final), sd(diagnostic_data$alpha_final)))
      cat(sprintf("  Beta mean: %.3f, sd: %.3f\n", mean(diagnostic_data$beta_final), sd(diagnostic_data$beta_final)))
      cat(sprintf("  Data scale range: [%.3f, %.3f]\n", min(diagnostic_data$data_scale), max(diagnostic_data$data_scale)))
      
      # Add method-specific metadata
      diagnostic_data$method <- method
      diagnostic_data$alpha_final_mean <- mean(diagnostic_data$alpha_final, na.rm = TRUE)
      diagnostic_data$beta_final_mean <- mean(diagnostic_data$beta_final, na.rm = TRUE)
      diagnostic_data$correction_magnitude <- mean(abs(diagnostic_data$alpha_final - 1.0), na.rm = TRUE)
      diagnostic_data$shift_magnitude <- mean(abs(diagnostic_data$beta_final), na.rm = TRUE)
      diagnostic_data$relative_correction_mean <- mean(diagnostic_data$relative_correction_magnitude, na.rm = TRUE)
      
      # Write diagnostic CSV
      write.csv(diagnostic_data, output_file, row.names = FALSE)
      
      cat(sprintf("✓ Diagnostic parameters saved: %s (%d genes)\n", output_file, n_genes))
      
    }, error = function(e) {
      cat(sprintf("Warning: Could not save diagnostic parameters for %s: %s\n", method, e$message))
    })
  }
  
  apply_batch_corrections <- function(dat, batch, group, dat_test, method) {
    if (method == "unadjusted") {
      # No correction - return original data
      return(list(
        dat_corrected = dat,
        dat_test_corrected = dat_test
      ))
    } else if (method == "naive") {
      # Naive correction: match means and variances across batches
      # This is a simple batch correction that standardizes each batch to have the same mean and variance
      
      cat(sprintf("[NAIVE CORRECTION] Applying naive batch correction (mean/variance matching)\n"))
      
      # Step 1: Correct training data by standardizing each batch
      dat_corrected <- dat
      unique_batches <- unique(batch)
      n_batches <- length(unique_batches)
      
      cat(sprintf("  Processing %d training batches: %s\n", n_batches, paste(unique_batches, collapse=", ")))
      
      # Calculate overall statistics across all training data
      overall_mean <- rowMeans(dat)
      overall_var <- apply(dat, 1, var)
      
      # For each batch, standardize to match overall statistics
      for (b in unique_batches) {
        batch_idx <- which(batch == b)
        batch_data <- dat[, batch_idx, drop = FALSE]
        
        # Calculate batch-specific statistics
        batch_mean <- rowMeans(batch_data)
        batch_var <- apply(batch_data, 1, var)
        
        # Avoid division by zero for genes with no variance
        batch_sd <- sqrt(pmax(batch_var, 1e-10))
        overall_sd <- sqrt(pmax(overall_var, 1e-10))
        
        # Standardize: (x - batch_mean) / batch_sd * overall_sd + overall_mean
        for (i in 1:nrow(dat)) {
          if (batch_sd[i] > 1e-10) {
            dat_corrected[i, batch_idx] <- (batch_data[i, ] - batch_mean[i]) / batch_sd[i] * overall_sd[i] + overall_mean[i]
          } else {
            # If no variance in batch, just shift to overall mean
            dat_corrected[i, batch_idx] <- overall_mean[i]
          }
        }
        
        cat(sprintf("    Batch %d: %d samples corrected\n", b, length(batch_idx)))
      }
      
      # Step 2: Correct test data to match training data distribution
      # Use the overall training statistics as the target
      test_mean <- rowMeans(dat_test)
      test_var <- apply(dat_test, 1, var)
      test_sd <- sqrt(pmax(test_var, 1e-10))
      
      dat_test_corrected <- dat_test
      for (i in 1:nrow(dat_test)) {
        if (test_sd[i] > 1e-10) {
          dat_test_corrected[i, ] <- (dat_test[i, ] - test_mean[i]) / test_sd[i] * overall_sd[i] + overall_mean[i]
        } else {
          # If no variance in test data, just shift to overall training mean
          dat_test_corrected[i, ] <- overall_mean[i]
        }
      }
      
      cat(sprintf("  Test data: %d samples corrected to match training distribution\n", ncol(dat_test)))
      
      # Save diagnostic parameters
      diagnostic_file <- file.path(dirname(dirname(dirname(dirname(output_file)))), 
                                  "diagnostics", "naive", 
                                  sprintf("naive_%s_%s_%s_diagnostics.csv", 
                                         classifier, num_datasets, test_study))
      save_correction_diagnostics("naive", dat, dat_corrected, rownames(dat), diagnostic_file)
      
      return(list(
        dat_corrected = dat_corrected,
        dat_test_corrected = dat_test_corrected
      ))
    } else if (method == "rank_samples") {
      # Rank samples within genes, batch-aware
      cat(sprintf("[RANK_SAMPLES ADJUSTMENT] Applying rank_samples batch correction\n"))
      
      # Step 1: Apply rank_samples adjustment to training data
      dat_corrected <- adjust_ranked_samples_with_batch_info(dat, batch, debug = FALSE)
      
      cat(sprintf("  Training data: %d samples corrected using rank_samples adjustment\n", ncol(dat)))
      
      # Step 2: Apply rank_samples adjustment to test data
      # Treat test data as a single batch for ranking
      test_batch <- rep(1, ncol(dat_test))
      dat_test_corrected <- adjust_ranked_samples_with_batch_info(dat_test, test_batch, debug = FALSE)
      
      cat(sprintf("  Test data: %d samples corrected using rank_samples adjustment\n", ncol(dat_test)))
      
      # Save diagnostic parameters
      diagnostic_file <- file.path(dirname(dirname(dirname(dirname(output_file)))), 
                                  "diagnostics", "rank_samples", 
                                  sprintf("rank_samples_%s_%s_%s_diagnostics.csv", 
                                         classifier, num_datasets, test_study))
      save_correction_diagnostics("rank_samples", dat, dat_corrected, rownames(dat), diagnostic_file)
      
      return(list(
        dat_corrected = dat_corrected,
        dat_test_corrected = dat_test_corrected
      ))
    } else if (method == "rank_twice") {
      # Double ranking: genes within samples, then samples within genes, batch-aware
      cat(sprintf("[RANK_TWICE ADJUSTMENT] Applying rank_twice batch correction\n"))
      
      # Step 1: Apply rank_twice adjustment to training data
      dat_corrected <- adjust_ranked_twice_with_batch_info(dat, batch, debug = FALSE)
      
      cat(sprintf("  Training data: %d samples corrected using rank_twice adjustment\n", ncol(dat)))
      
      # Step 2: Apply rank_twice adjustment to test data
      # Treat test data as a single batch for ranking
      test_batch <- rep(1, ncol(dat_test))
      dat_test_corrected <- adjust_ranked_twice_with_batch_info(dat_test, test_batch, debug = FALSE)
      
      cat(sprintf("  Test data: %d samples corrected using rank_twice adjustment\n", ncol(dat_test)))
      
      # Save diagnostic parameters
      diagnostic_file <- file.path(dirname(dirname(dirname(dirname(output_file)))), 
                                  "diagnostics", "rank_twice", 
                                  sprintf("rank_twice_%s_%s_%s_diagnostics.csv", 
                                         classifier, num_datasets, test_study))
      save_correction_diagnostics("rank_twice", dat, dat_corrected, rownames(dat), diagnostic_file)
      
      return(list(
        dat_corrected = dat_corrected,
        dat_test_corrected = dat_test_corrected
      ))
    } else if (method == "npn") {
      # NPN (Nonparanormal) quantile normalization
      cat(sprintf("[NPN ADJUSTMENT] Applying Nonparanormal quantile normalization\n"))
      
      # Load required library
      if (!requireNamespace("huge", quietly = TRUE)) {
        stop("Package 'huge' is required for NPN adjustment but is not installed.")
      }
      
      # Step 1: Apply NPN adjustment to training data
      dat_corrected <- adjust_npn(dat, batch, debug = FALSE)
      
      cat(sprintf("  Training data: %d samples corrected using NPN transformation\n", ncol(dat)))
      
      # Step 2: Apply NPN adjustment to test data
      # Treat test data as a single batch for NPN transformation
      test_batch <- rep(1, ncol(dat_test))
      dat_test_corrected <- adjust_npn(dat_test, test_batch, debug = FALSE)
      
      cat(sprintf("  Test data: %d samples corrected using NPN transformation\n", ncol(dat_test)))
      
      # Save diagnostic parameters
      diagnostic_file <- file.path(dirname(dirname(dirname(dirname(output_file)))), 
                                  "diagnostics", "npn", 
                                  sprintf("npn_%s_%s_%s_diagnostics.csv", 
                                         classifier, num_datasets, test_study))
      save_correction_diagnostics("npn", dat, dat_corrected, rownames(dat), diagnostic_file)
      
      return(list(
        dat_corrected = dat_corrected,
        dat_test_corrected = dat_test_corrected
      ))
    } else if (method == "combat") {
      # ComBat correction without labels (unsupervised)
      # Step 1: Correct batch effects within training data without using labels
      
      # DIAGNOSTIC: Check data before ComBat
      cat(sprintf("[PRE-COMBAT DIAGNOSTIC] Data before ComBat:\n"))
      cat(sprintf("  dat shape: %d x %d\n", nrow(dat), ncol(dat)))
      cat(sprintf("  dat range: [%.3f, %.3f]\n", min(dat), max(dat)))
      cat(sprintf("  dat mean: %.3f, sd: %.3f\n", mean(dat), sd(dat)))
      cat(sprintf("  Sample values: [%.3f, %.3f, %.3f]\n", dat[1,1], dat[1,2], dat[1,3]))
      
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
      
      # Save diagnostic parameters for timidity analysis
      diagnostic_file <- file.path(dirname(dirname(dirname(dirname(output_file)))), 
                                  "diagnostics", "combat", 
                                  sprintf("combat_%s_%s_%s_diagnostics.csv", 
                                         classifier, num_datasets, test_study))
      save_correction_diagnostics("combat", dat, dat_corrected, rownames(dat), diagnostic_file)
      
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
    } else if (method == "fast_mnn") {
      library(batchelor, quietly = TRUE)
      library(SummarizedExperiment, quietly = TRUE)

      # FastMNN correction using batchelor::fastMNN
      cat("Applying FastMNN correction...\n")
      
      # Ensure dat and dat_test have the same genes (rows)
      if (nrow(dat) != nrow(dat_test)) {
        stop(sprintf("Training data has %d genes but test data has %d genes", nrow(dat), nrow(dat_test)))
      }
      
      # Prepare data for fastMNN - split training data by batch and add test as separate matrix
      batch_names <- unique(batch)
      batch_matrices <- list()
      
      for (b in batch_names) {
        batch_idx <- which(batch == b)
        batch_matrices[[length(batch_matrices) + 1]] <- dat[, batch_idx, drop = FALSE]
      }
      
      # Add test data as the final matrix
      batch_matrices[[length(batch_matrices) + 1]] <- dat_test
      
      cat(sprintf("FastMNN: Processing %d training batches + 1 test batch\n", length(batch_names)))
      cat(sprintf("  Matrix dimensions:\n"))
      for (i in seq_along(batch_matrices)) {
        cat(sprintf("    Matrix %d: %d x %d\n", i, nrow(batch_matrices[[i]]), ncol(batch_matrices[[i]])))
      }
      
      # Apply fastMNN - pass matrices as separate arguments
      tryCatch({
        # Use do.call to pass matrices as separate arguments
        fastmnn_result <- do.call(batchelor::fastMNN, batch_matrices)
        
        # Get the reconstructed matrix (FastMNN uses "reconstructed" not "corrected")
        corrected_matrix <- SummarizedExperiment::assay(fastmnn_result, "reconstructed")
        
        # Ensure the result is a proper matrix
        corrected_matrix <- as.matrix(corrected_matrix)
        
        # Validate matrix dimensions
        if (ncol(corrected_matrix) != (ncol(dat) + ncol(dat_test))) {
          stop(sprintf("FastMNN output has wrong dimensions: expected %d columns, got %d", 
                      ncol(dat) + ncol(dat_test), ncol(corrected_matrix)))
        }
        
        # Split corrected data back into training and test
        n_train_samples <- ncol(dat)
        dat_corrected <- corrected_matrix[, 1:n_train_samples]
        dat_test_corrected <- corrected_matrix[, (n_train_samples + 1):ncol(corrected_matrix)]
        
        cat(sprintf("FastMNN correction completed successfully\n"))
        
        return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))
      }, error = function(e) {
        cat(sprintf("[ERROR] fastMNN failed: %s\n", e$message))
        stop(sprintf("FastMNN failed: %s", e$message))
      })
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
    } else if (startsWith(method, "posse_")) {
      # POSSE: Pathway Optimized Shift and Scale Ensemble with different parameter settings
      posse_variant <- sub("posse_", "", method)
      cat(sprintf("Applying POSSE adjustment (%s variant)...\n", posse_variant))
      
      # Load reticulate for Python integration
      if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("reticulate package is required for POSSE method. Please install it.")
      }
      library(reticulate, quietly = TRUE)
      
      # Source the POSSE wrapper
      posse_wrapper_path <- "scripts/posse_wrapper.py"
      if (!file.exists(posse_wrapper_path)) {
        stop(sprintf("POSSE wrapper not found at: %s", posse_wrapper_path))
      }
      
      tryCatch({
        # Import the POSSE wrapper module
        source_python(posse_wrapper_path)
        
        # Apply POSSE correction
        # Convert R matrices to format expected by Python (genes x samples)
        
        # Create diagnostic file path - ENSURE DIRECTORY EXISTS
        # Extract base output directory (remove /results/adjusters/individual part)
        base_output_dir <- dirname(dirname(dirname(dirname(output_file))))
        diagnostic_dir <- file.path(base_output_dir, "diagnostics", "posse")
        
        # Create diagnostic directory with error handling - FAIL FAST
        tryCatch({
          if (!dir.exists(diagnostic_dir)) {
            dir.create(diagnostic_dir, recursive = TRUE)
          }
          
          # Verify directory was created and is writable
          if (!dir.exists(diagnostic_dir)) {
            stop(sprintf("CRITICAL: Failed to create diagnostic directory: %s", diagnostic_dir))
          }
          
          if (file.access(diagnostic_dir, 2) != 0) {
            stop(sprintf("CRITICAL: Diagnostic directory is not writable: %s", diagnostic_dir))
          }
          
        }, error = function(e) {
          stop(sprintf("CRITICAL: Cannot create diagnostic directory %s: %s", diagnostic_dir, e$message))
        })
        
        diagnostic_file <- file.path(diagnostic_dir, sprintf("%s_%s_%s_%s_diagnostics.csv", 
                                                           method, classifier, num_datasets, test_study))
        
        cat(sprintf("Diagnostic file will be saved to: %s\n", diagnostic_file))
        
        posse_result <- apply_posse_correction(
          train_data = dat,
          test_data = dat_test,
          train_batch = batch,
          posse_variant = posse_variant,
          gene_names = rownames(dat),
          save_diagnostics = diagnostic_file
        )
        
        # CRITICAL VALIDATION: Ensure diagnostic file was created
        if (!file.exists(diagnostic_file)) {
          stop(sprintf("CRITICAL: POSSE diagnostic file was not created: %s", diagnostic_file))
        }
        
        file_size <- file.info(diagnostic_file)$size
        if (is.na(file_size) || file_size == 0) {
          stop(sprintf("CRITICAL: POSSE diagnostic file is empty: %s", diagnostic_file))
        }
        
        cat(sprintf("✓ POSSE diagnostic file validated: %s (%d bytes)\n", diagnostic_file, file_size))
        
        dat_corrected <- posse_result$train_corrected
        dat_test_corrected <- posse_result$test_corrected
        
        # Log detailed POSSE metrics if available
        if (!is.null(posse_result$metadata)) {
          metadata <- posse_result$metadata
          cat(sprintf("POSSE %s detailed metrics:\n", posse_variant))
          cat(sprintf("  S_HK=%.4f\n", ifelse(is.null(metadata$S_HK), 1, metadata$S_HK)))
          cat(sprintf("  C_null_final=%.4f\n", ifelse(is.null(metadata$C_null_final), 0, metadata$C_null_final)))
          cat(sprintf("  Iterations=%d\n", ifelse(is.null(metadata$Iterations), 0, metadata$Iterations)))
          
          # Extract parameter statistics if available
          if (!is.null(metadata$Parameters)) {
            alpha_params <- metadata$Parameters[[1]]
            beta_params <- metadata$Parameters[[2]]
            cat(sprintf("  Alpha: mean=%.4f, std=%.4f, min=%.4f, max=%.4f\n",
                        mean(alpha_params, na.rm=TRUE),
                        sd(alpha_params, na.rm=TRUE),
                        min(alpha_params, na.rm=TRUE),
                        max(alpha_params, na.rm=TRUE)))
            cat(sprintf("  Beta: mean=%.4f, std=%.4f\n",
                        mean(beta_params, na.rm=TRUE),
                        sd(beta_params, na.rm=TRUE)))
          }
          
          # Diagnostic interpretation
          s_hk <- ifelse(is.null(metadata$S_HK), 1, metadata$S_HK)
          if (s_hk < 0.5 || s_hk > 2.0) {
            cat("  [DIAGNOSTIC] Large housekeeping scale factor - significant technical differences\n")
          } else {
            cat("  [DIAGNOSTIC] Housekeeping scale factor within normal range\n")
          }
        }
        
        cat(sprintf("POSSE %s adjustment complete\n", posse_variant))
        
      }, error = function(e) {
        cat(sprintf("[WARNING] POSSE %s correction failed: %s\n", posse_variant, e$message))
        cat("[WARNING] Falling back to unadjusted data\n")
        dat_corrected <- dat
        dat_test_corrected <- dat_test
      })
      
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
    train_mean <- mean(dat_train, na.rm = TRUE)
    train_sd <- sd(as.vector(dat_train), na.rm = TRUE)
    
    # Check for invalid scaling parameters
    if (is.na(train_mean) || is.na(train_sd) || train_sd == 0 || !is.finite(train_sd) || train_sd < 1e-10) {
      cat(sprintf("Global scaling: mean=%s, sd=%s\n", 
                  ifelse(is.na(train_mean), "NA", sprintf("%.4f", train_mean)),
                  ifelse(is.na(train_sd), "NA", sprintf("%.4f", train_sd))))
      
      # If standard deviation is too small, use unit scaling (subtract mean only)
      if (!is.na(train_mean) && is.finite(train_mean)) {
        cat("[WARNING] Using unit scaling (mean centering only) due to low variance\n")
        dat_train_scaled <- dat_train - train_mean
        dat_test_scaled <- dat_test - train_mean
        
        return(list(
          dat_train = dat_train_scaled,
          dat_test = dat_test_scaled
        ))
      } else {
        stop("Scaling produced invalid values (NA, 0, or infinite standard deviation)")
      }
    }
    
    # Apply same transformation to both train and test
    dat_train_scaled <- (dat_train - train_mean) / train_sd
    dat_test_scaled <- (dat_test - train_mean) / train_sd
    
    # Check for NaN/Inf in results
    if (any(!is.finite(dat_train_scaled)) || any(!is.finite(dat_test_scaled))) {
      stop("Scaling produced non-finite values in output data")
    }
    
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
  
  # Check if batch correction produced valid data
  if(all(dat_corrected == dat_corrected[1,1]) || all(dat_test_corrected == dat_test_corrected[1,1])) {
    cat("[WARNING] Batch correction produced constant data - using original data instead\n")
    dat_corrected <- dat
    dat_test_corrected <- dat_test
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
  
  # ====================================================================
  # FEATURE REDUCTION FOR HIGH-DIMENSIONAL DATA (e.g., after POSSE)
  # ====================================================================
  # POSSE uses all genes for pathway analysis, but classifiers like nnet
  # can't handle 10k+ features. Reduce to top 1000 most variable genes.
  
  max_features_for_classifier <- 1000
  if (n_features > max_features_for_classifier) {
    cat(sprintf("[FEATURE REDUCTION] Reducing from %d to %d features for classifier\n", 
                n_features, max_features_for_classifier))
    
    # Select top variable genes from training data
    gene_vars <- rowVars(dat_train_norm)
    top_gene_indices <- order(gene_vars, decreasing = TRUE)[1:max_features_for_classifier]
    
    # Apply to both train and test
    dat_train_norm <- dat_train_norm[top_gene_indices, ]
    dat_test_norm <- dat_test_norm[top_gene_indices, ]
    
    n_features <- nrow(dat_train_norm)
    cat(sprintf("[FEATURE REDUCTION] New feature count: %d\n", n_features))
  }
  
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
