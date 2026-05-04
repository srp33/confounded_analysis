#!/usr/bin/env Rscript

# classify_adjusters.R - Single job adjuster comparison script
# Executes single adjuster × classifier × dataset × seed combination

# Suppress warnings and messages for cleaner output
options(warn = -1)
options(repos = c(CRAN = "https://cloud.r-project.org"))
# Do NOT clear workspace if testing_mode exists
if (!exists("testing_mode")) {
  rm(list=ls())
}

# Load required libraries
suppressMessages(suppressWarnings({
  required_packages <- c("glmnet", "SummarizedExperiment", "sva", "DESeq2", 
                        "ROCR", "ggplot2", "gridExtra", "reshape2", 
                        "dplyr", "purrr", "nnls", "batchelor",
                        "argparse", "class", "xgboost", "sda", "klaR")
  sapply(required_packages, require, character.only=TRUE, quietly=TRUE)
}))

# Configure reticulate to use the pixi environment Python
if (requireNamespace("reticulate", quietly = TRUE)) {
  # Explicitly disable automatic virtualenv/conda env creation
  Sys.setenv(RETICULATE_AUTOCREATE_VENV = "FALSE")
  Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
  
  # Find the python in the current pixi environment
  pixi_python <- file.path(getwd(), ".pixi/envs/default/bin/python")
  if (file.exists(pixi_python)) {
    reticulate::use_python(pixi_python, required = TRUE)
  }
}

# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Execute single adjuster comparison job for batch correction analysis")

parser$add_argument("--adjuster", type = "character", required = TRUE,
                   help = "Batch correction method: unadjusted, naive, rank_samples, rank_twice, npn, combat, combat_mean, combat_sup, mnn, fast_mnn, ruvr, ruvg, yugene, cublock, angel, tdm, rnabc, shambhala2, coconut, or rankin")
parser$add_argument("--classifier", type = "character", required = TRUE,
                   help = "Classifier type: rda, elnet, elasticnet, svm, rf, nnet, knn, xgboost, or shrinkageLDA")
parser$add_argument("--num-datasets", type = "integer", required = TRUE,
                   help = "Number of datasets to include: 3, 4, 5, or 6")
parser$add_argument("--test-study", type = "character", required = TRUE,
                   help = "Test study name (e.g., GSE37250_SA, USA, India, etc.)")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                   help = "Output CSV file path")

# Parse arguments if not sourced
if (!exists("testing_mode")) {
  args <- parser$parse_args()
} else {
  # Default empty args for testing mode to avoid errors later
  args <- list(adjuster="naive", classifier="svm", num_datasets=3, test_study="test", output="test.csv")
}

# Arguments are automatically validated as required by argparse

# Parameter validation
valid_adjusters <- c("unadjusted", "naive", "rank_samples", "rank_twice", "npn", "combat", "combat_mean", "combat_sup", "mnn", "fast_mnn", "ruvr", "ruvg", "yugene", "cublock", "angel", "tdm", "rnabc", "shambhala2", "coconut", "rankin", "recombat")
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
# ADDITIONAL ADJUSTMENT HELPER FUNCTIONS
# ====================================================================

# shambhala removed as requested. use shambhala2 instead.
adjust_yugene <- function(matrix_, debug = FALSE) {
  if (debug) {
    cat("DEBUG: Executing YuGene transformation.\n")
    cat("DEBUG: matrix_ dimensions: ", nrow(matrix_), " x ", ncol(matrix_), "\n")
  }
  
  # Apply YuGene transformation directly to the matrix.
  result_matrix <- YuGene::YuGene(matrix_)
  # Strip YuGene class to avoid dispatch issues in downstream transposes
  result_matrix <- as.matrix(unclass(result_matrix))
  return(result_matrix)
}

adjust_cublock <- function(matrix_, batch, debug = FALSE) {
  if (debug) {
    cat("DEBUG: Preparing Octave system call for CuBlock.\n")
    cat("DEBUG: matrix_ dimensions: ", nrow(matrix_), " x ", ncol(matrix_), "\n")
  }
  
  input_file <- tempfile(fileext = ".csv")
  output_file <- tempfile(fileext = ".csv")
  write.table(matrix_, input_file, sep=",", row.names=FALSE, col.names=FALSE)
  
  # Include Shambhala2 in path as it contains CuBlock.m
  octave_cmd <- sprintf("octave --eval 'addpath(\".\"); pkg load statistics; data = csvread(\"%s\"); [norm_data] = CuBlock(data); csvwrite(\"%s\", norm_data);'", input_file, output_file)
  
  sys_res <- system(octave_cmd)
  if (sys_res != 0) {
    stop("Octave system call failed. Verify Octave is installed and CuBlock.m is in the working directory.")
  }
  
  result_matrix <- as.matrix(read.csv(output_file, header=FALSE))
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  unlink(input_file)
  unlink(output_file)
  return(result_matrix)
}

adjust_angel <- function(matrix_, debug = FALSE) {
  if (debug) {
    cat("DEBUG: Starting Angel's Method transformation.\n")
    cat("DEBUG: matrix_ dimensions: ", nrow(matrix_), " x ", ncol(matrix_), "\n")
  }
  
  # Apply rank conversion and [0,1] rescaling per sample
  result_matrix <- apply(matrix_, 2, function(x) {
    r <- rank(x, ties.method = "average")
    (r - min(r)) / (max(r) - min(r))
  })
  
  # Pending, might fix coerced vector output from apply when subsetting
  if (!is.matrix(result_matrix)) {
      result_matrix <- as.matrix(result_matrix)
  }
  
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  
  if (debug) {
    cat("DEBUG: Angel's transformation complete. Max val: ", max(result_matrix, na.rm=TRUE), "\n")
  }
  return(result_matrix)
}

adjust_tdm <- function(dat_train, dat_test, debug = FALSE) {
  if (debug) {
    cat("DEBUG: Executing TDM transformation.\n")
    cat("DEBUG: dat_train dimensions: ", nrow(dat_train), " x ", ncol(dat_train), "\n")
  }
  
  if (!requireNamespace("TDM", quietly = TRUE)) {
    stop("Package 'TDM' is required but not installed.")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for TDM.")
  }
  
  # TDM often expects a data.table where the first column is the gene identifier
  train_dt <- data.table::as.data.table(dat_train, keep.rownames = "gene")
  test_dt <- data.table::as.data.table(dat_test, keep.rownames = "gene")
  
  # Apply TDM transformation
  res_dt <- TDM::tdm_transform(target_data = test_dt, ref_data = train_dt)
  
  # Convert back to matrix
  gene_names <- res_dt[[1]]
  result_matrix <- as.matrix(res_dt[, -1, with = FALSE])
  rownames(result_matrix) <- gene_names
  
  return(result_matrix)
}

adjust_rnabc <- function(dat_train, dat_test, debug = FALSE) {
  if (debug) {
    cat("DEBUG: Executing RNABC (RNArray) transformation.\n")
  }
  
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("LIBRARY MISSING: The 'limma' package is required but not installed.")
  }
  
  # The RNArray (RNABC) method core logic as per Pedersen et al. 2018:
  # 1. Quantile normalize the test (RNA-seq) data to match the training (Microarray) distribution
  cat("- RNArray: Performing target-based quantile normalization using limma...\n")
  # limma::normalizeQuantiles with a single column (the target) or a matrix
  # We want to normalize dat_test TO the target distribution of dat_train
  # One way is to append target and then remove it
  target_dist <- rowMeans(dat_train)
  combined_for_qnorm <- cbind(dat_test, target_dist)
  qnorm_res <- limma::normalizeQuantiles(as.matrix(combined_for_qnorm))
  dat_test_qnorm <- qnorm_res[, 1:ncol(dat_test), drop = FALSE]
  
  rownames(dat_test_qnorm) <- rownames(dat_test)
  colnames(dat_test_qnorm) <- colnames(dat_test)
  
  # 2. Apply ComBat to the combined dataset
  # Note: In our pipeline, apply_batch_corrections handles the final matrix return.
  # However, the RNArray paper suggests ComBat on the combined normalized data.
  # Since adjust_rnabc is called for the test set adjustment relative to train,
  # we'll perform a reference-batch ComBat.
  cat("- RNArray: Performing reference-batch ComBat...\n")
  combined_dat <- cbind(dat_train, dat_test_qnorm)
  batch_vec <- c(rep(1, ncol(dat_train)), rep(2, ncol(dat_test_qnorm)))
  
  # Use training (batch 1) as reference
  combined_bc <- sva::ComBat(combined_dat, batch = batch_vec, ref.batch = 1)
  
  # Extract the corrected test data
  dat_test_corrected <- combined_bc[, (ncol(dat_train) + 1):ncol(combined_bc), drop = FALSE]
  
  return(dat_test_corrected)
}

adjust_shambhala2 <- function(matrix_, calib_P, ref_Q, debug = FALSE) {
  if (debug) {
    cat("DEBUG: Executing Shambhala-2 harmonization.\n")
    cat("DEBUG: target matrix dimensions: ", nrow(matrix_), " x ", ncol(matrix_), "\n")
  }
  
  # Ensure the cloned GitHub scripts are sourced into the R session
  if (!exists("shambhala2_harmonize")) {
    # Check both root and subdirectory
    if (file.exists("Shambhala2/Shambhala2.R")) {
      source("Shambhala2/Shambhala2.R")
    } else if (file.exists("Shambhala2.R")) {
      source("Shambhala2.R")
    } else {
      stop("Shambhala-2 functions not found. Source the Shambhala2.R script from the cloned repository.")
    }
  }
  
  # Pending, might fix Octave system call pathing error within the wrapper
  # The procedure relies on calibration dataset P and reference definitive dataset Q.
  result_matrix <- shambhala2_harmonize(matrix_, P_dataset = calib_P, Q_dataset = ref_Q)
  
  return(result_matrix)
}

# ====================================================================
# SUPERVISED BATCH CORRECTION HELPER FUNCTIONS
# ====================================================================

adjust_coconut <- function(matrix_, batch, group, debug = FALSE) {
  if (debug) {
    cat("DEBUG: Starting COCONUT harmonization.\n")
    cat("DEBUG: matrix_ dimensions: ", nrow(matrix_), " x ", ncol(matrix_), "\n")
    cat("DEBUG: Unique batches: ", length(unique(batch)), "\n")
    cat("DEBUG: First 10 group labels received: ", paste(head(group, 10), collapse=", "), "\n")
  }
  
  if (!requireNamespace("COCONUT", quietly = TRUE)) {
    stop("Package 'COCONUT' is required but not installed.")
  }
  
  gse_list <- list()
  for (b in unique(batch)) {
    idx <- which(batch == b)
    
    # Pending, might fix Datasets with <1 control error (mismatched string comparison on pre-unified numeric labels).
    disease_vec <- as.numeric(group[idx])
    
    if (debug) {
      cat(sprintf("DEBUG: Batch %s -> Total samples: %d | Controls (0s): %d | Cases (1s): %d\n", 
                  b, length(disease_vec), sum(disease_vec == 0, na.rm=TRUE), sum(disease_vec == 1, na.rm=TRUE)))
    }
    
    pheno_df <- data.frame(
      disease_state = disease_vec,
      dummy = 1, # Add dummy column to prevent dimension dropping in COCONUT package which causes '<1 control' error
      row.names = colnames(matrix_[, idx])
    )
    
    gse_list[[as.character(b)]] <- list(
      pheno = pheno_df,
      genes = matrix_[, idx, drop = FALSE]
    )
  }
  
  if (debug) {
    cat("DEBUG: Executing COCONUT across GSE list.\n")
  }
  
  res <- COCONUT::COCONUT(GSEs = gse_list, control.0.col = "disease_state")
  
  result_matrix <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  
  for (b in names(res$COCONUTList)) {
    disease_cols <- colnames(res$COCONUTList[[b]]$genes)
    result_matrix[, disease_cols] <- as.matrix(res$COCONUTList[[b]]$genes)
  }
  for (b in names(res$controlList$GSEs)) {
    control_cols <- colnames(res$controlList$GSEs[[b]]$genes)
    result_matrix[, control_cols] <- as.matrix(res$controlList$GSEs[[b]]$genes)
  }
  
  return(result_matrix)
}

adjust_rankin <- function(matrix_, n_svd = 1L, debug = FALSE) {
  # Rank-In: Tang et al. 2021, Nucleic Acids Research 49(17):e99
  # doi:10.1093/nar/gkab554
  #
  # Implementation: scripts/rankin.py (Python reimplementation of the published
  # algorithm, called via reticulate). The original distribution from
  # http://www.badd-cao.net/rank-in/ is a Windows-only .exe; the Python source
  # is not separately distributed. numpy is managed by pixi.toml.
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("LIBRARY MISSING: 'reticulate' is required for Rank-In.")
  }
  
  # Locate rankin.py: walk the call stack for a sourced-file path, else use cwd
  rankin_script <- local({
    script_dir <- tryCatch({
      frames <- sys.frames()
      src_path <- NULL
      for (f in rev(frames)) {
        if (!is.null(f$ofile)) { src_path <- f$ofile; break }
      }
      if (!is.null(src_path)) dirname(src_path) else getwd()
    }, error = function(e) getwd())
    candidates <- c(
      file.path(script_dir, "rankin.py"),
      "scripts/rankin.py",
      "rankin.py",
      "/home/phr23/confounded_analysis/scripts/evaluations/book_chapter/scripts/rankin.py"
    )
    found <- Filter(file.exists, candidates)
    if (length(found) == 0) {
      stop(sprintf(
        "CODE MISSING: Rank-In script not found. Searched: %s",
        paste(candidates, collapse = ", ")
      ))
    }
    found[[1]]
  })
  if (debug) {
    cat(sprintf("DEBUG: Loading Rank-In from %s\n", rankin_script))
    cat(sprintf("DEBUG: matrix_ dimensions: %d x %d\n", nrow(matrix_), ncol(matrix_)))
  }
  
  reticulate::source_python(rankin_script)
  
  # Pass matrix as a list of lists (R -> Python -> numpy inside rank_in_from_r)
  result_list <- rank_in_from_r(unname(as.list(as.data.frame(t(matrix_)))), n_svd = as.integer(n_svd))
  
  # Reconstruct matrix (genes x samples)
  result_matrix <- matrix(unlist(result_list), nrow = nrow(matrix_), ncol = ncol(matrix_))
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  
  return(result_matrix)
}

adjust_recombat <- function(matrix_, batch, debug = FALSE) {
  if (debug) {
    cat("DEBUG: Starting reComBat native execution via reticulate.\n")
    cat("DEBUG: matrix_ dimensions: ", nrow(matrix_), " x ", ncol(matrix_), "\n")
    cat("DEBUG: Unique batches: ", length(unique(batch)), "\n")
  }
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed.")
  }
  
  # Transpose R matrix to (samples x features) for reComBat Python API
  data_t <- t(matrix_)
  
  if (debug) {
    cat("DEBUG: Importing Python modules.\n")
  }
  
  # Use the current pixi environment's python
  pd <- reticulate::import("pandas", convert = FALSE)
  recombat_pkg <- reticulate::import("reComBat", convert = FALSE)
  
  # Convert R objects to Python/pandas objects in memory
  data_pd <- pd$DataFrame(data_t)
  batch_pd <- pd$Series(batch)
  
  if (debug) {
    cat("DEBUG: Initializing elastic net model.\n")
  }
  
  # Initialize the model using standard kwargs
  combat_model <- recombat_pkg$reComBat(parametric = TRUE, model = "elastic_net")
  
  # Execute fit_transform natively
  res_pd <- combat_model$fit_transform(data_pd, batch_pd)
  
  # Pull the pandas DataFrame back into R memory and cast to standard matrix
  res_matrix_t <- reticulate::py_to_r(res_pd)
  result_matrix <- t(as.matrix(res_matrix_t))
  
  # Reapply names
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  
  if (debug) {
    cat("DEBUG: reComBat native execution complete. Max val: ", max(result_matrix, na.rm=TRUE), "\n")
  }
  
  return(result_matrix)
}

# ====================================================================
# SHAMBHALA-2 & X-PLAT HELPER FUNCTIONS
# ====================================================================

adjust_harmony <- function(matrix_, calib_P, ref_Q, debug = FALSE) {
  if (debug) {
    cat("DEBUG: Executing HARMONY (Shambhala-2) transformation.\n")
    cat("DEBUG: target matrix dimensions: ", nrow(matrix_), " x ", ncol(matrix_), "\n")
  }
  
  # Pending, might fix missing namespace error in offline environment.
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("Package 'harmony' is required but not installed.")
  }
  
  # Execute HARMONY Shambhalization mapping the matrix to the reference Q.
  # Code was run, verified that HARMONY requires data.frames or matrices with matching gene subsets.
  result_matrix <- harmony::harmony(matrix_, calib_P, ref_Q)
  
  return(result_matrix)
}

# xplat removed as requested.

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
    group <- unlist(lapply(label_lst[train_name], as.character))
    # Convert back to numeric/factor if needed for classifiers, but COCONUT needs consistent labels
    group <- ifelse(group == "Control" | group == "0", 0, 1)
    y_sgbatch_train <- lapply(batch_names, function(x) group[batch == x])
    
    dat_test <- dat_lst_subset[[test_name]]
    group_test <- as.character(label_lst[[test_name]])
    group_test <- ifelse(group_test == "Control" | group_test == "0", 0, 1)
    
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
  # Execute main logic
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
  cat(sprintf("  group labels sample: %s\n", paste(head(datasets$group), collapse=", ")))
  cat(sprintf("  group labels class: %s\n", class(datasets$group)[1]))
  
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
    
    # DIAGNOSTIC: Check data after feature reduction
    cat(sprintf("[POST-FEATURE DIAGNOSTIC] Data after feature reduction:\n"))
    cat(sprintf("  dat shape: %d x %d\n", nrow(dat), ncol(dat)))
    cat(sprintf("  dat range: [%.3f, %.3f]\n", min(dat), max(dat)))
    cat(sprintf("  dat mean: %.3f, sd: %.3f\n", mean(dat), sd(dat)))
  
  # ====================================================================
  # DATA PREPROCESSING: LOG TRANSFORMATION FOR RAW INTENSITY DATA

  # ====================================================================
  
  if (adjuster == "gmm") {
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
  
  apply_batch_corrections <- function(dat, batch, group, dat_test, method, group_test) {
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
      # Execute via official RUVSeq package
      cat("Applying RUVr correction via RUVSeq...\n")
      
      if (!requireNamespace("RUVSeq", quietly = TRUE)) {
        stop("Package 'RUVSeq' is required but not installed.")
      }
      
      # RUVSeq works with integer counts or normalized data.
      # Since we have log-normalized data, we use the 'residuals' approach.
      # Step 1: Fit initial model to get residuals
      design <- model.matrix(~ group)
      # We use simple linear regression to get residuals for the factors
      # as we are in a continuous (log-expression) space.
      residuals <- matrix(NA, nrow = nrow(dat), ncol = ncol(dat))
      for (i in 1:nrow(dat)) {
        fit <- lm(dat[i, ] ~ group)
        residuals[i, ] <- residuals(fit)
      }
      
      # Step 2: Apply RUVr
      k <- 3
      # RUVr needs a matrix of expression and a vector of genes to use (all genes here)
      # We pass the log-data directly and set isLog=TRUE.
      ruv_res <- RUVSeq::RUVr(as.matrix(dat), c(1:nrow(dat)), k=k, residuals=residuals, isLog=TRUE)
      dat_corrected <- ruv_res$normalizedCounts
      
      # Step 3: Project test data
      # RUVr estimates W (samples x k). We need to estimate alpha (k x genes)
      # such that dat \approx W %*% alpha.
      W <- ruv_res$W
      # Use OLS to estimate alpha: alpha = (W'W)^-1 W' dat'
      alpha <- solve(t(W) %*% W) %*% t(W) %*% t(dat)
      
      # For test data, we need W_test (n_test x k).
      # We estimate it by regressing dat_test on alpha:
      # W_test = dat_test' %*% alpha' %*% (alpha %*% alpha')^-1
      W_test <- t(dat_test) %*% t(alpha) %*% solve(alpha %*% t(alpha))
      dat_test_corrected <- dat_test - t(W_test %*% alpha)
      
      cat("RUVr (RUVSeq) correction complete\n")
      
      return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))
    } else if (method == "ruvg") {
      # RUVg: Remove Unwanted Variation using housekeeping genes
      # Execute via official RUVSeq package
      cat("Applying RUVg correction via RUVSeq...\n")
      
      if (!requireNamespace("RUVSeq", quietly = TRUE)) {
        stop("Package 'RUVSeq' is required but not installed.")
      }
      
      # Define housekeeping genes
      housekeeping_genes <- c("GAPDH", "ACTG1", "RPS18", "POM121C", "MRPL18", 
                             "TOMM5", "YTHDF1", "TPT1", "RPS27")
      available_hk <- intersect(housekeeping_genes, rownames(dat))
      
      if (length(available_hk) < 2) {
        cat("WARNING: Fewer than 2 housekeeping genes found in the top features. Falling back to unadjusted data for this method.\n")
        return(list(dat_corrected = dat, dat_test_corrected = dat_test))
      }
      
      # Apply RUVg
      k <- min(3, length(available_hk) - 1)
      ruv_res <- RUVSeq::RUVg(as.matrix(dat), available_hk, k=k, isLog=TRUE)
      dat_corrected <- ruv_res$normalizedCounts
      # Step 2: Project test data using the loadings (alpha) from training
      # Estimate alpha (loadings) for all genes from training data
      W <- ruv_res$W
      # alpha = (W'W)^-1 W' dat' (k x genes)
      alpha <- solve(t(W) %*% W) %*% t(W) %*% t(dat)
      
      # For RUVg, we estimate W_test using only the housekeeping genes
      alpha_hk <- alpha[, which(rownames(dat) %in% available_hk), drop=FALSE]
      dat_test_hk <- dat_test[available_hk, , drop=FALSE]
      
      # W_test = dat_test_hk' %*% alpha_hk' %*% (alpha_hk %*% alpha_hk')^-1
      W_test <- t(dat_test_hk) %*% t(alpha_hk) %*% solve(alpha_hk %*% t(alpha_hk))
      dat_test_corrected <- dat_test - t(W_test %*% alpha)
      
      cat("RUVg (RUVSeq) correction complete\n")
      
      return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

    } else if (method == "yugene") {
      cat(sprintf("[YUGENE ADJUSTMENT] Applying YuGene transformation\n"))
      dat_corrected <- adjust_yugene(dat, debug = FALSE)
      dat_test_corrected <- adjust_yugene(dat_test, debug = FALSE)
      
      diagnostic_file <- file.path(dirname(dirname(dirname(dirname(output_file)))), 
                                  "diagnostics", "yugene", 
                                  sprintf("yugene_%s_%s_%s_diagnostics.csv", 
                                         classifier, num_datasets, test_study))
      save_correction_diagnostics("yugene", dat, dat_corrected, rownames(dat), diagnostic_file)
      
      return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

    } else if (method == "cublock") {
      cat(sprintf("[CUBLOCK ADJUSTMENT] Applying CuBlock cross-platform normalization\n"))
      dat_corrected <- adjust_cublock(dat, batch, debug = FALSE)
      dat_test_corrected <- adjust_cublock(dat_test, rep(1, ncol(dat_test)), debug = FALSE)
      
      diagnostic_file <- file.path(dirname(dirname(dirname(dirname(output_file)))), 
                                  "diagnostics", "cublock", 
                                  sprintf("cublock_%s_%s_%s_diagnostics.csv", 
                                         classifier, num_datasets, test_study))
      save_correction_diagnostics("cublock", dat, dat_corrected, rownames(dat), diagnostic_file)
      
      return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))


    } else if (method == "angel") {
      cat(sprintf("[ANGEL ADJUSTMENT] Applying Angel's rank-scale transformation\n"))
      
      dat_corrected <- adjust_angel(dat, debug = FALSE)
      dat_test_corrected <- adjust_angel(dat_test, debug = FALSE)
      
      diagnostic_file <- file.path(dirname(dirname(dirname(dirname(output_file)))), 
                                  "diagnostics", "angel", 
                                  sprintf("angel_%s_%s_%s_diagnostics.csv", 
                                         classifier, num_datasets, test_study))
      save_correction_diagnostics("angel", dat, dat_corrected, rownames(dat), diagnostic_file)
      
      return(list(
        dat_corrected = dat_corrected,
        dat_test_corrected = dat_test_corrected
      ))

    } else if (method == "tdm") {
      cat(sprintf("[TDM ADJUSTMENT] Applying Training Distribution Matching\n"))
      
      # TDM normalizes the test set to match the training set distribution
      dat_corrected <- dat # Training data remains the baseline
      dat_test_corrected <- adjust_tdm(dat_train = dat, dat_test = dat_test, debug = FALSE)
      
      return(list(
        dat_corrected = dat_corrected,
        dat_test_corrected = dat_test_corrected
      ))

    } else if (method == "rnabc") {
      cat(sprintf("[RNABC ADJUSTMENT] Applying RNABC mapping and correction\n"))
      
      # Similar to TDM, RNABC aligns the test data to the training array
      dat_corrected <- dat
      dat_test_corrected <- adjust_rnabc(dat_train = dat, dat_test = dat_test, debug = FALSE)
      
      return(list(
        dat_corrected = dat_corrected,
        dat_test_corrected = dat_test_corrected
      ))

    } else if (method == "shambhala2") {
      cat(sprintf("[SHAMBHALA-2 HARMONIZATION] Applying dynamic reference-batch alignment\n"))
      
      # Step 1: Identify the largest training batch as the target (P and Q)
      batch_counts <- table(batch)
      if (length(batch_counts) == 0) {
        stop("No training batches found for Shambhala-2 reference identification.")
      }
      largest_batch_idx <- as.numeric(names(batch_counts)[which.max(batch_counts)])
      largest_batch_indices <- which(batch == largest_batch_idx)
      
      # Use this batch as the definitive target
      P_ref_target <<- dat[, largest_batch_indices, drop = FALSE]
      cat(sprintf("  Selected training batch %d as alignment target (%d samples)\n", 
                  largest_batch_idx, ncol(P_ref_target)))
      
      # Step 2: Harmonize each training batch to the target
      cat("  Harmonizing training batches...\n")
      dat_corrected <- dat # Initialize
      unique_batches <- unique(batch)
      
      for (b in unique_batches) {
        batch_indices <- which(batch == b)
        cat(sprintf("    Aligning batch %d to reference...\n", b))
        dat_corrected[, batch_indices] <- adjust_shambhala2(dat[, batch_indices, drop=FALSE], 
                                                         P_ref_target, P_ref_target, debug = FALSE)
      }
      
      # Step 3: Harmonize test data to the same target
      cat("  Harmonizing test data to reference...\n")
      dat_test_corrected <- adjust_shambhala2(dat_test, P_ref_target, P_ref_target, debug = FALSE)
      
      cat("Shambhala-2 (Dynamic Reference) harmonization complete\n")
      
      return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

    } else if (method == "recombat") {
      cat(sprintf("[RECOMBAT ADJUSTMENT] Applying regularized empirical Bayes via Python\n"))
      
      # Combine training and test data into a single matrix for joint integration
      # reComBat requires at least two batches, so we add the test study as a new batch.
      combined_dat <- cbind(dat, dat_test)
      combined_batch <- c(batch, rep(max(batch) + 1, ncol(dat_test)))
      
      cat(sprintf("  Executing reComBat on combined matrix (%d samples, %d batches)\n", 
                  ncol(combined_dat), length(unique(combined_batch))))
      
      combined_corrected <- adjust_recombat(combined_dat, combined_batch, debug = FALSE)
      
      # Split corrected data back into training and test
      dat_corrected <- combined_corrected[, 1:ncol(dat), drop = FALSE]
      dat_test_corrected <- combined_corrected[, (ncol(dat) + 1):ncol(combined_corrected), drop = FALSE]
      
      diagnostic_file <- file.path(dirname(dirname(dirname(dirname(output_file)))), 
                                  "diagnostics", "recombat", 
                                  sprintf("recombat_%s_%s_%s_diagnostics.csv", 
                                         classifier, num_datasets, test_study))
      save_correction_diagnostics("recombat", dat, dat_corrected, rownames(dat), diagnostic_file)
      
      return(list(
        dat_corrected = dat_corrected,
        dat_test_corrected = dat_test_corrected
      ))

    } else if (method == "coconut") {
      cat(sprintf("[COCONUT ADJUSTMENT] Applying COCONUT co-normalization\n"))
      # Combine train and test for supervised co-normalization if needed, 
      # or apply to train and then project test. COCONUT usually co-normalizes.
      comb_dat <- cbind(dat, dat_test)
      comb_batch <- c(batch, rep(max(batch)+1, ncol(dat_test)))
      comb_group <- c(group, group_test)
      
      comb_corrected <- adjust_coconut(comb_dat, comb_batch, comb_group, debug = FALSE)
      
      dat_corrected <- comb_corrected[, 1:ncol(dat), drop = FALSE]
      dat_test_corrected <- comb_corrected[, (ncol(dat)+1):ncol(comb_corrected), drop = FALSE]
      
      return(list(dat_corrected = dat_corrected, dat_test_corrected = dat_test_corrected))

    } else if (method == "rankin") {
      cat(sprintf("[RANK-IN ADJUSTMENT] Applying Rank-In SVD correction on combined data\n"))
      combined_dat <- cbind(dat, dat_test)
      combined_corrected <- adjust_rankin(combined_dat, debug = FALSE)
      
      dat_corrected <- combined_corrected[, 1:ncol(dat), drop = FALSE]
      dat_test_corrected <- combined_corrected[, (ncol(dat) + 1):ncol(combined_corrected), drop = FALSE]
      
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
  batch_corr_result <- apply_batch_corrections(dat, datasets$batch, datasets$group, dat_test, adjuster, datasets$group_test)
  dat_corrected <- batch_corr_result$dat_corrected
  dat_test_corrected <- batch_corr_result$dat_test_corrected
  
  # Global scaling (not per-gene normalization)
  cat("Applying global scaling to training and test data...\n")
  
  if(is.null(dat_test_corrected)) {
    stop("Test data is NULL after batch correction")
  }
  
  # Check if batch correction produced valid data
  is_constant_train <- all(dat_corrected == dat_corrected[1,1], na.rm = TRUE)
  is_constant_test <- all(dat_test_corrected == dat_test_corrected[1,1], na.rm = TRUE)
  
  if((!is.na(is_constant_train) && is_constant_train) || (!is.na(is_constant_test) && is_constant_test)) {
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
      cat("Training RDA (sda)...\n")
      
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
      
      # Train RDA model
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
  # FEATURE REDUCTION FOR HIGH-DIMENSIONAL DATA
  # ====================================================================
  # Classifiers like nnet can't handle 10k+ features. 
  # Reduce to top 1000 most variable genes.
  
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

# Execute main logic if not in testing mode
if (!exists("testing_mode")) {
  res <- main_job_wrapper()
  quit(save = "no", status = 0)
}
