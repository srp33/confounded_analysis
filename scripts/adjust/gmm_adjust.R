# GMM High-Level Interface and Legacy Compatibility
# This module provides user-facing functions and maintains backward compatibility

suppressPackageStartupMessages({
  library(tidyverse)
  library(mclust)
  library(data.table)
  library(foreach)
  library(doParallel)
  library(digest)
})

# Source core modules
source("scripts/adjust/gmm_parameters.R")
source("scripts/adjust/gmm_transforms.R")

# ============================================================================
# HIGH-LEVEL INTERFACE FUNCTIONS
# ============================================================================

#' Bimodal normalize with caching support (Main Interface)
#' 
#' @param data Input data matrix/data frame
#' @param cache_folder Optional cache folder path (cache file names generated deterministically)
#' @param force_recalculate Whether to force recalculation
#' @param debug Whether to enable debug logging
#' @param log_file Path to log file
#' @param adjustment_strategy Adjustment strategy to use
#' @param num_workers Number of workers to use. If NULL or 1, uses sequential processing.
#'                    If -1, uses all available cores. Otherwise uses minimum of specified number and available cores.
#' @return List with bimodal_data, recommended_modes, and performance info
bimodal_normalize_cached <- function(data, cache_folder = NULL, force_recalculate = FALSE,
                                   debug = FALSE, log_file = NULL, adjustment_strategy = "simple", num_workers = NULL) {
  
  # Get cached parameters
  param_result <- with_parameter_caching(
    data = data,
    batch = NULL,
    cache_folder = cache_folder,
    force_recalculate = force_recalculate,
    debug = debug,
    log_file = log_file,
    num_workers = num_workers
  )
  
  if (is.null(param_result) || is.null(param_result$gmm_params)) {
    if (debug) {
      message("Failed to get GMM parameters")
    }
    return(NULL)
  }

  if(debug) {
    message("DEBUG: Example GMM parameters retrieved:")
    message(param_result$gmm_params[1])
  }
  
  # Apply adjustments using cached parameters
  if (debug) {
    message("DEBUG: Applying adjustments using cached parameters")
  }
  adjustment_start <- Sys.time()
  
  adjusted_data <- apply_gmm_adjustment(data, param_result$gmm_params, adjustment_strategy, debug, log_file, num_workers = num_workers)
  
  adjustment_time <- as.numeric(difftime(Sys.time(), adjustment_start, units = "secs"))
  if (debug) {
    message("Adjustment completed in ", round(adjustment_time, 3), " seconds")
  }

  # Extract recommended modes from parameters
  recommended_modes <- sapply(param_result$gmm_params, function(p) p$recommended_modes)
  
  total_time <- param_result$cache_load_time_seconds + 
                param_result$extraction_time_seconds + 
                param_result$cahce_save_time_seconds + 
                adjustment_time
  
  return(list(
    bimodal_data = adjusted_data,
    recommended_modes = recommended_modes,
    cache_used = param_result$cache_used,
    performance_info = list(
      cache_load_time_seconds = param_result$cache_load_time_seconds,
      extraction_time_seconds = param_result$extraction_time_seconds,
      adjustment_time_seconds = adjustment_time,
      total_time_seconds = total_time
    )
  ))
}

#' GMM adjust with caching for batch processing (Main Interface)
#' 
#' @param data Input data matrix/data frame
#' @param batch Batch vector
#' @param debug Whether to enable debug logging
#' @param log_file Path to log file
#' @param adjustment_strategy Adjustment strategy
#' @param mixed_strategy Strategy for mixed scenarios
#' @param cache_folder Optional cache folder path (cache file names generated deterministically)
#' @param force_recalculate Whether to force recalculation
#' @param num_workers Number of workers to use. If NULL or 1, uses sequential processing.
#'                    If -1, uses all available cores. Otherwise uses minimum of specified number and available cores.
#' @return Adjusted data matrix
gmm_adjust_cached <- function(data, batch, debug = FALSE, log_file = NULL, 
                            adjustment_strategy = "simple", mixed_strategy = "unimodal",
                            cache_folder = NULL, force_recalculate = FALSE, num_workers = NULL) {
  
  if (debug) {
    message("DEBUG: Starting gmm_adjust_cached...")
    message("Input data dimensions:", nrow(data), "rows,", ncol(data), "columns")
    message("Input batch dimensions:", length(batch), "elements")
  }
  
  validate_inputs(data, batch, mixed_strategy)
  
  batch_factor <- as.factor(batch)
  batch_levels <- levels(batch_factor)
  
  n_samples <- nrow(data)
  n_genes <- ncol(data)
  
  unimodal_adjusted <- matrix(NA, nrow = n_samples, ncol = n_genes)
  bimodal_adjusted <- matrix(NA, nrow = n_samples, ncol = n_genes)
  
  recommended_modes_df <- data.frame(
    matrix(NA, nrow = length(batch_levels), ncol = n_genes),
    row.names = batch_levels
  )
  colnames(recommended_modes_df) <- colnames(data)
  
  for (b in batch_levels) {
    if (debug) {
      message("Processing batch ", b, " ...")
    }
    batch_indices <- which(batch == b)
    batch_data <- data[batch_indices, , drop = FALSE]
    
    # Use cached bimodal normalization with batch-specific cache folder
    bimodal_result <- bimodal_normalize_cached(
      batch_data, 
      cache_folder = cache_folder,
      force_recalculate = force_recalculate,
      debug = debug, 
      log_file = log_file, 
      adjustment_strategy = adjustment_strategy,
      num_workers = num_workers
    )
    
    if (is.null(bimodal_result)) {
      stop("Failed to process batch ", b)
    }
    
    batch_adjusted_unimodal <- unimodal_normalize(batch_data, debug = debug, num_workers = num_workers)
    batch_adjusted_bimodal <- bimodal_result$bimodal_data
    
    unimodal_adjusted[batch_indices, ] <- batch_adjusted_unimodal
    bimodal_adjusted[batch_indices, ] <- batch_adjusted_bimodal
    
    recommended_modes_df[b, ] <- bimodal_result$recommended_modes
  }
  
  # Apply final batch adjustment strategy
  adjusted_data <- gmm_batch_adjust(
    unimodal_adjusted, 
    bimodal_adjusted, 
    as.numeric(batch_factor), 
    recommended_modes_df, 
    debug, 
    log_file, 
    mixed_strategy
  )
  
  return(adjusted_data)
}

# ============================================================================
# LEGACY COMPATIBILITY FUNCTIONS
# ============================================================================

#' Logging wrapper that delegates to gmm_parameters.R
log_message <- function(..., log_file_path = NULL, iter = NULL, debug = TRUE) {
  if (!debug) return(invisible(NULL))
  
  if (!is.null(log_file_path)) {
    worker_log_message(..., iter = iter, log_file_path = log_file_path)
  } else {
    message(paste(...))
  }
}

#' Process single gene using parameter extraction from gmm_parameters.R
#' 
#' @param X Gene expression vector
#' @param iter Iteration number for logging
#' @param gene_id Gene identifier
#' @param debug Whether to enable debug logging
#' @param log_file Path to log file
#' @param strategy Adjustment strategy
#' @return List with unimodal, bimodal, and recommended_modes
process_single_gene <- function(X, iter, gene_id, debug, log_file, strategy) {
  # Use the parameter extraction function from gmm_parameters.R
  log_func <- function(..., iter = NULL) {
    worker_log_message(..., iter = iter, log_file_path = log_file)
  }
  
  gene_params <- process_single_gene_parameters(X, gene_id, debug, log_func, iter)
  
  # Generate unimodal fallback
  if (!is.numeric(X)) {
    X <- as.numeric(as.character(X))
  }
  
  if (all(is.na(X))) {
    return(list(unimodal = rep(0, length(X)), bimodal = rep(0, length(X)), recommended_modes = 1))
  }

  min_val <- min(X, na.rm = TRUE)
  X_transformed <- log(X - min_val + 1)
  
  unimodal_quantiles <- rank(X_transformed, na.last = "keep", ties.method = "average") / (sum(!is.na(X_transformed)) + 1)
  unimodal = qnorm(unimodal_quantiles)
  unimodal = unimodal/sd(unimodal, na.rm = TRUE)

  # If parameter extraction failed, return unimodal
  if (!gene_params$fit_successful) {
    return(list(unimodal = unimodal, bimodal = unimodal, recommended_modes = gene_params$recommended_modes))
  }
  
  # Apply strategy using cached parameters
  tryCatch({
    bimodal_result <- apply_adjustment_strategy(X, gene_params, strategy, debug, log_file)
    list(unimodal = unimodal, bimodal = bimodal_result, recommended_modes = gene_params$recommended_modes)
  }, error = function(e) {
    log_message(log_file_path=log_file, iter=iter, paste0("Error in process_single_gene for '", gene_id, "' with strategy '", strategy, "': ", e), debug = debug)
    return(list(unimodal = unimodal, bimodal = unimodal, recommended_modes = 1))
  })
}

#' Legacy bimodal normalize function (non-cached version)
#' 
#' This function provides backward compatibility for existing workflows.
#' For better performance, use bimodal_normalize_cached
#' 
#' @param data Input data matrix/data frame
#' @param debug Whether to enable debug logging
#' @param log_file Path to log file
#' @param adjustment_strategy Adjustment strategy
#' @param num_workers Number of workers to use. If NULL or 1, uses sequential processing.
#'                    If -1, uses all available cores. Otherwise uses minimum of specified number and available cores.
#' @return List with bimodal_data and recommended_modes
bimodal_normalize <- function(data, debug = FALSE, log_file = NULL, adjustment_strategy = "simple", num_workers = NULL) {
  validate_inputs(data)

  if (!is.null(log_file) && file.exists(log_file)) {
    file.remove(log_file)
  }

  cl <- setup_parallel(debug, num_workers)
  if (!is.null(cl)) {
    on.exit(stopCluster(cl), add = TRUE)
  }
  
  gene_names <- colnames(data)
  log_message(debug = debug, "Starting parallel processing...")
  start_time <- Sys.time()

  if (!is.null(cl)) {
    results_by_gene <- foreach(
      gene_name = gene_names,
      i = seq_along(gene_names)
    ) %dopar% {
      process_single_gene(
        data[, gene_name], i, gene_name, debug, log_file, adjustment_strategy
      )
    }
  } else {
    # Sequential processing
    results_by_gene <- list()
    for (i in seq_along(gene_names)) {
      gene_name <- gene_names[i]
      results_by_gene[[i]] <- process_single_gene(
        data[, gene_name], i, gene_name, debug, log_file, adjustment_strategy
      )
    }
  }

  end_time <- Sys.time()
  log_message(debug = debug, "Parallel processing took", round(difftime(end_time, start_time, units = "secs"), 1), "seconds")
  
  is_error <- function(x) inherits(x, "simpleError")
  
  bimodal_list <- lapply(results_by_gene, function(res) {
    if (is_error(res)) {
      log_message(debug = debug, "Error in process_single_gene:", res$message)
      return(rep(NA, nrow(data)))
    }
    res$bimodal
  })

  bimodal_data <- do.call(cbind, bimodal_list)
  
  recommended_modes <- sapply(results_by_gene, function(res) {
    if (is_error(res)) return(NA)
    res$recommended_modes
  })

  if (is.null(bimodal_data) || !is.matrix(bimodal_data) || !all(dim(bimodal_data) == dim(data))) {
    log_message(debug = debug, "Processing returned a NULL or incorrectly dimensioned object.")
    return(NULL)
  }

  colnames(bimodal_data) <- colnames(data)
  rownames(bimodal_data) <- rownames(data)
  
  return(list(
    bimodal_data = bimodal_data,
    recommended_modes = recommended_modes
  ))
}

#' Legacy multi-batch GMM adjustment function (non-cached version)
#' 
#' This function provides backward compatibility for existing workflows.
#' For better performance, use gmm_adjust_cached
#' 
#' @param data Input data matrix/data frame
#' @param batch Batch vector
#' @param debug Whether to enable debug logging
#' @param log_file Path to log file
#' @param adjustment_strategy Adjustment strategy
#' @param mixed_strategy Strategy for mixed scenarios
#' @param num_workers Number of workers to use. If NULL or 1, uses sequential processing.
#'                    If -1, uses all available cores. Otherwise uses minimum of specified number and available cores.
#' @return Adjusted data matrix
gmm_adjust <- function(data, batch, debug = FALSE, log_file = NULL, adjustment_strategy = "simple", mixed_strategy = "unimodal", num_workers = NULL) {
  log_message(debug = debug, "Starting gmm_adjust...")
  log_message(debug = debug, "Input data dimensions:", nrow(data), "rows,", ncol(data), "columns")
  log_message(debug = debug, "Input batch dimensions:", length(batch), "elements")

  validate_inputs(data, batch, mixed_strategy)

  batch_factor <- as.factor(batch)
  batch_levels <- levels(batch_factor)

  n_samples = nrow(data)
  n_genes = ncol(data)

  unimodal_adjusted = matrix(NA, nrow = n_samples, ncol = n_genes)
  bimodal_adjusted = matrix(NA, nrow = n_samples, ncol = n_genes)
  
  recommended_modes_df <- data.frame(
      matrix(NA, nrow = length(batch_levels), ncol = n_genes),
      row.names = batch_levels
  )
  colnames(recommended_modes_df) <- colnames(data)

  for (b in batch_levels) {
    log_message(debug = debug, "Processing batch ", b, " ...")
    batch_indices <- which(batch == b)
    batch_data <- data[batch_indices, , drop = FALSE]
    
    bimodal_result <- bimodal_normalize(batch_data, debug, log_file, adjustment_strategy, num_workers)
    batch_adjusted_unimodal <- unimodal_normalize(batch_data, debug = debug, num_workers = num_workers)
    batch_adjusted_bimodal <- bimodal_result$bimodal_data

    unimodal_adjusted[batch_indices, ] <- batch_adjusted_unimodal
    bimodal_adjusted[batch_indices, ] <- batch_adjusted_bimodal
    
    recommended_modes_df[b, ] <- bimodal_result$recommended_modes
    log_message(debug = debug, "For batch ", b, ", Number of 1 mode recommends: ", sum(recommended_modes_df[b, ] == 1), ", 2 mode recommends: ", sum(recommended_modes_df[b, ] == 2))
  }
  
  adjusted_data <- gmm_batch_adjust(
      unimodal_adjusted, 
      bimodal_adjusted, 
      as.numeric(batch_factor), 
      recommended_modes_df, 
      debug, 
      log_file, 
      mixed_strategy
  )

  return(adjusted_data)
}