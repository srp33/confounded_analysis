# gmm_adjust.R

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

#' Logging wrapper that delegates to gmm_parameters.R
log_message <- function(..., log_file_path = NULL, iter = NULL, debug = TRUE) {
  if (!debug) return(invisible(NULL))
  
  if (!is.null(log_file_path)) {
    worker_log_message(..., iter = iter, log_file_path = log_file_path)
  } else {
    message(paste(...))
  }
}

#' Process single gene using pre-extracted GMM parameters
#' 
#' Note: GMM parameters are now always pre-extracted alize via extract_gmm_parameter
#' This function no longer handles plogging extraction as a fall
#' 
#' @param X Gene expression vector
#' @param iter Iteration number for logging
#' @param gene_id Gene identifier
#' @param debug Whether to enable debug logging
#' @param log_file Path to log file
#' @param gmm_params Pre-extracted GMM parameters (required)
process_single_gene <- function(X, iter, gene_id, debug, log_file, strategy, gmm_params= NULL) {
  if (is.null(gmm_params)) {
    stop("GMM parameters are NULL for gene '", gene_id, "'. This should not happen with the new parameter extraction approach.")
  }
  
  if (!"gene_name" %in% colnames(gmm_params)) {
    stop("GMM parameters missing 'gene_name' column. Structure: ", paste(colnames(gmm_params), collapse = ", "))
  }
  
  if (!gene_id %in% gmm_params$gene_name) {
    stop("GMM parameters not found for gene '", gene_id, "'. Available genes: ", length(unique(gmm_params$gene_name)), ", First few: ", paste(head(gmm_params$gene_name, 5), collapse = ", "))
  }
  
  gene_params <- gmm_params[gmm_params$gene_name == gene_id, , drop = FALSE]
  
  if (nrow(gene_params) == 0) {
    stop("No parameters returned for gene '", gene_id, "' after filtering.")
  }

  if (!is.numeric(X)) {
    X <- as.numeric(as.character(X))
  }
  if (all(is.na(X))) {
    return(list(unimodal = rep(0, length(X)), bimodal = rep(0, length(X)), recommended_modes = 1, gene_params=gene_params))
  }

  min_val <- min(X, na.rm = TRUE)
  X_transformed <- log(X - min_val + 1)
  
  unimodal_quantiles <- rank(X_transformed, na.last = "keep", ties.method = "average") / (sum(!is.na(X_transformed)) + 1)
  unimodal = qnorm(unimodal_quantiles)
  unimodal = unimodal/sd(unimodal, na.rm = TRUE)

  # If parameter extraction failed, return unimodal
  fit_successful <- as.logical(gene_params$fit_successful)
  if (is.na(fit_successful) || !fit_successful) {
    recommended_modes <- as.numeric(gene_params$recommended_modes)
    return(list(unimodal = unimodal, bimodal = unimodal, recommended_modes = recommended_modes, gene_params=gene_params))
  }
  
  # Apply strategy using cached parameters
  tryCatch({
    bimodal_result <- apply_adjustment_strategy(X, gene_params, strategy, debug, log_file)
    recommended_modes <- as.numeric(gene_params$recommended_modes)
    list(unimodal = unimodal, bimodal = bimodal_result, recommended_modes = recommended_modes, gene_params=gene_params)
  }, error = function(e) {
    log_message(log_file_path=log_file, iter=iter, paste0("Error in process_single_gene for '", gene_id, "' with strategy '", strategy, "': ", e), debug = debug)
    return(list(unimodal = unimodal, bimodal = unimodal, recommended_modes = 1, gene_params=gene_params))
  })
}



#' Bimodal normalize function
#' 
#' @param data Input data matrix/data frame
#' @param gmm_parameters Dataframe where columns are genes, rows are parameters (optional - will be extracted if NULL)
#' @param batch_name Name of the batch for caching (optional, used when gmm_parameters is NULL)
#' @param cache_folder Path to cache directory (optional, used when gmm_parameters is NULL)
#' @param debug Whether to enable debug logging
#' @param log_file Path to log file
#' @param adjustment_strategy Adjustment strategy
#' @param num_workers Number of workers to use. If NULL or 1, uses sequential processing.
#'                    If -1, uses all available cores. Otherwise uses minimum of specified number and available cores.
#' @return List with bimodal_data and recommended_modes
bimodal_normalize <- function(data, gmm_parameters=NULL, batch_name=NULL, cache_folder=NULL, debug = FALSE, log_file = NULL, adjustment_strategy = "simple", num_workers = NULL) {
  validate_inputs(data)

  if (!is.null(log_file) && file.exists(log_file)) {
    file.remove(log_file)
  }

  gene_names <- colnames(data)
  log_message(debug = debug, "Bimodal output going to log file:", log_file)
  
  # Extract GMM parameters if not provided
  if (is.null(gmm_parameters)) {
    log_message(debug = debug, "Extracting GMM parameters using intelligent caching...")
    gmm_parameters <- extract_gmm_parameters(
      data, 
      batch_name = batch_name,
      cache_folder = cache_folder,
      debug = debug,
      log_file = log_file,
      num_workers = num_workers
    )
    log_message(debug = debug, "GMM parameters extracted for", if(is.null(gmm_parameters)) 0 else nrow(gmm_parameters), "genes")
  } else {
    log_message(debug = debug, "Using provided GMM parameters for", nrow(gmm_parameters), "genes")
  }
  
  start_time <- Sys.time()

  # Setup parallel processing if requested
  if (num_workers != 1) {
    available_cores <- detectCores()
    num_cores <- if (num_workers == -1) available_cores else min(num_workers, available_cores)
    
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl), add = TRUE)
    
    # Export necessary functions to workers
    clusterExport(cl, c('process_single_gene', 'apply_adjustment_strategy', 
                       'worker_log_message', 'inverse_cdf_gmm_R', 'bimodal_npn'), envir = .GlobalEnv)
    
    results_by_gene <- foreach(
      gene_name = gene_names,
      i = seq_along(gene_names),
      .packages = c('mclust'),
      .errorhandling = 'stop'
    ) %dopar% {
      process_single_gene(
        data[, gene_name], i, gene_name, debug, log_file, adjustment_strategy, gmm_parameters
      )
    }
  } else {
    # Sequential processing
    results_by_gene <- list()
    for (i in seq_along(gene_names)) {
      gene_name <- gene_names[i]
      results_by_gene[[i]] <- process_single_gene(
        data[, gene_name], i, gene_name, debug, log_file, adjustment_strategy, gmm_parameters
      )
    }
  }

  end_time <- Sys.time()
  log_message(debug = debug, "Bimodal normalize took", round(difftime(end_time, start_time, units = "secs"), 1), "seconds")
  
  is_error <- function(x) inherits(x, "simpleError")
  
  bimodal_list <- lapply(results_by_gene, function(res) {
    if (is_error(res)) {
      log_message(debug = debug, "Error in process_single_gene:", res$message)
      return(rep(NA, nrow(data)))
    }
    if (is.null(res$bimodal)) {
      return(rep(NA, nrow(data)))
    }
    res$bimodal
  })
  
  # Filter out NULL elements
  bimodal_list <- bimodal_list[!sapply(bimodal_list, is.null)]
  
  if (length(bimodal_list) == 0) {
    log_message(debug = debug, "No valid bimodal results found")
    return(NULL)
  }
  
  bimodal_data <- do.call(cbind, bimodal_list)
  
  recommended_modes <- sapply(results_by_gene, function(res) {
    if (is_error(res)) return(NULL)
    res$recommended_modes
  })

  if (is.null(gmm_parameters)) {
    gmm_params_list <- lapply(results_by_gene, function(res) {
      if (is_error(res)) return(NULL)
      res$gene_params
    })
    gmm_params_list <- gmm_params_list[!sapply(gmm_params_list, is.null)]
    
    if (length(gmm_params_list) > 0) {
      gmm_parameters <- do.call(rbind, gmm_params_list)
    }
  } else {
    # When using cached parameters, we might need to update with any newly computed ones
    new_params_list <- lapply(results_by_gene, function(res) {
      if (is_error(res)) return(NULL)
      res$gene_params
    })
    new_params_list <- new_params_list[!sapply(new_params_list, is.null)]
    
    if (length(new_params_list) > 0) {
      new_params <- do.call(rbind, new_params_list)
      
      # Update gmm_parameters with any newly computed parameters
      for (gene in new_params$gene_name) {
        if (gene %in% gmm_parameters$gene_name) {
          # Update existing gene
          gmm_parameters[gmm_parameters$gene_name == gene, ] <- new_params[new_params$gene_name == gene, ]
        } else {
          # Add new gene
          gmm_parameters <- rbind(gmm_parameters, new_params[new_params$gene_name == gene, ])
        }
      }
    }
  }

  if (is.null(bimodal_data)) {
    log_message(debug = debug, "Processing returned a NULL bimodal_data.")
    return(NULL)
  }
  
  if (!is.matrix(bimodal_data)) {
    log_message(debug = debug, "bimodal_data is not a matrix, attempting to convert.")
    bimodal_data <- as.matrix(bimodal_data)
  }
  
  if (!all(dim(bimodal_data) == dim(data))) {
    log_message(debug = debug, "Dimension mismatch: bimodal_data", dim(bimodal_data), "vs data", dim(data))
    return(NULL)
  }

  colnames(bimodal_data) <- colnames(data)
  rownames(bimodal_data) <- rownames(data)
  
  return(list(
    bimodal_data = bimodal_data,
    recommended_modes = recommended_modes,
    gmm_parameters = gmm_parameters
  ))
}


#' Multi-batch GMM adjustment function
#' 
#' @param data Input data matrix/data frame
#' @param batch Batch vector
#' @param gmm_parameters Named list with structure {batch_name = dataframe}. Dataframe columns are genes, rows are parameters
#' @param cache_folder Path to cache directory (optional, for intelligent caching)
#' @param debug Whether to enable debug logging
#' @param log_file Path to log file
#' @param adjustment_strategy Adjustment strategy
#' @param mixed_strategy Strategy for mixed scenarios
#' @param num_workers Number of workers to use. If NULL or 1, uses sequential processing.
#'                    If -1, uses all available cores. Otherwise uses minimum of specified number and available cores.
#' @return Adjusted data matrix
gmm_adjust <- function(data, batch, gmm_parameters = list(), cache_folder = NULL, debug = FALSE, log_file = NULL, return_gmm_parameters=FALSE, adjustment_strategy = "simple", mixed_strategy = "unimodal", num_workers = NULL) {
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
    batch_indices <- which(batch == b)
    batch_data <- data[batch_indices, , drop = FALSE]
    
    batch_gmm_params = gmm_parameters[[b]]
    log_message(debug = debug, "Starting bimodal_normalize for batch", b, "...")
    bimodal_result <- bimodal_normalize(
      batch_data, 
      gmm_parameters = batch_gmm_params, 
      batch_name = if(is.null(batch_gmm_params)) as.character(b) else NULL,
      cache_folder = if(is.null(batch_gmm_params)) cache_folder else NULL,
      debug = debug, 
      log_file = log_file, 
      adjustment_strategy = adjustment_strategy, 
      num_workers = num_workers
    )
    log_message(debug = debug, "Starting unimodal_normalize for batch ", b, " ...")
    batch_adjusted_unimodal <- unimodal_normalize(batch_data, debug = debug, num_workers = num_workers)
    batch_adjusted_bimodal <- bimodal_result$bimodal_data

    unimodal_adjusted[batch_indices, ] <- batch_adjusted_unimodal
    bimodal_adjusted[batch_indices, ] <- batch_adjusted_bimodal
    gmm_parameters[[b]] <- bimodal_result$gmm_parameters
    
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

  if (return_gmm_parameters) {
    return(list(
      adjusted_data = adjusted_data,
      gmm_parameters = gmm_parameters
    ))
  }
  return(adjusted_data)
}