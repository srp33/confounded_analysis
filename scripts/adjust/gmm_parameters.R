# GMM Parameter Extraction and Caching Infrastructure
# This module provides functions to extract, save, and load GMM parameters
# to eliminate redundant model fitting across different adjustment strategies.

suppressPackageStartupMessages({
  library(mclust)
  library(digest)
  library(foreach)
  library(doParallel)
})

# ============================================================================
# HELPER FUNCTIONS FOR PARAMETER EXTRACTION
# ============================================================================

#' Process a single gene for GMM parameter extraction
#' 
#' @param X Numeric vector of gene expression values
#' @param gene_name Character string, name of the gene
#' @param debug Logical, whether to enable debug logging
#' @param log_func Logging function for worker processes (defaults to message-based logging)
#' @param iter Optional iteration number for logging
#' @return List containing GMM parameters for the gene
process_single_gene_parameters <- function(X, gene_name, debug = FALSE, 
                                         log_func = function(..., iter = NULL) if(debug) message(...), 
                                         iter = NULL) {
  # Initialize default parameters
  gene_params <- list(
    gene_name = gene_name,
    n_components = 1,
    means = NA,
    variances = NA,
    weights = NA,
    recommended_modes = 1,
    boundary_coefficients = list(A = NA, B = NA, C = NA),
    fit_successful = FALSE,
    error_message = NULL
  )
  
  tryCatch({
    # Ensure X is numeric
    if (!is.numeric(X)) {
      X <- as.numeric(as.character(X))
      log_func("Converted non-numeric data to numeric for gene", gene_name, iter = iter)
    }
    
    # Check for all NA values
    if (all(is.na(X))) {
      gene_params$error_message <- "All values are NA"
      log_func("Gene", gene_name, "has all NA values", iter = iter)
      return(gene_params)
    }
    
    # Apply log-shift transformation
    min_val <- min(X, na.rm = TRUE)
    X_transformed <- log(X - min_val + 1)
    
    # Check for identical values
    if (all(X_transformed == X_transformed[1])) {
      gene_params$error_message <- "All values are identical"
      log_func("Gene", gene_name, "has all identical values", iter = iter)
      return(gene_params)
    }
    
    # Check minimum data points
    if (sum(is.finite(X_transformed)) < 10) {
      gene_params$error_message <- "Insufficient finite values (< 10)"
      log_func("Gene", gene_name, "has insufficient finite values (<10)", iter = iter)
      return(gene_params)
    }
    
    # Fit GMM with 1 or 2 components
    log_func("Fitting GMM for gene", gene_name, iter = iter)
    gmm <- Mclust(X_transformed, G = 1:2, modelNames = "V", verbose = FALSE)
    
    if (is.null(gmm)) {
      gene_params$error_message <- "Mclust failed"
      log_func("Mclust failed for gene", gene_name, iter = iter)
      return(gene_params)
    }
    
    # Extract basic parameters
    gene_params$recommended_modes <- gmm$G
    gene_params$fit_successful <- TRUE
    
    log_func("GMM fit successful for gene", gene_name, "- components:", gmm$G, iter = iter)
    
    # If only one component, force fit with 2 for bimodal parameters
    if (gmm$G == 1) {
      gmm_bimodal <- Mclust(X_transformed, G = 2, modelNames = "V", verbose = FALSE)
      if (!is.null(gmm_bimodal)) {
        gmm <- gmm_bimodal  # Use 2-component fit for parameter extraction
      }
    }
    
    # Extract and sort parameters
    params <- gmm$parameters
    sort_idx <- order(params$mean)
    
    gene_params$n_components <- length(params$mean)
    gene_params$means <- params$mean[sort_idx]
    
    # Handle variances (can be list or vector)
    covs <- if (is.list(params$variance$sigmasq)) {
      unlist(params$variance$sigmasq)
    } else {
      params$variance$sigmasq
    }
    gene_params$variances <- covs[sort_idx]
    gene_params$weights <- params$pro[sort_idx]
    
    # Calculate boundary coefficients for 2-component case
    if (length(gene_params$means) >= 2) {
      m0 <- gene_params$means[1]
      m1 <- gene_params$means[2]
      v0 <- gene_params$variances[1]
      v1 <- gene_params$variances[2]
      w0 <- gene_params$weights[1]
      w1 <- gene_params$weights[2]
      
      # Calculate boundary coefficients
      A <- 1 / (2 * v0) - 1 / (2 * v1)
      B <- m1 / v1 - m0 / v0
      C <- (m0^2 / (2 * v0)) - (m1^2 / (2 * v1)) - log(w0 * sqrt(v1) / (w1 * sqrt(v0)))
      
      gene_params$boundary_coefficients <- list(A = A, B = B, C = C)
    }
    
    return(gene_params)
    
  }, error = function(e) {
    gene_params$error_message <- as.character(e)
    return(gene_params)
  })
}

#' Worker logging function for parallel processing
#' 
#' @param ... Messages to log
#' @param iter Optional iteration number for logging
#' @param log_file_path Path to log file
worker_log_message <- function(..., iter = NULL, log_file_path = NULL) {
  if (is.null(log_file_path)) return(invisible(NULL))
  
  msg <- paste0(
    format(Sys.time(), "%m-%d %H:%M:%S"),
    " | PID ", Sys.getpid(),
    if (!is.null(iter)) paste0(" | Gene ", iter) else "",
    " | ",
    paste(..., collapse = " "),
    "\n"
  )
  
  # Write to parallel log file
  tryCatch({
    con <- file(log_file_path, "a")
    on.exit(close(con))
    cat(msg, file = con)
  }, error = function(e) {
    # Silently continue if logging fails
  })
}

#' Setup parallel processing for parameter extraction
#' 
#' @param debug Logical, whether to enable debug logging
#' @param num_workers Number of workers to use. If NULL or 1, uses sequential processing.
#'                    If -1, uses all available cores. Otherwise uses minimum of specified number and available cores.
#' @return Cluster object
setup_parallel <- function(debug = FALSE, num_workers = NULL) {
  # Determine number of workers
  available_cores <- detectCores()
  
  if (is.null(num_workers) || num_workers == 1) {
    if (debug) message("Since num_workers is: ", num_workers, " using sequential version.")
    return(NULL)  # Use sequential processing
  } else if (num_workers == -1) {
    num_cores <- available_cores
  } else {
    num_cores <- min(num_workers, available_cores)
  }
  
  if (debug) {
    message("DEBUG: Setting up parallel processing with ", num_cores, " cores (", available_cores, " available)")
  }
  
  tryCatch({
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    # Export required functions and libraries to workers
    cluster_result <- clusterEvalQ(cl, {
      library(mclust)
      return("mclust loaded")
    })
    
    # Export the parameter processing function and logging to workers
    clusterExport(cl, c("process_single_gene_parameters", "worker_log_message"), envir = .GlobalEnv)
    
    return(cl)
  }, error = function(e) {
    stop("Failed to setup parallel cluster: ", e$message)
  })
}

# ============================================================================
# MAIN PARAMETER EXTRACTION FUNCTIONS
# ============================================================================

#' Extract GMM parameters using parallel processing (helper function)
#' 
#' @param data Matrix or data frame with samples as rows, genes as columns
#' @param debug Logical, whether to enable debug logging
#' @param num_workers Number of workers to use for parallel processing
#' @return List containing GMM parameters for each gene, or NULL if parallel processing fails
extract_gmm_parameters_parallel <- function(data, debug = FALSE, log_file = NULL, num_workers = NULL) {
  # Create parallel-specific log file
  parallel_log_file <- if (!is.null(log_file)) {
    gsub("\\.log$", "_parallel.log", log_file)
  } else {
    "outputs/gmm_parallel.log"
  }
  
  # Clear parallel log file if it exists
  if (file.exists(parallel_log_file)) {
    file.remove(parallel_log_file)
  }
  gene_names <- colnames(data)
  n_genes <- length(gene_names)
  
  if (debug) {
    message("DEBUG: Extracting GMM parameters for ", n_genes, " genes using parallel processing")
  }
  
  # Setup parallel processing
  cl <- tryCatch({
    setup_parallel(debug, num_workers)
  }, error = function(e) {
    if (debug) {
      message("DEBUG: Failed to setup parallel processing: ", e$message)
      message("DEBUG: Error details: ", toString(e))
    }
    return(NULL)
  })
  
  if (is.null(cl)) {
    return(NULL)
  }
  
  on.exit({
    if (debug) {
      message("Stopping parallel cluster")
    }
    stopCluster(cl)
  }, add = TRUE)
  
  if (debug) {
    message("DEBUG: Starting parallel parameter extraction...")
  }
  start_time <- Sys.time()
  if (debug) {
    message("DEBUG: Start time: ", format(start_time, "%m-%d %H:%M:%S"))
  }
  
  # Process genes in parallel
  results_by_gene <- tryCatch({
    if (debug) {
      message("DEBUG: Starting parameter foreach loop for ", n_genes, " genes")
    }
    
    foreach_result <- foreach(
      gene_name = gene_names,
      i = seq_along(gene_names),
      .combine = 'c',  # Use 'c' to combine into a simple list
      .multicombine = TRUE,
      .errorhandling = 'remove',  # Remove failed tasks but continue with others
      .packages = c('mclust')
    ) %dopar% {
      tryCatch({
        # Log start of gene processing
        worker_log_message("Starting parameter processing for gene ", gene_name, iter = i, log_file_path = parallel_log_file)
        
        X <- data[, gene_name]
        
        # Create a logging function for this worker
        worker_log_func <- function(..., iter = NULL) {
          worker_log_message(..., iter = iter, log_file_path = parallel_log_file)
        }
        
        gene_params <- process_single_gene_parameters(X, gene_name, debug = FALSE, 
                                                    log_func = worker_log_func, iter = i)
        
        # Log completion
        if (gene_params$fit_successful) {
          worker_log_message("Successfully processed gene ", gene_name, " - modes: ", gene_params$recommended_modes, 
                           iter = i, log_file_path = parallel_log_file)
        } else {
          worker_log_message("Failed to process gene ", gene_name, " - error: ", gene_params$error_message, 
                           iter = i, log_file_path = parallel_log_file)
        }
        
        # Return as a named list element
        result <- list(gene_params)
        names(result) <- gene_name
        return(result)
        
      }, error = function(e) {
        # Log worker error
        worker_log_message("ERROR: Worker failed for gene ", gene_name, ": ", e$message, 
                         iter = i, log_file_path = parallel_log_file)
        
        # Create error result
        error_result <- list(
          gene_name = gene_name,
          n_components = 1,
          means = NA,
          variances = NA,
          weights = NA,
          recommended_modes = 1,
          boundary_coefficients = list(A = NA, B = NA, C = NA),
          fit_successful = FALSE,
          error_message = paste("Parallel worker error:", e$message)
        )
        
        result <- list(error_result)
        names(result) <- gene_name
        return(result)
      })
    }
    
    if (debug) {
      message("DEBUG: Parameter foreach completed, got ", length(foreach_result), " results")
    }
    
    return(foreach_result)
    
  }, error = function(e) {
    if (debug) {
      message("DEBUG: Parallel parameter processing failed with error: ", e$message)
      message("DEBUG: Full error: ", toString(e))
    }
    return(NULL)
  })




  
  end_time <- Sys.time()
  if (debug) {
    processing_time <- round(difftime(end_time, start_time, units = "secs"), 1)
    message("DEBUG: Parallel parameter extraction completed in ", processing_time, " seconds")
  }
  start_time <- Sys.time()
  
  # Handle partial results from parallel processing
  if (is.null(results_by_gene)) {
    if (debug) {
      message("ERROR: Parallel parameter processing returned NULL results")
    }
    return(NULL)
  }
  
  if (debug) {
    message("Parallel parameter processing completed, got ", length(results_by_gene), " results for ", n_genes, " genes")
    if (length(results_by_gene) > 0) {
      message("DEBUG: Result names (first 10): ", paste(head(names(results_by_gene), 10), collapse = ", "))
    }
  }
  
  # Check if we got results for all genes
  actual_gene_names <- names(results_by_gene)
  missing_genes <- setdiff(gene_names, actual_gene_names)
  
  if (length(missing_genes) > 0) {
    if (debug) {
      message("DEBUG: Missing ", length(missing_genes), " genes from parallel processing")
      message("DEBUG: First 10 missing genes: ", paste(head(missing_genes, 10), collapse = ", "))
    }
    return(NULL)
  }
  
  end_time <- Sys.time()
  processing_time <- round(difftime(end_time, start_time, units = "secs"), 1)
  
  if (debug) {
    message("DEBUG: Parallel parameter processing completed in ", processing_time, " seconds")
  }
  
  return(results_by_gene)
}

#' Extract GMM parameters using sequential processing (helper function)
#' 
#' @param data Matrix or data frame with samples as rows, genes as columns
#' @param debug Logical, whether to enable debug logging
#' @param log_file Path to log file for debug output
#' @return List containing GMM parameters for each gene
extract_gmm_parameters_sequential <- function(data, debug = FALSE, log_file = NULL) {
  gene_names <- colnames(data)
  n_genes <- length(gene_names)
  
  if (debug) {
    message("DEBUG: Extracting GMM parameters for ", n_genes, " genes using sequential processing")
  }
  
  start_time <- Sys.time()
  gmm_params <- list()
  
  for (i in seq_along(gene_names)) {
    gene_name <- gene_names[i]
    X <- data[, gene_name]
    
    if (debug && i %% 100 == 0) {
      message("DEBUG: Processing gene ", i, "/", n_genes, ": ", gene_name)
    }
    
    gene_params <- process_single_gene_parameters(X, gene_name, debug = debug, iter = i)
    gmm_params[[gene_name]] <- gene_params
  }
  
  end_time <- Sys.time()
  processing_time <- round(difftime(end_time, start_time, units = "secs"), 1)
  
  if (debug) {
    message("DEBUG: Sequential parameter extraction completed in ", processing_time, " seconds")
    successful_fits <- sum(sapply(gmm_params, function(x) x$fit_successful))
    message("DEBUG: Successfully fitted GMM for ", successful_fits, "/", n_genes, " genes")
  }
  
  return(gmm_params)
}

#' Extract GMM parameters from data (main function)
#' 
#' Fits GMM models once per gene and extracts all parameters needed for
#' different adjustment strategies. Automatically chooses between parallel
#' and sequential processing based on data size and availability.
#' 
#' @param data Matrix or data frame with samples as rows, genes as columns
#' @param debug Logical, whether to enable debug logging
#' @param log_file Path to log file for debug output
#' @param use_parallel Logical, whether to use parallel processing. If NULL (default), 
#'                     automatically decides based on number of genes (parallel for >10 genes)
#' @param num_workers Number of workers to use. If NULL or 1, uses sequential processing.
#'                    If -1, uses all available cores. Otherwise uses minimum of specified number and available cores.
#' @return List containing GMM parameters for each gene
extract_gmm_parameters <- function(data, debug = FALSE, log_file = NULL, use_parallel = NULL, num_workers = NULL) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input 'data' must be a data frame or matrix.")
  }
  
  # Clear log file if it exists
  if (!is.null(log_file) && file.exists(log_file)) {
    file.remove(log_file)
  }
  
  gene_names <- colnames(data)
  n_genes <- length(gene_names)
  
  # Handle num_workers parameter
  if (!is.null(num_workers) && num_workers == 1) {
    use_parallel <- FALSE
  } else if (!is.null(num_workers)) {
    use_parallel <- TRUE
  }
  
  # Automatically decide on parallelization if not specified
  if (is.null(use_parallel)) {
    use_parallel <- n_genes > 10  # Use parallel for datasets with >10 genes
    if (debug) {
      message("DEBUG: Auto-selecting ", 
              ifelse(use_parallel, "parallel", "sequential"), 
              " parameter processing for ", n_genes, " genes")
    }
  }
  
  # Try parallel processing first if requested
  if (use_parallel) {
    parallel_result <- extract_gmm_parameters_parallel(data, debug, log_file, num_workers)
    
    if (!is.null(parallel_result)) {
      return(parallel_result)
    } else {
      if (debug) {
        message("DEBUG: Parallel parameter processing failed, falling back to sequential processing")
      }
    }
  }
  
  # Use sequential processing (either requested or as fallback)
  return(extract_gmm_parameters_sequential(data, debug, log_file))
}

# ============================================================================
# PARAMETER CACHING UTILITIES
# ============================================================================

#' Generate cache key for GMM parameters
#' 
#' @param data Input data matrix/data frame
#' @param batch Optional batch vector
#' @return Character string cache key
generate_cache_key <- function(data, batch = NULL) {
  # Create a hash based on data dimensions and content sample
  # Deterministic due to seed
  set.seed(42)
  data_info <- list(
    nrow = nrow(data),
    ncol = ncol(data),
    colnames = colnames(data),
    sample_values = if (nrow(data) > 0 && ncol(data) > 0) {
      sample_rows <- sample(min(10, nrow(data)), min(10, nrow(data)))
      sample_cols <- sample(min(10, ncol(data)), min(10, ncol(data)))
      as.vector(data[sample_rows, sample_cols])
    } else {
      NULL
    },
    batch_info = if (!is.null(batch)) {
      list(length = length(batch), levels = unique(batch))
    } else {
      NULL
    }
  )
  
  return(digest(data_info, algo = "md5"))
}

#' Generate deterministic cache file path
#' 
#' @param data Input data matrix/data frame
#' @param cache_folder Cache folder path
#' @param batch Optional batch vector for batch-specific caching
#' @return Character string cache file path
generate_cache_file_path <- function(data, cache_folder, batch = NULL) {
  if (is.null(cache_folder)) return(NULL)
  
  # Generate cache key based on data characteristics
  cache_key <- generate_cache_key(data, batch)
  
  # Create deterministic filename
  filename <- paste0(cache_key, ".rds")
  
  # Return full path
  return(file.path(cache_folder, filename))
}

#' Save GMM parameters to cache file
#' 
#' @param gmm_params List of GMM parameters
#' @param cache_file Path to cache file
#' @param debug Whether to enable debug logging
save_gmm_cache <- function(gmm_params, cache_file, debug = FALSE) {
  if (is.null(cache_file)) return(invisible(NULL))
  
  tryCatch({
    # Create directory if it doesn't exist
    cache_dir <- dirname(cache_file)
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
    }
    
    # Save parameters with metadata
    cache_data <- list(
      gmm_params = gmm_params,
      timestamp = Sys.time(),
      version = "1.0"
    )
    
    saveRDS(cache_data, cache_file)
    
    if (debug) {
      message("Saved GMM parameters to cache: ", cache_file)
    }
  }, error = function(e) {
    if (debug) {
      message("Failed to save cache: ", e$message)
    }
  })
}

#' Load GMM parameters from cache file
#' 
#' @param cache_file Path to cache file
#' @param debug Whether to enable debug logging
#' @return List of GMM parameters or NULL if loading fails
load_gmm_cache <- function(cache_file, debug = FALSE) {
  if (is.null(cache_file) || !file.exists(cache_file)) {
    if (debug) {
      message("DEBUG: Cache file: ", cache_file, " does not exist.")
    }
    return(NULL)
  }
  
  tryCatch({
    cache_data <- readRDS(cache_file)
    
    if (debug) {
      message("Loaded GMM parameters from cache: ", cache_file)
    }
    return(cache_data$gmm_params)
  }, error = function(e) {
    if (debug) {
      message("Failed to load cache: ", e$message)
    }
    return(NULL)
  })
}

#' Extract GMM parameters with caching support
#' 
#' @param data Input data matrix/data frame
#' @param batch Optional batch vector
#' @param cache_folder Optional cache folder path (preferred method)
#' @param force_recalculate Whether to force recalculation
#' @param debug Whether to enable debug logging
#' @param log_file Path to log file
#' @param num_workers Number of workers to use for parallel processing
#' @return List with gmm_params, cache_used, and timing information
with_parameter_caching <- function(data, batch = NULL, cache_folder = NULL,
                                 force_recalculate = FALSE, debug = FALSE, log_file = NULL, num_workers = NULL) {
  
  cache_load_start <- Sys.time()
  
  # Determine cache file path
  if (!is.null(cache_folder)) {
    cache_file <- generate_cache_file_path(data, cache_folder, batch)
  }
  
  # Try to load from cache first
  gmm_params <- NULL
  cache_used <- FALSE

  if (force_recalculate) {
    if (debug) {
      message("DEBUG: Forced recalculation requested")
    }
  }
  else if (!is.null(cache_file)) {
    gmm_params <- load_gmm_cache(cache_file, debug)
    if (!is.null(gmm_params)) {
      cache_used <- TRUE
      if (debug) {
        message("DEBUG: Using cached GMM parameters")
      }
    }
    else if(debug) {
      message("DEBUG: No cached GMM parameters found.")
    }
  }
  
  cache_load_time <- as.numeric(difftime(Sys.time(), cache_load_start, units = "secs"))
  
  # Extract parameters if not cached
  extraction_time <- 0
  cache_save_time <- 0
  if (is.null(gmm_params)) {
    
    extraction_start <- Sys.time()
    gmm_params <- extract_gmm_parameters(data, debug, log_file, use_parallel = NULL, num_workers = num_workers)
    extraction_time <- as.numeric(difftime(Sys.time(), extraction_start, units = "secs"))
    
    # Save to cache
    cache_save_start = Sys.time()
    if (!is.null(cache_file) && !is.null(gmm_params)) {
      save_gmm_cache(gmm_params, cache_file, debug)
    }
    cache_save_time <- as.numeric(difftime(Sys.time(), cache_save_start, units = "secs"))
  }
  
  return(list(
    gmm_params = gmm_params,
    cache_used = cache_used,
    cache_load_time_seconds = cache_load_time,
    extraction_time_seconds = extraction_time,
    cache_save_time_seconds = cache_save_time
  ))
}