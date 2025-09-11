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
  # Initialize default parameters as dataframe with single column (flattened structure)
  gene_params <- data.frame(
    c(gene_name,1,NA,NA,NA,NA,NA,NA,1,FALSE,NA),
    row.names = c(
      "gene_name",
      "n_components", 
      "mean0",
      "mean1", 
      "variance0",
      "variance1",
      "weight0",
      "weight1",
      "recommended_modes",
      "fit_successful",
      "error_message"
    ),
    stringsAsFactors = FALSE
  )
  colnames(gene_params) <- gene_name
  
  tryCatch({
    # Ensure X is numeric
    if (!is.numeric(X)) {
      X <- as.numeric(as.character(X))
      log_func("Converted non-numeric data to numeric for gene", gene_name, iter = iter)
    }
    
    # Check for all NA values
    if (all(is.na(X))) {
      gene_params["error_message", gene_name] <- "All values are NA"
      log_func("Gene", gene_name, "has all NA values", iter = iter)
      return(gene_params)
    }
    
    # Apply log-shift transformation
    min_val <- min(X, na.rm = TRUE)
    X_transformed <- log(X - min_val + 1)
    
    # Check for identical values
    if (all(X_transformed == X_transformed[1])) {
      gene_params["error_message", gene_name] <- "All values are identical"
      log_func("Gene", gene_name, "has all identical values", iter = iter)
      return(gene_params)
    }
    
    # Check minimum data points
    if (sum(is.finite(X_transformed)) < 10) {
      gene_params["error_message", gene_name] <- "Insufficient finite values (< 10)"
      log_func("Gene", gene_name, "has insufficient finite values (<10)", iter = iter)
      return(gene_params)
    }
    
    # Fit GMM with 1 or 2 components
    log_func("Fitting GMM for gene", gene_name, iter = iter)
    gmm <- Mclust(X_transformed, G = 1:2, modelNames = "V", verbose = FALSE)
    
    if (is.null(gmm)) {
      gene_params["error_message", gene_name] <- "Mclust failed"
      log_func("Mclust failed for gene", gene_name, iter = iter)
      return(gene_params)
    }
    
    # Extract basic parameters
    gene_params["recommended_modes", gene_name] <- gmm$G
    gene_params["fit_successful", gene_name] <- TRUE
    
    log_func("GMM fit successful for gene", gene_name, "- components:", gmm$G, iter = iter)
    
    # If only one component, force fit with 2 for bimodal parameters
    if (gmm$G == 1) {
      gmm_bimodal <- Mclust(X_transformed, G = 2, modelNames = "V", verbose = FALSE)
      if (!is.null(gmm_bimodal)) {
        gmm <- gmm_bimodal  # Use 2-component fit for parameter extraction
      }
    }
    
    # Extract and sort parameters (flattened structure)
    params <- gmm$parameters
    sort_idx <- order(params$mean)
    
    gene_params["n_components", gene_name] <- length(params$mean)
    
    # Handle variances (can be list or vector)
    covs <- if (is.list(params$variance$sigmasq)) {
      unlist(params$variance$sigmasq)
    } else {
      params$variance$sigmasq
    }
    
    # Store flattened parameters
    sorted_means <- params$mean[sort_idx]
    sorted_variances <- covs[sort_idx]
    sorted_weights <- params$pro[sort_idx]
    
    gene_params["mean0", gene_name] <- sorted_means[1]
    gene_params["mean1", gene_name] <- if(length(sorted_means) > 1) sorted_means[2] else NA
    gene_params["variance0", gene_name] <- sorted_variances[1]
    gene_params["variance1", gene_name] <- if(length(sorted_variances) > 1) sorted_variances[2] else NA
    gene_params["weight0", gene_name] <- sorted_weights[1]
    gene_params["weight1", gene_name] <- if(length(sorted_weights) > 1) sorted_weights[2] else NA
    
    return(gene_params)
    
  }, error = function(e) {
    gene_params["error_message", gene_name] <- as.character(e)
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
# PARAMETER CACHING FUNCTIONS
# ============================================================================

#' Save GMM parameters to CSV file
#' 
#' @param gmm_params Dataframe with genes as columns and parameters as rows
#' @param batch_name Name of the batch (used in filename)
#' @param cache_folder Path to cache directory
save_gmm_parameters_to_csv <- function(gmm_params, batch_name, cache_folder) {
  if (!dir.exists(cache_folder)) {
    dir.create(cache_folder, recursive = TRUE)
  }
  
  cache_file <- file.path(cache_folder, paste0(batch_name, "_gmm_params.csv"))
  
  write.csv(gmm_params, cache_file, row.names = TRUE)
  message("Saved GMM parameters for batch '", batch_name, "' to: ", cache_file)
}

#' Load GMM parameters from CSV file
#' 
#' @param batch_name Name of the batch (used in filename)
#' @param cache_folder Path to cache directory
#' @return Dataframe with parameters or NULL if file doesn't exist
load_gmm_parameters_from_csv <- function(batch_name, cache_folder) {
  cache_file <- file.path(cache_folder, paste0(batch_name, "_gmm_params.csv"))
  
  if (!file.exists(cache_file)) {
    return(NULL)
  }
  
  params <- read.csv(cache_file, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Convert columns back to proper types
  if ("fit_successful" %in% rownames(params)) {
    params["fit_successful", ] <- as.logical(params["fit_successful", ])
  }
  
  numeric_rows <- c("n_components", "mean0", "mean1", "variance0", "variance1", "weight0", "weight1", "recommended_modes")
  for (row_name in numeric_rows) {
    if (row_name %in% rownames(params)) {
      params[row_name, ] <- as.numeric(params[row_name, ])
    }
  }
  
  return(params)
}

#' Validate cached parameters against current gene set
#' 
#' @param cached_params Cached parameter dataframe
#' @param current_genes Vector of current gene names
#' @return Logical indicating if cache is valid
validate_cached_parameters <- function(cached_params, current_genes) {
  if (is.null(cached_params)) {
    return(FALSE)
  }
  
  # Check if all current genes have cached parameters
  if (!all(current_genes %in% colnames(cached_params))) {
    return(FALSE)
  }
  
  return(TRUE)
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
      .combine = 'c',
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
        if (gene_params["fit_successful", gene_name]) {
          worker_log_message("Successfully processed gene ", gene_name, " - modes: ", gene_params["recommended_modes", gene_name], 
                           iter = i, log_file_path = parallel_log_file)
        } else {
          worker_log_message("Failed to process gene ", gene_name, " - error: ", gene_params["error_message", gene_name], 
                           iter = i, log_file_path = parallel_log_file)
        }
        
        # Return as a dataframe column
        result <- list(gene_params)
        names(result) <- gene_name
        return(result)
        
      }, error = function(e) {
        # Log worker error
        worker_log_message("ERROR: Worker failed for gene ", gene_name, ": ", e$message, 
                         iter = i, log_file_path = parallel_log_file)
        
        # Create error result (flattened structure)
        error_result <- data.frame(
          c(gene_name,1,NA,NA,NA,NA,NA,NA,1,FALSE,paste("Parallel worker error:", e$message)),
          row.names = c(
            "gene_name",
            "n_components", 
            "mean0",
            "mean1", 
            "variance0",
            "variance1",
            "weight0",
            "weight1",
            "recommended_modes",
            "fit_successful",
            "error_message"
          ),
          stringsAsFactors = FALSE
        )
        colnames(error_result) <- gene_name
        
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
    successful_fits <- sum(sapply(gmm_params, function(x) x["fit_successful", colnames(x)[1]]))
    message("DEBUG: Successfully fitted GMM for ", successful_fits, "/", n_genes, " genes")
  }
  
  return(gmm_params)
}

#' Extract GMM parameters from data (main function)
#' 
#' Fits GMM models once per gene and extracts all parameters needed for
#' different adjustment strategies. 
#' 
#' @param data Matrix or data frame with samples as rows, genes as columns
#' @param debug Logical, whether to enable debug logging
#' @param log_file Path to log file for debug output
#' @param num_workers Number of workers to use. If NULL or 1, uses sequential processing.
#'                    If -1, uses all available cores. Otherwise uses minimum of specified number and available cores.
#' @return List containing GMM parameters for each gene
extract_gmm_parameters <- function(data, debug = FALSE, log_file = NULL, num_workers = 1) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input 'data' must be a data frame or matrix.")
  }
  
  # Clear log file if it exists
  if (!is.null(log_file) && file.exists(log_file)) {
    file.remove(log_file)
  }
  
  gene_names <- colnames(data)
  n_genes <- length(gene_names)
  
  # Try parallel processing first if requested
  if (num_workers != 1) {
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
