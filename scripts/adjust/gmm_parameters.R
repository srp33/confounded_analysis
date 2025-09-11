# GMM Parameter Extraction and Caching Infrastructure
# This module provides functions to extract, save, and load GMM parameters
# to eliminate redundant model fitting across different adjustment strategies.
# 
# NEW STRUCTURE: Genes as rows, parameters as columns (more intuitive)

suppressPackageStartupMessages({
  library(mclust)
  library(digest)
  library(foreach)
  library(doParallel)
})

# ============================================================================
# CONSTANTS AND HELPER FUNCTIONS
# ============================================================================

# Define parameter structure once to avoid duplication
GMM_PARAM_COLS <- c(
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
)

# Define which columns should be numeric for type conversion
NUMERIC_PARAM_COLS <- c("n_components", "mean0", "mean1", "variance0", "variance1", "weight0", "weight1", "recommended_modes")

#' Create default parameter dataframe for a gene
#' 
#' @param gene_name Character string, name of the gene
#' @param error_msg Optional error message
#' @return Dataframe with default parameters (gene as row, parameters as columns)
create_default_gene_params <- function(gene_name, error_msg = NA) {
  gene_params <- data.frame(
    gene_name = gene_name,
    n_components = 1,
    mean0 = NA_real_,
    mean1 = NA_real_,
    variance0 = NA_real_,
    variance1 = NA_real_,
    weight0 = NA_real_,
    weight1 = NA_real_,
    recommended_modes = 1,
    fit_successful = FALSE,
    error_message = as.character(error_msg),
    stringsAsFactors = FALSE
  )
  rownames(gene_params) <- gene_name
  return(gene_params)
}

#' Convert parameter dataframe columns to proper types
#' 
#' @param params Parameter dataframe to convert
#' @return Dataframe with properly typed columns
convert_param_types <- function(params) {
  # Convert logical columns
  if ("fit_successful" %in% colnames(params)) {
    params$fit_successful <- as.logical(params$fit_successful)
  }
  
  # Convert numeric columns
  for (col_name in NUMERIC_PARAM_COLS) {
    if (col_name %in% colnames(params)) {
      params[[col_name]] <- as.numeric(params[[col_name]])
    }
  }
  
  return(params)
}

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
#' @return Dataframe with GMM parameters for the gene (1 row)
process_single_gene_parameters <- function(X, gene_name, debug = FALSE, 
                                         log_func = function(..., iter = NULL) if(debug) message(...), 
                                         iter = NULL) {
  # Initialize default parameters using helper function
  gene_params <- create_default_gene_params(gene_name)
  
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
    
    # Handle variances (can be list or vector)
    covs <- if (is.list(params$variance$sigmasq)) {
      unlist(params$variance$sigmasq)
    } else {
      params$variance$sigmasq
    }
    
    # Store sorted parameters
    sorted_means <- params$mean[sort_idx]
    sorted_variances <- covs[sort_idx]
    sorted_weights <- params$pro[sort_idx]
    
    gene_params$mean0 <- sorted_means[1]
    gene_params$mean1 <- if(length(sorted_means) > 1) sorted_means[2] else NA
    gene_params$variance0 <- sorted_variances[1]
    gene_params$variance1 <- if(length(sorted_variances) > 1) sorted_variances[2] else NA
    gene_params$weight0 <- sorted_weights[1]
    gene_params$weight1 <- if(length(sorted_weights) > 1) sorted_weights[2] else NA
    
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

# ============================================================================
# PARAMETER CACHING FUNCTIONS
# ============================================================================

#' Save GMM parameters to CSV file
#' 
#' @param gmm_params Dataframe with genes as rows and parameters as columns
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
  
  # Convert columns back to proper types using helper function
  return(convert_param_types(params))
}

#' Validate cached parameters against current gene set
#' 
#' @param cached_params Cached parameter dataframe
#' @param current_genes Vector of current gene names
#' @return List with 'valid' (logical), 'missing_genes' (character vector), 'cached_genes' (character vector)
validate_cached_parameters <- function(cached_params, current_genes) {
  if (is.null(cached_params)) {
    return(list(
      valid = FALSE,
      missing_genes = current_genes,
      cached_genes = character(0)
    ))
  }
  
  cached_genes <- intersect(current_genes, cached_params$gene_name)
  missing_genes <- setdiff(current_genes, cached_params$gene_name)
  
  return(list(
    valid = length(missing_genes) == 0,
    missing_genes = missing_genes,
    cached_genes = cached_genes
  ))
}

#' Merge cached parameters with newly computed parameters
#' 
#' @param cached_params Existing cached parameter dataframe (can be NULL)
#' @param new_params Newly computed parameter list or dataframe
#' @param gene_order Vector specifying the desired gene order in final result
#' @return Merged parameter dataframe with genes in specified order
merge_gmm_parameters <- function(cached_params, new_params, gene_order) {
  # Convert new_params list to dataframe if needed
  if (is.list(new_params) && !is.data.frame(new_params)) {
    if (length(new_params) > 0) {
      new_params_df <- do.call(rbind, new_params)
    } else {
      new_params_df <- NULL
    }
  } else {
    new_params_df <- new_params
  }
  
  # Handle case where we only have cached params
  if (is.null(new_params_df)) {
    if (is.null(cached_params)) {
      return(NULL)
    }
    return(cached_params[cached_params$gene_name %in% gene_order, , drop = FALSE])
  }
  
  # Handle case where we only have new params
  if (is.null(cached_params)) {
    return(new_params_df[new_params_df$gene_name %in% gene_order, , drop = FALSE])
  }
  
  # Merge both cached and new parameters
  # New parameters take precedence over cached ones
  all_genes <- union(cached_params$gene_name, new_params_df$gene_name)
  
  # Start with cached parameters
  merged_params <- cached_params
  
  # Update/add new parameters
  for (gene in new_params_df$gene_name) {
    if (gene %in% cached_params$gene_name) {
      # Update existing gene
      merged_params[merged_params$gene_name == gene, ] <- new_params_df[new_params_df$gene_name == gene, ]
    } else {
      # Add new gene
      merged_params <- rbind(merged_params, new_params_df[new_params_df$gene_name == gene, ])
    }
  }
  
  # Return only genes in the specified order
  final_genes <- intersect(gene_order, merged_params$gene_name)
  result <- merged_params[merged_params$gene_name %in% final_genes, , drop = FALSE]
  
  # Reorder to match gene_order
  result <- result[match(final_genes, result$gene_name), , drop = FALSE]
  
  return(result)
}

# ============================================================================
# MAIN PARAMETER EXTRACTION FUNCTIONS
# ============================================================================

#' Extract GMM parameters using sequential processing (helper function)
#' 
#' @param data Matrix or data frame with samples as rows, genes as columns
#' @param debug Logical, whether to enable debug logging
#' @param log_file Path to log file for debug output
#' @return Dataframe containing GMM parameters for each gene
extract_gmm_parameters_sequential <- function(data, debug = FALSE, log_file = NULL) {
  gene_names <- colnames(data)
  n_genes <- length(gene_names)
  
  if (debug) {
    message("DEBUG: Extracting GMM parameters for ", n_genes, " genes using sequential processing")
  }
  
  start_time <- Sys.time()
  gmm_params_list <- list()
  
  for (i in seq_along(gene_names)) {
    gene_name <- gene_names[i]
    X <- data[, gene_name]
    
    if (debug && i %% 100 == 0) {
      message("DEBUG: Processing gene ", i, "/", n_genes, ": ", gene_name)
    }
    
    gene_params <- process_single_gene_parameters(X, gene_name, debug = debug, iter = i)
    gmm_params_list[[i]] <- gene_params
  }
  
  # Combine all gene parameters into single dataframe
  gmm_params <- do.call(rbind, gmm_params_list)
  
  end_time <- Sys.time()
  processing_time <- round(difftime(end_time, start_time, units = "secs"), 1)
  
  if (debug) {
    message("DEBUG: Sequential parameter extraction completed in ", processing_time, " seconds")
    successful_fits <- sum(gmm_params$fit_successful)
    message("DEBUG: Successfully fitted GMM for ", successful_fits, "/", n_genes, " genes")
  }
  
  return(gmm_params)
}

#' Extract GMM parameters using parallel processing (helper function)
#' 
#' @param data Matrix or data frame with samples as rows, genes as columns
#' @param debug Logical, whether to enable debug logging
#' @param num_workers Number of workers to use for parallel processing
#' @return Dataframe containing GMM parameters for each gene, or NULL if parallel processing fails
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
  if (num_workers != 1) {
    available_cores <- detectCores()
    num_cores <- if (num_workers == -1) available_cores else min(num_workers, available_cores)
    
    if (debug) {
      message("Setting up parallel processing with ", num_cores, " cores")
    }
    
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl), add = TRUE)
  } else {
    if (debug) message("Using sequential processing")
    return(NULL)
  }
  
  start_time <- Sys.time()
  if (debug) {
    message("DEBUG: Starting parallel parameter extraction at ", format(start_time, "%m-%d %H:%M:%S"))
  }
  
  # Process genes in parallel
  results_by_gene <- tryCatch({
    if (debug) {
      message("DEBUG: Starting parameter foreach loop for ", n_genes, " genes")
    }
    
    foreach_result <- foreach(
      gene_name = gene_names,
      i = seq_along(gene_names),
      .combine = 'rbind',
      .multicombine = TRUE,
      .errorhandling = 'stop',
      .packages = c('mclust'),
      .export = c('process_single_gene_parameters', 'worker_log_message', 'create_default_gene_params')
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
        
        return(gene_params)
        
      }, error = function(e) {
        # Log worker error
        worker_log_message("ERROR: Worker failed for gene ", gene_name, ": ", e$message, 
                         iter = i, log_file_path = parallel_log_file)
        
        # Create error result using helper function
        error_result <- create_default_gene_params(gene_name, paste("Parallel worker error:", e$message))
        return(error_result)
      })
    }
    
    if (debug) {
      message("DEBUG: Parameter foreach completed, got ", nrow(foreach_result), " results")
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
  processing_time <- round(difftime(end_time, start_time, units = "secs"), 1)
  
  # Handle partial results from parallel processing
  if (is.null(results_by_gene)) {
    if (debug) {
      message("ERROR: Parallel parameter processing returned NULL results")
    }
    return(NULL)
  }
  
  if (debug) {
    message("DEBUG: Parallel parameter extraction completed in ", processing_time, " seconds")
    message("DEBUG: Got ", nrow(results_by_gene), " results for ", n_genes, " genes")
  }
  
  # Check if we got results for all genes
  if (nrow(results_by_gene) != n_genes) {
    if (debug) {
      message("DEBUG: Missing ", n_genes - nrow(results_by_gene), " genes from parallel processing")
    }
    return(NULL)
  }
  
  return(results_by_gene)
}

#' Extract GMM parameters from data with intelligent caching
#' 
#' Fits GMM models once per gene and extracts all parameters needed for
#' different adjustment strategies. Uses cached parameters when available
#' and only computes parameters for missing genes.
#' 
#' @param data Matrix or data frame with samples as rows, genes as columns
#' @param batch_name Name of the batch for caching (optional)
#' @param cache_folder Path to cache directory (optional)
#' @param debug Logical, whether to enable debug logging
#' @param log_file Path to log file for debug output
#' @param num_workers Number of workers to use. If 1, uses sequential processing.
#'                    If -1, uses all available cores. Otherwise uses minimum of specified number and available cores.
#' @return Dataframe containing GMM parameters for each gene (genes as rows, parameters as columns)
extract_gmm_parameters <- function(data, batch_name = NULL, cache_folder = NULL, debug = FALSE, log_file = NULL, num_workers = 1) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input 'data' must be a data frame or matrix.")
  }
  
  # Clear log file if it exists
  if (!is.null(log_file) && file.exists(log_file)) {
    file.remove(log_file)
  }
  
  current_genes <- colnames(data)
  cached_params <- NULL
  
  # Try to load cached parameters if cache info provided
  if (!is.null(batch_name) && !is.null(cache_folder)) {
    cached_params <- load_gmm_parameters_from_csv(batch_name, cache_folder)
    
    if (!is.null(cached_params) && debug) {
      message("DEBUG: Loaded cached parameters for ", nrow(cached_params), " genes")
    }
  }
  
  # Validate cache and identify missing genes
  cache_status <- validate_cached_parameters(cached_params, current_genes)
  
  if (cache_status$valid) {
    if (debug) {
      message("DEBUG: All ", length(current_genes), " genes found in cache")
    }
    return(cached_params[cached_params$gene_name %in% current_genes, , drop = FALSE])
  }
  
  # Determine which genes need computation
  genes_to_compute <- cache_status$missing_genes
  cached_genes <- cache_status$cached_genes
  
  if (debug) {
    message("DEBUG: Found ", length(cached_genes), " cached genes, need to compute ", length(genes_to_compute), " missing genes")
    if (length(genes_to_compute) > 0 && length(genes_to_compute) <= 10) {
      message("DEBUG: Missing genes: ", paste(genes_to_compute, collapse = ", "))
    } else if (length(genes_to_compute) > 10) {
      message("DEBUG: First 10 missing genes: ", paste(head(genes_to_compute, 10), collapse = ", "))
    }
  }
  
  # Extract parameters for missing genes only
  new_params <- NULL
  if (length(genes_to_compute) > 0) {
    data_subset <- data[, genes_to_compute, drop = FALSE]
    
    if (debug) {
      message("DEBUG: Computing parameters for ", length(genes_to_compute), " genes")
    }
    
    # Try parallel processing first if requested, fall back to sequential
    if (num_workers != 1) {
      new_params <- extract_gmm_parameters_parallel(data_subset, debug, log_file, num_workers)
      if (is.null(new_params) && debug) {
        message("DEBUG: Parallel processing failed, falling back to sequential")
      }
    }
    
    # Use sequential processing if parallel failed or not requested
    if (is.null(new_params)) {
      new_params <- extract_gmm_parameters_sequential(data_subset, debug, log_file)
    }
  }
  
  # Merge cached and new parameters
  merged_params <- merge_gmm_parameters(cached_params, new_params, current_genes)
  
  # Save updated parameters to cache if cache info provided
  if (!is.null(batch_name) && !is.null(cache_folder) && !is.null(merged_params)) {
    # Only save if we computed new parameters
    if (length(genes_to_compute) > 0) {
      save_gmm_parameters_to_csv(merged_params, batch_name, cache_folder)
      if (debug) {
        message("DEBUG: Saved updated cache with ", nrow(merged_params), " genes")
      }
    }
  }
  
  return(merged_params)
}