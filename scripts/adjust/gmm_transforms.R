# gmm_transforms.R
# GMM Data Transformation Functions
# This module provides functions for applying GMM-based transformations
# using pre-extracted parameters from gmm_parameters.R

suppressPackageStartupMessages({
  library(mclust)
  library(foreach)
  library(doParallel)
})

# Source parameter extraction functions
source("/scripts/adjust/gmm_parameters.R")

# ============================================================================
# INPUT VALIDATION
# ============================================================================

#' Validate inputs for GMM adjustment functions
#' 
#' @param data Input data matrix/data frame
#' @param batch Optional batch vector
#' @param mixed_strategy Optional mixed strategy parameter
validate_inputs <- function(data, batch = NULL, mixed_strategy = NULL) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input 'data' must be a data frame or matrix.")
  }
  
  if (!is.null(batch) && length(batch) != nrow(data)) {
    stop("Length of 'batch' must equal the number of rows in 'data'.")
  }
  
  if (!is.null(mixed_strategy) && !mixed_strategy %in% c("unimodal", "bimodal")) {
    stop("mixed_strategy must be either 'unimodal' or 'bimodal'.")
  }
}

# ============================================================================
# NPN TRANSFORMATION UTILITIES
# ============================================================================

#' Inverse CDF for GMM distribution
#' 
#' @param p Probability values
#' @param means Component means
#' @param variances Component variances
#' @param weights Component weights
#' @return Quantile values
inverse_cdf_gmm_R <- function(p, means, variances, weights, initial_guess = NULL) {
  # Input validation
  if (any(is.na(means)) || any(is.na(variances)) || any(is.na(weights))) {
    return(rep(NA, length(p)))
  }
  
  if (length(means) != length(variances) || length(means) != length(weights)) {
    return(rep(NA, length(p)))
  }
  
  # For large datasets, use approximation for speed
  if (length(p) > 100000) {
    # Use linear interpolation on a grid for speed
    return(inverse_cdf_gmm_fast(p, means, variances, weights))
  }
  
  gmm_cdf <- function(x) {
    sum(weights * pnorm(x, mean = means, sd = sqrt(pmax(variances, 1e-10))))
  }
  
  stds <- sqrt(pmax(variances, 1e-10))
  search_interval <- c(min(means - 10 * stds), max(means + 10 * stds))
  
  # Expand search interval if needed
  if (abs(search_interval[2] - search_interval[1]) < 1e-6) {
    search_interval <- c(search_interval[1] - 10, search_interval[2] + 10)
  }
  
  solve_for_single_p <- function(p_val, prev_result = NULL) {
    if (is.na(p_val)) return(NA)
    if (p_val <= 0.0) return(-Inf)
    if (p_val >= 1.0) return(Inf)
    
    root_function <- function(x) gmm_cdf(x) - p_val
    
    # Use previous result as starting point if available
    start_point <- if (!is.null(prev_result) && is.finite(prev_result)) {
      prev_result
    } else {
      # Use weighted mean as initial guess
      sum(weights * means)
    }
    
    tryCatch({
      # Try with initial guess first
      if (abs(root_function(start_point)) < 1e-6) {
        return(start_point)
      }
      
      # Adjust search interval around the starting point
      local_interval <- c(
        max(search_interval[1], start_point - 5),
        min(search_interval[2], start_point + 5)
      )
      
      # Check if root is in local interval
      if (root_function(local_interval[1]) * root_function(local_interval[2]) <= 0) {
        result <- uniroot(root_function, interval = local_interval, tol = 1e-4)
        return(result$root)
      } else {
        # Fall back to full interval
        result <- uniroot(root_function, interval = search_interval, tol = 1e-4)
        return(result$root)
      }
    }, error = function(e) {
      # Fallback to normal quantile
      return(qnorm(p_val))
    })
  }
  
  # Process with previous results as hints
  result <- numeric(length(p))
  prev_result <- NULL
  
  for (i in seq_along(p)) {
    result[i] <- solve_for_single_p(p[i], prev_result)
    if (is.finite(result[i])) {
      prev_result <- result[i]
    }
  }
  
  # Replace any remaining problematic values
  if (any(is.infinite(result) | is.na(result))) {
    bad_indices <- which(is.infinite(result) | is.na(result))
    result[bad_indices] <- qnorm(p[bad_indices])
  }
  
  return(result)
}

#' Fast approximation of inverse CDF for large datasets
#' Uses pre-computed grid and linear interpolation
inverse_cdf_gmm_fast <- function(p, means, variances, weights) {
  # Create a grid of x values
  stds <- sqrt(pmax(variances, 1e-10))
  x_min <- min(means - 8 * stds)
  x_max <- max(means + 8 * stds)
  
  # Use adaptive grid size based on data size
  grid_size <- min(10000, max(1000, length(p) / 100))
  x_grid <- seq(x_min, x_max, length.out = grid_size)
  
  # Compute CDF values on grid
  cdf_grid <- sapply(x_grid, function(x) {
    sum(weights * pnorm(x, mean = means, sd = sqrt(pmax(variances, 1e-10))))
  })
  
  # Remove duplicates and ensure monotonicity to avoid "collapsing" warning
  # Find unique CDF values and corresponding x values
  unique_indices <- !duplicated(cdf_grid)
  cdf_unique <- cdf_grid[unique_indices]
  x_unique <- x_grid[unique_indices]
  
  # Ensure strict monotonicity by adding small increments to duplicates
  if (length(cdf_unique) < length(cdf_grid)) {
    # If we lost too many points, create a more refined grid
    if (length(cdf_unique) < grid_size / 2) {
      # Use a different approach: create grid based on quantiles
      prob_grid <- seq(0.001, 0.999, length.out = grid_size)
      x_grid_new <- numeric(length(prob_grid))
      
      for (i in seq_along(prob_grid)) {
        # Simple bisection for each quantile
        x_grid_new[i] <- find_quantile_bisection(prob_grid[i], means, variances, weights, x_min, x_max)
      }
      
      cdf_unique <- prob_grid
      x_unique <- x_grid_new
    }
  }
  
  # Use linear interpolation with unique values
  result <- tryCatch({
    approx(cdf_unique, x_unique, xout = p, rule = 2, ties = "ordered")$y
  }, warning = function(w) {
    # If we still get warnings, fall back to simple normal approximation
    qnorm(p)
  })
  
  # Handle edge cases
  result[p <= 0] <- -Inf
  result[p >= 1] <- Inf
  
  # Fallback for any NAs
  na_indices <- is.na(result)
  if (any(na_indices)) {
    result[na_indices] <- qnorm(p[na_indices])
  }
  
  return(result)
}

#' Simple bisection method to find quantile
find_quantile_bisection <- function(prob, means, variances, weights, x_min, x_max, tol = 1e-4) {
  gmm_cdf <- function(x) {
    sum(weights * pnorm(x, mean = means, sd = sqrt(pmax(variances, 1e-10))))
  }
  
  # Bisection method
  left <- x_min
  right <- x_max
  
  for (i in 1:50) {  # Max 50 iterations
    mid <- (left + right) / 2
    cdf_mid <- gmm_cdf(mid)
    
    if (abs(cdf_mid - prob) < tol) {
      return(mid)
    }
    
    if (cdf_mid < prob) {
      left <- mid
    } else {
      right <- mid
    }
  }
  
  return((left + right) / 2)
}

#' Bimodal non-paranormal transformation
#' 
#' @param X Input data vector
#' @param m0 Mean of first component
#' @param m1 Mean of second component
#' @param v0 Variance of first component
#' @param v1 Variance of second component
#' @param w0 Weight of first component
#' @param w1 Weight of second component
#' @return Transformed data vector
bimodal_npn <- function(X, m0, m1, v0, v1, w0, w1) {
  new_mean1 <- -1
  new_mean2 <- 1
  
  old_mean_diff <- m1 - m0
  new_mean_diff <- new_mean2 - new_mean1
  
  new_sd1 <- sqrt(v0) * new_mean_diff / old_mean_diff
  new_sd2 <- sqrt(v1) * new_mean_diff / old_mean_diff
  
  new_var1 <- new_sd1^2
  new_var2 <- new_sd2^2
  
  ranks <- rank(X, na.last = "keep", ties.method = "average")
  quantiles <- ranks / (sum(!is.na(X)) + 1)
  
  inverse_cdf_gmm_R(
    quantiles,
    means = c(new_mean1, new_mean2),
    variances = c(new_var1, new_var2),
    weights = c(w0, w1)
  )
}

#' Process single gene for unimodal normalization
#' 
#' @param X Gene expression vector
#' @param gene_name Gene identifier
#' @return Normalized gene expression vector
process_unimodal_gene <- function(X, gene_name) {
  if (all(is.na(X))) {
    return(rep(0, length(X)))
  }
  
  if (all(X == X[1])) {
    return(X)
  }
  
  ranks <- rank(X, na.last = "keep", ties.method = "average")
  quantiles <- ranks / (sum(!is.na(X)) + 1)
  unimodal <- qnorm(quantiles)
  return(unimodal / sd(unimodal, na.rm = TRUE))
}

#' Unimodal normalize function with parallel processing
#' 
#' @param data Input data matrix/data frame
#' @param use_parallel Whether to use parallel processing (default: auto-detect based on gene count)
#' @param debug Whether to enable debug logging
#' @param num_workers Number of workers to use. If NULL or 1, uses sequential processing.
#'                    If -1, uses all available cores. Otherwise uses minimum of specified number and available cores.
#' @return Normalized data matrix
unimodal_normalize <- function(data, use_parallel = NULL, debug = FALSE, num_workers = NULL) {
  validate_inputs(data)
  
  gene_names <- colnames(data)
  n_samples <- nrow(data)
  n_genes <- ncol(data)
  
  # Handle num_workers parameter
  if (!is.null(num_workers) && num_workers == 1) {
    use_parallel <- FALSE
  } else if (!is.null(num_workers)) {
    use_parallel <- TRUE
  }
  
  # Auto-decide on parallelization
  if (is.null(use_parallel)) {
    use_parallel <- n_genes > 50  # Use parallel for datasets with >50 genes
  }
  
  if (use_parallel && n_genes > 1) {
    if (debug) {
      message("DEBUG: Using parallel processing for unimodal normalization of ", n_genes, " genes")
    }
    
    # Setup parallel processing
    available_cores <- detectCores()
    num_cores <- if (num_workers == -1) available_cores else min(num_workers, available_cores)
    
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl), add = TRUE)
    
    # Export function to workers
    clusterExport(cl, c("process_unimodal_gene"), envir = .GlobalEnv)
      
    # Process genes in parallel
    results_list <- tryCatch({
      foreach(
        gene_name = gene_names,
        .combine = 'cbind',
        .multicombine = TRUE,
        .errorhandling = 'stop'
      ) %dopar% {
        X <- data[, gene_name]
        process_unimodal_gene(X, gene_name)
      }
    }, error = function(e) {
      if (debug) {
        message("Parallel unimodal processing failed: ", e$message)
      }
      return(NULL)
    })
    
    if (!is.null(results_list)) {
      colnames(results_list) <- gene_names
      rownames(results_list) <- rownames(data)
      return(results_list)
    }
  }
  
  # Sequential processing (fallback or requested)
  if (debug) {
    message("DEBUG: Using sequential processing for unimodal normalization")
  }
  
  unimodal_data <- matrix(NA, nrow = n_samples, ncol = n_genes)
  colnames(unimodal_data) <- gene_names
  rownames(unimodal_data) <- rownames(data)
  
  for (i in seq_along(gene_names)) {
    X <- data[, gene_names[i]]
    unimodal_data[, i] <- process_unimodal_gene(X, gene_names[i])
  }
  
  return(unimodal_data)
}

# ==========================================================================
# ADJUSTMENT STRATEGY IMPLEMENTATIONS
# ============================================================================

#' Apply specific adjustment strategy using cached parameters
#' 
#' @param X Input data vector
#' @param gene_params GMM parameters for the gene
#' @param strategy Adjustment strategy ("simple", "npn")
#' @param debug Whether to enable debug logging
#' @param log_file Path to log file
#' @return Transformed data vector
apply_adjustment_strategy <- function(X, gene_params, strategy = "simple", debug = FALSE, log_file = NULL) {
  fit_successful <- as.logical(gene_params$fit_successful)
  n_components <- as.numeric(gene_params$n_components)
  if (is.na(fit_successful) || !fit_successful || n_components < 2) {
    # Fallback to unimodal
    ranks <- rank(X, na.last = "keep", ties.method = "average")
    quantiles <- ranks / (sum(!is.na(X)) + 1)
    unimodal <- qnorm(quantiles)
    return(unimodal / sd(unimodal, na.rm = TRUE))
  }
  
  # Extract parameters
  m0 <- as.numeric(gene_params$mean0)
  m1 <- as.numeric(gene_params$mean1)
  v0 <- as.numeric(gene_params$variance0)
  v1 <- as.numeric(gene_params$variance1)
  w0 <- as.numeric(gene_params$weight0)
  w1 <- as.numeric(gene_params$weight1)
  
  # Apply log-shift transformation (same as in parameter extraction)
  min_val <- min(X, na.rm = TRUE)
  X_transformed <- log(X - min_val + 1)
  
  if (strategy == "simple") {
    # Simple bimodal transformation
    ranks <- rank(X_transformed, na.last = "keep", ties.method = "average")
    quantiles <- ranks / (sum(!is.na(X_transformed)) + 1)
    
    # Transform to bimodal normal with means at -1 and 1
    new_mean1 <- -1
    new_mean2 <- 1
    
    old_mean_diff <- m1 - m0
    new_mean_diff <- new_mean2 - new_mean1
    
    new_sd1 <- sqrt(v0) * new_mean_diff / old_mean_diff
    new_sd2 <- sqrt(v1) * new_mean_diff / old_mean_diff
    
    new_var1 <- new_sd1^2
    new_var2 <- new_sd2^2
    
    result <- inverse_cdf_gmm_R(
      quantiles,
      means = c(new_mean1, new_mean2),
      variances = c(new_var1, new_var2),
      weights = c(w0, w1)
    )
    
    return(result)
    
  } else if (strategy == "npn") {
    # Non-paranormal transformation
    return(bimodal_npn(X_transformed, m0, m1, v0, v1, w0, w1))
    
  } else if (strategy == "npn_unit_std") {
    adjusted = bimodal_npn(X_transformed, m0, m1, v0, v1, w0, w1)
    new_mean = mean(adjusted, na.rm = TRUE)
    adjusted = adjusted - new_mean
    adjusted = adjusted / sd(adjusted, na.rm = TRUE)
    return(adjusted + new_mean)
  }
  
  else {
    stop("Unknown adjustment strategy: ", strategy)
  }
}

#' Process single gene for GMM adjustment
#' 
#' @param X Gene expression vector
#' @param gene_name Gene identifier
#' @param gene_params GMM parameters for the gene
#' @param strategy Adjustment strategy
#' @param debug Whether to enable debug logging
#' @param log_file Path to log file
#' @param iter Iteration number for logging
#' @return Adjusted gene expression vector
process_gmm_adjustment_gene <- function(X, gene_name, gene_params, strategy, debug, log_file, iter) {
  tryCatch({
    if (is.null(gene_params) || !gene_params$fit_successful) {
      # Unimodal fallback
      worker_log_message("No valid cached parameters for gene", gene_name, "- using unimodal fallback", iter=iter, log_file_path=log_file)
      ranks <- rank(X, na.last = "keep", ties.method = "average")
      quantiles <- ranks / (sum(!is.na(X)) + 1)
      unimodal <- qnorm(quantiles)
      return(unimodal / sd(unimodal, na.rm = TRUE))
    }
    
    # Edge case checks
    if (all(is.na(X))) {
      return(rep(0, length(X)))
    }
    
    if (all(X == X[1])) {
      return(X)
    }
    
    # Use cached strategy functions for fast adjustment
    if (gene_params$n_components >= 2 && gene_params$recommended_modes == 2) {
      result <- apply_adjustment_strategy(X, gene_params, strategy, debug, log_file)
      worker_log_message("Applied", strategy, "adjustment for gene", gene_name, iter=iter, log_file_path=log_file)
      return(result)
    } else {
      # Single component - unimodal fallback
      ranks <- rank(X, na.last = "keep", ties.method = "average")
      quantiles <- ranks / (sum(!is.na(X)) + 1)
      unimodal <- qnorm(quantiles)
      worker_log_message("Single component for gene", gene_name, "- used unimodal adjustment", iter=iter, log_file_path=log_file)
      return(unimodal / sd(unimodal, na.rm = TRUE))
    }
    
  }, error = function(e) {
    worker_log_message("Error applying", strategy, "adjustment for gene", gene_name, ":", e$message, iter=iter, log_file_path=log_file)
    # Unimodal fallback
    ranks <- rank(X, na.last = "keep", ties.method = "average")
    quantiles <- ranks / (sum(!is.na(X)) + 1)
    unimodal <- qnorm(quantiles)
    return(unimodal / sd(unimodal, na.rm = TRUE))
  })
}

#' Apply GMM adjustment using cached parameters with parallel processing
#' 
#' @param data Input data matrix/data frame
#' @param gmm_params GMM parameters from extract_gmm_parameters
#' @param strategy Adjustment strategy
#' @param debug Whether to enable debug logging
#' @param log_file Path to log file
#' @param use_parallel Whether to use parallel processing (default: auto-detect based on gene count)
#' @param num_workers Number of workers to use. If NULL or 1, uses sequential processing.
#'                    If -1, uses all available cores. Otherwise uses minimum of specified number and available cores.
#' @return Adjusted data matrix
apply_gmm_adjustment <- function(data, gmm_params, strategy = "simple", debug = FALSE, log_file = NULL, num_workers = NULL) {
  if (debug) {
    message("Starting apply_gmm_adjustment with strategy:", strategy)
  }
  
  validate_inputs(data)
  
  if (is.null(gmm_params) || !is.list(gmm_params)) {
    stop("gmm_params must be a valid list of GMM parameters")
  }
  
  gene_names <- colnames(data)
  n_samples <- nrow(data)
  n_genes <- ncol(data)
  
  # Handle num_workers parameter
  if (is.null(num_workers) || num_workers == 1) {
    use_parallel <- FALSE
  } else {
    use_parallel <- TRUE
  }
  
  # Auto-decide on parallelization
  if (is.null(use_parallel)) {
    use_parallel <- n_genes > 50  # Use parallel for datasets with >50 genes
  }
  
  if (use_parallel && n_genes > 1) {  
    # Setup parallel processing
    available_cores <- detectCores()
    num_cores <- if (num_workers == -1) available_cores else min(num_workers, available_cores)
    
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl), add = TRUE)
    
    # Export functions to workers
    clusterExport(cl, c("process_gmm_adjustment_gene", "apply_adjustment_strategy", 
                       "inverse_cdf_gmm_R", "bimodal_npn", "worker_log_message"), 
                 envir = .GlobalEnv)

      if (debug) {
        message("Starting parallel gmm adjustment")
      }
      
      # Process genes in parallel
      results_list <- tryCatch({
        foreach(
          gene_name = gene_names,
          i = seq_along(gene_names),
          .combine = 'cbind',
          .multicombine = TRUE,
          .errorhandling = 'stop'
        ) %dopar% {
          X <- data[, gene_name]
          gene_params <- gmm_params[[gene_name]]
          process_gmm_adjustment_gene(X, gene_name, gene_params, strategy, debug, log_file, i)
        }
      }, error = function(e) {
        if (debug) {
          message("Parallel GMM adjustment failed: ", e$message)
        }
        return(NULL)
      })
      
      if (!is.null(results_list)) {
        colnames(results_list) <- gene_names
        rownames(results_list) <- rownames(data)
        return(results_list)
      }
  }
  
  # Sequential processing (fallback or requested)
  if (debug) {
    message("DEBUG: Using sequential processing for GMM adjustment")
  }
  
  adjusted_data <- matrix(NA, nrow = n_samples, ncol = n_genes)
  colnames(adjusted_data) <- gene_names
  rownames(adjusted_data) <- rownames(data)
  
  for (i in seq_along(gene_names)) {
    gene_name <- gene_names[i]
    X <- data[, gene_name]
    gene_params <- gmm_params[[gene_name]]
    adjusted_data[, i] <- process_gmm_adjustment_gene(X, gene_name, gene_params, strategy, debug, log_file, i)
  }
  
  if (debug) {
    message("Completed apply_gmm_adjustment for ", n_genes, " genes")
  }
  return(adjusted_data)
}

# ============================================================================
# BATCH PROCESSING UTILITIES
# ============================================================================

#' Adjust data based on recommended modes and mixed strategy
#' 
#' @param unimodal_adjusted Unimodal adjusted data
#' @param bimodal_adjusted Bimodal adjusted data
#' @param batch_numeric Numeric batch vector
#' @param recommended_modes_df Data frame of recommended modes per batch
#' @param debug Whether to enable debug logging
#' @param log_file Path to log file
#' @param mixed_strategy Strategy for mixed scenarios
#' @return Final adjusted data matrix
gmm_batch_adjust <- function(unimodal_adjusted, bimodal_adjusted, batch_numeric, 
                           recommended_modes_df, debug, log_file, mixed_strategy) {
  
  n_samples <- nrow(unimodal_adjusted)
  n_genes <- ncol(unimodal_adjusted)
  adjusted_data <- matrix(NA, nrow = n_samples, ncol = n_genes)
  colnames(adjusted_data) <- colnames(unimodal_adjusted)
  rownames(adjusted_data) <- rownames(unimodal_adjusted)
  
  batch_levels <- rownames(recommended_modes_df)
  
  for (j in 1:n_genes) {
    gene_name <- colnames(unimodal_adjusted)[j]
    
    # Count modes across batches for this gene
    modes_for_gene <- recommended_modes_df[, j]
    unimodal_batches <- sum(modes_for_gene == 1, na.rm = TRUE)
    bimodal_batches <- sum(modes_for_gene == 2, na.rm = TRUE)
    
    # Determine strategy
    if (unimodal_batches > 0 && bimodal_batches > 0) {
      # Mixed case
      if (mixed_strategy == "unimodal") {
        adjusted_data[, j] <- unimodal_adjusted[, j]
        worker_log_message("Gene", gene_name, "- mixed modes, using unimodal strategy", log_file_path = log_file)
      } else {
        adjusted_data[, j] <- bimodal_adjusted[, j]
        worker_log_message("Gene", gene_name, "- mixed modes, using bimodal strategy", log_file_path = log_file)
      }
    } else if (bimodal_batches > 0) {
      # All bimodal
      adjusted_data[, j] <- bimodal_adjusted[, j]
      worker_log_message("Gene", gene_name, "- all batches bimodal", log_file_path = log_file)
    } else {
      # All unimodal or no valid modes
      adjusted_data[, j] <- unimodal_adjusted[, j]
      worker_log_message("Gene", gene_name, "- all batches unimodal", log_file_path = log_file)
    }
  }
  
  return(adjusted_data)
}