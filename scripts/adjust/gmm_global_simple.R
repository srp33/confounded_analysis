# gmm_global_simple.R
# Simple Gene-Global GMM Adjustment Strategy
# 
# This implements a very simple global adjustment approach:
# 1. Flatten all expression values to create a global distribution
# 2. Fit a two-component GMM to this global distribution in log space
# 3. Find the midpoint between the two modes
# 4. Apply transformation strategy using global GMM parameters

suppressPackageStartupMessages({
  library(mclust)
})

# Source functions from gmm_transforms.R for NPN functionality
source("/scripts/adjust/gmm_transforms.R")

#' Core function for gene-global GMM adjustment with multiple strategies
#' 
#' This function implements a global adjustment strategy that:
#' 1. Transforms all data to log space
#' 2. Flattens all expression values to create a global distribution
#' 3. Fits a two-component GMM to the log-transformed global distribution
#' 4. Applies the specified transformation strategy using global GMM parameters
#' 
#' @param data Input data matrix/data frame (samples x genes)
#' @param strategy Adjustment strategy: "simple" (centering/scaling) or "npn" (bimodal NPN)
#' @param debug Whether to enable debug logging
#' @return List with adjusted_data and adjustment_info
gmm_global_core <- function(data, strategy = "simple", debug = FALSE) {
  
  if (debug) {
    message("Starting gene-global GMM adjustment with strategy: ", strategy)
    message("Input data dimensions: ", nrow(data), " samples x ", ncol(data), " genes")
  }
  
  # Step 1: Transform all data to log space
  min_val <- min(data, na.rm = TRUE)
  data_log <- log(data - min_val + 1)
  
  # Step 2: Flatten all expression values to create a global distribution (in log space)
  global_distribution_log <- as.vector(data_log)
  
  if (debug) {
    message("Global distribution (log space) summary:")
    print(summary(global_distribution_log))
  }
  
  # Helper function for fallback values
  set_fallback_values <- function(reason = "GMM fitting failed") {
    if (debug) {
      message(reason, ", using median as adjustment")
    }
    return(list(
      adjustment_value_log = median(global_distribution_log, na.rm = TRUE),
      fit_successful = FALSE,
      component_info = list(scale_factor = 1.0, min_val = min_val)
    ))
  }
  
  # Step 3: Fit a two-component GMM to the global distribution (in log space)
  result <- tryCatch({
    if (debug) {
      message("Fitting GMM to all expression values in log space...")
    }
    
    # For very large datasets, subsample for GMM fitting to speed up
    if (length(global_distribution_log) > 1000000) {
      if (debug) {
        message("Large dataset detected, using subsample for GMM fitting...")
      }
      sample_size <- min(500000, length(global_distribution_log))
      sample_indices <- sample(length(global_distribution_log), sample_size)
      gmm_data <- global_distribution_log[sample_indices]
    } else {
      gmm_data <- global_distribution_log
    }
    
    # Fit GMM with exactly 2 components
    gmm_fit <- Mclust(gmm_data, G = 2, modelNames = "V", verbose = FALSE)
    
    if (is.null(gmm_fit)) {
      return(set_fallback_values("GMM fitting returned NULL"))
    }
    
    # Step 4: Find the midpoint between the two modes and calculate scaling
    means <- gmm_fit$parameters$mean
    weights <- gmm_fit$parameters$pro
    
    # Sort by mean to get consistent ordering
    sort_idx <- order(means)
    mean1 <- means[sort_idx[1]]
    mean2 <- means[sort_idx[2]]
    weight1 <- weights[sort_idx[1]]
    weight2 <- weights[sort_idx[2]]
    
    # Midpoint between the two modes (in log space)
    midpoint_log <- (mean1 + mean2) / 2
    
    # Calculate scaling factor to make modes end up at -1 and +1
    # Distance between modes in log space
    mode_distance_log <- abs(mean2 - mean1)
    # We want this distance to be 2 (from -1 to +1)
    scale_factor <- 2.0 / mode_distance_log
    
    component_info <- list(
      mean1 = mean1,
      mean2 = mean2,
      weight1 = weight1,
      weight2 = weight2,
      midpoint_log = midpoint_log,
      mode_distance_log = mode_distance_log,
      scale_factor = scale_factor,
      min_val = min_val
    )
    
    if (debug) {
      message("GMM fit successful:")
      message("  Component 1 - Mean: ", round(mean1, 3), ", Weight: ", round(weight1, 3))
      message("  Component 2 - Mean: ", round(mean2, 3), ", Weight: ", round(weight2, 3))
      message("  Midpoint (log): ", round(midpoint_log, 3))
      message("  Mode distance (log): ", round(mode_distance_log, 3))
      message("  Scale factor: ", round(scale_factor, 3))
    }
    
    list(
      adjustment_value_log = midpoint_log,
      fit_successful = TRUE,
      component_info = component_info
    )
    
  }, error = function(e) {
    return(set_fallback_values(paste("Error in GMM fitting:", e$message)))
  })
  
  # Extract results
  adjustment_value_log <- result$adjustment_value_log
  fit_successful <- result$fit_successful
  component_info <- result$component_info
  
  # Step 5: Apply the specified transformation strategy
  # Initialize adjusted_data with fallback
  adjusted_data <- data_log
  
  if (strategy == "simple") {
    # Simple centering and scaling in log space
    centered_data_log <- data_log - adjustment_value_log
    
    if (fit_successful && !is.null(component_info$scale_factor)) {
      adjusted_data_log <- centered_data_log * component_info$scale_factor
      if (debug) {
        message("Applied simple centering and scaling in log space")
        message("  Adjustment value (log): ", round(adjustment_value_log, 3))
        message("  Scale factor: ", round(component_info$scale_factor, 3))
      }
    } else {
      adjusted_data_log <- centered_data_log
      if (debug) {
        message("Applied centering only - no scaling due to GMM failure")
      }
    }
    
    adjusted_data <- adjusted_data_log
    
  } else if (strategy == "npn") {
    # Bimodal NPN transformation using global GMM parameters
    if (fit_successful) {
      if (debug) {
        message("Applying bimodal NPN transformation using global GMM parameters")
      }
      
      # Use the already flattened data for ranking (same as global_distribution_log)
      all_values <- global_distribution_log
      
      # Create target bimodal distribution with modes at -1 and +1
      target_means <- c(-1, 1)
      target_variances <- c(component_info$weight1, component_info$weight2) * 0.5  # Scale variances
      target_weights <- c(component_info$weight1, component_info$weight2)
      
      # Rank all values globally
      ranks <- rank(all_values, na.last = "keep", ties.method = "average")
      quantiles <- ranks / (sum(!is.na(all_values)) + 1)
      
      # Transform through inverse CDF of target bimodal distribution
      transformed_values <- tryCatch({
        inverse_cdf_gmm_R(
          quantiles,
          means = target_means,
          variances = target_variances,
          weights = target_weights
        )
      }, error = function(e) {
        if (debug) {
          message("Error in inverse_cdf_gmm_R: ", e$message)
        }
        return(NULL)
      })
      
      
      transformed_values <- qnorm(quantiles)
      
      # Reshape back to original matrix dimensions
      adjusted_data <- matrix(transformed_values, nrow = nrow(data_log), ncol = ncol(data_log))
      
    } else {
      # Fallback to simple unimodal NPN if GMM failed
      if (debug) {
        message("GMM failed, applying unimodal NPN transformation")
      }
      
      # Use the already flattened data (same as global_distribution_log)
      all_values <- global_distribution_log
      ranks <- rank(all_values, na.last = "keep", ties.method = "average")
      quantiles <- ranks / (sum(!is.na(all_values)) + 1)
      transformed_values <- qnorm(quantiles)
      
      # Reshape back and standardize
      adjusted_data <- matrix(transformed_values, nrow = nrow(data_log), ncol = ncol(data_log))
      adjusted_data <- adjusted_data / sd(adjusted_data, na.rm = TRUE)
    }
    
  } else {
    stop("Unknown strategy: ", strategy, ". Supported strategies are 'simple' and 'npn'.")
  }
  
  if (debug) {
    message("Adjustment completed successfully with strategy: ", strategy)
    if (!is.null(adjusted_data)) {
      message("Final data range: [", round(min(adjusted_data, na.rm = TRUE), 3), ", ", round(max(adjusted_data, na.rm = TRUE), 3), "]")
    }
  }
  
  return(list(
    adjusted_data = adjusted_data,
    adjustment_info = list(
      strategy = strategy,
      adjustment_value_log = adjustment_value_log,
      fit_successful = fit_successful,
      component_info = component_info,
      global_distribution_summary = summary(global_distribution_log),
      min_val_original = min_val
    )
  ))
}

#' Simple gene-global GMM adjustment (backward compatibility)
#' 
#' @param data Input data matrix/data frame (samples x genes)
#' @param debug Whether to enable debug logging
#' @return List with adjusted_data and adjustment_info
gmm_global_simple_core <- function(data, debug = FALSE) {
  return(gmm_global_core(data, strategy = "simple", debug = debug))
}

#' Bimodal NPN gene-global GMM adjustment
#' 
#' @param data Input data matrix/data frame (samples x genes)
#' @param debug Whether to enable debug logging
#' @return List with adjusted_data and adjustment_info
gmm_global_npn_core <- function(data, debug = FALSE) {
  return(gmm_global_core(data, strategy = "npn", debug = debug))
}