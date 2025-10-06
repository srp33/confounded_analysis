# gmm_adjust.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# ============================================================================
# GMM IMPLEMENTATION
# ============================================================================

#' 2-component Gaussian Mixture Model
GaussianMixture2D <- function(max_iter = 100, tol = 1e-4, alpha0 = 10.0) {
  structure(list(
    max_iter = max_iter,
    tol = tol,
    alpha0 = alpha0,
    means_ = NULL,
    variances_ = NULL,
    weights_ = NULL,
    resp_ = NULL
  ), class = "GaussianMixture2D")
}

#' Normal PDF calculation
normal_pdf <- function(x, mean, sd) {
  exp(-0.5 * ((x - mean) / sd)^2) / (sd * sqrt(2 * pi))
}

#' Fit 2-component GMM model with prior on means and weights
fit.GaussianMixture2D <- function(model, X, prior_alpha = 0.5) {
  X <- as.vector(X)
  n <- length(X)
  K <- 2  # Always 2 components
  eps <- 1e-12

  # Initialize means using percentiles (25th and 75th percentiles)
  percentiles <- c(25, 75)
  means <- quantile(X, percentiles / 100, na.rm = TRUE)

  variances <- rep(var(X, na.rm = TRUE), K)
  weights <- rep(0.5, K)  # Equal weights initially

  log_likelihood_old <- -Inf

  # Prior centers c_k from percentiles
  centers <- means
  # Prior variance s_k^2 = prior_alpha * data variance
  data_var <- var(X, na.rm = TRUE)
  prior_var <- prior_alpha * data_var + eps

  for (iter in 1:model$max_iter) {
    # E-step: responsibilities
    pdfs <- matrix(0, nrow = n, ncol = K)
    for (k in 1:K) {
      pdfs[, k] <- weights[k] * normal_pdf(X, means[k], sqrt(variances[k]))
    }
    responsibilities <- pdfs / (rowSums(pdfs) + eps)

    # M-step: update parameters
    Nk <- colSums(responsibilities)

    # Update weights with Dirichlet prior (MAP estimate)
    weights <- (Nk + model$alpha0 - 1) / (n + K * (model$alpha0 - 1))

    # Update means with prior (MAP)
    for (k in 1:K) {
      shrink_factor <- variances[k] / prior_var
      weighted_sum <- sum(responsibilities[, k] * X)
      means[k] <- (weighted_sum + shrink_factor * centers[k]) / (Nk[k] + shrink_factor + eps)
    }

    # Update variances
    for (k in 1:K) {
      variances[k] <- sum(responsibilities[, k] * (X - means[k])^2) / (Nk[k] + eps)
    }
    variances <- pmax(variances, 1e-6)  # Variance floor

    # Check convergence
    log_likelihood <- sum(log(rowSums(pdfs) + eps))
    if (abs(log_likelihood - log_likelihood_old) < model$tol) {
      break
    }
    log_likelihood_old <- log_likelihood
  }

  model$means_ <- means
  model$variances_ <- variances
  model$weights_ <- weights
  model$resp_ <- responsibilities

  return(model)
}


#' Predict probabilities
predict_proba.GaussianMixture2D <- function(model, X) {
  X <- as.vector(X)
  n <- length(X)
  K <- 2
  
  pdfs <- matrix(0, nrow = n, ncol = K)
  for (k in 1:K) {
    pdfs[, k] <- model$weights_[k] * normal_pdf(X, model$means_[k], sqrt(model$variances_[k]))
  }
  
  return(pdfs / (rowSums(pdfs) + 1e-12))
}

# ============================================================================
# INVERSE CDF TRANSFORMATION
# ============================================================================

#' Inverse CDF for GMM distribution with robust error handling
inverse_cdf_gmm <- function(p, means, variances, weights) {
  # Input validation
  if (any(is.na(means)) || any(is.na(variances)) || any(is.na(weights))) {
    return(rep(NA, length(p)))
  }
  
  # Check for degenerate cases
  if (length(unique(means)) == 1 || any(variances < 1e-12)) {
    # Fallback to normal quantiles if parameters are degenerate
    return(qnorm(p))
  }
  
  gmm_cdf <- function(x) {
    sum(weights * pnorm(x, mean = means, sd = sqrt(pmax(variances, 1e-10))))
  }
  
  # More robust search interval calculation
  stds <- sqrt(pmax(variances, 1e-10))
  min_bound <- min(means - 15 * stds)
  max_bound <- max(means + 15 * stds)
  
  # Ensure minimum interval width
  interval_width <- max_bound - min_bound
  if (interval_width < 20) {
    center <- (min_bound + max_bound) / 2
    min_bound <- center - 10
    max_bound <- center + 10
  }
  
  solve_for_single_p <- function(p_val) {
    if (is.na(p_val)) return(NA)
    if (p_val <= 1e-10) return(min_bound - 5)
    if (p_val >= 1 - 1e-10) return(max_bound + 5)
    
    root_function <- function(x) gmm_cdf(x) - p_val
    
    # Try with initial search interval
    search_interval <- c(min_bound, max_bound)
    
    # Check if root exists in interval by evaluating endpoints
    f_left <- root_function(search_interval[1])
    f_right <- root_function(search_interval[2])
    
    # If same sign, try expanding the interval
    expansion_attempts <- 0
    max_expansions <- 3
    
    while (sign(f_left) == sign(f_right) && expansion_attempts < max_expansions) {
      expansion_factor <- 2^(expansion_attempts + 1)
      search_interval[1] <- min_bound - 10 * expansion_factor
      search_interval[2] <- max_bound + 10 * expansion_factor
      
      f_left <- root_function(search_interval[1])
      f_right <- root_function(search_interval[2])
      expansion_attempts <- expansion_attempts + 1
    }
    
    # Try uniroot with error handling
    tryCatch({
      if (sign(f_left) != sign(f_right)) {
        return(uniroot(root_function, interval = search_interval, tol = 1e-4)$root)
      } else {
        # If still same sign, use normal quantile as fallback
        return(qnorm(p_val))
      }
    }, error = function(e) {
      # Fallback to normal quantile on any error
      return(qnorm(p_val))
    })
  }
  
  result <- sapply(p, solve_for_single_p)
  
  # Final cleanup for any remaining problematic values
  if (any(is.infinite(result) | is.na(result))) {
    bad_indices <- which(is.infinite(result) | is.na(result))
    result[bad_indices] <- qnorm(p[bad_indices])
  }
  
  return(result)
}

# ============================================================================
# ADJUSTMENT FUNCTIONS
# ============================================================================

get_gene_gmm_transform <- function(
    gene_exp,
    alpha0 = 10,
    nonlinear = TRUE,
    mean_mean_zero = TRUE,
    unit_var = TRUE,
    means_at_1 = FALSE
) {
  if (all(is.na(gene_exp))) return(gene_exp)

  # --- Log-transform ---
  min_val <- min(gene_exp, na.rm = TRUE)
  x_transformed <- log(gene_exp - min_val + 1)

  mean_shift_fallback <- scale(x_transformed, scale = FALSE)[, 1]

  # Check for small variance
  if (var(x_transformed, na.rm = TRUE) < 1e-8) return(mean_shift_fallback)

  # --- Quantiles ---
  ranks <- rank(x_transformed, na.last = "keep", ties.method = "average")
  n_valid <- sum(!is.na(x_transformed))
  quantiles <- ranks / (n_valid + 1)

  # --- Fit 2-component GMM ---
  gmm <- GaussianMixture2D(alpha0 = alpha0)
  gmm <- fit(gmm, x_transformed)

  qnorm_fallback <- qnorm(quantiles)

  # Validate GMM parameters
  if (any(is.na(gmm$means_)) || any(is.na(gmm$variances_)) || any(is.na(gmm$weights_))) {
    return(qnorm_fallback)
  }

  # Ensure lower component first
  if (gmm$means_[1] > gmm$means_[2]) {
    means <- c(gmm$means_[2], gmm$means_[1])
    variances <- c(gmm$variances_[2], gmm$variances_[1])
    weights <- c(gmm$weights_[2], gmm$weights_[1])
  } else {
    means <- gmm$means_
    variances <- gmm$variances_
    weights <- gmm$weights_
  }

  # Check for small variances, likely not truly bimodal
  if (any(variances < 1e-9)) return(qnorm_fallback)

  # --- Nonlinear mapping to original GMM distribution ---
  if (nonlinear) {
    mapped <- tryCatch({
      inverse_cdf_gmm(
        quantiles,
        means = means,
        variances = variances,
        weights = weights
      )
    }, error = function(e) {
      NA
    })

    if (any(is.na(mapped) | is.infinite(mapped))) {
      x_transformed <- mean_shift_fallback
    } else {
      x_transformed <- mapped
    }
  }

  # --- Affine corrections ---
  if (mean_mean_zero) {
    mean_center <- 0.5 * (means[1] + means[2])
    x_transformed <- x_transformed - mean_center
  }

  if (unit_var) {
    variance <- 0.5 * (variances[1] + variances[2]) + 0.25 * (means[2] - means[1])^2
    scale_factor <- sqrt(max(variance, 1e-9))
    x_transformed <- x_transformed / scale_factor
  }

  if (means_at_1) {
    if (unit_var) stop("Cannot have both means_at_1 and unit_var")
    if (!mean_mean_zero) stop("Cannot have means_at_1 without mean_mean_zero")

    # Compute safe scaling so component means map to ±1
    scale_factor <- (means[2] - means[1]) / 2
    if (abs(scale_factor) < 1e-6) {
      return(mean_shift_fallback)
    }
    x_transformed <- x_transformed / scale_factor
  }

  return(x_transformed)
}


#' Bimodal normalization using GMM with full inverse CDF
bimodal_normalize <- function(data, alpha0 = 10, mean_only = FALSE, debug = FALSE) {
  gene_names <- colnames(data)
  n_genes <- length(gene_names)
  n_samples <- nrow(data)
  
  # Initialize results
  bimodal_data <- matrix(NA, nrow = n_samples, ncol = n_genes)
  
  for (i in seq_along(gene_names)) {
    gene_name <- gene_names[i]
    gene_exp <- data[, i]
    
    if (debug && i %% 1000 == 0) {
      cat("Processed", i, "/", n_genes, "genes\n")
    }
    
    # Skip if all NA
    if (all(is.na(gene_exp))) {
      bimodal_data[, i] <- gene_exp
      next
    }
    
    # Skip if all identical values
    if (all(gene_exp == gene_exp[1], na.rm = TRUE)) {
      bimodal_data[, i] <- gene_exp
      next
    }
    
    tryCatch({
      # Apply GMM transformation (with or without inverse CDF)
      bimodal_result <- get_gene_gmm_transform(gene_exp, alpha0, mean_only)
      
      # Validate result
      if (all(is.na(bimodal_result)) || all(is.infinite(bimodal_result))) {
        stop("Invalid transformation result")
      }
      
      bimodal_data[, i] <- bimodal_result
      
    }, error = function(e) {
      if (debug) {
        cat("Error processing gene", gene_name, ":", e$message, "\n")
      }
      # Fallback to simple quantile normalization
      tryCatch({
        min_val <- min(gene_exp, na.rm = TRUE)
        x_transformed <- log(gene_exp - min_val + 1)
        n_valid <- sum(!is.na(x_transformed))
        
        if (n_valid > 1) {
          quantiles <- rank(x_transformed, na.last = "keep", ties.method = "average") / (n_valid + 1)
          fallback_result <- qnorm(quantiles)
          
          # Standardize if possible
          result_sd <- sd(fallback_result, na.rm = TRUE)
          if (!is.na(result_sd) && result_sd > 1e-10) {
            bimodal_data[, i] <- fallback_result / result_sd
          } else {
            bimodal_data[, i] <- fallback_result
          }
        } else {
          # If only one valid value or all NA, keep original
          bimodal_data[, i] <- gene_exp
        }
      }, error = function(e2) {
        # Ultimate fallback - keep original data
        bimodal_data[, i] <- gene_exp
      })
    })
  }
  
  colnames(bimodal_data) <- gene_names
  rownames(bimodal_data) <- rownames(data)
  
  return(bimodal_data)
}

#' GMM adjustment for multiple batches
#' 
#' @param data Matrix of gene expression data (samples x genes)
#' @param batch Vector of batch labels for each sample
#' @param alpha0 Dirichlet prior parameter for GMM weights
#' @param mean_only If TRUE, only adjust means without using inverse CDF transformation
#' @param debug If TRUE, print progress messages
#' 
#' @details When mean_only=TRUE, the function performs a simpler adjustment that
#' only shifts the means of the two GMM components to -1 and 1, without applying
#' the full inverse CDF transformation. This preserves more of the original data
#' structure while still achieving bimodal normalization.
gmm_adjust <- function(data, batch, alpha0 = 10, mean_only = FALSE, debug = FALSE) {
  if (debug) {
    cat("Starting GMM adjustment...\n")
    cat("Input data dimensions:", nrow(data), "rows,", ncol(data), "columns\n")
    cat("Input batch dimensions:", length(batch), "elements\n")
  }
  
  batch_factor <- as.factor(batch)
  batch_levels <- levels(batch_factor)
  
  n_samples <- nrow(data)
  n_genes <- ncol(data)
  
  adjusted_data <- matrix(NA, nrow = n_samples, ncol = n_genes)
  
  for (b in batch_levels) {
    if (debug) {
      cat("Processing batch:", b, "\n")
    }
    
    batch_indices <- which(batch == b)
    batch_data <- data[batch_indices, , drop = FALSE]
    
    # Apply bimodal normalization to all genes
    batch_adjusted <- bimodal_normalize(batch_data, alpha0, mean_only, debug)
    adjusted_data[batch_indices, ] <- batch_adjusted
  }
  
  colnames(adjusted_data) <- colnames(data)
  rownames(adjusted_data) <- rownames(data)
  
  return(adjusted_data)
}