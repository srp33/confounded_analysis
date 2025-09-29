# gmm_adjust_streamlined.R
# Streamlined GMM adjustment - bimodal adjustment for all genes

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# ============================================================================
# STREAMLINED GMM IMPLEMENTATION
# ============================================================================

#' Streamlined 2-component Gaussian Mixture Model
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

#' Inverse CDF for GMM distribution
inverse_cdf_gmm <- function(p, means, variances, weights) {
  # Input validation
  if (any(is.na(means)) || any(is.na(variances)) || any(is.na(weights))) {
    return(rep(NA, length(p)))
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
  
  solve_for_single_p <- function(p_val) {
    if (is.na(p_val)) return(NA)
    if (p_val <= 0.0) return(-Inf)
    if (p_val >= 1.0) return(Inf)
    
    root_function <- function(x) gmm_cdf(x) - p_val
    
    return(uniroot(root_function, interval = search_interval, tol = 1e-4)$root)
  }
  
  result <- sapply(p, solve_for_single_p)
  
  # Replace any remaining problematic values
  if (any(is.infinite(result) | is.na(result))) {
    bad_indices <- which(is.infinite(result) | is.na(result))
    result[bad_indices] <- qnorm(p[bad_indices])
  }
  
  return(result)
}

# ============================================================================
# STREAMLINED ADJUSTMENT FUNCTIONS
# ============================================================================
#' Get GMM parameters and apply mean-only centering and scaling transformation
get_gene_gmm_transform <- function(gene_exp, alpha0 = 10, mean_only = FALSE) {
  # Apply log transformation first
  min_val <- min(gene_exp, na.rm = TRUE)
  x_transformed <- log(gene_exp - min_val + 1)

  # Fit 2-component GMM
  gmm <- GaussianMixture2D(alpha0 = alpha0)
  gmm <- fit(gmm, x_transformed)

  # Extract parameters (ensure lower component is first)
  if (gmm$means_[1] > gmm$means_[2]) {
    means <- c(gmm$means_[2], gmm$means_[1])
    variances <- c(gmm$variances_[2], gmm$variances_[1])
    weights <- c(gmm$weights_[2], gmm$weights_[1])
  } else {
    means <- gmm$means_
    variances <- gmm$variances_
    weights <- gmm$weights_
  }

  if (mean_only) {
    # Move mean center to 0
    mean_center <- sum(means * weights)  # weighted mean of components
    x_centered <- x_transformed - mean_center

    # Compute variance as if weights were 0.5 each
    avg_variance <- sum(0.5 * (variances + (means - mean_center)^2))
    scale_factor <- sqrt(avg_variance)

    # Scale down the data
    result <- x_centered / scale_factor

  } else {
    # Transform to bimodal normal with means at -1 and 1
    new_mean1 <- -1
    new_mean2 <- 1

    old_mean_diff <- means[2] - means[1]
    new_mean_diff <- new_mean2 - new_mean1

    new_sd1 <- sqrt(variances[1]) * new_mean_diff / old_mean_diff
    new_sd2 <- sqrt(variances[2]) * new_mean_diff / old_mean_diff

    new_var1 <- new_sd1^2
    new_var2 <- new_sd2^2

    # Get quantiles from original data
    ranks <- rank(x_transformed, na.last = "keep", ties.method = "average")
    quantiles <- ranks / (sum(!is.na(x_transformed)) + 1)

    # Apply inverse CDF transformation
    result <- inverse_cdf_gmm(
      quantiles,
      means = c(new_mean1, new_mean2),
      variances = c(new_var1, new_var2),
      weights = weights
    )
  }

  return(result)
}


#' Bimodal normalization using streamlined GMM with full inverse CDF
bimodal_normalize_streamlined <- function(data, alpha0 = 10, mean_only = FALSE, debug = FALSE) {
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
      bimodal_data[, i] <- bimodal_result
      
    }, error = function(e) {
      if (debug) {
        cat("Error processing gene", gene_name, ":", e$message, "\n")
      }
      # Fallback to simple quantile normalization
      min_val <- min(gene_exp, na.rm = TRUE)
      x_transformed <- log(gene_exp - min_val + 1)
      quantiles <- rank(x_transformed, na.last = "keep", ties.method = "average") / (sum(!is.na(x_transformed)) + 1)
      fallback_result <- qnorm(quantiles)
      bimodal_data[, i] <- fallback_result / sd(fallback_result, na.rm = TRUE)
    })
  }
  
  colnames(bimodal_data) <- gene_names
  rownames(bimodal_data) <- rownames(data)
  
  return(bimodal_data)
}

#' Streamlined GMM adjustment for multiple batches
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
gmm_adjust_streamlined <- function(data, batch, alpha0 = 10, mean_only = FALSE, debug = FALSE) {
  if (debug) {
    cat("Starting streamlined GMM adjustment...\n")
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
    batch_adjusted <- bimodal_normalize_streamlined(batch_data, alpha0, mean_only, debug)
    adjusted_data[batch_indices, ] <- batch_adjusted
  }
  
  colnames(adjusted_data) <- colnames(data)
  rownames(adjusted_data) <- rownames(data)
  
  return(adjusted_data)
}
