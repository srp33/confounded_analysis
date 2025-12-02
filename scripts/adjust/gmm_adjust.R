# ============================================================================
# gmm_adjust.R - 1D Gaussian Mixture Model with posterior-mean priors
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(foreach)
  library(doParallel)
})

# -----------------------------------------------------------------------------
# GMM MODEL CONSTRUCTOR
# -----------------------------------------------------------------------------

#' 1D 2-component Gaussian Mixture Model
#'
#' @param max_iter Maximum number of EM iterations
#' @param tol Convergence tolerance for log-likelihood
#' @param weight_alpha Dirichlet prior pseudo-count for mixture weights (>= 1.0)
#' @param variance_alpha Inverse-Gamma prior shape parameter for variances (>1 for posterior mean to exist)
GaussianMixture1D <- function(
  max_iter = 100,
  tol = 1e-4,
  weight_alpha = NULL,
  variance_alpha = NULL
) {
  structure(list(
    max_iter = max_iter,
    tol = tol,
    weight_alpha = weight_alpha,
    variance_alpha = variance_alpha,
    means_ = NULL,
    variances_ = NULL,
    weights_ = NULL,
    resp_ = NULL
  ), class = "GaussianMixture1D")
}

# -----------------------------------------------------------------------------
# NORMAL PDF
# -----------------------------------------------------------------------------

normal_pdf <- function(x, mean, sd) {
  exp(-0.5 * ((x - mean) / sd)^2) / (sd * sqrt(2 * pi))
}

# -----------------------------------------------------------------------------
# FIT FUNCTION
# -----------------------------------------------------------------------------

fit.GaussianMixture1D <- function(model, X) {
  X <- as.vector(X)
  n <- length(X)
  K <- 2
  eps <- 1e-12

  # Initialize means at 25th and 75th percentiles
  means <- quantile(X, c(0.25, 0.75), na.rm = TRUE)
  variances <- rep(var(X, na.rm = TRUE) + 1e-6, K)
  weights <- rep(0.5, K)

  if (is.null(model$weight_alpha)) {
    model$weight_alpha <- 3.0 + n / 100.0
  }
  if (is.null(model$variance_alpha)) {
    model$variance_alpha <- 6.0 + n / 50.0
  }

  log_likelihood_old <- -Inf

  alpha_w <- model$weight_alpha
  alpha_v <- max(model$variance_alpha, 1.01)  # ensure posterior mean exists

  for (iter in 1:model$max_iter) {
    # -------------------------------
    # E-step: responsibilities
    # -------------------------------
    pdfs <- matrix(0, nrow = n, ncol = K)
    for (k in 1:K) {
      pdfs[, k] <- weights[k] * normal_pdf(X, means[k], sqrt(variances[k]))
    }
    resp <- pdfs / (rowSums(pdfs) + eps)

    # -------------------------------
    # M-step: posterior mean updates
    # -------------------------------
    Nk <- colSums(resp)
    variances_old <- variances

    # 1) Update weights: posterior mean of Dirichlet
    weights <- (Nk + alpha_w) / (n + K * alpha_w)

    # 2) Update means (weighted ML, no prior)
    for (k in 1:K) {
      means[k] <- sum(resp[, k] * X) / (Nk[k] + eps)
    }

    # 3) Update variances: inverse-gamma prior, coupled across components
    for (k in 1:K) {
      other_k <- 3 - k  # if k=1, other=2; if k=2, other=1
      if (Nk[k] < 1e-8) {
        variances[k] <- variances_old[k]
        next
      }

      S_k <- sum(resp[, k] * (X - means[k])^2)
      v_other <- variances_old[other_k]
      beta0 <- (alpha_v - 1) * v_other

      alpha_post <- alpha_v + 0.5 * Nk[k]
      beta_post <- beta0 + 0.5 * S_k

      if (alpha_post <= 1 + 1e-12) {
        variances[k] <- max(variances_old[k], 1e-6)
      } else {
        variances[k] <- beta_post / (alpha_post - 1)
      }
    }

    variances <- pmax(variances, 1e-6)

    # -------------------------------
    # Convergence check
    # -------------------------------
    log_likelihood <- sum(log(rowSums(pdfs) + eps))
    if (abs(log_likelihood - log_likelihood_old) < model$tol) break
    log_likelihood_old <- log_likelihood
  }

  model$means_ <- means
  model$variances_ <- variances
  model$weights_ <- weights
  model$resp_ <- resp

  return(model)
}


# -----------------------------------------------------------------------------
# SIMPLE FALLBACK
# -----------------------------------------------------------------------------

simple_fallback <- function(gene_exp) {
  min_val <- min(gene_exp, na.rm = TRUE)
  x_transformed <- log(gene_exp - min_val + 1)
  n_valid <- sum(!is.na(x_transformed))
  if (n_valid > 1) {
    quantiles <- rank(x_transformed, na.last = "keep", ties.method = "average") / (n_valid + 1)
    qnorm(quantiles)
  } else x_transformed
}


# -----------------------------------------------------------------------------
# GENE TRANSFORM / BIMODAL NORMALIZE / BATCH ADJUST
# -----------------------------------------------------------------------------

get_gene_gmm_transform <- function(
    gene_exp,
    weight_alpha = NULL,
    variance_alpha = NULL,
    mean_mean_zero = TRUE,
    mean1_zero = FALSE,
    unit_var = TRUE,
    means_at_1 = FALSE,
    diff_exp = FALSE,
    output_counts = FALSE,
    log_transform = TRUE
) {
  if (all(is.na(gene_exp))) return(gene_exp)

  # --- Log-transform (optional) ---
  if (log_transform) {
    min_val <- min(gene_exp, na.rm = TRUE)
    x_transformed <- log(gene_exp - min_val + 1)
  } else {
    x_transformed <- gene_exp
  }
  mean_shift_fallback <- scale(x_transformed, scale = FALSE)[, 1]

  # Check for small variance
  if (var(x_transformed, na.rm = TRUE) < 1e-8) return(mean_shift_fallback)

  # --- Quantiles ---
  ranks <- rank(x_transformed, na.last = "keep", ties.method = "average")
  n_valid <- sum(!is.na(x_transformed))
  quantiles <- ranks / (n_valid + 1)

  # --- Fit 2-component GMM ---
  gmm <- GaussianMixture1D(weight_alpha = weight_alpha, variance_alpha = variance_alpha)
  gmm <- fit.GaussianMixture1D(gmm, x_transformed)

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

  # --- Affine corrections ---
  if (mean1_zero) {
    # Adjusts the first mean to be 0
    x_transformed <- x_transformed - means[1]
  }

  if (mean_mean_zero) {
    # Centers the means on zero, one on either side
    mean_center <- 0.5 * (means[1] + means[2])
    x_transformed <- x_transformed - mean_center
  }

  if (unit_var) {
    # Scales the variance to be 1, assuming equal weights
    variance <- 0.5 * (variances[1] + variances[2]) + 0.25 * (means[2] - means[1])^2
    scale_factor <- sqrt(max(variance, 1e-9))
    x_transformed <- x_transformed / scale_factor
  }

  if (means_at_1) {
    # Compute safe scaling so component means map to ±1
    scale_factor <- (means[2] - means[1]) / 2
    if (abs(scale_factor) < 1e-6) {
      return(mean_shift_fallback)
    }
    x_transformed <- x_transformed / scale_factor
  }

  if (output_counts) {
    # Output integral data
    x_transformed = round(exp(x_transformed) * 250)
  }

  return(x_transformed)
}

#' Bimodal normalization using GMM with parallel processing
bimodal_normalize <- function(data, weight_alpha=NULL, variance_alpha=NULL, mean_mean_zero = TRUE, unit_var = TRUE, 
                             mean1_zero = FALSE, diff_exp = FALSE, means_at_1 = FALSE, 
                             output_counts = FALSE, log_transform = TRUE, debug = FALSE, num_workers = NULL) {
  gene_names <- colnames(data)
  if (is.null(gene_names)) gene_names <- paste0("Gene", 1:ncol(data))

  if (diff_exp && unit_var) stop("Unit variance not allowed for diff exp")
  if (diff_exp && means_at_1) stop("Means at 1 not allowed for diff exp")
  if (means_at_1 && unit_var) stop("Cannot have both means_at_1 and unit_var")
  if (means_at_1 && !mean_mean_zero) stop("Cannot have means_at_1 without mean_mean_zero")
  
  if (debug) message("%%%%%%% Num workers: ", num_workers)
  # Use parallel processing if num_workers > 1
  if (!is.null(num_workers) && (num_workers != 1)) {
    num_cores <- if (num_workers == -1) detectCores() else min(num_workers, detectCores())
    if (debug) message("%%%%%%% Num cores: ", num_cores)
    if (.Platform$OS.type == "unix") {
      # Prefer fork clusters on Unix systems (no port conflicts)
      cl <- tryCatch({
        if (debug) cat("Creating fork cluster\n")
        parallel::makeForkCluster(num_cores)
      }, error = function(e) {
        if (debug) cat("Fork cluster failed, using PSOCK with retry:", e$message, "\n")
        # Fallback to PSOCK with retry logic
        cl_temp <- NULL
        max_attempts <- 5
        for (attempt in 1:max_attempts) {
          tryCatch({
            cl_temp <- parallel::makePSOCKcluster(num_cores)
            return(cl_temp)
          }, error = function(e2) {
            if (debug) cat("PSOCK cluster attempt", attempt, "failed:", e2$message, "\n")
            if (attempt == max_attempts) {
              stop("Failed to create cluster after ", max_attempts, " attempts: ", e2$message)
            }
            Sys.sleep(runif(1, 0.5, 2))  # Random delay before retry
          })
        }
      })
    } else {
      if (debug) cat("Creating PSOCK cluster\n")
      # Try to create PSOCK cluster with port retry logic
      cl <- NULL
      max_attempts <- 5
      for (attempt in 1:max_attempts) {
        tryCatch({
          cl <- parallel::makePSOCKcluster(num_cores)
          break
        }, error = function(e) {
          if (debug) cat("Cluster creation attempt", attempt, "failed:", e$message, "\n")
          if (attempt == max_attempts) {
            stop("Failed to create cluster after ", max_attempts, " attempts: ", e$message)
          }
          Sys.sleep(runif(1, 0.5, 2))  # Random delay before retry
        })
      }
    }
    registerDoParallel(cl)
    on.exit(stopCluster(cl), add = TRUE)
    
    # Explicitly export all required functions to avoid namespace collisions
    results <- foreach(i = seq_along(gene_names), .combine = cbind, 
                      .packages = c("stats"),
                      .export = c("get_gene_gmm_transform", "simple_fallback", "GaussianMixture1D",
                                "fit.GaussianMixture1D", "normal_pdf"),
                      .errorhandling = 'remove') %dopar% {
      gene_exp <- data[, i]
      if (all(is.na(gene_exp)) || all(gene_exp == gene_exp[1], na.rm = TRUE)) {
        return(gene_exp)
      }
      
      tryCatch({
        get_gene_gmm_transform(gene_exp, weight_alpha, variance_alpha, mean_mean_zero, mean1_zero, 
                              unit_var, means_at_1, diff_exp, output_counts, log_transform)
      }, error = function(e) {
        simple_fallback(gene_exp)
      })
    }
    bimodal_data <- as.matrix(results)
  } else {
    # Sequential processing
    if (debug) message("Doing sequential processing!")
    bimodal_data <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
    for (i in seq_along(gene_names)) {
      gene_exp <- data[, i]
      if (all(is.na(gene_exp)) || all(gene_exp == gene_exp[1], na.rm = TRUE)) {
        if (debug) message("All NA for gene: ", gene_names[i])
        bimodal_data[, i] <- gene_exp
        next
      }
      
      tryCatch({
        bimodal_data[, i] <- get_gene_gmm_transform(gene_exp, weight_alpha, variance_alpha, mean_mean_zero, mean1_zero, 
                                                   unit_var, means_at_1, diff_exp, output_counts, log_transform)
      }, error = function(e) {
        if (debug) message("Got Error, simple fallback")
        if (debug) message(e)
        bimodal_data[, i] <- simple_fallback(gene_exp)
      })

      if (debug) {
        message("Min: ", min(bimodal_data[, i], na.rm = TRUE))
        message("Max: ", max(bimodal_data[, i], na.rm = TRUE))
        message("Num non-finite: ", sum(!is.finite(bimodal_data[, i])))
      }
    }
  }
  
  colnames(bimodal_data) <- gene_names
  rownames(bimodal_data) <- rownames(data)
  return(bimodal_data)
}

#' GMM adjustment for multiple batches
#' 
#' @param data Matrix of gene expression data (samples x genes)
#' @param batch Vector of batch labels for each sample
#' @param weight_alpha Dirichlet prior parameter for GMM weights
#' @param variance_alpha 
#' @param mean_mean_zero If TRUE, center means around zero (default TRUE)
#' @param unit_var If TRUE, scale to unit variance (default TRUE)
#' @param diff_exp If TRUE, adjust first mean to zero for differential expression preservation (default FALSE)
#' @param means_at_1 If TRUE, place means at ±1 (default FALSE)
#' @param output_counts If TRUE, attempt to preserve count structure (default FALSE)
#' @param log_transform If TRUE, apply log-transformation to data (default TRUE). Set to FALSE for pre-transformed data.
#' @param debug If TRUE, print progress messages
#' 
#' @details The function applies GMM-based bimodal normalization to each batch separately.
#' Various transformation options allow control over the final distribution properties.
gmm_adjust <- function(data, batch, genes_are_columns=TRUE, weight_alpha=NULL, variance_alpha=NULL, mean_mean_zero = TRUE,
                      mean1_zero = FALSE, unit_var = TRUE, diff_exp = FALSE, means_at_1 = FALSE, 
                      output_counts = FALSE, log_transform = TRUE, debug = FALSE, num_workers = NULL) {
  batch_factor <- as.factor(batch)
  batch_levels <- levels(batch_factor)
  if (!genes_are_columns) {
    if (debug) message("Transposing data")
    data <- t(data)
  }
  adjusted_data <- matrix(NA, nrow = nrow(data), ncol = ncol(data))

  if (debug) message("Batch: ",  batch)
  
  for (b in batch_levels) {
    batch_indices <- which(batch == b)
    if (debug) message("Batch indices: ", batch_indices)
    batch_data <- data[batch_indices, , drop = FALSE]
    
    batch_adjusted <- bimodal_normalize(
      batch_data, 
      weight_alpha=weight_alpha, variance_alpha=variance_alpha,
      mean_mean_zero=mean_mean_zero, unit_var=unit_var,
      mean1_zero=mean1_zero, diff_exp=diff_exp, means_at_1=means_at_1, 
      output_counts=output_counts, log_transform=log_transform,
      debug=debug, num_workers=num_workers
    )
    
    adjusted_data[batch_indices, ] <- batch_adjusted
  }
  
  colnames(adjusted_data) <- colnames(data)
  rownames(adjusted_data) <- rownames(data)

  if (!genes_are_columns) {
    if (debug) message("Transposing data back")
    adjusted_data <- t(adjusted_data)
  }
  return(adjusted_data)
}