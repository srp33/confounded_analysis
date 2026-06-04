# ComBat_nat: supervised ComBat with batch-only variance for reconstruction.
#
# Hypothesis: standard supervised ComBat compresses within-class variance to
# sigma^2_residual (full-model residual) because var.pooled is computed after
# removing BOTH batch AND class effects. When the test batch is then projected
# via unsupervised step-2, step-2 sees training data with tiny within-class
# spread relative to the test's natural variance, causing delta.star to crush
# the test class signal.
#
# This variant keeps everything identical EXCEPT the reconstruction scale:
# var.pooled.nat is computed from batch-only residuals (class signal still
# present), so within-class spread in the corrected training output is
# sigma_batchonly = sqrt(sigma^2_class + sigma^2_residual) — closer to the
# natural scale that step-2 expects.
#
# Usage: source("scripts/ComBat_nat.R"), then call ComBat_nat() in place of
# sva::ComBat() with the same arguments.

ComBat_nat <- function(dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE,
                       mean.only = FALSE, ref.batch = NULL,
                       BPPARAM = BiocParallel::bpparam("SerialParam")) {

  if (length(dim(batch)) > 1) stop("This version of ComBat only allows one batch variable")

  dat   <- as.matrix(dat)
  batch <- as.factor(batch)

  # Remove genes with zero within-batch variance
  zero.rows.lst <- lapply(levels(batch), function(bl) {
    if (sum(batch == bl) > 1)
      which(apply(dat[, batch == bl], 1, function(x) var(x) == 0))
    else
      integer(0)
  })
  zero.rows <- Reduce(union, zero.rows.lst)
  keep.rows <- setdiff(seq_len(nrow(dat)), zero.rows)
  if (length(zero.rows) > 0) {
    cat(sprintf("Found %d genes with uniform expression within a batch; not adjusted.\n",
                length(zero.rows)))
    dat.orig <- dat
    dat      <- dat[keep.rows, ]
  }

  if (any(table(batch) == 1)) mean.only <- TRUE
  if (mean.only) message("Using the 'mean only' version of ComBat")

  batchmod <- model.matrix(~-1 + batch)
  if (!is.null(ref.batch)) {
    if (!(ref.batch %in% levels(batch)))
      stop("ref.batch is not one of the levels of the batch variable")
    message("Using batch =", ref.batch, "as reference batch")
    ref <- which(levels(as.factor(batch)) == ref.batch)
    batchmod[, ref] <- 1
  } else {
    ref <- NULL
  }
  message("Found", nlevels(batch), "batches")

  n.batch   <- nlevels(batch)
  batches   <- lapply(levels(batch), function(bl) which(batch == bl))
  n.batches <- sapply(batches, length)
  if (any(n.batches == 1)) { mean.only <- TRUE; message("One batch has 1 sample; mean.only=TRUE") }
  n.array   <- sum(n.batches)

  design <- cbind(batchmod, mod)
  check  <- apply(design, 2, function(x) all(x == 1))
  if (!is.null(ref)) check[ref] <- FALSE
  design <- as.matrix(design[, !check])

  message("Adjusting for", ncol(design) - ncol(batchmod), "covariate(s)")

  if (qr(design)$rank < ncol(design)) {
    if (ncol(design) == (n.batch + 1))
      stop("The covariate is confounded with batch!")
    if (ncol(design) > (n.batch + 1)) {
      if (qr(design[, -c(1:n.batch)])$rank < ncol(design[, -c(1:n.batch)]))
        stop("The covariates are confounded!")
      else
        stop("At least one covariate is confounded with batch!")
    }
  }

  NAs <- any(is.na(dat))
  if (NAs) message("Found", sum(is.na(dat)), "missing values")

  message("Standardizing Data across genes")
  if (!NAs) {
    B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
  } else {
    B.hat <- apply(dat, 1, sva:::Beta.NA, design)
  }

  if (!is.null(ref.batch)) {
    grand.mean <- t(B.hat[ref, ])
  } else {
    grand.mean <- crossprod(n.batches / n.array, B.hat[1:n.batch, ])
  }

  # --- standard var.pooled (full-model residual, used for s.data and gamma/delta) ---
  if (!NAs) {
    if (!is.null(ref.batch)) {
      ref.dat    <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat - t(design[batches[[ref]], ] %*% B.hat))^2) %*%
                    rep(1/n.batches[ref], n.batches[ref])
    } else {
      var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array, n.array)
    }
  } else {
    if (!is.null(ref.batch)) {
      ref.dat    <- dat[, batches[[ref]]]
      var.pooled <- matrixStats::rowVars(ref.dat - t(design[batches[[ref]], ] %*% B.hat), na.rm = TRUE)
    } else {
      var.pooled <- matrixStats::rowVars(dat - t(design %*% B.hat), na.rm = TRUE)
    }
  }

  # --- natural-scale var.pooled (batch-only residual, used for reconstruction) ---
  # This retains the class signal in the residuals, giving larger pooled variance.
  # Reconstruction with this scale avoids compressing within-class spread to the
  # tiny full-model residual — the key change relative to standard supervised ComBat.
  batch.cols <- 1:n.batch
  batch.pred <- t(design[, batch.cols, drop = FALSE] %*% B.hat[batch.cols, , drop = FALSE])
  if (!NAs) {
    var.pooled.nat <- ((dat - batch.pred)^2) %*% rep(1/n.array, n.array)
  } else {
    var.pooled.nat <- matrixStats::rowVars(dat - batch.pred, na.rm = TRUE)
  }

  stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
  if (!is.null(design)) {
    tmp <- design
    tmp[, seq_len(n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% B.hat)
  }
  s.data <- (dat - stand.mean) / (sqrt(var.pooled) %*% t(rep(1, n.array)))

  message("Fitting L/S model and finding priors")
  batch.design <- design[, seq_len(n.batch)]
  if (!NAs) {
    gamma.hat <- solve(crossprod(batch.design),
                       tcrossprod(t(batch.design), as.matrix(s.data)))
  } else {
    gamma.hat <- apply(s.data, 1, sva:::Beta.NA, batch.design)
  }
  delta.hat <- NULL
  for (i in batches) {
    if (mean.only) {
      delta.hat <- rbind(delta.hat, rep(1, nrow(s.data)))
    } else {
      delta.hat <- rbind(delta.hat, matrixStats::rowVars(s.data[, i, drop = FALSE], na.rm = TRUE))
    }
  }

  gamma.bar <- rowMeans(gamma.hat)
  t2        <- matrixStats::rowVars(gamma.hat)
  a.prior   <- apply(delta.hat, 1, sva:::aprior)
  b.prior   <- apply(delta.hat, 1, sva:::bprior)

  gamma.star <- delta.star <- matrix(NA, nrow = n.batch, ncol = nrow(s.data))
  if (par.prior) {
    message("Finding parametric adjustments")
    results <- BiocParallel::bplapply(seq_len(n.batch), function(i) {
      if (mean.only) {
        gs <- sva:::postmean(gamma.hat[i, ], gamma.bar[i], 1, 1, t2[i])
        ds <- rep(1, nrow(s.data))
      } else {
        temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                             delta.hat[i, ], gamma.bar[i], t2[i],
                             a.prior[i], b.prior[i])
        gs <- temp[1, ]
        ds <- temp[2, ]
      }
      list(gamma.star = gs, delta.star = ds)
    }, BPPARAM = BPPARAM)
    for (i in seq_len(n.batch)) {
      gamma.star[i, ] <- results[[i]]$gamma.star
      delta.star[i, ] <- results[[i]]$delta.star
    }
  } else {
    message("Finding nonparametric adjustments")
    results <- BiocParallel::bplapply(seq_len(n.batch), function(i) {
      if (mean.only) delta.hat[i, ] <- 1
      temp <- sva:::int.eprior(as.matrix(s.data[, batches[[i]]]),
                               gamma.hat[i, ], delta.hat[i, ])
      list(gamma.star = temp[1, ], delta.star = temp[2, ])
    }, BPPARAM = BPPARAM)
    for (i in seq_len(n.batch)) {
      gamma.star[i, ] <- results[[i]]$gamma.star
      delta.star[i, ] <- results[[i]]$delta.star
    }
  }

  if (!is.null(ref.batch)) {
    gamma.star[ref, ] <- 0
    delta.star[ref, ] <- 1
  }

  message("Adjusting the Data\n")
  bayesdata <- s.data
  j <- 1
  for (i in batches) {
    bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i, ] %*% gamma.star)) /
                      (sqrt(delta.star[j, ]) %*% t(rep(1, n.batches[j])))
    j <- j + 1
  }

  # KEY CHANGE: reconstruct using batch-only var.pooled (natural scale)
  # rather than full-model var.pooled (compressed scale).
  bayesdata <- (bayesdata * (sqrt(var.pooled.nat) %*% t(rep(1, n.array)))) + stand.mean

  if (!is.null(ref.batch)) {
    bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }

  if (length(zero.rows) > 0) {
    dat.orig[keep.rows, ] <- bayesdata
    bayesdata <- dat.orig
  }

  return(bayesdata)
}
