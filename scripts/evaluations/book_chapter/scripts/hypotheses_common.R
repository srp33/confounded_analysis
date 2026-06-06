#!/usr/bin/env Rscript
# hypotheses_common.R
# Shared setup sourced by every per-hypothesis script (scripts/hypotheses/H*.R)
# and by the sequential driver evaluate_hypotheses.R.
#
# Provides: libraries, helpers, data load, primary-scenario ComBat fits,
# the H1 within-class variance ratio (reused by H4), and a file-writing add().
# Each hypothesis writes its own verdict CSV (verdict_H<k>.csv) so that the
# hypotheses can be evaluated in parallel (e.g. a SLURM job array) and merged
# afterwards by aggregate_verdicts.R.

suppressMessages({
  library(dplyr); library(tidyr); library(sva); library(limma); library(matrixStats)
})
source("scripts/ComBat_nat.R")

DATA_FILE <- "data/TB_real_data.RData"
ADJ_DIR   <- "outputs/adjusted_data/all_scenarios"
OUT_DIR   <- "outputs/diagnostics/hypothesis_tests"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
load(DATA_FILE)

ALL_STUDIES <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1", 1, "Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES != test][seq_len(n)]

# ── Generic helpers ───────────────────────────────────────────────────────────

knn_pur <- function(ref_mat, ref_lab, tst_mat, tst_lab, k = 5, n_pcs = 50) {
  cg  <- intersect(rownames(ref_mat), rownames(tst_mat))
  pca <- prcomp(t(ref_mat[cg, ]), center = TRUE, scale. = FALSE, rank. = n_pcs)
  rr  <- pca$x
  tr  <- predict(pca, t(tst_mat[cg, ]))
  d   <- as.matrix(dist(rbind(tr, rr)))[
    seq_len(nrow(tr)), (nrow(tr) + 1):(nrow(tr) + nrow(rr))]
  mean(sapply(seq_len(nrow(tr)), function(i)
    mean(ref_lab[order(d[i, ])[seq_len(k)]] == tst_lab[i])))
}

cls_sep_test <- function(ref_mat, ref_lab, tst_mat, tst_lab) {
  cg <- intersect(rownames(ref_mat), rownames(tst_mat))
  w  <- rowMeans(ref_mat[cg, ref_lab == 1, drop = FALSE]) -
        rowMeans(ref_mat[cg, ref_lab == 0, drop = FALSE])
  w  <- w / sqrt(sum(w^2))
  s  <- as.numeric(t(tst_mat[cg, ]) %*% w)
  mean(s[tst_lab == 1]) - mean(s[tst_lab == 0])
}

var_decomp <- function(mat, lab) {
  lab <- as.integer(lab)
  n0  <- sum(lab == 0); n1 <- sum(lab == 1)
  v0  <- rowVars(mat[, lab == 0, drop = FALSE])
  v1  <- rowVars(mat[, lab == 1, drop = FALSE])
  within  <- ((n0 - 1) * v0 + (n1 - 1) * v1) / (n0 + n1 - 2)
  total   <- rowVars(mat)
  sep     <- rowMeans(mat[, lab == 1, drop = FALSE]) -
             rowMeans(mat[, lab == 0, drop = FALSE])
  list(within = within, total = total, sep = sep)
}

prep <- function(n, test) {
  refs <- train_studies(n, test)
  cg   <- Reduce(intersect, lapply(c(refs, test),
                                   function(s) rownames(dat_lst[[s]])))
  dat  <- do.call(cbind, lapply(refs, function(s) dat_lst[[s]][cg, ]))
  bat  <- do.call(c, lapply(refs, function(s) rep(s, ncol(dat_lst[[s]]))))
  lab  <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))
  tst  <- dat_lst[[test]][cg, ]
  ltst <- bin(label_lst[[test]])
  # Log transform matching the main pipeline (shift negative values first)
  log_safe <- function(m) {
    if (max(m, na.rm = TRUE) > 100) {
      mn <- min(m, na.rm = TRUE)
      if (mn < 0) m <- m - mn
      m <- log2(m + 1)
    }
    m
  }
  dat <- log_safe(dat); tst <- log_safe(tst)
  list(dat = dat, bat = bat, lab = lab, tst = tst, ltst = ltst, genes = cg)
}

global_scale <- function(train, test) {
  m <- mean(train); s <- sd(as.vector(train))
  list(train = (train - m) / s, test = (test - m) / s)
}

# Inline combat_sup pipeline (replicates adjust_target_data_SA_UK.R)
run_combat_sup <- function(dat, bat, lab, tst, ltst) {
  dat_c <- if (length(unique(bat)) >= 2)
    suppressMessages(ComBat(dat, batch = bat, mod = model.matrix(~lab)))
  else dat
  comb  <- cbind(dat_c, tst)
  cb    <- c(rep(1L, ncol(dat_c)), rep(2L, ncol(tst)))
  clabs <- c(lab, ltst)
  out   <- suppressMessages(
    ComBat(comb, batch = cb, mod = model.matrix(~clabs), ref.batch = 1L))
  list(ref = dat_c,
       tgt = out[, (ncol(dat_c) + 1):ncol(out), drop = FALSE])
}

run_combat_sup_nat <- function(dat, bat, lab, tst) {
  dat_c <- if (length(unique(bat)) >= 2)
    limma::removeBatchEffect(dat, batch = bat, design = model.matrix(~lab))
  else dat
  comb  <- cbind(dat_c, tst)
  cb    <- c(rep(1L, ncol(dat_c)), rep(2L, ncol(tst)))
  out   <- suppressMessages(
    ComBat(comb, batch = cb, mod = NULL, ref.batch = 1L, mean.only = TRUE))
  list(ref = dat_c,
       tgt = out[, (ncol(dat_c) + 1):ncol(out), drop = FALSE])
}

# Primary scenario: n=3, test=USA
sc <- prep(3, "USA")
suppressMessages({
  tr_unsup <- ComBat(sc$dat, batch = sc$bat, mod = NULL)
  tr_sup   <- ComBat(sc$dat, batch = sc$bat, mod = model.matrix(~sc$lab))
  tr_mo_cls <- ComBat(sc$dat, batch = sc$bat, mod = model.matrix(~sc$lab),
                      mean.only = TRUE)
  tr_mo_raw <- ComBat(sc$dat, batch = sc$bat, mod = NULL, mean.only = TRUE)
  tr_nat    <- limma::removeBatchEffect(sc$dat, batch = sc$bat,
                                        design = model.matrix(~sc$lab))
  tr_nat_v2 <- tr_nat  # same step-1 as nat; step-2 difference only
})

# H1 within-class variance ratio — computed here because H4 also depends on it.
vd_u <- var_decomp(tr_unsup, sc$lab)
vd_s <- var_decomp(tr_sup,   sc$lab)
ratio_h1  <- vd_s$within / pmax(vd_u$within, 1e-10)
med_ratio <- median(ratio_h1, na.rm = TRUE)

# ── Verdict writer ────────────────────────────────────────────────────────────
# In parallel mode each hypothesis writes its own one-row verdict CSV, keyed by
# the leading token of the hypothesis label ("H1", "H4", ...). aggregate_verdicts.R
# merges them into verdict_table.csv.
add <- function(h, v, stat) {
  cat(sprintf("\n  VERDICT [%s] %s\n  %s\n", v, h, stat))
  tag <- gsub("[^A-Za-z0-9]+", "_", trimws(strsplit(h, " ")[[1]][1]))  # "H1","H4",...
  write.csv(data.frame(H = h, verdict = v, key_stat = stat, stringsAsFactors = FALSE),
            file.path(OUT_DIR, sprintf("verdict_%s.csv", tag)), row.names = FALSE)
}

cat("═══════════════════════════════════════════════════════════════════\n")
cat("  HYPOTHESIS EVALUATIONS (shared setup loaded)\n")
cat("═══════════════════════════════════════════════════════════════════\n")
