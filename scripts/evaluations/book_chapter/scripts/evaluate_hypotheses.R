#!/usr/bin/env Rscript
# evaluate_hypotheses.R
# Systematic evaluation of H1 (old) through H14 (current).
# Each section states the hypothesis, defines a decisive test, and prints
# a PASS/FAIL verdict.

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

# ── Verdict table ─────────────────────────────────────────────────────────────
VT <- data.frame(H = character(), verdict = character(),
                 key_stat = character(), stringsAsFactors = FALSE)
add <- function(h, v, stat) {
  cat(sprintf("\n  VERDICT [%s] %s\n  %s\n", v, h, stat))
  VT[nrow(VT) + 1, ] <<- list(h, v, stat)
}

cat("═══════════════════════════════════════════════════════════════════\n")
cat("  HYPOTHESIS EVALUATIONS\n")
cat("═══════════════════════════════════════════════════════════════════\n")

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H1 (old): Supervised step-1 compresses within-class variance ──────\n")
# Test: per-gene within-class variance ratio sup / unsup.
# H1 TRUE => ratio << 1.  H1 FALSE => ratio ≈ 1.

vd_u <- var_decomp(tr_unsup, sc$lab)
vd_s <- var_decomp(tr_sup,   sc$lab)
ratio_h1 <- vd_s$within / pmax(vd_u$within, 1e-10)
med_ratio <- median(ratio_h1, na.rm = TRUE)
wt_h1 <- wilcox.test(vd_s$within, vd_u$within, paired = TRUE, exact = FALSE)

cat(sprintf("  Median within-class var: unsup=%.4f  sup=%.4f\n",
            median(vd_u$within, na.rm = TRUE), median(vd_s$within, na.rm = TRUE)))
cat(sprintf("  Ratio (sup/unsup): median=%.4f  p=%.3g\n",
            med_ratio, wt_h1$p.value))

# Criterion: compress means ratio < 0.5; equal means ratio in [0.8, 1.2]
if (med_ratio < 0.5) {
  add("H1 (old)", "TRUE (unexpectedly)", sprintf("ratio=%.3f <0.5", med_ratio))
} else {
  add("H1 (old)", "FALSE", sprintf("ratio=%.3f ≈ 1 — no within-class compression", med_ratio))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H2 (old): Step-2 sees inflated test variance → large delta.hat ─────\n")
# Test: compute delta.hat for the test batch manually under supervised step-2,
# using (a) supervised step-1 training vs (b) unsupervised step-1 training.
# H2 TRUE => delta.hat much larger for (a). H2 FALSE => delta.hat ≈ 1 in both.

compute_delta_hat <- function(train_mat, train_lab, test_mat, test_lab) {
  cg   <- intersect(rownames(train_mat), rownames(test_mat))
  comb <- cbind(train_mat[cg, ], test_mat[cg, ])
  bat2 <- factor(c(rep(1, ncol(train_mat)), rep(2, ncol(test_mat))))
  cls2 <- c(train_lab, test_lab)
  design2 <- model.matrix(~0 + bat2 + cls2)
  B <- tryCatch(solve(crossprod(design2), tcrossprod(t(design2), comb)),
                error = function(e) NULL)
  if (is.null(B)) return(rep(NA, nrow(comb)))
  ref_idx <- seq_len(ncol(train_mat))
  vp  <- pmax(rowMeans((comb[, ref_idx] - t(design2[ref_idx, ] %*% B))^2), 1e-10)
  tst_idx <- (ncol(train_mat) + 1):ncol(comb)
  s_tst   <- (comb[, tst_idx] - t(design2[tst_idx, ] %*% B)) / sqrt(vp)
  rowVars(s_tst)
}

dh_sup   <- compute_delta_hat(tr_sup,   sc$lab, sc$tst, sc$ltst)
dh_unsup <- compute_delta_hat(tr_unsup, sc$lab, sc$tst, sc$ltst)

med_sup   <- median(dh_sup,   na.rm = TRUE)
med_unsup <- median(dh_unsup, na.rm = TRUE)
cat(sprintf("  Median delta.hat test batch: unsup_step1=%.3f  sup_step1=%.3f\n",
            med_unsup, med_sup))
# Fraction of genes with delta.hat > 2
f_sup   <- mean(dh_sup   > 2, na.rm = TRUE)
f_unsup <- mean(dh_unsup > 2, na.rm = TRUE)
cat(sprintf("  Fraction delta.hat > 2: unsup=%.3f  sup=%.3f\n", f_unsup, f_sup))

# H2 requires delta.hat to be substantially > 1 under supervised step-1
if (med_sup > 2 && med_sup > 2 * med_unsup) {
  add("H2 (old)", "TRUE (unexpectedly)", sprintf("delta.hat sup=%.2f >> unsup=%.2f", med_sup, med_unsup))
} else {
  add("H2 (old)", "FALSE", sprintf("delta.hat sup=%.3f vs unsup=%.3f — no inflation", med_sup, med_unsup))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H3 (old): ComBat_nat inflates training variance 5–7× vs test ───────\n")
# Test: per-gene SD ratio (ComBat_nat training output) / (raw test).
# H3 confirmed (FALSE for the fix) if ratio >> 1.

tr_nat_raw_scale <- ComBat_nat(sc$dat, batch = sc$bat, mod = model.matrix(~sc$lab))
sd_train_nat <- rowSds(tr_nat_raw_scale)
sd_test_raw  <- rowSds(sc$tst)
cg3 <- intersect(names(sd_train_nat), names(sd_test_raw))
scale_ratio  <- sd_train_nat[cg3] / pmax(sd_test_raw[cg3], 1e-10)
med_r3 <- median(scale_ratio, na.rm = TRUE)
q3     <- quantile(scale_ratio, c(0.25, 0.75), na.rm = TRUE)

cat(sprintf("  Per-gene SD ratio ComBat_nat_train / raw_test: median=%.2f [%.2f, %.2f]\n",
            med_r3, q3[1], q3[2]))

if (med_r3 > 3) {
  add("H3 (old)", "CONFIRMED", sprintf("Scale inflation ratio=%.2f (predicted 5-7x)", med_r3))
} else {
  add("H3 (old)", "FALSIFIED", sprintf("Scale inflation ratio=%.2f — not as inflated as claimed", med_r3))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H4: total=σ²_class+σ²_residual, within=σ²_residual in both methods ─\n")
# Test: confirm that (within + between) ≈ total for both methods,
# and that within/total ratios match theoretical prediction.

vd_s2 <- var_decomp(tr_sup, sc$lab)
vd_u2 <- var_decomp(tr_unsup, sc$lab)

# Variance decomposition identity: total = within + (sep² × p(1-p))
p <- mean(sc$lab)
between_s <- vd_s2$sep^2 * p * (1 - p)
between_u <- vd_u2$sep^2 * p * (1 - p)

resid_s <- median((vd_s2$within + between_s) / pmax(vd_s2$total, 1e-10), na.rm = TRUE)
resid_u <- median((vd_u2$within + between_u) / pmax(vd_u2$total, 1e-10), na.rm = TRUE)

cat(sprintf("  Decomp check (should = 1): unsup=%.4f  sup=%.4f\n", resid_u, resid_s))
cat(sprintf("  Median within/total: unsup=%.4f  sup=%.4f\n",
            median(vd_u2$within / pmax(vd_u2$total, 1e-10), na.rm = TRUE),
            median(vd_s2$within / pmax(vd_s2$total, 1e-10), na.rm = TRUE)))

if (abs(resid_u - 1) < 0.05 && abs(resid_s - 1) < 0.05 &&
    abs(med_ratio - 1) < 0.20) {
  add("H4", "TRUE",
      sprintf("Decomp checks: unsup=%.3f, sup=%.3f; within-var ratio=%.3f ≈ 1",
              resid_u, resid_s, med_ratio))
} else {
  add("H4", "FALSE",
      sprintf("Decomp checks: unsup=%.3f, sup=%.3f", resid_u, resid_s))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H5: Supervised step-1 preserves class separation; unsupervised reduces it ─\n")
# Test: compare class separation on the class axis across confounded scenarios.
# Compute ratio sep_sup / sep_unsup from the pre-corrected matrices in ADJ_DIR
# across all 18 scenarios.  H5: ratio > 1 (supervised preserves more).

h5_rows <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    ref_sup  <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_reference.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    ref_uns  <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_n%d_test%s_reference.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    if (is.null(ref_sup) || is.null(ref_uns)) next
    labs <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
    cg   <- intersect(rownames(ref_sup), rownames(ref_uns))
    sep_s <- {
      w <- rowMeans(ref_sup[cg, labs==1, drop=FALSE]) -
           rowMeans(ref_sup[cg, labs==0, drop=FALSE])
      w <- w / sqrt(sum(w^2))
      s <- as.numeric(t(ref_sup[cg, ]) %*% w)
      mean(s[labs==1]) - mean(s[labs==0])
    }
    sep_u <- {
      w <- rowMeans(ref_uns[cg, labs==1, drop=FALSE]) -
           rowMeans(ref_uns[cg, labs==0, drop=FALSE])
      w <- w / sqrt(sum(w^2))
      s <- as.numeric(t(ref_uns[cg, ]) %*% w)
      mean(s[labs==1]) - mean(s[labs==0])
    }
    h5_rows[[length(h5_rows)+1]] <- data.frame(n=n, test=test,
                                                sep_sup=sep_s, sep_unsup=sep_u,
                                                ratio=sep_s/max(sep_u, 0.01))
  }
}
h5 <- do.call(rbind, h5_rows)
cat("  Class separation on training class axis (supervised vs unsupervised):\n")
print(h5[, c("n","test","sep_sup","sep_unsup","ratio")], row.names = FALSE, digits = 3)

med_h5_ratio <- median(h5$ratio, na.rm = TRUE)
wt_h5 <- wilcox.test(h5$sep_sup, h5$sep_unsup, paired = TRUE, exact = FALSE)
cat(sprintf("  Median ratio sup/unsup=%.2f; Wilcoxon p=%.3g\n",
            med_h5_ratio, wt_h5$p.value))

if (med_h5_ratio > 1.05 && wt_h5$p.value < 0.05) {
  add("H5", "TRUE",
      sprintf("sup class-axis sep > unsup in %d/%d scenarios; median ratio=%.2f (p=%.3g)",
              sum(h5$sep_sup > h5$sep_unsup), nrow(h5), med_h5_ratio, wt_h5$p.value))
} else {
  add("H5", "FALSE",
      sprintf("No consistent advantage: median ratio=%.2f (p=%.3g)", med_h5_ratio, wt_h5$p.value))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H6: combat_sup fails due to step-2 batch-class confounding ──────────\n")
# Test A: correlation of (condition number of step-2 design) with purity loss.
# Test B: subsample test to 50/50 balance; rerun combat_sup; compare purity.

h6_rows <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    sc_i <- prep(n, test)
    # Condition number of step-2 supervised design matrix
    cls_all <- c(sc_i$lab, sc_i$ltst)
    bat_all <- factor(c(rep(1, ncol(sc_i$dat)), rep(2, ncol(sc_i$tst))))
    des     <- model.matrix(~0 + bat_all + cls_all)
    kappa_  <- kappa(des, exact = TRUE)
    # Test class proportion
    prop_active <- mean(sc_i$ltst)
    h6_rows[[length(h6_rows)+1]] <- data.frame(
      n=n, test=test, kappa=kappa_, prop_active=prop_active)
  }
}
h6 <- do.call(rbind, h6_rows)

# Load purity at k=5 from precomputed data (H6: purity loss ~ condition number)
purity_sup <- sapply(seq_len(nrow(h6)), function(i) {
  ref <- tryCatch(as.matrix(read.csv(
    file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_reference.csv", h6$n[i], h6$test[i])),
    row.names=1)), error=function(e) NULL)
  tgt <- tryCatch(as.matrix(read.csv(
    file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_target.csv",    h6$n[i], h6$test[i])),
    row.names=1)), error=function(e) NULL)
  if (is.null(ref) || is.null(tgt)) return(NA)
  labs_r <- do.call(c, lapply(train_studies(h6$n[i], h6$test[i]), function(s) bin(label_lst[[s]])))
  labs_t <- bin(label_lst[[h6$test[i]]])
  knn_pur(ref, labs_r, tgt, labs_t, k=5)
})
purity_nat <- sapply(seq_len(nrow(h6)), function(i) {
  ref <- tryCatch(as.matrix(read.csv(
    file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_reference.csv", h6$n[i], h6$test[i])),
    row.names=1)), error=function(e) NULL)
  tgt <- tryCatch(as.matrix(read.csv(
    file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_target.csv",    h6$n[i], h6$test[i])),
    row.names=1)), error=function(e) NULL)
  if (is.null(ref) || is.null(tgt)) return(NA)
  labs_r <- do.call(c, lapply(train_studies(h6$n[i], h6$test[i]), function(s) bin(label_lst[[s]])))
  labs_t <- bin(label_lst[[h6$test[i]]])
  knn_pur(ref, labs_r, tgt, labs_t, k=5)
})

h6$purity_sup  <- purity_sup
h6$purity_nat  <- purity_nat
h6$purity_loss <- purity_nat - purity_sup   # positive = sup is worse

# Correlation: condition number vs purity loss
cc_h6 <- cor(log(h6$kappa), h6$purity_loss, use="complete.obs", method="spearman")
cat(sprintf("  Spearman cor(log(kappa), purity_loss)=%.3f\n", cc_h6))

# Test B: re-run combat_sup on balanced test (50/50) for primary scenario
set.seed(42)
ltst0 <- sc$ltst
idx0  <- which(ltst0 == 0); idx1 <- which(ltst0 == 1)
n_bal <- min(length(idx0), length(idx1))
idx_b <- c(sample(idx0, n_bal), sample(idx1, n_bal))
tst_bal  <- sc$tst[, idx_b, drop=FALSE]
ltst_bal <- ltst0[idx_b]

res_nat_bal  <- run_combat_sup_nat(sc$dat, sc$bat, sc$lab, tst_bal)
res_nat_full <- run_combat_sup_nat(sc$dat, sc$bat, sc$lab, sc$tst)
res_sup_bal  <- run_combat_sup(sc$dat, sc$bat, sc$lab, tst_bal, ltst_bal)
res_sup_full <- run_combat_sup(sc$dat, sc$bat, sc$lab, sc$tst, sc$ltst)

pur_sup_natural  <- knn_pur(res_sup_full$ref,  sc$lab,  res_sup_full$tgt,  sc$ltst)
pur_sup_balanced <- knn_pur(res_sup_bal$ref,   sc$lab,  res_sup_bal$tgt,   ltst_bal)
pur_nat_natural  <- knn_pur(res_nat_full$ref,  sc$lab,  res_nat_full$tgt,  sc$ltst)
pur_nat_balanced <- knn_pur(res_nat_bal$ref,   sc$lab,  res_nat_bal$tgt,   ltst_bal)

cat(sprintf("  n3/USA combat_sup purity: natural=%.3f  balanced=%.3f (Δ=%.3f)\n",
            pur_sup_natural, pur_sup_balanced, pur_sup_balanced - pur_sup_natural))
cat(sprintf("  n3/USA combat_sup_nat purity: natural=%.3f  balanced=%.3f (Δ=%.3f)\n",
            pur_nat_natural, pur_nat_balanced, pur_nat_balanced - pur_nat_natural))
cat(sprintf("  Natural test class proportion (Active): %.2f\n", mean(sc$ltst)))

h6_pass <- (pur_sup_balanced - pur_sup_natural > 0.10) && (cc_h6 > 0.3)
if (h6_pass) {
  add("H6", "TRUE",
      sprintf("Balanced test: purity %.3f→%.3f (+%.3f); kappa corr=%.3f",
              pur_sup_natural, pur_sup_balanced,
              pur_sup_balanced - pur_sup_natural, cc_h6))
} else {
  add("H6", "FALSE or PARTIAL",
      sprintf("Balanced Δpurity=%.3f; kappa corr=%.3f",
              pur_sup_balanced - pur_sup_natural, cc_h6))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H7: Mean-only step-2 cannot invert class ordering ──────────────────\n")
# Test: check that class-axis separation for combat_sup_nat and combat_sup_mo
# is positive in ALL 18 scenarios (from precomputed data).

h7_rows <- list()
for (m in c("combat_sup_nat", "combat_sup_mo")) {
  for (n in 2:5) {
    for (test in ALL_STUDIES) {
      ref <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("%s_n%d_test%s_reference.csv", m, n, test)),
        row.names=1)), error=function(e) NULL)
      tgt <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("%s_n%d_test%s_target.csv", m, n, test)),
        row.names=1)), error=function(e) NULL)
      if (is.null(ref) || is.null(tgt)) next
      labs_r <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
      labs_t <- bin(label_lst[[test]])
      sep    <- cls_sep_test(ref, labs_r, tgt, labs_t)
      h7_rows[[length(h7_rows)+1]] <- data.frame(method=m, n=n, test=test, sep=sep)
    }
  }
}
h7 <- do.call(rbind, h7_rows)
n_neg   <- sum(h7$sep <= 0)
n_total <- nrow(h7)
cat(sprintf("  Negative class-axis separations: %d / %d\n", n_neg, n_total))
if (n_neg > 0) cat("  Negative cases:\n")
for (i in which(h7$sep <= 0))
  cat(sprintf("    %s n=%d test=%s sep=%.4f\n",
              h7$method[i], h7$n[i], h7$test[i], h7$sep[i]))

# Also verify mathematically: mean.only step applies a constant shift
# → cannot change ordering. Check that the shift is indeed constant across samples.
step2_shifts <- {
  comb  <- cbind(tr_nat, sc$tst)
  cb    <- factor(c(rep(1, ncol(tr_nat)), rep(2, ncol(sc$tst))))
  out   <- suppressMessages(
    ComBat(comb, batch=cb, mod=NULL, ref.batch=1L, mean.only=TRUE))
  tgt_out <- out[, (ncol(tr_nat)+1):ncol(out), drop=FALSE]
  shift   <- tgt_out - sc$tst[rownames(tgt_out), ]
  # Check: does shift vary across samples?
  col_sds <- apply(shift, 1, sd)
  median(col_sds)
}
cat(sprintf("  Median per-gene SD of sample-to-sample shift variation: %.6f\n",
            step2_shifts))

if (n_neg == 0 && step2_shifts < 0.001) {
  add("H7", "TRUE",
      sprintf("All %d mean-only method scores positive; shift is constant (SD=%.6f)",
              n_total, step2_shifts))
} else {
  add("H7", "FALSE or PARTIAL",
      sprintf("%d/%d negative; shift SD=%.6f", n_neg, n_total, step2_shifts))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H8: Supervised step-1 distorts gene-gene covariance vs test ─────────\n")
# Test: L2 distance of per-gene SD profiles:
# dist(sup_train, test_raw) vs dist(unsup_train, test_raw) vs dist(nat_train, test_raw)
# H8 TRUE: sup_train further from test than unsup_train and nat_train.

cg8    <- intersect(rownames(tr_sup), rownames(sc$tst))
sd_sup  <- rowSds(tr_sup[cg8, ])
sd_uns  <- rowSds(tr_unsup[cg8, ])
sd_nat  <- rowSds(tr_nat[cg8, ])
sd_mo   <- rowSds(tr_mo_cls[cg8, ])
sd_tst  <- rowSds(sc$tst[cg8, ])

d_sup <- sqrt(mean((sd_sup - sd_tst)^2))
d_uns <- sqrt(mean((sd_uns - sd_tst)^2))
d_nat <- sqrt(mean((sd_nat - sd_tst)^2))
d_mo  <- sqrt(mean((sd_mo  - sd_tst)^2))

cat(sprintf("  RMSE of per-gene SD vs raw test:\n"))
cat(sprintf("    unsup=%.4f  sup=%.4f  nat(mean-only)=%.4f  mo(mean-only)=%.4f\n",
            d_uns, d_sup, d_nat, d_mo))

# Also test variance covariance structure via first PC variance explained
pca_sup <- prcomp(t(tr_sup[cg8, ]), center=TRUE, scale.=FALSE, rank.=10)
pca_uns <- prcomp(t(tr_unsup[cg8, ]), center=TRUE, scale.=FALSE, rank.=10)
pca_nat <- prcomp(t(tr_nat[cg8, ]), center=TRUE, scale.=FALSE, rank.=10)
pca_tst <- prcomp(t(sc$tst[cg8, ]), center=TRUE, scale.=FALSE, rank.=10)

pct_sup <- pca_sup$sdev[1]^2 / sum(pca_sup$sdev^2)
pct_uns <- pca_uns$sdev[1]^2 / sum(pca_uns$sdev^2)
pct_nat <- pca_nat$sdev[1]^2 / sum(pca_nat$sdev^2)
pct_tst <- pca_tst$sdev[1]^2 / sum(pca_tst$sdev^2)

cat(sprintf("  PC1 variance explained: test=%.3f unsup=%.3f sup=%.3f nat=%.3f\n",
            pct_tst, pct_uns, pct_sup, pct_nat))

if (d_sup > d_uns * 1.2 && d_sup > d_nat * 1.2) {
  add("H8", "TRUE",
      sprintf("SD profile dist: sup=%.4f > unsup=%.4f, nat=%.4f", d_sup, d_uns, d_nat))
} else {
  add("H8", "FALSE or PARTIAL",
      sprintf("SD profile dist: sup=%.4f, unsup=%.4f, nat=%.4f", d_sup, d_uns, d_nat))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H9: Mean-only with class protection preserves class sep; without loses it ─\n")
# Test: compare class separation in training output for:
# (a) mean-only WITH class protection (tr_mo_cls)
# (b) mean-only WITHOUT class protection (tr_mo_raw)
# Across confounded (n=3/USA) and report.
# H9: sep_with > sep_without, especially in confounded scenarios.

h9_rows <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    sc_i  <- prep(n, test)
    labs_r <- sc_i$lab
    suppressMessages({
      t_with    <- ComBat(sc_i$dat, batch=sc_i$bat, mod=model.matrix(~labs_r), mean.only=TRUE)
      t_without <- ComBat(sc_i$dat, batch=sc_i$bat, mod=NULL, mean.only=TRUE)
    })
    sep_with    <- { w <- rowMeans(t_with[,labs_r==1,drop=FALSE]) - rowMeans(t_with[,labs_r==0,drop=FALSE])
                     w <- w/sqrt(sum(w^2)); s <- as.numeric(t(t_with)%*%w)
                     mean(s[labs_r==1])-mean(s[labs_r==0]) }
    sep_without <- { w <- rowMeans(t_without[,labs_r==1,drop=FALSE]) - rowMeans(t_without[,labs_r==0,drop=FALSE])
                     w <- w/sqrt(sum(w^2)); s <- as.numeric(t(t_without)%*%w)
                     mean(s[labs_r==1])-mean(s[labs_r==0]) }
    h9_rows[[length(h9_rows)+1]] <- data.frame(n=n, test=test,
      sep_with=sep_with, sep_without=sep_without, ratio=sep_with/max(sep_without,0.01))
  }
}
h9 <- do.call(rbind, h9_rows)
cat("  Training class separation — with vs without class protection:\n")
print(h9, row.names=FALSE, digits=3)
wt_h9 <- wilcox.test(h9$sep_with, h9$sep_without, paired=TRUE, exact=FALSE)
cat(sprintf("  Wilcoxon paired test: p=%.3g\n", wt_h9$p.value))

if (median(h9$ratio) > 1.5 && wt_h9$p.value < 0.05) {
  add("H9", "TRUE",
      sprintf("Class sep with protection: median=%.2f; without: %.2f (ratio=%.2f, p=%.3g)",
              median(h9$sep_with), median(h9$sep_without), median(h9$ratio), wt_h9$p.value))
} else {
  add("H9", "FALSE or PARTIAL",
      sprintf("Ratio=%.2f (p=%.3g)", median(h9$ratio), wt_h9$p.value))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H10: limma and ComBat mean-only step-1 are near-identical ───────────\n")
# Test: L2 norm of difference between nat and mo corrected matrices,
# and Pearson correlation, across all scenarios.

h10_rows <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    ref_nat <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_reference.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    ref_mo  <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_mo_n%d_test%s_reference.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    if (is.null(ref_nat) || is.null(ref_mo)) next
    cg <- intersect(rownames(ref_nat), rownames(ref_mo))
    # Align columns if needed
    sc_cols <- intersect(colnames(ref_nat), colnames(ref_mo))
    if (length(sc_cols) == 0) next
    diff_mat <- ref_nat[cg, sc_cols] - ref_mo[cg, sc_cols]
    l2_rel   <- sqrt(mean(diff_mat^2)) / sqrt(mean(ref_nat[cg, sc_cols]^2))
    r        <- cor(as.vector(ref_nat[cg, sc_cols]), as.vector(ref_mo[cg, sc_cols]))
    h10_rows[[length(h10_rows)+1]] <- data.frame(n=n, test=test, l2_rel=l2_rel, r=r)
  }
}
h10 <- do.call(rbind, h10_rows)
cat(sprintf("  Relative L2 difference nat vs mo: median=%.5f [%.5f, %.5f]\n",
            median(h10$l2_rel, na.rm=TRUE),
            quantile(h10$l2_rel, 0.25, na.rm=TRUE),
            quantile(h10$l2_rel, 0.75, na.rm=TRUE)))
cat(sprintf("  Pearson r between nat and mo matrices: median=%.6f\n",
            median(h10$r, na.rm=TRUE)))

if (median(h10$l2_rel, na.rm=TRUE) < 0.02 && median(h10$r, na.rm=TRUE) > 0.999) {
  add("H10", "TRUE",
      sprintf("Relative L2=%.5f, r=%.6f — near-identical outputs",
              median(h10$l2_rel, na.rm=TRUE), median(h10$r, na.rm=TRUE)))
} else {
  add("H10", "FALSE",
      sprintf("Relative L2=%.5f, r=%.6f — meaningful differences",
              median(h10$l2_rel, na.rm=TRUE), median(h10$r, na.rm=TRUE)))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H11: k=5 gap (combat vs nat) correlates with variance heterogeneity ─\n")
# Test: Spearman correlation between cross-batch variance heterogeneity
# (CV of per-batch median per-gene SD) and (combat purity - nat purity) at k=5.

h11_rows <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    sc_i <- prep(n, test)
    # Variance heterogeneity: CV of per-batch per-gene SD
    batch_sds <- sapply(unique(sc_i$bat), function(b) {
      median(rowSds(sc_i$dat[, sc_i$bat == b, drop=FALSE]), na.rm=TRUE)
    })
    cv_het <- sd(batch_sds) / mean(batch_sds)

    ref_c <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_n%d_test%s_reference.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    tgt_c <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_n%d_test%s_target.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    ref_nat <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_reference.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    tgt_nat <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_target.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    if (any(sapply(list(ref_c, tgt_c, ref_nat, tgt_nat), is.null))) next

    labs_r <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
    labs_t <- bin(label_lst[[test]])
    pur_c   <- knn_pur(ref_c,   labs_r, tgt_c,   labs_t)
    pur_nat <- knn_pur(ref_nat, labs_r, tgt_nat, labs_t)
    gap     <- pur_c - pur_nat   # positive = combat is better

    h11_rows[[length(h11_rows)+1]] <- data.frame(n=n, test=test,
                                                  cv_het=cv_het, gap=gap)
  }
}
h11 <- do.call(rbind, h11_rows)
cor_h11 <- cor(h11$cv_het, h11$gap, method="spearman", use="complete.obs")
cat("  Cross-batch variance heterogeneity vs purity gap (combat - nat) at k=5:\n")
print(h11, row.names=FALSE, digits=3)
cat(sprintf("  Spearman rho=%.3f\n", cor_h11))

if (cor_h11 > 0.35) {
  add("H11", "TRUE",
      sprintf("Spearman rho=%.3f — variance heterogeneity predicts gap", cor_h11))
} else {
  add("H11", "FALSE or PARTIAL",
      sprintf("Spearman rho=%.3f — no strong correlation", cor_h11))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H12: Larger class axis separation predicts k=11 advantage ───────────\n")
# Test: for each scenario × method, correlate training class-axis separation
# with purity at k=11 minus k=5.

h12_rows <- list()
for (m in c("combat", "combat_sup_nat", "combat_sup_mo")) {
  for (n in 2:5) {
    for (test in ALL_STUDIES) {
      ref <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("%s_n%d_test%s_reference.csv", m, n, test)),
        row.names=1)), error=function(e) NULL)
      tgt <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("%s_n%d_test%s_target.csv", m, n, test)),
        row.names=1)), error=function(e) NULL)
      if (is.null(ref) || is.null(tgt)) next
      labs_r <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
      labs_t <- bin(label_lst[[test]])
      sep    <- cls_sep_test(ref, labs_r, tgt, labs_t)
      pur5   <- knn_pur(ref, labs_r, tgt, labs_t, k=5)
      pur11  <- knn_pur(ref, labs_r, tgt, labs_t, k=11)
      h12_rows[[length(h12_rows)+1]] <- data.frame(
        method=m, n=n, test=test, sep=sep, pur5=pur5, pur11=pur11,
        delta_k=pur11 - pur5)
    }
  }
}
h12 <- do.call(rbind, h12_rows)
cor_h12 <- cor(h12$sep, h12$delta_k, method="spearman", use="complete.obs")
cat(sprintf("  Spearman rho(class-sep, purity[k=11]-purity[k=5])=%.3f\n", cor_h12))

# Also test: nat methods outperform combat at k=11 more when class sep is larger
h12_nat <- h12[h12$method == "combat_sup_nat", ]
h12_com <- h12[h12$method == "combat", ]
key <- paste(h12_nat$n, h12_nat$test)
key2 <- paste(h12_com$n, h12_com$test)
idx  <- match(key, key2)
h12_comp <- data.frame(
  n=h12_nat$n, test=h12_nat$test,
  sep_nat = h12_nat$sep,
  p11_nat = h12_nat$pur11,
  p11_com = h12_com$pur11[idx],
  adv_nat = h12_nat$pur11 - h12_com$pur11[idx])
cor_h12b <- cor(h12_comp$sep_nat, h12_comp$adv_nat, method="spearman")
cat(sprintf("  Spearman rho(nat class-sep, nat advantage at k=11)=%.3f\n", cor_h12b))

if (cor_h12 > 0.30 || cor_h12b > 0.35) {
  add("H12", "TRUE",
      sprintf("sep→(k11-k5) rho=%.3f; sep→nat_adv_k11 rho=%.3f", cor_h12, cor_h12b))
} else {
  add("H12", "FALSE or PARTIAL",
      sprintf("sep→(k11-k5) rho=%.3f; sep→nat_adv_k11 rho=%.3f", cor_h12, cor_h12b))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H13: v2 degradation is monotone in n ────────────────────────────────\n")
# Test: (nat_v1 - nat_v2) purity at k=5, grouped by n.
# H13 TRUE: degradation increases with n (Kruskal-Wallis + trend test).

h13_rows <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    for (v in c("nat","mo")) {
      r1 <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("combat_sup_%s_n%d_test%s_reference.csv", v, n, test)),
        row.names=1)), error=function(e) NULL)
      t1 <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("combat_sup_%s_n%d_test%s_target.csv", v, n, test)),
        row.names=1)), error=function(e) NULL)
      r2 <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("combat_sup_%s_v2_n%d_test%s_reference.csv", v, n, test)),
        row.names=1)), error=function(e) NULL)
      t2 <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("combat_sup_%s_v2_n%d_test%s_target.csv", v, n, test)),
        row.names=1)), error=function(e) NULL)
      if (any(sapply(list(r1,t1,r2,t2), is.null))) next
      labs_r <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
      labs_t <- bin(label_lst[[test]])
      p1 <- knn_pur(r1, labs_r, t1, labs_t)
      p2 <- knn_pur(r2, labs_r, t2, labs_t)
      h13_rows[[length(h13_rows)+1]] <- data.frame(n=n, test=test, variant=v,
                                                    v1=p1, v2=p2, delta=p1-p2)
    }
  }
}
h13 <- do.call(rbind, h13_rows)
cat("  Mean (v1 - v2) purity by n:\n")
h13_summ <- h13 %>% group_by(n) %>%
  summarise(mean_delta=mean(delta, na.rm=TRUE), n_obs=n(), .groups="drop")
print(h13_summ, digits=3)

# Jonckheere-Terpstra trend test (or just Spearman)
cor_h13 <- cor(h13$n, h13$delta, method="spearman")
cat(sprintf("  Spearman rho(n, v1-v2 delta)=%.3f\n", cor_h13))

if (cor_h13 > 0.40 && all(h13_summ$mean_delta >= 0)) {
  add("H13", "TRUE",
      sprintf("Spearman rho=%.3f; mean deltas: n2=%.3f n3=%.3f n4=%.3f n5=%.3f",
              cor_h13, h13_summ$mean_delta[1], h13_summ$mean_delta[2],
              h13_summ$mean_delta[3], h13_summ$mean_delta[4]))
} else {
  add("H13", "FALSE or PARTIAL",
      sprintf("Spearman rho=%.3f", cor_h13))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n── H14: Africa underperforms due to per-gene variance-scale difference ─\n")
# Test A: Compare per-gene SD distribution of Africa vs other studies (KS test).
# Test B: Pre-normalise all studies to equal per-gene SD, re-run combat_sup_nat,
#         check if Africa purity gap closes.

cg14 <- Reduce(intersect, lapply(ALL_STUDIES, function(s) rownames(dat_lst[[s]])))

study_sd <- sapply(ALL_STUDIES, function(s) {
  d <- dat_lst[[s]][cg14, ]
  if (max(d, na.rm=TRUE) > 100) d <- log2(d + 1)
  median(rowSds(d), na.rm=TRUE)
})
cat("  Median per-gene SD per study:\n")
print(round(study_sd, 4))

africa_sd <- rowSds(log2(dat_lst[["Africa"]][cg14, ] + 1))
other_sds <- unlist(lapply(setdiff(ALL_STUDIES, "Africa"), function(s) {
  d <- dat_lst[[s]][cg14, ]
  if (max(d, na.rm=TRUE) > 100) d <- log2(d+1)
  rowSds(d)
}))
ks_h14 <- ks.test(africa_sd, other_sds)
cat(sprintf("  KS test Africa vs others per-gene SD: D=%.4f, p=%.3g\n",
            ks_h14$statistic, ks_h14$p.value))

# Test B: normalise each study's variance to its global median SD, then re-run
# for the most discriminating scenario: n=5, test=Africa
sc_af <- prep(5, "Africa")
# Scale each batch to have the same median per-gene SD as the reference batch
ref_batch_sd <- median(rowSds(sc_af$dat[, sc_af$bat == ALL_STUDIES[1], drop=FALSE]), na.rm=TRUE)
dat_scaled <- sc_af$dat
for (b in unique(sc_af$bat)) {
  idx_b  <- which(sc_af$bat == b)
  b_sd   <- median(rowSds(sc_af$dat[, idx_b, drop=FALSE]), na.rm=TRUE)
  dat_scaled[, idx_b] <- sc_af$dat[, idx_b] * (ref_batch_sd / max(b_sd, 1e-6))
}
tst_scaled <- sc_af$tst * (ref_batch_sd / max(median(rowSds(sc_af$tst), na.rm=TRUE), 1e-6))

# Run nat on original and scaled
res_orig   <- run_combat_sup_nat(sc_af$dat, sc_af$bat, sc_af$lab, sc_af$tst)
res_scaled <- run_combat_sup_nat(dat_scaled, sc_af$bat, sc_af$lab, tst_scaled)
res_combat <- suppressMessages({
  dat_c <- ComBat(sc_af$dat, batch=sc_af$bat, mod=NULL)
  comb  <- cbind(dat_c, sc_af$tst)
  cb    <- c(rep(1L, ncol(dat_c)), rep(2L, ncol(sc_af$tst)))
  out   <- ComBat(comb, batch=cb, mod=NULL, ref.batch=1L)
  list(ref=dat_c, tgt=out[,(ncol(dat_c)+1):ncol(out),drop=FALSE])
})

labs5_r <- do.call(c, lapply(train_studies(5, "Africa"), function(s) bin(label_lst[[s]])))
labs5_t <- bin(label_lst[["Africa"]])
pur_nat_orig   <- knn_pur(res_orig$ref,   labs5_r, res_orig$tgt,   labs5_t)
pur_nat_scaled <- knn_pur(res_scaled$ref, labs5_r, res_scaled$tgt, labs5_t)
pur_combat     <- knn_pur(res_combat$ref, labs5_r, res_combat$tgt, labs5_t)

cat(sprintf("  n5/Africa KNN purity: combat=%.3f nat=%.3f nat_var-normalised=%.3f\n",
            pur_combat, pur_nat_orig, pur_nat_scaled))
gap_before <- pur_combat - pur_nat_orig
gap_after  <- pur_combat - pur_nat_scaled
cat(sprintf("  Gap (combat - nat): before normalisation=%.3f, after=%.3f\n",
            gap_before, gap_after))

h14_pass <- ks_h14$p.value < 0.001 && gap_after < gap_before * 0.5
if (h14_pass) {
  add("H14", "TRUE",
      sprintf("Africa SD differs from others (KS p=%.3g); variance normalisation closes gap %.3f→%.3f",
              ks_h14$p.value, gap_before, gap_after))
} else {
  add("H14", "FALSE or PARTIAL",
      sprintf("KS p=%.3g; gap before=%.3f after=%.3f", ks_h14$p.value, gap_before, gap_after))
}

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n\n═══════════════════════════════════════════════════════════════════\n")
cat("  SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n")
print(VT, row.names=FALSE)
write.csv(VT, file.path(OUT_DIR, "verdict_table.csv"), row.names=FALSE)
cat(sprintf("\nVerdict table written to %s/verdict_table.csv\n", OUT_DIR))
