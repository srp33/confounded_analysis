#!/usr/bin/env Rscript
# evaluate_mechanism_v2.R
# Efficient mechanism tests using cached reference matrices and
# precomputed purity from m1_pc1_alignment.csv.
#
# Key optimization: reference matrices are IDENTICAL within each n
# (step-1 is test-independent), so load once per n (4 files, not 18).
# Compression is derived theoretically from delta.hat + var_pooled,
# so target CSV files are never loaded.

suppressMessages({
  library(sva)
  library(limma)
  library(matrixStats)
})
source("scripts/ComBat_nat.R")

DATA_FILE <- "data/TB_real_data.RData"
ADJ_DIR   <- "outputs/adjusted_data/all_scenarios"
OUT_DIR   <- "outputs/diagnostics/hypothesis_tests"
load(DATA_FILE)

ALL_STUDIES <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1", 1, "Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES != test][seq_len(n)]

log_safe <- function(m) {
  if (max(m, na.rm = TRUE) > 100) {
    mn <- min(m, na.rm = TRUE)
    if (mn < 0) m <- m - mn
    m <- log2(m + 1)
  }
  m
}

prep_raw <- function(n, test) {
  refs <- train_studies(n, test)
  cg   <- Reduce(intersect, lapply(c(refs, test), function(s) rownames(dat_lst[[s]])))
  dat  <- do.call(cbind, lapply(refs, function(s) dat_lst[[s]][cg, ]))
  bat  <- do.call(c, lapply(refs, function(s) rep(s, ncol(dat_lst[[s]]))))
  lab  <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))
  tst  <- dat_lst[[test]][cg, ]
  ltst <- bin(label_lst[[test]])
  list(dat = log_safe(dat), bat = bat, lab = lab,
       tst = log_safe(tst), ltst = ltst)
}

# Returns both delta.hat (per-gene test-batch variance in standardized space)
# and var_pooled (per-gene pooled residual variance from training).
# These two quantities fully determine the ComBat step-2 variance correction:
#   SD_after ≈ sqrt(var_pooled / delta_hat)
compute_dh_vp <- function(train_mat, train_lab, test_mat, test_lab,
                           supervised = TRUE) {
  cg      <- intersect(rownames(train_mat), rownames(test_mat))
  comb    <- cbind(train_mat[cg, ], test_mat[cg, ])
  bat2    <- factor(c(rep(1, ncol(train_mat)), rep(2, ncol(test_mat))))
  if (supervised) {
    cls2    <- c(train_lab, test_lab)
    design2 <- model.matrix(~0 + bat2 + cls2)
  } else {
    design2 <- model.matrix(~0 + bat2)
  }
  B <- tryCatch(
    solve(crossprod(design2), tcrossprod(t(design2), comb)),
    error = function(e) NULL)
  if (is.null(B)) return(list(dh = rep(NA_real_, nrow(comb)),
                               vp = rep(NA_real_, nrow(comb))))
  ref_idx <- seq_len(ncol(train_mat))
  tst_idx <- (ncol(train_mat) + 1):ncol(comb)
  resid_r <- comb[, ref_idx] - t(design2[ref_idx, ] %*% B)
  vp      <- pmax(rowMeans(resid_r^2), 1e-10)
  s_tst   <- (comb[, tst_idx] - t(design2[tst_idx, ] %*% B)) / sqrt(vp)
  list(dh = rowVars(s_tst), vp = vp)
}

cat("═══════════════════════════════════════════════════════════════════\n")
cat("  MECHANISM EVALUATIONS (v2)\n")
cat("═══════════════════════════════════════════════════════════════════\n")

# ── M1: summarise pre-computed results ────────────────────────────────────
cat("\n── M1 (H8-causal): PC1 alignment summary ──────────────────────────────\n")

m1 <- read.csv(file.path(OUT_DIR, "m1_pc1_alignment.csv"))

cat("  Per-scenario results:\n")
print(m1[, c("n","test","cos_sup","cos_uns","frac_sup","frac_uns","pur_sup","purity_loss")],
      row.names = FALSE, digits = 3)

med_cos_sup  <- median(m1$cos_sup,  na.rm = TRUE)
med_cos_uns  <- median(m1$cos_uns,  na.rm = TRUE)
med_frac_sup <- median(m1$frac_sup, na.rm = TRUE)
med_frac_uns <- median(m1$frac_uns, na.rm = TRUE)

wt_cos <- wilcox.test(m1$cos_sup, m1$cos_uns, paired = TRUE, exact = FALSE)

cat(sprintf("\n  Median |cos(PC1, class_axis)|:  sup=%.3f  uns=%.3f  (paired Wilcoxon p=%.3g)\n",
            med_cos_sup, med_cos_uns, wt_cos$p.value))
cat(sprintf("  Median PC1 variance fraction:   sup=%.3f  uns=%.3f\n",
            med_frac_sup, med_frac_uns))

rho_frac_loss <- cor(m1$frac_sup, m1$purity_loss,
                     method = "spearman", use = "complete.obs")
rho_cos_loss  <- cor(m1$cos_sup,  m1$purity_loss,
                     method = "spearman", use = "complete.obs")
rho_frac_pur  <- cor(m1$frac_sup, m1$pur_sup,
                     method = "spearman", use = "complete.obs")

cat(sprintf("\n  Spearman rho(PC1_frac_sup, purity_loss) = %.3f  [+: more collapse -> more loss]\n",
            rho_frac_loss))
cat(sprintf("  Spearman rho(cos_sup,  purity_loss)     = %.3f  [+: more alignment -> more loss]\n",
            rho_cos_loss))
cat(sprintf("  Spearman rho(PC1_frac_sup, pur_sup)     = %.3f  [-: more collapse -> worse purity]\n",
            rho_frac_pur))

# ── M2: delta.hat -> compression -> purity ────────────────────────────────
cat("\n── M2 (H2->KNN): delta.hat -> test compression -> purity ──────────────\n")
cat("   Compression is THEORETICAL: SD_after ≈ sqrt(var_pooled / delta_hat)\n")
cat("   Mean-only step-2 (nat): compression = 1.0 exactly (no variance correction)\n")
cat("   Purity values reused from m1_pc1_alignment.csv\n\n")

m2_rows <- list()

for (n in 2:5) {
  for (test in ALL_STUDIES) {
    ref_sup <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_reference.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    if (is.null(ref_sup)) next

    sc_i  <- prep_raw(n, test)
    labs_r <- do.call(c, lapply(train_studies(n, test),
                                 function(s) bin(label_lst[[s]])))
    labs_t <- sc_i$ltst
    cg     <- intersect(rownames(ref_sup), rownames(sc_i$tst))

    # delta.hat and var_pooled from supervised step-2 design
    res_sup <- compute_dh_vp(ref_sup[cg, ], labs_r, sc_i$tst[cg, ], labs_t,
                              supervised = TRUE)
    # delta.hat from nat step-2 design (batch only)
    res_nat <- compute_dh_vp(ref_sup[cg, ], labs_r, sc_i$tst[cg, ], labs_t,
                              supervised = FALSE)

    dh_sup  <- res_sup$dh;   vp_sup <- res_sup$vp
    dh_nat  <- res_nat$dh
    sd_bef  <- rowSds(sc_i$tst[cg, ])

    # Theoretical compression ratio: SD_after_sup / SD_before
    # ComBat corrects: y_out = sqrt(vp) * (z - gamma) / sqrt(dh) + fitted_mean
    # => Var(y_out) ≈ vp / dh   => SD_out ≈ sqrt(vp/dh)
    # => compression = sqrt(vp/dh) / SD_before
    comp_sup <- median(sqrt(pmax(vp_sup, 1e-10) / pmax(dh_sup, 1e-10)) /
                       pmax(sd_bef, 1e-10), na.rm = TRUE)

    # nat is mean-only: Var unchanged => compression = 1.0
    comp_nat <- 1.0

    # Theoretical check: does 1/sqrt(dh) alone match comp × SD_before/sqrt(vp)?
    pred_simple <- median(1 / sqrt(pmax(dh_sup, 1e-10)), na.rm = TRUE)

    # Purity from M1 precomputed data
    m1_row <- m1[m1$n == n & m1$test == test, ]
    pur_sup <- if (nrow(m1_row) > 0) m1_row$pur_sup[1] else NA_real_
    pur_nat <- if (nrow(m1_row) > 0) m1_row$pur_nat[1] else NA_real_

    m2_rows[[length(m2_rows) + 1]] <- data.frame(
      n = n, test = test,
      med_dh_sup = median(dh_sup, na.rm = TRUE),
      med_dh_nat = median(dh_nat, na.rm = TRUE),
      frac_dh_gt2_sup = mean(dh_sup > 2, na.rm = TRUE),
      frac_dh_gt2_nat = mean(dh_nat > 2, na.rm = TRUE),
      comp_sup   = comp_sup,
      comp_nat   = comp_nat,
      pred_simple = pred_simple,
      pur_sup    = pur_sup,
      pur_nat    = pur_nat,
      purity_loss = pur_nat - pur_sup
    )
  }
}

m2 <- do.call(rbind, m2_rows)

cat("\n  Per-scenario results:\n")
print(m2[, c("n","test","med_dh_sup","med_dh_nat","frac_dh_gt2_sup",
             "comp_sup","pur_sup","purity_loss")],
      row.names = FALSE, digits = 3)

# Theoretical compression prediction vs delta.hat relationship
rho_pred <- cor(m2$comp_sup, 1/sqrt(m2$med_dh_sup),
                method = "spearman", use = "complete.obs")
cat(sprintf("\n  Predicted compression (1/sqrt(dh)) tracks theoretical: rho=%.3f\n",
            rho_pred))

# ── Chain correlations ─────────────────────────────────────────────────────
rho_dh_comp  <- cor(m2$med_dh_sup, m2$comp_sup,     method = "spearman", use = "complete.obs")
rho_comp_pur <- cor(m2$comp_sup,   m2$pur_sup,       method = "spearman", use = "complete.obs")
rho_dh_pur   <- cor(m2$med_dh_sup, m2$pur_sup,       method = "spearman", use = "complete.obs")
rho_dh_loss  <- cor(m2$med_dh_sup, m2$purity_loss,   method = "spearman", use = "complete.obs")

cat(sprintf("\n  Chain correlations:\n"))
cat(sprintf("  A: rho(delta.hat,   compression_sup) = %+.3f  [expect -: more dh -> more shrinkage]\n",
            rho_dh_comp))
cat(sprintf("  B: rho(compression, purity_sup)      = %+.3f  [expect +: less shrinkage -> better purity]\n",
            rho_comp_pur))
cat(sprintf("  C: rho(delta.hat,   purity_sup)      = %+.3f  [expect -: direct chain]\n",
            rho_dh_pur))
cat(sprintf("  D: rho(delta.hat,   purity_loss)     = %+.3f  [expect +: more dh -> more loss]\n",
            rho_dh_loss))

# Mediation: partial r(delta.hat, purity | compression)
# r_partial = (r_xy - r_xz*r_yz) / sqrt((1-r_xz^2)*(1-r_yz^2))
r_xy <- cor(m2$med_dh_sup, m2$pur_sup,  use = "complete.obs")
r_xz <- cor(m2$med_dh_sup, m2$comp_sup, use = "complete.obs")
r_yz <- cor(m2$comp_sup,   m2$pur_sup,  use = "complete.obs")
r_partial <- (r_xy - r_xz * r_yz) / sqrt(pmax((1 - r_xz^2) * (1 - r_yz^2), 1e-10))

cat(sprintf("\n  Mediation: partial r(delta.hat, purity | compression) = %.3f\n",
            r_partial))
cat("  (0 = compression fully mediates; unchanged = compression is irrelevant)\n")

cat(sprintf("\n  Summary: median delta.hat sup=%.2f, nat=%.2f\n",
            median(m2$med_dh_sup), median(m2$med_dh_nat)))
cat(sprintf("  Median compression sup=%.3f (nat=1.000 by definition)\n",
            median(m2$comp_sup)))
cat(sprintf("  Fraction genes with dh>2 (heavy compression): sup=%.3f, nat=%.3f\n",
            median(m2$frac_dh_gt2_sup), median(m2$frac_dh_gt2_nat)))

write.csv(m2, file.path(OUT_DIR, "m2_delta_compression.csv"), row.names = FALSE)

# ── Full mechanism summary ─────────────────────────────────────────────────
cat("\n\n═══════════════════════════════════════════════════════════════════\n")
cat("  MECHANISM SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("M1 — PC1 collapses onto class axis in supervised training:\n")
cat(sprintf("  |cos(PC1_sup, class_axis)| = %.3f  vs  |cos(PC1_uns)| = %.3f\n",
            med_cos_sup, med_cos_uns))
cat(sprintf("  PC1 fraction:  sup=%.3f  uns=%.3f  (Wilcoxon p=%.3g)\n",
            med_frac_sup, med_frac_uns, wt_cos$p.value))
cat(sprintf("  rho(PC1_frac, purity_loss) = %.3f\n", rho_frac_loss))
cat(sprintf("  rho(cos_sup,  purity_loss) = %.3f\n", rho_cos_loss))

cat("\nM2 — Elevated delta.hat compresses test variance:\n")
cat(sprintf("  rho(dh, comp)  = %.3f  rho(comp, pur) = %.3f  rho(dh, pur) = %.3f\n",
            rho_dh_comp, rho_comp_pur, rho_dh_pur))
cat(sprintf("  Partial r(dh, pur | comp) = %.3f\n", r_partial))

mechanism_supported <-
  med_cos_sup > 0.7 && med_cos_uns < med_cos_sup &&
  rho_frac_loss > 0.3 && rho_dh_comp < -0.3 && rho_comp_pur > 0.3

cat("\nOVERALL VERDICT:", if (mechanism_supported) "MECHANISM CONFIRMED" else "PARTIAL — see details above", "\n")
