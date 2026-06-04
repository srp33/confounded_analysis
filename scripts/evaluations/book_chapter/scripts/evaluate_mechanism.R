#!/usr/bin/env Rscript
# evaluate_mechanism.R
# Causal mechanism tests for the combat_sup failure story:
#   M1 (H8-causal): PC1 of supervised training aligns with class axis;
#                   PC1 concentration cross-scenario predicts purity loss.
#   M2 (H2->KNN):   delta.hat predicts test compression;
#                   compression predicts purity loss.
#
# Both tests use precomputed reference/target matrices where possible to
# avoid re-running full ComBat pipelines.

suppressMessages({
  library(dplyr)
  library(sva)
  library(limma)
  library(matrixStats)
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
  dat  <- log_safe(dat)
  tst  <- log_safe(tst)
  list(dat = dat, bat = bat, lab = lab, tst = tst, ltst = ltst)
}

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

# delta.hat: per-gene variance of standardized test-batch residuals under
# the step-2 design (batch + class for supervised, batch-only for unsupervised).
compute_delta_hat <- function(train_mat, train_lab, test_mat, test_lab,
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
  if (is.null(B)) return(rep(NA_real_, nrow(comb)))
  ref_idx <- seq_len(ncol(train_mat))
  vp      <- pmax(rowMeans((comb[, ref_idx] - t(design2[ref_idx, ] %*% B))^2), 1e-10)
  tst_idx <- (ncol(train_mat) + 1):ncol(comb)
  s_tst   <- (comb[, tst_idx] - t(design2[tst_idx, ] %*% B)) / sqrt(vp)
  rowVars(s_tst)
}

cat("═══════════════════════════════════════════════════════════════════\n")
cat("  MECHANISM EVALUATIONS\n")
cat("═══════════════════════════════════════════════════════════════════\n")

# ═══════════════════════════════════════════════════════════════════════════
cat("\n── M1 (H8-causal): PC1 alignment with class axis ─────────────────────\n")
cat("   Test A: |cos(PC1_sup, class_axis)| >> |cos(PC1_uns, class_axis)|\n")
cat("   Test B: cross-scenario rho(PC1_frac_sup, purity_loss) > 0\n\n")

m1_rows <- list()

for (n in 2:5) {
  for (test in ALL_STUDIES) {
    ref_sup <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_reference.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    tgt_sup <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_target.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    ref_nat <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_reference.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    tgt_nat <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_target.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    ref_uns <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_n%d_test%s_reference.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    if (any(sapply(list(ref_sup, tgt_sup, ref_nat, tgt_nat, ref_uns), is.null))) next

    labs_r <- do.call(c, lapply(train_studies(n, test),
                                function(s) bin(label_lst[[s]])))
    labs_t <- bin(label_lst[[test]])
    cg     <- rownames(ref_sup)

    # ── Class-axis unit vector from supervised reference
    w_sup <- rowMeans(ref_sup[cg, labs_r == 1, drop = FALSE]) -
             rowMeans(ref_sup[cg, labs_r == 0, drop = FALSE])
    w_sup <- w_sup / sqrt(sum(w_sup^2))

    # ── PC1 of supervised reference
    pca_sup   <- prcomp(t(ref_sup[cg, ]), center = TRUE, scale. = FALSE, rank. = 10)
    pc1_sup   <- pca_sup$rotation[, 1]
    cos_sup   <- abs(sum(pc1_sup[names(w_sup)] * w_sup))
    frac_sup  <- pca_sup$sdev[1]^2 / sum(pca_sup$sdev^2)

    # ── Class-axis unit vector from unsupervised reference
    w_uns <- rowMeans(ref_uns[cg, labs_r == 1, drop = FALSE]) -
             rowMeans(ref_uns[cg, labs_r == 0, drop = FALSE])
    w_uns <- w_uns / sqrt(sum(w_uns^2))

    # ── PC1 of unsupervised reference
    pca_uns   <- prcomp(t(ref_uns[cg, ]), center = TRUE, scale. = FALSE, rank. = 10)
    pc1_uns   <- pca_uns$rotation[, 1]
    cos_uns   <- abs(sum(pc1_uns[names(w_uns)] * w_uns))
    frac_uns  <- pca_uns$sdev[1]^2 / sum(pca_uns$sdev^2)

    # ── Purity
    pur_sup <- knn_pur(ref_sup, labs_r, tgt_sup, labs_t, k = 5)
    pur_nat <- knn_pur(ref_nat, labs_r, tgt_nat, labs_t, k = 5)

    m1_rows[[length(m1_rows) + 1]] <- data.frame(
      n = n, test = test,
      cos_sup  = cos_sup,  cos_uns  = cos_uns,
      frac_sup = frac_sup, frac_uns = frac_uns,
      pur_sup  = pur_sup,  pur_nat  = pur_nat,
      purity_loss = pur_nat - pur_sup
    )
  }
}

m1 <- do.call(rbind, m1_rows)

cat("  Per-scenario results:\n")
print(m1[, c("n","test","cos_sup","cos_uns","frac_sup","frac_uns","pur_sup","purity_loss")],
      row.names = FALSE, digits = 3)

med_cos_sup <- median(m1$cos_sup, na.rm = TRUE)
med_cos_uns <- median(m1$cos_uns, na.rm = TRUE)
med_frac_sup <- median(m1$frac_sup, na.rm = TRUE)
med_frac_uns <- median(m1$frac_uns, na.rm = TRUE)

cat(sprintf("\n  Median |cos(PC1, class_axis)|:  sup=%.3f  uns=%.3f\n",
            med_cos_sup, med_cos_uns))
cat(sprintf("  Median PC1 variance fraction:   sup=%.3f  uns=%.3f\n",
            med_frac_sup, med_frac_uns))

wt_cos <- wilcox.test(m1$cos_sup, m1$cos_uns, paired = TRUE, exact = FALSE)
cat(sprintf("  Paired Wilcoxon cos_sup vs cos_uns: p=%.3g\n", wt_cos$p.value))

# B: does cross-scenario PC1 fraction predict purity loss?
rho_frac_loss  <- cor(m1$frac_sup, m1$purity_loss, method = "spearman", use = "complete.obs")
rho_cos_loss   <- cor(m1$cos_sup,  m1$purity_loss, method = "spearman", use = "complete.obs")
rho_frac_pur   <- cor(m1$frac_sup, m1$pur_sup,     method = "spearman", use = "complete.obs")

cat(sprintf("\n  Spearman rho(PC1_frac_sup, purity_loss) = %.3f\n", rho_frac_loss))
cat(sprintf("  Spearman rho(cos_sup, purity_loss)       = %.3f\n", rho_cos_loss))
cat(sprintf("  Spearman rho(PC1_frac_sup, pur_sup)      = %.3f\n", rho_frac_pur))

# Save per-scenario table for plotting
write.csv(m1, file.path(OUT_DIR, "m1_pc1_alignment.csv"), row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════
cat("\n── M2 (H2->KNN): delta.hat -> test compression -> purity loss ─────────\n")
cat("   Test A: compression_sup << compression_nat (~1.0)\n")
cat("   Test B: rho(delta.hat, compression_sup) < 0\n")
cat("   Test C: rho(compression_sup, purity_sup) > 0 (more compression -> worse purity)\n")
cat("   Test D: mediation — direct delta.hat->purity weakens after conditioning on compression\n\n")

m2_rows <- list()

for (n in 2:5) {
  for (test in ALL_STUDIES) {
    sc_i <- prep_raw(n, test)

    ref_sup <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_reference.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    tgt_sup <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_target.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    ref_nat <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_reference.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    tgt_nat <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_target.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    if (any(sapply(list(ref_sup, tgt_sup, ref_nat, tgt_nat), is.null))) next

    labs_r <- do.call(c, lapply(train_studies(n, test),
                                function(s) bin(label_lst[[s]])))
    labs_t <- sc_i$ltst
    cg     <- intersect(rownames(ref_sup), rownames(sc_i$tst))

    # delta.hat for supervised step-2
    dh_sup <- compute_delta_hat(ref_sup[cg, ], labs_r, sc_i$tst[cg, ], labs_t,
                                supervised = TRUE)
    med_dh_sup <- median(dh_sup, na.rm = TRUE)

    # delta.hat for nat step-2 (batch-only design)
    dh_nat <- compute_delta_hat(ref_nat[cg, ], labs_r, sc_i$tst[cg, ], labs_t,
                                supervised = FALSE)
    med_dh_nat <- median(dh_nat, na.rm = TRUE)

    # Fraction of genes with delta.hat > 2 (heavy over-variance)
    frac_heavy_sup <- mean(dh_sup > 2, na.rm = TRUE)
    frac_heavy_nat <- mean(dh_nat > 2, na.rm = TRUE)

    # Compression: post-step2 SD / pre-step2 SD
    sd_before      <- rowSds(sc_i$tst[cg, ])
    sd_after_sup   <- rowSds(tgt_sup[cg, , drop = FALSE])
    sd_after_nat   <- rowSds(tgt_nat[cg, , drop = FALSE])

    comp_sup <- median(sd_after_sup / pmax(sd_before, 1e-10), na.rm = TRUE)
    comp_nat <- median(sd_after_nat / pmax(sd_before, 1e-10), na.rm = TRUE)

    # Theoretical prediction: if delta.hat drives compression,
    # compression_ratio ≈ 1 / sqrt(delta.hat)
    pred_comp <- 1 / sqrt(med_dh_sup)

    # KNN purity from precomputed matrices
    pur_sup <- knn_pur(ref_sup, labs_r, tgt_sup, labs_t, k = 5)
    pur_nat <- knn_pur(ref_nat, labs_r, tgt_nat, labs_t, k = 5)

    m2_rows[[length(m2_rows) + 1]] <- data.frame(
      n = n, test = test,
      med_dh_sup = med_dh_sup, med_dh_nat = med_dh_nat,
      frac_heavy_sup = frac_heavy_sup, frac_heavy_nat = frac_heavy_nat,
      comp_sup = comp_sup, comp_nat = comp_nat, pred_comp = pred_comp,
      pur_sup = pur_sup, pur_nat = pur_nat,
      purity_loss = pur_nat - pur_sup
    )
  }
}

m2 <- do.call(rbind, m2_rows)

cat("  Per-scenario results:\n")
print(m2[, c("n","test","med_dh_sup","med_dh_nat","comp_sup","comp_nat","pred_comp","pur_sup","purity_loss")],
      row.names = FALSE, digits = 3)

# Theoretical prediction check: comp_sup vs 1/sqrt(delta.hat)
rho_pred <- cor(m2$comp_sup, m2$pred_comp, method = "spearman", use = "complete.obs")
cat(sprintf("\n  Predicted vs observed compression (1/sqrt(dh) vs comp_sup): rho=%.3f\n",
            rho_pred))

# Chain correlations
rho_dh_comp   <- cor(m2$med_dh_sup, m2$comp_sup,  method = "spearman", use = "complete.obs")
rho_comp_pur  <- cor(m2$comp_sup,   m2$pur_sup,    method = "spearman", use = "complete.obs")
rho_dh_pur    <- cor(m2$med_dh_sup, m2$pur_sup,    method = "spearman", use = "complete.obs")
rho_dh_loss   <- cor(m2$med_dh_sup, m2$purity_loss, method = "spearman", use = "complete.obs")

cat(sprintf("\n  Chain A: rho(delta.hat_sup, comp_sup) = %.3f  (>0 means more delta = less compression)\n",
            rho_dh_comp))
cat(sprintf("  Chain B: rho(comp_sup, pur_sup)        = %.3f  (>0 means less compressed = better purity)\n",
            rho_comp_pur))
cat(sprintf("  Chain C: rho(delta.hat_sup, pur_sup)   = %.3f  (direct link, expect negative)\n",
            rho_dh_pur))
cat(sprintf("  Chain D: rho(delta.hat_sup, loss)      = %.3f  (expect positive: more dh = more loss)\n",
            rho_dh_loss))

# Mediation test: partial correlation of delta.hat and purity, controlling for compression
# r(X,Y | Z) = (r_xy - r_xz*r_yz) / sqrt((1-r_xz^2)*(1-r_yz^2))
r_xy <- cor(m2$med_dh_sup, m2$pur_sup, use = "complete.obs")   # delta.hat, purity
r_xz <- cor(m2$med_dh_sup, m2$comp_sup, use = "complete.obs")  # delta.hat, compression
r_yz <- cor(m2$comp_sup,   m2$pur_sup,  use = "complete.obs")  # compression, purity
partial_dh_pur_given_comp <-
  (r_xy - r_xz * r_yz) / sqrt(pmax((1 - r_xz^2) * (1 - r_yz^2), 1e-10))

cat(sprintf("\n  Mediation: partial r(delta.hat, purity | compression) = %.3f\n",
            partial_dh_pur_given_comp))
cat("  (If compression fully mediates, this should be near 0)\n")

# Summary compression comparison
cat(sprintf("\n  Median compression: sup=%.3f  nat=%.3f\n",
            median(m2$comp_sup, na.rm = TRUE),
            median(m2$comp_nat, na.rm = TRUE)))
cat(sprintf("  Fraction heavily over-dispersed (dh>2): sup=%.3f  nat=%.3f\n",
            median(m2$frac_heavy_sup, na.rm = TRUE),
            median(m2$frac_heavy_nat, na.rm = TRUE)))

write.csv(m2, file.path(OUT_DIR, "m2_delta_compression.csv"), row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════
cat("\n\n═══════════════════════════════════════════════════════════════════\n")
cat("  MECHANISM SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("M1 (PC1 alignment with class axis):\n")
cat(sprintf("  Supervised   |cos(PC1, class_axis)| = %.3f\n", med_cos_sup))
cat(sprintf("  Unsupervised |cos(PC1, class_axis)| = %.3f\n", med_cos_uns))
cat(sprintf("  Paired Wilcoxon p = %.3g\n", wt_cos$p.value))
cat(sprintf("  PC1 variance fraction: sup=%.3f, uns=%.3f\n", med_frac_sup, med_frac_uns))
cat(sprintf("  rho(PC1_frac, purity_loss) = %.3f  [> 0 supports mechanism]\n", rho_frac_loss))
cat(sprintf("  rho(cos_sup,  purity_loss) = %.3f  [> 0 supports mechanism]\n", rho_cos_loss))

cat("\nM2 (delta.hat chain):\n")
cat(sprintf("  rho(dh, compression)   = %.3f  [< 0: more dh -> more shrinkage]\n", rho_dh_comp))
cat(sprintf("  rho(compression, pur)  = %.3f  [> 0: less shrinkage -> better purity]\n", rho_comp_pur))
cat(sprintf("  rho(dh, purity)        = %.3f  [< 0: direct chain]\n", rho_dh_pur))
cat(sprintf("  partial r(dh, pur|comp)= %.3f  [near 0: mediated by compression]\n",
            partial_dh_pur_given_comp))

cat("\n")
if (med_cos_sup > 0.7 && wt_cos$p.value < 0.05 &&
    rho_frac_loss > 0.3 && rho_dh_comp < -0.3 && rho_comp_pur > 0.3) {
  cat("VERDICT: MECHANISM CONFIRMED\n")
  cat("  PC1 collapse onto class axis is real and predicts purity loss;\n")
  cat("  delta.hat chain to test compression is empirically supported.\n")
} else {
  cat("VERDICT: PARTIAL — inspect per-chain correlations above\n")
}
