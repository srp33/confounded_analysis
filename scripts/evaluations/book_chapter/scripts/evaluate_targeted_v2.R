#!/usr/bin/env Rscript
# evaluate_targeted_v2.R
# Follow-up on the T1 reversal and T3 n-dependence.
#
# T1 REVERSAL: high-SNR genes have LOWER dh_sup (2.41) not higher.
# Why? Supervised step-2 includes class in model -> class variance is removed
# from BOTH training AND test residuals by beta*class term. Residual delta.hat
# for high-SNR genes = Var(noise_test)/Var(noise_train) ≈ 1.
# Low-SNR genes have noise-dominated residuals; test noise > train noise
# (baseline inflation) -> higher delta.hat.
#
# T2 KEY FINDING: despite T1, class-axis separation still drops 65% after
# supervised step-2. Why? Even moderate compression (dh=2.41) substantially
# shrinks class-axis scores because test variance is CLASS-dominated:
#   compression_g = sqrt(vp_g/dh_g) / SD_before_g
#   For class-dominated genes: SD_before >> sqrt(vp_g), so compression << 1
#   even at moderate dh.
#
# T3: India purity 0.305->0.662 from n=3->n=4 despite dh staying ~3.9.
# Hypothesis: with n=4, more training batches -> lower PC1_frac -> KNN uses
# PC2-50 which still carry class signal. Measure: class signal in PC2-50.

suppressMessages({
  library(sva)
  library(limma)
  library(matrixStats)
  library(dplyr)
})
source("scripts/ComBat_nat.R")

DATA_FILE <- "data/TB_real_data.RData"
ADJ_DIR   <- "outputs/adjusted_data/all_scenarios"
OUT_DIR   <- "outputs/diagnostics/hypothesis_tests"
load(DATA_FILE)

ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES != test][seq_len(n)]
log_safe <- function(m) {
  if (max(m, na.rm=TRUE) > 100) {
    mn <- min(m, na.rm=TRUE); if (mn < 0) m <- m - mn; m <- log2(m+1)
  }; m
}

full_step2 <- function(train_mat, train_lab, test_mat, test_lab,
                       supervised=TRUE, mean_only=FALSE) {
  cg    <- intersect(rownames(train_mat), rownames(test_mat))
  comb  <- cbind(train_mat[cg,], test_mat[cg,])
  bat   <- factor(c(rep(1,ncol(train_mat)), rep(2,ncol(test_mat))))
  cls   <- c(train_lab, test_lab)
  des   <- if (supervised) model.matrix(~0+bat+cls) else model.matrix(~0+bat)
  B     <- tryCatch(solve(crossprod(des), tcrossprod(t(des),comb)), error=function(e) NULL)
  if (is.null(B)) return(NULL)
  ri    <- seq_len(ncol(train_mat))
  ti    <- (ncol(train_mat)+1):ncol(comb)
  vp    <- pmax(rowMeans((comb[,ri] - t(des[ri,]%*%B))^2), 1e-10)
  s_tst <- (comb[,ti] - t(des[ti,]%*%B)) / sqrt(vp)
  dh    <- rowVars(s_tst)
  gh    <- rowMeans(s_tst)
  if (mean_only) {
    corrected <- test_mat[cg,] - gh * sqrt(vp)
  } else {
    corrected <- sqrt(vp) * (s_tst - gh) / sqrt(pmax(dh,1e-10)) +
                 t(des[ti,] %*% B)
  }
  list(dh=dh, gh=gh, vp=vp, corrected=corrected, cg=cg)
}

knn_pur <- function(ref_mat, ref_lab, tst_mat, tst_lab, k=5, n_pcs=50) {
  cg  <- intersect(rownames(ref_mat), rownames(tst_mat))
  pca <- prcomp(t(ref_mat[cg,]), center=TRUE, scale.=FALSE, rank.=n_pcs)
  rr  <- pca$x
  tr  <- predict(pca, t(tst_mat[cg,]))
  d   <- as.matrix(dist(rbind(tr,rr)))[seq_len(nrow(tr)),(nrow(tr)+1):(nrow(tr)+nrow(rr))]
  mean(sapply(seq_len(nrow(tr)), function(i)
    mean(ref_lab[order(d[i,])[seq_len(k)]] == tst_lab[i])))
}

sep_fn <- function(s, l) mean(s[l==1]) - mean(s[l==0])

cat("═══════════════════════════════════════════════════════════════════\n")
cat("  MECHANISM FOLLOW-UP: T1 REVERSAL AND T3 N-DEPENDENCE\n")
cat("═══════════════════════════════════════════════════════════════════\n")

# ── T1b: Why do low-SNR genes have high delta.hat? ────────────────────────
cat("\n── T1b: Explaining the T1 reversal ────────────────────────────────────\n")
cat("   Sup step-2 includes class in model -> class variance removed from residuals\n")
cat("   Low-SNR genes: residuals dominated by noise -> noise heterogeneity drives dh\n\n")

ref3 <- as.matrix(read.csv(
  file.path(ADJ_DIR, "combat_sup_n3_testIndia_reference.csv"), row.names=1))
ru3  <- as.matrix(read.csv(
  file.path(ADJ_DIR, "combat_n3_testIndia_reference.csv"), row.names=1))
labs3 <- do.call(c, lapply(train_studies(3, "India"), function(s) bin(label_lst[[s]])))
india_raw  <- log_safe(dat_lst[["India"]])
labs_india <- bin(label_lst[["India"]])
cg3 <- intersect(rownames(ref3), rownames(india_raw))

# Per-gene quantities
gam3  <- rowMeans(ref3[cg3, labs3==1, drop=FALSE]) -
         rowMeans(ref3[cg3, labs3==0, drop=FALSE])
p3    <- mean(labs3)
s2cl  <- gam3^2 * p3*(1-p3)                                    # class variance
s2w   <- pmax(rowVars(ref3[cg3,]), 1e-10) - s2cl               # residual approx
s2w   <- pmax(s2w, 1e-10)
SNR3  <- s2cl / s2w

# Step-2 results (reuse)
r_sup <- full_step2(ref3[cg3,], labs3, india_raw[cg3,], labs_india, supervised=TRUE)

# Decompose test variance per gene: class + noise
gam_tst   <- rowMeans(india_raw[cg3, labs_india==1, drop=FALSE]) -
             rowMeans(india_raw[cg3, labs_india==0, drop=FALSE])
s2cl_tst  <- gam_tst^2 * mean(labs_india)*(1-mean(labs_india))
s2tot_tst <- pmax(rowVars(india_raw[cg3,]), 1e-10)
s2w_tst   <- pmax(s2tot_tst - s2cl_tst, 1e-10)

# For sup step-2: test residual variance AFTER model subtraction ≈ noise only
# (because class effect is removed by beta*class in model)
# Predicted dh_sup ≈ s2w_tst / vp_sup  (noise ratio)
vp_s      <- pmax(r_sup$vp, 1e-10)
dh_pred_noise <- s2w_tst / vp_s    # predicted dh if only noise contributes
rho_noise_pred <- cor(dh_pred_noise, r_sup$dh, method="spearman", use="complete.obs")

# For comparison: if class were NOT removed from test residuals
dh_pred_full  <- s2tot_tst / vp_s  # predicted dh with full test variance
rho_full_pred <- cor(dh_pred_full, r_sup$dh, method="spearman", use="complete.obs")

cat(sprintf("  rho(noise_ratio s2w_tst/vp_sup, dh_sup) = %.3f\n", rho_noise_pred))
cat(sprintf("  rho(total_ratio s2tot_tst/vp_sup, dh_sup) = %.3f\n", rho_full_pred))
cat("  (Higher rho for noise_ratio confirms class IS removed from test residuals)\n")

# Show that dh for high-SNR genes ≈ noise ratio, not class ratio
q90 <- quantile(SNR3, 0.90, na.rm=TRUE)
q10 <- quantile(SNR3, 0.10, na.rm=TRUE)
top <- names(SNR3)[SNR3 >= q90]
bot <- names(SNR3)[SNR3 <= q10]

noise_ratio_top <- median(s2w_tst[top] / vp_s[top], na.rm=TRUE)
noise_ratio_bot <- median(s2w_tst[bot] / vp_s[bot], na.rm=TRUE)
full_ratio_top  <- median(s2tot_tst[top] / vp_s[top], na.rm=TRUE)
full_ratio_bot  <- median(s2tot_tst[bot] / vp_s[bot], na.rm=TRUE)

cat(sprintf("\n  High-SNR genes (top 10%%): noise_ratio=%.2f  dh_sup=%.2f  full_ratio=%.2f\n",
            noise_ratio_top, median(r_sup$dh[top], na.rm=TRUE), full_ratio_top))
cat(sprintf("  Low-SNR  genes (bot 10%%): noise_ratio=%.2f  dh_sup=%.2f  full_ratio=%.2f\n",
            noise_ratio_bot, median(r_sup$dh[bot], na.rm=TRUE), full_ratio_bot))
cat("  -> dh_sup tracks noise_ratio, not full_ratio, confirming class is model-removed\n")

# But class axis compression still happens because compression_g = sqrt(vp/dh)/SD_total
# For high-SNR genes: SD_total >> sqrt(vp) -> large compression even at moderate dh
w3 <- gam3 / sqrt(sum(gam3^2))
comp_g <- sqrt(vp_s / pmax(r_sup$dh, 1e-10)) / sqrt(s2tot_tst)  # SD ratio after/before
cat(sprintf("\n  Median compression (sqrt(vp/dh)/SD_total): overall=%.3f\n",
            median(comp_g, na.rm=TRUE)))
cat(sprintf("    high-SNR genes: %.3f  (large sep -> large compression)\n",
            median(comp_g[top], na.rm=TRUE)))
cat(sprintf("    low-SNR  genes: %.3f\n", median(comp_g[bot], na.rm=TRUE)))
cat("  -> High-SNR genes are MORE compressed in absolute SD terms, explaining T2\n")

# ── T2b: Class-axis separation across all scenarios (fixed ratio) ──────────
cat("\n── T2b: Class-axis separation across all scenarios (corrected ratios) ─\n\n")
m1 <- read.csv(file.path(OUT_DIR, "m1_pc1_alignment.csv"))

sep_rows <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    rs <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_reference.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    if (is.null(rs)) next
    lr  <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
    tst  <- log_safe(dat_lst[[test]])
    lt   <- bin(label_lst[[test]])
    cg   <- intersect(rownames(rs), rownames(tst))
    gam  <- rowMeans(rs[cg, lr==1, drop=FALSE]) -
            rowMeans(rs[cg, lr==0, drop=FALSE])
    w    <- gam / sqrt(sum(gam^2))
    s_b  <- as.numeric(t(tst[cg,]) %*% w)
    sep_b <- sep_fn(s_b, lt)
    # Supervised step-2
    rs2  <- full_step2(rs[cg,], lr, tst[cg,], lt, supervised=TRUE)
    sep_s <- if (!is.null(rs2)) sep_fn(as.numeric(t(rs2$corrected) %*% w), lt) else NA
    # Nat step-2 (mean-only)
    rn2  <- full_step2(rs[cg,], lr, tst[cg,], lt, supervised=FALSE, mean_only=TRUE)
    sep_n <- if (!is.null(rn2)) sep_fn(as.numeric(t(rn2$corrected) %*% w), lt) else NA
    # Ratio: use sign-aware absolute reduction
    abs_ratio_sup <- abs(sep_s) / pmax(abs(sep_b), 1e-10)
    sign_flip_sup <- sign(sep_s) != sign(sep_b)
    abs_ratio_nat <- abs(sep_n) / pmax(abs(sep_b), 1e-10)
    # Purity
    pur_row <- m1[m1$n==n & m1$test==test, ]
    pur_sup <- if (nrow(pur_row)>0) pur_row$pur_sup[1] else NA
    sep_rows[[length(sep_rows)+1]] <- data.frame(
      n=n, test=test, sep_before=sep_b, sep_sup=sep_s, sep_nat=sep_n,
      abs_ratio_sup=abs_ratio_sup, sign_flip=sign_flip_sup,
      abs_ratio_nat=abs_ratio_nat, pur_sup=pur_sup
    )
  }
}
sep_all <- do.call(rbind, sep_rows)

cat("  |sep_sup|/|sep_before| < 1 means class signal reduced;\n")
cat("  sign_flip=TRUE means class order INVERTED after supervised step-2\n\n")
print(sep_all[, c("n","test","sep_before","sep_sup","abs_ratio_sup","sign_flip","pur_sup")],
      row.names=FALSE, digits=3)

# Summary
n_flip <- sum(sep_all$sign_flip, na.rm=TRUE)
cat(sprintf("\n  Class-axis sign flip (Active<->Control inverted): %d / %d scenarios\n",
            n_flip, nrow(sep_all)))
cat(sprintf("  Median |sep| ratio after sup step-2: %.3f  (nat: %.3f)\n",
            median(sep_all$abs_ratio_sup, na.rm=TRUE),
            median(sep_all$abs_ratio_nat, na.rm=TRUE)))

rho_sep_pur <- cor(sep_all$abs_ratio_sup, sep_all$pur_sup,
                   method="spearman", use="complete.obs")
cat(sprintf("  Spearman rho(sep_retention, purity_sup) = %.3f\n", rho_sep_pur))
write.csv(sep_all, file.path(OUT_DIR, "class_sep_ratios_v2.csv"), row.names=FALSE)

# ── T3: N-dependence: class signal in PC2-50 ────────────────────────────
cat("\n── T3: Why does n help? PC2-50 class signal as mediator ───────────────\n")
cat("   India test: purity 0.305->0.662 despite dh ~3.9-4.0\n\n")

t3_rows <- list()
for (n in 2:5) {
  rs  <- as.matrix(read.csv(
    file.path(ADJ_DIR, sprintf("combat_sup_n%d_testIndia_reference.csv", n)),
    row.names=1))
  lr  <- do.call(c, lapply(train_studies(n, "India"), function(s) bin(label_lst[[s]])))
  cg  <- intersect(rownames(rs), rownames(india_raw))

  # PCA of supervised training reference
  n_pcs <- min(50, ncol(rs)-1, nrow(rs)-1)
  pca   <- prcomp(t(rs[cg,]), center=TRUE, scale.=FALSE, rank.=n_pcs)
  frac1 <- pca$sdev[1]^2 / sum(pca$sdev^2)
  frac15 <- sum(pca$sdev[1:min(5,n_pcs)]^2) / sum(pca$sdev^2)

  # Class axis in gene space
  gam   <- rowMeans(rs[cg, lr==1, drop=FALSE]) -
           rowMeans(rs[cg, lr==0, drop=FALSE])
  w_cls <- gam / sqrt(sum(gam^2))
  # PC loadings
  V     <- pca$rotation      # genes x n_pcs
  # Class signal in each PC = |loading · class_direction| weighted by PC variance
  cls_in_pc <- abs(as.numeric(t(V[names(w_cls),]) %*% w_cls))
  # PC1 class loading and cumulative class signal in PC2-50
  cls_pc1  <- cls_in_pc[1]
  cls_pc25 <- sum(cls_in_pc[2:min(5,length(cls_in_pc))])   # PCs 2-5
  cls_pcall <- sum(cls_in_pc)                               # all PCs

  # Step-2 for India
  r_s   <- full_step2(rs[cg,], lr, india_raw[cg,], labs_india, supervised=TRUE)
  r_n   <- full_step2(rs[cg,], lr, india_raw[cg,], labs_india,
                      supervised=FALSE, mean_only=TRUE)
  dh_med <- median(r_s$dh, na.rm=TRUE)

  # Project corrected test onto PCA space; class separation per PC
  if (!is.null(r_s)) {
    tr_sup <- predict(pca, t(r_s$corrected))   # test in PCA space after sup step-2
    tr_nat <- predict(pca, t(r_n$corrected))   # test in PCA space after nat step-2
    tr_raw <- predict(pca, t(india_raw[cg,]))  # raw test in PCA space
    # Class separation per PC (Active mean - Control mean in each PC)
    pc_sep_sup <- sapply(seq_len(ncol(tr_sup)), function(j)
      sep_fn(tr_sup[,j], labs_india))
    pc_sep_nat <- sapply(seq_len(ncol(tr_nat)), function(j)
      sep_fn(tr_nat[,j], labs_india))
    pc_sep_raw <- sapply(seq_len(ncol(tr_raw)), function(j)
      sep_fn(tr_raw[,j], labs_india))
    # Total class signal = sum of |sep| weighted by PC variance (normalized)
    # (KNN uses Euclidean distance across all PCs equally scaled)
    tot_sep_raw <- sum(pc_sep_raw^2)
    tot_sep_sup <- sum(pc_sep_sup^2)
    tot_sep_nat <- sum(pc_sep_nat^2)
    # Signal in PC1 vs rest
    sep_pc1_sup  <- pc_sep_sup[1]^2
    sep_rest_sup <- sum(pc_sep_sup[2:length(pc_sep_sup)]^2)
    sep_pc1_nat  <- pc_sep_nat[1]^2
    sep_rest_nat <- sum(pc_sep_nat[2:length(pc_sep_nat)]^2)
  } else {
    tot_sep_sup <- sep_pc1_sup <- sep_rest_sup <- NA
    tot_sep_nat <- sep_pc1_nat <- sep_rest_nat <- NA
    tot_sep_raw <- NA
  }

  # Purity
  pur_row <- m1[m1$n==n & m1$test=="India", ]
  pur_sup <- if (nrow(pur_row)>0) pur_row$pur_sup[1] else NA

  t3_rows[[length(t3_rows)+1]] <- data.frame(
    n=n, n_train=ncol(rs),
    frac_pc1=frac1, frac_pc1_5=frac15,
    cls_in_pc1=cls_pc1, cls_in_pc2_5=cls_pc25, cls_in_all=cls_pcall,
    dh_med=dh_med,
    tot_sep_raw=tot_sep_raw,
    tot_sep_sup=tot_sep_sup, tot_sep_nat=tot_sep_nat,
    sep_pc1_sup=sep_pc1_sup, sep_rest_sup=sep_rest_sup,
    sep_frac_sup = tot_sep_sup / pmax(tot_sep_raw, 1e-10),
    sep_frac_nat = tot_sep_nat / pmax(tot_sep_raw, 1e-10),
    pur_sup=pur_sup
  )
}
t3 <- do.call(rbind, t3_rows)

cat("  Per-n training structure and class signal in PCA space (India test):\n")
print(t3[, c("n","frac_pc1","cls_in_pc1","cls_in_pc2_5",
             "dh_med","sep_frac_sup","sep_frac_nat","pur_sup")],
      row.names=FALSE, digits=3)

cat(sprintf("\n  n=3 -> n=4 transitions:\n"))
cat(sprintf("    PC1 fraction:     %.3f -> %.3f\n", t3$frac_pc1[1], t3$frac_pc1[2]))
cat(sprintf("    Class in PC1:     %.3f -> %.3f\n", t3$cls_in_pc1[1], t3$cls_in_pc1[2]))
cat(sprintf("    Class in PC2-5:   %.3f -> %.3f\n", t3$cls_in_pc2_5[1], t3$cls_in_pc2_5[2]))
cat(sprintf("    Total class signal retained (sup): %.3f -> %.3f\n",
            t3$sep_frac_sup[1], t3$sep_frac_sup[2]))
cat(sprintf("    Total class signal retained (nat): %.3f -> %.3f\n",
            t3$sep_frac_nat[1], t3$sep_frac_nat[2]))
cat(sprintf("    delta.hat:        %.2f -> %.2f\n", t3$dh_med[1], t3$dh_med[2]))
cat(sprintf("    KNN purity:       %.3f -> %.3f\n", t3$pur_sup[1], t3$pur_sup[2]))

# Predictor correlations (India, n=3,4,5,6)
cat("\n  Correlations with KNN purity across n for India:\n")
for (v in c("frac_pc1","frac_pc1_5","cls_in_pc1","cls_in_pc2_5",
            "dh_med","sep_frac_sup","tot_sep_sup")) {
  r <- suppressWarnings(cor(t3[[v]], t3$pur_sup,
                             method="spearman", use="complete.obs"))
  cat(sprintf("    rho(%-20s, purity) = %+.3f\n", v, r))
}

write.csv(t3, file.path(OUT_DIR, "n_dependence_india_v2.csv"), row.names=FALSE)

# ── Grand summary ─────────────────────────────────────────────────────────
cat("\n\n═══════════════════════════════════════════════════════════════════\n")
cat("  MECHANISM GRAND SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("Why compression destroys class signal (T1b):\n")
cat(sprintf("  sup step-2 removes class from test residuals (rho noise_pred/dh = %.3f)\n",
            rho_noise_pred))
cat(sprintf("  But compression_g = sqrt(vp/dh)/SD_total: high-SNR genes\n"))
cat(sprintf("    have large SD_total (class-dominated) -> large absolute compression\n"))
cat(sprintf("  High-SNR gene compression: %.3f  Low-SNR: %.3f\n",
            median(comp_g[top],na.rm=TRUE), median(comp_g[bot],na.rm=TRUE)))

cat("\nDirect class signal destruction (T2b):\n")
cat(sprintf("  Median |sep| retained by sup: %.3f\n",
            median(sep_all$abs_ratio_sup, na.rm=TRUE)))
cat(sprintf("  Median |sep| retained by nat: %.3f\n",
            median(sep_all$abs_ratio_nat, na.rm=TRUE)))
cat(sprintf("  Class axis INVERTED in %d/%d scenarios after sup step-2\n",
            n_flip, nrow(sep_all)))
cat(sprintf("  rho(sep_retention, purity) = %.3f\n", rho_sep_pur))

cat("\nN-dependence mechanism (T3):\n")
cat(sprintf("  PC1 fraction: n=3=%.3f, n=4=%.3f, n=5=%.3f, n=6=%.3f\n",
            t3$frac_pc1[1],t3$frac_pc1[2],t3$frac_pc1[3],t3$frac_pc1[4]))
cat(sprintf("  Class in PC2-5: n=3=%.3f, n=4=%.3f, n=5=%.3f, n=6=%.3f\n",
            t3$cls_in_pc2_5[1],t3$cls_in_pc2_5[2],t3$cls_in_pc2_5[3],t3$cls_in_pc2_5[4]))
cat(sprintf("  Total class signal retained (sup): n=3=%.3f, n=4=%.3f, n=5=%.3f, n=6=%.3f\n",
            t3$sep_frac_sup[1],t3$sep_frac_sup[2],t3$sep_frac_sup[3],t3$sep_frac_sup[4]))
