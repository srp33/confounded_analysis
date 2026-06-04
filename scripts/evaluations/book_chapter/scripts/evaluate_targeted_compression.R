#!/usr/bin/env Rscript
# evaluate_targeted_compression.R
# Three targeted tests of the compression mechanism:
#
# T1: Per-gene — delta.hat_sup is highest for the most class-informative genes
#     (supervised specifically destroys signal where it matters most)
#
# T2: Class-axis scores — supervised step-2 collapses the class-axis separation
#     of test samples to near-zero; nat step-2 preserves it
#     (direct explanation of KNN failure)
#
# T3: Why does n help despite constant delta.hat?
#     India purity 0.305 (n=3) -> 0.662 (n=4) while dh stays ~3.9
#     Hypothesis: more training batches -> lower PC1_frac -> compression
#     distributed across more axes -> class signal survives in other PCs

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

ALL_STUDIES <- c("GSE37250_SA","USA","India","GSE37250_M","Africa","GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1",1,"Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES != test][seq_len(n)]
log_safe <- function(m) {
  if (max(m, na.rm=TRUE) > 100) {
    mn <- min(m, na.rm=TRUE); if (mn < 0) m <- m - mn; m <- log2(m+1)
  }; m
}

# Returns delta.hat, gamma.hat, var_pooled, and corrected test matrix
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
  gh    <- rowMeans(s_tst)           # gamma.hat per gene (batch mean in std space)
  if (mean_only) {
    # Mean-only correction: shift only
    corrected <- test_mat[cg,] - gh * sqrt(vp)
  } else {
    # Full correction: shift + scale
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
  d   <- as.matrix(dist(rbind(tr,rr)))[
    seq_len(nrow(tr)), (nrow(tr)+1):(nrow(tr)+nrow(rr))]
  mean(sapply(seq_len(nrow(tr)), function(i)
    mean(ref_lab[order(d[i,])[seq_len(k)]] == tst_lab[i])))
}

cat("═══════════════════════════════════════════════════════════════════\n")
cat("  TARGETED COMPRESSION INVESTIGATION\n")
cat("═══════════════════════════════════════════════════════════════════\n")

# Primary scenario: n=3, test=India (most catastrophic)
# Also include n=4, test=India (recovery scenario) for T3
cat("\nLoading reference matrices for India analysis...\n")
load_ref <- function(method, n) as.matrix(read.csv(
  file.path(ADJ_DIR, sprintf("%s_n%d_testIndia_reference.csv", method, n)),
  row.names=1))

ref_sup3 <- load_ref("combat_sup", 3)
ref_uns3 <- load_ref("combat",     3)
ref_sup4 <- load_ref("combat_sup", 4)
ref_uns4 <- load_ref("combat",     4)

labs3 <- do.call(c, lapply(train_studies(3, "India"), function(s) bin(label_lst[[s]])))
labs4 <- do.call(c, lapply(train_studies(4, "India"), function(s) bin(label_lst[[s]])))

india_raw  <- log_safe(dat_lst[["India"]])
labs_india <- bin(label_lst[["India"]])
cg3 <- intersect(rownames(ref_sup3), rownames(india_raw))
cg4 <- intersect(rownames(ref_sup4), rownames(india_raw))

cat(sprintf("  n=3 genes: %d,  n=4 genes: %d\n", length(cg3), length(cg4)))

# ── T1: Per-gene targeted compression ─────────────────────────────────────
cat("\n── T1: Per-gene class SNR vs delta.hat ────────────────────────────────\n")
cat("   Prediction: rho(class_SNR, dh_sup) >> rho(class_SNR, dh_nat)\n\n")

# Class SNR per gene from supervised training reference
p3   <- mean(labs3)
gam3 <- rowMeans(ref_sup3[cg3, labs3==1, drop=FALSE]) -
        rowMeans(ref_sup3[cg3, labs3==0, drop=FALSE])
# Within-class variance = var_pooled from supervised training
s2w3 <- rowVars(ref_sup3[cg3, labs3==1, drop=FALSE]) * sum(labs3==1-1) +
        rowVars(ref_sup3[cg3, labs3==0, drop=FALSE]) * sum(labs3==0-1)
s2w3 <- s2w3 / (length(labs3) - 2)
s2w3 <- pmax(s2w3, 1e-10)
SNR3 <- gam3^2 * p3 * (1-p3) / s2w3   # class variance / within-class variance

# Step-2 delta.hat under three conditions (same test, n=3):
r_sup  <- full_step2(ref_sup3[cg3,], labs3, india_raw[cg3,], labs_india, supervised=TRUE)
r_nat  <- full_step2(ref_sup3[cg3,], labs3, india_raw[cg3,], labs_india, supervised=FALSE)
r_plain <- full_step2(ref_uns3[cg3,], labs3, india_raw[cg3,], labs_india, supervised=FALSE)

rho_snr_sup   <- cor(SNR3, r_sup$dh,   method="spearman", use="complete.obs")
rho_snr_nat   <- cor(SNR3, r_nat$dh,   method="spearman", use="complete.obs")
rho_snr_plain <- cor(SNR3, r_plain$dh, method="spearman", use="complete.obs")

cat(sprintf("  n=3, test=India:\n"))
cat(sprintf("    rho(class_SNR, dh_sup)   = %+.3f  [+: informative genes hit hardest]\n",
            rho_snr_sup))
cat(sprintf("    rho(class_SNR, dh_nat)   = %+.3f  [reference: same training, no class in design]\n",
            rho_snr_nat))
cat(sprintf("    rho(class_SNR, dh_plain) = %+.3f  [reference: unsupervised training]\n",
            rho_snr_plain))

# Top decile of class-SNR genes vs bottom decile delta.hat comparison
q90 <- quantile(SNR3, 0.90, na.rm=TRUE)
q10 <- quantile(SNR3, 0.10, na.rm=TRUE)
top_genes <- names(SNR3)[SNR3 >= q90]
bot_genes <- names(SNR3)[SNR3 <= q10]

cat(sprintf("\n  Top 10%% class-SNR genes: median dh_sup=%.2f  dh_nat=%.2f  dh_plain=%.2f\n",
            median(r_sup$dh[top_genes], na.rm=TRUE),
            median(r_nat$dh[top_genes], na.rm=TRUE),
            median(r_plain$dh[top_genes], na.rm=TRUE)))
cat(sprintf("  Bot 10%% class-SNR genes: median dh_sup=%.2f  dh_nat=%.2f  dh_plain=%.2f\n",
            median(r_sup$dh[bot_genes], na.rm=TRUE),
            median(r_nat$dh[bot_genes], na.rm=TRUE),
            median(r_plain$dh[bot_genes], na.rm=TRUE)))
cat(sprintf("  Ratio top/bot dh_sup=%.2f  dh_nat=%.2f  dh_plain=%.2f\n",
            median(r_sup$dh[top_genes],na.rm=TRUE)/median(r_sup$dh[bot_genes],na.rm=TRUE),
            median(r_nat$dh[top_genes],na.rm=TRUE)/median(r_nat$dh[bot_genes],na.rm=TRUE),
            median(r_plain$dh[top_genes],na.rm=TRUE)/median(r_plain$dh[bot_genes],na.rm=TRUE)))

# Save per-gene data for plotting
pg <- data.frame(gene=cg3, SNR=SNR3,
                 dh_sup=r_sup$dh, dh_nat=r_nat$dh, dh_plain=r_plain$dh,
                 gh_sup=r_sup$gh, gh_nat=r_nat$gh,
                 vp_sup=r_sup$vp, vp_nat=r_nat$vp)
write.csv(pg, file.path(OUT_DIR, "per_gene_targeted.csv"), row.names=FALSE)

# ── T2: Class-axis score compression ──────────────────────────────────────
cat("\n── T2: Class-axis score compression in test samples ───────────────────\n")
cat("   Prediction: sup step-2 collapses class-axis separation; nat preserves it\n\n")

# Class axis (unit vector) from supervised training reference
w3 <- gam3 / sqrt(sum(gam3^2))

# Class-axis scores: project each sample onto w
score_before <- as.numeric(t(india_raw[cg3,]) %*% w3)  # raw test
score_sup    <- as.numeric(t(r_sup$corrected)  %*% w3)  # after sup step-2
score_nat    <- as.numeric(t(full_step2(ref_sup3[cg3,], labs3,
                                        india_raw[cg3,], labs_india,
                                        supervised=FALSE, mean_only=TRUE)$corrected) %*% w3)

sep_fn <- function(s, l) mean(s[l==1]) - mean(s[l==0])

sep_before <- sep_fn(score_before, labs_india)
sep_sup    <- sep_fn(score_sup,    labs_india)
sep_nat    <- sep_fn(score_nat,    labs_india)

sd_before  <- sd(score_before)
sd_sup     <- sd(score_sup)
sd_nat     <- sd(score_nat)

cat(sprintf("  Class-axis separation (Active mean - Control mean):\n"))
cat(sprintf("    Before step-2 (raw test):  sep=%.3f  sd=%.3f\n", sep_before, sd_before))
cat(sprintf("    After sup step-2:          sep=%.3f  sd=%.3f  (ratio=%.3f)\n",
            sep_sup, sd_sup, sep_sup/max(sep_before,1e-10)))
cat(sprintf("    After nat step-2:          sep=%.3f  sd=%.3f  (ratio=%.3f)\n",
            sep_nat, sd_nat, sep_nat/max(sep_before,1e-10)))

# Also show effect on within-class spread
sd_act_before <- sd(score_before[labs_india==1])
sd_ctl_before <- sd(score_before[labs_india==0])
sd_act_sup    <- sd(score_sup[labs_india==1])
sd_ctl_sup    <- sd(score_sup[labs_india==0])

cat(sprintf("\n  Within-class spread on class axis:\n"))
cat(sprintf("    Before: Active sd=%.3f  Control sd=%.3f\n", sd_act_before, sd_ctl_before))
cat(sprintf("    After sup: Active sd=%.3f  Control sd=%.3f\n", sd_act_sup, sd_ctl_sup))

# Gamma.hat alignment with class axis
cos_gh_sup <- sum(r_sup$gh * w3) / sqrt(sum(r_sup$gh^2) * sum(w3^2))
cos_gh_nat <- sum(r_nat$gh * w3) / sqrt(sum(r_nat$gh^2) * sum(w3^2))
cat(sprintf("\n  Gamma.hat class-axis alignment:\n"))
cat(sprintf("    cos(gamma_hat_sup, w_class) = %.3f\n", cos_gh_sup))
cat(sprintf("    cos(gamma_hat_nat, w_class) = %.3f\n", cos_gh_nat))
cat("    (near 0 = neutral; +1 = pushes in class direction; -1 = opposes class)\n")

# Save class-axis scores for plotting
cs_df <- data.frame(
  sample=colnames(india_raw),
  class=labs_india,
  score_before=score_before,
  score_sup=score_sup,
  score_nat=score_nat
)
write.csv(cs_df, file.path(OUT_DIR, "class_axis_scores.csv"), row.names=FALSE)

# Replicate for all 18 scenarios — save class-axis sep ratios
cat("\n  Class-axis separation ratio (after/before) across all 18 scenarios:\n")
cat("  Method: sup step-2 vs nat step-2\n\n")
sep_rows <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    rs <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_reference.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    if (is.null(rs)) next
    lr <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
    tst_raw  <- log_safe(dat_lst[[test]])
    lt       <- bin(label_lst[[test]])
    cg       <- intersect(rownames(rs), rownames(tst_raw))
    gam      <- rowMeans(rs[cg, lr==1, drop=FALSE]) -
                rowMeans(rs[cg, lr==0, drop=FALSE])
    w        <- gam / sqrt(sum(gam^2))
    # Before
    s_bef    <- as.numeric(t(tst_raw[cg,]) %*% w)
    sep_bef  <- sep_fn(s_bef, lt)
    # Sup step-2
    r_s      <- full_step2(rs[cg,], lr, tst_raw[cg,], lt, supervised=TRUE)
    s_s      <- if(!is.null(r_s)) as.numeric(t(r_s$corrected) %*% w) else rep(NA, length(lt))
    sep_s    <- sep_fn(s_s, lt)
    # Nat step-2 (mean-only)
    r_n      <- full_step2(rs[cg,], lr, tst_raw[cg,], lt, supervised=FALSE, mean_only=TRUE)
    s_n      <- if(!is.null(r_n)) as.numeric(t(r_n$corrected) %*% w) else rep(NA, length(lt))
    sep_n    <- sep_fn(s_n, lt)
    sep_rows[[length(sep_rows)+1]] <- data.frame(
      n=n, test=test,
      sep_before=sep_bef, sep_sup=sep_s, sep_nat=sep_n,
      ratio_sup=sep_s/max(sep_bef,1e-10),
      ratio_nat=sep_n/max(sep_bef,1e-10)
    )
  }
}
sep_all <- do.call(rbind, sep_rows)
print(sep_all[, c("n","test","sep_before","sep_sup","sep_nat","ratio_sup","ratio_nat")],
      row.names=FALSE, digits=3)

cat(sprintf("\n  Median ratio sup=%.3f  nat=%.3f\n",
            median(sep_all$ratio_sup, na.rm=TRUE),
            median(sep_all$ratio_nat, na.rm=TRUE)))

rho_ratio_pur <- cor(sep_all$ratio_sup,
                     (read.csv(file.path(OUT_DIR,"m1_pc1_alignment.csv")) %>%
                       `[`(order(.$n,.$test),) %>% `[[`("pur_sup"))[seq_len(nrow(sep_all))],
                     method="spearman", use="complete.obs")
cat(sprintf("  Spearman rho(sep_ratio_sup, purity_sup) = %.3f\n", rho_ratio_pur))

write.csv(sep_all, file.path(OUT_DIR, "class_axis_sep_ratios.csv"), row.names=FALSE)

# ── T3: Why does n help despite constant delta.hat? ───────────────────────
cat("\n── T3: Why does n help (India: purity 0.305->0.662 with similar dh)? ──\n")
cat("   Hypothesis: PC1_frac decreases with n -> compression distributed\n")
cat("               across more axes -> class signal survives in other PCs\n\n")

t3_rows <- list()
for (n in 2:5) {
  rs  <- as.matrix(read.csv(
    file.path(ADJ_DIR, sprintf("combat_sup_n%d_testIndia_reference.csv", n)),
    row.names=1))
  lr  <- do.call(c, lapply(train_studies(n, "India"), function(s) bin(label_lst[[s]])))
  cg  <- intersect(rownames(rs), rownames(india_raw))

  # PCA of supervised reference
  pca    <- prcomp(t(rs[cg,]), center=TRUE, scale.=FALSE, rank.=20)
  frac1  <- pca$sdev[1]^2 / sum(pca$sdev^2)
  frac5  <- sum(pca$sdev[1:5]^2) / sum(pca$sdev^2)

  # Class separation in training reference on training class axis
  gam    <- rowMeans(rs[cg, lr==1, drop=FALSE]) -
            rowMeans(rs[cg, lr==0, drop=FALSE])
  w_cls  <- gam / sqrt(sum(gam^2))
  s_tr   <- as.numeric(t(rs[cg,]) %*% w_cls)
  sep_tr <- mean(s_tr[lr==1]) - mean(s_tr[lr==0])

  # Effective dimensionality (participation ratio of PCA spectrum)
  lambda <- pca$sdev^2
  eff_dim <- sum(lambda)^2 / sum(lambda^2)

  # delta.hat for India test
  r_s    <- full_step2(rs[cg,], lr, india_raw[cg,], labs_india, supervised=TRUE)
  dh_med <- median(r_s$dh, na.rm=TRUE)

  # Class-axis survival after sup step-2
  s_bef  <- as.numeric(t(india_raw[cg,]) %*% w_cls)
  s_aft  <- if(!is.null(r_s)) as.numeric(t(r_s$corrected) %*% w_cls) else rep(NA,ncol(india_raw))
  sep_bef <- mean(s_bef[labs_india==1]) - mean(s_bef[labs_india==0])
  sep_aft <- mean(s_aft[labs_india==1], na.rm=TRUE) - mean(s_aft[labs_india==0], na.rm=TRUE)
  sep_ratio <- sep_aft / max(sep_bef, 1e-10)

  # Within-space diversity: cos between PC1 and class axis
  cos_pc1 <- abs(sum(pca$rotation[,1] * w_cls[rownames(pca$rotation)]))

  # KNN purity from m1_pc1_alignment.csv (preloaded)
  m1      <- read.csv(file.path(OUT_DIR, "m1_pc1_alignment.csv"))
  pur_row <- m1[m1$n==n & m1$test=="India",]
  pur_sup <- if(nrow(pur_row)>0) pur_row$pur_sup[1] else NA

  t3_rows[[length(t3_rows)+1]] <- data.frame(
    n=n, n_train=ncol(rs),
    frac_pc1=frac1, frac_pc1_5=frac5, eff_dim=eff_dim,
    cos_pc1_class=cos_pc1,
    sep_tr=sep_tr, dh_med=dh_med,
    sep_ratio=sep_ratio, pur_sup=pur_sup
  )
}
t3 <- do.call(rbind, t3_rows)

cat("  Per-n results for test=India:\n")
print(t3, row.names=FALSE, digits=3)

cat(sprintf("\n  Key pattern: n=3 vs n=4\n"))
cat(sprintf("    PC1 fraction:  n=3=%.3f -> n=4=%.3f  (change: %+.3f)\n",
            t3$frac_pc1[1], t3$frac_pc1[2], t3$frac_pc1[2]-t3$frac_pc1[1]))
cat(sprintf("    PC1-5 fraction: n=3=%.3f -> n=4=%.3f\n",
            t3$frac_pc1_5[1], t3$frac_pc1_5[2]))
cat(sprintf("    Eff dimensionality: n=3=%.1f -> n=4=%.1f\n",
            t3$eff_dim[1], t3$eff_dim[2]))
cat(sprintf("    cos(PC1, class): n=3=%.3f -> n=4=%.3f\n",
            t3$cos_pc1_class[1], t3$cos_pc1_class[2]))
cat(sprintf("    Sep ratio (after/before sup): n=3=%.3f -> n=4=%.3f\n",
            t3$sep_ratio[1], t3$sep_ratio[2]))
cat(sprintf("    delta.hat:     n=3=%.2f -> n=4=%.2f\n",
            t3$dh_med[1], t3$dh_med[2]))
cat(sprintf("    KNN purity:    n=3=%.3f -> n=4=%.3f\n",
            t3$pur_sup[1], t3$pur_sup[2]))

# Cross-n correlation: what predicts purity?
if (nrow(t3) >= 3) {
  cat("\n  Correlations with KNN purity (n=3,4,5,6 for India):\n")
  for (v in c("frac_pc1","frac_pc1_5","eff_dim","cos_pc1_class","dh_med","sep_ratio")) {
    r <- cor(t3[[v]], t3$pur_sup, method="spearman", use="complete.obs")
    cat(sprintf("    rho(%s, purity) = %+.3f\n", v, r))
  }
}

write.csv(t3, file.path(OUT_DIR, "n_dependence_india.csv"), row.names=FALSE)

# ── Summary ───────────────────────────────────────────────────────────────
cat("\n\n═══════════════════════════════════════════════════════════════════\n")
cat("  TARGETED COMPRESSION SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("T1 — Is compression targeted at class-informative genes?\n")
cat(sprintf("  rho(SNR, dh_sup)   = %.3f  (supervised, targeted)\n", rho_snr_sup))
cat(sprintf("  rho(SNR, dh_nat)   = %.3f  (nat/batch-only design)\n", rho_snr_nat))
cat(sprintf("  rho(SNR, dh_plain) = %.3f  (plain combat)\n", rho_snr_plain))
cat(sprintf("  Top vs bot decile dh_sup ratio: %.2f (nat: %.2f)\n",
            median(r_sup$dh[top_genes],na.rm=TRUE)/median(r_sup$dh[bot_genes],na.rm=TRUE),
            median(r_nat$dh[top_genes],na.rm=TRUE)/median(r_nat$dh[bot_genes],na.rm=TRUE)))

cat("\nT2 — Does sup step-2 collapse class-axis separation?\n")
cat(sprintf("  n=3/India class-axis sep ratio: sup=%.3f  nat=%.3f\n",
            sep_sup/sep_before, sep_nat/sep_before))
cat(sprintf("  All-scenario median sep ratio:  sup=%.3f  nat=%.3f\n",
            median(sep_all$ratio_sup,na.rm=TRUE),
            median(sep_all$ratio_nat,na.rm=TRUE)))
cat(sprintf("  rho(sep_ratio, purity) = %.3f\n", rho_ratio_pur))

cat("\nT3 — What drives n-dependence of purity for India?\n")
best_v <- c("frac_pc1","frac_pc1_5","eff_dim","cos_pc1_class","dh_med","sep_ratio")[
  which.max(abs(sapply(c("frac_pc1","frac_pc1_5","eff_dim","cos_pc1_class","dh_med","sep_ratio"),
                       function(v) cor(t3[[v]], t3$pur_sup, method="spearman", use="complete.obs"))))]
cat(sprintf("  Best predictor across n: %s\n", best_v))
cat(sprintf("  sep_ratio (class-axis survival): %.3f -> %.3f -> %.3f -> %.3f\n",
            t3$sep_ratio[1], t3$sep_ratio[2], t3$sep_ratio[3], t3$sep_ratio[4]))
