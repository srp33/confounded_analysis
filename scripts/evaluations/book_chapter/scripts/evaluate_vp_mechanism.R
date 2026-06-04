#!/usr/bin/env Rscript
# evaluate_vp_mechanism.R
# Root-cause investigation: WHY is delta.hat elevated?
#
# Core hypothesis: supervised step-1 strips class variance from training
# residuals, producing a smaller var_pooled. The test data still contains
# class variance, so it looks inflated relative to training -> delta.hat >> 1.
#
# Quantitative prediction:
#   vp_uns / vp_sup  ≈  1 + sigma2_class / sigma2_residual   (per gene)
#   dh_sup / dh_uns  ≈  vp_uns / vp_sup                       (same test, different training)
#   plain combat step-2 delta.hat  ≈  1.0                     (uns training matches test scale)

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
    mn <- min(m, na.rm=TRUE)
    if (mn < 0) m <- m - mn
    m <- log2(m + 1)
  }
  m
}

prep_raw <- function(n, test) {
  refs <- train_studies(n, test)
  cg   <- Reduce(intersect, lapply(c(refs,test), function(s) rownames(dat_lst[[s]])))
  lab  <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))
  tst  <- dat_lst[[test]][cg,]
  ltst <- bin(label_lst[[test]])
  list(tst=log_safe(tst), ltst=ltst, lab=lab, cg=cg)
}

# Returns delta.hat AND var_pooled vectors for the test batch.
# Design: supervised (batch+class) or unsupervised (batch only).
dh_vp <- function(train_mat, train_lab, test_mat, test_lab, supervised=TRUE) {
  cg   <- intersect(rownames(train_mat), rownames(test_mat))
  comb <- cbind(train_mat[cg,], test_mat[cg,])
  bat  <- factor(c(rep(1,ncol(train_mat)), rep(2,ncol(test_mat))))
  des  <- if (supervised) model.matrix(~0+bat+c(train_lab,test_lab))
          else             model.matrix(~0+bat)
  B    <- tryCatch(solve(crossprod(des), tcrossprod(t(des),comb)), error=function(e) NULL)
  if (is.null(B)) return(list(dh=rep(NA,nrow(comb)), vp=rep(NA,nrow(comb))))
  ri   <- seq_len(ncol(train_mat))
  ti   <- (ncol(train_mat)+1):ncol(comb)
  res  <- comb[,ri] - t(des[ri,] %*% B)
  vp   <- pmax(rowMeans(res^2), 1e-10)
  s    <- (comb[,ti] - t(des[ti,] %*% B)) / sqrt(vp)
  list(dh=rowVars(s), vp=vp)
}

cat("═══════════════════════════════════════════════════════════════════\n")
cat("  VAR_POOLED MECHANISM INVESTIGATION\n")
cat("═══════════════════════════════════════════════════════════════════\n")

rows <- list()

for (n in 2:5) {
  for (test in ALL_STUDIES) {
    rs <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_reference.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    ru <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_n%d_test%s_reference.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    if (is.null(rs) || is.null(ru)) next

    sc  <- prep_raw(n, test)
    lr  <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
    cg  <- intersect(rownames(rs), sc$cg)
    lt  <- sc$ltst

    # A: step-2 using SUPERVISED step-1 training (sup design)
    a <- dh_vp(rs[cg,], lr, sc$tst[cg,], lt, supervised=TRUE)
    # B: step-2 using UNSUPERVISED step-1 training (sup design — same test)
    b <- dh_vp(ru[cg,], lr, sc$tst[cg,], lt, supervised=TRUE)
    # C: step-2 using UNSUPERVISED step-1 training (uns design — plain combat)
    cc <- dh_vp(ru[cg,], lr, sc$tst[cg,], lt, supervised=FALSE)

    # Summaries
    med_dh_sup  <- median(a$dh, na.rm=TRUE)
    med_dh_uns  <- median(b$dh, na.rm=TRUE)
    med_dh_plain<- median(cc$dh, na.rm=TRUE)
    med_vp_sup  <- median(a$vp, na.rm=TRUE)
    med_vp_uns  <- median(b$vp, na.rm=TRUE)
    med_vp_plain<- median(cc$vp, na.rm=TRUE)

    # Observed dh ratio vs predicted from vp ratio (per gene)
    # Prediction: dh_sup / dh_uns ≈ vp_uns / vp_sup  (same test, different training vp)
    vp_ratio_obs <- b$vp / pmax(a$vp, 1e-10)   # vp_uns / vp_sup per gene
    dh_ratio_obs <- a$dh / pmax(b$dh, 1e-10)   # dh_sup / dh_uns per gene
    rho_pred     <- cor(vp_ratio_obs, dh_ratio_obs, method="spearman", use="complete.obs")
    med_vp_ratio <- median(vp_ratio_obs, na.rm=TRUE)
    med_dh_ratio <- median(dh_ratio_obs, na.rm=TRUE)

    # Class variance decomposition from supervised reference
    # sigma2_class per gene = gamma^2 * p*(1-p)
    p      <- mean(lr)
    gamma  <- rowMeans(rs[cg, lr==1, drop=FALSE]) - rowMeans(rs[cg, lr==0, drop=FALSE])
    s2cls  <- gamma^2 * p * (1-p)
    # sigma2_residual ≈ var_pooled from supervised step-2 (a$vp)
    # Predicted vp ratio = 1 + sigma2_class / sigma2_residual
    pred_vp_ratio <- 1 + s2cls / pmax(a$vp, 1e-10)
    rho_cls_vp    <- cor(pred_vp_ratio, vp_ratio_obs, method="spearman", use="complete.obs")
    med_pred_vp   <- median(pred_vp_ratio, na.rm=TRUE)

    rows[[length(rows)+1]] <- data.frame(
      n=n, test=test,
      med_vp_sup=med_vp_sup, med_vp_uns=med_vp_uns, med_vp_plain=med_vp_plain,
      med_dh_sup=med_dh_sup, med_dh_uns=med_dh_uns, med_dh_plain=med_dh_plain,
      med_vp_ratio=med_vp_ratio, med_dh_ratio=med_dh_ratio,
      med_pred_vp=med_pred_vp,
      rho_pred=rho_pred, rho_cls_vp=rho_cls_vp,
      med_sd_test_raw=median(rowSds(sc$tst[cg,]), na.rm=TRUE)
    )
  }
}

res <- do.call(rbind, rows)

# ── Test 1: vp_sup vs vp_uns ──────────────────────────────────────────────
cat("\n── Test 1: var_pooled comparison (sup vs uns step-1) ─────────────────\n")
cat("   Prediction: vp_sup << vp_uns (class variance stripped from training)\n\n")

cat("  Per-scenario median var_pooled and delta.hat:\n")
print(res[, c("n","test","med_vp_sup","med_vp_uns","med_vp_plain",
              "med_dh_sup","med_dh_uns","med_dh_plain")],
      row.names=FALSE, digits=3)

cat(sprintf("\n  Overall median vp_sup=%.4f  vp_uns=%.4f  vp_plain=%.4f\n",
            median(res$med_vp_sup), median(res$med_vp_uns), median(res$med_vp_plain)))
cat(sprintf("  Overall median dh_sup=%.3f   dh_uns=%.3f   dh_plain=%.3f\n",
            median(res$med_dh_sup), median(res$med_dh_uns), median(res$med_dh_plain)))

wt_vp <- wilcox.test(res$med_vp_sup, res$med_vp_uns, paired=TRUE, exact=FALSE)
wt_dh <- wilcox.test(res$med_dh_sup, res$med_dh_uns, paired=TRUE, exact=FALSE)
cat(sprintf("  Wilcoxon vp_sup < vp_uns: p=%.3g\n", wt_vp$p.value))
cat(sprintf("  Wilcoxon dh_sup > dh_uns: p=%.3g\n", wt_dh$p.value))

# ── Test 2: dh_ratio ≈ vp_ratio ──────────────────────────────────────────
cat("\n── Test 2: dh_sup/dh_uns ≈ vp_uns/vp_sup? ────────────────────────────\n")
cat("   Prediction: delta.hat inflation = var_pooled compression ratio\n\n")

cat(sprintf("  Median per-scenario vp_uns/vp_sup = %.3f  (expected > 1)\n",
            median(res$med_vp_ratio)))
cat(sprintf("  Median per-scenario dh_sup/dh_uns = %.3f  (expected ≈ vp ratio)\n",
            median(res$med_dh_ratio)))
cat(sprintf("  Ratio of ratios (dh/vp): %.3f  (1.0 = perfect match)\n",
            median(res$med_dh_ratio) / median(res$med_vp_ratio)))
cat(sprintf("  Median per-gene rho(vp_ratio, dh_ratio): %.3f\n",
            median(res$rho_pred)))

# ── Test 3: class variance predicts vp_ratio ──────────────────────────────
cat("\n── Test 3: sigma2_class/sigma2_resid predicts vp_uns/vp_sup ──────────\n")
cat("   Prediction: rho(1+s2cls/vp_sup, vp_uns/vp_sup) ≈ 1 per scenario\n\n")

cat(sprintf("  Median predicted vp_ratio (1+s2cls/vp_sup) = %.3f\n",
            median(res$med_pred_vp)))
cat(sprintf("  Median observed  vp_ratio (vp_uns/vp_sup)  = %.3f\n",
            median(res$med_vp_ratio)))
cat(sprintf("  Per-gene Spearman rho(predicted, observed): median=%.3f  [IQR: %.3f-%.3f]\n",
            median(res$rho_cls_vp),
            quantile(res$rho_cls_vp, 0.25), quantile(res$rho_cls_vp, 0.75)))

# ── Test 4: plain combat delta.hat ≈ 1 ────────────────────────────────────
cat("\n── Test 4: plain combat (uns step-1, uns design) delta.hat ───────────\n")
cat("   Prediction: dh_plain ≈ 1.0 (training and test have same variance)\n\n")

cat(sprintf("  Median dh_plain = %.3f  (expected ≈ 1.0)\n", median(res$med_dh_plain)))
cat(sprintf("  Range: [%.3f, %.3f]\n", min(res$med_dh_plain), max(res$med_dh_plain)))
cat(sprintf("  Fraction scenarios with dh_plain < 1.5: %.3f\n",
            mean(res$med_dh_plain < 1.5)))

# ── Test 5: nat step-2 sees same vp_sup but ignores dh ────────────────────
cat("\n── Test 5: nat step-2 sees elevated dh but ignores it ────────────────\n")
cat("   (Confirming mean-only is the bypass, not a different delta.hat)\n\n")

# dh_nat from M2 is the dh when using sup step-1 + batch-only design
# We already have dh_sup (sup step-1 + class design) and dh_uns (uns step-1 + class design)
# The dh for nat step-2 = dh with sup step-1 + batch-only design = compute here
rows2 <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    rs <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_reference.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    if (is.null(rs)) next
    sc  <- prep_raw(n, test)
    lr  <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
    cg  <- intersect(rownames(rs), sc$cg)
    # nat step-2: sup step-1 training + batch-only design (no class)
    d   <- dh_vp(rs[cg,], lr, sc$tst[cg,], sc$ltst, supervised=FALSE)
    rows2[[length(rows2)+1]] <- data.frame(n=n, test=test,
                                            med_dh_nat_step2=median(d$dh, na.rm=TRUE))
  }
}
res2 <- do.call(rbind, rows2)
res  <- merge(res, res2, by=c("n","test"))

cat("  Step-2 delta.hat values (same test data, different step-1 and design):\n")
cat(sprintf("  Method                        | median delta.hat\n"))
cat(sprintf("  sup step-1 + class design     | %.3f  (combat_sup step-2)\n", median(res$med_dh_sup)))
cat(sprintf("  sup step-1 + batch-only design| %.3f  (nat step-2 — ignored by mean.only)\n", median(res$med_dh_nat_step2)))
cat(sprintf("  uns step-1 + class design     | %.3f  (same design, different step-1)\n", median(res$med_dh_uns)))
cat(sprintf("  uns step-1 + batch-only design| %.3f  (plain combat step-2)\n", median(res$med_dh_plain)))

# ── Causal chain summary ──────────────────────────────────────────────────
cat("\n\n═══════════════════════════════════════════════════════════════════\n")
cat("  CAUSAL CHAIN SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("Step-1 effect on var_pooled:\n")
cat(sprintf("  Supervised step-1 strips class variance from training residuals.\n"))
cat(sprintf("  vp_sup = %.4f  vs  vp_uns = %.4f  (ratio = %.2f×)\n",
            median(res$med_vp_sup), median(res$med_vp_uns),
            median(res$med_vp_uns) / median(res$med_vp_sup)))

cat("\nEffect on delta.hat (same test, different training vp):\n")
cat(sprintf("  dh_sup = %.2f  vs  dh_uns = %.2f  vs  dh_plain = %.2f\n",
            median(res$med_dh_sup), median(res$med_dh_uns), median(res$med_dh_plain)))
cat(sprintf("  dh_ratio = %.2f matches vp_ratio = %.2f?  ratio-of-ratios = %.3f\n",
            median(res$med_dh_ratio), median(res$med_vp_ratio),
            median(res$med_dh_ratio) / median(res$med_vp_ratio)))

cat("\nClass variance as the driver:\n")
cat(sprintf("  Per-gene: rho(1 + s2cls/s2resid, vp_uns/vp_sup) = %.3f\n",
            median(res$rho_cls_vp)))
cat(sprintf("  Predicted vp_ratio = %.3f  observed = %.3f\n",
            median(res$med_pred_vp), median(res$med_vp_ratio)))

cat("\nEffect on compression (median test SD ratio after step-2):\n")
comp_sup   <- median(sqrt(res$med_vp_sup / res$med_dh_sup) /
                     pmax(res$med_sd_test_raw, 1e-10), na.rm=TRUE)
cat(sprintf("  Sup compression ≈ sqrt(vp_sup/dh_sup) / SD_test_raw = %.3f (nat = 1.000)\n",
            comp_sup))

cat("\nConclusion:\n")
cat("  Supervised step-1 reduces var_pooled by removing class variance.\n")
cat("  Same test data then appears inflated (dh >> 1).\n")
cat("  Nat/mo bypasses the damage by applying mean-only step-2 (ignoring dh).\n")
cat("  Plain combat avoids the gap by using unsupervised step-1 (vp includes class).\n")

write.csv(res, file.path(OUT_DIR, "vp_mechanism.csv"), row.names=FALSE)
cat(sprintf("\nResults written to %s/vp_mechanism.csv\n", OUT_DIR))
