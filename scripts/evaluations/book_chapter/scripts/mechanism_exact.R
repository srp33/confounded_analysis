#!/usr/bin/env Rscript
# mechanism_exact.R
# EXACT mechanism of the combat_sup x KNN failure, in the PRODUCTION metric.
#
# Production KNN = top-1000 highly-variable genes (by training variance) + plain
# Euclidean class::knn (helper.R). No PCA, no whitening. We reuse the *actual*
# production functions predKNN_pp / predWrapper so the numbers match the pipeline.
#
# Hypothesised chain (from the combat_sup_mg note in adjust_target_data_SA_UK.R):
#   step1 supervised  -> reference within-class var_train_g shrinks on class genes
#   step2 ComBat      -> per-gene delta_g = sd_test_g / sd_train_g blows up there
#                        => test is divided by a large factor on class genes
#                        => class axis in the test is crushed / inverted
#   KNN (Euclidean, no whitening) then picks wrong-class neighbours -> MCC < 0
#   LDA/RDA (S_W^-1 whitening) undoes the per-gene rescaling -> survives.
#
# We verify every link in gene space, and show mean-only step2 (the fix) recovers.

suppressMessages({ library(sva); library(matrixStats); library(class) })
source("scripts/helper.R")          # predKNN_pp, predWrapper (the production KNN)

DATA_FILE <- "data/TB_real_data.RData"
OUT_DIR   <- "outputs/diagnostics/hypothesis_tests"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
load(DATA_FILE)

ALL_STUDIES <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1", 1, "Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES != test][seq_len(n)]
log_safe <- function(m) {
  if (max(m, na.rm = TRUE) > 100) { mn <- min(m, na.rm = TRUE); if (mn < 0) m <- m - mn; m <- log2(m + 1) }
  m
}

# ── Two-step ComBat with independently-set supervision and step-2 mode ─────────
# step2_mode: "full" (per-gene location+scale, mod by sup2), "meanonly" (location
# only, no per-gene delta), used to test the fix.
two_step <- function(dat, bat, lab, tst, ltst, sup1, sup2, step2_mode = "full") {
  if (length(unique(bat)) >= 2) {
    mod1  <- if (sup1) model.matrix(~lab) else NULL
    dat_c <- suppressMessages(ComBat(dat, batch = bat, mod = mod1))
  } else dat_c <- dat
  comb <- cbind(dat_c, tst)
  cb   <- c(rep(1L, ncol(dat_c)), rep(2L, ncol(tst)))
  if (step2_mode == "meanonly") {
    out <- suppressMessages(ComBat(comb, batch = cb, mod = NULL, ref.batch = 1L, mean.only = TRUE))
  } else {
    mod2 <- if (sup2) model.matrix(~ c(lab, ltst)) else NULL
    out  <- suppressMessages(ComBat(comb, batch = cb, mod = mod2, ref.batch = 1L))
  }
  list(ref = dat_c, tgt = out[, (ncol(dat_c) + 1):ncol(out), drop = FALSE])
}

# ── Production KNN: top-1000 HVG + class::knn, returns balanced_acc and MCC ─────
prod_knn <- function(ref, tgt, ref_lab, tgt_lab, n_genes = 1000) {
  cg  <- intersect(rownames(ref), rownames(tgt)); ref <- ref[cg, ]; tgt <- tgt[cg, ]
  hvg <- order(rowVars(ref), decreasing = TRUE)[seq_len(min(n_genes, nrow(ref)))]
  fit  <- predKNN_pp(ref[hvg, , drop = FALSE], ref_lab)               # CV-selects k
  pred <- predWrapper(fit$mod, tgt[hvg, , drop = FALSE], "knn")        # 0/1 vector
  tp <- sum(pred == 1 & tgt_lab == 1); tn <- sum(pred == 0 & tgt_lab == 0)
  fp <- sum(pred == 1 & tgt_lab == 0); fn <- sum(pred == 0 & tgt_lab == 1)
  tpr <- if ((tp + fn) > 0) tp / (tp + fn) else NA
  tnr <- if ((tn + fp) > 0) tn / (tn + fp) else NA
  den <- sqrt(as.double(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (den > 0) (tp * tn - fp * fn) / den else 0
  list(bacc = mean(c(tpr, tnr), na.rm = TRUE), mcc = mcc, k = fit$mod$k, hvg = cg[hvg])
}

# Test class-axis separation (d') in gene space, axis = reference class-mean diff
axis_dprime <- function(ref, tgt, ref_lab, tgt_lab) {
  cg <- intersect(rownames(ref), rownames(tgt)); ref <- ref[cg, ]; tgt <- tgt[cg, ]
  w  <- rowMeans(ref[, ref_lab == 1, drop = FALSE]) - rowMeans(ref[, ref_lab == 0, drop = FALSE])
  w  <- w / sqrt(sum(w^2))
  s  <- as.numeric(t(tgt) %*% w)
  sp <- sqrt(((sum(tgt_lab==1)-1)*var(s[tgt_lab==1]) + (sum(tgt_lab==0)-1)*var(s[tgt_lab==0])) /
             (length(s) - 2))
  (mean(s[tgt_lab == 1]) - mean(s[tgt_lab == 0])) / max(sp, 1e-8)
}

prep <- function(n, test) {
  refs <- train_studies(n, test)
  cg   <- Reduce(intersect, lapply(c(refs, test), function(s) rownames(dat_lst[[s]])))
  dat  <- log_safe(do.call(cbind, lapply(refs, function(s) dat_lst[[s]][cg, ])))
  bat  <- do.call(c, lapply(refs, function(s) rep(s, ncol(dat_lst[[s]]))))
  lab  <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))
  tst  <- log_safe(dat_lst[[test]][cg, ]); ltst <- bin(label_lst[[test]])
  list(dat = dat, bat = bat, lab = lab, tst = tst, ltst = ltst)
}

cat("══════════════════════════════════════════════════════════════════════\n")
cat("  EXACT mechanism: combat_sup x KNN in the production gene-space metric\n")
cat("══════════════════════════════════════════════════════════════════════\n")

# ══ SECTION 1: reproduce the failure with the real production KNN ══════════════
cat("\n── 1. Production KNN (1000-HVG, class::knn) per scenario ──────────────\n")
methods <- list(
  combat        = function(p) two_step(p$dat,p$bat,p$lab,p$tst,p$ltst, FALSE, FALSE, "full"),
  combat_sup    = function(p) two_step(p$dat,p$bat,p$lab,p$tst,p$ltst, TRUE,  TRUE,  "full"),
  sup1_only     = function(p) two_step(p$dat,p$bat,p$lab,p$tst,p$ltst, TRUE,  FALSE, "full"),
  sup_meanonly  = function(p) two_step(p$dat,p$bat,p$lab,p$tst,p$ltst, TRUE,  TRUE,  "meanonly"))

rows <- list()
for (n in 2:5) for (test in ALL_STUDIES) {
  p <- prep(n, test)
  for (mn in names(methods)) {
    fit <- tryCatch(methods[[mn]](p), error = function(e) NULL); if (is.null(fit)) next
    k   <- prod_knn(fit$ref, fit$tgt, p$lab, p$ltst)
    dp  <- axis_dprime(fit$ref, fit$tgt, p$lab, p$ltst)
    rows[[length(rows)+1]] <- data.frame(n=n, test=test, method=mn,
                                         bacc=k$bacc, mcc=k$mcc, k=k$k, axis_dp=dp)
  }
  cat(sprintf("  done n=%d test=%s\n", n, test))
}
R <- do.call(rbind, rows)
write.csv(R, file.path(OUT_DIR, "mechanism_exact_knn.csv"), row.names = FALSE)
cat("\n  Medians by method (24 scenarios):\n")
agg <- do.call(rbind, lapply(names(methods), function(mn) {
  s <- R[R$method == mn, ]
  data.frame(method=mn, bacc=round(median(s$bacc),3), mcc=round(median(s$mcc),3),
             axis_dp=round(median(s$axis_dp),3),
             n_inverted=sum(s$axis_dp < 0, na.rm=TRUE))
}))
print(agg, row.names = FALSE)
cat("  (axis_dp<0 = test class ordering INVERTED on the reference class axis)\n")

# ══ SECTION 2: the per-gene step-2 scaling, on a focal scenario ════════════════
cat("\n── 2. Per-gene step-2 scaling vs class-discriminability (n=3, USA) ────\n")
p <- prep(3, "USA")
fit_c <- two_step(p$dat,p$bat,p$lab,p$tst,p$ltst, FALSE, FALSE, "full")
fit_s <- two_step(p$dat,p$bat,p$lab,p$tst,p$ltst, TRUE,  TRUE,  "full")
cg <- Reduce(intersect, list(rownames(p$tst), rownames(fit_c$tgt), rownames(fit_s$tgt),
                             rownames(fit_c$ref), rownames(fit_s$ref)))
raw_t <- p$tst[cg, ]
# Net per-gene scaling the TEST underwent in step 2 = sd(out)/sd(raw test in)
s_combat <- rowSds(fit_c$tgt[cg, ]) / pmax(rowSds(raw_t), 1e-8)
s_sup    <- rowSds(fit_s$tgt[cg, ]) / pmax(rowSds(raw_t), 1e-8)
# Per-gene class discriminability, measured on the adjusted reference
disc <- function(ref) {
  m1 <- rowMeans(ref[cg, p$lab==1, drop=FALSE]); m0 <- rowMeans(ref[cg, p$lab==0, drop=FALSE])
  v1 <- rowVars(ref[cg, p$lab==1, drop=FALSE]);  v0 <- rowVars(ref[cg, p$lab==0, drop=FALSE])
  n1 <- sum(p$lab==1); n0 <- sum(p$lab==0)
  sw <- sqrt(((n1-1)*v1 + (n0-1)*v0)/(n1+n0-2))
  abs(m1 - m0) / pmax(sw, 1e-8)
}
disc_c <- disc(fit_c$ref); disc_s <- disc(fit_s$ref)
# Reference within-class SD per gene (the var_train_g that supervised step1 shrinks)
within_sd <- function(ref) {
  v1 <- rowVars(ref[cg, p$lab==1, drop=FALSE]); v0 <- rowVars(ref[cg, p$lab==0, drop=FALSE])
  n1 <- sum(p$lab==1); n0 <- sum(p$lab==0); sqrt(((n1-1)*v1 + (n0-1)*v0)/(n1+n0-2))
}
wsd_c <- within_sd(fit_c$ref); wsd_s <- within_sd(fit_s$ref)

top <- order(disc_c, decreasing = TRUE)[1:100]   # 100 most class-discriminative genes
cat(sprintf("  Reference within-class SD on top-100 class genes:  combat %.3f  combat_sup %.3f  (ratio %.2f)\n",
            median(wsd_c[top]), median(wsd_s[top]), median(wsd_s[top])/median(wsd_c[top])))
cat(sprintf("  Step-2 net scaling s_g on those genes:             combat %.3f  combat_sup %.3f\n",
            median(s_combat[top]), median(s_sup[top])))
cat(sprintf("  Step-2 net scaling s_g on all genes (median):      combat %.3f  combat_sup %.3f\n",
            median(s_combat), median(s_sup)))
cat(sprintf("  Spearman rho(class-discriminability, step-2 scaling): combat %+.3f  combat_sup %+.3f\n",
            cor(disc_c, s_combat, method="spearman"), cor(disc_s, s_sup, method="spearman")))
# HVG retention: do the class genes survive the 1000-HVG filter?
hvg_c <- order(rowVars(fit_c$ref[cg,]), decreasing=TRUE)[1:1000]
hvg_s <- order(rowVars(fit_s$ref[cg,]), decreasing=TRUE)[1:1000]
cat(sprintf("  Top-100 class genes retained in 1000-HVG set:      combat %d/100  combat_sup %d/100\n",
            sum(top %in% hvg_c), sum(top %in% hvg_s)))
write.csv(data.frame(gene=cg, disc_combat=disc_c, disc_sup=disc_s,
                     s_combat=s_combat, s_sup=s_sup, wsd_combat=wsd_c, wsd_sup=wsd_s),
          file.path(OUT_DIR, "mechanism_exact_pergene.csv"), row.names = FALSE)

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("  Read: if combat_sup crushes (small s_g) exactly the high-discriminability\n")
cat("  genes (negative rho) while combat does not, that IS the mechanism; the\n")
cat("  sup_meanonly row recovering bacc/mcc confirms step-2 per-gene delta is\n")
cat("  the executor, supervised step-1 the enabler.\n")
cat("══════════════════════════════════════════════════════════════════════\n")
