#!/usr/bin/env Rscript
# mechanism_why.R
# WHY does supervised step-1 ComBat produce a class axis that fails to generalise
# (collapses/inverts on the held-out test)?  Hypothesis: under batch-class
# confounding, protecting `class` in the model lets ComBat fold batch-specific
# structure into the class term beta. The reference class axis (mu1-mu0) is then
# overfit to the training studies and does not transfer to the test.
#
# Test, per scenario, on the reference class axis w = (mu1-mu0)/||.|| learned on
# the step-1 reference:
#   train_dp = resubstitution d' of the reference itself along w   (fit)
#   test_dp  = d' of the held-out test along the same w            (generalisation)
#   gap = train_dp - test_dp                                       (overfitting)
# and a confounding score = Cramer's V between study identity and class in train.
# Hypothesis: supervised step-1 -> large train_dp, small/negative test_dp (big gap),
# and the gap grows with confounding; unsupervised step-1 -> generalises.

suppressMessages({ library(sva); library(matrixStats) })
DATA_FILE <- "data/TB_real_data.RData"
OUT_DIR   <- "outputs/diagnostics/hypothesis_tests"
load(DATA_FILE)
ALL_STUDIES <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1", 1, "Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES != test][seq_len(n)]
log_safe <- function(m) { if (max(m,na.rm=TRUE)>100){mn<-min(m,na.rm=TRUE); if(mn<0)m<-m-mn; m<-log2(m+1)}; m }

step1 <- function(dat, bat, lab, sup) {
  if (length(unique(bat)) < 2) return(dat)
  suppressMessages(ComBat(dat, batch = bat, mod = if (sup) model.matrix(~lab) else NULL))
}
dprime_along <- function(mat, lab, w) {
  s <- as.numeric(t(mat) %*% w)
  sp <- sqrt(((sum(lab==1)-1)*var(s[lab==1]) + (sum(lab==0)-1)*var(s[lab==0]))/(length(s)-2))
  (mean(s[lab==1]) - mean(s[lab==0])) / max(sp, 1e-8)
}
cramers_v <- function(study, cls) {
  tab <- table(study, cls); ch <- suppressWarnings(chisq.test(tab)$statistic)
  sqrt(as.numeric(ch) / (sum(tab) * (min(dim(tab)) - 1)))
}

rows <- list()
for (n in 2:5) for (test in ALL_STUDIES) {
  refs <- train_studies(n, test)
  cg   <- Reduce(intersect, lapply(c(refs, test), function(s) rownames(dat_lst[[s]])))
  dat  <- log_safe(do.call(cbind, lapply(refs, function(s) dat_lst[[s]][cg, ])))
  bat  <- do.call(c, lapply(refs, function(s) rep(s, ncol(dat_lst[[s]]))))
  lab  <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))
  tst  <- log_safe(dat_lst[[test]][cg, ]); ltst <- bin(label_lst[[test]])
  cv   <- cramers_v(bat, lab)
  for (sup in c(FALSE, TRUE)) {
    ref <- step1(dat, bat, lab, sup)
    w   <- rowMeans(ref[, lab==1, drop=FALSE]) - rowMeans(ref[, lab==0, drop=FALSE])
    w   <- w / sqrt(sum(w^2))
    tr  <- dprime_along(ref, lab, w)
    te  <- dprime_along(tst, ltst, w)         # raw test on the reference-learned axis
    rows[[length(rows)+1]] <- data.frame(n=n, test=test, sup=sup, confound=cv,
                                         train_dp=tr, test_dp=te, gap=tr-te)
  }
  cat(sprintf("  done n=%d test=%s (confound V=%.2f)\n", n, test, cv))
}
D <- do.call(rbind, rows)
write.csv(D, file.path(OUT_DIR, "mechanism_why.csv"), row.names = FALSE)

cat("\n── Medians: reference class axis fit vs generalisation ───────────────\n")
agg <- do.call(rbind, lapply(c(FALSE, TRUE), function(s) {
  x <- D[D$sup == s, ]
  data.frame(step1 = if (s) "supervised" else "unsupervised",
             train_dp = round(median(x$train_dp), 2),
             test_dp  = round(median(x$test_dp),  2),
             gap      = round(median(x$gap),      2),
             n_test_inverted = sum(x$test_dp < 0))
}))
print(agg, row.names = FALSE)

sup <- D[D$sup == TRUE, ]
cat(sprintf("\n  Supervised step-1: Spearman rho(confounding, test_dp) = %+.3f\n",
            cor(sup$confound, sup$test_dp, method = "spearman")))
cat(sprintf("  Supervised step-1: Spearman rho(confounding, gap)     = %+.3f\n",
            cor(sup$confound, sup$gap, method = "spearman")))
cat("\n  Interpretation: large train_dp + small/negative test_dp under supervision\n")
cat("  = the protected class axis is overfit to the training studies; if rho(confound,\n")
cat("  test_dp) is negative, batch-class confounding is what breaks generalisation.\n")
