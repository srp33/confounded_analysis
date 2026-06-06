#!/usr/bin/env Rscript
# H9.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

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

