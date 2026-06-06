#!/usr/bin/env Rscript
# H1.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

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

