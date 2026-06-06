#!/usr/bin/env Rscript
# H4.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

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

