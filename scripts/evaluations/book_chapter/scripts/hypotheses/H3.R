#!/usr/bin/env Rscript
# H3.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

cat("\n── H3 (old): ComBat_nat inflates training variance 5–7× vs test ───────\n")
# Test: per-gene SD ratio (ComBat_nat training output) / (raw test).
# H3 confirmed (FALSE for the fix) if ratio >> 1.

tr_nat_raw_scale <- ComBat_nat(sc$dat, batch = sc$bat, mod = model.matrix(~sc$lab))
sd_train_nat <- rowSds(tr_nat_raw_scale)
sd_test_raw  <- rowSds(sc$tst)
cg3 <- intersect(names(sd_train_nat), names(sd_test_raw))
scale_ratio  <- sd_train_nat[cg3] / pmax(sd_test_raw[cg3], 1e-10)
med_r3 <- median(scale_ratio, na.rm = TRUE)
q3     <- quantile(scale_ratio, c(0.25, 0.75), na.rm = TRUE)

cat(sprintf("  Per-gene SD ratio ComBat_nat_train / raw_test: median=%.2f [%.2f, %.2f]\n",
            med_r3, q3[1], q3[2]))

if (med_r3 > 3) {
  add("H3 (old)", "CONFIRMED", sprintf("Scale inflation ratio=%.2f (predicted 5-7x)", med_r3))
} else {
  add("H3 (old)", "FALSIFIED", sprintf("Scale inflation ratio=%.2f — not as inflated as claimed", med_r3))
}

