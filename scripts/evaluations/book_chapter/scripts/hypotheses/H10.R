#!/usr/bin/env Rscript
# H10.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

cat("\n── H10: limma and ComBat mean-only step-1 are near-identical ───────────\n")
# Test: L2 norm of difference between nat and mo corrected matrices,
# and Pearson correlation, across all scenarios.

h10_rows <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    ref_nat <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_reference.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    ref_mo  <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_mo_n%d_test%s_reference.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    if (is.null(ref_nat) || is.null(ref_mo)) next
    cg <- intersect(rownames(ref_nat), rownames(ref_mo))
    # Align columns if needed
    sc_cols <- intersect(colnames(ref_nat), colnames(ref_mo))
    if (length(sc_cols) == 0) next
    diff_mat <- ref_nat[cg, sc_cols] - ref_mo[cg, sc_cols]
    l2_rel   <- sqrt(mean(diff_mat^2)) / sqrt(mean(ref_nat[cg, sc_cols]^2))
    r        <- cor(as.vector(ref_nat[cg, sc_cols]), as.vector(ref_mo[cg, sc_cols]))
    h10_rows[[length(h10_rows)+1]] <- data.frame(n=n, test=test, l2_rel=l2_rel, r=r)
  }
}
h10 <- do.call(rbind, h10_rows)
cat(sprintf("  Relative L2 difference nat vs mo: median=%.5f [%.5f, %.5f]\n",
            median(h10$l2_rel, na.rm=TRUE),
            quantile(h10$l2_rel, 0.25, na.rm=TRUE),
            quantile(h10$l2_rel, 0.75, na.rm=TRUE)))
cat(sprintf("  Pearson r between nat and mo matrices: median=%.6f\n",
            median(h10$r, na.rm=TRUE)))

if (median(h10$l2_rel, na.rm=TRUE) < 0.02 && median(h10$r, na.rm=TRUE) > 0.999) {
  add("H10", "TRUE",
      sprintf("Relative L2=%.5f, r=%.6f — near-identical outputs",
              median(h10$l2_rel, na.rm=TRUE), median(h10$r, na.rm=TRUE)))
} else {
  add("H10", "FALSE",
      sprintf("Relative L2=%.5f, r=%.6f — meaningful differences",
              median(h10$l2_rel, na.rm=TRUE), median(h10$r, na.rm=TRUE)))
}

