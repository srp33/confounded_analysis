#!/usr/bin/env Rscript
# H7.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

cat("\n── H7: Mean-only step-2 cannot invert class ordering ──────────────────\n")
# Test: check that class-axis separation for combat_sup_nat and combat_sup_mo
# is positive in ALL 18 scenarios (from precomputed data).

h7_rows <- list()
for (m in c("combat_sup_nat", "combat_sup_mo")) {
  for (n in 2:5) {
    for (test in ALL_STUDIES) {
      ref <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("%s_n%d_test%s_reference.csv", m, n, test)),
        row.names=1)), error=function(e) NULL)
      tgt <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("%s_n%d_test%s_target.csv", m, n, test)),
        row.names=1)), error=function(e) NULL)
      if (is.null(ref) || is.null(tgt)) next
      labs_r <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
      labs_t <- bin(label_lst[[test]])
      sep    <- cls_sep_test(ref, labs_r, tgt, labs_t)
      h7_rows[[length(h7_rows)+1]] <- data.frame(method=m, n=n, test=test, sep=sep)
    }
  }
}
h7 <- do.call(rbind, h7_rows)
n_neg   <- sum(h7$sep <= 0)
n_total <- nrow(h7)
cat(sprintf("  Negative class-axis separations: %d / %d\n", n_neg, n_total))
if (n_neg > 0) cat("  Negative cases:\n")
for (i in which(h7$sep <= 0))
  cat(sprintf("    %s n=%d test=%s sep=%.4f\n",
              h7$method[i], h7$n[i], h7$test[i], h7$sep[i]))

# Also verify mathematically: mean.only step applies a constant shift
# → cannot change ordering. Check that the shift is indeed constant across samples.
step2_shifts <- {
  comb  <- cbind(tr_nat, sc$tst)
  cb    <- factor(c(rep(1, ncol(tr_nat)), rep(2, ncol(sc$tst))))
  out   <- suppressMessages(
    ComBat(comb, batch=cb, mod=NULL, ref.batch=1L, mean.only=TRUE))
  tgt_out <- out[, (ncol(tr_nat)+1):ncol(out), drop=FALSE]
  shift   <- tgt_out - sc$tst[rownames(tgt_out), ]
  # Check: does shift vary across samples?
  col_sds <- apply(shift, 1, sd)
  median(col_sds)
}
cat(sprintf("  Median per-gene SD of sample-to-sample shift variation: %.6f\n",
            step2_shifts))

if (n_neg == 0 && step2_shifts < 0.001) {
  add("H7", "TRUE",
      sprintf("All %d mean-only method scores positive; shift is constant (SD=%.6f)",
              n_total, step2_shifts))
} else {
  add("H7", "FALSE or PARTIAL",
      sprintf("%d/%d negative; shift SD=%.6f", n_neg, n_total, step2_shifts))
}

