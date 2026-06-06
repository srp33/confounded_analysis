#!/usr/bin/env Rscript
# H2.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

cat("\n── H2 (old): Step-2 sees inflated test variance → large delta.hat ─────\n")
# Test: compute delta.hat for the test batch manually under supervised step-2,
# using (a) supervised step-1 training vs (b) unsupervised step-1 training.
# H2 TRUE => delta.hat much larger for (a). H2 FALSE => delta.hat ≈ 1 in both.

compute_delta_hat <- function(train_mat, train_lab, test_mat, test_lab) {
  cg   <- intersect(rownames(train_mat), rownames(test_mat))
  comb <- cbind(train_mat[cg, ], test_mat[cg, ])
  bat2 <- factor(c(rep(1, ncol(train_mat)), rep(2, ncol(test_mat))))
  cls2 <- c(train_lab, test_lab)
  design2 <- model.matrix(~0 + bat2 + cls2)
  B <- tryCatch(solve(crossprod(design2), tcrossprod(t(design2), comb)),
                error = function(e) NULL)
  if (is.null(B)) return(rep(NA, nrow(comb)))
  ref_idx <- seq_len(ncol(train_mat))
  vp  <- pmax(rowMeans((comb[, ref_idx] - t(design2[ref_idx, ] %*% B))^2), 1e-10)
  tst_idx <- (ncol(train_mat) + 1):ncol(comb)
  s_tst   <- (comb[, tst_idx] - t(design2[tst_idx, ] %*% B)) / sqrt(vp)
  rowVars(s_tst)
}

dh_sup   <- compute_delta_hat(tr_sup,   sc$lab, sc$tst, sc$ltst)
dh_unsup <- compute_delta_hat(tr_unsup, sc$lab, sc$tst, sc$ltst)

med_sup   <- median(dh_sup,   na.rm = TRUE)
med_unsup <- median(dh_unsup, na.rm = TRUE)
cat(sprintf("  Median delta.hat test batch: unsup_step1=%.3f  sup_step1=%.3f\n",
            med_unsup, med_sup))
# Fraction of genes with delta.hat > 2
f_sup   <- mean(dh_sup   > 2, na.rm = TRUE)
f_unsup <- mean(dh_unsup > 2, na.rm = TRUE)
cat(sprintf("  Fraction delta.hat > 2: unsup=%.3f  sup=%.3f\n", f_unsup, f_sup))

# H2 requires delta.hat to be substantially > 1 under supervised step-1
if (med_sup > 2 && med_sup > 2 * med_unsup) {
  add("H2 (old)", "TRUE (unexpectedly)", sprintf("delta.hat sup=%.2f >> unsup=%.2f", med_sup, med_unsup))
} else {
  add("H2 (old)", "FALSE", sprintf("delta.hat sup=%.3f vs unsup=%.3f — no inflation", med_sup, med_unsup))
}

