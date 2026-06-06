#!/usr/bin/env Rscript
# H6.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

cat("\n── H6: combat_sup fails due to step-2 batch-class confounding ──────────\n")
# Test A: correlation of (condition number of step-2 design) with purity loss.
# Test B: subsample test to 50/50 balance; rerun combat_sup; compare purity.

h6_rows <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    sc_i <- prep(n, test)
    # Condition number of step-2 supervised design matrix
    cls_all <- c(sc_i$lab, sc_i$ltst)
    bat_all <- factor(c(rep(1, ncol(sc_i$dat)), rep(2, ncol(sc_i$tst))))
    des     <- model.matrix(~0 + bat_all + cls_all)
    kappa_  <- kappa(des, exact = TRUE)
    # Test class proportion
    prop_active <- mean(sc_i$ltst)
    h6_rows[[length(h6_rows)+1]] <- data.frame(
      n=n, test=test, kappa=kappa_, prop_active=prop_active)
  }
}
h6 <- do.call(rbind, h6_rows)

# Load purity at k=5 from precomputed data (H6: purity loss ~ condition number)
purity_sup <- sapply(seq_len(nrow(h6)), function(i) {
  ref <- tryCatch(as.matrix(read.csv(
    file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_reference.csv", h6$n[i], h6$test[i])),
    row.names=1)), error=function(e) NULL)
  tgt <- tryCatch(as.matrix(read.csv(
    file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_target.csv",    h6$n[i], h6$test[i])),
    row.names=1)), error=function(e) NULL)
  if (is.null(ref) || is.null(tgt)) return(NA)
  labs_r <- do.call(c, lapply(train_studies(h6$n[i], h6$test[i]), function(s) bin(label_lst[[s]])))
  labs_t <- bin(label_lst[[h6$test[i]]])
  knn_pur(ref, labs_r, tgt, labs_t, k=5)
})
purity_nat <- sapply(seq_len(nrow(h6)), function(i) {
  ref <- tryCatch(as.matrix(read.csv(
    file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_reference.csv", h6$n[i], h6$test[i])),
    row.names=1)), error=function(e) NULL)
  tgt <- tryCatch(as.matrix(read.csv(
    file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_target.csv",    h6$n[i], h6$test[i])),
    row.names=1)), error=function(e) NULL)
  if (is.null(ref) || is.null(tgt)) return(NA)
  labs_r <- do.call(c, lapply(train_studies(h6$n[i], h6$test[i]), function(s) bin(label_lst[[s]])))
  labs_t <- bin(label_lst[[h6$test[i]]])
  knn_pur(ref, labs_r, tgt, labs_t, k=5)
})

h6$purity_sup  <- purity_sup
h6$purity_nat  <- purity_nat
h6$purity_loss <- purity_nat - purity_sup   # positive = sup is worse

# Correlation: condition number vs purity loss
cc_h6 <- cor(log(h6$kappa), h6$purity_loss, use="complete.obs", method="spearman")
cat(sprintf("  Spearman cor(log(kappa), purity_loss)=%.3f\n", cc_h6))

# Test B: re-run combat_sup on balanced test (50/50) for primary scenario
set.seed(42)
ltst0 <- sc$ltst
idx0  <- which(ltst0 == 0); idx1 <- which(ltst0 == 1)
n_bal <- min(length(idx0), length(idx1))
idx_b <- c(sample(idx0, n_bal), sample(idx1, n_bal))
tst_bal  <- sc$tst[, idx_b, drop=FALSE]
ltst_bal <- ltst0[idx_b]

res_nat_bal  <- run_combat_sup_nat(sc$dat, sc$bat, sc$lab, tst_bal)
res_nat_full <- run_combat_sup_nat(sc$dat, sc$bat, sc$lab, sc$tst)
res_sup_bal  <- run_combat_sup(sc$dat, sc$bat, sc$lab, tst_bal, ltst_bal)
res_sup_full <- run_combat_sup(sc$dat, sc$bat, sc$lab, sc$tst, sc$ltst)

pur_sup_natural  <- knn_pur(res_sup_full$ref,  sc$lab,  res_sup_full$tgt,  sc$ltst)
pur_sup_balanced <- knn_pur(res_sup_bal$ref,   sc$lab,  res_sup_bal$tgt,   ltst_bal)
pur_nat_natural  <- knn_pur(res_nat_full$ref,  sc$lab,  res_nat_full$tgt,  sc$ltst)
pur_nat_balanced <- knn_pur(res_nat_bal$ref,   sc$lab,  res_nat_bal$tgt,   ltst_bal)

cat(sprintf("  n3/USA combat_sup purity: natural=%.3f  balanced=%.3f (Δ=%.3f)\n",
            pur_sup_natural, pur_sup_balanced, pur_sup_balanced - pur_sup_natural))
cat(sprintf("  n3/USA combat_sup_nat purity: natural=%.3f  balanced=%.3f (Δ=%.3f)\n",
            pur_nat_natural, pur_nat_balanced, pur_nat_balanced - pur_nat_natural))
cat(sprintf("  Natural test class proportion (Active): %.2f\n", mean(sc$ltst)))

h6_pass <- (pur_sup_balanced - pur_sup_natural > 0.10) && (cc_h6 > 0.3)
if (h6_pass) {
  add("H6", "TRUE",
      sprintf("Balanced test: purity %.3f→%.3f (+%.3f); kappa corr=%.3f",
              pur_sup_natural, pur_sup_balanced,
              pur_sup_balanced - pur_sup_natural, cc_h6))
} else {
  add("H6", "FALSE or PARTIAL",
      sprintf("Balanced Δpurity=%.3f; kappa corr=%.3f",
              pur_sup_balanced - pur_sup_natural, cc_h6))
}

