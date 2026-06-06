#!/usr/bin/env Rscript
# H11.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

cat("\n── H11: k=5 gap (combat vs nat) correlates with variance heterogeneity ─\n")
# Test: Spearman correlation between cross-batch variance heterogeneity
# (CV of per-batch median per-gene SD) and (combat purity - nat purity) at k=5.

h11_rows <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    sc_i <- prep(n, test)
    # Variance heterogeneity: CV of per-batch per-gene SD
    batch_sds <- sapply(unique(sc_i$bat), function(b) {
      median(rowSds(sc_i$dat[, sc_i$bat == b, drop=FALSE]), na.rm=TRUE)
    })
    cv_het <- sd(batch_sds) / mean(batch_sds)

    ref_c <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_n%d_test%s_reference.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    tgt_c <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_n%d_test%s_target.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    ref_nat <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_reference.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    tgt_nat <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_nat_n%d_test%s_target.csv", n, test)),
      row.names=1)), error=function(e) NULL)
    if (any(sapply(list(ref_c, tgt_c, ref_nat, tgt_nat), is.null))) next

    labs_r <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
    labs_t <- bin(label_lst[[test]])
    pur_c   <- knn_pur(ref_c,   labs_r, tgt_c,   labs_t)
    pur_nat <- knn_pur(ref_nat, labs_r, tgt_nat, labs_t)
    gap     <- pur_c - pur_nat   # positive = combat is better

    h11_rows[[length(h11_rows)+1]] <- data.frame(n=n, test=test,
                                                  cv_het=cv_het, gap=gap)
  }
}
h11 <- do.call(rbind, h11_rows)
cor_h11 <- cor(h11$cv_het, h11$gap, method="spearman", use="complete.obs")
cat("  Cross-batch variance heterogeneity vs purity gap (combat - nat) at k=5:\n")
print(h11, row.names=FALSE, digits=3)
cat(sprintf("  Spearman rho=%.3f\n", cor_h11))

if (cor_h11 > 0.35) {
  add("H11", "TRUE",
      sprintf("Spearman rho=%.3f — variance heterogeneity predicts gap", cor_h11))
} else {
  add("H11", "FALSE or PARTIAL",
      sprintf("Spearman rho=%.3f — no strong correlation", cor_h11))
}

