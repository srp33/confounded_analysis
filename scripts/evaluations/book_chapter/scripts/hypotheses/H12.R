#!/usr/bin/env Rscript
# H12.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

cat("\n── H12: Larger class axis separation predicts k=11 advantage ───────────\n")
# Test: for each scenario × method, correlate training class-axis separation
# with purity at k=11 minus k=5.

h12_rows <- list()
for (m in c("combat", "combat_sup_nat", "combat_sup_mo")) {
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
      pur5   <- knn_pur(ref, labs_r, tgt, labs_t, k=5)
      pur11  <- knn_pur(ref, labs_r, tgt, labs_t, k=11)
      h12_rows[[length(h12_rows)+1]] <- data.frame(
        method=m, n=n, test=test, sep=sep, pur5=pur5, pur11=pur11,
        delta_k=pur11 - pur5)
    }
  }
}
h12 <- do.call(rbind, h12_rows)
cor_h12 <- cor(h12$sep, h12$delta_k, method="spearman", use="complete.obs")
cat(sprintf("  Spearman rho(class-sep, purity[k=11]-purity[k=5])=%.3f\n", cor_h12))

# Also test: nat methods outperform combat at k=11 more when class sep is larger
h12_nat <- h12[h12$method == "combat_sup_nat", ]
h12_com <- h12[h12$method == "combat", ]
key <- paste(h12_nat$n, h12_nat$test)
key2 <- paste(h12_com$n, h12_com$test)
idx  <- match(key, key2)
h12_comp <- data.frame(
  n=h12_nat$n, test=h12_nat$test,
  sep_nat = h12_nat$sep,
  p11_nat = h12_nat$pur11,
  p11_com = h12_com$pur11[idx],
  adv_nat = h12_nat$pur11 - h12_com$pur11[idx])
cor_h12b <- cor(h12_comp$sep_nat, h12_comp$adv_nat, method="spearman")
cat(sprintf("  Spearman rho(nat class-sep, nat advantage at k=11)=%.3f\n", cor_h12b))

if (cor_h12 > 0.30 || cor_h12b > 0.35) {
  add("H12", "TRUE",
      sprintf("sep→(k11-k5) rho=%.3f; sep→nat_adv_k11 rho=%.3f", cor_h12, cor_h12b))
} else {
  add("H12", "FALSE or PARTIAL",
      sprintf("sep→(k11-k5) rho=%.3f; sep→nat_adv_k11 rho=%.3f", cor_h12, cor_h12b))
}

