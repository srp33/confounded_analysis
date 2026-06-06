#!/usr/bin/env Rscript
# H13.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

cat("\n── H13: v2 degradation is monotone in n ────────────────────────────────\n")
# Test: (nat_v1 - nat_v2) purity at k=5, grouped by n.
# H13 TRUE: degradation increases with n (Kruskal-Wallis + trend test).

h13_rows <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    for (v in c("nat","mo")) {
      r1 <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("combat_sup_%s_n%d_test%s_reference.csv", v, n, test)),
        row.names=1)), error=function(e) NULL)
      t1 <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("combat_sup_%s_n%d_test%s_target.csv", v, n, test)),
        row.names=1)), error=function(e) NULL)
      r2 <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("combat_sup_%s_v2_n%d_test%s_reference.csv", v, n, test)),
        row.names=1)), error=function(e) NULL)
      t2 <- tryCatch(as.matrix(read.csv(
        file.path(ADJ_DIR, sprintf("combat_sup_%s_v2_n%d_test%s_target.csv", v, n, test)),
        row.names=1)), error=function(e) NULL)
      if (any(sapply(list(r1,t1,r2,t2), is.null))) next
      labs_r <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
      labs_t <- bin(label_lst[[test]])
      p1 <- knn_pur(r1, labs_r, t1, labs_t)
      p2 <- knn_pur(r2, labs_r, t2, labs_t)
      h13_rows[[length(h13_rows)+1]] <- data.frame(n=n, test=test, variant=v,
                                                    v1=p1, v2=p2, delta=p1-p2)
    }
  }
}
h13 <- do.call(rbind, h13_rows)
cat("  Mean (v1 - v2) purity by n:\n")
h13_summ <- h13 %>% group_by(n) %>%
  summarise(mean_delta=mean(delta, na.rm=TRUE), n_obs=n(), .groups="drop")
print(h13_summ, digits=3)

# Jonckheere-Terpstra trend test (or just Spearman)
cor_h13 <- cor(h13$n, h13$delta, method="spearman", use="complete.obs")
cat(sprintf("  Spearman rho(n, v1-v2 delta)=%.3f\n", cor_h13))

if (isTRUE(cor_h13 > 0.40) && all(h13_summ$mean_delta >= 0, na.rm = TRUE)) {
  add("H13", "TRUE",
      sprintf("Spearman rho=%.3f; mean deltas: n2=%.3f n3=%.3f n4=%.3f n5=%.3f",
              cor_h13, h13_summ$mean_delta[1], h13_summ$mean_delta[2],
              h13_summ$mean_delta[3], h13_summ$mean_delta[4]))
} else {
  add("H13", "FALSE or PARTIAL",
      sprintf("Spearman rho=%.3f", cor_h13))
}

