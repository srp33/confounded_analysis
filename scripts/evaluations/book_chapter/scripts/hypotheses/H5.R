#!/usr/bin/env Rscript
# H5.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

cat("\n── H5: Supervised step-1 preserves class separation; unsupervised reduces it ─\n")
# Test: compare class separation on the class axis across confounded scenarios.
# Compute ratio sep_sup / sep_unsup from the pre-corrected matrices in ADJ_DIR
# across all 18 scenarios.  H5: ratio > 1 (supervised preserves more).

h5_rows <- list()
for (n in 2:5) {
  for (test in ALL_STUDIES) {
    ref_sup  <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_sup_n%d_test%s_reference.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    ref_uns  <- tryCatch(as.matrix(read.csv(
      file.path(ADJ_DIR, sprintf("combat_n%d_test%s_reference.csv", n, test)),
      row.names = 1)), error = function(e) NULL)
    if (is.null(ref_sup) || is.null(ref_uns)) next
    labs <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
    cg   <- intersect(rownames(ref_sup), rownames(ref_uns))
    sep_s <- {
      w <- rowMeans(ref_sup[cg, labs==1, drop=FALSE]) -
           rowMeans(ref_sup[cg, labs==0, drop=FALSE])
      w <- w / sqrt(sum(w^2))
      s <- as.numeric(t(ref_sup[cg, ]) %*% w)
      mean(s[labs==1]) - mean(s[labs==0])
    }
    sep_u <- {
      w <- rowMeans(ref_uns[cg, labs==1, drop=FALSE]) -
           rowMeans(ref_uns[cg, labs==0, drop=FALSE])
      w <- w / sqrt(sum(w^2))
      s <- as.numeric(t(ref_uns[cg, ]) %*% w)
      mean(s[labs==1]) - mean(s[labs==0])
    }
    h5_rows[[length(h5_rows)+1]] <- data.frame(n=n, test=test,
                                                sep_sup=sep_s, sep_unsup=sep_u,
                                                ratio=sep_s/max(sep_u, 0.01))
  }
}
h5 <- do.call(rbind, h5_rows)
cat("  Class separation on training class axis (supervised vs unsupervised):\n")
print(h5[, c("n","test","sep_sup","sep_unsup","ratio")], row.names = FALSE, digits = 3)

med_h5_ratio <- median(h5$ratio, na.rm = TRUE)
wt_h5 <- wilcox.test(h5$sep_sup, h5$sep_unsup, paired = TRUE, exact = FALSE)
cat(sprintf("  Median ratio sup/unsup=%.2f; Wilcoxon p=%.3g\n",
            med_h5_ratio, wt_h5$p.value))

if (med_h5_ratio > 1.05 && wt_h5$p.value < 0.05) {
  add("H5", "TRUE",
      sprintf("sup class-axis sep > unsup in %d/%d scenarios; median ratio=%.2f (p=%.3g)",
              sum(h5$sep_sup > h5$sep_unsup), nrow(h5), med_h5_ratio, wt_h5$p.value))
} else {
  add("H5", "FALSE",
      sprintf("No consistent advantage: median ratio=%.2f (p=%.3g)", med_h5_ratio, wt_h5$p.value))
}

