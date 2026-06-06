#!/usr/bin/env Rscript
# H8.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

cat("\n── H8: Supervised step-1 distorts gene-gene covariance vs test ─────────\n")
# Test: L2 distance of per-gene SD profiles:
# dist(sup_train, test_raw) vs dist(unsup_train, test_raw) vs dist(nat_train, test_raw)
# H8 TRUE: sup_train further from test than unsup_train and nat_train.

cg8    <- intersect(rownames(tr_sup), rownames(sc$tst))
sd_sup  <- rowSds(tr_sup[cg8, ])
sd_uns  <- rowSds(tr_unsup[cg8, ])
sd_nat  <- rowSds(tr_nat[cg8, ])
sd_mo   <- rowSds(tr_mo_cls[cg8, ])
sd_tst  <- rowSds(sc$tst[cg8, ])

d_sup <- sqrt(mean((sd_sup - sd_tst)^2))
d_uns <- sqrt(mean((sd_uns - sd_tst)^2))
d_nat <- sqrt(mean((sd_nat - sd_tst)^2))
d_mo  <- sqrt(mean((sd_mo  - sd_tst)^2))

cat(sprintf("  RMSE of per-gene SD vs raw test:\n"))
cat(sprintf("    unsup=%.4f  sup=%.4f  nat(mean-only)=%.4f  mo(mean-only)=%.4f\n",
            d_uns, d_sup, d_nat, d_mo))

# Also test variance covariance structure via first PC variance explained
pca_sup <- prcomp(t(tr_sup[cg8, ]), center=TRUE, scale.=FALSE, rank.=10)
pca_uns <- prcomp(t(tr_unsup[cg8, ]), center=TRUE, scale.=FALSE, rank.=10)
pca_nat <- prcomp(t(tr_nat[cg8, ]), center=TRUE, scale.=FALSE, rank.=10)
pca_tst <- prcomp(t(sc$tst[cg8, ]), center=TRUE, scale.=FALSE, rank.=10)

pct_sup <- pca_sup$sdev[1]^2 / sum(pca_sup$sdev^2)
pct_uns <- pca_uns$sdev[1]^2 / sum(pca_uns$sdev^2)
pct_nat <- pca_nat$sdev[1]^2 / sum(pca_nat$sdev^2)
pct_tst <- pca_tst$sdev[1]^2 / sum(pca_tst$sdev^2)

cat(sprintf("  PC1 variance explained: test=%.3f unsup=%.3f sup=%.3f nat=%.3f\n",
            pct_tst, pct_uns, pct_sup, pct_nat))

if (d_sup > d_uns * 1.2 && d_sup > d_nat * 1.2) {
  add("H8", "TRUE",
      sprintf("SD profile dist: sup=%.4f > unsup=%.4f, nat=%.4f", d_sup, d_uns, d_nat))
} else {
  add("H8", "FALSE or PARTIAL",
      sprintf("SD profile dist: sup=%.4f, unsup=%.4f, nat=%.4f", d_sup, d_uns, d_nat))
}

