#!/usr/bin/env Rscript
# H14.R — auto-generated from evaluate_hypotheses.R; sources shared setup then runs one hypothesis.
if (!exists("sc")) source("scripts/hypotheses_common.R")

cat("\n── H14: Africa underperforms due to per-gene variance-scale difference ─\n")
# Test A: Compare per-gene SD distribution of Africa vs other studies (KS test).
# Test B: Pre-normalise all studies to equal per-gene SD, re-run combat_sup_nat,
#         check if Africa purity gap closes.

cg14 <- Reduce(intersect, lapply(ALL_STUDIES, function(s) rownames(dat_lst[[s]])))

study_sd <- sapply(ALL_STUDIES, function(s) {
  d <- dat_lst[[s]][cg14, ]
  if (max(d, na.rm=TRUE) > 100) d <- log2(d + 1)
  median(rowSds(d), na.rm=TRUE)
})
cat("  Median per-gene SD per study:\n")
print(round(study_sd, 4))

africa_sd <- rowSds(log2(dat_lst[["Africa"]][cg14, ] + 1))
other_sds <- unlist(lapply(setdiff(ALL_STUDIES, "Africa"), function(s) {
  d <- dat_lst[[s]][cg14, ]
  if (max(d, na.rm=TRUE) > 100) d <- log2(d+1)
  rowSds(d)
}))
ks_h14 <- ks.test(africa_sd, other_sds)
cat(sprintf("  KS test Africa vs others per-gene SD: D=%.4f, p=%.3g\n",
            ks_h14$statistic, ks_h14$p.value))

# Test B: normalise each study's variance to its global median SD, then re-run
# for the most discriminating scenario: n=5, test=Africa
sc_af <- prep(5, "Africa")
# Scale each batch to have the same median per-gene SD as the reference batch
ref_batch_sd <- median(rowSds(sc_af$dat[, sc_af$bat == ALL_STUDIES[1], drop=FALSE]), na.rm=TRUE)
dat_scaled <- sc_af$dat
for (b in unique(sc_af$bat)) {
  idx_b  <- which(sc_af$bat == b)
  b_sd   <- median(rowSds(sc_af$dat[, idx_b, drop=FALSE]), na.rm=TRUE)
  dat_scaled[, idx_b] <- sc_af$dat[, idx_b] * (ref_batch_sd / max(b_sd, 1e-6))
}
tst_scaled <- sc_af$tst * (ref_batch_sd / max(median(rowSds(sc_af$tst), na.rm=TRUE), 1e-6))

# Run nat on original and scaled
res_orig   <- run_combat_sup_nat(sc_af$dat, sc_af$bat, sc_af$lab, sc_af$tst)
res_scaled <- run_combat_sup_nat(dat_scaled, sc_af$bat, sc_af$lab, tst_scaled)
res_combat <- suppressMessages({
  dat_c <- ComBat(sc_af$dat, batch=sc_af$bat, mod=NULL)
  comb  <- cbind(dat_c, sc_af$tst)
  cb    <- c(rep(1L, ncol(dat_c)), rep(2L, ncol(sc_af$tst)))
  out   <- ComBat(comb, batch=cb, mod=NULL, ref.batch=1L)
  list(ref=dat_c, tgt=out[,(ncol(dat_c)+1):ncol(out),drop=FALSE])
})

labs5_r <- do.call(c, lapply(train_studies(5, "Africa"), function(s) bin(label_lst[[s]])))
labs5_t <- bin(label_lst[["Africa"]])
pur_nat_orig   <- knn_pur(res_orig$ref,   labs5_r, res_orig$tgt,   labs5_t)
pur_nat_scaled <- knn_pur(res_scaled$ref, labs5_r, res_scaled$tgt, labs5_t)
pur_combat     <- knn_pur(res_combat$ref, labs5_r, res_combat$tgt, labs5_t)

cat(sprintf("  n5/Africa KNN purity: combat=%.3f nat=%.3f nat_var-normalised=%.3f\n",
            pur_combat, pur_nat_orig, pur_nat_scaled))
gap_before <- pur_combat - pur_nat_orig
gap_after  <- pur_combat - pur_nat_scaled
cat(sprintf("  Gap (combat - nat): before normalisation=%.3f, after=%.3f\n",
            gap_before, gap_after))

h14_pass <- ks_h14$p.value < 0.001 && gap_after < gap_before * 0.5
if (h14_pass) {
  add("H14", "TRUE",
      sprintf("Africa SD differs from others (KS p=%.3g); variance normalisation closes gap %.3f→%.3f",
              ks_h14$p.value, gap_before, gap_after))
} else {
  add("H14", "FALSE or PARTIAL",
      sprintf("KS p=%.3g; gap before=%.3f after=%.3f", ks_h14$p.value, gap_before, gap_after))
}

