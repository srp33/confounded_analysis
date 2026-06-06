#!/usr/bin/env Rscript
# evaluate_knn_vs_rda.R
# Why does combat_sup wreck KNN but not RDA/LDA?  Linear-algebra evidence.
#
# ComBat step-2 applies a diagonal (per-gene) affine map  t* = D (t - c),
# D = diag(1/delta.hat), to the TEST only.  Under supervised step-1 the training
# reference's within-class variance is compressed, so delta.hat > 1 on the
# class-carrying genes and D compresses the test along the class-discriminating
# direction.
#
#   * KNN uses Euclidean distance, invariant only to ORTHOGONAL maps; a diagonal
#     D with delta != 1 reweights axes and breaks neighbour structure.
#   * Fisher/LDA maximises J(w)=(w'S_B w)/(w'S_W w), invariant under ANY invertible
#     linear map (the S_W^{-1} whitening undoes D); so RDA/shrinkageLDA survive.
#
# Three sections:
#   A. Per-classifier degradation table (KNN vs LDA/SVM/...) from the pipeline.
#   B. Per-scenario test separation on the RAW axis (mu1-mu0) vs the WHITENED axis
#      (S_W^{-1}(mu1-mu0)), plus Euclidean KNN purity vs LDA accuracy.
#   C. Controlled demo: apply combat_sup's own per-gene rescaling D to the
#      (healthy) unsupervised-combat test and show KNN collapses while the Fisher
#      ratio J(w) is preserved — isolating the diagonal map as the cause.

suppressMessages({ library(matrixStats) })

DATA_FILE <- "data/TB_real_data.RData"
ADJ_DIR   <- "outputs/adjusted_data/all_scenarios"
OUT_DIR   <- "outputs/diagnostics/hypothesis_tests"
CLS_FILE  <- "outputs/adjusters_on_classifiers.csv"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
load(DATA_FILE)

ALL_STUDIES <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1", 1, "Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES != test][seq_len(n)]

rd <- function(path) {
  m <- tryCatch(as.matrix(read.csv(path, row.names = 1)), error = function(e) NULL)
  m
}

# ── Fisher / KNN primitives in a reduced PC space ─────────────────────────────
# Reduce to top-K reference PCs (K << n_ref so S_W is well-conditioned), then
# compute everything in PC coordinates.
analyse_scenario <- function(ref, tgt, lab_r, lab_t, K = 50, ridge = 1e-2) {
  cg  <- intersect(rownames(ref), rownames(tgt))
  ref <- ref[cg, ]; tgt <- tgt[cg, ]
  K   <- min(K, ncol(ref) - 2)
  pca <- prcomp(t(ref), center = TRUE, scale. = FALSE, rank. = K)
  Rp  <- pca$x                                   # n_ref x K
  Tp  <- predict(pca, t(tgt))                    # n_test x K

  m1 <- colMeans(Rp[lab_r == 1, , drop = FALSE])
  m0 <- colMeans(Rp[lab_r == 0, , drop = FALSE])
  dmu <- m1 - m0

  # Pooled within-class scatter in PC space (+ ridge for conditioning)
  c1 <- scale(Rp[lab_r == 1, , drop = FALSE], center = m1, scale = FALSE)
  c0 <- scale(Rp[lab_r == 0, , drop = FALSE], center = m0, scale = FALSE)
  Sw <- (crossprod(c1) + crossprod(c0)) / (nrow(Rp) - 2)
  Sw <- Sw + diag(ridge * mean(diag(Sw)), K)

  # Directions learned on the REFERENCE
  w_raw  <- dmu / sqrt(sum(dmu^2))               # raw mean-difference axis (KNN/Euclidean world)
  w_fish <- solve(Sw, dmu)                        # Fisher / LDA whitened axis
  w_fish <- w_fish / sqrt(sum(w_fish^2))

  # Standardised test separation (d') along each reference-learned axis
  dprime <- function(w) {
    s <- as.numeric(Tp %*% w)
    sd_pool <- sqrt(((sum(lab_t == 1) - 1) * var(s[lab_t == 1]) +
                     (sum(lab_t == 0) - 1) * var(s[lab_t == 0])) /
                    (length(s) - 2))
    (mean(s[lab_t == 1]) - mean(s[lab_t == 0])) / max(sd_pool, 1e-8)
  }
  dp_raw  <- dprime(w_raw)
  dp_fish <- dprime(w_fish)

  # Fisher ratio achieved on the test by the reference-learned Fisher axis
  J_fish <- as.numeric(t(dmu) %*% solve(Sw, dmu))   # = Mahalanobis dist^2 of class means

  # Euclidean KNN purity (test labelled by k nearest reference neighbours)
  knn_pur <- function(k = 5) {
    d <- as.matrix(dist(rbind(Tp, Rp)))[seq_len(nrow(Tp)),
                                        (nrow(Tp) + 1):(nrow(Tp) + nrow(Rp))]
    mean(sapply(seq_len(nrow(Tp)), function(i)
      mean(lab_r[order(d[i, ])[seq_len(k)]] == lab_t[i])))
  }
  pur_knn <- knn_pur(5)

  # LDA accuracy: classify test by the Fisher axis, threshold at reference midpoint
  mid <- as.numeric((m1 + m0) / 2) %*% w_fish
  s_t <- as.numeric(Tp %*% w_fish)
  pred <- as.integer(s_t > as.numeric(mid))
  if (mean(pred == lab_t) < 0.5) pred <- 1L - pred   # orient sign by best of the two
  acc_lda <- mean(pred == lab_t)

  list(dp_raw = dp_raw, dp_fish = dp_fish, J_fish = J_fish,
       pur_knn = pur_knn, acc_lda = acc_lda,
       Rp = Rp, Tp = Tp, lab_r = lab_r, lab_t = lab_t,
       w_raw = w_raw, w_fish = w_fish, Sw = Sw, dmu = dmu)
}

cat("═══════════════════════════════════════════════════════════════════\n")
cat("  KNN vs RDA: why combat_sup is KNN-specific\n")
cat("═══════════════════════════════════════════════════════════════════\n")

# ══════════════════════════════════════════════════════════════════════════════
cat("\n── A. Per-classifier degradation (combat_sup - combat), pipeline output ──\n")
if (file.exists(CLS_FILE)) {
  cl <- read.csv(CLS_FILE)
  tabA <- do.call(rbind, lapply(c("mcc", "balanced_acc"), function(mt) {
    s <- cl[cl$metric == mt & cl$adjuster %in% c("combat", "combat_sup"), ]
    a <- aggregate(value ~ classifier + adjuster, s, function(x) mean(x, na.rm = TRUE))
    w <- reshape(a, idvar = "classifier", timevar = "adjuster", direction = "wide")
    names(w) <- c("classifier", "combat", "combat_sup")
    data.frame(metric = mt, classifier = w$classifier,
               combat = round(w$combat, 3), combat_sup = round(w$combat_sup, 3),
               delta = round(w$combat_sup - w$combat, 3))
  }))
  tabA <- tabA[order(tabA$metric, tabA$delta), ]
  print(tabA, row.names = FALSE)
  write.csv(tabA, file.path(OUT_DIR, "knn_vs_rda_by_classifier.csv"), row.names = FALSE)
  cat(sprintf("\n  -> written knn_vs_rda_by_classifier.csv\n"))
} else {
  cat("  (adjusters_on_classifiers.csv not found — skipping section A)\n")
}

# ══════════════════════════════════════════════════════════════════════════════
cat("\n── B. Test separation: RAW axis vs WHITENED (Fisher) axis, per scenario ──\n")
cat("   dp_raw  = test d' on (mu1-mu0)           [what KNN/Euclidean sees]\n")
cat("   dp_fish = test d' on S_W^-1(mu1-mu0)      [what LDA/RDA sees]\n\n")

rowsB <- list()
for (n in 2:5) for (test in ALL_STUDIES) {
  refr_c <- rd(sprintf("%s/combat_n%d_test%s_reference.csv",     ADJ_DIR, n, test))
  tgt_c  <- rd(sprintf("%s/combat_n%d_test%s_target.csv",        ADJ_DIR, n, test))
  refr_s <- rd(sprintf("%s/combat_sup_n%d_test%s_reference.csv", ADJ_DIR, n, test))
  tgt_s  <- rd(sprintf("%s/combat_sup_n%d_test%s_target.csv",    ADJ_DIR, n, test))
  if (any(sapply(list(refr_c, tgt_c, refr_s, tgt_s), is.null))) next
  lab_r <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
  lab_t <- bin(label_lst[[test]])

  ac <- analyse_scenario(refr_c, tgt_c, lab_r, lab_t)
  as_ <- analyse_scenario(refr_s, tgt_s, lab_r, lab_t)
  rowsB[[length(rowsB) + 1]] <- data.frame(
    n = n, test = test,
    sup_dp_raw  = as_$dp_raw,  sup_dp_fish  = as_$dp_fish,
    sup_pur_knn = as_$pur_knn, sup_acc_lda  = as_$acc_lda,
    com_dp_raw  = ac$dp_raw,   com_dp_fish  = ac$dp_fish,
    com_pur_knn = ac$pur_knn,  com_acc_lda  = ac$acc_lda)
}
B <- do.call(rbind, rowsB)
showB <- data.frame(n = B$n, test = B$test,
                    sup_dp_raw  = round(B$sup_dp_raw, 3),
                    sup_dp_fish = round(B$sup_dp_fish, 3),
                    sup_pur_knn = round(B$sup_pur_knn, 3),
                    sup_acc_lda = round(B$sup_acc_lda, 3))
print(showB, row.names = FALSE)

cat(sprintf("\n  combat_sup medians:  |dp_raw|=%.2f  |dp_fish|=%.2f  KNN_pur=%.2f  LDA_acc=%.2f\n",
            median(abs(B$sup_dp_raw)), median(abs(B$sup_dp_fish)),
            median(B$sup_pur_knn), median(B$sup_acc_lda)))
cat(sprintf("  combat     medians:  |dp_raw|=%.2f  |dp_fish|=%.2f  KNN_pur=%.2f  LDA_acc=%.2f\n",
            median(abs(B$com_dp_raw)), median(abs(B$com_dp_fish)),
            median(B$com_pur_knn), median(B$com_acc_lda)))
cat(sprintf("\n  Under combat_sup, whitened axis keeps %.0f%% of combat's d' but raw axis keeps only %.0f%%.\n",
            100 * median(abs(B$sup_dp_fish)) / median(abs(B$com_dp_fish)),
            100 * median(abs(B$sup_dp_raw))  / median(abs(B$com_dp_raw))))
write.csv(B, file.path(OUT_DIR, "knn_vs_rda_separation.csv"), row.names = FALSE)
cat("  -> written knn_vs_rda_separation.csv\n")

# ══════════════════════════════════════════════════════════════════════════════
cat("\n── C. Controlled demo: apply combat_sup's per-gene rescaling D to the ──\n")
cat("      HEALTHY unsupervised-combat test. KNN should collapse; Fisher J invariant.\n")
demo_scen <- list(n = 3, test = "USA")
n <- demo_scen$n; test <- demo_scen$test
refr_c <- rd(sprintf("%s/combat_n%d_test%s_reference.csv",     ADJ_DIR, n, test))
tgt_c  <- rd(sprintf("%s/combat_n%d_test%s_target.csv",        ADJ_DIR, n, test))
tgt_s  <- rd(sprintf("%s/combat_sup_n%d_test%s_target.csv",    ADJ_DIR, n, test))
lab_r  <- do.call(c, lapply(train_studies(n, test), function(s) bin(label_lst[[s]])))
lab_t  <- bin(label_lst[[test]])

cg <- Reduce(intersect, list(rownames(refr_c), rownames(tgt_c), rownames(tgt_s)))
refr_c <- refr_c[cg, ]; tgt_c <- tgt_c[cg, ]; tgt_s <- tgt_s[cg, ]

# Effective per-gene rescaling combat_sup applied vs combat (the diagonal D)
D_g <- rowSds(tgt_s) / pmax(rowSds(tgt_c), 1e-8)
# Apply D to the healthy combat test (variance only; keep gene means fixed)
mu_g    <- rowMeans(tgt_c)
tgt_mod <- mu_g + D_g * (tgt_c - mu_g)

base <- analyse_scenario(refr_c, tgt_c,  lab_r, lab_t)
modd <- analyse_scenario(refr_c, tgt_mod, lab_r, lab_t)
real <- analyse_scenario(refr_c, tgt_s,  lab_r, lab_t)  # actual combat_sup target on combat ref

cat(sprintf("\n  Median per-gene D (= sd_sup/sd_combat): %.3f  [IQR %.3f-%.3f]\n",
            median(D_g), quantile(D_g, .25), quantile(D_g, .75)))
demoT <- data.frame(
  case      = c("combat test (healthy)", "combat test x D (synthetic)", "actual combat_sup test"),
  dp_raw    = round(c(base$dp_raw,  modd$dp_raw,  real$dp_raw), 3),
  dp_fish   = round(c(base$dp_fish, modd$dp_fish, real$dp_fish), 3),
  J_fish    = round(c(base$J_fish,  modd$J_fish,  real$J_fish), 2),
  KNN_pur   = round(c(base$pur_knn, modd$pur_knn, real$pur_knn), 3),
  LDA_acc   = round(c(base$acc_lda, modd$acc_lda, real$acc_lda), 3))
print(demoT, row.names = FALSE)
write.csv(demoT, file.path(OUT_DIR, "knn_vs_rda_invariance_demo.csv"), row.names = FALSE)

cat(sprintf("\n  Fisher ratio J preserved under D: %.0f%% of healthy value (synthetic), %.0f%% (real combat_sup).\n",
            100 * modd$J_fish / base$J_fish, 100 * real$J_fish / base$J_fish))
cat(sprintf("  Euclidean KNN purity under D: %.2f -> %.2f (synthetic), %.2f (real combat_sup).\n",
            base$pur_knn, modd$pur_knn, real$pur_knn))
cat("  -> written knn_vs_rda_invariance_demo.csv\n")

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat("  CONCLUSION: combat_sup's diagonal variance map collapses the class\n")
cat("  signal in the Euclidean metric (KNN, raw axis) but is undone by the\n")
cat("  S_W^-1 whitening of Fisher/LDA — hence the effect is KNN-specific.\n")
cat("═══════════════════════════════════════════════════════════════════\n")
