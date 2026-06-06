#!/usr/bin/env Rscript
# decompose_steps.R
# Clean 2x2 attribution of the combat_sup KNN collapse to STEP 1 vs STEP 2.
#
# The two-step supervised-ComBat pipeline has two independent "supervision" knobs:
#   STEP 1 (harmonize training studies to each other):  mod = NULL  vs  mod = ~lab
#   STEP 2 (map raw test onto the corrected reference):  mod = NULL  vs  mod = ~clabs
# Section B of evaluate_knn_vs_rda.R confounded the two (both ref+test supervised);
# Section C only probed the step-2 per-gene rescaling. This script builds all four
# internally-consistent combinations from scratch and reports, on the TEST:
#   dp_raw  = d' on the raw mean-difference axis     (what Euclidean KNN sees)
#   dp_fish = d' on the whitened Fisher axis          (what LDA/RDA sees)
#   J_fish  = Mahalanobis dist^2 of class means       (Fisher separation)
#   KNN_pur = 5-NN purity of the test labelled by ref neighbours
#   LDA_acc = LDA accuracy on the test
# so we can see which knob moves which metric.

suppressMessages({ library(sva); library(matrixStats) })

DATA_FILE <- "data/TB_real_data.RData"
OUT_DIR   <- "outputs/diagnostics/hypothesis_tests"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
load(DATA_FILE)

ALL_STUDIES <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
bin <- function(l) as.integer(ifelse(l %in% c("1", 1, "Active"), 1L, 0L))
train_studies <- function(n, test) ALL_STUDIES[ALL_STUDIES != test][seq_len(n)]

log_safe <- function(m) {
  if (max(m, na.rm = TRUE) > 100) {
    mn <- min(m, na.rm = TRUE); if (mn < 0) m <- m - mn; m <- log2(m + 1)
  }
  m
}

# ── Build one (reference, test) pair for a chosen supervision of each step ─────
two_step <- function(dat, bat, lab, tst, ltst, sup1, sup2) {
  # STEP 1: harmonize training batches
  if (length(unique(bat)) >= 2) {
    mod1  <- if (sup1) model.matrix(~lab) else NULL
    dat_c <- suppressMessages(ComBat(dat, batch = bat, mod = mod1))
  } else dat_c <- dat
  # STEP 2: map raw test onto corrected reference (ref.batch = 1)
  comb  <- cbind(dat_c, tst)
  cb    <- c(rep(1L, ncol(dat_c)), rep(2L, ncol(tst)))
  mod2  <- if (sup2) model.matrix(~ c(lab, ltst)) else NULL
  out   <- suppressMessages(ComBat(comb, batch = cb, mod = mod2, ref.batch = 1L))
  list(ref = dat_c, tgt = out[, (ncol(dat_c) + 1):ncol(out), drop = FALSE])
}

# ── Metrics in a reduced reference-PC space ───────────────────────────────────
analyse <- function(ref, tgt, lab_r, lab_t, K = 50, ridge = 1e-2) {
  cg  <- intersect(rownames(ref), rownames(tgt)); ref <- ref[cg, ]; tgt <- tgt[cg, ]
  K   <- min(K, ncol(ref) - 2)
  pca <- prcomp(t(ref), center = TRUE, scale. = FALSE, rank. = K)
  Rp  <- pca$x; Tp <- predict(pca, t(tgt))
  m1  <- colMeans(Rp[lab_r == 1, , drop = FALSE])
  m0  <- colMeans(Rp[lab_r == 0, , drop = FALSE]); dmu <- m1 - m0
  c1  <- scale(Rp[lab_r == 1, , drop = FALSE], center = m1, scale = FALSE)
  c0  <- scale(Rp[lab_r == 0, , drop = FALSE], center = m0, scale = FALSE)
  Sw  <- (crossprod(c1) + crossprod(c0)) / (nrow(Rp) - 2)
  Sw  <- Sw + diag(ridge * mean(diag(Sw)), K)
  w_raw  <- dmu / sqrt(sum(dmu^2))
  w_fish <- solve(Sw, dmu); w_fish <- w_fish / sqrt(sum(w_fish^2))
  dprime <- function(w) {
    s <- as.numeric(Tp %*% w)
    sp <- sqrt(((sum(lab_t==1)-1)*var(s[lab_t==1]) + (sum(lab_t==0)-1)*var(s[lab_t==0])) /
               (length(s) - 2))
    (mean(s[lab_t==1]) - mean(s[lab_t==0])) / max(sp, 1e-8)
  }
  J_fish <- as.numeric(t(dmu) %*% solve(Sw, dmu))
  d <- as.matrix(dist(rbind(Tp, Rp)))[seq_len(nrow(Tp)), (nrow(Tp)+1):(nrow(Tp)+nrow(Rp))]
  pur <- mean(sapply(seq_len(nrow(Tp)), function(i) mean(lab_r[order(d[i,])[1:5]] == lab_t[i])))
  mid <- as.numeric((m1 + m0) / 2) %*% w_fish
  pred <- as.integer(as.numeric(Tp %*% w_fish) > as.numeric(mid))
  if (mean(pred == lab_t) < 0.5) pred <- 1L - pred
  list(dp_raw = dprime(w_raw), dp_fish = dprime(w_fish), J_fish = J_fish,
       pur_knn = pur, acc_lda = mean(pred == lab_t))
}

cat("══════════════════════════════════════════════════════════════════════\n")
cat("  2x2 decomposition: STEP-1 supervision x STEP-2 supervision\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

combos <- list(
  c(FALSE, FALSE),  # plain combat
  c(TRUE,  FALSE),  # supervised step 1 only
  c(FALSE, TRUE),   # supervised step 2 only
  c(TRUE,  TRUE))   # combat_sup
names(combos) <- c("unsup1_unsup2", "SUP1_unsup2", "unsup1_SUP2", "SUP1_SUP2")

rows <- list()
for (n in 2:5) for (test in ALL_STUDIES) {
  refs <- train_studies(n, test)
  cg   <- Reduce(intersect, lapply(c(refs, test), function(s) rownames(dat_lst[[s]])))
  dat  <- do.call(cbind, lapply(refs, function(s) dat_lst[[s]][cg, ]))
  bat  <- do.call(c, lapply(refs, function(s) rep(s, ncol(dat_lst[[s]]))))
  lab  <- do.call(c, lapply(refs, function(s) bin(label_lst[[s]])))
  tst  <- dat_lst[[test]][cg, ]; ltst <- bin(label_lst[[test]])
  dat  <- log_safe(dat); tst <- log_safe(tst)
  for (cn in names(combos)) {
    sp <- combos[[cn]]
    fit <- tryCatch(two_step(dat, bat, lab, tst, ltst, sp[1], sp[2]),
                    error = function(e) NULL)
    if (is.null(fit)) next
    a <- analyse(fit$ref, fit$tgt, lab, ltst)
    rows[[length(rows)+1]] <- data.frame(
      n = n, test = test, combo = cn,
      dp_raw = a$dp_raw, dp_fish = a$dp_fish, J_fish = a$J_fish,
      pur_knn = a$pur_knn, acc_lda = a$acc_lda)
  }
  cat(sprintf("  done n=%d test=%s\n", n, test))
}
D <- do.call(rbind, rows)
write.csv(D, file.path(OUT_DIR, "step_decomposition.csv"), row.names = FALSE)

cat("\n── Medians by combination (across all 24 scenarios) ──────────────────\n")
agg <- do.call(rbind, lapply(names(combos), function(cn) {
  s <- D[D$combo == cn, ]
  data.frame(combo = cn,
             dp_raw  = round(median(abs(s$dp_raw)),  3),
             dp_fish = round(median(abs(s$dp_fish)), 3),
             J_fish  = round(median(s$J_fish),       2),
             KNN_pur = round(median(s$pur_knn),      3),
             LDA_acc = round(median(s$acc_lda),      3))
}))
print(agg, row.names = FALSE)

cat("\n  Reading the table:\n")
cat("   * If KNN_pur / dp_raw collapse is driven by STEP 1, then SUP1_unsup2\n")
cat("     looks like SUP1_SUP2 (low), and unsup1_SUP2 looks like baseline (high).\n")
cat("   * If driven by STEP 2, the opposite pattern holds.\n")
cat("  -> written step_decomposition.csv\n")
