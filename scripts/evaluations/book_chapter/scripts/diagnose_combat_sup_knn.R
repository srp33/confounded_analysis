#!/usr/bin/env Rscript
# diagnose_combat_sup_knn.R
# Investigates why ComBat supervised has a strong negative interaction with KNN.
# Three angles:
#   1. PCA — do class clusters look different under combat vs combat_sup?
#   2. KNN neighbor purity — for each test sample, what fraction of its k nearest
#      training neighbors are the correct class?
#   3. Within- vs between-class distance ratio (class separability in feature space)

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})

DATA_FILE    <- "data/TB_real_data.RData"
ADJ_DIR      <- "outputs/adjusted_data/all_scenarios"
OUT_DIR      <- "outputs/diagnostics/combat_sup_knn"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ALL_STUDIES  <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
SCENARIOS    <- list(
  list(n = "3", test = "USA"),
  list(n = "3", test = "India"),
  list(n = "4", test = "GSE37250_M"),
  list(n = "4", test = "USA")
)
K_NEIGHBORS  <- c(1, 5, 11)
N_PCS        <- 50   # reduce to this many PCs before any distance computation

# ── helpers ──────────────────────────────────────────────────────────────────

load_scenario <- function(adjuster, n, test_study) {
  stem <- file.path(ADJ_DIR, sprintf("%s_n%s_test%s", adjuster, n, test_study))
  ref  <- as.matrix(read.csv(paste0(stem, "_reference.csv"), row.names = 1))
  tgt  <- as.matrix(read.csv(paste0(stem, "_target.csv"),    row.names = 1))
  list(ref = ref, tgt = tgt)
}

# Project train+test into shared PCA space fit on train only, return list of
# reduced matrices (samples × N_PCS) and the PCA object for variance reporting.
pca_project <- function(ref_t, tgt_t, n_pcs = N_PCS) {
  pca   <- prcomp(ref_t, center = TRUE, scale. = FALSE, rank. = n_pcs)
  ref_r <- pca$x[, seq_len(n_pcs), drop = FALSE]
  tgt_r <- predict(pca, tgt_t)[, seq_len(n_pcs), drop = FALSE]
  list(ref = ref_r, tgt = tgt_r, pca = pca)
}

knn_purity <- function(train_mat, train_labels, test_mat, test_labels, k) {
  # All in reduced PCA space (samples × PCs) — fast
  dists <- as.matrix(dist(rbind(test_mat, train_mat)))[
    seq_len(nrow(test_mat)),
    (nrow(test_mat) + 1):(nrow(test_mat) + nrow(train_mat))
  ]
  sapply(seq_len(nrow(test_mat)), function(i) {
    nn_idx    <- order(dists[i, ])[seq_len(k)]
    nn_labels <- train_labels[nn_idx]
    mean(nn_labels == test_labels[i])
  })
}

class_separability <- function(ref_r, tgt_r, train_labels, test_labels) {
  # ratio: mean between-class / mean within-class distance in PCA space
  mat    <- rbind(ref_r, tgt_r)
  labels <- c(train_labels, test_labels)
  d      <- as.matrix(dist(mat))
  within  <- d[outer(labels, labels, "==") & lower.tri(d)]
  between <- d[outer(labels, labels, "!=") & lower.tri(d)]
  list(
    within_mean  = mean(within),
    between_mean = mean(between),
    ratio        = mean(between) / mean(within)
  )
}

# ── main loop ────────────────────────────────────────────────────────────────

load(DATA_FILE)   # dat_lst, label_lst

purity_rows <- list()
sep_rows    <- list()
pca_plots   <- list()

for (sc in SCENARIOS) {
  n_str <- sc$n; test_study <- sc$test
  ref_studies <- ALL_STUDIES[seq_len(as.integer(n_str))]
  label_str   <- sprintf("n%s_%s", n_str, test_study)
  cat(sprintf("\n[%s] scenario %s\n", format(Sys.time(), "%H:%M:%S"), label_str))

  train_labels <- do.call(c, label_lst[ref_studies])
  test_labels  <- label_lst[[test_study]]
  train_labels <- as.integer(ifelse(train_labels %in% c("1", 1, "Active"), 1L, 0L))
  test_labels  <- as.integer(ifelse(test_labels  %in% c("1", 1, "Active"), 1L, 0L))

  for (adjuster in c("combat", "combat_sup")) {
    cat(sprintf("  %s ... loading\n", adjuster))
    d <- tryCatch(load_scenario(adjuster, n_str, test_study), error = function(e) NULL)
    if (is.null(d)) {
      message(sprintf("  Missing file for %s %s — skipping", adjuster, label_str)); next
    }

    common_genes <- intersect(rownames(d$ref), rownames(d$tgt))
    ref_t <- t(d$ref[common_genes, , drop = FALSE])   # samples × genes
    tgt_t <- t(d$tgt[common_genes, , drop = FALSE])

    # Project into shared PCA space (fit on train) — all distance work done here
    cat(sprintf("  %s ... PCA projection (%d genes → %d PCs)\n",
                adjuster, length(common_genes), N_PCS))
    proj <- pca_project(ref_t, tgt_t)
    pct  <- round(100 * proj$pca$sdev^2 / sum(proj$pca$sdev^2), 1)

    # ── 1. KNN purity (in PCA space) ─────────────────────────────────────────
    cat(sprintf("  %s ... KNN purity\n", adjuster))
    for (k in K_NEIGHBORS) {
      purity <- knn_purity(proj$ref, train_labels, proj$tgt, test_labels, k)
      purity_rows[[length(purity_rows) + 1]] <- data.frame(
        scenario    = label_str,
        adjuster    = adjuster,
        k           = k,
        test_sample = colnames(d$tgt),
        true_label  = test_labels,
        purity      = purity,
        stringsAsFactors = FALSE
      )
    }

    # ── 2. Class separability (in PCA space) ──────────────────────────────────
    cat(sprintf("  %s ... class separability\n", adjuster))
    sep <- class_separability(proj$ref, proj$tgt, train_labels, test_labels)
    sep_rows[[length(sep_rows) + 1]] <- data.frame(
      scenario     = label_str,
      adjuster     = adjuster,
      within_mean  = sep$within_mean,
      between_mean = sep$between_mean,
      ratio        = sep$ratio,
      stringsAsFactors = FALSE
    )

    # ── 3. PCA plot (PC1/PC2 of the same projection) ──────────────────────────
    all_labels <- c(train_labels, test_labels)
    pca_df <- data.frame(
      PC1   = c(proj$ref[, 1], proj$tgt[, 1]),
      PC2   = c(proj$ref[, 2], proj$tgt[, 2]),
      label = factor(all_labels, levels = c(0, 1), labels = c("Control", "Active TB")),
      split = c(rep("Train", nrow(proj$ref)), rep("Test", nrow(proj$tgt))),
      stringsAsFactors = FALSE
    )
    p <- ggplot(pca_df, aes(PC1, PC2, colour = label, shape = split)) +
      geom_point(alpha = 0.7, size = 2) +
      scale_shape_manual(values = c(Train = 16, Test = 17)) +
      labs(
        title    = sprintf("%s — %s", adjuster, label_str),
        subtitle = sprintf("PC1 %.1f%% | PC2 %.1f%%", pct[1], pct[2]),
        x = sprintf("PC1 (%.1f%%)", pct[1]),
        y = sprintf("PC2 (%.1f%%)", pct[2])
      ) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom")
    pca_plots[[paste(adjuster, label_str)]] <- p
  }
}

# ── collate and save ─────────────────────────────────────────────────────────

purity_df <- do.call(rbind, purity_rows)
sep_df    <- do.call(rbind, sep_rows)

write.csv(purity_df, file.path(OUT_DIR, "knn_purity.csv"),      row.names = FALSE)
write.csv(sep_df,    file.path(OUT_DIR, "class_separability.csv"), row.names = FALSE)

# ── plot 1: KNN purity distributions, combat vs combat_sup ───────────────────
purity_summary <- purity_df %>%
  group_by(scenario, adjuster, k) %>%
  summarise(
    mean_purity = mean(purity),
    se          = sd(purity) / sqrt(n()),
    .groups     = "drop"
  )

p_purity <- ggplot(purity_summary, aes(x = factor(k), y = mean_purity,
                                        colour = adjuster, group = adjuster)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_purity - se, ymax = mean_purity + se), width = 0.15) +
  facet_wrap(~scenario, nrow = 2) +
  scale_colour_manual(values = c(combat = "#4e79a7", combat_sup = "#e15759")) +
  labs(
    title = "KNN neighbor purity: ComBat vs ComBat Supervised",
    subtitle = "Fraction of k nearest training neighbors with the correct class label",
    x = "k (number of neighbors)", y = "Mean purity (correct class fraction)"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(OUT_DIR, "knn_purity.png"), p_purity,
       width = 10, height = 6, dpi = 150)

# ── plot 2: between/within distance ratio ────────────────────────────────────
p_sep <- ggplot(sep_df, aes(x = scenario, y = ratio, fill = adjuster)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(combat = "#4e79a7", combat_sup = "#e15759")) +
  labs(
    title    = "Class separability: between / within distance ratio",
    subtitle = "Higher = more separable. Covers train + test combined.",
    x = NULL, y = "Between / within distance ratio"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(OUT_DIR, "class_separability.png"), p_sep,
       width = 8, height = 4, dpi = 150)

# ── plot 3: per-sample purity violin, k=5, combat vs combat_sup ──────────────
purity_k5 <- purity_df %>% filter(k == 5)

p_violin <- ggplot(purity_k5, aes(x = adjuster, y = purity, fill = adjuster)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.5) +
  facet_grid(true_label ~ scenario,
             labeller = labeller(true_label = c(`0` = "Control", `1` = "Active TB"))) +
  scale_fill_manual(values = c(combat = "#4e79a7", combat_sup = "#e15759")) +
  labs(
    title    = "Per-sample KNN purity at k=5, by true class",
    subtitle = "Each dot = one test sample. Purity = fraction of 5 nearest neighbors with correct class.",
    x = NULL, y = "KNN purity (k=5)"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "knn_purity_by_class.png"), p_violin,
       width = 12, height = 6, dpi = 150)

# ── plot 4: PCA grid ──────────────────────────────────────────────────────────
if (length(pca_plots) > 0) {
  ncols    <- 4
  nrows    <- ceiling(length(pca_plots) / ncols)
  png(file.path(OUT_DIR, "pca_grid.png"), width = 18 * 150, height = nrows * 4.5 * 150, res = 150)
  grid.arrange(grobs = pca_plots, ncol = ncols)
  dev.off()
}

# ── print summary ─────────────────────────────────────────────────────────────
cat("\n=== Class separability summary ===\n")
print(sep_df %>% select(scenario, adjuster, ratio) %>%
  pivot_wider(names_from = adjuster, values_from = ratio) %>%
  mutate(combat_sup_vs_combat = combat_sup / combat))

cat("\n=== KNN purity summary (k=5) ===\n")
print(purity_summary %>% filter(k == 5) %>%
  select(scenario, adjuster, mean_purity) %>%
  pivot_wider(names_from = adjuster, values_from = mean_purity) %>%
  mutate(delta = combat_sup - combat))

cat(sprintf("\nOutputs written to: %s\n", OUT_DIR))
