#!/usr/bin/env Rscript
# check_rank_preservation.R
# For each gene, compute Spearman correlation between unadjusted and adjusted
# test samples. A healthy affine adjustment should give cor ~ +1 for all genes.
# A negative scale in ComBat's correction will show up as cor ~ -1.

suppressMessages({ library(ggplot2); library(dplyr) })

DATA_FILE <- "data/TB_real_data.RData"
ADJ_DIR   <- "outputs/adjusted_data/all_scenarios"
OUT_DIR   <- "outputs/diagnostics/combat_sup_knn"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SCENARIO  <- list(n = "3", test = "USA")   # single scenario for speed

load(DATA_FILE)

n_str      <- SCENARIO$n
test_study <- SCENARIO$test

# Raw (unadjusted) test data
raw_test <- dat_lst[[test_study]]   # genes × samples

results <- list()

for (adjuster in c("unadjusted", "combat", "combat_sup")) {

  if (adjuster == "unadjusted") {
    adj_test <- raw_test
  } else {
    f <- file.path(ADJ_DIR,
                   sprintf("%s_n%s_test%s_target.csv", adjuster, n_str, test_study))
    adj_test <- as.matrix(read.csv(f, row.names = 1))
  }

  common_genes <- intersect(rownames(raw_test), rownames(adj_test))
  common_samps <- intersect(colnames(raw_test), colnames(adj_test))

  cat(sprintf("%s: %d genes, %d common test samples\n",
              adjuster, length(common_genes), length(common_samps)))

  raw_sub <- raw_test[common_genes, common_samps, drop = FALSE]
  adj_sub <- adj_test[common_genes, common_samps, drop = FALSE]

  # Per-gene Spearman correlation between unadjusted and adjusted sample ordering
  cors <- vapply(seq_len(nrow(raw_sub)), function(i)
    cor(raw_sub[i, ], adj_sub[i, ], method = "spearman"),
    numeric(1))

  results[[adjuster]] <- data.frame(
    adjuster = adjuster,
    gene     = common_genes,
    spearman = cors
  )
}

df <- do.call(rbind, results)
df$adjuster <- factor(df$adjuster, levels = c("unadjusted", "combat", "combat_sup"))

# ── summary ──────────────────────────────────────────────────────────────────
cat("\n=== Spearman correlation summary (unadjusted vs adjusted, per gene) ===\n")
df %>%
  group_by(adjuster) %>%
  summarise(
    mean   = mean(spearman,          na.rm = TRUE),
    pct_negative = mean(spearman < 0, na.rm = TRUE),
    pct_below_neg_half = mean(spearman < -0.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

write.csv(df, file.path(OUT_DIR, "rank_preservation.csv"), row.names = FALSE)

# ── plot: density of per-gene Spearman correlations ──────────────────────────
p_rank <- ggplot(df, aes(x = spearman, fill = adjuster, colour = adjuster)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80,
                 alpha = 0.45, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values   = c(unadjusted = "#999999", combat = "#4e79a7", combat_sup = "#e15759")) +
  scale_colour_manual(values = c(unadjusted = "#999999", combat = "#4e79a7", combat_sup = "#e15759")) +
  labs(
    title    = sprintf("Per-gene rank preservation in test set (%s, n=%s)", test_study, n_str),
    subtitle = "Spearman correlation between unadjusted and adjusted sample ordering per gene\nHealthy = spike near +1; rank reversal = spike near -1",
    x = "Spearman correlation (unadjusted vs adjusted)", y = "Density"
  ) +
  theme_minimal(base_size = 12)

# ── distribution shift: where do test samples land relative to training? ─────
# Load adjusted data for train (reference) and test (target) for each adjuster
load(DATA_FILE)
ALL_STUDIES <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
ref_studies  <- ALL_STUDIES[seq_len(as.integer(n_str))]
train_labels <- do.call(c, lapply(label_lst[ref_studies], function(l)
  as.integer(ifelse(l %in% c("1", 1, "Active"), 1L, 0L))))
test_labels  <- as.integer(ifelse(label_lst[[test_study]] %in% c("1", 1, "Active"), 1L, 0L))

dist_rows <- list()
for (adjuster in c("combat", "combat_sup")) {
  stem    <- file.path(ADJ_DIR, sprintf("%s_n%s_test%s", adjuster, n_str, test_study))
  ref_mat <- as.matrix(read.csv(paste0(stem, "_reference.csv"), row.names = 1))
  tgt_mat <- as.matrix(read.csv(paste0(stem, "_target.csv"),    row.names = 1))

  # PCA on reference (train), project test into same space
  ref_t  <- t(ref_mat)
  tgt_t  <- t(tgt_mat[intersect(rownames(ref_mat), rownames(tgt_mat)), , drop = FALSE])
  pca    <- prcomp(ref_t, center = TRUE, scale. = FALSE, rank. = 2)
  ref_pc <- pca$x
  tgt_pc <- predict(pca, tgt_t)

  dist_rows[[paste(adjuster, "train")]] <- data.frame(
    adjuster = adjuster, split = "Train",
    PC1 = ref_pc[, 1], PC2 = ref_pc[, 2],
    label = factor(train_labels, levels = c(0,1), labels = c("Control","Active TB")))
  dist_rows[[paste(adjuster, "test")]] <- data.frame(
    adjuster = adjuster, split = "Test",
    PC1 = tgt_pc[, 1], PC2 = tgt_pc[, 2],
    label = factor(test_labels, levels = c(0,1), labels = c("Control","Active TB")))
}
dist_df <- do.call(rbind, dist_rows)

p_pca <- ggplot(dist_df, aes(PC1, PC2, colour = label, shape = split, size = split)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~adjuster) +
  scale_shape_manual(values = c(Train = 16, Test = 17)) +
  scale_size_manual(values  = c(Train = 1.5, Test = 3)) +
  scale_colour_manual(values = c(Control = "#4e79a7", "Active TB" = "#e15759")) +
  labs(
    title    = sprintf("PCA: where do test samples land in training space? (%s, n=%s)",
                       test_study, n_str),
    subtitle = "Triangles = test samples. If test triangles are inside wrong-class cluster → distribution shift.",
    x = "PC1 (train)", y = "PC2 (train)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

# save both plots
png(file.path(OUT_DIR, "rank_preservation.png"), width = 8*150, height = 4.5*150, res = 150)
print(p_rank)
dev.off()

png(file.path(OUT_DIR, "test_vs_train_pca.png"), width = 10*150, height = 5*150, res = 150)
print(p_pca)
dev.off()

cat(sprintf("\nPlots saved to %s/\n", OUT_DIR))
