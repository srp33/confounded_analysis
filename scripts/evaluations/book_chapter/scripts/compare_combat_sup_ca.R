#!/usr/bin/env Rscript
# compare_combat_sup_ca.R
# Compares KNN neighbor purity for combat vs combat_sup vs combat_sup_ca
# on the n3/testUSA scenario, to test whether class-agnostic projection
# recovers the KNN performance that combat_sup destroyed.

suppressMessages({ library(ggplot2); library(dplyr) })

DATA_FILE <- "data/TB_real_data.RData"
ADJ_DIR   <- "outputs/adjusted_data/all_scenarios"
OUT_DIR   <- "outputs/diagnostics/combat_sup_knn"

ALL_STUDIES <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
n_str       <- "3"
test_study  <- "USA"
N_PCS       <- 50
K_NEIGHBORS <- c(1, 5, 11)
METHODS     <- c("combat", "combat_sup", "combat_sup_ca", "combat_sup_mg")

load(DATA_FILE)
ref_studies  <- ALL_STUDIES[seq_len(as.integer(n_str))]
train_labels <- do.call(c, lapply(label_lst[ref_studies], function(l)
  as.integer(ifelse(l %in% c("1", 1, "Active"), 1L, 0L))))
test_labels  <- as.integer(ifelse(label_lst[[test_study]] %in% c("1", 1, "Active"), 1L, 0L))

knn_purity <- function(train_mat, train_labels, test_mat, test_labels, k) {
  dists <- as.matrix(dist(rbind(test_mat, train_mat)))[
    seq_len(nrow(test_mat)),
    (nrow(test_mat) + 1):(nrow(test_mat) + nrow(train_mat))]
  sapply(seq_len(nrow(test_mat)), function(i) {
    nn <- order(dists[i, ])[seq_len(k)]
    mean(train_labels[nn] == test_labels[i])
  })
}

purity_rows <- list()
pca_rows    <- list()

for (m in METHODS) {
  stem    <- file.path(ADJ_DIR, sprintf("%s_n%s_test%s", m, n_str, test_study))
  ref_mat <- as.matrix(read.csv(paste0(stem, "_reference.csv"), row.names = 1))
  tgt_mat <- as.matrix(read.csv(paste0(stem, "_target.csv"),    row.names = 1))
  cg      <- intersect(rownames(ref_mat), rownames(tgt_mat))
  ref_t   <- t(ref_mat[cg, , drop = FALSE])
  tgt_t   <- t(tgt_mat[cg, , drop = FALSE])

  pca   <- prcomp(ref_t, center = TRUE, scale. = FALSE, rank. = N_PCS)
  ref_r <- pca$x
  tgt_r <- predict(pca, tgt_t)
  pct   <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  for (k in K_NEIGHBORS) {
    purity_rows[[length(purity_rows)+1]] <- data.frame(
      method = m, k = k,
      purity = knn_purity(ref_r, train_labels, tgt_r, test_labels, k))
  }

  pca_rows[[length(pca_rows)+1]] <- data.frame(
    method = m,
    PC1 = c(ref_r[,1], tgt_r[,1]), PC2 = c(ref_r[,2], tgt_r[,2]),
    label = factor(c(train_labels, test_labels), levels=c(0,1),
                   labels=c("Control","Active TB")),
    split = c(rep("Train", nrow(ref_r)), rep("Test", nrow(tgt_r))),
    pc1v = pct[1], pc2v = pct[2])
}

purity_df <- do.call(rbind, purity_rows)
purity_df$method <- factor(purity_df$method, levels = METHODS)

cat("\n=== KNN purity by method (mean over test samples) ===\n")
purity_df %>%
  group_by(method, k) %>%
  summarise(mean_purity = mean(purity), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = k, values_from = mean_purity,
                     names_prefix = "k=") %>%
  print()

# ── purity plot ──
p1 <- purity_df %>%
  group_by(method, k) %>%
  summarise(mean_purity = mean(purity),
            se = sd(purity)/sqrt(n()), .groups = "drop") %>%
  ggplot(aes(factor(k), mean_purity, colour = method, group = method)) +
  geom_line() + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_purity - se, ymax = mean_purity + se), width = 0.15) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c(combat = "#4e79a7", combat_sup = "#e15759",
                                 combat_sup_ca = "#59a14f")) +
  labs(title = sprintf("KNN purity: class-agnostic projection (%s, n=%s)", test_study, n_str),
       subtitle = "Dashed line = random chance (0.5). Higher is better.",
       x = "k", y = "Mean KNN purity") +
  theme_minimal(base_size = 12)

png(file.path(OUT_DIR, "combat_sup_ca_purity.png"), width = 8*150, height = 5*150, res = 150)
print(p1); dev.off()

# ── PCA grid ──
pca_df <- do.call(rbind, pca_rows)
pca_df$method <- factor(pca_df$method, levels = METHODS)
p2 <- ggplot(pca_df, aes(PC1, PC2, colour = label, shape = split, size = split)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~method, scales = "free") +
  scale_shape_manual(values = c(Train = 16, Test = 17)) +
  scale_size_manual(values = c(Train = 1.5, Test = 3)) +
  scale_colour_manual(values = c(Control = "#4e79a7", "Active TB" = "#e15759")) +
  labs(title = sprintf("Test projection in training PCA space (%s, n=%s)", test_study, n_str),
       subtitle = "Triangles = test. combat_sup_ca should keep test samples in correct class regions.",
       x = "PC1 (train)", y = "PC2 (train)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

png(file.path(OUT_DIR, "combat_sup_ca_pca.png"), width = 14*150, height = 5*150, res = 150)
print(p2); dev.off()

cat(sprintf("\nPlots: %s/combat_sup_ca_purity.png, combat_sup_ca_pca.png\n", OUT_DIR))
