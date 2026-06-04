#!/usr/bin/env Rscript
# compare_combat_sup_nat.R
# Tests whether using batch-only var.pooled for reconstruction (ComBat_nat)
# fixes the test-class inversion that defeats KNN under supervised ComBat.
# Compares: combat, combat_sup, combat_sup_mg, combat_sup_nat
# Diagnostics: KNN purity, class-axis ordering, PCA grid.

suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

DATA_FILE   <- "data/TB_real_data.RData"
ADJ_DIR     <- "outputs/adjusted_data/all_scenarios"
OUT_DIR     <- "outputs/diagnostics/combat_sup_knn"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ALL_STUDIES <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
n_str       <- "3"
test_study  <- "USA"
N_PCS       <- 50
K_NEIGHBORS <- c(1, 5, 11)
METHODS     <- c("combat", "combat_sup", "combat_sup_mg", "combat_sup_nat", "combat_sup_mo")

load(DATA_FILE)
bin <- function(l) as.integer(ifelse(l %in% c("1", 1, "Active"), 1L, 0L))
ref_studies  <- ALL_STUDIES[seq_len(as.integer(n_str))]
train_labels <- do.call(c, lapply(ref_studies, function(s) bin(label_lst[[s]])))
test_labels  <- bin(label_lst[[test_study]])

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
axis_rows   <- list()

for (m in METHODS) {
  stem    <- file.path(ADJ_DIR, sprintf("%s_n%s_test%s", m, n_str, test_study))
  ref_mat <- tryCatch(as.matrix(read.csv(paste0(stem, "_reference.csv"), row.names = 1)),
                      error = function(e) { message("Missing: ", stem); NULL })
  tgt_mat <- tryCatch(as.matrix(read.csv(paste0(stem, "_target.csv"), row.names = 1)),
                      error = function(e) NULL)
  if (is.null(ref_mat) || is.null(tgt_mat)) next

  cg    <- intersect(rownames(ref_mat), rownames(tgt_mat))
  ref_t <- t(ref_mat[cg, , drop = FALSE])   # samples × genes
  tgt_t <- t(tgt_mat[cg, , drop = FALSE])

  pca   <- prcomp(ref_t, center = TRUE, scale. = FALSE, rank. = N_PCS)
  ref_r <- pca$x
  tgt_r <- predict(pca, tgt_t)
  pct   <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  # KNN purity
  for (k in K_NEIGHBORS) {
    purity_rows[[length(purity_rows)+1]] <- data.frame(
      method = m, k = k,
      purity = knn_purity(ref_r, train_labels, tgt_r, test_labels, k))
  }

  # Class-axis ordering (using raw gene space, not PCA)
  ref_g <- t(ref_mat[cg, ])
  tgt_g <- t(tgt_mat[cg, ])
  w <- colMeans(ref_g[train_labels == 1, , drop=FALSE]) -
       colMeans(ref_g[train_labels == 0, , drop=FALSE])
  w <- w / sqrt(sum(w^2))
  ref_scores <- as.numeric(ref_g %*% w)
  tgt_scores <- as.numeric(tgt_g %*% w)
  axis_rows[[length(axis_rows)+1]] <- data.frame(
    method   = m,
    split    = c(rep("Train", length(ref_scores)), rep("Test", length(tgt_scores))),
    class    = factor(c(train_labels, test_labels), levels=c(0,1), labels=c("Control","Active TB")),
    score    = c(ref_scores, tgt_scores))

  # PCA data
  pca_rows[[length(pca_rows)+1]] <- data.frame(
    method = m,
    PC1 = c(ref_r[,1], tgt_r[,1]), PC2 = c(ref_r[,2], tgt_r[,2]),
    label = factor(c(train_labels, test_labels), levels=c(0,1),
                   labels=c("Control","Active TB")),
    split = c(rep("Train", nrow(ref_r)), rep("Test", nrow(tgt_r))),
    pc1v = pct[1], pc2v = pct[2])
}

purity_df <- do.call(rbind, purity_rows)
pca_df    <- do.call(rbind, pca_rows)
axis_df   <- do.call(rbind, axis_rows)

for (df in list(purity_df, pca_df, axis_df)) df$method <- factor(df$method, levels = METHODS)

# ── 1. KNN purity summary ──────────────────────────────────────────────────────
cat("\n=== KNN purity (mean over test samples) ===\n")
purity_df %>%
  group_by(method, k) %>%
  summarise(mean_purity = mean(purity), .groups = "drop") %>%
  pivot_wider(names_from = k, values_from = mean_purity, names_prefix = "k=") %>%
  print()

# ── 2. Class-axis ordering ─────────────────────────────────────────────────────
cat("\n=== Class ordering on training class axis (Active - Control mean score) ===\n")
axis_df %>%
  group_by(method, split) %>%
  summarise(sep = mean(score[class == "Active TB"]) - mean(score[class == "Control"]),
            .groups = "drop") %>%
  pivot_wider(names_from = split, values_from = sep) %>%
  print()

# ── 3. Purity plot ─────────────────────────────────────────────────────────────
METHOD_COLOURS <- c(combat         = "#4e79a7",
                    combat_sup     = "#e15759",
                    combat_sup_mg  = "#f28e2b",
                    combat_sup_nat = "#59a14f",
                    combat_sup_mo  = "#b07aa1")

p_purity <- purity_df %>%
  group_by(method, k) %>%
  summarise(mean_purity = mean(purity), se = sd(purity)/sqrt(n()), .groups = "drop") %>%
  ggplot(aes(factor(k), mean_purity, colour = method, group = method)) +
  geom_line() + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_purity - se, ymax = mean_purity + se), width = 0.15) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = METHOD_COLOURS) +
  labs(title    = sprintf("KNN purity: natural-variance supervised ComBat (%s, n=%s)", test_study, n_str),
       subtitle = "Dashed = random chance. combat_sup_nat uses batch-only var.pooled for reconstruction.",
       x = "k", y = "Mean KNN purity") +
  theme_minimal(base_size = 12)

png(file.path(OUT_DIR, "combat_sup_nat_purity.png"), width = 8*150, height = 5*150, res = 150)
print(p_purity); dev.off()

# ── 4. Class axis density plot ─────────────────────────────────────────────────
p_axis <- ggplot(axis_df, aes(x = score, fill = class, linetype = split)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~method, scales = "free", ncol = 2) +
  scale_fill_manual(values = c(Control = "#4e79a7", "Active TB" = "#e15759")) +
  scale_linetype_manual(values = c(Train = "solid", Test = "dashed")) +
  labs(title    = sprintf("Class-axis projection (%s, n=%s)", test_study, n_str),
       subtitle = "Dashed = test. Correct = test Active (red dashed) right of test Control (blue dashed).",
       x = "Score on training class axis", y = "Density") +
  theme_minimal(base_size = 11)

png(file.path(OUT_DIR, "combat_sup_nat_class_axis.png"), width = 12*150, height = 8*150, res = 150)
print(p_axis); dev.off()

# ── 5. PCA grid ────────────────────────────────────────────────────────────────
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = label, shape = split, size = split)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~method, scales = "free", ncol = 2,
             labeller = labeller(method = function(x) {
               sprintf("%s (PC1=%.1f%%)", x, pca_df$pc1v[match(x, pca_df$method)])
             })) +
  scale_shape_manual(values = c(Train = 16, Test = 17)) +
  scale_size_manual(values = c(Train = 1.5, Test = 3)) +
  scale_colour_manual(values = c(Control = "#4e79a7", "Active TB" = "#e15759")) +
  labs(title    = sprintf("PCA: test projection into training space (%s, n=%s)", test_study, n_str),
       subtitle = "Triangles = test samples. Correct = red triangles in red training cloud.",
       x = "PC1 (train)", y = "PC2 (train)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

png(file.path(OUT_DIR, "combat_sup_nat_pca.png"), width = 14*150, height = 10*150, res = 150)
print(p_pca); dev.off()

cat(sprintf("\nPlots written to %s/\n", OUT_DIR))
