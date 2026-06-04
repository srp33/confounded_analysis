#!/usr/bin/env Rscript
# check_batch_class_split.R
# Replot combat_sup training space coloured by BATCH (not class) to verify
# whether the two sub-clouds correspond to the two training datasets.
# Also plots class within each batch to see if the split is batch-driven
# and whether the class direction is inverted between batches.

suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

DATA_FILE <- "data/TB_real_data.RData"
ADJ_DIR   <- "outputs/adjusted_data/all_scenarios"
OUT_DIR   <- "outputs/diagnostics/combat_sup_knn"

load(DATA_FILE)

ALL_STUDIES  <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
n_str        <- "3"
test_study   <- "USA"
ref_studies  <- ALL_STUDIES[seq_len(as.integer(n_str))]   # GSE37250_SA, USA, India
ref_studies  <- setdiff(ref_studies, test_study)           # GSE37250_SA, India

for (adjuster in c("combat", "combat_sup")) {
  stem    <- file.path(ADJ_DIR, sprintf("%s_n%s_test%s", adjuster, n_str, test_study))
  ref_mat <- as.matrix(read.csv(paste0(stem, "_reference.csv"), row.names = 1))
  tgt_mat <- as.matrix(read.csv(paste0(stem, "_target.csv"),    row.names = 1))

  common_genes <- intersect(rownames(ref_mat), rownames(tgt_mat))
  ref_t <- t(ref_mat[common_genes, , drop = FALSE])   # samples × genes
  tgt_t <- t(tgt_mat[common_genes, , drop = FALSE])

  # Derive batch and class labels by matching reference CSV sample names to dat_lst
  ref_samples  <- rownames(ref_t)   # sample names (ref_t is samples × genes)
  all_ref_studies <- ALL_STUDIES[seq_len(as.integer(n_str))]
  train_batch  <- sapply(ref_samples, function(s) {
    hit <- Filter(function(st) s %in% colnames(dat_lst[[st]]), all_ref_studies)
    if (length(hit)) hit[[1]] else NA_character_
  })
  label_lookup <- do.call(c, setNames(
    lapply(all_ref_studies, function(s) {
      l <- label_lst[[s]]
      setNames(as.integer(ifelse(l %in% c("1", 1, "Active"), 1L, 0L)),
               colnames(dat_lst[[s]]))
    }), NULL))
  train_labels <- label_lookup[ref_samples]
  test_labels  <- as.integer(ifelse(label_lst[[test_study]] %in% c("1", 1, "Active"), 1L, 0L))

  # PCA fit on training, project test
  pca    <- prcomp(ref_t, center = TRUE, scale. = FALSE, rank. = 2)
  ref_pc <- as.data.frame(pca$x)
  tgt_pc <- as.data.frame(predict(pca, tgt_t))
  pct    <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  ref_pc$batch <- train_batch
  ref_pc$class <- factor(train_labels, levels=c(0,1), labels=c("Control","Active TB"))
  ref_pc$split <- "Train"
  tgt_pc$batch <- test_study
  tgt_pc$class <- factor(test_labels,  levels=c(0,1), labels=c("Control","Active TB"))
  tgt_pc$split <- "Test"

  df <- rbind(ref_pc, tgt_pc)

  # ── plot 1: train coloured by BATCH, test grey ──────────────────────────────
  p1 <- ggplot(df, aes(PC1, PC2)) +
    geom_point(data = subset(df, split == "Train"),
               aes(colour = batch, shape = class), size = 2.5, alpha = 0.8) +
    geom_point(data = subset(df, split == "Test"),
               aes(shape = class), colour = "grey60", size = 2, alpha = 0.6) +
    labs(
      title    = sprintf("%s — training coloured by BATCH (n=%s, test=%s)", adjuster, n_str, test_study),
      subtitle = sprintf("PC1 %.1f%% | PC2 %.1f%%  — grey = test samples", pct[1], pct[2]),
      x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2])
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")

  # ── plot 2: class within each batch — does class direction match or flip? ───
  p2 <- ggplot(subset(df, split == "Train"), aes(PC1, PC2, colour = class)) +
    geom_point(size = 2, alpha = 0.7) +
    facet_wrap(~batch) +
    scale_colour_manual(values = c(Control = "#4e79a7", "Active TB" = "#e15759")) +
    labs(
      title    = sprintf("%s — class label within each training batch", adjuster),
      subtitle = "If class direction flips between batches, ComBat is confusing batch with biology",
      x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2])
    ) +
    theme_minimal(base_size = 12)

  png(file.path(OUT_DIR, sprintf("batch_split_%s.png", adjuster)),
      width = 9*150, height = 5*150, res = 150)
  print(p1)
  dev.off()

  png(file.path(OUT_DIR, sprintf("class_within_batch_%s.png", adjuster)),
      width = 10*150, height = 5*150, res = 150)
  print(p2)
  dev.off()

  # Print class balance per batch
  cat(sprintf("\n=== %s: class balance per training batch ===\n", adjuster))
  ref_pc %>%
    group_by(batch, class) %>%
    summarise(n = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = class, values_from = n) %>%
    mutate(pct_active = round(100 * `Active TB` / (Control + `Active TB`), 1)) %>%
    print()
}

cat(sprintf("\nOutputs in %s/\n", OUT_DIR))
