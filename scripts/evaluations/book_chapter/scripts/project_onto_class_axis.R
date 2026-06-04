#!/usr/bin/env Rscript
# project_onto_class_axis.R
# Direction is correct (cos(class_adj, class_USA)=0.79) but KNN still fails.
# Test whether it's a POSITIONING problem: project corrected train + test onto
# the class axis (centroid_Active - centroid_Control of corrected TRAINING) and
# see where train-Active / train-Control / test-Active / test-Control land.
#
# If test clusters are shifted/scaled off the training clusters along this axis,
# that's a translation/scale mismatch, not a rotation — and pinpoints the fix.

suppressMessages({ library(ggplot2); library(dplyr) })

DATA_FILE <- "data/TB_real_data.RData"
ADJ_DIR   <- "outputs/adjusted_data/all_scenarios"
OUT_DIR   <- "outputs/diagnostics/combat_sup_knn"

ALL_STUDIES <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
n_str       <- "3"
test_study  <- "USA"
METHODS     <- c("combat", "combat_sup", "combat_sup_ca", "combat_sup_mg")

load(DATA_FILE)
bin <- function(l) as.integer(ifelse(l %in% c("1", 1, "Active"), 1L, 0L))

# training label order matches adjust script: cbind over first-n studies in order
ref_studies  <- ALL_STUDIES[seq_len(as.integer(n_str))]   # SA, USA, India
train_labels <- do.call(c, lapply(ref_studies, function(s) bin(label_lst[[s]])))
test_labels  <- bin(label_lst[[test_study]])

rows <- list()
for (m in METHODS) {
  stem    <- file.path(ADJ_DIR, sprintf("%s_n%s_test%s", m, n_str, test_study))
  ref_mat <- as.matrix(read.csv(paste0(stem, "_reference.csv"), row.names = 1))
  tgt_mat <- as.matrix(read.csv(paste0(stem, "_target.csv"),    row.names = 1))
  cg      <- intersect(rownames(ref_mat), rownames(tgt_mat))
  ref     <- ref_mat[cg, , drop = FALSE]   # genes × train
  tgt     <- tgt_mat[cg, , drop = FALSE]   # genes × test

  # class axis from corrected TRAINING (unit vector)
  w <- rowMeans(ref[, train_labels == 1, drop = FALSE]) -
       rowMeans(ref[, train_labels == 0, drop = FALSE])
  w <- w / sqrt(sum(w^2))

  # project (genes are rows): score = t(w) %*% data
  ref_proj <- as.numeric(crossprod(w, ref))
  tgt_proj <- as.numeric(crossprod(w, tgt))

  rows[[length(rows)+1]] <- data.frame(
    method = m,
    score  = c(ref_proj, tgt_proj),
    class  = factor(c(train_labels, test_labels), levels = c(0,1),
                    labels = c("Control","Active TB")),
    split  = c(rep("Train", length(ref_proj)), rep("Test", length(tgt_proj))))
}
df <- do.call(rbind, rows)
df$method <- factor(df$method, levels = METHODS)
df$group  <- paste(df$split, df$class)

# ── numeric summary: centroid of each group along the class axis ──────────────
cat("\n=== Mean projection onto class axis, by method × split × class ===\n")
summ <- df %>%
  group_by(method, split, class) %>%
  summarise(mean_score = mean(score), sd_score = sd(score), .groups = "drop")
summ %>%
  tidyr::pivot_wider(names_from = c(split, class),
                     values_from = c(mean_score, sd_score)) %>%
  print(width = 200)

# Separation check: is test Active > test Control (correct ordering)?
cat("\n=== Class ordering along axis (Active mean - Control mean) ===\n")
df %>%
  group_by(method, split) %>%
  summarise(active_minus_control =
              mean(score[class=="Active TB"]) - mean(score[class=="Control"]),
            .groups = "drop") %>%
  tidyr::pivot_wider(names_from = split, values_from = active_minus_control) %>%
  print()

# ── plot: distribution along class axis ───────────────────────────────────────
p <- ggplot(df, aes(x = score, fill = class, linetype = split)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~method, scales = "free", ncol = 1) +
  scale_fill_manual(values = c(Control = "#4e79a7", "Active TB" = "#e15759")) +
  labs(title    = sprintf("Projection onto class axis (%s, n=%s)", test_study, n_str),
       subtitle = "Solid = train, dashed = test. If test (dashed) clusters don't sit over the matching-colour train cluster → positioning mismatch.",
       x = "Score along training class axis", y = "Density") +
  theme_minimal(base_size = 12)

png(file.path(OUT_DIR, "class_axis_projection.png"), width = 9*150, height = 9*150, res = 150)
print(p); dev.off()
cat(sprintf("\nPlot: %s/class_axis_projection.png\n", OUT_DIR))
