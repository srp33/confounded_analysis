#!/usr/bin/env Rscript
# make_mechanism_figure.R  (distance-first rebuild)
# Chapter figure for the combat_sup x KNN mechanism, stated in DISTANCES not axes.
#  A. Direction-free failure: nearest-centroid balanced accuracy (Euclidean), per
#     scenario, unsupervised vs supervised step 1. Supervised drops below chance.
#  B. Why: that distance-based accuracy degrades with batch-class confounding.
#  C. What changed: the per-gene class-difference DIRECTION. Unsupervised stays
#     aligned with the raw/test class signal; supervised anti-aligns (confounded
#     partial coefficient). Cosines across scenarios.

suppressMessages({ library(ggplot2); library(dplyr); library(tidyr); library(gridExtra) })
OUT_DIR <- "outputs/diagnostics/hypothesis_tests"

aud <- read.csv(file.path(OUT_DIR, "mechanism_audit.csv"))
why <- read.csv(file.path(OUT_DIR, "mechanism_why.csv"))
conf <- why %>% filter(sup == TRUE) %>% select(n, test, confound)
aud  <- left_join(aud, conf, by = c("n", "test"))

col_u <- "#4e79a7"; col_s <- "#e15759"

# ── Panel A: nearest-centroid balanced accuracy (pure distance) ───────────────
A <- aud %>% select(n, test, Unsupervised = cen_bacc_unsup, Supervised = cen_bacc_sup) %>%
  pivot_longer(c(Unsupervised, Supervised), names_to = "step1", values_to = "bacc")
A$step1 <- factor(A$step1, levels = c("Unsupervised", "Supervised"))
pA <- ggplot(A, aes(step1, bacc, fill = step1)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.5, fill = "grey92") +
  annotate("text", x = "Unsupervised", y = 0.34, label = "below chance", colour = "grey45",
           size = 3, fontface = "italic", hjust = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey45") +
  geom_line(aes(group = interaction(n, test)), colour = "grey75", linewidth = 0.3) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.45) +
  geom_jitter(width = 0.08, size = 1.4, alpha = 0.7) +
  scale_fill_manual(values = c(Unsupervised = col_u, Supervised = col_s), guide = "none") +
  labs(title = "A  KNN failure is purely about distance",
       subtitle = "Nearest-centroid balanced acc (raw Euclidean, no fitted axis), per scenario",
       x = "Step-1 ComBat", y = "Nearest-centroid balanced accuracy") +
  theme_minimal(base_size = 11)

# ── Panel B: distance accuracy vs batch-class confounding (supervised) ────────
rhoB <- cor(aud$confound, aud$cen_bacc_sup, method = "spearman", use = "complete.obs")
pB <- ggplot(aud, aes(confound, cen_bacc_sup)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey45") +
  geom_smooth(method = "lm", se = TRUE, colour = col_s, fill = "#e1575933", linewidth = 0.8) +
  geom_point(size = 2.3, alpha = 0.8, colour = col_s) +
  annotate("text", x = min(aud$confound, na.rm = TRUE), y = 0.2, hjust = 0, size = 3.6,
           label = sprintf("Spearman rho = %.2f", rhoB)) +
  labs(title = "B  Confounding drives the distance failure",
       subtitle = "Supervised step 1: more batch-class confounding -> worse",
       x = "Batch-class confounding (Cramer's V)",
       y = "Nearest-centroid balanced acc (supervised)") +
  theme_minimal(base_size = 11)

# ── Panel C: per-gene class-difference DIRECTION (cosines) ────────────────────
Cd <- aud %>% transmute(
  `raw vs\nunsup`   = cos_raw_unsup,
  `raw vs\nsup`     = cos_raw_sup,
  `unsup vs\ntest`  = cos_unsup_test,
  `sup vs\ntest`    = cos_sup_test) %>%
  pivot_longer(everything(), names_to = "pair", values_to = "cos")
Cd$pair <- factor(Cd$pair, levels = c("raw vs\nunsup","raw vs\nsup","unsup vs\ntest","sup vs\ntest"))
Cd$grp  <- ifelse(grepl("unsup", Cd$pair), "unsupervised", "supervised")
pC <- ggplot(Cd, aes(pair, cos, fill = grp)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.1, size = 1.1, alpha = 0.5) +
  scale_fill_manual(values = c(unsupervised = col_u, supervised = col_s), guide = "none") +
  labs(title = "C  The supervised class direction is the wrong one",
       subtitle = "Cosine of per-gene class-difference vectors (1 = same direction, <0 = opposed)",
       x = NULL, y = "Cosine similarity") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(size = 8.5))

top <- arrangeGrob(pA, pB, ncol = 2)
fig <- arrangeGrob(top, pC, ncol = 1, heights = c(1, 0.9))
ggsave(file.path(OUT_DIR, "mechanism_figure.png"), fig, width = 12, height = 10, dpi = 200)
ggsave(file.path(OUT_DIR, "mechanism_figure.pdf"), fig, width = 12, height = 10)
cat(sprintf("Wrote mechanism_figure.png/.pdf | Panel B rho=%.3f\n", rhoB))
cat(sprintf("medians: cen_unsup=%.3f cen_sup=%.3f | cos raw_unsup=%.2f raw_sup=%.2f unsup_test=%.2f sup_test=%.2f\n",
            median(aud$cen_bacc_unsup), median(aud$cen_bacc_sup),
            median(aud$cos_raw_unsup), median(aud$cos_raw_sup),
            median(aud$cos_unsup_test), median(aud$cos_sup_test)))
