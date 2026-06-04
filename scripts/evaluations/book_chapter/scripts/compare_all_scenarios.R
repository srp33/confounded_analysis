#!/usr/bin/env Rscript
# compare_all_scenarios.R
# Cross-scenario comparison of combat, combat_sup, combat_sup_mo, combat_sup_nat
# across all n × test_study combinations present in the adjusted-data directory.

suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(matrixStats)
})

DATA_FILE   <- "data/TB_real_data.RData"
ADJ_DIR     <- "outputs/adjusted_data/all_scenarios"
OUT_DIR     <- "outputs/diagnostics/combat_sup_knn"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ALL_STUDIES <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
METHODS     <- c("combat", "combat_sup", "combat_sup_mo", "combat_sup_nat",
                 "combat_sup_mo_v2", "combat_sup_nat_v2")
N_PCS       <- 50
K_NEIGHBORS <- c(1, 5, 11)

METHOD_COLOURS <- c(
  combat            = "#4e79a7",
  combat_sup        = "#e15759",
  combat_sup_mo     = "#b07aa1",
  combat_sup_nat    = "#59a14f",
  combat_sup_mo_v2  = "#76b7b2",
  combat_sup_nat_v2 = "#edc948"
)

load(DATA_FILE)
bin <- function(l) as.integer(ifelse(l %in% c("1", 1, "Active"), 1L, 0L))

knn_purity <- function(train_mat, train_labels, test_mat, test_labels, k) {
  dists <- as.matrix(dist(rbind(test_mat, train_mat)))[
    seq_len(nrow(test_mat)),
    (nrow(test_mat) + 1):(nrow(test_mat) + nrow(train_mat))]
  sapply(seq_len(nrow(test_mat)), function(i) {
    nn <- order(dists[i, ])[seq_len(k)]
    mean(train_labels[nn] == test_labels[i])
  })
}

# Discover all n × test_study combos present for ALL methods
find_scenarios <- function() {
  rows <- list()
  for (n in 3:6) {
    for (test in ALL_STUDIES[seq_len(n)]) {
      files_exist <- sapply(METHODS, function(m) {
        file.exists(file.path(ADJ_DIR,
          sprintf("%s_n%d_test%s_target.csv", m, n, test)))
      })
      if (all(files_exist))
        rows[[length(rows)+1]] <- list(n=n, test=test)
    }
  }
  rows
}

scenarios <- find_scenarios()
cat(sprintf("Found %d complete scenarios\n", length(scenarios)))

purity_rows <- list()
axis_rows   <- list()

for (sc in scenarios) {
  n <- sc$n; test_study <- sc$test
  ref_studies  <- ALL_STUDIES[seq_len(n)]
  sc_label     <- sprintf("n%d_%s", n, test_study)
  train_labels <- do.call(c, lapply(ref_studies, function(s) bin(label_lst[[s]])))
  test_labels  <- bin(label_lst[[test_study]])

  cat(sprintf("  %s ...\n", sc_label))

  for (m in METHODS) {
    stem    <- file.path(ADJ_DIR, sprintf("%s_n%d_test%s", m, n, test_study))
    ref_mat <- tryCatch(
      as.matrix(read.csv(paste0(stem, "_reference.csv"), row.names = 1)),
      error = function(e) NULL)
    tgt_mat <- tryCatch(
      as.matrix(read.csv(paste0(stem, "_target.csv"), row.names = 1)),
      error = function(e) NULL)
    if (is.null(ref_mat) || is.null(tgt_mat)) next

    cg    <- intersect(rownames(ref_mat), rownames(tgt_mat))
    ref_t <- t(ref_mat[cg, , drop=FALSE])
    tgt_t <- t(tgt_mat[cg, , drop=FALSE])

    pca   <- prcomp(ref_t, center=TRUE, scale.=FALSE, rank.=N_PCS)
    ref_r <- pca$x
    tgt_r <- predict(pca, tgt_t)

    for (k in K_NEIGHBORS) {
      purity_rows[[length(purity_rows)+1]] <- data.frame(
        scenario = sc_label, n = n, test_study = test_study,
        method = m, k = k,
        purity = knn_purity(ref_r, train_labels, tgt_r, test_labels, k),
        stringsAsFactors = FALSE)
    }

    # Class-axis ordering (gene space, not PCA)
    w <- colMeans(t(ref_mat[cg,])[train_labels==1,,drop=FALSE]) -
         colMeans(t(ref_mat[cg,])[train_labels==0,,drop=FALSE])
    w <- w / sqrt(sum(w^2))
    axis_rows[[length(axis_rows)+1]] <- data.frame(
      scenario   = sc_label, n = n, test_study = test_study, method = m,
      split      = c(rep("Train", ncol(ref_mat)), rep("Test", ncol(tgt_mat))),
      class      = factor(c(train_labels, test_labels), levels=c(0,1),
                          labels=c("Control","Active TB")),
      score      = c(as.numeric(t(ref_mat[cg,]) %*% w),
                     as.numeric(t(tgt_mat[cg,]) %*% w)),
      stringsAsFactors = FALSE)
  }
}

purity_df <- do.call(rbind, purity_rows)
axis_df   <- do.call(rbind, axis_rows)
purity_df$method <- factor(purity_df$method, levels=METHODS)
axis_df$method   <- factor(axis_df$method,   levels=METHODS)

# ── summary tables ─────────────────────────────────────────────────────────────
purity_summary <- purity_df %>%
  group_by(scenario, method, k) %>%
  summarise(mean_purity = mean(purity), .groups="drop")

cat("\n=== KNN purity at k=5, all scenarios ===\n")
purity_summary %>%
  filter(k == 5) %>%
  pivot_wider(names_from=method, values_from=mean_purity) %>%
  mutate(
    nat_vs_sup    = combat_sup_nat    - combat_sup,
    mo_vs_sup     = combat_sup_mo     - combat_sup,
    nat_v2_vs_sup = combat_sup_nat_v2 - combat_sup,
    mo_v2_vs_sup  = combat_sup_mo_v2  - combat_sup,
    nat_v2_vs_combat = combat_sup_nat_v2 - combat,
    mo_v2_vs_combat  = combat_sup_mo_v2  - combat
  ) %>%
  arrange(n, test_study) %>%
  print(width=200)

cat("\n=== Class-axis ordering (Test: Active - Control) at k=5 ===\n")
axis_df %>%
  filter(split == "Test") %>%
  group_by(scenario, method) %>%
  summarise(
    sep = mean(score[class=="Active TB"]) - mean(score[class=="Control"]),
    .groups="drop") %>%
  pivot_wider(names_from=method, values_from=sep) %>%
  arrange(scenario) %>%
  print(width=200)

# ── plot 1: purity by scenario (k=5) ──────────────────────────────────────────
p1 <- purity_summary %>%
  filter(k == 5) %>%
  ggplot(aes(x=method, y=mean_purity, fill=method)) +
  geom_col() +
  geom_hline(yintercept=0.5, linetype="dashed", colour="grey40") +
  facet_wrap(~scenario, ncol=4) +
  scale_fill_manual(values=METHOD_COLOURS) +
  scale_x_discrete(labels=c(
    combat         = "combat",
    combat_sup     = "sup",
    combat_sup_mo  = "mo",
    combat_sup_nat = "nat")) +
  labs(title="KNN purity (k=5) across all scenarios",
       subtitle="nat = limma mean-only + mean-only step2; mo = ComBat mean-only + mean-only step2",
       x=NULL, y="Mean KNN purity") +
  theme_minimal(base_size=10) +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=30, hjust=1))

ggsave(file.path(OUT_DIR, "all_scenarios_purity_k5.png"), p1,
       width=14, height=10, dpi=150)

# ── plot 2: purity improvement nat vs sup ─────────────────────────────────────
delta_df <- purity_summary %>%
  filter(k == 5) %>%
  pivot_wider(names_from=method, values_from=mean_purity) %>%
  mutate(
    `nat (mean-only step2)`    = combat_sup_nat    - combat_sup,
    `mo  (mean-only step2)`    = combat_sup_mo     - combat_sup,
    `nat_v2 (full step2)`      = combat_sup_nat_v2 - combat_sup,
    `mo_v2  (full step2)`      = combat_sup_mo_v2  - combat_sup
  ) %>%
  select(scenario, n,
         `nat (mean-only step2)`, `mo  (mean-only step2)`,
         `nat_v2 (full step2)`,   `mo_v2  (full step2)`) %>%
  pivot_longer(cols=c(`nat (mean-only step2)`, `mo  (mean-only step2)`,
                      `nat_v2 (full step2)`,   `mo_v2  (full step2)`),
               names_to="comparison", values_to="delta_purity")

p2 <- ggplot(delta_df, aes(x=scenario, y=delta_purity, fill=comparison)) +
  geom_col(position="dodge") +
  geom_hline(yintercept=0, colour="grey30") +
  scale_fill_manual(values=c(
    "nat (mean-only step2)"  = "#59a14f",
    "mo  (mean-only step2)"  = "#b07aa1",
    "nat_v2 (full step2)"    = "#edc948",
    "mo_v2  (full step2)"    = "#76b7b2")) +
  labs(title="Purity improvement over combat_sup (k=5)",
       subtitle="v2 = mean-only step-1, full variance correction step-2",
       x=NULL, y="Δ mean KNN purity") +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=35, hjust=1), legend.position="bottom")

ggsave(file.path(OUT_DIR, "all_scenarios_delta_purity.png"), p2,
       width=14, height=5, dpi=150)

# ── plot 3: purity line plot over k, faceted by scenario ──────────────────────
p3 <- purity_summary %>%
  ggplot(aes(x=factor(k), y=mean_purity, colour=method, group=method)) +
  geom_line() + geom_point(size=2) +
  geom_hline(yintercept=0.5, linetype="dashed", colour="grey50") +
  facet_wrap(~scenario, ncol=4) +
  scale_colour_manual(values=METHOD_COLOURS) +
  labs(title="KNN purity across k and scenarios",
       x="k", y="Mean KNN purity") +
  theme_minimal(base_size=10) +
  theme(legend.position="bottom")

ggsave(file.path(OUT_DIR, "all_scenarios_purity_by_k.png"), p3,
       width=14, height=10, dpi=150)

cat(sprintf("\nPlots written to %s/\n", OUT_DIR))
