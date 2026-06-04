#!/usr/bin/env Rscript
# compare_covariance.R
# Compares covariance structure of test data: unadjusted vs combat vs combat_sup.
# Uses SVD of the centered data matrix (equivalent to eigendecomposition of the
# covariance matrix, but numerically faster and avoids forming the full n_genes^2 matrix).
#
# Outputs:
#   1. eigenspectrum.png   — sorted eigenvalues (variance per PC)
#   2. gene_variance.png   — per-gene variance ratio (adjusted / unadjusted)
#   3. pc_alignment.png    — cosine similarity between top PCs of each pair

suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

DATA_FILE <- "data/TB_real_data.RData"
ADJ_DIR   <- "outputs/adjusted_data/all_scenarios"
OUT_DIR   <- "outputs/diagnostics/combat_sup_knn"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_SCENARIO <- list(n = "3", test = "USA")
N_TOP_PCS  <- 50    # how many PCs to compare for alignment

load(DATA_FILE)
n_str      <- N_SCENARIO$n
test_study <- N_SCENARIO$test

# ── load matrices (genes × samples, test set only) ───────────────────────────

raw_test <- dat_lst[[test_study]]

mats <- list(unadjusted = raw_test)
for (adj in c("combat", "combat_sup")) {
  f <- file.path(ADJ_DIR, sprintf("%s_n%s_test%s_target.csv", adj, n_str, test_study))
  mats[[adj]] <- as.matrix(read.csv(f, row.names = 1))
}

# align to common genes (rows) across all three
common_genes <- Reduce(intersect, lapply(mats, rownames))
mats <- lapply(mats, function(m) m[common_genes, , drop = FALSE])
cat(sprintf("Comparing covariance for %d genes × %d test samples\n",
            length(common_genes), ncol(mats[[1]])))

# ── SVD-based eigenspectrum ───────────────────────────────────────────────────
# Center by gene (row), then SVD. Eigenvalues of covariance = singular values² / (n-1)

loadings_list <- list()
eig_df_list   <- list()

for (nm in names(mats)) {
  cat(sprintf("  SVD: %s ... ", nm)); flush.console()
  X   <- mats[[nm]]
  X_c <- X - rowMeans(X)
  sv  <- svd(X_c, nu = 0, nv = N_TOP_PCS)
  cat("done\n")
  n_keep <- min(N_TOP_PCS, length(sv$d))
  loadings_list[[nm]] <- sv$v[, seq_len(n_keep), drop = FALSE]
  eig_df_list[[nm]]   <- data.frame(
    adjuster   = nm,
    pc         = seq_len(n_keep),
    eigenvalue = sv$d[seq_len(n_keep)]^2 / (ncol(X) - 1)
  )
}

eig_df <- do.call(rbind, eig_df_list)

# ── 1. Eigenspectrum plot ─────────────────────────────────────────────────────

p_eig <- ggplot(eig_df, aes(x = pc, y = eigenvalue, colour = adjuster)) +
  geom_line() +
  geom_point(size = 1) +
  scale_colour_manual(values = c(unadjusted = "#999999",
                                  combat     = "#4e79a7",
                                  combat_sup = "#e15759")) +
  scale_y_log10() +
  labs(
    title    = sprintf("Eigenspectrum of test covariance: %s (n=%s)", test_study, n_str),
    subtitle = sprintf("Top %d eigenvalues via SVD of centered data (log scale)", N_TOP_PCS),
    x = "Principal component", y = "Eigenvalue (log scale)"
  ) +
  theme_minimal(base_size = 12)

png(file.path(OUT_DIR, "eigenspectrum.png"), width = 8*150, height = 5*150, res = 150)
print(p_eig)
dev.off()
cat("Saved eigenspectrum.png\n")

# ── 2. Per-gene variance ratio ────────────────────────────────────────────────

var_unadj <- apply(mats[["unadjusted"]], 1, var)

gene_var_df <- lapply(c("combat", "combat_sup"), function(adj) {
  var_adj <- apply(mats[[adj]], 1, var)
  data.frame(
    adjuster = adj,
    gene     = common_genes,
    var_unadj = var_unadj,
    var_adj   = var_adj,
    ratio     = var_adj / pmax(var_unadj, 1e-10)
  )
})
gene_var_df <- do.call(rbind, gene_var_df)

cat("\n=== Per-gene variance ratio summary ===\n")
gene_var_df %>%
  group_by(adjuster) %>%
  summarise(
    mean_ratio    = mean(ratio,   na.rm = TRUE),
    median_ratio  = median(ratio, na.rm = TRUE),
    pct_increased = mean(ratio > 1, na.rm = TRUE),
    pct_gt2x      = mean(ratio > 2, na.rm = TRUE),
    pct_lt_half   = mean(ratio < 0.5, na.rm = TRUE),
    .groups = "drop"
  ) %>% print()

p_var <- ggplot(gene_var_df, aes(x = log2(ratio), fill = adjuster, colour = adjuster)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100,
                 alpha = 0.45, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values   = c(combat = "#4e79a7", combat_sup = "#e15759")) +
  scale_colour_manual(values = c(combat = "#4e79a7", combat_sup = "#e15759")) +
  labs(
    title    = sprintf("Per-gene variance change: %s (n=%s)", test_study, n_str),
    subtitle = "log2(var_adjusted / var_unadjusted) per gene. Right = increased, Left = decreased.",
    x = "log2 variance ratio (adjusted / unadjusted)", y = "Density"
  ) +
  theme_minimal(base_size = 12)

png(file.path(OUT_DIR, "gene_variance.png"), width = 8*150, height = 5*150, res = 150)
print(p_var)
dev.off()
cat("Saved gene_variance.png\n")

# ── 3. Top-PC alignment ───────────────────────────────────────────────────────
# Cosine similarity between top PCs of each pair.
# |V1^T V2| measures how much each PC in V1 aligns with any PC in V2.

pc_align <- function(V1, V2, label1, label2) {
  sim <- abs(t(V1) %*% V2)   # N_PCS × N_PCS absolute cosine similarity
  data.frame(
    pc1   = rep(seq_len(nrow(sim)), ncol(sim)),
    pc2   = rep(seq_len(ncol(sim)), each = nrow(sim)),
    sim   = as.vector(sim),
    pair  = sprintf("%s vs %s", label1, label2)
  )
}

align_df <- rbind(
  pc_align(loadings_list[["unadjusted"]], loadings_list[["combat"]],     "unadjusted", "combat"),
  pc_align(loadings_list[["unadjusted"]], loadings_list[["combat_sup"]], "unadjusted", "combat_sup"),
  pc_align(loadings_list[["combat"]],     loadings_list[["combat_sup"]], "combat",     "combat_sup")
)

# Summarise: for each PC of matrix 1, what is the max alignment to any PC of matrix 2?
cat("\n=== Top-PC alignment: max cosine similarity of each top PC ===\n")
align_df %>%
  group_by(pair, pc1) %>%
  summarise(max_sim = max(sim), .groups = "drop") %>%
  group_by(pair) %>%
  summarise(
    mean_max_sim   = mean(max_sim),
    pct_well_aligned = mean(max_sim > 0.9),
    .groups = "drop"
  ) %>% print()

# Plot: heatmap of top-20 PC alignment
p_align <- align_df %>%
  filter(pc1 <= 20, pc2 <= 20) %>%
  ggplot(aes(x = pc2, y = pc1, fill = sim)) +
  geom_tile() +
  facet_wrap(~pair) +
  scale_fill_gradient(low = "white", high = "#e15759", limits = c(0, 1),
                      name = "|cos sim|") +
  labs(
    title    = sprintf("Top-PC alignment (top 20 PCs): %s (n=%s)", test_study, n_str),
    subtitle = "Diagonal = identity (same PCs). Off-diagonal structure = geometric rotation.",
    x = "PC of matrix 2", y = "PC of matrix 1"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

png(file.path(OUT_DIR, "pc_alignment.png"), width = 12*150, height = 4.5*150, res = 150)
print(p_align)
dev.off()
cat("Saved pc_alignment.png\n")

cat(sprintf("\nAll outputs in: %s\n", OUT_DIR))
