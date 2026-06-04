#!/usr/bin/env Rscript
# check_axis_contamination.R
# Tests whether the supervised class axis is contaminated by batch.
#
# Setup: confounded training pair GSE37250_SA (mostly Active TB) + India (mostly
# Control), test study USA. All in shared log-gene space.
#
# Direction vectors (per gene):
#   class_adj    = batch-ADJUSTED class effect (ComBat's B.hat covariate row:
#                  group coefficient from design = [batch indicators, group])
#   class_naive  = mean(Active) - mean(Control), pooled, NOT batch-adjusted
#   batch_contr  = mean(SA) - mean(India)   (the dominant batch axis)
#   class_USA    = mean(USA Active) - mean(USA Control)  (test's own class axis)
#
# Cosines that matter:
#   cos(class_naive, batch_contr)  — raw confounding (expect HIGH: SA=active)
#   cos(class_adj,   batch_contr)  — residual contamination after batch adjust
#   cos(class_adj,   class_USA)    — does supervised axis generalize to test?
#   cos(class_naive, class_USA)    — does naive class axis generalize to test?

suppressMessages({ library(ggplot2) })

DATA_FILE <- "data/TB_real_data.RData"
OUT_DIR   <- "outputs/diagnostics/combat_sup_knn"

load(DATA_FILE)

bin <- function(l) as.integer(ifelse(l %in% c("1", 1, "Active"), 1L, 0L))

train_studies <- c("GSE37250_SA", "India")   # confounded pair (exclude USA overlap)
test_study    <- "USA"
all_used      <- c(train_studies, test_study)

# common genes across all involved studies
cg <- Reduce(intersect, lapply(dat_lst[all_used], rownames))
cat(sprintf("Common genes: %d\n", length(cg)))

# build combined training matrix + log transform (mirror pipeline heuristic)
maybe_log <- function(m) {
  if (max(m, na.rm = TRUE) > 100) {
    mn <- min(m, na.rm = TRUE); if (mn < 0) m <- m - mn
    log2(m + 1)
  } else m
}

dat_SA    <- maybe_log(dat_lst[["GSE37250_SA"]][cg, , drop = FALSE])
dat_India <- maybe_log(dat_lst[["India"]][cg, , drop = FALSE])
dat_USA   <- maybe_log(dat_lst[["USA"]][cg, , drop = FALSE])

lab_SA    <- bin(label_lst[["GSE37250_SA"]])
lab_India <- bin(label_lst[["India"]])
lab_USA   <- bin(label_lst[["USA"]])

cat(sprintf("SA:    %d Active / %d Control\n", sum(lab_SA==1),    sum(lab_SA==0)))
cat(sprintf("India: %d Active / %d Control\n", sum(lab_India==1), sum(lab_India==0)))
cat(sprintf("USA:   %d Active / %d Control\n", sum(lab_USA==1),   sum(lab_USA==0)))

# training matrix and metadata
dat_train   <- cbind(dat_SA, dat_India)
batch_train <- c(rep("SA", ncol(dat_SA)), rep("India", ncol(dat_India)))
group_train <- c(lab_SA, lab_India)

# ── direction vectors ─────────────────────────────────────────────────────────

# batch-adjusted class effect via least squares (ComBat's B.hat covariate row)
design <- cbind(model.matrix(~ -1 + factor(batch_train)), group = group_train)
Bhat   <- solve(crossprod(design), crossprod(design, t(dat_train)))  # params × genes
class_adj <- Bhat["group", ]

# naive (non-batch-adjusted) class effect
class_naive <- rowMeans(dat_train[, group_train == 1, drop = FALSE]) -
               rowMeans(dat_train[, group_train == 0, drop = FALSE])

# batch contrast
batch_contr <- rowMeans(dat_SA) - rowMeans(dat_India)

# test study's own class axis
class_USA <- rowMeans(dat_USA[, lab_USA == 1, drop = FALSE]) -
             rowMeans(dat_USA[, lab_USA == 0, drop = FALSE])

# ── cosine similarity ─────────────────────────────────────────────────────────
cosine <- function(a, b) sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))

vecs <- list(
  class_adj   = class_adj,
  class_naive = class_naive,
  batch_contr = batch_contr,
  class_USA   = class_USA
)

cat("\n=== Cosine similarity matrix between direction vectors ===\n")
M <- outer(names(vecs), names(vecs), Vectorize(function(i, j) cosine(vecs[[i]], vecs[[j]])))
dimnames(M) <- list(names(vecs), names(vecs))
print(round(M, 3))

cat("\n=== Key comparisons ===\n")
cat(sprintf("  cos(class_naive, batch_contr) = %+.3f   raw confounding (SA=active, India=control)\n",
            cosine(class_naive, batch_contr)))
cat(sprintf("  cos(class_adj,   batch_contr) = %+.3f   residual batch contamination after adjustment\n",
            cosine(class_adj, batch_contr)))
cat(sprintf("  cos(class_adj,   class_USA)   = %+.3f   does SUPERVISED axis generalize to test?\n",
            cosine(class_adj, class_USA)))
cat(sprintf("  cos(class_naive, class_USA)   = %+.3f   does NAIVE class axis generalize to test?\n",
            cosine(class_naive, class_USA)))

# ── plot: cosine matrix heatmap ───────────────────────────────────────────────
df <- expand.grid(v1 = names(vecs), v2 = names(vecs))
df$cos <- mapply(function(i, j) cosine(vecs[[as.character(i)]], vecs[[as.character(j)]]),
                 df$v1, df$v2)

p <- ggplot(df, aes(v1, v2, fill = cos)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%+.2f", cos)), size = 4) +
  scale_fill_gradient2(low = "#4e79a7", mid = "white", high = "#e15759",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title    = "Direction-vector alignment (cosine similarity)",
       subtitle = "SA+India training (confounded) vs USA test. Class axis vs batch axis vs test axis.",
       x = NULL, y = NULL, fill = "cosine") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

png(file.path(OUT_DIR, "axis_contamination.png"), width = 7*150, height = 6*150, res = 150)
print(p); dev.off()
cat(sprintf("\nPlot: %s/axis_contamination.png\n", OUT_DIR))
