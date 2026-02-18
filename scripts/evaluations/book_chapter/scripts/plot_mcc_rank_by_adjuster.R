#!/usr/bin/env Rscript

# plot_mcc_rank_spotlight.R
# Plot MCC ranks with adjuster on x-axis, rank on y-axis.

options(warn = -1)
suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(ggbeeswarm)
})

# Load utils if present
if (file.exists("scripts/adjuster_plot_utils.R")) {
  source("scripts/adjuster_plot_utils.R")
}

parser <- ArgumentParser(description = "Plot MCC rank by adjuster, highlighting outliers")
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Input CSV file with adjuster results")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output PNG file")
parser$add_argument("--width", type = "double", default = 14)
parser$add_argument("--height", type = "double", default = 9)
parser$add_argument("--dpi", type = "integer", default = 300)
parser$add_argument("--adjusters", type = "character", default = NULL,
                    help = "Comma-separated list of adjusters to include")
parser$add_argument("--n-datasets", type = "character", default = "all",
                    help = "Number of datasets to filter on, or 'all'")
parser$add_argument("--debug", action = "store_true", default = FALSE,
                    help = "Enable debug prints")

args <- parser$parse_args()

# ==============================================================================
# Setup Helpers & Constants
# ==============================================================================

# Define label mapping
classifier_labels_map <- c(
  "rda" = "RDA", "logistic" = "Logistic", "elasticnet" = "ElasticNet",
  "svm" = "SVM", "rf" = "Random Forest", "knn" = "KNN",
  "xgboost" = "XGBoost", "nnet" = "Neural Net", "shrinkageLDA" = "LDA"
)

# Parse dynamic filters
target_n <- if (tolower(args$n_datasets) != "all") as.integer(args$n_datasets) else NULL
target_adjusters <- if (!is.null(args$adjusters)) trimws(strsplit(args$adjusters, ",")[[1]]) else NULL

# Set plot title suffix
n_datasets_label <- if (is.null(target_n)) "All Studies" else sprintf("%d-Study", target_n)

if (args$debug) {
  cat("DEBUG: Starting pipeline with n_datasets:", n_datasets_label, "\n")
  cat("DEBUG: target_n:", target_n, "(NULL check:", is.null(target_n), ")\n")
}

# ==============================================================================
# Main Pipeline: Load -> Clean -> Rank -> Aggregate -> Format
# ==============================================================================

plot_data <- read.csv(args$input, stringsAsFactors = FALSE) %>%
  # Clean whitespace
  mutate(across(c(adjuster, classifier), trimws)) %>%

  # Filter: Metrics, Exclusions, and Dynamic Args
  filter(
    metric == "mcc",
    !is.na(value),
    !is.na(n_datasets),
    adjuster != "within_study_cv",
    !classifier %in% c("logistic", "rda"),
    if (is.null(target_n)) TRUE else n_datasets == target_n,
    is.null(target_adjusters) | adjuster %in% target_adjusters
  ) %>%

  # Calculate Rank (Rank 1 = Highest MCC)
  mutate(
    rank = rank(-value, ties.method = "average"),
    .by = c(classifier, n_datasets, test_study)
  ) %>%

  # Aggregate Stats (IQR and Mean)
  summarise(
    avg_rank = mean(rank, na.rm = TRUE),
    q25 = quantile(rank, 0.25, na.rm = TRUE),
    q75 = quantile(rank, 0.75, na.rm = TRUE),
    n_obs = n(),
    .by = c(classifier, adjuster)
  ) %>%

  # Apply Logic: Labels
  mutate(
    # Map classifier labels
    classifier_label = recode(classifier, !!!classifier_labels_map),

    # Apply external formatter if exists
    adjuster_label = if (exists("format_adjuster_label")) {
      sapply(adjuster, format_adjuster_label)
    } else {
      adjuster
    }
  ) %>%

  # Apply Logic: Outlier Stats and Ordering
  mutate(
    median_rank_group = median(avg_rank, na.rm = TRUE),
    .by = adjuster_label
  ) %>%
  mutate(
    dist_from_median = abs(avg_rank - median_rank_group),
    global_dist_iqr = IQR(dist_from_median, na.rm = TRUE),
    standardized_dev = dist_from_median / global_dist_iqr,
    dist_for_alpha = pmin(standardized_dev, 2.5)
  ) %>%
  # Sort to find top outliers
  arrange(desc(dist_from_median)) %>%
  mutate(
    rank_outlier = row_number(),
    label_text = ifelse(rank_outlier <= 6, as.character(classifier_label), ""),
    # Reorder Adjuster factor by performance (mean avg_rank)
    adjuster_label = reorder(adjuster_label, avg_rank, FUN = mean),
    # Ensure Classifier legend is consistent (alphabetical or specific order)
    classifier_label = factor(classifier_label, levels = sort(unique(classifier_label)))
  )

# Calculate axis limits
max_rank <- ceiling(max(plot_data$q75, na.rm = TRUE))
if (args$debug) cat("DEBUG: Calculated max_rank for y-axis:", max_rank, "\n")

# ==============================================================================
# Plotting
# ==============================================================================

pos_hex <- position_beeswarm(
  method = "hex",
  cex = 1.5,
  groupOnX = TRUE
)

p <- ggplot(plot_data, aes(x = adjuster_label, y = avg_rank)) +

  # A. Background Mass (Boxplot)
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.1,
    color = "grey80",
    fill = "grey90"
  ) +

  # B. Error Bars (IQR)
  geom_errorbar(
    aes(ymin = q25,
        ymax = q75,
        alpha = dist_for_alpha,
        group = adjuster_label),
    position = pos_hex,
    width = 0,
    linewidth = 0.4,
    color = "grey60"
  ) +

  # C. Points
  geom_point(
    aes(color = classifier_label,
        alpha = dist_for_alpha,
        group = adjuster_label),
    position = pos_hex,
    size = 5.0,
  ) +

  # D. Outlier Labels
  geom_text(
    aes(label = label_text,
        alpha = dist_for_alpha,
        group = adjuster_label),
    position = pos_hex,
    vjust = -0.8,
    hjust = -0.3,
    size = 3.5,
    fontface = "bold",
    color = "black",
    show.legend = FALSE
  ) +

  scale_alpha_continuous(range = c(0.4, 1.0), guide = "none") +

  scale_y_reverse(
    breaks = seq(1, max_rank, by = 1),
    limits = c(max_rank + 0.5, 0.5)
  ) +

  scale_x_discrete(expand = expansion(add = 1.0)) +
  coord_cartesian(clip = "off") +

  labs(
    title = sprintf("MCC Rank by Adjuster (%s)", n_datasets_label),
    subtitle = "Hex-packed by Rank. Bars show 25th-75th Percentiles. Labels show top 6 outliers.",
    x = NULL,
    y = "Average Rank (Lower is Better)",
    color = "Classifier"
  ) +

  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", color = "black"),
    panel.grid.major.x = element_line(color = "gray95"),
    panel.grid.major.y = element_line(color = "gray92", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  )

ggsave(args$output, p, width = args$width, height = args$height, dpi = args$dpi, bg = "white")
if (args$debug) cat("DEBUG: Saved plot to:", args$output, "\n")