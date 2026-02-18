#!/usr/bin/env Rscript

# plot_mcc_rank_spotlight.R
# Plot classifier rank distributions by adjuster, spotlighting significant outliers.

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
parser$add_argument("--n-labeled-outliers", type = "integer", default = 6,
                    help = "Number of top outliers to label")

args <- parser$parse_args()

# ==============================================================================
# Constants & Helpers
# ==============================================================================

# Map internal classifier names to display labels
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

# Helper to calculate Hodges-Lehmann estimator
calc_hodges_lehmann <- function(x) {
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA)
  # Calculate median of all pairwise averages (Walsh averages)
  # This is robust like median, but smooth like mean
  return(median(outer(x, x, "+") / 2))
}

# ==============================================================================
# Data Processing
# ==============================================================================

# 1. Load and Clean Raw Data
# --------------------------
raw_data <- read.csv(args$input, stringsAsFactors = FALSE) %>%
  mutate(across(c(adjuster, classifier), trimws)) %>%
  # Filter for MCC metric and valid entries
  filter(
    metric == "mcc",
    !is.na(value),
    !is.na(n_datasets),
    adjuster != "within_study_cv",
    !classifier %in% c("logistic", "rda"), # Exclude specific classifiers per requirement
    if (is.null(target_n)) TRUE else n_datasets == target_n,
    is.null(target_adjusters) | adjuster %in% target_adjusters
  )

# 2. Calculate Ranks per Study
# ----------------------------
ranked_data <- raw_data %>%
  # Rank descending (1 = Highest MCC)
  mutate(
    rank = rank(-value, ties.method = "average"),
    .by = c(classifier, n_datasets, test_study)
  ) %>%
  # Add display labels
  mutate(
    classifier_label = recode(classifier, !!!classifier_labels_map),
    adjuster_label = if (exists("format_adjuster_label")) {
      sapply(adjuster, format_adjuster_label)
    } else {
      adjuster
    }
  )

# Calculate max rank for plot limits
max_rank <- max(ranked_data$rank, na.rm = TRUE)

# 3. Aggregate Stats & Identify Outliers
# --------------------------------------
outlier_stats <- ranked_data %>%
  # UPDATE 1: Calculate group center using HL instead of simple median
  # This ensures we measure deviation from the "HL center"
  mutate(
    group_center = calc_hodges_lehmann(rank),
    .by = adjuster_label
  ) %>%
  summarise(
    avg_rank = calc_hodges_lehmann(rank),
    group_center = first(group_center),
    .by = c(classifier_label, adjuster_label)
  ) %>%
  # Compute deviation metrics for spotlight effect
  mutate(
    abs_dev = abs(avg_rank - group_center),
    global_iqr = IQR(abs_dev, na.rm = TRUE),
    z_score_dev = abs_dev / global_iqr,
    highlight_intensity = z_score_dev
  ) %>%
  arrange(desc(abs_dev)) %>%
  mutate(
    outlier_rank = row_number(),
    label_text = ifelse(outlier_rank <= args$n_labeled_outliers, as.character(classifier_label), "")
  )

# ==============================================================================
# 4. Synchronize Factor Levels
# ==============================================================================

# Calculate order: Primary = Median, Secondary = HL (was Mean)
# UPDATE 2: Use HL for secondary sort key
adjuster_order <- ranked_data %>%
  group_by(adjuster_label) %>%
  summarise(
    med = median(rank, na.rm = TRUE),
    hl_avg = calc_hodges_lehmann(rank)
  ) %>%
  arrange(med, hl_avg) %>%  # Sort by Median first, then HL
  pull(adjuster_label)

# Apply the factor levels
ranked_data$adjuster_label <- factor(ranked_data$adjuster_label, levels = adjuster_order)
outlier_stats$adjuster_label <- factor(outlier_stats$adjuster_label, levels = adjuster_order)

# Sort classifier legend alphabetically
outlier_stats$classifier_label <- factor(outlier_stats$classifier_label, 
                                         levels = sort(unique(outlier_stats$classifier_label)))

# ==============================================================================
# Visualization
# ==============================================================================

# Define beeswarm position strategy
pos_swarm <- position_beeswarm(
  method = "hex",
  cex = 1.5,
  groupOnX = TRUE
)

p <- ggplot(mapping = aes(x = adjuster_label)) +

  # Background Distribution (Violins)
  # Show full distribution of raw ranks
  geom_violin(
    data = ranked_data, 
    aes(y = rank),
    width = 0.7,
    alpha = 0.1,
    color = "grey80",
    fill = "grey90",
    trim = TRUE
  ) +

  # Median/HL Lines (Crossbars)
  # UPDATE 3: Use HL for the visual crossbar so it matches the data point logic
  stat_summary(
    data = ranked_data,
    aes(y = rank),
    fun = calc_hodges_lehmann,
    fun.min = calc_hodges_lehmann,
    fun.max = calc_hodges_lehmann,
    geom = "crossbar",
    width = 0.5,
    size = 0.4, 
    color = "grey40"
  ) +

  # HL Ranks (Points)
  # Highlight points based on deviation from group center
  geom_point(
    data = outlier_stats,
    aes(y = avg_rank,
        color = classifier_label,
        alpha = highlight_intensity,
        group = adjuster_label),
    position = pos_swarm,
    size = 5.0,
  ) +

  # Outlier Labels
  # Text for top N deviations
  geom_text(
    data = outlier_stats,
    aes(y = avg_rank,
        # 1. Add 2 spaces for a fixed horizontal buffer
        label = paste0("  ", label_text), 
        group = adjuster_label),
    position = pos_swarm,
    # 2. Use hjust = 0 so the "start" of the string (the spaces) anchors to the point
    hjust = 0,        
    vjust = 0.5,      # Center vertically
    size = 3.5,
    fontface = "bold",
    color = "black",
    show.legend = FALSE
  ) +

  # Scales & Theme
  scale_alpha_continuous(range = c(0.4, 1.0), guide = "none") +

  scale_y_reverse(
    breaks = seq(1, max_rank, by = 1)
  ) +

  scale_x_discrete(expand = expansion(add = 1.0)) +
  coord_cartesian(clip = "off") +

  labs(
    title = sprintf("Robustness of Adjustment Methods to Classifier Selection"),
    subtitle = "Violins: Distribution of Rank across Classification Scenarios\nPoints: Pseudomedian of Adjuster Rank by Classifier",
    x = NULL,
    y = "Adjuster Rank Among Other Adjusters",
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

# Save Output
ggsave(args$output, p, width = args$width, height = args$height, dpi = args$dpi, bg = "white")
cat("Saved plot to:", args$output, "\n")