#!/usr/bin/env Rscript

# plot_mcc_rank_spotlight.R
# Plot classifier rank distributions horizontally by adjuster with centered labels.

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

parser <- ArgumentParser(description = "Plot MCC rank by adjuster (Horizontal), highlighting outliers")
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Input CSV file with adjuster results")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output PNG file")
parser$add_argument("--width", type = "double", default = 11)
parser$add_argument("--height", type = "double", default = 9)
parser$add_argument("--dpi", type = "integer", default = 300)
parser$add_argument("--adjusters", type = "character", default = NULL,
                    help = "Comma-separated list of adjusters to include")
parser$add_argument("--n-datasets", type = "character", default = "all",
                    help = "Number of datasets to filter on, or 'all'")
parser$add_argument("--n-labeled-outliers", type = "integer", default = 8,
                    help = "Number of top outliers to label")

args <- parser$parse_args()

# ==============================================================================
# Constants & Helpers
# ==============================================================================

classifier_labels_map <- c(
  "rda" = "RDA", "logistic" = "Logistic Regression", "elasticnet" = "ElasticNet",
  "svm" = "SVM", "rf" = "Random Forest", "knn" = "KNN",
  "xgboost" = "XGBoost", "nnet" = "Neural Net", "shrinkageLDA" = "Shrinkage LDA"
)

shortened_labels_map <- c(
  "rda" = "RDA", "logistic" = "LR", "elasticnet" = "ENet",
  "svm" = "SVM", "rf" = "RF", "knn" = "KNN",
  "xgboost" = "XGB", "nnet" = "NN", "shrinkageLDA" = "LDA"
)

target_n <- if (tolower(args$n_datasets) != "all") as.integer(args$n_datasets) else NULL
target_adjusters <- if (!is.null(args$adjusters)) trimws(strsplit(args$adjusters, ",")[[1]]) else NULL

# Calculate Hodges-Lehmann estimator (Pseudomedian)
calc_hodges_lehmann <- function(x) {
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA)
  return(median(outer(x, x, "+") / 2))
}

# ==============================================================================
# Data Processing
# ==============================================================================

# Load and clean raw data
raw_data <- read.csv(args$input, stringsAsFactors = FALSE) %>%
  mutate(across(c(adjuster, classifier), trimws)) %>%
  filter(
    metric == "mcc",
    !is.na(value),
    !is.na(n_datasets),
    adjuster != "within_study_cv",
    !classifier %in% c("logistic", "rda"),
    if (is.null(target_n)) TRUE else n_datasets == target_n,
    is.null(target_adjusters) | adjuster %in% target_adjusters
  )

# Calculate ranks per study
ranked_data <- raw_data %>%
  mutate(
    rank = rank(-value, ties.method = "average"),
    .by = c(classifier, n_datasets, test_study)
  ) %>%
  mutate(
    classifier_label = recode(classifier, !!!classifier_labels_map),
    adjuster_label = if (exists("format_adjuster_label")) {
      sapply(adjuster, format_adjuster_label)
    } else {
      adjuster
    }
  )

max_rank <- max(ranked_data$rank, na.rm = TRUE)

# Aggregate stats and identify outliers
outlier_stats <- ranked_data %>%
  mutate(
    group_center = calc_hodges_lehmann(rank),
    .by = adjuster_label
  ) %>%
  summarise(
    avg_rank = calc_hodges_lehmann(rank),
    group_center = first(group_center),
    classifier = first(classifier),
    .by = c(classifier_label, adjuster_label)
  ) %>%
  mutate(
    abs_dev = abs(avg_rank - group_center),
    global_iqr = IQR(abs_dev, na.rm = TRUE),
    z_score_dev = abs_dev / global_iqr,
    highlight_intensity = z_score_dev
  ) %>%
  arrange(desc(abs_dev)) %>%
  mutate(
    outlier_rank = row_number(),
    label_text = ifelse(outlier_rank <= args$n_labeled_outliers, 
                       as.character(recode(classifier, !!!shortened_labels_map)), "")
  )

# ==============================================================================
# Synchronize Factor Levels
# ==============================================================================

adjuster_order <- ranked_data %>%
  group_by(adjuster_label) %>%
  summarise(
    med = median(rank, na.rm = TRUE),
    hl_avg = calc_hodges_lehmann(rank)
  ) %>%
  arrange(desc(med), desc(hl_avg)) %>% 
  pull(adjuster_label)

ranked_data$adjuster_label <- factor(ranked_data$adjuster_label, levels = adjuster_order)
outlier_stats$adjuster_label <- factor(outlier_stats$adjuster_label, levels = adjuster_order)
outlier_stats$classifier_label <- factor(outlier_stats$classifier_label, 
                                         levels = sort(unique(outlier_stats$classifier_label)))

# ==============================================================================
# Visualization
# ==============================================================================

# Shared beeswarm position logic
pos_swarm_h <- position_beeswarm(method = "hex", cex = 1.5, groupOnX = FALSE)

p <- ggplot(mapping = aes(y = adjuster_label)) +

  # Background distribution
  geom_violin(
    data = ranked_data, 
    aes(x = rank),
    width = 0.7,
    alpha = 0.1,
    color = "grey80",
    fill = "grey90",
    trim = TRUE
  ) +

  # HL center lines
  stat_summary(
    data = ranked_data,
    aes(x = rank),
    fun = calc_hodges_lehmann,
    fun.min = calc_hodges_lehmann,
    fun.max = calc_hodges_lehmann,
    geom = "crossbar",
    width = 0.5,
    size = 0.4, 
    color = "grey40"
  ) +

  # HL ranks per classifier
  geom_point(
    data = outlier_stats,
    aes(x = avg_rank,
        color = classifier_label,
        alpha = highlight_intensity,
        group = adjuster_label),
    position = pos_swarm_h,
    size = 5.0,
  ) +

  # Labels centered under points
  geom_text(
    data = filter(outlier_stats, label_text != ""),
    aes(x = avg_rank, label = label_text, group = adjuster_label),
    position = pos_swarm_h,
    size = 3.0,
    fontface = "bold",
    color = "black",
    vjust = 2.0,
    show.legend = FALSE
  ) +

  # Scales and theme
  scale_alpha_continuous(range = c(0.4, 1.0), guide = "none") +
  scale_x_continuous(breaks = seq(1, max_rank, by = 1), expand = expansion(add = 0.5)) +
  scale_y_discrete(expand = expansion(add = 1.0)) +
  coord_cartesian(xlim = c(1, max_rank), clip = "off") +
  
  labs(
    y = NULL,
    x = "Adjuster Rank Among Other Adjusters (1 = Best)",
    color = "Classifier"
  ) +

  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(face = "bold", color = "black"),
    panel.grid.major.y = element_line(color = "gray95"),
    panel.grid.major.x = element_line(color = "gray92", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 18, margin = margin(b = 5)),
    plot.subtitle = element_text(color = "grey40", size = 11, margin = margin(b = 15))
  )

# Save output
ggsave(args$output, p, width = args$width, height = args$height, dpi = args$dpi, bg = "white")
cat("Saved plot to:", args$output, "\n")