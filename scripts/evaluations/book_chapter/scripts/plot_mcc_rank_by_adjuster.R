#!/usr/bin/env Rscript

# plot_mcc_rank_by_adjuster.R
# Plot MCC ranks with adjuster on x-axis, rank on y-axis (rank 1 at top)
# Colored by classifier, averaged over test sets with standard error bars
# Classifiers staggered within adjuster columns

options(warn = -1)
suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("scripts/adjuster_plot_utils.R")

parser <- ArgumentParser(description = "Plot MCC rank by adjuster, colored by classifier")
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Input CSV file with adjuster results")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output PNG file")
parser$add_argument("--width", type = "double", default = 12)
parser$add_argument("--height", type = "double", default = 8)
parser$add_argument("--dpi", type = "integer", default = 300)
parser$add_argument("--adjusters", type = "character", default = NULL,
                    help = "Comma-separated list of adjusters to include")
parser$add_argument("--n-datasets", type = "character", default = "all",
                    help = "Number of datasets to filter on, or 'all' (default: all)")

args <- parser$parse_args()

# Load data
data <- read.csv(args$input, stringsAsFactors = FALSE)

# Filter to MCC metric, exclude Within-study CV and logistic classifier
data <- data %>%
  filter(
    metric == "mcc", 
    !is.na(value), 
    !is.na(n_datasets),
    adjuster != "within_study_cv",
    classifier != "logistic"
  )

# Filter to specified n_datasets unless "all"
n_datasets_label <- args$n_datasets
if (tolower(args$n_datasets) != "all") {
  n_val <- as.integer(args$n_datasets)
  data <- data %>% filter(n_datasets == n_val)
  n_datasets_label <- sprintf("%d-Study", n_val)
} else {
  n_datasets_label <- "All Studies"
}

# Filter adjusters if specified
if (!is.null(args$adjusters)) {
  adjusters_filter <- trimws(strsplit(args$adjusters, ",")[[1]])
  data <- data %>% filter(adjuster %in% adjusters_filter)
}

# For each classifier, n_datasets, and test_study, rank the adjusters by MCC
# Higher MCC is better, so rank 1 = best
ranked_data <- data %>%
  group_by(classifier, n_datasets, test_study) %>%
  mutate(rank = rank(-value, ties.method = "average")) %>%
  ungroup()

# Calculate average rank per adjuster per classifier (averaged over test studies)
avg_rank <- ranked_data %>%
  group_by(classifier, adjuster) %>%
  summarise(
    avg_rank = mean(rank, na.rm = TRUE),
    se_rank = sd(rank, na.rm = TRUE) / sqrt(n()),
    n_obs = n(),
    .groups = "drop"
  )

# Create nice classifier labels
classifier_labels <- c(
  "rda" = "RDA",
  "logistic" = "Logistic",
  "elasticnet" = "ElasticNet", 
  "svm" = "SVM",
  "rf" = "Random Forest",
  "knn" = "KNN",
  "xgboost" = "XGBoost",
  "nnet" = "Neural Net",
  "shrinkageLDA" = "Shrinkage LDA"
)

avg_rank$classifier_label <- classifier_labels[avg_rank$classifier]
avg_rank$classifier_label[is.na(avg_rank$classifier_label)] <- avg_rank$classifier[is.na(avg_rank$classifier_label)]

# Format adjuster labels
avg_rank$adjuster_label <- sapply(avg_rank$adjuster, format_adjuster_label)

# Order adjusters by overall average rank (best performers first)
adjuster_order <- avg_rank %>%
  group_by(adjuster, adjuster_label) %>%
  summarise(overall_avg = mean(avg_rank), .groups = "drop") %>%
  arrange(overall_avg)

avg_rank$adjuster_label <- factor(
  avg_rank$adjuster_label,
  levels = adjuster_order$adjuster_label
)

# Order classifiers consistently
classifier_order <- unique(avg_rank$classifier_label)
avg_rank$classifier_label <- factor(avg_rank$classifier_label, levels = classifier_order)

# Create numeric x positions for adjusters
adjuster_levels <- levels(avg_rank$adjuster_label)
n_adjusters <- length(adjuster_levels)
avg_rank$adjuster_num <- as.numeric(avg_rank$adjuster_label)

# Stagger classifiers within each adjuster column
n_classifiers <- length(unique(avg_rank$classifier_label))
stagger_width <- 0.3  # Total width for staggering
stagger_offsets <- seq(-stagger_width/2, stagger_width/2, length.out = n_classifiers)
names(stagger_offsets) <- levels(avg_rank$classifier_label)

avg_rank <- avg_rank %>%
  mutate(x_staggered = adjuster_num + stagger_offsets[as.character(classifier_label)])

# Get max rank for y-axis
max_rank <- ceiling(max(avg_rank$avg_rank + avg_rank$se_rank, na.rm = TRUE))

# Build the plot
p <- ggplot(avg_rank, aes(x = x_staggered, y = avg_rank, color = classifier_label)) +
  # Error bars
  geom_errorbar(
    aes(ymin = avg_rank - se_rank, ymax = avg_rank + se_rank),
    width = 0.1,
    linewidth = 0.5
  ) +
  # Points
  geom_point(size = 3) +
  # Reverse y-axis so rank 1 is at top
  scale_y_reverse(
    breaks = seq(1, max_rank, by = 1),
    limits = c(max_rank + 0.5, 0.5)
  ) +
  # X-axis with adjuster labels
  scale_x_continuous(
    breaks = 1:n_adjusters,
    labels = adjuster_levels,
    limits = c(0.5, n_adjusters + 0.5)
  ) +
  labs(
    title = sprintf("MCC Rank by Adjuster (%s)", n_datasets_label),
    subtitle = "Rank 1 (best) at top; error bars show ± SE across test sets",
    x = "Adjuster",
    y = "Average Rank",
    color = "Classifier"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save the plot
ggsave(args$output, p, width = args$width, height = args$height, dpi = args$dpi, bg = "white")
cat("Saved plot to:", args$output, "\n")
