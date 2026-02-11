#!/usr/bin/env Rscript

# plot_average_rank_by_classifier.R
# Generate a radar chart showing average rank of adjusters per classifier for the 4-study case
# Axes: classifiers (radial), Rank 1 on outside, error bands as semi-transparent areas

options(warn = -1)
suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("scripts/adjuster_plot_utils.R")

parser <- ArgumentParser(description = "Create average rank radar chart for 4-study case")
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Input CSV file with adjuster results")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output PNG file")
parser$add_argument("--width", type = "double", default = 10)
parser$add_argument("--height", type = "double", default = 10)
parser$add_argument("--dpi", type = "integer", default = 300)
parser$add_argument("--adjusters", type = "character", default = NULL,
                    help = "Comma-separated list of adjusters to include")
parser$add_argument("--n-datasets", type = "character", default = "4",
                    help = "Number of datasets to filter on, or 'all' to average across all (default: 4)")

args <- parser$parse_args()

# Load data
data <- read.csv(args$input, stringsAsFactors = FALSE)

# Filter to MCC metric
data <- data %>%
  filter(metric == "mcc", !is.na(value), !is.na(n_datasets))

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

# For each classifier, n_datasets, and test_study, rank the adjusters by MCC (higher is better, so rank 1 = best)
ranked_data <- data %>%
  group_by(classifier, n_datasets, test_study) %>%
  mutate(rank = rank(-value, ties.method = "average")) %>%
  ungroup()

# Calculate average rank per adjuster per classifier (averaged over all test studies and n_datasets)
avg_rank <- ranked_data %>%
  group_by(classifier, adjuster) %>%
  summarise(
    avg_rank = mean(rank, na.rm = TRUE),
    se_rank = sd(rank, na.rm = TRUE) / sqrt(n()),
    n_obs = n(),
    .groups = "drop"
  )

# Create nice labels
classifier_labels <- c(
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

# Get max rank for inversion (rank 1 = outside)
max_rank <- max(avg_rank$avg_rank + avg_rank$se_rank, na.rm = TRUE)

# Invert ranks so rank 1 is on outside
avg_rank <- avg_rank %>%
  mutate(
    inverted_rank = max_rank + 1 - avg_rank,
    inverted_lower = max_rank + 1 - (avg_rank + se_rank),
    inverted_upper = max_rank + 1 - (avg_rank - se_rank)
  )

# Prepare data for radar chart using coord_polar
# Need to close the polygon by repeating the first classifier
classifiers_ordered <- names(classifier_labels)
n_classifiers <- length(classifiers_ordered)

# Assign numeric positions for classifiers
classifier_pos <- setNames(1:n_classifiers, classifiers_ordered)
avg_rank$classifier_num <- classifier_pos[avg_rank$classifier]

# Close the loop: duplicate first point at the end for each adjuster
close_loop <- avg_rank %>%
  filter(classifier == classifiers_ordered[1]) %>%
  mutate(classifier_num = n_classifiers + 1)

radar_data <- bind_rows(avg_rank, close_loop)

# Create color palette
n_adjusters <- length(unique(radar_data$adjuster_label))
colors <- scales::hue_pal()(n_adjusters)

# Build the radar chart
p <- ggplot(radar_data, aes(x = classifier_num, group = adjuster_label, color = adjuster_label, fill = adjuster_label)) +
  # Error bands as ribbons
  geom_ribbon(aes(ymin = inverted_lower, ymax = inverted_upper), alpha = 0.15, color = NA) +
  # Lines connecting points
  geom_line(aes(y = inverted_rank), linewidth = 0.8) +
  # Points
  geom_point(aes(y = inverted_rank), size = 2) +
  # Convert to polar coordinates
  coord_polar(start = -pi / n_classifiers) +
  # Set axis breaks and labels - extend limits to include the closing segment
  scale_x_continuous(
    breaks = 1:n_classifiers,
    labels = classifier_labels[classifiers_ordered],
    limits = c(1, n_classifiers + 1)
  ) +
  scale_y_continuous(
    limits = c(0, max_rank + 1),
    breaks = seq(0, max_rank + 1, by = 2),
    labels = function(x) round(max_rank + 1 - x)  # Show original rank values
  ) +
  labs(
    title = sprintf("Average Adjuster Rank by Classifier (%s)", n_datasets_label),
    subtitle = "Rank 1 (best) on outside; shaded areas show ± SE",
    color = "Adjuster",
    fill = "Adjuster"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_line(color = "gray80"),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save the plot
ggsave(args$output, p, width = args$width, height = args$height, dpi = args$dpi, bg = "white")
cat("Saved plot to:", args$output, "\n")
