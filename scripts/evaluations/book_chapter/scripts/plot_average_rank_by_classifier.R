#!/usr/bin/env Rscript

# plot_average_rank_by_classifier.R
# Generate violin plots showing MCC distribution per classifier with adjuster means
# X-axis: classifier, Y-axis: MCC, Violin: full distribution, Points: adjuster means

options(warn = -1)
suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("scripts/adjuster_plot_utils.R")

parser <- ArgumentParser(description = "Create MCC violin plot by classifier with adjuster means")
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Input CSV file with adjuster results")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output PNG file")
parser$add_argument("--width", type = "double", default = 12)
parser$add_argument("--height", type = "double", default = 8)
parser$add_argument("--dpi", type = "integer", default = 300)
parser$add_argument("--adjusters", type = "character", default = NULL,
                    help = "Comma-separated list of adjusters to include")
parser$add_argument("--n-datasets", type = "character", default = "4",
                    help = "Number of datasets to filter on, or 'all' to average across all (default: 4)")

args <- parser$parse_args()

# Load data
data <- read.csv(args$input, stringsAsFactors = FALSE)

# Filter to MCC metric and remove NA classifiers
data <- data %>%
  filter(metric == "mcc", !is.na(value), !is.na(n_datasets), !is.na(classifier))

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

data$classifier_label <- classifier_labels[data$classifier]

# Format adjuster labels
data$adjuster_label <- sapply(data$adjuster, format_adjuster_label)

# Order classifiers by median MCC (best performers first)
classifier_order <- data %>%
  group_by(classifier, classifier_label) %>%
  summarise(median_mcc = median(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_mcc))

data$classifier_label <- factor(
  data$classifier_label,
  levels = classifier_order$classifier_label
)

# Calculate adjuster means per classifier for the points
adjuster_means <- data %>%
  group_by(classifier_label, adjuster_label) %>%
  summarise(
    mean_mcc = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

# Build the violin plot
p <- ggplot(data, aes(x = classifier_label, y = value)) +
  # Violin plot showing full distribution
  geom_violin(fill = "gray85", color = "gray50", alpha = 0.7, scale = "width") +
  # Add boxplot for quartiles
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.5, outlier.shape = NA) +
  # Add colored points for adjuster means
  geom_point(
    data = adjuster_means,
    aes(x = classifier_label, y = mean_mcc, color = adjuster_label),
    size = 4,
    position = position_dodge(width = 0.3)
  ) +
  # Styling
  scale_color_brewer(palette = "Set1", name = "Adjuster Mean") +
  labs(
    title = sprintf("MCC Distribution by Classifier (%s)", n_datasets_label),
    subtitle = "Violin shows full distribution | Box shows quartiles | Points show adjuster means",
    x = "Classifier",
    y = "MCC"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

# Save the plot
ggsave(args$output, p, width = args$width, height = args$height, dpi = args$dpi, bg = "white")
cat("Saved plot to:", args$output, "\n")
