#!/usr/bin/env Rscript

# plot_rank_ridgeline.R
# Generate a ridgeline plot showing rank distributions per adjuster
# Each ridge is an adjuster, x-axis shows rank (1 = best)

options(warn = -1)
suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(ggridges)
})

source("scripts/adjuster_plot_utils.R")

parser <- ArgumentParser(description = "Create ridgeline plot of adjuster rank distributions")
parser$add_argument("-i", "--input", type = "character", required = TRUE,
                    help = "Input CSV file with adjuster results")
parser$add_argument("-o", "--output", type = "character", required = TRUE,
                    help = "Output PNG file")
parser$add_argument("--width", type = "double", default = 10)
parser$add_argument("--height", type = "double", default = 10)
parser$add_argument("--dpi", type = "integer", default = 300)
parser$add_argument("--adjusters", type = "character", default = NULL,
                    help = "Comma-separated list of adjusters to include")
parser$add_argument("--n-datasets", type = "character", default = "all",
                    help = "Number of datasets to filter on, or 'all' (default: all)")

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

# For each classifier, n_datasets, and test_study, rank the adjusters by MCC
ranked_data <- data %>%
  group_by(classifier, n_datasets, test_study) %>%
  mutate(rank = rank(-value, ties.method = "average")) %>%
  ungroup()

# Format adjuster labels
ranked_data$adjuster_label <- sapply(ranked_data$adjuster, format_adjuster_label)

# Order adjusters by overall average rank (best performers at top)
adjuster_order <- ranked_data %>%
  group_by(adjuster, adjuster_label) %>%
  summarise(overall_avg = mean(rank, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall_avg))  # Reversed so best is at top in ridgeline

ranked_data$adjuster_label <- factor(
  ranked_data$adjuster_label,
  levels = adjuster_order$adjuster_label
)

# Get number of adjusters for x-axis limits
n_adjusters <- length(unique(ranked_data$adjuster))

# Build the ridgeline plot - x-axis is rank, y-axis is adjuster
p <- ggplot(ranked_data, aes(x = rank, y = adjuster_label, fill = stat(x))) +
  geom_density_ridges_gradient(
    scale = 2,
    rel_min_height = 0.01,
    gradient_lwd = 0.5
  ) +
  scale_x_continuous(
    breaks = 1:n_adjusters,
    limits = c(0.5, n_adjusters + 0.5),
    expand = c(0, 0)
  ) +
  scale_fill_viridis_c(name = "Rank", option = "plasma", direction = -1) +
  labs(
    x = "Rank (1 = best)",
    y = "Adjuster",
    title = sprintf("Rank Distribution by Adjuster (%s)", n_datasets_label),
    subtitle = "Distribution of ranks across all classifiers and test studies"
  ) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Save the plot
ggsave(args$output, p, width = args$width, height = args$height, dpi = args$dpi, bg = "white")
cat("Saved plot to:", args$output, "\n")
