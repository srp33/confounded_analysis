#!/usr/bin/env Rscript

# plot_mcc_rank_spotlight.R
# Plot MCC ranks with adjuster on x-axis, rank on y-axis (rank 1 at top)
# Highlights outliers:
# - Background boxplot shows the "consensus" mass.
# - Points are LARGE (size 5).
# - Outliers are highlighted via GLOBAL IQR SCALE.
# - Spacing: Uses 'ggbeeswarm' for robust hex/honeycomb packing.

options(warn = -1)
suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
  
  # Dynamic check for ggbeeswarm
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
    cat("Installing missing package: ggbeeswarm\n")
    install.packages("ggbeeswarm", repos = "http://cran.us.r-project.org")
  }
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
                    help = "Number of datasets to filter on, or 'all' (default: all)")

args <- parser$parse_args()

# ==============================================================================
# 1. Data Loading & Cleaning
# ==============================================================================

data <- read.csv(args$input, stringsAsFactors = FALSE)

# Clean whitespace
if("adjuster" %in% colnames(data)) data$adjuster <- trimws(data$adjuster)
if("classifier" %in% colnames(data)) data$classifier <- trimws(data$classifier)

# Filter Metrics & Exclusions
data <- data %>%
  filter(
    metric == "mcc", 
    !is.na(value), 
    !is.na(n_datasets),
    adjuster != "within_study_cv",
    classifier != "logistic",
    classifier != "rda"
  )

# Filter N-Datasets
n_datasets_label <- args$n_datasets
if (tolower(args$n_datasets) != "all") {
  n_val <- as.integer(args$n_datasets)
  data <- data %>% filter(n_datasets == n_val)
  n_datasets_label <- sprintf("%d-Study", n_val)
} else {
  n_datasets_label <- "All Studies"
}

# Filter Adjusters
if (!is.null(args$adjusters)) {
  adjusters_filter <- trimws(strsplit(args$adjusters, ",")[[1]])
  data <- data %>% filter(adjuster %in% adjusters_filter)
}

# ==============================================================================
# 2. Rank Calculation
# ==============================================================================

# Calculate rank WITHIN each (Classifier + Dataset + Study) group.
ranked_data <- data %>%
  group_by(classifier, n_datasets, test_study) %>%
  mutate(rank = rank(-value, ties.method = "average")) %>% # Rank 1 = Highest MCC
  ungroup()

# Aggregate to get Average Rank and Standard Error
avg_rank <- ranked_data %>%
  group_by(classifier, adjuster) %>%
  summarise(
    avg_rank = mean(rank, na.rm = TRUE),
    se_rank = sd(rank, na.rm = TRUE) / sqrt(n()),
    n_obs = n(),
    .groups = "drop"
  )

# ==============================================================================
# 3. Formatting & Logic
# ==============================================================================

# Classifier Labels
classifier_labels <- c(
  "rda" = "RDA", "logistic" = "Logistic", "elasticnet" = "ElasticNet", 
  "svm" = "SVM", "rf" = "Random Forest", "knn" = "KNN",
  "xgboost" = "XGBoost", "nnet" = "Neural Net", "shrinkageLDA" = "LDA"
)
avg_rank$classifier_label <- classifier_labels[avg_rank$classifier]
avg_rank$classifier_label[is.na(avg_rank$classifier_label)] <- avg_rank$classifier[is.na(avg_rank$classifier_label)]

# Adjuster Labels
if (exists("format_adjuster_label")) {
    avg_rank$adjuster_label <- sapply(avg_rank$adjuster, format_adjuster_label)
} else {
    avg_rank$adjuster_label <- avg_rank$adjuster
}

# Order Adjusters by Overall Performance (Best Average Rank on Left)
adjuster_order <- avg_rank %>%
  group_by(adjuster_label) %>%
  summarise(overall_avg = mean(avg_rank), .groups = "drop") %>%
  arrange(overall_avg)

avg_rank$adjuster_label <- factor(avg_rank$adjuster_label, levels = adjuster_order$adjuster_label)
# avg_rank$adjuster_num <- as.numeric(avg_rank$adjuster_label) # Not needed for ggbeeswarm

# Order Classifiers (Consistent Legend)
classifier_order <- sort(unique(avg_rank$classifier_label))
avg_rank$classifier_label <- factor(avg_rank$classifier_label, levels = classifier_order)

# ==============================================================================
# 4. Outlier Logic
# ==============================================================================

# Identify Outliers (Global Standardized Scale)
avg_rank <- avg_rank %>%
  group_by(adjuster_label) %>%
  mutate(median_rank_group = median(avg_rank, na.rm = TRUE)) %>%
  ungroup()

avg_rank$dist_from_median <- abs(avg_rank$avg_rank - avg_rank$median_rank_group)
global_dist_iqr <- IQR(avg_rank$dist_from_median, na.rm = TRUE)

avg_rank <- avg_rank %>%
  mutate(
    # Removed safety floor as requested (IQR assumed non-zero)
    standardized_dev = dist_from_median / global_dist_iqr, 
    dist_for_alpha = pmin(standardized_dev, 2.5)
  ) %>%
  # NEW: Sort by deviation to label top 6 outliers
  arrange(desc(dist_from_median)) %>%
  mutate(
    rank_outlier = row_number(),
    label_text = ifelse(rank_outlier <= 6, as.character(classifier_label), NA)
  )

max_rank <- ceiling(max(avg_rank$avg_rank + avg_rank$se_rank, na.rm = TRUE))

# ==============================================================================
# 5. Plotting
# ==============================================================================

# Define Shared Position Object
# 'method = "hex"' creates the honeycomb packing.
pos_hex <- position_beeswarm(
  method = "hex",
  cex = 1.5,
  groupOnX = TRUE
)

p <- ggplot(avg_rank, aes(x = adjuster_label, y = avg_rank)) +
  
  # A. The "Mass" Background (Boxplot)
  # Uses the factor x-axis directly
  geom_boxplot(
    width = 0.6, 
    outlier.shape = NA,
    alpha = 0.1, 
    color = "grey80", 
    fill = "grey90"
  ) +

  # B. Error Bars (Attached via shared position)
  # We MUST include y in aes() for position_beeswarm to work on errorbars
  geom_errorbar(
    aes(ymin = avg_rank - se_rank, 
        ymax = avg_rank + se_rank, 
        alpha = dist_for_alpha,
        group = classifier_label),
    position = pos_hex, # This syncs the jitter with the points
    width = 0,
    linewidth = 0.4, 
    color = "grey60" # Slightly darker to match increased point visibility
  ) +

  # C. The Points (Large, Attached via shared position)
  geom_point(
    aes(color = classifier_label, 
        group = classifier_label,
        alpha = dist_for_alpha),
    position = pos_hex,
    size = 5.0,
  ) +

  # D. Labels for Top 6 Outliers
  geom_text(
    aes(label = label_text, alpha = dist_for_alpha, group = classifier_label),
    position = pos_hex,
    # Move to top right:
    vjust = -0.8, # Slightly above (y is reversed)
    hjust = -0.3, # Start text to the right of the point
    size = 3.5,
    fontface = "bold",
    color = "black",
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  
  # Set baseline alpha to 0.6 to ensure points are definitely visible
  scale_alpha_continuous(range = c(0.6, 1.0), guide = "none") +

  scale_y_reverse(
    breaks = seq(1, max_rank, by = 1),
    limits = c(max_rank + 0.5, 0.5)
  ) +
  
  labs(
    title = sprintf("MCC Rank by Adjuster (%s)", n_datasets_label),
    subtitle = "Hex-packed by Rank. Labels show top 6 outliers. Opacity based on Global Outlier Scale.",
    x = NULL,
    y = "Average Rank (Lower is Better)",
    color = "Classifier"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face="bold", color="black"),
    panel.grid.major.x = element_line(color = "gray95"),
    panel.grid.major.y = element_line(color = "gray92", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face="bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40")
  )

ggsave(args$output, p, width = args$width, height = args$height, dpi = args$dpi, bg = "white")
cat("Saved plot to:", args$output, "\n")