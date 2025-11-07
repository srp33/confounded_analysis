#!/usr/bin/env Rscript

# plot_batch_on_classifiers.R
# Script to create batch effects on classifiers visualization
# Expected to be called from Snakemake workflow

# Suppress warnings and messages for cleaner output
options(warn = -1)
suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(RColorBrewer)
  library(gridExtra)
  library(ggpubr)
  library(scales)
})

# Define command line arguments
parser <- ArgumentParser(description = "Create batch effects on classifiers visualization")

parser$add_argument("-i", "--input", type = "character", required = TRUE,
                   help = "Input CSV file with batch effects on classifiers data")

parser$add_argument("-o", "--output", type = "character", default = "batch_on_classifiers.png",
                   help = "Output PNG file path (default: %(default)s)")

parser$add_argument("--width", type = "double", default = 16,
                   help = "Plot width in inches (default: %(default)s)")

parser$add_argument("--height", type = "double", default = 12,
                   help = "Plot height in inches (default: %(default)s)")

parser$add_argument("--dpi", type = "integer", default = 300,
                   help = "Plot resolution in DPI (default: %(default)s)")

# Parse arguments and input file
args <- parser$parse_args()

cat("Reading input data from:", args$input, "\n")
data <- read.csv(args$input, stringsAsFactors = FALSE)

cat("Data dimensions:", nrow(data), "rows,", ncol(data), "columns\n")
cat("Column names:", paste(colnames(data), collapse = ", "), "\n")

# Validate expected columns for batch effects data
expected_cols <- c("classifier", "mean", "variance", "seed", "metric", "value")
missing_cols <- setdiff(expected_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

#MXE = mean cross entropy
#MCC = Matthews Correlation Coefficient

# Filter to balanced accuracy for the main visualization
metric = "balanced_acc"
metric_name = "Balanced Accuracy"
metric_data <- data[data$metric == metric, ]
if (nrow(metric_data) == 0) {
  cat("Warning: No data found for metric 'balanced_acc', trying 'auc' instead\n")
  metric = "auc"
  metric_name = "AUC"
  metric_data <- data[data$metric == metric, ]
  if (nrow(metric_data) == 0) {
    cat("Available metrics:", paste(unique(data$metric), collapse = ", "), "\n")
    stop("No data found for metric:", metric)
  }
}

cat("Filtered to", nrow(metric_data), " observations\n")

# Create better labels and groupings
metric_data$classifier_label <- factor(metric_data$classifier,
  levels = c("logistic", "elasticnet", "svm", "rf", "knn", "xgboost", "nn", "lightgbm"),
  labels = c("Logistic", "ElasticNet", "SVM", "Random Forest", "KNN", "XGBoost", "Neural Net", "LightGBM"))

# Create batch effect intensity labels
metric_data$batch_intensity <- paste("Mean:", metric_data$mean, "Var:", metric_data$variance)
metric_data$batch_label <- factor(metric_data$batch_intensity,
  levels = unique(metric_data$batch_intensity[order(metric_data$mean, metric_data$variance)]))

# Create variance grouping for color coding
metric_data$variance_group <- factor(metric_data$variance,
  levels = sort(unique(metric_data$variance)),
  labels = paste("Variance =", sort(unique(metric_data$variance))))

cat("Creating figure\n")

# Calculate summary statistics for each combination
sumstats <- metric_data %>%
  group_by(classifier_label, mean, variance, variance_group, batch_label) %>%
  summarise(
    Avg = mean(value),
    Up = quantile(value, 0.975),
    Down = quantile(value, 0.025),
    .groups = "drop"
  )

# Create color palette for variance levels
n_variances <- length(unique(metric_data$variance))
variance_colors <- RColorBrewer::brewer.pal(max(3, min(n_variances, 9)), "Set1")[1:n_variances]
names(variance_colors) <- paste("Variance =", sort(unique(metric_data$variance)))

# Create the main plot
p_main <- ggplot(sumstats, aes(x = mean, y = Avg, color = variance_group)) +
  geom_errorbar(aes(ymin = Down, ymax = Up), width = 0.05, alpha = 0.7) +
  geom_line(aes(group = variance_group), size = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ classifier_label, scales = "free_y", ncol = 4) +
  scale_color_manual(values = variance_colors) +
  scale_x_continuous(breaks = sort(unique(metric_data$mean))) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey95"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey40")
  ) +
  labs(
    x = "Batch Effect Mean",
    y =  metric_name,
    title = "Impact of Batch Effects on Classifier Performance",
    subtitle = "Higher values indicate worse performance. Error bars show 95% confidence intervals."
  ) +
  guides(color = guide_legend(title = NULL, nrow = 1))

# Create a summary heatmap showing relative performance degradation
baseline_data <- sumstats %>%
  filter(mean == min(mean) & variance == min(variance)) %>%
  select(classifier_label, baseline_avg = Avg)

degradation_data <- sumstats %>%
  left_join(baseline_data, by = "classifier_label") %>%
  mutate(
    degradation = (Avg - baseline_avg) / baseline_avg * 100,
    degradation_label = sprintf("%.1f%%", degradation)
  )

p_heatmap <- ggplot(degradation_data, aes(x = factor(mean), y = factor(variance), fill = degradation)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = degradation_label), color = "white", size = 3, fontface = "bold") +
  facet_wrap(~ classifier_label, ncol = 4) +
  scale_fill_gradient2(
    low = "darkblue", mid = "white", high = "darkred",
    midpoint = 0, name = "Performance\nDegradation (%)"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "right",
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Batch Effect Mean",
    y = "Batch Effect Variance",
    title = "Relative Performance Degradation Heatmap"
  )

# Combine plots
combined_plot <- ggarrange(
  p_main, p_heatmap,
  ncol = 1, nrow = 2,
  heights = c(2, 1.2),
  labels = c("A", "B"),
  font.label = list(size = 16, face = "bold")
)

# Add overall title and caption
final_plot <- annotate_figure(
  combined_plot,
  top = text_grob("Batch Effects Impact on Machine Learning Classifiers", 
                  face = "bold", size = 18),
  bottom = text_grob("Panel A: Performance trends across batch effect intensities. Panel B: Relative degradation from baseline (no batch effects).",
                     size = 10, color = "grey40")
)

# Save the plot
cat("Saving plot to:", args$output, "\n")
ggsave(
  filename = args$output,
  plot = final_plot,
  width = args$width,
  height = args$height,
  dpi = args$dpi,
  units = "in"
)

cat("Plot saved successfully!\n")
cat("Output file:", args$output, "\n")
cat("Dimensions:", args$width, "x", args$height, "inches at", args$dpi, "DPI\n")