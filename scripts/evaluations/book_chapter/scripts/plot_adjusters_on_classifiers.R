#!/usr/bin/env Rscript

# plot_adjusters_on_classifiers.R
# Script to create adjuster effectiveness on classifiers visualization
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
parser <- ArgumentParser(description = "Create adjuster effectiveness on classifiers visualization")

parser$add_argument("-i", "--input", type = "character", required = TRUE,
                   help = "Input CSV file with adjusters on classifiers data")

parser$add_argument("-o", "--output", type = "character", default = "adjusters_on_classifiers.png",
                   help = "Output PNG file path (default: %(default)s)")

parser$add_argument("--width", type = "double", default = 14,
                   help = "Plot width in inches (default: %(default)s)")

parser$add_argument("--height", type = "double", default = 16,
                   help = "Plot height in inches (default: %(default)s)")

parser$add_argument("--dpi", type = "integer", default = 300,
                   help = "Plot resolution in DPI (default: %(default)s)")

# Parse arguments and input file
args <- parser$parse_args()

cat("Reading input data from:", args$input, "\n")
data <- read.csv(args$input, stringsAsFactors = FALSE)

cat("Data dimensions:", nrow(data), "rows,", ncol(data), "columns\n")
cat("Column names:", paste(colnames(data), collapse = ", "), "\n")

# Validate expected columns for adjusters data
expected_cols <- c("adjuster", "classifier", "n_datasets", "seed", "metric", "value")
missing_cols <- setdiff(expected_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Filter to MXE metric (cross entropy) for the main visualization
mxe_data <- data[data$metric == "mxe", ]
if (nrow(mxe_data) == 0) {
  stop("No MXE data found in input file")
}

cat("Filtered to", nrow(mxe_data), "MXE observations\n")

# Create better labels and groupings
mxe_data$classifier_label <- factor(mxe_data$classifier,
  levels = c("logistic", "elasticnet", "svm", "rf", "knn", "xgboost", "nn", "lightgbm"),
  labels = c("Logistic", "ElasticNet", "SVM", "Random Forest", "KNN", "XGBoost", "Neural Net", "LightGBM"))

mxe_data$adjuster_label <- factor(mxe_data$adjuster,
  levels = c("unadjusted", "combat", "mnn"),
  labels = c("Unadjusted", "ComBat", "MNN"))

mxe_data$adjuster_type <- "Batch Correction"
mxe_data$adjuster_type[mxe_data$adjuster == "unadjusted"] <- "Original Data"

mxe_data$dataset_label <- paste(mxe_data$n_datasets, "studies")

cat("Creating figure\n")

# Calculate summary statistics for each combination (following original pattern)
sumstats <- mxe_data %>%
  group_by(adjuster_label, classifier_label, dataset_label, adjuster_type) %>%
  summarise(
    Avg = mean(value),
    Up = quantile(value, 0.975),
    Down = quantile(value, 0.025),
    .groups = "drop"
  )

# Calculate frequency of best method for annotations (simplified version)
freq_data <- mxe_data %>%
  group_by(classifier_label, dataset_label, seed) %>%
  summarise(
    best_adjuster = adjuster_label[which.min(value)],
    .groups = "drop"
  ) %>%
  group_by(classifier_label, dataset_label, best_adjuster) %>%
  summarise(
    freq = n(),
    .groups = "drop"
  ) %>%
  group_by(classifier_label, dataset_label) %>%
  mutate(
    total = sum(freq),
    pct = freq / total
  ) %>%
  ungroup()

# Add frequency annotations to summary stats
sumstats <- sumstats %>%
  left_join(
    freq_data %>% 
      select(classifier_label, dataset_label, best_adjuster, pct) %>%
      rename(adjuster_label = best_adjuster),
    by = c("classifier_label", "dataset_label", "adjuster_label")
  ) %>%
  mutate(
    annot = ifelse(is.na(pct), "", percent(pct, accuracy = 1))
  )

# Create color scheme (following original pattern)
type_colors <- c("Original Data" = "#999999", "Batch Correction" = "#E69F00")

# Create individual plots for each classifier (grouped by classifier first)
plot_list <- list()
classifiers <- levels(mxe_data$classifier_label)

for (classifier in classifiers) {
  plot_data <- sumstats[sumstats$classifier_label == classifier, ]
  
  p <- ggplot(plot_data, aes(x = adjuster_label, y = Avg, color = adjuster_type)) +
    geom_errorbar(aes(ymin = Down, ymax = Up), width = 0.1) +
    geom_text(aes(label = annot, y = Up), color = "black", size = 3, vjust = -0.5) +
    geom_line(aes(group = 1), color = "grey", size = 0.5) +
    geom_point(size = 2) +
    facet_wrap(~ dataset_label, scales = "fixed", ncol = 4) +
    scale_color_manual(values = type_colors) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
      axis.title.y = element_text(size = 10),
      legend.title = element_blank(),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "grey90", size = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 9),
      plot.title = element_text(size = 12, hjust = 0.5)
    ) +
    labs(
      y = "Mean Cross-Entropy Loss",
      title = paste("Adjuster Effectiveness -", classifier)
    )
  
  plot_list[[classifier]] <- p
}

# Create legend from one of the plots
legend <- get_legend(
  plot_list[["Logistic"]] + 
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(color = guide_legend(title = NULL))
)

# Arrange plots in a grid (8 classifiers arranged in 4x2 grid)
combined_plots <- ggarrange(
  plot_list[["Logistic"]], plot_list[["ElasticNet"]],
  plot_list[["SVM"]], plot_list[["Random Forest"]],
  plot_list[["KNN"]], plot_list[["XGBoost"]],
  plot_list[["Neural Net"]], plot_list[["LightGBM"]],
  ncol = 2, nrow = 4,
  common.legend = TRUE,
  legend = "bottom"
)

# Add overall title
final_plot <- annotate_figure(
  combined_plots,
  top = text_grob("Batch Correction Method Effectiveness Across Classifiers", 
                  face = "bold", size = 16),
  bottom = text_grob("Lower values indicate better performance. Percentages show frequency of being the best method.",
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