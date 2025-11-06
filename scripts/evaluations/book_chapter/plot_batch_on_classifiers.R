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
})

# Define command line arguments
parser <- ArgumentParser(description = "Create batch effects on classifiers visualization")

parser$add_argument("-i", "--input", 
                   type = "character", 
                   required = TRUE,
                   help = "Input CSV file with batch effects on classifiers data")

parser$add_argument("-o", "--output", 
                   type = "character", 
                   default = "batch_on_classifiers.png",
                   help = "Output PNG file path (default: %(default)s)")

parser$add_argument("--width", 
                   type = "double", 
                   default = 12,
                   help = "Plot width in inches (default: %(default)s)")

parser$add_argument("--height", 
                   type = "double", 
                   default = 8,
                   help = "Plot height in inches (default: %(default)s)")

parser$add_argument("--dpi", 
                   type = "integer", 
                   default = 300,
                   help = "Plot resolution in DPI (default: %(default)s)")

# Parse arguments
opt <- parser$parse_args()

# Read and validate input data
cat("Reading input data from:", opt$input, "\n")
data <- read.csv(opt$input, stringsAsFactors = FALSE)

cat("Data dimensions:", nrow(data), "rows,", ncol(data), "columns\n")
cat("Column names:", paste(colnames(data), collapse = ", "), "\n")

# Validate expected columns for batch effects data
expected_cols <- c("classifier", "mean", "variance", "seed", "method", "metric", "value")
missing_cols <- setdiff(expected_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Filter to AUC metric for the main visualization (following original pattern)
auc_data <- data[data$metric == "auc", ]
if (nrow(auc_data) == 0) {
  stop("No AUC data found in input file")
}

cat("Filtered to", nrow(auc_data), "AUC observations\n")

# Create batch effect condition labels (following original pattern)
auc_data$batch_condition <- paste0("m", auc_data$mean, "_v", auc_data$variance)

# Create readable labels for batch conditions
auc_data$batch_label <- factor(auc_data$batch_condition,
  levels = c("m5_v1", "m5_v3", "m5_v5"),
  labels = c("Mean difference 5\nVariance fold change 1",
             "Mean difference 5\nVariance fold change 3", 
             "Mean difference 5\nVariance fold change 5"))

# Create readable method labels
auc_data$method_label <- factor(auc_data$method,
  levels = c("unadjusted", "combat", "mnn"),
  labels = c("No adjustment", "ComBat", "MNN"))

# Create method type groupings (following original pattern)
auc_data$method_type <- "Batch Correction"
auc_data$method_type[auc_data$method == "unadjusted"] <- "With batch effect,\nno adjustment"

# Create readable classifier labels
auc_data$classifier_label <- factor(auc_data$classifier,
  levels = c("logistic", "elasticnet", "svm", "rf", "knn", "xgboost", "nn", "lightgbm"),
  labels = c("Logistic", "ElasticNet", "SVM", "Random Forest", "KNN", "XGBoost", "Neural Net", "LightGBM"))

# Calculate baseline performance (no batch effects) for reference line
# This would typically come from a separate analysis, but we'll estimate from the data
baseline_auc <- median(auc_data$value[auc_data$method == "unadjusted"])

cat("Creating visualization...\n")

# Create the main plot (following Figure 1 pattern from original code)
p <- ggplot(auc_data, aes(x = classifier_label, y = value, fill = method_type)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = baseline_auc,
                 linetype = "Baseline AUC\n(estimated)"), 
             color = "red", size = 0.8) +
  scale_linetype_manual(name = "Reference", values = 2,
                        guide = guide_legend(override.aes = list(color = c("red")))) +
  scale_fill_manual(values = c("With batch effect,\nno adjustment" = "#999999", 
                               "Batch Correction" = "#E69F00")) +
  facet_wrap(~ batch_label, ncol = 3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    plot.margin = margin(0.1, 0.1, 0.1, 0.4, "cm"),
    strip.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    y = "AUC",
    title = "Classifier Response to Batch Effects",
    subtitle = "Higher AUC indicates better performance"
  )





# Save the plot
cat("Saving plot to:", opt$output, "\n")
ggsave(
  filename = opt$output,
  plot = p,
  width = opt$width,
  height = opt$height,
  dpi = opt$dpi,
  units = "in"
)

cat("Plot saved successfully!\n")
cat("Output file:", opt$output, "\n")
cat("Dimensions:", opt$width, "x", opt$height, "inches at", opt$dpi, "DPI\n")