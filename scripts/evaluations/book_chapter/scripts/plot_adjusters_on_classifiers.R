#!/usr/bin/env Rscript

# plot_adjusters_on_classifiers.R
# Script to create adjuster effectiveness on classifiers visualization
# Expected to be called from Snakemake workflow

SHARE_Y_AXIS <- TRUE

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
expected_cols <- c("adjuster", "classifier", "n_datasets", "test_study", "metric", "value")
missing_cols <- setdiff(expected_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("Available metrics:", paste(sort(unique(data$metric)), collapse = ", "), "\n")

# Filter to MCC metric for the main visualization
mxe_data <- data[data$metric == "mcc", ]
if (nrow(mxe_data) == 0) {
  stop("No MCC data found in input file")
}

cat("Filtered to", nrow(mxe_data), "MCC observations\n")

# Debug n_datasets values
cat("Unique n_datasets values:", paste(sort(unique(mxe_data$n_datasets)), collapse = ", "), "\n")
cat("n_datasets value counts:\n")
print(table(mxe_data$n_datasets, useNA = "always"))

# Check for and report missing n_datasets values
na_count <- sum(is.na(mxe_data$n_datasets))
if (na_count > 0) {
  cat("Warning: Found", na_count, "rows with missing n_datasets values. Removing them.\n")
  mxe_data <- mxe_data[!is.na(mxe_data$n_datasets), ]
  cat("After removing NA values:", nrow(mxe_data), "observations remain\n")
}

# Create better labels and groupings
mxe_data$classifier_label <- factor(mxe_data$classifier,
  levels = c("logistic", "elasticnet", "svm", "rf", "knn", "xgboost", "nnet", "rvc"),
  labels = c("Logistic", "ElasticNet", "SVM", "Random Forest", "KNN", "XGBoost", "Neural Net", "RVC"))

# Dynamically determine adjuster levels from the data
unique_adjusters <- sort(unique(mxe_data$adjuster))
cat("Unique adjusters in data:", paste(unique_adjusters, collapse = ", "), "\n")

# Create labels with proper capitalization
adjuster_labels <- sapply(unique_adjusters, function(x) {
  if (x == "unadjusted") return("Unadjusted")
  if (x == "combat") return("ComBat")
  if (x == "combat_sup") return("ComBat-Sup")
  if (x == "mnn") return("MNN")
  return(tools::toTitleCase(x))
})

mxe_data$adjuster_label <- factor(mxe_data$adjuster,
  levels = unique_adjusters,
  labels = adjuster_labels)

mxe_data$adjuster_type <- "Batch Correction"
mxe_data$adjuster_type[mxe_data$adjuster == "unadjusted"] <- "Original Data"

mxe_data$dataset_label <- paste(mxe_data$n_datasets, "studies")

# Convert to factor with explicit levels to avoid NA levels
mxe_data$dataset_label <- factor(mxe_data$dataset_label, 
                                levels = c("3 studies", "4 studies", "5 studies", "6 studies"))

# Debug dataset_label values
cat("Unique dataset_label values:", paste(sort(unique(mxe_data$dataset_label)), collapse = ", "), "\n")
cat("dataset_label value counts:\n")
print(table(mxe_data$dataset_label, useNA = "always"))

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

# Calculate frequency of best method for annotations
freq_data <- mxe_data %>%
  group_by(classifier_label, dataset_label, test_study) %>%
  summarise(
    best_adjuster = adjuster_label[which.max(value)],
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

# Debug final sumstats
cat("Unique dataset_label values in sumstats:", paste(sort(unique(sumstats$dataset_label)), collapse = ", "), "\n")
cat("Any NA values in sumstats dataset_label:", any(is.na(sumstats$dataset_label)), "\n")

# Create color scheme (following original pattern)
type_colors <- c("Original Data" = "#999999", "Batch Correction" = "#E69F00")

# Determine which classifiers actually have data (excluding NA classifiers)
classifiers_with_data <- sumstats %>%
  filter(!is.na(classifier_label)) %>%
  group_by(classifier_label) %>%
  summarise(has_data = n() > 0, .groups = "drop") %>%
  filter(has_data) %>%
  pull(classifier_label)

cat("Classifiers with data:", paste(classifiers_with_data, collapse = ", "), "\n")
cat("Y-axis sharing:", ifelse(SHARE_Y_AXIS, "enabled", "disabled"), "\n")

# Calculate global y-axis limits if sharing is enabled
if (SHARE_Y_AXIS) {
  global_y_min <- 0 #min(mxe_data$value, na.rm = TRUE)
  global_y_max <- max(mxe_data$value, na.rm = TRUE)
  global_y_range <- global_y_max - global_y_min
  # Add some padding for annotations
  global_y_limits <- c(global_y_min - 0.05 * global_y_range, 
                       global_y_max + 0.15 * global_y_range)
  cat("Using shared y-axis limits:", round(global_y_limits[1], 3), "to", round(global_y_limits[2], 3), "\n")
}

# Create individual plots for each classifier (grouped by classifier first)
plot_list <- list()

for (classifier in classifiers_with_data) {
  plot_data <- sumstats[sumstats$classifier_label == classifier & !is.na(sumstats$classifier_label), ]
  
  # Skip if no data after filtering
  if (nrow(plot_data) == 0) {
    cat("Skipping", classifier, "- no data after filtering\n")
    next
  }
  
  # Ensure consistent factor levels across all plots (only the valid ones)
  plot_data$dataset_label <- factor(plot_data$dataset_label, 
                                   levels = c("3 studies", "4 studies", "5 studies", "6 studies"))
  plot_data$adjuster_label <- factor(plot_data$adjuster_label,
                                    levels = adjuster_labels)
  
  # Debug plot_data for this classifier
  cat("Classifier:", classifier, "\n")
  cat("  Unique dataset_label values in plot_data:", paste(sort(unique(plot_data$dataset_label)), collapse = ", "), "\n")
  cat("  Factor levels of dataset_label:", paste(levels(plot_data$dataset_label), collapse = ", "), "\n")
  
  # Get raw data for this classifier for boxplots
  raw_data <- mxe_data[mxe_data$classifier_label == classifier & !is.na(mxe_data$classifier_label), ]
  
  # Calculate annotation position at the top of the plot
  if (SHARE_Y_AXIS) {
    annotation_y <- global_y_limits[2] * 0.98  # Position at 98% of max y-axis
  } else {
    local_y_max <- max(raw_data$value, na.rm = TRUE)
    local_y_min <- min(raw_data$value, na.rm = TRUE)
    local_y_range <- local_y_max - local_y_min
    annotation_y <- local_y_max + 0.12 * local_y_range  # Position at top with padding
  }
  
  p <- ggplot(raw_data, aes(x = adjuster_label, y = value, fill = adjuster_type)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 1, alpha = 0.7) +
    geom_text(data = plot_data, aes(x = adjuster_label, y = annotation_y, 
                                   label = annot, fill = NULL), 
              color = "black", size = 3, vjust = 1) +
    facet_wrap(~ dataset_label, scales = "fixed", ncol = 4) +
    {if (SHARE_Y_AXIS) {
      scale_y_continuous(limits = global_y_limits, expand = expansion(mult = c(0, 0)))
    } else {
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
    }} +
    scale_fill_manual(values = type_colors) +
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
      y = "Matthews Correlation Coefficient",
      title = classifier
    )
  
  plot_list[[as.character(classifier)]] <- p
}

# Debug plot_list structure
cat("Plot list names:", paste(names(plot_list), collapse = ", "), "\n")
cat("Plot list length:", length(plot_list), "\n")

# Create legend from the first available plot
first_plot_name <- names(plot_list)[1]
legend <- get_legend(
  plot_list[[first_plot_name]] + 
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(color = guide_legend(title = NULL))
)

# Try a simpler approach with grid.arrange
library(gridExtra)
plot_vector <- unname(plot_list)

# Arrange plots manually based on how many we have
# Fallback for other numbers of plots
  final_plot <- grid.arrange(
    grobs = c(plot_vector, list(legend)),
    ncol = 2
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