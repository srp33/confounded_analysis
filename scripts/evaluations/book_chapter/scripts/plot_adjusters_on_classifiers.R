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

parser$add_argument("--width", type = "double", default = 20,
                   help = "Plot width in inches (default: %(default)s)")

parser$add_argument("--height", type = "double", default = 16,
                   help = "Plot height in inches (default: %(default)s)")

parser$add_argument("--dpi", type = "integer", default = 300,
                   help = "Plot resolution in DPI (default: %(default)s)")

parser$add_argument("--adjusters", type = "character", default = NULL,
                   help = "Comma-separated list of adjusters to include (default: all)")

# Parse arguments and input file
args <- parser$parse_args()

cat("Reading input data from:", args$input, "\n")
data <- read.csv(args$input, stringsAsFactors = FALSE)

# Filter to specified adjusters if provided
if (!is.null(args$adjusters)) {
  adjusters_to_include <- trimws(strsplit(args$adjusters, ",")[[1]])
  cat("Filtering to adjusters:", paste(adjusters_to_include, collapse = ", "), "\n")
  data <- data[data$adjuster %in% adjusters_to_include, ]
  cat("After filtering:", nrow(data), "rows remain\n")
}

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

# Calculate overall performance for each adjuster to determine ordering
adjuster_performance <- mxe_data %>%
  group_by(adjuster) %>%
  summarise(overall_mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall_mean))

# Order adjusters by performance (best on the left)
unique_adjusters <- adjuster_performance$adjuster
cat("Adjusters ordered by overall performance (best first):\n")
for (i in 1:nrow(adjuster_performance)) {
  cat(sprintf("  %d. %s (mean MCC: %.4f)\n", 
              i, adjuster_performance$adjuster[i], adjuster_performance$overall_mean[i]))
}

# Calculate best adjuster per classifier
cat("\nBest adjuster per classifier:\n")
classifier_best <- mxe_data %>%
  group_by(classifier, adjuster) %>%
  summarise(mean_mcc = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(classifier) %>%
  arrange(desc(mean_mcc)) %>%
  slice(1) %>%
  ungroup()

for (i in 1:nrow(classifier_best)) {
  cat(sprintf("  %s: %s (mean MCC: %.4f)\n", 
              classifier_best$classifier[i], 
              classifier_best$adjuster[i], 
              classifier_best$mean_mcc[i]))
}

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
  # Get raw data for this classifier
  raw_data <- mxe_data[mxe_data$classifier_label == classifier & !is.na(mxe_data$classifier_label), ]
  
  # Skip if no data after filtering
  if (nrow(raw_data) == 0) {
    cat("Skipping", classifier, "- no data after filtering\n")
    next
  }
  
  # Calculate classifier-specific adjuster ordering (best on left)
  classifier_adjuster_order <- raw_data %>%
    group_by(adjuster, adjuster_label) %>%
    summarise(mean_mcc = mean(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_mcc))
  
  classifier_specific_levels <- classifier_adjuster_order$adjuster_label
  
  cat("Classifier:", classifier, "\n")
  cat("  Adjuster ordering for this classifier:\n")
  for (i in 1:nrow(classifier_adjuster_order)) {
    cat(sprintf("    %d. %s (mean MCC: %.4f)\n", 
                i, classifier_adjuster_order$adjuster_label[i], 
                classifier_adjuster_order$mean_mcc[i]))
  }
  
  # Apply classifier-specific ordering to raw data
  raw_data$adjuster_label <- factor(raw_data$adjuster_label,
                                    levels = classifier_specific_levels)
  
  # Ensure dataset_label is properly factored
  raw_data$dataset_label <- factor(raw_data$dataset_label, 
                                   levels = c("3 studies", "4 studies", "5 studies", "6 studies"))
  
  # Get plot_data with same ordering
  plot_data <- sumstats[sumstats$classifier_label == classifier & !is.na(sumstats$classifier_label), ]
  plot_data$adjuster_label <- factor(plot_data$adjuster_label,
                                    levels = classifier_specific_levels)
  plot_data$dataset_label <- factor(plot_data$dataset_label, 
                                   levels = c("3 studies", "4 studies", "5 studies", "6 studies"))
  
  # Calculate annotation position at the top of the plot
  if (SHARE_Y_AXIS) {
    annotation_y <- global_y_limits[2] * 0.98  # Position at 98% of max y-axis
  } else {
    local_y_max <- max(raw_data$value, na.rm = TRUE)
    local_y_min <- min(raw_data$value, na.rm = TRUE)
    local_y_range <- local_y_max - local_y_min
    annotation_y <- local_y_max + 0.12 * local_y_range  # Position at top with padding
  }
  
  # Calculate mean for each adjuster to show as horizontal line
  mean_data <- raw_data %>%
    group_by(adjuster_label, dataset_label, adjuster_type) %>%
    summarise(mean_value = mean(value), .groups = "drop")
  
  p <- ggplot(raw_data, aes(x = adjuster_label, y = value)) +
    # Add mean line for each adjuster
    geom_segment(data = mean_data, 
                aes(x = as.numeric(adjuster_label) - 0.3, 
                    xend = as.numeric(adjuster_label) + 0.3,
                    y = mean_value, yend = mean_value),
                color = "gray40", linewidth = 0.8) +
    # Add individual points colored by test study
    geom_point(aes(color = test_study, shape = test_study), 
              size = 3, alpha = 0.5) +
    geom_text(data = plot_data, aes(x = adjuster_label, y = annotation_y, 
                                   label = annot), 
              color = "black", size = 2.0, vjust = 1) +
    facet_wrap(~ dataset_label, scales = "fixed", ncol = 4) +
    scale_x_discrete() +
    {if (SHARE_Y_AXIS) {
      scale_y_continuous(limits = global_y_limits, expand = expansion(mult = c(0, 0)))
    } else {
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
    }} +
    scale_color_brewer(palette = "Set2", name = "Test Study") +
    scale_shape_manual(values = c(15, 16, 17, 25, 18, 19), name = "Test Study") +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 9, face = "bold"),
      legend.position = "right",
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

# Arrange plots in grid (legend is now part of each plot)
library(gridExtra)
plot_vector <- unname(plot_list)

final_plot <- grid.arrange(
  grobs = plot_vector,
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

# ====================================================================
# CREATE PERFORMANCE DIFFERENCE DISTRIBUTION PLOT
# ====================================================================

cat("\nCreating performance difference distribution plot...\n")

# Get the top adjuster (first in the ordered list)
top_adjuster <- unique_adjusters[1]
top_adjuster_label <- adjuster_labels[1]
cat("Top adjuster:", top_adjuster, "(", top_adjuster_label, ")\n")

# Calculate differences for each observation
# For each combination of classifier, n_datasets, and test_study,
# compute the difference between top adjuster and each other adjuster
difference_data <- mxe_data %>%
  # Create a unique identifier for each experimental condition
  mutate(condition_id = paste(classifier, n_datasets, test_study, sep = "_")) %>%
  # Get the top adjuster's value for each condition
  group_by(condition_id) %>%
  mutate(top_value = value[adjuster == top_adjuster]) %>%
  ungroup() %>%
  # Calculate difference (top - current)
  mutate(difference = top_value - value) %>%
  # Filter out the top adjuster itself (difference = 0)
  filter(adjuster != top_adjuster)

# Add labels
difference_data$adjuster_label <- factor(difference_data$adjuster,
  levels = unique_adjusters[-1],  # Exclude top adjuster
  labels = adjuster_labels[-1])

cat("Calculated", nrow(difference_data), "pairwise differences\n")

# ====================================================================
# STATISTICAL TESTING
# ====================================================================

cat("\nPerforming statistical tests...\n")
cat("Adjusters compared to", top_adjuster_label, ":\n")

# Perform paired t-test and Wilcoxon signed-rank test for each adjuster
stat_results <- data.frame()

for (adj in unique(difference_data$adjuster_label)) {
  adj_data <- difference_data[difference_data$adjuster_label == adj, ]
  n_obs <- nrow(adj_data)
  mean_diff <- mean(adj_data$difference, na.rm = TRUE)
  median_diff <- median(adj_data$difference, na.rm = TRUE)
  sd_diff <- sd(adj_data$difference, na.rm = TRUE)
  
  # One-sample t-test (testing if mean difference is significantly > 0)
  t_test <- t.test(adj_data$difference, mu = 0, alternative = "greater")
  
  # One-sample Wilcoxon signed-rank test (non-parametric alternative)
  wilcox_test <- wilcox.test(adj_data$difference, mu = 0, alternative = "greater")
  
  # Format p-values
  format_pval <- function(p) {
    if (p < 0.001) return("p < 0.001")
    if (p < 0.01) return(sprintf("p = %.3f", p))
    return(sprintf("p = %.2f", p))
  }
  
  # Significance stars
  get_stars <- function(p) {
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    return("ns")
  }
  
  cat(sprintf("  %s: n=%d, mean_diff=%.4f, median_diff=%.4f, sd=%.4f\n", 
              adj, n_obs, mean_diff, median_diff, sd_diff))
  cat(sprintf("    t-test: t=%.3f, %s %s\n", 
              t_test$statistic, format_pval(t_test$p.value), get_stars(t_test$p.value)))
  cat(sprintf("    Wilcoxon: %s %s\n", 
              format_pval(wilcox_test$p.value), get_stars(wilcox_test$p.value)))
  
  # Store results
  stat_results <- rbind(stat_results, data.frame(
    adjuster_label = adj,
    n = n_obs,
    mean_diff = mean_diff,
    median_diff = median_diff,
    sd_diff = sd_diff,
    t_statistic = t_test$statistic,
    t_pvalue = t_test$p.value,
    wilcox_pvalue = wilcox_test$p.value,
    significance = get_stars(t_test$p.value),
    stringsAsFactors = FALSE
  ))
}

# Add statistical annotations to the plot data
stat_results$pval_label <- sapply(stat_results$t_pvalue, function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
})

# Create violin/box plot showing distribution of differences
# Calculate y position for annotations (above the highest point)
y_max <- max(difference_data$difference, na.rm = TRUE)
y_min <- min(difference_data$difference, na.rm = TRUE)
y_range <- y_max - y_min
annotation_y <- y_max + 0.05 * y_range

diff_plot <- ggplot(difference_data, aes(x = adjuster_label, y = difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_violin(fill = "#E69F00", alpha = 0.3, color = NA) +
  geom_boxplot(width = 0.2, outlier.size = 1, outlier.alpha = 0.5, fill = "white") +
  geom_jitter(aes(color = classifier_label), width = 0.15, height = 0, 
              size = 1.5, alpha = 0.4) +
  # Add p-value annotations
  geom_text(data = stat_results, 
            aes(x = adjuster_label, y = annotation_y, label = pval_label),
            size = 5, fontface = "bold", vjust = 0) +
  scale_color_brewer(palette = "Set3", name = "Classifier") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    plot.caption = element_text(size = 9, hjust = 0, color = "gray40")
  ) +
  labs(
    x = "Adjuster",
    y = sprintf("Performance Difference\n(%s MCC - Adjuster MCC)", top_adjuster_label),
    title = "Distribution of Performance Differences Relative to Top Adjuster",
    subtitle = sprintf("Positive values indicate %s performs better", top_adjuster_label),
    caption = "Significance from one-sample t-test (H₁: difference > 0): *** p<0.001, ** p<0.01, * p<0.05, ns = not significant"
  )

# Save the difference plot
diff_output <- sub("\\.png$", "_differences.png", args$output)
cat("Saving difference plot to:", diff_output, "\n")
ggsave(
  filename = diff_output,
  plot = diff_plot,
  width = 12,
  height = 8,
  dpi = args$dpi,
  units = "in"
)

cat("Difference plot saved successfully!\n")
cat("Output file:", diff_output, "\n")

# Save statistical results to CSV
stat_output <- sub("\\.png$", "_statistics.csv", args$output)
write.csv(stat_results, stat_output, row.names = FALSE)
cat("\nStatistical results saved to:", stat_output, "\n")
cat("\nSummary of statistical tests:\n")
print(stat_results[, c("adjuster_label", "n", "mean_diff", "t_pvalue", "wilcox_pvalue", "significance")])

# Save classifier-specific best adjusters to CSV
classifier_best_output <- sub("\\.png$", "_best_per_classifier.csv", args$output)
write.csv(classifier_best, classifier_best_output, row.names = FALSE)
cat("\nBest adjuster per classifier saved to:", classifier_best_output, "\n")

# ====================================================================
# CREATE RELATIVE PERFORMANCE PLOT (SIMILAR TO MAIN FIGURE)
# ====================================================================

cat("\nCreating relative performance plot (mirroring main figure structure)...\n")

# Prepare data with differences from top adjuster for each condition
# Start fresh from mxe_data but keep the raw adjuster names
relative_data <- mxe_data %>%
  mutate(condition_id = paste(classifier, n_datasets, test_study, sep = "_")) %>%
  group_by(condition_id) %>%
  mutate(top_value = value[adjuster == top_adjuster]) %>%
  ungroup() %>%
  mutate(relative_value = value - top_value) %>%
  # Keep all adjusters including the top one (which will be at 0)
  select(-top_value, -condition_id, -adjuster_label)  # Remove old adjuster_label

# We'll set adjuster_label per classifier in the loop below

# Calculate mean relative performance for each combination
relative_mean_data <- relative_data %>%
  group_by(adjuster_label, classifier_label, dataset_label) %>%
  summarise(mean_relative = mean(relative_value, na.rm = TRUE), .groups = "drop")

# Create individual plots for each classifier
relative_plot_list <- list()

for (classifier in classifiers_with_data) {
  raw_data <- relative_data[relative_data$classifier_label == classifier & !is.na(relative_data$classifier_label), ]
  
  if (nrow(raw_data) == 0) {
    cat("Skipping", classifier, "- no data after filtering\n")
    next
  }
  
  # Calculate classifier-specific adjuster ordering based on original MCC values
  classifier_adjuster_order <- mxe_data %>%
    filter(classifier_label == classifier & !is.na(classifier_label)) %>%
    group_by(adjuster) %>%
    summarise(mean_mcc = mean(value, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_mcc))
  
  # Create labels for this ordering
  classifier_specific_labels <- sapply(classifier_adjuster_order$adjuster, function(x) {
    if (x == "unadjusted") return("Unadjusted")
    if (x == "combat") return("ComBat")
    if (x == "combat_sup") return("ComBat-Sup")
    if (x == "mnn") return("MNN")
    return(tools::toTitleCase(x))
  })
  
  # Apply classifier-specific ordering to raw_data
  raw_data$adjuster_label <- factor(raw_data$adjuster,
                                    levels = classifier_adjuster_order$adjuster,
                                    labels = classifier_specific_labels)
  raw_data$dataset_label <- factor(raw_data$dataset_label, 
                                   levels = c("3 studies", "4 studies", "5 studies", "6 studies"))
  
  # Calculate mean for each adjuster
  mean_data <- raw_data %>%
    group_by(adjuster_label, dataset_label) %>%
    summarise(mean_value = mean(relative_value, na.rm = TRUE), .groups = "drop")
  
  # Calculate y-axis limits
  local_y_max <- max(raw_data$relative_value, na.rm = TRUE)
  local_y_min <- min(raw_data$relative_value, na.rm = TRUE)
  local_y_range <- local_y_max - local_y_min
  y_limits <- c(local_y_min - 0.05 * local_y_range, 
                local_y_max + 0.1 * local_y_range)
  
  p <- ggplot(raw_data, aes(x = adjuster_label, y = relative_value)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.6, alpha = 0.7) +
    # Add mean line for each adjuster
    geom_segment(data = mean_data, 
                aes(x = as.numeric(adjuster_label) - 0.3, 
                    xend = as.numeric(adjuster_label) + 0.3,
                    y = mean_value, yend = mean_value),
                color = "gray40", linewidth = 0.8) +
    # Add individual points colored by test study
    geom_point(aes(color = test_study, shape = test_study), 
              size = 3, alpha = 0.5) +
    facet_wrap(~ dataset_label, scales = "fixed", ncol = 4) +
    scale_x_discrete() +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0, 0))) +
    scale_color_brewer(palette = "Set2", name = "Test Study") +
    scale_shape_manual(values = c(15, 16, 17, 25, 18, 19), name = "Test Study") +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 9, face = "bold"),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey90", size = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 9),
      plot.title = element_text(size = 12, hjust = 0.5)
    ) +
    labs(
      y = sprintf("Relative MCC\n(vs. %s)", top_adjuster_label),
      title = classifier
    )
  
  relative_plot_list[[as.character(classifier)]] <- p
}

# Arrange plots in grid
relative_plot_vector <- unname(relative_plot_list)

relative_final_plot <- grid.arrange(
  grobs = relative_plot_vector,
  ncol = 2
)

# Save the relative performance plot
relative_output <- sub("\\.png$", "_relative.png", args$output)
cat("Saving relative performance plot to:", relative_output, "\n")
ggsave(
  filename = relative_output,
  plot = relative_final_plot,
  width = args$width,
  height = args$height,
  dpi = args$dpi,
  units = "in"
)

cat("Relative performance plot saved successfully!\n")
cat("Output file:", relative_output, "\n")