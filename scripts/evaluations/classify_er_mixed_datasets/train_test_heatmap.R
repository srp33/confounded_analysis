# train_test_heatmap.R
#
# This script creates 2x2 heatmaps of AUC scores and Matthews Correlation Coefficient (MCC)
# for dataset combinations from the ER classification results CSV file.
# Usage: Rscript train_test_heatmap.R <adjuster>
# Example: Rscript train_test_heatmap.R unadjusted

# Install necessary packages if they are not already installed
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(readr)) install.packages("readr")
if (!require(dplyr)) install.packages("dplyr")
if (!require(tidyr)) install.packages("tidyr")

# Load libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)

# --- Parse Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  cat("Usage: Rscript train_test_heatmap.R <adjuster>\n")
  cat("Example: Rscript train_test_heatmap.R unadjusted\n")
  cat("Example: Rscript train_test_heatmap.R combat\n")
  quit(status = 1)
}

adjuster <- args[1]
cat("Processing adjuster:", adjuster, "\n")

# --- Configuration ---
CSV_FILE <- paste0("/outputs/metrics/er_classification_all_", adjuster, ".csv")
FIG_DIR <- "/outputs/figures"

# --- Helper Functions ---

# Function to calculate Matthews Correlation Coefficient (vectorized)
calculate_mcc <- function(tp, tn, fp, fn) {
  numerator <- (tp * tn) - (fp * fn)
  denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  
  # Handle division by zero (vectorized)
  result <- ifelse(denominator == 0, 0, numerator / denominator)
  
  return(result)
}

# Function to read and prepare data
read_and_prepare_data <- function(csv_file) {
  cat("Reading data from:", csv_file, "\n")
  data <- tryCatch({
    read_csv(csv_file, show_col_types = FALSE)
  }, error = function(e) {
    message(paste("Error reading file:", csv_file))
    message("Please ensure the CSV file exists in the current directory.")
    return(NULL)
  })
  
  if (is.null(data)) {
    stop("Aborting script due to missing or invalid input file.", call. = FALSE)
  }
  
  # Display data structure
  cat("Data dimensions:", nrow(data), "rows,", ncol(data), "columns\n")
  cat("Column names:", paste(colnames(data), collapse = ", "), "\n")
  
  return(data)
}

# Function to prepare metric data
prepare_metric_data <- function(data, metric_col) {
  if (metric_col == "MCC") {
    # For MCC, we need to calculate it from confusion matrix values
    metric_data <- data %>%
      # Remove entries where Train or Test contains "combined"
      filter(!grepl("combined", Train, ignore.case = TRUE)) %>%
      filter(!grepl("combined", Test, ignore.case = TRUE)) %>%
      # Remove entries where Train or Test contains a colon
      filter(!grepl(":", Train)) %>%
      filter(!grepl(":", Test)) %>%
      # Filter for rows with complete confusion matrix data
      filter(!is.na(`True Positive`) & !is.na(`True Negative`) & 
             !is.na(`False Positive`) & !is.na(`False Negative`)) %>%
      # Calculate MCC from confusion matrix values
      mutate(MCC = calculate_mcc(`True Positive`, `True Negative`, 
                                `False Positive`, `False Negative`)) %>%
      select(Train, Test, MCC) %>%
      group_by(Train, Test) %>%
      summarise(Mean_Metric = mean(MCC, na.rm = TRUE), .groups = 'drop')
  } else {
    # For other metrics like ROC AUC
    metric_data <- data %>%
      filter(!is.na(.data[[metric_col]])) %>%
      # Remove entries where Train or Test contains "combined"
      filter(!grepl("combined", Train, ignore.case = TRUE)) %>%
      filter(!grepl("combined", Test, ignore.case = TRUE)) %>%
      # Remove entries where Train or Test contains a colon
      filter(!grepl(":", Train)) %>%
      filter(!grepl(":", Test)) %>%
      select(Train, Test, all_of(metric_col)) %>%
      group_by(Train, Test) %>%
      summarise(Mean_Metric = mean(.data[[metric_col]], na.rm = TRUE), .groups = 'drop')
  }
  
  return(metric_data)
}

# Function to create heatmap
create_heatmap <- function(metric_data, title, subtitle, legend_name, 
                          midpoint = 0.5, limits = c(0, 1), is_mcc = FALSE) {
  
  base_plot <- ggplot(metric_data, aes(x = Test, y = Train, fill = Mean_Metric)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", Mean_Metric)), 
              color = "white", size = 3.5, fontface = "bold") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Test Dataset",
      y = "Train Dataset"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 0),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed(ratio = 1)
  
  # Use different color schemes for MCC vs other metrics
  if (is_mcc) {
    # For MCC: use a scheme where 0 is clearly visible (light grey)
    base_plot + 
      scale_fill_gradient2(name = legend_name, 
                           low = "darkred",     # negative MCC (poor performance)
                           mid = "lightgrey",   # MCC = 0 (no correlation)
                           high = "darkgreen",  # positive MCC (good performance)
                           midpoint = 0,
                           limits = limits)
  } else {
    # For other metrics like AUC: use the original blue-white-red scheme
    base_plot + 
      scale_fill_gradient2(name = legend_name, 
                           low = "darkred", 
                           mid = "lightgrey", 
                           high = "darkgreen",
                           midpoint = midpoint,
                           limits = limits)
  }
}

# Function to save plot and display summary
save_plot_and_summary <- function(plot, filename, metric_data, metric_name) {
  # Create output directory if it doesn't exist
  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  
  # Save the plot
  output_path <- file.path(FIG_DIR, filename)
  ggsave(
    output_path,
    plot = plot,
    device = "pdf",
    width = 10,
    height = 8,
    units = "in"
  )
  
  cat("Heatmap saved to:", output_path, "\n")
  
  # Display summary statistics
  cat(paste("\nSummary of", metric_name, "scores:\n"))
  print(summary(metric_data$Mean_Metric))
  
  # Display the data table
  cat(paste("\n", metric_name, "Matrix:\n"))
  metric_matrix <- metric_data %>%
    pivot_wider(names_from = Test, values_from = Mean_Metric)
  print(metric_matrix)
  
  cat("\n", rep("=", 50), "\n\n")
}

# --- Main Execution ---

# Read and prepare data
data <- read_and_prepare_data(CSV_FILE)

# Display unique train/test combinations (after filtering)
unique_trains <- unique(data$Train[!grepl("combined", data$Train, ignore.case = TRUE) & 
                                   !grepl(":", data$Train)])
unique_tests <- unique(data$Test[!grepl("combined", data$Test, ignore.case = TRUE) & 
                                 !grepl(":", data$Test)])
cat("Unique Train datasets:", paste(unique_trains, collapse = ", "), "\n")
cat("Unique Test datasets:", paste(unique_tests, collapse = ", "), "\n\n")

# Generate AUC Heatmap
cat("Generating ROC AUC Heatmap...\n")
auc_data <- prepare_metric_data(data, "ROC AUC")
auc_plot <- create_heatmap(
  auc_data, 
  paste("ROC AUC Heatmap: Dataset Combinations (", adjuster, ")", sep = ""),
  "Mean AUC scores for Train/Test dataset pairs",
  "Mean\nROC AUC",
  midpoint = 0.5,
  limits = c(0, 1)
)
auc_filename <- paste0("auc_heatmap_", adjuster, ".pdf")
save_plot_and_summary(auc_plot, auc_filename, auc_data, "ROC AUC")

# Generate MCC Heatmap
cat("Generating MCC Heatmap...\n")
mcc_data <- prepare_metric_data(data, "MCC")
mcc_plot <- create_heatmap(
  mcc_data,
  paste("Matthews Correlation Coefficient Heatmap: Dataset Combinations (", adjuster, ")", sep = ""), 
  "Mean MCC scores for Train/Test dataset pairs",
  "Mean\nMCC",
  midpoint = 0,
  limits = c(-1, 1),
  is_mcc = TRUE
)
mcc_filename <- paste0("mcc_heatmap_", adjuster, ".pdf")
save_plot_and_summary(mcc_plot, mcc_filename, mcc_data, "MCC")

cat("All heatmaps generated successfully for adjuster:", adjuster, "\n")