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
if (!require(stringr)) install.packages("stringr")

# Load libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

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
CSV_FILE <- paste0("/outputs/metrics/er_classification_", adjuster, ".csv")
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
  
  # Check if data is empty
  if (nrow(data) == 0) {
    warning("Input CSV file is empty. No data to process.")
    return(data.frame())  # Return empty data frame
  }
  
  # Convert to regular data.frame for easier processing
  data <- as.data.frame(data)
  
  return(data)
}

# Function to filter data (common filtering logic)
filter_datasets <- function(data) {
  data %>%
    filter(!str_detect(Train, regex("combined", ignore_case = TRUE)),
           !str_detect(Test, regex("combined", ignore_case = TRUE)),
           !str_detect(Train, ";"),
           !str_detect(Test, ";"))
}

# Function to prepare metric data
prepare_metric_data <- function(data, metric_col) {
  # Apply common filtering first
  filtered_data <- filter_datasets(data)
  
  if (metric_col == "MCC") {
    # For MCC, calculate from confusion matrix values
    filtered_data %>%
      filter(!is.na(`True Positive`) & !is.na(`True Negative`) & 
             !is.na(`False Positive`) & !is.na(`False Negative`)) %>%
      mutate(!!metric_col := calculate_mcc(`True Positive`, `True Negative`, 
                                          `False Positive`, `False Negative`)) %>%
      select(Train, Test, all_of(metric_col)) %>%
      group_by(Train, Test) %>%
      summarise(Mean_Metric = mean(.data[[metric_col]], na.rm = TRUE), .groups = 'drop')
  } else {
    # For existing metrics like ROC AUC
    filtered_data %>%
      filter(!is.na(.data[[metric_col]])) %>%
      select(Train, Test, all_of(metric_col)) %>%
      group_by(Train, Test) %>%
      summarise(Mean_Metric = mean(.data[[metric_col]], na.rm = TRUE), .groups = 'drop')
  }
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

# Function to generate heatmap for a specific metric
generate_metric_heatmap <- function(data, metric_name, metric_col, adjuster) {
  cat("Generating", metric_name, "Heatmap...\n")
  
  # Prepare data
  metric_data <- prepare_metric_data(data, metric_col)
  
  # Set parameters based on metric type
  if (metric_col == "MCC") {
    title <- paste("Matthews Correlation Coefficient Heatmap: Dataset Combinations (", adjuster, ")", sep = "")
    subtitle <- "Mean MCC scores for Train/Test dataset pairs"
    legend_name <- "Mean\nMCC"
    midpoint <- 0
    limits <- c(-1, 1)
    is_mcc <- TRUE
    filename <- paste0("mcc_heatmap_", adjuster, ".pdf")
  } else {
    title <- paste("ROC AUC Heatmap: Dataset Combinations (", adjuster, ")", sep = "")
    subtitle <- "Mean AUC scores for Train/Test dataset pairs"
    legend_name <- "Mean\nROC AUC"
    midpoint <- 0.5
    limits <- c(0, 1)
    is_mcc <- FALSE
    filename <- paste0("auc_heatmap_", adjuster, ".pdf")
  }
  
  # Create plot
  plot <- create_heatmap(metric_data, title, subtitle, legend_name, 
                        midpoint, limits, is_mcc)
  
  # Save and summarize
  save_plot_and_summary(plot, filename, metric_data, metric_name)
  
  return(metric_data)
}

# Display unique train/test combinations (after filtering)
filtered_data <- filter_datasets(data)
unique_trains <- unique(filtered_data$Train)
unique_tests <- unique(filtered_data$Test)
cat("Unique Train datasets:", paste(unique_trains, collapse = ", "), "\n")
cat("Unique Test datasets:", paste(unique_tests, collapse = ", "), "\n\n")

# Generate both heatmaps
auc_data <- generate_metric_heatmap(data, "ROC AUC", "ROC AUC", adjuster)
mcc_data <- generate_metric_heatmap(data, "MCC", "MCC", adjuster)

cat("All heatmaps generated successfully for adjuster:", adjuster, "\n")