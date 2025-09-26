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
if (!require(ComplexHeatmap)) install.packages("ComplexHeatmap")
if (!require(circlize)) install.packages("circlize")
if (!require(tibble)) install.packages("tibble")

# Load libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(tibble)

# --- Parse Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  cat("Usage: Rscript complex_heatmap.R <adjuster>\n")
  cat("Example: Rscript complex_heatmap.R unadjusted\n")
  cat("Example: Rscript complex_heatmap.R combat\n")
  quit(status = 1)
}

adjuster <- args[1]
cat("Processing adjuster:", adjuster, "\n")

# --- Configuration ---
CSV_FILE <- paste0("/outputs/metrics/er_classification_", adjuster, ".csv")
FIG_DIR <- "/outputs/figures"

platform_df <- read.csv("/scripts/evaluations/geo_metadata.csv")
platform_df$platform <- trimws(platform_df$platform)
dataset_to_platform <- setNames(platform_df$platform, platform_df$GSE_ID)

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
  input_data <- tryCatch({
    read_csv(csv_file, show_col_types = FALSE)
  }, error = function(e) {
    message(paste("Error reading file:", csv_file))
    message("Please ensure the CSV file exists in the current directory.")
    return(NULL)
  })
  
  if (is.null(input_data)) {
    stop("Aborting script due to missing or invalid input file.", call. = FALSE)
  }
  
  # Display data structure
  cat("Data dimensions:", nrow(input_data), "rows,", ncol(input_data), "columns\n")
  cat("Column names:", paste(colnames(input_data), collapse = ", "), "\n")
  
  # Check if data is empty
  if (nrow(input_data) == 0) {
    warning("Input CSV file is empty. No data to process.")
    return(data.frame())  # Return empty data frame
  }
  
  # Convert to regular data.frame for easier processing
  input_data <- as.data.frame(input_data)
  
  return(input_data)
}

# Function to filter data (common filtering logic)
filter_datasets <- function(input_data) {
  input_data %>%
    filter(!str_detect(Train, regex("combined", ignore_case = TRUE)),
           !str_detect(Test, regex("combined", ignore_case = TRUE)),
           !str_detect(Train, ";"),
           !str_detect(Test, ";"))
}

# Function to prepare metric data
prepare_metric_data <- function(input_data, metric_col) {
  # Apply common filtering first
  filtered_data <- filter_datasets(input_data)
  
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

# # --- Main Execution ---

# # Read and prepare data
input_data <- read_and_prepare_data(CSV_FILE)

generate_metric_difference_heatmap <- function(adjuster, metric_col, metric_display_name) {
  cat("\nGenerating", metric_display_name, "Difference Heatmap (", adjuster, " vs unadjusted)...\n")

  # Load both datasets
  adj_file <- paste0("/outputs/metrics/er_classification_", adjuster, ".csv")
  unadj_file <- "outputs/metrics/er_classification_ranked1.csv"

  adj_input <- read_and_prepare_data(adj_file)
  unadj_input <- read_and_prepare_data(unadj_file)

  adj_metric <- prepare_metric_data(adj_input, metric_col) %>%
    rename(Adj = Mean_Metric)
  unadj_metric <- prepare_metric_data(unadj_input, metric_col) %>%
    rename(Unadj = Mean_Metric)

  # Compute delta
  delta_data <- full_join(adj_metric, unadj_metric, by = c("Train", "Test")) %>%
    mutate(Delta = Adj - Unadj)

  # Create full matrix
  all_datasets <- sort(union(delta_data$Train, delta_data$Test))

  delta_matrix <- expand.grid(Train = all_datasets, Test = all_datasets) %>%
    left_join(delta_data, by = c("Train", "Test")) %>%
    pivot_wider(names_from = Test, values_from = Delta) %>%
    column_to_rownames("Train") %>%
    as.matrix()

  delta_matrix <- delta_matrix[all_datasets, all_datasets]

  # Platform annotations
  row_platform <- dataset_to_platform[rownames(delta_matrix)]
  col_platform <- dataset_to_platform[colnames(delta_matrix)]
  row_platform[is.na(row_platform)] <- "Unknown"
  col_platform[is.na(col_platform)] <- "Unknown"

  platform_type_map <- c(
    "Affymetrix Human Genome U133 Plus 2.0 Array" = "Microarray",
    "Affymetrix Human Genome U133A Array" = "Microarray",
    "Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]" = "Microarray",
    "Affymetrix Human Transcriptome Array 2.0 [probe set (exon) version] / Custom Affymetrix Human Transcriptome Array" = "Microarray",
    "Illumina HumanHT-12 V3.0 expression beadchip" = "Microarray",
    "Illumina Genome Analyzer (Homo sapiens)" = "RNAseq",
    "Illumina HiSeq 2000 (Homo sapiens)" = "RNAseq",
    "Illumina NextSeq 500 (Homo sapiens)" = "RNAseq",
    "Illumina HiSeq 2500 (Homo sapiens)" = "RNAseq"
  )

  row_platform_type <- platform_type_map[row_platform]
  col_platform_type <- platform_type_map[col_platform]
  row_platform_type[is.na(row_platform_type)] <- "Unknown"
  col_platform_type[is.na(col_platform_type)] <- "Unknown"

  row_split <- list(
    factor(row_platform_type, levels = c("Microarray", "RNAseq", "Unknown")),
    factor(row_platform, levels = unique(row_platform))
  )

  col_split <- list(
    factor(col_platform_type, levels = c("Microarray", "RNAseq", "Unknown")),
    factor(col_platform, levels = unique(col_platform))
  )

  platform_colors <- c(
    "Affymetrix Human Genome U133 Plus 2.0 Array" = "#56B4E9FF",
    "Affymetrix Human Genome U133A Array" = "#009E73FF",
    "Illumina HiSeq 2000 (Homo sapiens)" = "#E69F00FF",
    "Illumina NextSeq 500 (Homo sapiens)" = "#D55E00FF",
    "Illumina HiSeq 2500 (Homo sapiens)" = "#CC79A7FF",
    "Illumina Genome Analyzer (Homo sapiens)" = "#F0E442FF",
    "Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]" = "#0072B2FF",
    "Affymetrix Human Transcriptome Array 2.0 [probe set (exon) version] / Custom Affymetrix Human Transcriptome Array" = "#4682B4FF",
    "Illumina HumanHT-12 V3.0 expression beadchip" = "#6A9FB5FF",
    "Unknown" = "#000000"
  )

  row_ha <- rowAnnotation(
    Platform = row_platform,
    col = list(Platform = platform_colors),
    show_annotation_name = FALSE
  )
  col_ha <- HeatmapAnnotation(
    Platform = col_platform,
    col = list(Platform = platform_colors),
    show_annotation_name = FALSE
  )

  # Color scale for differences
  diff_range <- if (metric_col == "MCC") c(-1, 0, 1) else c(-0.5, 0, 0.5)
  col_fun <- circlize::colorRamp2(diff_range, c("#d73027", "#fad6b2ff", "#66c2a5"))

  # Title and legend
  title_text <- paste0("Δ ", metric_display_name, ": ", adjuster, " - unadjusted")
  legend_title <- paste0("Δ ", metric_display_name)

  ht <- Heatmap(delta_matrix,
    name = legend_title,
    col = col_fun,
    na_col = "white",
    row_split = row_split,
    column_split = col_split,
    top_annotation = col_ha,
    left_annotation = row_ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 10, fontface = "plain"),
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 10, fontface = "plain"),
    column_names_rot = 45,
    column_title = title_text,
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    row_title = "Train Dataset",
    row_title_gp = gpar(fontsize = 14),
    heatmap_legend_param = list(title = legend_title),
    heatmap_width = unit(1, "npc"),
    heatmap_height = unit(1, "npc"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      value <- delta_matrix[i, j]
      if (!is.na(value)) {
        grid.text(sprintf("%.2f", value), x, y,
          gp = gpar(fontsize = 10, col = ifelse(abs(value) > 0.2, "white", "black"))
        )
      }
    }
  )

  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  pdf_file <- file.path(FIG_DIR, paste0("delta_", tolower(metric_col), "_heatmap_", adjuster, ".pdf"))
  pdf(pdf_file, width = 14, height = 8)
  draw(ht, padding = unit(c(10, 10, 10, 10), "mm"), merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
  grid.text("Test Dataset", x = unit(0.3, "npc"), y = unit(0.02, "npc"), gp = gpar(fontsize = 14))
  dev.off()

  cat("Saved:", pdf_file, "\n")
}

# Function to generate heatmap for a specific metric
generate_metric_heatmap <- function(input_data, metric_name, metric_col, adjuster) {
    cat("Generating", metric_name, "Heatmap...\n")
    
    # Prepare data
    metric_data <- prepare_metric_data(input_data, metric_col)
    
    # Create list of all datasets to ensure consistent ordering
    all_datasets <- sort(union(unique(metric_data$Train), unique(metric_data$Test)))

    # Create matrix with all datasets, filling missing values with NA
    metric_matrix <- expand.grid(Train = all_datasets, Test = all_datasets) %>%
    left_join(metric_data, by = c("Train", "Test")) %>%
    pivot_wider(names_from = Test, values_from = Mean_Metric) %>%
    column_to_rownames("Train") %>%
    as.matrix()


    metric_matrix <- metric_matrix[all_datasets, all_datasets]

    # Order the datasets
    row_platform <- dataset_to_platform[rownames(metric_matrix)]
    col_platform <- dataset_to_platform[colnames(metric_matrix)]

    row_platform[is.na(row_platform)] <- "Unknown"
    col_platform[is.na(col_platform)] <- "Unknown"

    platform_type_map <- c(
        "Affymetrix Human Genome U133 Plus 2.0 Array" = "Microarray",
        "Affymetrix Human Genome U133A Array" = "Microarray",
        "Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]" = "Microarray",
        "Affymetrix Human Transcriptome Array 2.0 [probe set (exon) version] / Custom Affymetrix Human Transcriptome Array" = "Microarray",
        "Illumina HumanHT-12 V3.0 expression beadchip" = "Microarray",
        "Illumina Genome Analyzer (Homo sapiens)" = "RNAseq",
        "Illumina HiSeq 2000 (Homo sapiens)" = "RNAseq",
        "Illumina NextSeq 500 (Homo sapiens)" = "RNAseq",
        "Illumina HiSeq 2500 (Homo sapiens)" = "RNAseq"
    )

    row_platform_type <- platform_type_map[row_platform]
    col_platform_type <- platform_type_map[col_platform]

    row_platform_type[is.na(row_platform_type)] <- "Unknown"
    col_platform_type[is.na(col_platform_type)] <- "Unknown"

    # Group by platform type and platform name (no manual ordering)
    row_split <- list(
      factor(row_platform_type, levels = c("Microarray", "RNAseq", "Unknown")),
      factor(row_platform, levels = unique(row_platform))
    )

    col_split <- list(
      factor(col_platform_type, levels = c("Microarray", "RNAseq", "Unknown")),
      factor(col_platform, levels = unique(col_platform))
    )

    # Define colors
    platform_colors <- c(
        "Affymetrix Human Genome U133 Plus 2.0 Array" = "#56B4E9FF",
        "Affymetrix Human Genome U133A Array" = "#009E73FF",
        "Illumina HiSeq 2000 (Homo sapiens)" = "#E69F00FF",
        "Illumina NextSeq 500 (Homo sapiens)" = "#D55E00FF",
        "Illumina HiSeq 2500 (Homo sapiens)" = "#CC79A7FF",
        "Illumina Genome Analyzer (Homo sapiens)" = "#F0E442FF",
        "Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]" = "#0072B2FF",
        "Affymetrix Human Transcriptome Array 2.0 [probe set (exon) version] / Custom Affymetrix Human Transcriptome Array" = "#4682B4FF",
        "Illumina HumanHT-12 V3.0 expression beadchip" = "#6A9FB5FF",
        "Unknown" = "#000000"
    )

    # NA columns as white
    na_col <- "white"

    # Annotations for rows and columns
    row_ha <- rowAnnotation(
        Platform = row_platform,
        col = list(Platform = platform_colors),
        show_annotation_name = FALSE
    )
    col_ha <- HeatmapAnnotation(
        Platform = col_platform,
        col = list(Platform = platform_colors),
        show_annotation_name = FALSE
    )

    # Define color function for heatmap fill
    if (metric_col == "MCC") {
        col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#d73027", "#fad6b2ff", "#66c2a5"))
    } else {
        col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#d73027", "#fad6b2ff", "#66c2a5"))
    }

    

    # Update titles
    if (metric_col == "MCC") {
        title_text <- paste0("Matthews Correlation Coefficient Heatmap: Dataset Combinations (", adjuster, ")")
        legend_title <- "Mean MCC"
        col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#d73027", "#fad6b2ff", "#66c2a5"))
    }   else {
        title_text <- paste0("ROC AUC Heatmap: Dataset Combinations (", adjuster, ")")
        legend_title <- "Mean ROC AUC"
        col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#d73027", "#fad6b2ff", "#66c2a5"))
    }

    # Create heatmap object
    ht <- Heatmap(metric_matrix,
                    name = metric_name, 
                    col = col_fun,
                    na_col = na_col,
                    row_split = row_split,
                    column_split = col_split, 
                    top_annotation = col_ha,
                    left_annotation = row_ha,
                    cluster_rows = FALSE, 
                    cluster_columns = FALSE, 
                    show_row_names = TRUE,
                    row_names_gp = gpar(fontsize=10, fontface="plain"),
                    show_column_names = TRUE, 
                    column_names_gp = gpar(fontsize=10, fontface="plain"),
                    column_names_rot = 45,
                    heatmap_legend_param = list(title = legend_title),
                    column_title = title_text,
                    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
                    row_title = "Train Dataset",
                    row_title_gp = gpar(fontsize = 14),
                    heatmap_width = unit(1, "npc"),
                    heatmap_height = unit(1, "npc"),
                    cell_fun = function(j, i, x, y, width, height, fill) {
                value <- metric_matrix[i, j]
                if (!is.na(value)) {
                  grid.text(sprintf("%.2f", value), x, y, 
                            gp = gpar(fontsize = 10, col = ifelse(abs(value) > 0.5, "white", "black")))
                }
              }
    )


    # Save heatmap to PDF
    dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
    pdf_file <- file.path(FIG_DIR, paste0(tolower(gsub(" ", "_", metric_name)), "_heatmap_", adjuster, ".pdf"))

    pdf(pdf_file, width = 14, height = 8)
    draw(ht, padding = unit(c(10, 10, 10, 10), "mm"), merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")

    # Add bottom label 
    grid.text("Test Dataset",
        x = unit(0.3, "npc"),
        y = unit(0.02, "npc"),
        gp = gpar(fontsize = 14))

    dev.off()

    cat("Heatmap saved to: ", pdf_file, "\n")

    # Print summary
    cat(paste("\nSummary of", metric_name, "scores:\n"))
    print(summary(metric_data$Mean_Metric))

    return(metric_data)
  
}

# Display unique train/test combinations (after filtering)
filtered_data <- filter_datasets(input_data)
unique_trains <- unique(filtered_data$Train)
unique_tests <- unique(filtered_data$Test)
cat("Unique Train datasets:", paste(unique_trains, collapse = ", "), "\n")
cat("Unique Test datasets:", paste(unique_tests, collapse = ", "), "\n\n")

# Generate all heatmaps
auc_data <- generate_metric_heatmap(input_data, "ROC AUC", "ROC AUC", adjuster)
mcc_data <- generate_metric_heatmap(input_data, "MCC", "MCC", adjuster)
if (adjuster != "unadjusted") {
  generate_metric_difference_heatmap(adjuster, "ROC AUC", "ROC AUC")
  generate_metric_difference_heatmap(adjuster, "MCC", "MCC")
}

cat("All heatmaps generated successfully for adjuster:", adjuster, "\n")