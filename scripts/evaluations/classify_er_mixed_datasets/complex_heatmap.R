# complex_heatmap.R
#
# This script creates 2x2 heatmaps of AUC scores and Matthews Correlation Coefficient (MCC)
# for dataset combinations from the ER classification results CSV files.
# The script automatically detects all available adjusters from CSV files matching:
# /outputs/metrics/er_classification_*.csv
# Usage: Rscript complex_heatmap.R
# No arguments required - processes all found adjusters automatically

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

# --- Auto-detect Available Adjusters ---
# Find all CSV files matching the pattern er_classification_*.csv
csv_files <- list.files("/outputs/metrics", pattern = "^er_classification_.*\\.csv$", full.names = FALSE)

if (length(csv_files) == 0) {
  cat("Error: No CSV files found matching pattern 'er_classification_*.csv' in /outputs/metrics/\n")
  quit(status = 1)
}

# Extract adjuster names from filenames
adjusters <- gsub("^er_classification_(.+)\\.csv$", "\\1", csv_files)
cat("Found adjusters:", paste(adjusters, collapse = ", "), "\n")

# --- Configuration ---
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

  # Check if data is empty
  if (nrow(input_data) == 0) {
    warning("Input CSV file is empty. No data to process.")
    return(data.frame())
  }

  input_data <- as.data.frame(input_data)

  # --- 🔧 Add MCC calculation ---
  input_data$MCC <- calculate_mcc(
    tp = input_data$`True Positive`,
    tn = input_data$`True Negative`,
    fp = input_data$`False Positive`,
    fn = input_data$`False Negative`
  )

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
prepare_metric_data <- function(df, metric_col) {
  df %>%
    group_by(Train, Test) %>%
    summarise(Mean_Metric = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop")
}

prepare_delta_metric_data <- function(df_adj, df_unadj, metric_col) {
  data_adj <- prepare_metric_data(df_adj, metric_col) %>% rename(Adj = Mean_Metric)
  data_unadj <- prepare_metric_data(df_unadj, metric_col) %>% rename(Unadj = Mean_Metric)

  full_join(data_adj, data_unadj, by = c("Train", "Test")) %>%
    mutate(Mean_Metric = Adj - Unadj)
}

get_platform_annotations <- function(datasets) {
  platforms <- dataset_to_platform[datasets]
  platforms[is.na(platforms)] <- "Unknown"

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

  types <- platform_type_map[platforms]
  types[is.na(types)] <- "Unknown"

  splits <- list(
    factor(types, levels = c("Microarray", "RNAseq", "Unknown")),
    factor(platforms, levels = unique(platforms))
  )

  colors <- c(
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

  annotation <- list(
    splits = splits,
    colors = colors,
    platforms = platforms
  )

  return(annotation)
}

draw_heatmap <- function(data_matrix, metric_col, adjuster, is_difference = FALSE) {
  row_anno_info <- get_platform_annotations(rownames(data_matrix))
  col_anno_info <- get_platform_annotations(colnames(data_matrix))

  row_ha <- rowAnnotation(
    Platform = row_anno_info$platforms,
    col = list(Platform = row_anno_info$colors),
    show_annotation_name = FALSE
  )

  col_ha <- HeatmapAnnotation(
    Platform = col_anno_info$platforms,
    col = list(Platform = col_anno_info$colors),
    show_annotation_name = FALSE
  )

  # Fix scale to [-1, 1] for all heatmaps
  col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#D62728", "#FFFFFF", "#2CA02C"))

  legend_title <- if (is_difference) {
    paste0("Diff ", metric_col)
  } else {
    metric_col
  }

  title_text <- if (is_difference) {
    paste0("Diff ", metric_col, ": ", adjuster, " - unadjusted")
  } else {
    paste0(metric_col, ": Dataset Combinations (", adjuster, ")")
  }

  # Optional: If metric is AUC, rescale the matrix before plotting
  if (metric_col == "ROC AUC" && !is_difference) {
    # Rescale AUC from [0,1] to [-1,1]
    data_matrix <- 2 * (data_matrix - 0.5)
  }

  ht <- Heatmap(data_matrix,
                name = legend_title,
                col = col_fun,
                na_col = "white",
                row_split = row_anno_info$splits,
                column_split = col_anno_info$splits,
                top_annotation = col_ha,
                left_annotation = row_ha,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 9),
                show_column_names = TRUE,
                column_names_gp = gpar(fontsize = 9),
                column_names_rot = 45,
                column_title = title_text,
                column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                row_title = "Train Dataset",
                row_title_gp = gpar(fontsize = 12),
                heatmap_legend_param = list(title = legend_title),
                heatmap_width = unit(1, "npc"),
                heatmap_height = unit(1, "npc"),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  val <- data_matrix[i, j]
                  if (!is.na(val)) {
                    grid.text(sprintf("%.2f", val), x, y,
                              gp = gpar(fontsize = 8, col = ifelse(abs(val) > 0.5, "white", "black")))
                  }
                })

  return(ht)
}

prepare_metric_matrix <- function(metric_data, metric_col) {
  all_datasets <- sort(union(metric_data$Train, metric_data$Test))
  
  # Ensure metric_data is unique for (Train, Test)
  metric_data_unique <- metric_data %>%
    group_by(Train, Test) %>%
    summarise(Mean_Metric = mean(Mean_Metric, na.rm = TRUE), .groups = "drop")
  
  metric_matrix <- expand.grid(Train = all_datasets, Test = all_datasets) %>%
    left_join(metric_data_unique, by = c("Train", "Test")) %>%
    pivot_wider(names_from = Test, values_from = Mean_Metric) %>%
    column_to_rownames("Train") %>%
    as.matrix()
  
  return(metric_matrix)
}

generate_all_heatmaps_to_pdf <- function(adjuster, fig_dir = "/outputs/figures") {
  file_adjusted <- paste0("/outputs/metrics/er_classification_", adjuster, ".csv")
  file_unadjusted <- "/outputs/metrics/er_classification_unadjusted.csv"

  df_adj <- read_and_prepare_data(file_adjusted) %>% filter_datasets()
  df_unadj <- read_and_prepare_data(file_unadjusted) %>% filter_datasets()

  # Prepare data for each metric
  heatmap_list <- list()

  # Δ MCC
  delta_mcc_data <- prepare_delta_metric_data(df_adj, df_unadj, "MCC")
  delta_mcc_matrix <- prepare_metric_matrix(delta_mcc_data, "MCC")
  heatmap_list[["delta_mcc"]] <- draw_heatmap(delta_mcc_matrix, "MCC", adjuster, is_difference = TRUE)

  # Δ AUC
  delta_auc_data <- prepare_delta_metric_data(df_adj, df_unadj, "ROC AUC")
  delta_auc_matrix <- prepare_metric_matrix(delta_auc_data, "ROC AUC")
  heatmap_list[["delta_auc"]] <- draw_heatmap(delta_auc_matrix, "ROC AUC", adjuster, is_difference = TRUE)

  # MCC
  mcc_data <- prepare_metric_data(df_adj, "MCC")
  mcc_matrix <- prepare_metric_matrix(mcc_data, "MCC")
  heatmap_list[["mcc"]] <- draw_heatmap(mcc_matrix, "MCC", adjuster)

  # AUC
  auc_data <- prepare_metric_data(df_adj, "ROC AUC")
  auc_matrix <- prepare_metric_matrix(auc_data, "ROC AUC")
  heatmap_list[["auc"]] <- draw_heatmap(auc_matrix, "ROC AUC", adjuster)

  ## Save all to one PDF
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(fig_dir, paste0("combined_heatmaps_", adjuster, ".pdf"))
  pdf(pdf_file, width = 14, height = 8)

  for (ht in heatmap_list) {
    draw(ht, padding = unit(c(10, 10, 10, 10), "mm"),
         merge_legend = TRUE,
         heatmap_legend_side = "right",
         annotation_legend_side = "right")
    grid.text("Test Dataset", x = unit(0.3, "npc"), y = unit(0.02, "npc"), gp = gpar(fontsize = 14))
    grid.newpage()
  }

  dev.off()
  cat("All heatmaps saved to:", pdf_file, "\n")
}

# Process all found adjusters
for (adjuster in adjusters) {
  cat("\nProcessing adjuster:", adjuster, "\n")
  
  # --- Configuration ---
  CSV_FILE <- paste0("/outputs/metrics/er_classification_", adjuster, ".csv")
  
  # --- Main Execution ---
  
  # Read and prepare data
  input_data <- read_and_prepare_data(CSV_FILE)

  # Display unique train/test combinations (after filtering)
  filtered_data <- filter_datasets(input_data)
  unique_trains <- unique(filtered_data$Train)
  unique_tests <- unique(filtered_data$Test)
  cat("Unique Train datasets:", paste(unique_trains, collapse = ", "), "\n")
  cat("Unique Test datasets:", paste(unique_tests, collapse = ", "), "\n\n")

  # Generate all heatmaps
  generate_all_heatmaps_to_pdf(adjuster)

  cat("All heatmaps generated successfully for adjuster:", adjuster, "\n")
}

cat("\nProcessing complete for all adjusters.\n")