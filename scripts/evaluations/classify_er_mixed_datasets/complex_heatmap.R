# complex_heatmap.R
#
# This script creates 2x2 heatmaps of AUC scores and Matthews Correlation Coefficient (MCC)
# for dataset combinations from the ER classification results CSV file.
# Usage: Rscript train_test_heatmap.R 
# Example: Rscript train_test_heatmap.R 

# Install necessary packages if they are not already installed
# if (!require(ggplot2)) install.packages("ggplot2")
# if (!require(readr)) install.packages("readr")
# if (!require(dplyr)) install.packages("dplyr")
# if (!require(tidyr)) install.packages("tidyr")
# if (!require(stringr)) install.packages("stringr")
# if (!require(ComplexHeatmap)) install.packages("ComplexHeatmap")
# if (!require(circlize)) install.packages("circlize")
# if (!require(tibble)) install.packages("tibble")

# Load libraries
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(purrr)

# --- Auto-detect Available Adjusters ---
# Find all CSV files matching the pattern er_classification_*.csv
csv_files <- list.files("/outputs/metrics", pattern = "^er_classification_.*\\.csv$", full.names = FALSE)

# Extract adjuster names from filenames
all_adjusters <- gsub("^er_classification_(.+)\\.csv$", "\\1", csv_files)
cat("Found potential adjusters:", paste(all_adjusters, collapse = ", "), "\n")

# Function to check if a CSV file has meaningful data (more than just headers)
has_meaningful_data <- function(csv_file) {
  tryCatch({
    data <- read_csv(csv_file, show_col_types = FALSE)
    return(nrow(data) > 0)
  }, error = function(e) {
    return(FALSE)
  })
}

# Filter out adjusters with empty files
adjusters <- c()
for (adj in all_adjusters) {
  csv_path <- file.path("/outputs/metrics", paste0("er_classification_", adj, ".csv"))
  if (has_meaningful_data(csv_path)) {
    adjusters <- c(adjusters, adj)
  } else {
    cat("Skipping adjuster", adj, "- file is empty or contains only headers\n")
  }
}

cat("Adjusters with data:", paste(adjusters, collapse = ", "), "\n")

# Check if we have any adjusters with data
if (length(adjusters) == 0) {
  stop("No adjusters found with meaningful data. All CSV files appear to be empty or contain only headers.", call. = FALSE)
}

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

prepare_metric_matrix <- function(metric_data, metric_col) {
  # Check if metric_data is empty
  if (nrow(metric_data) == 0) {
    cat("Warning: No data available for", metric_col, "\n")
    return(NULL)
  }
  
  all_datasets <- sort(union(metric_data$Train, metric_data$Test))
  
  # Check if we have any datasets
  if (length(all_datasets) == 0) {
    cat("Warning: No datasets found for", metric_col, "\n")
    return(NULL)
  }
  
  # Ensure metric_data is unique for (Train, Test)
  metric_data_unique <- metric_data %>%
    group_by(Train, Test) %>%
    summarise(Mean_Metric = mean(Mean_Metric, na.rm = TRUE), .groups = "drop")
  
  metric_matrix <- expand.grid(Train = all_datasets, Test = all_datasets) %>%
    left_join(metric_data_unique, by = c("Train", "Test")) %>%
    pivot_wider(names_from = Test, values_from = Mean_Metric) %>%
    column_to_rownames("Train") %>%
    as.matrix()
  
  # Check if matrix is valid
  if (nrow(metric_matrix) == 0 || ncol(metric_matrix) == 0) {
    cat("Warning: Empty matrix created for", metric_col, "\n")
    return(NULL)
  }
  
  cat("Created matrix for", metric_col, "with dimensions:", nrow(metric_matrix), "x", ncol(metric_matrix), "\n")
  return(metric_matrix)
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

  # Check if data is empty or only contains headers
  if (nrow(input_data) == 0) {
    warning(paste("Input CSV file", csv_file, "is empty. No data to process."))
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
filter_datasets <- function(input_data, train_combined) {
  if (train_combined) {
    # Use train on both datasets
    result <- input_data %>%
      filter(str_detect(Train, ";")) %>%
      filter(!str_detect(Test, ";")) %>%
      mutate(
        Train = map2_chr(Train, Test, function(train_val, test_val) {
          train_parts <- str_split(train_val, ";")[[1]]
          # Return the part that's NOT the test value
          train_parts[train_parts != test_val][1]
        })
      )
    return(result)
  } else {
    # Use cross-training data
    return(input_data %>%
      filter(!str_detect(Train, ";"), !str_detect(Test, ";")))
  }
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

draw_heatmap <- function(data_matrix, metric_col, adjuster, train_combined, is_difference = FALSE) {
  # Check if matrix is NULL or empty
  if (is.null(data_matrix) || nrow(data_matrix) == 0 || ncol(data_matrix) == 0) {
    cat("Skipping heatmap for", metric_col, "- no valid data matrix\n")
    return(NULL)
  }
  
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

  train_text <- if (train_combined) {
    " (Trained on Combined)"
  } else {
    " (Cross-Trained)"
  }

  title_text <- if (is_difference) {
    paste0("Diff ", metric_col, ": ", adjuster, " - unadjusted", train_text)
  } else {
    paste0(metric_col, ": Dataset Combinations (", adjuster, ")", train_text)
  }

  # Optional: If metric is AUC, rescale the matrix before plotting
  if (metric_col == "ROC AUC" && !is_difference) {
    # Rescale AUC from [0,1] to [-1,1]
    data_matrix <- 2 * (data_matrix - 0.5)
  }

  row_title <- if (train_combined) {
    "Training Dataset, in Combination with Test"
  } else {
    "Train Dataset"
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
                row_title = row_title,
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

# # --- Main Execution ---
generate_jitter_plot <- function(all_diff_data, fig_dir) {
  p <- ggplot(all_diff_data, aes(x = Adjuster, y = Mean_Metric, color = Metric)) +
  geom_jitter(width = 0.25, height = 0, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("MCC" = "skyblue", "AUC" = "orange")) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", linewidth = 1, fatten = 1) +
  theme_minimal() +
  labs(
    title = "Distribution of Metric Differences for Adjusters",
    x = "Adjuster",
    y = "Difference in Metric (Adjusted - Unadjusted)",
    fill = "Metric"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12)
  ) +
  facet_wrap(~ Metric, scales = "free_y")

  ggsave(file.path(fig_dir, "jitter_plot_adjusters.png"), plot = p, width = 10, height = 6)
  cat("Jitter plot saved to:", file.path(fig_dir, "jitter_plot_adjusters.png"), "\n")
}

generate_all_heatmaps_to_pdf <- function(adjuster, train_combined, fig_dir = "/outputs/figures") {
  file_adjusted <- paste0("/outputs/metrics/er_classification_", adjuster, ".csv")
  file_unadjusted <- "/outputs/metrics/er_classification_unadjusted.csv"

  df_adj <- read_and_prepare_data(file_adjusted) %>% filter_datasets(train_combined)
  df_unadj <- read_and_prepare_data(file_unadjusted) %>% filter_datasets(train_combined)

  # Check if we have data after filtering
  if (nrow(df_adj) == 0 || nrow(df_unadj) == 0) {
    cat("Skipping", adjuster, "- no data after filtering for train_combined =", train_combined, "\n")
    return()
  }

  # Prepare data for each metric
  heatmap_list <- list()

  # Δ MCC
  delta_mcc_data <- prepare_delta_metric_data(df_adj, df_unadj, "MCC")
  delta_mcc_matrix <- prepare_metric_matrix(delta_mcc_data, "MCC")
  if (!is.null(delta_mcc_matrix)) {
    heatmap_list[["delta_mcc"]] <- draw_heatmap(delta_mcc_matrix, "MCC", adjuster, train_combined, is_difference = TRUE)
  }

  # Δ AUC
  delta_auc_data <- prepare_delta_metric_data(df_adj, df_unadj, "ROC AUC")
  delta_auc_matrix <- prepare_metric_matrix(delta_auc_data, "ROC AUC")
  if (!is.null(delta_auc_matrix)) {
    heatmap_list[["delta_auc"]] <- draw_heatmap(delta_auc_matrix, "ROC AUC", adjuster, train_combined, is_difference = TRUE)
  }

  # MCC
  mcc_data <- prepare_metric_data(df_adj, "MCC")
  mcc_matrix <- prepare_metric_matrix(mcc_data, "MCC")
  if (!is.null(mcc_matrix)) {
    heatmap_list[["mcc"]] <- draw_heatmap(mcc_matrix, "MCC", adjuster, train_combined)
  }

  # AUC
  auc_data <- prepare_metric_data(df_adj, "ROC AUC")
  auc_matrix <- prepare_metric_matrix(auc_data, "ROC AUC")
  if (!is.null(auc_matrix)) {
    heatmap_list[["auc"]] <- draw_heatmap(auc_matrix, "ROC AUC", adjuster, train_combined)
  }

  # Filter out NULL heatmaps
  heatmap_list <- heatmap_list[!sapply(heatmap_list, is.null)]
  
  if (length(heatmap_list) == 0) {
    cat("No valid heatmaps to save for", adjuster, "with train_combined =", train_combined, "\n")
    return()
  }

  ## Save all to one PDF
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  combined_suffix = ifelse(train_combined, "_train_combined", "")
  pdf_file <- file.path(fig_dir, paste0("combined_heatmaps_", adjuster, combined_suffix, ".pdf"))
  pdf(pdf_file, width = 14, height = 8)

  for (ht in heatmap_list) {
    if (!is.null(ht)) {
      draw(ht, padding = unit(c(10, 10, 10, 10), "mm"),
           merge_legend = TRUE,
           heatmap_legend_side = "right",
           annotation_legend_side = "right")
      grid.text("Test Dataset", x = unit(0.3, "npc"), y = unit(0.02, "npc"), gp = gpar(fontsize = 14))
      grid.newpage()
    }
  }

  dev.off()
  cat("All heatmaps saved to:", pdf_file, "\n")
}

all_adjuster_diffs <- data.frame()
for (adjuster in adjusters) {
  cat("\n=== Processing adjuster:", adjuster, "===\n")
  CSV_FILE <- paste0("/outputs/metrics/er_classification_", adjuster, ".csv")
  file_adjusted <- paste0("/outputs/metrics/er_classification_", adjuster, ".csv")
  file_unadjusted <- "/outputs/metrics/er_classification_unadjusted.csv"

  for(train_combined in c(TRUE, FALSE)) {
    cat("Processing train_combined =", train_combined, "\n")
    
    df_adj <- read_and_prepare_data(file_adjusted) %>% filter_datasets(train_combined)
    df_unadj <- read_and_prepare_data(file_unadjusted) %>% filter_datasets(train_combined)

    cat("After filtering - df_adj rows:", nrow(df_adj), ", df_unadj rows:", nrow(df_unadj), "\n")

    if (nrow(df_adj) == 0 || nrow(df_unadj) == 0) {
      cat("Skipping due to empty data after filtering\n")
      next
    }

    # Prepare the differences
    delta_mcc_data <- prepare_delta_metric_data(df_adj, df_unadj, "MCC")
    delta_auc_data <- prepare_delta_metric_data(df_adj, df_unadj, "ROC AUC")

    # Add Adjuster column to each
    delta_mcc_data$Adjuster <- adjuster
    delta_auc_data$Adjuster <- adjuster

    # Combine the data
    if (train_combined) {
      all_adjuster_diffs <- bind_rows(
        all_adjuster_diffs,
        delta_mcc_data %>% mutate(Metric = "MCC"),
        delta_auc_data %>% mutate(Metric = "AUC")
      )
    }

    # Generate Heatmaps for the current adjuster
    tryCatch({
      generate_all_heatmaps_to_pdf(adjuster, train_combined, FIG_DIR)
    }, error = function(e) {
      cat("Error generating heatmaps for", adjuster, "with train_combined =", train_combined, ":", e$message, "\n")
    })
  }
}
 
# Generate jitter Plot
generate_jitter_plot(all_adjuster_diffs, FIG_DIR)

cat("All heatmaps and scatter plot generated successfully for adjuster:", adjuster, "\n")

# Display unique train/test combinations (after filtering)
combined_data <- bind_rows(df_adj, df_unadj)

filtered_data <- filter_datasets(combined_data, FALSE)
unique_trains <- unique(filtered_data$Train)
unique_tests <- unique(filtered_data$Test)
cat("Unique Train datasets:", paste(unique_trains, collapse = ", "), "\n")
cat("Unique Test datasets:", paste(unique_tests, collapse = ", "), "\n\n")



cat("All heatmaps generated successfully.")