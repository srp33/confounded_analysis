# plot_performance.R
#
# This script combines multiple adjuster CSVs into a single dataframe,
# calculates MCC, computes deltas vs unadjusted, and generates scaling plots.

# --- Load Libraries ---
library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(patchwork)
library(argparse)

# --- Parse Arguments ---
parser <- ArgumentParser(description = "Generate performance plots from aggregated data.")

parser$add_argument('--metrics_file', required = TRUE,
    help = "Path to CSV containing aggregated metrics")
parser$add_argument('--figures_dir', required = TRUE, 
    help = "Directory to save figures")
parser$add_argument('--metadata_file', required = TRUE,
    help = "Path to CSV containing dataset metadata")
parser$add_argument('--cv_file', required = FALSE, default = NULL, 
    help = "Path to CSV containing cv metrics for test sets.")
parser$add_argument('--order_folder', required = TRUE, 
    help = "Path to folder containing the order files for each test set.")

args <- parser$parse_args()

metrics_file <- args$metrics_file
figures_dir <- args$figures_dir
metadata_file <- args$metadata_file
cv_file <- args$cv_file
order_folder <- args$order_folder

# Ensure figures directory exists
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load Data ---
all_metrics <- read_csv(metrics_file, show_col_types = FALSE)%>%
  mutate(
    n_studies = as.numeric(str_extract(subset_file, "(\\d+)(?=studies)"))
  )

# Treat log_transformed as baseline
df_unadj <- all_metrics %>% filter(adjuster == "log_transformed")

all_adjusters <- all_metrics %>% filter(adjuster != "log_transformed")

# --- Helper Functions ---

# Calculate MCC
calculate_mcc <- function(tp, tn, fp, fn) {
  numerator <- (tp * tn) - (fp * fn)
  denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  ifelse(denominator == 0, 0, numerator / denominator)
}

# Prepare delta metrics vs unadjusted
prepare_delta <- function(df_adj, df_unadj, metric_col) {
  adj <- df_adj %>%
    group_by(adjuster, n_studies, test_source) %>%
    summarise(Adj = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop")
  
  unadj <- df_unadj %>%
    group_by(n_studies, test_source) %>%
    summarise(Unadj = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop")
  
  full_join(adj, unadj, by = c("n_studies", "test_source")) %>%
    mutate(Mean_Metric = Adj - Unadj)
}

# given a test source and n_studies, return the ordered training source
get_train_source <- function(test_source, n_studies, order_folder) {
  order_file <- file.path(order_folder, paste0(test_source, "_order.csv"))

  if (!file.exists(order_file)) {
    stop(paste("Order file not found for test source:", test_source))
  }

  order_df <- read_csv(order_file, show_col_types = FALSE)

  if (!"train_source" %in% colnames(order_df)) {
    stop(paste("train_source column missing in", order_file))
  }

  if (n_studies == 0) {
    return(NA_character_)
  }

  if (n_studies > nrow(order_df)) {
    stop(paste(
      "n_studies =", n_studies,
      "exceeds rows in", order_file
    ))
  }

  order_df$train_source[n_studies]
}

prepare_training_order_df <- function(metrics_df, order_folder) {
  order_cache <- list()

  get_train_cached <- function(test_source, n_studies) {
    if (n_studies == 0) return(NA_character_)

    if (!test_source %in% names(order_cache)) {
      order_file <- file.path(order_folder, paste0(test_source, "_order.csv"))
      order_cache[[test_source]] <<- read_csv(order_file, show_col_types = FALSE)
    }

    order_cache[[test_source]]$train_source[n_studies]
  }

  metrics_df %>%
    mutate(
      train_source = mapply(get_train_cached, test_source, n_studies)
      ) %>%
    filter(!is.na(train_source)) %>%
    group_by(test_source) %>%
    mutate(
      train_source = factor(
        train_source,
        levels = unique(train_source[order(n_studies)])
      )
    ) %>%
    ungroup()
}

plot_scaling_performance <- function(
  metrics_df,
  order_folder,
  metric_col, 
  figures_dir,
  filename,
  y_label = NULL,
  cv_value
) {

  plot_df <- prepare_training_order_df(metrics_df, order_folder)

  if (is.null(y_label)) {
    y_label <- metric_col
  }

  p <- ggplot(
    plot_df, 
    aes(
      x = train_source,
      y = .data[[metric_col]],
      color = adjuster,
      group = interaction(adjuster, test_source)
    )
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~ test_source, scales = "free_x") +
    labs(
      x = "Training dataset added",
      y = y_label, 
      color = "Adjuster"
    ) +
    theme_bw() +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

    # Add cross-validation line if provided
    if (!is.null(cv_value)) {
      p <- p + geom_hline(yintercept = cv_value, linetype = "dashed", color = "black") +
        annotate("text", x=1, y = cv_value, label = paste0("cv = ", cv_value),
                  vjust = -0.5, hjust = 0)
  }
  
  ggsave(
    filename = file.path(figures_dir, filename),
    plot = p,
    width = 12,
    height = 8
  )
  
  return(p)
}

# generate_test_source_scaling_plot_absolute <- function(
#     df,
#     metric = "ROC_AUC",
#     gse_metadata_path = metadata_file,
#     order_folder = order_folder,
#     cv_value = NULL,
#     fig_dir = figures_dir
# ) {

#   message("Generating ABSOLUTE test-source scaling plot for metric: ", metric)

#   # Create study labels 
#   study_labels <- setNames(LETTERS[seq_along(unique(df$test_source))],
#                       unique(df$test_source))

#   # Compute mean & SE across replicates
#   df_sum <- df %>%
#     group_by(test_source, test_source_label, adjuster, n_studies, technology, sample_size) %>%
#     summarize(
#       mean_val = mean(Value, na.rm = TRUE),
#       se_val   = sd(Value, na.rm = TRUE) / sqrt(n()),
#       .groups = "drop"
#     )

#   p <- ggplot(df_sum, aes(
#     x = n_studies,
#     y = mean_val,
#     color = adjuster, 
#     shape = technology,
#     size = sample_size,
#     group = adjuster
#   )) +
#     geom_point(position = position_dodge(width = 0.5)) +
#     geom_errorbar(aes(
#       ymin = mean_val - se_val,
#       ymax = mean_val + se_val
#     ), width = 0.25, position = position_dodge(width = 0.5)) +
#     geom_line(aes(group = adjuster), linewidth = 1, position = position_dodge(width = 0.5)) +
#     facet_wrap(~ test_source) +
#     theme_minimal(base_size = 14) +
#     theme(
#       panel.grid.minor = element_blank(),
#       axis.text.x = element_text(angle = 45, hjust = 1)
#     ) +
#     scale_size_continuous(name = "Sample Size", range = c(2, 6)) +
#     labs(
#       title = paste0("Absolute Scaling Performance Across Training Sizes (", metric, ")"),
#       x = "Test Study",
#       y = metric,
#       color = "Adjuster",
#       shape = "Technology"
#     )

#   # Add cross-validation line if provided
#   if (!is.null(cv_value)) {
#     p <- p + geom_hline(yintercept = cv_value, linetype = "dashed", color = "black") +
#       annotate("text", x=1, y = cv_value, label = paste0("cv = ", cv_value),
#                 vjust = -0.5, hjust = 0, linewidth = 4)
#   }

#   file_path <- file.path(fig_dir, paste0("absolute_scaling_", metric, "_enhanced.png"))
#   ggsave(file_path, p, width = 12, height = 7)

#   message("Saved: ", file_path)
# }

get_cv_value <- function(cv_data, test_source, metric) {
  if (is.null(cv_data)) return(NULL)
  
  cv_row <- cv_data %>% filter(test_source == !!test_source)
  
  if (nrow(cv_row) == 0) {
    warning("Test source not found in cv_file: ", test_source)
    return(NULL)
  }
  
  if (!metric %in% colnames(cv_row)) {
    warning("Metric not found in cv_file: ", metric)
    return(NULL)
  }
  
  cv_value <- cv_row[[metric]][1]  # take first match
  return(cv_value)
}

# --- Main Processing ---

# Cross-Validation line
if (!is.null(cv_file)) {
  cv_data <- read_csv(args$cv_file, col_types = cols())
} else {
  cv_data <- NULL
} 

absolute_data <- all_metrics %>%
  select(adjuster, n_studies, test_source, `ROC_AUC`, MCC) %>%
  pivot_longer(cols = c(`ROC_AUC`, MCC), names_to = "Metric", values_to = "Value")

# Compute delta metrics for ROC_AUC and MCC
results <- data.frame()
for (metric in c("ROC_AUC", "MCC")) {
  # Compute delta metrics
  delta_df <- prepare_delta(all_adjusters, df_unadj, metric) %>%
    mutate(Metric = metric)
  
  results <- bind_rows(results, delta_df)
  
  # Compute absolute data for the metric
  abs_data <- absolute_data %>% filter(Metric == metric)
  
  # Example: use first test_source in the data for CV line
  first_test_source <- unique(abs_data$test_source)[1]
  
  cv_value <- get_cv_value(cv_data, test_source = first_test_source, metric = metric)

  # Generate other scaling plot DEBUG
  # generate_absolute_scaling_plot(
  #   df = abs_data, 
  #   order_folder = order_folder, 
  #   metric = metric, 
  #   gse_metadata_path = metadata_file, 
  #   fig_dir = figures_dir
  # )

  plot_metrics <- all_metrics %>% 
    group_by(adjuster, n_studies, test_source) %>%
    summarise(
      ROC_AUC = mean(ROC_AUC),
      MCC = mean(MCC), 
      .groups = "drop"
    )

  plot_scaling_performance(
    metrics_df = plot_metrics, 
    order_folder = order_folder, 
    metric_col = metric, 
    figures_dir = figures_dir, 
    filename = paste0("absolute_scaling_", metric, ".png"), 
    y_label = metric, cv_value = cv_value
  )
}

# Save combined CSV
write_csv(results, file.path(figures_dir, "scaling_comparison_results.csv"))


cat("✅ Scaling plots generated successfully in:", figures_dir, "\n")
