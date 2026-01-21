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
library(tidytext)

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
  ) %>%
  mutate(test_source = tolower(test_source))

gse_metadata <- read_csv(metadata_file, show_col_types = FALSE) %>%
  mutate(
    gse_id = tolower(trimws(gse_id)),
    technology = factor(technology)
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

prepare_training_order_df <- function(metrics_df, order_folder, gse_metadata) {
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
      train_source = mapply(get_train_cached, test_source, n_studies),
      train_order = n_studies
    ) %>%
    filter(!is.na(train_source)) %>%
    left_join(
      gse_metadata,
      by = c("train_source" = "gse_id")
    )  %>%
    group_by(test_source) %>%
    mutate(train_order = n_studies) %>%
    ungroup()
}

plot_scaling_performance <- function(
  metrics_df,
  order_folder,
  metric_col, 
  gse_metadata,
  figures_dir,
  filename,
  y_label = NULL,
  cv_value = NULL
) {

  plot_df <- prepare_training_order_df(
    metrics_df,
    order_folder,
    gse_metadata
  )

  plot_df <- plot_df %>%
    group_by(test_source, adjuster, n_studies, train_source, train_order) %>%
    summarize(
      mean_val = mean(.data[[metric_col]], na.rm = TRUE),
      se_val   = sd(.data[[metric_col]], na.rm = TRUE) / sqrt(n()),
      technology = first(technology),
      .groups = "drop"
    )

  if (is.null(y_label)) {
    y_label <- metric_col
  }

  all_train_sources <- plot_df %>%
    arrange(train_order) %>%
    pull(train_source) %>%
    unique()

  train_label_map <- setNames(LETTERS[seq_along(all_train_sources)], all_train_sources)

  plot_df <- plot_df %>%
    mutate(train_label = train_label_map[train_source])

  legend_df <- data.frame(
    train_label = names(train_label_map),
    letter = train_label_map, 
    gse_id = names(train_label_map)
  )

  facet_labels <- train_label_map[names(train_label_map) %in% unique(plot_df$test_source)]
  facet_labels <- setNames(
    paste0(names(facet_labels), " (", facet_labels, ")"),
    names(facet_labels)
  )

  p <- ggplot(
    plot_df, 
    aes(
      x = reorder_within(train_source, train_order, test_source),
      y = mean_val,
      color = adjuster,
      shape = technology,
      group = adjuster
    )
  ) +
    scale_x_reordered(
      labels = function(x) {
        # Remove the __<facet> part added by reorder_within
        gsub("__.*$", "", x) %>%
          { train_label_map[.] }
      }
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_errorbar(
      aes(
        ymin = mean_val - se_val,
        ymax = mean_val + se_val
      ),
      width = 0.2
    ) +
    facet_wrap(~ test_source, scales = "free_x",
      labeller = labeller(test_source = facet_labels)) +
    labs(
      x = "Training dataset added (ordered)",
      y = y_label, 
      color = "Adjuster",
      shape = "Technology"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(hjust = 1)
    )

  p <- p + geom_blank(
    data = legend_df,
    aes(x = 0, y = 0, fill = letter),
    inherit.aes = FALSE
    ) + 
    scale_fill_manual(
      name = "Training Dataset",
      values = setNames(rep("black", length(train_label_map)), train_label_map),
      labels = paste0(legend_df$letter, " = ", legend_df$gse_id)
    ) +
    guides(
      fill = guide_legend(override.aes = list(shape = NA))
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
    gse_metadata = gse_metadata,
    figures_dir = figures_dir, 
    filename = paste0("new_absolute_scaling_", metric, ".png"), 
    y_label = metric, 
    cv_value = cv_value
  )
}

# Save combined CSV
write_csv(results, file.path(figures_dir, "scaling_comparison_results.csv"))


cat("✅ Scaling plots generated successfully in:", figures_dir, "\n")
