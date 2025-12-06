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

parser$add_argument('--aggregated_metrics_file', required = TRUE,
    help = "Path to CSV containing aggregated metrics")
parser$add_argument('--figures_dir', required = TRUE, 
    help = "Directory to save figures")
parser$add_argument('--metadata_file', required = FALSE,
    help = "Path to CSV containing dataset metadata")

args <- parser$parse_args()

metrics_file <- args$aggregated_metrics_file
figures_dir <- args$figures_dir
metadata_file <- args$metadata_file

# Ensure figures directory exists
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load Data ---
all_metrics <- read_csv(metrics_file, show_col_types = FALSE)%>%
  mutate(
    n_studies = as.numeric(str_extract(subset_file, "(?<=subset)\\d+(?=studies)")),
    n_studies = ifelse(is.na(n_studies), 0, n_studies) # default to 0 if missing
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

# Generate scaling plots
generate_scaling_plot <- function(all_data, metric = "ROC_AUC", fig_dir = figures_dir) {
  if (nrow(all_data) == 0) return()
  
  p <- ggplot(all_data, aes(x = adjuster, y = Mean_Metric, fill = adjuster)) +
    geom_boxplot() +
    facet_wrap(~ n_studies, scales = "free_y", 
      labeller = labeller(n_studies = function(x) paste0(x, " Studies"))) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_blank(),   # remove x-axis text
      axis.ticks.x = element_blank()   # remove x-axis ticks
    ) +
    labs(
      title = paste0("Distribution of Delta Performance (", metric, ")"),
      x = "Adjuster",
      y = paste0("Delta vs log_transformed"),
      fill = "Adjuster"
    )
  
  file_path <- file.path(fig_dir, paste0("boxplot_", metric, "_by_nstudies.png"))
  ggsave(file_path, p, width = 12, height = 6)
  cat("Saved: ", file_path, "\n")
}

# Generate scaling plots
generate_test_source_scaling_plot <- function(all_data, metric = "ROC_AUC", fig_dir = figures_dir) {
  if (nrow(all_data) == 0) return()
  
  p <- ggplot(all_data, aes(x = n_studies, y = Mean_Metric, color = adjuster)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ test_source, scales = "free_y") +
    theme_minimal(base_size = 13) +
    labs(
      title = paste0("Scaling Performance (", metric, ")"),
      x = "Number of Training Studies",
      y = metric,
      color = "Adjuster"
    )
  
  file_path <- file.path(fig_dir, paste0("scaling_", metric, ".png"))
  ggsave(file_path, p, width = 10, height = 6)
  cat("Saved: ", file_path, "\n")
}

# Absolute scaling plot -- this is the one I am using
generate_test_source_scaling_plot_absolute <- function(
    df,
    metric = "ROC_AUC",
    gse_metadata_path=metadata_file,
    cv_value = NULL,
    fig_dir = figures_dir
) {

  message("Generating ABSOLUTE test-source scaling plot for metric: ", metric)

  # Load metadata
  gse_meta <- read_csv(gse_metadata_path, col_types = cols()) %>%
    select(gse_id, technology)

  # Filter for desired metric
  df <- df %>% filter(Metric == metric)

  # Merge sequencing technology info
  df <- df %>%
    left_join(gse_meta, by = c("test_source" = "gse_id"))

  # Create study labels 
  study_labels <- setNames(LETTERS[seq_along(unique(df$test_source))],
                      unique(df$test_source))

  df <- df %>%
    mutate(
      test_source_label = study_labels[test_source],
      adjuster = factor(adjuster),
      technology = factor(technology, levels = c("microarray", "rna-seq"))
    )

  # Compute mean & SE across replicates
  df_sum <- df %>%
    group_by(test_source, test_source_label, adjuster, n_studies) %>%
    summarize(
      mean_val = mean(Value, na.rm = TRUE),
      se_val   = sd(Value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  p <- ggplot(df_sum, aes(
    x = n_studies,
    y = mean_val,
    color = adjuster, 
    shape = technology,
    size = sample_size,
    group = adjuster
  )) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(
      ymin = mean_val - se_val,
      ymax = mean_val + se_val
    ), width = 0.25, position = position_dodge(width = 0.5)) +
    geom_line(aes(group = adjuster), linewidth = 1, position = position_dodge(width = 0.5)) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_size_continuous(name = "Sample Size") +
    labs(
      title = paste0("Absolute Scaling Performance Across Training Sizes (", metric, ")"),
      x = "Test Study",
      y = metric,
      color = "Adjuster",
      shape = "Technology"
    )

  # Add cross-validation line if provided
  if (!is.null(cv_value)) {
    p <- p + geom_hline(y_intercept = cv_value, linetype = "dashed", color = "black") +
      annotate("text", x=1, y = cv_value, label = paste0("cv = ", cv_value),
                vjust = -0.5, hjust = 0, size = 4)
  }

  file_path <- file.path(fig_dir, paste0("absolute_scaling_", metric, "_enhanced.png"))
  ggsave(file_path, p, width = 12, height = 7)

  message("Saved: ", file_path)
}


# --- Main Processing ---

# Compute delta metrics for ROC_AUC and MCC
results <- data.frame()
for (metric in c("ROC_AUC", "MCC")) {
  delta_df <- prepare_delta(all_adjusters, df_unadj, metric) %>%
    mutate(Metric = metric)
  results <- bind_rows(results, delta_df)
}

absolute_data <- all_metrics %>%
  select(adjuster, n_studies, test_source, `ROC_AUC`, MCC) %>%
  pivot_longer(cols = c(`ROC_AUC`, MCC), names_to = "Metric", values_to = "Value")

# Save combined CSV
write_csv(results, file.path(figures_dir, "scaling_comparison_results.csv"))

# Generate plots
generate_test_source_scaling_plot(results %>% filter(Metric == "ROC_AUC"), "ROC_AUC")
generate_test_source_scaling_plot(results %>% filter(Metric == "MCC"), "MCC")

generate_test_source_scaling_plot_absolute(absolute_data, "ROC_AUC")
generate_test_source_scaling_plot_absolute(absolute_data, "MCC")


cat("✅ Scaling plots generated successfully in:", figures_dir, "\n")
