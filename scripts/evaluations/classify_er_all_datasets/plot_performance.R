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

# --- Configuration ---
CONFIG <- list(
  aggregated_metrics_file = "/outputs/classify_er_all/all_metrics.csv",
  figures_dir = "/outputs/classify_all_figures"
)

# Read the aggregated CSV
all_metrics <- read_csv(CONFIG$aggregated_metrics_file, show_col_types = FALSE)%>%
  mutate(
    n_studies = as.numeric(str_extract(subset_file, "(?<=subset)\\d+(?=studies)")),
    n_studies = ifelse(is.na(n_studies), 0, n_studies) # default to 0 if missing
  )

# Treat log_transformed as baseline
df_unadj <- all_metrics %>% filter(adjuster == "log_transformed")

all_adjusters <- all_metrics %>% filter(adjuster != "log_transformed")

# Ensure figures directory exists
dir.create(CONFIG$figures_dir, showWarnings = FALSE, recursive = TRUE)

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
generate_scaling_plot <- function(all_data, metric = "ROC AUC", fig_dir = CONFIG$figures_dir) {
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

generate_absolute_boxplot <- function(all_data, metric_name = "ROC AUC", fig_dir = CONFIG$figures_dir) {
  # Filter for the metric we want
  data_plot <- all_data %>% filter(Metric == metric_name)
  
  if (nrow(data_plot) == 0) return()
  
  p <- ggplot(data_plot, aes(x = adjuster, y = Value, fill = adjuster)) +
    geom_boxplot() +
    facet_wrap(~ n_studies, scales = "free_y", 
      labeller = labeller(n_studies = function(x) paste0(x, " Studies"))) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_blank(),   # remove x-axis text
      axis.ticks.x = element_blank()   # remove x-axis ticks
    ) +
    labs(
      title = paste0("Absolute Performance (", metric_name, ")"),
      x = "Adjuster",
      y = metric_name,
      fill = "Adjuster"
    )
  
  file_path <- file.path(fig_dir, paste0("absolute_boxplot_", metric_name, "_by_nstudies.png"))
  ggsave(file_path, p, width = 12, height = 6)
  cat("Saved: ", file_path, "\n")
}


# --- Main Processing ---

# Compute delta metrics for ROC AUC and MCC
results <- data.frame()
for (metric in c("ROC AUC", "MCC")) {
  delta_df <- prepare_delta(all_adjusters, df_unadj, metric) %>%
    mutate(Metric = metric)
  results <- bind_rows(results, delta_df)
}

# results_avg <- results %>%
#   group_by(adjuster, n_studies, Metric) %>%
#   summarise(
#     Mean_Metric = mean(Mean_Metric, na.rm = TRUE),
#     .groups = "drop"
#   )

absolute_data <- all_metrics %>%
  select(adjuster, n_studies, test_source, `ROC AUC`, MCC) %>%
  pivot_longer(cols = c(`ROC AUC`, MCC), names_to = "Metric", values_to = "Value")

# Save combined CSV
write_csv(results, file.path(CONFIG$figures_dir, "scaling_comparison_results.csv"))

# Generate plots
generate_scaling_plot(results %>% filter(Metric == "ROC AUC"), "ROC AUC")
generate_scaling_plot(results %>% filter(Metric == "MCC"), "MCC")

generate_absolute_boxplot(absolute_data, "ROC AUC")
generate_absolute_boxplot(absolute_data, "MCC")


cat("✅ Scaling plots generated successfully in:", CONFIG$figures_dir, "\n")
