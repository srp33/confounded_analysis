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
generate_scaling_plot <- function(all_data, metric = "ROC_AUC", fig_dir = CONFIG$figures_dir) {
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
generate_test_source_scaling_plot <- function(all_data, metric = "ROC_AUC", fig_dir = CONFIG$figures_dir) {
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

generate_test_source_scaling_plot_absolute <- function(
    df,
    metric = "ROC_AUC",
    fig_dir = CONFIG$figures_dir
) {

  message("Generating ABSOLUTE test-source scaling plot for metric: ", metric)

  # Filter for desired metric
  df <- df %>% filter(Metric == metric)

  # Preserve test_source order in the raw data file
  test_source_order <- df %>%
    distinct(test_source) %>%
    pull(test_source)

  df <- df %>%
    mutate(
      test_source = factor(test_source, levels = test_source_order),
      adjuster = factor(adjuster),
      test_source_label = paste0("Test study: ", toupper(test_source))
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
    color = adjuster
  )) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(
      ymin = mean_val - se_val,
      ymax = mean_val + se_val
    ), width = 0.25) +
    facet_wrap(~ test_source_label, scales = "fixed") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = paste0("Absolute Scaling Performance Across Training Sizes (", metric, ")"),
      x = "Number of Training Studies",
      y = metric,
      color = "Adjuster"
    )

  file_path <- file.path(fig_dir, paste0("absolute_scaling_", metric, ".png"))
  ggsave(file_path, p, width = 12, height = 7)

  message("Saved: ", file_path)
}

generate_absolute_boxplot <- function(all_data, metric_name = "ROC_AUC", fig_dir = CONFIG$figures_dir) {
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

generate_adjuster_comparison_by_testset <- function(
    df,
    metric = "ROC_AUC",
    fig_dir = CONFIG$figures_dir
) {

  message("Generating adjuster comparison plots for metric: ", metric)

  # Filter to desired metric
  df <- df %>% filter(Metric == metric)

  # Preserve test_source order AS APPEARS in the data file
  test_source_order <- df %>%
    distinct(test_source) %>%
    pull(test_source)

  df <- df %>%
    mutate(
      test_source = factor(test_source, levels = test_source_order),
      adjuster = factor(adjuster)
    )

  # Loop over number of studies and make 1 figure per setting
  for (n in sort(unique(df$n_studies))) {

    df_n <- df %>% filter(n_studies == n)

    if (nrow(df_n) == 0) next

    p <- ggplot(df_n, aes(
      x = adjuster,
      y = Value,
      fill = adjuster
    )) +
      geom_boxplot(outlier.shape = NA, alpha = 0.8) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 1.8) +
      facet_wrap(~ test_source, scales = "free_y") +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = paste0("Adjuster Performance Across Test Sets – ", n, " Studies"),
        subtitle = paste("Metric:", metric),
        x = "Adjustment Method",
        y = metric,
        fill = "Adjuster"
      )

    # Save figure
    out_file <- file.path(fig_dir, paste0(
      "adjusters_by_testset_",
      metric, "_",
      n, "_studies.png"
    ))

    ggsave(out_file, p, width = 14, height = 8)
    message("Saved: ", out_file)
  }
}

generate_adjuster_lineplots_by_testset <- function(
    df,
    metric = "ROC_AUC",
    fig_dir = CONFIG$figures_dir
) {

  message("Generating mean ± SE line plots for metric: ", metric)

  # Filter metric
  df <- df %>% filter(Metric == metric)

  # Preserve test_source order in the dataset
  test_source_order <- df %>%
    distinct(test_source) %>%
    pull(test_source)

  # Standardize adjuster order globally
  adjuster_order <- df %>%
    distinct(adjuster) %>%
    pull(adjuster)

  df <- df %>%
    mutate(
      test_source = factor(test_source, levels = test_source_order),
      adjuster = factor(adjuster, levels = adjuster_order)
    )

  # Summary: mean + SE per (adjuster, test_source, n_studies)
  df_sum <- df %>%
    group_by(n_studies, test_source, adjuster) %>%
    summarize(
      mean_val = mean(Value, na.rm = TRUE),
      se_val   = sd(Value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      test_source_label = paste0("Test study: ", toupper(test_source))
    )

  # Loop through number of studies
  for (n in sort(unique(df_sum$n_studies))) {

    df_n <- df_sum %>% filter(n_studies == n)
    if (nrow(df_n) == 0) next

    p <- ggplot(df_n, aes(
      x = adjuster,
      y = mean_val,
      group = test_source_label,
      color = adjuster
    )) +
      geom_line(linewidth = 1.1, alpha = 0.8) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = mean_val - se_val,
                        ymax = mean_val + se_val),
                    width = 0.15) +
      facet_wrap(~ test_source_label, scales = "fixed") +  # standardized x-axis
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = "none"
      ) +
      labs(
        title = paste0("Performance Across Adjusters – ", n, " Studies"),
        subtitle = paste("Metric:", metric),
        x = "Adjustment Method",
        y = metric
      )

    # Output file
    out_file <- file.path(
      fig_dir,
      paste0("adjuster_lineplots_", metric, "_", n, "_studies.png")
    )

    ggsave(out_file, p, width = 13, height = 9)
    message("Saved: ", out_file)
  }
}


# --- Main Processing ---

# Compute delta metrics for ROC_AUC and MCC
results <- data.frame()
for (metric in c("ROC_AUC", "MCC")) {
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
  select(adjuster, n_studies, test_source, `ROC_AUC`, MCC) %>%
  pivot_longer(cols = c(`ROC_AUC`, MCC), names_to = "Metric", values_to = "Value")

# Save combined CSV
write_csv(results, file.path(CONFIG$figures_dir, "scaling_comparison_results.csv"))

# Generate plots
generate_test_source_scaling_plot(results %>% filter(Metric == "ROC_AUC"), "ROC_AUC")
generate_test_source_scaling_plot(results %>% filter(Metric == "MCC"), "MCC")

generate_test_source_scaling_plot_absolute(absolute_data, "ROC_AUC")
generate_test_source_scaling_plot_absolute(absolute_data, "MCC")

generate_absolute_boxplot(absolute_data, "ROC_AUC")
generate_absolute_boxplot(absolute_data, "MCC")

generate_adjuster_lineplots_by_testset(absolute_data, "ROC_AUC")
generate_adjuster_lineplots_by_testset(absolute_data, "MCC")


cat("✅ Scaling plots generated successfully in:", CONFIG$figures_dir, "\n")
