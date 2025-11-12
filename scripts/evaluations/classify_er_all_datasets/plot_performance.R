# plot_performance.R
#
# This script combines multiple adjuster CSVs into a single dataframe,
# calculates MCC, computes deltas vs unadjusted, and generates scaling plots.

# --- Load Libraries ---
library(ggplot2)
library(readr)
library(dplyr)
library(stringr)

# --- Configuration ---
CONFIG <- list(
  metrics_dir = "/outputs/classify_er_all",   # <-- change this to your folder
  figures_dir = "/outputs/figures",
  unadjusted_file = "unadjusted.csv"      # <-- put your baseline CSV here
)

# Ensure figures directory exists
dir.create(CONFIG$figures_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper Functions ---

# Calculate MCC
calculate_mcc <- function(tp, tn, fp, fn) {
  numerator <- (tp * tn) - (fp * fn)
  denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  ifelse(denominator == 0, 0, numerator / denominator)
}

# Read a CSV and add Adjuster column
read_adjuster_csv <- function(file_path) {
  adjuster <- str_remove(basename(file_path), "\\.csv$")
  df <- read_csv(file_path, show_col_types = FALSE)
  if (nrow(df) == 0) return(data.frame())
  
  df <- df %>%
    mutate(
      Adjuster = adjuster,
      n_studies = as.numeric(str_extract(subset_file, "(?<=subset)\\d+")),
      n_studies = ifelse(is.na(n_studies), 0, n_studies),
      MCC = calculate_mcc(`True Positive`, `True Negative`, `False Positive`, `False Negative`)
    )
  return(df)
}

# Prepare delta metrics vs unadjusted
prepare_delta <- function(df_adj, df_unadj, metric_col) {
  adj <- df_adj %>%
    group_by(Adjuster, n_studies, test_source) %>%
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
  
  p <- ggplot(all_data, aes(x = n_studies, y = Mean_Metric, color = Adjuster)) +
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

# --- Main Processing ---

# Load unadjusted baseline
file_unadjusted <- file.path(CONFIG$metrics_dir, CONFIG$unadjusted_file)
df_unadj <- read_adjuster_csv(file_unadjusted) %>% mutate(Adjuster = "unadjusted")

# Load all adjuster CSVs
adjuster_files <- list.files(CONFIG$metrics_dir, pattern = "\\.csv$", full.names = TRUE)
adjuster_files <- adjuster_files[!basename(adjuster_files) %in% CONFIG$unadjusted_file]

all_adjusters <- lapply(adjuster_files, read_adjuster_csv) %>% bind_rows()

# Compute delta metrics for ROC AUC and MCC
results <- data.frame()
for (metric in c("ROC AUC", "MCC")) {
  delta_df <- prepare_delta(all_adjusters, df_unadj, metric) %>%
    mutate(Metric = metric)
  results <- bind_rows(results, delta_df)
}

# Save combined CSV
write_csv(results, file.path(CONFIG$metrics_dir, "scaling_comparison_results.csv"))

# Generate plots
generate_scaling_plot(results %>% filter(Metric == "ROC AUC"), "ROC AUC")
generate_scaling_plot(results %>% filter(Metric == "MCC"), "MCC")

cat("✅ Scaling plots generated successfully in:", CONFIG$figures_dir, "\n")
