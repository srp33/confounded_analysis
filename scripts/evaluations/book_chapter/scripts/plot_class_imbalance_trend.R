#!/usr/bin/env Rscript

# plot_class_imbalance_trend.R - Plot adjuster performance trends across class imbalance levels
# Aggregation strategy:
# 1. Average MCC over training pairs and replicates -> (classifier, adjuster, test_set, imbalance)
# 2. Rank adjusters within (classifier, test_set, imbalance)
# 3. Average ranks over classifiers -> (adjuster, test_set, imbalance)
# 4. Plot mean rank over test sets with min/max uncertainty zone

library(argparse)
library(dplyr)
library(ggplot2)
library(tidyr)

# Parse command line arguments
parser <- ArgumentParser(description = "Create class imbalance trend plot with aggregated rankings")
parser$add_argument("--input-data", required = TRUE, help = "Path to class imbalanced results CSV")
parser$add_argument("-o", "--output", required = TRUE, help = "Output PNG file path")

args <- parser$parse_args()

# Load data
cat("Loading class imbalanced results from:", args$input_data, "\n")
data <- read.csv(args$input_data, stringsAsFactors = FALSE)

cat("Data dimensions:", nrow(data), "rows\n")
cat("Adjusters:", paste(unique(data$adjuster), collapse = ", "), "\n")
cat("Classifiers:", paste(unique(data$classifier), collapse = ", "), "\n")
cat("Imbalance levels:", paste(sort(unique(data$imbalance_pct)), collapse = ", "), "\n")

# Filter to target adjusters
target_adjusters <- c("unadjusted", "combat", "combat_sup")
data_filtered <- data %>%
  filter(adjuster %in% target_adjusters)

cat("Filtered to", nrow(data_filtered), "rows with target adjusters\n")

# Step 1: Average MCC over training pairs and replicates
# Result: one value per (classifier, adjuster, test_set, imbalance)
cat("\nStep 1: Averaging MCC over training pairs and replicates...\n")
mcc_averaged <- data_filtered %>%
  group_by(classifier, adjuster, test_dataset, imbalance_pct) %>%
  summarise(
    mean_mcc = mean(mcc, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

cat("  Result:", nrow(mcc_averaged), "unique (classifier, adjuster, test_set, imbalance) combinations\n")

# Step 2: Rank adjusters within each (classifier, test_set, imbalance) group
cat("\nStep 2: Ranking adjusters within (classifier, test_set, imbalance)...\n")
ranked_by_classifier <- mcc_averaged %>%
  group_by(classifier, test_dataset, imbalance_pct) %>%
  arrange(desc(mean_mcc)) %>%
  mutate(rank = rank(-mean_mcc, ties.method = "average")) %>%
  ungroup()

cat("  Ranks assigned\n")

# Step 3: Average ranks over classifiers
# Result: one value per (adjuster, test_set, imbalance)
cat("\nStep 3: Averaging ranks over classifiers...\n")
ranks_averaged_over_classifiers <- ranked_by_classifier %>%
  group_by(adjuster, test_dataset, imbalance_pct) %>%
  summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    n_classifiers = n(),
    .groups = "drop"
  )

cat("  Result:", nrow(ranks_averaged_over_classifiers), "unique (adjuster, test_set, imbalance) combinations\n")

# Step 4: Calculate mean, min, max over test sets for plotting
cat("\nStep 4: Computing statistics over test sets...\n")
plot_data <- ranks_averaged_over_classifiers %>%
  group_by(adjuster, imbalance_pct) %>%
  summarise(
    mean_rank = mean(mean_rank, na.rm = TRUE),
    min_rank = min(mean_rank, na.rm = TRUE),
    max_rank = max(mean_rank, na.rm = TRUE),
    n_test_sets = n(),
    .groups = "drop"
  )

cat("  Final plot data:", nrow(plot_data), "points\n")

# Create adjuster labels for plotting
plot_data <- plot_data %>%
  mutate(
    adjuster_label = case_when(
      adjuster == "unadjusted" ~ "Unadjusted",
      adjuster == "combat" ~ "ComBat (Unsupervised)",
      adjuster == "combat_sup" ~ "ComBat (Supervised)",
      TRUE ~ adjuster
    ),
    adjuster_label = factor(adjuster_label, 
                           levels = c("Unadjusted", "ComBat (Unsupervised)", "ComBat (Supervised)")),
    imbalance_pct_num = imbalance_pct * 100
  )

# Print summary
cat("\nSummary of plot data:\n")
print(plot_data %>% select(adjuster_label, imbalance_pct_num, mean_rank, min_rank, max_rank))

# Create the plot
p <- ggplot(plot_data, aes(x = imbalance_pct_num, y = mean_rank, color = adjuster_label, fill = adjuster_label)) +
  # Uncertainty ribbon (min/max over test sets)
  geom_ribbon(aes(ymin = min_rank, ymax = max_rank), alpha = 0.2, color = NA) +
  # Mean line
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  # Reverse y-axis so rank 1 is at top
  scale_y_reverse(
    breaks = c(1, 1.5, 2, 2.5, 3),
    labels = c("1st", "1.5", "2nd", "2.5", "3rd"),
    limits = c(3, 1)
  ) +
  scale_x_continuous(
    breaks = c(20, 30, 40, 50),
    labels = c("20%", "30%", "40%", "50%")
  ) +
  scale_color_manual(
    values = c(
      "Unadjusted" = "#E31A1C",
      "ComBat (Unsupervised)" = "#1F78B4",
      "ComBat (Supervised)" = "#33A02C"
    ),
    name = "Batch Adjuster"
  ) +
  scale_fill_manual(
    values = c(
      "Unadjusted" = "#E31A1C",
      "ComBat (Unsupervised)" = "#1F78B4",
      "ComBat (Supervised)" = "#33A02C"
    ),
    name = "Batch Adjuster"
  ) +
  labs(
    title = "Batch Correction Performance Across Class Imbalance Levels",
    subtitle = "Ranks averaged over classifiers | Shaded region shows min-max range across test sets",
    x = "Class Imbalance Level (% Active TB in High-Imbalance Training Set)",
    y = "Average Performance Rank",
    caption = "Lower rank = better performance | Each point represents mean over test sets"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

# Save the plot
ggsave(args$output, p, width = 10, height = 7, dpi = 300)

cat("\nPlot saved to:", args$output, "\n")

# Calculate and print trend statistics
cat("\nTrend analysis (linear regression of rank vs imbalance):\n")
trend_stats <- plot_data %>%
  group_by(adjuster_label) %>%
  do({
    model <- lm(mean_rank ~ imbalance_pct_num, data = .)
    data.frame(
      slope = coef(model)[2],
      intercept = coef(model)[1],
      p_value = summary(model)$coefficients[2, 4],
      r_squared = summary(model)$r.squared
    )
  }) %>%
  ungroup()

print(trend_stats)

cat("\nInterpretation:\n")
cat("- Positive slope: Performance degrades (higher rank) with increasing imbalance\n")
cat("- Negative slope: Performance improves (lower rank) with increasing imbalance\n")
cat("- p < 0.05: Statistically significant trend\n")

# Calculate rank stability (how much variation across test sets)
cat("\nRank stability across test sets (mean range):\n")
stability <- plot_data %>%
  mutate(range = max_rank - min_rank) %>%
  group_by(adjuster_label) %>%
  summarise(
    mean_range = mean(range),
    max_range = max(range),
    .groups = "drop"
  )

print(stability)
cat("Lower range = more consistent performance across different test sets\n")
