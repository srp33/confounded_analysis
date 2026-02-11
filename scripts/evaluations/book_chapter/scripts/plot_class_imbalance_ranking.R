#!/usr/bin/env Rscript

# plot_class_imbalance_ranking.R - Plot adjuster performance rankings across class imbalance levels
# Ranks adjusters within classifier/test_set/combination groups, then plots by imbalance level

library(argparse)
library(dplyr)
library(ggplot2)
library(tidyr)

# Parse command line arguments
parser <- ArgumentParser(description = "Create class imbalance ranking comparison plot")
parser$add_argument("--input-data", required = TRUE, help = "Path to class imbalanced results CSV")
parser$add_argument("-o", "--output", required = TRUE, help = "Output PNG file path")
parser$add_argument("--average-by", default = "none", 
                   help = "Averaging mode: 'none', 'test_sets', 'training_sets', 'classifiers', or 'combinations'")

args <- parser$parse_args()

# Load data
cat("Loading class imbalanced results from:", args$input_data, "\n")
data <- read.csv(args$input_data, stringsAsFactors = FALSE)

cat("Data dimensions:", nrow(data), "rows\n")
cat("Adjusters:", unique(data$adjuster), "\n")
cat("Classifiers:", unique(data$classifier), "\n")
cat("Imbalance levels:", sort(unique(data$imbalance_pct)), "\n")
if ("replicate" %in% colnames(data)) {
  cat("Replicates:", sort(unique(data$replicate)), "\n")
}

# Filter to only the adjusters we want to plot
target_adjusters <- c("unadjusted", "combat", "combat_sup")
data_filtered <- data %>%
  filter(adjuster %in% target_adjusters)

cat("Filtered to", nrow(data_filtered), "rows with target adjusters\n")

# Rank adjusters within each classifier/test_dataset/training_pair/imbalance_pct/replicate group
# This ensures fair comparison within each specific scenario
ranked_data <- data_filtered %>%
  group_by(classifier, test_dataset, training_pair, imbalance_pct, replicate) %>%
  arrange(desc(mcc)) %>%  # Higher MCC = better performance = lower rank number
  mutate(rank = rank(-mcc, ties.method = "average")) %>%  # Negative MCC so higher MCC gets rank 1
  ungroup()

cat("Ranking completed\n")

# Apply averaging based on user specification
if (args$average_by == "none") {
  # No averaging - treat each point as independent
  plot_data <- ranked_data %>%
    mutate(
      group_id = paste(classifier, test_dataset, training_pair, sep = "_"),
      point_id = row_number()
    )
  
} else if (args$average_by == "test_sets") {
  # Average ranks across test sets (within classifier/training_pair/imbalance)
  plot_data <- ranked_data %>%
    group_by(classifier, training_pair, imbalance_pct, adjuster) %>%
    summarise(
      rank = mean(rank, na.rm = TRUE),
      se_rank = sd(rank, na.rm = TRUE) / sqrt(n()),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    mutate(group_id = paste(classifier, training_pair, sep = "_"))
  
} else if (args$average_by == "training_sets") {
  # Average ranks across training pairs (within classifier/test_dataset/imbalance)
  plot_data <- ranked_data %>%
    group_by(classifier, test_dataset, imbalance_pct, adjuster) %>%
    summarise(
      rank = mean(rank, na.rm = TRUE),
      se_rank = sd(rank, na.rm = TRUE) / sqrt(n()),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    mutate(group_id = paste(classifier, test_dataset, sep = "_"))
  
} else if (args$average_by == "classifiers") {
  # Average ranks across classifiers (within test_dataset/training_pair/imbalance)
  plot_data <- ranked_data %>%
    group_by(test_dataset, training_pair, imbalance_pct, adjuster) %>%
    summarise(
      rank = mean(rank, na.rm = TRUE),
      se_rank = sd(rank, na.rm = TRUE) / sqrt(n()),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    mutate(group_id = paste(test_dataset, training_pair, sep = "_"))
  
} else if (args$average_by == "combinations") {
  # Average ranks across training pairs (within classifier/test_dataset/imbalance)
  plot_data <- ranked_data %>%
    group_by(classifier, test_dataset, imbalance_pct, adjuster) %>%
    summarise(
      rank = mean(rank, na.rm = TRUE),
      se_rank = sd(rank, na.rm = TRUE) / sqrt(n()),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    mutate(group_id = paste(classifier, test_dataset, sep = "_"))
  
} else {
  stop(sprintf("Invalid averaging mode: %s", args$average_by))
}

# Create adjuster labels for plotting
plot_data <- plot_data %>%
  mutate(
    adjuster_label = case_when(
      adjuster == "unadjusted" ~ "Unadjusted",
      adjuster == "combat" ~ "ComBat (Unsupervised)",
      adjuster == "combat_sup" ~ "ComBat (Supervised)",
      TRUE ~ adjuster
    ),
    adjuster_label = factor(adjuster_label, levels = c("Unadjusted", "ComBat (Unsupervised)", "ComBat (Supervised)")),
    imbalance_pct_label = sprintf("%.0f%%", imbalance_pct * 100)
  )

cat("Plot data prepared with", nrow(plot_data), "points\n")

# Create the base plot
if (args$average_by == "none") {
  # No error bars for individual points
  p <- ggplot(plot_data, aes(x = imbalance_pct * 100, y = rank, color = adjuster_label)) +
    geom_point(alpha = 0.6, size = 2, position = position_jitter(width = 0.5, height = 0)) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.3, linewidth = 1.2)
  
} else {
  # Include error bars for averaged data
  p <- ggplot(plot_data, aes(x = imbalance_pct * 100, y = rank, color = adjuster_label)) +
    geom_errorbar(aes(ymin = rank - se_rank, ymax = rank + se_rank), 
                  width = 1, alpha = 0.7, position = position_dodge(width = 1)) +
    geom_point(size = 3, position = position_dodge(width = 1)) +
    geom_line(aes(group = adjuster_label), linewidth = 1.2, position = position_dodge(width = 1))
}

# Complete the plot
p <- p +
  scale_y_reverse(breaks = c(1, 2, 3), labels = c("1st", "2nd", "3rd")) +
  scale_x_continuous(breaks = c(20, 30, 40, 50), labels = c("20%", "30%", "40%", "50%")) +
  scale_color_manual(
    values = c("Unadjusted" = "#E31A1C", "ComBat (Unsupervised)" = "#1F78B4", "ComBat (Supervised)" = "#33A02C"),
    name = "Batch Adjuster"
  ) +
  labs(
    title = "Adjuster Performance Rankings Across Class Imbalance Levels",
    subtitle = sprintf("Ranking method: %s | Lower rank = better performance", 
                      switch(args$average_by,
                             "none" = "Individual trials",
                             "test_sets" = "Averaged across test sets",
                             "training_sets" = "Averaged across training set combinations", 
                             "classifiers" = "Averaged across classifiers",
                             "combinations" = "Averaged across dataset combinations")),
    x = "Class Imbalance Level (% Active TB in High-Imbalance Training Set)",
    y = "Average Performance Rank",
    caption = ifelse(args$average_by == "none", 
                    "Points show individual trials with LOESS smoothing",
                    "Error bars show standard error of the mean rank")
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

# Add faceting if we have multiple groups and averaging is applied
if (args$average_by != "none" && length(unique(plot_data$group_id)) > 1 && length(unique(plot_data$group_id)) <= 12) {
  # Only facet if we have a reasonable number of groups
  if (args$average_by %in% c("test_sets", "training_sets")) {
    p <- p + facet_wrap(~ group_id, scales = "free_y", ncol = 3)
  }
}

# Save the plot
ggsave(args$output, p, width = 12, height = 8, dpi = 300)

cat("Plot saved to:", args$output, "\n")

# Print summary statistics
cat("\nSummary of ranking patterns:\n")

if (args$average_by != "none") {
  # Calculate trend statistics for averaged data
  trend_stats <- plot_data %>%
    group_by(adjuster_label) %>%
    do({
      if (nrow(.) >= 3) {
        model <- lm(rank ~ I(imbalance_pct * 100), data = .)
        data.frame(
          slope = coef(model)[2],
          p_value = summary(model)$coefficients[2, 4],
          r_squared = summary(model)$r.squared
        )
      } else {
        data.frame(slope = NA, p_value = NA, r_squared = NA)
      }
    }) %>%
    ungroup()
  
  print(trend_stats)
  
  cat("\nInterpretation:\n")
  cat("- Positive slope: Performance gets worse (higher rank) with increasing imbalance\n")
  cat("- Negative slope: Performance gets better (lower rank) with increasing imbalance\n")
  cat("- p < 0.05: Statistically significant trend\n")
  
} else {
  # For individual points, show distribution by imbalance level
  summary_stats <- plot_data %>%
    group_by(adjuster_label, imbalance_pct) %>%
    summarise(
      mean_rank = mean(rank, na.rm = TRUE),
      median_rank = median(rank, na.rm = TRUE),
      n_trials = n(),
      .groups = "drop"
    )
  
  print(summary_stats)
}

# Calculate and print ranking stability
cat("\nRanking stability (coefficient of variation of ranks):\n")
stability_stats <- plot_data %>%
  group_by(adjuster_label) %>%
  summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    sd_rank = sd(rank, na.rm = TRUE),
    cv_rank = sd_rank / mean_rank,
    .groups = "drop"
  )

print(stability_stats)