#!/usr/bin/env Rscript

# Plot ranking comparison between balanced and unbalanced TB analysis
# Ranks adjusters within each classifier/test_study/num_datasets combination
# Then averages ranks within classifier and plots with error bars

library(argparse)
library(dplyr)
library(ggplot2)
library(tidyr)

# Parse command line arguments
parser <- ArgumentParser(description = "Create ranking comparison plot between balanced and unbalanced results")
parser$add_argument("--balanced-data", required = TRUE, help = "Path to balanced results CSV")
parser$add_argument("--unbalanced-data", required = TRUE, help = "Path to unbalanced results CSV")
parser$add_argument("-o", "--output", required = TRUE, help = "Output PNG file path")

args <- parser$parse_args()

# Load data
cat("Loading balanced results from:", args$balanced_data, "\n")
balanced_data <- read.csv(args$balanced_data, stringsAsFactors = FALSE)

cat("Loading unbalanced results from:", args$unbalanced_data, "\n")
unbalanced_data <- read.csv(args$unbalanced_data, stringsAsFactors = FALSE)

# Filter balanced data to match unbalanced scenarios
# Unbalanced only tests 3 and 5 dataset scenarios with specific adjusters
unbalanced_adjusters <- c("unadjusted", "combat", "combat_sup")
test_scenarios <- c(3, 5)

balanced_filtered <- balanced_data %>%
  filter(
    adjuster %in% unbalanced_adjusters,
    n_datasets %in% test_scenarios,
    metric == "mcc"  # Only use MCC for ranking
  ) %>%
  mutate(
    condition = "balanced",
    mcc = value,
    num_datasets = n_datasets  # Standardize column name
  )

unbalanced_filtered <- unbalanced_data %>%
  filter(adjuster %in% unbalanced_adjusters) %>%
  mutate(
    condition = "unbalanced",
    n_datasets = num_datasets  # Standardize column name
  )

# Combine datasets
combined_data <- bind_rows(balanced_filtered, unbalanced_filtered)

cat("Combined data dimensions:", nrow(combined_data), "rows\n")
cat("Conditions:", unique(combined_data$condition), "\n")
cat("Adjusters:", unique(combined_data$adjuster), "\n")
cat("Classifiers:", unique(combined_data$classifier), "\n")

# Rank adjusters within each classifier/test_study/n_datasets/condition combination
ranked_data <- combined_data %>%
  group_by(classifier, test_study, n_datasets, condition) %>%
  arrange(desc(mcc)) %>%  # Higher MCC = better performance = lower rank number
  mutate(rank = rank(-mcc, ties.method = "average")) %>%  # Negative MCC so higher MCC gets rank 1
  ungroup()

# Calculate average ranks and standard errors within each classifier/condition
avg_ranks <- ranked_data %>%
  group_by(classifier, adjuster, condition) %>%
  summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    se_rank = sd(rank, na.rm = TRUE) / sqrt(n()),
    n_obs = n(),
    .groups = "drop"
  )

cat("Average ranks calculated for", nrow(avg_ranks), "combinations\n")

# Create adjuster labels for plotting
avg_ranks <- avg_ranks %>%
  mutate(
    adjuster_label = case_when(
      adjuster == "unadjusted" ~ "Unadjusted",
      adjuster == "combat" ~ "ComBat (Unsupervised)",
      adjuster == "combat_sup" ~ "ComBat (Supervised)",
      TRUE ~ adjuster
    ),
    adjuster_label = factor(adjuster_label, levels = c("Unadjusted", "ComBat (Unsupervised)", "ComBat (Supervised)")),
    classifier_label = case_when(
      classifier == "elasticnet" ~ "Elastic Net",
      classifier == "knn" ~ "k-NN",
      classifier == "logistic" ~ "Logistic Regression",
      classifier == "nnet" ~ "Neural Network",
      classifier == "rf" ~ "Random Forest",
      classifier == "shrinkageLDA" ~ "Shrinkage LDA",
      classifier == "svm" ~ "SVM",
      classifier == "xgboost" ~ "XGBoost",
      TRUE ~ classifier
    ),
    classifier_label = factor(classifier_label, levels = c("Elastic Net", "k-NN", "Logistic Regression", 
                                                          "Neural Network", "Random Forest", "Shrinkage LDA", 
                                                          "SVM", "XGBoost"))
  )

# Create the plot
p <- ggplot(avg_ranks, aes(x = classifier, y = mean_rank, color = adjuster_label, group = adjuster_label)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = mean_rank - se_rank, ymax = mean_rank + se_rank), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  geom_line(aes(linetype = condition), position = position_dodge(width = 0.3), alpha = 0.7) +
  facet_wrap(~ condition, ncol = 2, labeller = labeller(condition = c("balanced" = "Balanced", "unbalanced" = "Imbalanced"))) +
  scale_y_reverse(breaks = c(1, 2, 3), labels = c("1st", "2nd", "3rd")) +
  scale_color_manual(
    values = c("Unadjusted" = "#E31A1C", "ComBat (Unsupervised)" = "#1F78B4", "ComBat (Supervised)" = "#33A02C"),
    name = "Batch Adjuster"
  ) +
  scale_linetype_manual(values = c("balanced" = "solid", "unbalanced" = "dashed"), guide = "none") +
  labs(
    title = "Adjuster Performance Rankings: Balanced vs Imbalanced Data",
    subtitle = "Average rank across test scenarios (3 and 5 datasets)\nLower rank = better performance",
    x = "Classifier",
    y = "Average Rank",
    caption = "Error bars show standard error of the mean rank"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

# Add connecting lines between conditions for the same adjuster/classifier
# This requires reshaping the data to have balanced and unbalanced in separate columns
connection_data <- avg_ranks %>%
  select(classifier_label, adjuster_label, condition, mean_rank) %>%
  pivot_wider(names_from = condition, values_from = mean_rank, names_prefix = "rank_") %>%
  filter(!is.na(rank_balanced) & !is.na(rank_unbalanced))

# Create a version with connecting lines
dodge_width <- 0.3

# Create position mapping for dodging
pos_dodge <- position_dodge(width = dodge_width)

p_with_lines <- ggplot(avg_ranks, aes(x = classifier_label, y = mean_rank, color = adjuster_label, shape = condition)) +
  # Add connecting lines for each adjuster between balanced and unbalanced
  geom_line(aes(group = adjuster_label), alpha = 0.6, linewidth = 0.8, position = pos_dodge) +
  geom_errorbar(aes(ymin = mean_rank - se_rank, ymax = mean_rank + se_rank), 
                width = 0.2, position = pos_dodge) +
  geom_point(size = 3, position = pos_dodge) +
  scale_y_reverse(breaks = c(1, 2, 3), labels = c("1st", "2nd", "3rd")) +
  scale_color_manual(
    values = c("Unadjusted" = "#E31A1C", "ComBat (Unsupervised)" = "#1F78B4", "ComBat (Supervised)" = "#33A02C"),
    name = "Batch Adjuster"
  ) +
  scale_shape_manual(
    values = c("balanced" = 16, "unbalanced" = 17),
    labels = c("balanced" = "Balanced", "unbalanced" = "Imbalanced"),
    name = "Data Condition"
  ) +
  labs(
    title = "Adjuster Performance Rankings: Balanced vs Imbalanced Data",
    subtitle = "Average rank across test scenarios (3 and 5 datasets)\nLower rank = better performance",
    x = "Classifier",
    y = "Average Rank",
    caption = "Error bars show standard error of the mean rank"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11)
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 16)),
    shape = guide_legend(override.aes = list(color = "black"))
  )

# Save the plot
ggsave(args$output, p_with_lines, width = 12, height = 8, dpi = 300)

cat("Plot saved to:", args$output, "\n")

# Print summary statistics
cat("\nSummary of ranking changes:\n")
ranking_summary <- avg_ranks %>%
  select(classifier_label, adjuster_label, condition, mean_rank) %>%
  pivot_wider(names_from = condition, values_from = mean_rank, names_prefix = "rank_") %>%
  mutate(rank_change = rank_unbalanced - rank_balanced) %>%
  filter(!is.na(rank_change))

print(ranking_summary)

cat("\nAverage ranking change by adjuster:\n")
adjuster_changes <- ranking_summary %>%
  group_by(adjuster_label) %>%
  summarise(
    mean_change = mean(rank_change, na.rm = TRUE),
    se_change = sd(rank_change, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
print(adjuster_changes)