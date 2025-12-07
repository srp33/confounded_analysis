# generate_differences_plot.R
# Generate performance difference distribution plot with statistical tests

#' Generate differences plot
#' @param mxe_data Data frame with MCC values
#' @param top_adjuster Name of top performing adjuster
#' @param top_adjuster_label Formatted label for top adjuster
#' @param unique_adjusters Vector of all adjuster names
#' @param adjuster_labels Vector of formatted adjuster labels
#' @param output_file Path to save plot
#' @param dpi Plot resolution
#' @return Data frame with statistical test results
generate_differences_plot <- function(mxe_data, top_adjuster, top_adjuster_label, 
                                     unique_adjusters, adjuster_labels, 
                                     output_file, dpi = 300) {
  
  cat("\nCreating performance difference distribution plot...\n")
  cat("Top adjuster:", top_adjuster, "(", top_adjuster_label, ")\n")
  
  # Calculate differences
  difference_data <- mxe_data %>%
    mutate(condition_id = paste(classifier, n_datasets, test_study, sep = "_")) %>%
    group_by(condition_id) %>%
    mutate(top_value = value[adjuster == top_adjuster]) %>%
    ungroup() %>%
    mutate(difference = top_value - value) %>%
    filter(adjuster != top_adjuster)
  
  difference_data$adjuster_label <- factor(difference_data$adjuster,
    levels = unique_adjusters[-1],
    labels = adjuster_labels[-1])
  
  cat("Calculated", nrow(difference_data), "pairwise differences\n")
  
  # Statistical testing
  cat("\nPerforming statistical tests...\n")
  cat("Adjusters compared to", top_adjuster_label, ":\n")
  
  stat_results <- data.frame()
  
  for (adj in unique(difference_data$adjuster_label)) {
    adj_data <- difference_data[difference_data$adjuster_label == adj, ]
    n_obs <- nrow(adj_data)
    mean_diff <- mean(adj_data$difference, na.rm = TRUE)
    median_diff <- median(adj_data$difference, na.rm = TRUE)
    sd_diff <- sd(adj_data$difference, na.rm = TRUE)
    
    t_test <- t.test(adj_data$difference, mu = 0, alternative = "greater")
    wilcox_test <- wilcox.test(adj_data$difference, mu = 0, alternative = "greater")
    
    get_stars <- function(p) {
      if (p < 0.001) return("***")
      if (p < 0.01) return("**")
      if (p < 0.05) return("*")
      return("ns")
    }
    
    format_pval <- function(p) {
      if (p < 0.001) return("p < 0.001")
      if (p < 0.01) return(sprintf("p = %.3f", p))
      return(sprintf("p = %.2f", p))
    }
    
    cat(sprintf("  %s: n=%d, mean_diff=%.4f, median_diff=%.4f, sd=%.4f\n", 
                adj, n_obs, mean_diff, median_diff, sd_diff))
    cat(sprintf("    t-test: t=%.3f, %s %s\n", 
                t_test$statistic, format_pval(t_test$p.value), get_stars(t_test$p.value)))
    cat(sprintf("    Wilcoxon: %s %s\n", 
                format_pval(wilcox_test$p.value), get_stars(wilcox_test$p.value)))
    
    stat_results <- rbind(stat_results, data.frame(
      adjuster_label = adj,
      n = n_obs,
      mean_diff = mean_diff,
      median_diff = median_diff,
      sd_diff = sd_diff,
      t_statistic = t_test$statistic,
      t_pvalue = t_test$p.value,
      wilcox_pvalue = wilcox_test$p.value,
      significance = get_stars(t_test$p.value),
      stringsAsFactors = FALSE
    ))
  }
  
  stat_results$pval_label <- sapply(stat_results$t_pvalue, function(p) {
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    return("ns")
  })
  
  # Create plot
  y_max <- max(difference_data$difference, na.rm = TRUE)
  y_min <- min(difference_data$difference, na.rm = TRUE)
  y_range <- y_max - y_min
  annotation_y <- y_max + 0.05 * y_range
  
  diff_plot <- ggplot(difference_data, aes(x = adjuster_label, y = difference)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_violin(fill = "#E69F00", alpha = 0.3, color = NA) +
    geom_boxplot(width = 0.2, outlier.size = 1, outlier.alpha = 0.5, fill = "white") +
    geom_jitter(aes(color = classifier_label), width = 0.15, height = 0, 
                size = 1.5, alpha = 0.4) +
    geom_text(data = stat_results, 
              aes(x = adjuster_label, y = annotation_y, label = pval_label),
              size = 5, fontface = "bold", vjust = 0) +
    scale_color_brewer(palette = "Set3", name = "Classifier") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 9, hjust = 0, color = "gray40")
    ) +
    labs(
      x = "Adjuster",
      y = sprintf("Performance Difference\n(%s MCC - Adjuster MCC)", top_adjuster_label),
      title = "Distribution of Performance Differences Relative to Top Adjuster",
      subtitle = sprintf("Positive values indicate %s performs better", top_adjuster_label),
      caption = "Significance from one-sample t-test (H₁: difference > 0): *** p<0.001, ** p<0.01, * p<0.05, ns = not significant"
    )
  
  # Save plot
  cat("Saving difference plot to:", output_file, "\n")
  ggsave(filename = output_file, plot = diff_plot, width = 12, height = 8, dpi = dpi, units = "in")
  cat("Difference plot saved successfully!\n")
  
  # Return statistics for saving
  stat_results
}
