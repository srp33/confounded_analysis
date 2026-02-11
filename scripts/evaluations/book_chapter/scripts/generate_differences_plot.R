# generate_differences_plot.R
# Generate performance difference plot using estimation statistics (dabestr)

generate_differences_plot <- function(mxe_data, top_adjuster, top_adjuster_label, 
                                     unique_adjusters, adjuster_labels, 
                                     output_file, dpi = 300) {
  
  cat("\nCreating performance difference plot using estimation statistics (dabestr)...\n")
  cat("Top adjuster:", top_adjuster, "(", top_adjuster_label, ")\n")
  
  # Load dabestr package
  if (!require("dabestr", quietly = TRUE)) {
    stop("Package 'dabestr' is required. Install with: install.packages('dabestr')")
  }
  
  # Prepare data for dabestr
  plot_data <- mxe_data %>%
    mutate(condition_id = paste(classifier, n_datasets, test_study, sep = "_"))
  
  plot_data$adjuster_label <- factor(plot_data$adjuster,
    levels = unique_adjusters,
    labels = adjuster_labels)
  
  adjusters_to_compare <- c(top_adjuster_label, adjuster_labels[adjuster_labels != top_adjuster_label])
  
  cat("Comparing", top_adjuster_label, "against", length(adjusters_to_compare) - 1, "other adjusters\n")
  cat("Adjusters:", paste(adjusters_to_compare, collapse = ", "), "\n")
  
  plot_data <- plot_data %>%
    filter(adjuster_label %in% adjusters_to_compare)
  
  cat("\nCreating dabest object with paired analysis...\n")
  dabest_obj <- plot_data %>%
    dabestr::dabest(
      x = adjuster_label,
      y = value,
      idx = adjusters_to_compare,
      paired = TRUE,
      id.col = condition_id
    )
  
  cat("Computing mean differences using 5000 bootstrap resamples...\n")
  dabest_diff <- dabestr::mean_diff(dabest_obj, resamples = 5000)
  
  effect_results <- dabest_diff$result %>%
    select(control_group, test_group, difference, ci_low, ci_high, 
           pvalue_permutation, is_paired) %>%
    rename(
      reference = control_group,
      comparison = test_group,
      mean_difference = difference,
      ci_lower = ci_low,
      ci_upper = ci_high,
      p_value = pvalue_permutation
    )
  
  cat("\nEffect size summary:\n")
  for (i in 1:nrow(effect_results)) {
    cat(sprintf("  %s vs %s: Δ = %.4f [95%% CI: %.4f, %.4f], p = %.4f\n",
                effect_results$reference[i],
                effect_results$comparison[i],
                effect_results$mean_difference[i],
                effect_results$ci_lower[i],
                effect_results$ci_upper[i],
                effect_results$p_value[i]))
  }
  
  # Create and save plot
  cat("\nGenerating Gardner-Altman plot...\n")
  
  # Save plot
  png(filename = output_file, width = 12, height = 8, units = "in", res = dpi, bg = "white")
  
  plot(dabest_diff, 
       rawplot.ylabel = "Matthews Correlation Coefficient (MCC)",
       effsize.ylabel = sprintf("MCC Difference\n(vs. %s)", top_adjuster_label),
       theme = ggplot2::theme_bw(base_size = 12))
  
  dev.off()
  
  cat("Difference plot saved to:", output_file, "\n")
  cat("Plot saved successfully!\n")
  
  # Return effect size statistics for saving
  effect_results
}
