# generate_relative_plot.R
# Generate relative performance plot using estimation statistics (dabestr)

generate_relative_plot <- function(mxe_data, top_adjuster, top_adjuster_label, 
                                   output_file, width = 20, height = 16, dpi = 300) {
  
  cat("\nCreating relative performance plot using estimation statistics (dabestr)...\n")
  cat("Using classifier-specific best adjusters as reference\n")
  
  # Load required packages
  if (!require("dabestr", quietly = TRUE)) {
    stop("Package 'dabestr' is required. Install with: install.packages('dabestr')")
  }
  if (!require("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required. Install with: install.packages('cowplot')")
  }
  
  # Get classifiers with data
  classifiers_with_data <- unique(mxe_data$classifier_label[!is.na(mxe_data$classifier_label)])
  
  # Collect all effect size results for saving
  all_effect_results <- data.frame()
  
  # Create individual plots for each classifier
  relative_plot_list <- list()
  
  for (classifier in classifiers_with_data) {
    # Get data for this classifier
    classifier_data <- mxe_data[mxe_data$classifier_label == classifier, ]
    
    if (nrow(classifier_data) == 0) {
      cat(sprintf("Skipping %s - no data\n", classifier))
      next
    }
    
    # Get classifier-specific ordering and best adjuster
    classifier_adjuster_order <- get_classifier_ordering(mxe_data, classifier)
    classifier_specific_labels <- create_adjuster_labels(classifier_adjuster_order)
    classifier_best_adjuster <- classifier_adjuster_order$adjuster[1]
    classifier_best_label <- classifier_specific_labels[1]
    
    cat(sprintf("\nClassifier %s: best adjuster is %s\n", classifier, classifier_best_label))
    
    # Prepare data for dabestr
    plot_data <- classifier_data %>%
      select(adjuster, classifier, classifier_label, n_datasets, test_study, value) %>%
      mutate(condition_id = paste(classifier, n_datasets, test_study, sep = "_"))
    
    # Apply adjuster labels
    plot_data$adjuster_label <- factor(plot_data$adjuster,
                                       levels = classifier_adjuster_order$adjuster,
                                       labels = classifier_specific_labels)
    
    # Add dataset label
    plot_data$dataset_label <- factor(paste(plot_data$n_datasets, "studies"),
                                     levels = c("3 studies", "4 studies", "5 studies", "6 studies"))
    
    # Create separate dabest objects for each n_datasets group
    for (n_studies in c(3, 4, 5, 6)) {
      dataset_label <- paste(n_studies, "studies")
      subset_data <- plot_data %>% filter(dataset_label == !!dataset_label)
      
      if (nrow(subset_data) == 0) next
      
      # Define comparison groups (best vs all others)
      adjusters_in_subset <- unique(subset_data$adjuster_label)
      adjusters_to_compare <- c(classifier_best_label, 
                               classifier_specific_labels[classifier_specific_labels != classifier_best_label & 
                                                         classifier_specific_labels %in% adjusters_in_subset])
      
      if (length(adjusters_to_compare) < 2) next
      
      # Filter to only include adjusters we're comparing
      subset_data <- subset_data %>%
        filter(adjuster_label %in% adjusters_to_compare)
      
      # Create dabest object
      tryCatch({
        # Debug: check data structure
        cat(sprintf("    Attempting dabestr for %s with %d rows, %d adjusters\n", 
                   dataset_label, nrow(subset_data), length(adjusters_to_compare)))
        
        dabest_obj <- dabestr::load(
            data = subset_data,
            x = adjuster_label,
            y = value,
            idx = adjusters_to_compare,
            paired = "baseline",
            id_col = condition_id
          )
        
        if (is.null(dabest_obj)) {
          stop("dabestr::load returned NULL")
        }
        
        # Compute mean differences
        dabest_diff <- dabestr::mean_diff(dabest_obj, perm_count = 5000)
        
        if (is.null(dabest_diff)) {
          stop("mean_diff returned NULL")
        }
        
        if (is.null(dabest_diff$boot_result)) {
          stop("boot_result is NULL")
        }
        
        if (is.null(dabest_diff$permtest_pvals)) {
          stop("permtest_pvals is NULL")
        }
        
        # Extract effect size results from boot_result and permtest_pvals
        boot_data <- dabest_diff$boot_result %>%
          select(control_group, test_group, difference, bca_ci_low, bca_ci_high)
        
        perm_data <- dabest_diff$permtest_pvals %>%
          select(control_group, test_group, pval_permtest)
        
        effect_results <- boot_data %>%
          left_join(perm_data, by = c("control_group", "test_group")) %>%
          rename(
            reference = control_group,
            comparison = test_group,
            mean_difference = difference,
            ci_lower = bca_ci_low,
            ci_upper = bca_ci_high,
            p_value = pval_permtest
          ) %>%
          mutate(
            classifier = classifier,
            dataset_label = dataset_label,
            is_paired = dabest_diff$is_paired
          )
        
        # Collect for saving
        all_effect_results <- rbind(all_effect_results, effect_results)
        
        cat(sprintf("  %s: computed %d comparisons\n", dataset_label, nrow(effect_results)))
        
      }, error = function(e) {
        cat(sprintf("  Warning: Could not compute dabestr for %s %s: %s\n", 
                   classifier, dataset_label, e$message))
      })
    }
    
    # Create Gardner-Altman style plot with raw data and effect sizes
    # Prepare data for visualization
    plot_data <- classifier_data %>%
      select(adjuster, classifier, classifier_label, n_datasets, test_study, value) %>%
      mutate(dataset_label = factor(paste(n_datasets, "studies"),
                                    levels = c("3 studies", "4 studies", "5 studies", "6 studies")))
    
    # Apply ordering
    plot_data$adjuster_label <- factor(plot_data$adjuster,
                                       levels = classifier_adjuster_order$adjuster,
                                       labels = classifier_specific_labels)
    
    # Calculate means for raw data
    mean_data <- plot_data %>%
      group_by(adjuster_label, dataset_label) %>%
      summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")
    
    # Get effect size data for this classifier
    effect_data <- all_effect_results %>%
      filter(classifier == !!classifier) %>%
      mutate(dataset_label = factor(dataset_label,
                                    levels = c("3 studies", "4 studies", "5 studies", "6 studies")))
    
    # Create two-panel plot: raw data on top, effect sizes on bottom
    # Top panel: Raw MCC values
    p_raw <- ggplot(plot_data, aes(x = adjuster_label, y = value)) +
      geom_point(aes(color = test_study, shape = test_study), size = 2.5, alpha = 0.6, 
                position = position_jitter(width = 0.1, height = 0)) +
      geom_segment(data = mean_data, 
                  aes(x = as.numeric(adjuster_label) - 0.35, 
                      xend = as.numeric(adjuster_label) + 0.35,
                      y = mean_value, yend = mean_value),
                  color = "black", linewidth = 1) +
      facet_wrap(~ dataset_label, scales = "free_x", ncol = 4) +
      scale_color_brewer(palette = "Set2", name = "Test Study") +
      scale_shape_manual(values = c(15, 16, 17, 25, 18, 19), name = "Test Study") +
      theme_bw() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 9),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 9),
        plot.margin = margin(5, 5, 0, 5)
      ) +
      labs(y = "MCC", title = classifier)
    
    # Bottom panel: Effect sizes (mean differences with CI)
    if (nrow(effect_data) > 0) {
      # Prepare effect data with proper factor levels
      effect_data$comparison <- factor(effect_data$comparison,
                                       levels = classifier_specific_labels)
      
      p_effect <- ggplot(effect_data, aes(x = comparison, y = mean_difference)) +
        geom_hline(yintercept = 0, linetype = "solid", color = "gray60", linewidth = 0.5) +
        geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 0.8) +
        geom_point(size = 3, shape = 21, fill = "white", stroke = 1.5) +
        geom_text(aes(label = ifelse(p_value < 0.05, "*", "")), 
                 vjust = -1.5, size = 5, fontface = "bold") +
        facet_wrap(~ dataset_label, scales = "free_x", ncol = 4) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          axis.title.y = element_text(size = 9),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_blank(),
          plot.margin = margin(0, 5, 5, 5)
        ) +
        labs(x = "Adjuster", y = sprintf("Δ MCC\n(vs. %s)", classifier_best_label))
    } else {
      # If no effect data, create empty bottom panel
      p_effect <- ggplot() + 
        theme_void() +
        labs(caption = "Effect sizes not computed (dabestr failed)")
    }
    
    # Combine panels
    p <- cowplot::plot_grid(p_raw, p_effect, ncol = 1, rel_heights = c(1.2, 1), align = "v")
    
    relative_plot_list[[classifier]] <- p
  }
  
  # Arrange and save
  cat("\nArranging relative plots in grid\n")
  relative_final_plot <- gridExtra::grid.arrange(grobs = unname(relative_plot_list), ncol = 2)
  
  cat("Saving relative performance plot to:", output_file, "\n")
  ggplot2::ggsave(filename = output_file, plot = relative_final_plot, 
         width = width, height = height, dpi = dpi, units = "in")
  
  cat("Relative performance plot saved successfully!\n")
  
  # Save effect size results (always create the file, even if empty)
  effect_output <- sub("\\.png$", "_effect_sizes.csv", output_file)
  if (nrow(all_effect_results) > 0) {
    write.csv(all_effect_results, effect_output, row.names = FALSE)
    cat("Effect size results saved to:", effect_output, "\n")
  } else {
    # Create empty CSV with expected columns
    empty_results <- data.frame(
      reference = character(),
      comparison = character(),
      mean_difference = numeric(),
      ci_lower = numeric(),
      ci_upper = numeric(),
      p_value = numeric(),
      is_paired = logical(),
      classifier = character(),
      dataset_label = character()
    )
    write.csv(empty_results, effect_output, row.names = FALSE)
    cat("No effect size results computed (dabestr failed). Empty CSV saved to:", effect_output, "\n")
  }
}
