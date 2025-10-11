# GMM Debugging Helper Functions
debug_gmm_data <- function(data, context = "") {
  if (nrow(data) > 0) {
    methods <- unique(data$Method)
    gmm_present <- "GMM" %in% methods
    gmm_count <- sum(data$Method == "GMM", na.rm = TRUE)
    
    cat("DEBUG GMM", context, ":\n")
    cat("  - Methods found:", paste(methods, collapse = ", "), "\n")
    cat("  - GMM present:", gmm_present, "\n")
    cat("  - GMM rows:", gmm_count, "\n")
    
    if (gmm_present && gmm_count > 0) {
      gmm_sample <- data[data$Method == "GMM", ][1:min(3, gmm_count), ]
      cat("  - Sample GMM values:", paste(round(gmm_sample$value, 4), collapse = ", "), "\n")
    }
  } else {
    cat("DEBUG GMM", context, ": NO DATA\n")
  }
}

debug_figure_data <- function(plot_data, figure_name = "") {
  if ("Method" %in% colnames(plot_data)) {
    methods_in_plot <- unique(plot_data$Method)
    gmm_in_plot <- "GMM" %in% methods_in_plot
    
    cat("DEBUG FIGURE", figure_name, ":\n")
    cat("  - Methods in plot:", paste(methods_in_plot, collapse = ", "), "\n") 
    cat("  - GMM in final plot:", gmm_in_plot, "\n")
    cat("  - Total plot rows:", nrow(plot_data), "\n")
  }
}