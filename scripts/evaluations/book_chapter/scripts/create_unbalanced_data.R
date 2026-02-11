#!/usr/bin/env Rscript

# create_unbalanced_data.R - Create batch-specific imbalanced TB datasets
# Creates different class distributions in different training batches to test Combat's handling of confounded effects

suppressMessages(suppressWarnings({
  required_packages <- c("argparse")
  sapply(required_packages, require, character.only=TRUE, quietly=TRUE)
}))

# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Create batch-specific imbalanced TB dataset")

parser$add_argument("--input", type = "character", required = TRUE,
                   help = "Input RData file path (TB_real_data.RData)")
parser$add_argument("--output", type = "character", required = TRUE,
                   help = "Output RData file path for unbalanced data")
parser$add_argument("--num-datasets", type = "integer", required = TRUE,
                   help = "Number of datasets (3 or 5)")
parser$add_argument("--test-study", type = "character", required = TRUE,
                   help = "Test study name")
parser$add_argument("--seed", type = "integer", default = 123,
                   help = "Random seed for reproducible filtering (default: 123)")

args <- parser$parse_args()

# ====================================================================
# BATCH-SPECIFIC IMBALANCE PATTERNS
# ====================================================================

# Define imbalance patterns for different scenarios
get_imbalance_pattern <- function(num_datasets, test_study) {
  # All available studies in order
  all_studies <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
  
  # Get the selected studies for this num_datasets (same logic as classify_unbalanced.R)
  selected_studies <- all_studies[1:num_datasets]
  
  # Training studies are all selected studies except the test study
  train_studies <- selected_studies[selected_studies != test_study]
  
  if (num_datasets == 3) {
    # 3-dataset: Strong confounding between batch and class
    # Assign alternating high/low TB ratios to training studies
    train_ratios <- rep(c(0.8, 0.2), length.out = length(train_studies))
    return(list(
      train_ratios = train_ratios,
      train_names = train_studies,
      all_studies = selected_studies,
      test_study = test_study
    ))
  } else if (num_datasets == 5) {
    # 5-dataset: Gradient of confounding
    # Assign gradient of TB ratios to training studies
    train_ratios <- rep(c(0.8, 0.2, 0.6, 0.4), length.out = length(train_studies))
    return(list(
      train_ratios = train_ratios,
      train_names = train_studies,
      all_studies = selected_studies,
      test_study = test_study
    ))
  } else {
    stop(sprintf("Unsupported num_datasets: %d. Only 3 and 5 are supported.", num_datasets))
  }
}

# ====================================================================
# MAIN FILTERING FUNCTION
# ====================================================================

create_batch_imbalanced_data <- function(input_file, output_file, num_datasets, test_study, seed) {
  
  cat("=== BATCH-SPECIFIC IMBALANCED TB DATA CREATION ===\n")
  cat(sprintf("Input: %s\n", input_file))
  cat(sprintf("Output: %s\n", output_file))
  cat(sprintf("Num datasets: %d\n", num_datasets))
  cat(sprintf("Test study: %s\n", test_study))
  cat(sprintf("Seed: %d\n", seed))
  cat("==================================================\n\n")
  
  # Load original data
  if (!file.exists(input_file)) {
    stop(sprintf("Input file not found: %s", input_file))
  }
  
  load(input_file)  # Loads dat_lst and label_lst
  
  # Validate loaded data
  if (!exists("dat_lst") || !exists("label_lst")) {
    stop("Required objects 'dat_lst' and 'label_lst' not found in input file")
  }
  
  # Get imbalance pattern for this scenario
  pattern <- get_imbalance_pattern(num_datasets, test_study)
  train_studies <- pattern$train_names
  target_ratios <- pattern$train_ratios
  all_studies <- pattern$all_studies
  
  cat("Original data loaded:\n")
  cat(sprintf("  Studies: %s\n", paste(names(dat_lst), collapse=", ")))
  cat(sprintf("  Selected studies (%d total): %s\n", num_datasets, paste(all_studies, collapse=", ")))
  cat(sprintf("  Training studies: %s\n", paste(train_studies, collapse=", ")))
  cat(sprintf("  Test study: %s\n", test_study))
  cat(sprintf("  Target TB ratios for training: %s\n", paste(sprintf("%.1f%%", target_ratios*100), collapse=", ")))
  
  # Print original label distribution for all selected studies
  cat("\nOriginal label distribution (0=LTBI, 1=Active TB):\n")
  for (study in all_studies) {
    if (study %in% names(label_lst)) {
      labels <- label_lst[[study]]
      n_ltbi <- sum(labels == 0)
      n_tb <- sum(labels == 1)
      total <- length(labels)
      cat(sprintf("  %s: %d LTBI (%.1f%%), %d Active TB (%.1f%%), Total: %d\n", 
                  study, n_ltbi, 100*n_ltbi/total, n_tb, 100*n_tb/total, total))
    }
  }
  
  # Set seed for reproducible filtering
  set.seed(seed)
  
  # Create batch-imbalanced datasets
  dat_lst_imbalanced <- list()
  label_lst_imbalanced <- list()
  
  cat(sprintf("\nCreating batch-specific imbalanced datasets:\n"))
  
  # Process training studies with specific imbalance patterns
  for (i in seq_along(train_studies)) {
    study <- train_studies[i]
    target_tb_ratio <- target_ratios[i]
    
    if (!study %in% names(dat_lst)) {
      stop(sprintf("Training study '%s' not found in data", study))
    }
    
    labels <- label_lst[[study]]
    data <- dat_lst[[study]]
    
    # Identify LTBI (0) and Active TB (1) samples
    ltbi_indices <- which(labels == 0)
    tb_indices <- which(labels == 1)
    
    n_ltbi_orig <- length(ltbi_indices)
    n_tb_orig <- length(tb_indices)
    
    # Calculate target sample sizes to achieve desired TB ratio
    # We'll keep a reasonable total sample size (at least 30 samples)
    min_total_samples <- 30
    max_total_samples <- length(labels)
    
    # Try different total sample sizes to find one that works
    best_total <- NULL
    best_n_tb <- NULL
    best_n_ltbi <- NULL
    
    for (total_samples in seq(min_total_samples, max_total_samples, by = 5)) {
      n_tb_target <- round(total_samples * target_tb_ratio)
      n_ltbi_target <- total_samples - n_tb_target
      
      # Check if we have enough samples of each class
      if (n_tb_target <= n_tb_orig && n_ltbi_target <= n_ltbi_orig && 
          n_tb_target >= 3 && n_ltbi_target >= 3) {
        best_total <- total_samples
        best_n_tb <- n_tb_target
        best_n_ltbi <- n_ltbi_target
        break
      }
    }
    
    if (is.null(best_total)) {
      # Fallback: use maximum possible samples while maintaining ratio
      max_tb_for_ratio <- floor(n_tb_orig / target_tb_ratio)
      max_ltbi_for_ratio <- floor(n_ltbi_orig / (1 - target_tb_ratio))
      
      if (max_tb_for_ratio <= max_ltbi_for_ratio) {
        best_n_tb <- n_tb_orig
        best_n_ltbi <- round(n_tb_orig * (1 - target_tb_ratio) / target_tb_ratio)
      } else {
        best_n_ltbi <- n_ltbi_orig
        best_n_tb <- round(n_ltbi_orig * target_tb_ratio / (1 - target_tb_ratio))
      }
      best_total <- best_n_tb + best_n_ltbi
    }
    
    # Sample the specified numbers
    keep_tb <- sample(tb_indices, best_n_tb, replace = FALSE)
    keep_ltbi <- sample(ltbi_indices, best_n_ltbi, replace = FALSE)
    
    # Combine indices to keep
    keep_indices <- c(keep_ltbi, keep_tb)
    keep_indices <- sort(keep_indices)
    
    # Filter data and labels
    dat_lst_imbalanced[[study]] <- data[, keep_indices, drop = FALSE]
    label_lst_imbalanced[[study]] <- labels[keep_indices]
    
    # Report filtering results
    new_labels <- label_lst_imbalanced[[study]]
    n_ltbi_new <- sum(new_labels == 0)
    n_tb_new <- sum(new_labels == 1)
    total_new <- length(new_labels)
    actual_tb_ratio <- n_tb_new / total_new
    
    cat(sprintf("  %s (TRAIN): %d LTBI (%.1f%%), %d Active TB (%.1f%%), Total: %d (%.1f%% of original)\n", 
                study, n_ltbi_new, 100*n_ltbi_new/total_new, n_tb_new, 100*n_tb_new/total_new, 
                total_new, 100*total_new/length(labels)))
    cat(sprintf("    Target TB ratio: %.1f%%, Actual TB ratio: %.1f%%\n", 
                100*target_tb_ratio, 100*actual_tb_ratio))
  }
  
  # Keep test study unchanged
  dat_lst_imbalanced[[test_study]] <- dat_lst[[test_study]]
  label_lst_imbalanced[[test_study]] <- label_lst[[test_study]]
  
  test_labels <- label_lst_imbalanced[[test_study]]
  n_ltbi_test <- sum(test_labels == 0)
  n_tb_test <- sum(test_labels == 1)
  total_test <- length(test_labels)
  
  cat(sprintf("  %s (TEST): %d LTBI (%.1f%%), %d Active TB (%.1f%%), Total: %d (unchanged)\n", 
              test_study, n_ltbi_test, 100*n_ltbi_test/total_test, n_tb_test, 100*n_tb_test/total_test, total_test))
  
  # Add any remaining studies from all_studies that aren't training or test (shouldn't happen in practice)
  for (study in all_studies) {
    if (!study %in% names(dat_lst_imbalanced)) {
      dat_lst_imbalanced[[study]] <- dat_lst[[study]]
      label_lst_imbalanced[[study]] <- label_lst[[study]]
      cat(sprintf("  %s (OTHER): kept unchanged\n", study))
    }
  }
  
  # Validate that we have samples in all required studies
  for (study in all_studies) {
    if (!study %in% names(dat_lst_imbalanced)) {
      stop(sprintf("Study %s missing from imbalanced data", study))
    }
    if (ncol(dat_lst_imbalanced[[study]]) == 0) {
      stop(sprintf("Study %s has no samples after filtering", study))
    }
    if (length(unique(label_lst_imbalanced[[study]])) < 2) {
      warning(sprintf("Study %s has only one class after filtering", study))
    }
  }
  
  # Save imbalanced data with same variable names for compatibility
  dat_lst <- dat_lst_imbalanced
  label_lst <- label_lst_imbalanced
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  save(dat_lst, label_lst, file = output_file)
  
  cat(sprintf("\n✅ Batch-imbalanced data saved to: %s\n", output_file))
  
  # Final summary
  cat("\nFinal batch-imbalanced dataset summary:\n")
  total_ltbi <- sum(sapply(label_lst[train_studies], function(x) sum(x == 0)))
  total_tb <- sum(sapply(label_lst[train_studies], function(x) sum(x == 1)))
  total_samples <- total_ltbi + total_tb
  
  cat(sprintf("  Training total LTBI: %d (%.1f%%)\n", total_ltbi, 100*total_ltbi/total_samples))
  cat(sprintf("  Training total Active TB: %d (%.1f%%)\n", total_tb, 100*total_tb/total_samples))
  cat(sprintf("  Training total samples: %d\n", total_samples))
  
  # Calculate batch-class confounding strength
  batch_tb_ratios <- sapply(train_studies, function(study) {
    labels <- label_lst[[study]]
    sum(labels == 1) / length(labels)
  })
  confounding_strength <- sd(batch_tb_ratios)
  cat(sprintf("  Batch-class confounding strength (SD of TB ratios): %.3f\n", confounding_strength))
}

# ====================================================================
# EXECUTE MAIN FUNCTION
# ====================================================================

tryCatch({
  create_batch_imbalanced_data(args$input, args$output, args$num_datasets, args$test_study, args$seed)
}, error = function(e) {
  cat(sprintf("[ERROR] %s\n", e$message), file = stderr())
  quit(status = 1)
})