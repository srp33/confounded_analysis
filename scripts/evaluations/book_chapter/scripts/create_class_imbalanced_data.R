#!/usr/bin/env Rscript

# create_class_imbalanced_data.R - Create systematically class-imbalanced TB datasets
# For each pair of training datasets, determine optimal imbalance arrangement and test on remaining datasets

suppressMessages(suppressWarnings({
  required_packages <- c("argparse", "combinat")
  sapply(required_packages, require, character.only=TRUE, quietly=TRUE)
}))

# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Create systematically class-imbalanced TB datasets")

parser$add_argument("--input", type = "character", required = TRUE,
                   help = "Input RData file path (TB_real_data.RData)")
parser$add_argument("--output-dir", type = "character", required = TRUE,
                   help = "Output directory for imbalanced datasets")
parser$add_argument("--report-file", type = "character", required = TRUE,
                   help = "Output file for imbalance report")
parser$add_argument("--seed", type = "integer", default = 123,
                   help = "Random seed for reproducible sampling (default: 123)")
parser$add_argument("--n-replicates", type = "integer", default = 5,
                   help = "Number of replicate random samples per scenario (default: 5)")

args <- parser$parse_args()

# ====================================================================
# HELPER FUNCTIONS
# ====================================================================

# Generate all pairs of training datasets from 6 total datasets
generate_training_pairs <- function() {
  all_studies <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
  pairs <- combn(all_studies, 2, simplify = FALSE)
  return(pairs)
}

# Calculate optimal sample sizes for given imbalance target, keeping total constant
calculate_optimal_samples <- function(n_active, n_latent, target_active_pct, target_total) {
  
  # Calculate target sample sizes
  target_active <- round(target_total * target_active_pct)
  target_latent <- target_total - target_active
  
  # Check feasibility
  if (target_active > n_active || target_latent > n_latent || 
      target_active < 3 || target_latent < 3) {
    return(NULL)  # Not feasible
  }
  
  return(list(
    total = target_total,
    active = target_active,
    latent = target_latent,
    actual_active_pct = target_active / target_total
  ))
}

# Determine optimal arrangement and CONSISTENT sample size across all imbalance levels
determine_optimal_arrangement_consistent <- function(dataset1_stats, dataset2_stats, all_imbalance_pcts) {
  
  # For each arrangement, find the MINIMUM total samples across ALL imbalance levels
  # This ensures consistent sample size regardless of imbalance level
  
  # Arrangement 1: dataset1 high active, dataset2 low active
  max_totals_1 <- sapply(all_imbalance_pcts, function(imbalance_pct) {
    min(
      floor(dataset1_stats$n_active / imbalance_pct),
      floor(dataset1_stats$n_latent / (1 - imbalance_pct)),
      floor(dataset2_stats$n_active / (1 - imbalance_pct)),
      floor(dataset2_stats$n_latent / imbalance_pct)
    )
  })
  consistent_total_1 <- min(max_totals_1)
  
  # Arrangement 2: dataset1 low active, dataset2 high active  
  max_totals_2 <- sapply(all_imbalance_pcts, function(imbalance_pct) {
    min(
      floor(dataset1_stats$n_active / (1 - imbalance_pct)),
      floor(dataset1_stats$n_latent / imbalance_pct),
      floor(dataset2_stats$n_active / imbalance_pct),
      floor(dataset2_stats$n_latent / (1 - imbalance_pct))
    )
  })
  consistent_total_2 <- min(max_totals_2)
  
  # Choose arrangement that maximizes the CONSISTENT total samples
  if (consistent_total_1 >= consistent_total_2 && consistent_total_1 >= 60) {
    return(list(
      arrangement = 1,
      high_active_dataset = dataset1_stats$dataset,
      low_active_dataset = dataset2_stats$dataset,
      total_samples = consistent_total_1,  # Same for all imbalance levels
      feasible = TRUE
    ))
  } else if (consistent_total_2 >= 60) {
    return(list(
      arrangement = 2,
      high_active_dataset = dataset2_stats$dataset,
      low_active_dataset = dataset1_stats$dataset,
      total_samples = consistent_total_2,  # Same for all imbalance levels
      feasible = TRUE
    ))
  } else {
    return(list(feasible = FALSE))
  }
}

# Create imbalanced subset for a single dataset
create_imbalanced_subset <- function(data, labels, target_active_pct, target_total, seed_offset = 0) {
  set.seed(args$seed + seed_offset)
  
  active_indices <- which(labels == 1)
  latent_indices <- which(labels == 0)
  
  n_active <- length(active_indices)
  n_latent <- length(latent_indices)
  
  # Calculate target sample sizes
  optimal <- calculate_optimal_samples(n_active, n_latent, target_active_pct, target_total)
  
  if (is.null(optimal)) {
    cat(sprintf("        Cannot create subset: need %.0f active (have %d) and %.0f latent (have %d) for %.0f%% imbalance with %d total\n",
                round(target_total * target_active_pct), n_active,
                target_total - round(target_total * target_active_pct), n_latent,
                target_active_pct * 100, target_total))
    return(NULL)  # Cannot achieve this imbalance level
  }
  
  # Sample the specified numbers
  keep_active <- sample(active_indices, optimal$active, replace = FALSE)
  keep_latent <- sample(latent_indices, optimal$latent, replace = FALSE)
  
  # Combine and sort indices
  keep_indices <- sort(c(keep_active, keep_latent))
  
  return(list(
    data = data[, keep_indices, drop = FALSE],
    labels = labels[keep_indices],
    indices = keep_indices,
    stats = list(
      original_active = n_active,
      original_latent = n_latent,
      original_total = n_active + n_latent,
      final_active = optimal$active,
      final_latent = optimal$latent,
      final_total = optimal$total,
      target_active_pct = target_active_pct,
      actual_active_pct = optimal$actual_active_pct
    )
  ))
}

# ====================================================================
# MAIN FUNCTION
# ====================================================================

create_class_imbalanced_datasets <- function() {
  
  cat("=== CLASS IMBALANCE ANALYSIS SETUP ===\n")
  cat(sprintf("Input: %s\n", args$input))
  cat(sprintf("Output directory: %s\n", args$output_dir))
  cat(sprintf("Report file: %s\n", args$report_file))
  cat(sprintf("Seed: %d\n", args$seed))
  cat(sprintf("Replicates: %d\n", args$n_replicates))
  cat("=====================================\n\n")
  
  # Load original data
  if (!file.exists(args$input)) {
    stop(sprintf("Input file not found: %s", args$input))
  }
  
  load(args$input)  # Loads dat_lst and label_lst
  
  # Validate loaded data
  if (!exists("dat_lst") || !exists("label_lst")) {
    stop("Required objects 'dat_lst' and 'label_lst' not found in input file")
  }
  
  # Store original data in separate variables that won't be modified
  dat_lst_original <- dat_lst
  label_lst_original <- label_lst
  
  # Create output directory
  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
  }
  
  # Generate all training pairs
  training_pairs <- generate_training_pairs()
  cat(sprintf("Generated %d training pairs from 6 datasets\n", length(training_pairs)))
  
  # Define imbalance levels to test
  imbalance_levels <- c(0.20, 0.30, 0.40, 0.50)  # 20%, 30%, 40%, 50% active TB
  
  # Get all datasets for testing
  all_studies <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
  
  # Initialize report data
  report_data <- data.frame()
  
  # Analyze original class distributions for all datasets
  cat("Original class distributions:\n")
  original_stats <- list()
  for (dataset in all_studies) {
    if (dataset %in% names(label_lst_original)) {
      labels <- label_lst_original[[dataset]]
      n_active <- sum(labels == 1)
      n_latent <- sum(labels == 0)
      n_total <- length(labels)
      active_pct <- n_active / n_total
      
      original_stats[[dataset]] <- list(
        dataset = dataset,
        n_active = n_active,
        n_latent = n_latent,
        n_total = n_total,
        active_pct = active_pct
      )
      
      cat(sprintf("  %s: %d active (%.1f%%), %d latent (%.1f%%), total %d\n",
                  dataset, n_active, active_pct * 100, n_latent, (1 - active_pct) * 100, n_total))
    }
  }
  
  # Process each training pair
  for (pair_idx in seq_along(training_pairs)) {
    pair <- training_pairs[[pair_idx]]
    pair_name <- paste(pair, collapse = "-")
    
    cat(sprintf("\n--- Training Pair %d/%d: %s ---\n", pair_idx, length(training_pairs), pair_name))
    cat(sprintf("  Pair elements: '%s' and '%s'\n", pair[1], pair[2]))
    cat(sprintf("  Available datasets: %s\n", paste(names(dat_lst_original), collapse = ", ")))
    
    # Check that both datasets exist
    missing_datasets <- setdiff(pair, names(dat_lst_original))
    if (length(missing_datasets) > 0) {
      cat(sprintf("Skipping pair - missing datasets: %s\n", paste(missing_datasets, collapse = ", ")))
      next
    }
    
    # Get stats for both datasets in the pair
    dataset1_stats <- original_stats[[pair[1]]]
    dataset2_stats <- original_stats[[pair[2]]]
    
    # Determine optimal arrangement with CONSISTENT sample size across all imbalance levels
    # For this specific training pair
    arrangement <- determine_optimal_arrangement_consistent(dataset1_stats, dataset2_stats, imbalance_levels)
    
    if (!arrangement$feasible) {
      cat(sprintf("Skipping pair - cannot achieve sufficient samples across all imbalance levels\n"))
      next
    }
    
    high_dataset <- arrangement$high_active_dataset
    low_dataset <- arrangement$low_active_dataset
    consistent_total <- arrangement$total_samples
    total_samples_per_dataset <- consistent_total / 2  # Split equally between two training datasets
    
    cat(sprintf("  Optimal arrangement: %s=high active, %s=low active\n", high_dataset, low_dataset))
    cat(sprintf("  CONSISTENT sample size for this pair: %.0f per dataset (across all imbalance levels)\n", total_samples_per_dataset))
    cat(sprintf("  This gives %.0f total training samples\n", consistent_total))
    
    # Test on each remaining dataset
    test_datasets <- setdiff(all_studies, pair)
    
    # For each imbalance level
    for (imbalance_pct in imbalance_levels) {
      cat(sprintf("\n  Imbalance level: %d%% active TB\n", imbalance_pct * 100))
      
      # For each test dataset
      for (test_dataset in test_datasets) {
        
        # For each replicate
        for (replicate_idx in 1:args$n_replicates) {
          scenario_name <- sprintf("%s-imbal%.0f-test%s-rep%d", 
                                  pair_name, imbalance_pct * 100, test_dataset, replicate_idx)
          
          if (replicate_idx == 1) {
            cat(sprintf("    Creating scenarios for test=%s (%d replicates)\n", test_dataset, args$n_replicates))
          }
          
          # Create imbalanced subsets for training datasets with replicate-specific seed
          seed_offset <- pair_idx * 100000 + 
                        which(imbalance_levels == imbalance_pct) * 10000 + 
                        which(test_datasets == test_dataset) * 100 + 
                        replicate_idx
          
          high_subset <- create_imbalanced_subset(
            dat_lst_original[[high_dataset]], 
            label_lst_original[[high_dataset]], 
            imbalance_pct,
            total_samples_per_dataset,
            seed_offset = seed_offset
          )
          
          low_subset <- create_imbalanced_subset(
            dat_lst_original[[low_dataset]], 
            label_lst_original[[low_dataset]], 
            1 - imbalance_pct,  # Inverse imbalance
            total_samples_per_dataset,
            seed_offset = seed_offset + 50000
          )
          
          if (is.null(high_subset) || is.null(low_subset)) {
            cat(sprintf("      Failed to create subsets for %s\n", scenario_name))
            next
          }
          
          # Create final dataset
          dat_lst_imbalanced <- list()
          label_lst_imbalanced <- list()
          
          # Add imbalanced training datasets
          dat_lst_imbalanced[[high_dataset]] <- high_subset$data
          label_lst_imbalanced[[high_dataset]] <- high_subset$labels
          
          dat_lst_imbalanced[[low_dataset]] <- low_subset$data
          label_lst_imbalanced[[low_dataset]] <- low_subset$labels
          
          # Add unchanged test dataset from original data
          dat_lst_imbalanced[[test_dataset]] <- dat_lst_original[[test_dataset]]
          label_lst_imbalanced[[test_dataset]] <- label_lst_original[[test_dataset]]
          
          # Save dataset
          output_file <- file.path(args$output_dir, sprintf("%s.RData", scenario_name))
          # Save in a local environment to avoid overwriting loop variables
          local({
            dat_lst <- dat_lst_imbalanced
            label_lst <- label_lst_imbalanced
            save(dat_lst, label_lst, file = output_file)
          })
          
          # Record statistics for report
          test_labels <- label_lst_imbalanced[[test_dataset]]
          test_active <- sum(test_labels == 1)
          test_latent <- sum(test_labels == 0)
          
          report_row <- data.frame(
            training_pair = pair_name,
            test_dataset = test_dataset,
            replicate = replicate_idx,
            high_active_dataset = high_dataset,
            low_active_dataset = low_dataset,
            target_imbalance_pct = imbalance_pct,
            samples_per_training_dataset = total_samples_per_dataset,
            high_dataset_active = high_subset$stats$final_active,
            high_dataset_latent = high_subset$stats$final_latent,
            high_dataset_total = high_subset$stats$final_total,
            high_dataset_actual_pct = high_subset$stats$actual_active_pct,
            low_dataset_active = low_subset$stats$final_active,
            low_dataset_latent = low_subset$stats$final_latent,
            low_dataset_total = low_subset$stats$final_total,
            low_dataset_actual_pct = low_subset$stats$actual_active_pct,
            test_dataset_active = test_active,
            test_dataset_latent = test_latent,
            test_dataset_total = test_active + test_latent,
            total_train_samples = high_subset$stats$final_total + low_subset$stats$final_total,
            scenario_name = scenario_name,
            output_file = output_file,
            stringsAsFactors = FALSE
          )
          
          report_data <- rbind(report_data, report_row)
        }
        
        # Print summary for this test dataset (after all replicates)
        cat(sprintf("      %s: %s=%.0f/%.0f (%.1f%%), %s=%.0f/%.0f (%.1f%%), test=%s (%d samples) × %d replicates\n",
                    gsub("_rep.*", "", scenario_name),
                    high_dataset, high_subset$stats$final_active, high_subset$stats$final_total,
                    high_subset$stats$actual_active_pct * 100,
                    low_dataset, low_subset$stats$final_active, low_subset$stats$final_total,
                    low_subset$stats$actual_active_pct * 100,
                    test_dataset, test_active + test_latent,
                    args$n_replicates))
      }
    }
  }
  
  # Save report
  write.csv(report_data, args$report_file, row.names = FALSE)
  
  cat(sprintf("\n✅ Created %d imbalanced datasets\n", nrow(report_data)))
  cat(sprintf("✅ Report saved to: %s\n", args$report_file))
  cat(sprintf("✅ Datasets saved to: %s\n", args$output_dir))
  
  # Print summary statistics
  cat("\nSummary by imbalance level:\n")
  summary_stats <- aggregate(cbind(total_train_samples, test_dataset_total) ~ target_imbalance_pct, 
                           data = report_data, FUN = function(x) c(mean = mean(x), sd = sd(x)))
  print(summary_stats)
  
  cat("\nSummary by training pair:\n")
  pair_summary <- aggregate(cbind(total_train_samples) ~ training_pair, 
                          data = report_data, FUN = function(x) c(mean = mean(x), count = length(x)))
  print(pair_summary)
  
  return(report_data)
}

# ====================================================================
# EXECUTE MAIN FUNCTION
# ====================================================================

tryCatch({
  result <- create_class_imbalanced_datasets()
}, error = function(e) {
  cat(sprintf("[ERROR] %s\n", e$message), file = stderr())
  quit(status = 1)
})