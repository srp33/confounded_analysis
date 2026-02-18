#!/usr/bin/env Rscript

# create_class_imbalanced_data.R - Create systematically class-imbalanced TB datasets
# For each pair of training datasets, determine optimal imbalance arrangement and test on remaining datasets.
#
# IMPORTANT: This script FAILS HARD if any pair × imbalance level is infeasible.
# The imbalance levels must be chosen (in config.yaml) so that all C(6,2)=15 pairs work.

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
parser$add_argument("--n-replicates", type = "integer", default = 1,
                   help = "Number of replicate random samples per scenario (default: 1)")
parser$add_argument("--imbalance-levels", type = "character", required = TRUE,
                   help = "Comma-separated imbalance proportions, e.g. '0.10,0.20,0.30,0.40,0.50'")

args <- parser$parse_args()

# Parse imbalance levels from comma-separated string
imbalance_levels <- as.numeric(strsplit(args$imbalance_levels, ",")[[1]])
cat(sprintf("Imbalance levels: %s\n", paste(imbalance_levels, collapse = ", ")))

# ====================================================================
# HELPER FUNCTIONS
# ====================================================================

# Calculate optimal sample sizes for given imbalance target, keeping total constant
calculate_optimal_samples <- function(n_active, n_latent, target_active_pct, target_total) {
  target_active <- round(target_total * target_active_pct)
  target_latent <- target_total - target_active

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
  if (consistent_total_1 >= consistent_total_2 && consistent_total_1 >= 30) {
    return(list(
      arrangement = 1,
      high_active_dataset = dataset1_stats$dataset,
      low_active_dataset = dataset2_stats$dataset,
      total_samples = consistent_total_1,
      feasible = TRUE
    ))
  } else if (consistent_total_2 >= 30) {
    return(list(
      arrangement = 2,
      high_active_dataset = dataset2_stats$dataset,
      low_active_dataset = dataset1_stats$dataset,
      total_samples = consistent_total_2,
      feasible = TRUE
    ))
  } else {
    return(list(
      feasible = FALSE,
      best_total = max(consistent_total_1, consistent_total_2)
    ))
  }
}

# Create imbalanced subset for a single dataset
create_imbalanced_subset <- function(data, labels, target_active_pct, target_total, seed_offset = 0) {
  set.seed(args$seed + seed_offset)

  active_indices <- which(labels == 1)
  latent_indices <- which(labels == 0)

  n_active <- length(active_indices)
  n_latent <- length(latent_indices)

  optimal <- calculate_optimal_samples(n_active, n_latent, target_active_pct, target_total)

  if (is.null(optimal)) {
    stop(sprintf("Cannot create subset: need %.0f active (have %d) and %.0f latent (have %d) for %.0f%% imbalance with %.0f total",
                 round(target_total * target_active_pct), n_active,
                 target_total - round(target_total * target_active_pct), n_latent,
                 target_active_pct * 100, target_total))
  }

  keep_active <- sample(active_indices, optimal$active, replace = FALSE)
  keep_latent <- sample(latent_indices, optimal$latent, replace = FALSE)
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
  cat(sprintf("Imbalance levels: %s\n", paste(imbalance_levels, collapse = ", ")))
  cat("=====================================\n\n")

  # Load original data
  if (!file.exists(args$input)) {
    stop(sprintf("Input file not found: %s", args$input))
  }

  load(args$input)  # Loads dat_lst and label_lst

  if (!exists("dat_lst") || !exists("label_lst")) {
    stop("Required objects 'dat_lst' and 'label_lst' not found in input file")
  }

  dat_lst_original <- dat_lst
  label_lst_original <- label_lst

  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
  }

  all_studies <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
  training_pairs <- combn(all_studies, 2, simplify = FALSE)
  cat(sprintf("Generated %d training pairs from %d datasets\n", length(training_pairs), length(all_studies)))

  # Compute original stats
  cat("Original class distributions:\n")
  original_stats <- list()
  for (dataset in all_studies) {
    labels <- label_lst_original[[dataset]]
    n_active <- sum(labels == 1)
    n_latent <- sum(labels == 0)
    n_total <- length(labels)
    original_stats[[dataset]] <- list(
      dataset = dataset, n_active = n_active, n_latent = n_latent,
      n_total = n_total, active_pct = n_active / n_total
    )
    cat(sprintf("  %s: %d active (%.1f%%), %d latent, total %d\n",
                dataset, n_active, n_active / n_total * 100, n_latent, n_total))
  }

  # ---- FEASIBILITY CHECK: fail hard if any pair is infeasible ----
  cat("\nChecking feasibility for all pairs...\n")
  for (pair in training_pairs) {
    arrangement <- determine_optimal_arrangement_consistent(
      original_stats[[pair[1]]], original_stats[[pair[2]]], imbalance_levels
    )
    if (!arrangement$feasible) {
      stop(sprintf(
        "INFEASIBLE: pair %s-%s cannot achieve consistent samples (best total=%.0f, need >=30) for levels [%s]. Adjust imbalance_levels in config.yaml.",
        pair[1], pair[2], arrangement$best_total,
        paste(imbalance_levels * 100, collapse = ", ")
      ))
    }
    cat(sprintf("  %s — %s: OK (per_ds=%.0f)\n",
                pair[1], pair[2], arrangement$total_samples))
  }
  cat("All pairs feasible.\n")

  report_data <- data.frame()

  # Process each training pair
  for (pair_idx in seq_along(training_pairs)) {
    pair <- training_pairs[[pair_idx]]
    pair_name <- paste(pair, collapse = "-")

    cat(sprintf("\n--- Training Pair %d/%d: %s ---\n", pair_idx, length(training_pairs), pair_name))

    arrangement <- determine_optimal_arrangement_consistent(
      original_stats[[pair[1]]], original_stats[[pair[2]]], imbalance_levels
    )

    high_dataset <- arrangement$high_active_dataset
    low_dataset <- arrangement$low_active_dataset
    consistent_total <- arrangement$total_samples
    total_samples_per_dataset <- consistent_total  # total_samples is already per-dataset

    cat(sprintf("  Arrangement: %s=high active, %s=low active\n", high_dataset, low_dataset))
    cat(sprintf("  Consistent sample size: %.0f per dataset\n", total_samples_per_dataset))

    test_datasets <- setdiff(all_studies, pair)

    for (imbalance_pct in imbalance_levels) {
      cat(sprintf("\n  Imbalance level: %.0f%% active TB\n", imbalance_pct * 100))

      for (test_dataset in test_datasets) {
        for (replicate_idx in 1:args$n_replicates) {
          scenario_name <- sprintf("%s-imbal%.0f-test%s-rep%d",
                                  pair_name, imbalance_pct * 100, test_dataset, replicate_idx)

          seed_offset <- pair_idx * 100000 +
                        which(imbalance_levels == imbalance_pct) * 10000 +
                        which(test_datasets == test_dataset) * 100 +
                        replicate_idx

          high_subset <- create_imbalanced_subset(
            dat_lst_original[[high_dataset]],
            label_lst_original[[high_dataset]],
            imbalance_pct, total_samples_per_dataset,
            seed_offset = seed_offset
          )

          low_subset <- create_imbalanced_subset(
            dat_lst_original[[low_dataset]],
            label_lst_original[[low_dataset]],
            1 - imbalance_pct, total_samples_per_dataset,
            seed_offset = seed_offset + 50000
          )

          # Build the 3-dataset list (2 training + 1 test)
          dat_lst_imbalanced <- list()
          label_lst_imbalanced <- list()
          dat_lst_imbalanced[[high_dataset]] <- high_subset$data
          label_lst_imbalanced[[high_dataset]] <- high_subset$labels
          dat_lst_imbalanced[[low_dataset]] <- low_subset$data
          label_lst_imbalanced[[low_dataset]] <- low_subset$labels
          dat_lst_imbalanced[[test_dataset]] <- dat_lst_original[[test_dataset]]
          label_lst_imbalanced[[test_dataset]] <- label_lst_original[[test_dataset]]

          output_file <- file.path(args$output_dir, sprintf("%s.RData", scenario_name))
          local({
            dat_lst <- dat_lst_imbalanced
            label_lst <- label_lst_imbalanced
            save(dat_lst, label_lst, file = output_file)
          })

          test_labels <- label_lst_imbalanced[[test_dataset]]
          test_active <- sum(test_labels == 1)
          test_latent <- sum(test_labels == 0)

          report_row <- data.frame(
            training_pair = paste(sort(c(high_dataset, low_dataset)), collapse = "-"),
            test_dataset = test_dataset,
            replicate = replicate_idx,
            train_dataset_1 = high_dataset,
            train_dataset_2 = low_dataset,
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

        cat(sprintf("      %s: %s=%.0f/%.0f (%.1f%%), %s=%.0f/%.0f (%.1f%%), test=%s (%.0f samples)\n",
                    scenario_name,
                    high_dataset, high_subset$stats$final_active, high_subset$stats$final_total,
                    high_subset$stats$actual_active_pct * 100,
                    low_dataset, low_subset$stats$final_active, low_subset$stats$final_total,
                    low_subset$stats$actual_active_pct * 100,
                    test_dataset, test_active + test_latent))
      }
    }
  }

  write.csv(report_data, args$report_file, row.names = FALSE)

  cat(sprintf("\nCreated %d imbalanced datasets\n", nrow(report_data)))
  cat(sprintf("Report saved to: %s\n", args$report_file))
  cat(sprintf("Datasets saved to: %s\n", args$output_dir))
}

# ====================================================================
# EXECUTE
# ====================================================================

tryCatch({
  result <- create_class_imbalanced_datasets()
}, error = function(e) {
  cat(sprintf("[ERROR] %s\n", e$message), file = stderr())
  quit(status = 1)
})
