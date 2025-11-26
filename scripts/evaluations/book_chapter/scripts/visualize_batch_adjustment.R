#!/usr/bin/env Rscript

# visualize_batch_adjustment.R
# Generate PCA, LDA, and UMAP visualizations for batch adjustment methods
# Follows single-responsibility principle with modular functions

# Suppress warnings and messages
options(warn = -1)
suppressMessages(suppressWarnings({
  required_packages <- c("argparse", "ggplot2", "dplyr", "purrr", "umap", "MASS", 
                        "genefilter", "sva", "batchelor", "SummarizedExperiment")
  sapply(required_packages, require, character.only = TRUE, quietly = TRUE)
}))

# ====================================================================
# COMMAND-LINE ARGUMENT PARSING
# ====================================================================

parser <- ArgumentParser(description = "Visualize batch adjustment effects using PCA, LDA, and UMAP")

parser$add_argument("--adjuster", type = "character", required = TRUE,
                   help = "Batch correction method: unadjusted, combat, combat_sup, or mnn")
parser$add_argument("--num-datasets", type = "integer", required = TRUE,
                   help = "Number of datasets to include: 3, 4, 5, or 6")
parser$add_argument("--test-study", type = "character", required = TRUE,
                   help = "Study to use as test set (e.g., 'India', 'USA', 'Africa')")
parser$add_argument("--output-dir", type = "character", required = TRUE,
                   help = "Output directory for visualization files")
parser$add_argument("--reduce", type = "integer", default = 0,
                   help = "Number of dimensions to reduce to (default: 0, no reduction)")

args <- parser$parse_args()

# Validate arguments
valid_adjusters <- c("unadjusted", "combat", "combat_sup", "mnn")
valid_num_datasets <- c(3, 4, 5, 6)

if (!args$adjuster %in% valid_adjusters) {
  stop(sprintf("Invalid adjuster '%s'. Must be one of: %s", 
               args$adjuster, paste(valid_adjusters, collapse = ", ")))
}

if (!args$num_datasets %in% valid_num_datasets) {
  stop(sprintf("Invalid num-datasets '%d'. Must be one of: %s", 
               args$num_datasets, paste(valid_num_datasets, collapse = ", ")))
}

# Extract parameters
adjuster <- args$adjuster
num_datasets <- args$num_datasets
test_study <- args$test_study
output_dir <- args$output_dir
reduce <- args$reduce

# Create output directory structure
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "pca"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "lda"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "umap"), recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Starting visualization: adjuster=%s, num_datasets=%d, test_study=%s\n", 
            adjuster, num_datasets, test_study))

# Load helper functions
source("scripts/helper.R")


# ====================================================================
# DATA LOADING AND PREPARATION (Single Responsibility: Data I/O)
# ====================================================================

#' Load and filter TB real data
#' @param n_studies Number of studies to include
#' @return List with filtered data and labels
load_and_filter_data <- function(n_studies) {
  data_path <- "data/TB_real_data.RData"
  if (!file.exists(data_path)) {
    stop(sprintf("Data file not found: %s", data_path))
  }
  
  load(data_path)
  
  all_studies <- c("GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M")
  selected_studies <- all_studies[1:n_studies]
  
  dat_lst <- dat_lst[selected_studies]
  label_lst <- label_lst[selected_studies]
  
  cat(sprintf("Loaded %d studies: %s\n", 
              n_studies, paste(selected_studies, collapse = ", ")))
  
  list(
    dat_lst = dat_lst,
    label_lst = label_lst,
    study_names = selected_studies
  )
}

#' Validate test study selection
#' @param study_names Vector of study names
#' @param test_study Name of test study
#' @return Name of test study (validated)
validate_test_study <- function(study_names, test_study) {
  if (!test_study %in% study_names) {
    stop(sprintf("Test study '%s' not found in available studies: %s",
                 test_study, paste(study_names, collapse = ", ")))
  }
  
  cat(sprintf("Using test study: %s\n", test_study))
  cat(sprintf("Training studies: %s\n", 
              paste(setdiff(study_names, test_study), collapse = ", ")))
  
  test_study
}

#' Prepare training and test datasets
#' @param dat_lst List of data matrices
#' @param label_lst List of label vectors
#' @param test_name Name of test study
#' @param study_names All study names
#' @return List with prepared datasets and metadata
prepare_train_test_split <- function(dat_lst, label_lst, test_name, study_names) {
  train_names <- setdiff(study_names, test_name)
  
  # Combine training data
  dat <- do.call(cbind, dat_lst[train_names])
  
  # Create batch assignments using study names (not numeric IDs)
  batch <- rep(train_names, times = sapply(dat_lst[train_names], ncol))
  batches_ind <- lapply(train_names, function(name) which(batch == name))
  names(batches_ind) <- train_names
  
  group <- do.call(c, label_lst[train_names])
  
  # Test data
  dat_test <- dat_lst[[test_name]]
  group_test <- label_lst[[test_name]]
  
  cat(sprintf("Training: %d samples from %d batches\n", ncol(dat), length(train_names)))
  cat(sprintf("  Batches: %s\n", paste(train_names, collapse = ", ")))
  cat(sprintf("Test: %d samples from %s\n", ncol(dat_test), test_name))
  
  list(
    dat = dat,
    batch = batch,
    batches_ind = batches_ind,
    batch_names = train_names,
    group = group,
    dat_test = dat_test,
    group_test = group_test,
    test_name = test_name
  )
}

#' Reduce features to top N most variable genes
#' @param dat Training data matrix
#' @param dat_test Test data matrix
#' @param n_genes Number of genes to select
#' @return List with reduced datasets
reduce_features <- function(dat, dat_test, n_genes = 1000) {
  genes_sel_idx <- order(rowVars(dat), decreasing = TRUE)[1:n_genes]
  
  cat(sprintf("Selected top %d most variable genes\n", n_genes))
  
  list(
    dat = dat[genes_sel_idx, ],
    dat_test = dat_test[genes_sel_idx, ]
  )
}


# ====================================================================
# BATCH CORRECTION (Single Responsibility: Data Transformation)
# ====================================================================

#' Apply batch correction method
#' @param dat Training data matrix
#' @param batch Batch assignments (character vector with study names)
#' @param group Sample labels
#' @param dat_test Test data matrix
#' @param method Correction method: "unadjusted", "combat", or "mnn"
#' @return List with corrected training and test data
apply_batch_correction <- function(dat, batch, group, dat_test, method) {
  cat(sprintf("Applying batch correction: %s\n", method))
  
  if (method == "unadjusted") {
    return(list(
      dat_corrected = dat,
      dat_test_corrected = dat_test
    ))
  }
  
  if (method == "combat") {
    library(sva, quietly = TRUE)
    
    # ComBat correction without labels (unsupervised)
    # Step 1: Correct batch effects within training data without using labels
    dat_corrected <- ComBat(dat, batch=batch, mod=NULL)
    
    # Step 2: Adjust test data to match corrected training distribution
    combined_dat <- cbind(dat_corrected, dat_test)
    ref_batch_id <- 1
    test_batch_id <- 2
    combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), 
                       rep(test_batch_id, ncol(dat_test)))
    
    combat_combined <- ComBat(combined_dat, batch=combined_batch, 
                             mod=NULL, ref.batch=ref_batch_id)
    
    dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  if (method == "combat_sup") {
    library(sva, quietly = TRUE)
    
    # ComBat correction with labels (supervised)
    # Step 1: Correct batch effects within training data while preserving biological signal
    dat_corrected <- ComBat(dat, batch=batch, mod=model.matrix(~group))
    
    # Step 2: Adjust test data to match corrected training distribution
    combined_dat <- cbind(dat_corrected, dat_test)
    ref_batch_id <- 1
    test_batch_id <- 2
    combined_batch <- c(rep(ref_batch_id, ncol(dat_corrected)), 
                       rep(test_batch_id, ncol(dat_test)))
    
    combat_combined <- ComBat(combined_dat, batch=combined_batch, 
                             mod=NULL, ref.batch=ref_batch_id)
    
    dat_test_corrected <- combat_combined[, (ncol(dat_corrected) + 1):ncol(combat_combined)]
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  if (method == "mnn") {
    library(batchelor, quietly = TRUE)
    library(SummarizedExperiment, quietly = TRUE)
    
    # Combine data
    combined_dat <- cbind(dat, dat_test)
    # Test set gets a unique batch ID
    test_id <- "TEST_SET"
    combined_batch <- c(batch, rep(test_id, ncol(dat_test)))
    
    # Determine merge order: training batches by size, test last
    u_batches <- unique(batch)
    train_sizes <- table(batch)[u_batches]
    train_ord <- order(train_sizes, decreasing = TRUE)
    merge_ord <- c(u_batches[train_ord], test_id)
    
    # Apply MNN correction
    mnn_object <- batchelor::mnnCorrect(
      combined_dat, 
      batch = combined_batch, 
      merge.order = merge_ord
    )
    mnn_matrix <- SummarizedExperiment::assay(mnn_object)
    
    # Split back
    dat_corrected <- mnn_matrix[, 1:ncol(dat)]
    dat_test_corrected <- mnn_matrix[, (ncol(dat) + 1):ncol(mnn_matrix)]
    
    return(list(
      dat_corrected = dat_corrected,
      dat_test_corrected = dat_test_corrected
    ))
  }
  
  stop(sprintf("Unknown batch correction method: %s", method))
}

#' Normalize data within batches
#' @param dat Data matrix
#' @param batch Batch assignments (character vector)
#' @param batch_names Batch identifiers (character vector)
#' @return Normalized data matrix
normalize_within_batches <- function(dat, batch, batch_names) {
  # Create empty matrix with same dimensions
  dat_norm <- matrix(NA, nrow = nrow(dat), ncol = ncol(dat))
  
  # Normalize each batch separately
  for (batch_name in batch_names) {
    batch_indices <- batch == batch_name
    batch_data <- dat[, batch_indices, drop = FALSE]
    
    # Normalize this batch
    normalized_batch <- normalizeData(batch_data)
    
    # Remove dimnames to avoid assignment issues
    dimnames(normalized_batch) <- NULL
    
    # Assign to output matrix
    dat_norm[, batch_indices] <- normalized_batch
  }
  
  # Set dimnames after filling the matrix
  rownames(dat_norm) <- rownames(dat)
  colnames(dat_norm) <- colnames(dat)
  
  dat_norm
}


# ====================================================================
# DIMENSIONALITY REDUCTION (Single Responsibility: Transformation)
# ====================================================================

#' Perform PCA on training data and project test data
#' @param dat_train Training data matrix (genes x samples)
#' @param dat_test Test data matrix (genes x samples)
#' @return List with PCA results
compute_pca <- function(dat_train, dat_test) {
  cat("Computing PCA (fit on training, project test)...\n")
  
  # Fit PCA on training data only
  # Transpose: PCA expects samples x features
  pca_fit <- prcomp(t(dat_train), center = TRUE, scale. = TRUE)
  
  # Project training data
  train_coords <- pca_fit$x[, 1:2]
  
  # Project test data using the same transformation
  test_coords <- predict(pca_fit, newdata = t(dat_test))[, 1:2]
  
  # Combine coordinates
  combined_coords <- rbind(train_coords, test_coords)
  
  # Calculate variance explained
  var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
  
  list(
    coords = combined_coords,
    var_explained = var_explained[1:2],
    method = "PCA"
  )
}

#' Perform LDA on training data and project test data
#' @param dat_train Training data matrix (genes x samples)
#' @param dat_test Test data matrix (genes x samples)
#' @param labels_train Training sample labels
#' @param batch_train Training batch assignments
#' @return List with LDA results
compute_lda <- function(dat_train, dat_test, labels_train, batch_train) {
  cat("Computing LDA (fit on training, project test)...\n")
  
  library(MASS, quietly = TRUE)
  
  # Transpose: LDA expects samples x features
  dat_train_t <- t(dat_train)
  dat_test_t <- t(dat_test)
  
  # Create combined grouping for training: label + batch
  # This creates groups like "0_1", "0_2", "1_1", "1_2" etc.
  combined_group_train <- paste(labels_train, batch_train, sep = "_")
  
  # LDA requires at least 2 classes
  if (length(unique(combined_group_train)) < 2) {
    warning("LDA requires at least 2 classes. Skipping.")
    return(NULL)
  }
  
  # Fit LDA on training data only
  lda_fit <- lda(dat_train_t, grouping = as.factor(combined_group_train))
  
  # Project training data
  train_coords <- predict(lda_fit, dat_train_t)$x
  
  # Project test data using the same transformation
  test_coords <- predict(lda_fit, dat_test_t)$x
  
  # Combine coordinates
  combined_coords <- rbind(train_coords, test_coords)
  
  # Handle 1D case (only 1 discriminant)
  if (is.null(dim(combined_coords))) {
    combined_coords <- cbind(combined_coords, rep(0, length(combined_coords)))
    colnames(combined_coords) <- c("LD1", "LD2")
  } else if (ncol(combined_coords) == 1) {
    combined_coords <- cbind(combined_coords, rep(0, nrow(combined_coords)))
    colnames(combined_coords) <- c("LD1", "LD2")
  }
  
  list(
    coords = combined_coords[, 1:2],
    method = "LDA"
  )
}

#' Perform UMAP on data
#' @param dat Data matrix (genes x samples)
#' @return List with UMAP results
compute_umap <- function(dat) {
  cat("Computing UMAP...\n")
  
  library(umap, quietly = TRUE)
  
  # Use fixed seed for reproducibility
  set.seed(42)
  
  # Transpose: UMAP expects samples x features
  # Use custom config for stability
  custom_config <- umap.defaults
  custom_config$random_state <- 42
  custom_config$n_neighbors <- min(15, ncol(dat) - 1)
  
  umap_result <- umap(t(dat), config = custom_config)
  
  list(
    coords = umap_result$layout,
    method = "UMAP"
  )
}


# ====================================================================
# VISUALIZATION (Single Responsibility: Plotting)
# ====================================================================

#' Create a dimensionality reduction plot
#' @param coords 2D coordinates matrix
#' @param batch Batch assignments
#' @param labels Sample labels
#' @param method Method name (PCA, LDA, UMAP)
#' @param var_explained Variance explained (for PCA)
#' @param title Plot title
#' @return ggplot object
create_reduction_plot <- function(coords, batch, labels, dataset_type, method, 
                                 var_explained = NULL, title = NULL) {
  
  # Create data frame for plotting
  plot_df <- data.frame(
    Dim1 = coords[, 1],
    Dim2 = coords[, 2],
    Batch = as.factor(batch),
    Label = as.factor(labels),
    DatasetType = as.factor(dataset_type)
  )
  
  # Calculate sample counts per batch
  batch_counts <- table(batch)
  
  # Create new labels with sample counts: "StudyName (n=XX)"
  batch_labels <- sapply(levels(plot_df$Batch), function(b) {
    sprintf("%s (n=%d)", b, batch_counts[b])
  })
  names(batch_labels) <- levels(plot_df$Batch)
  
  # Determine axis labels
  if (method == "PCA" && !is.null(var_explained)) {
    xlab <- sprintf("PC1 (%.1f%%)", var_explained[1] * 100)
    ylab <- sprintf("PC2 (%.1f%%)", var_explained[2] * 100)
  } else if (method == "LDA") {
    xlab <- "LD1"
    ylab <- "LD2"
  } else if (method == "UMAP") {
    xlab <- "UMAP1"
    ylab <- "UMAP2"
  } else {
    xlab <- "Dimension 1"
    ylab <- "Dimension 2"
  }
  
  # Create base plot
  p <- ggplot(plot_df, aes(x = Dim1, y = Dim2))
  
  # Add points: color by batch, shape by label, border by dataset type
  # Use larger size for test set to make it stand out
  p <- p + geom_point(aes(color = Batch, shape = Label, size = DatasetType, alpha = DatasetType))
  p <- p + scale_size_manual(values = c("Training" = 2.5, "Test" = 3.2), guide = "none")
  p <- p + scale_alpha_manual(values = c("Training" = 0.6, "Test" = 0.75))
  
  # Styling
  p <- p + theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold")
    ) +
    labs(
      x = xlab,
      y = ylab,
      title = title,
      color = "Study",
      shape = "Status"
    ) +
    scale_shape_manual(values = c(16, 17), labels = c("Non-Progressing", "Progressing")) +
    # Use a colorblind-friendly palette with enough colors
    scale_color_brewer(palette = "Set2", labels = batch_labels)
  
  p
}

#' Save plot to file
#' @param plot ggplot object
#' @param filepath Output file path
#' @param width Plot width in inches
#' @param height Plot height in inches
save_plot <- function(plot, filepath, width = 8, height = 6) {
  ggsave(
    filename = filepath,
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
  cat(sprintf("Saved plot: %s\n", filepath))
}

#' Create and save all three reduction plots
#' @param dat_train Training data matrix
#' @param dat_test Test data matrix
#' @param batch_train Training batch assignments
#' @param labels_train Training sample labels
#' @param labels_test Test sample labels
#' @param test_study Test study name
#' @param adjuster Adjuster name
#' @param num_datasets Number of datasets
#' @param output_dir Output directory
create_all_visualizations <- function(dat_train, dat_test, batch_train, 
                                     labels_train, labels_test, test_study,
                                     adjuster, num_datasets, output_dir) {
  
  # Combine training and test data
  dat_combined <- cbind(dat_train, dat_test)
  
  # Create combined batch vector (test set gets its own label)
  batch_combined <- c(batch_train, rep(test_study, ncol(dat_test)))
  
  # Create combined labels
  labels_combined <- c(labels_train, labels_test)
  
  # Create dataset type indicator (train vs test)
  dataset_type <- c(rep("Training", ncol(dat_train)), rep("Test", ncol(dat_test)))
  
  base_name <- sprintf("%s_n%d_test%s", adjuster, num_datasets, test_study)
  
  # PCA (fit on training, project test)
  pca_result <- compute_pca(dat_train, dat_test)
  pca_plot <- create_reduction_plot(
    coords = pca_result$coords,
    batch = batch_combined,
    labels = labels_combined,
    dataset_type = dataset_type,
    method = "PCA",
    var_explained = pca_result$var_explained,
    title = sprintf("PCA: %s (n=%d, test=%s)", 
                   tools::toTitleCase(adjuster), num_datasets, test_study)
  )
  save_plot(
    pca_plot,
    file.path(output_dir, "pca", paste0(base_name, ".png"))
  )
  
  # LDA (fit on training, project test)
  lda_result <- compute_lda(dat_train, dat_test, labels_train, batch_train)
  if (!is.null(lda_result)) {
    lda_plot <- create_reduction_plot(
      coords = lda_result$coords,
      batch = batch_combined,
      labels = labels_combined,
      dataset_type = dataset_type,
      method = "LDA",
      title = sprintf("LDA: %s (n=%d, test=%s)", 
                     tools::toTitleCase(adjuster), num_datasets, test_study)
    )
    save_plot(
      lda_plot,
      file.path(output_dir, "lda", paste0(base_name, ".png"))
    )
  }
  
  # UMAP (still uses combined data - cannot project separately)
  dat_combined <- cbind(dat_train, dat_test)
  umap_result <- compute_umap(dat_combined)
  umap_plot <- create_reduction_plot(
    coords = umap_result$coords,
    batch = batch_combined,
    labels = labels_combined,
    dataset_type = dataset_type,
    method = "UMAP",
    title = sprintf("UMAP: %s (n=%d, test=%s)", 
                   tools::toTitleCase(adjuster), num_datasets, test_study)
  )
  save_plot(
    umap_plot,
    file.path(output_dir, "umap", paste0(base_name, ".png"))
  )
  
  cat("All visualizations created successfully\n")
}


# ====================================================================
# MAIN EXECUTION
# ====================================================================

main <- function() {
  tryCatch({
    cat("=== BATCH ADJUSTMENT VISUALIZATION ===\n")
    cat(sprintf("Adjuster: %s\n", adjuster))
    cat(sprintf("Num datasets: %d\n", num_datasets))
    cat(sprintf("Test study: %s\n", test_study))
    cat(sprintf("Output directory: %s\n", output_dir))
    cat("======================================\n\n")
    
    # Load and prepare data
    cat("Step 1: Loading data...\n")
    data <- load_and_filter_data(num_datasets)
    
    cat("\nStep 2: Validating test study...\n")
    test_name <- validate_test_study(data$study_names, test_study)
    
    cat("\nStep 3: Preparing train/test split...\n")
    datasets <- prepare_train_test_split(
      data$dat_lst, 
      data$label_lst, 
      test_name, 
      data$study_names
    )
    
    cat("\nStep 4: Reducing features...\n")
    if (reduce == 0) {
      reduced = list(dat = datasets$dat, dat_test = datasets$dat_test)
    }
    else {
      reduced <- reduce_features(datasets$dat, datasets$dat_test, n_genes = reduce)
    }
    
    cat("\nStep 5: Applying batch correction...\n")
    corrected <- apply_batch_correction(
      reduced$dat,
      datasets$batch,
      datasets$group,
      reduced$dat_test,
      adjuster
    )
    
    cat("\nStep 6: Normalizing data...\n")
    # Normalize training data within batches
    dat_train_norm <- normalize_within_batches(
      corrected$dat_corrected,
      datasets$batch,
      datasets$batch_names
    )
    
    # Normalize test data (as a single batch)
    dat_test_norm <- normalizeData(corrected$dat_test_corrected)
    
    cat("\nStep 7: Creating visualizations...\n")
    create_all_visualizations(
      dat_train = dat_train_norm,
      dat_test = dat_test_norm,
      batch_train = datasets$batch,
      labels_train = datasets$group,
      labels_test = datasets$group_test,
      test_study = test_study,
      adjuster = adjuster,
      num_datasets = num_datasets,
      output_dir = output_dir
    )
    
    cat("\n=== VISUALIZATION COMPLETE ===\n")
    cat(sprintf("Output directory: %s\n", output_dir))
    cat("Files organized by method: pca/, lda/, umap/\n")
    
    return(0)
    
  }, error = function(e) {
    cat(sprintf("\n[ERROR] Visualization failed: %s\n", e$message), file = stderr())
    cat(sprintf("[ERROR] Traceback:\n"), file = stderr())
    traceback_lines <- capture.output(traceback())
    for (line in traceback_lines) {
      cat(sprintf("[ERROR] %s\n", line), file = stderr())
    }
    return(1)
  })
}

# Execute main function
exit_code <- main()
quit(status = exit_code)
