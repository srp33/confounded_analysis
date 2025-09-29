# Load dependencies -----------------------------------------------------------

# Force OpenMP to run in single-threaded mode.
# This prevents errors when running multiple R processes in parallel.
Sys.setenv(OMP_NUM_THREADS = 1)

# Suppress package startup messages for a cleaner console output.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(vroom)
  library(stringr)
  library(argparse)
  library(sva)
  library(limma)
  library(Seurat)
  library(batchelor)
  library(rliger)
  library(fairadapt)
  library(future)
  library(huge)
  library(preprocessCore)
  library(data.table)
})

source("/scripts/adjust/gmm_adjust.R")
source("/scripts/adjust/gmm_global_simple.R")
source("/scripts/adjust/gmm_adjust_streamlined.R")

get_allocated_cores <- function() {
  # Check SLURM-provided environment variables first
  slurm_cpus_per_task <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
  slurm_ntasks <- as.integer(Sys.getenv("SLURM_NTASKS", unset = NA))
  slurm_cpus_on_node <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", unset = NA))

  if (!is.na(slurm_cpus_per_task) && slurm_cpus_per_task > 0) {
    return(slurm_cpus_per_task)
  }
  if (!is.na(slurm_ntasks) && slurm_ntasks > 0) {
    return(slurm_ntasks)
  }
  if (!is.na(slurm_cpus_on_node) && slurm_cpus_on_node > 0) {
    return(slurm_cpus_on_node)
  }

  # Fall back to full machine
  return(parallel::detectCores(logical = TRUE))
}

# Helper Functions ------------------------------------------------------------

transpose_essential <- function(gene_df) {
  #' Transpose a quantitative data frame robustly.
  #'
  #' This multi-step process ensures sample names (row names) are
  #' correctly preserved during the conversion from a data frame to a
  #' transposed matrix (features x samples).
  #'
  #' @param gene_df A data frame with samples as rows and features as columns.
  #' @return A transposed matrix with features as rows and samples as columns.
  
  df_version <- as.data.frame(gene_df)
  mat_untransposed <- as.matrix(df_version)
  rownames(mat_untransposed) <- rownames(df_version)
  mat_quantitative <- t(mat_untransposed)
  
  return(mat_quantitative)
}


ComBat_ignore_nonvariance <- function(matrix_, batch, design) {
  #' Wrap ComBat to handle features with zero variance.
  #'
  #' ComBat fails if any feature (row) has zero variance across all samples.
  #' This function identifies such features, excludes them from the ComBat
  #' adjustment, and then adds them back to the result.
  #'
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param design The design matrix.
  #' @return The adjusted matrix.
  
  varying_row_mask <- apply(matrix_, 1, function(x) { length(unique(x)) > 1 })
  
  if (sum(varying_row_mask) < nrow(matrix_)) {
    message(
      sprintf(
        "Found %d features with zero variance. These will be ignored by ComBat.",
        nrow(matrix_) - sum(varying_row_mask)
      )
    )
  }
  
  matrix_[varying_row_mask, ] <- ComBat(matrix_[varying_row_mask, ], batch, mod = design, prior.plots = FALSE)
  return(matrix_)
}


create_design_matrix <- function(categorical_df, columns_to_use = NULL, use_all = FALSE) {
  #' Create a design matrix from a data frame of categorical variables.
  #'
  #' @param categorical_df A data frame with samples as rows and categorical metadata as columns.
  #' @param columns_to_use A character vector of specific columns to include in the model.
  #' @param use_all A boolean indicating whether to use all columns in the data frame.
  #' @return A design matrix.
  
  # Ensure use_all is a proper logical value
  if (is.null(use_all) || length(use_all) == 0) {
    use_all <- FALSE
  }
  
  if (!is.null(columns_to_use)) {
    if (!all(columns_to_use %in% colnames(categorical_df))) {
      stop("One or more specified columns for the design matrix were not found in the metadata.")
    }
    design_df <- categorical_df[, columns_to_use, drop = FALSE]
    message("Creating design matrix from specified columns: ", paste(columns_to_use, collapse = ", "))
  } else if (use_all) {
    design_df <- categorical_df
    message("Creating design matrix from all available categorical variables.")
  } else {
    message("Creating design matrix with intercept only.")
    return(matrix(1, nrow = nrow(categorical_df), ncol = 1, dimnames = list(NULL, "Intercept")))
  }
  
  if (ncol(design_df) == 0) {
    message("No columns selected for design matrix. Returning intercept-only model.")
    return(matrix(1, nrow = nrow(categorical_df), ncol = 1, dimnames = list(NULL, "Intercept")))
  }
  
  design <- model.matrix(~ ., data = design_df)
  return(design)
}


prep_seurat_like <- function(df_, batch, data_are_counts) {
  #' Prepare a Seurat object for downstream adjustments.
  #'
  #' This helper function handles common preprocessing steps:
  #' 1. Stores original feature (gene) and sample names.
  #' 2. Sanitizes names to be compatible with Seurat (e.g., replaces underscores).
  #' 3. Creates a Seurat object.
  #' 4. Normalizes the data using appropriate methods.
  #'
  #' @param df_ The data matrix (features x samples).
  #' @param batch The batch variable vector.
  #' @param data_are_counts Logical, TRUE if data is raw counts.
  #' @return A list containing the Seurat object and original names.
  
  original_feature_names <- rownames(df_)
  original_sample_names <- colnames(df_)
  
  sanitized_feature_names <- make.unique(gsub("_", "-", original_feature_names))
  sanitized_sample_names <- make.unique(gsub("_", "-", original_sample_names))
  
  df_copy <- df_
  rownames(df_copy) <- sanitized_feature_names
  colnames(df_copy) <- sanitized_sample_names
  
  meta <- data.frame(Batch = batch)
  rownames(meta) <- sanitized_sample_names
  
  seurat_obj <- CreateSeuratObject(counts = df_copy, meta.data = meta)
  
  if (data_are_counts) {
    message("Data appears to be raw counts. Normalizing using LogNormalize.")
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  } else {
    message("Data appears pre-normalized. Setting 'data' layer directly from input.")
    seurat_obj <- SetAssayData(seurat_obj, layer = "data", new.data = df_copy)
  }
  
  return(list(
    obj = seurat_obj,
    orig_features = original_feature_names,
    orig_samples = original_sample_names
  ))
}


restore_names <- function(matrix, prep_list) {
  #' Restore original feature and sample names to a matrix.
  #' @param matrix The matrix with sanitized names.
  #' @param prep_list The list returned by `prep_seurat_like`.
  #' @return The matrix with original names.
  
  feature_name_map <- data.frame(original = prep_list$orig_features, sanitized = rownames(prep_list$obj))
  sample_name_map <- data.frame(original = prep_list$orig_samples, sanitized = colnames(prep_list$obj))
  
  original_rownames <- feature_name_map$original[match(rownames(matrix), feature_name_map$sanitized)]
  original_colnames <- sample_name_map$original[match(colnames(matrix), sample_name_map$sanitized)]
  
  rownames(matrix) <- original_rownames
  colnames(matrix) <- original_colnames
  
  return(matrix)
}


# Adjustment Functions --------------------------------------------------------

adjust_min_mean <- function(matrix_, batch, ..., debug = FALSE) {
  #' Adjust matrix by matching minimum and mean values across batches.
  #' Assumes the batch vector contains no NA values.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @return The adjusted matrix.
  
  message("Adjusting data by matching minimum and mean values across batches.")
  
  batch_levels <- unique(batch)
  
  global_mins <- apply(matrix_, 1, min)
  global_means <- apply(matrix_, 1, mean)
  
  adjusted_matrix <- matrix_
  
  for (b in batch_levels) {
    batch_indices <- which(batch == b)
    if (length(batch_indices) > 0) {
      batch_data <- matrix_[, batch_indices, drop = FALSE]
      
      batch_mins <- apply(batch_data, 1, min)
      batch_means <- rowMeans(batch_data)
      
      adjusted_matrix[, batch_indices] <- global_mins + (batch_data - batch_mins) * (global_means - global_mins) / (batch_means - batch_mins)
      
      if (debug) {
        message("DEBUG: Adjusted batch ", b, ": ", length(batch_indices), " samples")
        adjusted_batch_means <- rowMeans(adjusted_matrix[, batch_indices, drop = FALSE])
        mean_diff <- mean(abs(adjusted_batch_means - global_means))
        message("DEBUG: Mean absolute difference from global means: ", round(mean_diff, 6))
      }
    }
  }
  
  return(adjusted_matrix)
}


adjust_combat <- function(matrix_, batch, design, data_are_counts, debug = FALSE) {
  #' Adjust matrix using ComBat or ComBat_seq.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param design The design matrix.
  #' @param data_are_counts If TRUE, use ComBat_seq for count data.
  #' @return The adjusted matrix.
  
  if (data_are_counts) {
    message("Using ComBat_seq for count data.")
    return(ComBat_seq(matrix_, batch, covar_mod = design))
  } else {
    message("Using ComBat for continuous data.")
    return(ComBat_ignore_nonvariance(matrix_, batch, design))
  }
}


adjust_limma <- function(matrix_, batch, design, ..., debug = FALSE) {
  #' Adjust matrix using limma's removeBatchEffect.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param design The design matrix.
  #' @return The adjusted matrix.
  
  message("Adjusting data with limma::removeBatchEffect.")
  return(limma::removeBatchEffect(matrix_, batch = batch, design = design))
}


adjust_quantile <- function(matrix_, batch, debug = FALSE) {
  #' Adjust matrix using quantile normalization with a global reference.
  #' Assumes the batch vector contains no NA values.
  #' @param matrix_ The matrix to adjust (genes x samples).
  #' @param batch A vector identifying the batch for each sample (column) in matrix_.
  #' @param debug A logical flag to enable verbose debugging messages.
  #' @return The adjusted matrix.
  
  # Transpose because preprocessCore normalizes columns (we want to normalize genes).
  matrix_ <- t(matrix_)
  
  message("Determining global target distribution from the entire dataset.")
  global_target_distribution <- preprocessCore::normalize.quantiles.determine.target(matrix_)
  
  if (debug) {
    message("DEBUG: Global target distribution vector generated. Length: ", length(global_target_distribution))
  }
  
  result_matrix <- matrix_
  batch_levels <- unique(batch)
  
  for (b in batch_levels) {
    message("Normalizing batch: '", b, "'")
    rows_for_batch <- which(batch == b)
    
    normalized_batch_matrix <- preprocessCore::normalize.quantiles.use.target(
      x = matrix_[rows_for_batch, , drop = FALSE],
      target = global_target_distribution
    )
    
    result_matrix[rows_for_batch, ] <- normalized_batch_matrix
  }
  
  # Transpose back to original orientation (features x samples).
  return(t(result_matrix))
}


adjust_npn <- function(matrix_, batch, debug = FALSE) {
  #' Adjust matrix using Nonparanormal (NPN) transformation.
  #' If batch is NULL, the entire matrix is adjusted at once.
  #' Assumes the batch vector contains no NA values.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @return The adjusted matrix.
  
  if (is.null(batch)) {
    # If no batch is provided, adjust the whole matrix.
    message("Batch is NULL. Adjusting entire matrix with NPN transformation.")
    
    # Transpose to (samples x features) for huge.npn.
    matrix_t <- t(matrix_)
    
    npn_transformed_t <- huge::huge.npn(matrix_t, verbose = FALSE)
    
    return(t(npn_transformed_t))
    
  } else {
    message("Adjusting using Nonparanormal (NPN) transformation by batch.")
    
    # Split the matrix by batch.
    batch_levels <- unique(batch)
    matrix_by_batch <- list()
    
    for (b in batch_levels) {
      batch_indices <- which(batch == b)
      if (length(batch_indices) > 0) {
        matrix_by_batch[[as.character(b)]] <- matrix_[, batch_indices, drop = FALSE]
      }
    }
    
    # Apply NPN transformation to each batch.
    for (b in names(matrix_by_batch)) {
      matrix_t <- t(matrix_by_batch[[b]])
      npn_transformed_t <- huge::huge.npn(matrix_t, verbose = FALSE)
      matrix_by_batch[[b]] <- t(npn_transformed_t)
    }
    
    # Reassemble the matrix from the adjusted batches.
    result_matrix <- matrix_
    for (b in names(matrix_by_batch)) {
      batch_indices <- which(batch == as.character(b))
      result_matrix[, batch_indices] <- matrix_by_batch[[b]]
    }
    
    return(result_matrix)
  }
}


adjust_seurat_scaling <- function(df_, batch, data_are_counts, debug = FALSE) {
  #' Adjust using Seurat's ScaleData regression method.
  #' @param df_ The data matrix (features x samples).
  #' @param batch The batch variable vector.
  #' @param data_are_counts Logical, TRUE if data is raw counts.
  #' @return The adjusted matrix.
  
  message("Adjusting data with Seurat's ScaleData regression.")
  prep_list <- prep_seurat_like(df_, batch, data_are_counts)
  seurat_obj <- prep_list$obj
  
  all_features <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "Batch", features = all_features, verbose = FALSE)
  
  scaled_matrix <- as.matrix(GetAssayData(seurat_obj, layer = "scale.data"))
  
  return(restore_names(scaled_matrix, prep_list))
}


adjust_seurat_integration <- function(df_, batch, data_are_counts, debug = FALSE) {
  #' Adjust using Seurat's anchor-based integration workflow.
  #'
  #' This version integrates the default set of variable features and then
  #' combines the result with the original, unadjusted data for the remaining features.
  #'
  #' @param df_ The data matrix (features x samples).
  #' @param batch The batch variable vector.
  #' @param data_are_counts Logical, TRUE if data is raw counts.
  #' @return The adjusted matrix.
  
  message("Adjusting data with Seurat's RPCA integration (variable features only).")
  prep_list <- prep_seurat_like(df_, batch, data_are_counts)
  seurat_obj <- prep_list$obj
  
  seurat_obj.list <- SplitObject(seurat_obj, split.by = "Batch")
  
  # Pre-process each batch object.
  for (i in 1:length(seurat_obj.list)) {
    if (data_are_counts) {
      seurat_obj.list[[i]] <- FindVariableFeatures(seurat_obj.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    } else {
      if (debug) message("DEBUG: Finding variable features for batch ", i, " by variance ranking.")
      hvf_data <- GetAssayData(seurat_obj.list[[i]], layer = "data")
      variances <- apply(hvf_data, 1, var, na.rm = TRUE)
      top_features <- names(sort(variances, decreasing = TRUE)[1:2000])
      VariableFeatures(seurat_obj.list[[i]]) <- top_features
    }
    
    all_features_in_obj <- rownames(seurat_obj.list[[i]])
    seurat_obj.list[[i]] <- ScaleData(seurat_obj.list[[i]], features = all_features_in_obj, verbose = FALSE)
    
    num_cells <- ncol(seurat_obj.list[[i]])
    npcs_to_use <- min(50, num_cells - 1)
    seurat_obj.list[[i]] <- RunPCA(seurat_obj.list[[i]], npcs = npcs_to_use, features = VariableFeatures(seurat_obj.list[[i]]), verbose = FALSE)
  }
  
  # Find integration anchors.
  min_batch_size <- min(sapply(seurat_obj.list, ncol))
  k_anchor <- min(5, min_batch_size - 1)
  dims_to_use <- 1:min(30, min_batch_size - 1)
  
  anchors <- FindIntegrationAnchors(object.list = seurat_obj.list, reduction = "rpca", dims = dims_to_use, k.anchor = k_anchor, verbose = FALSE)
  
  if (nrow(anchors@anchors) < 2) {
    stop("Integration failed: Not enough anchors were found between datasets.")
  }
  
  # Integrate data.
  k_weight <- min(100, min_batch_size - 1) # Seurat default is 100
  if (debug) message("DEBUG: Setting k.weight to ", k_weight)
  
  seurat_obj.integrated <- IntegrateData(anchorset = anchors, k.weight = k_weight, verbose = FALSE)
  
  # Combine integrated and non-integrated features.
  integrated_matrix_sanitized <- as.matrix(GetAssayData(seurat_obj.integrated, assay = "integrated", layer = "data"))
  
  integrated_features_sanitized <- rownames(integrated_matrix_sanitized)
  all_features_sanitized <- rownames(seurat_obj)
  non_integrated_features_sanitized <- setdiff(all_features_sanitized, integrated_features_sanitized)
  
  feature_name_map <- data.frame(original = prep_list$orig_features, sanitized = rownames(seurat_obj))
  sample_name_map <- data.frame(original = prep_list$orig_samples, sanitized = colnames(seurat_obj))
  
  non_integrated_features_original <- feature_name_map$original[feature_name_map$sanitized %in% non_integrated_features_sanitized]
  
  original_non_integrated_data <- df_[non_integrated_features_original, , drop = FALSE]
  
  integrated_features_original <- feature_name_map$original[feature_name_map$sanitized %in% integrated_features_sanitized]
  rownames(integrated_matrix_sanitized) <- integrated_features_original
  colnames(integrated_matrix_sanitized) <- sample_name_map$original[match(colnames(integrated_matrix_sanitized), sample_name_map$sanitized)]
  
  combined_matrix <- rbind(integrated_matrix_sanitized, original_non_integrated_data)
  
  if (debug) message("DEBUG: Seurat integration - Combined matrix dimensions: ", nrow(combined_matrix), " rows, ", ncol(combined_matrix), " cols")
  
  # Ensure final matrix has original feature and sample order.
  final_matrix <- combined_matrix[prep_list$orig_features, prep_list$orig_samples, drop = FALSE]
  
  return(final_matrix)
}


adjust_fairadapt <- function(gene_df, batch, design, ..., debug = FALSE) {
  #' Adjust data using the fairadapt method.
  #' Note: `fairadapt` requires a design matrix with exactly one variable to preserve.
  #' @param gene_df The quantitative data (samples x features).
  #' @param batch The batch variable vector.
  #' @param design The design matrix.
  #' @return The adjusted matrix (features x samples).
  
  message("Adjusting using fairadapt.")
  
  if (ncol(design) != 2) {
    stop("fairadapt requires a design matrix with exactly one column to preserve (plus the intercept).")
  }
  
  design_col_name <- colnames(design)[colnames(design) != "(Intercept)"]
  if (length(design_col_name) != 1) {
    stop("Could not identify a unique column to preserve from the design matrix.")
  }
  
  if (debug) {
    message("DEBUG: gene_df dimensions: ", nrow(gene_df), " rows, ", ncol(gene_df), " cols")
    message("DEBUG: batch length: ", length(batch))
    message("DEBUG: design dimensions: ", nrow(design), " rows, ", ncol(design), " cols")
    message("DEBUG: design_col_name: ", design_col_name)
  }
  
  data_for_adj <- gene_df
  data_for_adj$batch <- batch
  data_for_adj[[design_col_name]] <- design[, design_col_name]
  
  adj.mat <- matrix(0, nrow = ncol(data_for_adj), ncol = ncol(data_for_adj))
  colnames(adj.mat) <- rownames(adj.mat) <- colnames(data_for_adj)
  
  batch_idx <- which(colnames(adj.mat) == "batch")
  design_idx <- which(colnames(adj.mat) == design_col_name)
  
  adj.mat[, batch_idx] <- 1
  adj.mat[design_idx, ] <- 0
  diag(adj.mat) <- 0
  
  formula <- as.formula(paste(design_col_name, "~ ."))
  
  mod <- fairadapt(formula,
    train.data = data_for_adj,
    prot.attr = "batch",
    adj.mat = adj.mat
  )
  
  adjusted_df <- mod$adapt.train[, colnames(gene_df), drop = FALSE]
  
  return(t(as.matrix(adjusted_df)))
}


adjust_liger <- function(df_, batch, data_are_counts, debug = FALSE) {
  #' Adjust using the LIGER method.
  message("Adjusting with LIGER.")
  prep_list <- prep_seurat_like(df_, batch, data_are_counts)
  
  liger_obj <- seuratToLiger(prep_list$obj)
  
  if (data_are_counts) {
    liger_obj <- normalize(liger_obj)
    liger_obj <- selectGenes(liger_obj)
  } else {
    message("Data is pre-normalized. Bypassing LIGER's normalize() and selectGenes().")
    hvf_data <- GetAssayData(prep_list$obj, layer = "data")
    variances <- apply(hvf_data, 1, var, na.rm = TRUE)
    top_features <- names(sort(variances, decreasing = TRUE)[1:2000])
    liger_obj@var.features <- top_features
  }
  
  liger_obj <- scaleNotCenter(liger_obj)
  liger_obj <- runIntegration(liger_obj, k = 20, verbose = FALSE)
  liger_obj <- quantileNorm(liger_obj, verbose = FALSE)
  
  if (debug) message("DEBUG: Available slots in liger object: ", paste(slotNames(liger_obj), collapse = ", "))
  
  corrected_matrix_sanitized <- liger_obj@W %*% t(liger_obj@H.norm)
  final_matrix_sanitized <- as.matrix(GetAssayData(prep_list$obj, layer = "data"))
  
  common_features <- intersect(rownames(corrected_matrix_sanitized), rownames(final_matrix_sanitized))
  if (debug) message("DEBUG: Liger - Found ", length(common_features), " common features between corrected and target matrices.")
  
  final_matrix_sanitized[common_features, ] <- corrected_matrix_sanitized[common_features, ]
  if (debug) message("DEBUG: Liger - Successfully overwrote variable features in the target matrix.")
  
  final_matrix <- restore_names(final_matrix_sanitized, prep_list)
  if (debug) message("DEBUG: Liger - Final matrix with restored names dimensions: ", nrow(final_matrix), " rows, ", ncol(final_matrix), " cols")
  
  return(final_matrix)
}


adjust_mnn <- function(df_, batch, data_are_counts, debug = FALSE) {
  #' Adjust using the mnnCorrect method from batchelor.
  #' @param df_ The data matrix (features x samples).
  #' @param batch The batch variable vector.
  #' @param data_are_counts Logical, TRUE if data is raw counts.
  #' @param debug Logical flag for debug output.
  #' @return The adjusted matrix.
  
  message("Adjusting with MNN.")
  prep_list <- prep_seurat_like(df_, batch, data_are_counts)
  sce_list <- lapply(unique(batch), function(b) as.SingleCellExperiment(prep_list$obj[, prep_list$obj$Batch == b]))
  sce_corrected <- do.call(batchelor::mnnCorrect, c(sce_list, list(assay.type = "logcounts")))
  corrected_matrix <- as.matrix(assay(sce_corrected, "corrected"))
  
  return(restore_names(corrected_matrix, prep_list))
}


adjust_gmm_common <- function(matrix_, batch, adjustment_strategy, strategy_name, debug = FALSE, meta_file = NULL) {
  #' Common implementation for all GMM-based adjustment methods.
  #' Wrapper for gmm_adjust to fit into the main adjustment pipeline with caching support.
  #' Handles the transposition of data to match gmm_adjust's input requirements.
  #' Parameter extraction and caching is now handled internally by bimodal_normalize and gmm_adjust.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param adjustment_strategy The strategy to pass to the underlying GMM functions.
  #' @param strategy_name Human-readable name for the strategy (for logging).
  #' @param debug Logical flag for debug output.
  #' @param meta_file Path to save the recommended modes for each gene.
  #' @return The adjusted matrix (features x samples).
  
  message("Adjusting with GMM-based adjustment (", strategy_name, ", with caching).")

  genes_df <- as.data.frame(t(matrix_))
  cache_folder <- "/data/.cache/gmm_cache"
  if (!dir.exists(cache_folder)) {
    dir.create(cache_folder, recursive = TRUE)
  }

  if (is.null(batch)) {
    # Single batch processing - bimodal_normalize will handle parameter extraction and caching
    batch_name <- "single_batch"
    
    result <- bimodal_normalize(
      genes_df, 
      gmm_parameters = NULL,  # Let bimodal_normalize extract parameters
      batch_name = batch_name, 
      cache_folder = cache_folder,
      debug = debug, 
      log_file = "/outputs/parallel/bimodal_parallel.log", 
      adjustment_strategy = adjustment_strategy, 
      num_workers = get_allocated_cores()
    )
    adjusted_genes_df <- result$bimodal_data
    
    if (!is.null(meta_file)) {
      recommended_modes <- result$recommended_modes
      write.csv(data.frame(single_batch = recommended_modes), file = meta_file, row.names = TRUE)
    }
    
  } else {
    # Multi-batch processing - gmm_adjust will coordinate per-batch parameter extraction
    adjusted_genes_df <- gmm_adjust(
      genes_df, 
      batch, 
      gmm_parameters = list(),  # Empty list - let gmm_adjust handle parameter extraction
      cache_folder = cache_folder,
      debug = debug, 
      adjustment_strategy = adjustment_strategy, 
      log_file = "/outputs/parallel/gmm_adjust_parallel.log", 
      num_workers = get_allocated_cores()
    )
    
    if (!is.null(meta_file)) {
      # For meta file, we need to extract recommended modes
      # This is a bit more complex now since we don't pre-extract parameters
      # We could either modify gmm_adjust to return this info, or extract it here
      # For now, let's create a placeholder - this could be improved
      batch_levels <- unique(batch)
      recommended_modes_df <- data.frame(
        matrix(1, nrow = length(batch_levels), ncol = ncol(genes_df)),  # Default to 1 mode
        row.names = batch_levels
      )
      colnames(recommended_modes_df) <- colnames(genes_df)
      write.csv(recommended_modes_df, file = meta_file, row.names = TRUE)
      
      if (debug) {
        message("DEBUG: Meta file created with default recommended modes. Consider enhancing gmm_adjust to return this information.")
      }
    }
  }

  return(t(as.matrix(adjusted_genes_df)))
}

adjust_gmm <- function(matrix_, batch, debug = FALSE, meta_file = NULL) {
  #' Wrapper for GMM adjustment with simple strategy.
  return(adjust_gmm_common(matrix_, batch, "simple", "Simple", debug, meta_file))
}

adjust_gmm_npn <- function(matrix_, batch, debug = FALSE, meta_file = NULL) {
  #' Wrapper for GMM adjustment with NPN strategy.
  return(adjust_gmm_common(matrix_, batch, "npn", "NPN", debug, meta_file))
}

adjust_gmm_npn_unit_std <- function(matrix_, batch, debug = FALSE, meta_file = NULL) {
  #' Wrapper for GMM adjustment with NPN unit standard deviation strategy.
  return(adjust_gmm_common(matrix_, batch, "npn_unit_std", "NPN, Unit Std", debug, meta_file))
}

adjust_gmm_global_simple <- function(matrix_, batch = NULL, debug = FALSE, meta_file = NULL) {
  #' Simple gene-global GMM adjustment strategy with scaling.
  #' Computes global distribution by averaging across genes, fits 2-component GMM,
  #' centers data at midpoint between modes, and scales so modes end up at -1 and +1.
  #' Works entirely in log space for mathematical consistency.
  
  message("Adjusting with simple gene-global GMM (log space).")
  data_transposed <- t(matrix_)
  
  # Apply the simple global adjustment
  result <- gmm_global_simple_core(data_transposed, debug = debug)
  
  return(t(result$adjusted_data))
}

adjust_gmm_global_npn <- function(matrix_, batch = NULL, debug = FALSE, meta_file = NULL) {
  #' Gene-global GMM adjustment with bimodal NPN transformation.
  #' Computes global distribution by averaging across genes, fits 2-component GMM,
  #' then applies bimodal NPN transformation to all expression values using global parameters.
  
  message("Adjusting with gene-global GMM bimodal NPN.")
  data_transposed <- t(matrix_)
  
  # Apply the global NPN adjustment
  result <- gmm_global_npn_core(data_transposed, debug = debug)
  
  # Ensure we have a matrix before transposing
  adjusted_data <- result$adjusted_data
  if (!is.matrix(adjusted_data)) {
    if (debug) {
      message("Converting to matrix before transpose")
    }
    adjusted_data <- as.matrix(adjusted_data)
  }
  
  return(t(adjusted_data))
}

adjust_gmm_streamlined <- function(matrix_, batch, debug = FALSE, meta_file = NULL, mean_only = FALSE) {
  #' Streamlined GMM adjustment using the fast implementation.
  #' Applies bimodal GMM transformation to all genes.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param debug Logical flag for debug output.
  #' @param meta_file Path to save metadata (currently unused in streamlined version).
  #' @param mean_only If TRUE, only adjust means without using inverse CDF transformation.
  #' @return The adjusted matrix (features x samples).
  
  if (mean_only) {
    message("Adjusting with streamlined GMM (mean-only adjustment, no inverse CDF).")
  } else {
    message("Adjusting with streamlined GMM (bimodal for all genes with inverse CDF).")
  }
  
  # Convert to the format expected by gmm_adjust_streamlined (samples x genes)
  genes_df <- as.data.frame(t(matrix_))
  
  if (is.null(batch)) {
    # Single batch - create dummy batch
    batch <- rep("batch1", nrow(genes_df))
  }
  
  # Apply streamlined adjustment with mean_only option
  adjusted_genes_df <- gmm_adjust_streamlined(genes_df, batch, alpha0 = 10, mean_only = mean_only, debug = debug)
  
  # Convert back to matrix format (features x samples)
  return(t(as.matrix(adjusted_genes_df)))
}

adjust_gmm_streamlined_mean_only <- function(matrix_, batch, debug = FALSE, meta_file = NULL) {
  #' Streamlined GMM adjustment with mean-only transformation.
  #' Applies bimodal GMM transformation to all genes but only adjusts means without inverse CDF.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param debug Logical flag for debug output.
  #' @param meta_file Path to save metadata (currently unused in streamlined version).
  #' @return The adjusted matrix (features x samples).
  
  return(adjust_gmm_streamlined(matrix_, batch, debug = debug, meta_file = meta_file, mean_only = TRUE))
}


rank_normalized <- function(matrix_, dim) {
  if (dim < 1 || dim > 2) {
    stop("Invalid dimension. Must be 1 for rows or 2 for columns.")
  }
  ranked = apply(matrix_, dim, rank, ties.method = "average")
  
  # apply() transposes the result when dim=1, so we need to transpose it back
  # When dim=1: apply ranks across columns (samples) for each row (feature)
  # When dim=2: apply ranks across rows (features) for each column (sample)
  if (dim == 1 && is.matrix(ranked)) {
    ranked = t(ranked)
  }
  
  return(ranked / max(ranked, na.rm = TRUE))
}

adjust_ranked <- function(matrix_, debug = FALSE) {
  #' Normalize sample-wise by ranking the genes within the sample.
  message("Adjusting with ranked.")
  # Data has genes as rows, so we need to rank along the rows.
  return(rank_normalized(matrix_, 1))
}

adjust_ranked_twice <- function(matrix_, debug = FALSE) {
  message("Adjusting with ranked twice.")
  #' Normalize sample-wise by ranking the genes within the sample.
  # Next, rank by sample. This tells us something about whether the gene is up or down regulated in the sample.
  return(rank_normalized(rank_normalized(matrix_, 1), 2))
}

adjust_ranked_projection <- function(matrix_, debug = FALSE) {
  message("Adjusting with ranked projection.")
  # This is equivalent to projecting the ranks onto a n-1 dimensional ball.
  ranked <- apply(matrix_, 1, rank, ties.method = "average")
  centered <- ranked - rowMeans(ranked)
  return(centered / sqrt(rowSums(centered^2)))
}

adjust_ranked_with_batch_info <- function(matrix_, batch, debug = FALSE) {
  #' Normalize sample-wise by ranking the genes within the sample, and then by batch.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param debug Logical flag for debug output.
  #' @return The adjusted matrix (features x samples).
  
  message("Adjusting with ranked with batch info.")
  ranked = rank_normalized(matrix_, 1)
  
  if (debug) {
    message("DEBUG: matrix_ dimensions: ", nrow(matrix_), " x ", ncol(matrix_))
    message("DEBUG: ranked dimensions: ", nrow(ranked), " x ", ncol(ranked))
  }
  
  batch_levels <- unique(batch)
  ranked2 <- matrix(NA, nrow = nrow(ranked), ncol = ncol(ranked))
  
  for (b in batch_levels) {
    # For each batch, we rank by sample.
    batch_indices <- which(batch == b)
    batch_data <- ranked[, batch_indices, drop = FALSE]
    
    if (debug) {
      message("DEBUG: Processing batch '", b, "' with ", length(batch_indices), " samples")
      message("DEBUG: batch_data dimensions: ", nrow(batch_data), " x ", ncol(batch_data))
    }
    
    # Only apply ranking if there's more than one sample in the batch
    if (ncol(batch_data) > 1) {
      batch_ranked <- rank_normalized(batch_data, 2)
      if (debug) {
        message("DEBUG: batch_ranked dimensions: ", nrow(batch_ranked), " x ", ncol(batch_ranked))
      }
      ranked2[, batch_indices] <- batch_ranked
    } else {
      # For single-sample batches, just use the original ranked values
      ranked2[, batch_indices] <- batch_data
    }
  }
  
  # Handle any remaining NA values
  if (any(is.na(ranked2))) {
    message("WARNING: Found NA values in ranked2 matrix. Replacing with original ranked values.")
    ranked2[is.na(ranked2)] <- ranked[is.na(ranked2)]
  }
  
  max_val <- max(ranked2, na.rm = TRUE)
  if (max_val == 0) {
    message("WARNING: Maximum value in ranked2 is 0. Using 1 as denominator.")
    max_val <- 1
  }
  
  return(ranked2 / max_val)
}


adjust_log <- function(matrix_, batch, debug = FALSE) {
  message("Adjusting with log.")
  # For each batch, subtract the minimum from every entry. Then take the log of the data.
  batch_levels <- unique(batch)
  adjusted <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  for (b in batch_levels) {
    batch_indices <- which(batch == b)
    batch_data <- matrix_[, batch_indices, drop = FALSE]
    min_val <- min(batch_data, na.rm = TRUE)
    adjusted[, batch_indices] <- log(batch_data - min_val + 1)
  }
  return(adjusted)
}

adjust_log_combat <- function(matrix_, batch, design, debug = FALSE) {
  message("Adjusting with log combat.")
  # For each batch, subtract the minimum from every entry. Then take the log of the data.
  # On this data, run ComBat
  batch_levels <- unique(batch)
  adjusted <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  for (b in batch_levels) {
    batch_indices <- which(batch == b)
    batch_data <- matrix_[, batch_indices, drop = FALSE]
    min_val <- min(batch_data, na.rm = TRUE)
    adjusted[, batch_indices] <- log(batch_data - min_val + 1)
  }
  
  return(ComBat_ignore_nonvariance(adjusted, batch, design))
}


# Main Orchestration Function -------------------------------------------------

batch_adjust_tidy <- function(df, input_file, adjuster, batch_col, column, full_design_matrix, reduce_rows=NULL, reduce_cols=NULL, debug=FALSE, meta_file=NULL) {
  #' Main function to orchestrate the batch adjustment process.
  #' @param df The input tidy data frame (samples x columns).
  #' @param adjuster The name of the adjustment method to use.
  #' @param batch_col The name of the batch column.
  #' @return A tidy data frame with adjusted values.
  
  # --- 0. Pre-processing: Handle NA values ---
  message("0. Pre-processing: Checking for NA values in the batch column.")
  if (!is.null(batch_col) && any(is.na(df[[batch_col]]))) {
    na_count <- sum(is.na(df[[batch_col]]))
    message("WARNING: Found ", na_count, " NA values in batch column ('", batch_col, "'). These samples (rows) will be removed before adjustment.")
    
    df <- df[!is.na(df[[batch_col]]), ]
    
    if (nrow(df) == 0) {
      stop("All samples were removed due to NA values in the batch column. Cannot proceed.")
    }
  }

  # --- 0.5. Reduce rows, if applicable ---
  if (!is.null(reduce_rows)) {
    message("Reducing rows to ", reduce_rows, " rows.")
    df <- df[1:reduce_rows, ]
  }
  
  # --- 1. Separate data components ---
  message("1. Separating metadata, batch, and gene data.")
  original_colnames <- colnames(df)

  # Handle first column containing sample IDs
  first_col_name <- colnames(df)[1]
  if (first_col_name %in% c("...1", "", "X", "meta_Sample_ID")) {
    if (first_col_name == "meta_Sample_ID") {
      message("Detected meta_Sample_ID column - setting as row names")
    } else {
      message("Detected unnamed first column with sample IDs - setting as row names")
    }
    # Convert to data.frame to avoid tibble deprecation warning
    df <- as.data.frame(df)
    rownames(df) <- df[[1]]
    df <- df[, -1]
  }

  if (!is.null(batch_col)) {
    batch <- df[[batch_col]]
    df[[batch_col]] <- NULL
  }
  else {
    batch <- NULL
  }
  
  # Metadata columns start with "meta_".

  meta_data_names <- colnames(df)[startsWith(colnames(df), "meta_")]
  metadata_cols <- df[, startsWith(colnames(df), "meta_"), drop = FALSE]
  message("Metadata cols: ", paste(colnames(metadata_cols), collapse = ", "))
  
  gene_col_names = colnames(df)[!startsWith(colnames(df), "meta_")]
  genes <- df[, !startsWith(colnames(df), "meta_")]
  
  # Verify all gene columns are numeric
  numeric_cols <- sapply(genes, is.numeric)
  if (!all(numeric_cols)) {
    non_numeric <- colnames(genes)[!numeric_cols]
    stop(sprintf("Non-numeric columns found in gene data: %s", paste(non_numeric, collapse = ", ")))
  }
  
  # --- 2. Create design matrix ---
  message("2. Creating design matrix.")
  design <- create_design_matrix(metadata_cols, column, full_design_matrix)
  
  # --- 3. Transpose gene data (with caching) ---
  transposed_cache_file <- sub("(\\.[^.]+)$", "_t\\1", input_file)
  
  if (file.exists(transposed_cache_file)) {
    message("Loading cached transposed data from '", transposed_cache_file, "'")
    mat_genes <- as.matrix(read.csv(transposed_cache_file, row.names = 1, check.names = FALSE))
    
    # Verify cached data matches current gene columns
    if (ncol(mat_genes) != nrow(genes) || !all(colnames(mat_genes) == rownames(genes))) {
      message("WARNING: Cached transposed data doesn't match current gene data. Regenerating...")
      mat_genes <- transpose_essential(genes)
      message("Updating cached transposed data to '", transposed_cache_file, "'")
      write.csv(mat_genes, transposed_cache_file, row.names = TRUE, quote = FALSE)
    }
  } else {
    message("3. Transposing gene data for adjustment (features x samples).")
    mat_genes <- transpose_essential(genes)
    message("Caching transposed data to '", transposed_cache_file, "'")
    write.csv(mat_genes, transposed_cache_file, row.names = TRUE, quote = FALSE)
  }

  if (!is.null(reduce_cols)) {
    message("Reducing gene expression columns to ", reduce_cols, " columns.")
    mat_genes <- mat_genes[1:reduce_cols, ]
  }
  if (!is.null(reduce_rows)) {
    message("Reducing gene expression rows to ", reduce_rows, " rows.")
    mat_genes <- mat_genes[, 1:reduce_rows]
  }

  # Determine if data are counts (all non-negative).
  data_are_counts <- !any(mat_genes < 0, na.rm = TRUE)
  
  # --- 4. Apply adjustment method ---
  message("4. Applying '", adjuster, "' adjustment method.")
  
  adjusted_matrix <- switch(adjuster,
    "gmm" = adjust_gmm(mat_genes, batch, debug=debug, meta_file=meta_file),
    "gmm_npn" = adjust_gmm_npn(mat_genes, batch, debug=debug, meta_file=meta_file),
    "gmm_npn_unit_std" = adjust_gmm_npn_unit_std(mat_genes, batch, debug=debug, meta_file=meta_file),
    "gmm_streamlined" = adjust_gmm_streamlined(mat_genes, batch, debug=debug, meta_file=meta_file),
    "gmm_streamlined_mean_only" = adjust_gmm_streamlined_mean_only(mat_genes, batch, debug=debug, meta_file=meta_file),
    "gmm_global_simple" = adjust_gmm_global_simple(mat_genes, batch, debug=TRUE, meta_file=meta_file),
    "gmm_global_npn" = adjust_gmm_global_npn(mat_genes, batch, debug=TRUE, meta_file=meta_file),
    "ranked1" = adjust_ranked(mat_genes, debug = debug),
    "ranked2" = adjust_ranked_twice(mat_genes, debug = debug),
    "ranked_batch" = adjust_ranked_with_batch_info(mat_genes, batch, debug = debug),
    "ranked_projection" = adjust_ranked_projection(mat_genes, debug = debug),
    "min_mean" = adjust_min_mean(mat_genes, batch, debug = debug),
    "combat" = adjust_combat(mat_genes, batch, design, data_are_counts, debug = debug),
    "log" = adjust_log(mat_genes, batch, debug = debug),
    "log_combat" = adjust_log_combat(mat_genes, batch, design, debug = debug),
    "limma" = adjust_limma(mat_genes, batch, design, debug = debug),
    "quantile" = adjust_quantile(mat_genes, batch, debug = debug),
    "npn" = adjust_npn(mat_genes, batch, debug = debug),
    "seurat_scaling" = adjust_seurat_scaling(mat_genes, batch, data_are_counts, debug = debug),
    "seurat_integration" = adjust_seurat_integration(mat_genes, batch, data_are_counts, debug = debug),
    "fairadapt" = adjust_fairadapt(genes, batch, design, debug = debug),
    "mnn" = adjust_mnn(mat_genes, batch, data_are_counts, debug = debug),
    "liger" = adjust_liger(mat_genes, batch, data_are_counts, debug = debug),
    stop(sprintf("Unknown adjuster '%s'", adjuster))
  )
  
  if (debug) message("DEBUG: Dimensions of final adjusted matrix before transposing: ", nrow(adjusted_matrix), " rows, ", ncol(adjusted_matrix), " cols")
  
  # --- 5. Reconstruct tidy data frame ---
  message("5. Reconstructing the tidy data frame.")
  adjusted_df <- as.data.frame(t(adjusted_matrix))

  # Ensure column names match - handle dimension mismatches
  if (ncol(adjusted_df) != length(gene_col_names)) {
    message("WARNING: Dimension mismatch between adjusted data (", ncol(adjusted_df), " cols) and original gene columns (", length(gene_col_names), " cols)")
    if (ncol(adjusted_df) < length(gene_col_names)) {
      # Use only the available columns
      gene_col_names <- gene_col_names[1:ncol(adjusted_df)]
      message("Using first ", ncol(adjusted_df), " gene column names")
    } else {
      # Pad with generic names if needed
      extra_names <- paste0("Gene_", (length(gene_col_names) + 1):ncol(adjusted_df))
      gene_col_names <- c(gene_col_names, extra_names)
      message("Added ", length(extra_names), " generic column names")
    }
  }
  
  colnames(adjusted_df) <- gene_col_names

  if (debug) message("DEBUG: Dimensions of final adjusted matrix after transposing: ", nrow(adjusted_df), " rows, ", ncol(adjusted_df), " cols")

  # Restore other columns if reduced
  if (!is.null(reduce_cols)) {
    message(" 5.1 Restoring skipped columns (reduced)")
    adjusted_df <- cbind(adjusted_df, genes[ ,reduce_cols:ncol(genes)])
  }
  
  if (is.null(batch)){
    message(" 5.2 Restoring metadata columns")
    final_df <- cbind(metadata_cols, adjusted_df)
  }
  else {
    message(" 5.2 Restoring batch and metadata columns")
    final_df <- cbind(batch, metadata_cols, adjusted_df)
    colnames(final_df)[1] <- batch_col
  }
  
  if (sum(startsWith(colnames(final_df), "meta_")) == 0 && length(meta_data_names) > 0) {
    stop("No metadata columns found in final data frame.")
  }
  
  return(final_df)
}


# Set future plan for parallel processing -------------------------------------
# Increase the maximum size for global variables.
options(future.globals.maxSize = 2000 * 1024^2)


# Parse command line arguments ------------------------------------------------

parser <- ArgumentParser(description = "A script to apply various batch correction methods to tidy data.")

parser$add_argument("input_file", help = "Path to the input CSV file. Rows are samples, columns are features/metadata.")
parser$add_argument("output_file", help = "Path for the output adjusted CSV file.")
parser$add_argument("-a", "--adjuster",
  default = "combat",
  help = "Batch adjustment method to use."
)
parser$add_argument("-b", "--batch-col", default = NULL, help = "Name of the column identifying the batch for each sample.")
parser$add_argument("-c", "--debug", action = "store_true", help = "Enable verbose debugging messages.")
parser$add_argument("-R", "--reduce-rows", default = NULL, help = "Number of rows to reduce to.")
parser$add_argument("-r", "--reduce-cols", default = NULL, help = "Number of columns to reduce to.")
parser$add_argument("-m", "--meta-file", default = NULL, help = "Path to save the recommended modes for each gene.")
parser$add_argument("--mean-only", action = "store_true", help = "For GMM streamlined methods, only adjust means without using inverse CDF transformation.")
parser$add_argument("--column", default = NULL, help = "Specific metadata columns to include in the design matrix.")
parser$add_argument("--full-design-matrix", action = "store_true", help = "Use all metadata columns in the design matrix.")


args <- parser$parse_args()

# Main Execution --------------------------------------------------------------

message("Reading input file '", args$input_file, "'")
suppressMessages(df <- vroom(args$input_file, show_col_types = FALSE))
message("Input file has ", nrow(df), " rows and ", ncol(df), " columns.")

if (!is.null(args$batch_col) && !(args$batch_col %in% names(df))) {
  available_meta_cols <- names(df)[startsWith(names(df), "meta_")]
  stop(sprintf(
    "The specified batch column ('%s') was not found in the input file. Please check the column name. Available metadata columns are: %s",
    args$batch_col, paste(available_meta_cols, collapse = ", ")
  ))
}

# Changing "Sample_ID" column to "meta_Sample_ID" to simplify later code.
if ("Sample_ID" %in% names(df)) {
  df <- df %>% dplyr::rename(meta_Sample_ID = Sample_ID)
}



message("Starting batch adjustment with method: '", args$adjuster, "'")

initial_start_time = Sys.time()
adjusted_data <- batch_adjust_tidy(
  df,
  input_file = args$input_file,
  adjuster = args$adjuster,
  batch_col = args$batch_col,
  column = args$column,
  full_design_matrix = args$`full-design-matrix`,
  debug = args$debug,
  meta_file = args$meta_file
)
message("Adjusting ", args$input_file, " took ", Sys.time() - initial_start_time, " seconds.")

start_time = Sys.time()
message("Writing adjusted data to '", args$output_file, "'")
# Round numeric columns and write to output file.
adjusted_data %>%
  mutate(across(where(is.numeric), ~ sprintf("%.6f", .))) %>%
  fwrite(args$output_file)

elapsed_time = Sys.time() - start_time
message("Successfully saved adjusted data to '", args$output_file, "'", " in ", elapsed_time, " seconds.")
