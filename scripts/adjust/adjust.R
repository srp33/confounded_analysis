# Load dependencies --------------------------------
# Suppress package startup messages for a cleaner console output
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
})

# Increase the maximum size for global variables in parallel processing
# This is necessary for Seurat integration with large datasets.
options(future.globals.maxSize = 2000 * 1024^2)

# Print Seurat version for debugging
message(paste("Using Seurat version:", packageVersion("Seurat")))

# Parse command line arguments --------------------------

parser <- ArgumentParser(description = "A script to apply various batch correction methods to tidy data.")

parser$add_argument("input_file", help = "Path to the input CSV file. Rows are samples, columns are features/metadata.")
parser$add_argument("output_file", help = "Path for the output adjusted CSV file.")
parser$add_argument("-a", "--adjuster", default = "combat",
                    choices = c("combat", "limma", "seurat_scaling", "seurat_integration", "harmony", "quantile", "fairadapt", "scanorama", "liger", "fastMNN"),
                    help = "Batch adjustment method to use.")
parser$add_argument("-b", "--batch-col", default = "Batch", help = "Name of the column identifying the batch for each sample.")
parser$add_argument("-c", "--column", nargs = "+", default = NULL, help = "Predictive columns to preserve. If specified, these are used to build the design matrix.")
parser$add_argument("-f", "--full-design-matrix", action = "store_true", help = "If set, the design matrix will include all categorical metadata variables.")

args <- parser$parse_args()

# Helper Functions ---------------------------------

transpose_essential <- function(quantitative_df) {
  #' Transposes a quantitative data frame robustly.
  #'
  #' This multi-step process ensures that sample names (row names) are
  #' correctly preserved during the conversion from a data frame to a
  #' transposed matrix (features x samples).
  #'
  #' @param quantitative_df A data frame with samples as rows and features as columns.
  #' @return A transposed matrix with features as rows and samples as columns.
  
  df_version <- as.data.frame(quantitative_df)
  mat_untransposed <- as.matrix(df_version)
  rownames(mat_untransposed) <- rownames(df_version)
  mat_quantitative <- t(mat_untransposed)
  
  return(mat_quantitative)
}

ComBat_ignore_nonvariance <- function(matrix_, batch, design) {
  #' A wrapper for ComBat that handles features with zero variance.
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
  
  if(sum(varying_row_mask) < nrow(matrix_)) {
    message(
      sprintf("Found %d features with zero variance. These will be ignored by ComBat.",
              nrow(matrix_) - sum(varying_row_mask))
    )
  }
  
  matrix_[varying_row_mask, ] <- ComBat(matrix_[varying_row_mask, ], batch, mod = design, prior.plots = FALSE)
  return(matrix_)
}

create_design_matrix <- function(categorical_df, columns_to_use = NULL, use_all = FALSE) {
  #' Creates a design matrix from a data frame of categorical variables.
  #'
  #' @param categorical_df A data frame with samples as rows and categorical metadata as columns.
  #' @param columns_to_use A character vector of specific columns to include in the model.
  #' @param use_all A boolean indicating whether to use all columns in the data frame.
  #' @return A design matrix.

  if (!is.null(columns_to_use)) {
    if (!all(columns_to_use %in% colnames(categorical_df))) {
      stop("One or more specified columns for the design matrix were not found in the metadata.")
    }
    design_df <- categorical_df[, columns_to_use, drop = FALSE]
    message(sprintf("Creating design matrix from specified columns: %s", paste(columns_to_use, collapse = ", ")))
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
  #' Prepares a Seurat object for downstream adjustments.
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
    #' Restores original feature and sample names to a matrix.
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

# Adjustment Functions ---------------------------------

adjust_combat <- function(matrix_, batch, design, data_are_counts) {
  #' Adjusts a matrix using ComBat or ComBat_seq.
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

adjust_limma <- function(matrix_, batch, design, ...) {
  #' Adjusts a matrix using limma's removeBatchEffect.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @param batch The batch variable vector.
  #' @param design The design matrix.
  #' @return The adjusted matrix.

  message("Adjusting data with limma::removeBatchEffect.")
  return(limma::removeBatchEffect(matrix_, batch = batch, design = design))
}

adjust_quantile <- function(matrix_, ...) {
  #' Adjusts a matrix using quantile normalization.
  #' @param matrix_ The matrix to adjust (features x samples).
  #' @return The adjusted matrix.

  message("Adjusting using quantile normalization.")
  return(t(normalizeQuantiles(as.matrix(t(matrix_)))))
}

adjust_seurat_scaling <- function(df_, batch, data_are_counts, ...) {
  #' Adjusts using Seurat's ScaleData regression method.
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

adjust_seurat_integration <- function(df_, batch, data_are_counts, ...) {
  #' Adjusts using Seurat's anchor-based integration workflow.
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
  
  # The RPCA workflow requires that PCA be run on each individual object in the list.
  for (i in 1:length(seurat_obj.list)) {
      if (data_are_counts) {
          seurat_obj.list[[i]] <- FindVariableFeatures(seurat_obj.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      } else {
          # For pre-normalized data, statistical models like 'vst' or 'disp' can fail.
          # A more robust method is to rank by variance directly.
          message(sprintf("Finding variable features for batch %d by variance ranking.", i))
          hvf_data <- GetAssayData(seurat_obj.list[[i]], layer="data")
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
  
  min_batch_size <- min(sapply(seurat_obj.list, ncol))
  k_anchor <- min(5, min_batch_size - 1)
  dims_to_use <- 1:min(30, min_batch_size - 1)
  
  anchors <- FindIntegrationAnchors(object.list = seurat_obj.list, reduction = "rpca", dims = dims_to_use, k.anchor = k_anchor, verbose = FALSE)
  
  if (nrow(anchors@anchors) < 2) {
    stop("Integration failed: Not enough anchors were found between datasets.")
  }
  
  k_weight <- min(67, min_batch_size - 1)
  message(sprintf("DEBUG: Setting k.weight to %d", k_weight))

  seurat_obj.integrated <- IntegrateData(anchorset = anchors, k.weight = k_weight, verbose = FALSE)
  
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
  
  message(sprintf("DEBUG: Seurat integration - Combined matrix dimensions: %d rows, %d cols", nrow(combined_matrix), ncol(combined_matrix)))
  
  final_matrix <- combined_matrix[prep_list$orig_features, prep_list$orig_samples, drop = FALSE]
  
  return(final_matrix)
}

adjust_fairadapt <- function(quantitative_df, batch, design, ...) {
    #' Adjusts data using the fairadapt method.
    #' Note: `fairadapt` requires a design matrix with exactly one variable to preserve.
    #' @param quantitative_df The quantitative data (samples x features).
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

    data_for_adj <- quantitative_df
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
                     adj.mat = adj.mat)

    adjusted_df <- mod$adapt.train[, colnames(quantitative_df), drop = FALSE]

    return(t(as.matrix(adjusted_df)))
}

adjust_liger <- function(df_, batch, data_are_counts) {
  #' Adjusts using the LIGER method.
  message("Adjusting with LIGER.")
  prep_list <- prep_seurat_like(df_, batch, data_are_counts)
  
  liger_obj <- seuratToLiger(prep_list$obj)

  liger_obj <- normalize(liger_obj)
  liger_obj <- selectGenes(liger_obj)
  liger_obj <- scaleNotCenter(liger_obj)

  # Use the modern rliger functions
  liger_obj <- runIntegration(liger_obj, k = 20, verbose = FALSE) 
  liger_obj <- quantileNorm(liger_obj, verbose = FALSE)
  
  # DEBUG: Print the available slots in the liger object to find the correct one.
  message("DEBUG: Available slots in liger object: ", paste(slotNames(liger_obj), collapse = ", "))
  
  # The corrected data is in the 'H.norm' slot after quantileNorm.
  # This is a matrix of cell factor loadings (cells x k).
  # We reconstruct the expression matrix by multiplying with W (genes x k).
  corrected_matrix_sanitized <- liger_obj@W %*% t(liger_obj@H.norm)
  
  message(sprintf("DEBUG: Liger - Reconstructed matrix dimensions: %d rows, %d cols", nrow(corrected_matrix_sanitized), ncol(corrected_matrix_sanitized)))
  
  # Create a final matrix with the original data. We will overwrite the variable features.
  final_matrix <- df_
  
  # Restore original names to the corrected matrix so we can match them
  corrected_matrix_original_names <- restore_names(corrected_matrix_sanitized, prep_list)

  message(sprintf("DEBUG: Liger - Corrected matrix with original names dimensions: %d rows, %d cols", nrow(corrected_matrix_original_names), ncol(corrected_matrix_original_names)))
  
  # Overwrite the original data with the corrected data for the variable features
  # This is robust to cases where LIGER may have dropped some features/cells
  final_matrix[rownames(corrected_matrix_original_names), colnames(corrected_matrix_original_names)] <- corrected_matrix_original_names
  
  message(sprintf("DEBUG: Liger - Final combined matrix dimensions: %d rows, %d cols", nrow(final_matrix), ncol(final_matrix)))

  return(final_matrix)
}


# Main Orchestration Function --------------------------------

batch_adjust_tidy <- function(df, adjuster, batch_col) {
  #' Main function to orchestrate the batch adjustment process.
  #' @param df The input tidy data frame (samples x columns).
  #' @param adjuster The name of the adjustment method to use.
  #' @param batch_col The name of the batch column.
  #' @return A tidy data frame with adjusted values.

  message("1. Separating metadata, batch, and quantitative data.")
  original_colnames <- colnames(df)
  batch <- df[[batch_col]]
  
  df[[batch_col]] <- NULL

  is_numeric_col <- sapply(df, is.numeric)
  quantitative <- df[, is_numeric_col, drop = FALSE]
  metadata_cols <- df[, !is_numeric_col, drop = FALSE]
  
  message("2. Creating design matrix.")
  design <- create_design_matrix(metadata_cols, args$column, args$full_design_matrix)

  message("3. Transposing quantitative data for adjustment (features x samples).")
  mat_quantitative <- transpose_essential(quantitative)
  
  data_are_counts <- !any(mat_quantitative < 0, na.rm = TRUE)
  
  message(sprintf("4. Applying '%s' adjustment method.", adjuster))
  
  adjusted_matrix <- switch(
    adjuster,
    combat = adjust_combat(mat_quantitative, batch, design, data_are_counts),
    limma = adjust_limma(mat_quantitative, batch, design),
    quantile = adjust_quantile(mat_quantitative),
    seurat_scaling = adjust_seurat_scaling(mat_quantitative, batch, data_are_counts),
    seurat_integration = adjust_seurat_integration(mat_quantitative, batch, data_are_counts),
    fairadapt = adjust_fairadapt(quantitative, batch, design),
    fastMNN = {
        message("Adjusting with fastMNN.")
        prep_list <- prep_seurat_like(mat_quantitative, batch, data_are_counts)
        sce_list <- lapply(unique(batch), function(b) as.SingleCellExperiment(prep_list$obj[, prep_list$obj$Batch == b]))
        # The 'scale.data is empty' warning is expected as fastMNN works on the logcounts layer.
        sce_corrected <- do.call(fastMNN, c(sce_list, list(assay.type = "logcounts")))
        corrected_matrix <- as.matrix(assay(sce_corrected, "reconstructed"))
        restore_names(corrected_matrix, prep_list)
    },
    scanorama = {
        message("Adjusting with Scanorama.")
        prep_list <- prep_seurat_like(mat_quantitative, batch, data_are_counts)
        sce_list <- lapply(unique(batch), function(b) as.SingleCellExperiment(prep_list$obj[, prep_list$obj$Batch == b]))
        sce_corrected <- do.call(runScanorama, c(sce_list, list(assay.type = "logcounts")))
        corrected_matrix <- as.matrix(assay(sce_corrected, "reconstructed"))
        restore_names(corrected_matrix, prep_list)
    },
    liger = adjust_liger(mat_quantitative, batch, data_are_counts),
    stop(sprintf("Unknown adjuster '%s'", adjuster))
  )
  
  message(sprintf("DEBUG: Dimensions of final adjusted matrix before transposing: %d rows, %d cols", nrow(adjusted_matrix), ncol(adjusted_matrix)))
  
  message("5. Reconstructing the tidy data frame.")
  adjusted_df <- as.data.frame(t(adjusted_matrix))
  
  final_df <- cbind(batch, metadata_cols, adjusted_df)
  colnames(final_df)[1] <- batch_col
  
  return(final_df[, original_colnames])
}


# Main Execution --------------------------------

message(sprintf("Reading input file '%s'", args$input_file))
suppressMessages(df <- vroom(args$input_file, show_col_types = FALSE))
message(sprintf("Input file has %d rows and %d columns.", nrow(df), ncol(df)))

if (!(args$batch_col %in% names(df))) {
  stop(sprintf(
    "The specified batch column ('%s') was not found in the input file. Please check the column name.",
    args$batch_col
  ))
}

message(sprintf("Starting batch adjustment with method: '%s'", args$adjuster))

adjusted_data <- batch_adjust_tidy(
  df,
  adjuster = args$adjuster,
  batch_col = args$batch_col
)

adjusted_data %>%
  mutate(across(where(is.numeric), ~round(., 6))) %>%
  write_csv(args$output_file)

message(sprintf("Successfully saved adjusted data to '%s'", args$output_file))
