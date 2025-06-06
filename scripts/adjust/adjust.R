#Debugging: 
options(error=function() { traceback(4); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

# Load tampor package
source("/opt/TAMPOR/TAMPOR.R")

# Load dependencies --------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(docstring)
  library(stringr)
  library(argparse)
  library(sva)
  library(doParallel)
  library(ggplot2)
  library(ggpubr)
  library(vsn)
  library(limma)
  library(vroom)
  library(Seurat)
  library(patchwork)
})

# Parse command line args --------------------------

parser <- ArgumentParser()

parser$add_argument("input_file", help = "path to the input file")
parser$add_argument("output_file", help = "path to the output file.")
parser$add_argument("-a", "--adjuster", default = "combat", choices = c("combat", "min_mean", "tampor", "limma"), help = "method to use for adjustment")
parser$add_argument("-b", "--batch-col", default = "Batch", help = "title of batch column to adjust for")
parser$add_argument("-c", "--column", nargs="+", default = NULL, help = "Predictive columns. If specified, batch correction will attempt to preserve predictive power for these columns.")
parser$add_argument("-f", "--full-design-matrix", action = "store_true", help = "If set, the design matrix will include all categorical variables")
args <- parser$parse_args()

# Define functions ---------------------------------
message_structure <- function(df, pretext = "Structure of data frame") {
  message(
    paste0(
      pretext, ": ",
      paste(utils::capture.output(str(df)), collapse = "\n")
    )
  )
}

adjust_combat <- function(matrix_, batch, design, data_are_counts = FALSE) {
  #' Adjusts a matrix using ComBat.
  #'
  #' @param matrix_ The matrix to batch adjust with ComBat. Columns are samples,
  #' rows are features.
  #' @param batch The per-sample batch assignments. 
  #' @param design The design matrix to use for adjustment.
  #' @param data_are_counts If TRUE, use ComBat_seq.
  #'
  #' @return The matrix_ after batch adjustment.
  
  if (data_are_counts) {
    message("Using ComBat_seq for count data.")
    return(ComBat_seq(matrix_, batch, covar_mod=design))
  } else {
    message("Using ComBat for continuous data.")
    return(ComBat_ignore_nonvariance(matrix_, batch, design))
  }
}

ComBat_ignore_nonvariance <- function(matrix_, batch, design) {
  #' Run ComBat and ignore nonvarying features.
  #'
  #' ComBat requires that all features have some variance (and probably assumes
  #' that all features are normally distributed). Since some features don't
  #' vary across samples, this function ignores nonvarying features before
  #' running ComBat.
  #'
  varying_row_mask <- apply(matrix_, 1, function(x) { length(unique(x)) > 1 })
  message("Varying rows: ", sum(varying_row_mask), " out of ", nrow(matrix_))
  matrix_[varying_row_mask,] <- ComBat(matrix_[varying_row_mask,], batch, mod=design, prior.plots=FALSE)
  matrix_
}

match_two_stats <- function(matrix_, batch, stat1, stat2) {
  #' Matches batches by scaling so that two statistics are equal.
  #'
  #' @param matrix_ The matrix to batch adjust by scaling. Columns are
  #' samples, rows are features.
  #' @param batch The per-sample batch assignments.
  #' @param stat1 The first statistic to match. (as a function, not a string)
  #' @param stat2 The second statistic to match.
  #'
  #' @return The matrix_ after batch adjustment.
  row_names = rownames(matrix_)
  # Get rowwise stats
  overall_stat1 <- apply(matrix_, 1, stat1)
  overall_stat2 <- apply(matrix_, 1, stat2)

  # Scale each batch individually so that stat1 is 0 and stat2 is 1
  for (b in levels(factor(batch))) {
    # Get columns for this batch
    # drop=F makes it return a matrix when you only grab one column
    batch_cols <- matrix_[, batch == b, drop = FALSE]

    if (ncol(batch_cols) == 0) {
      stop(sprintf("Can't scale rows: batch '%s' has <= 1 sample.", b))
    }

    adjusted = apply(batch_cols, 1, function(x) {
      if (all(x == 0))
        return(x)

      numerator = x - stat1(x) 
      denominator = stat2(x) - stat1(x)
      denominator[denominator == 0] <- 1

      return(numerator / denominator)
    })

    # Merge adjustment back in
    matrix_[, batch == b] <- adjusted
  }
  ## Scale back up to match overall
  matrix_ = sapply(1:row(matrix_), function(i) {
    x = matrix_[i]
    pre_stat1 = overall_stat1[i]
    pre_stat2 = overall_stat2[i]
    x * (pre_stat2 - pre_stat1) + pre_stat1
  })
  row_names(matrix_) = row_names
  matrix_
}

match_min_mean <- function(matrix_, batch) {
  #' Scales so the mins and means of each batch match
  #'
  #' @param matrix_ The matrix to batch adjust by scaling. Columns are
  #' samples, rows are features.
  #' @param batch The per-sample batch assignments.
  #' 
  #' @return The matrix_ after batch adjustment.
  #'
  #' @examples
  #' scale_adjust(data, c(rep(1, 5000), rep(2, 5000)))

  match_two_stats(
    matrix_,
    batch,
    function(x) { min(x) },
    function(x) { mean(x) }
  )
}

adjust_tampor <- function(df_, batch) {
  #' Adjusts using the tampor method.
  #'
  #' @param df_ The dataframe to batch adjust by scaling. Columns are
  #' samples, rows are features.
  #' @param batch The per-sample batch assignments.
  #' 
  #' @return The dataframe after batch adjustment.


  # From TAMPOR documentation:
  # Input is two data frames:
  # 1. (normalized or RAW) abundance. Columns are samples, rows are features.
  # 2. traits (metadata).	Columns are traits, rows are samples
  # Sample names (abundance columns, trait rows) must match exactly.

  sample_names = paste0("Sample_", seq_len(ncol(df_)))

  colnames(df_) = sample_names

  batch_as_df = data.frame(batch)
  colnames(batch_as_df) = "Batch"
  rownames(batch_as_df) = sample_names

  num_batches = length(unique(batch_as_df$Batch))

  # Save current working directory
  current_dir <- getwd()
  cat("Current working directory:", current_dir, "\n")
  # Cd into the output directory
  output_dir <- "/outputs/figures"
  setwd(output_dir)

  # If batches contain 0s, increment by 1 to avoid Tampor bug
  if (any(batch_as_df$Batch == 0)) {
    message("Incrementing batch values by 1 to avoid Tampor bug with batch values of 0.")
    batch_as_df$Batch <- batch_as_df$Batch + 1
  }

  # Print the number of negative values
  num_negatives <- sum(df_ < 0)
  message(sprintf("Number of negative values in data pre-exponentiation: %d", num_negatives))
  if (any(df_ < 0)) {
    # If there are negative values, we need to exponentiate
    message("Data contains negative values, exponentiating to make all values non-negative.")
    df_ = 2^df_
  }
  # Print the number of negative values
  num_negatives <- sum(df_ < 0)
  message(sprintf("Number of negative values in data post-exponentiation: %d", num_negatives))

  output <- TAMPOR(
    dat=df_,
    traits=batch_as_df,
    noGIS = TRUE,
    parallelThreads = num_batches,
    path = output_dir,
    iterations=500
  )$cleanRelAbun

  # Add back removed samples or features
  if (nrow(output) < nrow(df_)) {
    message(sprintf("Output has fewer rows (%d) than input (%d). Adding back removed samples.", nrow(output), nrow(df_)))
    missing_samples <- setdiff(rownames(df_), rownames(output))
    output <- rbind(output, df_[missing_samples, , drop = FALSE])
  }
  if (ncol(output) < ncol(df_)) {
    message(sprintf("Output has fewer columns (%d) than input (%d). Adding back removed features.", ncol(output), ncol(df_)))
    missing_features <- setdiff(colnames(df_), colnames(output))
    output <- cbind(output, df_[, missing_features, drop = FALSE])
  }

  # Set back to the original working directory
  setwd(current_dir)
  message(sprintf("Output has %d rows and %d columns", nrow(output), ncol(output)))
  output
}

create_design_matrix <- function(categorical) {
  #' Creates a design matrix from categorical data.
  #' @param categorical The data frame containing categorical variables.
  #' @return A design matrix suitable for use in batch adjustment.
  message_structure(categorical, "DEBUG: create_design_matrix - Initial categorical data frame")

  if (ncol(categorical) == 0) {
    stop("No categorical variables provided to create a design matrix.")
  }

  # Remove any columns with NA values in the categorical data
  na_col_mask <- colSums(is.na(categorical)) > 0
  if (any(na_col_mask)) {
    columns_with_na = colnames(categorical)[na_col_mask]
    message(sprintf("The following columns in the categorical data contain NA values. Removing these columns: %s", paste(columns_with_na, collapse = ", ")))
    categorical <- categorical[, !na_col_mask, drop = FALSE]
    message_structure(categorical, "DEBUG: create_design_matrix - Structure of 'categorical' after removing NA columns")
  }
  # Remove columns with more than 10 unique values (e.g. the column with sample names)
  columns_more_than_10 <- sapply(categorical, function(col) {length(unique(col)) > 10})
  if (any(columns_more_than_10)) {
    message(sprintf("Removing columns with more than 10 unique values: %s", paste(colnames(categorical)[columns_more_than_10], collapse = ", ")))
    categorical <- categorical[, !columns_more_than_10, drop = FALSE]
    message_structure(categorical, "DEBUG: create_design_matrix - Structure of 'categorical' after removing columns with more than 10 unique values")
  }

  design <- model.matrix(~ ., data = categorical)
  colnames(design) <- gsub("categorical", "", colnames(design))
  message_structure(design, "DEBUG: create_design_matrix - New design matrix structure")

  return(design)
}

adjust_limma <- function(x, batch, design) {
  # x is a data frame with samples as columns and probes as rows

  # Ensure the batch variable is a factor
  batch <- as.factor(batch)
  message_structure(batch, "DEBUG: adjust_limma - 'batch' vector details:")
  
  # Number of rows in batch and design should match.
  if (length(batch) != nrow(design)) {
    message("Batch and design matrix lengths do not match.")
    message(sprintf("Batch length: %d, Design rows: %d Original rows: %d", length(batch), nrow(design), ncol(x)))
    stop(sprintf("Batch length: %d, Design matrix rows: %d", length(batch), nrow(design)))
  }

  # Remove batch effects using limma's removeBatchEffect function
  limma::removeBatchEffect(x, batch = batch, design = design)
}

prep_seurat <- function(df_, batch, data_are_counts) {
  #' Prepares a Seurat object from a dataframe.
  #'
  #' @param df_ The dataframe to convert to a Seurat object. Columns are
  #' samples, rows are features.
  #' @param batch The per-sample batch assignments.
  #'
  #' @return A Seurat object ready for further processing.

  seurat_obj <- CreateSeuratObject(counts = df_, meta.data = data.frame(Batch = batch))
  
  if (data_are_counts) {
    # If the data are counts, we need to normalize it
    message("Data is not normalized. Normalizing using Seurat's NormalizeData function.")
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  } 
  
  # Find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  return(seurat_obj)
}


adjust_seurat <- function(df_, batch, data_are_counts) {
  #' Adjusts using the Seurat method.
  #'
  #' @param df_ The dataframe to batch adjust by scaling. Columns are
  #' samples, rows are features.
  #' @param batch The per-sample batch assignments.
  #' 
  #' @return The dataframe after batch adjustment.
  #' 
  #' @examples
  #' adjust_seurat(data, c(rep(1, 5000), rep(2, 5000)))

  seurat_obj <- prep_seurat(df_, batch, data_are_counts)
  
  # Scale the data
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "Batch")
  
  # Return the scaled data
  GetAssayData(seurat_obj, slot = "data")
}

adjust_harmony <- function(df_, batch, design, data_are_counts) {
  #' Adjusts using the Harmony method.
  #'
  #' @param df_ The dataframe to batch adjust by scaling. Columns are
  #' samples, rows are features.
  #' @param batch The per-sample batch assignments.
  #' @param design The design matrix to use for adjustment.
  #'
  #' @return The dataframe after batch adjustment.
  #' 
  #' @examples
  #' adjust_harmony(data, c(rep(1, 5000), rep(2, 5000)), design)
  seurat_obj <- prep_seurat(df_, batch, data_are_counts)  
  
  # Run Harmony
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "Batch", assay.use = "RNA", reduction = "pca", project.dim = FALSE)
  
  # Return the adjusted data
  GetAssayData(seurat_obj, slot = "data")
}


is.whole <- function(a, tol = 1e-7) { 
  all(abs(a - floor(a)) <= tol)
}

batch_adjust_tidy <- function(df, adjuster, batch_col = "Batch") {
  message("Separating batch column from data frame.")
  orig_col_names = colnames(df)
  batch = df[[batch_col]]
  df[[batch_col]] = NULL

  message("Separating quantitative and categorical columns.")
  is_categorical <- vapply(df, function(col) !is.numeric(col) || is.whole(col), logical(1))
  categorical <- df[, is_categorical, drop = FALSE]
  quantitative <- df[, !is_categorical, drop = FALSE]

  message_structure(categorical, "DEBUG: batch_adjust_tidy - Categorical data frame")

  if (args$full_design_matrix) {
    message("Creating full design matrix with all categorical variables.")
    design <- create_design_matrix(categorical)
  } else if (!is.null(args$column)) {
    message(sprintf("Creating design matrix with specified columns: %s", paste(args$column, collapse = ", ")))
    if (all(args$column %in% colnames(categorical))) {
      design <- model.matrix(~ ., data = categorical[, args$column, drop = FALSE])
      colnames(design) <- gsub("categorical", "", colnames(design))
    } else {
      stop(sprintf("Specified columns '%s' not found in categorical data.", paste(args$column, collapse = ", ")))
    }
  } else {
    message("Creating design matrix with intercept only.")
    design <- matrix(1, nrow = nrow(categorical), ncol = 1)
    colnames(design) <- "Intercept"
  }

  quantitative <- as.matrix(quantitative)
  quantitative <- t(quantitative)  # Transpose to have samples as columns and features as rows
  message(sprintf("Quantitative data frame transposed: %d rows, %d columns", nrow(quantitative), ncol(quantitative)))

  # Determine if data are counts or continuous
  data_are_counts = any(quantitative < 0)

  message(sprintf("Adjusting %d quantitative columns with %s method.", ncol(quantitative), adjuster))
  if (adjuster == "combat") {
    adjusted = adjust_combat(quantitative, batch, design, data_are_counts)
  } else if (adjuster == "min_mean") {
    adjusted = match_min_mean(quantitative, batch)
  } else if (adjuster == "tampor") {
    adjusted = adjust_tampor(quantitative, batch)
  } else if (adjuster == "limma") {
    adjusted = adjust_limma(quantitative, batch, design)
  } else if (adjuster == "seurat") {
    adjusted = adjust_seurat(quantitative, batch, data_are_counts)
  } else if (adjuster == "harmony") {
    adjusted = adjust_harmony(quantitative, batch, design, data_are_counts)
  }
  else {
    stop(sprintf("Unknown adjuster '%s'", adjuster))
  }
  adjusted = t(adjusted)  # Transpose back to columns as features and rows as samples
  message(sprintf("Adjusted columns: %d, rows: %d", ncol(adjusted), nrow(adjusted)))
  message_structure(adjusted, "DEBUG: batch_adjust_tidy - Adjusted quantitative data frame")
  
  message("Combining adjusted and categorical columns.")
  adjusted = cbind(batch, categorical, adjusted)
  colnames(adjusted)[1] = batch_col
  adjusted[,orig_col_names]
}


# Main execution --------------------------------

message(sprintf("\nReading input file '%s'", args$input_file))
suppressMessages(df <- vroom(args$input_file, show_col_types = FALSE))
message(sprintf("Input file has %d rows and %d columns", nrow(df), ncol(df)))

if (!(args$batch_col %in% names(df))) {
  discrete_col_names <- df %>%
    select_if(~!is.numeric(.) || is.whole(.)) %>%
    names()

  error_message <- sprintf(
    "--batch-col argument (default 'Batch', selected '%s') must be a column in 'input_path' csv. Options: [%s]",
    args$batch_col,
    paste(discrete_col_names, collapse = ", ")
  )
  stop(error_message)
}

message(sprintf("Batch adjust tidy"))
batch_adjust_tidy(
  df, 
  batch_col = args$batch_col,
  adjuster = args$adjuster
) %>% 
  # Make sure to only keep 6 decimal places for floating point numbers
  mutate(across(where(is.numeric), ~round(., 6))) %>%
  write_csv(args$output_file)

message(sprintf("Saved output to '%s'", args$output_file))
