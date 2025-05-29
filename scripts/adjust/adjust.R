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
})

# Parse command line args --------------------------

parser <- ArgumentParser()

parser$add_argument("input_file", help = "path to the input file")
parser$add_argument("output_file", help = "path to the output file.")
parser$add_argument("-a", "--adjuster", default = "combat", choices = c("combat", "min_mean", "tampor"), help = "method to use for adjustment")
parser$add_argument("-b", "--batch-col", default = "Batch", help = "title of batch column to adjust for")

args <- parser$parse_args()

# Define functions ---------------------------------

ComBat_ignore_nonvariance <- function(matrix_, batch) {
  #' Run ComBat and ignore nonvarying features.
  #'
  #' ComBat requires that all features have some variance (and probably assumes
  #' that all features are normally distributed). Since some features don't
  #' vary across samples, this function ignores nonvarying features before
  #' running ComBat.
  #'
  #' @param matrix_ The matrix to batch adjust with ComBat. Columns are features,
  #' rows are samples.
  #' @param batch The per-sample batch assignments. See the ComBat function for
  #' more information.
  #' 
  #' @return The matrix_ after batch adjustment.
  #'
  #' @examples
  #' ComBat_ignore_nonvariance(data, c(rep(1, 5000), rep(2, 5000)))
  matrix_ <- t(matrix_)

  varying_row_mask <- apply(matrix_, 1, function(x) { length(unique(x)) > 1 })

  matrix_[varying_row_mask,] <- ComBat(matrix_[varying_row_mask,], batch)

  t(matrix_)
}

match_two_stats <- function(matrix_, batch, stat1, stat2) {
  #' Matches batches by scaling so that two statistics are equal.
  #'
  #' @param matrix_ The matrix to batch adjust by scaling. Columns are
  #' features, rows are samples.
  #' @param batch The per-sample batch assignments.
  #' @param stat1 The first statistic to match. (as a function, not a string)
  #' @param stat2 The second statistic to match.
  #'
  #' @return The matrix_ after batch adjustment.
  #'
  #' @examples
  #' match_two_stats(data, c(rep(1, 5000), rep(2, 5000)), "mean", "sd")
  column_names = colnames(matrix_)
  # Get columnwise stats
  overall_stat1 <- apply(matrix_, 2, stat1)
  overall_stat2 <- apply(matrix_, 2, stat2)

  # Scale each batch individually so that stat1 is 0 and stat2 is 1
  for (b in levels(factor(batch))) {
    # drop=F makes it return a matrix when you only grab one row.
    batch_rows <- matrix_[batch == b, , drop = FALSE]

    if (nrow(batch_rows) <= 1) {
      stop(sprintf("Can't scale columns: batch '%s' has <= 1 sample.", b))
    }

    adjusted = apply(batch_rows, 2, function(x) {
      if (all(x == 0))
        return(x)

      numerator = x - stat1(x) 
      denominator = stat2(x) - stat1(x)
      denominator[denominator == 0] <- 1

      return(numerator / denominator)
    })

    # Merge adjustment back in
    matrix_[batch == b] = adjusted
  }
  ## Scale back up to match overall
  matrix_ = sapply(1:ncol(matrix_), function(i) {
    x = matrix_[,i]
    pre_stat1 = overall_stat1[i]
    pre_stat2 = overall_stat2[i]
    x * (pre_stat2 - pre_stat1) + pre_stat1
  })
  colnames(matrix_) = column_names
  matrix_
}

match_min_mean <- function(matrix_, batch) {
  #' Scales so the mins and means of each batch match
  #'
  #' @param matrix_ The matrix to batch adjust by scaling. Columns are
  #' features, rows are samples.
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
  #' features, rows are samples.
  #' @param batch The per-sample batch assignments.
  #' 
  #' @return The dataframe after batch adjustment.
  #'
  #' @examples
  #' adjust_tampor(data, c(rep(1, 5000), rep(2, 5000)))


  # From TAMPOR documentation:
  # Input is two data frames:
  # 1. (normalized or RAW) abundance. Columns are samples, rows are features.
  # 2. traits (metadata).	Columns are traits, rows are samples
  # Sample names (abundance columns, trait rows) must match exactly.

  sample_names = paste0("Sample_", seq_len(nrow(df_)))

  transposed = t(df_)
  colnames(transposed) = sample_names

  batch_as_df = data.frame(batch)
  colnames(batch_as_df) = "Batch"
  rownames(batch_as_df) = sample_names

  num_batches = length(unique(batch_as_df$Batch))


  if (any(transposed < 0)) {
    # If there are negative values, we need to exponentiate
    message("Data contains negative values, exponentiating to make all values non-negative.")
    transposed = 2^transposed
  }

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


  output <- TAMPOR(
    dat=transposed,
    traits=batch_as_df,
    noGIS = TRUE,
    parallelThreads = num_batches,
    path = output_dir,
    iterations=500
  )$cleanRelAbun

  # Set back to the original working directory
  setwd(current_dir)
  # Transpose back to original format and return
  t(output)
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
  # Check if the column is numeric and has no non-whole numbers
  # Using vapply for speed, logical(1) is the return type
  is_categorical <- vapply(df, function(col) !is.numeric(col) || is.whole(col), logical(1))
  categorical <- df[, is_categorical, drop = FALSE]
  quantitative <- df[, !is_categorical, drop = FALSE]

  message("Adjusting using the '%s' adjuster", adjuster)
  if (adjuster == "combat") {
    adjusted = ComBat_ignore_nonvariance(as.matrix(quantitative), batch)
  } else if (adjuster == "min_mean") {
    adjusted = match_min_mean(as.matrix(quantitative), batch)
  } else if (adjuster == "tampor") {
    adjusted = adjust_tampor(as.matrix(quantitative), batch)
  } else {
    stop(sprintf("Unknown adjuster '%s'", adjuster))
  }
  
  message("Combining adjusted and categorical columns.")
  adjusted = cbind(batch, categorical, adjusted)
  colnames(adjusted)[1] = batch_col
  adjusted[,orig_col_names]
}

message("Reading input file.")

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
) %>% write_csv(args$output_file)

message(sprintf("Saved output to '%s'", args$output_file))
