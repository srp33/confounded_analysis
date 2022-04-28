#Debugging: 
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

# Load dependencies --------------------------------
library(dplyr)
library(readr)
library(docstring)
library(stringr)
library(argparse)
library(sva)

# Parse command line args --------------------------

parser <- ArgumentParser()

parser$add_argument("input_file", help = "path to the input file")
parser$add_argument("output_file", help = "path to the output file.")
parser$add_argument("-a", "--adjuster", default = "combat", choices = c("combat", "scale"), help = "method to use for adjustment")
parser$add_argument("-b", "--batch-col", default = "Batch", help = "title of batch column to adjust for")

args <- parser$parse_args()

# Define functions ---------------------------------

is.whole <- function(a, tol = 1e-7) { 
  # Snatched from https://stat.ethz.ch/pipermail/r-help/2003-April/032471.html
  is.eq <- function(x,y) { 
    r <- all.equal(x,y, tol=tol)
    is.logical(r) && r 
  }
  (is.numeric(a) && is.eq(a, floor(a))) ||
    (is.complex(a) && {ri <- c(Re(a),Im(a)); is.eq(ri, floor(ri))})
}

ComBat_ignore_nonvariance <- function(matrix_, batch)
{
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

scale_adjust <- function(matrix_, batch)
{
  #' Run ComBat and ignore nonvarying features.
  #'
  #' @param matrix_ The matrix to batch adjust by scaling. Columns are
  #' features, rows are samples.
  #' @param batch The per-sample batch assignments.
  #' 
  #' @return The matrix_ after batch adjustment.
  #'
  #' @examples
  #' scale_adjust(data, c(rep(1, 5000), rep(2, 5000)))

  column_names = colnames(matrix_)
  
  # Get columnwise mins & maxes
  mins <- apply(matrix_, 2, min)
  maxes <- apply(matrix_, 2, max)

  # Scale each batch individually to [0, 1]
  for (b in levels(factor(batch))) {
    # drop=F makes it return a matrix when you only grab one row.
    batch_rows <- matrix_[batch == b, , drop = FALSE]

    if (nrow(batch_rows) <= 1) {
      stop(sprintf("Can't scale columns: batch '%s' has <= 1 sample.", b))
    }

    adjusted = apply(batch_rows, 2, function(x) {
      if (all(x == 0))
        return(x)

      numerator = x - min(x)
      denominator = max(x) - min(x)

      return(numerator / denominator)
    })

    # Merge adjustment back in
    matrix_[batch == b] = adjusted
  }

  ## Scale back up to [min, max]
  matrix_ = sapply(1:ncol(matrix_), function(i) {
    x = matrix_[,i]
    pre_min = mins[i]
    pre_max = maxes[i]
    x * (pre_max - pre_min) + pre_min
  })

  colnames(matrix_) = column_names
  matrix_
}

batch_adjust_tidy <- function(df, adjuster, batch_col = "Batch") {
  orig_col_names = colnames(df)
  batch = pull(df, batch_col)
  df = select(df, -batch_col)

  categorical <- df %>%
    select_if(~!is.numeric(.) || is.whole(.))
  quantitative <- df %>%
    select_if(~is.numeric(.) && !is.whole(.))

  if (adjuster == "combat") {
    adjusted = ComBat_ignore_nonvariance(as.matrix(quantitative), batch)
  } else {
    adjusted = scale_adjust(as.matrix(quantitative), batch)
  }
  
  adjusted = cbind(batch, categorical, adjusted)
  colnames(adjusted)[1] = batch_col
  adjusted[,orig_col_names]
}

message("Reading input file.")

suppressMessages(df <- read_csv(args$input_file))

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

message(sprintf("Adjusting using the '%s' adjuster", args$adjuster))
batch_adjust_tidy(
  df, 
  batch_col = args$batch_col,
  adjuster = args$adjuster
) %>% write_csv(args$output_file)

message(sprintf("Saved output to '%s'", args$output_file))
