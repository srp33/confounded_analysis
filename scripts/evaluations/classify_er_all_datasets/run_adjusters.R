#!/usr/bin/env Rscript

# ============================================================
# inter_source_analysis.R
# Train on one dataset, test on another using scikit-learn
# ============================================================

# ---- Libraries ----
suppressPackageStartupMessages({
  library(reticulate)
  library(readr)
  library(dplyr)
  library(optparse)
})

# ---- Parse Arguments ----
option_list <- list(
  make_option("--input", type="character", help="Path to combined data file (.csv)"),
  make_option("--output", type="character", help="Path to output results file (.csv)"),
  make_option("--adjuster", type="character", help="Adjustment method name"),
  make_option("--train", type="character", help="Training dataset meta_source name"),
  make_option("--test", type="character", help="Testing dataset meta_source name"),
  make_option("--n_repeats", type="integer", default=10, help="Number of repeats"),
  make_option("--python_env", type="character", default=NULL, help="Path to Python env with sklearn")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ---- Setup Python Environment ----
if (!is.null(opt$python_env)) {
  use_python(opt$python_env, required = TRUE)
}

# Import scikit-learn modules
sklearn <- import("sklearn")
np <- import("numpy")
metrics <- sklearn$metrics
ensemble <- sklearn$ensemble

# ---- Helper: run sklearn model ----
run_sklearn_model <- function(X_train, y_train, X_test, y_test, random_state = 42) {
  # Convert to numpy arrays
  X_train_np <- np$array(as.matrix(X_train))
  y_train_np <- np$array(as.integer(y_train))
  X_test_np <- np$array(as.matrix(X_test))
  y_test_np <- np$array(as.integer(y_test))

  # Model: HistGradientBoostingClassifier (same as your Python code)
  model <- ensemble$HistGradientBoostingClassifier(max_iter = as.integer(100), random_state = as.integer(random_state))
  model$fit(X_train_np, y_train_np)

  y_pred <- model$predict(X_test_np)
  y_proba <- model$predict_proba(X_test_np)[, 2]

  # Metrics (matching Python's calculate_metrics)
  tryCatch({
    cm <- metrics$confusion_matrix(y_test_np, y_pred, labels = np$array(c(0L, 1L)))
    tn <- cm[1, 1]; fp <- cm[1, 2]; fn <- cm[2, 1]; tp <- cm[2, 2]
    acc <- metrics$accuracy_score(y_test_np, y_pred)
    auc <- tryCatch(metrics$roc_auc_score(y_test_np, y_proba), error = function(e) NA)
    sens <- ifelse(tp + fn > 0, tp / (tp + fn), NA)
    spec <- ifelse(tn + fp > 0, tn / (tn + fp), NA)
    mcc <- metrics$matthews_corrcoef(y_test_np, y_pred)

    data.frame(
      Accuracy = acc,
      `ROC AUC` = auc,
      Sensitivity = sens,
      Specificity = spec,
      MCC = mcc,
      `True Negative` = tn,
      `False Positive` = fp,
      `False Negative` = fn,
      `True Positive` = tp,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning(paste("Metric calculation failed:", e))
    return(NULL)
  })
}

# ---- Load Data ----
cat("Loading data from:", opt$input, "\n")
df <- read.csv(opt$input, stringsAsFactors = FALSE)

if (!all(c("meta_source", "meta_er_status") %in% colnames(df))) {
  stop("Input data must contain columns: meta_source, meta_er_status")
}

# ---- Filter train/test ----
df_train <- df %>% filter(meta_source == opt$train)
df_test  <- df %>% filter(meta_source == opt$test)

if (nrow(df_train) == 0 || nrow(df_test) == 0) {
  stop(paste("No samples found for train =", opt$train, "or test =", opt$test))
}

# Drop non-numeric/meta columns
cols_to_drop <- c("meta_er_status", "meta_source")
X_train <- df_train %>% select(where(is.numeric))
X_test  <- df_test  %>% select(where(is.numeric))
y_train <- df_train$meta_er_status
y_test  <- df_test$meta_er_status

# ---- Run analysis for each repeat ----
results_list <- list()
cat("Running inter-source analysis for", opt$train, "→", opt$test, "with adjuster:", opt$adjuster, "\n")

for (i in seq_len(opt$n_repeats)) {
  seed <- 42 + i
  metrics_df <- run_sklearn_model(X_train, y_train, X_test, y_test, random_state = seed)
  metrics_df$train <- opt$train
  metrics_df$test <- opt$test
  metrics_df$classifier <- "HistGradientBoosting"
  metrics_df$adjuster <- opt$adjuster
  metrics_df$repeat <- i
  results_list[[i]] <- metrics_df
}

results_all <- bind_rows(results_list)

# ---- Save results ----
if (!dir.exists(dirname(opt$output))) dir.create(dirname(opt$output), recursive = TRUE)
write_csv(results_all, opt$output)

cat("Finished inter-source analysis\n")
cat("Output saved to:", opt$output, "\n")
