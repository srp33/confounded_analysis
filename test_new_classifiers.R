#!/usr/bin/env Rscript

# Simple test script to verify KNN and XGBoost functions work
cat("Testing new classifier functions...\n")

# Load helper functions
source("scripts/evaluations/book_chapter/scripts/helper.R")

# Test getPredFunctions
cat("Testing getPredFunctions...\n")
knn_func <- getPredFunctions("knn")
xgboost_func <- getPredFunctions("xgboost")
cat("✓ getPredFunctions works for both knn and xgboost\n")

# Create simple test data
set.seed(42)
n_genes <- 100
n_samples <- 20
test_data <- matrix(rnorm(n_genes * n_samples), nrow = n_genes, ncol = n_samples)
test_labels <- rep(c(0, 1), each = n_samples/2)

cat("Created test data: ", n_genes, " genes, ", n_samples, " samples\n")

# Test KNN function
cat("Testing KNN function...\n")
tryCatch({
  knn_result <- predKNN_pp(test_data[, 1:15], test_data[, 16:20], test_labels[1:15])
  cat("✓ KNN function works\n")
  cat("  Training predictions length:", length(knn_result$pred_trn_prob), "\n")
  cat("  Test predictions length:", length(knn_result$pred_tst_prob), "\n")
}, error = function(e) {
  cat("✗ KNN function failed:", e$message, "\n")
})

# Test XGBoost function
cat("Testing XGBoost function...\n")
tryCatch({
  xgb_result <- predXGBoost_pp(test_data[, 1:15], test_data[, 16:20], test_labels[1:15])
  cat("✓ XGBoost function works\n")
  cat("  Training predictions length:", length(xgb_result$pred_trn_prob), "\n")
  cat("  Test predictions length:", length(xgb_result$pred_tst_prob), "\n")
}, error = function(e) {
  cat("✗ XGBoost function failed:", e$message, "\n")
})

cat("Test completed!\n")