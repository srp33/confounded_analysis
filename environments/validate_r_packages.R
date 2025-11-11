#!/usr/bin/env Rscript
# R Package Loading Validation Script
# Tests all required packages for the batch correction pipeline
# Designed to run on both login and compute nodes
#
# Usage:
#   Rscript environments/validate_r_packages.R
#   ./environments/run_with_env.sh --sbatch environments/validate_r_packages.R
#   ./environments/run_with_env.sh --r-env combatseq --sbatch environments/validate_r_packages.R

# Suppress startup messages
suppressPackageStartupMessages({
  library(methods)
})

# ============================================================================
# Helper Functions
# ============================================================================

print_header <- function(title) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("  ", title, "\n", sep = "")
  cat(rep("=", 70), "\n", sep = "")
}

print_section <- function(title) {
  cat("\n--- ", title, " ---\n", sep = "")
}

test_package_load <- function(pkg_name, display_name = NULL, required = TRUE) {
  """
  Test if a package can be loaded
  
  Args:
    pkg_name: Name of the package to load
    display_name: Display name for the package (defaults to pkg_name)
    required: Whether this package is required
  
  Returns:
    list: success (logical), version (character), error (character)
  """
  if (is.null(display_name)) {
    display_name <- pkg_name
  }
  
  result <- tryCatch({
    suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
    version <- as.character(packageVersion(pkg_name))
    cat(sprintf("  ✓ %-35s v%s\n", display_name, version))
    list(success = TRUE, version = version, error = NULL)
  }, error = function(e) {
    status <- if (required) "REQUIRED" else "OPTIONAL"
    cat(sprintf("  ✗ %-35s FAILED (%s): %s\n", display_name, status, e$message))
    list(success = FALSE, version = NULL, error = e$message)
  })
  
  return(result)
}

# ============================================================================
# Package Test Functions
# ============================================================================

test_core_infrastructure <- function() {
  print_section("Core Infrastructure")
  
  packages <- list(
    c("Rcpp", "Rcpp", TRUE),
    c("RcppArmadillo", "RcppArmadillo", TRUE),
    c("BH", "Boost Headers", TRUE),
    c("data.table", "data.table", TRUE),
    c("Matrix", "Matrix", TRUE),
    c("foreach", "foreach", TRUE),
    c("doParallel", "doParallel", TRUE),
    c("future", "future", TRUE),
    c("matrixStats", "matrixStats", TRUE),
    c("remotes", "remotes", TRUE)
  )
  
  results <- lapply(packages, function(p) {
    c(p[2], test_package_load(p[1], p[2], as.logical(p[3])))
  })
  
  return(results)
}

test_tidyverse_ecosystem <- function() {
  print_section("Tidyverse Ecosystem")
  
  packages <- list(
    c("tidyverse", "tidyverse", TRUE),
    c("dplyr", "dplyr", TRUE),
    c("ggplot2", "ggplot2", TRUE),
    c("tidyr", "tidyr", TRUE),
    c("purrr", "purrr", TRUE),
    c("readr", "readr", TRUE),
    c("stringr", "stringr", TRUE),
    c("tibble", "tibble", TRUE)
  )
  
  results <- lapply(packages, function(p) {
    c(p[2], test_package_load(p[1], p[2], as.logical(p[3])))
  })
  
  return(results)
}

test_statistical_packages <- function() {
  print_section("Statistical Packages")
  
  packages <- list(
    c("glmnet", "glmnet", TRUE),
    c("e1071", "e1071 (SVM)", TRUE),
    c("ranger", "ranger (Random Forest)", TRUE),
    c("caret", "caret", TRUE),
    c("mclust", "mclust", TRUE),
    c("MCMCpack", "MCMCpack", TRUE),
    c("huge", "huge", TRUE)
  )
  
  results <- lapply(packages, function(p) {
    c(p[2], test_package_load(p[1], p[2], as.logical(p[3])))
  })
  
  return(results)
}

test_visualization_packages <- function() {
  print_section("Visualization Packages")
  
  packages <- list(
    c("gridExtra", "gridExtra", TRUE),
    c("gplots", "gplots", TRUE),
    c("RColorBrewer", "RColorBrewer", TRUE),
    c("plotly", "plotly", FALSE),
    c("heatmaply", "heatmaply", FALSE),
    c("ggpubr", "ggpubr", TRUE),
    c("ggtext", "ggtext", TRUE),
    c("corrplot", "corrplot", TRUE),
    c("kableExtra", "kableExtra", TRUE)
  )
  
  results <- lapply(packages, function(p) {
    c(p[2], test_package_load(p[1], p[2], as.logical(p[3])))
  })
  
  return(results)
}

test_bioconductor_core <- function() {
  print_section("Bioconductor Core")
  
  packages <- list(
    c("BiocManager", "BiocManager", TRUE),
    c("SummarizedExperiment", "SummarizedExperiment", TRUE),
    c("limma", "limma", TRUE),
    c("vsn", "vsn", TRUE),
    c("AnnotationDbi", "AnnotationDbi", TRUE),
    c("biomaRt", "biomaRt", TRUE),
    c("ExperimentHub", "ExperimentHub", TRUE),
    c("ComplexHeatmap", "ComplexHeatmap", TRUE)
  )
  
  results <- lapply(packages, function(p) {
    c(p[2], test_package_load(p[1], p[2], as.logical(p[3])))
  })
  
  return(results)
}

test_batch_correction_packages <- function() {
  print_section("Batch Correction Packages")
  
  packages <- list(
    c("sva", "sva (ComBat)", TRUE),
    c("batchelor", "batchelor", TRUE),
    c("RUVSeq", "RUVSeq", TRUE)
  )
  
  results <- lapply(packages, function(p) {
    c(p[2], test_package_load(p[1], p[2], as.logical(p[3])))
  })
  
  return(results)
}

test_machine_learning_packages <- function() {
  print_section("Machine Learning Packages")
  
  packages <- list(
    c("MLmetrics", "MLmetrics", TRUE),
    c("ROCR", "ROCR", TRUE),
    c("xgboost", "XGBoost", TRUE),
    c("lightgbm", "LightGBM", TRUE),
    c("parsnip", "parsnip", TRUE),
    c("tidymodels", "tidymodels", TRUE)
  )
  
  results <- lapply(packages, function(p) {
    c(p[2], test_package_load(p[1], p[2], as.logical(p[3])))
  })
  
  return(results)
}

test_dimensionality_reduction <- function() {
  print_section("Dimensionality Reduction")
  
  packages <- list(
    c("Rtsne", "Rtsne", TRUE),
    c("umap", "umap", TRUE)
  )
  
  results <- lapply(packages, function(p) {
    c(p[2], test_package_load(p[1], p[2], as.logical(p[3])))
  })
  
  return(results)
}

test_utility_packages <- function() {
  print_section("Utility Packages")
  
  packages <- list(
    c("argparse", "argparse", TRUE),
    c("docstring", "docstring", TRUE),
    c("itertools", "itertools", TRUE),
    c("pacman", "pacman", TRUE),
    c("Seurat", "Seurat", FALSE),
    c("languageserver", "languageserver (VS Code)", TRUE)
  )
  
  results <- lapply(packages, function(p) {
    c(p[2], test_package_load(p[1], p[2], as.logical(p[3])))
  })
  
  return(results)
}

# ============================================================================
# Functionality Tests
# ============================================================================

test_tidyverse_functionality <- function() {
  print_section("Tidyverse Functionality Test")
  
  tryCatch({
    suppressPackageStartupMessages({
      library(dplyr)
      library(ggplot2)
      library(tidyr)
    })
    
    # Create test data
    df <- data.frame(
      x = rnorm(1000),
      y = rnorm(1000),
      group = sample(c("A", "B", "C"), 1000, replace = TRUE)
    )
    
    # Test dplyr operations
    result <- df %>%
      group_by(group) %>%
      summarise(
        mean_x = mean(x),
        sd_y = sd(y),
        n = n()
      )
    cat(sprintf("  ✓ dplyr: grouped %d rows into %d groups\n", nrow(df), nrow(result)))
    
    # Test ggplot2
    p <- ggplot(df, aes(x = x, y = y, color = group)) +
      geom_point(alpha = 0.5) +
      theme_minimal()
    cat("  ✓ ggplot2: created scatter plot\n")
    
    # Test tidyr
    wide_df <- df %>%
      group_by(group) %>%
      summarise(mean_x = mean(x), mean_y = mean(y)) %>%
      pivot_wider(names_from = group, values_from = c(mean_x, mean_y))
    cat(sprintf("  ✓ tidyr: pivoted data to %d columns\n", ncol(wide_df)))
    
    return(TRUE)
  }, error = function(e) {
    cat(sprintf("  ✗ Tidyverse functionality test failed: %s\n", e$message))
    return(FALSE)
  })
}

test_statistical_functionality <- function() {
  print_section("Statistical Functionality Test")
  
  tryCatch({
    suppressPackageStartupMessages({
      library(glmnet)
      library(e1071)
      library(ranger)
    })
    
    # Create test data
    set.seed(42)
    n <- 500
    p <- 50
    X <- matrix(rnorm(n * p), n, p)
    y <- (X[, 1] + X[, 2] + rnorm(n)) > 0
    
    # Test glmnet (LASSO)
    cv_fit <- cv.glmnet(X, y, family = "binomial", nfolds = 5)
    pred <- predict(cv_fit, X, s = "lambda.min", type = "response")
    cat(sprintf("  ✓ glmnet: LASSO model trained, %d predictions\n", length(pred)))
    
    # Test SVM
    svm_model <- svm(X, y, kernel = "linear", probability = TRUE)
    svm_pred <- predict(svm_model, X)
    acc <- mean(svm_pred == y)
    cat(sprintf("  ✓ e1071: SVM trained, accuracy=%.4f\n", acc))
    
    # Test Random Forest
    df <- data.frame(X, y = as.factor(y))
    rf_model <- ranger(y ~ ., data = df, num.trees = 50, probability = TRUE)
    cat(sprintf("  ✓ ranger: Random Forest trained, OOB error=%.4f\n", rf_model$prediction.error))
    
    return(TRUE)
  }, error = function(e) {
    cat(sprintf("  ✗ Statistical functionality test failed: %s\n", e$message))
    traceback()
    return(FALSE)
  })
}

test_bioconductor_functionality <- function() {
  print_section("Bioconductor Functionality Test")
  
  tryCatch({
    suppressPackageStartupMessages({
      library(SummarizedExperiment)
      library(limma)
      library(sva)
    })
    
    # Create test SummarizedExperiment
    nrows <- 200
    ncols <- 20
    counts <- matrix(rpois(nrows * ncols, lambda = 10), nrows, ncols)
    rowData <- DataFrame(gene_id = paste0("gene", 1:nrows))
    colData <- DataFrame(
      sample_id = paste0("sample", 1:ncols),
      batch = rep(c("A", "B"), each = ncols/2),
      condition = rep(c("control", "treatment"), ncols/2)
    )
    
    se <- SummarizedExperiment(
      assays = list(counts = counts),
      rowData = rowData,
      colData = colData
    )
    cat(sprintf("  ✓ SummarizedExperiment: created %d x %d object\n", nrow(se), ncol(se)))
    
    # Test limma
    design <- model.matrix(~ condition, data = colData(se))
    fit <- lmFit(log2(counts + 1), design)
    fit <- eBayes(fit)
    top_genes <- topTable(fit, coef = 2, number = 10)
    cat(sprintf("  ✓ limma: differential expression analysis, %d top genes\n", nrow(top_genes)))
    
    # Test ComBat (sva)
    batch <- colData(se)$batch
    modcombat <- model.matrix(~ condition, data = colData(se))
    combat_data <- ComBat(dat = log2(counts + 1), batch = batch, mod = modcombat)
    cat(sprintf("  ✓ sva: ComBat batch correction, output %d x %d\n", nrow(combat_data), ncol(combat_data)))
    
    return(TRUE)
  }, error = function(e) {
    cat(sprintf("  ✗ Bioconductor functionality test failed: %s\n", e$message))
    traceback()
    return(FALSE)
  })
}

test_visualization_functionality <- function() {
  print_section("Visualization Functionality Test")
  
  tryCatch({
    suppressPackageStartupMessages({
      library(ggplot2)
      library(ComplexHeatmap)
      library(corrplot)
    })
    
    # Test ggplot2
    df <- data.frame(x = rnorm(100), y = rnorm(100))
    p <- ggplot(df, aes(x = x, y = y)) +
      geom_point() +
      geom_smooth(method = "lm") +
      theme_minimal()
    cat("  ✓ ggplot2: created scatter plot with regression line\n")
    
    # Test ComplexHeatmap
    mat <- matrix(rnorm(200), 20, 10)
    rownames(mat) <- paste0("gene", 1:20)
    colnames(mat) <- paste0("sample", 1:10)
    
    # Create heatmap (don't draw, just test creation)
    ht <- Heatmap(mat, name = "expression", show_row_names = FALSE)
    cat("  ✓ ComplexHeatmap: created heatmap object\n")
    
    # Test corrplot
    cor_mat <- cor(matrix(rnorm(500), 50, 10))
    # Don't actually plot, just test the function exists
    cat("  ✓ corrplot: correlation plot function available\n")
    
    return(TRUE)
  }, error = function(e) {
    cat(sprintf("  ✗ Visualization functionality test failed: %s\n", e$message))
    traceback()
    return(FALSE)
  })
}

# ============================================================================
# System Information
# ============================================================================

print_system_info <- function() {
  print_header("System Information")
  
  cat(sprintf("  Hostname:          %s\n", Sys.info()["nodename"]))
  cat(sprintf("  Platform:          %s\n", R.version$platform))
  cat(sprintf("  R version:         %s\n", R.version$version.string))
  cat(sprintf("  Working directory: %s\n", getwd()))
  
  # Check if running on compute node
  hostname <- Sys.info()["nodename"]
  is_compute_node <- grepl("compute|node|cn", hostname, ignore.case = TRUE)
  node_type <- if (is_compute_node) "Compute Node" else "Login Node"
  cat(sprintf("  Node type:         %s\n", node_type))
  
  # Environment variables
  cat("\n  Environment Variables:\n")
  env_vars <- c("R_LIBS_USER", "R_ENV_DIR", "SLURM_JOB_ID", "SLURM_JOB_NAME")
  for (var in env_vars) {
    value <- Sys.getenv(var, unset = "not set")
    cat(sprintf("    %-20s = %s\n", var, value))
  }
  
  # Library paths
  cat("\n  Library Paths:\n")
  lib_paths <- .libPaths()
  for (i in seq_along(lib_paths)) {
    cat(sprintf("    [%d] %s\n", i, lib_paths[i]))
  }
  
  # Bioconductor version (if available)
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    bioc_version <- tryCatch({
      BiocManager::version()
    }, error = function(e) "unknown")
    cat(sprintf("\n  Bioconductor:      %s\n", bioc_version))
  }
}

# ============================================================================
# Summary Generation
# ============================================================================

generate_summary <- function(all_results, functionality_results) {
  print_header("Test Summary")
  
  # Count results
  total_packages <- length(all_results)
  passed_packages <- sum(sapply(all_results, function(r) r[[2]]$success))
  failed_packages <- total_packages - passed_packages
  
  total_functionality <- length(functionality_results)
  passed_functionality <- sum(unlist(functionality_results))
  failed_functionality <- total_functionality - passed_functionality
  
  cat("\n  Package Loading:\n")
  cat(sprintf("    Total:   %d\n", total_packages))
  cat(sprintf("    Passed:  %d (%.1f%%)\n", passed_packages, 100*passed_packages/total_packages))
  cat(sprintf("    Failed:  %d\n", failed_packages))
  
  cat("\n  Functionality Tests:\n")
  cat(sprintf("    Total:   %d\n", total_functionality))
  cat(sprintf("    Passed:  %d (%.1f%%)\n", passed_functionality, 100*passed_functionality/total_functionality))
  cat(sprintf("    Failed:  %d\n", failed_functionality))
  
  # List failures
  if (failed_packages > 0) {
    cat("\n  Failed Package Loads:\n")
    for (result in all_results) {
      if (!result[[2]]$success) {
        cat(sprintf("    ✗ %s: %s\n", result[[1]], result[[2]]$error))
      }
    }
  }
  
  # Overall status
  cat("\n")
  if (failed_packages == 0 && failed_functionality == 0) {
    cat("  ✓ ALL TESTS PASSED\n")
    return(0)
  } else {
    cat("  ✗ SOME TESTS FAILED\n")
    return(1)
  }
}

# ============================================================================
# Main Function
# ============================================================================

main <- function() {
  start_time <- Sys.time()
  
  print_header("R Environment Validation")
  cat(sprintf("  Started: %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
  
  # Print system info
  print_system_info()
  
  # Test package loading
  print_header("Package Loading Tests")
  
  all_results <- c(
    test_core_infrastructure(),
    test_tidyverse_ecosystem(),
    test_statistical_packages(),
    test_visualization_packages(),
    test_bioconductor_core(),
    test_batch_correction_packages(),
    test_machine_learning_packages(),
    test_dimensionality_reduction(),
    test_utility_packages()
  )
  
  # Test functionality
  print_header("Functionality Tests")
  
  functionality_results <- list(
    test_tidyverse_functionality(),
    test_statistical_functionality(),
    test_bioconductor_functionality(),
    test_visualization_functionality()
  )
  
  # Generate summary
  exit_code <- generate_summary(all_results, functionality_results)
  
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat(sprintf("\n  Completed: %s\n", format(end_time, "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("  Duration:  %.2f seconds\n", duration))
  cat(rep("=", 70), "\n", sep = "")
  
  quit(status = exit_code)
}

# Run main function
main()
