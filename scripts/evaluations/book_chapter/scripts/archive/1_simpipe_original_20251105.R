# Simulation Pipeline

# Suppress all output and warnings
options(warn = -1)
suppressMessages(suppressWarnings({
  rm(list=ls())
  setwd("/scripts/evaluations/book_chapter")
  
  if(!dir.exists("/scripts/evaluations/book_chapter/results")){
    dir.create("/scripts/evaluations/robustifying/results")
  }
  
  required_packages <- c("SummarizedExperiment", "plyr", "sva", "MCMCpack", "ROCR", "ggplot2", 
                        "limma", "nnls", "glmnet", "rpart", "genefilter", "nnet", "e1071", 
                        "RcppArmadillo", "foreach", "parallel", "doParallel", "ranger", "scales",
                        "purrr", "dplyr", "lightgbm", "batchelor")
  
  package_results <- sapply(required_packages, require, character.only=TRUE, quietly=TRUE)
  failed_packages <- names(package_results)[!package_results]
  if(length(failed_packages) > 0) {
    cat("Warning: Failed to load packages:", paste(failed_packages, collapse=", "), "\n")
  }
  cat("Package loading results:", all(package_results), "\n")
}))

source("/scripts/evaluations/book_chapter/scripts/helper.R")
load("/scripts/evaluations/book_chapter/data/combined_sub.RData")

# Get allocated cores from SLURM or default to 1
get_allocated_cores <- function() {
  return(max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))))
}

####  Load and prepare data
# Select 1000 genes with largest variance in training set (Africa)
var_trn <- rowVars(train_expr)
genes_sel <- rownames(train_expr)[order(var_trn, decreasing=TRUE)[1:1000]]
train_expr <- train_expr[genes_sel, ]
test_expr <- test_expr[genes_sel, ]

####  Parse parameters 
command_args <- commandArgs(trailingOnly=TRUE)  
if(length(command_args)!=3){stop("Not enough input parameters!")}

## Degree of batch effect (strength of signal)
N_batch <- 3
N_sample_size <- as.numeric(command_args[1])   
max_batch_mean <- as.numeric(command_args[2]) 
max_batch_var <- as.numeric(command_args[3]) 

## Prediction model
learner_types <- c("logistic", "elnet", "svm", "rf", "lightgbm", "nnet")

hyper_pars <- list(
  hyper_mu=seq(from=-max_batch_mean, to=max_batch_mean, length.out=N_batch),  
  hyper_sd=sqrt(rep(0.01, N_batch)),
  hyper_alpha=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                    v=rep(0.01, N_batch))$alpha,
  hyper_beta=mv2ab(m=seq(from=1/max_batch_var, to=max_batch_var, length.out=N_batch), 
                   v=rep(0.01, N_batch))$beta
)

## Pipeline parameters
iterations <- 100
norm_data <- TRUE

exp_name <- sprintf('batchN%s_m%s_v%s', N_sample_size, 
                    gsub('.', '', max_batch_mean, fixed=T), gsub('.', '', max_batch_var, fixed=T))  
perf_measures <- c("mxe", "auc")    

####  Pipeline Functions  ####

#' Process single iteration of the simulation
#' @param iteration_id The iteration number
#' @param train_expr Training expression data
#' @param test_expr Test expression data  
#' @param y_train Training labels
#' @param y_test Test labels
#' @param learner_types Vector of learner types to test
#' @param perf_measures Vector of performance measures
#' @param other_params List of other parameters
process_single_iteration <- function(iteration_id, train_expr, test_expr, y_train, y_test,
                                   learner_types, perf_measures, other_params) {
  
  ## Subset training set in batches
  batches_ind <- subsetBatch(condition=y_train, N_sample_size=other_params$N_sample_size, 
                            N_batch=other_params$N_batch)
  y_sgbatch_train <- map(1:other_params$N_batch, ~y_train[batches_ind[[.x]]])
  
  batch <- rep(0, ncol(train_expr))
  iwalk(batches_ind, ~{batch[.x] <<- .y})
  
  curr_train_expr <- train_expr[, do.call(c, batches_ind)]
  curr_y_train <- y_train[do.call(c, batches_ind)]
  batch <- batch[do.call(c, batches_ind)]  
  batches_ind <- map(1:other_params$N_batch, ~which(batch == .x))
  
  ## Remove genes with only 0 values in any batch
  g_keep <- map(1:other_params$N_batch, ~{
    which(apply(curr_train_expr[, batch == .x], 1, function(x){!all(x==0)}))
  })
  g_keep <- Reduce(intersect, g_keep)  
  curr_train_expr <- curr_train_expr[g_keep, ]
  curr_test_expr <- test_expr[g_keep, ]
  
  ## Simulate batch effect 
  sim_batch_res <- simBatch(dat=curr_train_expr, condition=curr_y_train, 
                           batches_ind=batches_ind, batch=batch, 
                           hyper_pars=other_params$hyper_pars)
  train_expr_batch <- sim_batch_res$new_dat
  
  ## Normalize datasets before training
  if(other_params$norm_data){
    train_expr_norm <- normalizeData(curr_train_expr)
    test_expr_norm <- normalizeData(curr_test_expr)
    train_expr_batch_whole_norm <- normalizeData(train_expr_batch)
    
    # Normalize within each batch
    train_expr_batch_norm <- matrix(NA, nrow=nrow(train_expr_batch), ncol=ncol(train_expr_batch), 
                                   dimnames=dimnames(train_expr_batch))
    iwalk(batches_ind, ~{
      train_expr_batch_norm[, .x] <<- normalizeData(train_expr_batch[, .x])
    })
  } else {
    train_expr_norm <- curr_train_expr
    test_expr_norm <- curr_test_expr
    train_expr_batch_whole_norm <- train_expr_batch_norm <- train_expr_batch
  }
  
  ## Process all learner types using map
  learner_results <- map(learner_types, ~{
    l_type <- .x
    learner_fit <- getPredFunctions(l_type)
    
    tryCatch({
      ## Prediction from original training to test, without batch effect
      pred_base_res <- trainPipe(train_set=train_expr_norm, train_label=curr_y_train, 
                                test_set=test_expr_norm, lfit=learner_fit)
      
      ## Prediction from training WITH batch effect to test
      pred_batch_res <- trainPipe(train_set=train_expr_batch_whole_norm, train_label=curr_y_train, 
                                 test_set=test_expr_norm, lfit=learner_fit)
      
      ## Prediction from training after ComBat adjustment (with reference batch)
      # Apply ComBat to combined training and test data
      ref_batch <- min(batch)
      combined_dat <- cbind(train_expr_batch, curr_test_expr)
      combined_batch <- c(batch, rep(ref_batch, ncol(curr_test_expr)))  # Test set uses reference batch
      # Create model matrix for training samples only (ComBat will handle test samples without labels)
      combined_labels <- c(curr_y_train, rep(0, ncol(curr_test_expr)))  # Use dummy labels for test
      train_expr_combat <- ComBat(combined_dat, batch=combined_batch, mod=model.matrix(~combined_labels), ref.batch=ref_batch)
      
      # Split back into training and test
      train_expr_combat_adj <- train_expr_combat[, 1:ncol(train_expr_batch)]
      test_expr_combat_adj <- train_expr_combat[, (ncol(train_expr_batch) + 1):ncol(train_expr_combat)]
      
      if(other_params$norm_data){
        train_expr_combat_norm <- normalizeData(train_expr_combat_adj)
        test_expr_combat_norm <- normalizeData(test_expr_combat_adj)
      } else {
        train_expr_combat_norm <- train_expr_combat_adj
        test_expr_combat_norm <- test_expr_combat_adj
      }
      pred_combat_res <- trainPipe(train_set=train_expr_combat_norm, train_label=curr_y_train, 
                                  test_set=test_expr_combat_norm, lfit=learner_fit)
      
      ## Prediction from training after MNN correction
      library(batchelor, quietly = TRUE)
      # Combine training and test data for MNN correction
      combined_dat <- cbind(train_expr_batch, curr_test_expr)
      combined_batch <- c(batch, rep(max(batch) + 1, ncol(curr_test_expr)))  # Test set gets new batch ID
      
      # Create batch list for MNN (each batch as separate matrix)
      unique_batches <- sort(unique(combined_batch))
      batch_list <- map(unique_batches, ~{
        combined_dat[, combined_batch == .x]
      })
      
      # Apply MNN correction with test set last in merge order
      mnn_result <- do.call(fastMNN, c(batch_list, list(merge.order = seq_along(unique_batches))))
      corrected_combined <- assay(mnn_result, "corrected")
      
      # Split back into training and test
      train_expr_mnn_adj <- corrected_combined[, 1:ncol(train_expr_batch)]
      test_expr_mnn_adj <- corrected_combined[, (ncol(train_expr_batch) + 1):ncol(corrected_combined)]
      
      if(other_params$norm_data){
        train_expr_mnn_norm <- normalizeData(train_expr_mnn_adj)
        test_expr_mnn_norm <- normalizeData(test_expr_mnn_adj)
      } else {
        train_expr_mnn_norm <- train_expr_mnn_adj
        test_expr_mnn_norm <- test_expr_mnn_adj
      }
      pred_mnn_res <- trainPipe(train_set=train_expr_mnn_norm, train_label=curr_y_train, 
                               test_set=test_expr_mnn_norm, lfit=learner_fit)
      
      ## Evaluate performance
      tst_scores <- list(
        NoBatch = pred_base_res$pred_tst_prob, 
        Batch = pred_batch_res$pred_tst_prob,
        ComBat = pred_combat_res$pred_tst_prob, 
        MNNcorrect = pred_mnn_res$pred_tst_prob
      )
      
      perf_df <- map(perf_measures, ~{
        perf_name <- .x
        as.data.frame(t(map_dbl(tst_scores, ~{
          preds <- .x
          if(perf_name=="mxe"){preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)}
          rocr_pred <- prediction(preds, as.numeric(as.character(y_test)))
          if(perf_name %in% c("acc", "f")){
            curr_perf <- performance(rocr_pred, perf_name)  
            return(curr_perf@y.values[[1]][which.min(abs(curr_perf@x.values[[1]]-0.5))])
          } else {
            curr_perf <- performance(rocr_pred, perf_name)
            return(as.numeric(curr_perf@y.values))
          }
        })))
      })
      names(perf_df) <- perf_measures
      perf_df <- do.call(rbind, perf_df)
      
      list(learner = l_type, perf_df = perf_df, success = TRUE, error = NULL)
      
    }, error = function(e) {
      # PRESERVED: Detailed error reporting for batch correction failures
      if(grepl("MNN", e$message)) {
        cat("ERROR in MNN correction:", e$message, "\n")
        cat("Batch info - unique batches:", unique(batch), "\n")
        cat("Samples per batch:", table(batch), "\n")
        cat("Data summary:\n")
        print(summary(as.vector(train_expr_batch)))
      } else if(grepl("ComBat", e$message)) {
        cat("ERROR in ComBat correction:", e$message, "\n")
      }
      
      list(learner = l_type, perf_df = NULL, success = FALSE, error = e$message)
    })
  })
  
  names(learner_results) <- learner_types
  
  list(iteration = iteration_id, learner_results = learner_results)
}

#' Write results for a single learner and performance measure
write_learner_results <- function(learner_results, perf_measure, exp_name) {
  
  # Extract successful results for this performance measure
  successful_results <- keep(learner_results, "success")
  
  if(length(successful_results) == 0) return(NULL)
  
  # Process each learner type
  walk(names(successful_results), ~{
    l_type <- .x
    learner_data <- successful_results[[l_type]]
    
    if(!is.null(learner_data$perf_df) && perf_measure %in% rownames(learner_data$perf_df)) {
      first_file <- !file.exists(sprintf('results/%s_%s_%s.csv', l_type, perf_measure, exp_name))
      write.table(learner_data$perf_df[perf_measure, ], 
                 sprintf('results/%s_%s_%s.csv', l_type, perf_measure, exp_name),
                 append=!first_file, col.names=first_file, row.names=FALSE, sep=",")
    }
  })
}

####  Run Pipeline with Preserved Retry Logic  ####

# Create parameter list for passing to iteration function
other_params <- list(
  N_batch = N_batch,
  N_sample_size = N_sample_size,
  hyper_pars = hyper_pars,
  norm_data = norm_data,

)

cat("Starting simulation pipeline...\n")
start_time <- Sys.time()

# PRESERVED: Original retry logic for failed iterations
# This is critical for maintaining exact statistical behavior
ID <- 1
all_results <- list()
while(ID <= iterations) {
  if(ID %% 10 == 0) cat(sprintf("Processing iteration %d/%d\n", ID, iterations))
  
  result <- process_single_iteration(ID, train_expr, test_expr, y_train, y_test,
                                   learner_types, perf_measures, other_params)
  
  # Check if any learner failed - if so, retry this iteration
  any_failed <- any(map_lgl(result$learner_results, ~!.x$success))
  
  if(any_failed) {
    # PRESERVED: Decrement ID to retry failed iteration (original: ID <- ID - 1; break; next)
    cat(sprintf("Iteration %d failed, retrying...\n", ID))
    ID <- ID - 1
  } else {
    all_results[[ID]] <- result
  }
  
  ID <- ID + 1
}

# Extract and write results using functional programming
cat("Writing results...\n")

# Group results by iteration and learner
iteration_learner_results <- map(all_results, "learner_results")

# Write results for each performance measure
walk(perf_measures, ~{
  perf_measure <- .x
  cat(sprintf("Writing results for %s...\n", perf_measure))
  
  walk(iteration_learner_results, ~{
    write_learner_results(.x, perf_measure, exp_name)
  })
})

end_time <- Sys.time()
duration <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)

# Summary statistics
successful_iterations <- map_int(iteration_learner_results, ~{
  sum(map_lgl(.x, "success"))
})

total_successful <- sum(successful_iterations)
total_expected <- iterations * length(learner_types)

cat(sprintf("\nSimulation pipeline completed in %.2f minutes\n", duration))
cat(sprintf("Successful model fits: %d/%d (%.1f%%)\n", 
           total_successful, total_expected, 100 * total_successful / total_expected))

# Error summary
all_errors <- map(iteration_learner_results, ~{
  failed_results <- keep(.x, ~!.x$success)
  if(length(failed_results) > 0) {
    map_chr(failed_results, "error")
  } else {
    NULL
  }
}) %>% compact() %>% unlist()

if(length(all_errors) > 0) {
  cat("\nError summary:\n")
  error_counts <- table(all_errors)
  iwalk(error_counts, ~{
    cat(sprintf("  %s: %d occurrences\n", .y, .x))
  })
}