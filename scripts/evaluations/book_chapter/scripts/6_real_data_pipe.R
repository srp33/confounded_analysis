# 6-Study Real Data Analysis for Robustifying Evaluation
# Refactored with helper functions to reduce code duplication
# FIXED: All variables explicitly passed as parameters, no global scope dependencies

options(warn = -1)
suppressMessages(suppressWarnings({
  rm(list=ls())
  if(!dir.exists("/scripts/evaluations/robustifying/results_real_6studies")){
    dir.create("/scripts/evaluations/robustifying/results_real_6studies")
  }
  sapply(c("glmnet", "SummarizedExperiment", "sva", "DESeq2", "ROCR", "ggplot2", 
           "gridExtra", "reshape2", "dplyr", "nnls"), require, character.only=TRUE, quietly=TRUE)
}))

load("/scripts/evaluations/robustifying/data/TB_real_data.RData")
source("/scripts/evaluations/robustifying/code/helper.R")
source("/scripts/adjust/gmm_adjust.R")
set.seed(123)

####  Parameters  ####
command_args <- commandArgs(trailingOnly=TRUE)
debug_mode <- length(command_args) > 0 && command_args[1] == "debug"

norm_data <- TRUE
use_ref_combat <- FALSE
n_highvar_genes <- 1000
B <- if(debug_mode) 3 else 100

if(debug_mode) {
  cat("=== DEBUG MODE ENABLED ===\n")
  cat("Running with", B, "iterations instead of 100\n")
  cat("==========================\n\n")
}

learner_types <- c("lasso", "rf", "svm")
perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc")
perf_measures_names <- c("Mean cross-entropy loss", "AUC", "Root-mean-squared error", 
                         "F1 score", "Error rate", "Accuracy")
names(perf_measures_names) <- perf_measures

study_names <- names(dat_lst)
cat("Running 6-study analysis with studies:", paste(study_names, collapse=", "), "\n")

####  Helper Functions  ####

#' Prepare training and test data
prepare_datasets <- function(dat_lst, label_lst, test_name, study_names) {
  train_name <- setdiff(study_names, test_name)
  
  dat <- do.call(cbind, dat_lst[train_name])
  batch <- rep(1:length(train_name), times=sapply(dat_lst[train_name], ncol))
  batches_ind <- lapply(1:length(train_name), function(i){which(batch==i)})
  batch_names <- levels(factor(batch))
  group <- do.call(c, label_lst[train_name])
  y_sgbatch_train <- lapply(batch_names, function(k){group[batch==k]})
  
  dat_test <- dat_lst[[test_name]]
  group_test <- label_lst[[test_name]]
  
  list(dat=dat, batch=batch, batches_ind=batches_ind, batch_names=batch_names, 
       group=group, y_sgbatch_train=y_sgbatch_train, 
       dat_test=dat_test, group_test=group_test)
}

#' Select highly variable genes and reduce feature space
reduce_features <- function(dat, dat_test, n_genes=1000) {
  genes_sel_names <- order(rowVars(dat), decreasing=TRUE)[1:n_genes]
  list(dat=dat[genes_sel_names, ], 
       dat_test=dat_test[genes_sel_names, ])
}

#' Apply batch correction methods
apply_batch_corrections <- function(dat, batch, group) {
  dat_combat <- ComBat(dat, batch=batch, mod=model.matrix(~group))
  
  dat_gmm_adj <- gmm_adjust(data=t(dat), batch=batch, 
                            mean_mean_zero=TRUE, unit_var=TRUE, debug=FALSE)
  dat_gmm_adj <- t(dat_gmm_adj)
  
  list(dat_combat=dat_combat, dat_gmm_adj=dat_gmm_adj)
}

#' Normalize data matrices within batches
normalize_within_batches <- function(dat, batch, batch_names) {
  dat_whole_norm <- normalizeData(dat)
  dat_batch_norm <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat), dimnames=dimnames(dat))
  for(k in batch_names) {
    dat_batch_norm[, batch==k] <- normalizeData(dat[, batch==k])
  }
  list(whole_norm=dat_whole_norm, batch_norm=dat_batch_norm)
}

#' Train model for single learner type
train_single_learner <- function(l_type, dat_batch_whole_norm, dat_combat_whole_norm, 
                                 dat_gmm_whole_norm, dat_batch_norm, group, batch, batch_names,
                                 use_ref_combat) {
  learner_fit <- getPredFunctions(l_type)
  
  pred_unadj_res <- trainPipe(train_set=dat_batch_whole_norm, train_label=group, 
                              test_set=NULL, lfit=learner_fit, use_ref_combat=use_ref_combat)
  
  pred_combat_res <- trainPipe(train_set=dat_combat_whole_norm, train_label=group, 
                               test_set=NULL, lfit=learner_fit, use_ref_combat=use_ref_combat)
  
  pred_gmm_res <- trainPipe(train_set=dat_gmm_whole_norm, train_label=group, 
                            test_set=NULL, lfit=learner_fit, use_ref_combat=use_ref_combat)
  
  pred_sgbatch_res <- lapply(batch_names, function(k) {
    trainPipe(train_set=dat_batch_norm[, batch==k], train_label=group[batch==k], 
              test_set=NULL, lfit=learner_fit, use_ref_combat=use_ref_combat)
  })
  names(pred_sgbatch_res) <- paste0("Batch", batch_names)
  sgbatch_mod <- lapply(pred_sgbatch_res, function(res) res$mod)
  
  result <- list(unadj_mod=pred_unadj_res$mod, combat_mod=pred_combat_res$mod, 
                 gmm_mod=pred_gmm_res$mod, sgbatch_mod=sgbatch_mod)
  
  # Explicit cleanup of large intermediate objects to match original behavior
  rm(pred_unadj_res, pred_combat_res, pred_gmm_res, pred_sgbatch_res, learner_fit, envir=environment())
  
  result
}

#' Train models for all learner types
train_all_learners <- function(learner_types, dat_batch_whole_norm, dat_combat_whole_norm, 
                               dat_gmm_whole_norm, dat_batch_norm, group, batch, batch_names,
                               use_ref_combat) {
  unadj_mod_lst <- combat_mod_lst <- gmm_mod_lst <- sgbatch_mod_lst <- list()
  
  for(l_type in learner_types) {
    trained <- train_single_learner(l_type, dat_batch_whole_norm, dat_combat_whole_norm,
                                    dat_gmm_whole_norm, dat_batch_norm, group, batch, batch_names,
                                    use_ref_combat)
    unadj_mod_lst[[l_type]] <- trained$unadj_mod
    combat_mod_lst[[l_type]] <- trained$combat_mod
    gmm_mod_lst[[l_type]] <- trained$gmm_mod
    sgbatch_mod_lst[[l_type]] <- trained$sgbatch_mod
  }
  
  list(unadj_mod_lst=unadj_mod_lst, combat_mod_lst=combat_mod_lst, 
       gmm_mod_lst=gmm_mod_lst, sgbatch_mod_lst=sgbatch_mod_lst)
}

#' Calculate single-learner ensemble weights for one learner type
calculate_weights_single_learner <- function(l_type, train_lst, y_sgbatch_train, batch, 
                                             batch_names, use_ref_combat) {
  learner_fit <- getPredFunctions(l_type)
  
  cs_zmat <- CS_zmatrix(study_lst=train_lst, label_lst=y_sgbatch_train,
                        lfit=learner_fit, perf_name="mxe", 
                        use_ref_combat=use_ref_combat)
  cs_weights <- CS_weight(cs_zmat)
  
  reg_ssl_res <- Reg_SSL_pred(study_lst=train_lst, label_lst=y_sgbatch_train,
                              lfit=learner_fit, use_ref_combat=use_ref_combat)
  reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res$coef),
                             n_seq=table(batch)[batch_names])
  
  stacked_pred <- do.call(rbind, reg_ssl_res$pred)
  stacked_label <- do.call(c, lapply(y_sgbatch_train, as.character))
  reg_s_beta <- nnls(A=stacked_pred, b=as.numeric(stacked_label))$x
  reg_s_beta <- reg_s_beta / sum(reg_s_beta)
  
  list(cs_zmat=cs_zmat, cs_weights=cs_weights, reg_ssl_res=reg_ssl_res, 
       reg_a_beta=reg_a_beta, reg_s_beta=reg_s_beta)
}

#' Calculate single-learner ensemble weights for all learner types
calculate_single_learner_weights <- function(learner_types, train_lst, y_sgbatch_train, 
                                             batch, batch_names, use_ref_combat) {
  cs_zmat_lst <- cs_weights_seq <- reg_ssl_res <- reg_a_beta <- reg_s_beta <- list()
  
  for(l_type in learner_types) {
    weights <- calculate_weights_single_learner(l_type, train_lst, y_sgbatch_train, 
                                                batch, batch_names, use_ref_combat)
    cs_zmat_lst[[l_type]] <- weights$cs_zmat
    cs_weights_seq[[l_type]] <- weights$cs_weights
    reg_ssl_res[[l_type]] <- weights$reg_ssl_res
    reg_a_beta[[l_type]] <- weights$reg_a_beta
    reg_s_beta[[l_type]] <- weights$reg_s_beta
  }
  
  list(cs_zmat_lst=cs_zmat_lst, cs_weights_seq=cs_weights_seq, 
       reg_ssl_res=reg_ssl_res, reg_a_beta=reg_a_beta, reg_s_beta=reg_s_beta)
}

#' Calculate cross-learner ensemble weights
calculate_cross_learner_weights <- function(learner_types, train_lst, y_sgbatch_train, 
                                            batch, batches_ind, cs_zmat_lst, reg_ssl_res,
                                            use_ref_combat) {
  cm_navg_weights <- rep((as.matrix(sapply(batches_ind, length)) / sum(sapply(batches_ind, length))),
                         length(learner_types))
  cm_navg_weights <- cm_navg_weights / sum(cm_navg_weights)
  
  cm_cs_weights_seq <- CS_weight_crossmod(cs_zmat_lst)
  
  cm_reg_ssl_res <- crossmod_Reg_SSL_pred(study_lst=train_lst, label_lst=y_sgbatch_train,
                                          learner_lst=learner_types, use_ref_combat=use_ref_combat)
  cm_reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, cm_reg_ssl_res$coef),
                                n_seq=sapply(batches_ind, length))
  
  cm_stacked_pred <- do.call(rbind, cm_reg_ssl_res$pred)
  cm_stacked_label <- do.call(c, lapply(y_sgbatch_train, as.character))
  cm_reg_s_beta <- nnls(A=cm_stacked_pred, b=as.numeric(cm_stacked_label))$x
  cm_reg_s_beta <- cm_reg_s_beta / sum(cm_reg_s_beta)
  
  list(cm_navg_weights=cm_navg_weights, cm_cs_weights_seq=cm_cs_weights_seq,
       cm_reg_ssl_res=cm_reg_ssl_res, cm_reg_a_beta=cm_reg_a_beta, 
       cm_reg_s_beta=cm_reg_s_beta)
}

#' Generate predictions for single learner type in bootstrap iteration
generate_predictions_single_learner <- function(l_type, unadj_mod_lst, combat_mod_lst, 
                                                gmm_mod_lst, sgbatch_mod_lst, dat_test_norm,
                                                batch_names, navg_weights, cs_weights_seq, 
                                                reg_a_beta, reg_s_beta, perf_measures, 
                                                group_test, use_ref_combat) {
  unadj_tst_prob <- predWrapper(unadj_mod_lst[[l_type]], dat_test_norm, l_type)
  combat_tst_prob <- predWrapper(combat_mod_lst[[l_type]], dat_test_norm, l_type)
  gmm_tst_prob <- predWrapper(gmm_mod_lst[[l_type]], dat_test_norm, l_type)
  onestep_res <- ensemble_wrapper_realdata(sgbatch_mod_lst, l_type, dat_test_norm,
                                           navg_weights, cs_weights_seq, reg_a_beta, reg_s_beta)
  
  tst_scores <- c(list(Batch=unadj_tst_prob), onestep_res$pred_test_lst,
                  list(ComBat=combat_tst_prob, GMM=gmm_tst_prob,
                       Avg=onestep_res$pred_avg, n_Avg=onestep_res$pred_N_avg,
                       CS_Avg=onestep_res$pred_cs_avg, Reg_a=onestep_res$pred_reg_a,
                       Reg_s=onestep_res$pred_reg_s))
  perf_df <- perf_wrapper(perf_measures, tst_scores, group_test)
  
  result <- list(perf_df=perf_df, tst_scores=tst_scores, unadj_tst_prob=unadj_tst_prob,
                 combat_tst_prob=combat_tst_prob, gmm_tst_prob=gmm_tst_prob)
  
  # Explicit cleanup of intermediate objects
  rm(onestep_res, envir=environment())
  
  result
}

#' Generate predictions for all learner types in bootstrap iteration
generate_predictions <- function(learner_types, unadj_mod_lst, combat_mod_lst, 
                                gmm_mod_lst, sgbatch_mod_lst, dat_test_norm, batch_names,
                                navg_weights, cs_weights_seq, reg_a_beta, reg_s_beta,
                                perf_measures, group_test, use_ref_combat) {
  perf_df_lst <- tst_scores_modlst <- list()
  unadj_tst_prob <- combat_tst_prob <- gmm_tst_prob <- NULL
  
  for(l_type in learner_types) {
    result <- generate_predictions_single_learner(l_type, unadj_mod_lst, combat_mod_lst,
                                                  gmm_mod_lst, sgbatch_mod_lst, dat_test_norm,
                                                  batch_names, navg_weights, cs_weights_seq,
                                                  reg_a_beta, reg_s_beta, perf_measures,
                                                  group_test, use_ref_combat)
    perf_df_lst[[l_type]] <- result$perf_df
    tst_scores_modlst[[l_type]] <- result$tst_scores
    unadj_tst_prob <- result$unadj_tst_prob
    combat_tst_prob <- result$combat_tst_prob
    gmm_tst_prob <- result$gmm_tst_prob
  }
  
  list(perf_df_lst=perf_df_lst, tst_scores_modlst=tst_scores_modlst,
       unadj_tst_prob=unadj_tst_prob, combat_tst_prob=combat_tst_prob, 
       gmm_tst_prob=gmm_tst_prob)
}

#' Generate cross-learner predictions for bootstrap iteration
generate_crossmod_predictions <- function(tst_scores_modlst, batch_names, batches_ind,
                                          cm_navg_weights, cm_cs_weights_seq, cm_reg_a_beta,
                                          cm_reg_s_beta, unadj_tst_prob, combat_tst_prob,
                                          gmm_tst_prob, perf_measures, group_test) {
  preds_crossmod <- lapply(tst_scores_modlst, function(x) {
    do.call(cbind, x[paste0("Batch", batch_names)])
  })
  
  cm_onestep_res <- ensemble_crossmod_wrapperNew(preds_crossmod, length(batches_ind),
                                                 cm_navg_weights, cm_cs_weights_seq,
                                                 cm_reg_a_beta, cm_reg_s_beta)
  
  tst_cm_scores <- list(Avg=cm_onestep_res$cm_avg, n_Avg=cm_onestep_res$cm_N_avg,
                        CS_Avg=cm_onestep_res$cm_cs_avg, Reg_a=cm_onestep_res$cm_reg_a,
                        Reg_s=cm_onestep_res$cm_reg_s)
  perf_crossmod_df <- perf_wrapper(perf_measures, tst_cm_scores, group_test)
  
  tst_cm_scores_full <- c(list(Batch=unadj_tst_prob, ComBat=combat_tst_prob, GMM=gmm_tst_prob),
                          tst_cm_scores)
  perf_crossmod_full_df <- perf_wrapper(perf_measures, tst_cm_scores_full, group_test)
  
  list(perf_crossmod_full_df=perf_crossmod_full_df, tst_cm_scores_full=tst_cm_scores_full)
}

#' Extract performance metrics by model and measure
extract_performance_data <- function(perf_df_lst, perf_measure, n_batches) {
  perf_data <- list()
  
  for(model_name in names(perf_df_lst)) {
    perf_res <- perf_df_lst[[model_name]]
    if(perf_measure %in% rownames(perf_res)) {
      if(model_name == "crossmod") {
        perf_data[[model_name]] <- perf_res[perf_measure, , drop=FALSE]
      } else {
        start_col <- 2 + n_batches
        if(ncol(perf_res) >= start_col) {
          selected_cols <- c(1, start_col:ncol(perf_res))
          perf_data[[model_name]] <- perf_res[perf_measure, selected_cols]
        } else {
          perf_data[[model_name]] <- perf_res[perf_measure, 1, drop=FALSE]
        }
      }
    }
  }
  
  perf_data
}

#' Write performance results to CSV
write_performance_results <- function(perf_data, test_name, perf_measure, iteration) {
  if(length(perf_data) == 0) return(NULL)
  
  summary_df <- reshape2::melt(perf_data)
  summary_df$iteration <- iteration
  colnames(summary_df) <- c("Method", "value", "Model", "Iteration")
  
  filepath <- sprintf('/scripts/evaluations/robustifying/results_real_6studies/test%s_%s.csv',
                      test_name, perf_measure)
  first_file <- !file.exists(filepath)
  write.table(summary_df, filepath, append=!first_file, col.names=first_file,
              row.names=FALSE, sep=",")
  
  summary_df
}

####  Main Pipeline  ####
for(s in study_names) {
  ## Prepare datasets
  datasets <- prepare_datasets(dat_lst, label_lst, s, study_names)
  test_name <- s
  train_name <- setdiff(study_names, test_name)
  
  ## Feature reduction
  feat_reduced <- reduce_features(datasets$dat, datasets$dat_test, n_highvar_genes)
  dat <- feat_reduced$dat
  dat_test <- feat_reduced$dat_test
  
  ## Batch correction
  batch_corr <- apply_batch_corrections(dat, datasets$batch, datasets$group)
  dat_combat <- batch_corr$dat_combat
  dat_gmm_adj <- batch_corr$dat_gmm_adj
  
  ## Normalize data
  if(norm_data) {
    norm_raw <- normalize_within_batches(dat, datasets$batch, datasets$batch_names)
    norm_combat <- normalize_within_batches(dat_combat, datasets$batch, datasets$batch_names)
    norm_gmm <- normalize_within_batches(dat_gmm_adj, datasets$batch, datasets$batch_names)
  } else {
    norm_raw <- list(whole_norm=dat, batch_norm=dat)
    norm_combat <- list(whole_norm=dat_combat, batch_norm=dat_combat)
    norm_gmm <- list(whole_norm=dat_gmm_adj, batch_norm=dat_gmm_adj)
  }
  
  dat_batch_whole_norm <- norm_raw$whole_norm
  dat_batch_norm <- norm_raw$batch_norm
  dat_combat_whole_norm <- norm_combat$whole_norm
  dat_combat_norm <- norm_combat$batch_norm
  dat_gmm_whole_norm <- norm_gmm$whole_norm
  dat_gmm_norm <- norm_gmm$batch_norm
  
  train_lst <- lapply(datasets$batch_names, function(k) dat_batch_norm[, datasets$batch==k])
  
  ## Train all models
  trained_models <- train_all_learners(learner_types, dat_batch_whole_norm, 
                                       dat_combat_whole_norm, dat_gmm_whole_norm,
                                       dat_batch_norm, datasets$group, datasets$batch,
                                       datasets$batch_names, use_ref_combat)
  
  ## Calculate ensemble weights
  single_weights <- calculate_single_learner_weights(learner_types, train_lst, 
                                                     datasets$y_sgbatch_train, datasets$batch,
                                                     datasets$batch_names, use_ref_combat)
  
  navg_weights <- as.matrix(table(datasets$batch)[datasets$batch_names]) / sum(table(datasets$batch))
  
  cross_weights <- calculate_cross_learner_weights(learner_types, train_lst, 
                                                   datasets$y_sgbatch_train, datasets$batch,
                                                   datasets$batches_ind, single_weights$cs_zmat_lst,
                                                   single_weights$reg_ssl_res, use_ref_combat)
  
  ## Save weights
  save(navg_weights, 
       cs_zmat_lst=single_weights$cs_zmat_lst,
       cs_weights_seq=single_weights$cs_weights_seq,
       reg_ssl_res=single_weights$reg_ssl_res,
       reg_a_beta=single_weights$reg_a_beta,
       reg_s_beta=single_weights$reg_s_beta,
       cm_navg_weights=cross_weights$cm_navg_weights,
       cm_cs_weights_seq=cross_weights$cm_cs_weights_seq,
       cm_reg_ssl_res=cross_weights$cm_reg_ssl_res,
       cm_reg_a_beta=cross_weights$cm_reg_a_beta,
       cm_reg_s_beta=cross_weights$cm_reg_s_beta,
       file=sprintf('/scripts/evaluations/robustifying/results_real_6studies/test%s_weights.RData', test_name))
  
  ## Prediction & Ensemble loop
  dat_testOri <- datasets$dat_test
  group_testOri <- datasets$group_test
  
  b <- 1
  while(b <= B) {
    boot_ind <- sample(1:ncol(dat_testOri), ncol(dat_testOri), replace=TRUE)
    dat_test <- dat_testOri[, boot_ind]
    group_test <- group_testOri[boot_ind]
    
    if(norm_data) {
      dat_test_norm <- normalizeData(dat_test)
    } else {
      dat_test_norm <- dat_test
    }
    
    if(any(is.na(dat_test_norm))) {
      b <- b - 1
    } else {
      ## Generate single-learner predictions
      single_preds <- generate_predictions(learner_types, trained_models$unadj_mod_lst,
                                          trained_models$combat_mod_lst, 
                                          trained_models$gmm_mod_lst,
                                          trained_models$sgbatch_mod_lst, dat_test_norm,
                                          datasets$batch_names, navg_weights,
                                          single_weights$cs_weights_seq,
                                          single_weights$reg_a_beta,
                                          single_weights$reg_s_beta, perf_measures,
                                          group_test, use_ref_combat)
      
      perf_df_lst <- single_preds$perf_df_lst
      tst_scores_modlst <- single_preds$tst_scores_modlst
      
      ## Generate cross-learner predictions
      cross_preds <- generate_crossmod_predictions(tst_scores_modlst, datasets$batch_names,
                                                   datasets$batches_ind,
                                                   cross_weights$cm_navg_weights,
                                                   cross_weights$cm_cs_weights_seq,
                                                   cross_weights$cm_reg_a_beta,
                                                   cross_weights$cm_reg_s_beta,
                                                   single_preds$unadj_tst_prob,
                                                   single_preds$combat_tst_prob,
                                                   single_preds$gmm_tst_prob,
                                                   perf_measures, group_test)
      
      perf_df_lst[["crossmod"]] <- cross_preds$perf_crossmod_full_df
      tst_scores_modlst[["crossmod"]] <- cross_preds$tst_cm_scores_full
      
      ## Write performance results
      n_batches <- length(unique(datasets$batch))
      for(i in 1:length(perf_measures)) {
        perf_data <- extract_performance_data(perf_df_lst, perf_measures[i], n_batches)
        write_performance_results(perf_data, test_name, perf_measures[i], b)
      }
      
      # Explicit cleanup of large intermediate objects to match original behavior
      rm(boot_ind, dat_test, group_test, perf_df_lst, tst_scores_modlst, 
         single_preds, cross_preds, envir=environment())
    }
    
    b <- b + 1
  }
}

cat("6-study analysis completed successfully\n")