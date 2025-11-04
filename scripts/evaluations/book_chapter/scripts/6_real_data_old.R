# 6-Study Real Data Analysis for Robustifying Evaluation
# This script runs the same analysis as 3_real_data_pipe.R but uses all 6 studies
# instead of the subset of 4 studies

# Suppress all output and warnings
options(warn = -1)
suppressMessages(suppressWarnings({

rm(list=ls())
if(!dir.exists("/scripts/evaluations/robustifying/results_real_6studies")){dir.create("/scripts/evaluations/robustifying/results_real_6studies")}
sapply(c("glmnet", "SummarizedExperiment", "sva", "DESeq2", "ROCR", "ggplot2", 
         "gridExtra", "reshape2", "dplyr", "nnls"), require, character.only=TRUE, quietly=TRUE)

}))
load("/scripts/evaluations/robustifying/data/TB_real_data.RData")
source("/scripts/evaluations/robustifying/code/helper.R")
source("/scripts/adjust/gmm_adjust.R")
set.seed(123)

####  Parameters  ####
# Check for debug mode command line argument
command_args <- commandArgs(trailingOnly=TRUE)
debug_mode <- length(command_args) > 0 && command_args[1] == "debug"

norm_data <- TRUE  # whether to normalize datasets by features using z-score scaling
use_ref_combat <- FALSE  # whether to use ref combat to adjust test set against training set

n_highvar_genes <- 1000  # number of highly variable genes to use in feature reduction
B <- if(debug_mode) 3 else 100  # bootstrap samples (3 for debug, 100 for full run)

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

# Use ALL 6 studies instead of subset of 4
# All studies: "GSE37250_SA", "GSE37250_M", "GSE39941_M", "US", "Africa", "India"
study_names <- names(dat_lst)
cat("Running 6-study analysis with studies:", paste(study_names, collapse=", "), "\n")

####  Start pipeline  ####
#s = "India"
for(s in study_names){
  ## Get training & test set
  test_name <- s
  train_name <- setdiff(study_names, test_name)
  
  # training set
  dat <- do.call(cbind, dat_lst[train_name])
  batch <- rep(1:length(train_name), times=sapply(dat_lst[train_name], ncol))
  batches_ind <- lapply(1:length(train_name), function(i){which(batch==i)})
  batch_names <- levels(factor(batch))
  group <- do.call(c, label_lst[train_name])
  y_sgbatch_train <- lapply(batch_names, function(k){group[batch==k]})
  
  # test
  dat_test <- dat_lst[[test_name]]
  group_test <- label_lst[[test_name]]
  
  
  ####  Preprocess
  ## feature reduction - select highly variable genes in training data
  genes_sel_names <- order(rowVars(dat), decreasing=TRUE)[1:n_highvar_genes]
  dat <- dat[genes_sel_names, ]
  dat_test <- dat_test[genes_sel_names, ]
  
  ## batch correction again after feature
  dat_combat <- ComBat(dat, batch=batch, mod=model.matrix(~group))
  
  ## Apply gmm_adjust batch correction
  dat_gmm_adj <- gmm_adjust(data=t(dat), batch=batch, 
                           mean_mean_zero=TRUE, unit_var=TRUE, debug=FALSE)
  
  dat_gmm_adj <- t(dat_gmm_adj)
  
  ## normalize features
  if(norm_data){
    dat_batch_whole_norm <- normalizeData(dat)  # norm training set as a whole
    # norm training set within each batch
    dat_batch_norm <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat), dimnames=dimnames(dat))  
    for(k in batch_names){dat_batch_norm[, batch==k] <- normalizeData(dat[, batch==k])}
    
    dat_combat_whole_norm <- normalizeData(dat_combat)  # norm combat adjusted data as a whole
    # norm combat adjusted data within each batch
    dat_combat_norm <- matrix(NA, nrow=nrow(dat_combat), ncol=ncol(dat_combat), dimnames=dimnames(dat_combat))  
    for(k in batch_names){dat_combat_norm[, batch==k] <- normalizeData(dat_combat[, batch==k])}
    
    dat_gmm_whole_norm <- normalizeData(dat_gmm_adj)  # norm gmm adjusted data as a whole
    # norm gmm adjusted data within each batch
    dat_gmm_norm <- matrix(NA, nrow=nrow(dat_gmm_adj), ncol=ncol(dat_gmm_adj), dimnames=dimnames(dat_gmm_adj))  
    for(k in batch_names){dat_gmm_norm[, batch==k] <- normalizeData(dat_gmm_adj[, batch==k])}
  }else{
    dat_batch_whole_norm <- dat_batch_norm <- dat 
    dat_combat_whole_norm <- dat_combat_norm <- dat_combat
    dat_gmm_whole_norm <- dat_gmm_norm <- dat_gmm_adj
  }
  train_lst <- lapply(batch_names, function(k){dat_batch_norm[, batch==k]})
  
  
  
  ####  Training 
  #l_type="lasso"
  unadj_mod_lst <- combat_mod_lst <- gmm_mod_lst <- sgbatch_mod_lst <- sgbatch_mod_comb <- list()
  cs_zmat_lst <- cs_weights_seq <- reg_ssl_res <- reg_a_beta <- reg_s_beta <- list()
  for(l_type in learner_types){
    learner_fit <- getPredFunctions(l_type)
    
    ##  Training on original train set
    pred_unadj_res <- trainPipe(train_set=dat_batch_whole_norm, train_label=group, test_set=NULL, 
                                lfit=learner_fit, use_ref_combat=use_ref_combat)
    unadj_mod_lst[[l_type]] <- pred_unadj_res$mod
    
    ##  Training on train set after batch adjustment (Merged)
    pred_combat_res <- trainPipe(train_set=dat_combat_whole_norm, train_label=group, test_set=NULL, 
                                 lfit=learner_fit, use_ref_combat=use_ref_combat)
    combat_mod_lst[[l_type]] <- pred_combat_res$mod
    
    ##  Training on train set after GMM adjustment
    pred_gmm_res <- trainPipe(train_set=dat_gmm_whole_norm, train_label=group, test_set=NULL, 
                             lfit=learner_fit, use_ref_combat=use_ref_combat)
    gmm_mod_lst[[l_type]] <- pred_gmm_res$mod
    
    ##  Training within each batch from train set
    pred_sgbatch_res <- lapply(batch_names, function(k){
      trainPipe(train_set=dat_batch_norm[, batch==k], train_label=group[batch==k], test_set=NULL, 
                lfit=learner_fit, use_ref_combat=use_ref_combat)
    })
    names(pred_sgbatch_res) <- paste0("Batch", batch_names)
    sgbatch_mod_lst[[l_type]] <- lapply(pred_sgbatch_res, function(res){res$mod})
    
    
    ## Ensemble weights - single learner
    # cs
    cs_zmat_lst[[l_type]] <- CS_zmatrix(study_lst=train_lst, label_lst=y_sgbatch_train, 
                                        lfit=learner_fit, perf_name="mxe", 
                                        use_ref_combat=use_ref_combat)
    cs_weights_seq[[l_type]] <- CS_weight(cs_zmat_lst[[l_type]])
    # reg-a
    reg_ssl_res[[l_type]] <- Reg_SSL_pred(study_lst=train_lst, label_lst=y_sgbatch_train, 
                                          lfit=learner_fit, use_ref_combat=use_ref_combat)
    reg_a_beta[[l_type]] <- Reg_a_weight(coef_mat=do.call(rbind, reg_ssl_res[[l_type]]$coef), 
                                         n_seq=table(batch)[batch_names])
    # reg-s
    stacked_pred <- do.call(rbind, reg_ssl_res[[l_type]]$pred)
    stacked_label <- do.call(c, lapply(y_sgbatch_train, as.character))
    reg_s_beta[[l_type]] <- nnls(A=stacked_pred, b=as.numeric(stacked_label))$x
    reg_s_beta[[l_type]] <- reg_s_beta[[l_type]] / sum(reg_s_beta[[l_type]])
  }
  
  ## Ensemble weights - across learners
  navg_weights <- (as.matrix(table(batch)[batch_names]) / sum(table(batch)))  # single learner
  cm_navg_weights <- rep((as.matrix(sapply(batches_ind, length)) / sum(sapply(batches_ind, length))), 
                         length(learner_types))
  cm_navg_weights <- cm_navg_weights / sum(cm_navg_weights)  # across learners
  
  cm_cs_weights_seq <- CS_weight_crossmod(cs_zmat_lst)
  
  cm_reg_ssl_res <- crossmod_Reg_SSL_pred(study_lst=train_lst, label_lst=y_sgbatch_train,
                                          learner_lst=learner_types, use_ref_combat=use_ref_combat)
  cm_reg_a_beta <- Reg_a_weight(coef_mat=do.call(rbind, cm_reg_ssl_res$coef), 
                                n_seq=sapply(batches_ind, length))
  
  cm_stacked_pred <- do.call(rbind, cm_reg_ssl_res$pred)
  cm_stacked_label <- do.call(c, lapply(y_sgbatch_train, as.character))
  cm_reg_s_beta <- nnls(A=cm_stacked_pred, b=as.numeric(cm_stacked_label))$x
  cm_reg_s_beta <- cm_reg_s_beta / sum(cm_reg_s_beta)
  
  save(navg_weights, cs_zmat_lst, cs_weights_seq, reg_ssl_res, reg_a_beta, reg_s_beta, 
       cm_navg_weights, cm_cs_weights_seq, cm_reg_ssl_res, cm_reg_a_beta, cm_reg_s_beta,
       file=sprintf('/scripts/evaluations/robustifying/results_real_6studies/test%s_weights.RData', test_name))
  rm(pred_unadj_res, pred_combat_res, pred_gmm_res, pred_sgbatch_res)
  
  
  ####  Prediction & Ensemble
  b=1
  perf_df_lst <- tst_scores_modlst <- list()
  # save original test set
  dat_testOri <- dat_test; group_testOri <- group_test
  rm(dat_test, group_test)  # save original data for bootstrap
  
  while(b<=B){
    
    ## draw bootstrap sample
    boot_ind <- sample(1:ncol(dat_testOri), ncol(dat_testOri), replace=TRUE)
    dat_test <- dat_testOri[, boot_ind]
    group_test <- group_testOri[boot_ind]
    
    ## normalize
    if(norm_data){dat_test_norm <- normalizeData(dat_test)  # norm test set
    }else{dat_test_norm <- dat_test}
    if(any(is.na(dat_test_norm))){
      b <- b - 1
    }else{
      #l_type="lasso"
      for(l_type in learner_types){
        ##  Obtain prediction on each test set
        unadj_tst_prob <- predWrapper(unadj_mod_lst[[l_type]], dat_test_norm, l_type)
        combat_tst_prob <-  predWrapper(combat_mod_lst[[l_type]], dat_test_norm, l_type)
        gmm_tst_prob <- predWrapper(gmm_mod_lst[[l_type]], dat_test_norm, l_type)
        onestep_res <- ensemble_wrapper_realdata(sgbatch_mod_lst, l_type, dat_test_norm, 
                                                 navg_weights, cs_weights_seq, reg_a_beta, reg_s_beta)
        
        
        ##  Evaluate performance
        tst_scores <- c(list(Batch=unadj_tst_prob), onestep_res$pred_test_lst, 
                        list(ComBat=combat_tst_prob, GMM=gmm_tst_prob,
                             Avg=onestep_res$pred_avg, n_Avg=onestep_res$pred_N_avg, 
                             CS_Avg=onestep_res$pred_cs_avg, Reg_a=onestep_res$pred_reg_a, 
                             Reg_s=onestep_res$pred_reg_s))
        perf_df <- perf_wrapper(perf_measures, tst_scores, group_test)
        
        ##  Cache results
        perf_df_lst[[l_type]] <- perf_df
        tst_scores_modlst[[l_type]] <- tst_scores
      }
      
      
      ####  Ensemble across models
      
      preds_crossmod <- lapply(tst_scores_modlst, function(x){
        do.call(cbind, x[paste0("Batch", batch_names)])
      })
      cm_onestep_res <- ensemble_crossmod_wrapperNew(preds_crossmod, length(batches_ind), 
                                                     cm_navg_weights, cm_cs_weights_seq, cm_reg_a_beta, cm_reg_s_beta)
      
      ## calculate performance
      tst_cm_scores <- list(Avg=cm_onestep_res$cm_avg, n_Avg=cm_onestep_res$cm_N_avg, 
                            CS_Avg=cm_onestep_res$cm_cs_avg, Reg_a=cm_onestep_res$cm_reg_a, 
                            Reg_s=cm_onestep_res$cm_reg_s)
      perf_crossmod_df <- perf_wrapper(perf_measures, tst_cm_scores, group_test)
      
      # Create crossmod results by combining preprocessing methods with ensemble methods
      # Get preprocessing methods (Batch, ComBat, GMM) from individual model results
      # Use RF as the source since all individual models should have the same preprocessing results
      
      # Create tst_scores for crossmod that includes preprocessing + ensemble methods
      tst_cm_scores_full <- c(list(Batch=unadj_tst_prob, ComBat=combat_tst_prob, GMM=gmm_tst_prob),
                              tst_cm_scores)  # tst_cm_scores contains Avg, n_Avg, CS_Avg, Reg_a, Reg_s
      
      # Generate performance metrics for the full crossmod scores
      perf_crossmod_full_df <- perf_wrapper(perf_measures, tst_cm_scores_full, group_test)
      
      perf_df_lst[["crossmod"]] <- perf_crossmod_full_df
      tst_scores_modlst[["crossmod"]] <- tst_cm_scores_full
      
      for(i in 1:length(perf_measures)){
        # Extract performance data more robustly
        perf_data <- list()
        for(model_name in names(perf_df_lst)){
          perf_res <- perf_df_lst[[model_name]]
          if(perf_measures[i] %in% rownames(perf_res)){
            if(model_name == "crossmod"){
              # For crossmod, include all methods (no column filtering needed)
              perf_data[[model_name]] <- perf_res[perf_measures[i], , drop=FALSE]
            } else {
              # For individual models, use the original column selection logic
              # Get all columns except batch-specific ones (typically columns 2-4 for 3 batches)
              # Keep columns: 1 (first method), and from (2+n_batches) onwards
              n_batches <- length(unique(batch))
              start_col <- 2 + n_batches
              if(ncol(perf_res) >= start_col){
                selected_cols <- c(1, start_col:ncol(perf_res))
                perf_data[[model_name]] <- perf_res[perf_measures[i], selected_cols]
              } else {
                # Fallback: just take the first column if structure is unexpected
                perf_data[[model_name]] <- perf_res[perf_measures[i], 1, drop=FALSE]
              }
            }
          }
        }
        
        if(length(perf_data) > 0){
          summary_df <- reshape2::melt(perf_data)
          summary_df$iteration <- b
          # Rename columns to match what the figure generation script expects
          colnames(summary_df) <- c("Method", "value", "Model", "Iteration")
        } else {
          next  # Skip if no data available
        }
        
        ## write out performances - NOTE: Using results_real_6studies directory
        first_file <- !file.exists(sprintf('/scripts/evaluations/robustifying/results_real_6studies/test%s_%s.csv', test_name, perf_measures[i]))
        write.table(summary_df, sprintf('/scripts/evaluations/robustifying/results_real_6studies/test%s_%s.csv', test_name, perf_measures[i]),
                    append=!first_file, col.names=first_file, row.names=FALSE, sep=",")
      }
      
      rm(boot_ind, dat_test, group_test)
    }
    b <- b + 1
  }
}

cat("6-study analysis completed successfully\n")