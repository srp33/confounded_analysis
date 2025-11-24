# Suppress warnings and messages for cleaner output
options(warn = -1)


perf_wrapper <- function(perf_names, tst_scores, ytest){
  perf_df <- lapply(perf_names, function(perf_name){
    as.data.frame(t(sapply(tst_scores, function(preds){
      if(perf_name=="mxe"){preds <- pmax(pmin(preds, 1 - 1e-15), 1e-15)}  
      rocr_pred <- prediction(preds, as.numeric(as.character(ytest)))
      if(perf_name %in% c("acc", "f", "err")){
        curr_perf <- performance(rocr_pred, perf_name)  
        return(curr_perf@y.values[[1]][which.min(abs(curr_perf@x.values[[1]]-0.5))])
      }else{
        curr_perf <- performance(rocr_pred, perf_name)
        return(as.numeric(curr_perf@y.values))
      }
    })))
  })
  names(perf_df) <- perf_names
  perf_df <- do.call(rbind, perf_df)
  return(perf_df)
}

# Calculate confusion matrix elements and derived metrics
confusion_matrix_wrapper <- function(tst_scores, ytest, threshold = 0.5) {
  # Convert test labels to numeric
  ytest_numeric <- as.numeric(as.character(ytest))
  
  # Calculate confusion matrix elements for each method
  confusion_metrics <- sapply(tst_scores, function(preds) {
    # Convert predictions to binary using threshold
    pred_binary <- as.numeric(preds >= threshold)
    
    # Calculate confusion matrix elements
    tp <- sum(pred_binary == 1 & ytest_numeric == 1)
    fp <- sum(pred_binary == 1 & ytest_numeric == 0)
    tn <- sum(pred_binary == 0 & ytest_numeric == 0)
    fn <- sum(pred_binary == 0 & ytest_numeric == 1)
    
    # Calculate derived metrics
    # Matthews Correlation Coefficient
    mcc_denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc <- if(mcc_denom == 0) 0 else (tp * tn - fp * fn) / mcc_denom
    
    # Precision and Recall
    precision <- if(tp + fp == 0) 0 else tp / (tp + fp)
    recall <- if(tp + fn == 0) 0 else tp / (tp + fn)
    
    # Specificity
    specificity <- if(tn + fp == 0) 0 else tn / (tn + fp)
    
    # Balanced Accuracy
    sensitivity <- recall  # Same as recall
    balanced_acc <- (sensitivity + specificity) / 2
    
    return(c(
      tp = tp, fp = fp, tn = tn, fn = fn,
      mcc = mcc, precision = precision, recall = recall, 
      specificity = specificity, balanced_acc = balanced_acc
    ))
  })
  
  # Convert to data frame with proper structure
  # confusion_metrics should be a matrix with 9 rows (metrics) and n columns (methods)
  # After transpose: n rows (methods) and 9 columns (metrics)
  confusion_transposed <- t(confusion_metrics)
  confusion_df <- as.data.frame(confusion_transposed)
  
  # The rownames should be the metric names, but we transposed so they become column names
  expected_colnames <- c("tp", "fp", "tn", "fn", "mcc", "precision", "recall", "specificity", "balanced_acc")
  if(ncol(confusion_df) == length(expected_colnames)) {
    colnames(confusion_df) <- expected_colnames
  }
  
  # Now transpose back to get metrics as rows and methods as columns (like perf_df)
  confusion_df <- as.data.frame(t(confusion_df))
  
  # Remove any problematic row names
  if(!is.null(rownames(confusion_df))) {
    # Ensure row names are valid
    rownames(confusion_df) <- make.names(rownames(confusion_df), unique = TRUE)
  }
  
  return(confusion_df)
}



####  Helpers related to inverse gamma distribution
mv2ab <- function(m, v){
  a <- 2 + m^2/v
  b <- m * (a-1)
  return(list(alpha=a, beta=b))
}

ab2mv <- function(a, b){
  m <- b / (a-1)
  v <- b^2 / ((a-1)^2*(a-2))
  return(list(mean=m, var=v))
}


####  Take subset of dataset
reduceSize <- function(dat, y, N){
  reduced_ctrl <- sample(which(y==0)) 
  reduced_case <- sample(which(y==1))
  #identical(sort(reduced_ctrl), which(y==0))
  
  reduced_ctrl <- reduced_ctrl[1:N] 
  reduced_case <- reduced_case[1:N]
  
  reduced_indices <- c(reduced_ctrl, reduced_case)
  dat <- dat[, reduced_indices]
  y <- y[reduced_indices]
  return(list(dat=dat, y=y))
}


####  Split a dataset into man-made batches
splitBatch <- function(condition, N_batch){
  # split samples into case / control groups
  case_ind <- which(condition==1)
  ctrl_ind <- which(condition==0)
  
  # split each condition group into N_batch batches
  batches_ind_case <- split(case_ind, sample(N_batch,length(case_ind),replace=TRUE))
  batches_ind_ctrl <- split(ctrl_ind, sample(N_batch,length(ctrl_ind),replace=TRUE))
  #print((sum(sapply(batches_ind_case,length))==length(case_ind)) & (sum(sapply(batches_ind_ctrl,length))==length(ctrl_ind)))
  
  # combine case / control samples in each batch
  batches_ind <- list()
  for(i in 1:N_batch){
    batches_ind[[i]] <- sort(c(batches_ind_case[[i]], batches_ind_ctrl[[i]]))
  }
  return(batches_ind)
}


####  Take subset from dataset
subsetBatch <- function(condition, N_sample_size, N_batch){
  # split samples into case / control groups
  case_ind <- which(condition==1)
  ctrl_ind <- which(condition==0)
  
  # number of controls and cases to take
  N_ctrl <- N_sample_size / 2
  N_case <- N_sample_size / 2
  if(N_ctrl*N_batch > length(ctrl_ind) | N_case*N_batch > length(case_ind)){stop("Not enough samples to subset!")}
  
  # split each condition group into N_batch batches
  batches_ind_case <- split(case_ind, sample(N_batch,length(case_ind),replace=TRUE))
  while(any(sapply(batches_ind_case, length) < N_case)){batches_ind_case <- split(case_ind, sample(N_batch,length(case_ind),replace=TRUE))}
  batches_ind_case <- lapply(batches_ind_case, function(x){x[1:N_case]})
  
  batches_ind_ctrl <- split(ctrl_ind, sample(N_batch,length(ctrl_ind),replace=TRUE))
  while(any(sapply(batches_ind_ctrl, length) < N_ctrl)){batches_ind_ctrl <- split(ctrl_ind, sample(N_batch,length(ctrl_ind),replace=TRUE))}
  batches_ind_ctrl <- lapply(batches_ind_ctrl, function(x){x[1:N_ctrl]})
  
  # combine case / control samples in each batch
  batches_ind <- list()
  for(i in 1:N_batch){
    batches_ind[[i]] <- sort(c(batches_ind_case[[i]], batches_ind_ctrl[[i]]))
  }
  return(batches_ind)
}


####  Simulate batch effect based on ComBat assumption
simBatch <- function(dat, condition, batches_ind, batch, hyper_pars){
  n_batches <- sapply(batches_ind, length) # number of samples in each batch
  n_genes <- nrow(dat)
  
  ## Organize hyper batch parameters
  batch_par <- list()
  for(i in 1:length(n_batches)){
    batch_par[[i]] <- sapply(hyper_pars, function(item){item[i]}) # mean, sd of gaussian; alpha, beta for InvGamma
  }
    
  ## Simulate batch parameters from hyper-pars
  gamma <- delta2 <- list()
  for(i in 1:length(n_batches)){
    gamma[[i]] <- rnorm(n_genes, mean=batch_par[[i]]["hyper_mu"], sd=batch_par[[i]]["hyper_sd"])
            delta2[[i]] <- MCMCpack::rinvgamma(n_genes, shape=batch_par[[i]]["hyper_alpha"], scale=batch_par[[i]]["hyper_beta"])
  }
    
  ## Simulate batch effect
  # fit linear model to data with no batch parameters, calculate residual variance
  X <- model.matrix(~Condition, data=data.frame(Condition=condition))
  beta <- solve(t(X) %*% X) %*% t(X) %*% t(dat)
  resid <- dat - t(X %*% beta)
  #range(apply(resid,1,mean)); range(apply(resid,1,var))
  
  # spike-in batch variance: multiply by condition adjusted data with delta
  resid_varbatch <- matrix(NA, nrow=nrow(dat), ncol=ncol(dat), dimnames=dimnames(dat))
  for(j in 1:length(n_batches)){
    curr_resid <- resid[, batches_ind[[j]]]
    spikein_var <- lapply(1:n_batches[j], function(col_ind){curr_resid[, col_ind] * sqrt(delta2[[j]])})
    resid_varbatch[, batches_ind[[j]]] <- do.call(cbind, spikein_var) 
  }
  #sapply(1:5, function(k){mean(apply(resid[, batches_ind[[k]]],1,var))})
  #sapply(1:5, function(k){mean(apply(resid_varbatch[, batches_ind[[k]]],1,var))})
  
  # construct mean batch parameter design matrix using gamma
  X_batch <- model.matrix(~-1+Batch, data=data.frame(Batch=factor(batch)))
  gamma_vec <- do.call(rbind, gamma)  #apply(gamma_vec,1,mean)
  
  # new data with added batch effect
  new_dat <- t(cbind(X, X_batch) %*% rbind(beta, gamma_vec)) + resid_varbatch
  if(!identical(rownames(new_dat), rownames(dat))){stop("BUG in simBatch function!")
  }else{colnames(new_dat) <- colnames(dat)}
  
  res <- list(new_dat=new_dat, batch_par=batch_par)
  return(res)
}


####  Gene-wise normalize datasets (z-score scaling)
normalizeData <- function(dat){
  dat_norm <- t(apply(dat, 1, scale, center=TRUE, scale=TRUE))
  dimnames(dat_norm) <- dimnames(dat)
  return(dat_norm)
}


#### Train pipeline: ref combat test with train, and fit learner
#train_set=train_expr_norm; train_label=y_train; test_set=test_expr_norm; lfit=learner_fit
trainPipe <- function(train_set, test_set, train_label, lfit=learner_fit){
  pred_res <- lfit(trn_set=train_set, y_trn=train_label, tst_set=test_set)
  return(pred_res)
}


####  Get functions corresponding to learner type
getPredFunctions <- function(learner_type){
  if(learner_type=="logistic"){return(predLogistic_pp)
  }else if(learner_type=="lasso"){return(predLasso_pp)  #return(predLasso)
  }else if(learner_type=="elnet"){return(predElnet_pp)  # Keep backward compatibility
  }else if(learner_type=="elasticnet"){return(predElnet_pp)  # New name
  }else if(learner_type=="svm"){return(predSVM)
  }else if(learner_type=="rf"){return(predRF_pp)  #return(predRF)
  }else if(learner_type=="nnet"){return(predNnet_pp)  # Keep backward compatibility
  }else if(learner_type=="nn"){return(predNnet_pp)  # New name
  }else if(learner_type=="naivebayes"){return(predNB)
  }else if(learner_type=="knn"){return(predKNN_pp) 
  }else if(learner_type=="xgboost"){return(predXGBoost_pp)
  }else if(learner_type=="rf_fs"){return(predRF_fs_pp) 
  }else if(learner_type=="plusminus"){return(predMas)
  }else{stop("Method not supported!")}
}


####  Some performance metrics 
LogLossBinary <- function(actual, predicted, eps=1e-15) {
  if(class(actual)=="factor"){actual <- as.numeric(as.character(actual))}
  predicted = pmin(pmax(predicted, eps), 1-eps)
  - (sum(actual * log(predicted) + (1 - actual) * log(1 - predicted))) / length(actual)
}

AccuracyBinary <- function(actual, predicted) {
  if(class(actual)=="factor"){actual <- as.numeric(as.character(actual))}
  predicted_bi <- as.numeric(predicted >= 0.5)
  sum(predicted_bi==actual) / length(actual)
}


####  Prediction functions
# lasso
predLasso <- function(
  trn_set,
  # gene-by-sample expression matrix for training
  tst_set,
  # gene-by-sample expression matrix for test
  y_trn
  # response of training set, binary & numeric
){
  library(glmnet, quietly = TRUE)
  obj <- cv.glmnet(x=t(trn_set), y=factor(y_trn), family="binomial", alpha=1,
                   lambda=10^(seq(from=-4, to=4, by=0.1)),
                   type.measure="mse", nfolds=10)#, intercept=FALSE)
  best_lambda <- obj$lambda.1se
  mod_logit <- glmnet(x=t(trn_set), y=factor(y_trn), family="binomial", alpha=1,
                      lambda=best_lambda)#, intercept=FALSE)

  # predictions
  pred_train_prob <- predict(mod_logit, t(trn_set), type="response")[,1]
  pred_test_prob <-  predict(mod_logit, t(tst_set), type="response")[,1]
  #type "response" gives the fitted probabilities for "binomial",
  #type "class" produces the class label corresponding to the maximum probability

  res <- list(beta=c(mod_logit$a0, mod_logit$beta[,1]),
              pred_trn_prob=pred_train_prob, pred_tst_prob=pred_test_prob)#,
              #pred_trn_class=pred_train_class, pred_tst_class=pred_test_class)
  return(res)
}

# elastic net
# predElnet <- function(
#   trn_set,
#   # gene-by-sample expression matrix for training
#   tst_set, 
#   # gene-by-sample expression matrix for test
#   y_trn 
#   # response of training set, binary & numeric
# ){
#   library(caret)
#   parGrid <- expand.grid(lambda=exp(seq(from=-10, to=10, by=1)),
#                          alpha=seq(from=0, to=1, by=0.1))
#   ctrl <- trainControl(method = "cv", number=4)
#   mod_elnet <- train(x=t(trn_set), y=as.factor(y_trn), family = "binomial",
#                      method="glmnet",
#                      trControl=ctrl,
#                      tuneGrid=parGrid)
#   pred_train_elnet <- predict(mod_elnet, t(trn_set), type="prob")[,"1"] 
#   pred_test_elnet <- predict(mod_elnet, t(tst_set), type="prob")[,"1"] 
#   # either "raw" or "prob", for the number/class predictions or class probabilities
#   res <- list(pred_trn=pred_train_elnet, pred_tst=pred_test_elnet)
#   return(res)
# }

# SVM
predSVM <- function(
  trn_set,
  # gene-by-sample expression matrix for training
  tst_set = NULL, 
  # gene-by-sample expression matrix for test
  y_trn
  # response of training set, binary & numeric
){
  library(e1071, quietly = TRUE)
  tune_ctrl <- tune.control(sampling="cross", cross=4)#, error.fun=LogLossBinary)
  obj <- tune(svm, train.x=t(trn_set), train.y=as.factor(y_trn),
              tunecontrol=tune_ctrl,
              ranges=list(type="C-classification",
                          kernel="linear",
                          cost=exp(seq(from=-10, to=10, by=2))))
  best_cost <- obj$best.parameters[,"cost"]
  mod_svm <- svm(x=t(trn_set), y=as.factor(y_trn),
                 type="C-classification", kernel="linear",
                 cost=best_cost, probability=TRUE)
  pred_train_svm <- predict(mod_svm, t(trn_set), probability=TRUE)
  pred_train_svm <- attr(pred_train_svm, "probabilities")[,"1"]
  if(!is.null(tst_set)){
    pred_test_svm <- predict(mod_svm, t(tst_set), probability=TRUE)
    pred_test_svm <- attr(pred_test_svm, "probabilities")[,"1"]
  }else{
    pred_test_svm <- NULL
  }
  res <- list(mod=mod_svm, pred_trn_prob=pred_train_svm, pred_tst_prob=pred_test_svm)
  return(res)
}

# random forest
predRF <- function(
  trn_set,
  # gene-by-sample expression matrix for training
  tst_set, 
  # gene-by-sample expression matrix for test
  y_trn
  # response of training set, binary & numeric
){
  library(caret, quietly = TRUE)
  
  trn_transposed <- t(trn_set)
  rownames(trn_transposed) <- NULL
  training_df <- data.frame(trn_transposed, as.factor(y_trn))
  colnames(training_df) <- c(paste("gene", 1:nrow(trn_set), sep=""), "response")
  rownames(training_df) <- 1:ncol(trn_set)
  
  tst_transposed <- t(tst_set)
  rownames(tst_transposed) <- NULL
  test_df <- data.frame(tst_transposed)
  colnames(test_df) <- paste("gene", 1:nrow(tst_set), sep="")
  rownames(test_df) <- 1:ncol(tst_set)
  
  ctrl <- trainControl(method="cv", number=10)
  f <- as.formula(paste("response ~ ", paste(colnames(training_df)[-ncol(training_df)], collapse= "+",sep="")))
  mod_rf <- train(form=f, data=training_df, 
                  method="rf", metric="Accuracy", 
                  tuneLength=5, trControl=ctrl)
  
  pred_train_rf <- predict(mod_rf, training_df, type="prob")[,"1"]
  pred_test_rf <- predict(mod_rf, test_df, type="prob")[,"1"]
  
  res <- list(pred_trn_prob=pred_train_rf, pred_tst_prob=pred_test_rf)
  return(res)
}

# neural net 
predNnet <- function(
  trn_set,
  # gene-by-sample expression matrix for training
  tst_set, 
  # gene-by-sample expression matrix for test
  y_trn
  # response of training set, binary & numeric
){
  library(caret, quietly = TRUE)
  parGrid <- expand.grid(size=seq(from=2, to=3, by=1),
                         decay=10^seq(from=-4, to=-1, by=0.5))
  ctrl <- trainControl(method = "cv", number=10)
  
  # Fix row names for training
  trn_transposed <- t(trn_set)
  rownames(trn_transposed) <- NULL
  
  mod_nnet <- train(x=trn_transposed, y=as.factor(y_trn), maxit=1000, MaxNWts=50000,
                    method="nnet", trace=FALSE,
                    trControl=ctrl, softmat=TRUE,
                    tuneGrid=parGrid)
  pred_train_nnet <- predict(mod_nnet, trn_transposed, type="prob")[,"1"]
  
  # Fix row names for test
  tst_transposed <- t(tst_set)
  rownames(tst_transposed) <- NULL
  pred_test_nnet <- predict(mod_nnet, tst_transposed, type="prob")[,"1"] 
  
  res <- list(pred_trn_prob=pred_train_nnet, pred_tst_prob=pred_test_nnet)
  return(res)
}

# mas-o-menos 
predMas <- function(
  trn_set,
  # gene-by-sample expression matrix for training
  tst_set, 
  # gene-by-sample expression matrix for test
  y_trn
  # response of training set, binary & numeric
){
  trn_set_norm <- t(scale(t(trn_set), center=TRUE, scale=TRUE))
  tst_set_norm <- t(scale(t(tst_set), center=TRUE, scale=TRUE))
  
  trn_norm_transposed <- t(trn_set_norm)
  rownames(trn_norm_transposed) <- NULL
  training_df_norm <- data.frame(trn_norm_transposed, as.factor(y_trn))
  colnames(training_df_norm) <- c(paste("gene", 1:nrow(trn_set_norm), sep=""), "response")
  rownames(training_df_norm) <- 1:ncol(trn_set_norm)
  
  alpha <- rep(0, nrow(trn_set_norm))
  for(j in 1:nrow(trn_set_norm)){
    f <- as.formula(paste("response ~ 0 +", paste(colnames(training_df_norm)[j], 
                                                  collapse= "+",sep="")))
    ctr <- glm.control(maxit=1000)
    mod_tmp <- glm(f, data=training_df_norm, family=binomial, control=ctr)
    alpha[j] <- coef(mod_tmp)
  }
  v <- (2*(alpha>0)-1)/sqrt(nrow(trn_set))
  
  pred_train_plusminus <- as.numeric(t(trn_set_norm) %*% as.matrix(v))
  pred_train_plusminus <- 1/(1+exp(- pred_train_plusminus))
  pred_test_plusminus <- as.numeric(t(tst_set_norm) %*% as.matrix(v))
  pred_test_plusminus <- 1/(1+exp(- pred_test_plusminus))
  
  res <- list(pred_trn_prob=pred_train_plusminus, pred_tst_prob=pred_test_plusminus)
  return(res)
}


####  Prediction functions slightly organized from Prasad Patil's github
predLasso_pp <- function(trn_set, tst_set=NULL, y_trn, ...){
  mod <- glmnet::cv.glmnet(x=t(trn_set), y=as.numeric(as.character(y_trn)), family="binomial", ...)
  pred_trn_prob <- as.vector(predict(mod, newx=t(trn_set), s="lambda.1se", type="response"))
  pred_trn_class <- as.vector(predict(mod, newx=t(trn_set), s="lambda.1se", type="class"))
  if(!is.null(tst_set)){
    pred_tst_prob <- as.vector(predict(mod, newx=t(tst_set), s="lambda.1se", type="response"))
    pred_tst_class <- as.vector(predict(mod, newx=t(tst_set), s="lambda.1se", type="class"))
  }else{
    pred_tst_prob <- NULL; pred_tst_class <- NULL
  }
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=pred_tst_prob,
              pred_trn_class=pred_trn_class, pred_tst_class=pred_tst_class))
}

#trn_set=train_set; y_trn=train_label; tst_set=test_refadj
predRF_pp <- function(trn_set, tst_set=NULL, y_trn){
  # Create training data with proper row names
  trn_transposed <- t(trn_set)
  rownames(trn_transposed) <- NULL  # Remove potentially problematic row names
  data <- data.frame(y=as.factor(y_trn), trn_transposed)
  mod <- ranger::ranger(y ~ ., data = data, write.forest=TRUE)
  mod_prob <- ranger::ranger(y ~ ., data = data, write.forest=TRUE, probability=TRUE)
  
  # For predictions, also fix row names
  trn_pred_data <- t(trn_set)
  rownames(trn_pred_data) <- NULL
  pred_trn_prob <- predict(mod_prob, data = data.frame(trn_pred_data))$predictions[, "1"]
  pred_trn_class <- predict(mod, data = data.frame(trn_pred_data))$predictions
  if(!is.null(tst_set)){
    tst_transposed <- t(tst_set)
    rownames(tst_transposed) <- NULL  # Remove potentially problematic row names
    newdata <- data.frame(tst_transposed)
    pred_tst_prob <- predict(mod_prob, data = newdata)$predictions[, "1"]
    pred_tst_class <- predict(mod, data = newdata)$predictions
  }else{
    pred_tst_prob <- NULL; pred_tst_class <- NULL
  }
  return(list(mod = mod_prob, pred_trn_prob=pred_trn_prob, pred_tst_prob=pred_tst_prob,
              pred_trn_class=pred_trn_class, pred_tst_class=pred_tst_class))
}

predRF_fs_pp <- function(trn_set, tst_set, y_trn){
  data <- data.frame(y=as.factor(y_trn), t(trn_set))
  newdata <- data.frame(t(tst_set))
  tout <- colttests(as.matrix(data[,-1]), data[,1], tstatOnly = TRUE)
  data <- data[, c("y",rownames(tout)[order(abs(tout[,1]), decreasing = T)[1:20]])]
  mod <- ranger(y ~ ., data = data, write.forest = TRUE, probability = T)
  predict(mod, data=newdata)$predictions[,2]
}

predNnet_pp <- function(trn_set, tst_set=NULL, y_trn){
  library(nnet, quietly = TRUE)
  
  # Simple validation
  if(is.null(trn_set) || is.null(y_trn)) {
    stop("Training set or labels are NULL")
  }
  
  # Create training data with proper row names
  trn_transposed <- t(trn_set)
  rownames(trn_transposed) <- NULL  # Remove potentially problematic row names
  data <- data.frame(y=as.factor(y_trn), trn_transposed)
  
  # Adjust network size for dataset
  n_samples <- nrow(data)
  network_size <- min(10, max(2, n_samples %/% 3))
  
  # Train neural network
  mod <- nnet(y ~ ., data = data, size = network_size, MaxNWts = 10000, 
              linout = FALSE, trace = FALSE, maxit = 200)
  
  pred_trn_prob <- as.vector(predict(mod, newdata = data[,-1]))
  pred_trn_class <- as.vector(predict(mod, newdata = data[,-1], type="class"))
  
  # Test predictions only if test set provided
  if(!is.null(tst_set)) {
    tst_transposed <- t(tst_set)
    rownames(tst_transposed) <- NULL  # Remove potentially problematic row names
    newdata <- data.frame(tst_transposed)
    pred_tst_prob <- as.vector(predict(mod, newdata = newdata))
    pred_tst_class <- as.vector(predict(mod, newdata = newdata, type="class"))
  } else {
    pred_tst_prob <- NULL
    pred_tst_class <- NULL
  }
  
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=pred_tst_prob,
              pred_trn_class=pred_trn_class, pred_tst_class=pred_tst_class))
}

# Logistic Regression with no regularization
predLogistic_pp <- function(trn_set, tst_set=NULL, y_trn){
  # Create training data with proper row names
  trn_transposed <- t(trn_set)
  rownames(trn_transposed) <- NULL  # Remove potentially problematic row names
  data <- data.frame(y=as.factor(y_trn), trn_transposed)
  
  mod <- glm(y ~ ., data = data, family = binomial())
  
  pred_trn_prob <- as.vector(predict(mod, newdata = data[,-1], type="response"))
  pred_trn_class <- as.vector(ifelse(pred_trn_prob >= 0.5, "1", "0"))
  
  if(!is.null(tst_set)){
    tst_transposed <- t(tst_set)
    rownames(tst_transposed) <- NULL  # Remove potentially problematic row names
    newdata <- data.frame(tst_transposed)
    pred_tst_prob <- as.vector(predict(mod, newdata = newdata, type="response"))
    pred_tst_class <- as.vector(ifelse(pred_tst_prob >= 0.5, "1", "0"))
  } else {
    pred_tst_prob <- NULL
    pred_tst_class <- NULL
  }
  
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=pred_tst_prob,
              pred_trn_class=pred_trn_class, pred_tst_class=pred_tst_class))
}

# ElasticNet (L1 + L2 regularization)
predElnet_pp <- function(trn_set, tst_set=NULL, y_trn){
  library(glmnet, quietly = TRUE)
  # Use alpha=0.5 for equal mix of L1 and L2 regularization
  mod <- cv.glmnet(x=t(trn_set), y=as.numeric(as.character(y_trn)), family="binomial", alpha=0.5)
  pred_trn_prob <- as.vector(predict(mod, newx=t(trn_set), s="lambda.1se", type="response"))
  pred_trn_class <- as.vector(predict(mod, newx=t(trn_set), s="lambda.1se", type="class"))
  if(!is.null(tst_set)){
    pred_tst_prob <- as.vector(predict(mod, newx=t(tst_set), s="lambda.1se", type="response"))
    pred_tst_class <- as.vector(predict(mod, newx=t(tst_set), s="lambda.1se", type="class"))
  }else{
    pred_tst_prob <- NULL; pred_tst_class <- NULL
  }
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=pred_tst_prob,
              pred_trn_class=pred_trn_class, pred_tst_class=pred_tst_class))
}

# K-Nearest Neighbors
predKNN_pp <- function(trn_set, tst_set=NULL, y_trn){
  library(class, quietly = TRUE)
  
  # Determine optimal k using cross-validation
  n_samples <- ncol(trn_set)
  k_values <- c(3, 5, 7, 9, 11)
  k_values <- k_values[k_values < n_samples]
  
  if(length(k_values) == 0) {
    k_opt <- min(3, n_samples - 1)
  } else {
    # Simple cross-validation to find best k
    cv_accuracy <- sapply(k_values, function(k) {
      # 5-fold CV or leave-one-out if too few samples
      if(n_samples < 10) {
        # Leave-one-out CV
        correct <- 0
        for(i in 1:n_samples) {
          train_idx <- setdiff(1:n_samples, i)
          trn_cv <- t(trn_set[, train_idx])
          tst_cv <- t(trn_set[, i, drop=FALSE])
          rownames(trn_cv) <- NULL
          rownames(tst_cv) <- NULL
          pred <- knn(train = trn_cv, 
                     test = tst_cv, 
                     cl = y_trn[train_idx], k = k)
          if(as.character(pred) == as.character(y_trn[i])) correct <- correct + 1
        }
        return(correct / n_samples)
      } else {
        # 5-fold CV
        folds <- cut(seq(1, n_samples), breaks = 5, labels = FALSE)
        correct <- 0
        total <- 0
        for(fold in 1:5) {
          test_idx <- which(folds == fold)
          train_idx <- which(folds != fold)
          if(length(test_idx) > 0 && length(train_idx) > 0) {
            trn_fold <- t(trn_set[, train_idx])
            tst_fold <- t(trn_set[, test_idx])
            rownames(trn_fold) <- NULL
            rownames(tst_fold) <- NULL
            pred <- knn(train = trn_fold, 
                       test = tst_fold, 
                       cl = y_trn[train_idx], k = k)
            correct <- correct + sum(as.character(pred) == as.character(y_trn[test_idx]))
            total <- total + length(test_idx)
          }
        }
        return(if(total > 0) correct / total else 0)
      }
    })
    k_opt <- k_values[which.max(cv_accuracy)]
  }
  
  # Train predictions (using leave-one-out to avoid overfitting)
  pred_trn_class <- character(n_samples)
  for(i in 1:n_samples) {
    train_idx <- setdiff(1:n_samples, i)
    trn_subset <- t(trn_set[, train_idx])
    tst_subset <- t(trn_set[, i, drop=FALSE])
    rownames(trn_subset) <- NULL
    rownames(tst_subset) <- NULL
    pred_trn_class[i] <- as.character(knn(train = trn_subset, 
                                         test = tst_subset, 
                                         cl = y_trn[train_idx], k = k_opt))
  }
  
  # Convert to probabilities (simple approach)
  pred_trn_prob <- as.numeric(pred_trn_class == "1")
  
  # Test predictions
  if(!is.null(tst_set)) {
    # Ensure proper matrix structure for KNN
    trn_for_knn <- t(trn_set)
    tst_for_knn <- t(tst_set)
    rownames(trn_for_knn) <- NULL
    rownames(tst_for_knn) <- NULL
    
    pred_tst_class <- as.character(knn(train = trn_for_knn, 
                                      test = tst_for_knn, 
                                      cl = y_trn, k = k_opt))
    pred_tst_prob <- as.numeric(pred_tst_class == "1")
  } else {
    pred_tst_prob <- NULL
    pred_tst_class <- NULL
  }
  
  # Store model info (k value used)
  mod <- list(k = k_opt, train_data = t(trn_set), train_labels = y_trn)
  
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=pred_tst_prob,
              pred_trn_class=pred_trn_class, pred_tst_class=pred_tst_class))
}

# XGBoost
predXGBoost_pp <- function(trn_set, tst_set=NULL, y_trn){
  library(xgboost, quietly = TRUE)
  
  # Prepare data - ensure labels are 0/1
  y_numeric <- as.numeric(as.factor(y_trn)) - 1  # Convert to 0/1 regardless of input format
  train_matrix <- xgb.DMatrix(data = t(trn_set), label = y_numeric)
  
  # Parameters for binary classification
  params <- list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    max_depth = 3,
    eta = 0.1,
    subsample = 0.8,
    colsample_bytree = 0.8,
    verbose = 0
  )
  
  # Adjust number of rounds based on dataset size
  n_samples <- ncol(trn_set)
  nrounds <- if(n_samples < 50) 50 else if(n_samples < 200) 100 else 200
  
  # Train model with early stopping if enough samples
  if(n_samples >= 20) {
    # Use cross-validation to determine optimal rounds
    cv_result <- xgb.cv(
      params = params,
      data = train_matrix,
      nrounds = nrounds,
      nfold = min(5, n_samples %/% 4),
      early_stopping_rounds = 10,
      verbose = 0,
      showsd = FALSE
    )
    # Check if best_iteration exists and is valid
    if(!is.null(cv_result$best_iteration) && length(cv_result$best_iteration) > 0 && !is.na(cv_result$best_iteration)) {
      best_nrounds <- cv_result$best_iteration
    } else {
      # Fallback: use the iteration with minimum test error
      if(!is.null(cv_result$evaluation_log) && nrow(cv_result$evaluation_log) > 0) {
        best_nrounds <- which.min(cv_result$evaluation_log$test_logloss_mean)
      } else {
        best_nrounds <- min(50, nrounds)
      }
    }
  } else {
    best_nrounds <- min(30, nrounds)
  }
  
  # Train final model
  mod <- xgb.train(
    params = params,
    data = train_matrix,
    nrounds = best_nrounds,
    verbose = 0
  )
  
  # Training predictions
  pred_trn_prob <- predict(mod, train_matrix)
  pred_trn_class <- as.character(as.numeric(pred_trn_prob >= 0.5))
  
  # Test predictions
  if(!is.null(tst_set)) {
    test_matrix <- xgb.DMatrix(data = t(tst_set))
    pred_tst_prob <- predict(mod, test_matrix)
    pred_tst_class <- as.character(as.numeric(pred_tst_prob >= 0.5))
  } else {
    pred_tst_prob <- NULL
    pred_tst_class <- NULL
  }
  
  return(list(mod=mod, pred_trn_prob=pred_trn_prob, pred_tst_prob=pred_tst_prob,
              pred_trn_class=pred_trn_class, pred_tst_class=pred_tst_class))
}

predWrapper <- function(mod, tst_set, function_name){
  if(function_name=='logistic'){
    tst_transposed <- t(tst_set)
    rownames(tst_transposed) <- NULL
    newdata <- data.frame(tst_transposed)
    res <- as.vector(predict(mod, newdata = newdata, type="response"))
  }else if(function_name=='lasso'){
    res <- as.vector(predict(mod, newx=t(tst_set), s="lambda.1se", type="response"))
  }else if(function_name=='elnet' || function_name=='elasticnet'){
    res <- as.vector(predict(mod, newx=t(tst_set), s="lambda.1se", type="response"))
  }else if(function_name=='svm'){
    res <- predict(mod, t(tst_set), probability=TRUE)
    res <- attr(res, "probabilities")[,"1"]
  }else if(function_name=='rf'){
    tst_transposed <- t(tst_set)
    rownames(tst_transposed) <- NULL
    newdata <- data.frame(tst_transposed)
    res <- predict(mod, data = newdata)$predictions[, "1"]
  }else if(function_name=='nnet' || function_name=='nn'){
    tst_transposed <- t(tst_set)
    rownames(tst_transposed) <- NULL
    newdata <- data.frame(tst_transposed)
    res <- as.vector(predict(mod, newdata = newdata))
  }else if(function_name=='knn'){
    tst_for_knn <- t(tst_set)
    rownames(tst_for_knn) <- NULL
    res <- knn(train = mod$train_data, test = tst_for_knn, cl = mod$train_labels, k = mod$k)
    res <- as.numeric(as.character(res) == "1")
  }else if(function_name=='xgboost'){
    test_matrix <- xgb.DMatrix(data = t(tst_set))
    res <- predict(mod, test_matrix)
  }
  return(res)
}

####  Quiet melt wrapper
# Suppress "No id variables; using all as measure variables" messages
quiet_melt <- function(...) {
  suppressMessages(melt(...))
}

####  Smart Results Loader
# Automatically finds results in the correct directory structure
load_results <- function(study_type = c("simulation", "real_4studies", "real_6studies"), 
                        study_name = NULL, metric = NULL, base_dir = "/scripts/evaluations/robustifying") {
  
  study_type <- match.arg(study_type)
  
  # Define possible result directories in order of preference
  possible_dirs <- switch(study_type,
    "simulation" = c(
      file.path(base_dir, "results")
    ),
    "real_4studies" = c(
      file.path(base_dir, "results", "real_4studies"),  # Unified structure (preferred)
      file.path(base_dir, "results"),                   # Copied files
      file.path(base_dir, "results_real_4studies")      # Original structure
    ),
    "real_6studies" = c(
      file.path(base_dir, "results", "real_6studies"),  # Unified structure (preferred)
      file.path(base_dir, "results"),                   # Copied files  
      file.path(base_dir, "results_real_6studies")      # Original structure
    )
  )
  
  # For real data, construct filename
  if (study_type != "simulation" && !is.null(study_name) && !is.null(metric)) {
    filename <- sprintf("test%s_%s.csv", study_name, metric)
    cat("DEBUG: Looking for specific file:", filename, "\n")
    
    # Try each directory until we find the file
    for (dir in possible_dirs) {
      filepath <- file.path(dir, filename)
      cat("DEBUG: Checking filepath:", filepath, "\n")
      cat("DEBUG: File exists:", file.exists(filepath), "\n")
      if (file.exists(filepath)) {
        cat("DEBUG: Found file at:", filepath, "\n")
        return(list(dir = dir, file = filepath))
      }
    }
    
    stop(sprintf("Could not find %s in any of the expected directories: %s", 
                filename, paste(possible_dirs, collapse = ", ")))
  }
  
  # For simulation data or when just getting directory
  for (dir in possible_dirs) {
    if (dir.exists(dir)) {
      csv_files <- list.files(dir, pattern = "\\.csv$")
      
      # Check if this directory contains the right type of files
      if (study_type == "simulation") {
        # Simulation files have patterns like "lasso_auc_batchN20_m0_v1.csv"
        sim_pattern_files <- grep("_batchN[0-9]+_m[0-9]+_v[0-9]+\\.csv$", csv_files)
        if (length(sim_pattern_files) >= 1) {
          dir <- sub("/$", "", dir)
          return(list(dir = dir, file = NULL))
        }
      } else {
        # Real data files have patterns like "testGSE37250_SA_auc.csv"
        real_pattern_files <- grep("^test[A-Za-z0-9_]+_[a-z]+\\.csv$", csv_files)
        if (length(real_pattern_files) >= 10) {
          dir <- sub("/$", "", dir)
          return(list(dir = dir, file = NULL))
        }
      }
    }
  }
  
  stop(sprintf("Could not find results directory for %s in: %s", 
              study_type, paste(possible_dirs, collapse = ", ")))
}
