suppressPackageStartupMessages({
  library(reticulate)
})

# ---- Setup Python Environment ----
# Optionally move to your main script if you always use the same env
# use_python("/path/to/python_env", required = TRUE)

# Import scikit-learn modules once
sklearn <- import("sklearn")
np <- import("numpy")
metrics <- sklearn$metrics
ensemble <- sklearn$ensemble

# ---- Helper: run sklearn model ----
run_sklearn_model <- function(X_train, y_train, X_test, y_test, random_state = 42) {
  X_train_np <- np$array(as.matrix(X_train))
  y_train_np <- np$array(as.integer(y_train))
  X_test_np <- np$array(as.matrix(X_test))
  y_test_np <- np$array(as.integer(y_test))

  model <- ensemble$HistGradientBoostingClassifier(
    max_iter = as.integer(100),
    random_state = as.integer(random_state)
  )
  model$fit(X_train_np, y_train_np)

  y_pred <- model$predict(X_test_np)
  y_proba <- model$predict_proba(X_test_np)[, 2]

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
}
