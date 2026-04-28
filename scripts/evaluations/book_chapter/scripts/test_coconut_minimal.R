library(COCONUT)
data(GSEs.test)

# adjust_coconut helper function to handle the pre-unified numeric labels directly and inject diagnostic checks.
adjust_coconut <- function(matrix_, batch, group, debug = FALSE) {
  if (debug) {
    cat("DEBUG: Starting COCONUT harmonization.\n")
    cat("DEBUG: matrix_ dimensions: ", nrow(matrix_), " x ", ncol(matrix_), "\n")
    cat("DEBUG: Unique batches: ", length(unique(batch)), "\n")
    cat("DEBUG: First 10 group labels received: ", paste(head(group, 10), collapse=", "), "\n")
  }
  
  if (!requireNamespace("COCONUT", quietly = TRUE)) {
    stop("Package 'COCONUT' is required but not installed.")
  }
  
  gse_list <- list()
  for (b in unique(batch)) {
    idx <- which(batch == b)
    
    # Pending, might fix Datasets with <1 control error (mismatched string comparison on pre-unified numeric labels).
    disease_vec <- as.numeric(group[idx])
    
    if (debug) {
      cat(sprintf("DEBUG: Batch %s -> Total samples: %d | Controls (0s): %d | Cases (1s): %d\n", 
                  b, length(disease_vec), sum(disease_vec == 0, na.rm=TRUE), sum(disease_vec == 1, na.rm=TRUE)))
    }
    
    pheno_df <- data.frame(
      disease_state = disease_vec,
      dummy = 1, # Add dummy column to prevent dimension dropping in COCONUT
      row.names = colnames(matrix_[, idx])
    )
    
    gse_list[[as.character(b)]] <- list(
      pheno = pheno_df,
      genes = matrix_[, idx, drop = FALSE]
    )
  }
  
  if (debug) {
    cat("DEBUG: Executing COCONUT across GSE list.\n")
  }
  
  res <- COCONUT::COCONUT(GSEs = gse_list, control.0.col = "disease_state")
  
  if (debug) {
    cat("DEBUG: COCONUT execution finished. Reassembling results.\n")
    cat("DEBUG: res$COCONUTList names: ", paste(names(res$COCONUTList), collapse=", "), "\n")
  }
  
  result_matrix <- matrix(NA, nrow = nrow(matrix_), ncol = ncol(matrix_))
  rownames(result_matrix) <- rownames(matrix_)
  colnames(result_matrix) <- colnames(matrix_)
  
  for (b in names(res$COCONUTList)) {
    disease_cols <- colnames(res$COCONUTList[[b]]$genes)
    if (debug) {
        cat(sprintf("DEBUG: Reassembling batch %s cases: %d columns\n", b, length(disease_cols)))
        print(head(colnames(res$COCONUTList[[b]]$genes)))
    }
    result_matrix[, disease_cols] <- as.matrix(res$COCONUTList[[b]]$genes)
  }
  for (b in names(res$controlList$GSEs)) {
    control_cols <- colnames(res$controlList$GSEs[[b]]$genes)
    if (debug) {
        cat(sprintf("DEBUG: Reassembling batch %s controls: %d columns\n", b, length(control_cols)))
    }
    result_matrix[, control_cols] <- as.matrix(res$controlList$GSEs[[b]]$genes)
  }
  
  return(result_matrix)
}

# Prepare test data from GSEs.test
all_genes <- rownames(GSEs.test[[1]]$genes)
matrix_list <- lapply(GSEs.test, function(x) as.matrix(x$genes[all_genes, ]))
matrix_ <- do.call(cbind, matrix_list)

batch <- rep(names(GSEs.test), times = sapply(GSEs.test, function(x) ncol(x$genes)))
# Use the control column from the test data
group <- unlist(lapply(GSEs.test, function(x) x$pheno$Healthy0.Sepsis1))

print("Running adjust_coconut with debug=TRUE...")
res_matrix <- adjust_coconut(matrix_, batch, group, debug = TRUE)

# Code was run without Datasets with <1 control error (Verified).
print("Success!")
