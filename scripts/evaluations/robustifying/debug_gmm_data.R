#!/usr/bin/env Rscript

# Debug script to check GMM data in results files
setwd("/scripts/evaluations/robustifying")
source("code/helper.R")

# Check what's in the US test file
cat("=== Checking testUS_mxe.csv ===\n")
res <- read.csv("results/testUS_mxe.csv")
colnames(res) <- c("Method", "value", "Model", "Iteration")

cat("Unique methods:", paste(unique(res$Method), collapse=", "), "\n")
cat("GMM present:", "GMM" %in% res$Method, "\n")
cat("Number of GMM rows:", sum(res$Method == "GMM"), "\n")

# Check a few GMM values
gmm_data <- res[res$Method == "GMM" & res$Model == "crossmod", ]
cat("First 5 GMM values:\n")
print(head(gmm_data, 5))

# Check if filtering works
selected_methods <- c("Batch", "ComBat", "GMM", "n_Avg", "CS_Avg", "Reg_s")
filtered <- res[res$Method %in% selected_methods & res$Model == "crossmod", ]
cat("Methods after filtering:", paste(unique(filtered$Method), collapse=", "), "\n")
cat("GMM in filtered data:", "GMM" %in% filtered$Method, "\n")