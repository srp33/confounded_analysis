library(COCONUT)
# Create dummy data
set.seed(42)
genes <- 10
samples_per_batch <- 5
batch1_dat <- matrix(rnorm(genes * samples_per_batch), nrow=genes)
batch2_dat <- matrix(rnorm(genes * samples_per_batch), nrow=genes)
colnames(batch1_dat) <- paste0("B1_", 1:samples_per_batch)
colnames(batch2_dat) <- paste0("B2_", 1:samples_per_batch)

# Batch 1 has controls and cases
pheno1 <- c(0, 0, 1, 1, 1)
# Batch 2 has ONLY samples treated as controls
pheno2 <- c(0, 0, 0, 0, 0)

gse_list <- list(
  batch1 = list(pheno = data.frame(disease_state = pheno1, dummy=1, row.names=colnames(batch1_dat)), genes = batch1_dat),
  batch2 = list(pheno = data.frame(disease_state = pheno2, dummy=1, row.names=colnames(batch2_dat)), genes = batch2_dat)
)

tryCatch({
  res <- COCONUT::COCONUT(GSEs = gse_list, control.0.col = "disease_state")
  print("Success!")
}, error = function(e) {
  print(paste("Error:", e$message))
})
