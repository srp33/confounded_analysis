library(COCONUT)
# Create dummy data
set.seed(42)
genes <- 10
samples_per_batch <- 5
batch1_dat <- matrix(rnorm(genes * samples_per_batch), nrow=genes)
batch2_dat <- matrix(rnorm(genes * samples_per_batch), nrow=genes)

# Batch 1 has controls and cases
pheno1 <- c(0, 0, 1, 1, 1)
# Batch 2 has ONLY cases
pheno2 <- c(1, 1, 1, 1, 1)

gse_list <- list(
  batch1 = list(pheno = data.frame(disease_state = pheno1, dummy=1), genes = batch1_dat),
  batch2 = list(pheno = data.frame(disease_state = pheno2, dummy=1), genes = batch2_dat)
)

tryCatch({
  res <- COCONUT::COCONUT(GSEs = gse_list, control.0.col = "disease_state")
  print("Success!")
}, error = function(e) {
  print(paste("Error:", e$message))
})
