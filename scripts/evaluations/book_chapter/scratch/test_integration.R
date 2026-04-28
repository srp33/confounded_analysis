testing_mode <- TRUE
source("scripts/classify_adjusters.R")

# Mock data
dat <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(dat) <- paste0("Gene", 1:10)
colnames(dat) <- paste0("Sample", 1:10)

dat_test <- matrix(rnorm(50), nrow = 10, ncol = 5)
rownames(dat_test) <- paste0("Gene", 1:10)
colnames(dat_test) <- paste0("Test", 1:5)

# Test Angel
cat("Testing Angel's method...\n")
res_angel <- adjust_angel(dat, debug = TRUE)
print(dim(res_angel))
print(range(res_angel))

# Test TDM
cat("\nTesting TDM method...\n")
tryCatch({
  res_tdm <- adjust_tdm(dat, dat_test, debug = TRUE)
  print(dim(res_tdm))
}, error = function(e) cat("TDM Error:", e$message, "\n"))

# Test RNABC
cat("\nTesting RNABC method...\n")
tryCatch({
  res_rnabc <- adjust_rnabc(dat, dat_test, debug = TRUE)
  print(dim(res_rnabc))
}, error = function(e) cat("RNABC Error:", e$message, "\n"))

# Test Shambhala2 (will fail if Octave/P0/Q0 missing, but we check syntax)
cat("\nTesting Shambhala2 method (syntax/setup)...\n")
# Create dummy P0 and Q0 for syntax check if needed, but here we just check if the function exists
if (exists("adjust_shambhala2")) {
  cat("adjust_shambhala2 exists.\n")
}
