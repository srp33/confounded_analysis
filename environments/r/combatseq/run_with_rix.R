# run_with_rix.R
# Helper script to install rix and run generate_env.R

# Set up a temporary user library for rix installation
temp_lib <- file.path(tempdir(), "R_libs")
dir.create(temp_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(temp_lib, .libPaths()))

cat("R library paths:\n")
cat(paste("  ", .libPaths(), collapse = "\n"), "\n\n")

# Install rix if not already installed
if (!requireNamespace("rix", quietly = TRUE)) {
    cat("Installing rix package from CRAN to temporary library...\n")
    install.packages("rix", repos = "https://cloud.r-project.org", lib = temp_lib)
}

# Load rix
library(rix)
cat("rix version:", as.character(packageVersion("rix")), "\n\n")

# Source the generator script
cat("Running generate_env.R...\n")
source("generate_env.R")
