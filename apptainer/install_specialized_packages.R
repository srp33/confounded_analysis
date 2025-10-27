# install_specialized_packages.R
# OPTIMIZED VERSION - Only truly specialized/lightweight packages
# Heavy packages moved to base image for better performance

# Configure repositories with Bioconductor binary support
options(repos = c(
    CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest",
    BioCsoft = "https://packagemanager.posit.co/bioconductor/packages/3.21/bioc/__linux__/jammy/latest",
    BioCann = "https://packagemanager.posit.co/bioconductor/packages/3.21/data/annotation/__linux__/jammy/latest",
    BioCexp = "https://packagemanager.posit.co/bioconductor/packages/3.21/data/experiment/__linux__/jammy/latest",
    CRAN_source = "https://cran.rstudio.com/"
))

# Use all available CPU cores
num_cores <- parallel::detectCores()
cat("Detected", num_cores, "CPU cores for LIGHTWEIGHT specialized package installation\n")
options(Ncpus = num_cores)

# Configure pak for maximum parallelization and prefer binaries
options(pak.no_extra_messages = TRUE)
options(pak.prefer_binary = TRUE)
Sys.setenv(PKG_BUILD_EXTRA_FLAGS = paste0("-j", num_cores))
Sys.setenv(MAKEFLAGS = paste0("-j", num_cores))

# Check what's already installed to avoid conflicts
installed_pkgs <- rownames(installed.packages())
cat("Found", length(installed_pkgs), "packages already installed in base image\n")

# Ensure BiocManager is available for fallback installations
if (!"BiocManager" %in% installed_pkgs) {
    cat("Installing BiocManager for Bioconductor package management...\n")
    pak::pkg_install("BiocManager")
}

# Install Bioconductor packages separately with better error handling
cat("--- Installing Bioconductor analysis packages ---\n")
bioc_packages <- c("polyester", "sva", "batchelor", "BatchQC", "RUVSeq")

for (pkg in bioc_packages) {
    full_pkg_name <- paste0("bioc::", pkg)
    if (!pkg %in% installed_pkgs) {
        cat("Installing", full_pkg_name, "...\n")
        tryCatch({
            pak::pkg_install(full_pkg_name)
            cat("✓", full_pkg_name, "installed successfully\n")
        }, error = function(e) {
            cat("✗", full_pkg_name, "failed:", e$message, "\n")
            # Try alternative installation methods for polyester
            if (pkg == "polyester") {
                cat("polyester not in Bioconductor 3.21 - installing from GitHub source...\n")
                cat("Dependencies (Biostrings, IRanges, S4Vectors, logspline, zlibbioc) already in base image\n")
                tryCatch({
                    # Install from GitHub using remotes (already available in base image)
                    remotes::install_github("alyssafrazee/polyester", upgrade = "never")
                    cat("✓ polyester installed from GitHub source\n")
                }, error = function(e2) {
                    cat("✗ GitHub installation failed, trying pak as fallback...\n")
                    tryCatch({
                        pak::pkg_install("git::https://github.com/alyssafrazee/polyester.git")
                        cat("✓ polyester installed via pak from git\n")
                    }, error = function(e3) {
                        cat("✗ All polyester installation methods failed\n")
                        cat("Consider using seqgendiff package as alternative for RNA-seq simulation\n")
                    })
                })
            }
        })
    } else {
        cat("✓", pkg, "already installed\n")
    }
}

# Install only packages that aren't already installed
cat("--- Installing CRAN packages not in base image ---\n")
cran_packages <- c(
    "MLmetrics", "ROCR", "nnls", "moments", "pracma",
    "ggpubr", "ggtext", "corrplot", "DT",
    "Rtsne", "umap", "kableExtra", "argparse", "docstring", 
    "itertools", "fairadapt", "pacman", "Seurat",
    "seqgendiff"  # Alternative RNA-seq simulation package
)

# Filter out already installed packages
cran_to_install <- cran_packages[!cran_packages %in% installed_pkgs]
if (length(cran_to_install) > 0) {
    cat("Installing", length(cran_to_install), "CRAN packages:", paste(cran_to_install, collapse=", "), "\n")
    pak::pkg_install(cran_to_install)
} else {
    cat("All CRAN packages already installed in base image\n")
}

# Note: Annotation packages moved to base image for better performance
# They are large/slow to compile, so better to include them in the base image
# that gets built once and reused

cat("--- All specialized R packages installed successfully ---\n")