# install_annotation_packages.R
# ANNOTATION STAGE - Custom annotation packages from external sources
# This stage is separate to allow better debugging of annotation package issues

# Configure repositories
options(repos = c(
    CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest",
    BioCsoft = "https://packagemanager.posit.co/bioconductor/packages/3.21/bioc/__linux__/jammy/latest",
    BioCann = "https://packagemanager.posit.co/bioconductor/packages/3.21/data/annotation/__linux__/jammy/latest",
    BioCexp = "https://packagemanager.posit.co/bioconductor/packages/3.21/data/experiment/__linux__/jammy/latest",
    CRAN_source = "https://cran.rstudio.com/"
))

# Use all available CPU cores
num_cores <- parallel::detectCores()
cat("Detected", num_cores, "CPU cores for annotation package installation\n")
options(Ncpus = num_cores)

# Configure pak for maximum parallelization
options(pak.no_extra_messages = TRUE)
options(pak.prefer_binary = TRUE)
Sys.setenv(PKG_BUILD_EXTRA_FLAGS = paste0("-j", num_cores))
Sys.setenv(MAKEFLAGS = paste0("-j", num_cores))

# Check what's already installed
installed_pkgs <- rownames(installed.packages())
cat("Found", length(installed_pkgs), "packages already installed\n")

# Install custom annotation packages from mbni.org
cat("--- Installing custom annotation packages from mbni.org ---\n")
cat("These packages provide custom CDF environments for Affymetrix arrays\n")

custom_packages <- list(
    list(
        name = "hgu133ahsentrezgprobe", 
        url = "http://brainarray.mbni.med.umich.edu/bioc/src/contrib/hgu133ahsentrezgprobe_25.0.0.tar.gz",
        description = "Custom CDF for HG-U133A arrays with Entrez Gene mappings"
    ),
    list(
        name = "hgu133plus2hsentrezgprobe", 
        url = "http://brainarray.mbni.med.umich.edu/bioc/src/contrib/hgu133plus2hsentrezgprobe_25.0.0.tar.gz",
        description = "Custom CDF for HG-U133 Plus 2.0 arrays with Entrez Gene mappings"
    )
)

tmpDir <- tempdir()
success_count <- 0
total_count <- length(custom_packages)

for (pkg_info in custom_packages) {
    cat("\n--- Processing", pkg_info$name, "---\n")
    cat("Description:", pkg_info$description, "\n")
    cat("URL:", pkg_info$url, "\n")
    
    if (pkg_info$name %in% installed_pkgs) {
        cat("✓", pkg_info$name, "already installed\n")
        success_count <- success_count + 1
        next
    }
    
    tryCatch({
        cat("Downloading", pkg_info$name, "...\n")
        pkgFilePath <- file.path(tmpDir, basename(pkg_info$url))
        
        # Download with timeout and progress
        download.file(pkg_info$url, pkgFilePath, timeout = 60, mode = "wb")
        
        # Verify download
        if (!file.exists(pkgFilePath) || file.size(pkgFilePath) == 0) {
            stop("Download failed or file is empty")
        }
        
        cat("Installing", pkg_info$name, "from downloaded file...\n")
        pak::pkg_install(pkgFilePath)
        
        # Verify installation
        if (pkg_info$name %in% rownames(installed.packages())) {
            cat("✓", pkg_info$name, "installed and verified successfully\n")
            success_count <- success_count + 1
        } else {
            cat("✗", pkg_info$name, "installation completed but package not found\n")
        }
        
    }, error = function(e) {
        cat("✗", pkg_info$name, "failed:", e$message, "\n")
        
        # Provide helpful debugging information
        if (grepl("connect|timeout|resolve", e$message, ignore.case = TRUE)) {
            cat("  → Network issue: mbni.org server may be down or unreachable\n")
            cat("  → Try again later or check server status\n")
        } else if (grepl("download|404|not found", e$message, ignore.case = TRUE)) {
            cat("  → Download issue: file may have moved or been removed\n")
            cat("  → Check if URL is still valid:", pkg_info$url, "\n")
        } else {
            cat("  → Installation issue: package may have dependency problems\n")
        }
        
        cat("  → Continuing without", pkg_info$name, "\n")
        cat("  → Some microarray analyses may not work without this package\n")
    })
}

# Summary
cat("\n--- Annotation Package Installation Summary ---\n")
cat("Successfully installed:", success_count, "out of", total_count, "custom packages\n")

if (success_count == total_count) {
    cat("✓ All custom annotation packages installed successfully\n")
} else if (success_count > 0) {
    cat("⚠ Partial success - some packages failed but build can continue\n")
} else {
    cat("✗ No custom annotation packages installed - check network connectivity\n")
    cat("  → Standard Bioconductor annotation packages are still available\n")
    cat("  → Most analyses will work fine without custom CDFs\n")
}

cat("--- Annotation stage completed ---\n")