#!/usr/bin/env Rscript
# install_github_packages.R
# Install R packages from GitHub that aren't available via CRAN/Bioconductor
# Run after: source environments/load_envs.sh book_chapter

cat("Installing GitHub-only R packages...\n\n")

# Ensure remotes is available
if (!requireNamespace("remotes", quietly = TRUE)) {
    cat("Installing remotes package...\n")
    install.packages("remotes", repos = "https://cloud.r-project.org")
}

# GitHub packages to install
github_packages <- list(
    polyester = "alyssafrazee/polyester",
    # Add more GitHub packages here as needed
    # package_name = "github_user/repo"
)

for (pkg_name in names(github_packages)) {
    repo <- github_packages[[pkg_name]]
    
    if (requireNamespace(pkg_name, quietly = TRUE)) {
        cat("✓", pkg_name, "already installed\n")
    } else {
        cat("Installing", pkg_name, "from", repo, "...\n")
        tryCatch({
            remotes::install_github(repo, upgrade = "never")
            cat("✓", pkg_name, "installed successfully\n\n")
        }, error = function(e) {
            cat("✗", pkg_name, "installation failed:", e$message, "\n\n")
        })
    }
}

cat("GitHub package installation complete!\n")
