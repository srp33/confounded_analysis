if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if  (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

pkgs <- list(
  polyester = "alyssafrazee/polyester"
)

# GitHub packages (install only if missing)
for (pkg in names(pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing GitHub package ", pkg, " from ", pkgs[[pkg]])
    tryCatch(
      remotes::install_github(pkgs[[pkg]], upgrade = "never"),
      error = function(e) message("Failed to install ", pkg, ": ", e$message)
    )
  } else {
    message("Package ", pkg, " already installed")
  }
}
