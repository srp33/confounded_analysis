if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#packages = c("tidyverse", "doParallel", "readxl", "sva", "SCAN.UPC")
packages = c("pacman", "dplyr", "readr", "stringr", "tidyr", "tibble", "ggplot2", "gridExtra", "png", "magick", "colorspace", "pracma", "kableExtra", "Rtsne", "argparse", "docstring", "R.devices", "doParallel", "readxl", "sva", "SCAN.UPC")

tampor_packages = c("limma", "vsn", "doParallel", "ggplot2", "ggpubr")

BiocManager::install(packages)
BiocManager::install(tampor_packages)