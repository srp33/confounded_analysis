lines <- readLines("scripts/classify_adjusters.R")
fs_start <- which(grepl("# Step 1.4: Feature selection", lines))[1]
fs_end <- which(grepl("cat\\(sprintf\\(\"Data preparation completed:", lines))[1] - 1

fs_block <- lines[fs_start:fs_end]

# Remove the block from its current position
lines_no_fs <- lines[-(fs_start:fs_end)]

# Find where batch correction ends
bc_end <- which(grepl("# ====================================================================", lines_no_fs) & 
                grepl("CLASSIFIER TRAINING", lines_no_fs))[1] - 1

# Insert the block before classifier training
final_lines <- c(lines_no_fs[1:bc_end], 
                 "  # Step 1.4: Feature selection (1000 highly variable genes)",
                 "  # Moved here to allow batch correction to use all genes",
                 "  cat(sprintf(\"[FEATURE SELECTION] Selecting 1000 most variable genes...\\n\"))",
                 "  reduced <- reduce_features(dat_corrected, dat_test_corrected, n_genes=1000)",
                 "  dat_corrected <- reduced$dat",
                 "  dat_test_corrected <- reduced$dat_test",
                 "",
                 lines_no_fs[(bc_end+1):length(lines_no_fs)])

writeLines(final_lines, "scripts/classify_adjusters.R")
