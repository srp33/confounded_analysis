metadata <- read.delim("/scripts/evaluations/combat_seq/real_data_application/sample_metadata.tsv", stringsAsFactors = FALSE)
head(metadata)

library(SummarizedExperiment)

# Path to RDS files
rds_files <- list.files("/data/raw_counts/rds_output", pattern = "counts_only.*\\.rds$", full.names = TRUE)
count_list <- lapply(rds_files, readRDS)
output_file <- "/data/raw_counts/rds_output/signature_data.rds"

# Convert SummarizedExperiment to matrix if needed
count_list <- lapply(count_list, function(se) {
  if (inherits(se, "SummarizedExperiment")) {
    assay(se, "counts")
  } else {
    se
  }
})

# Example: combine and check sample names
all_sample_names <- unlist(lapply(count_list, colnames))
print(head(all_sample_names))

# Check overlap
sum(all_sample_names %in% metadata$geo_accession)

# Intersect gene names
common_genes <- Reduce(intersect, lapply(count_list, rownames))
count_list <- lapply(count_list, function(mat) mat[common_genes, ])

# Combine
combined_counts <- do.call(cbind, count_list)

# Ensure sample order matches metadata
metadata <- metadata[match(colnames(combined_counts), metadata$geo_accession), ]
stopifnot(all(metadata$geo_accession == colnames(combined_counts)))

# Keep only samples with batch info (remove NA)
metadata_filtered <- metadata[!is.na(metadata$batch), ]

# Subset combined_counts
combined_counts_filtered <- combined_counts[, colnames(combined_counts) %in% metadata_filtered$geo_accession]

# Check matching order
stopifnot(all(metadata_filtered$geo_accession == colnames(combined_counts_filtered)))

# Convert batch and group to factors
metadata_filtered$batch <- factor(metadata_filtered$batch)
metadata_filtered$group <- factor(metadata_filtered$group)

# Create SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = combined_counts_filtered),
  colData = metadata_filtered
)

# Save the combined object
cat("Saving filtered combined .rds file to ", output_file, "\n")
saveRDS(se, output_file)
