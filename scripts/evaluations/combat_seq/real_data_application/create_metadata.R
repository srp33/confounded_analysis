metadata <- read.delim("sample_metadata.tsv", stringsAsFactors = FALSE)
head(metadata)

library(SummarizedExperiment)

# Path to RDS files
rds_files <- list.files("/data/raw_counts/rds_output", pattern = "\\.rds$", full.names = TRUE)
count_list <- lapply(rds_files, readRDS)

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

col_data <- metadata
rownames(col_data) <- col_data$geo_accession

se <- SummarizedExperiment(
  assays = list(counts = combined_counts),
  colData = col_data
)

# Save the combined object
saveRDS(se, "/data/raw_counts/rds_output/signature_data.rds")
