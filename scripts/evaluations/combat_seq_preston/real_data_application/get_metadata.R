# On your local computer (with internet)
library(GEOquery)

gse83083 <- getGEO("GSE83083", GSEMatrix = TRUE)[[1]]
gse59765 <- getGEO("GSE59765", GSEMatrix = TRUE)[[1]]

# Extract useful metadata
metadata_83083 <- pData(gse83083)[, c("geo_accession", "title", "source_name_ch1")]
metadata_59765 <- pData(gse59765)[, c("geo_accession", "title", "source_name_ch1")]

# Add batch and group based on title/source_name
metadata_83083$batch <- ifelse(grepl("HER2", metadata_83083$title), 1,
                         ifelse(grepl("GFP18", metadata_83083$title), 1,
                         ifelse(grepl("KRAS", metadata_83083$title), 3, 
                         ifelse(grepl("GFP30", metadata_83083$title), 3, NA))))
metadata_83083$group <- ifelse(grepl("GFP18", metadata_83083$title), "gfp18",
                         ifelse(grepl("GFP30", metadata_83083$title), "gfp30", 
                         ifelse(grepl("HER2", metadata_83083$title), "her2",
                         ifelse(grepl("KRAS", metadata_83083$title), "kraswt", NA))))

metadata_59765$batch <- ifelse(grepl("Control", metadata_59765$title), 2,
                         ifelse(grepl("EGFR", metadata_59765$title), 2, NA))
metadata_59765$group <- ifelse(grepl("GFP", metadata_59765$title), "gfp_for_egfr",
                         ifelse(grepl("EGFR", metadata_59765$title), "egfr", NA))

# Combine and write to TSV
combined_metadata <- rbind(metadata_83083, metadata_59765)
write.table(combined_metadata, file = "sample_metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)