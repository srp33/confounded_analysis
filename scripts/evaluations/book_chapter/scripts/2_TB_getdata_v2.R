# ==============================================================================
# TB Data Retrieval & Harmonization Script - Version 2 (Strict Replication)
# ==============================================================================
# Purpose: Integrates heterogeneous TB gene expression datasets for batch effect 
#          analysis using the EXACT data topology from Zhang et al. (2020).
#
# COHORT DEMOGRAPHICS & SOURCES (Strict Replication):
# 1. Batch 1 (South Africa): GSE79362 - Zak et al. (2016).
#    - Population: South African Adolescents (Western Cape). Predominantly Black African/Mixed.
#    - Platform: Illumina RNA-seq.
#
# 2. Batch 2 (India): GSE119370 - Sweeney et al. (2016).
#    - Population: Indian adults (actual India cohort, not UK proxy).
#    - Platform: Illumina HumanHT-12 v4.
#
# 3. Validation (USA): GSE73408 - Walter et al. (2016).
#    - Population: US Adults (Denver, CO). Mixed demographics.
#    - Platform: Affymetrix HuGene 1.1 ST.
#
# 4. Validation (SA/Malawi): GSE37250 & GSE39941 - Kaforou/Berry.
#    - Population: South African & Malawian adults.
#    - Platform: Illumina HumanHT-12 v4.
#
# DIFFERENCE FROM V1:
# - V1 used GSE107994 (Leicester/UK, South Asian ethnicity) as India proxy
# - V1 used GSE73408 (USA) as training batch
# - V2 uses GSE119370 (actual India cohort) for strict replication
# - V2 uses GSE73408 (USA) for validation only (as per original paper)
# ==============================================================================

rm(list=ls())

# --- Dependencies ---
required_packages <- c("GEOquery", "annotate", "hugene11sttranscriptcluster.db",
                       "SummarizedExperiment", "limma", "BatchQC", "ggplot2", 
                       "matrixStats", "illuminaHumanv4.db")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(paste("Package missing:", pkg, "- Please install via BiocManager::install()"))
  }
}

set.seed(123)
options(warn = 1) # Print warnings as they occur

# ==============================================================================
# 1. Process Batch 1: South Africa (GSE79362 - Zak et al.)
# ==============================================================================
cat("\n--- Processing Batch 1: South Africa (GSE79362) ---\n")

if (!file.exists("data/combined_sub.RData")) {
  stop("Critical Error: 'data/combined_sub.RData' is missing. This file contains the processed RNA-seq data.")
}

load("data/combined_sub.RData")

# Extract Training set (Zak et al.)
dat_africa <- train_expr
group_africa <- as.numeric(y_train) - 1 # Convert factor levels 1/2 -> 0/1

# Ensure matrix format
if (is(dat_africa, "SummarizedExperiment")) dat_africa <- assay(dat_africa, 1)

cat(sprintf("  > Loaded %d samples (%d TB, %d Control)\n", 
            ncol(dat_africa), sum(group_africa==1), sum(group_africa==0)))

# ==============================================================================
# 2. Process Batch 2: TRUE India (GSE119370 - Sweeney et al.)
# ==============================================================================
cat("\n--- Processing Batch 2: TRUE India (GSE119370) ---\n")

# Download Series Matrix
gse_ind <- getGEO("GSE119370", destdir="data", GSEMatrix=TRUE)

# Handle list return
if(length(gse_ind) > 0) gse_ind <- gse_ind[[1]]

# Extract Data
dat_india_raw <- exprs(gse_ind)

# Annotate Probes (Illumina IDs -> Gene Symbols)
# This is critical because Sweeney is Illumina, but Walter (USA) was Affymetrix
mapped_probes <- mapIds(illuminaHumanv4.db, keys=rownames(dat_india_raw), 
                        column="SYMBOL", keytype="PROBEID", multiVals="first")

# Filter NAs
keep_probes <- !is.na(mapped_probes)
dat_india <- dat_india_raw[keep_probes, ]
rownames(dat_india) <- mapped_probes[keep_probes]

# Metadata Cleaning
pheno_ind <- pData(gse_ind)
# Inspect columns to find disease group (usually 'characteristics_ch1.1')
# Sweeney labels: "Group: Active TB", "Group: LTBI", "Group: Healthy"
dis_stat <- pheno_ind$characteristics_ch1.1 

is_active <- grepl("Active TB", dis_stat, ignore.case=TRUE)
is_control <- grepl("LTBI|Healthy", dis_stat, ignore.case=TRUE)

# Filter Samples
keep_ind <- is_active | is_control
dat_india <- dat_india[, keep_ind]
group_india <- as.numeric(is_active[keep_ind])

# Log Transform Check
# Illumina data is often not logged in GEO Series Matrix
if(max(dat_india, na.rm=TRUE) > 50) {
  cat("  > Detected linear scale. Applying Log2(x+1)...\n")
  dat_india <- log2(dat_india + 1)
}

cat(sprintf("  > Retained %d samples (%d TB, %d Control)\n", 
            ncol(dat_india), sum(group_india==1), sum(group_india==0)))

# ==============================================================================
# 3. Harmonization & Intersection (Training Batches Only)
# ==============================================================================
cat("\n--- Harmonizing Gene Sets (Training) ---\n")

common_genes <- intersect(rownames(dat_africa), rownames(dat_india))
cat(sprintf("  > Overlapping genes across 2 training batches: %d\n", length(common_genes)))

# Subset
dat_africa <- dat_africa[common_genes, ]
dat_india <- dat_india[common_genes, ]

# Variance Filtering (Quality Control)
# Keep genes with variance > 0 and detected in > 2 samples per batch
pass_afr <- rowVars(dat_africa) > 0 & rowSums(dat_africa != 0) > 2
pass_ind <- rowVars(dat_india) > 0 & rowSums(dat_india != 0) > 2

keep_final <- pass_afr & pass_ind
cat(sprintf("  > Genes passing variance filter in BOTH batches: %d\n", sum(keep_final)))

# Build Training List
dat_lst <- list(
  Africa = dat_africa[keep_final, ], 
  India = dat_india[keep_final, ]
)
label_lst <- list(
  Africa = group_africa, 
  India = group_india
)

# ==============================================================================
# 4. Process Validation Set: USA (GSE73408 - Walter et al.)
# ==============================================================================
cat("\n--- Processing Validation: USA (GSE73408) ---\n")

# Download Series Matrix
gse <- getGEO("GSE73408", destdir="data", GSEMatrix=TRUE)[[1]]

# Annotate Probes to Gene Symbols
x <- hugene11sttranscriptclusterSYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])

# Map and filter NA symbols
gene_symbols <- sapply(featureNames(gse), function(s) ifelse(s %in% names(xx), xx[[s]], NA))
gse <- gse[!is.na(gene_symbols), ]
gene_symbols <- gene_symbols[!is.na(gene_symbols)]

dat_usa <- exprs(gse)
rownames(dat_usa) <- gene_symbols

# Metadata Cleaning
pheno <- pData(gse)
disease_status <- pheno$characteristics_ch1.2 # "clinical group: TB/LTBI/PNA"

# Filter: Keep Active TB (1) and LTBI (0), Remove Pneumonia
is_tb <- grepl("clinical group: TB$", disease_status, ignore.case=TRUE)
is_control <- grepl("clinical group: LTBI", disease_status, ignore.case=TRUE)
is_pna <- grepl("clinical group: PNA", disease_status, ignore.case=TRUE)

keep_usa <- !is_pna & (is_tb | is_control)
dat_usa <- dat_usa[, keep_usa]
group_usa <- as.numeric(is_tb[keep_usa]) # 1=TB, 0=Control

cat(sprintf("  > Retained %d samples (%d TB, %d Control)\n", 
            ncol(dat_usa), sum(group_usa==1), sum(group_usa==0)))

# ==============================================================================
# 5. Process Validation Sets: SA/Malawi (GSE37250 / GSE39941)
# ==============================================================================
cat("\n--- Processing Validation: SA/Malawi ---\n")

# Helper to annotate Illumina V4
annotate_illumina <- function(dat) {
  mapped <- mapIds(illuminaHumanv4.db, keys=rownames(dat), column="SYMBOL", keytype="PROBEID", multiVals="first")
  dat <- dat[!is.na(mapped), ]
  rownames(dat) <- mapped[!is.na(mapped)]
  return(dat)
}

# Download GSE37250
gse37 <- getGEO("GSE37250", destdir="data", GSEMatrix=TRUE)[[1]]
p37 <- pData(gse37)

# Filter SA & Malawi
idx_sa <- p37$`hiv status:ch1`=="HIV negative" & 
  p37$`geographical region:ch1`=="South Africa" & 
  p37$`disease state:ch1` %in% c("active tuberculosis", "latent TB infection")

idx_mw <- p37$`hiv status:ch1`=="HIV negative" & 
  p37$`geographical region:ch1`=="Malawi" & 
  p37$`disease state:ch1` %in% c("active tuberculosis", "latent TB infection")

dat_sa_val <- annotate_illumina(exprs(gse37[, idx_sa]))
dat_mw_val <- annotate_illumina(exprs(gse37[, idx_mw]))

grp_sa_val <- as.numeric(p37$`disease state:ch1`[idx_sa] == "active tuberculosis")
grp_mw_val <- as.numeric(p37$`disease state:ch1`[idx_mw] == "active tuberculosis")

# Download GSE39941
gse39 <- getGEO("GSE39941", destdir="data", GSEMatrix=TRUE)[[1]]
p39 <- pData(gse39)

idx_mw2 <- p39$`hiv status:ch1`=="HIV negative" & 
  p39$`geographical region:ch1`=="Malawi" & 
  p39$`disease status:ch1` %in% c("active tuberculosis", "latent TB infection")

dat_mw_val2 <- annotate_illumina(exprs(gse39[, idx_mw2]))
grp_mw_val2 <- as.numeric(p39$`disease status:ch1`[idx_mw2] == "active tuberculosis")

# ==============================================================================
# 6. Final Merge and Save
# ==============================================================================
cat("\n--- Finalizing Data Object ---\n")

# Align Validation sets to the training genes
train_genes <- rownames(dat_lst$Africa)
dat_usa <- dat_usa[intersect(rownames(dat_usa), train_genes), ]
dat_sa_val <- dat_sa_val[intersect(rownames(dat_sa_val), train_genes), ]
dat_mw_val <- dat_mw_val[intersect(rownames(dat_mw_val), train_genes), ]
dat_mw_val2 <- dat_mw_val2[intersect(rownames(dat_mw_val2), train_genes), ]

# Note: This effectively drops genes not found in validation sets. 
# A robust pipeline might impute or handle NAs, but for now we intersect.
final_genes <- Reduce(intersect, list(train_genes, rownames(dat_usa), 
                                      rownames(dat_sa_val), rownames(dat_mw_val), 
                                      rownames(dat_mw_val2)))

# Subset everything one last time
dat_lst <- lapply(dat_lst, function(x) x[final_genes, ])
dat_lst$USA <- dat_usa[final_genes, ]
dat_lst$GSE37250_SA <- dat_sa_val[final_genes, ]
dat_lst$GSE37250_M  <- dat_mw_val[final_genes, ]
dat_lst$GSE39941_M  <- dat_mw_val2[final_genes, ]

label_lst$USA <- group_usa
label_lst$GSE37250_SA <- grp_sa_val
label_lst$GSE37250_M  <- grp_mw_val
label_lst$GSE39941_M  <- grp_mw_val2

cat(sprintf("Final gene count: %d\n", length(final_genes)))
cat("Sample distribution:\n")
print(sapply(label_lst, table))

save(dat_lst, label_lst, file="data/TB_real_data_v2.RData")
cat("\n✅ Success! Saved to 'data/TB_real_data_v2.RData'\n")
cat("\nData Topology (Strict Replication):\n")
cat("  Training: Africa (GSE79362) + India (GSE119370)\n")
cat("  Validation: USA (GSE73408) + SA/Malawi (GSE37250, GSE39941)\n")
