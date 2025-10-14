rm(list=ls())
sapply(c("sva", "dplyr", "DESeq2", "ggplot2", "reshape2", "gridExtra", "scales", 
         "RUVSeq", "ggpubr", "BatchQC"), require, character.only=TRUE)

## Parameters (change paths when necessary)
data_dir <- "/scripts/evaluations/ComBat-seq/real_data_application"  # path to the signature data (.rds)
source("/scripts/evaluations/ComBat-seq/real_data_application/gfrn_helpers.R")  # path to gfrn_helpers.R
#source("/scripts/evaluations/ComBat-seq/ComBat_seq.R"); source("/scripts/evaluations/ComBat-seq/helper_seq.R")   
# path to the combat-seq scripts (or use the sva package on github, in which case comment out the above line)

pathway_regex <- c("her2", "^egfr", "kraswt")  
set.seed(1)


## Load data
sigdata <- readRDS(file.path(data_dir, "signature_data.rds"))
cts_mat <- assay(sigdata, "counts")   # count matrix (also have tpm and fpkm in there)
rownames(cts_mat) <- paste0("gene", 1:nrow(cts_mat))
batch <- colData(sigdata)$batch
group <- colData(sigdata)$group

# Take the subset (all controls & 1 condition /batch specified by pathway_regex)
pathway_condition_ind <- grep(paste(pathway_regex, collapse="|"), group)
ctrl_ind <- which(group %in% c("gfp_for_egfr", "gfp18", "gfp30"))
subset_ind <- sort(c(ctrl_ind, pathway_condition_ind))  
cts_sub <- cts_mat[, subset_ind]
batch_sub <- batch[subset_ind]
group_sub <- group[subset_ind]

message("Printing table to find non-gfp condition per batch.")
table(group_sub, batch_sub)
# remove genes with only 0 counts in the subset & in any batch
keep1 <- apply(cts_sub[, batch_sub==1],1,function(x){!all(x==0)})
keep2 <- apply(cts_sub[, batch_sub==2],1,function(x){!all(x==0)})
keep3 <- apply(cts_sub[, batch_sub==3],1,function(x){!all(x==0)})
cts_sub <- cts_sub[keep1 & keep2 & keep3, ]


## Use ComBatSeq to adjust data
group_sub <- factor(as.character(group_sub), levels=c("gfp_for_egfr", "gfp18", "gfp30",  gsub("^", "", pathway_regex, fixed=T)))
group_sub <- plyr::revalue(group_sub, c("gfp_for_egfr"="gfp", "gfp18"="gfp", "gfp30"="gfp"))

start_time <- Sys.time()
combatseq_sub <- ComBat_seq(counts=cts_sub, batch=batch_sub, group=group_sub, shrink=FALSE)
end_time <- Sys.time()
print(end_time - start_time)

message("After ComBatSeq adjusts data")
table(group_sub, batch_sub)


## Use original ComBat on logCPM
combat_sub <- sva::ComBat(cpm(cts_sub, log=TRUE), batch=batch_sub, mod=model.matrix(~group_sub))


## RUVseq
group1 <- plyr::revalue(as.factor(as.character(group_sub[batch_sub==1])), c("gfp"="0", "her2"="1"))
deres1 <- edgeR_DEpipe(cts_sub[, batch_sub==1], batch=NULL, group=group1,
                       include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)
group2 <- plyr::revalue(as.factor(as.character(group_sub[batch_sub==2])), c("gfp"="0", "egfr"="1"))
deres2 <- edgeR_DEpipe(cts_sub[, batch_sub==2], batch=NULL, group=group2,
                       include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)
group3 <- plyr::revalue(as.factor(as.character(group_sub[batch_sub==3])), c("gfp"="0", "kraswt"="1"))
deres3 <- edgeR_DEpipe(cts_sub[, batch_sub==3], batch=NULL, group=group3,
                       include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)
null_genes <- Reduce(intersect, list(which(deres1$de_res$FDR>0.95), which(deres2$de_res$FDR>0.95), which(deres3$de_res$FDR>0.95)))
group_obj <- makeGroups(group_sub)
ruvseq_sub <- RUVs(cts_sub, cIdx=null_genes, scIdx=group_obj, k=1)$normalizedCounts


## Normalize library size
cts_norm <- apply(cts_sub, 2, function(x){x/sum(x)})
cts_adj_norm <- apply(combatseq_sub, 2, function(x){x/sum(x)})
cts_adjori_norm <- apply(combat_sub, 2, function(x){x/sum(x)})
cts_ruvseq_norm <- apply(ruvseq_sub, 2, function(x){x/sum(x)})


## PCA 
col_data <- data.frame(Batch=factor(batch_sub), Group=group_sub) 
rownames(col_data) <- colnames(cts_sub)

col_data_obj <- col_data[colnames(cts_norm), , drop=FALSE]
seobj <- SummarizedExperiment(assays=cts_norm, colData=col_data_obj)
pca_obj <- plotPCA(DESeqTransform(seobj), intgroup=c("Batch", "Group"))
plt <- ggplot(pca_obj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj$plot_env$percentVar[2])),
       title="Unadjusted") 

col_data_adj <- col_data[colnames(cts_adj_norm), , drop=FALSE]
seobj_adj <- SummarizedExperiment(assays=cts_adj_norm, colData=col_data_adj)
pca_obj_adj <- plotPCA(DESeqTransform(seobj_adj), intgroup=c("Batch", "Group"))
plt_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="ComBat-Seq") 

col_data_adjori <- col_data[colnames(cts_adjori_norm), , drop=FALSE]
seobj_adjori <- SummarizedExperiment(assays=cts_adjori_norm, colData=col_data_adjori)
pca_obj_adjori <- plotPCA(DESeqTransform(seobj_adjori), intgroup=c("Batch", "Group"))
plt_adjori <- ggplot(pca_obj_adjori$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) + 
  geom_point() + 
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[2])),
       title="Original ComBat") 

## ruvseq
col_data_ruv <- col_data[colnames(cts_ruvseq_norm), , drop=FALSE]
seobj_ruvseq <- SummarizedExperiment(assays=cts_ruvseq_norm, colData=col_data_ruv)
pca_obj_ruvseq <- plotPCA(DESeqTransform(seobj_ruvseq), intgroup=c("Batch", "Group"))
plt_ruvseq <- ggplot(pca_obj_ruvseq$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) +
  geom_point() +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_ruvseq$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_ruvseq$plot_env$percentVar[2])),
       title="RUV-Seq")

plt_PCA_full <- ggarrange(plt, plt_ruvseq, plt_adjori, plt_adj, ncol=1, nrow=4, common.legend=TRUE, legend="right")

## Create SummarizedExperiment objects for BatchQC compatibility
create_se_for_batchqc <- function(data_matrix, col_data) {
  SummarizedExperiment(assays = list(counts = data_matrix), colData = col_data)
}

# Create SE objects for each dataset with proper column matching
se_unadjusted <- create_se_for_batchqc(cpm(cts_sub, log=TRUE), col_data[colnames(cts_sub), , drop=FALSE])
se_combatseq <- create_se_for_batchqc(cpm(combatseq_sub, log=TRUE), col_data[colnames(combatseq_sub), , drop=FALSE])
se_combat <- create_se_for_batchqc(combat_sub, col_data[colnames(combat_sub), , drop=FALSE])
se_ruvseq <- create_se_for_batchqc(cpm(ruvseq_sub, log=TRUE), col_data[colnames(ruvseq_sub), , drop=FALSE])

# Function to safely get explained variation with error handling
safe_batchqc_explained_variation <- function(se, batch_col, condition_col, assay_name, method_name) {
  tryCatch({
    result <- batchqc_explained_variation(se, batch=batch_col, condition=condition_col, assay_name=assay_name)
    message("Result")
    str(result)
    # Check if result is valid and has explained_variation component
    if (is.null(result$EV_table_ind) || is.null(result$EV_table_ind$Explained) || length(result$EV_table_ind$Explained) == 0) {
      warning(paste("batchqc_explained_variation returned empty results for", method_name))
      # Return a dummy result with appropriate structure
      return(list(explained_variation = matrix(0, nrow=1, ncol=3, 
                                             dimnames=list(NULL, c("Full (Condition+Batch)", "Condition", "Batch")))))
    }
    return(result)
  }, error = function(e) {
    warning(paste("Error in batchqc_explained_variation for", method_name, ":", e$message))
    # Return a dummy result with appropriate structure  
    return(list(explained_variation = matrix(0, nrow=1, ncol=3,
                                           dimnames=list(NULL, c("Full (Condition+Batch)", "Condition", "Batch")))))
  })
}

# Updated function calls with error handling
varexp_full <- list(
  unadjusted = safe_batchqc_explained_variation(se_unadjusted, "Batch", "Group", "counts", "unadjusted")$EV_table_ind$Explained,
  combatseq = safe_batchqc_explained_variation(se_combatseq, "Batch", "Group", "counts", "combatseq")$EV_table_ind$Explained,
  combat = safe_batchqc_explained_variation(se_combat, "Batch", "Group", "counts", "combat")$EV_table_ind$Explained,
  ruvseq = safe_batchqc_explained_variation(se_ruvseq, "Batch", "Group", "counts", "ruvseq")$EV_table_ind$Explained
)

# Check if varexp_full has valid data before proceeding
if (all(sapply(varexp_full, function(x) all(x == 0)))) {
  warning("All batchqc_explained_variation results are empty or zero. Skipping explained variation plot.")
  # Create a simple dummy plot or skip this section
  plt_varexp_full <- ggplot() + 
    geom_text(aes(x=0.5, y=0.5, label="Explained variation data unavailable"), size=5) +
    theme_void() +
    labs(title="Explained Variation (Data Unavailable)")
} else {
  # Proceed with normal processing if data is available
  varexp_full_df <- reshape2::melt(varexp_full)
  varexp_full_df$L1 <- factor(varexp_full_df$L1, levels=c("unadjusted", "ruvseq", "combat","combatseq"))
  varexp_full_df$L1 <- plyr::revalue(varexp_full_df$L1, c("unadjusted"="Unadjusted", "combatseq"="ComBat-Seq",
                                                          "ruvseq"="RUV-Seq", "combat"="Original ComBat"))
  varexp_full_df$Var2 <- plyr::revalue(varexp_full_df$Var2, c("Full (Condition+Batch)"="Condition+Batch"))
  
  plt_varexp_full <- ggplot(varexp_full_df, aes(x=Var2, y=value)) +
    geom_boxplot() +
    facet_wrap(~L1, nrow=4, ncol=1) +
    labs(y="Explained variation") +
    theme(axis.title.x = element_blank())
}

# Final plot generation
ggarrange(plt_PCA_full, plt_varexp_full, ncol=2, widths=c(0.55, 0.45))
