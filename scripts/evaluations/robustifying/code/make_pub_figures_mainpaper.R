# Suppress all output and warnings
options(warn = -1)
suppressMessages(suppressWarnings({

rm(list = ls())
setwd("/scripts/evaluations/robustifying")
if (!dir.exists("/scripts/evaluations/robustifying/figures")) {
  dir.create("/scripts/evaluations/robustifying/figures")
}
sapply(c("ggplot2", "gridExtra", "reshape2", "DelayedMatrixStats", "plyr", "ggpubr", 
         "SummarizedExperiment", "MCMCpack", "dplyr", "scales"), require, character.only = TRUE, quietly = TRUE)

}))



########  Figure 1: simulation studies  ########
source("/scripts/evaluations/robustifying/code/helper.R")
results_dir <- "/scripts/evaluations/robustifying/results/"
file_lst <- grep(".csv", dir(results_dir), fixed = TRUE, value = TRUE)  # all results files
method_names <- c("lasso", "rf", "svm")
method_names_plt <- c("Lasso logistic regression", "Random Forest", "Support Vector Machines")
names(method_names_plt) <- method_names
Nbatch <- 3
N_sample_size <- 20
batch_mean_vec <- c(0, 3, 5)
batch_var_vec <- c(1, 3, 5)
perf_measures <- c("mxe", "auc")  
perf_measures_plt <- c("Mean cross-entropy loss", "AUC")  
names(perf_measures_plt) <- perf_measures
subset_colnames_single <- c("Batch", "ComBat", "GMM", "n_Avg", "CS_Avg", "Reg_s")
batch_levels <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)),
                           1, paste, collapse = "_"))
curr_files_mv_suffix <- sort(apply(expand.grid(paste0("m", batch_mean_vec), paste0("v", batch_var_vec)),
                                   1, paste, collapse = "_"))
curr_mod <- "rf"
plt_lst <- list()

for (curr_perf in perf_measures) {
  curr_file_lst <- sort(apply(expand.grid(method_names, paste0(curr_perf, "_batchN", N_sample_size, "_",
                                                               curr_files_mv_suffix, ".csv")),
                              1, paste, collapse = "_"))
  # print(curr_file_lst)  # Suppressed for cleaner output
  sgmod_res_lst <- base_lst <- list()
  
  for (bl in batch_levels) {
    curr_file_lst_bl <- grep(bl, curr_file_lst, value = TRUE)
    curr_res <- read.csv(paste0(results_dir, grep(curr_mod, curr_file_lst_bl, value = TRUE)), header = TRUE)
    tmp <- curr_res[, subset_colnames_single]
    tmp <- data.frame(quiet_melt(tmp),
                      Type = rep(c("With batch effect,\nno adjustment", "Merging", "Merging", rep("Ensemble", 3)),
                                 each = nrow(tmp)))
    sgmod_res_lst[[bl]] <- tmp
    base_lst[[bl]] <- curr_res[, "NoBatch"]
  }
  plt_df <- quiet_melt(sgmod_res_lst)

  plt_df$variable <- factor(plt_df$variable, levels=subset_colnames_single)
  plt_df$variable <- revalue(plt_df$variable, c("Batch" = "No adjustment",
                                                      "ComBat" = "Merge + ComBat",
                                                      "GMM" = "GMM adjust",
                                                      "n_Avg" = "Batch size",
                                                      "CS_Avg" = "Cross-study",
                                                      "Reg_s" = "Stacking regression"))
  plt_df$L1 <- revalue(plt_df$L1, c("m0_v1"="Mean difference 0\nVariance fold change 1",
                                    "m0_v3"="Mean difference 0\nVariance fold change 3",
                                    "m0_v5"="Mean difference 0\nVariance fold change 5",
                                    "m5_v1"="Mean difference 5\nVariance fold change 1",
                                    "m5_v3"="Mean difference 5\nVariance fold change 3",
                                    "m5_v5"="Mean difference 5\nVariance fold change 5"))
  plt_df$Type <- factor(plt_df$Type, levels=c("With batch effect,\nno adjustment", "Merging", "Ensemble"))

  plt_lst[[curr_perf]] <- ggplot(plt_df, aes(x=variable, y=value, fill=Type)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=median(do.call(c, base_lst)),
                   linetype="Average AUC on data\nwithout simulated batch\neffect"), color="red") +
    scale_linetype_manual(name="limit", values=2,
                          guide=guide_legend(override.aes=list(color=c("red")))) +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
    facet_wrap(~L1, ncol=length(batch_var_vec)) +
    theme(axis.text.x=element_text(angle=30, hjust=1),
          axis.title.x=element_blank(),
          legend.title=element_blank(),
          legend.position="top",
          legend.direction="vertical",
          plot.margin = margin(0.1, 0.1, 0.1, 0.4, "cm")) +
    labs(y=perf_measures_plt[curr_perf])
}

png("/scripts/evaluations/robustifying/figures/Fig1.png", width=5, height=8, units="in", res=300)
plt_lst[["auc"]]
invisible(dev.off())

########  Figure 1b: simulation studies - Mean shift 5 only  ########
# Create a version of Figure 1 showing only mean shift 5 conditions
plt_lst_mean5 <- list()

for (curr_perf in perf_measures) {
  curr_file_lst <- sort(apply(expand.grid(method_names, paste0(curr_perf, "_batchN", N_sample_size, "_",
                                                               curr_files_mv_suffix, ".csv")),
                              1, paste, collapse = "_"))
  
  sgmod_res_lst <- base_lst <- list()
  
  # Filter to only mean shift 5 conditions
  batch_levels_mean5 <- batch_levels[grepl("^m5_", batch_levels)]
  
  for (bl in batch_levels_mean5) {
    curr_file_lst_bl <- grep(bl, curr_file_lst, value = TRUE)
    curr_res <- read.csv(paste0(results_dir, grep(curr_mod, curr_file_lst_bl, value = TRUE)), header = TRUE)
    tmp <- curr_res[, subset_colnames_single]
    tmp <- data.frame(quiet_melt(tmp),
                      Type = rep(c("With batch effect,\nno adjustment", "Merging", "Merging", rep("Ensemble", 3)),
                                 each = nrow(tmp)))
    sgmod_res_lst[[bl]] <- tmp
    base_lst[[bl]] <- curr_res[, "NoBatch"]
  }
  plt_df <- quiet_melt(sgmod_res_lst)

  plt_df$variable <- factor(plt_df$variable, levels=subset_colnames_single)
  plt_df$variable <- revalue(plt_df$variable, c("Batch" = "No adjustment",
                                                      "ComBat" = "Merge + ComBat",
                                                      "GMM" = "GMM adjust",
                                                      "n_Avg" = "Batch size",
                                                      "CS_Avg" = "Cross-study",
                                                      "Reg_s" = "Stacking regression"))
  plt_df$L1 <- revalue(plt_df$L1, c("m5_v1"="Mean difference 5\nVariance fold change 1",
                                    "m5_v3"="Mean difference 5\nVariance fold change 3",
                                    "m5_v5"="Mean difference 5\nVariance fold change 5"))
  plt_df$Type <- factor(plt_df$Type, levels=c("With batch effect,\nno adjustment", "Merging", "Ensemble"))

  plt_lst_mean5[[curr_perf]] <- ggplot(plt_df, aes(x=variable, y=value, fill=Type)) +
    geom_boxplot() +
    geom_hline(aes(yintercept=median(do.call(c, base_lst)),
                   linetype="Average AUC on data\nwithout simulated batch\neffect"), color="red") +
    scale_linetype_manual(name="limit", values=2,
                          guide=guide_legend(override.aes=list(color=c("red")))) +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
    facet_wrap(~L1, ncol=length(batch_var_vec)) +
    theme(axis.text.x=element_text(angle=30, hjust=1),
          axis.title.x=element_blank(),
          legend.title=element_blank(),
          legend.position="top",
          legend.direction="vertical",
          plot.margin = margin(0.1, 0.1, 0.1, 0.4, "cm")) +
    labs(y=perf_measures_plt[curr_perf])
}

png("/scripts/evaluations/robustifying/figures/Fig1b_mean5only.png", width=8, height=5, units="in", res=300)
plt_lst_mean5[["auc"]]
invisible(dev.off())




########  Figure 2: real data application  ########
rm(list=ls())
source("/scripts/evaluations/robustifying/code/helper.R")

#### 4 studies
results_dir <- "/scripts/evaluations/robustifying/results_real_4studies/"
cat("=== DEBUG: Fig2 4-study section ===\n")
cat("DEBUG: Using results directory:", results_dir, "\n")
cat("DEBUG: Directory exists:", dir.exists(results_dir), "\n")
study_names <- c("GSE37250_SA", "GSE37250_M", "US", "India")
study_label <- c("F", "G", "E", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc")  
perf_measures_names <- c("Mean cross-entropy loss", "AUC")  
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "GMM", "n_Avg", "CS_Avg", "Reg_s")
MODEL <- "crossmod"

## calculate frequency table
#j=1; b=1
curr_perf <- "mxe"
freq_lst <- list()
for(j in 1:length(study_names)){
  curr_testset <- study_names[j]
  file_prefix <- sprintf("test%s", curr_testset)
  res <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
  colnames(res) <- c("Method", "value", "Model", "Iteration")
  res$Study <- curr_testset
  
  # Get unique methods and create frequency table
  methods_of_interest <- c("Batch", "ComBat", "GMM", "n_Avg", "CS_Avg", "Reg_s")
  total_count <- rep(0, length(methods_of_interest))
  names(total_count) <- methods_of_interest
  
  for(b in 1:100){
    curr_iter <- dplyr::filter(res, Iteration==b, Model=="crossmod", 
                               Method %in% methods_of_interest)
    
    # Handle duplicates by taking the first occurrence of each method
    curr_iter_unique <- curr_iter[!duplicated(curr_iter$Method), ]
    
    if(nrow(curr_iter_unique) == length(methods_of_interest)) {
      best_method <- curr_iter_unique$Method[which.min(curr_iter_unique$value)]
      total_count[best_method] <- total_count[best_method] + 1
    }
  }
  freq_lst[[study_names[j]]] <- total_count / 100
}
annot_tb <- do.call(rbind, freq_lst)
rownames(annot_tb) <- study_label


plt_4s <- res_lst <- list()
i=1;j=1
for(i in 1:length(perf_measures)){
  curr_perf <- perf_measures[i]
  res_lst[[curr_perf]] <- plt_4s[[curr_perf]] <- list()
  for(j in 1:length(study_names)){
    curr_testset <- study_names[j]
    file_prefix <- sprintf("test%s", curr_testset)
    res_lst[[curr_perf]][[curr_testset]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]][[curr_testset]]) <- c("Method", "value", "Model", "Iteration")
    
    filtered_data <- res_lst[[curr_perf]][[curr_testset]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "GMM", "n_Avg", "CS_Avg", "Reg_s"))
    cat("DEBUG: Methods after filtering for", curr_testset, curr_perf, ":", paste(unique(filtered_data$Method), collapse=", "), "\n")
    sumstats <- filtered_data %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025), .groups = "drop")
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method %in% c("ComBat", "GMM")] <- "Merging"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "Merging", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "GMM", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "GMM"="GMM adjust", "n_Avg"="Batch size", "CS_Avg"="Cross-study",
                                                        "Reg_s"="Stacking regression"))
    
    sumstats_crossmod <- dplyr::filter(sumstats, Model==MODEL)
    sumstats_crossmod$Annot <- ""
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Original data"] <- scales::percent(annot_tb[study_label[curr_testset], "Batch"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Merge + ComBat"] <- scales::percent(annot_tb[study_label[curr_testset], "ComBat"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="GMM adjust"] <- scales::percent(annot_tb[study_label[curr_testset], "GMM"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Batch size"] <- scales::percent(annot_tb[study_label[curr_testset], "n_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Cross-study"] <- scales::percent(annot_tb[study_label[curr_testset], "CS_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Stacking regression"] <- scales::percent(annot_tb[study_label[curr_testset], "Reg_s"])
    
    plt_4s[[curr_perf]][[curr_testset]] <- ggplot(sumstats_crossmod, aes(x=Method, y=Avg, color=Type)) +
      geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
      geom_text(aes(label=Annot, y=0.6), color="black", size=3.5) + 
      geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
      geom_point() +
      labs(y=perf_measures_names[perf_measures[i]],
           title=ifelse(curr_perf=="mxe", paste("Test set:", study_label[curr_testset]), "")) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=30, hjust=1),
            legend.title=element_blank(),
            panel.grid.major=element_blank())
    
    res_lst[[curr_perf]][[curr_testset]]$Study <- curr_testset
   }
}

#### 6 studies
results_dir <- "/scripts/evaluations/robustifying/results_real_6studies/"
cat("=== DEBUG: Fig2 6-study section ===\n")
cat("DEBUG: Using results directory:", results_dir, "\n")
cat("DEBUG: Directory exists:", dir.exists(results_dir), "\n")
study_names <- c("GSE37250_SA", "GSE37250_M", "GSE39941_M", "US", "Africa", "India")
study_label <- c("F", "G", "C", "E", "A", "D")
names(study_label) <- study_names
norm_data <- TRUE 
use_ref_combat <- FALSE 
match_preval <- FALSE
perf_measures <- c("mxe", "auc", "rmse", "f", "err", "acc") 
perf_measures_names <- c("Mean cross-entropy loss", "AUC", "Root-mean-squared error", 
                         "F1 score", "Error rate", "Accuracy") 
names(perf_measures_names) <- perf_measures
selected_method <- c("Batch", "ComBat", "GMM", "n_Avg", "CS_Avg", "Reg_s")
MODEL <- "crossmod"

curr_perf <- "mxe"
freq_lst <- list()
for(j in 1:length(study_names)){
    curr_testset <- study_names[j]
    file_prefix <- sprintf("test%s", curr_testset)
    res <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
  colnames(res) <- c("Method", "value", "Model", "Iteration")
  res$Study <- curr_testset
    
  # Get unique methods and create frequency table
  methods_of_interest <- c("Batch", "ComBat", "GMM", "n_Avg", "CS_Avg", "Reg_s")
  total_count <- rep(0, length(methods_of_interest))
  names(total_count) <- methods_of_interest
  
  for(b in 1:100){
    curr_iter <- dplyr::filter(res, Iteration==b, Model=="crossmod", 
                               Method %in% methods_of_interest)
    
    # Handle duplicates by taking the first occurrence of each method
    curr_iter_unique <- curr_iter[!duplicated(curr_iter$Method), ]
    
    if(nrow(curr_iter_unique) == length(methods_of_interest)) {
      best_method <- curr_iter_unique$Method[which.min(curr_iter_unique$value)]
      total_count[best_method] <- total_count[best_method] + 1
    }
  }
  freq_lst[[study_names[j]]] <- total_count / 100
}
annot_tb <- do.call(rbind, freq_lst)
rownames(annot_tb) <- study_label

plt_6s <- res_lst <- list()
for(i in 1:length(perf_measures)){
  curr_perf <- perf_measures[i]
  res_lst[[curr_perf]] <- plt_6s[[curr_perf]] <- list()
  for(j in 1:length(study_names)){
    curr_testset <- study_names[j]
    file_prefix <- sprintf("test%s", curr_testset)
    res_lst[[curr_perf]][[curr_testset]] <- read.csv(paste0(results_dir, "/", file_prefix, "_", curr_perf, ".csv"))
    colnames(res_lst[[curr_perf]][[curr_testset]]) <- c("Method", "value", "Model", "Iteration")
    
    filtered_data <- res_lst[[curr_perf]][[curr_testset]] %>%
      dplyr::filter(Method %in% c("Batch", "ComBat", "GMM", "n_Avg", "CS_Avg", "Reg_s"))
    cat("DEBUG: Methods after filtering for", curr_testset, curr_perf, ":", paste(unique(filtered_data$Method), collapse=", "), "\n")
    
    sumstats <- filtered_data %>%
      dplyr::group_by(Method, Model) %>%
      dplyr::summarise(Avg=mean(value), Up=quantile(value, 0.975), Down=quantile(value, 0.025), .groups = 'drop')
    
    sumstats$Type <- rep("Original Data", nrow(sumstats))
    sumstats$Type[sumstats$Method %in% c("ComBat", "GMM")] <- "Merging"
    sumstats$Type[sumstats$Method %in% c("n_Avg","CS_Avg","Reg_s")] <- "Ensemble"
    sumstats$Type <- factor(sumstats$Type, levels=c("Original Data", "Merging", "Ensemble"))
    sumstats$Method <- factor(sumstats$Method, levels=c("Batch", "ComBat", "GMM", "n_Avg", "CS_Avg", "Reg_s"))
    sumstats$Method <- plyr::revalue(sumstats$Method, c("Batch"="Original data", "ComBat"="Merge + ComBat",
                                                        "GMM"="GMM adjust", "n_Avg"="Batch size",
                                                        "CS_Avg"="Cross-study", "Reg_s"="Stacking regression"))
    
    sumstats_crossmod <- dplyr::filter(sumstats, Model==MODEL)
    sumstats_crossmod$Annot <- ""
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Original data"] <- scales::percent(annot_tb[study_label[curr_testset], "Batch"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Merge + ComBat"] <- scales::percent(annot_tb[study_label[curr_testset], "ComBat"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="GMM adjust"] <- scales::percent(annot_tb[study_label[curr_testset], "GMM"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Batch size"] <- scales::percent(annot_tb[study_label[curr_testset], "n_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Cross-study"] <- scales::percent(annot_tb[study_label[curr_testset], "CS_Avg"])
    sumstats_crossmod$Annot[sumstats_crossmod$Method=="Stacking regression"] <- scales::percent(annot_tb[study_label[curr_testset], "Reg_s"])
    
    plt_6s[[curr_perf]][[curr_testset]] <- ggplot(sumstats_crossmod, aes(x=Method, y=Avg, color=Type)) +
      geom_errorbar(aes(ymin=Down, ymax=Up), width=.1) +
      geom_text(aes(label=Annot, y=1), color="black", size=3.5) + 
      geom_line(aes(x=Method, y=Avg, group=1), color="grey") +
      geom_point() +
      labs(y=perf_measures_names[perf_measures[i]],
           title=ifelse(curr_perf=="mxe", paste("Test set:", study_label[curr_testset]), "")) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_text(angle=30, hjust=1),
            legend.title=element_blank(),
            panel.grid.major=element_blank())
    
    res_lst[[curr_perf]][[curr_testset]]$Study <- curr_testset
  }
}

lgd <- get_legend(plt_4s[["mxe"]][["GSE37250_SA"]])
plt_4studies <- ggarrange(plt_4s[["mxe"]][["India"]], plt_4s[["mxe"]][["US"]], 
                          plt_4s[["mxe"]][["GSE37250_SA"]], plt_4s[["mxe"]][["GSE37250_M"]], 
                          ncol=4, nrow=1, common.legend=TRUE, legend="none")
plt_4studies <- annotate_figure(plt_4studies, 
                                top=text_grob("4 studies", face="bold", size=18, hjust=1, x=0.15))
plt_6studies <- ggarrange(plt_6s[["mxe"]][["India"]], plt_6s[["mxe"]][["US"]], 
                          plt_6s[["mxe"]][["GSE37250_SA"]], plt_6s[["mxe"]][["GSE37250_M"]], 
                          plt_6s[["mxe"]][["Africa"]], plt_6s[["mxe"]][["GSE39941_M"]],
                          ncol=6, nrow=1, common.legend=TRUE, legend="none")
plt_6studies <- annotate_figure(plt_6studies, 
                                top=text_grob("6 studies", face="bold", size=18, hjust=1, x=0.1))
p <- arrangeGrob(plt_4studies, plt_6studies, lgd,                              
                 ncol=3, nrow=3, layout_matrix=rbind(c(2,2,2), c(NA,NA,NA), c(1,1,3)),
                 heights=c(1,0.1,1))

png("/scripts/evaluations/robustifying/figures/Fig2.png", width=13, height=6, units="in", res=300)
print(as_ggplot(p))
invisible(dev.off())

