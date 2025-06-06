# Setup ------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(vroom)
  library(tibble) 
})


# Load data -------
IN_DIR = "/outputs/metrics/"
FIG_DIR = "/outputs/figures/"
SAMPLE_DIR = "/../data/"

cbp2 <- c("#E69F00", "#56B4E9","#009E73","#F0E442", 
          "#0072B2", "#D55E00", "#CC79A7", "#999999")

all_info_df <- tribble(
  ~dataset,    ~title,                            ~batch_label,      ~true_label,
  "gse20194", "GSE 20194 ER",                   "meta_batch",      "meta_er_status",
  "gse20194", "GSE 20194 HER2",                 "meta_batch",      "meta_her2_status",
  "gse20194", "GSE 20194 PR",                   "meta_batch",      "meta_pr_status",
  "gse24080",  "GSE 24080 Eventfree Survival",  "meta_batch",      "meta_efs_outcome_label",
  "gse24080",  "GSE 24080 Overall Survival",    "meta_batch",      "meta_os_outcome_label",
  "gse49711",  "GSE 49711 Stage",               "meta_Class",      "meta_INSS_Stage"
)

order <- c("unadjusted", "min_mean", "limma", "limma_target", "combat", "combat_target", "tampor")

score_functions <- c("roc_auc_score", "mutual_info_score", "accuracy_score")


#---Function to determine what the random chance that the largest label would be chosen --- 
calculate_random_accuracy <- function(df, column_name) {
  values <- df[[column_name]]
  table(values) %>%
    max() / length(values)
}


#---Function to calculate the baseline for a given column and score function---
baseline_for_column <- function(df, column_name, score_function) {
  if (score_function == "accuracy_score") {
    # For accuracy score, we can use the random accuracy as the baseline
    return(calculate_random_accuracy(df, column_name))
  } else if (score_function == "roc_auc_score") {
    # For ROC AUC score, we use 50% as the baseline
    return(0.5)
  } else if (score_function == "mutual_info_score") {
    # For mutual information score, 0 is the baseline, but the maximum is the entropy of the column
    return(calculate_entropy(df, column_name))
  } else {
    stop("Unknown score function: ", score_function)
  }
}

#---Function to get random accuracy for a given dataset and column name, given a cache of accuracies---
get_baseline <- function(file_cache, dataset_name, column_name, accuracy_cache, score_function = "accuracy_score") {
  key <- paste(dataset_name, column_name, score_function, sep = "_")
  if (!key %in% names(accuracy_cache)) {
    message(paste("Calculating random accuracy for dataset:", dataset_name, "and column:", column_name))
    value <- baseline_for_column(file_cache[[dataset_name]], column_name, score_function)
    accuracy_cache[[key]] <- value
  }
  return(accuracy_cache[[key]])
}

#---Function for debugging, prints whether a column exists, its unique values, and whether certain values appear---
check_column <- function(data, column_name, values_to_check) {
  if (!column_name %in% colnames(data)) {
    stop("Column '", column_name, "' does not exist in the data.")
  }
  unique_values <- unique(data[[column_name]])
  message("Unique values in '", column_name, "' column: ", paste(unique_values, collapse = ", "))
  if (length(values_to_check) == 0) {
    message("No values to check provided.")
    return()
  }
  message("Checking for values: ", paste(values_to_check, collapse = ", "))
  values_missing <- values_to_check[!values_to_check %in% unique_values]
  if (length(values_missing) > 0) {
    message("The following values are missing in '", column_name, "': ", paste(values_missing, collapse = ", "))
  }
  # For each value to check, print how many times it appears
  for (value in values_to_check) {
    count <- sum(data[[column_name]] == value, na.rm = TRUE)
    message("Value '", value, "' appears ", count, " times in '", column_name, "'.")
  }
}

batchall <- vroom(file.path(IN_DIR, "batch_classification.csv"))
trueall  <- vroom(file.path(IN_DIR, "true_classification.csv"))
problems(batchall)
problems(trueall)

# Debugging checks
check_column(batchall, "dataset", all_info_df$dataset)
check_column(trueall, "dataset", all_info_df$dataset)

file_cache <- list() # Key: dataset_name, Value: the loaded data frame/tibble
random_accuracies_cache <- list() # Key: "{dataset_name}_{column_name}_{score_function}", Value: the random accuracy value

pdf(NULL) # Suppress Rplots.pdf generation

for (i in seq_len(nrow(all_info_df))) {
  dataset_name  <- all_info_df$dataset[i]
  title <- all_info_df$title[i]
  batch_col_name <- all_info_df$batch_label[i]
  true_col_name  <- all_info_df$true_label[i]

  if (!dataset_name %in% names(file_cache)) {
    message(paste("Loading data for dataset:", dataset_name))
    current_data <- vroom(file.path(SAMPLE_DIR, dataset_name, "unadjusted.csv"))
    file_cache[[dataset_name]] <- current_data
  }
  current_data <- file_cache[[dataset_name]]

  check_column(trueall, "column", true_col_name)

  batch <- batchall %>% 
    filter(dataset == dataset_name)
  true <- trueall %>% 
    filter(dataset == dataset_name) %>%
    filter(column == true_col_name)

  together <- rbind(batch, true)
  check_column(together, "metric", c())
  check_column(together, "score", score_functions)

  for(score_function in score_functions) {
    message(paste("Creating figure for dataset:", dataset_name, "with score function:", score_function))
    together_score <- filter(together, score == score_function)

    ran_true <- get_baseline(file_cache, dataset_name, true_col_name, random_accuracies_cache, score_function)
    ran_batch <- get_baseline(file_cache, dataset_name, batch_col_name, random_accuracies_cache, score_function)

    # Save figures ------------------

    ggplot() +
      geom_boxplot(data = together_score, mapping = aes(x = factor(adjuster, order), y = value, color = column)) + 
      geom_jitter(data = together_score, mapping = aes(x = factor(adjuster, order), y = value, color = column), position=position_jitterdodge()) +
      geom_hline(yintercept = ran_true, color = "#56B4E9") + 
      geom_hline(yintercept = ran_batch, color = "#E69F00") + 
      ggtitle(title) +
      theme_bw(base_size = 18) + theme(axis.title.x=element_blank(), legend.title=element_blank(), axis.text.x = element_text(angle = 90)) + 
      scale_y_continuous(name = score_function, limits = c(0.0, 1.0)) +
      facet_wrap(vars(metric), strip.position = "top") + 
      scale_colour_manual(values=cbp2)
    filename_for_plot <- paste(c(FIG_DIR, title, "_", score_function, ".pdf"), collapse = "")
    ggsave(file.path(FIG_DIR, filename_for_plot), width = 11, height = 8.5, units = 'in')
  }
}

print(batchall)

#metric comparison --------------
metriccomp <- read_csv(paste(c(IN_DIR, "singlemetriccomparison_minus.csv"), collapse = ""))
metricdata <- filter(metriccomp, dataset %in% datasets)

score_averages = group_by(metricdata, adjuster, metric, dataset) %>%
  summarize(ave=median(score), 
	    interact = interaction(adjuster, metric))
score_averages <- group_by(score_averages, adjuster, metric)

# Save figure --------------------
ggplot() +
  geom_jitter(data = metricdata, mapping = aes(x = factor(dataset, datasets), y = score,  color = adjuster), position=position_jitterdodge()) + 
  geom_boxplot(data = metricdata, mapping = aes(x = factor(dataset, datasets), y = score, color = adjuster)) +
  
  ggtitle("Single Metric Comparison")+
  facet_wrap(vars(metric))+
  theme_bw(base_size = 12) + 
  scale_y_continuous(name = "Score: (true - trueRandom) - abs(batchRandom - batch)", limits = c(-1.0, 1.0)) +
  theme(axis.title.x=element_blank(), legend.title=element_blank(), axis.text.x = element_text(angle = 90)) +
  scale_colour_manual(values=cbp2)
ggsave(paste(c(FIG_DIR, "singlemetriccomparison_minus.pdf"), collapse = ""), width = 11, height = 8.5, units = 'in')

ggplot() +
	geom_line(data = score_averages, mapping = aes(x = factor(dataset, datasets), y = ave, color = adjuster, linetype = metric, group = interact))+
	geom_point(data = score_averages, mapping = aes(x = factor(dataset, datasets), y = ave, color = adjuster, pch = metric)) +
	ggtitle("Single Metric Score Averages")+
  theme_bw(base_size = 12) +
  scale_y_continuous(name = "Average Score: (true - trueRandom) - abs(batchRandom - batch)", limits = c(-1.0, 1.0)) +
  theme(axis.title.x=element_blank(), legend.title=element_blank(), axis.text.x = element_text(angle = 90)) +
  scale_colour_manual(values=cbp2)
ggsave(paste(c(FIG_DIR, "singleMetricScoreAverages_minus.pdf"), collapse = ""), width = 11, height = 8.5, units = 'in')




# Classification accuracy --------
#df <- read_csv(paste0(IN_DIR, "/classification.csv")) %>%
#    mutate(adjuster = factor(adjuster, levels = order))


#df <- df %>% filter(model != "MLPClassifier")
#df <- df %>% filter(!str_detect(dataset, "pretrain"))

#df %>% filter(col_type == "batch_col") %>%
#  ggplot(aes(x = model, y = accuracy, fill = model)) +
#  geom_boxplot() +
#  facet_grid(dataset ~ adjuster, scales = "free_y") +
#  theme(axis.text.x = element_blank(),
#        axis.ticks.x = element_blank(),
#        axis.title.x = element_blank()) +
#  labs(y = "Batch Classification Accuracy")
#ggsave(paste0(FIG_DIR, "/batch_accuracy.pdf"))

#df %>% filter(col_type == "true_class_col") %>%
#  ggplot(aes(x = model, y = accuracy, fill = model)) +
#  geom_boxplot() +
#  facet_grid(dataset ~ adjuster, scales = "free_y") +
#  labs(y = "True Class Classification Accuracy") +
#  theme_bw(base_size = 12) +
#  theme(axis.text.x = element_blank(),
#        axis.ticks.x = element_blank(),
#        axis.title.x = element_blank())
#ggsave(paste0(FIG_DIR, "/true_class_accuracy.pdf"))


