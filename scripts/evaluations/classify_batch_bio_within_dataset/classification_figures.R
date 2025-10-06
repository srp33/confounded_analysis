# Setup ------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(vroom)
  library(tibble) 
  library(argparse)
})

parser <- ArgumentParser()

# Load data -------
parser$add_argument("-i", "--in_dir", type="character",
              help="Directory containing input files")
parser$add_argument("-f", "--fig_dir", type="character",
              help="Directory where figures will be saved")
parser$add_argument("-s", "--sample_dir", type="character", 
              help="Directory containing sample data")
parser$add_argument("-b", "--batch_filename", type="character",
              help="Filename for batch data")
parser$add_argument("-t", "--true_filename", type="character",
              help="Filename for true labels data")
                                        
args <- parser$parse_args()

IN_DIR <- args$in_dir
FIG_DIR <- args$fig_dir
SAMPLE_DIR <- args$sample_dir
BATCH_FILENAME <- args$batch_filename
TRUE_FILENAME <- args$true_filename

cbp2 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

all_info_df <- tribble(
  ~dataset,    ~title,                            ~batch_label,      ~true_label,
  "gse20194", "GSE 20194 ER",                   "meta_batch",      "meta_er_status",
  "gse20194", "GSE 20194 HER2",                 "meta_batch",      "meta_her2_status",
  "gse20194", "GSE 20194 PR",                   "meta_batch",      "meta_pr_status",
  "gse24080",  "GSE 24080 Cytogenetic Abnormality",    "meta_batch",      "meta_cytogenetic_abnormality",
  "gse49711",  "GSE 49711 Stage 3 4",           "meta_Sex",      "meta_INSS_Stage_Split_3_4",

  "gse_20194_62944", "GSE 20194 62944 ER",       "meta_source",      "meta_er_status",

  "2_dims_no_bio_no_batch",   "2 dims no bio no batch",   "meta_batch", "meta_bio",
  "2_dims_no_bio_yes_batch",  "2 dims no bio yes batch",  "meta_batch", "meta_bio",
  "2_dims_yes_bio_no_batch",  "2 dims yes bio no batch",  "meta_batch", "meta_bio",
  "2_dims_yes_bio_yes_batch",  "2 dims yes bio yes batch",  "meta_batch", "meta_bio",

  "400_dims_no_bio_no_batch",   "400 dims no bio no batch",   "meta_batch", "meta_bio",
  "400_dims_no_bio_yes_batch",  "400 dims no bio yes batch",  "meta_batch", "meta_bio",
  "400_dims_yes_bio_no_batch",  "400 dims yes bio no batch",  "meta_batch", "meta_bio",
  "400_dims_yes_bio_yes_batch",  "400 dims yes bio yes batch",  "meta_batch", "meta_bio",

  "1000_dims_no_bio_no_batch",   "1000 dims no bio no batch",   "meta_batch", "meta_bio",
  "1000_dims_no_bio_yes_batch",  "1000 dims no bio yes batch",  "meta_batch", "meta_bio",
  "1000_dims_yes_bio_no_batch",  "1000 dims yes bio no batch",  "meta_batch", "meta_bio",
  "1000_dims_yes_bio_yes_batch",  "1000 dims yes bio yes batch",  "meta_batch", "meta_bio",

  "structured_synthetic", "Structured Synthetic", "meta_batch", "meta_bio"
)

order <- c("unadjusted", "seurat_scaling", "liger", "seurat_integration", "monotonic", "non_monotonic", "wasserstein", "autoclass", "min_mean", "combat_target", "mnn", "combat", "quantile", "npn", "simple") #, "icvae", "fair_adapt", "limma_target", "limma", "harmony", 

score_functions <- c("roc_auc_score", "mutual_info_score", "accuracy_score")


#---Function to determine what the random chance that the largest label would be chosen --- 
calculate_random_accuracy <- function(df, column_name) {
  values <- df[[column_name]]
  table(values) %>%
    max() / length(values)
}

#---Function to calculate the entropy of a column in a data frame ---
calculate_entropy <- function(df, column_name) {
  values <- df[[column_name]]
  value_counts <- table(values)
  probabilities <- value_counts / sum(value_counts)
  entropy <- -sum(probabilities * log2(probabilities + 1e-10)) # Adding a small constant to avoid log(0)
  return(entropy)
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
    return(0)
  }
}

#---Function to get random accuracy for a given dataset and column name, given a cache of accuracies---
get_baseline <- function(file_cache, dataset_name, column_name, accuracy_cache, score_function = "accuracy_score") {
  key <- paste(dataset_name, column_name, score_function, sep = "_")
  if (!key %in% names(accuracy_cache)) {
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

batchall <- vroom(file.path(IN_DIR, BATCH_FILENAME))
trueall  <- vroom(file.path(IN_DIR, TRUE_FILENAME))
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

  check_column(trueall, "p_column", true_col_name)

  batch <- batchall %>% 
    filter(dataset == dataset_name) %>%
    filter(adjuster %in% order)
  true <- trueall %>% 
    filter(dataset == dataset_name) %>%
    filter(p_column == true_col_name) %>%
    filter(adjuster %in% order)

  together <- rbind(batch, true)
  check_column(together, "classifier", c())
  check_column(together, "metric", score_functions)

  for(score_function in score_functions) {
    message(paste("Creating figure for dataset:", dataset_name, "with score function:", score_function))
    together_score <- filter(together, metric == score_function)

    ran_true <- get_baseline(file_cache, dataset_name, true_col_name, random_accuracies_cache, score_function)
    ran_batch <- get_baseline(file_cache, dataset_name, batch_col_name, random_accuracies_cache, score_function)

    # Save figures ------------------

    ggplot() +
      geom_boxplot(data = together_score, mapping = aes(x = factor(adjuster, order), y = value, color = p_column, shape=c_value)) + 
      geom_jitter(data = together_score, mapping = aes(x = factor(adjuster, order), y = value, color = p_column, shape=c_value), position=position_jitterdodge()) +
      geom_hline(yintercept = ran_true, color = "#56B4E9") + 
      geom_hline(yintercept = ran_batch, color = "#E69F00") + 
      ggtitle(title) +
      theme_bw(base_size = 18) + theme(axis.title.x=element_blank(), legend.title=element_blank(), axis.text.x = element_text(angle = 90)) + 
      scale_y_continuous(name = score_function, limits = c(0.0, 1.0)) +
      facet_wrap(vars(classifier), strip.position = "top") + 
      scale_colour_manual(values=cbp2)
    filename_for_plot <- paste(c(title, "_", score_function, ".pdf"), collapse = "")
    ggsave(file.path(FIG_DIR, filename_for_plot), width = 11, height = 8.5, units = 'in')
  }
}
