# Setup ------
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load("tidyverse", "kableExtra")
source("functions.R")

# Load data -------
IN_DIR = "../data/metrics/"
FIG_DIR = "../data/output/"
TAB_DIR = "../data/output/"

mse <- read_csv(paste(c(IN_DIR, "/mse.csv"), collapse = ""))
mmd <- read_csv(paste(c(IN_DIR, "/mmd.csv"), collapse = ""))

order <- c("unadjusted", "scale", "combat", "confounded")

metrics <- rbind(mse, mmd) %>%
    filter(!str_detect(adjuster, "pretrained")) %>%
    mutate(adjuster = factor(adjuster, levels = order))

# This is a workaround for a ggplot2 bug.
pdf(NULL)

# Save figures ------------------
metrics %>%
  filter(metric == "MSE") %>%
  ggplot(aes(x = metric, y = value)) +
  geom_bar(stat = "identity") +
  facet_grid(dataset ~ adjuster, scales = "free_y", drop = TRUE) +
  labs(y = "MSE", x = "")
ggsave(paste(c(FIG_DIR, "/mse.pdf"), collapse = ""))

metrics %>%
  filter(metric == "MMD") %>%
  ggplot(aes(x = metric, y = value)) +
  geom_bar(stat = "identity") +
  facet_grid(dataset ~ adjuster, scales = "free_y", drop = TRUE) +
  labs(y = "MMD", x = "")
ggsave(paste(c(FIG_DIR, "/mmd.pdf"), collapse = ""))

# Save tables --------------
metrics_table <- function(metrics, metric_name) {
  return(
    metrics %>%
      mutate(value = formatC(value, format = "g", digits = 3)) %>%
      spread(adjuster, value) %>%
      filter(metric == metric_name) %>%
      select(-metric) %>%
      rename(Dataset = dataset)
  )
}

save_table <- function(table, path) {
  return(
    table %>%
      kable("latex") %>%
      kable_styling(latex_options = "striped") %>%
      as.character() %>%
      str_split("\n") %>%
      unlist() %>%
      head(-1) %>% tail(-1) %>% paste(collapse = "\n") %>%
      cat(file=path, sep="")
  )
}


metrics %>%
  metrics_table("MSE") %>%
  save_table(paste(c(TAB_DIR, "/mse.tex"), collapse = ""))
metrics %>%
  metrics_table("MMD") %>%
  save_table(paste(c(TAB_DIR, "/mmd.tex"), collapse = ""))

# Classification accuracy --------
df <- read_csv(paste0(IN_DIR, "/classification.csv")) %>%
    mutate(adjuster = factor(adjuster, levels = order))


df <- df %>% filter(model != "MLPClassifier")
df <- df %>% filter(!str_detect(dataset, "pretrain"))

df %>% filter(col_type == "batch_col") %>%
  ggplot(aes(x = model, y = accuracy, fill = model)) +
  geom_boxplot() +
  facet_grid(dataset ~ adjuster, scales = "free_y") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = "Batch Classification Accuracy")
ggsave(paste0(FIG_DIR, "/batch_accuracy.pdf"))

df %>% filter(col_type == "true_class_col") %>%
  ggplot(aes(x = model, y = accuracy, fill = model)) +
  geom_boxplot() +
  facet_grid(dataset ~ adjuster, scales = "free_y") +
  labs(y = "True Class Classification Accuracy") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
ggsave(paste0(FIG_DIR, "/true_class_accuracy.pdf"))

# Make the table

for (col in c("batch_col", "true_class_col")) {
  my_table <- df %>%
    mutate(model = str_replace_all(model, "Classifier", "")) %>%
    group_by(adjuster, model, dataset, col_type) %>%
    summarize(accuracy = round(mean(accuracy), 3), baseline = round(mean(baseline), 3)) %>%
    filter(col_type == col) %>%
    select(-col_type) %>%
    spread(model, accuracy) %>%
    arrange(dataset, adjuster) %>%
    select(dataset, adjuster, baseline, everything()) %>%
    rename(Dataset = dataset, Adjustment = adjuster, Baseline = baseline) %>%
    kable("latex") %>%
    kable_styling(latex_options = "striped") %>%
    collapse_rows(columns = 1:3, valign = "top")
#   my_table %>% print()
  # We need to remove the first and last line for it to work in the latex file...
  filename = paste(c(
    TAB_DIR,
    str_replace(col, "_col", ""),
    ".tex"
  ), collapse = "")
  my_table %>%
    as.character() %>%
    str_split("\n") %>%
    unlist() %>%
    head(-1) %>% tail(-1) %>% paste(collapse = "\n") %>%
    cat(file=filename, sep="")
}
