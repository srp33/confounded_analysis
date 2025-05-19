# Setup ------
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)


# Load data -------
IN_DIR = "/outputs/metrics/"
FIG_DIR = "/outputs/figures/"
SAMPLE_DIR = "/../data/"

cbp2 <- c("#E69F00", "#56B4E9","#009E73","#F0E442", 
          "#0072B2", "#D55E00", "#CC79A7", "#999999")


batchall <- read_csv(paste(c(IN_DIR, "batch_classification.csv"), collapse = ""))
trueall <- read_csv(paste(c(IN_DIR, "true_classification.csv"), collapse = ""))

#---Function to determine what the random chance that the largest label would be chosen --- 
random_accuracy <- function(data, label) { 
  samples <- read_csv(paste(c(SAMPLE_DIR, "/", data, "/unadjusted.csv"), collapse = ""))
  samples <- select(samples, label)
  max = 0
  total = 0
  output <- table(samples)
  for (sample in output) {
    total = total + sample
    if(max < sample) { 
      max = sample
    }
  }
  return(max / total)
}

#                     dataset,    title,                            batch_label,  true_label
gse20194_er_info =  c("gse20194", "GSE 20194 ER",                   "batch",      "er_status")
gse20194_her2_info =c("gse20194", "GSE 20194 HER2",                 "batch",      "her2_status")
gse20194_pr_info =  c("gse20194", "GSE 20194 PR",                   "batch",      "pr_status")

gse24080_efs_info = c("gse24080",  "GSE 24080 Eventfree Survival",  "batch",      "efs_outcome_label")
gse24080_os_info =  c("gse24080",  "GSE 24080 Overall Survival",    "batch",      "os_outcome_label")

gse49711_stage_info=c("gse49711",  "GSE 49711 Stage",               "Class",      "INSS_Stage")

# Combine into a list
all_info <- list(
  gse20194_er_info,
  gse20194_her2_info,
  gse20194_pr_info,
  gse24080_efs_info,
  gse24080_os_info,
  gse49711_stage_info
)

# Extract components into separate vectors
datasets     <- sapply(all_info, function(x) x[1])
titles       <- sapply(all_info, function(x) x[2])
batch_labels <- sapply(all_info, function(x) x[3])
true_labels  <- sapply(all_info, function(x) x[4])

order <- c("unadjusted", "scaled", "combat", "confounded")

pdf(NULL)

for(x in 1:length(datasets)) {
    batchall <- batchall %>% mutate(valType = batch_labels[x])
    trueall <- trueall %>% mutate(valType = true_labels[x])
    data = datasets[x]
    title = titles[x]
    ran_batch <- random_accuracy(data, batch_labels[x])
    ran_true <- random_accuracy(data, true_labels[x])

    batchs <- filter(batchall, dataset == data)
    trues <- filter(trueall, dataset == data)
    together <- rbind(batchs, trues)

    # This is a workaround for a ggplot2 bug.
    pdf(NULL)

    # Save figures ------------------

    ggplot() +
      geom_boxplot(data = together, mapping = aes(x = factor(adjuster, order), y = value, color = valType)) + 
      geom_jitter(data = together, mapping = aes(x = factor(adjuster, order), y = value, color = valType), position=position_jitterdodge()) +
      geom_hline(yintercept = ran_true, color = "#56B4E9") + 
      geom_hline(yintercept = ran_batch, color = "#E69F00") + 
      ggtitle(title) +
      theme_bw(base_size = 18) + theme(axis.title.x=element_blank(), legend.title=element_blank(), axis.text.x = element_text(angle = 90)) + 
      scale_y_continuous(name = "Accuracy", limits = c(0.0, 1.0)) +
      facet_wrap(vars(metric), strip.position = "top") + 
      scale_colour_manual(values=cbp2)
    ggsave(paste(c(FIG_DIR, data, ".pdf"), collapse = ""), width = 11, height = 8.5, units = 'in')

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


