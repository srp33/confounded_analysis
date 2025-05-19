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

datasets = c("simulated_expression", "bladderbatch", "gse37199", 
                "tcga_small", "tcga_medium", "tcga"
                )
titles = c("Simulated Data", "Bladderbatch", "GSE 37199", 
                "TCGA Small", "TCGA Medium", "TCGA"
                )
batch_labels = c("Batch", "batch", "plate",
              "CancerType", "CancerType", "CancerType"
	      )
true_labels = c("Class", "batch", "Stage",
                "TP53_Mutated", "TP53_Mutated", "TP53_Mutated"
                )

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
#-- TCGA dataset comparisons
title = "TCGA Dataset Comparisons"
batchall <- batchall %>% mutate(valType = batch_labels[x])
trueall <- trueall %>% mutate(valType = true_labels[x])
tcga_sets <- c("tcga", "tcga_medium", "tcga_small")
tcga_batch <- filter(batchall, dataset %in% tcga_sets)
tcga_true <-filter(trueall, dataset %in% tcga_sets)

#-- find the average for each
batch_averages = group_by(tcga_batch, valType, dataset, adjuster, metric) %>%
  summarize(ave=median(value))
true_averages = group_by(tcga_true, valType, dataset, adjuster, metric) %>%
  summarize(ave=median(value))
averages = rbind(true_averages, batch_averages)
print(averages)
averages = group_by(averages, valType, adjuster, metric)
level_order <- c("tcga_small", "tcga_medium", "tcga")

# Save figures ------------------

ggplot() +
  geom_line(data = averages, mapping = aes(x = factor(dataset, level_order), y = ave, color = adjuster, linetype = valType, group = interaction(adjuster, valType)), size = 1.5) +
  geom_point(data = averages, mapping = aes(x = factor(dataset, level_order), y = ave, color = adjuster), size = 1.5) +
  geom_hline(yintercept = random_accuracy("tcga", "TP53_Mutated"), color = "#000000", linetype = "dotted") +
  geom_hline(yintercept = random_accuracy("tcga", "CancerType"), color = "#000000", linetype = "solid") +
  ggtitle(title) +
  theme_bw(base_size = 18) + theme(axis.title.x=element_blank(), legend.title=element_blank(), axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(name = "Average Accuracy", limits = c(0.0, 1.0)) +
  scale_x_discrete(labels = c("small", "medium", "full")) + 
  facet_wrap(vars(metric), strip.position = "top") +
  scale_colour_manual(values=cbp2)
ggsave(paste(c(FIG_DIR, "tcga_comparisons.pdf"), collapse = ""), width = 11, height = 8.5, units = 'in')


#metric comparison --------------
metriccomp <- read_csv(paste(c(IN_DIR, "singlemetriccomparison_minus.csv"), collapse = ""))
tcga_only <- filter(metriccomp, dataset %in% tcga_sets)
metricdata <- filter(metriccomp, dataset %in% c("bladderbatch", "gse37199", "tcga", "simulated_expression"))

score_averages = group_by(metricdata, adjuster, metric, dataset) %>%
  summarize(ave=median(score), 
	    interact = interaction(adjuster, metric))
score_averages <- group_by(score_averages, adjuster, metric)

# Save figure --------------------
ggplot() +
  geom_jitter(data = metricdata, mapping = aes(x = factor(dataset, c("bladderbatch", "gse37199","tcga", "simulated_expression")), y = score,  color = adjuster), position=position_jitterdodge()) + 
  geom_boxplot(data = metricdata, mapping = aes(x = factor(dataset, c("bladderbatch", "gse37199","tcga", "simulated_expression")), y = score, color = adjuster)) +
  
  ggtitle("Single Metric Comparison")+
  facet_wrap(vars(metric))+
  theme_bw(base_size = 12) + 
  scale_y_continuous(name = "Score: (true - trueRandom) - abs(batchRandom - batch)", limits = c(-1.0, 1.0)) +
  theme(axis.title.x=element_blank(), legend.title=element_blank(), axis.text.x = element_text(angle = 90)) +
  scale_colour_manual(values=cbp2)
ggsave(paste(c(FIG_DIR, "singlemetriccomparison_minus.pdf"), collapse = ""), width = 11, height = 8.5, units = 'in')


ggplot() +
  geom_jitter(data = tcga_only, mapping = aes(x = factor(dataset, c("tcga_small", "tcga_medium", "tcga")), y = score,  color = adjuster), position=position_jitterdodge()) +
  geom_boxplot(data = tcga_only, mapping = aes(x = factor(dataset, c("tcga_small", "tcga_medium", "tcga")), y = score, color = adjuster)) +

  ggtitle("Single Metric Comparison - TCGA")+
  facet_wrap(vars(metric))+
  theme_bw(base_size = 12) +
  scale_y_continuous(name = "Score: (true - trueRandom) - abs(batchRandom - batch)", limits = c(-1.0, 1.0)) +
  theme(axis.title.x=element_blank(), legend.title=element_blank(), axis.text.x = element_text(angle = 90)) +
  scale_colour_manual(values=cbp2)
ggsave(paste(c(FIG_DIR, "singlemetriccomparison_tcga_minus.pdf"), collapse = ""), width = 11, height = 8.5, units = 'in')

ggplot() +
	geom_line(data = score_averages, mapping = aes(x = factor(dataset, c("bladderbatch", "gse37199","tcga", "simulated_expression")), y = ave, color = adjuster, linetype = metric, group = interact))+
	geom_point(data = score_averages, mapping = aes(x = factor(dataset, c("bladderbatch", "gse37199","tcga", "simulated_expression")), y = ave, color = adjuster, pch = metric)) +
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


