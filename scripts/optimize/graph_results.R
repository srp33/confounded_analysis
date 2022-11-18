library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

in_file_path = commandArgs(trailingOnly = TRUE)[1]

data = read_tsv(in_file_path)

plot_data = select(data, minibatch_size:learning_rate, combined_rank) %>%
  t() %>% # We do these steps so that all columns are the same type so we don't get an error when we pivot.
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  pivot_longer(minibatch_size:learning_rate, names_to = "parameter", values_to = "parameter_value") %>%
  select(parameter, parameter_value, combined_rank) %>%
  mutate(combined_rank = as.double(combined_rank))

set.seed(0)

# Prevents Rplots.pdf from being created
pdf(NULL)

#We are getting a permissions issue with Rplots.pdf when trying to run this code...
ggplot(plot_data, aes(x = parameter_value, y = combined_rank)) +
  geom_boxplot() +
  facet_wrap(vars(parameter), scales = "free_x") +
  xlab("Parameter") +
  ylab("Combined rank (higher is better)") +
  theme_bw()

ggsave("/outputs/optimizations/gse20194/summarized_results.pdf", height=6, width=8)
