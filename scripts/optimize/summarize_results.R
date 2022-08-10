library(dplyr)
library(readr)
library(tidyr)

results_file_path = commandArgs(trailingOnly = TRUE)[1]

data = read_tsv(results_file_path)

unadjusted = filter(data, dataset == "unadjusted") %>%
    select(-dataset) %>%
    group_by(algorithm, column) %>%
    summarize(unadjusted_value = mean(value)) %>%
    ungroup()

adjusted = filter(data, dataset != "unadjusted") %>%
    select(-dataset) %>%
    group_by(across(!value)) %>%
    summarize(value = mean(value)) %>%
    ungroup()

batch_data = inner_join(adjusted, unadjusted, by=c("algorithm", "column")) %>%
    filter(column == "Batch") %>%
    select(-column) %>%
    mutate(batch_score = unadjusted_value - value) %>%
    dplyr::rename(batch_value = value) %>%
    dplyr::rename(batch_unadjusted_value = unadjusted_value) %>%
    arrange(desc(batch_score))

class_data = inner_join(adjusted, unadjusted, by=c("algorithm", "column")) %>%
    filter(column == "Class") %>%
    select(-column) %>%
    mutate(class_score = value - unadjusted_value) %>%
    dplyr::rename(class_value = value) %>%
    dplyr::rename(class_unadjusted_value = unadjusted_value) %>%
    arrange(desc(class_score))

# Higher ranks are better
data = inner_join(batch_data, class_data) %>%
    mutate(batch_rank = rank(batch_score)) %>%
    mutate(class_rank = rank(class_score)) %>%
    mutate(combined_rank = batch_rank + class_rank) %>%
#    filter(batch_unadjusted_value < 0.3 | batch_unadjusted_value > 0.7) %>%
    filter(batch_value > 0.4 & batch_value < 0.6) %>%
#    filter(class_value > 0.8) %>%
    filter(class_value > (class_unadjusted_value - 0.02)) %>%
    arrange(desc(combined_rank))

#print(batch_data, n=5, width=Inf)
#print(class_data, n=5, width=Inf)
print(data, n=20, width=Inf)

#minibatch_size  ds_layers   ae_layers   code_size   scaling loss_weight minibatch_iterations    learning_rate         dataset    algorithm   column

#print(unadjusted, n=Inf, width=Inf)

#                score = (float(line[3]) - true_random) - abs(batch_random - float(line[4]))
