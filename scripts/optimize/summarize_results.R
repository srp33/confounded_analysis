library(dplyr)
library(readr)
library(tidyr)

results_file_path = commandArgs(trailingOnly = TRUE)[1]
unadjusted_name = commandArgs(trailingOnly = TRUE)[2]
out_file_path = commandArgs(trailingOnly = TRUE)[3]

data = read_tsv(results_file_path)

# Summarize results across iterations first, then across algorithms
unadjusted = filter(data, dataset == unadjusted_name) %>%
    select(algorithm, iteration, column, value) %>%
    group_by(algorithm, column) %>%
    summarize(unadjusted_value = mean(value)) %>%
    group_by(column) %>%
    summarize(unadjusted_value = mean(unadjusted_value)) %>%
    ungroup()

adjusted = filter(data, dataset != "unadjusted") %>%
    select(-dataset) %>%
    group_by(minibatch_size, ds_layers, ae_layers, code_size, scaling, loss_weight, minibatch_iterations, learning_rate, algorithm, column) %>%
    summarize(value = mean(value)) %>%
    group_by(minibatch_size, ds_layers, ae_layers, code_size, scaling, loss_weight, minibatch_iterations, learning_rate, column) %>%
    summarize(value = mean(value)) %>%
    ungroup()

# Higher batch scores are better
batch_data = inner_join(adjusted, unadjusted, by=c("column")) %>%
    filter(column == "Batch") %>%
    select(-column) %>%
    mutate(batch_score = unadjusted_value - value) %>%
    dplyr::rename(batch_value = value) %>%
    dplyr::rename(batch_unadjusted_value = unadjusted_value)

# Higher class scores are better
class_data = inner_join(adjusted, unadjusted, by=c("column")) %>%
    filter(column == "Class") %>%
    select(-column) %>%
    mutate(class_score = value - unadjusted_value) %>%
    dplyr::rename(class_value = value) %>%
    dplyr::rename(class_unadjusted_value = unadjusted_value)

combined_data = inner_join(batch_data, class_data)

# Convert numbers to ranks so they are easier to compare
# Higher ranks are better
combined_data %>%
    mutate(batch_rank = rank(batch_score) / nrow(combined_data)) %>%
    mutate(class_rank = rank(class_score) / nrow(combined_data)) %>%
    mutate(combined_rank = (batch_rank + class_rank) / 2) %>%
    arrange(desc(combined_rank)) %>%
    write_tsv(out_file_path)
