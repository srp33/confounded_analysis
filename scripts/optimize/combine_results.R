library(dplyr)
library(ggplot2)
library(readr)

args = commandArgs(trailingOnly = TRUE)

in_file_path1 = args[1]
in_file_path2 = args[2]
out_file_path_tsv = args[3]
out_file_path_pdf = args[4]

data1 = read_tsv(in_file_path1)
data2 = read_tsv(in_file_path2)

data = inner_join(data1, data2, by=c("minibatch_size", "ds_layers", "ae_layers", "code_size", "scaling", "loss_weight", "minibatch_iterations", "learning_rate")) %>%
    filter(!is.na(class_value.x)) %>%
    filter(!is.na(class_value.y)) %>%
    mutate(overall_rank = combined_rank.x + combined_rank.y) %>%
    arrange(desc(overall_rank))

write_tsv(data, out_file_path_tsv)
 
x = pull(data, combined_rank.x)
y = pull(data, combined_rank.y)
#print(shapiro.test(x))
#print(shapiro.test(y))
cor_coef = cor(x, y)
cor_p = cor.test(x, y, method="spearman")$p.value
title = paste0("rho = ", round(cor_coef, 2), " (p = ", round(cor_p, 2), ")")

# Prevents Rplots.pdf from being created
pdf(NULL)

ggplot(data, aes(x = combined_rank.x, y = combined_rank.y)) +
    geom_point() +
    ggtitle(title) +
    xlab("GSE20194 (microarray)") +
    ylab("GSE49711 (RNA-Seq)") +
    xlim(0, 0.5) +
    ylim(0, 0.5) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

ggsave(out_file_path_pdf, width=8, height=6)
