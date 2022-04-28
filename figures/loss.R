# Load packages -----
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

# Load data ---------
df <- read_csv("/output/metrics/mnist_confounded_log.csv") %>% 
  select(ae_loss, disc_loss, iteration) %>% 
  gather("Function", "Loss Function", -iteration) %>% 
  mutate(Function = str_replace(Function, "_loss", "")) %>%
  mutate(Function = str_replace(Function, "ae", "Autoencoder")) %>%
  mutate(Function = str_replace(Function, "disc", "Discriminator")) %>%
  rename(Iteration = iteration)

# This is a workaround for a ggplot2 bug.
pdf(NULL)

#filter(df, Iteration > 1000) %>% 
  ggplot(df, aes(x = Iteration, y = `Loss Function`)) + 
  geom_line() + 
  geom_smooth() + 
#  geom_vline(xintercept = 9600, color = "Red", linetype = "longdash") +
  facet_wrap(vars(Function), nrow=2, scales = "free_y") +
  theme_bw(base_size = 18) +
  ylab("Loss value")

ggsave("/output/figures/loss_chart.pdf")
