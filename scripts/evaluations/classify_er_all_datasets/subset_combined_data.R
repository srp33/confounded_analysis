# Compare all adjustment methods by their performance on classifying ER status
# Using the combined .csv file for all breast cancer datasets in /data/all_combined_data
# Will probably want to run each adjuster in parallel on the HPC using a job array

# Libraries
library(readr)
source("~/confounded_analysis/scripts/adjust/adjust.R")

# Load data and define splits
combined_data <- read_csv("/data/all_combined_data/combined.csv") # %>% column_to_rownames(var = "gene") if the first column is gene neames
glimpse(combined_data)


