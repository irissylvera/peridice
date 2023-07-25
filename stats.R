library(tidyverse)
library(cluster)
library(factoextra)

all_peaks_stats <- read_csv("csvs/all_peaks_stats.csv")

all_peaks_numeric <- all_peaks_stats %>% 
  select(-c(filename)) %>% 
  pivot_wider(names_from = )

all_peaks_scaled <- scale(all_peaks_numeric)
  