library(tidyverse)
library(cluster)
library(factoextra)

# load in targeted and untargeted peak data
all_peaks_stats <- read_csv("csvs/all_peaks_stats.csv")

# all_peaks_stats %>% 
#   group_by(replicate, mz, rt, treatment, date, compound) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   filter(n > 1L) 

all_peaks_stats <- all_peaks_stats %>% 
  mutate(filename = str_remove(filename, "_\\d$")) %>% 
  group_by(compound, filename) %>% 
  mutate(avg_area_per_L = mean(area_per_L)) %>% 
  distinct(avg_area_per_L) %>% 
  group_by(compound, filename) %>% 
  arrange(desc(avg_area_per_L)) %>% 
  slice(1)

stats_cluster <- all_peaks_stats %>% 
  ungroup() %>% 
  # filter(is.na(area_per_L)) %>% 
  complete(compound, filename) %>%
  mutate(avg_area_per_L = replace_na(avg_area_per_L, 0))
  # mutate(avg_area_per_L = as.character(avg_area_per_L))

all_peaks_numeric <- all_peaks_stats %>% 
  select(-c(filename)) %>% 
  pivot_wider(names_from = compound, values_from = avg_area_per_L)

all_peaks_scaled <- scale(all_peaks_numeric)
