library(tidyverse)

# load in untargeted peak list of chosen feats
untargeted_peak_list <- read_csv("csvs/untargeted_peak_list.csv")

# load in targeted peak list
targeted_peak_list <- read_csv("csvs/targeted_peak_list.csv")

# load in filtering metadata generated in dataframe setup
filt_data <- read_csv("csvs/filt_data.csv")

untarg_peak_tidy <- untargeted_peak_list %>% 
  select(-c(timepoint)) %>% 
  mutate(filename = str_remove(filename, ".mzML")) %>% 
  mutate(replicate = as.character(replicate))

targ_peak_list <- targeted_peak_list %>% 
  dplyr::rename(compound = Peptide) %>% 
  mutate(timepoint = str_extract(Replicate, "27June|30June|14July|21July|27July")) %>% 
  mutate(Replicate = str_remove(Replicate, "230616_")) %>% 
  mutate(timepoint = str_replace(timepoint, "27June", "6-27")) %>% 
  mutate(timepoint = str_replace(timepoint, "30June", "6-30")) %>% 
  mutate(timepoint = str_replace(timepoint, "14July", "7-14")) %>% 
  mutate(timepoint = str_replace(timepoint, "21July", "7-21")) %>% 
  mutate(timepoint = str_replace(timepoint, "27July", "7-27")) %>% 
  mutate(date = as.Date(timepoint, format="%m-%d")) %>% 
  select(-c(Protein, `Precursor Charge`, `Product Charge`, `Product Mz`, `Fragment Ion`, `Peak Rank`)) %>% 
  dplyr::rename(mz = `Precursor Mz`) %>% 
  dplyr::rename(rt = `Retention Time`) %>% 
  mutate(rt = as.numeric(rt)) %>% 
  mutate(Area = as.numeric(Area)) %>% 
  mutate(Background = as.numeric(Background)) %>% 
  mutate(rt = rt * 60) %>% 
  # na.omit() %>% 
  rename(filename = Replicate) %>% 
  rename(area = Area) %>% 
  mutate(treatment = str_extract(filename, "Std|Poo|Blk|C|ZF|ZL|ZH|RL|RH|LL|LH|Tote")) %>% 
  group_by(compound, filename, treatment) %>% 
  mutate(thresh = area/Background) %>% 
  filter(thresh > 1) %>% 
  select(-c(thresh)) %>%
  mutate(replicate = str_extract(filename, "\\d$"))
  # mutate(replicate = as.character(replicate))

filt_data_tidy <- filt_data %>% 
  mutate(timepoint = str_replace(timepoint, "27June", "6-27")) %>% 
  mutate(timepoint = str_replace(timepoint, "30June", "6-30")) %>% 
  mutate(timepoint = str_replace(timepoint, "14July", "7-14")) %>% 
  mutate(timepoint = str_replace(timepoint, "21July", "7-21")) %>% 
  mutate(timepoint = str_replace(timepoint, "27July", "7-27")) %>% 
  mutate(date = as.Date(timepoint, format="%m-%d")) %>% 
  na.omit() %>% 
  mutate(filename = str_remove(filename, ".mzML")) %>% 
  mutate(replicate = as.character(replicate))
  
  
targ_data_qc <- targ_peak_list %>% 
  left_join(filt_data_tidy, by = join_by(treatment, replicate, timepoint)) %>% 
  mutate(area_per_L = area/vol_L_filt) %>% 
  na.omit() %>% 
  select(c(compound, filename.x, replicate, mz, rt, treatment, date.x, area_per_L)) %>% 
  rename(filename = filename.x) %>% 
  rename(date = date.x)

write.csv(targ_data_qc, file = "csvs/targ_data_qc.csv", row.names = FALSE)

targ_data_qc <- read_csv("csvs/targ_data_qc.csv")

all_peaks_stats <- targ_data_qc %>% 
  rbind(untarg_peak_tidy)

write.csv(all_peaks_stats, file = "csvs/all_peaks_stats.csv", row.names = FALSE)
