# XCMS output and extract peak_data
msnexp_filled <- readRDS("./xcms_output/msnexp_filled2.rds")

peak_data_long <- msnexp_filled %>%
  chromPeaks() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(peakidx=row_number())

peak_data <- msnexp_filled %>%
  featureDefinitions() %>%
  as.data.frame() %>%
  select(mzmed, rtmed, npeaks, peakidx) %>%
  rownames_to_column("id") %>%
  unnest_longer(peakidx) %>%
  rename_with(~paste0("feat_", .x), .cols = -peakidx) %>%
  left_join(peak_data_long) %>%
  mutate(filename=basename(fileNames(msnexp_filled))[sample]) %>% 
  mutate(filename, str_extract(filename, "230616_.*")) %>% 
  mutate(`str_extract(filename, "230616_.*")`, str_remove(`str_extract(filename, "230616_.*")`, "230616_")) %>% 
  select(-c("filename", 'str_extract(filename, "230616_.*")')) %>% 
  rename(filename = 'str_remove(`str_extract(filename, "230616_.*")`, "230616_")') %>% 
  left_join(metadataframe %>% mutate(filename=basename(filename))) %>%
  filter(samp_type=="Smp") %>% 
  rename(compound = feat_id) %>% 
  rename(area = into)

peak_data_clean <- peak_data %>% 
  filter(samp_type=="Smp") %>% 
  select(compound, filename, area, mz, rt, treatment) %>% 
  group_by(compound, filename) %>%
  # mutate(n=n()) %>%
  group_by(compound) %>%
  # filter(all(n==1)) %>%
  ungroup() %>%
  # complete(compound, filename, fill=list(area=0)) %>%
  # select(-n) %>%
  # Extract treatment information from the filename
  # Longest one needs to go first!
  mutate(treatment=str_extract(filename, "Std|Poo|Blk|C|ZF|ZL|ZH|RL|RH|LL|LH|Tote1|Tote2|Tote3|Tote3|Tote4")) %>% 
  mutate(timepoint = str_extract(filename, "27June|30June|14July|21July|27July")) %>% 
  mutate(replicate = str_extract(filename, "\\d.mzML")) %>% 
  mutate(replicate = str_remove(replicate, ".mzML")) %>% 
  group_by(timepoint, treatment, replicate)

