# generate metadataframe, env_data, samp_data, filt_data
## metadataframe = binding to msnexp object for metadata associated with xcms output
## env_data and samp_data = create filt_data used to calculate area_per_L for some quality control/checking

raw_data <- list.files('HILICpos/', full.names = TRUE)
raw_data <- as.data.frame(raw_data) %>% 
  filter(!str_detect(raw_data, "DDA"))

# extract metadata from filenames
metadataframe <- raw_data %>%
  # Grab just the unique filenames and rename the column
  distinct(filename=`raw_data`) %>%
  # Extract treatment information from the filename
  # Longest one needs to go first!")) %>%
  # Create a new column with timepoint information
  mutate(treatment = str_extract(filename, "Std|Poo|Blk|C|ZF|ZL|ZH|RL|RH|LL|LH|Tote")) %>% 
  mutate(timepoint=str_extract(filename, "27June|30June|14July|21July|27July")) %>%
  # Create a new column with sample type information
  mutate(samp_type=str_extract(filename, "Blk|Smp|Std|Poo")) %>%
  #mutate(treatment=ifelse(samp_type=="Poo", NA, treatment)) %>% 
  mutate(colid = treatment) %>% 
  mutate(filename = str_remove(filename, "HILICpos//230616_"))
  

env_data <- raw_data %>%
  # Grab just the unique filenames and rename the column
  distinct(filename=`raw_data`) %>%
  # Create a new column with sample type information
  mutate(samp_type=str_extract(filename, "Blk|Smp|Std|Poo")) %>%
  # Create a new column with timepoint information
  mutate(timepoint=str_extract(filename, "27June|30June|14July|21July|27July")) %>%
  # Extract treatment information from the filename
  # Longest one needs to go first!
  mutate(treatment=str_extract(filename, "Std|Poo|Blk|-C|ZF|ZL|ZH|RL|RH|LL|LH|Tote1|Tote2|Tote3|Tote3|Tote4")) %>%
  mutate(treatment = str_remove(treatment, "-")) %>% 
  # Replace accidental "P" treaments from "Pooled" with NAs
  # mutate(filename = str_remove(filename,"27June-|30June-|14July-|21July-|27July-")) %>% 
  mutate(replicate = str_extract(filename, "\\d.mzML")) %>% 
  mutate(replicate = str_remove(replicate, ".mzML")) %>% 
  mutate(filename = str_remove(filename, "HILICpos//230616_"))

samp_notes <- read_xlsx("PERIDICE_sample_notes.xlsx") %>% 
  mutate(treatment = str_extract(treatment, "Std|Poo|Blk|C|ZF|ZL|ZH|RL|RH|LL|LH|Tote1|Tote2|Tote3|Tote3|Tote4")) %>% 
  select(-c(init, notes, sample_collection_date)) %>%  
  dplyr::rename(timepoint = date) %>% 
  mutate(replicate = as.character(replicate)) %>% 
  group_by(timepoint, treatment, replicate)

filt_data <- env_data %>% 
  full_join(samp_notes)

write.csv(filt_data, file = "csvs/filt_data.csv", row.names = FALSE)

