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
  mutate(filename = str_remove(filename, "HILICpos//230616_")) %>% 
  # Extract treatment information from the filename
  # Longest one needs to go first!")) %>%
  # Create a new column with timepoint information
  mutate(treatment = str_extract(filename, "Std|Poo|Blk|C|ZF|ZL|ZH|RL|RH|LL|LH|Tote1|Tote2|Tote3|Tote4")) %>% 
  mutate(timepoint=str_extract(filename, "27June|30June|14July|21July|27July")) %>%
  mutate(replicate = str_extract(filename, "\\d.mzML")) %>% 
  mutate(replicate = str_remove(replicate, ".mzML")) %>% 
  # Create a new column with sample type information
  mutate(samp_type=str_extract(filename, "Blk|Smp|Std|Poo")) %>%
  #mutate(treatment=ifelse(samp_type=="Poo", NA, treatment)) %>% 
  mutate(colid = levels(treatment)) %>% 
  group_by(timepoint, treatment, replicate)
  

# env_data <- metadataframe
# 
# samp_notes <- read_xlsx("xls/PERIDICE_sample_notes.xlsx") %>% 
#   mutate(treatment = str_extract(treatment, "Std|Poo|Blk|C|ZF|ZL|ZH|RL|RH|LL|LH|Tote1|Tote2|Tote3|Tote4")) %>% 
#   select(-c(init, notes, sample_collection_date)) %>%  
#   dplyr::rename(timepoint = date) %>% 
#   mutate(replicate = as.character(replicate)) 
# 
# filt_data <- env_data %>% 
#   left_join(samp_notes) 
# 
# write.csv(filt_data, file = "csvs/filt_data.csv", row.names = FALSE)

#filt data 
filt_data <- read_csv("csvs/filt_data.csv") %>% 
  mutate(timepoint = str_replace_all(timepoint, "14July", "7-14-22")) %>% 
  mutate(timepoint = str_replace_all(timepoint, "21July", "7-21-22")) %>% 
  mutate(timepoint = str_replace_all(timepoint, "30June", "6-30-22")) %>% 
  mutate(timepoint = str_replace_all(timepoint, "27July", "7-27-22")) %>% 
  mutate(timepoint = str_replace_all(timepoint, "27June", "6-27-22")) %>% 
  mutate(date = as.Date(timepoint, format = "%m-%d-%y")) %>% 
  select(-c("timepoint", "samp_type"))

# pcpn
pcpn <- read_xlsx("metadata/peridice_pcpn.xlsx") %>% 
  select(Tank, Treatment2, Date, `PN (uM)`, `PC (uM)`, Cnratio, AccN, AccC, AddNn, AddN) %>% 
  rename(tank = Tank, treatment = Treatment2, date = Date, pn_um = `PN (uM)`, pc_um = `PC (uM)`, cn_ratio = Cnratio, 
         acc_n = AccN, acc_c = AccC, add_nn = AddNn, add_n = AddN) %>% 
  mutate(tank = str_remove(tank, "/(2|3)")) %>% 
  mutate(replicate = str_extract(tank, "\\d$")) %>% 
  mutate(replicate = as.numeric(replicate)) %>% 
  select(-c("tank"))

# add_n is nitrate add_nn is nitrate_nitrite

# setup metadata
metadata_targ <- filt_data %>% 
  left_join(pcpn, by = c("date", "replicate", "treatment")) %>% 
  na.omit() %>% 
  select(-c("replicate")) %>% 
  rename(replicate = filename) %>% 
  mutate(replicate = str_remove(replicate, ".mzML"))

