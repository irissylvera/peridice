# trends with function

library(tidyverse)
library(ggplot2)
library(readr)
library(RColorBrewer)

peri <- read_csv("PERIDICE_metabolite_data.csv") %>% 
  group_by(treatment, date, triplicate) %>% 
  mutate(bulk_metab = sum(nmol)) %>% 
  mutate(nmol_per_bulk = nmol/bulk_metab)

peri <- read_csv("csvs/quant_data.csv") %>% 
  select(metabolite, filename = replicate, nmol = nM) %>% 
  filter(str_detect(filename, "Smp")) %>% 
  mutate(treatment = str_extract(filename, "C|ZL|ZF|ZH|LL|LH|RL|RH|Tote")) %>% 
  mutate(timepoint = str_extract(filename, "27June|30June|14July|21July|27July")) %>% 
  mutate(timepoint = str_replace_all(timepoint, "14July", "7-14-22")) %>% 
  mutate(timepoint = str_replace_all(timepoint, "21July", "7-21-22")) %>% 
  mutate(timepoint = str_replace_all(timepoint, "30June", "6-30-22")) %>% 
  mutate(timepoint = str_replace_all(timepoint, "27July", "7-27-22")) %>% 
  mutate(timepoint = str_replace_all(timepoint, "27June", "6-27-22")) %>% 
  mutate(date = as.Date(timepoint, format = "%m-%d-%y")) %>% 
  select(-c("timepoint")) %>% 
  group_by(filename) %>% 
  mutate(bulk_metab = sum(nmol)) %>% 
  ungroup() %>% 
  mutate(nmol_per_bulk = nmol/bulk_metab)

metab_groups <- read.csv(
  "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
  stringsAsFactors = FALSE, header = TRUE) %>% 
  select(metabolite = Compound_Name, emp_form = Empirical_Formula, metab_type = Compound_Type) %>% 
  distinct()

peri_groups <- peri %>% 
  left_join(metab_groups) %>% 
  group_by(metab_type, treatment, date) %>% 
  mutate(nmol_group = mean(nmol_per_bulk)) 

# %>% 
#   distinct(nmol_group, .keep_all = TRUE)

groups_of_int <- data.frame(metab_type = c("Amino Acid", "Amino Acid - degraded", "Amino Acid derivative", 
                                           "Amino Acid synthesis",
                                           "Nucleic Acid", "Nucleoside", "Nucleoside derivative", "Osmolyte", 
                                           "Sulfur", "Urea cycle"))

nb.cols <- 21
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)

sulf <- metab_groups %>% 
  filter(metab_type=="Sulfur")

peri_groups %>% 
  filter(str_detect(metab_type, "Amino Acid")) %>% 
  filter(metab_type != "Amino Acid") %>% 
  filter(!str_detect(filename, "Tote")) %>% 
  # filter(metab_type %in% groups_of_int$metab_type) %>% 
  ggplot() + 
  geom_col(aes(x = factor(treatment, levels = c("C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
               y = nmol_per_bulk, 
               fill = metab_type), 
           position = "fill") + 
  scale_fill_manual(values = mycolors) +
  facet_wrap(~date)

peri_groups %>% 
  mutate(ratio = str_extract(treatment, "H$|L$|C|Tote")) %>% 
  mutate(rate = str_extract(treatment, "R|L|Z|C|Tote")) %>% 
  filter(metab_type %in% groups_of_int$metab_type) %>% 
  # filter(treatment == "RH") %>% 
  filter(str_detect(filename, "Tote|30June|27July")) %>% 
  filter(treatment != "ZF") %>% 
  ggplot() + 
  geom_col(aes(x = factor(treatment, 
                          levels = c("Tote", "C", "ZL", "ZH", "LL", "LH", "RL", "RH")), 
               y = nmol_per_bulk, fill = metab_type), position = "fill") +
  scale_fill_manual(values = mycolors) +
  facet_wrap(~date)

peri_groups %>% 
  # filter(metab_type == "Sulfur") %>% 
  mutate(ratio = str_extract(treatment, "H$|L$|C|Tote")) %>% 
  mutate(rate = str_extract(treatment, "R|L|Z|C|Tote")) %>% 
  # filter(metab_type %in% groups_of_int$metab_type) %>% 
  # filter(treatment == "RH") %>% 
  # filter(str_detect(filename, "Tote|30June|27July")) %>% 
  filter(treatment != "ZF") %>% 
  ggplot() + 
  geom_col(aes(x = factor(date), 
               y = nmol_per_bulk, fill = metab_type), position = "stack") +
  scale_fill_manual(values = mycolors) +
  facet_wrap(~factor(treatment, 
                     levels = c("Tote", "C", "ZL", "ZH", "LL", "LH", "RL", "RH")))

# sulfur metabs accumulating in RL because euks are blooming - make sulfur metabs


peri_groups %>% 
  filter(!str_detect(filename, "Tote")) %>% 
  mutate(rate = str_extract(treatment, "H$|L$|C|Tote")) %>% 
  mutate(ratio = str_extract(treatment, "R|L|Z|C|Tote")) %>% 
  group_by(treatment, date) %>% 
  mutate(nmol = sum(bulk_metab)) %>% 
  # filter(treatment == "RH") %>% 
  # filter(str_detect(filename, "Tote|30June|27July")) %>% 
  filter(treatment != "ZF") %>% 
  ggplot() + 
  geom_col(aes(x = factor(date), 
               y = nmol, fill = treatment)) +
  # scale_fill_manual(values = mycolors) +
  facet_grid(~ratio~rate)
  

  
  
  
