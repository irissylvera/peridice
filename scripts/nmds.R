# nmds
library(tidyverse)
library(readr) 
library(ggplot2)
library(vegan)

# load data
peri <- read_csv("PERIDICE_metabolite_data.csv")%>% 
  group_by(treatment, date, triplicate) %>% 
  mutate(bulk_metab = sum(nmol)) %>% 
  mutate(nmol_per_bulk = nmol/bulk_metab)%>% 
  filter(!str_detect(filename, "Tote"))

metab_groups <- read.csv(
  "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
  stringsAsFactors = FALSE, header = TRUE) %>% 
  select(metabolite = Compound_Name, emp_form = Empirical_Formula, metab_type = Compound_Type) %>% 
  distinct()

metadataframe <- peri %>% 
  select(filename, treatment, date) %>% 
  mutate(samp_type = str_extract(filename, "Smp|Blk|Poo")) %>% 
  distinct(filename, .keep_all = TRUE)  %>% 
  filter(!str_detect(filename, "Tote"))

quant_nmds <- peri %>% 
  group_by(metabolite, treatment, date) %>%
  summarise(avg_nmol = mean(nmol_per_bulk)) %>%
  ungroup() %>%
  mutate(id = paste(date, treatment))

quant_mat <- quant_nmds %>%
  group_by(metabolite) %>%
  mutate(norm_conc = rank(avg_nmol)) %>%
  select(metabolite, norm_conc, id) %>% 
  pivot_wider(names_from = "metabolite", values_from = "norm_conc", values_fill = 0) %>%
  column_to_rownames("id") %>%
  data.matrix() 

mdsout <- quant_mat %>%
  metaMDS(k = 2, autotransform = FALSE)

mds_data <- metadataframe %>% 
  ungroup() %>% 
  #mutate(date = as.Date(timepoint, format = "%d-%m-%y")) %>% 
  filter(str_detect(filename, "Smp")) %>% 
  mutate(treatment = str_remove(treatment, "\\d")) %>% 
  mutate(id = paste(date, treatment))

mdsout$points %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  left_join(mds_data) %>%
  mutate(date_fct = factor(date)) %>% 
  ggplot() +
  geom_point(aes(x=MDS1, y=MDS2, color = factor(treatment, 
                                                      levels = c("C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), shape = date_fct), size=4) +
  scale_color_manual(name = "Treatment", values = c("plum2", "mediumpurple", "darkslategray2",
                                                    "lightseagreen", "gold", "lightcoral",
                                                    "orange", "brown3")) + 
  scale_shape_discrete(name = "Date") + 
  theme_bw() + 
  ggtitle("Variation of Metabolite Concentration (nM) Across Treatments")

# combine with culture data
aas <- metab_groups %>% 
  filter(metab_type == "Amino Acid")

sulfur <- metab_groups %>% 
  filter(metab_type == "Sulfur")

osmo <- metab_groups %>% 
  filter(metab_type == "Osmolyte")

nucleo <- metab_groups %>% 
  filter(metab_type == "Nucleoside")

quant_nmds <- peri %>% 
  group_by(metabolite, treatment, date) %>%
  filter(metabolite %in% nucleo$metabolite) %>% 
  summarise(avg_nmol = mean(nmol_per_bulk)) %>%
  ungroup() %>%
  mutate(id = paste(date, treatment))

quant_mat <- quant_nmds %>%
  group_by(metabolite) %>%
  mutate(norm_conc = rank(avg_nmol)) %>%
  select(metabolite, norm_conc, id) %>% 
  pivot_wider(names_from = "metabolite", values_from = "norm_conc", values_fill = 0) %>%
  column_to_rownames("id") %>%
  data.matrix() 

mdsout <- quant_mat %>%
  metaMDS(k = 2, autotransform = FALSE)

mds_data <- metadataframe %>% 
  ungroup() %>% 
  #mutate(date = as.Date(timepoint, format = "%d-%m-%y")) %>% 
  filter(str_detect(filename, "Smp")) %>% 
  mutate(treatment = str_remove(treatment, "\\d")) %>% 
  mutate(id = paste(date, treatment))

mdsout$points %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  left_join(mds_data) %>%
  mutate(date_fct = factor(date)) %>% 
  ggplot() +
  geom_point(aes(x=MDS1, y=MDS2, color = factor(treatment, 
                                                levels = c("C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), shape = date_fct), size=4) +
  scale_color_manual(name = "Treatment", values = c("plum2", "mediumpurple", "darkslategray2",
                                                    "lightseagreen", "gold", "lightcoral",
                                                    "orange", "brown3")) + 
  scale_shape_discrete(name = "Date") + 
  theme_bw() + 
  ggtitle("Variation of Metabolite Concentration (nM) Across Treatments")

