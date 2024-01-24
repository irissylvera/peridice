# nmds
library(tidyverse)
library(readr) 
library(ggplot2)
library(vegan)

# load data
peri <- read_csv("PERIDICE_metabolite_data.csv")
metadataframe <- peri %>% 
  select(filename, treatment, date) %>% 
  mutate(samp_type = str_extract(filename, "Smp|Blk|Poo")) %>% 
  distinct(filename, .keep_all = TRUE)

quant_nmds <- peri %>% 
  group_by(metabolite, treatment, date) %>%
  summarise(avg_nmol = mean(nmol_per_pc)) %>%
  ungroup() %>%
  mutate(id = paste(date, treatment)) %>% 
  filter(!str_detect(treatment, "Tote"))

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
  mutate(date = as.Date(timepoint, format = "%d-%m-%y")) %>% 
  filter(str_detect(filename, "Smp")) %>% 
  mutate(treatment = str_remove(treatment, "\\d")) %>% 
  mutate(id = paste(date, treatment))


mdsout$points %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  left_join(mds_data) %>%
  mutate(date_fct = factor(date)) %>% 
  ggplot() +
  geom_point(aes(x=MDS1, y=MDS2, color=factor(treatment, levels = level_order), shape = date_fct), size=4) +
  scale_color_manual(values = c("plum2","mediumpurple3","darkslategray2", "lightseagreen","gold", "lightcoral","orange", "brown3"),
                     labels = c("Control","0N:1P:1Fe, high dose P","0N:1P, low dose added", "0N:1P, high dose added",
                                "6N:1P, low dose added", "6N:1P, high dose added",
                                "16N:1P, low dose added", "16N:1P, high dose added"),
                     name = "Treatment") +
  scale_shape_discrete(name = "Date") + 
  theme_bw() + 
  ggtitle("Variation of Metabolite Concentration (nM) Across Treatments")