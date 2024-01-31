# nmol C
library(tidyverse)
library(ggplot2)
library(readr)

peri <- read_csv("PERIDICE_metabolite_data.csv") %>% 
  group_by(treatment, date, triplicate) %>% 
  mutate(bulk_metab = sum(nmol)) %>% 
  mutate(nmol_per_bulk = nmol/bulk_metab)

peri <- peri_pcpn %>% 
  group_by(treatment, date, triplicate) %>% 
  mutate(bulk_metab = sum(nmol)) %>% 
  mutate(nmol_per_bulk = nmol/bulk_metab)

peri %>% 
  filter(metabolite == "Sucrose") %>% 
  filter(str_detect(filename, "27July-RH_1"))

metab_groups <- read.csv(
  "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
  stringsAsFactors = FALSE, header = TRUE) %>% 
  select(metabolite = Compound_Name, emp_form = Empirical_Formula, metab_type = Compound_Type) %>% 
  distinct()

peri_c <- peri %>% 
  ungroup() %>% 
  left_join(metab_groups) %>% 
  mutate(carbon = str_remove(emp_form, "C\\d{1,3}")) %>% 
  mutate(carbon = str_extract(carbon, "\\d{1,3}")) %>% 
  mutate(carbon = as.numeric(carbon)) %>% 
  mutate(nmol_c_per_pc = nmol*carbon / (pc*1000))

peri_c %>% 
  arrange(desc(nmol_c_per_pc)) %>% 
  filter(metabolite %in% head(unique(metabolite), 10)) %>% 
  ggplot() + 
  geom_col(aes(x = treatment, 
               y = nmol_c_per_pc, fill = metabolite)) + 
  scale_y_continuous(labels = scales::percent_format()) + 
  facet_wrap(~date)
