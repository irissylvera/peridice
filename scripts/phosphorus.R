# phosphorus stress
library(tidyverse)
library(ggplot2)
library(readr)

peri <- read_csv("PERIDICE_metabolite_data.csv") %>% 
  group_by(treatment, date, triplicate) %>% 
  mutate(bulk_metab = sum(nmol)) %>% 
  mutate(nmol_per_bulk = nmol/bulk_metab)

metab_groups <- read.csv(
  "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
  stringsAsFactors = FALSE, header = TRUE) %>% 
  select(metabolite = Compound_Name, emp_form = Empirical_Formula, metab_type = Compound_Type)

peri_phospho <- peri %>% 
  left_join(metab_groups) %>% 
  distinct(nmol_per_bulk, .keep_all = TRUE) %>% 
  filter(str_detect(emp_form, "P"))

peri_sulf <- peri %>% 
  left_join(metab_groups) %>% 
  distinct(nmol_per_bulk, .keep_all = TRUE) %>% 
  filter(str_detect(emp_form, "S"))

# phosphorus stress - Zs have highest P (ZH, ZF > ZL), Rs have middle P (RL > RH), Ls have highest (LH > LL)

peri_phospho %>% 
  filter(metabolite == "Glycerophosphocholine") %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
                   y = nmol_per_bulk, color = metabolite)) + 
  facet_wrap(~date) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + 
  theme(legend.position = "none") + 
  ggtitle("Glycerophosphocholine") + 
  theme_bw()

peri_sulf %>% 
  filter(metabolite == "Dimethylsulfoniopropionate") %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
                   y = nmol_per_bulk, color = metabolite)) + 
  facet_wrap(~date) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

peri %>% 
  filter(metabolite == "O-Acetylcarnitine") %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
                   y = nmol_per_bulk, color = metabolite)) + 
  facet_wrap(~date) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
