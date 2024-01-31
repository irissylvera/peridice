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

peri_nitro <- peri %>% 
  left_join(metab_groups) %>% 
  distinct(nmol_per_bulk, .keep_all = TRUE) %>% 
  filter(str_detect(emp_form, "N"))

# phosphorus stress - Zs have highest P (ZH, ZF > ZL), Rs have middle P (RL > RH), Ls have highest (LH > LL)

peri_phospho %>% 
  # filter(str_detect(metabolite,"(R)|Dimethyl|G|I|L-")) %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(treatment, levels = c("Tote", "C", "RL", "RH", "LL", "LH", "ZL", "ZF", "ZH")), 
                   y = nmol_per_bulk, color = metabolite)) + 
  facet_grid(~metabolite~date, scales = "free") + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + 
  theme(legend.position = "none") + 
  # ggtitle("Glycerophosphocholine") + 
  theme_bw()

# (R)-2,3-Dihydroxypropane-1-sulfonate, Dimethylsulfonioacetate
# Dimethylsulfoniopropionate, Glutathione disulfide, Gonyol, Isethionic acid
# L-Cystathionine, L-Cysteic acid, L-Methionine, L-Methionine S-oxide
# Taurine

peri %>% 
  filter(metabolite == "Dimethylsulfoniopropionate") %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
                   y = nmol_per_bulk, color = metabolite)) + 
  facet_wrap(~date) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

peri %>% 
  filter(metabolite == "Dimethylsulfonioacetate") %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
                   y = nmol_per_bulk, color = metabolite)) + 
  facet_wrap(~date) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

peri_phospho %>% 
  # filter(metabolite == "Gonyol") %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
                   y = nmol_per_bulk, color = metabolite)) + 
  facet_grid(~metabolite~date, scales = "free") + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

peri_sulf %>% 
  filter(str_detect(metabolite, "Dimethylsulfonioacetate|L-Methionine")) %>% 
  filter(!str_detect(metabolite, "S")) %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
                   y = nmol_per_bulk, color = metabolite)) + 
  facet_wrap(~date, scales = "free") + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

peri_nitro %>% 
  filter(metabolite == "O-Propionylcarnitine") %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
                   y = nmol_per_bulk, color = metabolite)) + 
  facet_grid(~metabolite~date) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + 
  theme(legend.position = "none") + 
  # ggtitle("Glycerophosphocholine") + 
  theme_bw()
