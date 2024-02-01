# amino acids
library(tidyverse)
library(readr)
library(readxl)

# load in data
## peridice
peri <- read_csv("PERIDICE_metabolite_data.csv") %>% 
  group_by(treatment, date, triplicate) %>% 
  mutate(bulk_metab = sum(nmol)) %>% 
  mutate(nmol_per_bulk = nmol/bulk_metab)

## pcpn data
pcpn <- read_xlsx("metadata/peridice_pcpn.xlsx") %>% 
  select(tank = Tank, treatment = Treatment2, date = Date, pn = `PN (uM)`, 
         pc = `PC (uM)`, cn = Cnratio,added_N_uM = AddN) %>% 
  mutate(triplicate = str_extract(tank, "\\d")) %>% 
  mutate(tank =  str_remove (tank, "1/2|1/3|2/2"))  %>% 
  group_by(tank, date) %>% 
  mutate(mean_pc = mean(pc)) %>% 
  distinct(mean_pc, .keep_all = TRUE) %>% 
  mutate(replicate = paste0("230616_Smp_", ""))

metab_groups <- read.csv(
  "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
  stringsAsFactors = FALSE, header = TRUE) %>% 
  select(metabolite = Compound_Name, emp_form = Empirical_Formula, metab_type = Compound_Type)


aas <- metab_groups %>% 
  # filter(metab_type == "Amino Acid") %>% 
  filter(metabolite == "L-Alanine" | metabolite == "L-Proline" | metabolite == "L-Threonine" | 
           metabolite == "L-Leucine" | metabolite == "L-Asparagine" | metabolite == "L-Glutamine" |
           metabolite=="L-Methionine"|metabolite=="L-Histidine"|metabolite=="L-Arginine"|metabolite=="L-Tyrosine")
# different trends: "L-Aspartic acid"

peri %>% 
  # filter(metabolite %in% aas$metabolite) %>% 
  filter(!str_detect(treatment, "Tote")) %>% 
  distinct(nmol, .keep_all = TRUE) %>% 
  filter(metabolite == "L-Glutamic acid") %>% 
  # filter(nmol <= 1) %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(treatment, levels = c("C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
                   y = nmol_per_bulk, color = metabolite)) + 
  facet_wrap(~date, scales = "free") + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

# L-Alanine, L-Threonine
# Rate sensitive: L-Proline. L-Methionine
# L-Leucine, L-Glutamine, L-Histidine
# Ratio sensitive: L-Asparagine. L-Arginine(?), 
# L-Tyrosine is ratio at 2 weeks, rate at 3, ratio at 4

peri %>% 
  filter(metabolite %in% aas$metabolite) %>% 
  filter(!str_detect(treatment, "Tote")) %>% 
  distinct(nmol, .keep_all = TRUE) %>% 
  group_by(treatment, date) %>% 
  mutate(nmol = sum(nmol)) %>% 
  ggplot() + 
  geom_col(aes(x = factor(treatment, levels = c("C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), y = nmol)) + 
  facet_wrap(~date) + 
  ggtitle("Bulk Amino Acids") + 
  xlab("Treatment")
