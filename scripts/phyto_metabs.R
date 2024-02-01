# Goal: Trends of metabs associated with phytoplankton (Bryn's paper)
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

diatom_metabs <- data.frame(metabolite = c("Sarcosine", "Choline", "L-Cysteine",
                                           "Hydroxyproline",  "L-Ornithine","Trigonelline", "Betonicine", 
                                           "Turicine", "L-Cysteinylglycine", "Homarine",
                                           "Phosphocholine", "Butyrylcarnitine", 
                                           "Glycerophosphocholine"),
                            spec = "diatom",
                            type = "phytoplankton")

diatom_metabs <- data.frame(metabolite = c("L-Ornithine", "Homarine", "Trigonelline", ""),
                            spec = "diatom",
                            type = "phytoplankton")

dino_metabs <- data.frame(metabolite = c("4-Aminobutyric acid", "Arsenobetaine",
                                         "Gonyol"),
                          spec = "dino",
                          type = "phytoplankton")

phyto_metabs <- diatom_metabs %>% 
  rbind(dino_metabs)

peri %>% 
  filter(metabolite %in% diatom_metabs$metabolite) %>%
  # filter(metabolite == "Glycerophosphocholine") %>% 
  # filter(nmol <= 1) %>% 
  filter(!str_detect(treatment, "Tote")) %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(date), 
                   y = nmol_per_bulk, 
                   color = factor(treatment, 
                                  levels = c("C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")))) + 
  scale_color_manual(name = "Treatment", values = c("plum2", "mediumpurple", "darkslategray2",
                                                    "lightseagreen", "gold", "lightcoral",
                                                    "orange", "brown3")) + 
  facet_wrap(~metabolite, scales = "free") + 
  theme_bw()

peri %>% 
  filter(metabolite %in% dino_metabs$metabolite) %>%
  # filter(!str_detect(metabolite, "L-Ornithine|Sarcosine")) %>% 
  filter(!str_detect(treatment, "Tote")) %>% 
  distinct(nmol, .keep_all = TRUE) %>% 
  # filter(nmol <= 1) %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(treatment, levels = c("C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
                   y = nmol_per_bulk, color = metabolite)) + 
  facet_wrap(~date, scales = "free") + 
  theme_bw() + 
  ggtitle("diatom")

# Betonicine, Butyrylcarnitine, Choline, Glycerophosphocholine, 
# Homarine, L-Ornithine, Sarcosine, Trigonelline

peri %>% 
  filter(metabolite == "Arsenobetaine") %>% 
  filter(!str_detect(treatment, "Tote")) %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
                   y = nmol_per_bulk)) + 
  facet_wrap(~date) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + 
  ggtitle("Arsenobetaine") + 
  xlab("Treatment") + 
  theme_bw()

peri %>% 
  mutate(ratio = str_extract(treatment, "H$|L$|C|Tote")) %>% 
  mutate(rate = str_extract(treatment, "R|L|Z|C|Tote")) %>% 
  filter(metabolite %in% phyto_metabs$metabolite) %>%
  # filter(treatment == "RH") %>% 
  filter(metabolite != "Sarcosine") %>% 
  # filter(str_detect(filename, "Tote|30June|27July")) %>% 
  # filter(treatment != "ZF") %>% 
  ggplot() + 
  geom_col(aes(x = factor(treatment, 
                          levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
               y = nmol_per_bulk, fill = metabolite), position = "stack") +
  scale_fill_manual(values = mycolors) +
  facet_wrap(~date)
