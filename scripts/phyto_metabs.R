# Goal: Trends of metabs associated with phytoplankton (Bryn's paper)
library(tidyverse)
library(readr)
library(readxl)

# load in data
## peridice
peri <- read_csv("PERIDICE_metabolite_data.csv")
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
  filter(nmol <= 10) %>% 
  ggplot() +
  geom_col(aes(x = factor(date), y = nmol, fill = metabolite), position = "stack") + 
  facet_wrap(~treatment) + 
  theme_bw() + 
  ggtitle("diatoms")

phospho <- data.frame(metabolite = c("Glycerophosphocholine", "Choline", "Phosphocholine"))

peri %>% 
  filter(str_detect(metabolite, "Glycerophosphocholine")) %>%
  # mutate(trend = ifelse(treatment %in% c("LH|LL|ZL" > "RH|RL|ZH"),"red","blue")) %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(treatment, levels = c("RH", "LH", "RL", "LL", "ZH", "ZF", "ZL", "C")), 
                   y = nmol, fill = metabolite)) + 
  facet_wrap(~date, scales = "free") + 
  theme_bw() + 
  ggtitle("phytos")

sulfur <- metab_groups %>% 
  filter(metab_type == "Sulfur")

non_phospho <- metab_groups %>% 
  filter(!str_detect(emp_form, "P"))

peri %>% 
  filter(metabolite %in% non_phospho$metabolite) %>%
  group_by(metabolite, treatment, date, triplicate) %>% 
  mutate(nmol = sum(nmol)) %>% 
  filter(!str_detect(treatment, "Tote")) %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(treatment, levels = c("RH", "LH", "RL", "LL", "ZH", "ZF", "ZL", "C")), 
                   y = nmol_per_pc)) + 
  facet_grid(~metabolite~date, scales = "free") + 
  theme_bw() + 
  ggtitle("phytos")
