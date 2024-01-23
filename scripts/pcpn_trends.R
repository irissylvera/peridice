library(tidyverse)
library(readr)
library(readxl)
library(ggplot2)
library(ggpubr)

# Goal: correspondence between metabolites and poc or pon (grab specific metabs)

# load in data
## peridice
peri <- read_csv("PERIDICE_metabolite_data.csv")
## pcpn data
pcpn <- read_xlsx("metadata/peridice_pcpn.xlsx") %>% 
  select(tank = Tank, treatment = Treatment2, date = Date, pn = `PN (uM)`, 
         pc = `PC (uM)`, cn = Cnratio,added_N_uM = AddN) %>% 
  mutate(triplicate = str_extract(tank, "\\d")) %>% 
  mutate(tank = str_remove (tank, "1/2|1/3|2/2"))  %>% 
  group_by(tank, date) %>% 
  mutate(mean_pc = mean(pc)) %>% 
  distinct(mean_pc, .keep_all = TRUE) %>% 
  mutate(replicate = paste0("230616_Smp_", ""))

# plot general metab trends
peri %>% 
  group_by(treatment, date) %>% 
  mutate(nmol = mean(nmol)) %>% 
  ggplot() + 
  geom_col(aes(x = date, y = nmol)) + 
  facet_wrap(~treatment)
#### Observations: RH LH have highest week 2 and 3 values, 
#### 4 week mark RH and LL have highest values
#### two different reasons for that, RH has enough nutrients but hasn't come to equilibrium? 
#### potentially about to die off like RL did and LH did, whereas LL didn't have that sharp 
#### growth like a bloom, just grew slowly and adjusted.

# plot general pc/pn trends
## pc
pcpn %>%  
  filter(!str_detect(treatment, "Tote")) %>% 
  group_by(triplicate, treatment, date) %>%  
  ggplot(aes(x = date, y = pc, group = interaction(treatment, triplicate), color = treatment)) +
  geom_line() + 
  geom_point() + 
  facet_wrap(~treatment) + 
  theme_bw() +  
  ylab("Particulate Carbon") + 
  xlab("Date")

#### Observations: RH RL Lh all grow consistently, LL doesn't have significant change

## pn 
pcpn %>%  
  filter(!str_detect(treatment, "Tote")) %>% 
  group_by(triplicate, treatment, date) %>%  
  ggplot(aes(x = date, y = pn, group = interaction(treatment, triplicate), color = treatment)) +
  geom_line() + 
  geom_point() + 
  facet_wrap(~treatment) + 
  theme_bw() +  
  ylab("Particulate Nitrogen") + 
  xlab("Date")

pcpn %>%  
  filter(!str_detect(treatment, "Tote")) %>% 
  group_by(triplicate, treatment, date) %>%  
  ggplot(aes(x = date, y = cn, group = interaction(treatment, triplicate), color = treatment)) +
  geom_line() + 
  geom_point() + 
  facet_wrap(~treatment) + 
  theme_bw() +  
  ylab("C:N Ratio") + 
  xlab("Date")

#### Observations: RH and LH consistently grow with particulate N
#### this goes into why do all of the Zs have no PN accumulation? What's happening in those systems?

# Correlation
## use k means? cluster that follows pc trend or pn trend?

## k means
peri_km <- peri %>% 
  filter(!str_detect(metabolite, "Amino hydroxypropanesulfonate")) %>% 
  ungroup() %>% 
  distinct(nmol, .keep_all = TRUE) %>% 
  select(metabolite, nmol, filename)  %>% 
  pivot_wider(names_from = metabolite, values_from = nmol) %>% 
  column_to_rownames("filename") %>% 
  data.matrix() %>%
  `[<-`(is.na(.), 0) %>%
  scale() %>%
  t() %>%
  kmeans(centers = 10) %>% 
  pluck("cluster") %>%
  as.data.frame() %>%
  set_names("cluster") %>%
  rownames_to_column("metabolite")

peri_cluster <- peri %>% 
  left_join(peri_km, by = c("metabolite")) %>% 
  mutate(cluster = factor(cluster)) %>% 
  na.omit()

## compare to PC
peri_cluster %>% 
  group_by(cluster, treatment, date) %>% 
  mutate(nmol = mean(nmol)) %>% 
  ggplot() + 
  geom_line(aes(x = date, y = nmol, color = cluster)) + 
  facet_wrap(~treatment, scales = "free")+
  theme_bw()

peri_cluster %>% 
  filter(cluster == 2) %>% 
  group_by(metabolite, treatment, date) %>% 
  mutate(nmol = mean(nmol)) %>% 
  ggplot() + 
  geom_line(aes(x = date, y = nmol, color = metabolite)) + 
  facet_wrap(~treatment)+
  theme_bw()

peri_cluster %>% 
  group_by(metabolite, treatment, date) %>% 
  mutate(nmol = mean(nmol)) %>% 
  distinct(nmol, .keep_all = TRUE) %>% 
  ggplot() + 
  geom_line(aes(x = factor(date), y = nmol, group = metabolite, color = metabolite)) + 
  facet_wrap(~treatment)+
  theme_bw() 

peri %>% 
  ggplot() + 
  geom_point(aes(x = nmol, y = pn)) + 
  facet_wrap(~metabolite)

pc_metabs <- data.frame(metabolite = c("2-O-alpha-D-Glucosylglycerol",
                                       "Dimethylsulfoniopropionate", 
                                       "L-Alanine", "L-Glutamine"))
peri %>% 
  filter(metabolite %in% pc_metabs$metabolite) %>% 
  ggplot() + 
  geom_point(aes(x = nmol, y = pc, color = treatment)) + 
  facet_wrap(~metabolite, scales = "free") + 
  theme_bw()

peri %>% 
  filter(metabolite %in% pc_metabs$metabolite) %>% 
  ggscatter(x = "nmol", y = "pc", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = "nmol", ylab = "pc", color = "treatment",
            palette = c(C = "plum2", ZF = "mediumpurple3", ZL = "darkslategray2", 
                        ZH = "lightseagreen",
                       LL = "gold", LH = "lightcoral",
                       RL = "orange", RH = "brown3")) + 
  facet_wrap(~metabolite, scales = "free")

pn_metabs <- data.frame(metabolite = c("2-O-alpha-D-Glucosylglycerol",
                                       " Dimethylsulfoniopropionate",
                                       "L-Alanine", "L-Aspartic acid",
                                       "L-Glutamine", "L-Leucine"))

peri %>% 
  filter(metabolite %in% pn_metabs$metabolite) %>% 
  ggplot() + 
  geom_point(aes(x = nmol, y = pn, color = treatment)) + 
  facet_wrap(~metabolite, scales = "free") + 
  theme_bw()

peri %>% 
  group_by(metabolite) %>% 
  distinct(nmol, .keep_all = TRUE) %>% 
  distinct(pc, .keep_all = TRUE) %>% 
  mutate(rank_nmol = rank(nmol)) %>% 
  mutate(rank_pc = rank(pc)) %>% 
  filter(metabolite %in% pc_metabs$metabolite) %>% 
  ggscatter(x = "nmol", y = "pc", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = "nmol", ylab = "pc", color = "treatment",
            palette = c(C = "plum2", ZF = "mediumpurple3", ZL = "darkslategray2", 
                        ZH = "lightseagreen",
                        LL = "gold", LH = "lightcoral",
                        RL = "orange", RH = "brown3")) + 
  facet_wrap(~metabolite, scales = "free")

peri %>%
  group_by(metabolite) %>%
  filter(metabolite %in% pc_metabs$metabolite) %>%
  distinct(nmol, .keep_all = TRUE) %>% 
  distinct(pc, .keep_all = TRUE) %>% 
  summarize(spear_val=cor(pc, nmol, method="spearman", use = "pairwise"),
            pear_val=cor(pc, nmol, method="pearson", use = "pairwise"),
            rank_pear=cor(rank(pc), rank(nmol), method="pearson"),
            rank_spear=cor(rank(pc), rank(nmol), method="spearman"))

sulfur <- metab_groups %>% 
  filter(metab_type == "Sulfur")

peri %>% 
  filter(metabolite %in% sulfur$metabolite) %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(date), y = nmol_per_pc, color = treatment)) + 
  facet_wrap(~metabolite, scales = "free") + 
  theme_bw()


sugar <- metab_groups %>% 
  filter(metab_type == "Sugar")

peri %>% 
  filter(metabolite %in% sugar$metabolite) %>% 
  ggplot() + 
  geom_boxplot(aes(x = factor(date), y = nmol_per_pc, color = treatment)) + 
  facet_wrap(~metabolite, scales = "free") + 
  theme_bw()
