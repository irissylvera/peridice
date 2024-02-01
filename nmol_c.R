# nmol C
library(tidyverse)
library(ggplot2)
library(readr)
library(readxl)

peri <- read_csv("PERIDICE_metabolite_data.csv") %>% 
  group_by(treatment, date, triplicate) %>% 
  mutate(bulk_metab = sum(nmol)) %>% 
  mutate(nmol_per_bulk = nmol/bulk_metab)

pcpn <- read_xlsx("metadata/peridice_pcpn.xlsx") %>% 
  select(tank = Tank, treatment = Treatment2, date = Date, pn = `PN (uM)`, 
         pc = `PC (uM)`, cn = Cnratio,added_N_uM = AddN) %>% 
  mutate(triplicate = str_extract(tank, "\\d")) %>% 
  mutate(tank = str_remove (tank, "1/2|1/3|2/2"))  %>% 
  mutate(treatment = str_remove(treatment, "\\d")) %>% 
  group_by(tank, date) %>% 
  mutate(pc = mean(pc)) %>% 
  distinct(pc, .keep_all=TRUE)

peri_pcpn <- peri %>% 
  ungroup() %>% 
  select(metabolite, filename, treatment, date, nmol, triplicate) %>% 
  mutate(triplicate = str_extract(filename, "(?<=Tote)\\d|\\d$")) %>% 
  left_join(pcpn, by = c("treatment", "date", "triplicate"))

peri <- peri_pcpn %>% 
  group_by(treatment, date, triplicate) %>% 
  mutate(bulk_metab = sum(nmol)) %>% 
  mutate(nmol_per_bulk = nmol/bulk_metab)


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
  group_by(metabolite, treatment, date) %>% 
  mutate(nmol_c_per_pc = mean(nmol_c_per_pc)) %>%
  distinct(nmol_c_per_pc, .keep_all = TRUE) %>% 
  filter(metabolite != "Hydroxyisoleucine") %>% 
  # arrange(desc(nmol_c_per_pc)) %>%
  # filter(metabolite %in% head(unique(metabolite), 10)) %>%
  ggplot() + 
  geom_col(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
               y = nmol_c_per_pc, fill = metabolite)) + 
  scale_y_continuous(labels = scales::percent_format()) + 
  facet_wrap(~date)

+
  theme(legend.position = "none")

# LL has higher living biomass per detritus
# Rs have most detritus build up

# Check hydroxyisoleucine, l-leucine
# sucrose - die off of cyano
