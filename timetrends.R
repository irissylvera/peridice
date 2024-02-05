# t0 compared to day 3

peri <- read_csv("PERIDICE_metabolite_data.csv") %>% 
  group_by(treatment, date, triplicate) %>% 
  mutate(bulk_metab = sum(nmol)) %>% 
  mutate(nmol_per_bulk = nmol/bulk_metab) 

metab_groups <- read.csv(
  "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards.csv",
  stringsAsFactors = FALSE, header = TRUE) %>% 
  select(metabolite = Compound_Name, emp_form = Empirical_Formula, metab_type = Compound_Type) %>% 
  distinct()

peri_groups <- peri %>% 
  left_join(metab_groups) %>% 
  group_by(metab_type, treatment, date) %>% 
  mutate(nmol_group = mean(nmol_per_bulk)) %>% 
  distinct(nmol_group, .keep_all = TRUE) 

# %>% 
#   filter(str_detect(filename, "Tote|30June|27July"))

peri_groups %>% 
  filter(str_detect(metab_type, "Amino Acid|Osmolyte|Sugar|Sulfur|Urea")) %>% 
  ggplot() + 
  geom_col(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
               y = nmol_group, 
               fill = metab_type), 
           position = "fill") + 
  facet_wrap(~date)

aas <- metab_groups %>% 
  filter(str_detect(metab_type, "Amino Acid"))

# aa derivative is unchanged in terms of relative abundance except Ls, highest in LL
# it's hydroxyisoleucine :/
peri %>% 
  filter(metabolite %in% aas$metabolite) %>% 
  ggplot() + 
  geom_col(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
               y = nmol_per_bulk, 
               fill = metabolite), 
           position = "fill") + 
  facet_wrap(~date)

sugar <- metab_groups %>% 
  filter(str_detect(metab_type, "Sugar"))

peri %>% 
  filter(metabolite %in% sugar$metabolite) %>% 
  ggplot() + 
  geom_col(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
               y = nmol_per_bulk, 
               fill = metabolite, color = factor(triplicate)), 
           position = "stack") + 
  facet_wrap(~date)

# ectoine

# Nucleoside, Organic Acid (- nitrogenous), Osmolyte, Sugar, Sulfur, Urea cycle

nucleo <- metab_groups %>% 
  filter(str_detect(metab_type, "Nucleoside"))

peri %>% 
  filter(metabolite %in% nucleo$metabolite) %>% 
  ggplot() + 
  geom_col(aes(x = factor(treatment, levels = c("Tote", "C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), 
               y = nmol_per_bulk, 
               fill = metabolite), 
           position = "stack") + 
  facet_wrap(~date)
