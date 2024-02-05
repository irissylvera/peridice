# correlation with rate

# load in data
## peridice
peri <- read_csv("PERIDICE_metabolite_data.csv") %>% 
  ungroup() %>% 
  group_by(treatment, date) %>% 
  mutate(bulk_metab = sum(nmol)) %>% 
  mutate(nmol_per_bulk = nmol/bulk_metab)
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

peri %>% 
  ggplot() + 
  geom_point(aes(x = nmol_per_bulk, y = pc)) + 
  facet_wrap(~metabolite)

peri %>% 
  #filter(metabolite %in% pc_metabs$metabolite) %>% 
  ggscatter(x = "nmol_per_bulk", y = "pc", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "nmol", ylab = "pc", color = "treatment",
            palette = c(C = "plum2", ZF = "mediumpurple3", ZL = "darkslategray2", 
                        ZH = "lightseagreen",
                        LL = "gold", LH = "lightcoral",
                        RL = "orange", RH = "brown3")) + 
  facet_wrap(~metabolite, scales = "free")

peri_corr <- peri %>%
  group_by(metabolite) %>%
  distinct(bulk_metab, .keep_all = TRUE) %>% 
  distinct(pc, .keep_all = TRUE) %>% 
  summarize(spear_val=cor(pc, bulk_metab, method="spearman", use = "pairwise"),
            pear_val=cor(pc, bulk_metab, method="pearson", use = "pairwise"),
            rank_pear=cor(rank(pc), rank(bulk_metab), method="pearson"),
            rank_spear=cor(rank(pc), rank(bulk_metab), method="spearman")) %>% 
  left_join(peri)

peri_corr %>% 
  filter(!str_detect(treatment, "Tote")) %>% 
  # filter(spear_val > 0.5) %>% 
  ggplot(aes(x = bulk_metab, 
             y = pc)) + 
  geom_point(aes(color = factor(treatment, 
                                levels = c("C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")), size = spear_val)) + 
  scale_color_manual(name = "Treatment", values = c("plum2", "mediumpurple", "darkslategray2",
                                                    "lightseagreen", "gold", "lightcoral",
                                                    "orange", "brown3")) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") + 
  xlab("Bulk Metabolite Pool (nmol)") + 
  ylab("Particulate Organic Carbon (umol)") + 
  theme_bw()

