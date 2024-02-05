# k means and foldchange

set.seed(20)

peri_km <- peri %>% 
  # filter(!str_detect(metabolite, "Amino hydroxypropanesulfonate")) %>% 
  ungroup() %>% 
  distinct(nmol, .keep_all = TRUE) %>% 
  select(metabolite, nmol_per_bulk, filename)  %>% 
  pivot_wider(names_from = metabolite, values_from = nmol_per_bulk) %>% 
  column_to_rownames("filename") %>% 
  data.matrix() %>%
  `[<-`(is.na(.), 0) %>%
  scale() %>%
  t() %>%
  kmeans(centers = 6) %>% 
  pluck("cluster") %>%
  as.data.frame() %>%
  set_names("cluster") %>%
  rownames_to_column("metabolite")

peri_cluster <- peri %>% 
  left_join(peri_km, by ="metabolite")

peri_fc <- peri %>% 
  filter(!str_detect(treatment, "Tote")) %>% 
  select(metabolite, treatment, date, nmol_per_bulk) %>% 
  pivot_wider(names_from = treatment, values_from = nmol_per_bulk, values_fn = mean) %>% 
  pivot_longer(cols = c(ZF, ZL, ZH, LL, LH, RL, RH)) %>% 
  mutate(fold_change = value/C)

peri_cluster <- peri_fc %>% 
  left_join(peri_km, by ="metabolite") %>% 
  rename(control = C, treatment = name, nmol_per_bulk = value) %>% 
  mutate(cluster = factor(cluster)) %>% 
  na.omit()

peri_cluster %>% 
  group_by(cluster, treatment, date) %>% 
  mutate(nmol = mean(nmol_per_bulk)) %>% 
  ggplot(aes(x = date, y = nmol, 
             color = factor(treatment, levels = c("C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")))) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(name = "Treatment", 
                     values = c("plum2","mediumpurple3","darkslategray2", "lightseagreen","gold", "lightcoral","orange", "brown3")) + 
  facet_wrap(~cluster, scales = "free")
