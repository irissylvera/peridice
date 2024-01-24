# k means and foldchange

peri_km <- peri %>% 
  filter(!str_detect(metabolite, "Amino hydroxypropanesulfonate")) %>% 
  ungroup() %>% 
  distinct(nmol, .keep_all = TRUE) %>% 
  select(metabolite, nmol_per_pc, filename)  %>% 
  pivot_wider(names_from = metabolite, values_from = nmol_per_pc) %>% 
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
  left_join(peri_km, by ="metabolite")

peri_cluster %>% 
  ggplot() + 
  
  geom_boxplot(aes(x = factor(date), y = nmol, color = treatment)) + 
  facet_wrap(~cluster, scales = "free")

peri_fc <- peri %>% 
  filter(!str_detect(treatment, "Tote")) %>% 
  pivot_wider(names_from = treatment, values_from = nmol, values_fn = mean) %>% 
  pivot_longer(cols = c(ZF, ZL, ZH, LL, LH, RL, RH)) %>% 
  mutate(fold_change = value/C)

peri_cluster <- peri_fc %>% 
  left_join(peri_km, by ="metabolite") %>% 
  rename(control = C, treatment = name, nmol = value) %>% 
  mutate(cluster = factor(cluster)) %>% 
  na.omit()

peri_cluster %>% 
  group_by(cluster, treatment, date) %>% 
  mutate(nmol = mean(nmol)) %>% 
  ggplot(aes(x = date, y = nmol, 
             color = factor(treatment, levels = c("C", "ZL", "ZF", "ZH", "LL", "LH", "RL", "RH")))) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(name = "Treatment", 
                     values = c("plum2","mediumpurple3","darkslategray2", "lightseagreen","gold", "lightcoral","orange", "brown3")) + 
  facet_wrap(~cluster)

# i think i should do a k means clustering based on rank ordered foldchange?
