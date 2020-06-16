
ccnd1kd_norm <- ccnd1kd %>% 
  #mutate(KELLY = log2(KELLY) ) %>% 
  #mutate(SKNAS = log2(SKNAS) ) %>% 
  mutate(tf = substr(siRNA, 1, nchar(siRNA)-2)) %>% 
  mutate(Selected = if_else(tf == "EBF1", 1L , Selected)) %>% 
  filter(Selected == 1) %>% 
  #filter(siRNA %in% selected_tf) %>% 
  
  mutate(tf = if_else(tf == "MEIS", "MEIS2", tf))

# 
# table(x$tf %in% df$symb)
# x$tf[!x$tf %in% df$symb]



# y <- data.frame(id = names(TTT),
#                 sym = unlist(sapply(TTT, "[", 3)))

# 
# head(footprint_list)
# 
# 
# footprint_list[[1]]
# 
# lapply(head(footprint_list), function(t)unique(t[[1]]$ID))
# 
# 
# y <- data.frame(id = names(footprint_list),
#                 sym = unlist(sapply(footprint_list, function(t) unique(t[[1]]$ID))),
#                 sym2 = unlist(sapply(footprint_list, function(t) unique(t[[2]]$ID))))
# 


#
#df <- readRDS("data/cells/atacseq/footprint/footprints_SKNASvsKELLY_ttest.RDS")
MES_activity <- readRDS("analysis/tumor/VIPER/MES_TFactivity.RDS")
names(MES_activity$es$p.value)[MES_activity$es$p.value < 0.05]
MES_activity <- MES_activity$es$nes



footprint_ttest <- readRDS("data/cells/atacseq/footprint/footprints_SKNASvsKELLY_ttest_flipped.RDS")
head(footprint_ttest)
footprint_ttest$name[!ccnd1kd_norm$tf %in% footprint_ttest$name]

ccnd1kd_norm %>% 
  mutate(oppening = -1*footprint_ttest$val[match(tf, footprint_ttest$name)]) %>% 
  mutate(MES_activity =  MES_activity[match( tf, names(MES_activity))]) %>% 
  #mutate(mes =  1) %>% 
  arrange(oppening) %>% 
  mutate(pos = 1:n()) %>% 
  #mutate(ccnd1rel = KELLY/SKNAS) %>% 
  mutate(ccnd1rel = SKNAS - KELLY) %>% 
  ggplot(aes(x = oppening, y = ccnd1rel, color = MES_activity)) +
  #ggplot(aes(x = pos, y = ccnd1rel)) +
  geom_point() +
  geom_vline(xintercept = 0) +
  geom_text(aes(label = siRNA)) +
  #scale_color_gradientn(colours = brewer.pal(11,"PRGn")) +
  scale_color_viridis_c() +
  ylab("CCND1 expression after KD (SK-N-AS - KELLY)") +
  xlab("relative chromatin oppening \n open in KELLY to the left \n open in SK-N-AS to the right") +
  theme_cowplot()
  
ccnd1kd_norm


library(patchwork)
ggk <- ccnd1kd_norm %>% 
  mutate(oppening = -1*footprint_ttest$val[match(tf, footprint_ttest$name)]) %>% 
  mutate(mes =  MES_activity[match( tf, names(MES_activity))]) %>% 
  arrange(oppening) %>% 
  mutate(pos = 1:n()) %>% 
  ggplot(aes(x = oppening, y = KELLY, color = mes)) +
  geom_point() +
  geom_text(aes(label = tf)) +
  scale_color_viridis_c() +
  theme_cowplot()


ggs <- ccnd1kd_norm %>% 
  mutate(oppening = -1*footprint_ttest$val[match(tf, footprint_ttest$name)]) %>% 
  mutate(mes =  MES_activity[match( tf, names(MES_activity))]) %>% 
  arrange(oppening) %>% 
  mutate(pos = 1:n()) %>% 
  ggplot(aes(x = oppening, y = SKNAS, color = mes)) +
  geom_point() +
  geom_text(aes(label = tf)) +
  scale_color_viridis_c() +
  theme_cowplot()


ggk+ggs


ccnd1kd_norm %>% 
  # mutate(KELLY = log2(KELLY) ) %>%
  # mutate(SKNAS = log2(SKNAS) ) %>%
  mutate(oppening = -1*footprint_ttest$val[match(tf, footprint_ttest$name)]) %>% 
  mutate(mes =  MES_activity[match( tf, names(MES_activity))]) %>% 
  arrange(oppening) %>% 
  mutate(pos = 1:n()) %>% 
  ggplot(aes(x = SKNAS, y = KELLY, color = mes)) +
  geom_point() +
  geom_text(aes(label = tf)) +
  scale_color_viridis_c() +
  theme_cowplot()

