library(tidyverse)
library(cowplot)
ccnd1kd_val <- read_delim("data/cells/CRCsiRNAknockdown/KELLY_SKNAS_validation_long.txt", delim = "\t")

selected_tf <- c("CREB5 #1", "ETS1 #2", "ETV6 #1", "FOSL2 #1", "EBF1 #2", 
                 "GATA3 #1", "GATA2 #1", "MEIS2 #1", "GATA2 #1" )


#------------------------------------------------------------------------------#
#                                 Scatter plot                                 #
#------------------------------------------------------------------------------#

gg_sel <- ccnd1kd_val %>% 
  group_by(siRNA, Cell) %>% 
  filter(n() > 1) %>% 
  ungroup() %>% 
  mutate(siRNA = sub(" ", " #", siRNA)) %>% 
  filter(siRNA %in% selected_tf) %>% 
  #as.data.frame()
  ggplot(aes(x = siRNA, y = Expression, color = Cell)) +
  geom_point() +
  geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=.2,
                position=position_dodge(0.05))+
  #facet_wrap(.~Cell, scales = "free") +
  facet_grid(Cell~., scales = "free") +
  coord_flip() +
  scale_color_manual(values = c("#2FB47C", "#420A68")) +
  theme_cowplot()


gg_all <- ccnd1kd_val %>% 
  group_by(siRNA, Cell) %>% 
  filter(n() > 1) %>% 
  ungroup() %>% 
  mutate(siRNA = sub(" ", " #", siRNA)) %>% 
  #filter(siRNA %in% selected_tf) %>% 
  #as.data.frame()
  ggplot(aes(x = siRNA, y = Expression, color = Cell)) +
  geom_point() +
  geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=.2,
                position=position_dodge(0.05))+
  #facet_wrap(.~Cell, scales = "free") +
  facet_grid(Cell~., scales = "free") +
  coord_flip() +
  scale_color_manual(values = c("#2FB47C", "#420A68")) +
  theme_cowplot()


gg_all + gg_sel + plot_layout(ncol = 1, heights = c(0.7, 0.3))
ggsave("results/figures_revision/figure_CCND1_KD_Validation.pdf", width = 5, height = 8)






#------------------------------------------------------------------------------#
#                                 Bar plot                                     #
#------------------------------------------------------------------------------#

#ccnd1kd_val %>% 
gg_sel <- ccnd1kd_val %>% 
  dplyr::select(-rep) %>% 
  group_by(siRNA, Cell) %>% 
  filter(n() > 1) %>% 
  mutate(var = sd^2) %>% 
  summarise_all(mean) %>% 
  #mutate(var = log2(var)) %>%
  mutate(sd = sqrt(var)) %>% 
  ungroup() %>% 
  mutate(siRNA = sub(" ", " #", siRNA)) %>% 
  filter(siRNA %in% selected_tf) %>% 
  
  mutate(Expression = log2(Expression)) %>%
  #mutate(sd = log2(sd)) %>%
  
  #as.data.frame()
  ggplot(aes(x = siRNA, y = Expression, color = Cell)) +
  geom_point() +
  #geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=.2,
                position=position_dodge(0.05))+
  #facet_wrap(.~Cell, scales = "free") +
  facet_grid(Cell~., scales = "free") +
  coord_flip() +
  scale_color_manual(values = c("#2FB47C", "#420A68")) +
  theme_cowplot()


gg_all <- ccnd1kd_val %>% 
  dplyr::select(-rep) %>% 
  group_by(siRNA, Cell) %>% 
  filter(n() > 1) %>% 
  mutate(var = sd^2) %>% 
  summarise_all(mean) %>% 
  mutate(sd = sqrt(var)) %>% 
  
  ungroup() %>% 
  mutate(siRNA = sub(" ", " #", siRNA)) %>% 
  #filter(siRNA %in% selected_tf) %>% 
  #as.data.frame()
  ggplot(aes(x = siRNA, y = Expression, color = Cell)) +
  geom_point() +
  #geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=.2,
                position=position_dodge(0.05))+
  #facet_wrap(.~Cell, scales = "free") +
  facet_grid(Cell~., scales = "free") +
  coord_flip() +
  scale_color_manual(values = c("#2FB47C", "#420A68")) +
  theme_cowplot()


gg_all + gg_sel + plot_layout(ncol = 1, heights = c(0.7, 0.3))
ggsave("results/figures_revision/figure_CCND1_KD_Validation_mean.pdf", width = 5, height = 8)





gg_sel <- ccnd1kd_val %>% 
  dplyr::select(-rep) %>% 
  group_by(siRNA, Cell) %>% 
  filter(n() > 1) %>% 
  mutate(var = sd^2) %>% 
  summarise_all(mean) %>% 
  #mutate(var = log2(var)) %>%
  mutate(sd = sqrt(var)) %>% 
  ungroup() %>% 
  mutate(siRNA = sub(" ", " #", siRNA)) %>% 
  filter(siRNA %in% selected_tf) %>% 
  
  #mutate(Expression = log2(Expression)) %>%
  #mutate(sd = log2(sd)) %>%
  
  #as.data.frame()
  ggplot(aes(x = siRNA, y = Expression, color = Cell)) +
  geom_point(alpha = 0.99) +
  #geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=.2,
                position=position_dodge(0.05))+
  #facet_wrap(.~Cell, scales = "free") +
  facet_grid(Cell~., scales = "free") +
  coord_flip() +
  scale_y_continuous(breaks=c(0 , 0.25, 0.5, 0.75, 1), limits = c(0,1.1)) +
  scale_color_manual(values = c("#2FB47C", "#420A68")) +
  theme_cowplot()
gg_sel 
ggsave("results/figures_revision/figure_CCND1_KD_Validation_mean_selected.pdf", width = 5, height = 3.5)



ccnd1kd_val %>% 
  dplyr::select(-rep) %>% 
  group_by(siRNA, Cell) %>% 
  filter(n() > 1) %>% 
  mutate(var = sd^2) %>% 
  summarise_all(mean) %>% 
  #mutate(var = log2(var)) %>%
  mutate(sd = sqrt(var)) %>% 
  ungroup() %>% 
  mutate(siRNA = sub(" ", " #", siRNA)) %>% 
  filter(siRNA %in% selected_tf)

