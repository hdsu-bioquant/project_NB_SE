library(tidyverse)
library(readxl)
library(patchwork)

kd_df <- read_xlsx("~/Downloads/Fig 4e_RTqPCR_CCND1 expr_raw data.xlsx", col_names = TRUE, sheet = "Both_long")
kd_df


kd_df %>% 
  ggplot(aes(x = CCND1_exprs, y = siRNA, color = Biological_replicate)) +
  geom_jitter(height = .1) +
  #geom_point() +
  facet_grid(cell_line~., scales = "free") +
  theme_cowplot()
kd_df

kd_df %>% 
  # Mean CCND1 expression for all technical replicates
  group_by(cell_line, Biological_replicate, Plate, siRNA_class, siRNA) %>% 
  summarise(CCND1_exprs = mean(CCND1_exprs)) %>%
  # Select negative: mean both negatives
  mutate(siRNA = if_else(grepl("^neg ", siRNA), "neg", siRNA)) %>% 
  group_by(cell_line, Biological_replicate, Plate, siRNA_class, siRNA) %>% 
  summarise(CCND1_exprs = mean(CCND1_exprs)) %>% 
  ungroup() %>% 
  # Create negative col
  pivot_wider(names_from = "siRNA_class", values_from = "CCND1_exprs") %>% 
  group_by(cell_line, Biological_replicate, Plate) %>% 
  mutate(neg = na.omit(neg)) %>% 
  filter(!is.na(gene)) %>% # remove neg mock and  untr from gene exprs col
  # Normalize by neg
  mutate(CCND1_exprs = gene/neg) %>% 
  
  ggplot(aes(x = CCND1_exprs, y = siRNA, shape = Biological_replicate, color = cell_line, fill = cell_line)) +
  geom_jitter(height = .1, alpha = 0.95, size = 2.5) +
  
  # Error bars
  geom_errorbar(data = . %>% group_by(cell_line, siRNA) %>% 
                  summarise(sd = sd(CCND1_exprs), 
                            CCND1_exprs = mean(CCND1_exprs),
                            Biological_replicate = NA),
                aes(xmin=CCND1_exprs-sd, xmax=CCND1_exprs+sd), width=.2, color = "Black") +
  stat_summary(fun = mean, geom = "errorbar", 
               aes(xmax = ..x.., xmin = ..x.., group = siRNA),
               width = .3, linetype = "solid", color = "Black") +
  
  facet_grid(cell_line~., scales = "free") +

  xlim(c(0, 1.2)) +
  ggtitle("Mean of technical replicates\n Normalized by dividing mean of both controls ") +
  scale_color_manual(values = c("#2FB47C", "#420A68")) +
  scale_fill_manual(values = c("#2FB47C", "#420A68")) +
  scale_shape_manual(values = c(21,22,23,24)) +
  theme_cowplot()

ggsave("kd_figures/Figure_4e_KD_CCND1.pdf", width = 6, height = 4)
  
