annot_t <- readRDS("annotation/annotation_tumor_phenoData_NMFassign.RDS")

#Hchip_t  <- as.data.frame(t(readRDS("analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K4_Hmatrix_wnorm.RDS"))) %>% 
Hchip_t  <- as.data.frame(t(readRDS("analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS"))) %>% 
  #Hchip_t  <- as.data.frame(t(readRDS(params$Hchip_t))) %>% 
  rename_all(paste0, "_ChIP-seq") %>% 
  rownames_to_column("ProjectID") 

gg_chip <- Hchip_t %>% 
  left_join(annot_t, by = "ProjectID") %>% 
  mutate(Signature_ChIPseq = if_else(Signature_ChIPseq %in% c("MYCN", "MES"),
                                    "Other", Signature_ChIPseq)) %>% 
  mutate(Signature_ChIPseq = factor(Signature_ChIPseq, levels = c("MNA-LR" ,"MNA-HR", "Other"))) %>% 
  #head()  
  ggplot(aes(x = `MNA-LR_ChIP-seq`, `MNA-HR_ChIP-seq`, color = Signature_ChIPseq)) +
  geom_point(size = 1, alpha = 0.9) +
  scale_color_manual(values = c("firebrick", "#00afaf", "grey80")) +
  theme_cowplot()
gg_chip

#Hrna_tt  <- as.data.frame(t(readRDS("analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_Hmatrix_wnorm.RDS"))) %>% 
Hrna_tt  <- as.data.frame(t(readRDS("analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS"))) %>% 
  #Hrna_tt  <- as.data.frame(t(readRDS(params$Hrna_tt))) %>% 
  rename_all(paste0, "_RNA-seq") %>% 
  rownames_to_column("ProjectID") 



gg_rna <- Hrna_tt %>% 
  left_join(annot_t, by = "ProjectID") %>% 
  mutate(Signature_RNAseq = if_else(Signature_RNAseq %in% c("MYCN", "MES"),
                                    "Other", Signature_RNAseq)) %>% 
  mutate(Signature_RNAseq = factor(Signature_RNAseq, levels = c("MNA-LR" ,"MNA-HR", "Other"))) %>% 
  #Hrna_tt %>% 
  #head() %>% 
  ggplot(aes(x = `MNA-LR_RNA-seq`, `MNA-HR_RNA-seq`, color = Signature_RNAseq)) +
  geom_point(size = 1, alpha = 0.9) +
  scale_color_manual(values = c("firebrick", "#00afaf", "grey80")) +
  theme_cowplot()
gg_rna


gg_chip + gg_rna
ggsave("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB_highPurity/results/comparison_highPurity/MNA-LR_vs_MNA-HR.pdf", width=10, height=3)


sd <- gg_chip$data %>% 
  select(ProjectID, `MNA-LR_ChIP-seq`, `MNA-HR_ChIP-seq`, Signature_ChIPseq)
write_xlsx(list(`Extended Data figure 5c` = sd), 
           path = "results/figure_source_data/Extended_Data_figure_5c.xlsx")

sd <- gg_rna$data %>% 
  select(ProjectID, `MNA-LR_RNA-seq`, `MNA-HR_RNA-seq`, Signature_RNAseq)
write_xlsx(list(`Extended Data figure 5d` = sd), 
           path = "results/figure_source_data/Extended_Data_figure_5d.xlsx")



