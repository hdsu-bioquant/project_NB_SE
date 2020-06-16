setwd("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB")
library(patchwork)
library(tidyverse)
library(cowplot)
library(ggpubr)
ccnd1kd <- read.delim("data/cells/CRCsiRNAknockdown/CCND1_Regulator_siRNA_knockdown.txt", 
                      header = TRUE, stringsAsFactors = FALSE, dec = ",")

ccnd1kd <- ccnd1kd %>% 
  mutate(Kelly_rep1 = log2(Kelly_rep1) ) %>% 
  mutate(Kelly_rep2 = log2(Kelly_rep2) ) %>% 
  mutate(SKNAS_rep1 = log2(SKNAS_rep1) ) %>% 
  mutate(SKNAS_rep2 = log2(SKNAS_rep2) )

ccnd1kd$KELLY_sd <- apply(ccnd1kd, 1, function(x) sd(c(x["Kelly_rep1"], x["Kelly_rep2"]))) 
ccnd1kd$SKNAS_sd <- apply(ccnd1kd, 1, function(x) sd(c(x["SKNAS_rep1"], x["SKNAS_rep2"]))) 



ccnd1kd <- ccnd1kd %>% 
  mutate(KELLY = (Kelly_rep1 + Kelly_rep2)/2) %>% 
  mutate(SKNAS = (SKNAS_rep1 + SKNAS_rep2)/2) %>% 
  select(siRNA, Selected, KELLY, SKNAS, KELLY_sd, SKNAS_sd)



#pivot_longer(cols = -c("siRNA", "Selected"), names_to = "Replicate", values_to = "Expression") %>% 
#group_by(Replicate) %>% 
#summarise(Expression = mean(Expression))

#selected_tf <- c("CREB5_1", "ETS1_2", "ETV6_1", "FOSL2_1", "EBF1_2", "GATA3_1", "GATA2_1")
selected_tf <- c("CREB5 #1", "ETS1 #2", "ETV6 #1", "FOSL2 #1", "EBF1 #2", 
                 "GATA3 #1", "GATA2 #1", "MEIS2 #1", "GATA2 #1" )



validated_tf <- c("CREB5 #1" = "S-K-N-AS & KELLY", 
                  "ETS1 #2"  = "S-K-N-AS", 
                  "ETV6 #1"  = "S-K-N-AS & KELLY", 
                  "FOSL2 #1" = "S-K-N-AS",  
                  "EBF1 #2"  = "S-K-N-AS", 
                  "GATA3 #1" = "KELLY", 
                  "MEIS2 #1" = "KELLY", 
                  "GATA2 #1" = "KELLY" )
offset <- 0.7


ccnd1kd %>% 
  # mutate(KELLY = log2(KELLY) ) %>% 
  # mutate(SKNAS = log2(SKNAS) ) %>% 
  
  #filter(siRNA %in% selected_tf) %>% 
  
  mutate(tf = substr(siRNA, 1, nchar(siRNA)-2)) %>% 
  mutate(tf = if_else(tf == "MEIS", "MEIS2", tf)) %>% 
  mutate(siRNA_corr = substr(siRNA, nchar(siRNA), nchar(siRNA))) %>% 
  mutate(siRNA_corr = paste0(tf, " #", siRNA_corr)) %>%
  mutate(Validated = siRNA_corr %in% selected_tf) %>% 
  mutate(Validated = validated_tf[match(siRNA_corr, names(validated_tf))]) %>% 
  mutate(Validated = if_else(is.na(Validated), "No", Validated)) %>% 
  mutate(Validated = factor(Validated, levels = c("S-K-N-AS & KELLY", 
                                                  "KELLY", "S-K-N-AS","No"))) %>% 
  mutate(siRNA_corr = sub(" #", "-", siRNA_corr)) %>% 
  # filter(Validated)
  #as.data.frame()
  
  mutate(Selected = if_else(tf == "EBF1", 1L , Selected)) %>% 
  #filter(Selected == 1) %>% 
  
  
  group_by(tf) %>% 
  # Keep only siRNA with effect in the same direction
  mutate(KELLY_direc = sum(KELLY > 0)) %>% 
  mutate(KELLY_direc = if_else(KELLY_direc == n() | KELLY_direc == 0 , TRUE, FALSE)) %>% 
  mutate(SKNAS_direc = sum(SKNAS > 0)) %>% 
  mutate(SKNAS_direc = if_else(SKNAS_direc == n() | SKNAS_direc == 0 , TRUE, FALSE)) %>% 
  #filter(KELLY_direc & SKNAS_direc ) %>% #do filter
  ungroup() %>% 
  
  mutate(label_side = if_else(SKNAS > median(SKNAS), "r", "l")) %>% 
  mutate(midpoint = mean(c(min(SKNAS), max(SKNAS)))) %>% 
  arrange(-SKNAS) %>% 
  group_by(label_side) %>% 
  
  #mutate(labelpos = rev(seq(min(SKNAS), max(SKNAS), length.out = n()))) %>% 
  mutate(labelpos = if_else(label_side == "r",
                            rev(seq(min(midpoint)-offset, max(SKNAS), length.out = n())),
                            rev(seq(min(SKNAS), max(midpoint)+offset, length.out = n())))) %>% 
  
  
  ggplot(aes(x = KELLY, y = SKNAS)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline()+
  geom_point(aes(color = Validated), size = 2.5, alpha = 0.98) +
  
  #stat_smooth(method = "lm", se = FALSE, color = "black") +
  # stat_regline_equation(label.y = 1.5,
  #                       aes(label =  paste(..rr.label..))) +
  
  stat_cor(method = "pearson", label.x = -2, label.y = 1.5) +
  # geom_errorbar(aes(ymin=SKNAS-SKNAS_sd, ymax=SKNAS+SKNAS_sd),
  #               width=.2,
  #               position=position_dodge(0.05)) +
  # 
  # geom_errorbar(aes(xmin=KELLY-KELLY_sd, xmax=KELLY+KELLY_sd),
  #               width=.2,
  #               position=position_dodge(0.05)) +
  
  
  geom_text(data = function(x){x %>% filter(label_side == "r")},
            mapping = aes(x = max(KELLY) + 0.2, y = labelpos, label = siRNA_corr),
            parse = TRUE, hjust = 0, color = "black") +
  geom_text(data = function(x){x %>% filter(label_side == "l")},
            mapping = aes(x = min(KELLY) - 2, y = labelpos, label = siRNA_corr),
            parse = TRUE, hjust = 0, color = "black") +
  
  geom_segment(data = function(x){x %>% filter(label_side == "r")},
               mapping = aes(x = KELLY, xend = max(KELLY) + 0.1, y = SKNAS, yend = labelpos),
               color = "grey70", size = 0.1) +
  geom_segment(data = function(x){x %>% filter(label_side == "l")},
               mapping = aes(x = KELLY, xend = min(KELLY) - 0.1, y = SKNAS, yend = labelpos),
               color = "grey70", size = 0.1) +
  xlim(c(-3.5,4.5)) +
  
  #
  scale_color_manual(values = c("Red", "#2FB47C", "#420A68", "Black") ) +
  theme_cowplot()
ggsave("results/figures_revision/figure_CCND1_KD_Screening_all.pdf", width = 10, height = 6.5)
#









ccnd1kd %>% 
  # mutate(KELLY = log2(KELLY) ) %>% 
  # mutate(SKNAS = log2(SKNAS) ) %>% 
  
  #filter(siRNA %in% selected_tf) %>% 
  
  mutate(tf = substr(siRNA, 1, nchar(siRNA)-2)) %>% 
  mutate(tf = if_else(tf == "MEIS", "MEIS2", tf)) %>% 
  mutate(siRNA_corr = substr(siRNA, nchar(siRNA), nchar(siRNA))) %>% 
  mutate(siRNA_corr = paste0(tf, " #", siRNA_corr)) %>%
  mutate(Validated = siRNA_corr %in% selected_tf) %>% 
  mutate(siRNA_corr = sub(" #", "-", siRNA_corr)) %>% 
  # filter(Validated)
  #as.data.frame()
  
  mutate(Selected = if_else(tf == "EBF1", 1L , Selected)) %>% 
  filter(Selected == 1) %>% 
  
  
  group_by(tf) %>% 
  # Keep only siRNA with effect in the same direction
  mutate(KELLY_direc = sum(KELLY > 0)) %>% 
  mutate(KELLY_direc = if_else(KELLY_direc == n() | KELLY_direc == 0 , TRUE, FALSE)) %>% 
  mutate(SKNAS_direc = sum(SKNAS > 0)) %>% 
  mutate(SKNAS_direc = if_else(SKNAS_direc == n() | SKNAS_direc == 0 , TRUE, FALSE)) %>% 
  #filter(KELLY_direc & SKNAS_direc ) %>% #do filter
  ungroup() %>% 
  
  mutate(label_side = if_else(SKNAS > median(SKNAS), "r", "l")) %>% 
  mutate(midpoint = mean(c(min(SKNAS), max(SKNAS)))) %>% 
  arrange(-SKNAS) %>% 
  group_by(label_side) %>% 
  
  #mutate(labelpos = rev(seq(min(SKNAS), max(SKNAS), length.out = n()))) %>% 
  mutate(labelpos = if_else(label_side == "r",
                            rev(seq(min(midpoint), max(SKNAS), length.out = n())),
                            rev(seq(min(SKNAS), max(midpoint), length.out = n())))) %>% 
  
  
  ggplot(aes(x = KELLY, y = SKNAS, color = Validated)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline()+
  geom_point(size = 2.5, alpha = 0.98) +
  
  # geom_errorbar(aes(ymin=SKNAS-SKNAS_sd, ymax=SKNAS+SKNAS_sd),
  #               width=.2,
  #               position=position_dodge(0.05)) +
  # 
  # geom_errorbar(aes(xmin=KELLY-KELLY_sd, xmax=KELLY+KELLY_sd),
  #               width=.2,
  #               position=position_dodge(0.05)) +
  
  
  geom_text(data = function(x){x %>% filter(label_side == "r")},
            mapping = aes(x = max(KELLY) + 0.2, y = labelpos, label = siRNA_corr),
            parse = TRUE, hjust = 0, color = "black") +
  geom_text(data = function(x){x %>% filter(label_side == "l")},
            mapping = aes(x = min(KELLY) - 1, y = labelpos, label = siRNA_corr),
            parse = TRUE, hjust = 0, color = "black") +
  
  geom_segment(data = function(x){x %>% filter(label_side == "r")},
               mapping = aes(x = KELLY, xend = max(KELLY) + 0.1, y = SKNAS, yend = labelpos),
               color = "black", size = 0.1) +
  geom_segment(data = function(x){x %>% filter(label_side == "l")},
               mapping = aes(x = KELLY, xend = min(KELLY) - 0.1, y = SKNAS, yend = labelpos),
               color = "black", size = 0.1) +
  xlim(c(-2.5,3)) +
  
  #
  scale_color_manual(values = c("Black", "Red")) +
  theme_cowplot()
ggsave("results/figures_revision/figure_CCND1_KD_Screening_selected.pdf", width = 8, height = 6)
#



















ccnd1kd_summary <- ccnd1kd %>% 
  # mutate(KELLY = log2(KELLY) ) %>% 
  # mutate(SKNAS = log2(SKNAS) ) %>% 
  
  #filter(siRNA %in% selected_tf) %>% 
  
  mutate(tf = substr(siRNA, 1, nchar(siRNA)-2)) %>% 
  mutate(tf = if_else(tf == "MEIS", "MEIS2", tf)) %>% 
  mutate(siRNA_corr = substr(siRNA, nchar(siRNA), nchar(siRNA))) %>% 
  mutate(siRNA_corr = paste0(tf, " #", siRNA_corr)) %>%
  mutate(Validated = siRNA_corr %in% selected_tf) %>% 
  mutate(siRNA_corr = sub(" #", "-", siRNA_corr)) 
unique(ccnd1kd_summary$tf)
length(unique(ccnd1kd_summary$tf))
dim(ccnd1kd_summary)


cor(ccnd1kd$KELLY, ccnd1kd$SKNAS)




SE_annot <- readRDS("analysis/tumor/SE_annot/tumor_consensusSE_target_annotation_df.RDS")


library(rtracklayer)
SE <- readRDS("analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS")

foot_ranges_list <- lapply(c(
  "SK-N-AS" = "data/cells/atacseq/footprint/SK-N-AS_footprints.bed",
  "CLB-GA" = "data/cells/atacseq/footprint/CLB-GA_footprints.bed",
  "KELLY" = "data/cells/atacseq/footprint/KELLY_footprints.bed"), function(path){
  foot_ranges <- read.table(path, stringsAsFactors=FALSE)
  #highlight_ranges <- read.table("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/src/project_NB_SE/scripts/figure3/figure3g_highlight_ranges.bed", stringsAsFactors=FALSE)
  colnames(foot_ranges) <- c("chr", "start", "end", "id")
  foot_ranges <- makeGRangesFromDataFrame(foot_ranges, keep.extra.column = TRUE)
  foot_ranges <- foot_ranges[foot_ranges$id %in% ccnd1kd_summary$tf]
  foot_ranges$CellID <- sub("_footprints.bed", "", basename(path))
  #length(unique(foot_ranges$id))
  foot_ranges
  
})
# foot_ranges_list
# x <- as.data.frame(foot_ranges_list$`SK-N-AS`)

foot_gg_list <- lapply(foot_ranges_list, function(foot_ranges){
  foot_df <- bind_rows(lapply(c("PHOX2A", "PHOX2B", "LDB1", "GATA3", "ISL1", "IRS2", "CCND1"), function(tfid){
    tfse <- SE[SE$target_SYMBOL %in% tfid]
    foottf <- foot_ranges
    foottf$SE_target <- tfid
    foottf <- subsetByOverlaps(foottf, tfse)
    #print(foottf)
    as.data.frame(foottf)
  })) 
  
  foot_df %>% 
    ggplot(aes(x = SE_target, fill = id, label = id))+
    geom_bar() +
    #geom_text(size = 3, position = position_stack(vjust = 0.5)) +
    ylab(paste0("Footprints count in ", unique(foot_ranges$CellID), " inside SE")) +
    xlab("SE TARGET") +
    theme_cowplot()
  
})

foot_gg_list$KELLY + foot_gg_list$`SK-N-AS`


ccnd1kd_summary$tf[!ccnd1kd_summary$tf %in% foot_df_both$id ]

foot_df_both <- bind_rows(lapply(foot_ranges_list, function(foot_ranges){
  foot_df <- bind_rows(lapply(c("PHOX2A", "PHOX2B", "LDB1", "GATA3", "ISL1", "IRS2", "CCND1"), function(tfid){
    tfse <- SE[SE$target_SYMBOL %in% tfid]
    foottf <- foot_ranges
    foottf$SE_target <- tfid
    foottf <- subsetByOverlaps(foottf, tfse)
    #print(foottf)
    as.data.frame(foottf)
  }))
  print(ccnd1kd_summary$tf[!ccnd1kd_summary$tf %in% foot_df$id ])
  foot_df
}))

foot_df_both %>%
  group_by(SE_target, CellID, id) %>% 
  summarise(Count = n()) %>% 
  arrange(id) %>% 
  group_by(id) %>% 
  mutate(idnum = group_indices()) %>% 
  ungroup() %>% 
  mutate(id = paste0(idnum, ". ", id)) %>% 
  mutate(id = factor(id, levels = unique(id))) %>% 
  ggplot(aes(x = SE_target, y = Count, fill = id, label = idnum))+
  geom_bar(stat = "identity") +
  facet_grid(.~CellID) +
  ylab("Footprints count inside SEs") +
  xlab("SE TARGET") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90))
#

ggsave("results/figures_revision/figure_footprint_count_KD_TFs.pdf", width = 10, height = 11)
# 
# 
# x <- as.data.frame(SE)
# SE[SE$target_SYMBOL == "PHOX2B"]
# 
# 
# foot_gg_list$`SK-N-AS` + foot_gg_list$KELLY
# 
# 
# 
# 
# 
# 
# foot_ranges <- read.table("data/cells/atacseq/footprint/KELLY_footprints.bed", stringsAsFactors=FALSE)
# #highlight_ranges <- read.table("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/src/project_NB_SE/scripts/figure3/figure3g_highlight_ranges.bed", stringsAsFactors=FALSE)
# colnames(foot_ranges) <- c("chr", "start", "end", "id")
# foot_ranges <- makeGRangesFromDataFrame(foot_ranges, keep.extra.column = TRUE)
# foot_ranges <- foot_ranges[foot_ranges$id %in% x$tf]
# length(unique(foot_ranges$id))
# foot_ranges
# 
# 
# foot_df <- bind_rows(lapply(c("PHOX2A", "LDB1", "GATA3", "ISL1", "IRS2", "CCND1"), function(tfid){
#   tfse <- SE[SE$target_SYMBOL %in% tfid]
#   foottf <- foot_ranges
#   foottf$SE_target <- tfid
#   foottf <- subsetByOverlaps(foottf, tfse)
#   as.data.frame(foottf)
# }))
# foot_df %>% 
#   ggplot(aes(x = SE_target, fill = id))+
#   geom_bar() +
#   ylab("Footprints count in KELLY inside SE") +
#   xlab("SE TARGET") +
#   theme_cowplot()
# 
# 
# 
# 
# 
# 
# 
# 
