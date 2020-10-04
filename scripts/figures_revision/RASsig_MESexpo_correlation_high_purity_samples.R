library(tidyverse)
library(reshape2)
library(rio)
library(ggplot2)
library(gridExtra)
library(ggrepel)

#-------------------------------------------------------------------------------------------------
DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/"

NBexprs = paste0(DATAPATH, 'analysis/tumor/rnaseq/exprs/tumor_RNAseq_TPM_Matrix_filt_log.RDS')
tumoNMF = paste0(DATAPATH, 'analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS')
NBreg   = paste0(DATAPATH, 'analysis/tumor/ARACNe/network.txt')
NBmut   = paste0(DATAPATH, 'annotation/NB_mutation_matrix.RDS')

outpath = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB_revision/"

rasSigr_Loboda = readRDS(paste0(DATAPATH, 'db/publicGeneSigs/ras_target_genes.RDS'))
rasSigr = import("https://www.cell.com/cms/10.1016/j.celrep.2018.03.046/attachment/7ddfc075-bc33-4b25-832e-729c9653f7a4/mmc2.xlsx")
colnames(rasSigr) = rasSigr[2,]
rasSigr = rasSigr[- c(1:2),]
rasSigr = rasSigr[,c(1,2,3,5)]
rasSigr[,2:ncol(rasSigr)] = apply(rasSigr[,2:ncol(rasSigr)], 2, as.numeric)

tmp = list(exp_rasUp = unique(rasSigr$gene[rasSigr$fold_change > 0.15 & rasSigr$p < 0.001]),
           exp_rasDw = unique(rasSigr$gene[rasSigr$fold_change < -0.15 & rasSigr$p < 0.001]),
           mod_rasUp = unique(rasSigr$gene[rasSigr$classifier_weight > 0]),
           mod_rasDw = unique(rasSigr$gene[rasSigr$classifier_weight < 0]))

rasSigr = tmp
rm(tmp)

## NB expression data
dat = readRDS(NBexprs)
sym = sapply(strsplit(rownames(dat),"|",fixed=T),function(x)x[2])

expo = readRDS(tumoNMF)
identical(colnames(expo), colnames(dat))

## RAS mutation status
NBmuts = readRDS(NBmut)
NBmuts = NBmuts[rownames(NBmuts) %in% c("NRAS", "KRAS", "HRAS", "BRAF", "RAF1", "NF1"),]
NBmuts = NBmuts[, -which(apply(NBmuts, 2, function(x)sum(is.na(x))) == nrow(NBmuts))]

## RAS metagene signature
ras = rasSigr$exp_rasUp
ras1 = rasSigr_Loboda$upreg

## Median expression of RAS target genes in Way et.al and Loboda et.al signatures and their comparison
vals1 = apply(dat,2,function(x){
  x = median(x[sym %in% ras1])
})

vals2 = apply(dat,2,function(x){
  x = median(x[sym %in% ras])
})


if(identical(names(vals1), names(vals2)) & identical(names(vals1), colnames(expo)))
{
  purity = read.csv(paste0(DATAPATH,"annotation/purity_tumor.csv"), stringsAsFactors = F, header=T)
  
  mutInfo = droplevels(na.omit(melt(t(NBmuts))))
  mutInfo = data.frame(PatID = as.character(mutInfo$Var1),
                       Mutation = paste(mutInfo$Var2, mutInfo$value, sep="|"),
                       stringsAsFactors = F)
  
  mutInfo = t(sapply(split(mutInfo, mutInfo$PatID), function(x){
    if(nrow(x) > 1){
      x = data.frame(PatID = unique(x[,1]), 
                     Mutation = paste(x[,2], collapse="|"),
                     stringsAsFactors = F) 
    }
    return(x)
  }))
  
  mutInfo = data.frame(mutInfo, stringsAsFactors = F)
  mutInfo$PatID = unlist(mutInfo$PatID)
  mutInfo$Mutation = unlist(mutInfo$Mutation)
  
  df = data.frame(PatID = names(vals1),
                  Loboda_et_al = as.numeric(vals1),
                  Way_et_al = as.numeric(vals2),
                  MES_Exposure = expo["MES",],
                  stringsAsFactors = F)
  
  sum(df$PatID %in% mutInfo$PatID) == nrow(mutInfo)
  
  df = merge(df, mutInfo, all = T)
  df = df[df$PatID %in% purity$ProjectID[purity$Purity >= 0.5],]
  
  corr1 = cor.test(df$Loboda_et_al, df$MES_Exposure, use="pairwise.complete.obs", method="spearman")
  corr2 = cor.test(df$Way_et_al, df$MES_Exposure, use="pairwise.complete.obs", method="spearman")
  
  p1 = ggplot(df, aes(x = MES_Exposure, y = Loboda_et_al)) + theme_bw(base_size = 9) +
    labs(subtitle = paste(paste0("Samples with purity scores available = ", nrow(purity)),
                          paste0("Samples with high purity scores (>= 50%) = ", nrow(df)),
                          sep="\n"), x="Mesenchymal exposure", y="RAS signature (Loboda et.al)") + 
    geom_point(color=ifelse(is.na(df$Mutation), "#636363", "firebrick"), size = 0.5) + 
    geom_smooth(method = "lm", se = FALSE, colour="grey30") +
    geom_text_repel(data = df[!is.na(df$Mutation),], aes(label = Mutation), size = 3, color = "firebrick") +
    annotate("text", x = 0.0035, y = 3, 
             label = paste(paste0("p=",ifelse(corr1$p.value == 0 , "< 2.2e-16", signif(corr1$p.value,3))),
                          paste0("rho=",round(corr1$estimate,2)), sep="\n"), 
             size=2.8) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour="black"), 
          axis.line = element_blank())
  
  p2 = ggplot(df, aes(x = MES_Exposure, y = Way_et_al)) + theme_bw(base_size = 9) +
    
    labs(subtitle = paste(paste0("Samples with purity scores available = ", nrow(purity)),
                          paste0("Samples with high purity scores (>= 50%) = ", nrow(df)),
                          sep="\n"), x="Mesenchymal exposure", y="RAS signature (Way et.al)") + 
    
    geom_point(color=ifelse(is.na(df$Mutation), "#636363", "firebrick"), size = 0.5) + 
    geom_smooth(method = "lm", se = FALSE, colour="grey30") +
    geom_text_repel(data = df[!is.na(df$Mutation),], aes(label = Mutation), size = 3, color = "firebrick") +
    annotate("text", x = 0.0035, y = 4.25, 
             label = paste(paste0("p=",ifelse(corr2$p.value == 0 , "< 2.2e-16", signif(corr2$p.value,3))),
                           paste0("rho=",round(corr2$estimate,2)), sep="\n"), 
             size=2.8) +
    #ylim(4,5.2)+
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour="black"), 
          axis.line = element_blank())
  
  p = grid.arrange(p1, p2, ncol = 2)
  ggsave(filename = paste0(outpath, "results/RASsig_MESexpo_correlation_high_purity_samples.pdf"), plot = p, width = 7, height = 5)
  
  
}






write_xlsx(list(`Extended Data figure 10g` = p1$data %>% select(PatID, MES_Exposure, Loboda_et_al, Mutation)), 
           path = "results/figure_source_data/Extended_Data_figure_10g.xlsx")
write_xlsx(list(`Extended Data figure 10h` = p2$data %>% select(PatID, MES_Exposure, Way_et_al, Mutation)), 
           path = "results/figure_source_data/Extended_Data_figure_10h.xlsx")
