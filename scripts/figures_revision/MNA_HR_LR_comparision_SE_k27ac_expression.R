library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)

#-------------------------------------------------------------------------------
DATAPATH  = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/"
outpath   = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB_revision/"

NBexprs   = paste0(DATAPATH, 'analysis/tumor/rnaseq/exprs/tumor_RNAseq_TPM_Matrix_filt_log.RDS')
annoSamp  = paste0(DATAPATH, 'annotation/annotation_tumor_phenoData_NMFassign.RDS')

SEmat     = paste0(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE_SignalScore.RDS')
annoSE    = paste0(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_NMF_assignment.RDS')
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
output1 = paste0(outpath, "results/MNA_HR_MNA_LR_comparision_SE_K27ac_and_target_expression.pdf")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Consensus super enahncer signal matrix and the signature assignments
#-------------------------------------------------------------------------------

SEmat = readRDS(SEmat)
annoSE = readRDS(annoSE)[,c("SE_ID", "target_SYMBOL", "Signature_ChIPseq_MaxExposure")]

annoSamp = readRDS(annoSamp)
annoSamp = annoSamp[annoSamp$avail.ChIPseq,]
annoSamp = annoSamp[,c("ProjectID", "Signature_ChIPseq")]

#Sanity check !!
identical(colnames(SEmat), annoSamp$ProjectID)
#[1] TRUE

#----
# a = split(annoSE, annoSE$target_SYMBOL)
# a = a[sapply(a, nrow) > 1]
# 
# b = a[sapply(a, function(x) length(unique(x$Signature_ChIPseq_MaxExposure))) == 1]
# c = a[sapply(a, function(x) length(unique(x$Signature_ChIPseq_MaxExposure))) > 1]
# rm(a,b,c)
#---

#-------------------------------------------------------------------------------
# NB expression data and subsetting it to samples from CHiPseq cohort only
#-------------------------------------------------------------------------------

# NB expression data
NBexprs = readRDS(NBexprs)
symbols = sapply(strsplit(rownames(NBexprs), "|", fixed=T), function(x)x[2])

#-------------------------------------------------------------------------------
# MNA-HR/LR specific SE H3K27ac signal intensity
#-------------------------------------------------------------------------------

sel = SEmat[rownames(SEmat) %in% annoSE$SE_ID, 
            colnames(SEmat) %in% annoSamp$ProjectID]

#Sanity check !! 
identical(colnames(sel), annoSamp$ProjectID)

sel = sapply(split(data.frame(t(sel)), annoSamp$Signature_ChIPseq), function(x) apply(x, 2, mean))
sel = data.frame(SE_ID = rownames(sel), sel, stringsAsFactors = F, check.names = F)

sel = merge(sel, annoSE, by = "SE_ID")

#-------
dat1 = sel[sel$Signature_ChIPseq_MaxExposure %in% c("MNA-HR", "MNA-LR"),]
dat1 = droplevels(dat1)
dat1$HRvsLR = dat1$`MNA-HR`/dat1$`MNA-LR`
dat1 = dat1[,c(6,1,7,3,4,8)]
colnames(dat1)[1:3] = c("SE", "Gene", "SE_signature")
dat1 = dat1[order(dat1$SE_signature, dat1$HRvsLR, decreasing = T),]

write.table(dat1, paste0(outpath, "results/MNA_HR_LR_Specific_SE_H3K27ac_signal_in_MNA_HR_LR_samples.txt"), row.names = F, quote = F, sep="\t")
rm(dat1)
#-------

sel = melt(sel)
colnames(sel) = c("SE_ID","Gene", "Signature_specific_SE", "Signature_specific_samples", "H3K27ac_Signal")

lab = rep(0, nrow(sel))
lab[sel$Gene == "CCND1" & sel$Signature_specific_SE =="MNA-HR" ] = 1

# a = sel[sel$Signature_specific_SE == "MNA-HR" & sel$Signature_specific_samples == "MNA-HR",]
# a = a[order(a$H3K27ac_Signal),]
# #a = rbind(head(a, 20), tail(a, 20))
# a = tail(a, 10)
# a = apply(a[,1:4], 1, function(x) paste(x, collapse="|"))
# 
# b = sel[sel$Signature_specific_SE == "MNA-HR" & sel$Signature_specific_samples == "MNA-LR",]
# b = b[order(b$H3K27ac_Signal),]
# #b = rbind(head(b, 20), tail(b, 20))
# b = tail(b, 10)
# b = apply(b[,1:4], 1, function(x) paste(x, collapse="|"))
# 
# c = sel[sel$Signature_specific_SE == "MNA-LR" & sel$Signature_specific_samples == "MNA-HR",]
# c = c[order(c$H3K27ac_Signal),]
# #c = rbind(head(c, 20), tail(c, 20))
# c = tail(c, 10)
# c = apply(c[,1:4], 1, function(x) paste(x, collapse="|"))
# 
# d = sel[sel$Signature_specific_SE == "MNA-LR" & sel$Signature_specific_samples == "MNA-LR",]
# d = d[order(d$H3K27ac_Signal),]
# #d = rbind(head(d, 20), tail(d, 20))
# d = tail(d,10)
# d = apply(d[,1:4], 1, function(x) paste(x, collapse="|"))
# 
# lab = apply(sel[,1:4], 1, function(x) paste(x, collapse="|"))
# lab[lab %in% c(a,b,c,d)] = 1
# lab[lab != 1] = 0
# 
# sel$Plot = lab
# rm(a,b,c,d,lab)

sel$Plot = lab
rm(lab)

sel$Signature_specific_SE = factor(paste0(sel$Signature_specific_SE, " signature specific SE"))
sel$H3K27ac_Signal = log2(sel$H3K27ac_Signal)

p1 = ggplot(sel, aes(x = Signature_specific_samples, y = H3K27ac_Signal, color = Signature_specific_samples)) + 
     theme_bw(base_size = 9) + labs(x = "", y = "Summarized H3K27ac signal over SE (log2)") +
     geom_violin(draw_quantiles = c(0.5, 0.9)) +
     #geom_point(position = position_jitterdodge(jitter.width = 0.1, seed = 1),  alpha = 0.25, size = 0.5) +
     scale_color_manual(values=c("#397FB9", "#994F9F", "#50AF49", "#E21E25")) +
     facet_wrap(~ Signature_specific_SE, ncol = 2) + 
     ggsignif::geom_signif(comparisons = list(c(2,3)), textsize = 2) +
     geom_text_repel(data = subset(sel, Plot == 1), aes(label = Gene), size=1.5) + #position  = position_jitterdodge(jitter.width = 0.1, seed = 1)
     theme(panel.grid = element_blank(),
        axis.text = element_text(colour="black"), 
        axis.line = element_blank(),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.5),
        legend.position = "right")

rm(sel)
#-------------------------------------------------------------------------------
# MNA-HR/LR specific SE target gene expression
#-------------------------------------------------------------------------------

# Subsetting the expression matrix and sample annotations to keep common samples
cmn = intersect(annoSamp$ProjectID, colnames(NBexprs))
NBexprs = NBexprs[symbols %in% annoSE$target_SYMBOL, cmn]
annoSamp = annoSamp[annoSamp$ProjectID %in% cmn,]
rm(cmn)

#Sanity check !!
identical(colnames(NBexprs), annoSamp$ProjectID)
#[1] TRUE

sel = sapply(split(data.frame(t(NBexprs)), annoSamp$Signature_ChIPseq), function(x) apply(x, 2, mean))
sel = data.frame(target_SYMBOL = sapply(strsplit(rownames(sel), ".", fixed=T), function(x)x[2]), 
                 sel, stringsAsFactors = F, check.names = F)

sel = merge(sel, annoSE, by = "target_SYMBOL")

#-------
dat1 = sel[sel$Signature_ChIPseq_MaxExposure %in% c("MNA-HR", "MNA-LR"),]
dat1 = droplevels(dat1)
dat1$HRvsLR = dat1$`MNA-HR`/dat1$`MNA-LR`
dat1 = dat1[,c(6,1,7,3,4,8)]
colnames(dat1)[1:3] = c("SE", "Gene", "SE_signature")
dat1 = dat1[order(dat1$SE_signature, dat1$HRvsLR, decreasing = T),]

write.table(dat1, paste0(outpath, "results/MNA_HR_LR_Specific_SE_target_gene_expression_in_MNA_HR_LR_samples.txt"), row.names = F, quote = F, sep="\t")
rm(dat1)
#-------

sel = melt(sel)
colnames(sel) = c("Gene", "SE_ID", "Signature_specific_SE", "Signature_specific_samples", "Expression")

lab = rep(0, nrow(sel))
lab[sel$Gene == "CCND1" & sel$Signature_specific_SE =="MNA-HR" ] = 1
lab[sel$Gene %in% c("RGS5", "GAP43", "MAOA", "DDC", "KRT19", "CD44", "NRCAM", "CADM1", "EBF1", "LMO1", "DPP6", "CREB5", "ZNF423", "CASZ1", "SOX6", "PAX4") & 
    sel$Signature_specific_SE =="MNA-LR" & (sel$Signature_specific_samples == "MNA-HR" | sel$Signature_specific_samples == "MNA-LR") ] = 1

# a = sel[sel$Signature_specific_SE == "MNA-HR" & sel$Signature_specific_samples == "MNA-HR",]
# a = a[order(a$Expression),]
# #a = rbind(head(a, 20), tail(a, 20))
# a = tail(a, 10)
# a = apply(a[,1:4], 1, function(x) paste(x, collapse="|"))
# 
# b = sel[sel$Signature_specific_SE == "MNA-HR" & sel$Signature_specific_samples == "MNA-LR",]
# b = b[order(b$Expression),]
# #b = rbind(head(b, 20), tail(b, 20))
# b = tail(b, 10)
# b = apply(b[,1:4], 1, function(x) paste(x, collapse="|"))
# 
# c = sel[sel$Signature_specific_SE == "MNA-LR" & sel$Signature_specific_samples == "MNA-HR",]
# c = c[order(c$Expression),]
# #c = rbind(head(c, 20), tail(c, 20))
# c = tail(c, 10)
# c = apply(c[,1:4], 1, function(x) paste(x, collapse="|"))
# 
# d = sel[sel$Signature_specific_SE == "MNA-LR" & sel$Signature_specific_samples == "MNA-LR",]
# d = d[order(d$Expression),]
# #d = rbind(head(d, 20), tail(d, 20))
# d = tail(d,10)
# d = apply(d[,1:4], 1, function(x) paste(x, collapse="|"))
# 
# lab = apply(sel[,1:4], 1, function(x) paste(x, collapse="|"))
# lab[lab %in% c(a,b,c,d)] = 1
# lab[lab != 1] = 0
# rm(a,b,c,d,lab)

sel$Plot = lab
rm(lab)

sel$Signature_specific_SE = factor(paste0(sel$Signature_specific_SE, " signature specific SE target genes"))

p2 = ggplot(sel, aes(x = Signature_specific_samples, y = Expression, color = Signature_specific_samples)) + 
    theme_bw(base_size = 9) + labs(x = "", y = "Gene expression (TPM)") +
    geom_violin(draw_quantiles = c(0.5, 0.9)) +
    #geom_point(position = position_jitterdodge(jitter.width = 0.1, seed = 1),  alpha = 0.25, size = 0.5) +
    scale_color_manual(values=c("#397FB9", "#994F9F", "#50AF49", "#E21E25")) +
    facet_wrap(~ Signature_specific_SE, ncol = 2) + 
    ggsignif::geom_signif(comparisons = list(c(2,3)), textsize = 2) +
    geom_text_repel(data = subset(sel, Plot == 1), aes(label = Gene), size=1.5) + #position  = position_jitterdodge(jitter.width = 0.1, seed = 1)
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour="black"), 
          axis.line = element_blank(),
          strip.background = element_blank(),
          axis.ticks.x = element_blank(),
          #axis.text = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.5),
          legend.position = "right")

ggarrange(p1, p2, common.legend = T, ncol = 2) %>% ggexport(filename = output1, width=10, height=5)

rm(sel)


p1$data
p2$data
write_xlsx(list(`Extended Data figure 5e chip` = p1$data,
                `Extended Data figure 5e rna` = p2$data), 
           path = "results/figure_source_data/Extended_Data_figure_5e.xlsx")
