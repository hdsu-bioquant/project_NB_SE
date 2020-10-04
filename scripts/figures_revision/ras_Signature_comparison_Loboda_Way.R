library(reshape2)
library(readxl)
library(rio)
library(RColorBrewer)
library(parallel)
library(gplots)
library(tidyverse)
library(writexl)
#-------------------------------------------------------------------------------------------------
DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/"

NBexprs = paste0(DATAPATH, 'analysis/tumor/rnaseq/exprs/tumor_RNAseq_TPM_Matrix_filt_log.RDS')
tumoNMF = paste0(DATAPATH, 'analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_Hmatrix_wnorm.RDS')
NBreg   = paste0(DATAPATH, 'analysis/tumor/ARACNe/network.txt')
NBmut   = paste0(DATAPATH, 'annotation/NB_mutation_matrix.RDS')

outpath = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB_revision/results/"

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

#-------------------------------------------------------------------------------------------------

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
#---------------------------------
# Permutation testing
#---------------------------------

## Median expression of RAS target genes using Way et.al signature genes
medexp_way = apply(dat,2,function(x){
  x = median(x[sym %in% ras])
})

## Correlation between median expression of RAS target genes using Way et.al signature genes and Mesenchymal exposure
obs.cor = round(cor(medexp_way, expo[2,], method="spearman"),3)

## Calculating the correlation between Mesenchymal exposure and 1000 random gene sets of length equal to Way et.al
rho_random = mclapply(1:1000, function(i){
    set.seed(i)
    random = sample(sym, length(ras))
    
    ## Median expression of random geneset
    medexp = apply(dat,2,function(x){
      x = median(x[sym %in% random])
    })
    
    corr = round(cor(medexp, expo[2,], method="spearman"),3)
    return(corr)
  }, mc.cores = 60)

rho_random = unlist(rho_random)

## Calculating the permutation p value
perm_pval = sum(rho_random > obs.cor)/length(rho_random)
if(perm_pval == 0){
  perm_pval = "p<0.001"
}else{
  perm_pval = paste0("p=", perm_pval)
}

## Median expression of RAS target genes in Way et.al and Loboda et.al signatures and their comparison
vals1 = apply(dat,2,function(x){
  x = median(x[sym %in% ras1])
})

vals2 = apply(dat,2,function(x){
  x = median(x[sym %in% ras])
})

#----------
# RAS plot
#----------

pdf(paste0(outpath,"RAS_validation.pdf"), width=5.5, height=4)

    colGr <- colorRampPalette(brewer.pal(9, "YlOrRd"))(10)[as.numeric(cut(expo[2,],breaks = 10))]
    
    layout(matrix(c(1:4), ncol=2, byrow=T), width=1, height=1)
    par(cex=0.7) # xaxs="i", yaxs="i",
    
    #Plot 1 - Correlation Way et.al RAS signature and MES exposure
    par(mar = c(3.5,3.5,0.5,0.5), mgp=c(1.5,0.5,0))
    plot(medexp_way, expo[2,], las=2, xlim=c(2,6), ylim=c(-100,4000), pch=20, col = colGr,
         xlab = "RAS signature gene expression (Diff exp genes up)", ylab = "Mesenchymal exposure")
    legend("topleft", legend = paste0("rho=", obs.cor), bty="n")
    
    #Plot 2 - Test significance of Way et.al correlation value by permutation testing
    par(mar = c(4,3.5,0.5,0.5), mgp=c(2.5,0.5,0), yaxs="i")
    plot(density(rho_random, na.rm=T), xlim = c(-1,1), main="", 
         xlab = paste("Spearman correlations of MES exposure", "with random genes sets (n=1000)", sep="\n"))
    abline(v = obs.cor, lty=2, col="firebrick")
    legend("topleft", legend = perm_pval, bty="n")
    legend("topright", legend = paste0("observed rho=", obs.cor), bty="n", text.col = "firebrick")
    
    #Plot 3 - sSimilarity between Way et.al and Loboda et.al
    sig_cor = round(cor(vals1, vals2, method="spearman"), 3)
    par(mar = c(3,3,0.5,0.5), mgp=c(1.5,0.5,0))
    plot(x = vals1, y = vals2, pch=20, xlim = c(2.5,6.5), ylim = c(1.5,5.5), col = colGr,
         xlab = "RAS activity (NB) - Loboda et.al. BMC 2010", 
         ylab = "RAS activity (NB) - Way et.al. Cell.Rep. 2018")
    abline(lm(vals2 ~ vals1))
    abline(v = median(vals1), h = median(vals2), lty=2)
    legend("bottomright", legend = paste0("rho=", sig_cor), bty="n")
    
    xmark = vals1[which(names(vals1) %in% colnames(NBmuts))]
    ymark = vals2[which(names(vals2) %in% colnames(NBmuts))]
    points(x = xmark, y = ymark, pch = 4, cex=1, col = "black")
    rm(xmark, ymark)
    
    #Plot4 - Intersection between Way et.al and Loboda et.al
    par(mar = c(0,0,0,0))
    venn(list('Way et.al Cell Rep. 2018' = rasSigr$exp_rasUp, 'Loboda et.al BMC 2010' = rasSigr_Loboda$upreg))
    
dev.off()





write_xlsx(list(`Extended Data figure 9c` = tibble(ProjectID = names(medexp_way), RAS_signature_expression = medexp_way, Mesenchymal_score = expo[2,])), 
           path = "results/figure_source_data/Extended_Data_figure_9c.xlsx")

write_xlsx(list(`Extended Data figure 9d` = tibble(Assay = paste0("Assay", 1:1000), rho_random = rho_random)), 
           path = "results/figure_source_data/Extended_Data_figure_9d.xlsx")

write_xlsx(list(`Extended Data figure 9e` = tibble(ProjectID = names(vals1), RAS_activity_Loboda = vals1, RAS_activity_Way = vals2)), 
           path = "results/figure_source_data/Extended_Data_figure_9e.xlsx")

length(rasSigr$exp_rasUp)
length(rasSigr_Loboda$upreg)
length(c(rasSigr_Loboda$upreg, rep(NA, length(rasSigr$exp_rasUp)-length(rasSigr_Loboda$upreg))))

write_xlsx(list(`Extended Data figure 9f` = tibble('Way et.al Cell Rep. 2018' = rasSigr$exp_rasUp, 'Loboda et.al BMC 2010' = c(rasSigr_Loboda$upreg, rep(NA, length(rasSigr$exp_rasUp)-length(rasSigr_Loboda$upreg))))), 
           path = "results/figure_source_data/Extended_Data_figure_9f.xlsx")

