options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

SEtarget = as.character(args[1])
KDvals   = as.character(args[2])
cellExpr = as.character(args[3])
outpath  = as.character(args[4])

#---------------------------------------------------------------------------------------------------------
# DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"
# 
# SEtarget = paste0(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS')
# KDvals   = paste0(DATAPATH, 'db/DeepMap19Q2/cellKnockdownCERES.RDS')
# cellExpr = paste0(DATAPATH, 'db/DeepMap19Q2/cellExpression.RDS')
# outpath  = paste0(DATAPATH, "results/")
#---------------------------------------------------------------------------------------------------------

library(ggplot2)
library(ggrepel)
library(beeswarm)
library(GenomicRanges)

# Consensus SE target genes
consSEtar = readRDS(SEtarget)
SEtar = unique(mcols(consSEtar)$target_SYMBOL)
SEtar = SEtar[!is.na(SEtar)]
rm(consSEtar)

#SE target genes KD comparison NBcell vs rest

# DepMap CERES data
KDdat = readRDS(KDvals)
diffKD = KDdat$essentiality
class = rep("Others", ncol(diffKD))
identical(colnames(KDdat$essentiality), KDdat$annotation$DepMap_ID)
class[KDdat$annotation$disease_sutype == "neuroblastoma"] = "NB"
class = factor(class, levels = c("NB", "Others"))
diffKD = diffKD[rownames(diffKD) %in% SEtar,]

# Selecting only NB
nbcells.kd = diffKD[,class == "NB"]
medrnk = apply(nbcells.kd, 1, median)

# Diff KD analysis of SE target genes in NB cells vs remaining other cells
diffKD.res  = t(apply(diffKD, 1, function(x){
              w = wilcox.test(x~class)
              med.diff = median(x[class == "NB"]) - median(x[class == "Others"])
              w.pval = w$p.value
              return(c(Diff_NB_vs_Rest = med.diff, pval = w.pval))
              }))

diffKD.res = as.data.frame(diffKD.res)
diffKD.res = data.frame(Gene = rownames(diffKD.res), diffKD.res, FDR = p.adjust(diffKD.res$pval), medianKD_NB = medrnk)

# top hits selection
topHits = rep("No", nrow(diffKD.res))
topHits[diffKD.res$FDR < 0.01 & diffKD.res$medianKD_NB < -0.5 & diffKD.res$Diff_NB_vs_Rest < -0.25] = "Yes"
diffKD.res$topHits = topHits

p = ggplot(diffKD.res, aes(x = Diff_NB_vs_Rest, y = medianKD_NB)) + theme_bw(base_size = 9) + geom_point(aes(color = topHits)) + 
    labs(x = "NB cells vs Other cells (KD difference)", y = "Median KD effect in NB cells") +
    scale_color_manual(values=c("grey80", "black")) +
    #geom_hline(yintercept = -0.5) + geom_vline(xintercept = -0.25) + 
    geom_text_repel(data = subset(diffKD.res, diffKD.res$topHits == "Yes"), aes(label = Gene), size=2.5) +
    theme(panel.grid = element_blank(),
          legend.position = "none")

pdf(paste0(outpath,"figure3/SEtargetGenes_DiffKDprofileNBcellsVsRest.pdf"), width=3.5, height=3.5)
  print(p)
dev.off()

# Plotting the top hits from the above plot in KELLY and SKNAS cells
nbcells.kd.anno = KDdat$annotation[KDdat$annotation$disease_sutype == "neuroblastoma",]
identical(nbcells.kd.anno$DepMap_ID, colnames(nbcells.kd))
colnames(nbcells.kd) = nbcells.kd.anno$stripped_cell_line_name
nbcells.kd = nbcells.kd[,colnames(nbcells.kd) %in% c("KELLY","SKNAS")]
nbcells.kd = nbcells.kd[rownames(nbcells.kd) %in% as.character(diffKD.res$Gene[diffKD.res$topHits == "Yes"]),]
colnames(nbcells.kd) =  c("Kelly (ADRN)", "SK-N-AS (MES)")

pdf(paste0(outpath,"figure3/Kelly_SKNA_KDprofile_topHits.pdf"), width=2, height=2)
par(mar=c(4,3.5,0.25,0.25), mgp=c(2.2,0.5,0), cex=0.8, xaxs="i")
  barplot(as.matrix(t(nbcells.kd)), beside=T, ylab = "CERES KD score", las =2,
          legend.text = TRUE, args.legend = list(x = 15, y = -1, bty="n", cex=0.7, border=NA, x.intersp=0.5))
dev.off()

rm(diffKD, diffKD.res, p, medrnk, class, topHits, nbcells.kd.anno, nbcells.kd)

## CCND1 KD comparison across tumor types
ccnd1.tumorType = as.numeric(KDdat$essentiality[rownames(KDdat$essentiality) == "CCND1",])
ccnd1.tumorType = split(ccnd1.tumorType, KDdat$annotation$disease)
ccnd1.tumorType = ccnd1.tumorType[sapply(ccnd1.tumorType, length) >= 5]
ccnd1.tumorType = ccnd1.tumorType[order(sapply(ccnd1.tumorType, median))]

pdf(paste0(outpath,"sup_figure3/CCND1_KDprofile_across_tumorTypes_cells.pdf"), width=3.7, height=3)
par(mar=c(11,3.5,0.25,0.25), mgp=c(2.2,0.5,0), cex=0.8, xaxs="i")
boxplot(ccnd1.tumorType, boxlty=1, whisklty=1, staplelty=1, medlwd=1, boxwex=0.7, 
        ylim=c(-3,0), ylab = "CCND1 knockdown", xlab= "", las=2, outline=T)
abline(h = -1, lty=2)
dev.off()
rm(ccnd1.tumorType)

## CCND1 expression comparison across tumor types
cellsExp = readRDS(cellExpr)
ccnd1Exp = unlist(cellsExp$expression[rownames(cellsExp$expression) == "CCND1",])
ccnd1Exp = split(ccnd1Exp, cellsExp$annotation$disease)
ccnd1Exp = ccnd1Exp[sapply(ccnd1Exp, length) >= 10]
ccnd1Exp = ccnd1Exp[order(sapply(ccnd1Exp, median))]

pdf(paste0(outpath,"sup_figure3/CCND1_ExpressionProfile_across_tumorTypes_cells.pdf"), width=4.2, height=3.4)
par(mar=c(11,3.5,0.25,0.25), mgp=c(2.2,0.5,0), cex=0.8, xaxs="i")
boxplot(ccnd1Exp, boxlty=1, whisklty=1, staplelty=1, medlwd=1, boxwex=0.7, pch=20, ylim =c(0,12),
        ylab = "CCND1 expression (log2 TPM)", xlab= "", las=2, outline=T)
dev.off()
rm(ccnd1Exp)

#CCND1 expression comparision NBcells vs Rest
ccnd1Exp = t(cellsExp$expression[rownames(cellsExp$expression) == "CCND1",])
identical(rownames(ccnd1Exp), cellsExp$annotation$DepMap_ID)
class = rep("Others", length(ccnd1Exp))
class[cellsExp$annotation$disease_sutype == "neuroblastoma"] = "NBcells"
ccnd1Exp = data.frame(CCND1=ccnd1Exp, Type=class)

pdf(paste0(outpath,"sup_figure3/CCND1_ExpressionProfile_NBcellsVsRest.pdf"), width=1.5, height=2.5)
  par(mar=c(3.5,2.5,0.25,0.25), mgp=c(1.5,0.5,0), cex=0.8, xaxs="i")
  boxplot(CCND1~Type, data=ccnd1Exp, boxlty=1, whisklty=1, staplelty=1, medlwd=1, boxwex=0.7, outline=F,
          ylim=c(-0.2,12),ylab = "CCND1 expression (log2 TPM)", xlab= "",las=2)
  beeswarm(CCND1 ~ Type, data = ccnd1Exp, method="hex", pch=20, cex=0.12,  col="firebrick", add=T)
  text(x = 1, y = -0.3, labels = paste0("n=", sum(ccnd1Exp$Type == "NBcells")), cex=0.75)
  text(x = 2, y = -0.3, labels = paste0("n=", sum(ccnd1Exp$Type == "Others")), cex=0.75)
  text(x = 1.5, y = 12, labels = paste0("p=", signif(wilcox.test(ccnd1Exp$CCND1 ~ ccnd1Exp$Type)$p.value, digits = 2)), cex=0.75)
dev.off()
rm(cellsExp, ccnd1Exp, class, KDdat, SEtar)
