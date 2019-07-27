options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

NBexprs  = as.character(args[1])
rasSigr  = as.character(args[2])
NBreg    = as.character(args[3])
crcList  = as.character(args[4])
mesTFact = as.character(args[5])
outpath1 = as.character(args[6])
outpath2 = as.character(args[7])

#-------------------------------------------------------------------------------------------------
# DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"
# 
# NBexprs  = paste0(DATAPATH, 'data/tumor/rnaseq/exprs/tumor_RNAseq_Counts_Matrix.RDS')
# rasSigr  = paste0(DATAPATH, 'db/publicGeneSigs/ras_target_genes.RDS')
# NBreg    = paste0(DATAPATH, 'analysis/tumor/ARACNe/network.txt')
# crcList  = paste0(DATAPATH, 'results/supptables/crcTF_fractionObserved_combined.txt')
# mesTFact = paste0(DATAPATH, "analysis/tumor/VIPER/MES_TFactivity.RDS")
# outpath1 = paste0(DATAPATH, "analysis/")
# outpath2 = paste0(DATAPATH, "results/")
#-------------------------------------------------------------------------------------------------

library(readxl)
library(edgeR)
library(limma)
library(gridExtra)
library(ggrepel)
library(statmod)
library(WriteXLS)
library(fgsea)

# Matched Primary and relapse pairs
pri = c("NSP032-PT01", "NSP090-PT01", "NSP460-PT01")
rel = c("NSP032-RT03", "NSP090-RT01", "NSP460-RM01")

# Load the RNAseq dataset (counts) and subset for matched Primary vs Relapse pairs
nbcnt = readRDS(NBexprs)

a = nbcnt[,colnames(nbcnt) %in% pri] 
a = a[,match(pri, colnames(a))]

b = nbcnt[,colnames(nbcnt) %in% rel]
b = b[,match(rel, colnames(b))]

infdat = cbind(a,b)
infcli = data.frame(Samples=colnames(infdat),Type=c(rep("Primary",3),rep("Relapse",3)), Pairs=rep(c("PR1","PR2","PR3"),2), stringsAsFactors=F)
rm(a,b,nbcnt)

# RNAseq read-counts normalization
infNorm = DGEList(counts=infdat, group=c(1,1,1,0,0,0), samples = infcli, remove.zeros = T)
isexpr  = rowSums(cpm(infNorm) > 1) > 1
#table(isexpr)
infNorm = infNorm[isexpr, ,keep.lib.sizes=FALSE]
infNorm = calcNormFactors(infNorm, method="TMM")
rm(infdat,infcli,isexpr)

#pdf("~/Documents/tmp/primary_relapseMDS.pdf", width=3,height=3)
# par(mar=c(3,3,0.5,0.5), mgp=c(2,1,0), cex=0.75)
# plotMDS(infNorm, top=1000, xlim=c(-4,4), ylim=c(-4,4), col=as.numeric(factor(infNorm$samples$Pairs)))
#dev.off()

design = model.matrix(~ 0+Pairs+Type, data=infNorm$samples)
colnames(design) = c("PR1","PR2","PR3","Relapse")
infNorm <- estimateDisp(infNorm, design, robust=TRUE)

#infNorm$common.dispersion
#plotBCV(infNorm)

fit = glmFit(infNorm, design)
lrt = glmLRT(fit)
deGen = topTags(lrt, n=Inf)$table

topHits = deGen[deGen$FDR < 0.01 & abs(deGen$logFC) > 1,]
write.table(data.frame(Genes= rownames(topHits), topHits, stringsAsFactors = F),
            paste0(outpath1, "tumor/Rel_vs_Pri/RelapseVsPrimary_topDiffExpGenes.txt"),
            row.names = F, quote = F, sep="\t")

## Barcode-plot with selected pathways

#-- RAS metagene signature
rasUp = readRDS(rasSigr)$upreg

#-- JUN FOS target genes
net = read.table(NBreg, header=T, stringsAsFactors=F, sep="\t")[,1:2]
net =  net[net$Regulator %in% c("FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND"),]
net = net[order(net$Regulator),]
junfos = unique(net$Target)
rm(net)

# -- assign index
tmplist = list(RasUp=rasUp, JunFos=junfos)
idx = ids2indices(tmplist, id = sapply(strsplit(rownames(infNorm),"|",fixed=T),function(x) x[2]))
cam = camera(infNorm, idx, design)

pdf(paste0(outpath2,"figure4/JunFos_Ras_enrichment_barcodeplot_relapsevsPrimary.pdf"), width=7, height=1)
par(mfrow=c(1,2), mar=c(3,1,1,0.5), mgp=c(1.5,0.5,0), cex=0.7)
 barcodeplot(lrt$table$logFC,
            index  = idx[["RasUp"]],
            quantiles = quantile(lrt$table$logFC)[c(2,4)],
            labels = c(paste0("High in", "Primary", sep="\n"), paste0("High in", "Relapse", sep="\n")),
            main = paste0("Ras activated target genes | fdr=", signif(cam$FDR, 1)[1]),
            alpha=0.1, worm=F, col.bars="firebrick")

 barcodeplot(lrt$table$logFC,
            index  = idx[["JunFos"]],
            quantiles = quantile(lrt$table$logFC)[c(2,4)],
            labels = c(paste0("High in", "Primary", sep="\n"), paste0("High in", "Relapse", sep="\n")),
            main = paste0("Jun Fos target genes | fdr=", signif(cam$FDR, 1)[2]),
            alpha=0.1, worm=F, col.bars="firebrick")
dev.off()
rm(rasUp, tmplist, idx, cam)

# pdf("~/Documents/tmp/primary_relapseVolcano.pdf", width=3, height=3)
# par(mar=c(3,3,0.5,0.5), mgp=c(2,1,0), cex=0.75, xaxs="i", yaxs="i")
# plot(y = -log10(deGen$FDR), x = deGen$logFC, pch=20,cex=0.5,
#      xlab = "logFC", ylab = "-log10(fdr)", frame=F)
# abline(v = c(-1,1))
# abline(h = -log10(0.05))
# dev.off()

# Create TF regulome genesets
net = read.table(NBreg, header=T, stringsAsFactors=F, sep="\t")[,1:2]
net = split(net$Target, net$Regulator)

# CRCs and Mes TF activity
crc = read.table(crcList, header=T, sep="\t", stringsAsFactors = F)
mesTFact = readRDS(mesTFact)
mesTFact = data.frame(TFs = names(mesTFact$es$nes), TFact = mesTFact$es$nes, stringsAsFactors = F)
mesTFact$isCRC = "No"
mesTFact$isCRC[mesTFact$TFs %in% crc$TFs] = "Yes"

# TF activity kmeans
set.seed(123); nclust = 4 
km = data.frame(TFact = mesTFact$TFact, k = kmeans(mesTFact$TFact, nclust)$cluster)
km = split(km$TFact, km$k)
km = sort(sapply(km, max))
km = as.numeric(head(km, c(nclust-1)))
rm(nclust)

# CAMERA analysis
idx = ids2indices(net, id = sapply(strsplit(rownames(infNorm),"|",fixed=T),function(x) x[2]))
cam = camera(infNorm, idx, design)
cam = data.frame(TFs = rownames(cam), cam, stringsAsFactors = F)
cam = cam[order(cam$Direction, cam$FDR),]
cam$logFDR = -log10(cam$FDR)
cam$Direction[cam$FDR > 0.01] = "Not significant"
cam$Direction[cam$Direction == "Up"] = "High in Relapse"
cam$Direction[cam$Direction == "Down"] = "High in Primary"
rm(idx)

# Write the enrichment analysis results
write.table(cam, paste0(outpath1, "tumor/Rel_vs_Pri/RelapseVsPrimary_EnrichTFregulons.txt"),
            row.names = F, quote = F, sep="\t")

# Combine Camera results, Mes activity and CRC TF status
resAllTF = merge(cam, mesTFact)

# Plot all TF's
p.allTF = ggplot(resAllTF, aes(x=TFact, y=logFDR)) + theme_bw(base_size=10) + 
  labs(x="MES TF activity", y = paste("-log10 FDR", "Enrichment of target genes of TFs", "(Relapse vs Primary)", sep="\n")) +
  geom_point(aes(colour=Direction), size=1, alpha=0.5) + 
  scale_colour_manual(name = paste("Expression of", "TF target genes", sep="\n"), values = c("#1a9850", "#542788", "grey80")) +
  geom_text_repel(data=subset(resAllTF, FDR < 0.01 & isCRC == "Yes"), aes(label=TFs), size=2) + 
  
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = km, linetype = "dashed", colour = "grey50") +
  
  theme(legend.position="right",
        legend.direction="vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf(paste0(outpath2,"sup_figure4/all_TF_in_relapse_pri_enriched.pdf"), width=5, height=3.5)
print(p.allTF)
dev.off()

# Plot only CRCs
resOnlyCRC = resAllTF[resAllTF$isCRC == "Yes",]

p.crc = ggplot(resOnlyCRC, aes(x=TFact, y=logFDR)) + theme_bw(base_size=10) + 
  labs(x="MES TF activity", y = paste("-log10 FDR", "Enrichment of target genes of CRC TFs", "(Relapse vs Primary)", sep="\n")) +
  geom_point(aes(colour=Direction), size=2, alpha=0.5) + 
  scale_colour_manual(name = paste("Expression of CRC", "TF target genes", sep="\n"), values = c("#1a9850", "#542788", "grey80")) +
  geom_text_repel(data=subset(resOnlyCRC, FDR < 0.01 & isCRC == "Yes"), aes(label=TFs), size=2) + 
  
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = km, linetype = "dashed", colour = "grey50") +
  
  theme(legend.position="right",
        legend.direction="vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf(paste0(outpath2,"figure4/crcTF_in_relapse_pri_enriched.pdf"), width=5, height=3.5)
print(p.crc)
dev.off()
