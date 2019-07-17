options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

NBexprs = as.character(args[1])
tumoNMF = as.character(args[2])
rasSigr = as.character(args[3])
NBreg   = as.character(args[4])
NBmut   = as.character(args[5])
outpath = as.character(args[6])

#-------------------------------------------------------------------------------------------------
# DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"
# 
# NBexprs = paste0(DATAPATH, 'analysis/tumor/rnaseq/exprs/tumor_RNAseq_TPM_Matrix_filt_log.RDS')
# tumoNMF = paste0(DATAPATH, 'analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS')
# rasSigr = paste0(DATAPATH, 'db/publicGeneSigs/ras_target_genes.RDS')
# NBreg   = paste0(DATAPATH, 'analysis/tumor/ARACNe/network.txt')
# NBmut   = paste0(DATAPATH, 'annotation/NB_mutation_matrix.RDS')
# outpath = paste0(DATAPATH, "results/")
#-------------------------------------------------------------------------------------------------

library(reshape2)
library(readxl)

dat = readRDS(NBexprs)
sym = sapply(strsplit(rownames(dat),"|",fixed=T),function(x)x[2])

expo = readRDS(tumoNMF)
identical(colnames(expo), colnames(dat))

## RAS metagene signature
ras = readRDS(rasSigr)
ras = ras$upreg

## JUN FOS (AP1 complex) target genes
net = read.table(NBreg, header=T, stringsAsFactors=F, sep="\t")[,1:2]
net =  net[net$Regulator %in% c("FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND"),]
net = net[order(net$Regulator),]
junfos = unique(net$Target)
rm(net)

## Median expression of JUN/FOS target genes
vals1 = apply(dat,2,function(x){
x = median(x[sym %in% junfos])
})

## Median expression of RAS target genes
vals2 = apply(dat,2,function(x){
x = median(x[sym %in% ras])
})

## RAS mutation status
NBmuts = readRDS(NBmut)
NBmuts = NBmuts[rownames(NBmuts) %in% c("NRAS", "KRAS", "HRAS", "BRAF", "RAF1", "NF1"),]
NBmuts = NBmuts[, -which(apply(NBmuts, 2, function(x)sum(is.na(x))) == nrow(NBmuts))]


length(vals1)
length(vals2)
dim(expo)
#--------------
# JUN FOS plot
#--------------

pdf(paste0(outpath,"figure4/junfos_corr_exposures.pdf"), width=4.3, height=3.8)
layout(matrix(c(1:4), ncol=2,byrow=T), width=1,height=1)
par(oma=c(1,0,0,0),mar=c(2,3.5,0.25,0.5), mgp=c(2.5,0.5,0), cex=0.7) # xaxs="i", yaxs="i"

#Sig1
plot(x=vals1, y=expo[1,], xlab="",ylab=paste(rownames(expo)[1], "Exposure"),las=2, frame=F, 
     cex=1, pch=20,col=rgb(0, 0, 1, 0.2))
abline(lm(expo[1,]~vals1))
corr = cor(vals1, expo[1,], method="spearman")
text(x=3, y=5000, label=paste0("rho=",round(corr,2)))

#Sig2
plot(x=vals1, y=expo[2,], xlab="", ylab=paste(rownames(expo)[2], "Exposure"),las=2, frame=F, 
      cex=1, pch=20,col=rgb(0, 0, 1, 0.2))
abline(lm(expo[2,]~vals1))
corr = cor(vals1, expo[2,], method="spearman")
text(x=3, y=3500,label=paste0("rho=",round(corr,2)))

#Sig3
plot(x=vals1, y=expo[3,], xlab="",ylab=paste(rownames(expo)[3], "Exposure"),las=2, frame=F,
     cex=1, pch=20,col=rgb(0, 0, 1, 0.2))
abline(lm(expo[3,]~vals1))
corr = cor(vals1, expo[3,], method="spearman")
text(x=3, y=4000,label=paste0("rho=",round(corr,2)))

#Sig4
plot(x=vals1, y=expo[4,], xlab="", ylab=paste(rownames(expo)[4], "Exposure"), las=2, frame=F,
     cex=1, pch=20,col=rgb(0, 0, 1, 0.2)) 
abline(lm(expo[4,]~vals1))
corr = cor(vals1, expo[4,], method="spearman")
text(x=3, y=5000,label=paste0("rho=",round(corr,2)))

mtext(text="JUN/FOS TF family target genes median expression per sample (log2 TPM)", side=1,outer=T, line=0, cex=0.72)
dev.off()

#----------
# RAS plot
#----------

pdf(paste0(outpath,"figure4/ras_corr_exposures.pdf"), width=4.3, height=3.8)
layout(matrix(c(1:4), ncol=2,byrow=T), width=1,height=1)
par(oma=c(1,0,0,0),mar=c(2,3.5,0.25,0.5), mgp=c(2.5,0.5,0), cex=0.7) # xaxs="i", yaxs="i",

#Sig1
plot(x=vals2, y=expo[1,], xlab="",ylab=paste(rownames(expo)[1], "Exposure"),las=2, frame=F, 
   cex=1, pch=20,col=rgb(0, 0, 1, 0.2))
abline(lm(expo[1,]~vals2))
corr = cor(vals2, expo[1,], method="spearman")
text(x=3, y=5000,label=paste0("rho=",round(corr,2)))

xmark = vals2[which(names(vals2) %in% colnames(NBmuts))]
ymark = expo[1,][which(colnames(expo) %in% colnames(NBmuts))]
points(x = xmark, y = ymark, pch = 4, col = "firebrick")
rm(xmark, ymark)

#Sig2
plot(x=vals2, y=expo[2,], xlab="", ylab=paste(rownames(expo)[2], "Exposure"),las=2, frame=F, 
     cex=1, pch=20,col=rgb(0, 0, 1, 0.2))
abline(lm(expo[2,]~vals2))
corr = cor(vals2, expo[2,], method="spearman")
text(x=3, y=3500,label=paste0("rho=",round(corr,2)))

xmark = vals2[which(names(vals2) %in% colnames(NBmuts))]
ymark = expo[2,][which(colnames(expo) %in% colnames(NBmuts))]
points(x = xmark, y = ymark, pch = 4, col = "firebrick")
rm(xmark, ymark)

#Sig3
plot(x=vals2, y=expo[3,], xlab="",ylab=paste(rownames(expo)[3], "Exposure"),las=2, frame=F,
     cex=1, pch=20,col=rgb(0, 0, 1, 0.2))
abline(lm(expo[3,]~vals2))
corr = cor(vals2, expo[3,], method="spearman")
text(x=3, y=4000,label=paste0("rho=",round(corr,2)))

xmark = vals2[which(names(vals2) %in% colnames(NBmuts))]
ymark = expo[3,][which(colnames(expo) %in% colnames(NBmuts))]
points(x = xmark, y = ymark, pch = 4, col = "firebrick")
rm(xmark, ymark)

#Sig4
plot(x=vals2, y=expo[4,], xlab="", ylab=paste(rownames(expo)[4], "Exposure"),las=2, frame=F,
     cex=1, pch=20,col=rgb(0, 0, 1, 0.2))
abline(lm(expo[4,]~vals2))
corr = cor(vals2, expo[4,], method="spearman")
text(x=3, y=5000,label=paste0("rho=",round(corr,2)))

xmark = vals2[which(names(vals2) %in% colnames(NBmuts))]
ymark = expo[4,][which(colnames(expo) %in% colnames(NBmuts))]
points(x = xmark, y = ymark, pch = 4, col = "firebrick")
rm(xmark, ymark)

mtext(text="RAS target genes median expression per sample (log2 TPM)", side=1,outer=T, line=0, cex=0.72)
dev.off()

#----------------------------------
# Correlation to only RAS proteins
#----------------------------------

rasp = dat[sym %in% c("HRAS","KRAS","NRAS"),]
pdf(paste0(outpath,"sup_figure4/ras_proteins_corr_exposures.pdf"), width=7, height=7)
layout(matrix(c(1:12), ncol=3, nrow=4, byrow=T), width=1,height=1)
par(mar=c(3.5,3.5,1,0.25), mgp=c(2.5,0.5,0),cex=0.7) # xaxs="i", yaxs="i", 

for(i in 1:nrow(expo))
{
  for(j in 1: nrow(rasp))
  {
    x = as.numeric(rasp[j,])
    y = as.numeric(expo[i,])
    
    if(j==1){yl=rownames(expo)[i]; xl=""}else{yl="";xl=""}
    if(i==4){xl=rownames(rasp)[j]; yl=rownames(expo)[i]}

    corr = round(cor(x,y,method="spearman"),2)
    plot(x=x, y=y, xlab=xl, ylab=yl, cex=0.5, pch=20, las=2, frame=F, col=rgb(0, 0, 1, 0.2), main=paste0("rho=",corr), cex.main=0.8)
    abline(lm(y~x))
    
    xmark = x[which(colnames(rasp) %in% colnames(NBmuts))]
    ymark = y[which(colnames(expo) %in% colnames(NBmuts))]
    points(x = xmark, y = ymark, pch = 4, col = "firebrick")
    rm(xmark, ymark)
  }
}
dev.off()

#-------------------------------------------
# Correlation to only AP1 complex  proteins
#-------------------------------------------

ap1 = dat[sym %in% c("FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND"),]
rownames(ap1) = sapply(strsplit(rownames(ap1),"|",fixed=T), function(x)x[2])
ap1 = ap1[match(c("FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND"),rownames(ap1)),]

pdf(paste0(outpath,"sup_figure4/AP1_complex_proteins_corr_exposures.pdf"), width=9.5, height=6.5)
layout(matrix(c(1:28), ncol=7, nrow=4, byrow=T), width=1,height=1)
par(mar=c(3.5,3.5,1,0.25), mgp=c(2.5,0.5,0), cex=0.7) #xaxs="i", yaxs="i", 

for(i in 1:nrow(expo))
{
  for(j in 1: nrow(ap1))
  {
    x = as.numeric(ap1[j,])
    y = as.numeric(expo[i,])
    
    if(j==1){yl=rownames(expo)[i];xl=""}else{yl="";xl=""}
    if(i==4){xl=rownames(ap1)[j]; yl=rownames(expo)[i]}
    
    corr = round(cor(x,y,method="spearman"),2)
    plot(x=x, y=y, xlab=xl, ylab=yl, cex=0.5, pch=20, las=2, frame=F, col=rgb(0, 0, 1, 0.2), main=paste0("rho=",corr), cex.main=0.8)
    abline(lm(y~x))
  }
}
dev.off()


