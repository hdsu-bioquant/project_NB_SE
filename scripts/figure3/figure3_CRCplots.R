options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

KDdata   = as.character(args[1])
TFact    = as.character(args[2])
tumorNMF = as.character(args[3])
cellNMF  = as.character(args[4])
tumorCRC = as.character(args[5])
cellCRC  = as.character(args[6])
outpath  = as.character(args[7])

#-------------------------------------------------------------------------------------------------
 # DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"
 # 
 # KDdata   = paste0(DATAPATH, 'analysis/cells/crcGIEMSAkd/nbKDinhouse.RDS')
 # TFact    = paste0(DATAPATH, 'analysis/tumor/VIPER/')
 # tumorNMF = paste0(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS')
 # cellNMF  = paste0(DATAPATH, 'analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K3_Hmatrix_hnorm.RDS')
 # tumorCRC = paste0(DATAPATH, 'data/tumor/chipseq/H3K27ac/CRC')
 # cellCRC  = paste0(DATAPATH, 'data/cells/chipseq/H3K27ac/CRC')
 # outpath  = paste0(DATAPATH, "results/")
#-------------------------------------------------------------------------------------------------

## Load required packages
library(RColorBrewer)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(ggrepel)

## Load required files
# (1) Our knockdown analysis data (internal validation) 
# (2) TF activity file
# (3) NMF Signature specific annotation
# (4) List all CRC files for loading

## (1) Our knockdown analysis data (internal validation)
kd.NB.fc = readRDS(KDdata)
kd.NB.fc = kd.NB.fc[order(kd.NB.fc$log2FC),,drop=F]

## (2) TF activity file
a = readRDS(paste0(TFact, "MYCN_TFactivity.RDS"))
b = readRDS(paste0(TFact, "MES_TFactivity.RDS"))
c = readRDS(paste0(TFact, "MNA-LR_TFactivity.RDS"))
d = readRDS(paste0(TFact, "MNA-HR_TFactivity.RDS"))

viperTFactivity = data.frame(MYCN = a$es$nes, 
                             MES = b$es$nes,
                            `MNA-LR` = c$es$nes,
                            `MNA-HR` = d$es$nes)

viperTFactivity = apply(viperTFactivity, 2, function(x){c(x-median(x))/mad(x)})

write.table(data.frame(TFs = rownames(viperTFactivity), round(viperTFactivity,4), stringsAsFactors = F),
            paste0(outpath, "supptables/TFactivity_across_all_signatures_ZnormPerSig.txt"), row.names=F, quote=F, sep="\t")

pheatmap(viperTFactivity, clustering_method = "ward.D2", scale = "none", 
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100), 
         cutree_rows=1, show_rownames=F, treeheight_row = 15, treeheight_col = 15,
         filename = paste0(outpath,"sup_figure3/TFactivity_heatmap_across_all_signatures_ZnormPerSig.pdf"), width = 1.4, height = 3.7)
#dev.off()
rm(a,b,c,d)

## (3)  NMF Signature specific annotation

# Tumors
primsig = readRDS(tumorNMF)
primsig = apply(primsig, 2, which.max)
primsig = split(names(primsig), primsig)
names(primsig) = c("MYCN-Amp-T","Mesenchymal-T","LR-MYCN-nonAmp-T","HR-MYCN-nonAmp-T")

#Cells
cellsig = readRDS(cellNMF)
cellsig = apply(cellsig, 2, which.max)
cellsig = split(names(cellsig), cellsig)
names(cellsig) = c("MYCN-Amp-C", "Mesenchymal-C", "MYCN-nonAmp-C")

## (4) List all CRC files for loading

# Tumor CRC file list
files1  = list.files(tumorCRC, full.names = T)
sampId1 = sapply(strsplit(sapply(strsplit(files1,"/"),function(x){x[length(x)]}), "_"), function(y)y[1])

# Cell-line CRCs file list
files2  = list.files(cellCRC, full.names = T)
sampId2 = sapply(strsplit(sapply(strsplit(files2,"/"),function(x){x[length(x)]}), "_"), function(y)y[1])

# Tumor + Cell-line CRcs file list
files = c(files1, files2)
sampId = c(sampId1,sampId2)
anno = c(rep("Primary Tumors",length(sampId1)),rep("Cell-lines",length(sampId2)))
rm(files1,files2,sampId1,sampId2)

#table(anno)
#Cell-lines Primary Tumors 
#         25             60 

# Reading all the CRC files
crcTF = vector("list",length(sampId))
names(crcTF) = sampId

for( i in 1:length(files))
{
    tmp=readLines(files[i])
    if(length(tmp)>0)
    {
       tmp = sapply(strsplit(tmp, "\t"), function(x)x[1])
       pattern = c("[", "]", "'", " ")
       for(j in pattern){tmp = gsub(j, "", tmp, fixed=T)}
       tmp = strsplit(tmp,",")
       crcTF[[i]] = tmp
       rm(j, pattern, tmp)
    }
}
rm(i)

# CRC TF occurance matrix
cmn = unique(unlist(crcTF))
mat = matrix(0,ncol=length(sampId), nrow=length(cmn))
colnames(mat) = sampId
rownames(mat) = cmn

for( i in 1:ncol(mat))
{
  tmp = unique(unlist(crcTF[[colnames(mat)[i]]]))
  if(!is.null(tmp))
  {
    mat[which(rownames(mat) %in% tmp),i] = 1
  }
 rm(tmp)
}
rm(i,cmn)

# Sample filtering criteria - remove all samples where no CRCs were called
mat = mat[, -which(colSums(mat) == 0)]
rm(sampId, files)

# TF activity corresponding to each CRC TF
tfact = matrix(NA, ncol=1, nrow=nrow(mat))
rownames(tfact) = rownames(mat)
tfact = data.frame(Genes= rownames(tfact), ActivityNA=tfact[,1], stringsAsFactors=F)

viperTFactivity = data.frame(Genes= rownames(viperTFactivity), Activity=viperTFactivity, stringsAsFactors=F)

tfact = merge(x=tfact, y=viperTFactivity, all.x=T)
tfact = tfact[match(rownames(mat), tfact$Genes),]
rownames(tfact)= tfact$Genes
tfact = tfact[,- c(1:2), drop=F]

# TF filtering criteria - remove all CRC TF observed in less than 5% of all samples
selCRC = mat[which(rowSums(mat) >= round(0.05*ncol(mat))),]
selTF  = tfact[which(rowSums(mat) >= round(0.05*ncol(mat))),,drop=F]

selCRC = selCRC[order(selTF$Activity.MES, decreasing=T),]
selTF  = selTF[order(selTF$Activity.MES, decreasing=T),,drop=F]
rm(tfact,mat)

## Remove NA valued TF activity
rmv = which(is.na(selTF$Activity.MES))
selCRC = selCRC[-rmv,]
selTF = selTF[-rmv,,drop=F]
rm(rmv)

## Fraction of CRC TF observed per signature in tumors and cells
primsig = melt(lapply(primsig,function(x){x[x %in% colnames(selCRC)]}))
cellsig = melt(lapply(cellsig,function(x){x[x %in% colnames(selCRC)]}))
anno = data.frame(rbind(primsig,cellsig), Type=c(rep("tumor",nrow(primsig)),rep("cell",nrow(cellsig))))
colnames(anno) = c("Samples","Signature","Type")
anno$Samples = as.character(anno$Samples)
anno$Signature = factor(anno$Signature,levels=c("MYCN-Amp-T","Mesenchymal-T","LR-MYCN-nonAmp-T","HR-MYCN-nonAmp-T","MYCN-Amp-C","Mesenchymal-C","MYCN-nonAmp-C"))
anno$Type = factor(anno$Type, levels=c("tumor","cell"))
anno = anno[match(colnames(selCRC), anno$Samples),]
rm(primsig,cellsig)

#dim(selCRC)
# 65 76 
# 65 CRC TF and 76 samples (51 tumors and 25 cells)

#summary(anno)

#Samples                     Signature     Type   
#Length:76          MYCN-Amp-T      :13   tumor:51  
#Class :character   Mesenchymal-T   : 4   cell :25  
#Mode  :character   LR-MYCN-nonAmp-T:21             
#                   HR-MYCN-nonAmp-T:13             
#                   MYCN-Amp-C      : 8             
#                   Mesenchymal-C   : 6             
#                   MYCN-nonAmp-C   :11         

# Counting  number of CRC TF in tumors and  cells separately
#tmp1 = split(data.frame(t(selCRC)),anno$Type)
#sapply(tmp1, function(x)sum(colSums(x) > 0))
#tumor  cell 
#61    62 

#tmp2 = split(data.frame(t(selCRC)),anno$Signature)
#sapply(tmp2, function(x)sum(colSums(x) > 0))
#MYCN-Amp-T    Mesenchymal-T L    R-MYCN-nonAmp-T    HR-MYCN-nonAmp-T   MYCN-Amp-C    Mesenchymal-C    MYCN-nonAmp-C 
#        25               25                  47                   35           27               39               42
#rm(tmp1,tmp2)

## Computing Fraction observerd and Proportion observed in CRC TF
fracObs = split(data.frame(t(selCRC)),anno$Signature)
fracObs = round(sapply(fracObs,function(x){colSums(x)/nrow(x)})*100,2)

propObs = cbind(prop.table(fracObs[,1:4],1), prop.table(fracObs[,5:ncol(fracObs)],1))*100
propObs[is.nan(propObs)]=0

write.table(data.frame(TFs = rownames(fracObs), fracObs, stringsAsFactors = F), 
            paste0(outpath,"supptables/crcTF_fractionObserved.txt"), quote=F, sep="\t")

## Heatmap of fraction observed CRC across signatures
modules = pheatmap(fracObs, scale = "none", clustering_method = "ward.D2",  clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", silent = T)
crcClust = data.frame('Major class' = cutree(modules$tree_row, k=2), 
                      'Minor class' = as.character(cutree(modules$tree_row, k=5)), stringsAsFactors = F, check.names = F)

mes.marker =  crcClust$`Major class`[rownames(crcClust) == "MYC"]
crcClust$`Major class` = ifelse(crcClust$`Major class` == mes.marker, "MES", "ADRN")

labs = rep("Cell",ncol(fracObs))
labs[grep("-T",colnames(fracObs))] = "Tumor"
names(labs) = colnames(fracObs)
labs =  data.frame(Type = labs)

colabs = gsub("-T","", colnames(fracObs))
colabs = gsub("-C","", colabs)

colabs[colabs %in% "Mesenchymal"] = "MES"
colabs[colabs %in% "MYCN-Amp"] = "MYCN"
colabs[colabs %in% "MYCN-nonAmp"] = "MNA"
colabs[colabs %in% "LR-MYCN-nonAmp"] = "MNA-LR"
colabs[colabs %in% "HR-MYCN-nonAmp"] = "MNA-HR"

color = list(Type = c('Tumor' = "grey20", 'Cell' = "grey80"),
             'Major class' = c('MES' = "#420A68FF", 'ADRN' = "#2FB47CFF"),
             'Minor class' = c('1' = "#fb8072",'2' = "#80b1d3", '3' = "#fdb462",'4' = "#8dd3c7", '5' = "#fccde5")
             )

pheatmap(fracObs, color = c("white",colorRampPalette(brewer.pal(9,"Purples")[2:9])(49)),  breaks = c(0,seq(2,100,2)), border_color = "white",
scale = "none", clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", 
annotation_col = labs, annotation_row = crcClust, annotation_colors = color , 
clustering_method = "ward.D2",  treeheight_row = 30, treeheight_col = 10, fontsize = 6, labels_col = colabs, #cutree_rows = 5, cutree_cols=4,
filename = paste0(outpath,"figure3/crcTF_fractionObserved_tumors_cells.pdf"), 
#width = 4.3, height = 8.5) #if uncommenting selTFtrans above
width = 3, height = 7)
#dev.off()

rm(mes.marker, labs, colabs, color)

# Distribution of TF activities for the Major and Minor groups
# identified in the above plot across signatures
crcClust = crcClust[modules$tree_row$order,]
crcClust = data.frame(Genes = rownames(crcClust), crcClust, stringsAsFactors = F)
write.table(crcClust, paste0(outpath,"supptables/crcTF_modules.txt"), row.names = F, quote=F, sep="\t")

crcClust = merge(viperTFactivity, crcClust)
colnames(crcClust) = gsub("Activity.", "", colnames(crcClust), fixed=T)

#Plot
crcClust = melt(crcClust)
crcClust$variable = factor(as.character(crcClust$variable), levels = c("MES", "MYCN", "MNA.HR", "MNA.LR"))
crcClust$Major.class = factor(crcClust$Major.class, levels=c("MES","ADRN"))
crcClust$Minor.class = factor(crcClust$Minor.class, levels=c("1","5","3","2","4"))

# All groups
my_comparisons <- list( c("MES", "MYCN"), c("MES", "MNA.HR"), c("MES", "MNA.LR") )
p <- ggboxplot(crcClust, x = "variable", y = "value", xlab = "", ylab = "TF activity (Z normalized)",
               color = "variable", palette = "npg",
               add = "jitter", outlier.shape = NA) + 
  stat_compare_means(label = "p.format", method="wilcox", label.x=2, size=3, comparisons = my_comparisons) +
  #stat_compare_means(label = "p.format", method="anova", label.x=2, size=3) +
  geom_text_repel(aes(label=Genes), size=2) + 
  theme_pubr(base_size = 7)
p = facet(p, facet.by = "Minor.class", scales="free")

pdf(paste0(outpath,"sup_figure3/crcTF_TFactivity_per_Minor_class_and_Signature.pdf"), width = 8, height=8)
print(p)
dev.off()
rm(p)

# Only major groups
p <- ggboxplot(crcClust, x = "variable", y = "value", xlab = "", ylab = "TF activity (Z normalized)",
               color = "variable", palette = "npg",
               add = "jitter", outlier.shape = NA) + 
  stat_compare_means(label = "p.format", method="wilcox", label.x=2, size=3, comparisons = my_comparisons) +
  #geom_text_repel(aes(label=CRC), size=2) + 
  theme_pubr(base_size = 7)
p = facet(p, facet.by = "Major.class", scales="free")

pdf(paste0(outpath,"figure3/crcTF_TFactivity_per_Major_class_and_Signature.pdf"), width = 5, height = 5.5)
print(p)
dev.off()
rm(p)

rm(modules, crcClust, my_comparison)

## Selecting the Mes Activity
selTF.Mes = selTF[,"Activity.MES", drop=F]
    
## KD vs TF activity correlation (internal KD)
pdf(paste0(outpath,"figure3/crcTF_correlation_TFactivity_vs_KD.pdf"), width=3.5, height=3)
  cmn = intersect(rownames(kd.NB.fc), rownames(selTF.Mes))
  a = kd.NB.fc[rownames(kd.NB.fc) %in% cmn ,,drop=F]
  b = selTF.Mes[rownames(selTF.Mes) %in% cmn ,,drop=F]
  a = a[match(rownames(b),rownames(a)),,drop=F]
  
  corr = cor.test(a$log2FC, b$Activity.MES, use="pairwise.complete.obs", method="spearman")
  df = data.frame(a,b)
    
  ggplot(df, aes(x=log2FC, y=Activity.MES)) + theme_bw(base_size=9) + 
  labs(x="Knockdown sensitivity log2(GIMEN/NMB)", y="TF activity") + geom_point(size=1.5) +
  # geom_point(aes(color=log2FC),size=1.5) + scale_color_gradient(high="yellow", low="red") +
  geom_smooth(method = "lm", se = FALSE, colour="grey30") +
  geom_text_repel(aes(label=rownames(df)),size=2.5) + 
  annotate("text", x = -0.43, y = -3.3, label = paste(paste0("p=",round(corr$p.value,3)),paste0("rho=",round(corr$estimate,2)),sep="\n"), size=2.8)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  rm(a,b,cmn,corr,df)
dev.off()

# Multi-dimentional and Integrative CRC plot
pdf(paste0(outpath,"figure3/crcTF_oncoprints_SignatureSpecific.pdf"),width=7.4,height=5.8)

layout(rbind(
             c(1,0,0,0,2),
             c(3:7),
             c(0,8,8,9,9)
            ),
         width=c(2.04,0.35,0.1,0.35,1.1),
         height=c(0.015,0.7,0.08)
        )
  

# Layout while including the KD analysis and its legend
# layout(rbind(
#               c(1,0,0,0,2),
#               c(3:7),
#               rep(8,5), 
#               c(9,10,10,11,11)
#             ),
#        width=c(2.04,0.35,0.1,0.35,1.1),
#        height=c(0.015,0.7,0.08,0.08)
#       )

  # Tumors signatures
  tmpAnnoTumor = anno[anno$Type == "tumor",]
  tmpAnnoTumor = tmpAnnoTumor[order(tmpAnnoTumor$Signature),]
  tmpAnnoTumor = droplevels(tmpAnnoTumor)
  tmpCRCTumor  = selCRC[,colnames(selCRC) %in% tmpAnnoTumor$Samples]
  tmpCRCTumor  = tmpCRCTumor[,match(tmpAnnoTumor$Samples,colnames(tmpCRCTumor))]
  colnames(tmpCRCTumor) = gsub("B087_","",colnames(tmpCRCTumor))
  colnames(tmpCRCTumor) = gsub("_","",colnames(tmpCRCTumor))
  colnames(tmpCRCTumor) = gsub("-","",colnames(tmpCRCTumor))

  # Tumor sample annotation
  idT = rep(NA,length(tmpAnnoTumor$Signature))
  idT[tmpAnnoTumor$Signature == "MYCN-Amp-T"]=1
  idT[tmpAnnoTumor$Signature == "Mesenchymal-T"]=2
  idT[tmpAnnoTumor$Signature == "LR-MYCN-nonAmp-T"]=3
  idT[tmpAnnoTumor$Signature == "HR-MYCN-nonAmp-T"]=4
  idT = as.numeric(idT)
  idT = matrix(idT,ncol=1)
  sigBordersTumor = table(tmpAnnoTumor$Signature)

  # Cell line signatures
  tmpAnnoCell = anno[anno$Type == "cell",]
  tmpAnnoCell = tmpAnnoCell[order(tmpAnnoCell$Signature),]
  tmpAnnoCell = droplevels(tmpAnnoCell)
  tmpCRCcell  = selCRC[,colnames(selCRC) %in% tmpAnnoCell$Samples]
  tmpCRCcell  = tmpCRCcell[,match(tmpAnnoCell$Samples,colnames(tmpCRCcell))]

  # Cell sample annotation
  idC = rep(NA,length(tmpAnnoCell$Signature))
  idC[tmpAnnoCell$Signature == "MYCN-Amp-C"]=1
  idC[tmpAnnoCell$Signature == "Mesenchymal-C"]=2
  idC[tmpAnnoCell$Signature == "MYCN-nonAmp-C"]=3
  idC = as.numeric(idC)
  idC = matrix(idC,ncol=1)
  sigBordersCell = table(tmpAnnoCell$Signature)

  #Sample coloring Tumors
  par(mar=c(0,2.5,0.25,0.25),xaxs="i",yaxs="i")
  image(idT,col=c("#006837","#40004B","#31a354","#78c679"),xaxt="n", yaxt="n")
  grid(col="white",nx=nrow(idT), ny=0, lty=1,lwd=1); box()
  abline(v=seq(0,1,length.out=ncol(tmpCRCTumor))[c(sigBordersTumor[1], sum(sigBordersTumor[1:2]), sum(sigBordersTumor[1:3]))]+0.01, lwd=1)

  #Sample coloring cells
  par(mar=c(0,0.25,0.25,2.5),xaxs="i",yaxs="i")
  image(idC,col=c("#006837","#40004B","#31a354"),xaxt="n", yaxt="n")
  grid(col="white",nx=nrow(idC), ny=0, lty=1,lwd=1); box()
  abline(v=seq(0,1,length.out=ncol(tmpCRCTumor))[c(sigBordersCell[1], sum(sigBordersCell[1:2]), sum(sigBordersCell[1:3]))]+0.01, lwd=1)
  
  # Tumor CRC oncoprints
  par(mar=c(4,2.5,0.25,0.25),mgp=c(0,0.2,0),xaxs="i",yaxs="i")
  image(t(tmpCRCTumor),col=c("#f0f0f0","black"),xaxt="n",yaxt="n")
  axis(side=2, at=seq(0,1,length.out=nrow(tmpCRCTumor)), labels=rownames(tmpCRCTumor), las=2, cex.axis=0.6, tick=F)
  grid(nx=ncol(tmpCRCTumor),ny=nrow(tmpCRCTumor),col="white", lty=1,lwd=1); box()
  axis(side=1, at=seq(0,1,length.out=ncol(tmpCRCTumor)), labels=colnames(tmpCRCTumor), cex.axis=0.6, tick=F,las=2)
  abline(v=seq(0,1,length.out=ncol(tmpCRCTumor))[c(sigBordersTumor[1], sum(sigBordersTumor[1:2]), sum(sigBordersTumor[1:3]))]+0.01, lwd=1)
  #axis(side=3, at=0.5, labels="CRC frequency (Tumors)", cex.axis=1, tick=F)

  # Tumor CRC proportion barplots
  par(mar=c(4,0.1,0.25,0.25),mgp=c(0,0.2,0),xaxs="i",yaxs="i")
  barplot(t(propObs[,1:4]),horiz=T,las=2,cex.axis=0.9,cex.lab=0.9,cex.names=0.8,col=c("#006837","#40004B","#31a354","#78c679"),border=NA,yaxt="n",xaxt="n")
  axis(side=1, at=c(0,50,100), labels=c(0,50,100), cex.axis=0.9, tick=F,las=2)
  axis(side=3, at=50, labels="Proportion(%)", cex.axis=0.8, line= -0.3, tick=F)
  box(col="black")

  # TF activity heatmap
  par(mar=c(4,0.25,0.25,0.25),mgp=c(0,0.2,0))
  col = rev(colorRampPalette(brewer.pal(11,"PRGn"))(100))
  image(t(as.matrix(selTF.Mes)), col=col, xaxt="n", yaxt="n")
  grid(col="white",ny=nrow(selTF.Mes), nx=0, lty=1,lwd=1)
  box(col="black")
  axis(side=1, at=0, labels="TF activity", cex.axis=0.7, tick=F,las=2)

  # Cell line CRC proportion barplots
  par(mar=c(4,0.25,0.25,0.1),mgp=c(0,0.2,0),xaxs="i",yaxs="i")
  barplot(t(propObs[,5:ncol(propObs)]),horiz=T,las=2,cex.axis=0.9,cex.lab=0.9,cex.names=0.8,col=c("#006837","#40004B","#31a354"),border=NA,yaxt="n",xaxt="n")
  axis(side=1, at=c(0,50,100), labels=c(0,50,100), cex.axis=0.9, tick=F,las=2)
  axis(side=3, at=50, labels="Proportion(%)", cex.axis=0.8, line= -0.3, tick=F)
  box(col="black")

  # Cell CRC oncoprints
  par(mar=c(4,0.25,0.25,2.5),mgp=c(0,0.2,0),xaxs="i",yaxs="i")
  image(t(tmpCRCcell),col=c("#f0f0f0","black"),xaxt="n",yaxt="n")
  axis(side=4, at=seq(0,1,length.out=nrow(tmpCRCcell)), labels=rownames(tmpCRCcell), las=2, cex.axis=0.6, tick=F)
  grid(nx=ncol(tmpCRCcell),ny=nrow(tmpCRCcell),col="white", lty=1,lwd=1); box()
  axis(side=1, at=seq(0,1,length.out=ncol(tmpCRCcell)), labels=colnames(tmpCRCcell), cex.axis=0.6, tick=F,las=2)
  abline(v=seq(0,1,length.out=ncol(tmpCRCcell))[c(sigBordersCell[1], sum(sigBordersCell[1:2]))]+0.02,lwd=1)
  #axis(side=3, at=0.5, labels="CRC frequency (Cells)", cex.axis=1, tick=F)

  # #Internal KD analysis
  # par(mar=c(2.9,8,1,8),mgp=c(0,0.2,0))
  # col = rev(colorRampPalette(brewer.pal(9,"YlOrRd"))(100))
  # image(as.matrix(kd.NB.fc), col=col, xaxt="n", yaxt="n")
  # grid(col="white",nx=nrow(kd.NB.fc), ny=0, lty=1,lwd=1)
  # box(col="black")
  # axis(side=2, at=0, labels=paste("KD","validation",sep="\n"), cex.axis=0.8, tick=F,las=2)
  # valMes = seq(0,1,length.out=nrow(kd.NB.fc))[which(kd.NB.fc$log2FC < 0)]
  # valNonMes = seq(0,1,length.out=nrow(kd.NB.fc))[which(kd.NB.fc$log2FC > 0)]
  # axis(side=1, at=valMes, labels=rownames(kd.NB.fc)[which(kd.NB.fc$log2FC < 0)], cex.axis=0.75, tick=F,las=2,col.axis="#40004B")
  # axis(side=1, at=valNonMes, labels=rownames(kd.NB.fc)[which(kd.NB.fc$log2FC > 0)], cex.axis=0.75, tick=F,las=2,col.axis="#006837")
  # rm(valMes,valNonMes)
  # 
  # ## Legends for KD sensitivty
  # par(mar=c(1.5,4.5,2.5,17.5),mgp=c(1.5,0.5,0))
  # key = sort(runif(100, min = min(kd.NB.fc, na.rm=T), max = max(kd.NB.fc,na.rm=T)),decreasing = T)
  # image(as.matrix(key), col=colorRampPalette(brewer.pal(9,"YlOrRd"))(100), xaxt="n", yaxt="n");box()
  # title("Cell sensitivity", line=0.3, cex.main=1)
  # axis(side=1,at=c(0,.5,1),labels=rev(round(key[c(1,length(key)/2,length(key))])))
  # axis(side=2,at=0,labels=paste("High Mes","cell death",sep="\n"),tick=F,las=1)
  # axis(side=4,at=0,labels=paste("High NonMes","cell death",sep="\n"),tick=F,las=1)
  # rm(key)

  # Signature legend
  par(mar=c(0,0,1,0),xaxs="i",yaxs="i", xpd=F)
  plot.new(); legend("left", legend=c("MYCN-Amp","Mesenchymal","LR/MYCN-nonAmp","HR-MYCN-nonAmp"), bty="n", fill=c("#006837","#40004B","#31a354","#78c679"),
  x.intersp=0.3, cex=0.9, border=c("#006837","#40004B","#31a354","#78c679"))

  ## Legends for TF activity
  par(mar=c(1.5,9,2.5,3.8),mgp=c(1.5,0.5,0),xpd=F)
  key = sort(runif(100, min= min(selTF.Mes,na.rm=T), max=max(selTF.Mes, na.rm=T)),decreasing = T)
  image(as.matrix(key), col= colorRampPalette(brewer.pal(11,"PRGn"))(100), xaxt="n", yaxt="n");box()
  title("TF activity", line=0.3, cex.main=1)
  axis(side=1,at=c(0,.5,1),labels=rev(round(key[c(1,length(key)/2,length(key))])))
  axis(side=4,at=0,labels=paste("High in","Mes",sep="\n"),tick=F,las=1)
  axis(side=2,at=0,labels=paste("High in","NonMes",sep="\n"),tick=F,las=1)
dev.off()  
