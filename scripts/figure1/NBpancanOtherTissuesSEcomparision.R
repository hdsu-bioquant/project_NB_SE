options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

expdata = as.character(args[1])
sampdesp = as.character(args[2])
nbSE = as.character(args[3])
otherSE = as.character(args[4])
outpath = as.character(args[5])

#-------------------------------------------------------------------------------------------------
# DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"
# expdata  = paste0(DATAPATH,"db/TCGA_TARGET_GTex/TcgaTargetGtex_log2_fpkm.RDS")
# sampdesp = paste0(DATAPATH,"db/TCGA_TARGET_GTex/TcgaTargetGtex_sample_information.RDS")
# nbSE     = paste0(DATAPATH,"analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS")
# otherSE  = paste0(DATAPATH,"db/SEmultiTisuues")
# outpath  = paste0(DATAPATH,"results/")
#-------------------------------------------------------------------------------------------------

library(readxl)
library(pheatmap)
library(RColorBrewer)
library(GenomicRanges)

#---------------------------------------------------------------------------------------------------
# Comparision of NB super enhancer target gene expression across multiple tumor and normal tissues
#---------------------------------------------------------------------------------------------------

# Load the RNAseq dataset (tpm values) for TARGET TCGA and GTex
allexp = readRDS(expdata)
sampInfo = readRDS(sampdesp)

# Super enhancers
SE = readRDS(nbSE)
SE = unique(mcols(SE)$target_SYMBOL)

# Subset
symbols = sapply(strsplit(rownames(allexp), "|", fixed=T),function(x){x[2]})
allexpSE = allexp[symbols %in% SE,]
medallexpSE = apply(allexpSE, 2, function(x){median(x, na.rm=T)})

rm(symbols,SE)

# Plot expression of NB target SE across multiple tumors and tissues
pdf(paste0(outpath,"sup_figure1/TCGA_TARGET_GTex_NB_SE_target_genes_expression_comparision.pdf"),width=9.3,height=4)
par(mar=c(17.3,2.5,0.5,0.25), mgp=c(1.5,0.5,0), yaxs="i", xaxs="i", cex=0.65)

  id = paste(sampInfo$Primary_Disease_or_Tissue, sampInfo$SampleType, sep="|")
  vals = medallexpSE

  vals = split(vals,id)
  tmp = sapply(vals,median)
  vals = vals[order(tmp)]
  l = which(tmp > quantile(tmp,0.75))
  l = sort(tmp[names(tmp) %in% names(l)])
  l = names(l[1])

  a = boxplot(vals,plot=F)
  cols = sapply(strsplit(a$names,"|",fixed=T),function(x)x[2])
  ids  = sapply(strsplit(a$names,"|",fixed=T),function(x)x[1])
  
  boxplot(vals, ylab = "log2 (TPM + 0.001)", las=2, outline=F, frame=F, xaxt="n", ylim = c(-10,4),
           col = ifelse(cols == "Normal", "forestgreen", "firebrick"),
         boxlty=1, whisklty=1, staplelty=1, medlwd=1, boxwex=0.6, boxlwd=1, whisklwd=1, staplelwd=1)
  axis(1,at=1:max(a$group),labels=ids,line= -0.3, las=2, tick=F)
  legend("bottomright", legend=c("Normal Tissues","Primary Tumors"), fill=c("forestgreen", "firebrick"), border = c("forestgreen", "firebrick"), x.intersp=0.3, y.intersp=0.8, bty="n")
  abline(v=which(l == a$names)-0.5, h=0, lwd=1, col= "grey30")
  text(x=which(l == a$names)-0.5, y = -8, labels = "top 25%", srt = 90, pos = 4, col = "grey30")
  rm(id,vals,tmp,a,cols,ids,l)
dev.off()

rm(allexp, allexpSE, sampInfo, medallexpSE)

#-------------------------------------------------------------------------
# Comparision of NB super enhancer with SE identified from other tissues
# data from Richard Young's paper on SE cell identity
#-------------------------------------------------------------------------

# Read data and corres[ponding annotation
files = list.files(otherSE, full.names = T)
id    = gsub(".csv","",sapply(strsplit(files, "/"),function(x){x[length(x)]}))
anno  = as.data.frame(read_xlsx(paste0(otherSE,"/mmc2.xlsx"), col_names = TRUE))
anno  = anno[which(anno$`Rep. Dataset` == "X" | anno$`Is a cancer` == "Cancer"),]
anno  = anno[order(anno$`Rep. Dataset`),]
anno$Subtype[is.na(anno$Subtype)] = anno$`Tissue Type`[is.na(anno$Subtype)]
tmp = rep(NA,nrow(anno))
tmp[anno$`Rep. Dataset` == "X"] = "Normal"
tmp[anno$`Is a cancer` == "Cancer"] = "Tumor"
anno = anno[,c(1,3)]
anno = cbind(anno,Class=tmp)
rm(tmp)

# anno$`File Name`[!anno$`File Name` %in% id] 
#[1] "UCSD_Psoas" "NHDF_Ad"    "VACO_9m"    "Panc1" 
# These files are missing in the annotation file, hence not using them

cmn = intersect(id, anno$`File Name`)
files = files[which(id %in% cmn)]
anno = anno[which(anno$`File Name` %in% cmn),]
id = id[which(id %in% cmn)]
rm(cmn)

files = files[match(anno$`File Name`,id)]
id = id[match(anno$`File Name`,id)]

# View(cbind(anno, files, id))

#table(anno$Class)
#Normal  Tumor 
#24     15 

# Making a list with all the tissue specific super enhancers and Neuroblastoma super enhancers
SElist = vector("list", length(id))
names(SElist) = id
for(i in 1: length(SElist))
{
  dat = read.csv(files[i], skip=1, header=F, stringsAsFactors = F)
  dat = dat[dat$V7 == 1,] #Selecting only super enhancers
  SElist[[i]] = GRanges(seqnames=dat$V2, ranges=IRanges(start=dat$V3, end=dat$V4), Gene= dat$V5)
  print(paste("At step: ",i," of ",length(SElist),sep=""))
  rm(dat)
}
rm(i,id,files)
SElist$Neuroblastoma = readRDS(nbSE)

# Find the fraction overlap of SE from other tissues with NB super enhancers (at least 50%)
# Code adapted from https://support.bioconductor.org/p/72656/

SEoverlap = rep(NA, length(SElist))
for( i in 1: length(SElist))
{
    print(paste("At step: ",i," of ",length(SElist),sep=""))
  
    refGR  = SElist[[i]]
    testGR = SElist$Neuroblastoma
    
    # overlap with atleast 50% coverage
    hits = findOverlaps(refGR, testGR)
    overlaps = pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
    percentOverlap = width(overlaps) / width(testGR[subjectHits(hits)])
    hits = hits[percentOverlap > 0.5]

    # Compute jaccard's coefficient i.e intersection (at least 50%) divided by union
    SEoverlap[i] = round(length(hits) / length(reduce(c(refGR, testGR))), 3)
    
    # SOME OTHER METRIC TO FIND  FRACTION OVERLAP (NOT USED !!)
    #SEoverlap[i] = round(length(hits) / c(length(refGR) + length(testGR)), 3)
    #SEoverlap[i] = round(c(length(hits) / length(testGR)) * 100, 2)

    rm(refGR, testGR, hits ,overlaps, percentOverlap)
}
rm(i)

collab = anno$Subtype
collab[which(anno$Class == "Tumor")] = paste(collab[which(anno$Class == "Tumor")],anno$`File Name`[which(anno$Class == "Tumor")],sep="|")
collab = c(collab,"Neuroblastoma")
names(SEoverlap) = collab
SEoverlap = sort(SEoverlap)
rm(collab, anno, SElist)

# Main figure 
pdf(paste0(outpath,"figure1/SE_multiple_tissue_overlap_to_NB_SE_MainFig.pdf"),width=1.5, height=1.5)
  par(mar=c(2.5,2.5,0.25,0.25), mgp=c(1.6,0.5,0),xaxs="i", yaxs="i", cex=0.6)
  colr = rep("#5aa02c", length(SEoverlap))
  colr[which(names(SEoverlap) %in% grep("|",names(SEoverlap),value=T,fixed=T))]="#ff7f2a"
  colr[length(colr)] = "black"
  
  plot(SEoverlap, pch=20, cex=1.2, ylab= "Fraction overlap", xlab="Samples", ylim = c(-0.03,1.03),xlim=c(0,41),las=2,col=colr,frame=F)
  legend("topleft", legend=c("Primary tissue","Cell-lines"), fill=c("#5aa02c","#ff7f2a"), border=c("#5aa02c","#ff7f2a"),bty="n",x.intersp=0.3, y.intersp=0.8)
  abline(h=min(tail(SEoverlap,2)),lty=1,col="grey30")
  text(x=length(SEoverlap)/2, y=min(tail(SEoverlap,2)), paste0("max overlap ~",round(min(tail(SEoverlap,2))*100),"%"),cex=0.8,pos=3)
  rm(colr)
dev.off()

# Supplementary figure with sample names
pdf(paste0(outpath,"sup_figure1/SE_multiple_tissue_overlap_to_NB_SE_SupplFig.pdf"),width=4.5, height=2.5)
  par(mar=c(11,2.5,0.25,0.25), mgp=c(1.5,0.5,0),xaxs="i", yaxs="i", cex=0.6)
  colr = rep("#5aa02c", length(SEoverlap))
  colr[which(names(SEoverlap) %in% grep("|",names(SEoverlap),value=T,fixed=T))]="#ff7f2a"
  colr[length(colr)] = "grey"
  barplot(SEoverlap, ylab= "Fraction overlap", ylim = c(-0.03,1.03),las=2,col=colr,border=NA)
  legend("topleft", legend=c("Primary tissue","Cell-lines"), fill=c("#5aa02c","#ff7f2a"), border=c("#5aa02c","#ff7f2a"),bty="n",x.intersp=0.3, y.intersp=0.8)
  abline(h=min(tail(SEoverlap,2)),lty=1,col="grey30")
  text(x=length(SEoverlap)/2, y=min(tail(SEoverlap,2)), paste0("max overlap ~",round(min(tail(SEoverlap,2))*100),"%"),cex=0.8,pos=3)
  rm(colr)
dev.off()

rm(SEoverlap)
