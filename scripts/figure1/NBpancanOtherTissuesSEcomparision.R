options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

expdata = as.character(args[1])
sampdesp = as.character(args[2])
nbSE = as.character(args[3])
otherSE = as.character(args[4])
outpath = as.character(args[5])

#-------------------------------------------------------------------------------------------------
# DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/"
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
library(rio)
library(ggplot2)
library(ggrepel)
library(writexl)
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
dev.off()

library(tidyverse)
sd <- lapply(names(vals), function(id){
  tibble(type = id,
         Sample = names(vals[[id]]),
         val = vals[[id]])
})
sd <- bind_rows(sd)
write_xlsx(list(`Extended Data figure 1f` = sd), 
           path = "results/figure_source_data/Extended_Data_figure_1f.xlsx")

rm(id,vals,tmp,a,cols,ids,l)
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
#identical(anno$`File Name`, id)
#identical(anno$`File Name`, gsub(".csv","", basename(files)))

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

######################################
# NB cell line SE from Groningen et.al 
######################################

se2 = import(file = "https://media.nature.com/original/nature-assets/ng/journal/v49/n8/extref/ng.3899-S4.xlsx")
se2 = se2[,2:4]
se2 = se2[-1,]
colnames(se2) = se2[1,]
se2 = se2[-1,]
se2 = makeGRangesFromDataFrame(se2)

####################################
# NB cell line SE from Boeva et.al
####################################
se3 = import(file = "https://media.nature.com/original/nature-assets/ng/journal/v49/n9/extref/ng.3921-S3.xlsx")
se3 = se3[-1,]
colnames(se3) = se3[1,]
se3 = se3[-1,]
se3 = se3[, c(1:3,5,14,15,19)]

# Chosing a SE specific gene (single gene) based on highest correlation
se3$`Gene(s)` = sapply(strsplit(se3$`Gene(s)`, "|", fixed=T), function(x)x[1])
se3$`FC Score Group I over Group II` = round(log2(as.numeric(se3$`FC Score Group I over Group II`)),3)
se3$`Wilcoxon p-value (two sided test)` = signif(as.numeric(se3$`Wilcoxon p-value (two sided test)`),3)
se3 = se3[, -ncol(se3)]
tmp = makeGRangesFromDataFrame(se3[,1:3])
mcols(tmp) = se3[4:ncol(se3)]
se3 = tmp
rm(tmp)

###################################################
# Add NB specific SEs to multiple tissue derived SE
###################################################

SElist$Neuroblastoma = readRDS(nbSE)
SElist$Groningen_et_al = se2
SElist$Boeva_et_al = se3
rm(se2, se3)

colnames(anno) = c("Name", "Subtype", "Class")
NB_anno = data.frame(Name = c("Neuroblastoma", "Groningen_et_al","Boeva_et_al"),
                     Subtype = "Neuroblastoma",
                     Class = "Tumor", stringsAsFactors = F)
anno = rbind(anno, NB_anno)

anno$Type = rep("Primary", nrow(anno))
anno$Type[anno$Name %in% c("H1", "HMEC", "HUVEC", "DND41","GM12878", 
                          "Jurkat", "MM1S", "RPMI-8402", "K562", "u87", 
                          "HCC1954", "MCF-7", "HeLa", "HCT-116", "VACO_400", 
                          "VACO_503", "H2171", "LnCAP", "Groningen_et_al", "Boeva_et_al")] = "Cell-line"

rm(NB_anno)

# Find the fraction overlap of SE from other tissues with NB super enhancers (at least 50%)
# Code adapted from https://support.bioconductor.org/p/72656/

SEoverlap = rep(NA, length(SElist))

if(identical(names(SElist), anno$Name))
{
  names(SEoverlap) = anno$Subtype
  names(SEoverlap)[which(anno$Type == "Cell-line")] = paste(names(SEoverlap)[which(anno$Type == "Cell-line")],
                                                            anno$Name[which(anno$Type == "Cell-line")], sep=" | ")
}

overlap_threshold = 0.25

for( i in 1: length(SElist))
{
    print(paste("At step: ",i," of ",length(SElist),sep=""))
  
    refGR  = SElist$Neuroblastoma
    testGR = SElist[[i]]
      
    # overlap with atleast 50% coverage
    hits = findOverlaps(query = refGR, subject = testGR)
    overlaps = pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
    percentOverlap = width(overlaps) / width(testGR[subjectHits(hits)])
    hits = hits[percentOverlap > overlap_threshold]

    # % of neuroblastoma super enhancers found in trhe overlap
    SEoverlap[i] = round(c(length(unique(queryHits(hits)))/length(refGR)) * 100, 2)

    rm(refGR, testGR, hits ,overlaps, percentOverlap)
}
rm(i)

SEoverlap = data.frame(Label = names(SEoverlap), Overlap = as.numeric(SEoverlap))
SEoverlap = cbind(anno, SEoverlap)

# Main figure 
g <- ggplot(SEoverlap, aes(x = reorder(Label, Overlap), y = Overlap)) + 
  theme_bw(base_size = 9) +
  labs(x = "", y = "% of all SE in neuroblastoma (n=1973)") + 
  geom_point() +
  #scale_fill_manual(values=c("#397FB9", "#994F9F")) +
  geom_text_repel(data = subset(SEoverlap, Overlap > 10.5), aes(label = Label), size=2) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour="black"), 
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.8,0.8),
        axis.text.x = element_blank())

write_xlsx(list(`Figure 1b` = g$data), 
           path = "results/figure_source_data/Figure_1b.xlsx")
ggsave(paste0(outpath,"figure1/SE_multiple_tissue_overlap_to_NB_SE_MainFig.pdf"), g, width=3, height=3)

# Supplementary figure with sample names
g <- ggplot(SEoverlap, aes(x = reorder(Label, Overlap), y = Overlap)) + 
  theme_bw(base_size = 9) + ylim(0,106) + coord_flip() +
  labs(x = "", y = "% of all SE in neuroblastoma (n=1973)") + 
  geom_bar(stat="identity") +
  #scale_fill_manual(values=c("#397FB9", "#994F9F")) +
  geom_text(aes(label = paste0(round(Overlap),"%"), y = Overlap+6), 
            vjust = 0.5, hjust = 1, angle = 0, size=2) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour="black"), 
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.8,0.8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) 

write_xlsx(list(`Extended Data figure 1d` = g$data), 
           path = "results/figure_source_data/Extended_Data_figure_1d.xlsx")

ggsave(paste0(outpath,"sup_figure1/SE_multiple_tissue_overlap_to_NB_SE_SupplFig.pdf"), g, width=4.2, height=5)

rm(SEoverlap, overlap_threshold, anno, SElist)
