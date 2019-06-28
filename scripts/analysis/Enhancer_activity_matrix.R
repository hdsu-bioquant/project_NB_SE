library(GenomicRanges)
library(rtracklayer)
library(limma)
library(parallel)

DATAPATH  = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"

annot = paste0(DATAPATH, "annotation/annotation_tumor.RDS")
annot = readRDS(annot)
annot <- annot[annot$avail.ChIPseq,]
bw_paths  = paste0(DATAPATH, "data/tumor/chipseq/H3K27ac/bw/", annot$ProjectID,"_H3K27ac.bw")
names(bw_paths) <- annot$ProjectID
#names(bw_paths)

SE_target = paste0(DATAPATH, "analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS")
SE <- readRDS(SE_target)

enhancers = paste0(DATAPATH, "analysis/tumor/chipseq/H3K27ac/consensusEnhancers/tumor_H3K27ac_noH3K4me3_consensusEnhancers.bed")
enhancers <- read.table(enhancers, as.is=TRUE)
#head(enhancers)
enhancers <- enhancers[,1:4]
colnames(enhancers) = c('chr','start','end','ID')
enhancers = makeGRangesFromDataFrame(enhancers, keep.extra.column = TRUE)
enhancers <- enhancers[!enhancers %over% SE]

#----------------------------------------------
detectCores()
#[1] 32

length(bw_paths)
#[1] 60

# Splitting the jobs in 30 cores each twice
#-----------------------------------------------

enhancers_total_Signal_1_30 <- do.call(cbind, mclapply(bw_paths[1:30], function(bw_path){
  bw         <- import(bw_path)
  bw         <- subsetByOverlaps(bw, enhancers)
  bw_en_over <- findOverlaps(bw, enhancers)
  bw         <- bw[from(bw_en_over)]
  bw$ID      <- enhancers$ID[to(bw_en_over)]
  
  sapply(split(bw, bw$ID), function(SE_sig){
    sum(SE_sig$score)})
}, mc.cores = 30))

enhancers_total_Signal_31_60 <- do.call(cbind, mclapply(bw_paths[31:60], function(bw_path){
  bw         <- import(bw_path)
  bw         <- subsetByOverlaps(bw, enhancers)
  bw_en_over <- findOverlaps(bw, enhancers)
  bw         <- bw[from(bw_en_over)]
  bw$ID      <- enhancers$ID[to(bw_en_over)]
  
  sapply(split(bw, bw$ID), function(SE_sig){
    sum(SE_sig$score)})
}, mc.cores = 30))

identical(rownames(enhancers_total_Signal_1_30), rownames(enhancers_total_Signal_31_60))
#[1] TRUE

enhancers_total_Signal = cbind(enhancers_total_Signal_1_30, enhancers_total_Signal_31_60)
rm(enhancers_total_Signal_1_30, enhancers_total_Signal_31_60)

# Creating enhancer activity matrix
df1 = enhancers_total_Signal
df1[df1 < 0] = 0
df1 = log2(df1+1)
#plot(density(as.matrix(df1), na.rm=T))
df1 = normalizeQuantiles(df1)
df1 = round(df1, 3)
#plot(density(as.matrix(df1), na.rm=T))
#boxplot(df1, outline=F, las=2)

df1 = data.frame(ID = rownames(df1), df1, stringsAsFactors = F, check.names = F)
df2 = data.frame(enhancers, stringsAsFactors = F)
enh = merge(df2, df1)
enh = enh[, -c(5:6)]

saveRDS(enh, paste0(DATAPATH, "analysis/tumor/chipseq/H3K27ac/consensusEnhancerstumor_H3K27ac_noH3K4me3_consensusEnhancers_SignalScore.RDS"))
