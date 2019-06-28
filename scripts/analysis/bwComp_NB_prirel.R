library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(WriteXLS)

DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"

consSE           = paste0(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS')
tumorNMF         = paste0(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K4_Wmatrix_Wnorm.RDS')
mesTFact         = paste0(DATAPATH, "analysis/tumor/VIPER/MES_TFactivity.RDS")
relvspri_Expdiff = paste0(DATAPATH, 'analysis/tumor/Rel_vs_Pri/GSEA/relapseVSprimary.rnk')
relvspri_BWdiff  = paste0(DATAPATH, 'analysis/tumor/Rel_vs_Pri/bwComp/')
outpath          = paste0(DATAPATH, "analysis/tumor/Rel_vs_Pri/bwComp/")

#--------

mesTFact = readRDS(mesTFact)
mesTFact = data.frame(Genes = names(mesTFact$es$nes), MesTFactivity= mesTFact$es$nes,  stringsAsFactors = F)
summary(mesTFact)

relvspri_Expdiff = read.delim(relvspri_Expdiff, header=T, stringsAsFactors = F)
summary(relvspri_Expdiff)

consSE = readRDS(consSE)

SEtumor_MES = readRDS(tumorNMF)
SEtumor_MES = rownames(SEtumor_MES)[SEtumor_MES[,'MES'] > quantile(SEtumor_MES[,'MES'], 0.8)]
SEtumor_MES = unique(mcols(consSE[mcols(consSE)$ID %in% SEtumor_MES])$target_SYMBOL)
SEtumor_MES = data.frame(Genes = SEtumor_MES, Type="Mes specific", stringsAsFactors = F)

bw_files = list.files(relvspri_BWdiff, full.names = T)
bw_files = grep(".bw", bw_files, value=T, fixed=T)
relvspri_BWdiff = mclapply(bw_files, function(x){import(x)}, mc.cores = 4)
names(relvspri_BWdiff) = sapply(strsplit(bw_files, "/"), function(x) gsub("_bigwig_log2Diff.bw", "", x[length(x)]))

#-------

topHits = mclapply(relvspri_BWdiff, function(x){
  a = findOverlaps(x, consSE)
  
  df = data.frame(Genes = mcols(consSE)$target_SYMBOL[subjectHits(a)], logFC = score(x)[queryHits(a)])
  df = sapply(split(df$logFC, df$Genes), median)
  
  df = data.frame(Genes = names(df), logFC = as.numeric(df), stringsAsFactors = F)
  
  df = merge(x = df, y = relvspri_Expdiff)
  df = merge(x = df, y = mesTFact, all.x = T)
  df = merge( x= df, y = SEtumor_MES, all.x = T)
  
  df = df[which(df$logFC < -0.5 & df$Score < -1  | 
                  df$logFC > 0.5 & df$Score > 1),]
  df = df[order(df$logFC),]
  
  return(df)
}, mc.cores = 4)

WriteXLS(topHits, paste0(outpath,"bestCandidates_Rel_vs_Pri.xls"), row.names = F, AdjWidth = T, BoldHeaderRow = T, FreezeRow = 1)

# All overlap
set1 = Reduce(intersect,list(topHits$NSP032_RM01_PT01$Genes, 
                       topHits$NSP032_RT01_PT01$Genes, 
                       topHits$NSP032_RT02_PT01$Genes, 
                       topHits$NSP032_RT03_PT01$Genes))


set2 = Reduce(intersect, list(topHits$NSP032_RT01_PT01$Genes, 
                       topHits$NSP032_RT02_PT01$Genes, 
                       topHits$NSP032_RT03_PT01$Genes))

# Removin RT-02 as its of poor quality
set3 = Reduce(intersect, list(topHits$NSP032_RM01_PT01$Genes, 
                       topHits$NSP032_RT01_PT01$Genes, 
                       topHits$NSP032_RT03_PT01$Genes))

# Removing the relapse metastasis
set4 = Reduce(intersect, list(topHits$NSP032_RT01_PT01$Genes, 
                       topHits$NSP032_RT03_PT01$Genes))
