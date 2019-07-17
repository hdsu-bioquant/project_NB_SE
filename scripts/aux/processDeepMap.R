options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

projDesp = as.character(args[1])
cellAnno = as.character(args[2])
cellExpr = as.character(args[3])
valsKD   = as.character(args[4])
probsKD  = as.character(args[5])
outpath  = as.character(args[6])

#-------------------------------------------------------------------------------------------------
# DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"

# projDesp = "https://ndownloader.figshare.com/files/15023474"
# cellAnno = "https://ndownloader.figshare.com/files/15023525"
# cellExpr = "https://ndownloader.figshare.com/files/15023486"
# valsKD   = "https://ndownloader.figshare.com/files/15023465"
# probsKD  = "https://ndownloader.figshare.com/files/15023459"
# outpath  = paste0(DATAPATH, "db/")
#-------------------------------------------------------------------------------------------------

####################################################################################################
# AK Sharma - 19th June, 2019
# DepMap Version - 19Q2 (see README.txt)

# This script downloads the following DepMap 19Q2 data -

# CERES KD scores
# CERES KD probability
# CCLE cell line gene expression

# The script then performs various data wrangling and creates a
# consensus set of samples across different genomics platforms. It generate ready to use .RDS files
# The entire script should run ideally in < 5 min

####################################################################################################

# Required libraries and root directory definition
library(data.table)
library(dplyr)

#--------------------------------------------------------------------------------------------------#
# Cells - Download the README.txt file for the 19Q1 version                                        #
#--------------------------------------------------------------------------------------------------#

download.file(url = projDesp, destfile = paste0(outpath, "DeepMap19Q2/README.txt"))

#--------------------------------------------------------------------------------------------------#
# Cells - Download cell line annotations
#--------------------------------------------------------------------------------------------------#

cells.anno = fread(cellAnno)
cells.anno = cells.anno[,c(1:3,6:9)]
cells.anno = data.frame(cells.anno[,1:3], apply(cells.anno[,4:ncol(cells.anno)], 2, as.factor))
saveRDS(cells.anno, file = paste0(outpath, "DeepMap19Q2/cellAnnotation.RDS"))

#--------------------------------------------------------------------------------------------------#
# Cells - Download expression data (log(TPM) + 1) (only protein coding genes)
#--------------------------------------------------------------------------------------------------#
cells.exp = fread(cellExpr)

sampID = c(as.matrix(cells.exp[,1]))
genes = sapply(strsplit(colnames(cells.exp), " (", fixed=T), function(x)x[1])
genes = genes[- 1]
cells.exp = as.data.frame(t(cells.exp[, -1]))
colnames(cells.exp) = sampID
cells.exp = cbind(Gene = genes, cells.exp)
cells.exp$Gene = as.character(cells.exp$Gene)

# Remove only zeros containing rows
sum(rowSums(cells.exp[,2:ncol(cells.exp)]) == 0)
#[1] 21
cells.exp = cells.exp[rowSums(cells.exp[,2:ncol(cells.exp)]) > 0,]

# Checking if there are 'NAa in gene names or if the gene names are duplicated
sum(is.na(cells.exp$Gene))
# [1] 0
sum(duplicated(cells.exp$Gene))
#[1] 0

# Summarizing duplicated genes (Since no duplicates in this version, commmenting out)
# dupIDs = cells.exp$Gene[duplicated(cells.exp$Gene)]
# dupExp = cells.exp %>% filter(Gene %in% dupIDs)
# dupExp = dupExp %>% group_by(Gene) %>% summarise_all(median)
# nondup = cells.exp %>% filter(!Gene %in% dupIDs)
# cells.exp = rbind(nondup,dupExp)

# Cleaning and orderingt he expression matrix
rownames(cells.exp) = cells.exp$Gene
cells.exp =  cells.exp[, -1]

cmn = intersect(cells.anno$DepMap_ID, colnames(cells.exp))
cells.anno.exp = cells.anno[cells.anno$DepMap_ID %in% cmn,]
cells.exp = cells.exp[,colnames(cells.exp) %in% cmn]
cells.anno.exp = cells.anno.exp[match(colnames(cells.exp), cells.anno.exp$DepMap_ID),]
identical(colnames(cells.exp), cells.anno.exp$DepMap_ID)
cells.exp = list(expression = cells.exp, annotation = droplevels(cells.anno.exp))

saveRDS(cells.exp, file = paste0(outpath, "DeepMap19Q2/cellExpression.RDS"))
rm(dupIDs, dupExp, nondup, genes, sampID, cells.anno.exp, cmn)

#--------------------------------------------------------------------------------------------------------------#
# Cells - CRISPR KD CERES scores
#
# CERES data with principle components strongly related to known batch effects removed, then shifted
# and scaled per cell line so the median nonessential KO effect is 0 and the median essential KO effect is -1.
#---------------------------------------------------------------------------------------------------------------#

cells.kd.ceres = fread(valsKD)

sampID = c(as.matrix(cells.kd.ceres[,1]))
genes = sapply(strsplit(colnames(cells.kd.ceres), " (", fixed=T), function(x)x[1])
genes = genes[- 1]
cells.kd.ceres = as.data.frame(t(cells.kd.ceres[, -1]))
colnames(cells.kd.ceres) = sampID
cells.kd.ceres = cbind(Gene = genes, cells.kd.ceres)
rownames(cells.kd.ceres) = as.character(cells.kd.ceres$Gene)
cells.kd.ceres =  cells.kd.ceres[, -1]

cmn = intersect(cells.anno$DepMap_ID, colnames(cells.kd.ceres))
cells.anno.kd = cells.anno[cells.anno$DepMap_ID %in% cmn,]
cells.kd.ceres = cells.kd.ceres[,colnames(cells.kd.ceres) %in% cmn]
cells.anno.kd = cells.anno.kd[match(colnames(cells.kd.ceres), cells.anno.kd$DepMap_ID),]
identical(colnames(cells.kd.ceres), cells.anno.kd$DepMap_ID)
cells.kd.ceres = list(essentiality = cells.kd.ceres, annotation = droplevels(cells.anno.kd))

saveRDS(cells.kd.ceres, file = paste0(outpath, "DeepMap19Q2/cellKnockdownCERES.RDS"))
rm(sampID, genes, cells.anno.kd, cmn)

#--------------------------------------------------------------------------------------------------#
# Cells - CRISPR KD Probability
#
# Probability that knocking out the gene has a real depletion effect using gene_effect.
#--------------------------------------------------------------------------------------------------#

cells.kd.prob = fread(probsKD)

sampID = c(as.matrix(cells.kd.prob[,1]))
genes = sapply(strsplit(colnames(cells.kd.prob), " (", fixed=T), function(x)x[1])
genes = genes[- 1]
cells.kd.prob = as.data.frame(t(cells.kd.prob[, -1]))
colnames(cells.kd.prob) = sampID
cells.kd.prob = cbind(Gene = genes, cells.kd.prob)
rownames(cells.kd.prob) = as.character(cells.kd.prob$Gene)
cells.kd.prob = cells.kd.prob[, -1]

cmn = intersect(cells.anno$DepMap_ID, colnames(cells.kd.prob))
cells.anno.kd = cells.anno[cells.anno$DepMap_ID %in% cmn,]
cells.kd.prob = cells.kd.prob[,colnames(cells.kd.prob) %in% cmn]
cells.anno.kd = cells.anno.kd[match(colnames(cells.kd.prob), cells.anno.kd$DepMap_ID),]
identical(colnames(cells.kd.prob), cells.anno.kd$DepMap_ID)
cells.kd.prob = list(essentialityP = cells.kd.prob, annotation = droplevels(cells.anno.kd))

saveRDS(cells.kd.prob, file = paste0(outpath, "DeepMap19Q2/cellKnockdownProb.RDS"))
rm(sampID, genes, cells.anno.kd, cmn)
