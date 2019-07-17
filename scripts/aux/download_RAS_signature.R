options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

url = as.character(args[1])
outpath = as.character(args[2])

#-------------------------------------------------------------------------------------------------
# DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"
# 
# url = "https://static-content.springer.com/esm/art%3A10.1186%2F1755-8794-3-26/MediaObjects/12920_2010_161_MOESM3_ESM.XLS"
# outpath  = paste0(DATAPATH, "db/")
#-------------------------------------------------------------------------------------------------

library(readxl)
library(utils)

# This script downloads the RAS target genes identified in the article - 
# Loboda, A., Nebozhyn, M., Klinghoffer, R., Frazier, J., Chastain, M., Arthur, W., Roberts, B., Zhang, T., Chenard, M., Haines, B. and Andersen, J., 2010.
# A gene expression signature of RAS pathway dependence predicts response to PI3K and RAS pathway inhibitors and expands the population of RAS pathway activated tumors. 
# BMC medical genomics, 3(1), p.26. (https://doi.org/10.1186/1755-8794-3-26)

download.file(url, destfile = paste0(outpath,"publicGeneSigs/12920_2010_161_MOESM3_ESM.XLS"))
ras = as.data.frame(read_excel(paste0(outpath,"publicGeneSigs/12920_2010_161_MOESM3_ESM.XLS")), sheet = 1)

ras = list(upreg = ras$`RAS Signature UP`, dowreg = ras$`RAS Signature DOWN`)

ras = lapply(ras, function(x){
  x = x[! is.na(x)]
  x = unique(as.character(x))
})

saveRDS(ras, paste0(outpath,"publicGeneSigs/ras_target_genes.RDS"))
