options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

KDdata = as.character(args[1])
KDres = as.character(args[2])

#-------------------------------------------------------------------------------------------------
# DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"

# KDdata = paste0(DATAPATH,"data/cells/CRCsiRNAknockdown/")
# KDres = paste0(DATAPATH,"analysis/cells/crcGIEMSAkd/nbKDinhouse.RDS")
#-------------------------------------------------------------------------------------------------

# Load require packages
library("readxl")

#------------------------------
# In house KD data compilation
#------------------------------

# GIMEN
giemsaGIMEN = as.data.frame(read_xlsx(path=paste0(KDdata,"GIEMSA_analysis_T72_T96_GIMEN.xlsx"),sheet = 1, col_names = T))
giemsaGIMEN = giemsaGIMEN[-c(1:2),]
rownames(giemsaGIMEN) = gsub(" ", "", giemsaGIMEN[,1], fixed=T)
colnames(giemsaGIMEN) = c("Genes",  "96hr_Rep1", "96hr_Rep2")
rownames(giemsaGIMEN)[2] = "CEBPB#2" # was a typo in te file
giemsaGIMEN = giemsaGIMEN[,-1]
giemsaGIMEN[,1] = as.numeric(giemsaGIMEN[,1])
giemsaGIMEN[,2] = as.numeric(giemsaGIMEN[,2])

round(cor(giemsaGIMEN, method="spearman"),2)
#            96hr_Rep1 96hr_Rep2
#96hr_Rep1      1.00      0.56
#96hr_Rep2      0.56      1.00

id = sapply(strsplit(rownames(giemsaGIMEN), "#", fixed=T),function(x)x[1])
id = split(giemsaGIMEN, id)
giemsaGIMEN = t(sapply(id, function(x){apply(x,2,median)}))
colnames(giemsaGIMEN) = paste0("GIMEN_",colnames(giemsaGIMEN))
giemsaGIMENsum = data.frame(GIMEN_96hr= apply(giemsaGIMEN[,1:2],1,mean))
rm(id)

# NMB
giemsaNMB = as.data.frame(read_xlsx(path=paste0(KDdata, "GIEMSA_analysis_T72_T96_NMB.xlsx"), sheet = 1, col_names = T))
giemsaNMB = giemsaNMB[-c(1:2),]
giemsaNMB = giemsaNMB[,-c(2:3)]
rownames(giemsaNMB) = gsub(" ", "", giemsaNMB[,1], fixed=T)
colnames(giemsaNMB) = c("Genes",  "96hr_Rep1", "96hr_Rep2")
rownames(giemsaNMB)[2] = "CEBPB#2" # was a typo in te file
giemsaNMB = giemsaNMB[,-1]
giemsaNMB[,1] = as.numeric(giemsaNMB[,1])
giemsaNMB[,2] = as.numeric(giemsaNMB[,2])

round(cor(giemsaNMB, method="spearman"),2)
#            96hr_Rep1 96hr_Rep2
#96hr_Rep1      1.00      0.64
#96hr_Rep2      0.64      1.00

id = sapply(strsplit(rownames(giemsaNMB), "#", fixed=T),function(x)x[1])
id = split(giemsaNMB, id)
giemsaNMB = t(sapply(id, function(x){apply(x,2,median)}))
colnames(giemsaNMB) = paste0("NMB_",colnames(giemsaNMB))
giemsaNMBsum = data.frame(NMB_96hr= apply(giemsaNMB[,1:2],1,mean))
rm(id)

# Mes vs NonMes knockdown fold change
inhouse_nbKD_GIMENvsNMB = data.frame(KDlogFC = log2(giemsaGIMENsum) - log2(giemsaNMBsum))
colnames(inhouse_nbKD_GIMENvsNMB ) = "log2FC"
saveRDS(inhouse_nbKD_GIMENvsNMB, file = KDres)


# Not used because siRNA quality was poor for these cell-lines
if(FALSE)
{
      # SHEP
      giemsaSHEP = as.data.frame(read_xlsx(path=paste0(KDdata,"GIEMSA_analysis_T72_T96_SHEP.xlsx"),sheet = 1, col_names = T))
      rownames(giemsaSHEP) = gsub(" ", "", giemsaSHEP$Genes, fixed=T)
      colnames(giemsaSHEP) = c("Genes", "96hr_Rep1", "96hr_Rep2")
      rownames(giemsaSHEP)[2] = "CEBPB#2" # was a typo in te file
      giemsaSHEP = giemsaSHEP[,-1]

      #round(cor(giemsaSHEP, method="spearman"),2)
      #      96hr_Rep1 96hr_Rep2
      #96hr_Rep1  1.00  0.22
      #96hr_Rep2  0.22  1.00

      id = sapply(strsplit(rownames(giemsaSHEP), "#", fixed=T),function(x)x[1])
      id = split(giemsaSHEP,id)
      giemsaSHEP = t(sapply(id, function(x){apply(x,2,median)}))
      colnames(giemsaSHEP) = paste0("SHEP_",colnames(giemsaSHEP))
      giemsaSHEPsum = data.frame(SHEP_96hr= apply(giemsaSHEP[,1:2],1,mean))
      rm(id)

      # SKNAS
      giemsaSKNAS = as.data.frame(read_xlsx(path=paste0(KDdata,"GIEMSA_analysis_T72_T96_SKNAS.xlsx"),sheet = 1, col_names = T))
      rownames(giemsaSKNAS) = gsub(" ", "", giemsaSKNAS[,1], fixed=T)
      colnames(giemsaSKNAS) = c("Genes",  "96hr_Rep1", "96hr_Rep2")
      rownames(giemsaSKNAS)[2] = "CEBPB#2" # was a typo in te file
      giemsaSKNAS = giemsaSKNAS[,-1]

      # round(cor(giemsaSKNAS, method="spearman"),2)
      #            96hr_Rep1 96hr_Rep2
      #96hr_Rep1      1.00      0.12
      #96hr_Rep2      0.12      1.00

      id = sapply(strsplit(rownames(giemsaSKNAS), "#", fixed=T),function(x)x[1])
      id = split(giemsaSKNAS, id)
      giemsaSKNAS = t(sapply(id, function(x){apply(x,2,median)}))
      colnames(giemsaSKNAS) = paste0("SKNAS_",colnames(giemsaSKNAS))
      giemsaSKNASsum = data.frame(SKNAS_96hr= apply(giemsaSKNAS[,1:2],1,mean))
      rm(id)

      # MHH
      giemsaMHH = as.data.frame(read_xlsx(path=paste0(KDdata,"GIEMSA_analysis_T72_T96_MHH.xlsx"),sheet = 1, col_names = T))
      rownames(giemsaMHH) = gsub(" ", "", giemsaMHH[,1], fixed=T)
      colnames(giemsaMHH) = c("Genes",  "96hr_Rep1", "96hr_Rep2")
      rownames(giemsaMHH)[2] = "CEBPB#2" # was a typo in te file
      giemsaMHH = giemsaMHH[,-1]

      # round(cor(giemsaMHH, method="spearman"),2)
      #            96hr_Rep1 96hr_Rep2
      #96hr_Rep1      1.00     -0.11
      #96hr_Rep2     -0.11      1.00

      id = sapply(strsplit(rownames(giemsaMHH), "#", fixed=T),function(x)x[1])
      id = split(giemsaMHH, id)
      giemsaMHH = t(sapply(id, function(x){apply(x,2,median)}))
      colnames(giemsaMHH) = paste0("MHH_",colnames(giemsaMHH))
      giemsaMHHsum = data.frame(MHH_96hr= apply(giemsaMHH[,1:2],1,mean))
      rm(id)

      # NBLS
      giemsaNBLS = as.data.frame(read_xlsx(path=paste0(KDdata,"GIEMSA_analysis_T72_T96_NBLS.xlsx"),sheet = 1, col_names = T))
      rownames(giemsaNBLS) = gsub(" ", "", giemsaNBLS[,1], fixed=T)
      colnames(giemsaNBLS) = c("Genes",  "96hr_Rep1", "96hr_Rep2")
      rownames(giemsaNBLS)[2] = "CEBPB#2" # was a typo in te file
      giemsaNBLS = giemsaNBLS[,-1]

      #round(cor(giemsaNBLS, method="spearman"),2)
      #            96hr_Rep1 96hr_Rep2
      #96hr_Rep1      1.00      0.14
      #96hr_Rep2      0.14      1.00

      id = sapply(strsplit(rownames(giemsaNBLS), "#", fixed=T),function(x)x[1])
      id = split(giemsaNBLS, id)
      giemsaNBLS = t(sapply(id, function(x){apply(x,2,median)}))
      colnames(giemsaNBLS) = paste0("NBLS_",colnames(giemsaNBLS))
      giemsaNBLSsum = data.frame(NBLS_96hr= apply(giemsaNBLS[,1:2],1,mean))
      rm(id)

#####

# Read SHEP CRC list file
shepCRC = read.delim("/icgc/dkfzlsdf/analysis/B080/herrmanc/Projects/Neuroblastoma_Frank/data_new/celllines_for_analysis/SH-EP/H3K27ac/CRC/SE2GENE/CRC_peaks/rose_noH3K4me3/sp/500000/B087_SH-EP_none_H3K27ac_CRC_SCORES_SE2GENE.txt", stringsAsFactors=F, header=F)[,1]
pattern = c("[", "]", "'", " ")
for(j in pattern){shepCRC = gsub(j, "", shepCRC, fixed=T)}
rm(j, pattern)
crcLen = length(shepCRC)
shepCRC = unlist(strsplit(shepCRC,","))
cnt = (sort(table(shepCRC))/crcLen)*100
shepCRC = unique(shepCRC)

df = data.frame(SHEP = giemsaSHEPsum, GIMEN=giemsaGIMENsum, SKNAS=giemsaSKNASsum, NMB=giemsaNMBsum, NBLS=giemsaNBLSsum, MHH=giemsaMHHsum)
round(cor(df, method="spearman"),2)

df = as.matrix(df[,c(1,4)])
df = df[rownames(df) %in% shepCRC,]

pdf("~/Documents/tmp/SHEP_CRC_KD_sensitivity_comparision.pdf", width=3.5, height=1.5)
par(mar=c(4,3,0.5,0.5), mgp=c(1.5,0.5,0), cex=0.7, xaxs="i", yaxs="i")
barplot(t(df),beside=T, las=2, ylab = "Sensitivity", legend=T)
dev.off()

}
