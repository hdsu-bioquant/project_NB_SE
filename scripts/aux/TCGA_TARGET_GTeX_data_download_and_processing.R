options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

url1 = as.character(args[1])
url2 = as.character(args[2])
url3 = as.character(args[3])
outpath = as.character(args[4])

#-------------------------------------------------------------------------------------------------
# DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"
# url1 = "https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz"
# url2 = "https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz"
# url3 = "https://toil.xenahubs.net/download/probeMap/gencode.v23.annotation.gene.probemap"
# outpath = paste0(DATAPATH,"db/")
#-------------------------------------------------------------------------------------------------

library(data.table)
library(R.utils)

#--------------------------------------------------------------------------------------
# The compiled data from TCGA, TARGET and GTeX was downloaded from
# The "gene expression RNAseq" data named "RSEM tpm" was downloaded from below
# https://xenabrowser.net/datapages/?cohort=TCGA%20TARGET%20GTEx&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443 
# https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_rsem_gene_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#
# Tumor metadata
#   {
#    "cohort": "TCGA TARGET GTEx", 
#
#    "label": "RSEM tpm", 
#    "wrangling_procedure": "Data (file names: *.rsem_genes.results) are downloaded, tpm values are extracted, log2(x+0.001) transformed, and combined.", 
#    "version": "2016-09-03",
#    "dataProducer": "UCSC TOIL RNA-seq recompute", 
#    "type": "genomicMatrix", 
#    "unit": "log2(tpm+0.001)", 
#    "dataSubType": "gene expression RNAseq",
#    ":probeMap":"/probeMap/gencode.v23.annotation.gene.probemap"
#    }
#
# Sample information metadata
#   {
#    "type":"clinicalMatrix",
#    "version":"2016-09-15",
#    "cohort":"TCGA TARGET GTEx",
#    "dataSubType":"phenotype",
#    ":clinicalFeature":"TcgaTargetGTEX_phenotype_clinFeature",
#    "label":"TCGA TARGET GTEX selected phenotypes"
#    }
#
# Gene annotation metadata
#  {
#    "type":"probeMap",
#   "version":"2018-07-13",
#    "assembly":"hg38",
#    "xdb":"genecode Release 23 (GRCh38.p3)",
#    "url": "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_23/gencode.v23.annotation.gtf.gz",
#    "note":"comprehensive annotation, gene only",
#    "label":"Gene: GENCODE 23 comprehensive (hg38) e.g. ENSG00000223972.5",
#    "userLevel":"basic",
#    "idType":"gene"
#  }
#--------------------------------------------------------------------------------------
url1
allexp = data.frame(fread(url1), stringsAsFactors = F)
allexp[1:5,1:5]
rownames(allexp) = allexp$sample
allexp = allexp[,-1]

#dim(allexp)
#[1] 60498 19131

# Clean sample assignment
url2
sampInfo = data.frame(fread(url2), stringsAsFactors = F)
sampInfo = sampInfo[,-2]
colnames(sampInfo) = c("Sample", "Primary_Disease_or_Tissue", "PrimarySite", "SampleType","Gender", "Study")
sampInfo = sampInfo[sampInfo$SampleType %in% c("Normal Tissue", "Solid Tissue Normal", "Primary Tumor", 
                                               "Primary Solid Tumor", "Additional - New Primary",
                                               "Primary Blood Derived Cancer - Bone Marrow", "Primary Blood Derived Cancer - Peripheral Blood"),]
sampInfo$SampleType[sampInfo$SampleType %in% c("Normal Tissue", "Solid Tissue Normal")]="Normal"
sampInfo$SampleType[sampInfo$SampleType %in% c("Primary Tumor","Primary Solid Tumor", "Additional - New Primary","Primary Blood Derived Cancer - Bone Marrow", "Primary Blood Derived Cancer - Peripheral Blood")]="Tumor"

# Correct for long names, names with typos and missing values
sampInfo$PrimarySite[sampInfo$PrimarySite == "Sympathetic\xcaNervous System"] = "Sympathetic nervous system"
sampInfo$PrimarySite[sampInfo$PrimarySite == "Adrenal Gland"] = "Adrenal gland"
sampInfo$PrimarySite[sampInfo$PrimarySite == "Head and Neck region"] = "Head and Neck"
sampInfo$PrimarySite[sampInfo$PrimarySite == "Soft tissue,Bone"] = "Bone soft tissue"
sampInfo$PrimarySite[sampInfo$PrimarySite == "Lining of body cavities"] = "Body cavity lining"
sampInfo$PrimarySite[sampInfo$PrimarySite == ""] = sapply(strsplit(sampInfo$Primary_Disease_or_Tissue[sampInfo$PrimarySite == ""]," -",fixed=T),function(x)x[1])

# Keep samples (tumors or normals) from each tissue type with > 50 samples
a = paste(sampInfo$Study, sampInfo$Primary_Disease_or_Tissue, sampInfo$SampleType, sep="|")
df = data.frame(counts = table(a))
df = df[order(df[,1]),]
df = df[df[,2] >= 50,]
sampInfo = sampInfo[a %in% as.character(df[,1]),]
rm(a,df)

# Making expression data and clinical data compatible 
colnames(allexp) = gsub(".", "-", colnames(allexp), fixed=T)
cmn = intersect(colnames(allexp), sampInfo$Sample)
allexp = allexp[, colnames(allexp) %in% cmn]
sampInfo = sampInfo[sampInfo$Sample %in% cmn,]
sampInfo = sampInfo[match(colnames(allexp), sampInfo$Sample),]
rm(cmn)

# Gene mapping
gencode = data.frame(fread(url3), stringsAsFactors=F)
cmn = intersect(gencode$id, rownames(allexp))
allexp = allexp[rownames(allexp) %in% cmn,]
gencode = gencode[gencode$id %in% cmn,]
gencode = gencode[match(rownames(allexp), gencode$id),]
rownames(allexp) = paste(gencode$id,gencode$gene,sep="|")
rm(cmn,gencode)

# Saving RDS objects
saveRDS(allexp,   paste0(outpath, "TCGA_TARGET_GTex/TcgaTargetGtex_log2_fpkm.RDS"))
saveRDS(sampInfo, paste0(outpath, "TCGA_TARGET_GTex/TcgaTargetGtex_sample_information.RDS"))


