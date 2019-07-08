##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##            Download misc data and install missing R packages               ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# Parse Cmd line.
options(echo=TRUE) # if you want to see commands in output file
args <- commandArgs(TRUE)

#from pipeline
datapath <- args[1]
datapath
setwd(datapath)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                     paths of misc data to download                         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# hsapiens_genes_path <- "db/misc/EnsDb_Hsapiens_v75_genes.RDS"
# path_hic <- "db/hic/GSE63525_K562_HiCCUPS_looplist.txt"
# path_TADs <- "db/TADs/hESC_domains_hg19.RDS"

hsapiens_genes_path <- args[2]
path_hic            <- args[3]
path_TADs           <- args[4]

hsapiens_genes_path
path_hic
path_TADs


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                        install missing R packages                          ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
devtools::install_github('andquintero/bratwurst', ref='dev_hdsu', upgrade = "never")
devtools::install_github("thomasp85/patchwork", upgrade = "never")

#BiocManager::install("viper")


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", version = "3.8")
# BiocManager::install("Homo.sapiens", version = "3.8")
#BiocManager::install("EnsDb.Hsapiens.v75", version = "3.8")
#install.packages("pkgconfig")
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(Homo.sapiens)
#install.packages("pheatmap")
#install.packages("R.utils")


#BiocManager::install("GenomicInteractions", version = "3.8")
# devtools::install_github("cran/amap", ref = "40fe5e59fa3d84790f568438b752356d52ad15e0", upgrade = "never")
# BiocManager::install("DiffBind")

library(EnsDb.Hsapiens.v75)
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          Granges of human genes                            ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
edb <- EnsDb.Hsapiens.v75
organism(edb)
# supportedFilters(edb)
## Change the seqlevels style form Ensembl (default) to UCSC:
seqlevelsStyle(edb) <- "UCSC"
hsapiens_genes <- genes(edb)
hsapiens_genes
table(hsapiens_genes$gene_biotype)

saveRDS(hsapiens_genes, file = hsapiens_genes_path)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            Download HiC data                               ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# set ftp url to RNA-seq data
# download data
ftp.url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_K562_HiCCUPS_looplist.txt.gz"
download.file(url = ftp.url, path_hic)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                   Download TADs and write as Granges                       ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# set ftp url to RNA-seq data
# download data
data_url <- "http://compbio.med.harvard.edu/modencode/webpage/hic/hESC_domains_hg19.bed"
path_TADs_bed <- file.path(dirname(path_TADs), basename(data_url))
download.file(url = data_url, path_TADs_bed)

# read bed and write Granges
TADs <- read.table(path_TADs_bed, stringsAsFactors=FALSE)
colnames(TADs) <- c("chr", "start", "end")
TADs <- makeGRangesFromDataFrame(TADs)
TADs
saveRDS(TADs, path_TADs)


