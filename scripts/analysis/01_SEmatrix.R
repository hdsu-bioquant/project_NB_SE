##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##              Computes SE signal over consensus regions                     ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
#library(GenomicRanges)
# Parse Cmd line.
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)
#read paths of out files
out.RDS <- args[1]
out.RDS
out.matrix <- args[2]
out.matrix
#read paths of consensus SE
SE.path <- args[3]
SE.path
#read paths of files with signal of SE over consensus regions
averageOverBed.path <- args[4:length(args)]
averageOverBed.path

library(limma)
#
#setwd("/icgc/dkfzlsdf/analysis/B080/quintera/Projects/neuroblastoma/")
# averageOverBed.tumors <- list.files(path="data/tumors" ,pattern="*_bigWigAverageOverBed.txt", full.names = T, recursive = T)
# averageOverBed.cells <- list.files(path="data/cells" ,pattern="*_bigWigAverageOverBed.txt", full.names = T, recursive = T)
# averageOverBed.path <- c(averageOverBed.tumors, averageOverBed.cells)
# SE.path <- "results/tumors/H3K27ac/SE/rose_noH3K4me3/tumors_consensus_SE_list_filtered.bed"

#extract sample names

sample.names <- basename(averageOverBed.path)
sample.names <- sub("_H3K27ac_bigWigAverageOverBed.txt$","", sample.names)

names(averageOverBed.path) <- sample.names

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##              Read signal of consensus SE for each sample                   ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
#read consensus SE
consensusSE <- read.table(SE.path, stringsAsFactors=FALSE)
colnames(consensusSE) <- c("chr", "start", "end", "ID", "nMerged", "IDsMerged")


#read signal over each sample
sample.signalSE <- lapply(names(averageOverBed.path), function(sample){
  sample.df <- read.table(averageOverBed.path[sample], header=FALSE, stringsAsFactors=FALSE)
  colnames(sample.df) <- c("name", "size", "covered", "sum", "meanZeroes", "mean")
  #extract the mean contribution of each base to the total signal in this region
  sample.df <- sample.df[match(consensusSE$ID, sample.df$name),5, drop=FALSE]
  rownames(sample.df) <- consensusSE$ID
  colnames(sample.df) <- sample
  return(sample.df)
})
sample.signalSE <- do.call(cbind, sample.signalSE)

#remove negative values
sample.signalSE[sample.signalSE < 0] <- 0
SE.matrix <- as.matrix(sample.signalSE)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                           Quantile normalization                           ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
#perform quantile normalization with limma package
SE.quantnorm.matrix <- normalizeQuantiles(SE.matrix)
colnames(SE.quantnorm.matrix)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                                Save matrix                                 ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
#save as RDS
saveRDS(SE.quantnorm.matrix, file=out.RDS)

#save as txt
write.table(SE.quantnorm.matrix, file=out.matrix)

