library(GenomicInteractions)
library(tidyverse)
library(edgeR)

setwd("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/")

# Paths data
path_tumor_annot <- "annotation/annotation_tumor.RDS"
path_SE <- "analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed"
path_hichip <- c("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/data/cells/hichip/mango/SK-N-AS_HiChIP_mango.all",
                 "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/data/cells/hichip/mango/CLB-GA_HiChIP_mango.all")
names(path_hichip) <- sub("_HiChIP_mango.all", "", basename(path_hichip))

path_SE_signal  <- "analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE_SignalScore.RDS"
path_gene_exprs <- "data/tumor/rnaseq/exprs/tumor_RNAseq_Counts_Matrix.RDS"

# Paths external data
path_hic  <- "db/hic/GSE63525_K562_HiCCUPS_looplist.txt"
path_TADs <- "db/TADs/hESC_domains_hg19.RDS"

# Paths misc data
path_hsapiens_genes <- "db/misc/EnsDb_Hsapiens_v75_genes.RDS"

# path_hichip_clbga <- ""
# path_hichip_sknas


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                                Params                                      ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
NCORES <- 10
# For each SE find genes around this window size 
SE_neighbor_genes_window <- 500000

hichip_mango_FDR   <- 0.05


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                Read human genes Granges and find promoters                 ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

hsapiens_genes <- readRDS(path_hsapiens_genes)
hsapiens_genes$id <- hsapiens_genes$symbol
# expand Units:
proms <- promoters(hsapiens_genes, upstream = 5000, downstream = 1000) 



##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          Read Consensus SEs                                ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
#read consensus SE
consensusSE <- read.table(path_SE, stringsAsFactors=FALSE)
colnames(consensusSE) <- c("chr", "start", "end", "ID", "nMerged", "IDsMerged")
SE <- makeGRangesFromDataFrame(consensusSE, keep.extra.column = TRUE)
SE$id <- SE$ID

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##  Assign target gene based on correlation of signal and gene expression     ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##Read data - SE signal - Gene expression - TADs - Genes##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# Read data - Tumor annotation
tumor_annot <- readRDS(path_tumor_annot)
# Read data - SE signal
SE_signal  <- readRDS(path_SE_signal)
# Read data - Gene expression
gene_exprs <- readRDS(path_gene_exprs)
# Read data - TADs
hESC_domain_hg19 <- readRDS(path_TADs)
# Read data - Genes
hsapiens_genes <- readRDS(path_hsapiens_genes)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## Keep samples with available RNAseq and ChIPseq data  ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# Tumor annotation
tumor_annot <- tumor_annot[which(tumor_annot$avail.ChIPseq & tumor_annot$avail.RNAseq),]
dim(tumor_annot)
# SE signal
SE_signal   <- SE_signal[,tumor_annot$ProjectID]
dim(SE_signal)
# Gene expression
gene_exprs  <- gene_exprs[,tumor_annot$ProjectID]
dim(gene_exprs)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## Gene expression - Normalize - filter genes wo annot  ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# RNAseq read-counts normalization
gene_exprs_norm <- DGEList(counts  = gene_exprs, 
                           samples = colnames(gene_exprs), 
                           remove.zeros = TRUE)
gene_exprs_norm$samples$group <- rep(1,ncol(gene_exprs_norm))

#identify the cpm value that corresponds to atleast 20 reads in the smallest library
cutoff <- (20/min(gene_exprs_norm$samples$lib.size)) * 10^6  
#cutoff cpm should be present in atleast 25% of samples
isexpr <- rowSums(cpm(gene_exprs_norm) > cutoff) > round(0.25 * ncol(gene_exprs_norm)) 
table(isexpr)
gene_exprs_norm <- gene_exprs_norm[isexpr, ,keep.lib.sizes=FALSE]
# normalize using voom
gene_exprs_norm <- calcNormFactors(gene_exprs_norm, method="TMM")
gene_exprs_norm <- voom(gene_exprs_norm, plot=TRUE)
gene_exprs_norm <- gene_exprs_norm$E
hist(gene_exprs_norm)
dim(gene_exprs_norm)

# use only Ensembl ID to match gene names
rownames(gene_exprs_norm) <- sapply(strsplit(rownames(gene_exprs_norm), "\\."), "[[", 1)
# Filter H.sapiens Granges and expression matrix
# Keep only shared genes
gene_exprs_norm <- gene_exprs_norm[rownames(gene_exprs_norm) %in% hsapiens_genes$gene_id, ]
dim(gene_exprs_norm)
hsapiens_genes  <- hsapiens_genes[hsapiens_genes$gene_id %in% rownames(gene_exprs_norm)]
length(hsapiens_genes)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## Overlap SE and genes regions - Corr. Signal-Exprs    ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# find overlap between SE regions and gene regions, 
# filter out the regions which overlap 
# (duplicates will be present because some SE overlap with multiple genes and vice versa)
overlap <- findOverlaps(resize(SE, 
                               width = (2*SE_neighbor_genes_window) + width(SE),
                               fix   = "center"), 
                        hsapiens_genes)

SE_signal_exprs_corr <- mclapply(1:length(overlap), function(i){
  hit <- overlap[i]
  
  SE_hit   <- SE[queryHits(hit)]
  gene_hit <- hsapiens_genes[subjectHits(hit)]
  
  ct <- cor.test(SE_signal[SE_hit$ID,], gene_exprs_norm[gene_hit$gene_id,], method = "spearman")
  
  data.frame(SE      = SE_hit$ID,
             ENSEMBL = gene_hit$gene_id, 
             SYMBOL  = gene_hit$symbol,
             method  = "correlation",
             pval    = ct$p.value,
             corr    = ct$estimate,
             row.names = NULL, 
             stringsAsFactors = FALSE)
  
}, mc.cores = NCORES)
SE_signal_exprs_corr <- do.call(rbind, SE_signal_exprs_corr)
head(SE_signal_exprs_corr)
length(unique(SE_signal_exprs_corr$SE))

##––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##      Find if SE and target are in the same TAD       ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## annotate with TADS
load(TADS.path)
hESC_domain_hg19
# produce vector which tells for all regions for which correlations where computed if the se region and the prospective gene are located in the same TAD
ov.se.tad <- findOverlaps(SE,hESC_domain_hg19)
ov.genes.tad <- findOverlaps(promoters(GENES,upstream=1,downstream=1),hESC_domain_hg19)
ov.se.genes.tad <- merge(as.data.frame(ov.se.tad),as.data.frame(ov.genes.tad),by=1,all=TRUE)
tmp <- apply(ov.se.genes.tad,1,function(x) {x[2]==x[3]})
tmp[is.na(tmp)] <- TRUE
is.same.tad <- rep(TRUE,length(SE))
is.same.tad[ov.se.genes.tad[,1]] <- tmp


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                      Read HiChIP mango results                             ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

hichip <- do.call(rbind, lapply(path_hichip, function(path_mango){
  read_tsv(path_mango, col_names = F)  %>% 
    filter(X8 <= hichip_mango_FDR)
}))

# a <- read_tsv(paste(datadir, "/HiChIP-SK-N-AS.mango.all", sep = "", collapse = ""), col_names = F) %>% filter(X8 <= FDR)
# b <- read_tsv(paste(datadir, "/HiChIP-CLB-GA.mango.all", sep = "", collapse = ""), col_names = F) %>% filter(X8 <= FDR)
# hichip <- rbind(a,b)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                             Read HiC data                                  ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
loops <- read_tsv(gzfile(path_hic))

# loops <- read.table(gzfile(path_hic))
loops <- na.omit(loops)
loops$chr1 <- paste0("chr", loops$chr1)
loops$chr2 <- paste0("chr", loops$chr2)
colnames(hichip) <- colnames(loops)[1:8]
colnames(hichip) <- colnames(loops)[1:ncol(hichip)]

interA <- GRanges(loops$chr1, IRanges(loops$x1, loops$x2))
interB <- GRanges( loops$chr2, IRanges(loops$y1, loops$y2))
HiC <- GenomicInteractions(interA, interB,counts =  loops$o)


interA <- GRanges(hichip$chr1, IRanges(hichip$x1, hichip$x2))
interB <- GRanges( hichip$chr2, IRanges(hichip$y1, hichip$y2))
hichipI <- GenomicInteractions(interA, interB,counts =  hichip$o)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                             Read HiC data                                  ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## Calling targets

#names(genes.hg19) <- genes.hg19@elementMetadata$gid
annotation.features = list(SE = SE,
                           #SE5kb = SE.5kb,
                           promotor = proms
                           #gene5kb = genes.hg19.5kb
)


annotateInteractions(HiC, annotation.features)
annotateInteractions(hichipI, annotation.features)
#categoriseInteractions(loopsI)


# Annotated using interactions:
  

plotInteractionAnnotations(HiC)



x <- readRDS("/icgc/dkfzlsdf/analysis/B080/quintera/Projects/neuroblastoma/results/tumors/H3K27ac/SE/rose_noH3K4me3/tumors_consensusSE_GRanges_noH3K4me3.RDS")





