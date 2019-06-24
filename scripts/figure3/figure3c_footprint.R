library(limma)
library(parallel)
library(GenomicRanges)
library(cowplot)
library(readr)
library(tidyverse)
library(ggrepel)
library(tibble)
library(preprocessCore)


# filemap <- read.table("/icgc/dkfzlsdf/analysis/B080/quintera/Projects/neuroblastoma/src/fig3c/filemap.txt", header = TRUE, stringsAsFactors = FALSE)
# 
# apply(filemap, 1, function(x){
#   file.copy(from = x[3], to = x[2])
# })

setwd("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/")
# params <- list(SE = "analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS",
#                cellline1_foot = "data/cells/atacseq/footprint/CLB-GA_footprints_calls_GrangesList.RDS",
#                cellline2_foot = "data/cells/atacseq/footprint/SK-N-AS_footprints_calls_GrangesList.RDS",
#                MES_activity   = "analysis/tumor/VIPER/MES_TFactivity.RDS")

params <- list(SE = "analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS",
               cellline1_foot = "data/cells/atacseq/footprint/KELLY_footprints_calls_GrangesList.RDS",
               cellline2_foot = "data/cells/atacseq/footprint/SK-N-AS_footprints_calls_GrangesList.RDS",
               MES_activity   = "analysis/tumor/VIPER/MES_TFactivity.RDS")


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                             Read data                                      ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
SE <- readRDS(params$SE)
MES_activity <- readRDS(params$MES_activity)
names(MES_activity$es$p.value)[MES_activity$es$p.value < 0.05]

MES_activity <- MES_activity$es$nes


footprint_path <- setNames(c(params$cellline1_foot, params$cellline2_foot),
         sub("_.*", "", basename(c(params$cellline1_foot, params$cellline2_foot))))

footprint_list <- lapply(footprint_path, function(fppath){
  readRDS(fppath)
})

# use only common footprints
footIDs <- Reduce(intersect, lapply(footprint_list, names))
names(footIDs) <- footIDs

footprint_list <- lapply(footIDs, function(footID){
  lapply(footprint_list, '[[', footID)
})


gc()

footprint_list[[1]]

#================================================================================#
#                         Read CRC list                                          #
#================================================================================#
crcpath <- c(list.files("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/data/tumor/chipseq/H3K27ac/CRC", pattern = ".*_H3K27ac_ROSE_noH3K4me3_500Kb_CRC.txt$", full.names = TRUE),
             list.files("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/data/cells/chipseq/H3K27ac/CRC", pattern = ".*_H3K27ac_ROSE_noH3K4me3_500Kb_CRC.txt$", full.names = TRUE))

CRCvector <- unique(do.call(c, lapply(crcpath, function(path){
  x <- lapply(strsplit(readLines(path), split = '\t'), '[[', 1)
  do.call(c, strsplit(gsub('\\[\'|\'\\]', '', x), "\\', \\'"))
})))
sort(CRCvector)


# cellline1_foot <- readRDS(params$cellline1_foot)
# cellline2_foot <- readRDS(params$cellline2_foot)
#================================================================================#
#                                   Params                                       #
#================================================================================#

filter      <-  TRUE
#loadSEfile  = T
#ncores      <- 20
#Tumor       = F
quantileNormal <-  FALSE



#================================================================================#
#                           T test purity score                                  #
#================================================================================#

min.mot = ifelse(filter,10,100)

TTT <- mclapply(footprint_list,function(t) {
  #print(names(t))
  cl1 = t[[1]]
  cl2 = t[[2]]
  ##
  #A =normalizeQuantiles(data.frame(cl1=cl1$mfrow.purity,cl2=cl2$mfrow.purity))
  A = data.frame(cl1=cl1$mfrow.purity,cl2=cl2$mfrow.purity)
  ##
  cl1$score = A$cl1
  cl2$score = A$cl2
  ##
  if (filter) {
    ov1 = findOverlaps(cl1,SE)
    ov2 = findOverlaps(cl2,SE)
    ##
    if (length(queryHits(ov1))>0 & length(queryHits(ov2))>0) {
      cl1 = cl1[queryHits(ov1)]
      cl2 = cl2[queryHits(ov2)]
    }
  }
  ##                      
  ## get top 10% but not less than 0.75
  ##combined_purity <- c(cl1$mfrow.purity, cl2$mfrow.purity)
  ##cutOff <- min(head(rev(sort(combined_purity)), n = round(0.1*length(cl1))))
  ##
  combined_score <- c(cl1$score, cl2$score)
  cutOff <- min(head(rev(sort(combined_score)), n = round(0.1*length(cl1))))
  ##           
  ##                       cutOff <- min(head(rev(sort(combined_purity)), n = 500))
  ##
  if(cutOff < 0.7 ){
    cutOff <- 0.7
  }else{
    # check if the large cutoff reduces the number of hits to much:
    n <- length(cl1[cl1$score>=  cutOff | cl2$score >= cutOff]$score)
    iterated <- F
    
    while(n < min.mot & cutOff >= 0.705 & is.finite(cutOff)){
      cutOff   <- cutOff - 0.05
      n        <- length(cl1[cl1$score>=  cutOff | cl2$score >= cutOff]$score)
      iterated <- T
      #message(paste(n, cutOff))
    }
    if(n >= min.mot){
      if(iterated) message("iterativly defined cutoff")
    }
    if(cutOff < 0.7 ){
      cutOff <- 0.7
      
    }    
  }
  
  x     = cl1[cl1$score>=  cutOff | cl2$score >= cutOff]$score
  y     = cl2[cl1$score>= cutOff | cl2$score >= cutOff]$score
  ##
  n = length(x)
  ##
  m <- matrix(c(x,y), ncol = 2, byrow = F)
  if (length(x) >= min.mot) { # at least 10 motives
    if(quantileNormal){
      m <- normalize.quantiles(m)
    }
    #print(m)
    res <- t.test(m[,1],m[,2],paired=TRUE)
    
    # print(res)
    # print(res$p.value)
    # print(res$statistic)
    
    res <- data.frame(statistic.t = res$statistic,
                      p.value     = res$p.value)
    
    #res <- lapply(res, function(x) x)
    
    #print(lapply(as.list((res)), class))
    
    #print((res))
    
    #print(data.frame(list((res))))
    #res <- as.data.frame(as.list(unlist(res))) 
    
    res$symb <- unique(t[[1]]$ID)
    return(res)
  } else {
    return(NA)
  }
},mc.cores=20)
message("finished t-testing")
#TTT
head(TTT)
#tf = gsub('-calls.all.csv','',gsub('.*- ','',basename(cl1.files)))
#names(TTT) = tf

lapply((TTT[[1]]), class)


#==============================================================================#
#                         organize by t statistic                              #
#==============================================================================#

# Remove NAs
length(TTT)
TTT = TTT[!is.na(TTT)]
length(TTT)
# make data frame
df <- do.call(rbind, TTT)
dim(df)

head(df)

#==============================================================================#
#                              Format TF SYMBOL                                #
#==============================================================================#
# remove string after _
df$symb <- sapply(strsplit(df$symb, "_"), "[[", 1)
# to upper
df$symb <- toupper(df$symb)
# remove duplicated strings
df$symb <- sub(pattern = "^(.+?)\\1+$", "\\1", df$symb)
df$symb

# sort by t stat
df <- df[order(df$statistic.t),]
df$pos <-1:nrow(df)

x <- df[df$symb %in% CRCvector,]

# split by TF and find the max
df.perTF = split(df,df$symb)
lapply(df.perTF[1:2], as.data.frame)
as.data.frame(df.perTF$IRF8)
as.data.frame(df.perTF$MAFK)
#Keep only TFs that are consistently different
max.T = sapply(df.perTF,function(x) {
  if(sum(x$statistic.t>0)*sum(x$statistic.t<0) ==0) {
    max.t = sign(x$statistic.t[1])*max(abs(x$statistic.t))
  }else {
    max.t = NA
  }
  return(max.t)
})


med.T = sapply(df.perTF,function(x) {
  return(mean(x$statistic.t,  na.rm = T))
})


#lookup <- read_tsv(file.path('data','lookup_motif_TF.txt'), col_names = F)
#d <- sapply(df$tf, function(n){
#  n <- strsplit(n, split = " ")[[1]][2]
#  m <- match(n, lookup$X1)
#  print(m)
#  m <- m[1]
#  lookup[m,]$X2
#  })
#df$names <- d

df$logt = sign(df$statistic.t)*log(abs(df$statistic.t))

ddf = data.frame(pos=1:length(sort(max.T)),val=sort(max.T))
if(quantileNormal){
  ddf = data.frame(pos=1:length(sort(med.T)),val=sort(med.T))
}

head(ddf)


# flip around:
ddf.flipped = data.frame(pos=1:length(sort(max.T)),val=sort(max.T))
head(ddf.flipped)
identical(ddf.flipped,ddf)


if(quantileNormal){
  ddf.flipped = data.frame(pos=1:length(sort(med.T)),val=sort(med.T))
}
ddf.flipped$pos <- rev(ddf.flipped$pos)
ddf.flipped$val <- -ddf.flipped$val

m <- ggplot(ddf, aes(pos, val)) + geom_point() + 
  geom_text_repel(aes(label=as.character(rownames(ddf)))) +
  xlab("")+  ylab("t.statistic")

m.flipped <- ggplot(ddf.flipped, aes(pos, val)) + geom_point() + 
  geom_text_repel(aes(label=as.character(rownames(ddf)))) +
  xlab("")+
  ylab("t.statistic")


ddf$name <- rownames(ddf)
fName <- ""
if(filter){
  message("filtered")
  fName <- ".seFiltered"
}
if(quantileNormal){
  fName <- paste0(fName, "_quantile")
}

# write.table(ddf, file = paste0('/icgc/dkfzlsdf/analysis/B080/quintera/Projects/neuroblastoma/src/fig3c/atac.tstat_',cl1.name,'_',cl2.name,fName,'.csv'), row.names = F)
# ggsave(paste0('/icgc/dkfzlsdf/analysis/B080/quintera/Projects/neuroblastoma/src/fig3c/atac.tstat_',cl1.name,'_',cl2.name,fName,'.pdf'), m, width = 7, height = 3.5)
# ggsave(paste0('/icgc/dkfzlsdf/analysis/B080/quintera/Projects/neuroblastoma/src/fig3c/atac.tstat_',cl1.name,'_',cl2.name,fName,'_flipped.pdf'), m.flipped, width = 7, height = 3.5)
# 
# dfSubset <- ddf[rownames(ddf) %in% TFList,]
# dim(dfSubset)
# m2 <- ggplot(dfSubset, aes(pos, val)) + geom_point() + 
#   geom_text_repel(aes(label=as.character(rownames(dfSubset)))) +
#   xlab("")
# 
# ggsave(paste0('/icgc/dkfzlsdf/analysis/B080/quintera/Projects/neuroblastoma/src/fig3c/atac.tstat_',cl1.name,'_',cl2.name,fName,'_CRC.pdf'), m2, width = 15, height = 8)


#================================================================================#
#                         plot for figure 3C                                     #
#================================================================================#

filter      = T
loadSEfile  = T
ncores      <- 1
Tumor       = F
quantileNormal = F


# flip around:
ddf.flipped = ddf

ddf.flipped$pos <- rev(ddf.flipped$pos)
ddf.flipped$val <- -ddf.flipped$val



library(RColorBrewer)
col = rev(colorRampPalette(brewer.pal(11,"PRGn"))(100)) # spectrum of green(nonMes) to purple(Mes)
ddf.flipped$name <- sub("_.*", "", ddf.flipped$name)
idx <- match(ddf.flipped$name, names(MES_activity))
table(is.na(idx))



ddf.flipped$mes  <- MES_activity[match( ddf.flipped$name, names(MES_activity))]





co.v <- cor(ddf.flipped$val, ddf.flipped$mes, method = "spearman", use = "pairwise.complete.obs") 

cor.df <- data.frame(x = 35, y = -400, text = paste0("Spearman = ", round(co.v, 2)))

m <- ggplot(ddf.flipped, aes(pos, val, color = mes)) + 
  geom_text_repel(aes(label=ddf.flipped$name), color = "black") +
  geom_point(size = 2) + 
  xlab("")+
  geom_text(data = cor.df, aes(x,y, label = text), color  = "black", size = 5)+ 
  ylab("t.statistic") +  expand_limits(x = 50, y = 500) + scale_colour_gradientn(colours = col)#+ scale_colour_brewer(palette = "PRGn", type = "seq")
m
#figure 4
#ggsave("/icgc/dkfzlsdf/analysis/B080/quintera/Projects/neuroblastoma/src/fig3c/SKNAS-CLBGA.pdf",m,width = 7, height = 3.5)





#CRCvector <- readLines("/icgc/dkfzlsdf/analysis/B080/quintera/Projects/neuroblastoma/src/fig3c/CRClist.txt")


dfSubset <- ddf.flipped[rownames(ddf.flipped) %in% CRCvector,]
dfSubset <- df[rownames(df) %in% CRCvector,]

dim(dfSubset)
summary(dfSubset)

dfSubset <- dfSubset[order(dfSubset$val, decreasing = TRUE),]
dfSubset <- dfSubset[order(dfSubset$val, decreasing = FALSE),]


co.v <- cor.test(dfSubset$val, dfSubset$mes, method = "spearman", use = "pairwise.complete.obs") 

m2 <- ggplot(dfSubset, aes(pos, val, color = mes, label = name)) + 
  #geom_text_repel(data = subset(dfSubset, mes>0), color = "black", nudge_y = 50) +
  #geom_text_repel(data = subset(dfSubset, mes<=0), color = "black", nudge_y = -50) +
  
  geom_text_repel(aes(label=dfSubset$name), color = "black") +
  #geom_text_repel(data = subset(dfSubset, mes>0), color = "black", nudge_y = 15 - subset(dfSubset, mes>0)$val) +
  #geom_text_repel(data = subset(dfSubset, mes>0), color = "black", nudge_y = c(seq(from=-100, by=5, length.out = sum(dfSubset$mes>0)) - subset(dfSubset, mes>0)$val) ) +
  #geom_text_repel(data = subset(dfSubset, mes>0), color = "black",  direction="x", nudge_y = c(seq(from=3, by=3, length.out = sum(dfSubset$mes>0))) ) +
  #geom_text_repel(data = subset(dfSubset, mes>0), color = "black", nudge_y = c(seq(from=5, by=4, length.out = sum(dfSubset$mes>0)) ) ) +
  #geom_text_repel(data = subset(dfSubset, mes<=0), color = "black", nudge_y = c(seq(from=-100, by=5, length.out = sum(dfSubset$mes<=0))) ) +
  #geom_text_repel(data = subset(dfSubset, mes>0), color = "black", nudge_y = c(15:(15 +sum(dfSubset$mes>0) -1)) - subset(dfSubset, mes>0)$val ) +
  
  #geom_text_repel(data = subset(dfSubset, mes<=0), color = "black", nudge_y = 30 + subset(dfSubset, mes<=0)$val) +
  geom_point(size = 2) + 
  xlab("")+
  #geom_text(data = cor.df, aes(x,y, label = text), color  = "black", size = 5)+ 
  ylab("t.statistic") +  
  #ylim(-250, 50) +
  #expand_limits(x = 50, y = 200) + 
  scale_colour_gradientn(colours = col)#+ scale_colour_brewer(palette = "PRGn", type = "seq")
m2
ggsave("results/figure3/figure3c_footprint.pdf", m2, width = 12, height = 8)





