#------------------------------------------------------------------------------#
#                                 Read Annotation                              #
#------------------------------------------------------------------------------#

annot <- readRDS("annotation/annotation_tumor_phenoData_NMFassign.RDS")
head(annot)


#------------------------------------------------------------------------------#
#                  Add sig assignment based on top 10%                         #
#------------------------------------------------------------------------------#
top_10perc_assing <- function(wmatrix){
  sig_assign <- lapply(setNames(colnames(wmatrix), colnames(wmatrix)), function(sigID){
    selec_wmatrix <- do.call(cbind, lapply(as.data.frame(wmatrix), function(sign_expo){
      sign_expo[sign_expo < quantile(sign_expo, 0.90)] <- NA
      sign_expo
    }))
    rownames(selec_wmatrix) <- rownames(wmatrix)
    selec_wmatrix <- selec_wmatrix[!is.na(selec_wmatrix[,sigID]),,drop=FALSE]
    # Keep only the top feature if there's an overlap
    sig_SE_IDs <- rownames(selec_wmatrix[rowMaxs(selec_wmatrix, na.rm = TRUE) == selec_wmatrix[,sigID],])
    #print(sig_SE_IDs)
    setNames(rep(sigID, length(sig_SE_IDs)), sig_SE_IDs)
    
  })
  #print(unlist(sig_assign))
  #do.call(c, sig_assign)
  unlist(unname(sig_assign))
}

# Add signature assignment 
tumor_hchip <- readRDS("analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS")
#top_10perc_assing(t(tumor_hchip))
tumor_hchip <- top_10perc_assing(t(tumor_hchip))

tumor_hrna <- readRDS("analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS")
tumor_hrna <- top_10perc_assing(t(tumor_hrna))


annot <- annot %>% 
  mutate(Signature_ChIPseq = tumor_hchip[match(ProjectID, names(tumor_hchip))]) %>% 
  mutate(Signature_RNAseq = tumor_hrna[match(ProjectID, names(tumor_hrna))])

length(tumor_hchip)
table(is.na(annot$Signature_ChIPseq))
length(tumor_hrna)
table(is.na(annot$Signature_RNAseq))

saveRDS(annot, "annotation/annotation_tumor_phenoData_NMFassign_top10perc.RDS")





#------------------------------------------------------------------------------#
#                  Add sig assignment based on cluster                         #
#------------------------------------------------------------------------------#

# Run k means over all rows and assign features to the clusters
SignFeatures <- function(wmatrix){
  ssf <- apply(wmatrix, 1, function(x){
    #print(x)
    x <- sigmoidTransform(x)
    k <- kmeans(x, 2)
    max_idx <- which.max(k$centers)
    #paste(if_else(k$cluster == max_idx, "1", "0"), collapse = "")
    setNames(if_else(k$cluster == max_idx, 1, 0), names(x))
  })
  ssf
  
  ssf <- ssf[,colSums(ssf) == 1]
  setNames(rownames(ssf)[apply(ssf, 2, which.max)], colnames(ssf))
}
#SignFeatures(t(tumor_hchip))

# Add signature assignment 
tumor_hchip <- readRDS("analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS")
tumor_hchip <- SignFeatures(t(tumor_hchip))

tumor_hrna <- readRDS("analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS")
tumor_hrna <- SignFeatures(t(tumor_hrna))


annot <- annot %>% 
  mutate(Signature_ChIPseq = tumor_hchip[match(ProjectID, names(tumor_hchip))]) %>% 
  mutate(Signature_RNAseq = tumor_hrna[match(ProjectID, names(tumor_hrna))])

length(tumor_hchip)
table(is.na(annot$Signature_ChIPseq))
length(tumor_hrna)
table(is.na(annot$Signature_RNAseq))

saveRDS(annot, "annotation/annotation_tumor_phenoData_NMFassign_kmeansk2hmat.RDS")

