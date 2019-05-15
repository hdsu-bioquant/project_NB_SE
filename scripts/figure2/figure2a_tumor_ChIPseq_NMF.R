##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                  SE signal H matrix heatmap                                ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# Parse Cmd line.
options(echo=TRUE) # if you want to see commands in output file
args <- commandArgs(TRUE)


library(Bratwurst)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(knitr)


params <- list(K         = as.numeric(args[1]),
               annot     = args[2],
               norm_nmfW = args[3],
               norm_nmfH = args[4],
               hmatrix_wnorm = args[5],
               wmatrix_wnorm = args[6],
               nmf_features  = args[7],
               hmatrix_hnorm = args[8],
               figure2a      = args[9])
params
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          Read NMF objects                                  ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# read W normalized NMF
norm.nmf.exp <- readRDS(params$norm_nmfW)
# read H normalized NMF
norm.nmf.exp.H <- readRDS(params$norm_nmfH)

# read annotation 
annot <- readRDS(params$annot)
annot <- annot[match(colnames(norm.nmf.exp), annot$ProjectID),]
rownames(annot) <- annot$ProjectID


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                Signature annotation enrichment                             ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

auc <- function(rnk.list,max=NULL) {
  aux = sapply(rnk.list,function(rnk) {
    if (is.null(max)) {max = max(rnk)} 
    rnk = sort(rnk)
    X = 0
    i = 1
    ngenes = length(rnk)
    while ((rnk[i] <= max) && (i <= length(rnk))) {X = X + max -rnk[i];i = i+1}
    rauc = X/(i-1)/max
    return(rauc)
  })
  return(aux)
}

recovery_annot <- function(h, annot, annotID){
  annot.factor <- annot[,annotID]
  n.samples = nrow(annot)
  
  ALL.RNKS = lapply(levels(annot.factor),function(l) {
    RNKS=lapply(1:nrow(h),function(i) {
      exp = sort(h[i,],decreasing=TRUE)
      i.rnk = match(rownames(annot)[annot.factor==l],names(exp))
      i.rnk = sort(i.rnk[!is.na(i.rnk)])
      return(i.rnk)
    })
    names(RNKS) = paste0('Sig ',1:length(RNKS))
    return(RNKS)
  })
  names(ALL.RNKS) = levels(annot.factor)
  
  AUC.RAND = lapply(ALL.RNKS,function(r) {
    do.call('rbind',lapply(r, function(x) {
      ##
      l = lapply(1:500,function(i) {
        sample(1:n.samples,length(x))
      })
      aux = auc(l,max=n.samples)
      return(c(mean(aux),sd(aux)))
    }))
  })
  
  AUC = lapply(ALL.RNKS,auc,max=n.samples)
  
  
  PVAL = lapply(1:length(AUC),function(i) {
    x = data.frame(AUC.RAND[[i]],AUC[[i]])
    colnames(x) = c('mean','sd','val')
    z = (x[,3]-x[,1])/x[,2]
    p = ifelse(z>0,pnorm(z,lower.tail=FALSE),pnorm(z))
    x$z = z
    x$p = p
    return(x)
  })
  names(PVAL) = names(AUC)
  return(PVAL)
  
}
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                           H matrix                                         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
hmatrix <- HMatrix(norm.nmf.exp, k = params$K)
colnames(hmatrix) <- colnames(norm.nmf.exp)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##          Assign Signature Names according to annotation                    ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
signatureIDs <- c("MYCN", "MES", "MNA-LR", "MNA-HR")
names(signatureIDs) <- signatureIDs
signature_idx <- lapply(signatureIDs,function(x) 0)
signature_idx

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                        Find MYCN Signature                                 ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
rec_MYCN <- recovery_annot(hmatrix, annot, "MYCN")
signature_idx$MYCN <- which(rec_MYCN$Amp$p < 0.05 & rec_MYCN$Amp$z > 0)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                 Find High and Low  Risk Signatures                         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
rec_risk <- recovery_annot(hmatrix, annot, "Stage")
signature_idx$`MNA-HR` <- which(rec_risk$`4`$p < 0.05 & rec_risk$`4`$z > 0)
signature_idx$`MNA-LR` <- which(rec_MYCN$Amp$p < 0.05 & rec_MYCN$Amp$z < 0)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                        Assign MES Signature                                ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
signature_idx <- unlist(signature_idx)
signature_idx[!signature_idx %in% 1:nrow(hmatrix)] <- setdiff(1:nrow(hmatrix), signature_idx)
signature_idx

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##               Select opt K and compute feature stats                     ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
K <- params$K
norm.nmf.exp <- setOptK(norm.nmf.exp, K)
norm.nmf.exp@OptK

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                             NMH H Matrix                                 ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
norm.nmf.hmatrix <- HMatrix(norm.nmf.exp, k = K)[signature_idx,]
colnames(norm.nmf.hmatrix) <- colnames(norm.nmf.exp)
rownames(norm.nmf.hmatrix) <- names(signature_idx)
saveRDS(norm.nmf.hmatrix, params$hmatrix_wnorm)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                             NMH W Matrix                                 ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
norm.nmf.wmatrix <- WMatrix(norm.nmf.exp, k = K)[,signature_idx]
colnames(norm.nmf.wmatrix) <- names(signature_idx)
rownames(norm.nmf.wmatrix) <- rownames(norm.nmf.exp)
saveRDS(norm.nmf.wmatrix, params$wmatrix_wnorm)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                      NMF Signature features                              ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
norm.nmf.exp <- computeFeatureStats(norm.nmf.exp)
norm.nmf.exp <- computeSignatureSpecificFeatures(norm.nmf.exp, 1)
SignaFeat <- do.call(cbind, lapply(1:norm.nmf.exp@OptK, function(i){
  as.numeric(substr(norm.nmf.exp@FeatureStats$cluster, i, i))
}))
SignaFeat <- as.data.frame(SignaFeat)[,signature_idx]
#rownames are features names
rownames(SignaFeat) <- rownames(norm.nmf.exp)
#colnames are signatures
colnames(SignaFeat) <- names(signature_idx)
saveRDS(SignaFeat, params$nmf_features)


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                             NMH H Matrix                                 ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
norm.nmf.hmatrix.hnorm <- HMatrix(norm.nmf.exp.H, k = K)[signature_idx,]
colnames(norm.nmf.hmatrix.hnorm) <- colnames(norm.nmf.exp.H)
rownames(norm.nmf.hmatrix.hnorm) <- names(signature_idx)
saveRDS(norm.nmf.hmatrix, params$hmatrix_hnorm)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                        H matrix heatmap annotation                         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
#Annotation for H matrix heatmap
type.colVector <- list(Type    = c("Tumor"="grey80", "Metastasis"="grey60", "NA"="white"),
                       MYCN    = c("Amp"="grey20", "NonAmp"="grey80", "NA"="white"),
                       Stage   = c("4"="grey20", "1-3;4S"="grey80", "NA"="white"),
                       Age     = c(">18mo"="grey20", "<18mo"="grey80", "NA"="white"),
                       Relapse = c("Yes"="grey20", "No"="grey80", "NA"="grey100"))

#Heatmap annotation
heat.anno <- HeatmapAnnotation(df = data.frame(Type    = annot$Type,
                                               MYCN    = annot$MYCN,
                                               Stage   = annot$Stage,
                                               Age     = annot$Age,
                                               Relapse = annot$Relapse),
                               col = type.colVector,
                               show_annotation_name = TRUE, na_col = "white")


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##              Generate H matrix heatmap, W normalized                       ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
#plot H matrix
h.heatmap <- Heatmap(norm.nmf.hmatrix,
                     #col  = colorRamp2(seq(min(tmp.hmatrix), max(tmp.hmatrix), length.out = 100), viridis(n=100)),
                     col  = viridis(n=100),
                     name = "Exposure",
                     clustering_distance_columns = 'pearson',
                     show_column_dend = FALSE,
                     heatmap_legend_param = 
                       list(color_bar = "continuous", legend_height=unit(2, "cm")),
                     top_annotation = heat.anno,
                     show_column_names = FALSE,
                     show_row_names = FALSE,
                     cluster_rows = FALSE)
#print(h.heatmap)

#save H matrix heatmap
pdf(file = params$figure2a, width=10, height=7)
#print(h.heatmap)
#draw(h.heatmap)
h.heatmap
dev.off()