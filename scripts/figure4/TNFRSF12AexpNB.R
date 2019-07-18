options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

NBcells  = as.character(args[1])
NBanno   = as.character(args[2])
NBexprs  = as.character(args[3])
tumorNMF = as.character(args[4])
outpath  = as.character(args[5])

#-------------------------------------------------------------------------------------------------
# DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"
# 
# NBcells  = paste0(DATAPATH, 'analysis/cells/rnaseq/exprs/cells_RNAseq_TPM_Matrix_filt_log.RDS')
# NBanno   = paste0(DATAPATH, 'annotation/annotation_cells.RDS')
# NBexprs  = paste0(DATAPATH, 'analysis/tumor/rnaseq/exprs/tumor_RNAseq_TPM_Matrix_filt_log.RDS')
# tumorNMF = paste0(DATAPATH, 'analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS')
# outpath  = paste0(DATAPATH, 'results/')
#-------------------------------------------------------------------------------------------------

# TNFRSF12A expression in MES vs ADRN in NB primary tumors and cells
# TNFRSF12A is a marker gene for NB MES cells/tumors

NBcells.anno = readRDS(NBanno)
NBcells.exp  = readRDS(NBcells)
NBcells.exp = NBcells.exp[grep("TNFRSF12A", rownames(NBcells.exp)),]
NBcells.exp = sort(unlist(NBcells.exp))

#Primary tumors RNAseq NMF derived exposures values per signature (k=4) and expression data
NBtumors.sig = readRDS(tumorNMF)
NBtumors.sig = apply(NBtumors.sig, 2, function(x) rownames(NBtumors.sig)[which.max(x)])
NBtumors.sig[NBtumors.sig != "MES"] = "ADRN"

NBtumors.exp = readRDS(NBexprs)
identical(colnames(NBtumors.exp), names(NBtumors.sig))
NBtumors.exp =  as.numeric(NBtumors.exp[grep("TNFRSF12A", rownames(NBtumors.exp)),])

pdf(paste0(outpath,"figure4/TNFRSF12A_expression_NB.pdf"), width=3.5, height=2.5)
layout(matrix(c(1,2), byrow = T, ncol=2), width=c(1,0.5))
par( mar=c(3.3, 3.3,0.5,0.5), mgp = c(2,0.5,0), cex=0.7, las=2, xaxs="i", yaxs="i")
  cutoff = quantile(NBcells.exp, 0.75)
  plot(NBcells.exp, ylab = "TNFRSF12A expression (log2 TPM)", col=ifelse(NBcells.exp > cutoff, "black", "grey80"),
       pch=20, xlim = c(0, 25), ylim = c(0,8))
  abline(h=cutoff)
  text(x = 3, y = cutoff, labels = "Top 20%", cex=1, pos = 3)
  
  xpos = which(NBcells.exp > cutoff)
  ypos = NBcells.exp[xpos]
  text(xpos, ypos, labels = names(NBcells.exp)[xpos], pos=2, cex=0.7)
  
  boxplot(NBtumors.exp ~ NBtumors.sig, pch=20,  ylab = "TNFRSF12A expression (log2 TPM)", ylim=c(-0.5,9),
          boxlty=1, whisklty=1, staplelty=1, medlwd=0.9, boxwex=0.7, boxlwd=0.9, whisklwd=0.9, staplelwd=0.9)
  text(x = 1.5, y = 8.5, labels = paste0("p=", signif(wilcox.test(NBtumors.exp ~ NBtumors.sig)$p.value, 3)), cex=0.7)
  text(x = c(1,2) , y = 0, labels = c(paste0("n=", table(NBtumors.sig)[1]),paste0("n=", table(NBtumors.sig)[2])), cex=0.7)
dev.off()

