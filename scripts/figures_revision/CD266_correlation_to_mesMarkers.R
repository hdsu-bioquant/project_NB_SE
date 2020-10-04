#-------------------------------------------------------------------------------------------------
DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/"
outpath  = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB_revision/"
NBexprs  = paste0(DATAPATH, 'analysis/tumor/rnaseq/exprs/tumor_RNAseq_TPM_Matrix_filt_log.RDS')
tumorNMF = paste0(DATAPATH, 'analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS')
#-------------------------------------------------------------------------------------------------

#Primary tumors RNAseq NMF derived exposures values per signature (k=4) and expression data
NBtumors.sig = readRDS(tumorNMF)
mesSig = NBtumors.sig["MES",]
mesSig = c(mesSig - median(mesSig))/mad(mesSig)

# Expression of TNFRSF12A
NBtumors.exp = readRDS(NBexprs)
symbols = sapply(strsplit(rownames(NBtumors.exp), "|", fixed=T), function(x)x[2])
identical(colnames(NBtumors.exp), colnames(NBtumors.sig))

# Expression of CD266 and other well known Mes marker genes
cd266 = NBtumors.exp[symbols == "TNFRSF12A",]
vim   = NBtumors.exp[symbols == "VIM"]
fn1   = NBtumors.exp[symbols == "FN1"]
snail = NBtumors.exp[symbols == "SNAI1"]
#slug  = NBtumors.exp[symbols == "SNAI2"] # This is also highly correlated

# Emt marker genes from DOI: 10.1158/1078-0432.CCR-15-0876
# A Patient-Derived, Pan-Cancer EMT Signature Identifies Global Molecular Alterations and Immune Target Enrichment Following Epithelial-to-Mesenchymal Transition | Clinical Cancer Research | Feb 2016
emtSig = read.table(paste0(outpath, "data/emt_marker_genes.txt"), skip = 2, header=T, stringsAsFactors = F)
emtSig = split(emtSig$Genes, emtSig$Class)

# Subset NB expression data with Mesenchymal signature genes
NBemt.M = NBtumors.exp[symbols %in% emtSig$M,]

# Correlation of Mes signature genes with Mes exposures
cor_emtM_mesSig = apply(NBemt.M, 1, function(x) round(cor(mesSig, x, method="spearman"),2))

# Correlation of Mes signature genes with CD266 expression
cor_emtM_cd266 = apply(NBemt.M, 1, function(x) round(cor(cd266, x, method="spearman"),2))

# Correlation of 1000 random genes sets with Mes exposures
cor_rand_mesSig = t(sapply(seq(1:1000), function(x){
    rand_vals = apply(NBtumors.exp[sample(1:nrow(NBtumors.exp), length(emtSig$M)),], 1, function(y){
                round(cor(mesSig, y, method = "spearman"),2)})
}))

# Correlation of 1000 random genes sets with CD266 expression
cor_rand_cd266 = t(sapply(seq(1:1000), function(x){
  rand_vals = apply(NBtumors.exp[sample(1:nrow(NBtumors.exp), length(emtSig$M)),], 1, function(y){
    round(cor(cd266, y, method = "spearman"),2)})
}))

pdf(paste0(outpath, "results/CD266_correlation_to_MesMarkers.pdf"), width = 7, height = 4.5)
par(mfrow = c(2,3), mar = c(4,4,1,0.5), mgp = c(3,0.5,0), cex = 0.6, las = 2)

  plot(density(cor_emtM_mesSig), xlim=c(-1,1), ylim = c(0, 3), main = "", xlab= paste("Expression correlation between", "Mes. markers(Mak et.al) and Mes. exposure", sep="\n"), lwd=0)
  apply(cor_rand_mesSig, 2, function(x) lines(density(x),  col="grey90", lwd= 0.1))
  lines(density(cor_emtM_mesSig),  col="firebrick"); abline(v=0)
  legend("topleft", legend = paste("1000 random", "genesets", sep="\n"), lty=1, col = "grey90", bty="n", x.intersp = 0.5)
  
  plot(density(cor_emtM_cd266), xlim=c(-1,1), ylim = c(0, 3), main = "", xlab=paste("Expression correlation between", "Mes. markers(Mak et.al) and TNFRSF12A", sep="\n"), lwd=0)
  apply(cor_rand_cd266, 2, function(x) lines(density(x),  col="grey90", lwd= 0.1))
  lines(density(cor_emtM_cd266), col="firebrick"); abline(v=0)
  legend("topleft", legend = paste("1000 random", "genesets", sep="\n"), lty=1, col = "grey90", bty="n", x.intersp = 0.5)
  
  par(mgp = c(1.5,0.5,0))
  
  plot(cd266, mesSig, col="grey30", pch=20, xlab = "TNFRSF12A expression (log2 TPM)", ylab = "Mesenchymal score")
  legend("bottomright", 
         legend = c(paste0("rho=", round(cor(cd266, mesSig, method="spearman"),2)),
                    paste0("p=", signif(cor.test(cd266, mesSig, method="spearman")$p.value,2))),
         bty="n",y.intersp = 1, cex=1)
  abline(lm(mesSig~cd266))
  
  plot(cd266, vim, col="grey30", pch=20, xlab = "TNFRSF12A expression (log2 TPM)", ylab = "VIM expression (log2 TPM)")
  legend("bottomright", 
         legend = c(paste0("rho=", round(cor(cd266, vim, method="spearman"),2)),
                    paste0("p<2.2e-16")),
         bty="n",y.intersp = 1, cex=1)
  abline(lm(vim~cd266))
  
  plot(cd266, fn1, col="grey30", pch=20, xlab = "TNFRSF12A expression (log2 TPM)", ylab = "FN1 expression (log2 TPM)")
  legend("bottomright", 
         legend = c(paste0("rho=", round(cor(cd266, fn1, method="spearman"),2)),
                    paste0("p<2.2e-16")),
         bty="n",y.intersp = 1, cex = 1)
  abline(lm(fn1~cd266))
  
  plot(cd266, snail, col="grey30", pch=20, xlab = "TNFRSF12A expression (log2 TPM)", ylab = "SNAI1 expression (log2 TPM)")
  legend("bottomright", 
         legend = c(paste0("rho=", round(cor(cd266, snail, method="spearman"),2)),
                    paste0("p<2.2e-16")),
         bty="n",y.intersp = 1, cex = 1)
  abline(lm(snail~cd266))
dev.off()

sd <- rbind(cor_rand_cd266, cor_emtM_cd266)
rownames(sd) <- c(paste0("Random", 1:1000), "cor_emtM_cd266")
dim(sd)
write_xlsx(list(`Extended Data figure 9j` = as.data.frame(sd) %>% rownames_to_column("Assay")), 
           path = "results/figure_source_data/Extended_Data_figure_9j.xlsx")

write_xlsx(list(`Extended Data figure 9g` = tibble(ProjectID = names(cd266), TNFRSF12A_expression = cd266, VIM_expression = vim)), 
           path = "results/figure_source_data/Extended_Data_figure_9g.xlsx")

write_xlsx(list(`Extended Data figure 9h` = tibble(ProjectID = names(cd266), TNFRSF12A_expression = cd266, FN1_expression = fn1)), 
           path = "results/figure_source_data/Extended_Data_figure_9h.xlsx")

write_xlsx(list(`Extended Data figure 9i` = tibble(ProjectID = names(cd266), TNFRSF12A_expression = cd266, Mesenchymal_score = mesSig)), 
           path = "results/figure_source_data/Extended_Data_figure_9i.xlsx")
