library(ggplot2)
library(cowplot)
library(reshape2)
library(ggrepel)
library(gridExtra)

#-------------------------------------------------------------------------------
DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/"
outpath  = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB_revision/"
NBexprs  = paste0(DATAPATH, 'analysis/tumor/rnaseq/exprs/tumor_RNAseq_TPM_Matrix_filt_log.RDS')
tumorNMF = paste0(DATAPATH, 'analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS')
#-------------------------------------------------------------------------------

#Staining samples #NB-S-325, I032_004, and I040_001
patID = setNames(c("NSP050-PT01", "NSP254-RM01", NA), c("NB-S-325", "I032_004", "I040_001"))

#Primary tumors RNAseq NMF derived exposures values per signature (k=4)
NBtumors.sig = readRDS(tumorNMF)
NBtumors.sig = apply(NBtumors.sig, 2, function(x)rownames(NBtumors.sig)[which.max(x)])

# Expression of TNFRSF12A
NBtumors.exp = readRDS(NBexprs)
symbols = sapply(strsplit(rownames(NBtumors.exp), "|", fixed=T), function(x)x[2])

# Expression of  marker genes RNAseq cohort
#stain = t(NBtumors.exp[symbols %in% c("NCAM1", "VIM", "ALK", "FLI1", "SOX6", "CCND1", "ERG", "PHOX2B", "MAML3"),])
stain = t(NBtumors.exp[symbols %in% c("NCAM1", "VIM", "PHOX2B"),])

if(identical(names(NBtumors.sig), rownames(stain))){
  colnames(stain) = sapply(strsplit(colnames(stain), "|", fixed=T), function(x)x[2])
  stain = data.frame(SampleID = rownames(stain), stain, Signature = NBtumors.sig)
  stain = melt(stain)
  stain$label = "No"
  stain$label[stain$SampleID %in% patID] = "Yes"
  colnames(stain) = c("SampleID", "Signature", "Marker_gene", "Expression", "Label")
  stain$Signature = factor(as.character(stain$Signature), 
                           levels = c("MES", "MYCN", "MNA-HR", "MNA-LR"))
}

p1 = ggplot(stain, aes(x=Signature, y=Expression)) + theme_bw(base_size = 9) + labs(x="") +
    geom_violin(trim=FALSE, aes(fill = Signature), color="grey95", lwd=0.2, draw_quantiles = 0.5) + #geom_boxplot(width=0.1, outlier.shape = NA) +
    facet_wrap(~ Marker_gene) +
    scale_fill_manual(values=c("#397FB9", "#E21E25", "#994F9F", "#50AF49")) +
    #geom_point(data = subset(stain, Label == "Yes")) +
    #geom_text_repel(data = subset(stain, Label == "Yes"), aes(label = SampleID), size=3) +
    theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right")
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2 = ggplot(stain, aes(x=Marker_gene, y=Expression)) + theme_bw(base_size = 9) + labs(x="") +
  geom_violin(trim=FALSE, aes(fill = Marker_gene), color="grey95", lwd=0.2, draw_quantiles = 0.5) + #geom_boxplot(width=0.1, outlier.shape = NA) +
  facet_wrap(~ Signature) +
  #geom_point(data = subset(stain, Label == "Yes")) +
  #geom_text_repel(data = subset(stain, Label == "Yes"), aes(label = SampleID), size=3) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right")

p = plot_grid(p1, p2,  nrow = 2, rel_heights  = c(1/2, 1))

save_plot(paste0(outpath,"results/staining_marker_gene_expression_RNAseq_cohort.pdf"), 
          plot = p, base_width = 3.5, base_height = 4)


write_xlsx(list(`Figure 2h` = p1$data), 
           path = "results/figure_source_data/Figure_2h.xlsx")
