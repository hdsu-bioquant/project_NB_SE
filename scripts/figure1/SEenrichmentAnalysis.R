options(echo=TRUE) # if you want to see commands in output file
args <- commandArgs(TRUE)

library(enrichR)
library(ggplot2)
library(ggrepel)
library(ggbeeswarm)

#path = as.character(args[1])

path = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"

# Super enhancers target genes
se = readRDS(paste0(path, "analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS"))
se = as.character(unique(se$target_SYMBOL))

# Enrichment analysis
enr =  enrichr(se, "GO_Biological_Process_2018")
GOenrich = enr$GO_Biological_Process_2018[,c(1:4,9)]
GOenrich = GOenrich[GOenrich$Adjusted.P.value < 0.05,]
GOenrich$P.value = signif(GOenrich$Adjusted.P.value, digits = 4)
GOenrich$Adjusted.P.value = signif(GOenrich$Adjusted.P.value, digits = 4)
rm(enr,se)

# Going through the enriched geneset lists and manually assigning metageneset tags
tmp = rep(NA,nrow(GOenrich))

tmp[grep("signal", GOenrich$Term)]="Transcriptional regulation and signalling"
tmp[grep("transcription", GOenrich$Term)]="Transcriptional regulation and signalling"
tmp[grep("response", GOenrich$Term)]="Transcriptional regulation and signalling"
tmp[grep("insulin", GOenrich$Term)]="Transcriptional regulation and signalling"
tmp[grep("calcium", GOenrich$Term)]="Transcriptional regulation and signalling"
tmp[grep("expression", GOenrich$Term)]="Transcriptional regulation and signalling" 
tmp[grep("ERK1", GOenrich$Term)]="Transcriptional regulation and signalling" 
tmp[grep("MAPK", GOenrich$Term)]="Transcriptional regulation and signalling" 
tmp[grep("kinase", GOenrich$Term)]="Transcriptional regulation and signalling"
tmp[grep("phosphorylation", GOenrich$Term)]="Transcriptional regulation and signalling"

tmp[grep("mesenchymal", GOenrich$Term)]="Cell migration and EMT"
tmp[grep("migration", GOenrich$Term)]="Cell migration and EMT"
tmp[grep("actin", GOenrich$Term)]="Cell migration and EMT"
tmp[grep("adhesion", GOenrich$Term)]="Cell migration and EMT"
tmp[grep("motility", GOenrich$Term)]="Cell migration and EMT"

tmp[grep("macromolecule", GOenrich$Term)]="Metabolism"
tmp[grep("compound", GOenrich$Term)]="Metabolism"

tmp[grep("dendrite", GOenrich$Term)]="Neuronal developmental"
tmp[grep("axon", GOenrich$Term)]="Neuronal developmental"
tmp[grep("neuron", GOenrich$Term)]="Neuronal developmental"
tmp[grep("nervous", GOenrich$Term)]="Neuronal developmental"

tmp[is.na(tmp)] = "Developmental processes"
GOenrich = cbind(GOenrich, Class=tmp); rm(tmp)

# Arranging the enrichment results in a clean format
GOenrich = GOenrich[order(GOenrich$Class),]
GOenrich = GOenrich[,c(1:4,6,5)]
GOenrich = split(GOenrich,GOenrich$Class)

GOenrich = lapply(GOenrich, function(x){
  x = x[order(x$Adjusted.P.value),]
  x = cbind(x, Rank = 1:nrow(x))
}) # ordering by fdr

avgFDR = sapply(GOenrich, function(x) median(x$Adjusted.P.value))
GOenrich = GOenrich[order(avgFDR)] # ordering by average fdr of each broad group
avgFDR = names(sort(avgFDR))
GOenrich = do.call("rbind", GOenrich)
GOenrich$Class = factor(as.character(GOenrich$Class), levels = avgFDR)
rownames(GOenrich) = c(1:nrow(GOenrich))
rm(avgFDR)

GOenrich$Term = sapply(strsplit(GOenrich$Term, " (", fixed=TRUE), function(x)x[1])
GOenrich$Value = round(-log10(GOenrich$Adjusted.P.value), 4)

# Writing the results
write.table(GOenrich[,1:6], paste0(path,"results/suppltables/GO_BP_enrichment_SE_target_genes.txt"), row.names=F, quote=F, sep="\t")

# Plotting the enrichment analysis results
toShow = rep("No",nrow(GOenrich))
toShow[GOenrich$Term %in% c("positive regulation of epithelial cell migration",
                       "positive regulation of epithelial to mesenchymal transition",
                       "nervous system development",
                       "noradrenergic neuron differentiation",
                       "positive regulation of signal transduction",
                       "positive regulation of transcription, DNA-templated",
                       "positive regulation of cell differentiation",
                       "positive regulation of macromolecule metabolic process"
                       )] = "Yes"

p = ggplot(GOenrich, aes(x = Class, y = Rank, size = Value, color = Class)) + theme_void(base_size = 9) + scale_y_reverse() +
    labs(size = expression(paste(-log[10],"(FDR)")), colour = "Meta-genesets") +
    geom_quasirandom(shape=21) + 
    geom_text_repel(data=subset(GOenrich, toShow == "Yes"), nudge_y = 0.5, aes(label=Term), size=2.5) + 
    guides(colour = guide_legend(nrow = 3, ncol=2, byrow = F, title.position = "top"), 
           size = guide_legend(nrow = 3, ncol=1, byrow = F, title.position = "top")) +
    theme(legend.position = "bottom") + scale_colour_brewer(palette = "Set1") 

ggsave(filename = paste0(path,"results/figure1/GO_BP_enrichment_SE_target_genes.pdf"), 
       plot = p, width=4.2, height=3.5)
