library(rio)
library(ggplot2)
library(reshape2)
library(GenomicRanges)

path = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB_revision/"
#path = "/Users/ashwin/Documents/Projects/NB/00_superNB_revision/"

###################################
# Gartlgruber et.al patient samples
###################################

se1p = readRDS("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS")
#se1p = readRDS(paste0(path,"data/tumor_consensusSE_target_GRanges.RDS"))

###############################
# Gartlgruber et.al cell-lines
###############################

se1c = read.delim("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/analysis/cells/chipseq/H3K27ac/consensusSE/cells_H3K27ac_noH3K4me3_consensusSE.bed", header=F, stringsAsFactors = F)
#se1c = read.delim(paste0(path, "data/cells_H3K27ac_noH3K4me3_consensusSE.bed"), header=F, stringsAsFactors = F)
se1c = se1c[,1:3]
colnames(se1c) = c("chr", "start", "end")
se1c = makeGRangesFromDataFrame(se1c)

###################
# Groningen et.al 
###################

se2 = import(file = "https://media.nature.com/original/nature-assets/ng/journal/v49/n8/extref/ng.3899-S4.xlsx")
se2 = se2[,2:4]
se2 = se2[-1,]
colnames(se2) = se2[1,]
se2 = se2[-1,]
se2 = makeGRangesFromDataFrame(se2)

##################
# Boeva et.al
##################
se3 = import(file = "https://media.nature.com/original/nature-assets/ng/journal/v49/n9/extref/ng.3921-S3.xlsx")
se3 = se3[-1,]
colnames(se3) = se3[1,]
se3 = se3[-1,]
se3 = se3[, c(1:3,5,14,15,19)]

# Chosing a SE specific gene (single gene) based on highest correlation
se3$`Gene(s)` = sapply(strsplit(se3$`Gene(s)`, "|", fixed=T), function(x)x[1])
se3$`FC Score Group I over Group II` = round(log2(as.numeric(se3$`FC Score Group I over Group II`)),3)
se3$`Wilcoxon p-value (two sided test)` = signif(as.numeric(se3$`Wilcoxon p-value (two sided test)`),3)
se3 = se3[, -ncol(se3)]
tmp = makeGRangesFromDataFrame(se3[,1:3])
mcols(tmp) = se3[4:ncol(se3)]
se3 = tmp
rm(tmp)

################################
# Function to identify overlaps
################################

find_specific_fraction_overlap = function(refGR, testGR, fraction)
{
  hits <- findOverlaps(query = refGR, subject = testGR)
  overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
  
  percentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
  hits <- hits[percentOverlap > fraction]
  
  hits = round(c(length(unique(queryHits(hits)))/length(refGR)) * 100, 2)
  
  return(hits)
}

#######################################################################
# Overlap between SE identified from Gartlgruber et.al in NB patients
# and cell lines and SEs identified in Groningen et.al and Boeva et.al
# at 25% genomic overlap
#######################################################################

overlap_threshold = 0.25
res = data.frame(Class = c("Tumor_Gartlgruber_Groningen", "Tumor_Gartlgruber_Boeva", "Cell_Gartlgruber_Groningen", "Cell_Gartlgruber_Boeva"),
                 OverlapFraction = NA,
                 NonOverlapFraction = NA,
                 stringsAsFactors = F)

res$OverlapFraction[res$Class == "Tumor_Gartlgruber_Groningen"] = find_specific_fraction_overlap(refGR = se1p, testGR = se2, fraction = overlap_threshold)
res$OverlapFraction[res$Class == "Tumor_Gartlgruber_Boeva"]     = find_specific_fraction_overlap(refGR = se1p, testGR = se3, fraction = overlap_threshold)

res$OverlapFraction[res$Class == "Cell_Gartlgruber_Groningen"]  = find_specific_fraction_overlap(refGR = se1c, testGR = se2, fraction = overlap_threshold)
res$OverlapFraction[res$Class == "Cell_Gartlgruber_Boeva"]      = find_specific_fraction_overlap(refGR = se1c, testGR = se3, fraction = overlap_threshold)

res$NonOverlapFraction = 100 -  res$OverlapFraction

res = melt(res)

g <- ggplot(res, aes(fill=variable, y=value, x=Class)) + theme_bw(base_size = 9) +
labs(fill = "", x = "", y="Fraction") + coord_flip() +
geom_bar(position="fill", stat="identity") +
scale_fill_manual(values=c("#fc9272", "#a1d99b")) +
theme(panel.grid = element_blank(),
      axis.text = element_text(colour="black"), 
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave(paste0(path,"results/Overlap_SE_regions_different_studies_at_25_percent.pdf"), g, width=4, height=1)


library(writexl)
write_xlsx(list(`Extended Data figure 1e` = g$data), 
           path = "results/figure_source_data/Extended_Data_figure_1e.xlsx")


#######################################################################
# Overlap between SE identified from Gartlgruber et.al in NB patients
# and cell lines and SEs identified in Groningen et.al and Boeva et.al
# at 0 - 100% genomic overlap
#######################################################################

# overlap_threshold = seq(0, 1, 0.05)
# 
# ovl1 = rep(NA, length(overlap_threshold))
# ovl2 = rep(NA, length(overlap_threshold))
# ovl3 = rep(NA, length(overlap_threshold))
# ovl4 = rep(NA, length(overlap_threshold))
# 
# for(i in 1:length(overlap_threshold))
# {
#   ovl1[i] = find_specific_fraction_overlap(refGR = se1p, testGR = se2, fraction = overlap_threshold[i])
#   ovl2[i] = find_specific_fraction_overlap(refGR = se1p, testGR = se3, fraction = overlap_threshold[i])
#   
#   ovl3[i] = find_specific_fraction_overlap(refGR = se1c, testGR = se2, fraction = overlap_threshold[i])
#   ovl4[i] = find_specific_fraction_overlap(refGR = se1c, testGR = se3, fraction = overlap_threshold[i])
# }
# 
# res = cbind(ovl1, ovl2, ovl3, ovl4)
# rm(ovl1, ovl2, ovl3, ovl4, i)
# 
# pdf(paste0(path,"results/Overlap_SE_regions_different_studies.pdf"), width=3, height=3)
# 
#   par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2,0.5,0), cex=0.7, xaxt="n", las=2)
#  
#    plot(x = 1:21, y = res[,1], type="b", col="firebrick", ylim = c(0,100), pch=20,
#        xlab = "Minimum genomic coverage overlap threshold (in %)", 
#        ylab = "Super enhancers from Gartlgruber et.al found (in %)")
#   lines(x = 1:21, y = res[,2], type="b", col="forestgreen", pch=20)
#   
#   lines(x = 1:21, y = res[,3], type="b", col="firebrick", pch= 4)
#   lines(x = 1:21, y = res[,4], type="b", col="forestgreen", pch = 4)
#   
#   legend(x = 6, y = 100,
#          legend = c("Groningen et.al(n=1662)", 
#                     "Boeva et.al(n=5975)"), 
#          fill = c("firebrick", "forestgreen"),
#          bty="n", x.intersp = 0.5, cex=0.8, border = NA)
#   
#   legend(x = 6, y = 90,
#          legend = c("Gartlgruber et.al Tumors(n = 1973)",
#                     "Gartlgruber et.al Cell-lines(n = 2511)"), 
#          pch = c(20, 4),
#          bty="n", x.intersp = 0.5, cex=0.8, border = NA)
#   
#   mtext(text = paste0(overlap_threshold * 100, "%"), side = 1, line = 0.2, at = 1:21, cex = 0.7, las=2)
#   abline(v = 6, lwd=0.7, lty=2)
# dev.off()

