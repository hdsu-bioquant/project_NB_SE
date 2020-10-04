options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

tqcpath = as.character(args[1])
cqcpath = as.character(args[2])
tNMFval = as.character(args[3])
cNMFval = as.character(args[4])
outpath = as.character(args[5])

#-------------------------------------------------------------------------------------------------------------------------
# DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/"
# tqcpath = paste0(DATAPATH, "data/tumor/chipseq/H3K27ac/quality_control_metrics_tumors.RDS")
# cqcpath = paste0(DATAPATH, "data/cells/chipseq/H3K27ac/quality_control_metrics_cells.RDS")
# tNMFval = paste0(DATAPATH, "analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS")
# cNMFval = paste0(DATAPATH, "analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K3_Hmatrix_hnorm.RDS")
# outpath = paste0(DATAPATH, "results/")
#---------------------------------------------------------------------------------------------------------------------------

library(corrplot)
library(RColorBrewer)

## QC values tumors and cells
tqc = readRDS(tqcpath)
cqc = readRDS(cqcpath)

## Signature exposures tumors and cells
sigT = t(readRDS(tNMFval))
sigT = sigT[match(rownames(tqc),rownames(sigT)),]

sigC = t(readRDS(cNMFval))
sigC = sigC[match(rownames(cqc),rownames(sigC)),]

# Combined annotation of tumors and cells
tmp = matrix(ncol=3, nrow=nrow(cqc))
colnames(tmp) = colnames(tqc)[2:4]
tmp = cbind(cqc[,1, drop=F], tmp)
anno = rbind(cbind(Type="Tumor", tqc[,1:4]), 
             cbind(Type="Cell", tmp))

rowids = rownames(anno)
anno = as.data.frame(apply(anno,2,as.character),stringsAsFactors=F)

anno$Type[anno$Type == "Tumor"]=1
anno$Type[anno$Type == "Cell"]=2

anno$MYCN[anno$MYCN == "Amp"]=3
anno$MYCN[anno$MYCN == "NonAmp"]=4

anno$Stage[anno$Stage == "4"]=5
anno$Stage[anno$Stage == "1-3;4S"]=6

anno$Age[anno$Age == ">18mo"]=7
anno$Age[anno$Age == "<18mo"]=8

anno$Relapse[anno$Relapse == "Yes"]=9
anno$Relapse[anno$Relapse == "No"]=10

anno = apply(anno,2,as.numeric)
anno = anno[,rev(1:ncol(anno))]
rownames(anno) = rowids

rm(tmp, rowids)

## Combine QC and exposures
tqc = cbind(tqc,sigT)
cqc = cbind(cqc, sigC)
qc  = rbind(tqc[,5:14], cqc[,2:11])

## Correlation among QC parameters and exposures
tcor = round(cor(tqc[,5:ncol(tqc)], method="spearman"),2)
ccor = round(cor(cqc[,2:ncol(cqc)], method="spearman"),2)

tcorsig = cor.mtest(tqc[,5:ncol(tqc)], method="spearman")$p
ccorsig = cor.mtest(cqc[,2:ncol(cqc)], method="spearman")$p

## Clean variables
rm(tqc, sigT, cqc, sigC)

#-----------------------
# Final plotting
#-----------------------

# Integrative plot for cohort description and QC metrics
pdf(paste0(outpath,"sup_figure1/CHIPqc.pdf"), width=13,height=8)
layout(cbind(1:7,c(8,8,8,9,9,9,10)),  height=c(rep(0.6,6),0.8), width=c(1,0.25))

  par(mar=c(0.5,4,2.2,0.5), mgp=c(2,0.5,0), cex=0.7, xaxs="i", yaxs="i", xaxt="n")

  #SE frequency
  a = t(as.matrix(qc[,1]))
  x = barplot(a, las=2, ylab=paste("SE", "frequency",sep="\n"), border=NA, ylim=c(0,1000), col="#bfd3e6")
  abline(v=x[60]+0.75,lty=2)

  #CRC frequency
  b = t(as.matrix(qc[,2]))+1
  x = barplot(b, log="y", las=2, ylab=paste("CRC count", "(log scale)",sep="\n"), border=NA, col="#ccebc5", ylim=c(1,200))
  abline(v=x[60]+0.75,lty=2)

  #Library size
  c = t(as.matrix(qc[,3]))/10^6
  x = barplot(c, las=2, ylab=paste("Library size", "(in millions)",sep="\n"), border=NA, ylim=c(0,50), col="#f1b6da")
  abline(v=x[60]+0.75,lty=2)

  # Broad and Narrow peaks
  d = t(as.matrix(qc[,c(4:5)]))/1000
  x = barplot(d, beside=T, las=2, legend.text=T, args.legend=list(x=57,y=115, horiz=T, bty="n", border=NA, x.intersp=0.3, y.intersp=0.8), ylab=paste("# Peaks","(x1000)",sep="\n"), border=NA, ylim=c(0,100),col=c("#fb6a4a", "#fcbba1"))
  abline(v=x[,60][2]+0.75,lty=2)

  # FRiP scores
  e = t(as.matrix(qc[,c(6:7)]))
  x = barplot(e, beside=T, las=2, legend.text=T, args.legend=list(x=60,y=0.9, horiz=T, bty="n", border=NA, x.intersp=0.3, y.intersp=0.8), ylab=paste("FRiP", "scores",sep="\n"), border=NA, ylim=c(0,0.8),col=c('#fdbf6f','#ff7f00'))
  abline(v=x[,60][2]+0.75,lty=2)

  # NSC, RSC and PBC quality metrics
  f = t(as.matrix(qc[,c(8:10)]))
  x=barplot(f, beside=T, las=2, legend.text=T, args.legend=list(x=59.5,y=2.2, horiz=T, bty="n", border=NA, x.intersp=0.3, y.intersp=0.8), ylab=paste("QC", "values",sep="\n"), border=NA, ylim=c(0,3), col=c("#c7e9c0","#74c476","#238b45"), xaxt="n")
  abline(h=c(1.05,0.8,0.5), col=c("#c7e9c0","#74c476","#238b45"))
  abline(v=x[,60][3]+0.75,lty=2)

  par(mar=c(5.8,4,0,0.5), mgp=c(0,0.2,0), xaxt="s")

  # Sample annotation
  col = rev(brewer.pal(10,"Paired"))
  image(anno, col=col, xaxt="n", yaxt="n")
  grid(col="grey10",nx=nrow(anno), ny=ncol(anno), lty=1,lwd=0.5);box()
  axis(side=2, at=seq(0,1,length.out=ncol(anno)), labels=colnames(anno), cex.axis=0.9, tick=F,las=2)
  axis(side=1, at=seq(0,1,length.out=nrow(anno)), labels=rownames(anno), col.axis="grey30", las=2, line=0, tick=F)

  par(mar=c(0,0,0,0), mgp=c(0,0,0), cex=0.7)

  # QC correlation plots
  corrplot(tcor, type = "upper", tl.col = "black", tl.cex=0.8, diag=F, p.mat=tcorsig, sig.level = 10^-3, insig="pch", pch=4, pch.cex=0.8, title = "Tumor samples")
  corrplot(ccor, type = "upper", tl.col = "black", tl.cex=0.8,  diag=F, p.mat=ccorsig, sig.level = 10^-3, insig="pch", pch=4, pch.cex=0.8, title = "Cell samples")

  # Annotation legend
  par(mar=c(0,0,0,0))
  plot.new(); legend("top", legend=c("Tumor","Cell", "MYCN Amp","MYCN NonAmp", "Stage-4","Stage-1-3,4s", "Age>18mo","Age<18mo", "Relapse","Non Relapse"), fill=col,  border=NA, bty="n", x.intersp=0.3, y.intersp=0.9, cex=0.9)

  #rm(a,b,c,d,e,f,x)
dev.off()


source_data <- tibble(ProjectID = rownames(anno),
                      SE_frequency = a,
                      CRC_frequency = b,
                      Library_size = c,
                      Broad_peaks = d[1,],
                      Narrow_peaks = d[2,],
                      FRiP_score_MACS2 = e[1,],
                      FRiP_score_SICER = e[2,],
                      NSC_quality_metric = f[1,],
                      RSC_quality_metric = f[2,],
                      PBC_quality_metric = f[3,])

write_xlsx(list(`Extended Data figure 1a` = source_data), 
           path = "results/figure_source_data/Extended_Data_figure_1a.xlsx")

write_xlsx(list(`Extended Data figure 1b Tumor` = as.data.frame(tcor) %>% rownames_to_column("Metric"),
                `Extended Data figure 1b Cell` = as.data.frame(ccor) %>% rownames_to_column("Metric")), 
           path = "results/figure_source_data/Extended_Data_figure_1b.xlsx")





#Only cohort description
pdf(paste0(outpath,"figure1/cohortDespFig.pdf"),width=4.8,height=1.8)
  layout(matrix(c(rep(1,5),2:6),byrow=T,ncol=5), width=rep(1,5), height=c(1,0.8))
  par(mar=c(3,3,0.25,0.25),mgp=c(0,0.2,0))
  
  anno = anno[anno[,"Type"] == 1,] ## Selecting only tumors
  anno = anno[,-5] ## Deleting the Type (tumor or cell) column
  anno = anno[order(anno[,4],anno[,3],anno[,2],anno[,1], decreasing=T),]
  
  col = rev(brewer.pal(8,"Paired"))
  image(anno, col=col, xaxt="n", yaxt="n")
  grid(col="white",nx=nrow(anno), ny=ncol(anno), lty=1,lwd=0.8);# box()
  axis(side=2, at=seq(0,1,length.out=ncol(anno)), labels=colnames(anno), cex.axis=0.7, tick=F,las=2)
  axis(side=1, at=seq(0,1,length.out=nrow(anno)), labels=rownames(anno), col.axis="grey30", las=2, line=0, tick=F, cex.axis=0.7)
  
  par(mar=c(0.1,3.5,0,0), xpd=T)
  plot.new(); legend("center", legend=c("Amp","NonAmp"), title="MYCN", title.adj=0,  fill=col[1:2], border=NA, bty="n", x.intersp=0.3, y.intersp=0.8, cex=0.7)
  plot.new(); legend("center", legend=c("4","1-3,4s"), title="Stage", title.adj=0, fill=col[3:4], border=NA, bty="n", x.intersp=0.3, y.intersp=0.8, cex=0.7)
  plot.new(); legend("center", legend=c(">18mo","<18mo"), title="Age", title.adj=0, fill=col[5:6], border=NA, bty="n", x.intersp=0.3, y.intersp=0.8, cex=0.7)
  plot.new(); legend("center", legend=c("Yes","No"), title="Relapse", title.adj=0, fill=col[7:8], border=NA, bty="n", x.intersp=0.3, y.intersp=0.8, cex=0.7)
dev.off()

