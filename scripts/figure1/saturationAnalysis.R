# walltime=10:00:00,mem=1g

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

SEbed = as.character(args[1])
SEsaturation = as.character(args[2])
SEsaturationPlot = as.character(args[3])

#-------------------------------------------------------------------------------------------------------------------------
 # DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/"
 # SEbed = paste0(DATAPATH, "data/tumor/chipseq/H3K27ac/SE")
 # SEsaturation = paste0(DATAPATH, "analysis/tumor/SE_saturation/saturation.RDS")
 # SEsaturationPlot = paste0(DATAPATH, "results/figure1/SEsaturationAnalysis_coverage_with_extrapolation.pdf")
#---------------------------------------------------------------------------------------------------------------------------

# Load these packages
library(GenomicRanges)
library(rtracklayer)

# Paths to SE analyzed samples
goodsamples = list.files(SEbed, full.names = T)
sampID = gsub(SEbed, "", sapply(strsplit(goodsamples, "-", fixed=T), function(x)x[1]))
sampID = gsub("/", "", sampID)

# Saturation analysis of SE regions
allSElist = vector(mode="list", length=length(sampID))
names(allSElist) = sampID

for(i in 1: length(allSElist))
{
  allSElist[[i]] = import(con=goodsamples[i], format="BED")
}
rm(i)

perm_cnt = 1000
cmat = matrix(ncol=length(goodsamples), nrow=perm_cnt)

for(i in 1:length(goodsamples))
{
  for(p in 1:perm_cnt)
  {
    print(paste0("At file - ", i, " of ", length(goodsamples), " | Permutation step-", p))
    fileind = sample(sampID, i)
    tmp = allSElist[names(allSElist) %in% fileind]
    
    all = tmp[[1]]
    for(k in 1 : length(tmp))
    {
      if(k>1){all = suppressWarnings(c(all, tmp[[k]]))}
    }
    all = reduce(all)
    cmat[p,i] = (sum(width(all)))/10^6

    rm(fileind, tmp, all, k)
  }
}
rm(p, i, perm_cnt)
saveRDS(cmat, file = SEsaturation)

cmat <- readRDS(file = SEsaturation)

#Plotting
pdf(SEsaturationPlot, width=2.3, height=2.3)
par(mar=c(3,3.5,0.3,0.3), mgp=c(1.8,0.5,0),cex=1)

 df = data.frame(
                  upperQuart = apply(cmat,2,function(x){quantile(x,0.75)}),
                  median     = apply(cmat,2,function(x){quantile(x,0.50)}),
                  lowerQuart = apply(cmat,2,function(x){quantile(x,0.25)})
                )

 # Predicting the data
 x  = seq(1:nrow(df))
 x1 = c(nrow(df)+1):1000
 
 ymed  = df$median
 m.med = nls(ymed~a*x/(b+x), start=list(a=1, b=1))

 #Get some estimation of goodness of fit
 cor(ymed, predict(m.med))
 
 #Median only plot
 plot(x=x, ymed, col="#67000d", type="l", xlab="Samples", ylab=paste("Median coverage of","SE regions (Mb)",sep="\n"),
       pch=20, las=2, cex.axis=0.75, cex.lab=0.8, lwd=2, xlim=c(0,1000), ylim=c(0,210))

 # extrapolation
 a  = summary(m.med)$parameters[,1][1]
 b  = summary(m.med)$parameters[,1][2]
 y1 = (a*x1) / (b + x1)
 lines(x1, y1, lty=2, col="#fb6a4a", lwd=2)
 legend(x=330, y=200, legend=c("Observed", "Predicted"), lty=c(1,2), col = c("#67000d","#fb6a4a"), lwd=2, bty="n", cex=0.8, x.intersp=0.2,y.intersp=0.8)

 par("plt" = c(0.54, 0.95, 0.43, 0.77),mgp=c(1,0.5,0),cex=0.9)
 par(new = TRUE)

 #SE coverage
 plot(apply(cmat,2,median),type="l",col="#67000d", ylab="", xlab=paste0("Samples(n=",ncol(cmat),")"),las=2, cex.axis=0.75, cex.lab=0.8, lwd=2, ylim=c(0,130))
dev.off()

#################


source_data_saturation <- bind_rows(tibble(Sample = x, Median_coverage = ymed, Type = "Observed"),
          tibble(Sample = x1, Median_coverage = y1, Type = "Predicted"))

source_data_coverage <- tibble(Sample = 1:60, SE_coverage = apply(cmat,2,median))


write_xlsx(list(`Extended Data figure 1c satu` = source_data_saturation,
                `Extended Data figure 1c cov` = source_data_coverage), 
           path = "results/figure_source_data/Extended_Data_figure_1c.xlsx")


