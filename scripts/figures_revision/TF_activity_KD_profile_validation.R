library(ggplot2)
library(reshape2)
library(ggpubr)
library(parallel)
library(ggrepel)

DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/"
KDdata   = paste0(DATAPATH, 'analysis/cells/crcGIEMSAkd/nbKDinhouse.RDS')
TFact    = paste0(DATAPATH, 'analysis/tumor/VIPER/')
outpath  = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB_revision/"

## CRC TFs list
crcTF = read.table(paste0(DATAPATH, "results/supptables/crcTF_modules.txt"), sep="\t", header=T, stringsAsFactors = F)

# CRCs selected for KD analysis
usedForKD = c("NFKB2", "ETV6", "FOSL2", "SMAD3",
              "RUNX2", "ETS1","RUNX1", "RUNX3",
              "NR3C1", "RARB", "MYC", "SOX13",
              "SOX9", "GLI2", "SOX6", "NFIB", 
              "MYCN", "EBF1", "ZBTB7C", "GATA2",
              "ZNF423", "TBX2","POU2F2","FOXO3", "HMX1",
              "SOX11","SREBF2", "PBX1", "HAND2",
              "KLF7", "GATA3", "PHOX2B")

used_class = rep("No", nrow(crcTF))
used_class[crcTF$TF %in% usedForKD] = "Yes"

crcTF$used_for_KD_experiments = used_class
rm(used_class, usedForKD)

#----------------------------------------------------------------
# SANITY CHECKS 
#---------------------------------------------------------------

# Is there any bias in distribution of TFs selected for KD experiments among the groups ADRN and MES?
# The answer is NO

table(crcTF$used_for_KD_experiments, crcTF$Module)
chisq.test(table(crcTF$used_for_KD_experiments, crcTF$Module))

#      ADRN MES
# No    22  11
# Yes   19  13
# 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  table(crcTF$used_for_KD_experiments, crcTF$Module)
# X-squared = 0.12387, df = 1, p-value = 0.7249

# What mean counts of ADRN and MES TF would you expect if 50% of all CRCs were sampled randomly 1000 times
# Compare this expectation to the number of ARDN and MES selected for the KD experiments
# Is there a bias ? Did we inadvertantly chose more ADRN or MES TFs compared to random ?
# The answer is NO

cnts = mclapply(1:1000, function(i, drop_frac=0.5){
  set.seed(i)
  random = sample(crcTF$Module, round(drop_frac*nrow(crcTF)))
  random = c(table(random))
  return(random)
}, mc.cores = 60)

cnts = do.call("rbind", cnts)

# Expected mean count of ADRN and MES TFs
round(colMeans(cnts))

# ADRN  MES 
# 20    12 

# These numbers above can nicely be calculated using a hyper geometric distribution

#ADRN counts
#qhyper(m = 41, n = 24, k = 32, p = 0.5)
#[1] 20

# MES counts
#qhyper(m = 24, n = 41, k = 32, p = 0.5)
#[1] 12

# Actual counts of ADRN and MES TFs selected 
table(crcTF$Module[crcTF$used_for_KD_experiments == "Yes"])

# ADRN  MES 
# 19   13 


# cnts = melt(cnts)[,2:3]
# ggdensity(cnts, x = "value",
#           add = "mean", 
#           color = "Var2", palette = c("#00AFBB", "#E7B800")) + geom_vline(xintercept = 19, color = "red")

#---------------------------------------------------------------

## Our knockdown analysis data (internal validation)
kd.NB.fc = readRDS(KDdata)
kd.NB.fc = data.frame(TF = rownames(kd.NB.fc), log2FC_GIMEN_vs_NMB = kd.NB.fc$log2FC, stringsAsFactors = F)

## TF activity file
mesTFact = readRDS(paste0(TFact, "MES_TFactivity.RDS"))
mesTFact = data.frame(TF = names(mesTFact$es$nes), MES_TF_activity = mesTFact$es$nes, stringsAsFactors = F)

# Combine all data
crcTF = merge(x = crcTF, y = kd.NB.fc, all.x = T)
crcTF = merge(x = crcTF, y = mesTFact, all.x = T)

rm(kd.NB.fc, mesTFact)

# Check if the correlation matches to the one in the figure 
corr = cor.test(crcTF$MES_TF_activity, crcTF$log2FC_GIMEN_vs_NMB, method="spearman", use="pairwise.complete.obs")

ggplot(crcTF, aes(x=log2FC_GIMEN_vs_NMB, y=MES_TF_activity)) + theme_bw(base_size=9) + 
  labs(x="Knockdown sensitivity log2(GIMEN/NMB)", y="TF activity") + geom_point(size=1.5) +
  # geom_point(aes(color=log2FC),size=1.5) + scale_color_gradient(high="yellow", low="red") +
  geom_smooth(method = "lm", se = FALSE, colour="grey30") +
  geom_text_repel(aes(label=TF),size=2.5) + 
  annotate("text", x = -0.43, y = -3.3, label = paste(paste0("p=",round(corr$p.value,3)),paste0("rho=",round(corr$estimate,2)),sep="\n"), size=2.8)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Downsampling data and computing correlations to check effect of outliers
df = crcTF[, 4:5]
df = df[!is.na(df$log2FC_GIMEN_vs_NMB),]

rho_random = mclapply(1:1000, function(i, drop_frac=0.15){
  set.seed(i)
  
  #remove 25% of the data
  random = df[- sample(1:nrow(df), round(drop_frac*nrow(df))),]
  cor.subset = cor(random$MES_TF_activity, random$log2FC_GIMEN_vs_NMB, method="spearman", use="pairwise.complete.obs")
  return(cor.subset)
}, mc.cores = 60)

rho_random = unlist(rho_random)

p1 = ggdensity(rho_random, ylab = "Density", add = "mean", fill="grey90",
               xlab = paste("1000 correlations between TF activity & KD sensitivity", 
                            "with 15% (n=5) of data removed randomly for each correlation", sep="\n")) +
    geom_vline(xintercept = corr$estimate , color="red")

# TF activity between CRCs selected and non selected for the knockdown experiments
tmp = melt(crcTF[,c(3,5)])

p2 <- ggboxplot(tmp, x = "used_for_KD_experiments", y = "value", 
               xlab = "", ylab="MES TF activity",
               color = "used_for_KD_experiments", palette = "jco",
               add = "jitter") +
  stat_compare_means(label = "p.format")
p2 = p2 + stat_compare_means(label = "p.format")

pdf(paste0(outpath, "results/TF_activity_KD_profile_validation.pdf"), width=8, height=3)
ggarrange(p1,p2)
dev.off()



write_xlsx(list(`Extended Data figure 7e` = crcTF[,c(1,2,3,5)]), 
           path = "results/figure_source_data/Extended_Data_figure_7e.xlsx")
