library(glmnet)
library(survminer)
library(survival)
library(rio)
library(reshape2)

outpath = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB_revision/"
anno    = readRDS("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/annotation/annotation_tumor_phenoData_NMFassign.RDS")

#-------------------------------------------------------------------------------
#anno = readRDS("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/annotation/annotation_tumor_phenoData_NMFassign_top10perc.RDS")
#anno = readRDS("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/annotation/annotation_tumor_phenoData_NMFassign_kmeansk2hmat.RDS")
#-------------------------------------------------------------------------------

anno = anno[colnames(anno) %in% 
              c("ProjectID",
                #"Sex",
                "MYCN",
                "Stage", 
                "Age",
                'a_efs_all',
                "efs_day",
                "b_os_all",
                "os_day",
                "Signature_RNAseq",
                "avail.RNAseq")]

anno = anno[anno$avail.RNAseq,]
anno = anno[,c(1:4,6,8,7,9,10)]
#colnames(anno) = c("SampleID", "MYCN", "Stage", "Age", "Sex", "EFS_status", "EFS_days", "OS_status", "OS_days", "Signature")
colnames(anno) = c("SampleID", "MYCN", "Stage", "Age", "EFS_status", "EFS_days", "OS_status", "OS_days", "Signature")

#--------------------------------------------------
# anno$Signature[is.na(anno$Signature)]="Unassigned"
#--------------------------------------------------

anno$MYCN  = factor(anno$MYCN, levels = c("NonAmp", "Amp"))
anno$Stage = factor(anno$Stage, levels = c("1-3;4S", "4"))
anno$Age  = factor(anno$Age, levels = c("<18mo", ">18mo"))
#anno$Sex = factor(anno$Sex, levels = c("M", "F"))

anno$EFS_status = as.numeric(anno$EFS_status)
anno$EFS_days = as.numeric(anno$EFS_days)/365

anno$OS_status = as.numeric(anno$OS_status)
anno$OS_days = as.numeric(anno$OS_days)/365

anno$Signature = factor(anno$Signature, levels = c("MNA-LR", "MYCN", "MNA-HR", "MES"))

#-------------------------------------------------------------------------------
#anno$Signature = factor(anno$Signature, levels = c("MNA-LR", "MYCN", "MNA-HR", "MES", "Unassigned"))
#-------------------------------------------------------------------------------

anno_OS = droplevels(anno[!is.na(anno$OS_status),])
anno_EF = droplevels(anno[!is.na(anno$EFS_status),])

# Function to extract pvalues from posthoc test

posthocP = function(data) {
  data = reshape2::melt(data$p.value)
  data = data[!is.na(data$value),]
  data$value = signif(data$value,2)
  data = data[order(data$value),]
  id = paste0(paste(data$Var1, data$Var2, sep = " vs ")," = ", data$value)
  return(id)
}

#------------------------------------------------------------------------------#
# Survival analysis - Overall survival
#------------------------------------------------------------------------------#

fitOS = survfit(Surv(OS_days, OS_status) ~ Signature, data = anno_OS)
posthoc_OS = pairwise_survdiff(Surv(OS_days, OS_status) ~ Signature, data = anno_OS)

pOS = ggsurvplot( 
                  fit = fitOS,                       
                  data = anno_OS,              
                  pval.method = TRUE,
                  pval = F, 
                  conf.int = TRUE,
                  palette = c("#50AF49", "#E21E25", "#994F9F", "#397FB9"),
                  
                  xlim = c(0,18), 
                  break.time.by = 1, 
                 
                  ggtheme = theme_bw(base_size = 9),  
                  
                  risk.table = "absolute",
                  risk.table.y.text.col = T,
                  risk.table.y.text = F, 
                  
                  #surv.median.line = "hv",
                  #cumevents = T,
                  #cumcensor = T,
                  
                  fontsize = 3
              )

pOS$plot = pOS$plot + ggplot2::labs(x = "Time (in years)") +
          ggplot2::annotate("text", x = 0, 
                            y = seq(0.25, 0, length.out = 6), hjust = 0,
                            label = posthocP(posthoc_OS), size = 3)
  
pdf(paste0(outpath, "results/NB_signatures_OS.pdf"),width=6, height=6)
print(pOS, newpage=F)
dev.off()

pOS
write_xlsx(list(`Figure 2f` = pOS$data.survplot), 
           path = "results/figure_source_data/Figure_2f.xlsx")




# fit.coxph.OS <- coxph(Surv(OS_days, OS_status) ~ MYCN + Stage + Age + Sex + Signature, data = anno_OS)
# pOSX = ggforest(fit.coxph.OS, data = anno_OS, main = "Hazard ratio - Overall survival")
# ggsave(filename = paste0(outpath, "results/NB_signatures_OS_Cox.pdf"),
#        plot = pOSX, width=7, height=7, units ="in")

fit.coxph.OS <- coxph(Surv(OS_days, OS_status) ~ MYCN + Stage + Age + Signature, data = anno_OS)
pOSX = ggforest(fit.coxph.OS, data = anno_OS, main = "Hazard ratio - Overall survival")
pOSX
ggsave(filename = paste0(outpath, "results/NB_signatures_OS_Cox_without_sex.pdf"),
       plot = pOSX, width=7, height=7, units ="in")

fit.coxph.OS$score

sd <- summary(fit.coxph.OS)$coefficients %>% 
  as.data.frame() %>% 
  rownames_to_column("Levels")
write_xlsx(list(`Extended Data figure 5a` = sd), 
           path = "results/figure_source_data/Extended_Data_figure_5a.xlsx")

#------------------------------------------------------------------------------#
# Survival analysis - Event free survival
#------------------------------------------------------------------------------#

fitEFS = survfit(Surv(EFS_days, EFS_status) ~ Signature, data = anno_EF)
posthoc_EFS <- pairwise_survdiff(Surv(EFS_days, EFS_status) ~ Signature, data = anno_EF)

pEFS = ggsurvplot( 
        fitEFS,
        data = anno_EF,
        pval.method = TRUE,
        pval = F,
        conf.int = TRUE,
        palette = c("#50AF49", "#E21E25", "#994F9F", "#397FB9"),
        
        xlim = c(0,18),
        break.time.by = 1, 
        
        ggtheme = theme_bw(base_size = 9),
        
        risk.table = "absolute",
        risk.table.y.text.col = T,
        risk.table.y.text = F,
        
        #surv.median.line = "hv",
        #cumevents = T,
        #cumcensor = T,
        
        fontsize = 3
      )

pEFS$plot = pEFS$plot + ggplot2::labs(x = "Time (in years)") +
  ggplot2::annotate("text", x = 0, 
                    y = seq(0.25, 0, length.out = 6), hjust = 0,
                    label = posthocP(posthoc_EFS), size = 3)


pdf(paste0(outpath, "results/NB_signatures_EFS.pdf"),width=6, height=6)
print(pEFS, newpage=F)
dev.off()

fit.coxph.EFS <- coxph(Surv(EFS_days, EFS_status) ~ MYCN + Stage + Age + Sex + Signature, data = anno_EF)
pEFX = ggforest(fit.coxph.EFS, data = anno_EF, main = "Hazard ratio - Event free survival")
ggsave(filename = paste0(outpath, "results/NB_signatures_EFS_Cox.pdf"), 
       plot = pEFX, width=7, height=7, units ="in")

#------------------------------------------------------------------------------#
# glmnet - COx model
#------------------------------------------------------------------------------#

# cvfit = cv.glmnet(x = t(expDat.pri),
#                   y = Surv(time = cliDat.pri$OS, event =  cliDat.pri$OS.status), 
#                   family = "cox",  alpha=0, nlambda = 500)
# 
# pdf(paste0(path,"analysis/10foldCV_ridgeRegression.pdf"), width=1.5, height=1.5)
# par(mar=c(3,3,0.5,0.5), mgp=c(1.5,0.5,0), cex=0.6)
# plot(cvfit, col="grey80", xlim=c(-5,5), ylim=c(10,30), las=1.5, sign.lambda = -1)
# dev.off()
# 
# # The vertical lines represent these
# cvfit$lambda.min ; -log(cvfit$lambda.min)
# cvfit$lambda.1se ; -log(cvfit$lambda.1se)
# 
# coef.sel = coef(cvfit, s = cvfit$lambda.min)
# coef.sel = unlist(coef.sel[coef.sel[,1] != 0 ,])
# coef.sel = data.frame(do.call("rbind",strsplit(rownames(coef.sel),"|", fixed=T)),
#                       round(as.numeric(coef.sel),4), stringsAsFactors = F)
# colnames(coef.sel) = c("HGNC","Entrez","Coef")
# coef.sel = coef.sel[order(coef.sel$Coef, decreasing = T),]
# write.table(coef.sel, paste0(path, "analysis/ridgeRegression_coefficients.tsv"), 
#             sep="\t", quote=F, row.names = F)
