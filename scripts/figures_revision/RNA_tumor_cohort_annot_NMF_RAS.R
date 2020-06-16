rna_annotation <- readRDS("annotation/annotation_tumor_phenoData_NMFassign.RDS")
colnames(rna_annotation)
table(rna_annotation$Risk)
rna_annotation <- rna_annotation %>% 
  # Drop blacklisted
  filter(avail.RNAseq) %>% 
  # select final columns
  select(ProjectID, Signature_RNAseq,  Type, MYCN, Stage.class, Stage, Age.days, Age, Relapse, 
         TimeToRelapse, PrimaryID, Sex, Location, Purity, Risk, Status, 
         TERT, ALT, FOXR1, FOXR2, CCircle, GenomicSubtype, 
         chr1p, chr3p, chr11q, chr17q, avail.ChIPseq, avail.RNAseq, 
         geo_accession.x, a_efs_all, b_os_all, c_sex_all, d_fav_all,
         e_efs_hr, f_os_hr, efs_bin, efs_day, os_bin, os_day, 
         death_from_disease, progression, Ploidy) 

hmatrix <- readRDS("analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_Hmatrix_wnorm.RDS")
rownames(hmatrix) <- paste0("NMF_", rownames(hmatrix))

# table(ras %in% rownames(tumor_exprs))
# table(junfos %in% rownames(tumor_exprs))


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            RAS targets Loboda                              ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## RAS metagene signature
ras = readRDS("db/publicGeneSigs/ras_target_genes.RDS")
ras = ras$upreg
ras

tumor_exprs <- readRDS("analysis/tumor/rnaseq/exprs/tumor_RNAseq_TPM_Matrix_symbol.RDS")
#colnames(tumor_exprs) <- toupper(colnames(tumor_exprs))
table(colnames(hmatrix) %in% colnames(tumor_exprs) )
tumor_exprs <- tumor_exprs[,colnames(tumor_exprs) %in% colnames(hmatrix) ]
table(ras %in% rownames(tumor_exprs) )

tumor_exprs_ras <- tumor_exprs[rownames(tumor_exprs) %in% ras,]
dim(tumor_exprs_ras)

ras_loboda <- as.data.frame(tumor_exprs_ras) %>% 
  rownames_to_column("geneID") %>% 
  pivot_longer(cols = -geneID, names_to = "ProjectID", values_to = "Expression") %>% 
  group_by(ProjectID) %>% 
  summarise(RAS_Loboda = median(Expression))

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            RAS targets Way                                 ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## RAS metagene signature
ras = readRDS("db/publicGeneSigs/Way_RAS_sig.RDS")

tumor_exprs <- readRDS("analysis/tumor/rnaseq/exprs/tumor_RNAseq_TPM_Matrix_symbol.RDS")
tumor_exprs <- tumor_exprs[,colnames(tumor_exprs) %in% colnames(hmatrix) ]
table(ras %in% rownames(tumor_exprs) )

tumor_exprs_ras <- tumor_exprs[rownames(tumor_exprs) %in% ras,]
dim(tumor_exprs_ras)

ras_way <- as.data.frame(tumor_exprs_ras) %>% 
  rownames_to_column("geneID") %>% 
  pivot_longer(cols = -geneID, names_to = "ProjectID", values_to = "Expression") %>% 
  group_by(ProjectID) %>% 
  summarise(RAS_Way = median(Expression))







annot_NMF_RAS <- as.data.frame(t(hmatrix))
annot_NMF_RAS <- annot_NMF_RAS %>% 
  #rename_all(list(~ str_replace(., "Sig", "NMF_Sig"))) %>% 
  rownames_to_column("ProjectID") %>% 
  #mutate(Signature_assign = sig_assign[match(Sample, names(sig_assign))]) %>% 
  left_join(ras_loboda, by = "ProjectID") %>% 
  left_join(ras_way, by = "ProjectID") %>% 
  left_join(rna_annotation, by = "ProjectID")




saveRDS(annot_NMF_RAS, "annotation/In-house_annot_NMF_RAS.RDS")
library(writexl)
write_xlsx(list(`TARGET NB - NMF and RAS score` = annot_NMF_RAS), 
           path = "annotation/In-house_annot_NMF_RAS.xlsx")











