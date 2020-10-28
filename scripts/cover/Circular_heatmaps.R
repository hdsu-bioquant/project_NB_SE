library(circlize)
library(ComplexHeatmap)
library(viridis)

hmatrix_wnorm<- readRDS("analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K4_Hmatrix_wnorm.RDS")



Heatmap(hmatrix_wnorm,
        col  = viridis(n=100),
        name = "Exposure",
        clustering_distance_columns = 'pearson',
        show_column_dend = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        cluster_rows = FALSE)

col_fun1 <- colorRamp2(seq(min(hmatrix_wnorm), max(hmatrix_wnorm), length.out = 100), 
                       viridis(n=100))

circos.par


circos.clear()
circos.par$gap.after <- 0
circos.heatmap(t(hmatrix_wnorm),
               col = col_fun1,
               track.height = .8)



pdf(file = "results/figures_cover/tumor_SE_hmatrix.pdf", width=10, height=10)
circos.clear()
circos.par$gap.after <- 0
circos.heatmap(t(hmatrix_wnorm),
               col = col_fun1,
               track.height = .8)

dev.off()



dir.create("results/figures_cover")
