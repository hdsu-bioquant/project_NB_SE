setwd("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/")
setwd("/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/superNB/src/project_NB_SE/scripts")

window <- 1000

rmarkdown::render("src/project_NB_SE/scripts/figures_revision/MNA_LR_HR_IGV_plot.Rmd", 
                  params = list(width_window = window,
                                chr = "chr11",
                                start = 2310000,
                                end   = 2450000,
                                name  = "TSPAN32",
                                figure= "results/figures_revision/figure_MNA-HR_MNA-LR_IGV_TSPAN32.pdf",
                                width = 11,
                                height= 6))


rmarkdown::render("src/project_NB_SE/scripts/figures_revision/MNA_LR_HR_IGV_plot.Rmd", 
                  params = list(width_window = window,
                                chr = "chr6",
                                start = 145487585,
                                end   = 146127445,
                                name  = "EMP2A",
                                figure= "results/figures_revision/figure_MNA-HR_MNA-LR_IGV_EMP2A.pdf",
                                width = 11,
                                height= 6))


rmarkdown::render("src/project_NB_SE/scripts/figures_revision/MNA_LR_HR_IGV_plot.Rmd", 
                  params = list(width_window = window,
                                chr = "chr17",
                                start = 48260810,
                                end   = 48345708,
                                name  = "COL1A1",
                                figure= "results/figures_revision/figure_MNA-HR_MNA-LR_IGV_COL1A1.pdf",
                                width = 11,
                                height= 6))


rmarkdown::render("src/project_NB_SE/scripts/figures_revision/MNA_LR_HR_IGV_plot.Rmd", 
                  params = list(width_window = window,
                                chr = "chr5",
                                start = 51506328,
                                end   = 52263001,
                                name  = "ITGA1|PELO",
                                figure= "results/figures_revision/figure_MNA-HR_MNA-LR_IGV_ITGA1_PELO.pdf",
                                width = 11,
                                height= 6))




summbw_type <- "sumbw_bar"
nbars <- 150

rmarkdown::render("src/project_NB_SE/scripts/figures_revision/MNA_LR_HR_IGV_plot.Rmd", 
                  params = list(width_window = ceiling(abs((2310000 - 2450000)/nbars)),
                                chr = "chr11",
                                start = 2310000,
                                end   = 2450000,
                                name  = "TSPAN32",
                                figure= "results/figures_revision/figure_MNA-HR_MNA-LR_IGV_TSPAN32_bar.pdf",
                                width = 11,
                                height= 6,
                                summbw = summbw_type))


rmarkdown::render("src/project_NB_SE/scripts/figures_revision/MNA_LR_HR_IGV_plot.Rmd", 
                  params = list(width_window = ceiling(abs((145487585 - 146127445)/nbars)),
                                chr = "chr6",
                                start = 145487585,
                                end   = 146127445,
                                name  = "EMP2A",
                                figure= "results/figures_revision/figure_MNA-HR_MNA-LR_IGV_EMP2A_bar.pdf",
                                width = 11,
                                height= 6,
                                summbw = summbw_type))


rmarkdown::render("src/project_NB_SE/scripts/figures_revision/MNA_LR_HR_IGV_plot.Rmd", 
                  params = list(width_window = ceiling(abs((48260810 - 48345708)/nbars)),
                                chr = "chr17",
                                start = 48260810,
                                end   = 48345708,
                                name  = "COL1A1",
                                figure= "results/figures_revision/figure_MNA-HR_MNA-LR_IGV_COL1A1_bar.pdf",
                                width = 11,
                                height= 6,
                                summbw = summbw_type))


rmarkdown::render("src/project_NB_SE/scripts/figures_revision/MNA_LR_HR_IGV_plot.Rmd", 
                  params = list(width_window = ceiling(abs((51506328 - 52263001)/nbars)),
                                chr = "chr5",
                                start = 51506328,
                                end   = 52263001,
                                name  = "ITGA1|PELO",
                                figure= "results/figures_revision/figure_MNA-HR_MNA-LR_IGV_ITGA1_PELO_bar.pdf",
                                width = 11,
                                height= 6,
                                summbw = summbw_type))




