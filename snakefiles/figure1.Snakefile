#==============================================================================#
#            Analysis and figures included in manuscript figure 1              #
#==============================================================================#


rule compile_figure1:
    input:
        figure1a = join(DATAPATH, 'results/figure1/cohortDespFig.pdf'),
        figure1b = join(DATAPATH, 'results/figure1/SE_multiple_tissue_overlap_to_NB_SE_MainFig.pdf'),
        figure1c = join(DATAPATH, 'results/figure1/figure1c_HockeyStick_plot.pdf'),
        figure1d = join(DATAPATH, 'results/figure1/GO_BP_enrichment_SE_target_genes.pdf'),
        figure1e = join(DATAPATH, 'results/figure1/figure1e_IGV_plot.pdf'),
        figure1f = join(DATAPATH, 'results/figure1/SEsaturationAnalysis_coverage_with_extrapolation.pdf')
    output: join(DATAPATH, 'results/figure1/figure1_paths.txt')
    shell:
        """
        touch {output}
        echo 'Figure 1a {input.figure1a}' >> {output}
        echo 'Figure 1b {input.figure1b}' >> {output}
        echo 'Figure 1c {input.figure1c}' >> {output}
        echo 'Figure 1d {input.figure1d}' >> {output}
        echo 'Figure 1e {input.figure1e}' >> {output}
        echo 'Figure 1f {input.figure1f}' >> {output}
        
        """

#================================================================================#
# Neuroblastoma study cohort description and Chipseq QC metrics                  #
#================================================================================#

rule fig1_cohortDespChipQC:
    input:
        tumorQCvals = join(DATAPATH, 'data/tumor/chipseq/H3K27ac/quality_control_metrics_tumors.RDS'),
        cellQCvals  = join(DATAPATH, 'data/cells/chipseq/H3K27ac/quality_control_metrics_cells.RDS'),
        tumorNMF    = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS'),
        cellNMF     = join(DATAPATH, 'analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K3_Hmatrix_hnorm.RDS')
    output:
        figMain = join(DATAPATH, 'results/figure1/cohortDespFig.pdf'),
        figSupp = join(DATAPATH, 'results/sup_figure1/CHIPqc.pdf')
    params:
        script  = 'scripts/figure1/chipQCandCohortDesp.R',
        outpath = join(DATAPATH, 'results/')
    conda: '../envs/R3.5.yaml'
    shell:
        """
        Rscript {params.script} {input.tumorQCvals} {input.cellQCvals} {input.tumorNMF} {input.cellNMF} {params.outpath}
        """

#================================================================================#
# Neuroblastoma super enhancer target gene expression across multiple tissues    #
# Neuroblastoma super enhancer overlap with multiple other tissues               #
#================================================================================#

rule fig1_NBSEcomparison:
    input:
        expdata  = join(DATAPATH, 'db/TCGA_TARGET_GTex/TcgaTargetGtex_log2_fpkm.RDS'),
        sampInfo = join(DATAPATH, 'db/TCGA_TARGET_GTex/TcgaTargetGtex_sample_information.RDS'),
        SEnbs    = join(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS'),
        SEothers = join(DATAPATH, 'db/SEmultiTisuues/')
    output:
        figMain  = join(DATAPATH, 'results/figure1/SE_multiple_tissue_overlap_to_NB_SE_MainFig.pdf'),
        figSupp1 = join(DATAPATH, 'results/sup_figure1/TCGA_TARGET_GTex_NB_SE_target_genes_expression_comparision.pdf'),
        figSupp2 = join(DATAPATH, 'results/sup_figure1/SE_multiple_tissue_overlap_to_NB_SE_SupplFig.pdf')
    params:
        script  = 'scripts/figure1/NBpancanOtherTissuesSEcomparision.R',
        outpath = join(DATAPATH, 'results/')
    conda: '../envs/R3.5.yaml'

    shell:
        """
          Rscript {params.script} {input.expdata} {input.sampInfo} {input.SEnbs} {input.SEothers} {params.outpath}
        """




#================================================================================#
# Neuroblastoma Super enhancers saturation analysis                              #
#================================================================================#

rule fig1_NBsaturationAnalysis:
    input:
        SEbedPath = join(DATAPATH, 'data/tumor/chipseq/H3K27ac/SE')
    output:
        SEsaturation = join(DATAPATH, 'analysis/tumor/SE_saturation/saturation.RDS'),
        figMain      = join(DATAPATH, 'results/figure1/SEsaturationAnalysis_coverage_with_extrapolation.pdf')
    params:
        script  = 'scripts/figure1/saturationAnalysis.R'
    conda: '../envs/R3.5.yaml'
    shell:
        """
          Rscript {params.script} {input.SEbedPath} {output.SEsaturation} {output.figMain}
        """



#================================================================================#
#                      Figure 1e - MAML3 loci                                    #
#================================================================================#
rule fig1e_IGV:
    input:
        consensusSE = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed')
    output:
        report = join(DATAPATH, 'reports/figure1e_IGV_plot.html'),
        rmd    = temp(join(DATAPATH, 'reports/figure1e_IGV_plot.Rmd')),
        figure = join(DATAPATH, 'results/figure1/figure1e_IGV_plot.pdf')
    params:
        script       = 'scripts/figure1/figure1e_IGV_plot.Rmd',
        work_dir     = DATAPATH,
        path_config  = join(os.getcwd(), 'scripts/figure1/figure1e_tracks.txt'),
        window   = config['igv_plot']['figure1e']['window'],
        ymax     = config['igv_plot']['figure1e']['ymax'],
        gr_chr   = config['igv_plot']['figure1e']['chr'],
        gr_start = config['igv_plot']['figure1e']['start'],
        gr_end   = config['igv_plot']['figure1e']['end'],
        gr_name  = config['igv_plot']['figure1e']['name']
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  work_dir     = '{params.work_dir}', \
                  path_config  = '{params.path_config}', \
                  width_window = {params.window}, \
                  ymax  = '{params.ymax}', \
                  chr   = '{params.gr_chr}', \
                  start = {params.gr_start}, \
                  end   = {params.gr_end}, \
                  name  = '{params.gr_name}', \
                  figure = '{output.figure}' \
                ))"


        """

#================================================================================#
#                      Figure 1 - Enrichment analysis                            #
#================================================================================#

rule fig1_SEenrichmentAnalysis:
    input:
        SE_target = join(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS')
    output:
        table  = join(DATAPATH, 'results/supptables/GO_BP_enrichment_SE_target_genes.txt'),
        figure = join(DATAPATH, 'results/figure1/GO_BP_enrichment_SE_target_genes.pdf')
    params:
        script   = 'scripts/figure1/SEenrichmentAnalysis.R',
        work_dir = DATAPATH
    conda: '../envs/R3.5.yaml'
    shell:
        """
        Rscript {params.script} {params.work_dir} {input.SE_target}
        """



#================================================================================#
#                      Figure 1 - Hockey Stick plot                              #
#================================================================================#
rule fig1_Hockey:
    input:
        tumor_annot = join(DATAPATH, 'annotation/annotation_tumor.RDS'),
        SE_target   = join(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS'),
        enhancers   = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusEnhancers/tumor_H3K27ac_noH3K4me3_consensusEnhancers.bed'),
        bw = expand(join(DATAPATH, 'data/tumor/chipseq/H3K27ac/bw/{sample}_H3K27ac.bw'), zip, sample=TUMOR_SAMPLES_CHIP)
    output:
        report = join(DATAPATH, 'reports/figure1c_HockeyStick_plot.html'),
        rmd    = temp(join(DATAPATH, 'reports/figure1c_HockeyStick_plot.Rmd')),
        figure = join(DATAPATH, 'results/figure1/figure1c_HockeyStick_plot.pdf')
    params:
        script       = 'scripts/figure1/figure1c_HockeyStick_plot.Rmd',
        work_dir     = DATAPATH
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  work_dir  = '{params.work_dir}', \
                  annot     = '{input.tumor_annot}', \
                  SE_target = '{input.SE_target}', \
                  enhancers = '{input.enhancers}', \
                  figure    = '{output.figure}' \
                ))"


        """

