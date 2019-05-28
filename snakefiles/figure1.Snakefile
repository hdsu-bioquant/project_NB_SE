#==============================================================================#
#            Analysis and figures included in manuscript figure 1              #
#==============================================================================#


rule compile_figure1:
    input:
        #figure1a = join(DATAPATH, 'results/figure1/figure1a_tumor_SE_hmatrix.pdf'),
        #figure1b = join(DATAPATH, 'results/figure1/figure1b_cells_SE_hmatrix.pdf'),
        figure1c = join(DATAPATH, 'results/figure1/figure1c_HockeyStick_plot.pdf'),
        figure1d = join(DATAPATH, 'results/figure1/GO_BP_enrichment_SE_target_genes.pdf'),
        figure1e = join(DATAPATH, 'results/figure1/figure1e_IGV_plot.pdf')
    output: join(DATAPATH, 'results/figure1/figure1_paths.txt')
    shell:
        """
        touch {output}
        #echo 'Figure 1a {input.figure1e}' >> {output}
        #echo 'Figure 1b {input.figure1e}' >> {output}
        #echo 'Figure 1c {input.figure1c}' >> {output}
        #echo 'Figure 1d {input.figure1d}' >> {output}
        echo 'Figure 1e {input.figure1e}' >> {output}
        
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
        table  = join(DATAPATH, 'results/suppltables/GO_BP_enrichment_SE_target_genes.txt'),
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

