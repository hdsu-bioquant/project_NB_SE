#==============================================================================#
#            Analysis and figures included in manuscript figure 4              #
#==============================================================================#

rule compile_figure4:
    input:
        figure4a = join(DATAPATH,"results/figure4/junfos_corr_exposures.pdf"),
        figure4b = join(DATAPATH,"results/figure4/ras_corr_exposures.pdf"),
        #figure4c = join(DATAPATH, 'results/figure2/figure2c_tumor_cells_SE_UMAP.pdf'),
        #figure4d = join(DATAPATH, 'results/figure2/figure2d_01_tumor_riverplot.pdf'),
        #figure4e = join(DATAPATH, 'results/figure2/figure2e_tumor_SE_targets_hmatrix.pdf'),
        #figure4f = join(DATAPATH, 'results/figure2/figure2f_tumor_cells_density.pdf'),
        #figure4i = join(DATAPATH, 'results/figure4/figure4i_IGV_plot.pdf')
        figure4g = join(DATAPATH, 'results/figure4/tumors_RNAseq_PrimaryVsRelapse.pdf')
    output: join(DATAPATH, 'results/figure4/figure4_paths.txt')
    shell:
        """
        touch {output}
        #echo 'Figure 4a {input.figure4a}' >> {output}
        #echo 'Figure 4b {input.figure4b}' >> {output}
        #echo 'Figure 4c {input.figure4a}' >> {output}
        #echo 'Figure 4d {input.figure4a}' >> {output}
        #echo 'Figure 4e {input.figure4a}' >> {output}
        #echo 'Figure 4f {input.figure4a}' >> {output}
        #echo 'Figure 4g {input.figure4g}' >> {output}
        echo 'Figure 4i {input.figure4a}' >> {output}
        
        """

#================================================================================#
#   Correlation of RAS and JUN/FOS signature expression to Mesenchymal exposure  #
#================================================================================#
rule fig4_RAS_JUN_FOS:
    input:
        NBexprs = join(DATAPATH, 'analysis/tumor/rnaseq/exprs/tumor_RNAseq_TPM_Matrix_filt_log.RDS'),
        tumoNMF = join(DATAPATH, 'analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS'),
        rasSigr = join(DATAPATH, 'db/publicGeneSigs/ras_target_genes.RDS'),
        NBreg   = join(DATAPATH, 'analysis/tumor/ARACNe/network.txt'),
        NBmut   = join(DATAPATH, 'annotation/NB_mutation_matrix.RDS')
    output:
        figMain_1 = join(DATAPATH,"results/figure4/junfos_corr_exposures.pdf"),
        figMain_2 = join(DATAPATH,"results/figure4/ras_corr_exposures.pdf"),
        figSupp_1 = join(DATAPATH,"results/sup_figure4/ras_proteins_corr_exposures.pdf"),
        figSupp_2 = join(DATAPATH,"results/sup_figure4/AP1_complex_proteins_corr_exposures.pdf")
    params:
        script  = 'scripts/figure4/RasJunFosAnalysis.R',
        outpath = join(DATAPATH, 'results/')
    conda: '../envs/R3.5.yaml'
    shell:
        """
        Rscript {params.script} {input.NBexprs} {input.tumoNMF} {input.rasSigr} {input.NBreg} {input.NBmut} {params.outpath}
        """

#================================================================================#
#             Figure 4 - tumors RNAseq Primary vs. Relapse                       #
#================================================================================#
optK_tr = str(config['NMFparams']['tumor']['optimalK']['rnaseq'])
rule fig4_tumors_RNAseq_PrimaryVsRelapse:
    input:
        tumor_h = join(DATAPATH, ('analysis/tumor/rnaseq/NMF/tumor_consensusSE_K' + optK_tr + '_Hmatrix_wnorm.RDS'))
    output:
        report    = join(DATAPATH, 'reports/tumors_RNAseq_PrimaryVsRelapse.html'),
        rmd       = temp(join(DATAPATH, 'reports/tumors_RNAseq_PrimaryVsRelapse.Rmd')),
        figure    = join(DATAPATH, 'results/figure4/tumors_RNAseq_PrimaryVsRelapse.pdf')
    params:
        script   = 'scripts/figure4/tumors_RNAseq_PrimaryVsRelapse.Rmd',
        work_dir = DATAPATH
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  work_dir     = '{params.work_dir}', \
                  tumor_h      = '{input.tumor_h}', \
                  figure       = '{output.figure}' \
                ))"


        """


#================================================================================#
#                      Figure 4i - VSNL1 loci                                    #
#================================================================================#
rule fig4i_IGV:
    input:
        consensusSE = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed')
    output:
        report = join(DATAPATH, 'reports/figure4i_IGV_plot.html'),
        rmd    = temp(join(DATAPATH, 'reports/figure4i_IGV_plot.Rmd')),
        figure = join(DATAPATH, 'results/figure4/figure4i_IGV_plot.pdf')
    params:
        script       = 'scripts/figure4/figure4i_IGV_plot.Rmd',
        work_dir     = DATAPATH,
        path_config  = join(os.getcwd(), 'scripts/figure4/figure4i_tracks.txt'),
        width_window = config['igv_plot']['figure4i']['window'],
        ymax     = config['igv_plot']['figure4i']['ymax'],
        gr_chr   = config['igv_plot']['figure4i']['chr'],
        gr_start = config['igv_plot']['figure4i']['start'],
        gr_end   = config['igv_plot']['figure4i']['end'],
        gr_name  = config['igv_plot']['figure4i']['name']
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  work_dir     = '{params.work_dir}', \
                  path_config  = '{params.path_config}', \
                  width_window = {params.width_window}, \
                  ymax  = '{params.ymax}', \
                  chr   = '{params.gr_chr}', \
                  start = {params.gr_start}, \
                  end   = {params.gr_end}, \
                  name  = '{params.gr_name}', \
                  figure = '{output.figure}' \
                ))"


        """

