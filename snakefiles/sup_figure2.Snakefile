#==============================================================================#
#            Analysis and figures included in manuscript figure 2              #
#==============================================================================#
#workdir: "../../"
#workdir: "~/"

rule compile_sup_figure2:
    input:
        #sup_figure2a = join(DATAPATH, 'results/figures/figure2/figure2a_tumor_SE_hmatrix.pdf'),
        #sup_figure2b = join(DATAPATH, 'results/figures/figure2/figure2b_cells_SE_hmatrix.pdf')
        sup_figure2c = join(DATAPATH, 'results/sup_figure2/sup_figure2c_tumor_cells_SE_hmatrix.pdf'),
        sup_figure2d = join(DATAPATH, 'results/sup_figure2/sup_figure2d_tumor_corr_SEvsExprs_exposure.pdf'),
        sup_figure2e = join(DATAPATH, 'results/sup_figure2/sup_figure2e_tumor_mostVariablehmatrix.pdf')
    output: join(DATAPATH, 'results/sup_figure2/sup_figure2_paths.txt')
    shell:
        """
        touch {output}
        #echo 'Sup. Figure 2a {input.sup_figure2c}' >> {output}
        #echo 'Sup. Figure 2b {input.sup_figure2c}' >> {output}
        echo 'Sup. Figure 2c {input.sup_figure2c}' >> {output}
        echo 'Sup. Figure 2d {input.sup_figure2d}' >> {output}
        echo 'Sup. Figure 2e {input.sup_figure2e}' >> {output}
        
        """
    
#================================================================================#
#          Sup. Figure 2e - Tumor Most Variable Genes RNAseq NMF                 #
#================================================================================#
optK_tr = str(config['NMFparams']['tumor']['optimalK']['rnaseq'])
rule sup_fig2e_tumor_mostVarible_heatmap:
    input:
        annotation = join(DATAPATH, 'annotation/annotation_tumor.RDS'),
        norm_nmfW = join(DATAPATH, 'analysis/tumor/rnaseq/NMF_mostVariable/tumor_mostVariable_normNMF_W.RDS'),
        norm_nmfH = join(DATAPATH, 'analysis/tumor/rnaseq/NMF_mostVariable/tumor_mostVariable_normNMF_H.RDS')
    output:
        report    = join(DATAPATH, 'reports/sup_figure2e_tumor_mostVariable_Hmatrix.html'),
        rmd       = temp(join(DATAPATH, 'reports/sup_figure2e_tumor_mostVariable_Hmatrix.Rmd')),
        hmatrix_wnorm = join(DATAPATH, ('analysis/tumor/rnaseq/NMF_mostVariable/tumor_mostVariable_K' + optK_tr + '_Hmatrix_wnorm.RDS')),
        wmatrix_wnorm = join(DATAPATH, ('analysis/tumor/rnaseq/NMF_mostVariable/tumor_mostVariable_K' + optK_tr + '_Wmatrix_Wnorm.RDS')),
        nmf_features  = join(DATAPATH, ('analysis/tumor/rnaseq/NMF_mostVariable/tumor_mostVariable_K' + optK_tr + '_NMF_features.RDS')),
        hmatrix_hnorm = join(DATAPATH, ('analysis/tumor/rnaseq/NMF_mostVariable/tumor_mostVariable_K' + optK_tr + '_Hmatrix_hnorm.RDS')),
        sup_figure2e  = join(DATAPATH, 'results/sup_figure2/sup_figure2e_tumor_mostVariablehmatrix.pdf')
    params:
        script   = 'scripts/sup_figure2/sup_figure2e_tumor_mostVariable_Hmatrix.Rmd',
        optimalK = optK_tr
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  K         = {params.optimalK}, \
                  annot     = '{input.annotation}', \
                  norm_nmfW = '{input.norm_nmfW}', \
                  norm_nmfH = '{input.norm_nmfH}', \
                  hmatrix_wnorm = '{output.hmatrix_wnorm}', \
                  wmatrix_wnorm = '{output.wmatrix_wnorm}', \
                  nmf_features  = '{output.nmf_features}', \
                  hmatrix_hnorm = '{output.hmatrix_hnorm}', \
                  sup_figure2e  = '{output.sup_figure2e}' \
                ))"


        """

#================================================================================#
#    Sup. Figure 2d - Correlation of SE signal exposure to SE target exposure    #
#================================================================================#
optK_tc = str(config['NMFparams']['tumor']['optimalK']['chipseq'])
optK_tr = str(config['NMFparams']['tumor']['optimalK']['rnaseq'])
rule sup_fig2d_tumor_corr_SEvsExprs_exposure:
    input:
        h_SEsig = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_Hmatrix_wnorm.RDS')),
        h_exprs = join(DATAPATH, ('analysis/tumor/rnaseq/NMF/tumor_consensusSE_K' + optK_tr + '_Hmatrix_wnorm.RDS'))
    output:
        report = join(DATAPATH, 'reports/sup_figure2d_tumor_corr_SEvsExprs_exposure.html'),
        rmd    = temp(join(DATAPATH, 'reports/sup_figure2d_tumor_corr_SEvsExprs_exposure.Rmd')),
        figure = join(DATAPATH, 'results/sup_figure2/sup_figure2d_tumor_corr_SEvsExprs_exposure.pdf')
    params:
        script   = 'scripts/sup_figure2/sup_figure2d_tumor_corr_SEvsExprs_exposure.Rmd'
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  h_SEsig = '{input.h_SEsig}', \
                  h_exprs = '{input.h_exprs}', \
                  figure  = '{output.figure}' \
                ))"


        """

#================================================================================#
#           Sup. Figure 2c - Tumor and Cell lines SE signal H matrix             #
#================================================================================#
optK_tcc = str(config['NMFparams']['tumor_cells']['optimalK']['chipseq'])
rule sup_fig2c_tumor_cells_SE_Hmatrix:
    input:
        annot_tumor   = join(DATAPATH, 'annotation/annotation_tumor.RDS'),
        annot_cells   = join(DATAPATH, 'annotation/annotation_cells.RDS'),
        hmatrix_wnorm = join(DATAPATH, ('analysis/tumor_cells/chipseq/H3K27ac/NMF/tumor_cells_consensusSE_K' + optK_tcc + '_Hmatrix_wnorm.RDS'))
    output:
        report        = join(DATAPATH, 'reports/sup_figure2c_tumor_cells_Hmatrix.html'),
        rmd           = temp(join(DATAPATH, 'reports/sup_figure2c_tumor_cells_Hmatrix.Rmd')),
        sup_figure2c  = join(DATAPATH, 'results/sup_figure2/sup_figure2c_tumor_cells_SE_hmatrix.pdf')
    params:
        script   = 'scripts/sup_figure2/sup_figure2c_tumor_cells_Hmatrix.Rmd'
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  annot_tumor   = '{input.annot_tumor}', \
                  annot_cells   = '{input.annot_cells}', \
                  hmatrix_wnorm = '{input.hmatrix_wnorm}', \
                  sup_figure2c  = '{output.sup_figure2c}' \
                ))"


        """


