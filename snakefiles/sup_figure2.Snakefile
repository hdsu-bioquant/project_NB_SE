#==============================================================================#
#            Analysis and figures included in manuscript figure 2              #
#==============================================================================#
#workdir: "../../"
#workdir: "~/"

rule compile_sup_figure2:
    input:
        #sup_figure2a = join(DATAPATH, 'results/figures/figure2/figure2a_tumor_SE_hmatrix.pdf'),
        #sup_figure2b = join(DATAPATH, 'results/figures/figure2/figure2b_cells_SE_hmatrix.pdf')
        sup_figure2c = join(DATAPATH, 'results/figures/sup_figure2/figure2c_tumor_cells_SE_hmatrix.pdf')
    output: join(DATAPATH, 'results/figures/sup_figure2/sup_figure2_paths.txt')
    shell:
        """
        touch {output}
        #echo 'Sup. Figure 2a {input.sup_figure2c}' >> {output}
        #echo 'Sup. Figure 2b {input.sup_figure2c}' >> {output}
        echo 'Sup. Figure 2c {input.sup_figure2c}' >> {output}
        
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
        sup_figure2c  = join(DATAPATH, 'results/figures/sup_figure2/figure2c_tumor_cells_SE_hmatrix.pdf')
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


