#==============================================================================#
#            Analysis and figures included in manuscript figure 2              #
#==============================================================================#
#workdir: "../../"
#workdir: "~/"

rule compile_figure2:
    input:
        figure2a = join(DATAPATH, 'results/figures/figure2/figure2a_tumor_SE_hmatrix.pdf'),
        figure2b = join(DATAPATH, 'results/figures/figure2/figure2b_cells_SE_hmatrix.pdf')
    output: join(DATAPATH, 'results/figures/figure2/figure2_paths.txt')
    shell:
        """
        touch {output}
        echo 'Figure 2a {input.figure2a}' >> {output}
        echo 'Figure 2b {input.figure2b}' >> {output}
        
        """

#================================================================================#
#                      Figure 2b - Tumor SE signal NMF                           #
#================================================================================#
optK_cc = str(config['NMFparams']['cells']['optimalK']['chipseq'])
rule fig2b_cells_SE_heatmap:
    input:
        annotation = join(DATAPATH, 'annotation/annotation_cells.RDS'),
        norm_nmfW = join(DATAPATH, 'analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_SignalScore_normNMF_W.RDS'),
        norm_nmfH = join(DATAPATH, 'analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_SignalScore_normNMF_H.RDS')
    output:
        report    = join(DATAPATH, 'reports/figure2b_cells_ChIPseq_NMF_heatmap.html'),
        rmd       = temp(join(DATAPATH, 'reports/figure2b_cells_ChIPseq_NMF_heatmap.Rmd')),
        hmatrix_wnorm = join(DATAPATH, ('analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_cc + '_Hmatrix_wnorm.RDS')),
        wmatrix_wnorm = join(DATAPATH, ('analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_cc + '_Wmatrix_Wnorm.RDS')),
        nmf_features = join(DATAPATH, ('analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_cc + '_NMF_features.RDS')),
        hmatrix_hnorm = join(DATAPATH, ('analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_cc + '_Hmatrix_hnorm.RDS')),
        figure2b      = join(DATAPATH, 'results/figures/figure2/figure2b_cells_SE_hmatrix.pdf')
    params:
        script   = 'scripts/figure2/figure2b_cells_ChIPseq_NMF.Rmd',
        optimalK = optK_cc
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
                  figure2b      = '{output.figure2b}' \
                ))"


        """

#================================================================================#
#                      Figure 2a - Tumor SE signal NMF                           #
#================================================================================#
optK_tc = str(config['NMFparams']['tumor']['optimalK']['chipseq'])
rule fig2a_tumor_SE_heatmap:
    input:
        annotation = join(DATAPATH, 'annotation/annotation_tumor.RDS'),
        norm_nmfW = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_SignalScore_normNMF_W.RDS'),
        norm_nmfH = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_SignalScore_normNMF_H.RDS')
    output:
        report    = join(DATAPATH, 'reports/figure2a_tumor_ChIPseq_NMF_heatmap.html'),
        rmd       = temp(join(DATAPATH, 'reports/figure2a_tumor_ChIPseq_NMF_heatmap.Rmd')),
        hmatrix_wnorm = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_Hmatrix_wnorm.RDS')),
        wmatrix_wnorm = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_Wmatrix_Wnorm.RDS')),
        nmf_features = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_NMF_features.RDS')),
        hmatrix_hnorm = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_Hmatrix_hnorm.RDS')),
        figure2a      = join(DATAPATH, 'results/figures/figure2/figure2a_tumor_SE_hmatrix.pdf')
    params:
        script   = 'scripts/figure2/figure2a_tumor_ChIPseq_NMF.Rmd',
        optimalK = optK_tc
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
                  figure2a      = '{output.figure2a}' \
                ))"


        """


# rule fig2a_tumor_SE_heatmap:
#     input:
#         annotation = join(DATAPATH, 'annotation/annotation_tumor.RDS'),
#         norm_nmfW = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_SignalScore_normNMF_W.RDS'),
#         norm_nmfH = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_SignalScore_normNMF_H.RDS')
#     output:
#         #report    = join(DATAPATH, 'reports/figure2a_tumor_ChIPseq_NMF_heatmap.html'),
#         #rmd       = temp(join(DATAPATH, 'reports/figure2a_tumor_ChIPseq_NMF_heatmap.Rmd')),
#         hmatrix_wnorm = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_Hmatrix_wnorm.RDS')),
#         wmatrix_wnorm = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_Wmatrix_Wnorm.RDS')),
#         nmf_features = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_NMF_features.RDS')),
#         hmatrix_hnorm = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_Hmatrix_hnorm.RDS')),
#         figure2a      = join(DATAPATH, 'results/figures/figure2/figure2a_tumor_SE_hmatrix.pdf')
#     params:
#         script   = 'scripts/figure2/figure2a_tumor_ChIPseq_NMF.R',
#         optimalK = optK_tc
#     conda: '../envs/R3.5.yaml'
#     shell:
#         """
#
#         Rscript {params.script} {params.optimalK} \
#                                 {input.annotation} \
#                                 {input.norm_nmfW} \
#                                 {input.norm_nmfH} \
#                                 {output.hmatrix_wnorm} \
#                                 {output.hmatrix_wnorm} \
#                                 {output.wmatrix_wnorm} \
#                                 {output.nmf_features} \
#                                 {output.hmatrix_hnorm}
#
#
#         """
