#==============================================================================#
#            Analysis and figures included in manuscript figure 2              #
#==============================================================================#
#workdir: "../../"
#workdir: "~/"

rule compile_figure2:
    input:
        figure2a   = join(DATAPATH, 'results/figure2/figure2a_tumor_SE_hmatrix.pdf'),
        figure2b   = join(DATAPATH, 'results/figure2/figure2b_cells_SE_hmatrix.pdf'),
        figure2c   = join(DATAPATH, 'results/figure2/figure2c_tumor_cells_SE_UMAP.pdf'),
        figure2d01 = join(DATAPATH, 'results/figure2/figure2d_01_tumor_riverplot.pdf'),
        figure2d02 = join(DATAPATH, 'results/figure2/figure2d_02_cells_riverplot.pdf'),
        figure2e   = join(DATAPATH, 'results/figure2/figure2e_tumor_SE_targets_hmatrix.pdf'),
        figure2f   = join(DATAPATH, 'results/figure2/figure2f_tumor_cells_density.pdf'),
        figure2g   = join(DATAPATH, 'results/figure2/figure2g_tumor_SE_targets_NMF_recovery.pdf'),
        figure2h   = join(DATAPATH, 'results/figure2/figure2h_SEtargetExprs_mouseE12.5.pdf'),
        figure_pan = join(DATAPATH, 'results/figure2/figure2panels_signature_enrichment.pdf')
    output: join(DATAPATH, 'results/figure2/figure2_paths.txt')
    shell:
        """
        touch {output}
        echo 'Figure 2a {input.figure2a}' >> {output}
        echo 'Figure 2b {input.figure2b}' >> {output}
        echo 'Figure 2c {input.figure2c}' >> {output}
        echo 'Figure 2d {input.figure2d01}' >> {output}
        echo 'Figure 2d {input.figure2d02}' >> {output}
        echo 'Figure 2e {input.figure2e}' >> {output}
        echo 'Figure 2f {input.figure2f}' >> {output}
        echo 'Figure 2g {input.figure2g}' >> {output}
        echo 'Figure 2h {input.figure2h}' >> {output}
        echo 'Figure pan {input.figure_pan}' >> {output}
        
        """

#================================================================================#
#           Figure 2h - expression mapped to Mouse GSE99933 E12.5                #
#================================================================================#
optK_tc = str(config['NMFparams']['tumor']['optimalK']['chipseq'])
optK_tr = str(config['NMFparams']['tumor']['optimalK']['rnaseq'])
rule fig2h_SEexprs_mouse:
    input:
        SE_target    = join(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS'),
        wc_tumor     = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_Wmatrix_Wnorm.RDS')),
        wr_tumor     = join(DATAPATH, ('analysis/tumor/rnaseq/NMF/tumor_consensusSE_K' + optK_tr + '_Wmatrix_Wnorm.RDS')),
        mouse_pstime = join(DATAPATH, 'db/GSE99933_E12.5/GSE99933_E12.5.txt'),
        mouse_exprs  = join(DATAPATH, 'db/GSE99933_E12.5/GSE99933_E12.5_exprs_Zscore.txt')
    output:
        report = join(DATAPATH, 'reports/figure2h_SEtargetExprs_mouseE12.5.html'),
        rmd    = temp(join(DATAPATH, 'reports/figure2h_SEtargetExprs_mouseE12.5.Rmd')),
        figure = join(DATAPATH, 'results/figure2/figure2h_SEtargetExprs_mouseE12.5.pdf')
    params:
        script   = 'scripts/figure2/figure2h_SEtargetExprs_mouseE12.5.Rmd'
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  SE_target    = '{input.SE_target}', \
                  wc_tumor     = '{input.wc_tumor}', \
                  wr_tumor     = '{input.wr_tumor}', \
                  mouse_pstime = '{input.mouse_pstime}', \
                  mouse_exprs  = '{input.mouse_exprs}', \
                  figure = '{output.figure}' \
                ))"


        """


#================================================================================#
#          Figure 2g - Tumor SE Target RNAseq NMF recovery plots                 #
#================================================================================#
optK_tr = str(config['NMFparams']['tumor']['optimalK']['rnaseq'])
rule fig2g_tumor_SE_Targets_recovery:
    input:
        annotation = join(DATAPATH, 'annotation/annotation_tumor.RDS'),
        hmatrix_wnorm = join(DATAPATH, ('analysis/tumor/rnaseq/NMF/tumor_consensusSE_K' + optK_tr + '_Hmatrix_wnorm.RDS'))
    output:
        report    = join(DATAPATH, 'reports/figure2g_tumor_RNAseq_recovery_plots.html'),
        rmd       = temp(join(DATAPATH, 'reports/figure2g_tumor_RNAseq_recovery_plots.Rmd')),
        figure2g  = join(DATAPATH, 'results/figure2/figure2g_tumor_SE_targets_NMF_recovery.pdf')
    params:
        script   = 'scripts/figure2/figure2g_tumor_RNAseq_recovery_plots.Rmd'
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  annot         = '{input.annotation}', \
                  hmatrix_wnorm = '{input.hmatrix_wnorm}', \
                  figure2g      = '{output.figure2g}' \
                ))"


        """

#================================================================================#
#              Figure 2f - Tumor and Cell MES exposure Density                   #
#================================================================================#
optK_tc = str(config['NMFparams']['tumor']['optimalK']['chipseq'])
optK_cc = str(config['NMFparams']['cells']['optimalK']['chipseq'])
rule fig2f_tumor_cells_density:
    input:
        h_tumor = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_Hmatrix_wnorm.RDS')),
        h_cells = join(DATAPATH, ('analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_cc + '_Hmatrix_wnorm.RDS'))
    output:
        report   = join(DATAPATH, 'reports/figure2f_tumor_cells_density.html'),
        rmd      = temp(join(DATAPATH, 'reports/figure2f_tumor_cells_density.Rmd')),
        figure2f = join(DATAPATH, 'results/figure2/figure2f_tumor_cells_density.pdf')
    params:
        script   = 'scripts/figure2/figure2f_tumor_cells_density.Rmd'
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  h_tumor  = '{input.h_tumor}', \
                  h_cells  = '{input.h_cells}', \
                  figure2f = '{output.figure2f}' \
                ))"


        """

#================================================================================#
#            Figure 2 panels - Top NMF Signature's feature enrichment            #
#================================================================================#
optK_tc = str(config['NMFparams']['tumor']['optimalK']['chipseq'])
optK_cc = str(config['NMFparams']['cells']['optimalK']['chipseq'])
optK_tr = str(config['NMFparams']['tumor']['optimalK']['rnaseq'])


rule fig2panels_signature_enrichment:
    input:
        SE_target   = join(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS'),
        SE_target_df= join(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_annotation_df.RDS'),
        terms       = join(DATAPATH, 'results/supptables/GO_BP_enrichment_SE_target_genes.txt'),
        wChIP_tumor = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_Wmatrix_Wnorm.RDS')),
        wChIP_cells = join(DATAPATH, ('analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_cc + '_Wmatrix_Wnorm.RDS')),
        wRNAs_tumor = join(DATAPATH, ('analysis/tumor/rnaseq/NMF/tumor_consensusSE_K' + optK_tr + '_Wmatrix_Wnorm.RDS')),
        wRNAv_tumor = join(DATAPATH, ('analysis/tumor/rnaseq/NMF_mostVariable/tumor_mostVariable_K' + optK_tr + '_Wmatrix_Wnorm.RDS'))
    output:
        report      = join(DATAPATH, 'reports/figure2panels_signature_enrichment.html'),
        rmd         = temp(join(DATAPATH, 'reports/figure2panels_signature_enrichment.Rmd')),
        eChIP_tumor = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K4_GO_BP_enrichment.RDS'),
        eChIP_cells = join(DATAPATH, 'analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K3_GO_BP_enrichment.RDS'),
        eRNAs_tumor = join(DATAPATH, 'analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_GO_BP_enrichment.RDS'),
        eRNAv_tumor = join(DATAPATH, 'analysis/tumor/rnaseq/NMF_mostVariable/tumor_mostVariable_K4_GO_BP_enrichment.RDS'),
        figure      = join(DATAPATH, 'results/figure2/figure2panels_signature_enrichment.pdf'),
        sup_figure  = join(DATAPATH, 'results/sup_figure2/sup_figure2_MostVariable_signature_enrichment.pdf')
    params:
        script   = 'scripts/figure2/figure2panels_signature_enrichment.Rmd',
        optimalK = optK_tr
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  SE             = '{input.SE_target}', \
                  targetannot_df = '{input.SE_target_df}', \
                  terms          = '{input.terms}', \
                  wChIP_tumor    = '{input.wChIP_tumor}', \
                  wChIP_cells    = '{input.wChIP_cells}', \
                  wRNAs_tumor    = '{input.wRNAs_tumor}', \
                  wRNAv_tumor    = '{input.wRNAv_tumor}', \
                  eChIP_tumor= '{output.eChIP_tumor}', \
                  eChIP_cells= '{output.eChIP_cells}', \
                  eRNAs_tumor= '{output.eRNAs_tumor}', \
                  eRNAv_tumor= '{output.eRNAv_tumor}', \
                  figure         = '{output.figure}', \
                  sup_figure     = '{output.sup_figure}' \
                ))"


        """


#================================================================================#
#                    Figure 2e - Tumor SE Target RNAseq NMF                      #
#================================================================================#
optK_tr = str(config['NMFparams']['tumor']['optimalK']['rnaseq'])
rule fig2e_tumor_SE_Targets_heatmap:
    input:
        annotation = join(DATAPATH, 'annotation/annotation_tumor.RDS'),
        norm_nmfW = join(DATAPATH, 'analysis/tumor/rnaseq/NMF/tumor_consensusSE_targetExprs_normNMF_W.RDS'),
        norm_nmfH = join(DATAPATH, 'analysis/tumor/rnaseq/NMF/tumor_consensusSE_targetExprs_normNMF_H.RDS')        
    output:
        report    = join(DATAPATH, 'reports/figure2e_tumor_RNAseq_NMF.html'),
        rmd       = temp(join(DATAPATH, 'reports/figure2e_tumor_RNAseq_NMF.Rmd')),
        hmatrix_wnorm = join(DATAPATH, ('analysis/tumor/rnaseq/NMF/tumor_consensusSE_K' + optK_tr + '_Hmatrix_wnorm.RDS')),
        wmatrix_wnorm = join(DATAPATH, ('analysis/tumor/rnaseq/NMF/tumor_consensusSE_K' + optK_tr + '_Wmatrix_Wnorm.RDS')),
        nmf_features  = join(DATAPATH, ('analysis/tumor/rnaseq/NMF/tumor_consensusSE_K' + optK_tr + '_NMF_features.RDS')),
        hmatrix_hnorm = join(DATAPATH, ('analysis/tumor/rnaseq/NMF/tumor_consensusSE_K' + optK_tr + '_Hmatrix_hnorm.RDS')),
        figure2e      = join(DATAPATH, 'results/figure2/figure2e_tumor_SE_targets_hmatrix.pdf')
    params:
        script   = 'scripts/figure2/figure2e_tumor_RNAseq_NMF.Rmd',
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
                  figure2e      = '{output.figure2e}' \
                ))"


        """


#================================================================================#
#              Figure 2d - Tumor and Cell lines Riverplots.                      #
#================================================================================#
rule fig2d_tumor_cells_riverplots:
    input:
        nmf_tumor = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_SignalScore_normNMF_W.RDS'),
        nmf_cells = join(DATAPATH, 'analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_SignalScore_normNMF_W.RDS')
    output:
        report     = join(DATAPATH, 'reports/figure2d_tumor_cells_riverplots.html'),
        rmd        = temp(join(DATAPATH, 'reports/figure2d_tumor_cells_riverplots.Rmd')),
        figure2d01 = join(DATAPATH, 'results/figure2/figure2d_01_tumor_riverplot.pdf'),
        figure2d02 = join(DATAPATH, 'results/figure2/figure2d_02_cells_riverplot.pdf')
    params:
        script   = 'scripts/figure2/figure2d_tumor_cells_riverplots.Rmd'
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  nmf_tumor      = '{input.nmf_tumor}', \
                  nmf_cells      = '{input.nmf_cells}', \
                  figure2d_tumor = '{output.figure2d01}', \
                  figure2d_cells = '{output.figure2d02}' \
                ))"


        """



#================================================================================#
#              Figure 2c - Tumor and Cell lines SE signal UMAP                   #
#================================================================================#
optK_tcc = str(config['NMFparams']['tumor_cells']['optimalK']['chipseq'])
rule fig2c_tumor_cells_SE_umap:
    input:
        annot_tumor = join(DATAPATH, 'annotation/annotation_tumor.RDS'),
        annot_cells = join(DATAPATH, 'annotation/annotation_cells.RDS'),
        matrix      = join(DATAPATH, 'analysis/tumor_cells/chipseq/H3K27ac/consensusSE/tumor_cells_H3K27ac_noH3K4me3_consensusSE_SignalScore.RDS'),
        hmatrix_wnorm = join(DATAPATH, ('analysis/tumor_cells/chipseq/H3K27ac/NMF/tumor_cells_consensusSE_K' + optK_tcc + '_Hmatrix_wnorm.RDS'))
    output:
        report    = join(DATAPATH, 'reports/figure2c_tumor_cells_UMAP.html'),
        rmd       = temp(join(DATAPATH, 'reports/figure2c_tumor_cells_UMAP.Rmd')),
        figure2c  = join(DATAPATH, 'results/figure2/figure2c_tumor_cells_SE_UMAP.pdf')
    params:
        script   = 'scripts/figure2/figure2c_tumor_cells_UMAP.Rmd'
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  annot_tumor   = '{input.annot_tumor}', \
                  annot_cells   = '{input.annot_cells}', \
                  SE_signal     = '{input.matrix}', \
                  hmatrix_wnorm = '{input.hmatrix_wnorm}', \
                  figure2c      = '{output.figure2c}' \
                ))"


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
        nmf_features  = join(DATAPATH, ('analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_cc + '_NMF_features.RDS')),
        hmatrix_hnorm = join(DATAPATH, ('analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_cc + '_Hmatrix_hnorm.RDS')),
        figure2b      = join(DATAPATH, 'results/figure2/figure2b_cells_SE_hmatrix.pdf')
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
        nmf_features  = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_NMF_features.RDS')),
        hmatrix_hnorm = join(DATAPATH, ('analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_Hmatrix_hnorm.RDS')),
        figure2a      = join(DATAPATH, 'results/figure2/figure2a_tumor_SE_hmatrix.pdf')
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
