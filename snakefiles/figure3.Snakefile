#==============================================================================#
#            Analysis and figures included in manuscript figure 3              #
#==============================================================================#


rule compile_figure3:
    input:
        figure3_CRC1 = join(DATAPATH, 'results/figure3/crcTF_fraction_observed_tumor_and cells_signatures_combined.pdf'),
        figure3_CRC2 = join(DATAPATH, 'results/figure3/crcTF_TF_activity_tumor_and cells_signatures_combined.pdf'),
        figure3_CRC3 = join(DATAPATH, 'results/figure3/crcTF_correlation_TFactivity_vs_KD.pdf'),
        figure3_CCND11 = join(DATAPATH,"results/figure3/SEtargetGenes_DiffKDprofileNBcellsVsRest.pdf"),
        figure3_CCND12 = join(DATAPATH,"results/figure3/Kelly_SKNA_KDprofile_topHits.pdf"),
        figure3_foot   = join(DATAPATH, 'results/figure3/figure3_MES_vs_ADRN_footprint.pdf'),
        figure_enhan= join(DATAPATH, 'results/figure3/Super_Enhancer_interaction.pdf'),
        figure3g    = join(DATAPATH, 'results/figure3/figure3g_IGV_plot.pdf')
    output: join(DATAPATH, 'results/figure3/figure3_paths.txt')
    shell:
        """
        touch {output}
        echo 'Figure 3 CRC1 {input.figure3_CRC1}' >> {output}
        echo 'Figure 3 CRC2 {input.figure3_CRC2}' >> {output}
        echo 'Figure 3 CRC3 {input.figure3_CRC3}' >> {output}
        echo 'Figure 3 CCND1 1 {input.figure3_CCND11}' >> {output}
        echo 'Figure 3 CCND1 2 {input.figure3_CCND12}' >> {output}
        echo 'Figure 3 footprint {input.figure3_foot}' >> {output}
        echo 'Figure 3g {input.figure3g}' >> {output}
        
        """

#================================================================================#
#             Figure 3 - Enhancer Comparison using HiChIP                        #
#================================================================================#
optK_tc = str(config['NMFparams']['tumor']['optimalK']['chipseq'])
optK_cc = str(config['NMFparams']['cells']['optimalK']['chipseq'])
rule fig3_Enhancer_Comparison_using_HiChIP:
    input:
        consensusSE  = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed'),
        cells_h      = join(DATAPATH, 'analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_cc + '_Hmatrix_hnorm.RDS'),
        cells_w      = join(DATAPATH, 'analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_cc + '_Wmatrix_Wnorm.RDS'),
        tumor_h      = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K' + optK_tc + '_Hmatrix_hnorm.RDS'),
        SKNAS_HiChIP = join(DATAPATH, 'data/cells/hichip/mango/SK-N-AS_HiChIP_mango.all'),
        CLBGA_HiChIP = join(DATAPATH, 'data/cells/hichip/mango/CLB-GA_HiChIP_mango.all')
    output:
        report    = join(DATAPATH, 'reports/Enhancer_Comparison_using_HiChIP.html'),
        rmd       = temp(join(DATAPATH, 'reports/Enhancer_Comparison_using_HiChIP.Rmd')),
        figure    = join(DATAPATH, 'results/figure3/Super_Enhancer_interaction.pdf'),
        supfigure = join(DATAPATH, 'results/sup_figure3/Enhancer_interaction.pdf')
    params:
        script   = 'scripts/figure3/Enhancer_Comparison_using_HiChIP.Rmd',
        work_dir = DATAPATH
    conda: '../envs/R3.5_diffbind.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  work_dir     = '{params.work_dir}', \
                  SE           = '{input.consensusSE}', \
                  cells_h      = '{input.cells_h}', \
                  cells_w      = '{input.cells_w}', \
                  tumor_h      = '{input.tumor_h}', \
                  SKNAS_HiChIP = '{input.SKNAS_HiChIP}', \
                  CLBGA_HiChIP = '{input.CLBGA_HiChIP}', \
                  figure       = '{output.figure}', \
                  sup_figure   = '{output.supfigure}' \
                ))"


        """

#================================================================================#
#                      Figure 3 - MES vs. ADRN footprint                         #
#================================================================================#
rule fig3_MES_vs_ADRN_footprint:
    input:
        consensusSE    = join(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS'),
        cellline1_foot = join(DATAPATH, 'data/cells/atacseq/footprint/KELLY_footprints_calls_GrangesList.RDS'),
        cellline2_foot = join(DATAPATH, 'data/cells/atacseq/footprint/SK-N-AS_footprints_calls_GrangesList.RDS'),
        MES_activity   = join(DATAPATH, 'analysis/tumor/VIPER/MES_TFactivity.RDS'),
        tumor_CRCs     = join(DATAPATH, 'data/tumor/chipseq/H3K27ac/CRC/'),
        cells_CRCs     = join(DATAPATH, 'data/cells/chipseq/H3K27ac/CRC/')
    output:
        report = join(DATAPATH, 'reports/figure3_MES_vs_ADRN_footprint.html'),
        rmd    = temp(join(DATAPATH, 'reports/figure3_MES_vs_ADRN_footprint.Rmd')),
        figure = join(DATAPATH, 'results/figure3/figure3_MES_vs_ADRN_footprint.pdf')
    params:
        script   = 'scripts/figure3/figure3_MES_vs_ADRN_footprint.Rmd',
        work_dir = DATAPATH
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  SE             = '{input.consensusSE}', \
                  cellline1_foot = '{input.cellline1_foot}', \
                  cellline2_foot = '{input.cellline2_foot}', \
                  MES_activity   = '{input.MES_activity}', \
                  tumor_CRCs     = '{input.tumor_CRCs}', \
                  cells_CRCs     = '{input.cells_CRCs}', \
                  figure         = '{output.figure}' \
                ))"


        """



#================================================================================#
#                      Figure 3 - CCND1 loci                                     #
#================================================================================#
rule fig3_IGV:
    input:
        consensusSE = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed')
    output:
        report = join(DATAPATH, 'reports/figure3_IGV_plot.html'),
        rmd    = temp(join(DATAPATH, 'reports/figure3_IGV_plot.Rmd')),
        figure = join(DATAPATH, 'results/figure3/figure3g_IGV_plot.pdf')
    params:
        script   = 'scripts/figure3/figure3_IGV_plot.Rmd',
        work_dir = DATAPATH,
        highbed  = join(os.getcwd(), 'scripts/figure3/figure3_IGV_highlight_ranges.bed'),
        window   = config['igv_plot']['figure3_CCND1']['window'],
        gr_chr   = config['igv_plot']['figure3_CCND1']['chr'],
        gr_start = config['igv_plot']['figure3_CCND1']['start'],
        gr_end   = config['igv_plot']['figure3_CCND1']['end'],
        gr_name  = config['igv_plot']['figure3_CCND1']['name'],
        width    = config['igv_plot']['figure3_CCND1']['width'],
        height   = config['igv_plot']['figure3_CCND1']['height']
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  work_dir     = '{params.work_dir}', \
                  path_high_bed= '{params.highbed}', \
                  width_window = {params.window}, \
                  chr    = '{params.gr_chr}', \
                  start  = {params.gr_start}, \
                  end    = {params.gr_end}, \
                  name   = '{params.gr_name}', \
                  figure = '{output.figure}', \
                  width  = {params.width}, \
                  height = {params.height} \
                ))"


        """

#================================================================================#
# CRCs related plots -                                                           #  
# (1) CRC fraction observed heatmap                                              #
# (2) CRC occurrence per sample oncoprints                                       #
# (3) CRC knockdown vs TF activity correlation (Mesenchymal vs Adrenergic)       #
#================================================================================#

rule fig3_CRCplots:
    input:
        KDdata   = join(DATAPATH, 'analysis/cells/crcGIEMSAkd/nbKDinhouse.RDS'),
        network  = join(DATAPATH, 'analysis/tumor/VIPER/networkViper.txt'),
        tumorNMF = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS'),
        cellNMF  = join(DATAPATH, 'analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K3_Hmatrix_hnorm.RDS'),
        tumorCRC = join(DATAPATH, 'data/tumor/chipseq/H3K27ac/CRC/'),
        cellCRC  = join(DATAPATH, 'data/cells/chipseq/H3K27ac/CRC/'),
        t_annot  = join(DATAPATH, 'annotation/annotation_tumor.RDS'),
        c_annot  = join(DATAPATH, 'annotation/annotation_cells.RDS')
    output:
        crcList   = join(DATAPATH, 'results/supptables/crcTF_fractionObserved_combined.txt'),
        tabSupp_1 = join(DATAPATH, 'results/supptables/TFactivity_across_all_signatures_ZnormPerSig.txt'),
        tabSupp_2 = join(DATAPATH, 'results/supptables/crcTF_modules.txt'),

        figSupp_1 = join(DATAPATH, 'results/sup_figure3/TFactivity_heatmap_across_all_signatures_ZnormPerSig.pdf'),
        figSupp_2 = join(DATAPATH, 'results/sup_figure3/crcTF_TFactivity_per_crcModule_and_Signature_labelled.pdf'),
        figSupp_3 = join(DATAPATH, 'results/sup_figure3/crcTF_fraction_observed_tumors_only.pdf'),
        figSupp_4 = join(DATAPATH, 'results/sup_figure3/crcTF_fraction_observed_cells_only.pdf'),
        figSupp_5 = join(DATAPATH, 'results/sup_figure3/crcTF_oncoprints_Signature_Specific.pdf'),

        figMain_1 = join(DATAPATH, 'results/figure3/crcTF_fraction_observed_tumor_and cells_signatures_combined.pdf'),
        figMain_2 = join(DATAPATH, 'results/figure3/crcTF_TF_activity_tumor_and cells_signatures_combined.pdf'),
        figMain_3 = join(DATAPATH, 'results/figure3/crcTF_correlation_TFactivity_vs_KD.pdf')
    params:
        TFact    = join(DATAPATH, 'analysis/tumor/VIPER/'),
        script  = 'scripts/figure3/figure3_CRCplots.R',
        outpath = join(DATAPATH, 'results/')
    conda: '../envs/R3.5.yaml'
    shell:
        """
        Rscript {params.script} {input.KDdata} {params.TFact} {input.tumorNMF} \
                {input.cellNMF} {input.tumorCRC} {input.cellCRC} \
                {input.t_annot} {input.c_annot} {params.outpath} \
                
        """


#================================================================================#
# In depth analysis of CCND1 target                                              #
#================================================================================#

rule fig3_CRCknockdown:
    input:
        SEtarget = join(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS'),
        KDvals   = join(DATAPATH, 'db/DeepMap19Q2/cellKnockdownCERES.RDS'),
        cellExpr = join(DATAPATH, 'db/DeepMap19Q2/cellExpression.RDS')
    output:
        figMain_1 = join(DATAPATH,"results/figure3/SEtargetGenes_DiffKDprofileNBcellsVsRest.pdf"),
        figMain_2 = join(DATAPATH,"results/figure3/Kelly_SKNA_KDprofile_topHits.pdf"),

        figSupp_1 = join(DATAPATH,"results/sup_figure3/CCND1_KDprofile_across_tumorTypes_cells.pdf"),
        figSupp_2 = join(DATAPATH,"results/sup_figure3/CCND1_ExpressionProfile_across_tumorTypes_cells.pdf"),
        figSupp_3 = join(DATAPATH,"results/sup_figure3/CCND1_ExpressionProfile_NBcellsVsRest.pdf")
    params:
        script  = 'scripts/figure3/CCND1plots.R',
        outpath = join(DATAPATH, 'results/')
    conda: '../envs/R3.5.yaml'
    shell:
        """
        Rscript {params.script} {input.SEtarget} {input.KDvals} {input.cellExpr} {params.outpath}
        """

