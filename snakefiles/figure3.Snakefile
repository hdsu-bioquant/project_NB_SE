#==============================================================================#
#            Analysis and figures included in manuscript figure 3              #
#==============================================================================#


rule compile_figure3:
    input:
        figure3_CRC1 = join(DATAPATH, 'results/figure3/crcTF_fractionObserved_tumors_cells.pdf'),
        figure3_CRC2 = join(DATAPATH, 'results/figure3/crcTF_TFactivity_per_Major_class_and_Signature.pdf'),
        figure3_CRC3 = join(DATAPATH, 'results/figure3/crcTF_correlation_TFactivity_vs_KD.pdf'),
        figure3_CRC4 = join(DATAPATH, 'results/figure3/crcTF_oncoprints_SignatureSpecific.pdf'),
        figure3_CCND11 = join(DATAPATH,"results/figure3/SEtargetGenes_DiffKDprofileNBcellsVsRest.pdf"),
        figure3_CCND12 = join(DATAPATH,"results/figure3/Kelly_SKNA_KDprofile_topHits.pdf"),
        figure3g    = join(DATAPATH, 'results/figure3/figure3g_IGV_plot.pdf')
    output: join(DATAPATH, 'results/figure3/figure3_paths.txt')
    shell:
        """
        touch {output}
        echo 'Figure 3 CRC1 {input.figure3_CRC1}' >> {output}
        echo 'Figure 3 CRC2 {input.figure3_CRC2}' >> {output}
        echo 'Figure 3 CRC3 {input.figure3_CRC3}' >> {output}
        echo 'Figure 3 CRC4 {input.figure3_CRC4}' >> {output}
        echo 'Figure 3 CCND1 1 {input.figure3_CCND11}' >> {output}
        echo 'Figure 3 CCND1 2 {input.figure3_CCND12}' >> {output}
        echo 'Figure 3g {input.figure3g}' >> {output}
        
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
        TFact    = join(DATAPATH, 'analysis/tumor/VIPER/'),
        tumorNMF = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS'),
        cellNMF  = join(DATAPATH, 'analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K3_Hmatrix_hnorm.RDS'),
        tumorCRC = join(DATAPATH, 'data/tumor/chipseq/H3K27ac/CRC/'),
        cellCRC  = join(DATAPATH, 'data/cells/chipseq/H3K27ac/CRC/')
    output:
        tabSupp_1 = join(DATAPATH, 'results/supptables/TFactivity_across_all_signatures_ZnormPerSig.txt'),
        tabSupp_2 = join(DATAPATH, 'results/supptables/crcTF_modules.txt'),

        figSupp_1 = join(DATAPATH, 'results/sup_figure3/TFactivity_heatmap_across_all_signatures_ZnormPerSig.pdf'),
        figSupp_2 = join(DATAPATH, 'results/sup_figure3/crcTF_TFactivity_per_Minor_class_and_Signature.pdf'),

        figMain_1 = join(DATAPATH, 'results/figure3/crcTF_fractionObserved_tumors_cells.pdf'),
        figMain_2 = join(DATAPATH, 'results/figure3/crcTF_TFactivity_per_Major_class_and_Signature.pdf'),
        figMain_3 = join(DATAPATH, 'results/figure3/crcTF_correlation_TFactivity_vs_KD.pdf'),
        figMain_4 = join(DATAPATH, 'results/figure3/crcTF_oncoprints_SignatureSpecific.pdf')
    params:
        script  = 'scripts/figure3/figure3_CRCplots.R',
        outpath = join(DATAPATH, 'results/')
    conda: '../envs/R3.5.yaml'
    shell:
        """
        Rscript {params.script} {input.KDdata} {input.TFact} {input.tumorNMF} \
                {input.cellNMF} {input.tumorCRC} {input.cellCRC} {params.outpath}
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

