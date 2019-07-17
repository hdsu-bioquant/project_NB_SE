#==============================================================================#
#            Analysis and figures included in manuscript figure 4              #
#==============================================================================#

rule compile_figure4:
    input:
        figure4a = join(DATAPATH,"results/figure4/junfos_corr_exposures.pdf"),
        figure4b = join(DATAPATH,"results/figure4/ras_corr_exposures.pdf"),
        figure4cd = join(DATAPATH, 'results/figure4/RasJunFos_Exprs_mouseE12.5.pdf'),
        #figure4e = join(DATAPATH, 'results/figure2/figure2e_tumor_SE_targets_hmatrix.pdf'),
        #figure4f = join(DATAPATH, 'results/figure2/figure2f_tumor_cells_density.pdf'),
        #figure4i = join(DATAPATH, 'results/figure4/figure4i_IGV_plot.pdf')
        figure4ij = expand(join(DATAPATH, 'results/figure4/figure4_{ExamplegeneID}_loci.pdf'), zip, ExamplegeneID = config['igv_plot']['figure4']['name']),
        figure4g = join(DATAPATH, 'results/figure4/tumors_RNAseq_PrimaryVsRelapse.pdf')
    output: join(DATAPATH, 'results/figure4/figure4_paths.txt')
    shell:
        """
        touch {output}
        #echo 'Figure 4a {input.figure4a}' >> {output}
        #echo 'Figure 4b {input.figure4b}' >> {output}
        #echo 'Figure 4c-d {input.figure4cd}' >> {output}
        #echo 'Figure 4e {input.figure4a}' >> {output}
        #echo 'Figure 4f {input.figure4a}' >> {output}
        #echo 'Figure 4g {input.figure4g}' >> {output}
        echo 'Figure 4i-j {input.figure4ij}' >> {output}
        
        """

#================================================================================#
#         Primary versus Relapse differential gene expression analysis           #
#================================================================================#
rule fig4_primary_vs_relapse:
    input:
        NBexprs  = join(DATAPATH, 'data/tumor/rnaseq/exprs/tumor_RNAseq_Counts_Matrix.RDS'),
        rasSigr  = join(DATAPATH, 'db/publicGeneSigs/ras_target_genes.RDS'),
        NBreg    = join(DATAPATH, 'analysis/tumor/ARACNe/network.txt'),
        crcList  = join(DATAPATH, 'results/supptables/crcTF_fractionObserved.txt'),
        mesTFact = join(DATAPATH, "analysis/tumor/VIPER/MES_TFactivity.RDS")

    output:
        diffTab  = join(DATAPATH, 'analysis/tumor/Rel_vs_Pri/RelapseVsPrimary_topDiffExpGenes.txt'),
        mainFig1 = join(DATAPATH, 'results/figure4/JunFos_Ras_enrichment_barcodeplot_relapsevsPrimary.pdf'),
        mainFig2 = join(DATAPATH, 'results/figure4/crcTF_in_relapse_pri_enriched.pdf'),
        suppFig1 = join(DATAPATH, 'results/sup_figure4/all_TF_in_relapse_pri_enriched.pdf')

    params:
        script   = 'scripts/figure4/diffAnalysisPrimaryVSRelapse.R',
        outpath1 = join(DATAPATH, "analysis/"),
        outpath2 = join(DATAPATH, "results/")
    conda: '../envs/R3.5.yaml'
    shell:
        """
         Rscript {params.script} {input.NBexprs} {input.rasSigr} {input.NBreg} \
          {input.crcList} {input.mesTFact} {params.outpath1} {params.outpath2}
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
#                             TNFRSF12A related analysis                         #
#================================================================================#
rule fig4_TNFRSF12A_analysis:
    input:
        NBcells  = join(DATAPATH, 'analysis/cells/rnaseq/exprs/cells_RNAseq_TPM_Matrix_filt_log.RDS'),
        NBanno   = join(DATAPATH, 'annotation/annotation_cells.RDS'),
        NBexprs  = join(DATAPATH, 'analysis/tumor/rnaseq/exprs/tumor_RNAseq_TPM_Matrix_filt_log.RDS'),
        tumorNMF = join(DATAPATH, 'analysis/tumor/rnaseq/NMF/tumor_consensusSE_K4_Hmatrix_hnorm.RDS')
    output:
        mainFig1 = join(DATAPATH, 'results/figure4/TNFRSF12A_expression_NB.pdf')
    params:
        script  = 'scripts/figure4/TNFRSF12AexpNB.R',
        outpath = join(DATAPATH, 'results/')
    conda: '../envs/R3.5.yaml'

    shell:
        """
          Rscript {params.script} {input.NBcells} {input.NBanno} {input.NBexprs} {input.tumorNMF} {params.outpath}
        """


#================================================================================#
#        RAS JUN/FOS targets Expression mapped to Mouse GSE99933 E12.5           #
#================================================================================#
optK_tc = str(config['NMFparams']['tumor']['optimalK']['chipseq'])
rule fig4_RasJunFos_mouse:
    input:
        mouse_pstime = join(DATAPATH, 'db/GSE99933_E12.5/GSE99933_E12.5.txt'),
        mouse_exprs  = join(DATAPATH, 'db/GSE99933_E12.5/GSE99933_E12.5_exprs_Zscore.txt'),
        NBreg        = join(DATAPATH, 'analysis/tumor/ARACNe/network.txt'),
        rasSigr      = join(DATAPATH, 'db/publicGeneSigs/ras_target_genes.RDS')
    output:
        report = join(DATAPATH, 'reports/RasJunFos_Exprs_mouseE12.5.html'),
        rmd    = temp(join(DATAPATH, 'reports/RasJunFos_Exprs_mouseE12.5.Rmd')),
        figure = join(DATAPATH, 'results/figure4/RasJunFos_Exprs_mouseE12.5.pdf')
    params:
        script   = 'scripts/figure4/RasJunFos_Exprs_mouseE12.5.Rmd'
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  mouse_pstime = '{input.mouse_pstime}', \
                  mouse_exprs  = '{input.mouse_exprs}', \
                  NBreg        = '{input.NBreg}', \
                  rasSigr      = '{input.rasSigr}', \
                  figure       = '{output.figure}' \
                ))"


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
#                      Figure 4  - Example loci IGV                              #
#================================================================================#
optK_cc = str(config['NMFparams']['cells']['optimalK']['chipseq'])
rule fig4_IGV:
    input:
        consensusSE    = join(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS'),
        hsapiens_genes = join(DATAPATH, 'db/misc/EnsDb_Hsapiens_v75_genes.RDS'),
        cells_h        = join(DATAPATH, 'analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_cc + '_Hmatrix_hnorm.RDS')
    output:
        report = join(DATAPATH, 'reports/figure4_{ExamplegeneID}_IGV_plot.html'),
        rmd    = temp(join(DATAPATH, 'reports/figure4_{ExamplegeneID}_IGV_plot.Rmd')),
        figure = join(DATAPATH, 'results/figure4/figure4_{ExamplegeneID}_loci.pdf')
    params:
        script       = 'scripts/figure4/figure4_IGV_plot.Rmd',
        work_dir     = DATAPATH,
        width_window = config['igv_plot']['figure4']['window'],
        width        = config['igv_plot']['figure4']['width'],
        height       = config['igv_plot']['figure4']['height'],
        #gr_name  = config['igv_plot']['figure4i']['name']
        gr_name  = lambda wildcards: wildcards.ExamplegeneID
    conda: '../envs/R3.5.yaml'
    shell:
        """
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  work_dir       = '{params.work_dir}', \
                  SE             = '{input.consensusSE}', \
                  hsapiens_genes = '{input.hsapiens_genes}', \
                  cells_h        = '{input.cells_h}', \
                  width_window   = {params.width_window}, \
                  width          = {params.width}, \
                  height         = {params.height}, \
                  name           = '{params.gr_name}', \
                  figure         = '{output.figure}' \
                ))"


        """

