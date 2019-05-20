#==============================================================================#
#            Analysis and figures included in manuscript figure 2              #
#==============================================================================#
#workdir: "../../"
#workdir: "~/"

rule compile_figure4:
    input:
        #figure4a = join(DATAPATH, 'results/figure2/figure2a_tumor_SE_hmatrix.pdf'),
        #figure4b = join(DATAPATH, 'results/figure2/figure2b_cells_SE_hmatrix.pdf'),
        #figure4c = join(DATAPATH, 'results/figure2/figure2c_tumor_cells_SE_UMAP.pdf'),
        #figure4d = join(DATAPATH, 'results/figure2/figure2d_01_tumor_riverplot.pdf'),
        #figure4e = join(DATAPATH, 'results/figure2/figure2e_tumor_SE_targets_hmatrix.pdf'),
        #figure4f = join(DATAPATH, 'results/figure2/figure2f_tumor_cells_density.pdf'),
        #figure4g = join(DATAPATH, 'results/figure2/figure2g_tumor_SE_targets_NMF_recovery.pdf')
        figure4i = join(DATAPATH, 'results/figure4/figure4i_IGV_plot.pdf')
    output: join(DATAPATH, 'results/figure4/figure4_paths.txt')
    shell:
        """
        touch {output}
        #echo 'Figure 4a {input.figure4i}' >> {output}
        #echo 'Figure 4b {input.figure4i}' >> {output}
        #echo 'Figure 4c {input.figure4i}' >> {output}
        #echo 'Figure 4d {input.figure4i}' >> {output}
        #echo 'Figure 4e {input.figure4i}' >> {output}
        #echo 'Figure 4f {input.figure4i}' >> {output}
        #echo 'Figure 4g {input.figure4i}' >> {output}
        echo 'Figure 4i {input.figure4i}' >> {output}
        
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
        path_config  = join(os.getcwd(), 'scripts/figure4/primary_metastasis_relapse_summarized_dropBadSample.txt'),
        width_window = 1000,
        ymax     = 35,
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

