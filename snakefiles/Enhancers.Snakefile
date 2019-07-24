#================================================================================#
#                            Enhancers target genes                              #
#================================================================================#
# # Finds enhancers target gene, also saves RNAseq expression matrix (TPMs) for tumor and cells
# rule enhancers_target_genes:
#     input:
#         tumor_annot    = join(DATAPATH, 'annotation/annotation_tumor.RDS'),
#         SE_bed         = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed'),
#         hichip_SK_N_AS = join(DATAPATH, 'data/cells/hichip/mango/SK-N-AS_HiChIP_mango.all'),
#         hichip_CLB_GA  = join(DATAPATH, 'data/cells/hichip/mango/CLB-GA_HiChIP_mango.all'),
#         SE_signal      = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE_SignalScore.RDS'),
#         gene_exprs     = join(DATAPATH, 'data/tumor/rnaseq/exprs/tumor_RNAseq_TPM_Matrix.RDS'),
#         gene_exprs_cl  = join(DATAPATH, 'data/cells/rnaseq/exprs/cells_RNAseq_TPM_Matrix.RDS'),
#         hic            = join(DATAPATH, 'db/hic/GSE63525_K562_HiCCUPS_looplist.txt'),
#         TADs           = join(DATAPATH, 'db/TADs/hESC_domains_hg19.RDS'),
#         hsapiens_genes = join(DATAPATH, 'db/misc/EnsDb_Hsapiens_v75_genes.RDS')
#     output:
#         report          = join(DATAPATH, 'reports/02_SE_target_genes_report.html'),
#         rmd             = temp(join(DATAPATH, 'reports/02_SE_target_genes_report.Rmd')),
#         SE_target       = join(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS'),
#         SE_target_df    = join(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_annotation_df.RDS'),
#         tumor_exprs_fil = join(DATAPATH, 'analysis/tumor/rnaseq/exprs/tumor_RNAseq_TPM_Matrix_filt_log.RDS'),
#         cells_exprs_fil = join(DATAPATH, 'analysis/cells/rnaseq/exprs/cells_RNAseq_TPM_Matrix_filt_log.RDS')
#     params:
#         script   = 'scripts/analysis/02_SE_target_genes.Rmd',
#         workdir  = DATAPATH
#     conda: 'envs/R3.5.yaml'
#     shell:
#         """
#         #unset LD_LIBRARY_PATH
#         #export PATH="/usr/local/cuda/bin:$PATH"
#         #nvidia-smi
#
#         cp {params.script} {output.rmd}
#
#         Rscript -e "rmarkdown::render( '{output.rmd}', \
#             params = list( \
#                 work_dir       = '{params.workdir}', \
#                 tumor_annot    = '{input.tumor_annot}', \
#                 SE             = '{input.SE_bed}', \
#                 hichip_SK_N_AS = '{input.hichip_SK_N_AS}', \
#                 hichip_CLB_GA  = '{input.hichip_CLB_GA}', \
#                 SE_signal      = '{input.SE_signal}', \
#                 gene_exprs     = '{input.gene_exprs}', \
#                 hic            = '{input.hic}', \
#                 TADs           = '{input.TADs}', \
#                 hsapiens_genes = '{input.hsapiens_genes}', \
#                 SE_target_gr   = '{output.SE_target}', \
#                 gene_exprs_cells     = '{input.gene_exprs_cl}', \
#                 tumor_exprs_filtered = '{output.tumor_exprs_fil}', \
#                 cells_exprs_filtered = '{output.cells_exprs_fil}' \
#                 ))"
#
#
#         """


#================================================================================#
#                  Enhancers signal bigWig AVERAGE OVER BED                      #
#================================================================================#
### Computes the SE average score signal for tumors and cells
rule enhancers_SignalMatrix_combined:
    input:
        averageOverBed_tumor = expand((DATAPATH + 'analysis/tumor/chipseq/H3K27ac/consensusEnhancers/{sample}_H3K27ac_bigWigAverageOverBed.txt'), zip, sample = TUMOR_SAMPLES_CHIP),
        averageOverBed_cells = expand((DATAPATH + 'analysis/cells/chipseq/H3K27ac/consensusEnhancers/{sample}_H3K27ac_bigWigAverageOverBed.txt'), zip, sample = CELLS_SAMPLES_CHIP),
        consensusEnhan       = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusEnhancers/tumor_H3K27ac_noH3K4me3_consensusEnhancers.bed')
    output: 
        matrix_rds = join(DATAPATH, 'analysis/tumor_cells/chipseq/H3K27ac/consensusEnhancers/tumor_cells_H3K27ac_noH3K4me3_consensusEnhancers_SignalScore.RDS'),
        matrix_txt = join(DATAPATH, 'analysis/tumor_cells/chipseq/H3K27ac/consensusEnhancers/tumor_cells_H3K27ac_noH3K4me3_consensusEnhancers_SignalScore.txt')
    params:
        script = 'scripts/analysis/01_SEmatrix.R'
    conda:
        '../envs/R3.5.yaml'
    shell:
        """
        Rscript {params.script} {output.matrix_rds} {output.matrix_txt} {input.consensusEnhan} {input.averageOverBed_tumor} {input.averageOverBed_cells}
        """


### Computes the SE average score signal for tumors or cells
def find_bwAverage(wildcards):
    SAMPLES = TUMOR_SAMPLES_CHIP if wildcards.type == "tumor" else CELLS_SAMPLES_CHIP
    averageOverBed = expand((DATAPATH + 'analysis/' + wildcards.type + '/chipseq/H3K27ac/consensusEnhancers/{sample}_H3K27ac_bigWigAverageOverBed.txt'), zip, sample = SAMPLES)
    #print(wildcards.type)
    #print(averageOverBed)
    return averageOverBed
    
rule enhancers_SignalMatrix:
    input:
        averageOverBed_path = find_bwAverage,
        consensusEnhan      = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusEnhancers/tumor_H3K27ac_noH3K4me3_consensusEnhancers.bed')
    output: 
        matrix_rds = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/consensusEnhancers/{type}_H3K27ac_noH3K4me3_consensusEnhancers_SignalScore.RDS'),
        matrix_txt = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/consensusEnhancers/{type}_H3K27ac_noH3K4me3_consensusEnhancers_SignalScore.txt')
    params:
        script = 'scripts/analysis/01_SEmatrix.R'
    wildcard_constraints:
        type = "[a-z]+"
    conda:
        '../envs/R3.5.yaml'
    shell:
        """
        Rscript {params.script} {output.matrix_rds} {output.matrix_txt} {input.consensusEnhan} {input.averageOverBed_path}
        """


### Computes the average score over each bed for tumors
rule enhancers_bigwigaverageoverbed:
    input:
        hsapiens_genes = join(DATAPATH, 'db/misc/EnsDb_Hsapiens_v75_genes.RDS'),
        bw             = join(DATAPATH, 'data/{type}/chipseq/H3K27ac/bw/{sample}_H3K27ac.bw'), 
        consensusEnhan = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusEnhancers/tumor_H3K27ac_noH3K4me3_consensusEnhancers.bed')
    output:
        bw_over_bed    = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/consensusEnhancers/{sample}_H3K27ac_bigWigAverageOverBed.txt')
    conda:
        '../envs/generaltools.yaml'
    shell:
        """
        # Compute the average score of the SES_substract.bw bigWig over the noH3K4me3 consensus SE
        bigWigAverageOverBed {input.bw} {input.consensusEnhan} {output.bw_over_bed}
        """


#================================================================================#
#                     TUMORS CONSENSUS ENHANCERS  FILTER H3K4me3                 #
#================================================================================#
### Compute consensus enhancer list from enhancers called by rose for each sample after H3K4me3 filtering
rule tumors_consensus_enhancers_noH3K4me3:
    input:
        enH3K27ac_noH3K4me3 = expand(join(DATAPATH, 'data/tumor/chipseq/H3K27ac/enhancers/{sample}_H3K27ac_ROSE_noH3K4me3_Enhancers.bed'), zip, sample=TUMOR_SAMPLES_CHIP)
    output:
        consensusbed        = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusEnhancers/tumor_H3K27ac_noH3K4me3_consensusEnhancers.bed')
    conda:
        '../envs/generaltools.yaml'
    shell:
        """
        # Merge all SE
        cat {input.enH3K27ac_noH3K4me3}| sortBed | bedtools merge -c 4,4 -o distinct,count_distinct | 
        awk '$5 > 1' |sed -e 's/$/\tEnhancer_/' | sed -n 'p;=' |
        paste -d "" - - | awk 'BEGIN{{FS="\t";OFS="\t"}} {{ t = $4; $4 = $6; $6 = t; print;}} ' > {output.consensusbed}
        """


