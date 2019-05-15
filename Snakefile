#==============================================================================#
#                       Parse config and sample information                    #
#==============================================================================#
# IMPORT python libraries
from os.path import join
import csv
import re

# Import config file & parameters
configfile: 'config.yaml'

# Read Annotation CSV and find samples with ChIPseq or/and RNAseq data
TUMOR_SAMPLES_CHIP = []
TUMOR_SAMPLES_RNA  = []
TUMOR_SAMPLES_CHIP_RNA = []
with open(config['tumor_annotation_csv']) as f:
    reader = csv.DictReader(f, delimiter=',')
    for row in reader:
        if row['avail.ChIPseq'] == 'TRUE':
            TUMOR_SAMPLES_CHIP.append(row['ProjectID'])
        if row['avail.RNAseq'] == 'TRUE':
            TUMOR_SAMPLES_RNA.append(row['ProjectID'])
            if row['avail.ChIPseq'] == 'TRUE':
                TUMOR_SAMPLES_CHIP_RNA.append(row['ProjectID'])

CELLS_SAMPLES_CHIP = []
CELLS_SAMPLES_RNA  = []
CELLS_SAMPLES_CHIP_RNA = []
with open(config['cells_annotation_csv']) as f:
    reader = csv.DictReader(f, delimiter=',')
    for row in reader:
        if row['avail.ChIPseq'] == 'TRUE':
            CELLS_SAMPLES_CHIP.append(row['ProjectID'])
        if row['avail.RNAseq'] == 'TRUE':
            CELLS_SAMPLES_RNA.append(row['ProjectID'])
            if row['avail.ChIPseq'] == 'TRUE':
                CELLS_SAMPLES_CHIP_RNA.append(row['ProjectID'])


# TUMOR_SAMPLES_CHIP
# TUMOR_SAMPLES_RNA
# TUMOR_SAMPLES_CHIP_RNA
# len(TUMOR_SAMPLES_CHIP)
# len(TUMOR_SAMPLES_RNA)
# len(TUMOR_SAMPLES_CHIP_RNA)


#==============================================================================#
#                  Main path to store results and tmp data                     #
#==============================================================================#
# Import paths from config file
DATAPATH = config['main_working_directory']


#==============================================================================#
#                               Include extra rules                            #
#==============================================================================#
# Include snakefiles containing figure rules
include: "snakefiles/figure2.Snakefile"


        
#==============================================================================#
#               Print sample data at the pipeline's start.                     #
#==============================================================================#
    
def printExp():
  print("-------------------------------")
  print(str(len(TUMOR_SAMPLES_CHIP)) + " Tumor samples with ChIPseq data available")
  #print(TUMOR_SAMPLES_CHIP)
  print("-------------------------------")
  print(str(len(TUMOR_SAMPLES_RNA)) + " Tumor samples with RNAseq data available ")
  #print(TUMOR_SAMPLES_RNA)
  print("-------------------------------")
  print(str(len(TUMOR_SAMPLES_CHIP_RNA)) + " Tumor samples with ChIPseq and RNAseq data available")
  #print(TUMOR_SAMPLES_CHIP_RNA)
  print("")
  print("-------------------------------")
  print(str(len(CELLS_SAMPLES_CHIP)) + " Cell lines samples with ChIPseq data available")
  #print(CELLS_SAMPLES_CHIP)
  print("-------------------------------")
  print(str(len(CELLS_SAMPLES_RNA)) + " Cell lines samples with RNAseq data available ")
  #print(CELLS_SAMPLES_RNA)
  print("-------------------------------")
  print(str(len(CELLS_SAMPLES_CHIP_RNA)) + " Cell lines samples with ChIPseq and RNAseq data available")
  #print(CELLS_SAMPLES_CHIP_RNA)
  print("")
printExp()



#==============================================================================#
#                         Function to collect final files                      #
#==============================================================================#
#helper function to collect final files from pipeline
def inputall(wilcards):
    collectfiles = []
    if config["compileFigs"]["figure2"]:
        collectfiles.append(join(DATAPATH, 'results/figures/figure2/figure2_paths.txt'))
    if config["phase02_NMF"]["NMF_chipseq_tumor"]:
        collectfiles.append(join(DATAPATH, 'reports/03_tumor_chipseq_NMF_report.html'))
    if config["phase02_NMF"]["NMF_chipseq_cells"]:
        collectfiles.append(join(DATAPATH, 'reports/03_cells_chipseq_NMF_report.html'))
        #collectfiles.append(join(DATAPATH, 'reports/01_{type}_{omics}_NMF_report.html'))
    if config["phase01_consensusSE"]["SE_target_gene"]:
        collectfiles.append(join(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS'))
    #Consensus SE
    if config["phase01_consensusSE"]["consensus_tumor_SE"]:
        collectfiles.append(join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed'))
        collectfiles.append(join(DATAPATH, 'analysis/tumor_cells/chipseq/H3K27ac/consensusSE/tumor_cells_H3K27ac_noH3K4me3_consensusSE_SignalScore.RDS'))
        collectfiles.extend(expand(join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/consensusSE/{type}_H3K27ac_noH3K4me3_consensusSE_SignalScore.txt'), zip, type = ["tumor", "cells"]))
        #collectfiles.append('.snakemake/completeLibrary.txt')
        #collectfiles.append(join(DATAPATH, 'tmp.txt'))
    #return final list of all files to collect from the pipeline
    return collectfiles

# Collect pipeline result files
rule all:
    input: inputall


#"envs/R3.5.yaml"

rule placeh:
    input:
        consensusSE = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed')
    output: 
        outtmp = join(DATAPATH, 'tmp.txt')
    params:
        script='scripts/analysis/01_SEmatrix.R',
    conda:
        'envs/cuda_R3.5.yaml'
    shell:
        """
        Rscript {params.script} {output.outtmp} {input.consensusSE} 
        """


optK_tcc = str(config['NMFparams']['tumor_cells']['optimalK']['chipseq'])
rule NMF_report_chipseq_tumor_cells:
    input:
        matrix     = join(DATAPATH, 'analysis/tumor_cells_/chipseq/H3K27ac/consensusSE/tumor_cells__H3K27ac_noH3K4me3_consensusSE_SignalScore.RDS'),
        annotation = join(DATAPATH, 'annotation/annotation_{type}.RDS')
    output:
        report    = join(DATAPATH, 'reports/03_tumor_cells__chipseq_NMF_report.html'),
        rmd       = temp(join(DATAPATH, 'reports/03_tumor_cells_chipseq_NMF_report.Rmd')),
        nmf       = join(DATAPATH, 'analysis/tumor_cells/chipseq/H3K27ac/NMF/tumor_cells_consensusSE_SignalScore_NMF.RDS'),
        norm_nmfW = join(DATAPATH, 'analysis/tumor_cells/chipseq/H3K27ac/NMF/tumor_cells_consensusSE_SignalScore_normNMF_W.RDS'),
        norm_nmfH = join(DATAPATH, 'analysis/tumor_cells/chipseq/H3K27ac/NMF/tumor_cells_consensusSE_SignalScore_normNMF_H.RDS'),
        hmatrix_wnorm = join(DATAPATH, ('analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_tcc + '_Hmatrix_wnorm.RDS')),
        wmatrix_wnorm = join(DATAPATH, ('analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_tcc + '_Wmatrix_Wnorm.RDS')),
        nmf_features  = join(DATAPATH, ('analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_tcc + '_NMF_features.RDS')),
        hmatrix_hnorm = join(DATAPATH, ('analysis/cells/chipseq/H3K27ac/NMF/cells_consensusSE_K' + optK_tcc + '_Hmatrix_hnorm.RDS'))
    params:
        script   = 'scripts/analysis/03_tumor_cells_chipseq_NMF_report.Rmd',
        assayID  = '{type}_chipseq',
        workdir  = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/NMF/'),
        nmf_kmin = lambda wildcards: config['NMFparams'][wildcards.type]['k.min'],
        nmf_kmax = lambda wildcards: config['NMFparams'][wildcards.type]['k.max'],
        nmf_iter = lambda wildcards: config['NMFparams'][wildcards.type]['iterations']
    conda: 'envs/R3.5.yaml'
    shell:
        """
    
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  assayID   = '{params.assayID}', \
                  work_dir  = '{params.workdir}', \
                  nmf_kmin  = '{params.nmf_kmin}', \
                  nmf_kmax  = '{params.nmf_kmax}', \
                  nmf_iter  = '{params.nmf_iter}', \
                  nmf       = '{output.nmf}', \
                  norm_nmfW = '{output.norm_nmfW}', \
                  norm_nmfH = '{output.norm_nmfH}', \
                  matrix    = '{input.matrix}', \
                  metadata  = '{input.annotation}' \
                  hmatrix_wnorm = '{output.hmatrix_wnorm}', \
                  wmatrix_wnorm = '{output.wmatrix_wnorm}', \
                  nmf_features  = '{output.nmf_features}', \
                  hmatrix_hnorm = '{output.hmatrix_hnorm}', \
                ))"
        
        
        """


rule NMF_report_chipseq:
    input:
        matrix     = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/consensusSE/{type}_H3K27ac_noH3K4me3_consensusSE_SignalScore.RDS'),
        annotation = join(DATAPATH, 'annotation/annotation_{type}.RDS')
    output:
        report    = join(DATAPATH, 'reports/03_{type}_chipseq_NMF_report.html'),
        rmd       = temp(join(DATAPATH, 'reports/03_{type}_chipseq_NMF_report.Rmd')),
        nmf       = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/NMF/{type}_consensusSE_SignalScore_NMF.RDS'),
        norm_nmfW = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/NMF/{type}_consensusSE_SignalScore_normNMF_W.RDS'),
        norm_nmfH = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/NMF/{type}_consensusSE_SignalScore_normNMF_H.RDS')
    params:
        script   = 'scripts/analysis/03_chipseq_NMF_report.Rmd',
        assayID  = '{type}_chipseq',
        workdir  = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/NMF/'),
        nmf_kmin = lambda wildcards: config['NMFparams'][wildcards.type]['k.min'],
        nmf_kmax = lambda wildcards: config['NMFparams'][wildcards.type]['k.max'],
        nmf_iter = lambda wildcards: config['NMFparams'][wildcards.type]['iterations']
    conda: 'envs/R3.5.yaml'
    shell:
        """
        #unset LD_LIBRARY_PATH
        #export PATH="/usr/local/cuda/bin:$PATH"
        #nvidia-smi
    
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  assayID   = '{params.assayID}', \
                  work_dir  = '{params.workdir}', \
                  nmf_kmin  = '{params.nmf_kmin}', \
                  nmf_kmax  = '{params.nmf_kmax}', \
                  nmf_iter  = '{params.nmf_iter}', \
                  nmf       = '{output.nmf}', \
                  norm_nmfW = '{output.norm_nmfW}', \
                  norm_nmfH = '{output.norm_nmfH}', \
                  matrix    = '{input.matrix}', \
                  metadata  = '{input.annotation}' \
                ))"
        
        
        """


#================================================================================#
#                                  SE target genes                               #
#================================================================================#
rule SE_target_genes:
    input:
        tumor_annot = join(DATAPATH, 'annotation/annotation_tumor.RDS'),
        SE_bed      = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed'),
        hichip_SK_N_AS = join(DATAPATH, 'data/cells/hichip/mango/SK-N-AS_HiChIP_mango.all'),
        hichip_CLB_GA = join(DATAPATH, 'data/cells/hichip/mango/CLB-GA_HiChIP_mango.all'),
        SE_signal = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE_SignalScore.RDS'),
        gene_exprs = join(DATAPATH, 'data/tumor/rnaseq/exprs/tumor_RNAseq_Counts_Matrix.RDS'),
        hic = join(DATAPATH, 'db/hic/GSE63525_K562_HiCCUPS_looplist.txt'),
        TADs = join(DATAPATH, 'db/TADs/hESC_domains_hg19.RDS'),
        hsapiens_genes = join(DATAPATH, 'db/misc/EnsDb_Hsapiens_v75_genes.RDS')
    output:
        report    = join(DATAPATH, 'reports/02_SE_target_genes_report.html'),
        rmd       = temp(join(DATAPATH, 'reports/02_SE_target_genes_report.Rmd')),
        SE_target = join(DATAPATH, 'analysis/tumor/SE_annot/tumor_consensusSE_target_GRanges.RDS')
    params:
        script   = 'scripts/analysis/02_SE_target_genes.Rmd',
        workdir  = DATAPATH
    conda: 'envs/R3.5.yaml'
    shell:
        """
        #unset LD_LIBRARY_PATH
        #export PATH="/usr/local/cuda/bin:$PATH"
        #nvidia-smi
    
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
            params = list( \
                work_dir       = '{params.workdir}', \
                tumor_annot    = '{input.tumor_annot}', \
                SE             = '{input.SE_bed}', \
                hichip_SK_N_AS = '{input.hichip_SK_N_AS}', \
                hichip_CLB_GA  = '{input.hichip_CLB_GA}', \
                SE_signal      = '{input.SE_signal}', \
                gene_exprs     = '{input.gene_exprs}', \
                hic            = '{input.hic}', \
                TADs           = '{input.TADs}', \
                hsapiens_genes = '{input.hsapiens_genes}', \
                SE_target_gr   = '{output.SE_target}' \
                ))"
        
        
        """


#================================================================================#
#                     SE signal bigWig AVERAGE OVER BED                          #
#================================================================================#
### Computes the SE average score signal for tumors and cells
rule SE_SignalMatrix_combined:
    input:
        averageOverBed_tumor = expand((DATAPATH + 'analysis/tumor/chipseq/H3K27ac/consensusSE/{sample}_H3K27ac_bigWigAverageOverBed.txt'), zip, sample = TUMOR_SAMPLES_CHIP),
        averageOverBed_cells = expand((DATAPATH + 'analysis/cells/chipseq/H3K27ac/consensusSE/{sample}_H3K27ac_bigWigAverageOverBed.txt'), zip, sample = CELLS_SAMPLES_CHIP),
        consensusSE          = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed')
    output: 
        matrix_rds = join(DATAPATH, 'analysis/tumor_cells/chipseq/H3K27ac/consensusSE/tumor_cells_H3K27ac_noH3K4me3_consensusSE_SignalScore.RDS'),
        matrix_txt = join(DATAPATH, 'analysis/tumor_cells/chipseq/H3K27ac/consensusSE/tumor_cells_H3K27ac_noH3K4me3_consensusSE_SignalScore.txt')
    params:
        script='scripts/analysis/01_SEmatrix.R'
    conda:
        'envs/R3.5.yaml'
    shell:
        """
        Rscript {params.script} {output.matrix_rds} {output.matrix_txt} {input.consensusSE} {input.averageOverBed_tumor} {input.averageOverBed_cells}
        """


### Computes the SE average score signal for tumors or cells
def find_bwAverage(wildcards):
    SAMPLES = TUMOR_SAMPLES_CHIP if wildcards.type == "tumor" else CELLS_SAMPLES_CHIP
    averageOverBed = expand((DATAPATH + 'analysis/' + wildcards.type + '/chipseq/H3K27ac/consensusSE/{sample}_H3K27ac_bigWigAverageOverBed.txt'), zip, sample = SAMPLES)
    #print(wildcards.type)
    return averageOverBed
    
rule SE_SignalMatrix:
    input:
        averageOverBed_path = find_bwAverage,
        consensusSE         = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed')
    output: 
        matrix_rds = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/consensusSE/{type}_H3K27ac_noH3K4me3_consensusSE_SignalScore.RDS'),
        matrix_txt = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/consensusSE/{type}_H3K27ac_noH3K4me3_consensusSE_SignalScore.txt')
    params:
        script='scripts/analysis/01_SEmatrix.R'
    wildcard_constraints:
        type = "a-z"
    conda:
        'envs/R3.5.yaml'
    shell:
        """
        Rscript {params.script} {output.matrix_rds} {output.matrix_txt} {input.consensusSE} {input.averageOverBed_path}
        """


### Computes the average score over each bed for tumors
rule SE_bigwigaverageoverbed:
    input:
        bw = join(DATAPATH, 'data/{type}/chipseq/H3K27ac/bw/{sample}_H3K27ac.bw'), 
        consensusSE = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed')
    output:
        bw_over_bed=temp(join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/consensusSE/{sample}_H3K27ac_bigWigAverageOverBed.txt'))
    conda:
        'envs/generaltools.yaml'
    shell:
        """
        # Compute the average score of the SES_substract.bw bigWig over the noH3K4me3 consensus SE
        bigWigAverageOverBed {input.bw} {input.consensusSE} {output.bw_over_bed}
        """



#================================================================================#
#                TUMORS CONSENSUS SUPER ENHANCERS  FILTER H3K4me3                #
#================================================================================#
### Compute consensus SE list from SE called by rose for each sample ater H3K4me3 filtering
rule tumors_consensus_SE_noH3K4me3:
    input:
        seH3K27ac_noH3K4me3 = expand(join(DATAPATH, 'data/tumor/chipseq/H3K27ac/SE/{sample}_H3K27ac_ROSE_noH3K4me3_SuperEnhancers.bed'), zip, sample=TUMOR_SAMPLES_CHIP)
    output:
        consensusbed = join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed')
    conda:
        'envs/generaltools.yaml'
    shell:
        """
        # Merge all SE
        cat {input.seH3K27ac_noH3K4me3}| sortBed | bedtools merge -c 4,4 -o distinct,count_distinct | 
        awk '$5 > 1' |sed -e 's/$/\tSE_/' | sed -n 'p;=' |
        paste -d "" - - | awk 'BEGIN{{FS="\t";OFS="\t"}} {{ t = $4; $4 = $6; $6 = t; print;}} ' > {output.consensusbed}
        """



# Download auxiliary data and install missing R packages in conda env
rule down_misc_install_missing_R:
    output:
        hsapiens_genes = join(DATAPATH, 'db/misc/EnsDb_Hsapiens_v75_genes.RDS'),
        hic_K562       = join(DATAPATH, 'db/hic/GSE63525_K562_HiCCUPS_looplist.txt'),
        tads           = join(DATAPATH, 'db/TADs/hESC_domains_hg19.RDS')
    params:
        script  = 'scripts/aux/missing_packages_and_aux_data.R',
        workdir = DATAPATH
    conda: 'envs/R3.5.yaml'
    shell:
        """
        
        #unset LD_LIBRARY_PATH
        #export PATH="/usr/local/cuda/bin:$PATH"
        #nvidia-smi
        #conda install openssl=1.0
        
        Rscript {params.script} {params.workdir} {output.hsapiens_genes} {output.hic_K562} {output.tads}
        #git clone https://github.com/cudamat/cudamat.git
        #pip install cudamat/
        #rm -rf cudamat
        
        """


