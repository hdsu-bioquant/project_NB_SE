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
    #Consensus SE
    if config["consensusSE"]["call_consensus_tumor_SE"]:
        collectfiles.append(join(DATAPATH, 'analysis/tumor/chipseq/H3K27ac/consensusSE/tumor_H3K27ac_noH3K4me3_consensusSE.bed'))
        collectfiles.extend(expand(join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/consensusSE/{type}_H3K27ac_noH3K4me3_consensusSE_SignalScore.txt'), zip, type = ["tumor", "cells"]))
        collectfiles.extend(expand(join(DATAPATH, 'reports/01_{type}_chipseq_NMF_report.html'), zip, type = ["tumor", "cells"]))
        #collectfiles.append(join(DATAPATH, 'reports/01_{type}_{omics}_NMF_report.html'))
        collectfiles.append('.snakemake/completeLibrary.txt')
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



rule NMF_report_chipseq:
    input:
        matrix     = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/consensusSE/{type}_H3K27ac_noH3K4me3_consensusSE_SignalScore.RDS'),
        annotation = join(DATAPATH, 'annotation/annotation_{type}.RDS')
    output:
        report    = join(DATAPATH, 'reports/01_{type}_chipseq_NMF_report.html'),
        rmd       = temp(join(DATAPATH, 'reports/01_{type}_chipseq_NMF_report.Rmd')),
        nmf       = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/NMF/{type}_consensusSE_SignalScore_NMF.RDS'),
        norm_nmfW = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/NMF/{type}_consensusSE_SignalScore_normNMF_W.RDS'),
        norm_nmfH = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/NMF/{type}_consensusSE_SignalScore_normNMF_H.RDS')
    params:
        script   = 'scripts/analysis/99_NMF_report.Rmd',
        assayID  = '{type}_chipseq',
        workdir  = join(DATAPATH, 'analysis/{type}/chipseq/H3K27ac/NMF/'),
        nmf_kmin = lambda wildcards: config['NMFparams'][wildcards.type]['k.min'],
        nmf_kmax = lambda wildcards: config['NMFparams'][wildcards.type]['k.max'],
        nmf_iter = lambda wildcards: config['NMFparams'][wildcards.type]['iterations']
    conda: 'envs/cuda_R3.5.yaml'
    shell:
        """
        unset LD_LIBRARY_PATH
        export PATH="/usr/local/cuda/bin:$PATH"
        nvidia-smi
    
        cp {params.script} {output.rmd}

        Rscript -e "rmarkdown::render( '{output.rmd}', \
                params = list( \
                  assayID   = '{params.assayID}', \
                  work_dir  = '{params.workdir}', \
                  nmf_kmin  = '{params.nmf_kmin}', \
                  nmf_kmax  = '{params.nmf_kmax}', \
                  nmf_iter  = '{params.nmf_iter}', \
                  nmf       = '{output.nmf}' \
                  norm_nmfW = '{output.norm_nmfW}' \
                  norm_nmfH = '{output.norm_nmfH}' \
                  matrix    = '{input.matrix}', \
                  metadata  = '{input.annotation}' \
                ))"
        
        
        """



#================================================================================#
#                     SE signal bigWig AVERAGE OVER BED                          #
#================================================================================#
### Computes the SE average score signal for tumors
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
        script='scripts/analysis/01_SEmatrix.R',
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



# Install missing R packages in conda env cuda_R3.4
rule install_missing_R_01:
    output: '.snakemake/completeLibrary.txt'
    params:
        script  = 'scripts/aux/install_R_packages01.R',
    conda: 'envs/cuda_R3.5.yaml'
    shell:
        """
        
        unset LD_LIBRARY_PATH
        export PATH="/usr/local/cuda/bin:$PATH"
        nvidia-smi
        conda install openssl=1.0
        
        Rscript {params.script}
        git clone https://github.com/cudamat/cudamat.git
        pip install cudamat/
        rm -rf cudamat
    
        touch {output}
        
        """
