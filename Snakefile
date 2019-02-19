#==============================================================================#
#                       Parse config and sample information                    #
#==============================================================================#
import csv

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

