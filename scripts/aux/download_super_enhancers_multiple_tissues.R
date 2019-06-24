options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

url1 = as.character(args[1])
url2 = as.character(args[2])
outpath = as.character(args[3])

#-------------------------------------------------------------------------------------------------
# DATAPATH = "/icgc/dkfzlsdf/analysis/B080/crg/B087_Neuroblastoma/publication_GEO/"
# url1 = "https://www.cell.com/cms/10.1016/j.cell.2013.09.053/attachment/c44ace85-27a5-4f4f-b7e4-db375a76f583/mmc7.zip"
# url2 = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867413012270-mmc2.xlsx"
# outpath = paste0(DATAPATH,"db/")
#-------------------------------------------------------------------------------------------------

# This script downloads the super enhancers identified across multiple tissues and cells in the paper - 
# Hnisz, D., Abraham, B.J., Lee, T.I., Lau, A., Saint-André, V., Sigova, A.A., Hoke, H.A. and Young, R.A., 2013. 
# Super-enhancers in the control of cell identity and disease. Cell, 155(4), pp.934-947.

library(utils)

#Download the super enahancers and uzip
download.file(url = url1, destfile = paste0(outpath, "SEmultiTisuues/mmc7.zip"))
unzip(paste0(outpath, "SEmultiTisuues/mmc7.zip"), exdir = paste0(outpath,"SEmultiTisuues/"))
unlink(exdir = paste0(outpath,"SEmultiTisuues/mmc7.zip"))

# Download the super enhancers description file
download.file(url2, destfile = paste0(outpath,"SEmultiTisuues/mmc2.xlsx"))
