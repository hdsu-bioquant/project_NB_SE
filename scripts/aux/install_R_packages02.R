##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          Extra requirements                                ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

if(!require(wCorr)){
  install.packages('wCorr', repos = 'http://cran.us.r-project.org')
}
library(wCorr)

if(!require(EnsDb.Hsapiens.v75)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("EnsDb.Hsapiens.v75", suppressUpdates = T)  
}
library(EnsDb.Hsapiens.v75)  

if(!require(InteractionSet)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("InteractionSet", suppressUpdates = T)  
}
library(InteractionSet)  

if(!require(DT)){
  install.packages('DT', repos = 'http://cran.us.r-project.org')
}
library(DT)

if(!require(UpSetR)){
  install.packages('UpSetR', repos = 'http://cran.us.r-project.org')
}
library(UpSetR)

if(!require(eulerr)){
  install.packages('eulerr', repos = 'http://cran.us.r-project.org')
}
library(eulerr)