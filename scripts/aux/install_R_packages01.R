##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            Bratwurst requirements                          ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

# if(!require(kableExtra)){
#   install.packages('kableExtra', repos = 'http://cran.us.r-project.org')
# }
#
# if(!require(RcppCNPy)){
#   install.packages('RcppCNPy', repos = 'http://cran.us.r-project.org')
# }

if(!require(riverplot)){
  install.packages('riverplot', repos = 'http://cran.us.r-project.org')
}


library(devtools)
options(unzip = "internal")
devtools::install_github("wurst-theke/bratwurst")


library(Bratwurst)


