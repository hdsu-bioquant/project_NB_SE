options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(TRUE)

exppath = as.character(args[1])
regpath = as.character(args[2])
NMFexpo = as.character(args[3])
outpath = as.character(args[4])

#print(selsig)

library(viper)
library(parallel)

# NB normalized expression matrix
NBexp = read.table(exppath, header=T, stringsAsFactors=F, sep="\t", row.names = 1, check.names = FALSE)
NBexp = as.matrix(NBexp)

# Regulome object
regulon = aracne2regulon(afile = regpath, eset=NBexp, format="3col", verbose=T)

# Signature matrix from NMF
sig = readRDS(NMFexpo)

# Names of signatures
selsig = c(rownames(sig),"all")

if(identical(colnames(NBexp), colnames(sig)))
{
    # Conditions
    type = apply(sig, 2, function(x){rownames(sig)[which.max(x)]})

    lapply(selsig, function(x){

      if(x != "all")
      {
        #--------------------
        # Multi sample viper
        #--------------------

        # Assignment of signature groups
        sigN = which(type == x)
        rest = which(type != x)

        # Signature object
        signature = bootstrapTtest(NBexp[,sigN], NBexp[,rest], per=1000, seed=123, verbose=T, cores=40)
        print("Step 1")
        
        # Null model object
        nullmodel = ttestNull(NBexp[,sigN], NBexp[,rest], per=1000, repos=T, seed=123, verbose=T, cores=40)
        print("Step 2")
        
        # msVIPER analysis
        PerSigTFactivity = msviper(signature, regulon, nullmodel, pleiotropy = T, synergy = 0, cores=40, verbose = T)
        print("Step 3")
        
        # Compiling bootstrapped results
        PerSigTFactivity = bootstrapmsviper(PerSigTFactivity, "median")
        print("Step 4")
        
        # Leading edge analysis
        PerSigTFactivity = ledge(PerSigTFactivity)
        print("Step 5")
        
        # Saving results
        save(list=c("PerSigTFactivity","regulon", "nullmodel", "signature"), file=paste0(outpath, "/", x, "_viperAnalysis.Rdata"))
        saveRDS(PerSigTFactivity, file=paste0(outpath,"/", x, "_TFactivity.RDS"))
        print("Step save done")
        
      }else{

        #-----------------------
        # Sample specific Viper
        #-----------------------

        PerSampleTFactivity = viper(NBexp, regulon, pleiotropy = F, nes = TRUE,  method = "scale", bootstraps = 1000, cores = 10, verbose = TRUE)
        saveRDS(PerSampleTFactivity, file=paste0(outpath,"/NBpersampleTFactivity.RDS"))
      }
      return(NULL)
    })
 }
