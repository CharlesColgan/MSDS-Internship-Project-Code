#function for installation and library initialization
function(){
  #Install nessecary Packages
  Base=c("tidyverse",
         "metap",
         "quantreg",
         "statmod",
         "mvtnorm",
         "reshape2",
         "mice")
  for(i in 1:length(Base)){
    if(!requireNamespace(Base[i],quietly=TRUE)){
      install.packages(Base[i])
    }
  }
  
  if(!requireNamespace("BiocManager",quietly=TRUE)){
    install.packages("BiocManager")
    }
  BioC=c("topGO",
         "limma",
         "pcaMethods",
         "impute",
         "UniProt.ws",
         "multtest",
         "Biostrings",
         "org.Hs.eg.db",
         "protViz")
  for(i in 1:length(BioC)){
    if(!requireNamespace(BioC[i],quietly=TRUE)){
      BiocManager::install(BioC[i])
      }
  }

  if(!require(NAguideR)){
    devtools::install_github("wangshisheng/NAguideR")
  }
  
  if (!requireNamespace("reticulate",quietly=TRUE)){
    install.packages("reticulate")
  }
  
  #initialise libarries
  library(tidyverse)
  library(BiocManager)
  library(devtools)
  library(metap)
  library(mice)
  library(statmod)
  library(mvtnorm)
  library(quantreg)
  library(reshape2)
  library(MSstats)
  library(topGO)
  library(limma)
  library(Biostrings)
  library(BiocVersion)
  library(pcaMethods)
  library(NAguideR)
  library(impute)
  library(protViz)
  library(topGO)
  library(UniProt.ws)
  library(proteusLabelFree)
  library(proteusTMT)
  library(proteusSILAC)
  library(reticulate)
  library(readxl)
  library(readr)
  
  #cleaning Environmetn
  rm(Base,BioC,i)
}
