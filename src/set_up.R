## Packages
#CRAN
cran_packages <- c( "survival", "tidyverse", "magrittr", "Hmisc", "BiocManager",
                                 "smoothHR", "reshape2", "ggsci", "openxlsx", 
                                "randomForestSRC", "survivalROC")
# Load packages
for(p in cran_packages) {
  
  if(!require(p, character.only = TRUE)) {
    
    install.packages(p)
    
  }
  
}

#Bioconductor
bioc_packages <- c("SIAMCAT")

# Load packages
for(p in bioc_packages) {
  
  if(!require(p, character.only = TRUE)) {
    
    BiocManager::install(p)
    
  }
  
}




