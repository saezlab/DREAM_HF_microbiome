# Read the data
PARAM <- list()
PARAM$folder.R <- paste0(getwd(), "/")
PARAM$folder <- gsub("Rmd/", "", PARAM$folder.R)

# With event labels; for parameter optimization
S1 <- read.csv(file = paste0(PARAM$folder, "data/train/phenodata.csv"), row.names=1)
O1 <- read.csv(file = paste0(PARAM$folder, "data/train/readcounts.csv"), row.names=1)
T1 <- read.csv(file = paste0(PARAM$folder, "data/train/taxtable.csv"))

# Without event labels; for final scoring
S2 <- read.csv(file = paste0(PARAM$folder, "data/test/phenodata.csv"), row.names=1)
O2 <- read.csv(file = paste0(PARAM$folder, "data/test/readcounts.csv"), row.names=1)
T2 <- read.csv(file = paste0(PARAM$folder, "data/test/taxtable.csv"))



