# load libraries 
library(tidyverse)
library(ResourceSelection) # for hoslem test
set.seed(198657)

# create output folder
mainDir <- getwd()
subDir <- "output"

if (file.exists(subDir)){
  setwd(file.path(mainDir))
} else {
  dir.create(file.path(mainDir, subDir))
  setwd(file.path(mainDir))
}
# set parameters
PARAM <- list()
PARAM$folder.R <- paste0(getwd(), "/")
PARAM$folder.result <- paste0(PARAM$folder.R, "output/")
PARAM$folder.data <- paste0(PARAM$folder.R, "/")

# load data
print("Load data")
S.test <- read.csv(file = paste0(PARAM$folder.data, "test/pheno_test.csv"), row.names=1,)
S.train <- read.csv(file = paste0(PARAM$folder.data, "train/pheno_training.csv"), row.names=1,)
S.val <- read.csv(file = paste0(PARAM$folder.data, "validation/pheno_validation.csv"), row.names=1,)

O.test <- read.csv(file = paste0(PARAM$folder.data, "test/readcounts_test.csv"), row.names=1,)
O.train <- read.csv(file = paste0(PARAM$folder.data, "train/readcounts_training.csv"), row.names=1,)
O.val <- read.csv(file = paste0(PARAM$folder.data, "validation/readcounts_validation.csv"), row.names=1,)
# OPTIONAL btw line 31-59

# change rownames with fake ids for convenience
rows <- rownames(O.train)  
seq=seq(1:c(length(rows)))
fakename<- sapply(seq,
                  function(x) paste('bacteria', x))
fakename = gsub(" ","",as.character(fakename))
rownames(O.train) <- fakename

data <- cbind(meta(S.train), t(O.train))
data$Event_time<- as.numeric(data$Event_time)

# define potential predictors
predictors <- c(colnames(S.train), rownames(O.train))

###############################################################################
# Fix and combine train and test datasets
# change rownames with fake ids for convenience
rownames(O.test) <- fakename
rownames(O.val) <- fakename

df.test <- cbind(meta(S.test), t(O.test))
df.test$Event_time <- as.numeric(df.test$Event_time)
df.val <- cbind(meta(S.val), t(O.val))
df.val$Event_time <- as.numeric(df.val$Event_time)
rm(O.train, O.test, O.val, S.test, S.val, S.train,fakename )

df.total= rbind(data, df.val)
# choose only HF cases and eliminate "NA"
data <- subset(data, select=predictors, PrevalentHFAIL==0&!is.na(Event))

###############################################################################
# FOR VALIDATION PHASE ONLY ****************************************************
###############################################################################

# Read scores from participant
model=cox.mod.subset
scores=read.csv(file = paste0(PARAM$folder.result, "scores.csv")) 
collect.scores <- data.frame()
df=df.val
# HarrellC Index ***************************************************************
# Read the real labels that were hidden from that algorithm
labels <- df
labels$SampleID<- rownames(df)

# Align the user provided scores with the true event times
true.scores <- as.numeric(labels[scores$SampleID,"Event_time"])

# Calculate Harrell's C statistics
HarrellC <- Hmisc::rcorr.cens(scores$Score, true.scores, outx=FALSE)
print(HarrellC["C Index"])
collect.scores <- rbind(collect.scores, HarrellC["C Index"])

# HOSLEM ***********************************************************************
test <- hoslem.test(true.scores, scores$Score) # x numeric observations, y expected values.
p.val <- test$p.value
collect.scores <- cbind(collect.scores, p.val)
colnames(collect.scores) <- c("HarrellsC","Hoslem.p.val")

write.csv(collect.scores, file=paste0(PARAM$folder.result, 
                                    Sys.Date(),"_",
                                    "real.data.scores.csv"))



# g	number of bins to use to calculate quantiles.
