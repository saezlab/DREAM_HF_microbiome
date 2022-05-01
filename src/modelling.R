# Read data
# source("set_up.R")
#source("read_data.R")
setwd("~/Documents/SAEZ/projects/DREAM_HF_microbiome")

library(tidyverse)
library(SIAMCAT)
# This document contains the HF dream challenge modelling
# Modelling with siamcat-lasso regression 
set.seed(5497659)

# load the data
PARAM <- list()
PARAM$folder.R <- paste0(getwd(), "/")
PARAM$folder.testdata <- paste0(PARAM$folder.R, "data/test/")
PARAM$folder.traindata <- paste0(PARAM$folder.R, "data/train/")

S1 <- read.csv(file = paste0(PARAM$folder.R, "data/train/phenodata.csv"), row.names=1)
O1 <- read.csv(file = paste0(PARAM$folder.R, "data/train/readcounts.csv"), row.names=1)
T1 <- read.csv(file = paste0(PARAM$folder.R, "data/train/taxtable.csv"))


# metadata
meta.train <- S1

# choose only HF cases and eliminate "NA"
meta.train <- meta.train %>%
  filter(Event!="NA")
table(meta.train$Event)
dim(meta.train)

# read files for training data
feat.train <- O1
# remove unmapped
toremove = c("ssRNA_positive-strand_viruses", "1", "0")
feat.train <- feat.train[!rownames(feat.train) %in% toremove, ]
rows=rownames(feat.train)
# Convert all variable types to numeric
feat.train <- as.data.frame(apply(feat.train, 2, as.numeric))  
rownames(feat.train) <- rows

# choose sample data
data.sample.tmp <- feat.train[match(rownames(meta.train), colnames(feat.train))]
dim(data.sample.tmp)
# relative abundance
feat.train.rel <- prop.table(as.matrix(data.sample.tmp), 2)
feat.train.rel <- feat.train.rel[rowSums(feat.train.rel >= 10^-5) >= 2, ]
dim(feat.train.rel)

# modelling
runsiamcat <- function(featTable, metaTable, fileName, label_new, case, ml, norm, cutoff.p){
  dim(featTable)
  fileName=paste(fileName, ml, norm, cutoff.p)
  print(fileName)
  
  # create SIAMCAT object and classify
  siamcat <- siamcat(feat=featTable, 
                     meta=metaTable, 
                     label=label_new, 
                     case=case)
  # abundance and prevelance filtering
  siamcat <- filter.features(siamcat,
                             filter.method = 'abundance',
                             cutoff=0.0001,
                             verbose=3)
  # 
  siamcat <- filter.features(siamcat,
                             filter.method = 'prevalence',
                             cutoff = cutoff.p,
                             feature.type = 'filtered',
                             verbose=3)
  
  # normalize with log.clr
  siamcat <- normalize.features(siamcat, 
                                norm.method = norm, 
                                feature.type = 'filtered',
                                norm.param = list(log.n0=1e-05, sd.min.q=0.1))
  
  # compute associations 
  siamcat <- check.associations(siamcat, feature.type = 'normalized',
                                detect.lim = 10^-5, 
                                plot.type = "quantile.box",
                                fn.plot = paste0(PARAM$folder.results,
                                                 Sys.Date(), fileName, 'assoc.plot.pdf'))
  
  # train model
  siamcat <- create.data.split(siamcat, num.folds =2, num.resample = 2)  
  siamcat <- train.model(siamcat, method = ml, verbose = 3)
  siamcat <- make.predictions(siamcat)
  siamcat <- evaluate.predictions(siamcat)    
  print(siamcat@eval_data$auroc)
  # evaluation plot
  model.evaluation.plot(siamcat, fn.plot = paste0(PARAM$folder.results, Sys.Date(), '.',
                                                  fileName, 'eval.plot.pdf'))
  # interpretation plot
  model.interpretation.plot(siamcat, fn.plot = paste0(PARAM$folder.results, 
                                                      Sys.Date(), '.', 
                                                      fileName,
                                                      'interpret.plot.pdf'),
                            consens.thres = 0.5,
                            detect.lim = 1e-05,
                            heatmap.type = 'zscore')
  
  # save siamcat object
  save(siamcat, file = paste0(PARAM$folder.results, Sys.Date(), '.', fileName, 'siamcat.Rdata'))
  return(siamcat) 
}

# choose variables for testing confounder
metatest <- c("Age","BMI", "Smoking", "BPtreatment" , "PrevalentDiabetes",
              "PrevalentCHD", "Event", "SystolicBP", "Sex", "NonHDLcholesterol")
# run modelling
siamcat.main <- runsiamcat(feat.train.rel, meta.train, "siamcat", "Event","1", "lasso", "log.clr", "0.03")
