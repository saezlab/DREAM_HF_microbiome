# load libraries 
library(tidyverse)
library(ggplot2)
library(car)

# This document contains the code to compare the INT results with real Finrisk data 
set.seed(1978)
# set parameters
PARAM <- list()
setwd("/Users/ecekartal/Documents/SAEZ/projects/DREAM_HF_microbiome/")
PARAM$folder.R <- paste0(getwd(), "/")
PARAM$folder.result <- paste0(PARAM$folder.R, "results/")

# Finrisk cohort data
S1 <- read.csv(file = paste0(PARAM$folder.R, "data/train/phenodata.csv"), row.names=1, )
# INT-simulated data example
S2 <- read.csv(file = paste0(PARAM$folder.R, "data/test/phenodata.csv"), row.names=1)

# Finrisk cohort data split
F.test <- read.csv(file = paste0(PARAM$folder.R, "data/finrisk_dataset/test/pheno_test.csv"), row.names=1,)
F.train <- read.csv(file = paste0(PARAM$folder.R, "data/finrisk_dataset/train/pheno_training.csv"), row.names=1,)
F.val <- read.csv(file = paste0(PARAM$folder.R, "data/finrisk_dataset/validation/pheno_validation.csv"), row.names=1,)

# INT-simulated data split
S.test <- read.csv(file = paste0(PARAM$folder.R, "data/synthetic_dataset/test/pheno_test.csv"), row.names=1,)
S.train <- read.csv(file = paste0(PARAM$folder.R, "data/synthetic_dataset/train/pheno_training.csv"), row.names=1,)
S.val <- read.csv(file = paste0(PARAM$folder.R, "data/synthetic_dataset/validation/pheno_validation.csv"), row.names=1,)

# bind 2 meta datasets
S1$dataset<- paste0(0)
S2$dataset<- paste0(1)
S1=S1[,-c(7:8)] ##REMOVE THIS SINCE YOU HAVE EVENT INFO FOR BOTH DATASETS
metatotal<- rbind(S1, S2)
metatotal$dataset <- as.numeric(metatotal$dataset)

###############################################################################
# first we compare the overall synthetic data vs Finrisk data
###############################################################################
# QQ plots
# subset continues metadata 
idx = c( "Age", "BMI", "SystolicBP","NonHDLcholesterol", "Event_time")# add eventtime since you have this info

for (i in idx){
  x <- S1[[i]]
  y <- S2[[i]]
  p.qq <- ggplot(mapping = aes(x = sort(x), y = sort(y))) +
    geom_point() +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    xlab(paste0("Finrisk Dataset ",i))+
    ylab(paste0("Synthetic Dataset ", i))
  ggsave(filename=paste0(PARAM$folder.result, Sys.Date(), i, "_QQplot_FinriskvsSynthetic.pdf"),
         plot=p.qq)
}
###############################################################################
# subset categorical metadata for bar plot
idx = c("Smoking", "BPtreatment", "PrevalentDiabetes", "PrevalentCHD", "Sex", "dataset")
         #"Event",  include this since you have this info
df=melt(id="dataset", metatotal[idx]) %>% 
  group_by(variable,dataset) %>% 
  count(value)

p.bar <- ggplot(df,aes(x=as.factor(value), y=n, fill=as.factor(dataset))) +
  geom_bar(stat="identity", position = "dodge", alpha=.6) +
  facet_wrap(variable~., scales = "free") +
  labs(y = "Count", color = "Dataset", 
       fill = "Dataset",
       title = "Comparison of two datasets:0 for Finrisk and 1 for INN-simulated")
ggsave(filename=paste0(PARAM$folder.result, Sys.Date(),"_FinriskvsSynthetic_barplot.pdf"),        
        plot=p.bar) 
 
###############################################################################
# we compare the train-test-validation synthetic data vs Finrisk data
###############################################################################
# QQ plots
# subset continues metadata 
idx = c( "Age", "BodyMassIndex","SystolicBP", "NonHDLcholesterol","Event_time")
# HERE I COPY PASTE JUST FOR BEING LAZY FOR 3 DIFFIRENT COMBINATIONS
for (i in idx){
  x <- F.test[[i]]
  y <- S.test[[i]]
  p.qq <- ggplot(mapping = aes(x = sort(x), y = sort(y))) +
    geom_point() +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    xlab(paste0("Finrisk Dataset Test ",i))+
    ylab(paste0("Synthetic Dataset Test ", i))
  ggsave(filename=paste0(PARAM$folder.result, Sys.Date(),"_TEST_", i, "_QQplot_FinriskvsSynthetic.pdf"),
         width = 5, height = 5, plot=p.qq) 
}

for (i in idx){
  x <- F.train[[i]]
  y <- S.train[[i]]
  p.qq <- ggplot(mapping = aes(x = sort(x), y = sort(y))) +
    geom_point() +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    xlab(paste0("Finrisk Dataset Train ",i))+
    ylab(paste0("Synthetic Dataset Train ", i))
  ggsave(filename=paste0(PARAM$folder.result, Sys.Date(),"_TRAIN_", i, "_QQplot_FinriskvsSynthetic.pdf"),
         width = 5, height = 5, plot=p.qq) 
}
for (i in idx){
  x <- F.val[[i]]
  y <- S.val[[i]]
  p.qq <- ggplot(mapping = aes(x = sort(x), y = sort(y))) +
    geom_point() +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    xlab(paste0("Finrisk Dataset Validation ",i))+
    ylab(paste0("Synthetic Dataset Validation ", i))
  ggsave(filename=paste0(PARAM$folder.result, Sys.Date(),"_VALIDATION_", i, "_QQplot_FinriskvsSynthetic.pdf"),
         width = 5, height = 5, plot=p.qq) 
}

###############################################################################
# subset categorical metadata for bar plot
idx = c("Smoking", "BPTreatment", "PrevalentDiabetes", "PrevalentCHD", "PrevalentHFAIL", "Sex","Event", "dataset")
# bind meta datasets
F.val$dataset<- paste0("Fval")
F.test$dataset<- paste0("Ftest")
F.train$dataset<- paste0("Ftrain")

S.val$dataset<- paste0("Sval")
S.test$dataset<- paste0("Stest")
S.train$dataset<- paste0("Strain")

total<- rbind(F.train, F.test, F.val, S.train, S.test, S.val)

df=melt(id="dataset", total[idx]) %>% 
  group_by(variable,dataset) %>% 
  count(value)

p.bar <- ggplot(df,aes(x=as.factor(value), y=n, fill=as.factor(dataset))) +
  geom_bar(stat="identity", position = "dodge", alpha=.6) +
  facet_wrap(variable~., scales = "free") +
  labs(y = "Count", color = "Dataset", 
       fill = "Dataset",
       title = "Comparison of two datasets:0 for Finrisk and 1 for INN-simulated")
ggsave(filename=paste0(PARAM$folder.result, Sys.Date(),"_FinriskvsSynthetic_barplot_subgroups.pdf"),        
       plot=p.bar) 
