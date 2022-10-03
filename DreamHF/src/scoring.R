#example codes that can be run both on synthetic and real dataset
# load libraries 
library(tidyverse)
library(survminer)
library(survival)
library(microbiome)
library(magrittr)
library(randomForestSRC)
library(ResourceSelection)
library(GGally)
set.seed(198657)

args=(commandArgs(TRUE))
PARAM <- list()

#this argument is to specify the path of input folder
#the input folder structure is similar to DreamHF.zip both for synthetic and real dataset
PARAM$folder.R <- paste0(args[1]) 
#one output file (score.csv) has to be created in Team_Name_Submission_Number folder
#in bellow example your submission name is : FINRISK_TEST_1, please change to your submission name
#please avoid using (.) in your team name
dir.create(file.path(PARAM$folder.R, "FINRISK_TEST_1","output"))
PARAM$folder.data <- paste0(PARAM$folder.R, "/")
PARAM$folder.result <- paste0(PARAM$folder.data, "FINRISK_TEST_1/output/")

#load data

# load data
print("Load data")
S.test <- read.csv(file = paste0(PARAM$folder.data, 
                                 "test/pheno_test.csv"), 
                   row.names=1,)
S.train <- read.csv(file = paste0(PARAM$folder.data, 
                                  "train/pheno_training.csv"),
                   row.names=1,)

O.test <- read.csv(file = paste0(PARAM$folder.data, 
                                 "test/readcounts_test.csv"), 
                   row.names=1,)
O.train <- read.csv(file = paste0(PARAM$folder.data, 
                                  "train/readcounts_training.csv"), 
                    row.names=1,)

endpoints <- c("Event_time", "Event")

# change rownames with fake ids for convenience
rows <- rownames(O.train)  
seq=seq(1:c(length(rows)))
fakename<- sapply(seq,
                  function(x) paste('bacteria', x))
fakename = gsub(" ","",as.character(fakename))
rownames(O.train) <- fakename

df.train <- cbind(meta(S.train), t(O.train))
df.train$Event_time<- as.numeric(df.train$Event_time)
# exclude any minus value
df.train <- subset(df.train, Event_time > 0 & Event_time < 17)
df.train <- subset(df.train, !is.na(Event_time))
df.train <- subset(df.train, !is.na(Event))  

# define potential predictors
predictors <- c(colnames(S.train), rownames(O.train))

# choose only HF cases and eliminate "NA"
df.train <- subset(df.train, select=predictors, PrevalentHFAIL==0&!is.na(Event))

# remove PrevalentHFAIL from all matrix
predictors <- predictors[! predictors %in% "PrevalentHFAIL"]
df.train = df.train %>% 
  select(!PrevalentHFAIL) 
###############################################################################
# Fix and combine train and test datasets
# change rownames with fake ids for convenience
rownames(O.test) <- fakename

df.test <- cbind(meta(S.test), t(O.test))
df.test$Event_time <- as.numeric(df.test$Event_time)
# exclude any minus value
df.test <- subset(df.test, Event_time > 0 & Event_time < 17)
df.test = df.test %>% 
  select(!PrevalentHFAIL) 

# remove unnecessary files and save the main files
rm(O.train, O.test, S.test, S.train, fakename)
predictors <- setdiff(predictors, endpoints)
# Random survival forest********************************************************
# Random survival forest********************************************************
set.seed(9876)

# Labels are given for parameter optimization
dff1 <- as.data.frame(df.train[, c(endpoints, predictors)])
rownames(dff1)<-NULL
dff1=dff1[1:100,1:200]

# Labels are hidden for the second data subset 
dff2 <- as.data.frame(df.test[, predictors])
rownames(dff2)<-NULL
# Train Random Survival Model with the labeled training data
rf1 <- rfsrc(Surv(Event_time, Event) ~ .,
             data = dff1,
             importance = FALSE)

# Predictions for the test data
samples <- which(rowSums(is.na(dff2))==0)

pred <- predict(rf1, dff2[samples, ])
scores <- pred$predicted; names(scores) <- rownames(df.test)[samples]
res2 <- data.frame(SampleID=rownames(df.test), 
                  Score=sample(scores[rownames(df.test)]))  # Fully random
res <- res %>%  drop_na()

# Write the scoring table for evaluation
write.csv(res, file=paste0(PARAM$folder.result, 
                           "scores.csv"), 
          quote=FALSE, row.names=FALSE)

# read score file for test data
scores <- read.csv(file = paste0(PARAM$folder.result, 
                                  "scores.csv"))

# Harrells C *******************************************************************
labels=df.test[samples]
labels$SampleID<- rownames(df.test)

# normalize the scores btw 0-1
range01 <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}

# range the real values
labels$Event_time = range01(labels$Event_time)

# range the predicted values
# only apply 0-1 when it is not all 0s
if (length(unique(scores$Score))==1){
  scores$Score=scores$Score
} else {
  scores$Score = range01(scores$Score)
}

# Align the user provided scores with the true event times
true.scores <- as.numeric(labels[scores$SampleID,"Event_time"])

# Calculate Harrell's C statistics
C <- Hmisc::rcorr.cens(scores$Score, true.scores, outx=FALSE)
print(C)

# Hoslem test ***************************************************************** 
# https://stackoverflow.com/questions/44368326/hosmer-lemeshow-statistic-in-r

hosmerlem <- function(y, yhat, g=10)
       {cutyhat <- cut(yhat, breaks = quantile(yhat, 
                                               probs=seq(0,1, 1/g)), 
                       include.lowest=TRUE)
       obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
       expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
       chisq <- sum((obs - expect)^2/expect)
       P <- 1 - pchisq(chisq, g - 2)
       return(list(chisq=chisq,p.value=P))}

# remove NAs  for us, not valid for participants
scores <- scores %>%  drop_na()
hoslem <- hosmerlem(y=true.scores, yhat=scores$Score)
hoslem$harrelC <- C["C Index"][[1]]

# run with event info not event time
eventinfo <- as.numeric(labels[scores$SampleID,"Event"])
test <- hoslem.test(eventinfo, scores$Score) # x numeric observations, y expected values.
hoslem$resourcePval <- test$p.value

# Write the scoring table for evaluation
write.csv(hoslem, file=paste0(PARAM$folder.result, 
                         "stats.csv"), 
          quote=FALSE, row.names=FALSE)
