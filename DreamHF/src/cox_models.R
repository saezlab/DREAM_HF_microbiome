# load libraries 
library(tidyverse)
library(survminer)
library(survival)
library(microbiome)
library(magrittr)
set.seed(198657)

mainDir <- getwd()
subDir <- "output"

# create output folder
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
###############################################################################
# set parameters
alpha_level = 1
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

# choose only HF cases and eliminate "NA"
data <- subset(data, select=predictors, PrevalentHFAIL==0&!is.na(Event))

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
###############################################################################
# Univariate Cox regression analysis
###############################################################################
# A Cox regression of time to death on the time-constant covariates is specified as follow:
# cox ph model: allows both categorical and numeric variables 
# with categorical variable
# we cant estimate survival/hazard  but estimate hazard ratio
print("Cox survival univariate model building")
# the baseline survival function
cox.mod.subset <- coxph(Surv(Event_time, Event) ~ Sex + Age, 
                        data=data)
summary(cox.mod.subset)
print(cox.zph(cox.mod.subset))
# forest plot, graphical summary of a Cox model
pdf(file=paste0(PARAM$folder.result, 
                Sys.Date(), "_",
                "ggforest_cox.subset.pdf"))
ggforest(cox.mod.subset, data=data)
dev.off() # Turn the PDF device off

# All meta variables *****************************
cox.mod.total.meta <- coxph(Surv(Event_time, Event) ~ BodyMassIndex + Age + Smoking + 
                              PrevalentDiabetes + SystolicBP + BPTreatment + 
                              NonHDLcholesterol + PrevalentCHD + Sex + PrevalentHFAIL, 
                            data=data)
summary(cox.mod.total.meta)
print(cox.zph(cox.mod.total.meta))
# forest plot, graphical summary of a Cox model
pdf(file=paste0(PARAM$folder.result, 
                Sys.Date(), "_",
                "ggforest_cox.total.pdf"))
ggforest(cox.mod.total.meta, data=data)
dev.off() # Turn the PDF device off

###############################################################################
# univariate models for all predictors
univ_formulas <- sapply(predictors,
                        function(x) as.formula(paste('Surv(Event_time, Event) ~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = data)})
# Results *****************************
print("Univariate Model Results")
results <- lapply(predictors, function(x) {
  df <- summary(univ_models[[x]])$coefficients %>% as.data.frame()
  df <- df[nrow(df), ] %>% 
    dplyr::select(coef, "se(coef)", "z", "Pr(>|z|)") %>% 
    set_colnames(c("coef", "se_coef", "west_stat_value", "p")) %>% 
    dplyr::mutate(west_stat = "Wald")
  df <- df %>% 
    mutate(predictor = x)
}) %>% 
  do.call(rbind, .) 

# Multiple testing correction and results to better present
results <- results %>%
  dplyr::mutate(P_adjusted = p.adjust(p, "BH")) %>% 
  ungroup %>% 
  dplyr::mutate(PH = exp(coef)) %>% 
  dplyr::mutate(p_adj = P_adjusted) %>% 
  dplyr::mutate(direction = ifelse(coef < 0, "negative", "positive")) %>% 
  dplyr::group_by(predictor) 

# Results in neat form for presentation
neat_results <- results %>% 
  # filter(p == min(p)) %>% 
  dplyr:: ungroup() %>% 
  dplyr::mutate(HR = round(exp(coef),3)) %>% 
  dplyr::mutate(HR_lower_95 = round(exp(coef - 1.96*se_coef), 3),
         HR_upper_95 = round(exp(coef + 1.96*se_coef), 3),
         P = round(p, 5),
         Coefficient = round(coef, 3),
         "Coefficient SE" = round(se_coef, 3)) %>% 
  dplyr:: mutate(HR = paste0(HR, " (95% CI, ", HR_lower_95, "-", HR_upper_95, ")")) %>% 
  dplyr::select(Predictor = predictor, Coefficient, "Coefficient SE", HR, "p","P_adjusted", "west_stat_value", "west_stat") %>%
  dplyr::mutate(HR = ifelse(is.na(Coefficient), NA, HR))  %>% 
  dplyr::filter(P_adjusted < alpha_level) %>% 
  dplyr::arrange(P_adjusted) %>% 
  set_colnames(c("Predictor", "Coefficient", "Coefficient SE", "HR","P-value" ,"P (adjusted)", "west Statistic Value", "west Statistic"))
# export file
write.csv(neat_results, file=paste0(PARAM$folder.result, 
                              Sys.Date(), "_",
                              "univariate_results.csv"))

###############################################################################
# Get validation scores for all univariate and calculate Harrellc
###############################################################################

predict.uni <- function(univ_models, test.df, filename){
  print("Validation scores of Cox model in test or validation data")
  # Validation *****************************
  for (i in predictors){
    scores=as.data.frame(predict(univ_models[[i]], 
                                 newdata=test.df,     
                                 # the expected number of events given the covariates and follow-up time
                                 se.type="expected")) 
    scores = scores %>%
      set_colnames(c("Score")) %>%
      rownames_to_column(var = "SampleID")
    
    univ_models[[i]]$scores <- scores
    
    # Validation *****************************
    # Read the real labels that were hidden from that algorithm
    labels <- test.df
    labels$SampleID<- rownames(test.df)
    
    # Align the user provided scores with the true event times
    true.scores <- as.numeric(labels[scores$SampleID,"Event_time"])
    
    # Calculate Harrell's C statistics
    HarrellC <- Hmisc::rcorr.cens(scores$Score, true.scores, outx=FALSE)
    
    print(HarrellC["C Index"])
    univ_models[[i]]$HarrellC <- HarrellC["C Index"]
  }
  return(univ_models)
}

# Validation of univariate cox model in TEST and VALIDATION cohorts ************
univ_models_ontest <- predict.uni(univ_models, df.test, "ontest")
univ_models_onval <- predict.uni(univ_models, df.val, "onvalidation")

# Multivariable model **********************************************************
print("Cox survival multivariate model building")
cox.multivar <- coxph(Surv(Event_time, Event)~., data=data)

predict.multi <- function(multi_model, test.df, filename){
  print("Validation scores of Cox model")
  
# Validation *******************************************************************
  scores=as.data.frame(predict(multi_model, 
                               newdata=test.df, 
                               se.type="expected"))
  scores = scores %>%
    set_colnames(c("Score")) %>%
    rownames_to_column(var = "SampleID")
  
  multi_model$scores <- scores
  
  # HarrellC index *************************************************************
  print("HarrellC index of Cox model")
  # Read the real labels that were hidden from that algorithm
  labels <- test.df
  labels$SampleID<- rownames(test.df)
  
  # Align the user provided scores with the true event times
  true.scores <- as.numeric(labels[scores$SampleID,"Event_time"])
  
  # Calculate Harrell's C statistics
  HarrellC <- Hmisc::rcorr.cens(scores$Score, true.scores, outx=FALSE)
  multi_model$HarrellC <- HarrellC["C Index"]
  print(multi_model$HarrellC)
  return(multi_model)
}

# Validation of multivariate cox model in TEST and VALIDATION cohorts **********
multicox_ontest <- predict.multi(cox.multivar, df.test, "cox_multivar_ontest_")
multicox_onval <- predict.multi(cox.multivar, df.val, "cox_multivar_onvalidation_")

###############################################################################
# Hosmer lemesrow test
###############################################################################
source("src/Hosmer.test_etc.R")
# Multivariate models **********************************************************
model.list <- list(cox_subset = cox.mod.subset, 
                   cox_meta_total = cox.mod.total.meta, 
                   cox_multivar = cox.multivar)
collect.hoslem <- data.frame()

for (model in model.list){
  print(model$formula)
  pred <- Coxar(model, years=16.9, bh=NULL)
  test <- HosLem.test(model$y, pred, plot=FALSE)
  p.val <-  test$pval
  collect.hoslem <- rbind(collect.hoslem, p.val)
  colnames(collect.hoslem) <- "p.val"
}
rownames(collect.hoslem) <- names(model.list)
collect.hoslem$p.adj <- p.adjust(collect.hoslem$p.val)

# Univariate models ************************************************************
collect.hoslem.u <- data.frame()

for (mName in predictors){
  print(univ_models[[mName]]$formula)
  pred <- Coxar(univ_models[[mName]], years=16.9)
  test <- HosLem.test(univ_models[[mName]]$y, pred, plot=FALSE)
  p.val <-  test$pval
  collect.hoslem.u <- rbind(collect.hoslem.u, p.val)
  colnames(collect.hoslem.u) <- "p.val"
  
}
rownames(collect.hoslem.u) <- predictors
collect.hoslem.u$p.adj <- p.adjust(collect.hoslem.u$p.val)

# export file
hoslem.total <- rbind(collect.hoslem, collect.hoslem.u)
write.csv(hoslem.total, file=paste0(PARAM$folder.result, 
                                      Sys.Date(),"_",
                                      "Hosmer_lemesrow_scores.csv"))
# save.image(paste0(PARAM$folder.result, Sys.Date(),"cox.results.RData"))
