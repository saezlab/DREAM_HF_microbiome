# load libraries 
library(tidyverse)
library(survminer)
library(survival)
library(microbiome)
library(magrittr)
library(randomForestSRC)
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
source(paste0(mainDir,"/src/hosmerTest.R"))

# load data
print("Load data")

S.train <- read.csv(file = paste0(PARAM$folder.data, 
                                  "train/pheno_training.csv"),
                   row.names=1,)
S.test <- read.csv(file = paste0(PARAM$folder.data, 
                                 "test/pheno_test.csv"), 
                   row.names=1,)
S.val <- read.csv(file = paste0(PARAM$folder.data, 
                                "validation/pheno_validation.csv"), 
                  row.names=1,)

O.test <- read.csv(file = paste0(PARAM$folder.data, 
                                 "test/readcounts_test.csv"), 
                   row.names=1,)
O.train <- read.csv(file = paste0(PARAM$folder.data, 
                                  "train/readcounts_training.csv"), 
                    row.names=1,)
O.val <- read.csv(file = paste0(PARAM$folder.data, 
                                "validation/readcounts_validation.csv"), 
                  row.names=1,)
###############################################################################
# set parameters
alpha_level = 1
endpoints <- c("Event_time", "Event")

# change rownames with fake ids for convenience
rows <- rownames(O.train)  
seq=seq(1:c(length(rows)))
fakename<- sapply(seq,
                  function(x) paste('bacteria', x))
fakename = gsub(" ","",as.character(fakename))
rownames(O.train) <- fakename

data <- cbind(meta(S.train), t(O.train))
data$Event_time<- as.numeric(data$Event_time)
# exclude any minus value
data <- subset(data, Event_time > 0 & Event_time < 17)
data <- subset(data, !is.na(Event_time))
data <- subset(data, !is.na(Event))  

# define potential predictors
predictors <- c(colnames(S.train), rownames(O.train))

# choose only HF cases and eliminate "NA"
data <- subset(data, select=predictors, PrevalentHFAIL==0&!is.na(Event))

# remove PrevalentHFAIL from all matrix
predictors <- predictors[! predictors %in% "PrevalentHFAIL"]
data = data %>% 
  select(!PrevalentHFAIL) 
###############################################################################
# Fix and combine train and test datasets
# change rownames with fake ids for convenience
rownames(O.test) <- fakename
rownames(O.val) <- fakename

df.test <- cbind(meta(S.test), t(O.test))
df.test$Event_time <- as.numeric(df.test$Event_time)
# exclude any minus value
df.test <- subset(df.test, Event_time > 0 & Event_time < 17)
df.test = df.test %>% 
  select(!PrevalentHFAIL) 

df.val <- cbind(meta(S.val), t(O.val))
df.val$Event_time <- as.numeric(df.val$Event_time)
# exclude any minus value
df.val <- subset(df.val, Event_time > 0 & Event_time < 17)
df.val = df.val %>% 
  select(!PrevalentHFAIL) 

# remove unnecessary files and save the main files
rm(O.train, O.test, O.val, S.test, S.val, S.train,fakename )
save(data, df.test, df.val, file=paste0(PARAM$folder.result,"dream_synthetic_files.RData"))

###############################################################################
#SUBSET FOR TESTING
# data=data[1:300,1:20]
# predictors=predictors[1:20]

predictors <- setdiff(predictors, endpoints)
###############################################################################
# Univariate Cox regression analysis*******************************************

# A Cox regression of time to death on the time-constant covariates is specified as follow:
# cox ph model: allows both categorical and numeric variables 
# with categorical variable
# we cant estimate survival/hazard  but estimate hazard ratio
# the baseline survival function
# 
# data$age_fc <- cut(data$Age, c(0,54, 59,64, 69,74,79, 89, 110),
# labels = c(paste(c(50,55,60,65,70,75,80),
# c(54,59,64,69,74,79,89), sep='-'), "90+"))

cox.mod.sex <- coxph(Surv(Event_time, Event) ~ Sex+Age,
                     data=data) 

summary(cox.mod.sex)
print(cox.zph(cox.mod.sex))
# forest plot, graphical summary of a Cox model
pdf(file=paste0(PARAM$folder.result, 
                Sys.Date(), "_",
                "ggforest_cox.subset.pdf"))
ggforest(cox.mod.sex, data=data)
dev.off() # Turn the PDF device off

# All meta variables *****************************
cox.mod.total.meta <- coxph(Surv(Event_time, Event) ~ BodyMassIndex + Smoking + 
                              PrevalentDiabetes + SystolicBP + BPTreatment + 
                              NonHDLcholesterol + PrevalentCHD + Sex  + Age, 
                            data=data)
summary(cox.mod.total.meta)
print(cox.zph(cox.mod.total.meta))
# forest plot, graphical summary of a Cox model
pdf(file=paste0(PARAM$folder.result, 
                Sys.Date(), "_",
                "ggforest_cox.total.pdf"))
ggforest(cox.mod.total.meta, data=data)
dev.off() # Turn the PDF device off

# All univariate variables *****************************
# univariate models for all predictors
univ_formulas <- sapply(predictors,
                        function(x) as.formula(paste('Surv(Event_time, Event) ~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = data)})
univ_models$age_sex = cox.mod.sex
univ_models$total_cov = cox.mod.total.meta

###############################################################################
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
  dplyr::select(Predictor = predictor, Coefficient, "Coefficient SE", HR, "p",
                "P_adjusted", "west_stat_value", "west_stat") %>%
  dplyr::mutate(HR = ifelse(is.na(Coefficient), NA, HR))  %>% 
  dplyr::filter(P_adjusted < alpha_level) %>% 
  dplyr::arrange(P_adjusted) %>% 
  set_colnames(c("Predictor", "Coefficient", "Coefficient SE", "HR","P-value",
                 "P (adjusted)", "west Statistic Value", "west Statistic"))
# export file
write.csv(neat_results, file=paste0(PARAM$folder.result, 
                                    Sys.Date(), "_",
                                    "univariate_results.csv"))

###############################################################################
# Get validation scores for univariate models and calculate Harrellc
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

save(cox.mod.sex, cox.mod.total.meta, univ_models_ontest, univ_models_onval, 
     file=paste0(PARAM$folder.result, "univ_models.RData"))

# Multivariate models **********************************************************
# Calculate HarrelC
model.list <- list(cox_subset = cox.mod.sex, 
                   cox_meta_total = cox.mod.total.meta)

collect.harrelc <- data.frame()
collect.hoslem <- data.frame()

for (model in model.list){  
  print(model$formula)
  scores=as.data.frame(predict(model, 
                               newdata=df.test,     
                               # the expected number of events given the covariates and follow-up time
                               se.type="expected")) 
  scores = scores %>%
    set_colnames(c("Score")) %>%
    rownames_to_column(var = "SampleID")
  labels <- df.test
  labels$SampleID<- rownames(df.test)
  
  # Align the user provided scores with the true event times
  true.scores <- as.numeric(labels[scores$SampleID,"Event_time"])
  
  # Calculate Harrell's C statistics
  # Calculate Harrell's C statistics
  HarrellC <- Hmisc::rcorr.cens(scores$Score, true.scores, outx=FALSE)
  collect.harrelc <- rbind(collect.harrelc, HarrellC["C Index"])
  colnames(collect.harrelc) <- "HarrellC"

  # Calculate Hoslem
  pred <- Coxar(model, years=10, bh=NULL)
  test <- HosLem.test(model$y, pred, plot=FALSE)
  p.val <-  test$pval
  collect.hoslem <- rbind(collect.hoslem, p.val)
  colnames(collect.hoslem) <- "p.val"
}

rownames(collect.hoslem) <- names(model.list)
collect.hoslem$p.adj <- p.adjust(collect.hoslem$p.val)
rownames(collect.harrelc) <- names(model.list)
totalscores <-cbind(collect.harrelc, collect.hoslem)

# export file
write.csv(totalscores, file=paste0(PARAM$folder.result, 
                                    Sys.Date(),"_",
                                    "cox_scores.csv"))

# Hosmer lemesrow test based on library(ResourceSelection)**********************
library(ResourceSelection)
library(GGally)
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
      rownames_to_column(var = "SampleID") %>%
      drop_na()
    # normalize range from 0-1
    univ_models[[i]]$scores <- scores
    print(i)
    # Validation *****************************
    # Read the real labels that were hidden from that algorithm
    labels <- test.df
    labels$SampleID <- rownames(test.df)
    
    # normalize the scores btw 0-1
    range01 <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}
    labels$Event_time = range01(labels$Event_time)
   
     # only apply 0-1 when it is not all 0s
    if (length(unique(scores$Score))==1){
      scores$Score=scores$Score
    } else {
      scores$Score = range01(scores$Score)
    }

    # Align the user provided scores with the true event times
    true.scores <- as.numeric(labels[scores$SampleID,"Event_time"])
    
    # Calculate Harrell's C statistics
    HarrellC <- Hmisc::rcorr.cens(scores$Score, true.scores, outx=FALSE)
    
    print(HarrellC["C Index"])
    univ_models[[i]]$HarrellC <- HarrellC["C Index"]
    
    # hoslem test from 
    if (length(unique(scores$Score))==1){
      univ_models[[i]]$Hoslem.pval.rs <- "untested"
    } else {
      test <- hoslem.test(true.scores, scores$Score) # x numeric observations, y expected values.
      print(test)
      univ_models[[i]]$Hoslem.pval.rs <- test$p.value
    }
    
    # test <- hoslem.test(true.scores, scores$Score) # x numeric observations, y expected values.
    # print(test)
    # univ_models[[i]]$Hoslem.pval.rs <- test$p.value
    uni.models.rs=univ_models
  }
  return(uni.models.rs)
}

# Validation of univariate cox model in TEST and VALIDATION cohorts ************
univ_models_ontest <- predict.uni(univ_models, df.test, "ontest_rs")
univ_models_onval <- predict.uni(univ_models, df.val, "onvalidation_rs")


# plot cox 
# Create the new data  
sex_df <- with(data,
               data.frame(Sex = c(0, 1), 
                          Age = rep(mean(Age, na.rm = TRUE), 2)))

fit <- survfit(cox.mod.sex, newdata = sex_df)
plot.agesex=ggsurvplot(fit, conf.int = FALSE, data = sex_df, palette = "Dark2", 
                 censor = FALSE, surv.median.line = "hv")

ggsave(file = paste0(PARAM$folder.result, 
                     Sys.Date(), "_",
                     "ggsurvplot.age_sex.pdf"), plot.agesex$plot)


###############################################################################
# Random survival forest********************************************************
set.seed(9876)

# Labels are given for parameter optimization
dff1 <- as.data.frame(data[, c(endpoints, predictors)])
rownames(dff1)<-NULL

# Labels are hidden for the second data subset 
dff2 <- as.data.frame(df.test[, predictors])
rownames(dff2)<-NULL

# Train Random Survival Model with the labeled training data
rf1 <- rfsrc(Surv(Event_time, Event) ~ .,
             data = dff1,
             importance = FALSE)

# Predictions for the test data
samples <- which(rowSums(is.na(dff2))==0)
pred2 <- predict(rf1, dff2[samples, ])
scores <- pred2$predicted; names(scores) <- rownames(df.test)[samples]
res <- data.frame(SampleID=rownames(df.test), 
                  Score=sample(scores[rownames(df.test)]))  # Fully random
res <- res %>%  drop_na()
# Write the scoring table for evaluation
write.csv(res, file=paste0(PARAM$folder.result, 
                           Sys.Date(), 
                           "random_forest_scores.csv"), 
          quote=FALSE, row.names=FALSE)

# Harrells C *******************************************************************
labels=data
labels$SampleID<- rownames(data)

# Align the user provided scores with the true event times
true.scores <- as.numeric(labels[res$SampleID,"Event_time"])

# Calculate Harrell's C statistics
C <- Hmisc::rcorr.cens(res$Score, true.scores, outx=FALSE)
print(C)

Cindex.train <- get.cindex(rf1$yvar[,1], rf1$yvar[,2], rf1$predicted.oob)

# Hoslem test ***************************************************************** 
test.rs <- hoslem.test(true.scores, res$Score) # x numeric observations, y expected values.
#with aki code
pred <- Coxar(model, years=10, bh=NULL)
test <- HosLem.test(model$y, pred, plot=FALSE)

oo <- subsample(rf1)

#combine all
stats_rfs <- data.frame(harrellc.test = C[[1]], 
                        harrellc.train = Cindex.train,
                        hosLem.pval = test$pval,
                        hoslemrs.pval= test.rs$p.value) 
# export stats
write.csv(stats_rfs, file=paste0(PARAM$folder.result, 
                         Sys.Date(), 
                         "random_forest_stats.csv"))


pdf(paste0(PARAM$folder.result, "RSF_overview.pdf"), 
    width = 10, 
    height = 8)
par(oma = c(0.5, 0.5, 0.5, 0.5))
plot.survival(rf1, cens.model = "rfsrc",  show.plots = TRUE)
dev.off()

# plot rfs
pdf(paste0(PARAM$folder.result, "Variable Importance RSF.pdf"), 
    width = 10, 
    height = 15)
par(oma = c(0.5, 0.5, 0.5, 0.5))
par(cex.axis = 2.0,
    cex.lab = 2.0,
    cex.main = 2.0,
    mar = c(6.0,17,1,1),
    mgp = c(4, 1, 0))

plot.subsample(oo, alpha = .01, standardize = TRUE, normal = TRUE, pmax = 20)

dev.off()
