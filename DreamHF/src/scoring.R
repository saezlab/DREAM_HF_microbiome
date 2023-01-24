#example codes that can be run both on synthetic and real dataset
# load libraries 
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
df.train = df.train %>% 
  dplyr::select(!PrevalentHFAIL) 
# remove PrevalentHFAIL from all matrix
predictors <- predictors[! predictors %in% "PrevalentHFAIL"]

###############################################################################
# Fix and combine train and test datasets
# change rownames with fake ids for convenience
rownames(O.test) <- fakename

df.test <- cbind(meta(S.test), t(O.test))
df.test$Event_time <- as.numeric(df.test$Event_time)
df.test = df.test %>% 
  dplyr::select(!PrevalentHFAIL) 
# exclude any minus value
df.test <- subset(df.test, Event_time > 0 & Event_time < 17)

# remove unnecessary files and save the main files
rm(O.train, O.test, S.test, S.train, fakename)
predictors <- setdiff(predictors, endpoints)


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
res <- data.frame(SampleID=rownames(df.test), 
                  Score=sample(scores[rownames(df.test)]))  # Fully random
res <- res %>%  drop_na()

# Write the scoring table for evaluation
write.csv(res, file=paste0(PARAM$folder.result, 
                           "scores.csv"), 
          quote=FALSE, row.names=FALSE)

# SCORING ONLY

# read score file for test data
scores <- read.csv(file = paste0(PARAM$folder.result, 
                                  "scores.csv"))
collect.rsf <- data.frame()

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
# remove NAs  for us, not valid for participants
scores <- scores %>%  drop_na()
eventinfo <- as.numeric(labels[scores$SampleID,"Event"])

# Calculate Harrell's C statistics
C <- Hmisc::rcorr.cens(scores$Score, true.scores, outx=FALSE)
print(C)

# Hoslem test ***************************************************************** 
# run with event info not event time
hoslem.test <- hoslem.test(eventinfo, scores$Score) 
collect.rsf <- rbind(collect.rsf, hoslem.test$p.value)

collect.rsf <- cbind(C["C Index"][[1]], collect.rsf)

# Write the scoring table for evaluation
write.csv(collect.rsf, file=paste0(PARAM$folder.result, 
                         "stats.rsf.csv"), 
          quote=FALSE, row.names=FALSE)


# Cox ********************************************************
cox.mod.sex <- coxph(Surv(Event_time, Event) ~ Sex+Age,
                     data=df.train) 

cox.mod.sex.s <- coxph(Surv(Event_time, Event) ~  Age + strata(Sex),
                       data=df.train) 

cox.mod.total.meta <- coxph(Surv(Event_time, Event) ~ BodyMassIndex + Smoking + 
                              PrevalentDiabetes + SystolicBP + BPTreatment + 
                              NonHDLcholesterol + PrevalentCHD + Sex  + Age, 
                            data=df.train)

cox.mod.total.meta.s <- coxph(Surv(Event_time, Event) ~ BodyMassIndex + Smoking + 
                                PrevalentDiabetes + SystolicBP + BPTreatment + 
                                NonHDLcholesterol + PrevalentCHD +  Age + strata(Sex), 
                              data=df.train)


model.list <- list(cox_sex = cox.mod.sex, 
                   cox_sex.s = cox.mod.sex.s,
                   cox_meta_total = cox.mod.total.meta,
                   cox_meta_total.s = cox.mod.total.meta.s)

collect.stats <- data.frame()
collect.harrel <- data.frame()

for (model in model.list){  
  print(model$formula)
  scores=as.data.frame(predict(model, 
                               newdata=df.test,     
                               se.type="expected")) 
  
  scores = scores %>%
    set_colnames(c("Score")) %>%
    rownames_to_column(var = "SampleID")
  scores <- scores %>%  drop_na()
  # range the predicted values
  # only apply 0-1 when it is not all 0s
  scores$Score = range01(scores$Score)
  
  
  # Align the user provided scores with the true event times
  true.scores <- as.numeric(labels[scores$SampleID,"Event_time"])
  
  # Calculate Harrell's C statistics
  C <- Hmisc::rcorr.cens(scores$Score, true.scores, outx=FALSE)
  print(C)
  collect.harrel<- rbind(collect.harrel, C["C Index"])
  
  # Calculate Hoslem
  # 1
  eventinfo <- as.numeric(labels[scores$SampleID,"Event"])
  hoslem.test <- hoslem.test(eventinfo, scores$Score) # x numeric observations, y expected values.
  collect.stats <- rbind(collect.stats, hoslem.test$p.value)

}

collect.stats <- cbind(collect.harrel, collect.stats)
collect.stats[5,] <-  collect.rsf

colnames(collect.stats) <- c("HarrellC", "hoslem.test")
rownames (collect.stats)<- c("cox_sex",
                             "cox_sex.s",
                             "cox_meta_total",
                             "cox_meta_total.s",
                             "random survival")

write.csv(collect.stats, file=paste0(PARAM$folder.result, 
                                     "stats.total.csv"), 
          quote=FALSE, row.names=FALSE)

# plotting
scores$dataset <- factor(scores$dataset,levels = c("train", "test", "scoring"))
scores$model <- factor(scores$model,levels = c("Cox_subset", "Cox_allcovariates", "RSF"))

scores$harrelC=as.numeric(scores$harrelC)

p= ggplot(scores, aes(x=dataset, y=harrelC,  fill=dataset)) +    
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  xlab("") + ylab("Harrell's C index")+
  geom_bar(stat="identity", color="black")  + theme_classic() +
  theme(axis.text.x = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 90, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"))+
  facet_grid(rows=.~scores$model)
ggsave(p, file=paste0(PARAM$folder.results, 
                               Sys.Date(),'_', 
                               "barplot_harrelc.pdf"), width = 5, height=4)
# Set up for ggplot
kmi <- rep("KM",length(km_fit$time))
km_df <- data.frame(km_fit$time,km_fit$surv,kmi)
names(km_df) <- c("Time","Surv","Model")

coxi <- rep("Cox",length(cox_fit$time))
cox_df <- data.frame(cox_fit$time,cox_fit$surv,coxi)
names(cox_df) <- c("Time","Surv","Model")

rfi <- rep("RF",length(r_fit$unique.death.times))
rf_df <- data.frame(r_fit$unique.death.times,avg_prob,rfi)
names(rf_df) <- c("Time","Surv","Model")

plot_df <- rbind(km_df,cox_df,rf_df)

p <- ggplot(plot_df, aes(x = Time, y = Surv, color = Model))
p + geom_line()



