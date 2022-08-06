test_algorithm <- function (S1, O1, S2, O2, file) {

  # Prepare features by combining phenodata and abundance data; add some filtering
  D1 <- apply(cbind(S1, t(O1)), 2, as.numeric) # Train
  D2 <- apply(cbind(S2, t(O2)), 2, as.numeric) # Test

  set.seed(36636333)
  endpoints <- c("Event_time", "Event")
  covariates <- setdiff(colnames(S1), c("SampleID", endpoints))
  taxa <- sample(rownames(O1), 100)
  predictors <- c(taxa, covariates)

  # Labels are given for parameter optimization
  dff1 <- as.data.frame(D1[, c(endpoints, predictors)])
  dff1 <- subset(dff1, Event_time > 0 & Event_time < 18)
  dff1 <- subset(dff1, !is.na(Event_time))
  dff1 <- subset(dff1, !is.na(Event))  

  # Labels are hidden for the second data subset 
  dff2 <- as.data.frame(D2[, predictors])

  # Train Random Survival Model with the labeled training data
  rf1 <- rfsrc(Surv(Event_time, Event) ~ .,
             data = dff1,
             importance = FALSE)
  
  # Predictions for the test data
  samples <- which(rowSums(is.na(dff2))==0)
  pred2 <- predict(rf1, dff2[samples, ])
  scores <- pred2$predicted; names(scores) <- rownames(S2)[samples]

  res <- data.frame(SampleID=rownames(S2), Score=sample(scores[rownames(S2)]))  # Fully random

  # Write the scoring table for evaluation
  write.csv(res, file=file, quote=FALSE, row.names=FALSE)
  return(res)

}

# Write predictions to a file
scores <- test_algorithm( S.train, O.train,S.test, O.test, file = paste0(PARAM$folder.result, "scores.csv"))
###############################################################################
# Harralls C 
labels=S.train
labels$SampleID<- rownames(S.train)

# Align the user provided scores with the true event times
true.scores <- as.numeric(labels[scores$SampleID,"Event_time"])

# Calculate Harrell's C statistics
C <- Hmisc::rcorr.cens(scores$Score, true.scores, outx=FALSE)
#print(C["C Index"])
write.csv(C, file=paste0(PARAM$folder,"HarrellC.csv"))
###############################################################################
# hosmer lemesrow test
options(width=200) # characters per line
source("Hosmer.test_etc.R")

dat <- subset(data,select=predictors,PrevalentHFAIL==0&!is.na(Event))
var=predictors

res <- coxph(Surv(Event_time, Event)~.,data=dat)

pred <- Coxar(res,years=16)
test <- HosLem.test(res$y,pred,plot=FALSE)
test$pval
