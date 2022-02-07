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
  write.csv(res, file=file,quote=FALSE,row.names=FALSE)

}


# Write predictions to a file
tmp <- test_algorithm(S1, O1, S2, O2, file = "scores.csv")

