test_algorithm <- function (S1, O1, S2, O2, file) {
S1=S.train
S2=S.test
O1=O.train[1:20,]
O2=O.test[1:20,]
  # Prepare features by combining phenodata and abundance data; add some filtering
  D1 <- apply(cbind(S1, t(O1)), 2, as.numeric) # Train
  D2 <- apply(cbind(S2, t(O2)), 2, as.numeric) # Test

  set.seed(36636333)
  endpoints <- c("Event_time", "Event")
  covariates <- setdiff(colnames(S1), c("SampleID", endpoints))
  taxa <- sample(rownames(O1), 10)
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
  scores <- pred2$predicted; names(scores) <- rownames(dff2)[samples]

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















surv=univ_models$Sex
pred=Coxar(surv, years=10)
test <- HosLem.test(surv, pred, plot=FALSE)


HosLem.test <- function(surv,pred,plot=TRUE,DF.reduce=2) {
  # Cook-Ridker version test uses DF.reduce=2 (default)
  # Cook NR, Ridker PM. (2009) Ann Intern Med. 150(11): 795-802
  # D'Agostino-Nam version uses DF.reduce=1
  # About the differences: http://hdl.handle.net/1773/22648
  # (Guffey D., Hosmer-Lemeshow goodness-of-fit test: Translations to the Cox Proportional Hazards Model)
  # if plot is a name, a jpg plot is saved with this name (.jpg will be added)
  if(!DF.reduce%in%c(1,2)) stop('Please specify DF.reduce = 1 or 2')
  if(!is.Surv(surv$y)) stop('Please use a survival object')
  version <- c('D\'Agostino-Nam','Cook-Ridker')[DF.reduce]
  if(names(surv$means) %in% list.num){
    grp <- as.numeric(cut2(pred, g=10))
    pj <- as.numeric(by(pred,grp,mean))
    idx <- TRUE
  }
  if(!(names(surv$means) %in% list.num)){
    grp <- as.numeric(cut2(pred, g=10))
    pj <- as.numeric(by(pred,grp,mean))
    idx <- TRUE
  }
  while(any(idx)&length(pj)>3) {
    kmj <- 1-grpkm(grp,surv$y)[,3] # failure, not survival
    pj <- as.numeric(by(pred,grp,mean))
    nj <- table(grp)
    idx <- grp%in%which(nj<5|pj*nj<5|nj*kmj<5)&grp<max(grp)
    if(any(idx)) {
      grp[idx] <- grp[idx]+1
      grp <- match(grp,sort(unique(grp)))
    }
  }
  dan_chi2 <- sum(nj*(kmj-pj)^2/(pj*(1-pj)))
  pval <- pchisq(dan_chi2,length(pj)-DF.reduce,lower.tail=FALSE)
  plotname <- ''
  if(plot!=FALSE) {
    if(is.character(plot)) { plot <- paste0(plot,'.jpg'); jpeg(filename=plot); plotname <- plot }
    barplot(height=rbind(pj*100,kmj*100),beside=T,names.arg=1:length(pj),
            main=paste('Hosmer-Lemeshow (',version,') test\np-value =',format.pval(pval,digits=4),sep=''),
            xlab=xlab <- paste('Risk', tiles[length(pj)]),ylab='Risk, %',col=c('black','lightgray'))
    if(is.character(plot)) dev.off()
  }
  list(version=version,quantiles=tiles[length(pj)],nj=nj,kmj=kmj,pj=pj,chi2=dan_chi2,pval=pval,plotname=plotname)
}


