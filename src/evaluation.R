# Read the survival scores provided by the submitted algorithm
scores <- read.csv(file=paste0(PARAM$folder.testdata,"scores.csv"))

# Read the real labels that were hidden from that algorithm
labels <- read.csv(file=paste0(PARAM$folder,"test_phenodata_labels.csv"), row.names=1)

# Align the user provided scores with the true event times
true.scores <- labels[scores$SampleID,"Event_time"]

# Calculate Harrell's C statistics
C <- Hmisc::rcorr.cens(scores$Score, true.scores, outx=FALSE)
#print(C["C Index"])
write.csv(C, file=paste0(PARAM$folder,"HarrellC.csv"))

