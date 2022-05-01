# load libraries and 
library(tidyverse)
library(gtools)
library(Hmisc)
library(corrplot)
library(reshape2)
library(car)
library(vegan)
require(foreach)
require(bigmemory)
require(plyr)
library(doMC)

source("https://raw.githubusercontent.com/defleury/Toolbox_16S/master/R/function.alpha_diversity.R")
source("https://raw.githubusercontent.com/defleury/Toolbox_16S/master/R/function.rarefaction.R")

# This document contains the code to compare the INN results with real Finrisk data 
set.seed(594)# load the data
PARAM <- list()
setwd("/Users/ecekartal/Documents/SAEZ/projects/DREAM_HF_microbiome/")
PARAM$folder.R <- paste0(getwd(), "/")
PARAM$folder.finrisk.data <- paste0(PARAM$folder.R, "data/test/")
PARAM$folder.simulated.data <- paste0(PARAM$folder.R, "data/train/")

# Finrisk cohort data
S1 <- read.csv(file = paste0(PARAM$folder.R, "data/train/phenodata.csv"), row.names=1, )
O1 <- read.csv(file = paste0(PARAM$folder.R, "data/train/readcounts.csv"), row.names=1)

# INN-simulated data example
S2 <- read.csv(file = paste0(PARAM$folder.R, "data/test/phenodata.csv"), row.names=1)
O2 <- read.csv(file = paste0(PARAM$folder.R, "data/test/readcounts.csv"), row.names=1)

str(O1)
str(O2)
# both have a non-numeric value "_no_DNA_stage (Viruses)
which(stringr::str_starts(O1$Sample5854, "_no"  ))
which(stringr::str_starts(O2$Sample5916, "_no"  ))
# assign Na to those cells
O1$Sample5854[4774] = 0
O2$Sample5916[4774] = NA
# convert to numeric
O1$Sample5854=as.integer(O1$Sample5854)
O2$Sample5916=as.integer(O2$Sample5916)
toremove = c("ssRNA_positive-strand_viruses", "1", "0")
O1 <- O1[!rownames(O1) %in% toremove, ]
O2 <- O2[!rownames(O2) %in% toremove, ]
# bind otu tables
otu.table <-cbind(O1, O2)

# bind 2 meta datasets
S1$dataset<- paste0(0)
S2$dataset<- paste0(1)
S1=S1[,-c(7:8)] ##REMOVE THISIF YOU HAVE EVENT INFO FOR BOTH DATASETS
metatotal<- rbind(S1, S2)
metatotal$dataset <- as.numeric(metatotal$dataset)
###############################################################################
# subset metadata
idx = c( "Age", "BMI","Smoking", "BPtreatment", "PrevalentDiabetes", "PrevalentCHD",
        # "Event", "Event_time", 
         "SystolicBP", "NonHDLcholesterol", "Sex", "dataset")

testlist <-metatotal[idx]
# Correlation matrix with significance levels (p-value)
res2 <- rcorr(as.matrix(testlist))

# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
corr= flattenCorrMatrix(res2$r, res2$P)
# Insignificant correlation are crossed
# Positive correlations are displayed in blue and negative correlations in red color. 
corrplot(res2$r, type="upper",p.mat = res2$P, sig.level = 0.01)
###############################################################################
df=melt(metatotal)

p.bar <- df %>%
  ggplot() + 
  geom_histogram(mapping = aes(x = value, y = ..count.., 
                               color = dataset, fill = dataset), 
                 #stat = "binwidth",
                 alpha = 0.7,
                 position = "identity") + 
  facet_wrap(variable~., scales = "free")+
  labs(y = "Count", color = "Dataset", 
       fill = "Dataset",
       title = "Comparison of two datasets:0 for Finrisk and 1 for INN-simulated")
ggsave(filename=paste0(PARAM$folder.R, "barplot.pdf"), plot=p.bar) 

# same but density plot
p.density <- df %>%
  ggplot() + 
  geom_density(mapping = aes(x = value, y = ..count.., 
                               color = dataset, fill = dataset), 
               alpha = 0.5,
               position = "identity") + 
  facet_wrap(variable~., scales = "free")+
  labs(y = "Count", color = "Dataset", 
       fill = "Dataset",
       title = "Comparing two datasets")
ggsave(filename=paste0(PARAM$folder.R, "density.pdf"), plot=p.density) 

###############################################################################
# Calculate alpha diversities (Hill Diversities) for count tables rarefied to 1000 
#  For each sample, diversities at Hill numbers 
# q=0 (richness), 
# q=1 (exp(Shannon))  
# q=2 (inverse Simpson) 
# Hill-based evenness for q=1  
# Hill-based evenness for q=2 are computed.
quantile(colSums(O1), seq(0,1,0.05))
quantile(colSums(O2), seq(0,1,0.05))

calculate.hill <- function(df.rar, size.n){
  Hill <- Hill_Diversity.rarefied(df.rar, size=size.n, iterations=100, q.H=c(0, 1, 2))
  # get some statistics
  print(summary(Hill$q.0))
  print(summary(Hill$q.1))
  print(summary(Hill$q.2))
  colnames(Hill) = c("Sample", "Richness", "exp(Shannon)", "inv(Simpson)") 
  Hill[, "evenness (q1)"] <- Hill$`exp(Shannon)`/ Hill$Richness
  return(Hill)
}
Hill<- calculate.hill(otu.table, 5000)
# combine meta data and Hill
Hill <- merge(Hill, metatotal, by.x="Sample", by.y="row.names", all.x=TRUE)
##############################################################################
collect.confounders <- data.frame()
for (i in idx) {
  # calculate anova  
  lm<- lm(
    substitute(Richness~ as.factor(i), list(i = as.name(i))),
    data = Hill, na.action=na.omit)
  aov <- Anova(lm) %>% broom::tidy() %>% 
    mutate(metric="richness")
  aov <- aov[1,]
  aov$term <- i
  # collect data
  collect.confounders <- rbind(
    collect.confounders,
    aov
  )
}
print(collect.confounders)

###############################################################################
# are there significant community-level compositional shifts between datasets
# To approach these questions, the following analyses steps are performed:
# calculate pairwise sample compositional (dis)similarities according to various 
# beta div indices
###############################################################################

#Traditional indices of community similarity
#=> formulated as distances/dissimilarities (not similarities)
community.similarity.par <- function(o.t, distance="jaccard", use.cores=detectCores()) {
  #Register cluster
  registerDoMC(cores=use.cores);
  ########################
  samples <- colnames(o.t);
  #Turn OTU table into list format
  ot.occ.list <- apply((o.t > 0), 2, which);
  #ot.count.list <- mclapply(seq(1,length(ot.occ.list)), function(i) {o.t[ot.occ.list[[i]],i]}, mc.cores=use.cores);
  ot.count.list <- alply(o.t, .margins=2, .fun=function(o.vec) {o.vec[o.vec > 0]}, .parallel=T)
  names(ot.occ.list) <- names(ot.count.list) <- samples;
  ########################
  
  ########################
  #Compute community distance matrix
  ########################
  #Classic, unweighted Jaccard
  if (distance == "jaccard") {
    cs.list <- mclapply(ot.occ.list, function(a.list) {unlist(lapply(ot.occ.list, function(b.list) {1 - (length(intersect(a.list, b.list)) / length(union(a.list, b.list)))}))}, mc.cores=use.cores);
  }
  ########################
  #Classic, weighted Jaccard
  if (distance == "jaccard.abd") {
    cs.list <- mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count); a.list <- names(a.count);
      unlist(lapply(ot.count.list, function(b.count) {ab.shared <- intersect(a.list, names(b.count)); 1 - (sum(a.count[ab.shared]) + sum(b.count[ab.shared])) / (a.N + sum(b.count))}))
    }, mc.cores=use.cores);
  }
  ########################
  #Classic, weighted Jaccard, based on fractions (more balancing for uneven sample sizes)
  if (distance == "jaccard.abd.frac") {
    cs.list <- mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count); a.list <- names(a.count); a.rel <- a.count / a.N;
      unlist(lapply(ot.count.list, function(b.count) {
        b.N <- sum(b.count); ab.shared <- intersect(a.list, names(b.count)); b.rel <- b.count / b.N;
        1 - (0.5 * (sum(a.rel[ab.shared]) + sum(b.rel[ab.shared])))
      }))
    }, mc.cores=use.cores);
  }
  ########################
  #Classic, weighted Jaccard, based on fractions (more balancing for uneven sample sizes), Chao's version
  #=> after Chao et al., 2005, Ecology Letters
  if (distance == "jaccard.abd.frac_chao") {
    cs.list <- mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count); a.list <- names(a.count); a.rel <- a.count / a.N;
      unlist(lapply(ot.count.list, function(b.count) {
        b.N <- sum(b.count); ab.shared <- intersect(a.list, names(b.count)); b.rel <- b.count / b.N;
        U <- sum(a.rel[ab.shared]); V <- sum(b.rel[ab.shared]);
        if (U == 0 & V == 0) {return(1)} else {return(1 - ((U*V) / (U + V - (U*V))))}
      }))
    }, mc.cores=use.cores);
  }
  ########################
  #Classic, weighted Jaccard, alternative formulation (Bray-Curtis-like)
  if (distance == "jaccard.abd.alt") {
    cs.list <- mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count); a.list <- names(a.count);
      unlist(lapply(ot.count.list, function(b.count) {
        b.list <- names(b.count);
        ab.shared <- intersect(a.list, b.list); ab.union <- union(a.list, b.list);
        a.freq <- b.freq <- numeric(length=length(ab.union)); names(a.freq) <- names(b.freq) <- ab.union;
        a.freq[a.list] <- a.count; b.freq[b.list] <- b.count;
        1 - (sum(pmin(a.freq, b.freq)) / sum(pmax(a.freq, b.freq)));
      }))
    }, mc.cores=use.cores);
  }
  ########################
  #Classic, weighted Jaccard, Chao's version
  #=> after Chao et al., 2005, Ecology Letters
  if (distance == "jaccard.abd.chao") {
    cs.list <- mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count); a.list <- names(a.count);
      unlist(lapply(ot.count.list, function(b.count) {
        b.list <- names(b.count); b.N <- sum(b.count);
        #Get shared taxa & union of taxa
        ab.shared <- intersect(a.list, b.list); ab.union <- sort(union(a.list, b.list));
        a.freq <- b.freq <- numeric(length=length(ab.union)); names(a.freq) <- names(b.freq) <- ab.union;
        a.freq[a.list] <- a.count; b.freq[b.list] <- b.count;
        #If no taxa are shared, return trivial (distance = 1), otherwise compute Chao distance
        if (length(ab.shared) == 0) {d.chao <- 1} else {
          #Which taxa observed in sample a are singletons in sample b and vice versa?
          a_shared.b_singl <- which(a.freq >= 1 & b.freq == 1);
          b_shared.a_singl <- which(a.freq == 1 & b.freq >= 1);
          #How many shared taxa are singletons / doubletons in samples a & b?
          f.a.singl <- length(b_shared.a_singl);
          f.a.doubl <- length(which(a.freq == 2 & b.freq >= 1));
          f.b.singl <- length(a_shared.b_singl);
          f.b.doubl <- length(which(a.freq >= 1 & b.freq == 2));
          if (f.a.doubl == 0) {f.a.doubl = 1}; if (f.b.doubl == 0) {f.b.doubl = 1};
          #Get taxa abundance estimates for samples a & b
          a.adjust <- ((b.N-1)/b.N)*f.b.singl/(2*f.b.doubl)*sum(a.freq[a_shared.b_singl] / a.N);
          a.estimate <- sum(a.freq[ab.shared] / a.N) + a.adjust;
          b.adjust <- ((a.N-1)/a.N)*f.a.singl/(2*f.a.doubl)*sum(b.freq[b_shared.a_singl] / b.N);
          b.estimate <- sum(b.freq[ab.shared] / b.N) + b.adjust;
          if (a.estimate > 1) {a.estimate = 1}; if (b.estimate > 1) {b.estimate == 1};
          d.chao <- 1 - ((a.estimate * b.estimate) / (a.estimate + b.estimate - (a.estimate * b.estimate)));
          if (d.chao < 0) {d.chao = 0}
        }
        #Return
        d.chao
      }))
    }, mc.cores=use.cores);
  }
  ########################
  #Classic Bray-Curtis dissimilarity
  if (distance == "bray_curtis") {
    cs.list <- mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count); a.list <- names(a.count);
      unlist(lapply(ot.count.list, function(b.count) {
        b.N <- sum(b.count); b.list <- names(b.count);
        ab.shared <- intersect(a.list, b.list); ab.union <- union(a.list, b.list);
        a.freq <- b.freq <- numeric(length=length(ab.union)); names(a.freq) <- names(b.freq) <- ab.union;
        a.freq[a.list] <- a.count; b.freq[b.list] <- b.count;
        1 - (2 * (sum(pmin(a.freq, b.freq)) / (a.N + b.N)));
      }))
    }, mc.cores=use.cores);
  }
  ########################
  #Classic Morisita-Horn dissimilarity
  if (distance == "morisita_horn") {
    cs.list <- mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count); a.n <- length(a.count); a.list <- names(a.count);
      a.simpson <- sum(a.count ^ 2) / (a.N ^ 2);
      unlist(lapply(ot.count.list, function(b.count) {
        #Get current summary statistics for B
        b.N <- sum(b.count); b.n <- length(b.count); b.list <- names(b.count);
        b.simpson <- sum(b.count ^ 2) / (b.N ^ 2);
        #Get current freq vectors from intersect and union of a & b
        ab.shared <- intersect(a.list, b.list); ab.union <- union(a.list, b.list);
        a.freq <- b.freq <- numeric(length=length(ab.union)); names(a.freq) <- names(b.freq) <- ab.union;
        a.freq[a.list] <- a.count; b.freq[b.list] <- b.count;
        #Get summed abundance product
        ab.prod <- sum(a.freq * b.freq);
        #Return
        1 - ((2 * ab.prod) / ((a.simpson + b.simpson) * a.N * b.N));
      }))
    }, mc.cores=use.cores);
  }
  ########################
  
  #Return
  do.call("rbind", cs.list);
}

###############################################################################
otu.table.rel <- prop.table(as.matrix(otu.table), 2)
core = 3
# First, pre-allocate a results collector list.
betadiv <- list()
# Calculate Bray-Curtis dissimilarities on non-transformed data.
# Bray-Curtis
betadiv[["Bray_Curtis"]] <- community.similarity.par(otu.table.rel, 
                                                      distance="bray_curtis", 
                                                      use.cores=core)
  
# Calculate Bray-Curtis dissimilarities on sqrt-transformed data.
betadiv[["Bray_Curtis.sqrt"]] <- community.similarity.par(sqrt(otu.table.rel), 
                                                          distance="bray_curtis", 
                                                          use.cores=core)
# Calculate the abundance-weighted Jaccard dissimilarities
# Weighted Jaccard
betadiv[["Jaccard_w"]] <- community.similarity.par(otu.table.rel, 
                                                    distance="jaccard.abd.frac", 
                                                    use.cores=core)
# convert to distance matrix the betadiv file
dist.betadiv <- lapply(betadiv, as.dist)

# Perform some data quality checks and data transformations needed to run the analyses
names(betadiv)
dim(betadiv[[1]])
# check if ids are identical 
identical(rownames(as.matrix(betadiv[["Bray_Curtis"]])), 
          rownames(as.matrix(betadiv[["Jaccard_w"]]))) #needs to be true

curr.samples <- rownames(metatotal)
curr.meta <- metatotal[row.names(metatotal) %in% colnames(Bray),]

Bray <- as.matrix(dist.betadiv[["Bray_Curtis"]])[curr.samples, curr.samples]
BrayCurtis.sqrt <- as.matrix(dist.betadiv[["Bray_Curtis.sqrt"]])[curr.samples, curr.samples]
Jaccard <- as.matrix(dist.betadiv[["Jaccard_w"]])[curr.samples, curr.samples]

################################################################################

pco.analysis <- function(analysisName, type, data.sample) {
  #analysisName =  c("Bray", "Bray.sqrt", "Jaccard_w")
  data_analysisName <- dist.betadivEleg[[type]]
  # check the order of sample id 
  identical(rownames(as.matrix(data_analysisName)), as.character(rownames(data.sample)))
  # pcoa analysis
  
  pcoa_analysisName <- pcoa(data_analysisName, rn=rownames(dist.betadivEleg))
  curr.pcoa <- pcoa_analysisName
  # Prepare data for plotting
  plot.data <- data.frame(data.sample, 
                          Axis.1=curr.pcoa$vectors[, "Axis.1"], 
                          Axis.2=curr.pcoa$vectors[, "Axis.2"], 
                          Group=interaction(data.sample$dataset))
  # Get percent variation explained
  if (curr.pcoa$correction[1] == "none") {
    plot.var_explained <- round(100*curr.pcoa$values[1:2, "Relative_eig"]/
                                  sum(curr.pcoa$values[, "Relative_eig"]), digits=1)
  } else {
    plot.var_explained <- round(100*curr.pcoa$values[1:2, "Rel_corr_eig"]/
                                  sum(curr.pcoa$values[, "Rel_corr_eig"]), digits=1)
  }
  # Get convex hulls around points
  plot.hulls <- ddply(plot.data, "Group", 
                      function(df) df[chull(df$Axis.1, df$Axis.2), ])
  # Get centroids
  plot.centroids <- ddply(plot.data, "Group", 
                          function(df) c(mean(df$Axis.1), mean(df$Axis.2)))
  curr.plot.df <- merge(plot.data, plot.centroids, by="Group")
  
  # Plot ordination
  curr.plot <- ggplot(curr.plot.df, aes_string(x="Axis.1", y="Axis.2", colour="Group")) +
    # Add group centroids
    geom_point(data=curr.plot.df, aes(x=V1, y=V2, color=Group),  shape=15, size=8, alpha=0.6) +
    # Add segments of points to centroids
    geom_segment(data=curr.plot.df, aes(x=Axis.1, y=Axis.2, xend=V1, yend=V2), alpha=0.6) +
    # Add points (per sample)
    #geom_point(aes(color = curr.plot.df$status), size = 4) +
    geom_point(size=3, alpha=0.6) +
    # set color
    scale_color_manual(values=c("#0571b0","#f4a582"))+
    # Add ellipses (at 95%)
    stat_ellipse(level=0.95, segments=101, alpha=0.5) +
    # Add convex hull around all points per group
    geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.1) +
    # Add x & y axis labels
    xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
    ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep="")) +
    # add title
    ggtitle(list(title = paste0(type, ' ', " PCoA Analysis"))) +
    theme(plot.title = element_text(size=20)) +
    theme_classic()
  
  # save the plot
  ggsave(filename=paste0(PARAM$folder.R, type,' ', " PCoA Analysis.pdf"), plot=curr.plot) 
}

################################################################################

# choose right data set according to sites!!!
data.sample <- curr.meta
rownames(data.sample) <- as.character(rownames(data.sample))
dist.betadivEleg = lapply(dist.betadiv, function(x) 
{as.dist(as.matrix(x)[rownames(data.sample), rownames(data.sample)])}) 
# herewith we extract from the list all distance pairs of eligible cases and controls

# run the function
pco.analysis(Bray, "Bray_Curtis", data.sample)
pco.analysis(Jaccard_w, "Jaccard_w",  data.sample)
pco.analysis(Bray.sqrt, "Bray_Curtis.sqrt", data.sample)
