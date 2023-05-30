# load data
library(readr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggupset)
library(tidyverse)
library(patchwork)
# set parameters
mainDir <- getwd()
PARAM <- list()
PARAM$folder.R <- paste0(getwd(), "/")
PARAM$folder.result <- paste0(PARAM$folder.R, "output/")
PARAM$folder.data <- paste0(PARAM$folder.R, "final_results")
PARAM$folder.synthetic <- paste0(PARAM$folder.R) #change here accordingly
PARAM$folder.finrisk <- paste0(PARAM$folder.R) #change here accordingly

bootstrap_Harrell_c <- read_csv(paste0(PARAM$folder.data, "/bootstrap_Harrell_c_baseline3.csv"))[,-1]

bootstrap_hoslem <- read_csv(paste0(PARAM$folder.data,"/bootstrap_hoslem_baseline3.csv"))[,-1]

# prepare dataset
df.harrelc=as.data.frame(t(bootstrap_Harrell_c))
df.harrelc$Teams_Participants <- rownames(df.harrelc)
df.harrelc=melt(df.harrelc)
colnames(df.harrelc)  <- c("Teams/Participants", 
                           "variable", 
                           "Bootstrapped Harrell's C Index (N=1000)")

df.hoslem = as.data.frame(t(bootstrap_hoslem))
df.hoslem$Teams_Participants <- rownames(df.hoslem)
df.hoslem=melt(df.hoslem)
colnames(df.hoslem)  <- c("Teams/Participants", 
                          "variable", 
                          "Bootstrapped Hosmer-Lemeshow Test (N=1000)")

# winner model barplot
# load data
Final_leaderboard <- as.data.frame(read_csv(paste0(PARAM$folder.data, "/Final_leaderboard.csv")))

Final_leaderboard$`Teams/Participants` <- c("Metformin_121_6", "Yuanfang_Guan_and_Hanrui_Zhang", 
                                            "DenverFINRISKHacky","TristanF","SB2", "Pasolli", 
                                            "UTK_Bioinformatics","Baseline Model Age-Sex",
                                            "Baseline Model All Covariates", "Baseline Total")
Final_leaderboard$harrell_c.r <- round(Final_leaderboard$harrell_c, digits = 3)
library(scales)
Final_leaderboard$hoslem_sci<- scientific_format()(Final_leaderboard$hoslem_test)

#################################################################################
pharrelc.bar <- ggplot(Final_leaderboard, aes(x=`Teams/Participants`, 
                                              y=harrell_c,  
                                              color=`Teams/Participants`,
                                              fill=`Teams/Participants`,
                                              label=`Teams/Participants`)) +
  geom_segment(aes(y = 0, 
                   x = `Teams/Participants`, 
                   yend = harrell_c, 
                   xend = `Teams/Participants`), 
                   color = "grey") +
  geom_point(stat='identity',  size=12)  +
  geom_text(aes(label=harrell_c.r), color="black", size=3) +
  labs(y ="Harrell's C Index") +
  theme_bw() +
  theme(#axis.title.y=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank(),
    axis.text.x = element_text(color = "grey20", size = 10, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(color = "grey20", size = 10,  hjust = 1, vjust = 0, face = "plain"),
    axis.title.x = element_text(color = "grey20", size = 12, hjust = .5, vjust = 0, face = "bold"),
    axis.title.y = element_text(color = "grey20", size = 12, hjust = .5, vjust = .5, face = "bold"),
    legend.position="none") +  
  scale_x_discrete(limits=c("Pasolli", "UTK_Bioinformatics","Metformin_121_6",
                            "Yuanfang_Guan_and_Hanrui_Zhang",
                            "TristanF","Baseline Total", "Baseline Model All Covariates",
                            "Baseline Model Age-Sex", "DenverFINRISKHacky", "SB2"),
                   labels=c("Team 7 (9731666)", "Team 6 (9731713)","Team 5 (9731454)",
                            "Team 4 (9731490)","Team 3 (9731636)", "Baseline Total", 
                            "Baseline Model All Covariates",
                            "Baseline Model Age-Sex", "DenverFINRISKHacky", "SB2")) +
  scale_color_manual(values=c("#7E57C2", "#7E57C2", "#7E57C2", "#1976D2", "#bdbdbd", 
                              "#bdbdbd", "#009688", "#bdbdbd", "#bdbdbd", "#bdbdbd")) +
  scale_fill_manual(values=c("#7E57C2", "#7E57C2", "#7E57C2", "#1976D2", "#bdbdbd", 
                             "#bdbdbd", "#009688", "#bdbdbd", "#bdbdbd", "#bdbdbd")) +
  rremove("x.grid") +
  coord_flip() + 
  ggtitle("A")

pharrelc.bar
#################################################################################
phoslem.bar <- ggplot(Final_leaderboard, aes(x=`Teams/Participants`, 
                                             y=hoslem_test,  
                                             color=`Teams/Participants`,
                                             fill=`Teams/Participants`,
                                             label=`Teams/Participants`)) +
  geom_segment(aes(y = 0, 
                   x = `Teams/Participants`, 
                   yend = hoslem_test, 
                   xend = `Teams/Participants`), 
               color = "grey") +
  geom_point(stat='identity',  size=12) +
  geom_text(aes(label=hoslem_sci), color="black", size=3) +
  scale_y_log10()+
  theme_bw() +   
  labs(y ="Hosmer-Lemeshow Test") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(color = "grey20", size = 10, hjust = .5, vjust = .5, face = "plain"),
        #axis.text.y = element_text(color = "grey20", size = 10,  hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size =12, hjust = .5, vjust = 0, face = "bold"),
        #axis.title.y = element_text(color = "grey20", size =12, hjust = .5, vjust = .5, face = "bold"),
        legend.position="none") +  
  scale_x_discrete(limits=c("Pasolli", "UTK_Bioinformatics","Metformin_121_6",
                            "Yuanfang_Guan_and_Hanrui_Zhang",
                            "TristanF","Baseline Total", "Baseline Model All Covariates",
                            "Baseline Model Age-Sex", "DenverFINRISKHacky", "SB2"),
                   labels=c("Team 7 (9731666)", "Team 6 (9731713)","Team 5 (9731454)",
                            "Team 4 (9731490)","Team 3 (9731636)", "Baseline Total", 
                            "Baseline Model All Covariates",
                            "Baseline Model Age-Sex", "DenverFINRISKHacky", "SB2")) +
  scale_color_manual(values=c("#7E57C2", "#7E57C2", "#7E57C2", "#1976D2", "#bdbdbd", 
                              "#bdbdbd", "#009688", "#bdbdbd", "#bdbdbd", "#bdbdbd")) +
  scale_fill_manual(values=c("#7E57C2", "#7E57C2", "#7E57C2", "#1976D2", "#bdbdbd", 
                             "#bdbdbd", "#009688", "#bdbdbd", "#bdbdbd", "#bdbdbd")) + 
  rremove("x.grid") +
  coord_flip() + 
  ggtitle("B")

phoslem.bar
#################################################################################
# plot hoslem and harrelc bootstrep
pharrelc <- ggplot(df.harrelc, aes(x=`Teams/Participants`, 
                                y=`Bootstrapped Harrell's C Index (N=1000)`,  
                                color=`Teams/Participants`,
                                fill=`Teams/Participants`)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.2, fill="white", color="black")   +  
  scale_color_manual(values=c("#7E57C2", "#7E57C2", "#7E57C2", "#1976D2", "#bdbdbd", 
                              "#bdbdbd", "#009688", "#bdbdbd", "#bdbdbd", "#bdbdbd")) +
  scale_fill_manual(values=c("#7E57C2", "#7E57C2", "#7E57C2", "#1976D2", "#bdbdbd", 
                             "#bdbdbd", "#009688", "#bdbdbd", "#bdbdbd", "#bdbdbd")) + #009688 4DB6ACGREEN  1976D2 BLUE 
  theme_bw() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10,  hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "grey20", size = 12, hjust = .5, vjust = .5, face = "bold"),
        legend.position="none") +    
  scale_x_discrete(limits=c("Pasolli", "UTK_Bioinformatics","Metformin_121_6",
                            "Yuanfang_Guan_and_Hanrui_Zhang",
                            "TristanF","Baseline Total", "Baseline Model AllCovariates",
                            "Baseline Model Age-Sex", "DenverFINRISKHacky", "SB2"),
                   labels=c("Team 7 (9731666)", "Team 6 (9731713)","Team 5 (9731454)",
                            "Team 4 (9731490)","Team 3 (9731636)", "Baseline Total", "Baseline Model All Covariates",
                            "Baseline Model Age-Sex", "DenverFINRISKHacky", "SB2")) +
  coord_flip() + 
  ggtitle("C")
pharrelc
#################################################################################
phoslem <- ggplot(df.hoslem, aes(x=`Teams/Participants`, 
                                   y=`Bootstrapped Hosmer-Lemeshow Test (N=1000)`,  
                                   color=`Teams/Participants`)) +
  geom_boxplot()   +
  scale_color_manual(values=c("#7E57C2", "#7E57C2", "#7E57C2", "#1976D2", "#bdbdbd", 
                              "#bdbdbd", "#009688", "#bdbdbd", "#bdbdbd", "#bdbdbd")) +
  theme_bw() + rremove("x.grid") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(color = "grey20", size = 10, hjust = .5, vjust = .5, face = "plain"),
        #axis.text.y = element_text(color = "grey20", size = 10,  hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size =12, hjust = .5, vjust = 0, face = "bold"),
        #axis.title.y = element_text(color = "grey20", size =12, hjust = .5, vjust = .5, face = "bold"),
        legend.position="none") +  
  scale_x_discrete(limits=c("Pasolli", "UTK_Bioinformatics","Metformin_121_6",
                            "Yuanfang_Guan_and_Hanrui_Zhang",
                            "TristanF","Baseline Total", "Baseline Model AllCovariates",
                            "Baseline Model Age-Sex", "DenverFINRISKHacky", "SB2"),
                   labels=c("Team 7 (9731666)", "Team 6 (9731713)","Team 5 (9731454)",
                            "Team 4 (9731490)","Team 3 (9731636)", "Baseline Total", "Baseline Model AllCovariates",
                            "Baseline Model Age-Sex", "DenverFINRISKHacky", "SB2")) +
  coord_flip() + ggtitle("D")
phoslem
#################################################################################
figure2 <- pharrelc.bar + phoslem.bar + pharrelc+ phoslem + 
                     plot_layout(widths = c(1, 1))

figure2
ggsave(figure2, filename = paste0(PARAM$folder.result, 
                                  "Figure2_harrellc_hoslem.pdf"),
       width = 15, height = 10)

#################################################################################
# pair comparisons
#################################################################################

pair_model <- as.data.frame(read_csv(paste0(PARAM$folder.data, "/pair_average_eval.csv")))
pair_model$harrell_c.r <- round(pair_model$harrell_c, digits = 3)
pair_model=as_tibble(pair_model)
# add matches
list_name <-list(c("SB2", "DenverFINRISKHacky"), 
     c("SB2", "DenverFINRISKHacky", "Yuanfang_Guan_and_Hanrui_Zhang"),
     c("SB2", "DenverFINRISKHacky", "Metformin_121_6"), 
     c("SB2", "DenverFINRISKHacky","TristanF"),
     c("SB2", "DenverFINRISKHacky", "Yuanfang_Guan_and_Hanrui_Zhang", "Metformin_121_6"),
     c("SB2", "DenverFINRISKHacky", "Yuanfang_Guan_and_Hanrui_Zhang", "Metformin_121_6", "TristanF"),
     c("SB2", "DenverFINRISKHacky", "Yuanfang_Guan_and_Hanrui_Zhang", "Metformin_121_6", 
       "TristanF", "UTK_Bioinformatics"),
     c("SB2", "DenverFINRISKHacky", "Yuanfang_Guan_and_Hanrui_Zhang", "Metformin_121_6", 
       "TristanF", "UTK_Bioinformatics", "Pasolli"),
     c("SB2"), c("DenverFINRISKHacky"))

pair_model$teams <- list_name
values=c("#1976D2", "#009688", "#FFDB6D", "#C4961A", "#CC79A7", 
         "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
# plot
p.pairs <- pair_model %>%
  distinct(harrell_c, .keep_all=TRUE) %>%
  arrange(order) %>%
  ggplot(aes(x=teams, y=harrell_c, color=`Team Pairs`)) +
  geom_path(group=1, color="grey")+
  geom_point(alpha=0.8, color=values, size=7)+
  scale_x_upset(n_intersections = 20, order_by = "degree") +
  labs(y="Harrell's C Index") +
  theme_bw() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10,  face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, face = "bold"),
        axis.title.y = element_text(color = "grey20", size = 12, face = "bold"),
        legend.position="left") +  
  theme_combmatrix(combmatrix.panel.point.color.fill = "grey20",
                   combmatrix.panel.line.size = 1,
                   combmatrix.label.make_space = FALSE,
                   combmatrix.panel.point.size = 5)

p.pairs
ggsave(p.pairs, filename = paste0(PARAM$folder.result, 
                                  "Figure_pairs.pdf"),
       width = 15, height = 5)

###############################################################################
# comparing synthetic vs real in different groups
###############################################################################
# finrisk
S.test.F <- read.csv(file = paste0(PARAM$folder.finrisk, 
                                 "test/pheno_test.csv"), 
                   row.names=1,)
S.train.F <- read.csv(file = paste0(PARAM$folder.finrisk, 
                                  "train/pheno_training.csv"),
                    row.names=1,)
S.score.F <- read.csv(file = paste0(PARAM$folder.finrisk, 
                                "score/pheno_validation.csv"), 
                  row.names=1,)
S.test.F$data <- "FINRISK"
S.train.F$data <- "FINRISK"
S.score.F$data <- "FINRISK"

S.test.F$Group <- "Test"
S.train.F$Group <- "Train"
S.score.F$Group <- "Score"
# synthetic
S.test.S <- read.csv(file = paste0(PARAM$folder.synthetic, 
                                   "test/pheno_test.csv"), 
                     row.names=1,)
S.train.S <- read.csv(file = paste0(PARAM$folder.synthetic, 
                                    "train/pheno_training.csv"),
                      row.names=1,)
S.score.S <- read.csv(file = paste0(PARAM$folder.synthetic, 
                                  "score/pheno_validation.csv"), 
                    row.names=1,)

S.test.S$data <- "Synthetic"
S.train.S$data <- "Synthetic"
S.score.S$data <- "Synthetic"

S.test.S$Group <- "Test"
S.train.S$Group <- "Train"
S.score.S$Group <- "Score"

# bind all tables
phe.table <-rbind(S.test.S, S.train.S, S.score.S, S.test.F, S.train.F, S.score.F)
df=melt(phe.table)
###############################################################################
# plot categorical variables
library("viridis")
p.bar <- phe.table %>%
  select(Group, data, Sex, Smoking, BPTreatment, PrevalentCHD,
         PrevalentHFAIL, PrevalentDiabetes, Event) %>%
  melt() %>%
  count(Group, data, variable, value) %>%
  mutate(Groupset = paste0(Group, " ", data)) %>%
  ggplot(aes(x = Group, 
             y = n, 
             fill = Groupset, 
             color=as.factor(value))) +
  geom_bar(stat="identity", 
           position=position_dodge()) +
  facet_wrap(variable~., scales = "free")+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  labs(y = "", 
       x="", 
       color = "Group",
       fill = "Dataset",
       title = "Comparison of FINRISK and Synthetic Daatset in subgroupos")+
  theme_classic() 


p.bar
ggsave(filename=paste0(PARAM$folder.result, "barplot_categorical.pdf"),
       plot=p.bar) 

###############################################################################
# same but density plot
p.density <- phe.table %>%
  select(Group, data, Age, Event_time, BodyMassIndex, SystolicBP, NonHDLcholesterol) %>%
  melt() %>%  
  mutate(Group = paste0(data, " ", Group)) %>%
  ggplot() + 
  #geom_histogram(aes(y=..density..), colour="black", fill="white") +
  geom_density(mapping = aes(x = value, 
                             y = ..count.., 
                             color = as.character(Group)),
                             alpha = 0.8,
                             position = "identity",
               size=1) + 
  facet_wrap(.~variable, scales = "free", nrow = 4) +
  labs(y = "Total # of Individuals", 
       color = "Dataset", 
       fill = "Dataset",
       title = "Comparing two datasets:0 for FINRISK and 1 for Synthetic") +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_classic() +  theme(legend.position="bottom") 

p.density
ggsave(filename=paste0(PARAM$folder.result, "density_continous.pdf"),
       plot=p.density) 
###############################################################################
#scatter for continous 
df = phe.table %>% pivot_wider(names_from = data, values_from = Age)

p.scatter <- df %>% 
  ggplot(aes(x=Synthetic, y=FINRISK)) + 
  geom_point(color=) +
  #facet_wrap(.~variable, scales = "free", nrow = 4) +
  labs(y = "Total # of Individuals", 
       color = "Dataset", 
       fill = "Dataset",
       title = "Comparing FINRISK and Synthetic Datasets") +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_classic() +  theme(legend.position="bottom") 

p.scatter
