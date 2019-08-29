#list of files to transform - trasform into right format for analysis
list <- c("alpha_1_exons.csv", "alpha_2_exons.csv", "alpha_3_exons.csv", "alpha_4_exons.csv", "alpha_5_exons.csv", "alpha_6_exons.csv",
"alpha_7_exons.csv", "alpha_8_exons.csv", "alpha_9_exons.csv","beta_1_exons.csv", "beta_2_exons.csv")
for (i in list){
  exon <- read.csv("i", header = TRUE) 
  str(exon) #check it has loaded in properly
  attach(exon)  #attach the data columns
  exon <- t(exon)
  colnames(exon) <- exon[1, ]
  exon <- exon[-1, ]
  exon<-as.data.frame(exon)
  write.csv(exon, paste(i,".csv"))
}

  rm(list = ls())
  
  
  #load packages
  library("lattice")
library("RColorBrewer")
library("scales")
library("colorspace"))
library("rlang")
library("ggplot2")
library("magrittr")
library("permute")
library("vegan")
library("ggplot2")

#create dataframe and create plot
#read and attach data
exon <- read.csv("apis_variation.csv", header = TRUE) 
str(exon) #check it has loaded in properly
attach(exon)  #attach the data columns
colnames(exon)[1] <- "Sample"

#Sample C1s
d1 <- log(colSums(exon[which(exon$Sample =="DMSO_1"), 3:length(exon)])+1)
d2 <- log(colSums(exon[which(exon$Sample =="DMSO_2"), 3:length(exon)])+1)
d3 <- log(colSums(exon[which(exon$Sample =="DMSO_3"), 3:length(exon)])+1)
d4 <- log(colSums(exon[which(exon$Sample =="DMSO_4"), 3:length(exon)])+1)
d5 <- log(colSums(exon[which(exon$Sample =="DMSO_5"), 3:length(exon)])+1)

#Sample C2
s1 <- log(colSums(exon[which(exon$Sample =="SUGAR_1"), 3:length(exon)])+1)
s2 <- log(colSums(exon[which(exon$Sample =="SUGAR_2"), 3:length(exon)])+1)
s3 <- log(colSums(exon[which(exon$Sample =="SUGAR_3"), 3:length(exon)])+1)
s4 <- log(colSums(exon[which(exon$Sample =="SUGAR_4"), 3:length(exon)])+1)
s5 <- log(colSums(exon[which(exon$Sample =="SUGAR_5"), 3:length(exon)])+1)

#Sample C3
C31 <- log(colSums(exon[which(exon$Sample =="CLO_3_1"), 3:length(exon)])+1)
C32 <- log(colSums(exon[which(exon$Sample =="CLO_3_2"), 3:length(exon)])+1)
C33 <- log(colSums(exon[which(exon$Sample =="CLO_3_3"), 3:length(exon)])+1)
C34 <- log(colSums(exon[which(exon$Sample =="CLO_3_4"), 3:length(exon)])+1)
C35 <- log(colSums(exon[which(exon$Sample =="CLO_3_5"), 3:length(exon)])+1)

#Sample C3
I31 <- log(colSums(exon[which(exon$Sample =="IMD_3_1"), 3:length(exon)])+1)
I32 <- log(colSums(exon[which(exon$Sample =="IMD_3_2"), 3:length(exon)])+1)
I33 <- log(colSums(exon[which(exon$Sample =="IMD_3_3"), 3:length(exon)])+1)
I34 <- log(colSums(exon[which(exon$Sample =="IMD_3_4"), 3:length(exon)])+1)
I35 <- log(colSums(exon[which(exon$Sample =="IMD_3_5"), 3:length(exon)])+1)

#Sample C3
I031 <- log(colSums(exon[which(exon$Sample =="IMD_0_3_1"), 3:length(exon)])+1)
I032 <- log(colSums(exon[which(exon$Sample =="IMD_0_3_2"), 3:length(exon)])+1)
I033 <- log(colSums(exon[which(exon$Sample =="IMD_0_3_3"), 3:length(exon)])+1)
I034 <- log(colSums(exon[which(exon$Sample =="IMD_0_3_4"), 3:length(exon)])+1)
I035 <- log(colSums(exon[which(exon$Sample =="IMD_0_3_5"), 3:length(exon)])+1)


#create data frame
totals_by_caste <- data.frame(Subunits = rep(as.character(names(exon[, 3:length(exon)])), times = 5),
                              Treatment = rep(c("DMSO", "SUGAR", "CLO_3", "IMD_3",
                                                  "IMD_0_3"), each = 30), 
                              Count = c(d1, d2, d3, d4, d5, 
                                        s1, s2, s3, s4, s5, 
                                        C31,C32, C33, C34, C35, 
                                        I31, I32, I33, I34, I35, 
                                        I031, I032, I033, I034, I035))

#plot data
p <- ggplot(totals_by_caste, aes(x = Subunits, y = Count, fill = Treatment)) +
  geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = -90, hjust = 0, face = "italic")) 
  + ylab("RNA Exon Read Counts (/log)")

print(p + ggtitle("nAChR Receptor Subunit Usage in Pesticide Treated Apis Mellifera"))
