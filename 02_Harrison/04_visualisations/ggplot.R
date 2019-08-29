#load packages
library("lattice")
library("RColorBrewer")
library("scales")
library("colorspace")
library("ggplot2")
library("magrittr")
library("permute")
library("vegan")
library("ggplot2")

library("scales")

#read and attach data
exon <- read.csv("life_history_variation.csv", header = TRUE) 
str(exon) #check it has loaded in properly
attach(exon)  #attach the data columns
colnames(exon)[1] <- "Sample"

#Sample C1s
c1_MotherQueen <- (colSums(exon[which(exon$Life_Stage =="Queen" & exon$Individual == "C1"), 5:length(exon)]))
c1_WorkerLarvae <- (colSums(exon[which(exon$Life_Stage =="Worker Larvae" & exon$Individual == "C1"), 5:length(exon)]))
c1_ReproductiveAdultWorker <- (colSums(exon[which(exon$Life_Stage =="Worker Rep. Adult" & exon$Individual == "C1"), 5:length(exon)]))
c1_MaleAdult <- (colSums(exon[which(exon$Life_Stage =="Male Adult" & exon$Individual == "C1"), 5:length(exon)]))
c1_AdultWorkerwithUndeterminedReproductiveStatus <- (colSums(exon[which(exon$Life_Stage =="Worker Adult with Und. RS" & exon$Individual == "C1"), 5:length(exon)]))
c1_WorkerPupae <- (colSums(exon[which(exon$Life_Stage =="Worker Pupae" & exon$Individual == "C1"), 5:length(exon)]))
c1_MaleLarvae <- (colSums(exon[which(exon$Life_Stage =="Male Larvae" & exon$Individual == "C1"), 5:length(exon)]))
c1_MalePupae <- (colSums(exon[which(exon$Life_Stage =="Male Pupae" & exon$Individual == "C1"), 5:length(exon)]))

#Sample C2
C2_MotherQueen <- (colSums(exon[which(exon$Life_Stage =="Queen" & exon$Individual == "C2"), 5:length(exon)]))
C2_WorkerLarvae <- (colSums(exon[which(exon$Life_Stage =="Worker Larvae" & exon$Individual == "C2"), 5:length(exon)]))
C2_ReproductiveAdultWorker <- (colSums(exon[which(exon$Life_Stage =="Worker Rep. Adult" & exon$Individual == "C2"), 5:length(exon)]))
C2_MaleAdult <- (colSums(exon[which(exon$Life_Stage =="Male Adult" & exon$Individual == "C2"), 5:length(exon)]))
C2_AdultWorkerwithUndeterminedReproductiveStatus <- (colSums(exon[which(exon$Life_Stage =="Worker Adult with Und. RS" & exon$Individual == "C2"), 5:length(exon)]))
C2_WorkerPupae <- (colSums(exon[which(exon$Life_Stage =="Worker Pupae" & exon$Individual == "C2"), 5:length(exon)]))
C2_MaleLarvae <- (colSums(exon[which(exon$Life_Stage =="Male Larvae" & exon$Individual == "C2"), 5:length(exon)]))
C2_MalePupae <- (colSums(exon[which(exon$Life_Stage =="Male Pupae" & exon$Individual == "C2"), 5:length(exon)]))


#Sample C3
C3_MotherQueen <- (colSums(exon[which(exon$Life_Stage =="Queen" & exon$Individual == "C3"), 5:length(exon)]))
C3_WorkerLarvae <- (colSums(exon[which(exon$Life_Stage =="Worker Larvae" & exon$Individual == "C3"), 5:length(exon)]))
C3_ReproductiveAdultWorker <- (colSums(exon[which(exon$Life_Stage =="Worker Rep. Adult" & exon$Individual == "C3"), 5:length(exon)]))
C3_MaleAdult <- (colSums(exon[which(exon$Life_Stage =="Male Adult" & exon$Individual == "C3"), 5:length(exon)]))
C3_AdultWorkerwithUndeterminedReproductiveStatus <- (colSums(exon[which(exon$Life_Stage =="Worker Adult with Und. RS" & exon$Individual == "C3"), 5:length(exon)]))
C3_WorkerPupae <- (colSums(exon[which(exon$Life_Stage =="Worker Pupae" & exon$Individual == "C3"), 5:length(exon)]))
C3_MaleLarvae <- (colSums(exon[which(exon$Life_Stage =="Male Larvae" & exon$Individual == "C3"), 5:length(exon)]))
C3_MalePupae <- (colSums(exon[which(exon$Life_Stage =="Male Pupae" & exon$Individual == "C3"), 5:length(exon)]))


#create data frame
totals_by_caste <- data.frame(Subunits = rep(as.character(names(exon[, 5:length(exon)])), 
                                             times = 8),
                              Development = rep(c("Queen",
                                                  "Worker Larvae",
                                                  "Worker Rep. Adult",
                                                  "Male Adult",
                                                  "Worker Adult with Und. RS", 
                                                  "Worker Pupae", 
                                                  "Male Larvae", 
                                                  "Male Pupae"), each = 30), 
                              Count = c(c1_MotherQueen, C2_MotherQueen, 
                                        C3_MotherQueen, c1_WorkerLarvae,
                                        C2_WorkerLarvae, C3_WorkerLarvae,
                                        c1_ReproductiveAdultWorker, 
                                        C2_ReproductiveAdultWorker, 
                                        C3_ReproductiveAdultWorker, 
                                        c1_MaleAdult, C2_MaleAdult, 
                                        C3_MaleAdult,
                                        c1_AdultWorkerwithUndeterminedReproductiveStatus, 
                                        C2_AdultWorkerwithUndeterminedReproductiveStatus, 
                                        C3_AdultWorkerwithUndeterminedReproductiveStatus, 
                                        c1_WorkerPupae, C2_WorkerPupae, 
                                        C3_WorkerPupae, c1_MaleLarvae, 
                                        C2_MaleLarvae, C3_MaleLarvae, 
                                        c1_MalePupae, C2_MalePupae, C3_MalePupae))


#plot line graph of overall expression of subunit use across life stages, using mean, upper and lower counts for wisker boxplots


#####Figure 4
# create a vector for each subunit
subunit_vec <- totals_by_caste$Subunits

# create a vector for each stage
stage_vec <- totals_by_caste$Development

# create a result dataframe (empty now, will be filled by the loop)
result_df <- as.data.frame(matrix(NA, ncol = 5))

#name dataframe columns
colnames(result_df) <- c ("subunit", "stage", "count_mean", "lower_count", "upper_count")

#loop through criteria
for(stage_position in 1:length(stage_vec)){
  for(subunit_position in 1:length(subunit_vec)){
    #subset for stage and subunit
    totals_by_caste_subset <- subset(totals_by_caste, subset = Development == stage_vec[stage_position] & Subunits == subunit_vec[subunit_position])
    #obtain count mean
    count_mean <- mean(totals_by_caste_subset$Count)
    #obtain lower count
    lower_count <- min(totals_by_caste_subset$Count)
    #obtain higher count
    higher_count <- max(totals_by_caste_subset$Count)
    #combine all into vector
    result_vec <- c(subunit_vec[subunit_position],
                    stage_vec[stage_position],
                    count_mean,
                    lower_count,
                    higher_count)
    #add to the results df
    result_df <- rbind(result_df, result_vec)
  }
}

# remove MA
result_df <- result_df[complete.cases(result_df), ]

# change character for numeric
result_df$count_mean <- as.numeric(result_df$count_mean)
result_df$lower_count <- as.numeric(result_df$lower_count)
result_df$upper_count <- as.numeric(result_df$upper_count)

#check order of stages
stage_vec[stage_position]

#rename stages
result_df$stage <- replace(as.character(result_df$stage), result_df$stage == "1", "Male Adult")
result_df$stage <- replace(as.character(result_df$stage), result_df$stage == "2", "Male Larvae")
result_df$stage <- replace(as.character(result_df$stage), result_df$stage == "3", "Male Pupae")
result_df$stage <- replace(as.character(result_df$stage), result_df$stage == "4", "Queen")
result_df$stage <- replace(as.character(result_df$stage), result_df$stage == "5", "Worker Adult with Und. RS.")
result_df$stage <- replace(as.character(result_df$stage), result_df$stage == "6", "Worker Larvae")
result_df$stage <- replace(as.character(result_df$stage), result_df$stage == "7", "Worker Pupae")
result_df$stage <- replace(as.character(result_df$stage), result_df$stage == "8", "Worker Rep. Adult")

#rename subunits
result_df$subunit <- replace(as.character(result_df$subunit), result_df$subunit == "1", "alpha1")
result_df$subunit <- replace(as.character(result_df$subunit), result_df$subunit == "2", "alpha2")
result_df$subunit <- replace(as.character(result_df$subunit), result_df$subunit == "3", "alpha3")
result_df$subunit <- replace(as.character(result_df$subunit), result_df$subunit == "4", "alpha4")
result_df$subunit <- replace(as.character(result_df$subunit), result_df$subunit == "5", "alpha5")
result_df$subunit <- replace(as.character(result_df$subunit), result_df$subunit == "6", "alpha6")
result_df$subunit <- replace(as.character(result_df$subunit), result_df$subunit == "7", "alpha7")
result_df$subunit <- replace(as.character(result_df$subunit), result_df$subunit == "8", "alpha8")
result_df$subunit <- replace(as.character(result_df$subunit), result_df$subunit == "9", "alpha9beta2")
result_df$subunit <- replace(as.character(result_df$subunit), result_df$subunit == "10", "beta1")

# check result df (uncomment to check)
head(result_df)

# plot results
p <- ggplot(result_df, aes(x = stage, y = count_mean, Fill = subunit, group = subunit, colour = subunit))  + 
  scale_colour_manual(values=c("darkred", "red", "darkorange", "purple", "blue", "darkblue", "pink", "turquoise", "limegreen", "darkgreen")) + 
  geom_errorbar(data = result_df,
                mapping = aes(x = stage, ymin = upper_count, ymax = lower_count),
                width = 0.2,
                size = 1) +
  geom_point(data = result_df,
             mapping = aes(x = stage, y = count_mean),
             size = 4,
             shape = 21,
             fill = "white") +
  geom_line(data = result_df,
            mapping = aes(x = stage, y = count_mean)) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, face = "italic")) + 
  ylab("Mean Exon Read Counts (/log)") + 
  scale_y_continuous(trans=log2_trans()) 
  
print(p + ggtitle("Overall Subunit Expression in Different Life Stages in Bombus Terrestris"))


#plot boxplots of subunit expression at each life stage
##Figure 5
#plot data
p <- ggplot(totals_by_caste, aes(x = Subunits, y = Count, fill = Development)) +
  geom_boxplot(width=1) + theme_bw() + theme(axis.text.x = element_text(angle = -90, hjust = 0, face = "italic")) 
  + ylab("RNA Exon Read Counts (/log)")

p <- p + scale_y_continuous(trans=log2_trans())

p <- p + scale_fill_manual(values=c("black", "white", "skyblue", "limegreen", "purple", "yellow", "blue", "red",
                                    "brown", "electricblue"))

print(p + ggtitle("nAChR Subunit Usage in Different Life Stages in Bombus Terrestris"))

#test for significance
model2 <- glm(totals_by_caste$Count ~ totals_by_caste$Subunits*totals_by_caste$Development, family = "quasipoisson", data = exon)
summary(model2)


