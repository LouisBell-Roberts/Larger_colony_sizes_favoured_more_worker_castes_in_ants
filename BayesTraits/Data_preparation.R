########################Preparing files for BayesTraits######################
####Converting .csv data sets to .txt files
###Filtering and transforming variables to correct format
##Varying binary thresholds for continuous traits
#Louis Bell-Roberts
#17/01/2024

################################################################################################
####### IMPORTANT NOTE #########
#When choosing CS thresholds, should I be using all CS data that is available to calculate quantiles, or should I be using CS data for species present in that analysis?
#When I was preparing datasets to make the binary transformation thresholds for CS, I used threshold values based off all species for which we had CS data for (not including species with weird biology). Therefore, the same CS thresholds used are 40% = 188.600; 50 = 300; 60% = 501.800; 70% = 1000 (This explains why I had to go up to 70/40 threshold for CS_MF analysis because generally species in this analysis had quite large colony sizes, therefore too many species were being classed as having small colonies)
################################################################################################

#Caste, MF, CS and PG to be made discrete
##Using the 40% quantile and 60% quantile thresholds for the sensitivity analysis

#Filtering data - this will create separate data frames
##Filter for pairs of variables (CS and Caste, MF and Caste)

#Prune trees and prune databases

#Convert data frame to txt file

#Convert pruned phylogenies to Nexus format: if using .tre file then it is likely in Newick format
##BayesTraits requires Nexus format

library(tidyverse)
library(ape)

setwd("/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/BayesTraits/Data/")

#Read in the data
pre_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Following_review/ant_data.csv", header = T)

#Select columns of interest
ant_data <- pre_data %>% dplyr::select(species, colony.size, queen.mating.frequency, queen.mating.frequency.categorical, queen.number.continuous, queen.number.binary, queen.number.categorical, caste.number)

#Create binary caste.number columns
ant_data$CasteBin <- ifelse(ant_data$caste.number > 1, 1, 0)
ant_data$CasteBin <- as.factor(ant_data$CasteBin)

#Read in phylogenetic trees
ant_trees <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/Economo_2018_400.tre")


##################################
# CREATE DATA FRAMES FOR CS_Caste #
##################################

#Filter for species with colony size and caste number data
CS_Caste_data <- ant_data %>% dplyr::filter(complete.cases(colony.size), complete.cases(caste.number))

#Prune all 400 trees
CS_Caste_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, CS_Caste_data$species))

#Prune database to match species in the tree
CS_Caste_data_pruned <- dplyr::filter(CS_Caste_data, species %in% CS_Caste_trees_pruned[[1]]$tip.label)

#Calculate CS quantiles using all of the CS data that is available (not just the CS data available for that specific analysis)
quantile(ant_data$colony.size, probs = seq(0, 1, 0.1), na.rm = T) 

#Quantiles using all species which we have CS data for (excluding species with weird biology): 30% = 100; 40% = 172.00; median = 300; 60% = 509.00; 70% = 1000

#Transform data into binary variables
##Threshold = low (172)
CS_Caste_data_pruned_low_thresh <- CS_Caste_data_pruned
CS_Caste_data_pruned_low_thresh$CSBin <- cut(CS_Caste_data_pruned_low_thresh$colony.size,
                                     breaks=c(-Inf, 172, Inf),
                                     labels=c("0", "1")) #Make CS into a binary variable
CS_Caste_data_pruned_low_thresh$CSBin <- as.factor(CS_Caste_data_pruned_low_thresh$CSBin)
##Select columns of interest
CS_Caste_data_pruned_low_thresh <- CS_Caste_data_pruned_low_thresh %>% dplyr::select(species, CSBin, CasteBin)
table(CS_Caste_data_pruned_low_thresh$CSBin, CS_Caste_data_pruned_low_thresh$CasteBin)

##Threshold = median (300)
CS_Caste_data_pruned_median_thresh <- CS_Caste_data_pruned
CS_Caste_data_pruned_median_thresh$CSBin <- cut(CS_Caste_data_pruned_median_thresh$colony.size,
                                     breaks=c(-Inf, 300, Inf),
                                     labels=c("0", "1")) #Make CS into a binary variable
CS_Caste_data_pruned_median_thresh$CSBin <- as.factor(CS_Caste_data_pruned_median_thresh$CSBin)
##Select columns of interest
CS_Caste_data_pruned_median_thresh <- CS_Caste_data_pruned_median_thresh %>% dplyr::select(species, CSBin, CasteBin)
table(CS_Caste_data_pruned_median_thresh$CSBin, CS_Caste_data_pruned_median_thresh$CasteBin)

##Threshold = high (509)
CS_Caste_data_pruned_high_thresh <- CS_Caste_data_pruned
CS_Caste_data_pruned_high_thresh$CSBin <- cut(CS_Caste_data_pruned_high_thresh$colony.size,
                                     breaks=c(-Inf, 509, Inf),
                                     labels=c("0", "1")) #Make CS into a binary variable
CS_Caste_data_pruned_high_thresh$CSBin <- as.factor(CS_Caste_data_pruned_high_thresh$CSBin)
##Select columns of interest
CS_Caste_data_pruned_high_thresh <- CS_Caste_data_pruned_high_thresh %>% dplyr::select(species, CSBin, CasteBin)
table(CS_Caste_data_pruned_high_thresh$CSBin, CS_Caste_data_pruned_high_thresh$CasteBin)

##Threshold = very high (1000) ------ Do I need to include this analysis?


##################################
# CREATE DATA FRAMES FOR MF_Caste #
##################################

#####
#Policing threshold for queen mating frequency

#Filter for species with mating frequency and caste number data
MF_Caste_data_police_thresh <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(caste.number))

#Prune all 400 trees
MF_Caste_police_thresh_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, MF_Caste_data_police_thresh$species))

#Prune database to match species in the tree
MF_Caste_data_police_thresh_pruned <- dplyr::filter(MF_Caste_data_police_thresh, species %in% MF_Caste_police_thresh_trees_pruned[[1]]$tip.label)

#Transform data into binary variables
##Threshold = 2 (predicted threshold at which policing occurs)
MF_Caste_data_police_thresh_pruned$MFBin <- cut(MF_Caste_data_police_thresh_pruned$queen.mating.frequency,
                                             breaks=c(-Inf, 2, Inf),
                                             labels=c("0", "1")) #Make MF into a binary variable
MF_Caste_data_police_thresh_pruned$MFBin <- as.factor(MF_Caste_data_police_thresh_pruned$MFBin)
##Select columns of interest
MF_Caste_data_police_thresh_pruned <- MF_Caste_data_police_thresh_pruned %>% dplyr::select(species, MFBin, CasteBin)
table(MF_Caste_data_police_thresh_pruned$MFBin, MF_Caste_data_police_thresh_pruned$CasteBin)


#####
#Binary thresholds based on queen.mating.frequency.categorical data

#Filter for species with queen.mating.frequency.categorical and caste number data
MFcat_Caste_data <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency.categorical), complete.cases(caste.number))

##Monandry vs facultative polyandry
#Exclude species that are obligately polyandrous
MFcat_fac_Caste_data <- MFcat_Caste_data %>% dplyr::filter(queen.mating.frequency.categorical < 2)
#Prune all 400 trees
MFcat_fac_Caste_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, MFcat_fac_Caste_data$species))
#Prune database to match species in the tree
MFcat_fac_Caste_data_pruned <- dplyr::filter(MFcat_fac_Caste_data, species %in% MFcat_fac_Caste_trees_pruned[[1]]$tip.label)
##Select columns of interest
MFcat_fac_Caste_data_pruned <- MFcat_fac_Caste_data_pruned %>% dplyr::select(species, queen.mating.frequency.categorical, CasteBin)


##Monandry vs obligate polyandry
#Exclude species that are facultatively polyandrous
MFcat_ob_Caste_data <- MFcat_Caste_data %>% dplyr::filter(queen.mating.frequency.categorical != 1)
#Prune all 400 trees
MFcat_ob_Caste_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, MFcat_ob_Caste_data$species))
#Prune database to match species in the tree
MFcat_ob_Caste_data_pruned <- dplyr::filter(MFcat_ob_Caste_data, species %in% MFcat_ob_Caste_trees_pruned[[1]]$tip.label)
#Assign 1 to species that are obligately polyandrous
MFcat_ob_Caste_data_pruned$queen.mating.frequency.categorical[MFcat_ob_Caste_data_pruned$queen.mating.frequency.categorical == 2] <- 1
##Select columns of interest
MFcat_ob_Caste_data_pruned <- MFcat_ob_Caste_data_pruned %>% dplyr::select(species, queen.mating.frequency.categorical, CasteBin)
table(MFcat_ob_Caste_data_pruned$queen.mating.frequency.categorical, MFcat_ob_Caste_data_pruned$CasteBin)

##Monandry+facultative polyandry vs obligate polyandry
MFcat_fac_ob_Caste_data <- MFcat_Caste_data
#Assign 0 to facultatively polyandrous species
MFcat_fac_ob_Caste_data$queen.mating.frequency.categorical[MFcat_fac_ob_Caste_data$queen.mating.frequency.categorical == 1] <- 0 
#Assign 1 to obligately polyandrous species
MFcat_fac_ob_Caste_data$queen.mating.frequency.categorical[MFcat_fac_ob_Caste_data$queen.mating.frequency.categorical == 2] <- 1
#Prune all 400 trees
MFcat_fac_ob_Caste_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, MFcat_fac_ob_Caste_data$species))
#Prune database to match species in the tree
MFcat_fac_ob_Caste_data_pruned <- dplyr::filter(MFcat_fac_ob_Caste_data, species %in% MFcat_fac_ob_Caste_trees_pruned[[1]]$tip.label)
##Select columns of interest
MFcat_fac_ob_Caste_data_pruned <- MFcat_fac_ob_Caste_data_pruned %>% dplyr::select(species, queen.mating.frequency.categorical, CasteBin)
table(MFcat_fac_ob_Caste_data_pruned$queen.mating.frequency.categorical, MFcat_fac_ob_Caste_data_pruned$CasteBin)

##################################
# CREATE DATA FRAMES FOR MF_PG #
##################################

#Binary thresholds based on categorical mating frequency and categorical queen number data

#Filter for species with queen.mating.frequency.categorical and caste number data
MFcat_PGcat_data <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency.categorical), complete.cases(queen.number.categorical))

##Monandry vs facultative polyandry & Monogyny vs facultative polygyny
#Exclude species with obligate polyandry and obligate polygyny
MFcat_fac_PGcat_fac_data <- MFcat_PGcat_data %>% dplyr::filter(queen.mating.frequency.categorical < 2, queen.number.categorical < 2)
#Prune all 400 trees
MFcat_fac_PGcat_fac_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, MFcat_fac_PGcat_fac_data$species))
#Prune database to match species in the tree
MFcat_fac_PGcat_fac_data_pruned <- dplyr::filter(MFcat_fac_PGcat_fac_data, species %in% MFcat_fac_PGcat_fac_trees_pruned[[1]]$tip.label)
##Select columns of interest
MFcat_fac_PGcat_fac_data_pruned <- MFcat_fac_PGcat_fac_data_pruned %>% dplyr::select(species, queen.mating.frequency.categorical, queen.number.categorical)

##Monandry vs obligate polyandry & Monogyny vs facultative polygyny
#Exclude species with facultative polyandry and obligate polygyny
MFcat_ob_PGcat_fac_data <- MFcat_PGcat_data %>% dplyr::filter(queen.mating.frequency.categorical != 1, queen.number.categorical < 2)
#Assign 1 to species with obligate polyandry
MFcat_ob_PGcat_fac_data$queen.mating.frequency.categorical[MFcat_ob_PGcat_fac_data$queen.mating.frequency.categorical == 2] <- 1
#Prune all 400 trees
MFcat_ob_PGcat_fac_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, MFcat_ob_PGcat_fac_data$species))
#Prune database to match species in the tree
MFcat_ob_PGcat_fac_data_pruned <- dplyr::filter(MFcat_ob_PGcat_fac_data, species %in% MFcat_ob_PGcat_fac_trees_pruned[[1]]$tip.label)
##Select columns of interest
MFcat_ob_PGcat_fac_data_pruned <- MFcat_ob_PGcat_fac_data_pruned %>% dplyr::select(species, queen.mating.frequency.categorical, queen.number.categorical)

##Monandry+facultative polyandry vs obligate polyandry & Monogyny vs facultative polygyny
#Exclude species with obligate polygyny
MFcat_fac_ob_PGcat_fac_data <- MFcat_PGcat_data %>% dplyr::filter(queen.number.categorical < 2)
#Assign 0 to species with facultative polyandry
MFcat_fac_ob_PGcat_fac_data$queen.mating.frequency.categorical[MFcat_fac_ob_PGcat_fac_data$queen.mating.frequency.categorical == 1] <- 0
#Assign 1 to species with obligate polyandry
MFcat_fac_ob_PGcat_fac_data$queen.mating.frequency.categorical[MFcat_fac_ob_PGcat_fac_data$queen.mating.frequency.categorical == 2] <- 1
#Prune all 400 trees
MFcat_fac_ob_PGcat_fac_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, MFcat_fac_ob_PGcat_fac_data$species))
#Prune database to match species in the tree
MFcat_fac_ob_PGcat_fac_data_pruned <- dplyr::filter(MFcat_fac_ob_PGcat_fac_data, species %in% MFcat_fac_ob_PGcat_fac_trees_pruned[[1]]$tip.label)
##Select columns of interest
MFcat_fac_ob_PGcat_fac_data_pruned <- MFcat_fac_ob_PGcat_fac_data_pruned %>% dplyr::select(species, queen.mating.frequency.categorical, queen.number.categorical)

##Monandry vs facultative polyandry & Monogyny vs obligate polygyny
#Exclude species with obligate polyandry and facultative polygyny
MFcat_fac_PGcat_ob_data <- MFcat_PGcat_data %>% dplyr::filter(queen.mating.frequency.categorical < 2, queen.number.categorical != 1)
#Assign 1 to species with obligate polygyny
MFcat_fac_PGcat_ob_data$queen.number.categorical[MFcat_fac_PGcat_ob_data$queen.number.categorical == 2] <- 1
#Prune all 400 trees
MFcat_fac_PGcat_ob_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, MFcat_fac_PGcat_ob_data$species))
#Prune database to match species in the tree
MFcat_fac_PGcat_ob_data_pruned <- dplyr::filter(MFcat_fac_PGcat_ob_data, species %in% MFcat_fac_PGcat_ob_trees_pruned[[1]]$tip.label)
##Select columns of interest
MFcat_fac_PGcat_ob_data_pruned <- MFcat_fac_PGcat_ob_data_pruned %>% dplyr::select(species, queen.mating.frequency.categorical, queen.number.categorical)

##Monandry vs obligate polyandry & Monogyny vs obligate polygyny
#Exclude species with facultative polyandry and facultative polygyny
MFcat_ob_PGcat_ob_data <- MFcat_PGcat_data %>% dplyr::filter(queen.mating.frequency.categorical != 1, queen.number.categorical != 1)
#Assign 1 to species with obligate polyandry
MFcat_ob_PGcat_ob_data$queen.mating.frequency.categorical[MFcat_ob_PGcat_ob_data$queen.mating.frequency.categorical == 2] <- 1
#Assign 1 to species with obligate polygyny
MFcat_ob_PGcat_ob_data$queen.number.categorical[MFcat_ob_PGcat_ob_data$queen.number.categorical == 2] <- 1
#Prune all 400 trees
MFcat_ob_PGcat_ob_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, MFcat_ob_PGcat_ob_data$species))
#Prune database to match species in the tree
MFcat_ob_PGcat_ob_data_pruned <- dplyr::filter(MFcat_ob_PGcat_ob_data, species %in% MFcat_ob_PGcat_ob_trees_pruned[[1]]$tip.label)
##Select columns of interest
MFcat_ob_PGcat_ob_data_pruned <- MFcat_ob_PGcat_ob_data_pruned %>% dplyr::select(species, queen.mating.frequency.categorical, queen.number.categorical)

##Monandry+facultative polyandry vs obligate polyandry & Monogyny vs obligate polygyny
#Exclude species with facultative polygyny
MFcat_fac_ob_PGcat_ob_data <- MFcat_PGcat_data %>% dplyr::filter(queen.number.categorical != 1)
#Assign 0 to species with facultative polyandry
MFcat_fac_ob_PGcat_ob_data$queen.mating.frequency.categorical[MFcat_fac_ob_PGcat_ob_data$queen.mating.frequency.categorical == 1] <- 0
#Assign 1 to species with obligate polyandry
MFcat_fac_ob_PGcat_ob_data$queen.mating.frequency.categorical[MFcat_fac_ob_PGcat_ob_data$queen.mating.frequency.categorical == 2] <- 1
#Assign 1 to species with obligate polygyny
MFcat_fac_ob_PGcat_ob_data$queen.number.categorical[MFcat_fac_ob_PGcat_ob_data$queen.number.categorical == 2] <- 1
#Prune all 400 trees
MFcat_fac_ob_PGcat_ob_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, MFcat_fac_ob_PGcat_ob_data$species))
#Prune database to match species in the tree
MFcat_fac_ob_PGcat_ob_data_pruned <- dplyr::filter(MFcat_fac_ob_PGcat_ob_data, species %in% MFcat_fac_ob_PGcat_ob_trees_pruned[[1]]$tip.label)
##Select columns of interest
MFcat_fac_ob_PGcat_ob_data_pruned <- MFcat_fac_ob_PGcat_ob_data_pruned %>% dplyr::select(species, queen.mating.frequency.categorical, queen.number.categorical)
table(MFcat_fac_ob_PGcat_ob_data_pruned$queen.mating.frequency.categorical, MFcat_fac_ob_PGcat_ob_data_pruned$queen.number.categorical)

##Monandry vs facultative polyandry & Monogyny+facultative polygyny vs obligate polygyny
#Exclude species with obligate polyandry
MFcat_fac_PGcat_fac_ob_data <- MFcat_PGcat_data %>% dplyr::filter(queen.mating.frequency.categorical < 2)
#Assign 0 to species with facultative polygyny
MFcat_fac_PGcat_fac_ob_data$queen.number.categorical[MFcat_fac_PGcat_fac_ob_data$queen.number.categorical == 1] <- 0
#Assign 1 to species with obligate polygyny
MFcat_fac_PGcat_fac_ob_data$queen.number.categorical[MFcat_fac_PGcat_fac_ob_data$queen.number.categorical == 2] <- 1
#Prune all 400 trees
MFcat_fac_PGcat_fac_ob_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, MFcat_fac_PGcat_fac_ob_data$species))
#Prune database to match species in the tree
MFcat_fac_PGcat_fac_ob_data_pruned <- dplyr::filter(MFcat_fac_PGcat_fac_ob_data, species %in% MFcat_fac_PGcat_fac_ob_trees_pruned[[1]]$tip.label)
##Select columns of interest
MFcat_fac_PGcat_fac_ob_data_pruned <- MFcat_fac_PGcat_fac_ob_data_pruned %>% dplyr::select(species, queen.mating.frequency.categorical, queen.number.categorical)

##Monandry vs obligate polyandry & Monogyny+facultative polygyny vs obligate polygyny
#Exclude species with facultative polyandry
MFcat_ob_PGcat_fac_ob_data <- MFcat_PGcat_data %>% dplyr::filter(queen.mating.frequency.categorical != 1)
#Assign 1 to species with obligate polyandry
MFcat_ob_PGcat_fac_ob_data$queen.mating.frequency.categorical[MFcat_ob_PGcat_fac_ob_data$queen.mating.frequency.categorical == 2] <- 1
#Assign 0 to species with facultative polygyny
MFcat_ob_PGcat_fac_ob_data$queen.number.categorical[MFcat_ob_PGcat_fac_ob_data$queen.number.categorical == 1] <- 0
#Assign 1 to species with obligate polygyny
MFcat_ob_PGcat_fac_ob_data$queen.number.categorical[MFcat_ob_PGcat_fac_ob_data$queen.number.categorical == 2] <- 1
#Prune all 400 trees
MFcat_ob_PGcat_fac_ob_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, MFcat_ob_PGcat_fac_ob_data$species))
#Prune database to match species in the tree
MFcat_ob_PGcat_fac_ob_data_pruned <- dplyr::filter(MFcat_ob_PGcat_fac_ob_data, species %in% MFcat_ob_PGcat_fac_ob_trees_pruned[[1]]$tip.label)
##Select columns of interest
MFcat_ob_PGcat_fac_ob_data_pruned <- MFcat_ob_PGcat_fac_ob_data_pruned %>% dplyr::select(species, queen.mating.frequency.categorical, queen.number.categorical)

##Monandry+facultative polyandry vs obligate polyandry & Monogyny+facultative polygyny vs obligate polygyny
#No further filtering required
MFcat_fac_ob_PGcat_fac_ob_data <- MFcat_PGcat_data
#Assign 0 to species with facultative polyandry
MFcat_fac_ob_PGcat_fac_ob_data$queen.mating.frequency.categorical[MFcat_fac_ob_PGcat_fac_ob_data$queen.mating.frequency.categorical == 1] <- 0
#Assign 1 to species with obligate polyandry
MFcat_fac_ob_PGcat_fac_ob_data$queen.mating.frequency.categorical[MFcat_fac_ob_PGcat_fac_ob_data$queen.mating.frequency.categorical == 2] <- 1
#Assign 0 to species with facultative polygyny
MFcat_fac_ob_PGcat_fac_ob_data$queen.number.categorical[MFcat_fac_ob_PGcat_fac_ob_data$queen.number.categorical == 1] <- 0
#Assign 1 to species with obligate polygyny
MFcat_fac_ob_PGcat_fac_ob_data$queen.number.categorical[MFcat_fac_ob_PGcat_fac_ob_data$queen.number.categorical == 2] <- 1
#Prune all 400 trees
MFcat_fac_ob_PGcat_fac_ob_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, MFcat_fac_ob_PGcat_fac_ob_data$species))
#Prune database to match species in the tree
MFcat_fac_ob_PGcat_fac_ob_data_pruned <- dplyr::filter(MFcat_fac_ob_PGcat_fac_ob_data, species %in% MFcat_fac_ob_PGcat_fac_ob_trees_pruned[[1]]$tip.label)
##Select columns of interest
MFcat_fac_ob_PGcat_fac_ob_data_pruned <- MFcat_fac_ob_PGcat_fac_ob_data_pruned %>% dplyr::select(species, queen.mating.frequency.categorical, queen.number.categorical)
table(MFcat_fac_ob_PGcat_fac_ob_data_pruned$queen.mating.frequency.categorical, MFcat_fac_ob_PGcat_fac_ob_data_pruned$queen.number.categorical)


########
#Write dataframes to file as .txt files that are space delimited
##################################
# CREATE DATA FRAMES FOR CS_Caste #
##################################
write.table(CS_Caste_data_pruned_low_thresh, file = "CS_Caste/ant_data_BayesTraits_CS_low_thresh_Caste.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(CS_Caste_data_pruned_median_thresh, file = "CS_Caste/ant_data_BayesTraits_CS_median_thresh_Caste.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(CS_Caste_data_pruned_high_thresh, file = "CS_Caste/ant_data_BayesTraits_CS_high_thresh_Caste.txt", sep = "\t", row.names = F, col.names = F, quote = F)

##################################
# CREATE DATA FRAMES FOR MF_Caste #
##################################
write.table(MF_Caste_data_police_thresh_pruned, file = "MF_Caste/ant_data_BayesTraits_MF_police_thresh_Caste.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(MFcat_fac_Caste_data_pruned, file = "MF_Caste/ant_data_BayesTraits_MFcat_fac_Caste.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(MFcat_ob_Caste_data_pruned, file = "MF_Caste/ant_data_BayesTraits_MFcat_ob_Caste.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(MFcat_fac_ob_Caste_data_pruned, file = "MF_Caste/ant_data_BayesTraits_MFcat_fac_ob_Caste.txt", sep = "\t", row.names = F, col.names = F, quote = F)

##################################
# CREATE DATA FRAMES FOR MF_PG #
##################################
write.table(MFcat_fac_PGcat_fac_data_pruned, file = "MF_PG/ant_data_BayesTraits_MFcat_fac_PGcat_fac.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(MFcat_ob_PGcat_fac_data_pruned, file = "MF_PG/ant_data_BayesTraits_MFcat_ob_PGcat_fac.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(MFcat_fac_ob_PGcat_fac_data_pruned, file = "MF_PG/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac.txt", sep = "\t", row.names = F, col.names = F, quote = F)

write.table(MFcat_fac_PGcat_ob_data_pruned, file = "MF_PG/ant_data_BayesTraits_MFcat_fac_PGcat_ob.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(MFcat_ob_PGcat_ob_data_pruned, file = "MF_PG/ant_data_BayesTraits_MFcat_ob_PGcat_ob.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(MFcat_fac_ob_PGcat_ob_data_pruned, file = "MF_PG/ant_data_BayesTraits_MFcat_fac_ob_PGcat_ob.txt", sep = "\t", row.names = F, col.names = F, quote = F)

write.table(MFcat_fac_PGcat_fac_ob_data_pruned, file = "MF_PG/ant_data_BayesTraits_MFcat_fac_PGcat_fac_ob.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(MFcat_ob_PGcat_fac_ob_data_pruned, file = "MF_PG/ant_data_BayesTraits_MFcat_ob_PGcat_fac_ob.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(MFcat_fac_ob_PGcat_fac_ob_data_pruned, file = "MF_PG/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_ob.txt", sep = "\t", row.names = F, col.names = F, quote = F)


########
#Write 400 pruned trees to nexus file (necessary for BayesTraits) for each dataset
##Save as a single multiphylo object

##################################
# CREATE DATA FRAMES FOR CS_Caste #
##################################
setwd("/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/BayesTraits/Trees/")

write.nexus(CS_Caste_trees_pruned, file = "CS_Caste/ant_trees_BayesTraits_CS_Caste_pruned.nex")

##################################
# CREATE DATA FRAMES FOR MF_Caste #
##################################
write.nexus(MF_Caste_police_thresh_trees_pruned, file = "MF_Caste/ant_trees_BayesTraits_MF_police_thresh_Caste_pruned.nex")
write.nexus(MFcat_fac_Caste_trees_pruned, file = "MF_Caste/ant_trees_BayesTraits_MFcat_fac_Caste_pruned.nex")
write.nexus(MFcat_ob_Caste_trees_pruned, file = "MF_Caste/ant_trees_BayesTraits_MFcat_ob_Caste_pruned.nex")
write.nexus(MFcat_fac_ob_Caste_trees_pruned, file = "MF_Caste/ant_trees_BayesTraits_MFcat_fac_ob_Caste_pruned.nex")

##################################
# CREATE DATA FRAMES FOR MF_PG #
##################################
write.nexus(MFcat_fac_PGcat_fac_trees_pruned, file = "MF_PG/ant_trees_BayesTraits_MFcat_fac_PGcat_fac_pruned.nex")
write.nexus(MFcat_ob_PGcat_fac_trees_pruned, file = "MF_PG/ant_trees_BayesTraits_MFcat_ob_PGcat_fac_pruned.nex")
write.nexus(MFcat_fac_ob_PGcat_fac_trees_pruned, file = "MF_PG/ant_trees_BayesTraits_MFcat_fac_ob_PGcat_fac_pruned.nex")

write.nexus(MFcat_fac_PGcat_ob_trees_pruned, file = "MF_PG/ant_trees_BayesTraits_MFcat_fac_PGcat_ob_pruned.nex")
write.nexus(MFcat_ob_PGcat_ob_trees_pruned, file = "MF_PG/ant_trees_BayesTraits_MFcat_ob_PGcat_ob_pruned.nex")
write.nexus(MFcat_fac_ob_PGcat_ob_trees_pruned, file = "MF_PG/ant_trees_BayesTraits_MFcat_fac_ob_PGcat_ob_pruned.nex")

write.nexus(MFcat_fac_PGcat_fac_ob_trees_pruned, file = "MF_PG/ant_trees_BayesTraits_MFcat_fac_PGcat_fac_ob_pruned.nex")
write.nexus(MFcat_ob_PGcat_fac_ob_trees_pruned, file = "MF_PG/ant_trees_BayesTraits_MFcat_ob_PGcat_fac_ob_pruned.nex")
write.nexus(MFcat_fac_ob_PGcat_fac_ob_trees_pruned, file = "MF_PG/ant_trees_BayesTraits_MFcat_fac_ob_PGcat_fac_ob_pruned.nex")





