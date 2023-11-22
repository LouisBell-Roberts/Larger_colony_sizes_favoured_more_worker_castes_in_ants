########################ASR with corHMM - analysis across 400 trees######################
###Runs ASR for caste number as a binary variable (single/multiple castes)
##Analysing all species for which data is available on the number of worker castes
#Uses AICc values to identify the best ASR model
#Louis Bell-Roberts
#22/11/2023

#Load packages
library(dplyr)
library(ape)
library(phytools)
library(ggplot2)
library(corHMM)
library(readxl)

#Load custom functions for running corHMM models
corHMM_multiphylo_1rat_ER <- function(multiphylo, trait_data) {
  
  # Run corHMM ancestral state reconstruction for each phylogeny
  for (i in 1:length(multiphylo)) {
    phylo <- multiphylo[[i]]
    output <- corHMM(phy = phylo, data = trait_data, rate.cat = 1, model = "ER")
    
    # Save output to a file
    output_file <- paste0("corHMM_ASR_Caste_only_1rat_ER_", i, ".rda")
    saveRDS(output, file = output_file)
    cat("Saved output to", output_file, "\n")
    
  }
}

corHMM_multiphylo_1rat_ARD <- function(multiphylo, trait_data) {
  
  # Run corHMM ancestral state reconstruction for each phylogeny
  for (i in 1:length(multiphylo)) {
    phylo <- multiphylo[[i]]
    output <- corHMM(phy = phylo, data = trait_data, rate.cat = 1, model = "ARD")
    
    # Save output to a file
    output_file <- paste0("corHMM_ASR_Caste_only_1rat_ARD_", i, ".rda")
    saveRDS(output, file = output_file)
    cat("Saved output to", output_file, "\n")
    
  }
}

corHMM_multiphylo_2rat_ER <- function(multiphylo, trait_data) {
  
  # Run corHMM ancestral state reconstruction for each phylogeny
  for (i in 1:length(multiphylo)) {
    phylo <- multiphylo[[i]]
    output <- corHMM(phy = phylo, data = trait_data, rate.cat = 2, model = "ER")
    
    # Save output to a file
    output_file <- paste0("corHMM_ASR_Caste_only_2rat_ER_", i, ".rda")
    saveRDS(output, file = output_file)
    cat("Saved output to", output_file, "\n")
    
  }
}

corHMM_multiphylo_2rat_ARD <- function(multiphylo, trait_data) {
  
  # Run corHMM ancestral state reconstruction for each phylogeny
  for (i in 1:length(multiphylo)) {
    phylo <- multiphylo[[i]]
    output <- corHMM(phy = phylo, data = trait_data, rate.cat = 2, model = "ARD")
    
    # Save output to a file
    output_file <- paste0("corHMM_ASR_Caste_only_2rat_ARD_", i, ".rda")
    saveRDS(output, file = output_file)
    cat("Saved output to", output_file, "\n")
    
  }
}

corHMM_multiphylo_3rat_ER <- function(multiphylo, trait_data) {
  
  # Run corHMM ancestral state reconstruction for each phylogeny
  for (i in 1:length(multiphylo)) {
    phylo <- multiphylo[[i]]
    output <- corHMM(phy = phylo, data = trait_data, rate.cat = 3, model = "ER")
    
    # Save output to a file
    output_file <- paste0("corHMM_ASR_Caste_only_3rat_ER_", i, ".rda")
    saveRDS(output, file = output_file)
    cat("Saved output to", output_file, "\n")
    
  }
}

corHMM_multiphylo_3rat_ARD <- function(multiphylo, trait_data) {
  
  # Run corHMM ancestral state reconstruction for each phylogeny
  for (i in 1:length(multiphylo)) {
    phylo <- multiphylo[[i]]
    output <- corHMM(phy = phylo, data = trait_data, rate.cat = 3, model = "ARD")
    
    # Save output to a file
    output_file <- paste0("corHMM_ASR_Caste_only_3rat_ARD_", i, ".rda")
    saveRDS(output, file = output_file)
    cat("Saved output to", output_file, "\n")
    
  }
}

#Function that reads in 400 corHMM models, in order of their tree number, and then creates a vector of their AICc values
corHMM_AIC <- function(folder_path) {
  
  # get the file names
  file_names <- list.files(folder_path)
  
  # extract the numbers from the file names using regular expressions
  file_numbers <- as.numeric(gsub("[^0-9]", "", file_names))
  
  # sort the file names based on the corresponding numbers
  sorted_file_names <- file_names[order(file_numbers)]
  
  # read in the files and save them to a list
  corHMM_results <- lapply(sorted_file_names, function(x) readRDS(file.path(folder_path, x)))
  
  #Extracting all of the AICc values from the list of corHMM outputs
  # create an empty vector to hold the AIC values
  AIC_values <- numeric()
  
  # loop through each element of the list
  for (i in seq_along(corHMM_results)) {
    # extract the AIC value from the ith element and append it to the AIC_values vector
    AIC_values <- c(AIC_values, corHMM_results[[i]]$AICc)
  }
  return(AIC_values)
}

#########

#Read in data file
ant_data <- as.data.frame(read_xlsx("ant_data.xlsx", col_names = T))

#Set variables so that they're in the correct structure
ant_data[ant_data == ""] <- NA #Replace blank by NA
class(ant_data$caste.number) #numeric

#Make caste.number a binary variable
ant_data$CasteBin <- ifelse(ant_data$caste.number > 1,"polymorphic","monomorphic")
ant_data$CasteBin <- as.factor(ant_data$CasteBin)

#Rename 'species' column as 'animal'
ant_data <- ant_data %>% dplyr::rename(animal = species)

#Read in sample of 400 phylogenetic trees
ant.trees <- read.tree(file = "Economo_2018_400.tre")

#Filter data
ASRdata <- dplyr::filter(ant_data, caste.number >=1)
ASRdata <- dplyr::select(ASRdata, animal, CasteBin)

#Prune trees
ant_trees_pruned <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, ASRdata$animal))

#Prune database
ASRdata <- dplyr::filter(ASRdata, animal %in% ant_trees_pruned[[1]]$tip.label)

# check that species names in the tree and the data match
ASRdata$animal[which((ASRdata$animal %in% ant_trees_pruned[[1]]$tip.label) == FALSE)]	# 0 means all match
ant_trees_pruned[[1]]$tip.label[which((ant_trees_pruned[[1]]$tip.label %in% ASRdata$animal) == FALSE)]	# 0 means all match


### ESTIMATING ANCESTRAL STATES ###
## Models to estimate the number of worker castes at each node on the phylogeny ##

#Assign first column as species names
ASRdata <- data.frame(ASRdata$animal, ASRdata$CasteBin)

setwd("") #Assign file path for models to be saved
corHMM_multiphylo_1rat_ER(ant_trees_pruned, ASRdata)
setwd("") #Assign file path for models to be saved
corHMM_multiphylo_1rat_ARD(ant_trees_pruned, ASRdata)

setwd("") #Assign file path for models to be saved
corHMM_multiphylo_2rat_ER(ant_trees_pruned, ASRdata)
setwd("") #Assign file path for models to be saved
corHMM_multiphylo_2rat_ARD(ant_trees_pruned, ASRdata)

setwd("") #Assign file path for models to be saved
corHMM_multiphylo_3rat_ER(ant_trees_pruned, ASRdata)
setwd("") #Assign file path for models to be saved
corHMM_multiphylo_3rat_ARD(ant_trees_pruned, ASRdata)




###############
####Read in each of the different corHMM ASR models
###############
folder_path_1rat_ER <- "" #File path to location where single rate, equal rates models are saved
folder_path_1rat_ARD <- "" #File path to location where single rate, all rates different models are saved
folder_path_2rat_ER <- "" #File path to location where two rate, equal rates models are saved
folder_path_2rat_ARD <- "" #File path to location where two rate, all rates different models are saved
folder_path_3rat_ER <- "" #File path to location where three rate, equal rates models are saved
folder_path_3rat_ARD <- "" #File path to location where three rate, all rates different models are saved

#Obtain AICc scores for each type of model over the 400 trees
CS_Caste_1rat_ER_AICvec <- corHMM_AIC(folder_path_1rat_ER)
CS_Caste_1rat_ARD_AICvec <- corHMM_AIC(folder_path_1rat_ARD)
CS_Caste_2rat_ER_AICvec <- corHMM_AIC(folder_path_2rat_ER)
CS_Caste_2rat_ARD_AICvec <- corHMM_AIC(folder_path_2rat_ARD)
CS_Caste_3rat_ER_AICvec <- corHMM_AIC(folder_path_3rat_ER)
CS_Caste_3rat_ARD_AICvec <- corHMM_AIC(folder_path_3rat_ARD)

# Create data frame
AIC_df <- data.frame(AIC_1rat_ER = CS_Caste_1rat_ER_AICvec,
                     AIC_1rat_ARD = CS_Caste_1rat_ARD_AICvec,
                     AIC_2rat_ER = CS_Caste_2rat_ER_AICvec,
                     AIC_2rat_ARD = CS_Caste_2rat_ARD_AICvec,
                     AIC_3rat_ER = CS_Caste_3rat_ER_AICvec,
                     AIC_3rat_ARD = CS_Caste_3rat_ARD_AICvec)
apply(AIC_df, 2, mean)

#Identify the best model for each individual tree
# Identify column with smallest value for each row
smallest_column <- apply(AIC_df, 1, function(row) {
  col_index <- which.min(row)
  colnames(AIC_df)[col_index]
})
table(smallest_column)




########END OF SCRIPT##########
