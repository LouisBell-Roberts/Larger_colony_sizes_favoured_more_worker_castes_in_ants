########################ASR with corHMM estimating ancestral numbers of worker castes (1 to 4 castes)######################
##Runs ASR with data for just Caste number
#Louis Bell-Roberts
#01/02/2024


#Load packages
library(dplyr)
library(ape)
library(phytools)
library(ggplot2)
library(corHMM)
library(doParallel)
library(gtools)

#Run in parallel using 50 cores
registerDoParallel(50)

#Load custom functions for running corHMM models in parallel with adjustable model settings
corHMM_multiphylo <- function(multiphylo, trait_data, model_type, rate_cat) {
  
  # Run corHMM ancestral state reconstruction for each phylogeny
  foreach(i = 1:length(multiphylo)) %dopar% {
    output <- corHMM(phy = multiphylo[[i]], data = trait_data, rate.cat = rate_cat, model = model_type)
    
    # Define dynamic output file name
    output_file <- paste0("corHMM_ASR_Caste_only_discrete_", rate_cat, "rat_", model_type, "_", i, ".rda")
    
    # Save output to a file
    saveRDS(output, file = output_file)
    cat("Saved output to", output_file, "\n")
  }
}


#Function that reads in 400 corHMM models, in order of their tree number, and then creates a vector of their AICc values
corHMM_AIC <- function(folder_path) {
  
  # get the file names
  file_names <- list.files(folder_path)
  
  # sort the file names based on the corresponding numbers
  sorted_file_names <- mixedsort(file_names)
  
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
ant_data <- read.csv("ant_data.csv")

#Assign caste.number as numeric
ant_data$caste.number <- as.numeric(ant_data$caste.number)

#Rename 'species' column as 'animal'
ant_data <- ant_data %>% dplyr::rename(animal = species)

#Read in sample of 400 phylogenetic trees
ant.trees <- read.tree(file ="Economo_2018_400.tre")

#Filter data
ASRdata <- dplyr::filter(ant_data, caste.number >=1)
ASRdata <- dplyr::select(ASRdata, animal, caste.number)

#Prune trees
ant_trees_pruned <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, ASRdata$animal))

#Prune database
ASRdata <- dplyr::filter(ASRdata, animal %in% ant_trees_pruned[[1]]$tip.label)

#Check that species names in the tree and the data match
ASRdata$animal[which((ASRdata$animal %in% ant_trees_pruned[[1]]$tip.label) == FALSE)]	# 0 means all match
ant_trees_pruned[[1]]$tip.label[which((ant_trees_pruned[[1]]$tip.label %in% ASRdata$animal) == FALSE)]	# 0 means all match


### ESTIMATING ANCESTRAL STATES ###
## Models to estimate the number of worker castes at each node on the phylogeny ##

####### 1 rate
setwd("") #Assign file path for models to be saved
corHMM_multiphylo(multiphylo = ant_trees_pruned, trait_data = ASRdata, model_type = "ER", rate_cat = 1)
setwd("") #Assign file path for models to be saved
corHMM_multiphylo(multiphylo = ant_trees_pruned, trait_data = ASRdata, model_type = "SYM", rate_cat = 1)
setwd("") #Assign file path for models to be saved
corHMM_multiphylo(multiphylo = ant_trees_pruned, trait_data = ASRdata, model_type = "ARD", rate_cat = 1)


####### 2 rates
setwd("") #Assign file path for models to be saved
corHMM_multiphylo(multiphylo = ant_trees_pruned, trait_data = ASRdata, model_type = "ER", rate_cat = 2)
setwd("") #Assign file path for models to be saved
corHMM_multiphylo(multiphylo = ant_trees_pruned, trait_data = ASRdata, model_type = "SYM", rate_cat = 2)
setwd("") #Assign file path for models to be saved
corHMM_multiphylo(multiphylo = ant_trees_pruned, trait_data = ASRdata, model_type = "ARD", rate_cat = 2)

####### 3 rates
setwd("") #Assign file path for models to be saved
corHMM_multiphylo(multiphylo = ant_trees_pruned, trait_data = ASRdata, model_type = "ER", rate_cat = 3)
setwd("") #Assign file path for models to be saved
corHMM_multiphylo(multiphylo = ant_trees_pruned, trait_data = ASRdata, model_type = "SYM", rate_cat = 3)
setwd("") #Assign file path for models to be saved
corHMM_multiphylo(multiphylo = ant_trees_pruned, trait_data = ASRdata, model_type = "ARD", rate_cat = 3)


###############
####Read in each of the different corHMM ASR models
###############

folder_path_1rat_ER <- "" #File path to location where single rate, equal rates models are saved
folder_path_1rat_SYM <- "" #File path to location where single rate, symmetrical rates models are saved
folder_path_1rat_ARD <- "" #File path to location where single rate, all rates different models are saved
folder_path_2rat_ER <- "" #File path to location where two rate, equal rates models are saved
folder_path_2rat_SYM <- "" #File path to location where two rate, symmetrical rates models are saved
folder_path_2rat_ARD <- "" #File path to location where two rate, all rates different models are saved
folder_path_3rat_ER <- "" #File path to location where three rate, equal rates models are saved
folder_path_3rat_SYM <- "" #File path to location where three rate, symmetrical rates models are saved
folder_path_3rat_ARD <- "" #File path to location where three rate, all rates different models are saved

#Obtain AICc scores for each type of model over the 400 trees
CS_Caste_1rat_ER_AICvec <- corHMM_AIC(folder_path_1rat_ER)
CS_Caste_1rat_SYM_AICvec <- corHMM_AIC(folder_path_1rat_SYM)
CS_Caste_1rat_ARD_AICvec <- corHMM_AIC(folder_path_1rat_ARD)
CS_Caste_2rat_ER_AICvec <- corHMM_AIC(folder_path_2rat_ER)
CS_Caste_2rat_SYM_AICvec <- corHMM_AIC(folder_path_2rat_SYM)
CS_Caste_2rat_ARD_AICvec <- corHMM_AIC(folder_path_2rat_ARD)
CS_Caste_3rat_ER_AICvec <- corHMM_AIC(folder_path_3rat_ER)
CS_Caste_3rat_SYM_AICvec <- corHMM_AIC(folder_path_3rat_SYM)
CS_Caste_3rat_ARD_AICvec <- corHMM_AIC(folder_path_3rat_ARD)

# Create data frame
AIC_df <- data.frame(AIC_1rat_ER = CS_Caste_1rat_ER_AICvec,
                     AIC_1rat_SYM = CS_Caste_1rat_SYM_AICvec,
                     AIC_1rat_ARD = CS_Caste_1rat_ARD_AICvec,
                     AIC_2rat_ER = CS_Caste_2rat_ER_AICvec,
                     AIC_2rat_SYM = CS_Caste_2rat_SYM_AICvec,
                     AIC_2rat_ARD = CS_Caste_2rat_ARD_AICvec,
                     AIC_3rat_ER = CS_Caste_3rat_ER_AICvec,
                     AIC_3rat_SYM = CS_Caste_3rat_SYM_AICvec,
                     AIC_3rat_ARD = CS_Caste_3rat_ARD_AICvec)
apply(AIC_df, 2, mean)

#Identify the best model for each individual tree
# Identify column with smallest value for each row
smallest_column <- apply(AIC_df, 1, function(row) {
  col_index <- which.min(row)
  colnames(AIC_df)[col_index]
})
table(smallest_column)

#Best model: AIC_2rat_SYM - 58.5% of trees

#################

