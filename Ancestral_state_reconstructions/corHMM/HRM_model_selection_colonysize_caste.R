########################ASR with corHMM - analysis across 400 trees######################
###Runs ASR for caste number as a binary variable (single/multiple castes)
##Analysing all species for which data is available for both the number of worker castes as well as colony size
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
corHMM_multiphylo_2rat_ER <- function(multiphylo, trait_data) {
  
  # Run corHMM ancestral state reconstruction for each phylogeny
  for (i in 1:length(multiphylo)) {
    phylo <- multiphylo[[i]]
    output <- corHMM(phy = phylo, data = trait_data, rate.cat = 2, model = "ER")
    
    # Save output to a file
    output_file <- paste0("corHMM_ASR_Caste_CS_2rat_ER_", i, ".rda")
    saveRDS(output, file = output_file)
    cat("Saved output to", output_file, "\n")
    
  }
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
ASRdata <- dplyr::filter(ant_data, caste.number >=1, colony.size >=1)
ASRdata <- dplyr::select(ASRdata, animal, colony.size, CasteBin)

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
corHMM_multiphylo_2rat_ER(ant_trees_pruned, ASRdata)


########END OF SCRIPT##########










