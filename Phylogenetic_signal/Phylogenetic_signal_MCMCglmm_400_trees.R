########################Calculate measure of phylogenetic signal (heritability) - Colony size and Queen mating frequency######################
####Calculate over 400 trees using MCMCglmm
###Run MCMCglmm models in parallel
##Designed to run on West group sever
#Louis Bell-Roberts
#21/02/2024

.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))


#Load packages
library(tidyverse)
library(ape)
library(phytools)
library(ggplot2)
library(geiger)
library(doParallel)
library(MCMCglmm)
library(coda)

registerDoParallel(40)


#Read in data file
ant_data <- read.csv("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/Data/Trait_database/ant_data.csv")

#Set variables so that they're in the correct structure and apply transformations
class(ant_data$colony.size) #numeric
class(ant_data$queen.mating.frequency) #numeric
ant_data$colony.size <- log10(ant_data$colony.size)
ant_data$queen.mating.frequency <- log10(ant_data$queen.mating.frequency)

#Rename 'species' column as 'animal'
ant_data <- ant_data %>% dplyr::rename(animal = species)

#Read in sample of 400 phylogenetic trees
ant_trees <- read.tree(file ="/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/Data/15k_Economo_trees/Economo_2018_400.tre")

########################################## 

###################### 
#Signal in colony size
###################### 

########################################## 


#Filter data
sxData_CS <- ant_data %>% dplyr::filter(complete.cases(colony.size), complete.cases(caste.number)) %>% dplyr::select(animal, colony.size) #469 species

#Prune tree
ant_trees_pruned_CS <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, sxData_CS$animal))

#Prune database
sxData_CS <- dplyr::filter(sxData_CS, animal %in% ant_trees_pruned_CS[[1]]$tip.label) # 436 species

# check that species names in the tree and the data match
sxData_CS$animal[which((sxData_CS$animal %in% ant_trees_pruned_CS[[1]]$tip.label) == FALSE)]	# 0 means all match
ant_trees_pruned_CS[[1]]$tip.label[which((ant_trees_pruned_CS[[1]]$tip.label %in% sxData_CS$animal) == FALSE)]	# 0 means all match

#########
#MCMCglmm models
#########

##Set prior
# B is for the fixed effect to help with mixing
# R is for residual variance
# G is the phylogenetic/additive genetic variance
prior1 <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)

##############################
#Running models in parallel
##############################

# # Parallel models with outputs saved as .rds files
# foreach(i = 1:400) %dopar% {
#   # Create the MCMCglmm
#   ##1st model
#   model1 <- MCMCglmm(colony.size ~ 1, random = ~animal, family = "gaussian", prior = prior1, pedigree = ant_trees_pruned_CS[[i]], data = sxData_CS, nitt = 1100000, burnin = 100000, thin = 1000, verbose = F)
#   
#   #Run 2nd model
#   model2 <- MCMCglmm(colony.size ~ 1, random = ~animal, family = "gaussian", prior = prior1, pedigree = ant_trees_pruned_CS[[i]], data = sxData_CS, nitt = 1100000, burnin = 100000, thin = 1000, verbose = F)
#   
#   # Save the model as an .rds file for 1st and 2nd chain
#   saveRDS(model1, file = file.path("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/Phylogenetic_signal/Model_outputs/CS/1st_chain", paste0("CS_phylo_sig_1M_100k_1k_", i, "_1stRun.rds")))
#   saveRDS(model2, file = file.path("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/Phylogenetic_signal/Model_outputs/CS/2nd_chain", paste0("CS_phylo_sig_1M_100k_1k_", i, "_2ndRun.rds")))
# }

#Read in the MCMCglmm objects
setwd("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/Phylogenetic_signal/Model_outputs/CS/1st_chain/")
# get list of all RDS files in the working directory
model_files_CS <- list.files(pattern = "\\.rds$")
# read in all models using lapply()
mcmc_list_CS <- lapply(model_files_CS, readRDS)


##############################
#Combine posterior distribution of the random effects - VCV
##############################

# extract the VCV objects from each model output in the list and create a list of the mcmc objects
vcv_matrices_CS <- lapply(mcmc_list_CS, function(model) model$VCV)
#Combine the mcmc objects together
combined_vcv_CS <- do.call(rbind, vcv_matrices_CS)
#Assign as an mcmc object
MCMC_combined_CS <- as.mcmc(combined_vcv_CS)


#Calculate phylogenetic signal (heritability)
herit_CS <- MCMC_combined_CS[ , "animal"] / (MCMC_combined_CS[ , "animal"] + MCMC_combined_CS[ , "units"])

effectiveSize(herit_CS)

posterior.mode(herit_CS) #1st run: 0.7260537; CI = 0.5967257 to 0.8402467

posterior.mode(MCMC_combined_CS[ , "animal"]) #1.207741

posterior.mode(MCMC_combined_CS[ , "units"]) #0.4736025

#Compute 95% credible interval
HPDinterval(herit_CS) #0.5967257 to 0.8402467

HPDinterval(MCMC_combined_CS[ , "animal"]) #0.7972802 to 1.798957

HPDinterval(MCMC_combined_CS[ , "units"]) #0.3359746 to 0.6142431









########################################## 

#################################
#Signal in queen mating frequency
#################################

########################################## 

#Filter data
sxData_MF <- ant_data %>% dplyr::filter(complete.cases(queen.mating.frequency), complete.cases(caste.number)) %>% dplyr::select(animal, queen.mating.frequency) #110 species

#Prune tree
ant_trees_pruned_MF <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, sxData_MF$animal))

#Prune database
sxData_MF <- dplyr::filter(sxData_MF, animal %in% ant_trees_pruned_MF[[1]]$tip.label) #104 species

# check that species names in the tree and the data match
sxData_MF$animal[which((sxData_MF$animal %in% ant_trees_pruned_MF[[1]]$tip.label) == FALSE)]	# 0 means all match
ant_trees_pruned_MF[[1]]$tip.label[which((ant_trees_pruned_MF[[1]]$tip.label %in% sxData_MF$animal) == FALSE)]	# 0 means all match

#########
#MCMCglmm models
#########

##Set prior
# B is for the fixed effect to help with mixing
# R is for residual variance
# G is the phylogenetic/additive genetic variance
prior1 <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)


##############################
#Running models in parallel
##############################

# # Parallel models with outputs saved as .rds files
# foreach(i = 1:400) %dopar% {
#   # Create the MCMCglmm
#   ##1st model
#   model1 <- MCMCglmm(queen.mating.frequency ~ 1, random = ~animal, family = "gaussian", prior = prior1, pedigree = ant_trees_pruned_MF[[i]], data = sxData_MF, nitt = 1100000, burnin = 100000, thin = 1000, verbose = F)
#   
#   #Run 2nd model
#   model2 <- MCMCglmm(queen.mating.frequency ~ 1, random = ~animal, family = "gaussian", prior = prior1, pedigree = ant_trees_pruned_MF[[i]], data = sxData_MF, nitt = 1100000, burnin = 100000, thin = 1000, verbose = F)
#   
#   # Save the model as an .rds file for 1st and 2nd chain
#   saveRDS(model1, file = file.path("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/Phylogenetic_signal/Model_outputs/MF/1st_chain", paste0("MF_phylo_sig_1M_100k_1k_", i, "_1stRun.rds")))
#   saveRDS(model2, file = file.path("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/Phylogenetic_signal/Model_outputs/MF/2nd_chain", paste0("MF_phylo_sig_1M_100k_1k_", i, "_2ndRun.rds")))
# }

#Read in the MCMCglmm objects
setwd("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/Phylogenetic_signal/Model_outputs/MF/1st_chain/")
# get list of all RDS files in the working directory
model_files_MF <- list.files(pattern = "\\.rds$")
# read in all models using lapply()
mcmc_list_MF <- lapply(model_files_MF, readRDS)


##############################
#Combine posterior distribution of the random effects - VCV
##############################

# extract the VCV objects from each model output in the list and create a list of the mcmc objects
vcv_matrices_MF <- lapply(mcmc_list_MF, function(model) model$VCV)
#Combine the mcmc objects together
combined_vcv_MF <- do.call(rbind, vcv_matrices_MF)
#Assign as an mcmc object
MCMC_combined_MF <- as.mcmc(combined_vcv_MF)


#Calculate phylogenetic signal (heritability)
herit_MF <- MCMC_combined_MF[ , "animal"] / (MCMC_combined_MF[ , "animal"] + MCMC_combined_MF[ , "units"])

posterior.mode(herit_MF) #1st run: 0.8894263; CI = 0.6787336 to 0.9973513

posterior.mode(MCMC_combined_MF[ , "animal"]) #1st run: 0.108262

posterior.mode(MCMC_combined_MF[ , "units"]) #1st run: 0.01700555

#Compute 95% credible interval
HPDinterval(herit_MF) #0.7401037 to 0.9985775

HPDinterval(MCMC_combined_MF[ , "animal"]) #0.06341655 to 0.1966299

HPDinterval(MCMC_combined_MF[ , "units"]) #0.0003969985 to 0.03021275






