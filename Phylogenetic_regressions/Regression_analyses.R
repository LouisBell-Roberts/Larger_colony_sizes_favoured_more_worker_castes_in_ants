############################### BPMM analyses using BPMMs ########################################
#Louis Bell-Roberts
#14/01/2024

#############             #############               #############             

#Packages
library(tidyverse)
library(ape)
library(phytools)
library(MCMCglmm)
library(doParallel)
library(coda)

#Use 50 cores for parallel computation
registerDoParallel(50)

#Read in data file
ant_data <- read.csv("ant_data.csv")

#Set variables so that they're in the correct structure and apply transformations
class(ant_data$colony.size) #numeric
class(ant_data$queen.mating.frequency) #numeric
ant_data$queen.mating.frequency.categorical <- as.factor(ant_data$queen.mating.frequency.categorical)
class(ant_data$queen.mating.frequency.categorical) #factor
class(ant_data$queen.number.continuous) #numeric
ant_data$queen.number.binary <- as.factor(ant_data$queen.number.binary) #Assign as a factor
class(ant_data$queen.number.binary) #factor
ant_data$queen.number.categorical <- as.factor(ant_data$queen.number.categorical) #Assign as a factor
class(ant_data$queen.number.categorical) #factor
class(ant_data$worker.size.variation) #numeric
class(ant_data$caste.number) #numeric

ant_data$colony.size <- log10(ant_data$colony.size)
ant_data$queen.mating.frequency <- log10(ant_data$queen.mating.frequency)
ant_data$queen.number.continuous <- log10(ant_data$queen.number.continuous)
ant_data$worker.size.variation <- sqrt(ant_data$worker.size.variation)

#Rename 'species' column as 'animal'
ant_data <- ant_data %>% dplyr::rename(animal = species)

#Read in sample of 400 phylogenetic trees
ant.trees <- read.tree(file ="Economo_2018_400.tre")

###Subsets of the variables
##Analysis predicting number of worker castes
data_MF_caste <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(caste.number))
data_MFcat_caste <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency.categorical), complete.cases(caste.number))
data_CS_caste <- dplyr::filter(ant_data, complete.cases(colony.size), complete.cases(caste.number))
data_PG_caste <- dplyr::filter(ant_data, complete.cases(queen.number.continuous), complete.cases(caste.number))
data_PGbinary_caste <- dplyr::filter(ant_data, complete.cases(queen.number.binary), complete.cases(caste.number))
data_PGcat_caste <- dplyr::filter(ant_data, complete.cases(queen.number.categorical), complete.cases(caste.number))

##Analysis predicting variation in worker size
data_MF_siz_var <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(worker.size.variation))
data_MFcat_siz_var <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency.categorical), complete.cases(worker.size.variation))
data_CS_siz_var <- dplyr::filter(ant_data, complete.cases(colony.size), complete.cases(worker.size.variation))
data_PG_siz_var <- dplyr::filter(ant_data, complete.cases(queen.number.continuous), complete.cases(worker.size.variation))
data_PGbinary_siz_var <- dplyr::filter(ant_data, complete.cases(queen.number.binary), complete.cases(worker.size.variation))
data_PGcat_siz_var <- dplyr::filter(ant_data, complete.cases(queen.number.categorical), complete.cases(worker.size.variation))

##Analysis predicting variation in worker size in species with a single worker caste
data_MF_mono_siz_var <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(worker.size.variation), caste.number <2)
data_MFcat_mono_siz_var <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency.categorical), complete.cases(worker.size.variation), caste.number <2)
data_CS_mono_siz_var <- dplyr::filter(ant_data, complete.cases(colony.size), complete.cases(worker.size.variation), caste.number <2)
data_PG_mono_siz_var <- dplyr::filter(ant_data, complete.cases(queen.number.continuous), complete.cases(worker.size.variation), caste.number <2)
data_PGbinary_mono_siz_var <- dplyr::filter(ant_data, complete.cases(queen.number.binary), complete.cases(worker.size.variation), caste.number <2)
data_PGcat_mono_siz_var <- dplyr::filter(ant_data, complete.cases(queen.number.categorical), complete.cases(worker.size.variation), caste.number <2)

##Pairwise analysis among predictor variables
data_MF_CS <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(colony.size))
data_MFcat_CS <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency.categorical), complete.cases(colony.size))
data_MF_PG <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(queen.number.continuous))
data_MF_PGbinary <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(queen.number.binary))
data_MF_PGcat <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(queen.number.categorical))
data_MFcat_PG <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency.categorical), complete.cases(queen.number.continuous))
data_CS_PG <- dplyr::filter(ant_data, complete.cases(colony.size), complete.cases(queen.number.continuous))
data_CS_PGbinary <- dplyr::filter(ant_data, complete.cases(colony.size), complete.cases(queen.number.binary))
data_CS_PGcat <- dplyr::filter(ant_data, complete.cases(colony.size), complete.cases(queen.number.categorical))

##Prune multiphylo objects for each of the different sets of predictor variables
#Analyses predicting caste number
PT_data_MF_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_caste$animal))
PT_data_MFcat_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MFcat_caste$animal))
PT_data_CS_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_caste$animal))
PT_data_PG_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PG_caste$animal))
PT_data_PGbinary_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PGbinary_caste$animal))
PT_data_PGcat_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PGcat_caste$animal))

#Analyses predicting variation in worker size
PT_data_MF_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_siz_var$animal))
PT_data_MFcat_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MFcat_siz_var$animal))
PT_data_CS_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_siz_var$animal))
PT_data_PG_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PG_siz_var$animal))
PT_data_PGbinary_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PGbinary_siz_var$animal))
PT_data_PGcat_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PGcat_siz_var$animal))

#Analyses predicting variation in worker size in species with a single worker caste
PT_data_MF_mono_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_mono_siz_var$animal))
PT_data_MFcat_mono_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MFcat_mono_siz_var$animal))
PT_data_CS_mono_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_mono_siz_var$animal))
PT_data_PG_mono_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PG_mono_siz_var$animal))
PT_data_PGbinary_mono_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PGbinary_mono_siz_var$animal))
PT_data_PGcat_mono_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PGcat_mono_siz_var$animal))

#Pairwise analyses among the predictor variables
PT_data_MF_CS <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_CS$animal))
PT_data_MFcat_CS <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MFcat_CS$animal))
PT_data_MF_PG <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_PG$animal))
PT_data_MF_PGbinary <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_PGbinary$animal))
PT_data_MF_PGcat <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_PGcat$animal))
PT_data_MFcat_PG <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MFcat_PG$animal))
PT_data_CS_PG <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_PG$animal))
PT_data_CS_PGbinary <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_PGbinary$animal))
PT_data_CS_PGcat <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_PGcat$animal))

##Prune database to match the tree
#Filter through dataframe and select only the rows that match the tips of the tree

#Analyses predicting caste number
data_MF_caste <- filter(data_MF_caste, animal %in% PT_data_MF_caste[[1]]$tip.label)
data_MFcat_caste <- filter(data_MFcat_caste, animal %in% PT_data_MFcat_caste[[1]]$tip.label)
data_CS_caste <- filter(data_CS_caste, animal %in% PT_data_CS_caste[[1]]$tip.label)
data_PG_caste <- filter(data_PG_caste, animal %in% PT_data_PG_caste[[1]]$tip.label)
data_PGbinary_caste <- filter(data_PGbinary_caste, animal %in% PT_data_PGbinary_caste[[1]]$tip.label)
data_PGcat_caste <- filter(data_PGcat_caste, animal %in% PT_data_PGcat_caste[[1]]$tip.label)

#Analyses predicting variation in worker size
data_MF_siz_var <- filter(data_MF_siz_var, animal %in% PT_data_MF_siz_var[[1]]$tip.label)
data_MFcat_siz_var <- filter(data_MFcat_siz_var, animal %in% PT_data_MFcat_siz_var[[1]]$tip.label)
data_CS_siz_var <- filter(data_CS_siz_var, animal %in% PT_data_CS_siz_var[[1]]$tip.label)
data_PG_siz_var <- filter(data_PG_siz_var, animal %in% PT_data_PG_siz_var[[1]]$tip.label)
data_PGbinary_siz_var <- filter(data_PGbinary_siz_var, animal %in% PT_data_PGbinary_siz_var[[1]]$tip.label)
data_PGcat_siz_var <- filter(data_PGcat_siz_var, animal %in% PT_data_PGcat_siz_var[[1]]$tip.label)

#Analyses predicting variation in worker size in species with a single worker caste
data_MF_mono_siz_var <- filter(data_MF_mono_siz_var, animal %in% PT_data_MF_mono_siz_var[[1]]$tip.label)
data_MFcat_mono_siz_var <- filter(data_MFcat_mono_siz_var, animal %in% PT_data_MFcat_mono_siz_var[[1]]$tip.label)
data_CS_mono_siz_var <- filter(data_CS_mono_siz_var, animal %in% PT_data_CS_mono_siz_var[[1]]$tip.label)
data_PG_mono_siz_var <- filter(data_PG_mono_siz_var, animal %in% PT_data_PG_mono_siz_var[[1]]$tip.label)
data_PGbinary_mono_siz_var <- filter(data_PGbinary_mono_siz_var, animal %in% PT_data_PGbinary_mono_siz_var[[1]]$tip.label)
data_PGcat_mono_siz_var <- filter(data_PGcat_mono_siz_var, animal %in% PT_data_PGcat_mono_siz_var[[1]]$tip.label)

#Pairwise analyses among the predictor variables
data_MF_CS <- filter(data_MF_CS, animal %in% PT_data_MF_CS[[1]]$tip.label)
data_MFcat_CS <- filter(data_MFcat_CS, animal %in% PT_data_MFcat_CS[[1]]$tip.label)
data_MF_PG <- filter(data_MF_PG, animal %in% PT_data_MF_PG[[1]]$tip.label)
data_MF_PGbinary <- filter(data_MF_PGbinary, animal %in% PT_data_MF_PGbinary[[1]]$tip.label)
data_MF_PGcat <- filter(data_MF_PGcat, animal %in% PT_data_MF_PGcat[[1]]$tip.label)
data_MFcat_PG <- filter(data_MFcat_PG, animal %in% PT_data_MFcat_PG[[1]]$tip.label)
data_CS_PG <- filter(data_CS_PG, animal %in% PT_data_CS_PG[[1]]$tip.label)
data_CS_PGbinary <- filter(data_CS_PGbinary, animal %in% PT_data_CS_PGbinary[[1]]$tip.label)
data_CS_PGcat <- filter(data_CS_PGcat, animal %in% PT_data_CS_PGcat[[1]]$tip.label)


############################################################
#Phylogenetic regression analyses
############################################################
#Model setup

#Define function for running each model
##Parallel models where the outputs are saved as .rds files
###Loop through each element in the list of trees and perform parallel computation
runMCMCglmm <- function(response, predictor, family, prior, multiphylo, data, path) {
  
  # Parallel models with outputs saved as .rds files
  foreach(i = 1:400) %dopar% {
    # Create the MCMCglmm
    model <- MCMCglmm(fixed = as.formula(paste(response, "~", predictor)), random = ~animal, family = family, prior = prior, pedigree = multiphylo[[i]], data = data, nitt = 1100000, burnin = 100000, thin = 1000, verbose = F)
    
    # Save the model as an .rds file for 1st and 2nd chain
    saveRDS(model, file = file.path(path, "1st_chain", paste0(path, "_1M_100k_1k_", i, "_1stRun.rds")))
    saveRDS(model, file = file.path(path, "2nd_chain", paste0(path, "_1M_100k_1k_", i, "_2ndRun.rds")))
  }
}

#Set prior
#Default prior
prior1 <- list(R = list(V = 1, nu = 0.002),
               G = list(G1 = list(V = 1, nu = 0.002))
)
#Parameter expanded prior
prior.par.exp <- list(R = list(V = 1, nu = 1),
                      G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

############################################################
#Pairwise analyses predicting caste number using MCMCglmm
############################################################

##############################################################################
#Analyses predicting number of worker castes

############
#MF_Caste#
############
runMCMCglmm(response = "caste.number", predictor = "queen.mating.frequency", family = "poisson", prior = prior.par.exp, multiphylo = PT_data_MF_caste, data = data_MF_caste, path = "MF_Caste")

############
#MFcat_Caste#
############
runMCMCglmm(response = "caste.number", predictor = "queen.mating.frequency.categorical", family = "poisson", prior = prior.par.exp, multiphylo = PT_data_MFcat_caste, data = data_MFcat_caste, path = "MFcat_Caste")

############
#CS_Caste#
############
runMCMCglmm(response = "caste.number", predictor = "colony.size", family = "poisson", prior = prior.par.exp, multiphylo = PT_data_CS_caste, data = data_CS_caste, path = "CS_Caste")

############
#PG_Caste#
############
runMCMCglmm(response = "caste.number", predictor = "queen.number.continuous", family = "poisson", prior = prior.par.exp, multiphylo = PT_data_PG_caste, data = data_PG_caste, path = "PG_Caste")

############
#PGbinary_Caste#
############
runMCMCglmm(response = "caste.number", predictor = "queen.number.binary", family = "poisson", prior = prior.par.exp, multiphylo = PT_data_PGbinary_caste, data = data_PGbinary_caste, path = "PGbinary_Caste")

############
#PGcat_Caste#
#############
runMCMCglmm(response = "caste.number", predictor = "queen.number.categorical", family = "poisson", prior = prior.par.exp, multiphylo = PT_data_PGcat_caste, data = data_PGcat_caste, path = "PGcat_Caste")

##############################################################################
#Analyses predicting variation in worker size

############
#MF_siz_var#
############
runMCMCglmm(response = "worker.size.variation", predictor = "queen.mating.frequency", family = "gaussian", prior = prior1, multiphylo = PT_data_MF_siz_var, data = data_MF_siz_var, path = "MF_siz_var")

############
#MFcat_siz_var#
############
runMCMCglmm(response = "worker.size.variation", predictor = "queen.mating.frequency.categorical", family = "gaussian", prior = prior1, multiphylo = PT_data_MFcat_siz_var, data = data_MFcat_siz_var, path = "MFcat_siz_var")

############
#CS_siz_var#
############
runMCMCglmm(response = "worker.size.variation", predictor = "colony.size", family = "gaussian", prior = prior1, multiphylo = PT_data_CS_siz_var, data = data_CS_siz_var, path = "CS_siz_var")

############
#PG_siz_var#
############
runMCMCglmm(response = "worker.size.variation", predictor = "queen.number.continuous", family = "gaussian", prior = prior1, multiphylo = PT_data_PG_siz_var, data = data_PG_siz_var, path = "PG_siz_var")

############
#PGbinary_siz_var#
############
runMCMCglmm(response = "worker.size.variation", predictor = "queen.number.binary", family = "gaussian", prior = prior1, multiphylo = PT_data_PGbinary_siz_var, data = data_PGbinary_siz_var, path = "PGbinary_siz_var")

############
#PGcat_siz_var#
############
runMCMCglmm(response = "worker.size.variation", predictor = "queen.number.categorical", family = "gaussian", prior = prior1, multiphylo = PT_data_PGcat_siz_var, data = data_PGcat_siz_var, path = "PGcat_siz_var")

##############################################################################
#Analyses predicting variation in worker size in species with a single worker caste

############
#MF_mono_siz_var#
############
runMCMCglmm(response = "worker.size.variation", predictor = "queen.mating.frequency", family = "gaussian", prior = prior1, multiphylo = PT_data_MF_mono_siz_var, data = data_MF_mono_siz_var, path = "MF_mono_siz_var")

############
#MFcat_mono_siz_var#
############
runMCMCglmm(response = "worker.size.variation", predictor = "queen.mating.frequency.categorical", family = "gaussian", prior = prior1, multiphylo = PT_data_MFcat_mono_siz_var, data = data_MFcat_mono_siz_var, path = "MFcat_mono_siz_var")

############
#CS_mono_siz_var#
############
runMCMCglmm(response = "worker.size.variation", predictor = "colony.size", family = "gaussian", prior = prior1, multiphylo = PT_data_CS_mono_siz_var, data = data_CS_mono_siz_var, path = "CS_mono_siz_var")

############
#PG_mono_siz_var#
############
runMCMCglmm(response = "worker.size.variation", predictor = "queen.number.continuous", family = "gaussian", prior = prior1, multiphylo = PT_data_PG_mono_siz_var, data = data_PG_mono_siz_var, path = "PG_mono_siz_var")

############
#PGbinary_mono_siz_var#
############
runMCMCglmm(response = "worker.size.variation", predictor = "queen.number.binary", family = "gaussian", prior = prior1, multiphylo = PT_data_PGbinary_mono_siz_var, data = data_PGbinary_mono_siz_var, path = "PGbinary_mono_siz_var")

############
#PGcat_mono_siz_var#
############
runMCMCglmm(response = "worker.size.variation", predictor = "queen.number.categorical", family = "gaussian", prior = prior1, multiphylo = PT_data_PGcat_mono_siz_var, data = data_PGcat_mono_siz_var, path = "PGcat_mono_siz_var")

##############################################################################
##Pairwise analysis among predictor variables

############
#MF_CS#
############
runMCMCglmm(response = "colony.size", predictor = "queen.mating.frequency", family = "gaussian", prior = prior1, multiphylo = PT_data_MF_CS, data = data_MF_CS, path = "MF_CS")

############
#MFcat_CS#
############
runMCMCglmm(response = "colony.size", predictor = "queen.mating.frequency.categorical", family = "gaussian", prior = prior1, multiphylo = PT_data_MFcat_CS, data = data_MFcat_CS, path = "MFcat_CS")

############
#MF_PG#
############
runMCMCglmm(response = "queen.mating.frequency", predictor = "queen.number.continuous", family = "gaussian", prior = prior1, multiphylo = PT_data_MF_PG, data = data_MF_PG, path = "MF_PG")

############
#MF_PGbinary#
############
runMCMCglmm(response = "queen.mating.frequency", predictor = "queen.number.binary", family = "gaussian", prior = prior1, multiphylo = PT_data_MF_PGbinary, data = data_MF_PGbinary, path = "MF_PGbinary")

############
#MF_PGcat#
############
runMCMCglmm(response = "queen.mating.frequency", predictor = "queen.number.categorical", family = "gaussian", prior = prior1, multiphylo = PT_data_MF_PGcat, data = data_MF_PGcat, path = "MF_PGcat")

############
#MFcat_PG#
############
runMCMCglmm(response = "queen.number.continuous", predictor = "queen.mating.frequency.categorical", family = "gaussian", prior = prior1, multiphylo = PT_data_MFcat_PG, data = data_MFcat_PG, path = "MFcat_PG")

############
#CS_PG#
############
runMCMCglmm(response = "colony.size", predictor = "queen.number.continuous", family = "gaussian", prior = prior1, multiphylo = PT_data_CS_PG, data = data_CS_PG, path = "CS_PG")

############
#CS_PGbinary#
############
runMCMCglmm(response = "colony.size", predictor = "queen.number.binary", family = "gaussian", prior = prior1, multiphylo = PT_data_CS_PGbinary, data = data_CS_PGbinary, path = "CS_PGbinary")

############
#CS_PGcat#
############
runMCMCglmm(response = "colony.size", predictor = "queen.number.categorical", family = "gaussian", prior = prior1, multiphylo = PT_data_CS_PGcat, data = data_CS_PGcat, path = "CS_PGcat")



############ END ##########################


