############################### Regression analyses using BPMMs ########################################
#Louis Bell-Roberts
#15/11/2023

#Packages
library(tidyverse)
library(ape)
library(phytools)
library(MCMCglmm)
library(mulTree)
library(readxl)

#Read in data file
ant_data <- as.data.frame(read_xlsx("ant_data.xlsx", col_names = T))

#Set variables so that they're in the correct structure and apply transformations
ant_data[ant_data == ""] <- NA #Replace blank by NA
class(ant_data$colony.size) #numeric
class(ant_data$queen.mating.frequency) #numeric
ant_data$queen.number <- as.factor(ant_data$queen.number) #Assign as a factor
class(ant_data$queen.number) #factor
class(ant_data$worker.size.variation) #numeric
class(ant_data$caste.number) #numeric

ant_data$colony.size <- log10(ant_data$colony.size)
ant_data$queen.mating.frequency <- log10(ant_data$queen.mating.frequency)
ant_data$worker.size.variation <- sqrt(ant_data$worker.size.variation)

#Rename 'species' column as 'animal'
ant_data <- ant_data %>% dplyr::rename(animal = species)

#Read in sample of 400 phylogenetic trees
ant.trees <- read.tree(file = "Economo_2018_400.tre")

###Subsets of the variables
##Analysis predicting number of worker castes
data_MF_caste <- dplyr::filter(ant_data, queen.mating.frequency >=0, caste.number >=1)
data_CS_caste <- dplyr::filter(ant_data, colony.size >=0, caste.number >=1)
data_PGbinary_caste <- dplyr::filter(ant_data, complete.cases(queen.number), caste.number >=1)

##Analysis predicting variation in worker size
data_MF_siz_var <- dplyr::filter(ant_data, queen.mating.frequency >=0, worker.size.variation >=0)
data_CS_siz_var <- dplyr::filter(ant_data, colony.size >=0, worker.size.variation >=0)
data_PGbinary_siz_var <- dplyr::filter(ant_data, complete.cases(queen.number), worker.size.variation >=0)

##Analysis predicting variation in worker size in species with a single worker caste
data_MF_mono_siz_var <- dplyr::filter(ant_data, queen.mating.frequency >=0, worker.size.variation >=0, caste.number <2)
data_CS_mono_siz_var <- dplyr::filter(ant_data, colony.size >=0, worker.size.variation >=0, caste.number <2)
data_PGbinary_mono_siz_var <- dplyr::filter(ant_data, complete.cases(queen.number), worker.size.variation >=0, caste.number <2)

##Pairwise analysis among predictor variables
data_MF_CS <- dplyr::filter(ant_data, queen.mating.frequency >=0, colony.size >=0)
data_MF_PGbinary <- dplyr::filter(ant_data, queen.mating.frequency >=0, complete.cases(queen.number))
data_CS_PGbinary <- dplyr::filter(ant_data, colony.size >=0, complete.cases(queen.number))

##Prune multiphylo objects for each of the different sets of predictor variables
#Analyses predicting caste number
PT_data_MF_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_caste$animal))
PT_data_CS_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_caste$animal))
PT_data_PGbinary_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PGbinary_caste$animal))
#Analyses predicting variation in worker size
PT_data_MF_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_siz_var$animal))
PT_data_CS_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_siz_var$animal))
PT_data_PGbinary_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PGbinary_siz_var$animal))
#Analyses predicting variation in worker size in species with a single worker caste
PT_data_MF_mono_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_mono_siz_var$animal))
PT_data_CS_mono_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_mono_siz_var$animal))
PT_data_PGbinary_mono_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PGbinary_mono_siz_var$animal))
#Pairwise analyses among the predictor variables
PT_data_MF_CS <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_CS$animal))
PT_data_MF_PGbinary <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_PGbinary$animal))
PT_data_CS_PGbinary <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_PGbinary$animal))

##Prune database to match the tree
#Filter through my dataframe and select only the rows that match the tips of my tree

#Analyses predicting caste number
data_MF_caste <- filter(data_MF_caste, animal %in% PT_data_MF_caste[[1]]$tip.label)
data_CS_caste <- filter(data_CS_caste, animal %in% PT_data_CS_caste[[1]]$tip.label)
data_PGbinary_caste <- filter(data_PGbinary_caste, animal %in% PT_data_PGbinary_caste[[1]]$tip.label)
#Analyses predicting variation in worker size
data_MF_siz_var <- filter(data_MF_siz_var, animal %in% PT_data_MF_siz_var[[1]]$tip.label)
data_CS_siz_var <- filter(data_CS_siz_var, animal %in% PT_data_CS_siz_var[[1]]$tip.label)
data_PGbinary_siz_var <- filter(data_PGbinary_siz_var, animal %in% PT_data_PGbinary_siz_var[[1]]$tip.label)
#Analyses predicting variation in worker size in species with a single worker caste
data_MF_mono_siz_var <- filter(data_MF_mono_siz_var, animal %in% PT_data_MF_mono_siz_var[[1]]$tip.label)
data_CS_mono_siz_var <- filter(data_CS_mono_siz_var, animal %in% PT_data_CS_mono_siz_var[[1]]$tip.label)
data_PGbinary_mono_siz_var <- filter(data_PGbinary_mono_siz_var, animal %in% PT_data_PGbinary_mono_siz_var[[1]]$tip.label)
#Pairwise analyses among the predictor variables
data_MF_CS <- filter(data_MF_CS, animal %in% PT_data_MF_CS[[1]]$tip.label)
data_MF_PGbinary <- filter(data_MF_PGbinary, animal %in% PT_data_MF_PGbinary[[1]]$tip.label)
data_CS_PGbinary <- filter(data_CS_PGbinary, animal %in% PT_data_CS_PGbinary[[1]]$tip.label)



############################################################
#Phylogenetic regression analyses
############################################################

##Model settings
#Set prior
prior1 <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)
#Parameter expanded prior
prior.par.exp <- list(R = list(V = 1, nu = 1),
                    G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

#Set up model parameters
# The MCMC parameters (iterations, thin, burnin)
##1 million
nitt1M <- 1100000
burnin100K = 100000
thin1K = 1000
mul_parameters1M <- c(nitt1M, thin1K, burnin100K)

##400k
nitt400k <- 400000
burnin20k = 20000
thin200 = 240
mul_parameters400k <- c(nitt400k, thin200, burnin20k)



############################################################
#Pairwise analyses predicting caste number using MCMCglmm
############################################################

############
#Caste ~ CS#
############
#Prepare the data for the mulTree function
mulTree_data_CS_Caste <- as.mulTree(data = data_CS_caste, tree = PT_data_CS_caste,
                                    taxa = "animal")

CS_Caste <- caste.number ~ colony.size
mulTree(mulTree.data = mulTree_data_CS_Caste, formula = CS_Caste, priors = prior.par.exp,
        parameters = mul_parameters1M, output = "CS_Caste_400k", ESS = 200,
        chains = 2, family = "poisson")


############
#Caste ~ MF#
############
mulTree_data_MF_Caste <- as.mulTree(data = data_MF_caste, tree = PT_data_MF_caste,
                                    taxa = "animal")

MF_Caste <- caste.number ~ queen.mating.frequency
mulTree(mulTree.data = mulTree_data_MF_Caste, formula = MF_Caste, priors = prior.par.exp,
        parameters = mul_parameters1M, output = "MF_Caste_1M", ESS = 200,
        chains = 2, family = "poisson")


###############
#Caste ~ PGbinary#
###############
mulTree_data_PGbin_Caste <- as.mulTree(data = data_PGbinary_caste, tree = PT_data_PGbinary_caste,
                                       taxa = "animal")

PGbin_Caste <- caste.number ~ queen.number
mulTree(mulTree.data = mulTree_data_PGbin_Caste, formula = PGbin_Caste, priors = prior.par.exp,
        parameters = mul_parameters1M, output = "PGbin_Caste_1M", ESS = 200,
        chains = 2, family = "poisson")



############################################################
#Pairwise analyses predicting variation in worker size using MCMCglmm
############################################################

##########
#size var ~ CS#
##########
mulTree_data_CS_siz_var <- as.mulTree(data = data_CS_siz_var, tree = PT_data_CS_siz_var,
                                   taxa = "animal")

CS_siz_var <- worker.size.variation ~ colony.size
mulTree(mulTree.data = mulTree_data_CS_siz_var, formula = CS_siz_var, priors = prior1,
        parameters = mul_parameters1M, output = "CS_siz_var_1M", ESS = 200,
        chains = 2, family = "gaussian")



##########
#size var ~ MF#
##########
mulTree_data_MF_siz_var <- as.mulTree(data = data_MF_siz_var, tree = PT_data_MF_siz_var,
                                   taxa = "animal")

MF_siz_var <- worker.size.variation ~ queen.mating.frequency
mulTree(mulTree.data = mulTree_data_MF_poly, formula = MF_poly, priors = prior1,
        parameters = mul_parameters1M, output = "MF_siz_var_1M", ESS = 200,
        chains = 2, family = "gaussian")


##########
#size var ~ PGbinary#
##########
mulTree_data_PGbin_siz_var <- as.mulTree(data = data_PGbinary_siz_var, tree = PT_data_PGbinary_siz_var,
                                      taxa = "animal")

PGbin_siz_var <- worker.size.variation ~ queen.number
mulTree(mulTree.data = mulTree_data_PGbin_siz_var, formula = PGbin_siz_var, priors = prior1,
        parameters = mul_parameters1M, output = "PGbin_siz_var_1M", ESS = 200,
        chains = 2, family = "gaussian")




############################################################
#Pairwise analyses predicting variation in worker size in species with a single worker caste using MCMCglmm
############################################################

##########
#size var ~ CS#
##########
mulTree_data_CS_mono_siz_var <- as.mulTree(data = data_CS_mono_siz_var, tree = PT_data_CS_mono_siz_var,
                                   taxa = "animal")

CS_mono_siz_var <- worker.size.variation ~ colony.size
mulTree(mulTree.data = mulTree_data_CS_mono_siz_var, formula = CS_mono_siz_var, priors = prior1,
        parameters = mul_parameters1M, output = "CS_mono_siz_var_1M", ESS = 200,
        chains = 2, family = "gaussian")

##########
#size var ~ MF#
##########
mulTree_data_MF_mono_siz_var <- as.mulTree(data = data_MF_mono_siz_var, tree = PT_data_MF_mono_siz_var,
                                   taxa = "animal")

MF_mono_siz_var <- worker.size.variation ~ queen.mating.frequency
mulTree(mulTree.data = mulTree_data_MF_mono_siz_var, formula = MF_mono_siz_var, priors = prior1,
        parameters = mul_parameters1M, output = "MF_mono_siz_var_1M", ESS = 200,
        chains = 2, family = "gaussian")


##########
#size var ~ PGbinary#
##########
mulTree_data_PGbin_mono_siz_var <- as.mulTree(data = data_PGbinary_mono_siz_var, tree = PT_data_PGbinary_mono_siz_var,
                                      taxa = "animal")

PGbin_mono_siz_var <- worker.size.variation ~ queen.number
mulTree(mulTree.data = mulTree_data_PGbin_mono_siz_var, formula = PGbin_mono_siz_var, priors = prior1,
        parameters = mul_parameters1M, output = "PG_mono_siz_var_1M", ESS = 200,
        chains = 2, family = "gaussian")



################################################################
#Pairwise analyses among the predictor variables using MCMCglmm
################################################################

##########
#MF ~ CS#
##########
mulTree_data_MF_CS <- as.mulTree(data = data_MF_CS, tree = PT_data_MF_CS,
                                 taxa = "animal")

MF_CS <- colony.size ~ queen.mating.frequency
mulTree(mulTree.data = mulTree_data_MF_CS, formula = MF_CS, priors = prior1,
        parameters = mul_parameters1M, output = "MF_CS_1M", ESS = 200,
        chains = 2, family = "gaussian")


##########
#PGbin ~ CS#
##########
mulTree_data_CS_PGbin <- as.mulTree(data = data_CS_PGbinary, tree = PT_data_CS_PGbinary,
                                    taxa = "animal")

PGbin_CS <- colony.size ~ queen.number
mulTree(mulTree.data = mulTree_data_CS_PGbin, formula = PGbin_CS, priors = prior1,
        parameters = mul_parameters1M, output = "CS_PGbin_400k", ESS = 200,
        chains = 2, family = "gaussian")



##########
#PGbin ~ MF#
##########
mulTree_data_MF_PGbin <- as.mulTree(data = data_MF_PGbinary, tree = PT_data_MF_PGbinary,
                                    taxa = "animal")

PGbin_MF <- queen.mating.frequency ~ queen.number
mulTree(mulTree.data = mulTree_data_MF_PGbin, formula = PGbin_MF, priors = prior1,
        parameters = mul_parameters1M, output = "MF_PGbin_400k", ESS = 200,
        chains = 2, family = "gaussian")


############ END ##########################


