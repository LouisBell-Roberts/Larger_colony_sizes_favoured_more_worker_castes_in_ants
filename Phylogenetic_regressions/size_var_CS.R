############################### PGLMM: Variation in worker size ########################################
###Testing for a relationship between colony size and variation in worker size
#Louis Bell-Roberts
#15/06/2022


library(tidyverse)
library(ape)
library(phytools)
library(MCMCglmm)
library(mulTree)


#Read in data file - ensure that caste variable is set a numeric variable
d <- read.csv("data", header = T)

#Set variables so that they're in the correct structure
d[d == ""] <- NA                     # Replace blank by NA
d$Caste3 <- as.numeric(as.character(d$Caste3))
d$number.queens.MEAN <- as.numeric(as.character(d$number.queens.MEAN))

#Change underscores in species names to full stops so that it matches the Economo tree notation
d$animal <- gsub("_", ".", d$animal, fixed = T) #need to run if using the Economo tree

#Remove species that aren't ants, not on Antcat, that are supercolonial, social parasites, hybridisation caste systems, clonal and queenless
## Filter species in the analysis
data <- dplyr::filter(d, type == 'ant', Taxonomic_info.Antcat.registered_undescribed.species_Antcat.unregistered_incertae.sedis._subspecies <1, Supercolonial.Helantera.review_no.stringent.nonstringent <1, social.parasite.Antwiki_no.yes <1, hybridisation <1, clonal <1, GamergatesAntwiki_yes_no <1, SD_all_workers >=0, All_workers_HW_mean >=0)

#Create worker_CoV variable and square root transform it: 
data$worker_CoV <- c(data$SD_all_workers / data$All_workers_HW_mean)
data$worker_CoV <- sqrt(data$worker_CoV)

#Tree files - 'Genus_polytomy_tree.tre' is the most up-to-date tree and should be used
ant_trees_NCuniform_stem_posterior <- read.tree(file = "/drives/4tb/Louis/Data/15k_Economo_trees/15k_NCuniform_stem_posterior.tre")

ant_trees_NCuniform_crown_posterior <- read.tree(file = "/drives/4tb/Louis/Data/15k_Economo_trees/15k_NCuniform_crown_posterior.tre")

ant_trees_FBD_stem_posterior <- read.tree(file = "/drives/4tb/Louis/Data/15k_Economo_trees/15K_FBD_stem_posterior.tre")

ant_trees_FBD_crown_posterior <- read.tree(file = "/drives/4tb/Louis/Data/15k_Economo_trees/15K_FBD_crown_posterior.tre")

ant.trees <- c(ant_trees_NCuniform_stem_posterior, ant_trees_NCuniform_crown_posterior, ant_trees_FBD_stem_posterior, ant_trees_FBD_crown_posterior)
# ant.trees <- ant.trees[1]

########
#IMPORTANT - log transfrom all of the variables that need it at this stage as the mulTree function does not like these function in the formula ###
########

##Multiple regression analyses##
all_variables_PGcont <- dplyr::filter(data, Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1, number.queens.MEAN >=1)
all_variables_PGcont$eff.mating.freq.MEAN.harmonic <- log10(all_variables_PGcont$eff.mating.freq.MEAN.harmonic)
all_variables_PGcont$colony.size <- log10(all_variables_PGcont$colony.size)
all_variables_PGcont$number.queens.MEAN <- log10(all_variables_PGcont$number.queens.MEAN)

all_variables_PGbinary <- dplyr::filter(data, Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1, polygyny.clean.new >=0)
all_variables_PGbinary$eff.mating.freq.MEAN.harmonic <- log10(all_variables_PGbinary$eff.mating.freq.MEAN.harmonic)
all_variables_PGbinary$colony.size <- log10(all_variables_PGbinary$colony.size)

#all_variables_PGbinary <- dplyr::filter(data, Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1, polygyny.clean.new >=0, Worker_Sterility_0_1 >=0)
all_variables_PGbinary$Worker_Sterility_0_1 <- as.factor(all_variables_PGbinary$Worker_Sterility_0_1)

##Subsets of the variables
#Pairwise analyses predicting caste number
data_MF_caste <- dplyr::filter(data, eff.mating.freq.MEAN.harmonic >=1, Caste3 >=1)
data_MF_caste$eff.mating.freq.MEAN.harmonic <- log10(data_MF_caste$eff.mating.freq.MEAN.harmonic)


data_CS_poly <- dplyr::filter(data, colony.size >=1, worker_CoV >=0)
data_CS_poly$colony.size <- log10(data_CS_poly$colony.size)
#Select only the columns of interest
data_CS_poly <- dplyr::select(data_CS_poly, animal, colony.size, SD_all_workers, All_workers_HW_mean, Sample_size_all_workers, worker_CoV, Caste3)


data_PGcont_caste <- dplyr::filter(data, Caste3 >=1, number.queens.MEAN >=1)
data_PGcont_caste$number.queens.MEAN <- log10(data_PGcont_caste$number.queens.MEAN)

data_PGbinary_caste <- dplyr::filter(data, Caste3 >=1, polygyny.clean.new >=0)

#Pairwise analyses among the predictor variables
data_MF_CS <- dplyr::filter(data, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1)
data_MF_CS$eff.mating.freq.MEAN.harmonic <- log10(data_MF_CS$eff.mating.freq.MEAN.harmonic)
data_MF_CS$colony.size <- log10(data_MF_CS$colony.size)

data_MF_PGcont <- dplyr::filter(data, eff.mating.freq.MEAN.harmonic >=1, number.queens.MEAN >=1)
data_MF_PGcont$eff.mating.freq.MEAN.harmonic <- log10(data_MF_PGcont$eff.mating.freq.MEAN.harmonic)
data_MF_PGcont$number.queens.MEAN <- log10(data_MF_PGcont$number.queens.MEAN)

data_MF_PGbinary <- dplyr::filter(data, eff.mating.freq.MEAN.harmonic >=1, polygyny.clean.new >=0)
data_MF_PGbinary$eff.mating.freq.MEAN.harmonic <- log10(data_MF_PGbinary$eff.mating.freq.MEAN.harmonic)

data_CS_PGcont <- dplyr::filter(data, colony.size >=1, number.queens.MEAN >=1)
data_CS_PGcont$colony.size <- log10(data_CS_PGcont$colony.size)
data_CS_PGcont$number.queens.MEAN <- log10(data_CS_PGcont$number.queens.MEAN)

data_CS_PGbinary <- dplyr::filter(data, colony.size >=1, polygyny.clean.new >=0)
data_CS_PGbinary$colony.size <- log10(data_CS_PGbinary$colony.size)

##Prune trees - for each of the different sets of predictor variables
#These will now be multiphylo objects that have been pruned
#Pairwise analyses predicting caste number
PT_data_MF_caste<-drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_caste$animal))
PT_data_CS_poly<-drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_poly$animal))
PT_data_PGcont_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PGcont_caste$animal))
PT_data_PGbinary_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PGbinary_caste$animal))
#Pairwise analyses among the predictor variables
PT_data_MF_CS<-drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_CS$animal))
PT_data_MF_PGcont<-drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_PGcont$animal))
PT_data_MF_PGbinary<-drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_PGbinary$animal))
PT_data_CS_PGcont<-drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_PGcont$animal))
PT_data_CS_PGbinary<-drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_PGbinary$animal))
#Multiple regression analyses
PT_all_variables_PGcont<-drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, all_variables_PGcont$animal))
PT_all_variables_PGbinary<-drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, all_variables_PGbinary$animal))


##Prune database to match the tree
#Filter through my dataframe and select only the rows that match the tips of my tree
#Pairwise analyses predicting caste number
data_MF_caste <- filter(data_MF_caste, animal %in% PT_data_MF_caste[[1]]$tip.label)
data_CS_poly <- filter(data_CS_poly, animal %in% PT_data_CS_poly[[1]]$tip.label)
data_PGcont_caste <- filter(data_PGcont_caste, animal %in% PT_data_PGcont_caste[[1]]$tip.label)
data_PGbinary_caste <- filter(data_PGbinary_caste, animal %in% PT_data_PGbinary_caste[[1]]$tip.label)
#Pairwise analyses among the predictor variables
data_MF_CS<-filter(data_MF_CS, animal %in% PT_data_MF_CS[[1]]$tip.label)
data_MF_PGcont<-filter(data_MF_PGcont, animal %in% PT_data_MF_PGcont[[1]]$tip.label)
data_MF_PGbinary<-filter(data_MF_PGbinary, animal %in% PT_data_MF_PGbinary[[1]]$tip.label)
data_CS_PGcont <- filter(data_CS_PGcont, animal %in% PT_data_CS_PGcont[[1]]$tip.label)
data_CS_PGbinary <- filter(data_CS_PGbinary, animal %in% PT_data_CS_PGbinary[[1]]$tip.label)
#Multiple regression analyses
all_variables_PGcont<-filter(all_variables_PGcont, animal %in% PT_all_variables_PGcont[[1]]$tip.label)
all_variables_PGbinary<-filter(all_variables_PGbinary, animal %in% PT_all_variables_PGbinary[[1]]$tip.label)

##Phylogenetic regression analyses

############################################################
#Pairwise analyses predicting caste number using MCMCglmm
############################################################

##Set prior
# B is for the fixed effect to help with mixing
# R is for residual variance
# G is the phylogenetic/additive genetic variance
prior1 <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002)
)
#Parameter expanded prior - when variance components (VCV of the posterior distribution) are close to 0
#prior.exp <- list(R = list(V = 1, fix = 1), 
#                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
#Alternative parameter expanded prior
prior.Ville <- list(R = list(V = 1, nu = 1),
                    G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
#prior.exp.shinichi <- list(R = list(V = 1, nu = 0.002), 
#                           G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))


#Set up model parameters
# The formula
CS_poly <- worker_CoV ~ colony.size
# The MCMC parameters (iterations, thining, burnin)
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

################################################################
#Pairwise analyses among the predictor variables using MCMCglmm
################################################################

##########
#MF ~ CS#
##########
mulTree_data_CS_poly <- as.mulTree(data = data_CS_poly, tree = PT_data_CS_poly,
                                   taxa = "animal")

CS_poly <- worker_CoV ~ colony.size
mulTree(mulTree.data = mulTree_data_CS_poly, formula = CS_poly, priors = prior1,
        parameters = mul_parameters1M, output = "CS_poly_1M", ESS = 1000,
        chains = 2, family = "gaussian")



############ END ##########################


