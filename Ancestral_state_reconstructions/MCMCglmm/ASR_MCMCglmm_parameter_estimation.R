############################### Parameter estimation for MCMCglmm ASR analyses ########################################
###To quantifying phylogenetic uncertainty in parameter estimates, analyses run over a sample of 400 phylogenetic trees
##Louis Bell-Roberts
#01/02/2024

library(MCMCglmm)
library(coda)
library(gtools)

##########
#MF_Caste#
##########

##Multiple models
# setwd("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/ASR/Model_outputs/MCMCglmm/MF_Caste/1st_run/") 
setwd("/Volumes/ADATA SE800/DOL_worker_castes/ASR/MCMCglmm/MF_Caste/1st_run/") 

# get list of all RDS files in the working directory and sort alphanumerically
model_files <- mixedsort(list.files(pattern = "\\.rds$"))

# read in all models using lapply()
mcmc_list <- lapply(model_files, readRDS)

#Combine all 400 MCMCglmm Sol model outputs
# extract the Sol objects from each model output in the list
chain_list <- lapply(mcmc_list, function(x) x$Sol)

#Check the number of columns in each model output
num_cols <- lapply(chain_list, ncol)
num_cols==210

#Select only the columns of interest before combining posterior distributions
##Required because the "polymorphic.both" transition is not in every model output. Therefore, the number of columns is not the same in all model outputs
###Columns of interest: "monomorphic.both" and "monomorphic.onlymonomorphic"
chain_list_sub <- lapply(chain_list, function(x) x[, 1:2])

# combine the chains using rbind
MCMC_combined <- do.call(rbind, chain_list_sub)
MCMC_combined <- as.mcmc(MCMC_combined)

#Obtain point estimates - 10^ transformation used as log10 transformation was used prior to the analysis
10^(posterior.mode(MCMC_combined[,1:2])) #CAT2monomorphic.both = 2.462114; CAT2monomorphic.only monomorphic = 1.994259
10^(HPDinterval(MCMC_combined[,1:2])) #CAT2monomorphic.both = 1.1040772, 5.551794; CAT2monomorphic.only monomorphic = 0.9450067, 3.904074


# Compare monomorphic.both (GAIN OF POLYMORPHISM) with monomorphic.only monomorphic (NO CHANGE)
# to test if high colony size makes the origin of worker polymorphism more likely
table(MCMC_combined[,1] > MCMC_combined[,2]) / length(MCMC_combined[,2])	#TRUE = 0.841985. Therefore, Queen mating frequency is not significantly higher in ancestors where transitions to multiple castes occurred



################
#CS_Caste#
################

#Load in models
# setwd("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/ASR/Model_outputs/MCMCglmm/CS_Caste/1st_run/")
setwd("/Volumes/ADATA SE800/DOL_worker_castes/ASR/MCMCglmm/CS_Caste/1st_run/") 

# get list of all RDS files in the working directory and sort alphanumerically
model_files <- mixedsort(list.files(pattern = "\\.rds$"))
# read in all models using lapply()
mcmc_list <- lapply(model_files, readRDS)

#Combine all 400 MCMCglmm Sol model outputs
# extract the Sol objects from each model output in the list
chain_list <- lapply(mcmc_list, function(x) x$Sol)

#Check the number of columns in each model output
num_cols <- lapply(chain_list, ncol)
num_cols==874

#Select only the columns of interest before combining posterior distributions
##Required because the "polymorphic.both" transition is not in every model output. Therefore, the number of columns is not the same in all model outputs
###Columns of interest: "monomorphic.both" and "monomorphic.onlymonomorphic"
chain_list_sub <- lapply(chain_list, function(x) x[, 1:4])

# combine the chains using rbind
MCMC_combined <- do.call(rbind, chain_list_sub)
MCMC_combined <- as.mcmc(MCMC_combined)

#Obtain point estimates - 10^ transformation used as log10 transformation was used prior to the analysis
10^(posterior.mode(MCMC_combined[,1:4])) #CAT2monomorphic.both = 1025.0506; CAT2monomorphic.only monomorphic = 281.6716
10^(HPDinterval(MCMC_combined[,1:4])) #CAT2monomorphic.only monomorphic = 40.03965, 1590.898; CAT2monomorphic.both = 137.16739, 9064.821

# Compare monomorphic.both (GAIN OF POLYMORPHISM) with monomorphic.only monomorphic (NO CHANGE)
# to test if high colony size makes the origin of worker polymorphism more likely
table(MCMC_combined[,1] > MCMC_combined[,2]) / length(MCMC_combined[,2])	# TRUE = 0.9966. Therefore, colony size is significantly higher in ancestors where transitions to multiple castes occurred

