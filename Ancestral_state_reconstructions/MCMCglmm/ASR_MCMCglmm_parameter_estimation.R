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
#Assign file path for where MF_Caste MCMCglmm ASR models are saved
setwd("")

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
10^(posterior.mode(MCMC_combined[,1:2]))
10^(HPDinterval(MCMC_combined[,1:2]))


# Compare monomorphic.both (GAIN OF POLYMORPHISM) with monomorphic.only monomorphic (NO CHANGE)
# to test if high colony size makes the origin of worker polymorphism more likely
table(MCMC_combined[,1] > MCMC_combined[,2]) / length(MCMC_combined[,2])

################
#CS_Caste#
################

##Load in models
#Assign file path for where CS_Caste MCMCglmm ASR models are saved
setwd("") 

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
10^(posterior.mode(MCMC_combined[,1:4]))
10^(HPDinterval(MCMC_combined[,1:4]))

# Compare monomorphic.both (GAIN OF POLYMORPHISM) with monomorphic.only monomorphic (NO CHANGE)
# to test if high colony size makes the origin of worker polymorphism more likely
table(MCMC_combined[,1] > MCMC_combined[,2]) / length(MCMC_combined[,2])
