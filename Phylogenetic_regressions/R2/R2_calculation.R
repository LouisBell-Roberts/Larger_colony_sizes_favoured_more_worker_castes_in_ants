########## R2 calculations for regression analysis using MCMCglmm ###########
####Across 400 trees
#07/03/2024
#Louis Bell-Roberts



#Loading essential packages ####
library(data.table)
library(ape)
library(phytools)
library(MCMCglmm)
library(tidyverse)
library(gtools)
library(arm)
library(purrr)
library(gtools)
#

#########################################################

##
#####Custom functions
##


###
#
#Function that reads in all of the individual models for a given pair of variables analysed. Models for each pair of variables to be read in from their own separate folders.
model_reader_function <- function(pair) {
  
  # Create an empty list to store the loaded data
  mcmc_list <- list()
  
  # Get list of all rda files in the specified directory
  model_files <- list.files(path = file.path(pair), pattern = ".rds", full.names = TRUE)
  
  ## Sort the file names alphanumerically
  sorted_files <- mixedsort(model_files)
  
  ## Read in all models
  # Read in each MCMCglmm model
  mcmc_list <- lapply(sorted_files, function(file) {
    readRDS(file)
  }) # mcmc_list is the list of all of the models read in
  
  return(mcmc_list)
}




###
#
#Function for calculating R2_marginal along with associated CI for Gaussian response variables. Works on models with 1 random effect (=phylogeny)
MCMCglmm_R2_marginal_with_CI_gaussian <- function(mcmc_object) {
  mmF_list <- split(mcmc_object$Sol, seq(nrow(mcmc_object$Sol)))
  vmVarF <- purrr::map(mmF_list, ~.x %*% t(mcmc_object$X)) %>% map_dbl(~var(as.vector(.x)))
  
  # getting R2_marginal
  R2_marginal <- vmVarF / (vmVarF + mcmc_object$VCV[, 1] + mcmc_object$VCV[, 2])
  result <- list(
    mean = mean(R2_marginal),
    posterior_mode = posterior.mode(R2_marginal),
    credible_interval = HPDinterval(R2_marginal)
  )
  return(result)
}




###
#
#Function that converts marginal R2 list result to a dataframe
convert_results_list_to_dataframe <- function(results_list) {
  # Use lapply to extract values from the results_list
  mean_R2 <- unlist(lapply(results_list, function(x) x$mean))
  post_mode_R2 <- unlist(lapply(results_list, function(x) x$posterior_mode))
  lower_CI <- unlist(lapply(results_list, function(x) x$credible_interval[1]))
  upper_CI <- unlist(lapply(results_list, function(x) x$credible_interval[2]))
  
  # Create a data frame from the extracted values
  R2_distribution <- data.frame(mean_R2, post_mode_R2, lower_CI, upper_CI)
  
  return(R2_distribution)
}

#########################################################




#########################################################
#Continuous response variables
#########################################################

#
############
#MF_siz_var
############

#Read in all of the individual models
mcmc_list_MF_siz_var <- model_reader_function(pair = "MF_siz_var")

####
#Marginal
#Calculate R2_marginal values across the 400 trees along with their associated CI
R2_marginal_results_list_MF_siz_var <- lapply(mcmc_list_MF_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
#Create R2 marginal distribution data frame
R2_marginal_distribution_MF_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_MF_siz_var)
quantile(R2_marginal_distribution_MF_siz_var$post_mode_R2)


############
#MFcat_siz_var
############

#Read in all of the individual models
mcmc_list_MFcat_siz_var <- model_reader_function(pair = "MFcat_siz_var")

####
#Marginal
R2_marginal_results_list_MFcat_siz_var <- lapply(mcmc_list_MFcat_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MFcat_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_MFcat_siz_var)
quantile(R2_marginal_distribution_MFcat_siz_var$post_mode_R2)


#
############
#CS_siz_var
############

#Read in all of the individual models
mcmc_list_CS_siz_var <- model_reader_function(pair = "CS_siz_var")

####
#Marginal
R2_marginal_results_list_CS_siz_var <- lapply(mcmc_list_CS_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_CS_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_CS_siz_var)
quantile(R2_marginal_distribution_CS_siz_var$post_mode_R2)


#
############
#PG_siz_var
############

#Read in all of the individual models
mcmc_list_PG_siz_var <- model_reader_function(pair = "PG_siz_var")

####
#Marginal
R2_marginal_results_list_PG_siz_var <- lapply(mcmc_list_PG_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_PG_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_PG_siz_var)
quantile(R2_marginal_distribution_PG_siz_var$post_mode_R2)


############
#PGcat_siz_var
############

#Read in all of the individual models
mcmc_list_PGcat_siz_var <- model_reader_function(pair = "PGcat_siz_var")

####
#Marginal
R2_marginal_results_list_PGcat_siz_var <- lapply(mcmc_list_PGcat_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_PGcat_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_PGcat_siz_var)
quantile(R2_marginal_distribution_PGcat_siz_var$post_mode_R2)


#
############
#MF_mono_siz_var
############

#Read in all of the individual models
mcmc_list_MF_mono_siz_var <- model_reader_function(pair = "MF_mono_siz_var")

####
#Marginal
R2_marginal_results_list_MF_mono_siz_var <- lapply(mcmc_list_MF_mono_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MF_mono_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_MF_mono_siz_var)
quantile(R2_marginal_distribution_MF_mono_siz_var$post_mode_R2)



############
#MFcat_mono_siz_var
############

#Read in all of the individual models
mcmc_list_MFcat_mono_siz_var <- model_reader_function(pair = "MFcat_mono_siz_var")

####
#Marginal
R2_marginal_results_list_MFcat_mono_siz_var <- lapply(mcmc_list_MFcat_mono_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MFcat_mono_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_MFcat_mono_siz_var)
quantile(R2_marginal_distribution_MFcat_mono_siz_var$post_mode_R2)


#
############
#CS_mono_siz_var
############

#Read in all of the individual models
mcmc_list_CS_mono_siz_var <- model_reader_function(pair = "CS_mono_siz_var")

####
#Marginal
R2_marginal_results_list_CS_mono_siz_var <- lapply(mcmc_list_CS_mono_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_CS_mono_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_CS_mono_siz_var)
quantile(R2_marginal_distribution_CS_mono_siz_var$post_mode_R2)


#
############
#PG_mono_siz_var
############

#Read in all of the individual models
mcmc_list_PG_mono_siz_var <- model_reader_function(pair = "PG_mono_siz_var")

####
#Marginal
R2_marginal_results_list_PG_mono_siz_var <- lapply(mcmc_list_PG_mono_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_PG_mono_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_PG_mono_siz_var)
quantile(R2_marginal_distribution_PG_mono_siz_var$post_mode_R2)


############
#PGcat_mono_siz_var
############

#Read in all of the individual models
mcmc_list_PGcat_mono_siz_var <- model_reader_function(pair = "PGcat_mono_siz_var")

####
#Marginal
R2_marginal_results_list_PGcat_mono_siz_var <- lapply(mcmc_list_PGcat_mono_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_PGcat_mono_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_PGcat_mono_siz_var)
quantile(R2_marginal_distribution_PGcat_mono_siz_var$post_mode_R2)


#
############
#MF_PG
############

#Read in all of the individual models
mcmc_list_MF_PG <- model_reader_function(pair = "MF_PG")

####
#Marginal
R2_marginal_results_list_MF_PG <- lapply(mcmc_list_MF_PG, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MF_PG <- convert_results_list_to_dataframe(R2_marginal_results_list_MF_PG)
quantile(R2_marginal_distribution_MF_PG$post_mode_R2)


#
############
#MFcat_PG
############

#Read in all of the individual models
mcmc_list_MFcat_PG <- model_reader_function(pair = "MFcat_PG")

####
#Marginal
R2_marginal_results_list_MFcat_PG <- lapply(mcmc_list_MFcat_PG, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MFcat_PG <- convert_results_list_to_dataframe(R2_marginal_results_list_MFcat_PG)
quantile(R2_marginal_distribution_MFcat_PG$post_mode_R2)


#
############
#MF_CS
############

#Read in all of the individual models
mcmc_list_MF_CS <- model_reader_function(pair = "MF_CS")

####
#Marginal
R2_marginal_results_list_MF_CS <- lapply(mcmc_list_MF_CS, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MF_CS <- convert_results_list_to_dataframe(R2_marginal_results_list_MF_CS)
quantile(R2_marginal_distribution_MF_CS$post_mode_R2)


#
############
#MFcat_CS
############

#Read in all of the individual models
mcmc_list_MFcat_CS <- model_reader_function(pair = "MFcat_CS")

####
#Marginal
R2_marginal_results_list_MFcat_CS <- lapply(mcmc_list_MFcat_CS, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MFcat_CS <- convert_results_list_to_dataframe(R2_marginal_results_list_MFcat_CS)
quantile(R2_marginal_distribution_MFcat_CS$post_mode_R2)


#
############
#CS_PG
############

#Read in all of the individual models
mcmc_list_CS_PG <- model_reader_function(pair = "CS_PG")

####
#Marginal
R2_marginal_results_list_CS_PG <- lapply(mcmc_list_CS_PG, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_CS_PG <- convert_results_list_to_dataframe(R2_marginal_results_list_CS_PG)
quantile(R2_marginal_distribution_CS_PG$post_mode_R2)


#
############
#CS_PGcat
############

#Read in all of the individual models
mcmc_list_CS_PGcat <- model_reader_function(pair = "CS_PGcat")

####
#Marginal
R2_marginal_results_list_CS_PGcat <- lapply(mcmc_list_CS_PGcat, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_CS_PGcat <- convert_results_list_to_dataframe(R2_marginal_results_list_CS_PGcat)
quantile(R2_marginal_distribution_CS_PGcat$post_mode_R2)


########################################################################

#Create table of the R2 results
R2_marginal_vec <- c(
                     median(R2_marginal_distribution_MF_siz_var$post_mode_R2),
                     median(R2_marginal_distribution_MFcat_siz_var$post_mode_R2),
                     median(R2_marginal_distribution_CS_siz_var$post_mode_R2),
                     median(R2_marginal_distribution_PG_siz_var$post_mode_R2),
                     median(R2_marginal_distribution_PGcat_siz_var$post_mode_R2),
                     median(R2_marginal_distribution_MF_mono_siz_var$post_mode_R2),
                     median(R2_marginal_distribution_MFcat_mono_siz_var$post_mode_R2),
                     median(R2_marginal_distribution_CS_mono_siz_var$post_mode_R2),
                     median(R2_marginal_distribution_PG_mono_siz_var$post_mode_R2),
                     median(R2_marginal_distribution_PGcat_mono_siz_var$post_mode_R2),
                     median(R2_marginal_distribution_MF_CS$post_mode_R2),
                     median(R2_marginal_distribution_MFcat_CS$post_mode_R2),
                     median(R2_marginal_distribution_MF_PG$post_mode_R2),
                     median(R2_marginal_distribution_MFcat_PG$post_mode_R2),
                     median(R2_marginal_distribution_CS_PG$post_mode_R2),
                     median(R2_marginal_distribution_CS_PGcat$post_mode_R2))


analysis_vec <- c("Queen mating frequency - Variation in worker size", "Queen mating frequency (categorical) - Variation in worker size", "Colony size - Variation in worker size", "Number of queens - Variation in worker size", "Number of queens (categorical) - Variation in worker size", "Queen mating frequency - Variation in worker size (Monomorphic)", "Queen mating frequency (categorical) - Variation in worker size (Monomorphic)", "Colony size - Variation in worker size (Monomorphic)", "Number of queens - Variation in worker size (Monomorphic)", "Number of queens (categorical) - Variation in worker size (Monomorphic)", "Queen mating frequency - Colony size", "Queen mating frequency (categorical) - Colony size", "Queen mating frequency - Number of queens", "Queen mating frequency (categorical) - Number of queens", "Colony size - Number of queens", "Colony size - Number of queens (categorical)")

# Create dataframe
R2_df <- data.frame(Analysis = analysis_vec, R2_marginal = R2_marginal_vec)

R2_df$R2_marginal <- format(R2_df$R2_marginal, scientific = FALSE)
# write.csv(R2_df, file = "R2_result.csv", row.names = F)








