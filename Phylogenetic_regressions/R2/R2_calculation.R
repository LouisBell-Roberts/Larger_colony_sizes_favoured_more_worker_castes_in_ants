########## R2 calculations for regression analysis using MCMCglmm ###########
####Across 400 trees
#07/03/2024
#Louis Bell-Roberts

.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))


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
#Write a function that reads in all of the individual models
model_reader_function <- function(pair, type) {
  
  # Create an empty list to store the loaded data
  mcmc_list <- list()
  
  # Check the value of the 'type' argument
  if (type == "Int") {
    # Get list of all rda files in the specified directory for Intercept type
    model_files <- list.files(path = file.path("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/MCMC_regressions/R2/Model_outputs/Intercept_only",pair,"1st_chain"), pattern = ".rds", full.names = TRUE)
  } else if (type == "Fix") {
    # Get list of all rda files in the specified directory for Fix type
    model_files <- list.files(path = file.path("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/MCMC_regressions/Model_outputs",pair,"1st_chain"), pattern = ".rds", full.names = TRUE)
  } else {
    # Return an error message if 'type' argument is not valid
    stop("Invalid value for 'type'. Please specify 'Int' or 'Fix'.")
  }
  
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
#Function for calculating R2_conditional along with associated CI for Gaussian response variables. Works on models with 1 random effect (=phylogeny)
MCMCglmm_R2_conditional_with_CI_gaussian <- function(mcmc_object) {
  mmF_list <- split(mcmc_object$Sol, seq(nrow(mcmc_object$Sol)))
  vmVarF <- purrr::map(mmF_list, ~.x %*% t(mcmc_object$X)) %>% map_dbl(~var(as.vector(.x)))
  
  # getting R2_conditional
  R2_conditional <- (vmVarF + mcmc_object$VCV[, 1]) / (vmVarF + mcmc_object$VCV[, 1] + mcmc_object$VCV[, 2])
  result <- list(
    mean = mean(R2_conditional),
    posterior_mode = posterior.mode(R2_conditional),
    credible_interval = HPDinterval(R2_conditional)
  )
  return(result)
}



###
#
#Function for calculating R2_marginal along with associated CI for Poisson response variables. Works on models with 1 random effect (=phylogeny)
MCMCglmm_R2_marginal_with_CI_poisson <- function(mmF, mm0) { #mmF is model including fixed effects, while mm0 is the intercept only model
  mVdis <- log(1 + 1/exp(mm0$Sol[,1] + 0.5*(rowSums(mmF$VCV))))
  
  mmF_list <- split(mmF$Sol, seq(nrow(mmF$Sol)))
  vmVarF1 <- map(mmF_list, ~.x %*% t(mmF$X)) %>% map_dbl(~var(as.vector(.x)))
  
  R2m <- vmVarF1/(vmVarF1 + mmF$VCV[,1] + mmF$VCV[,2] + mean(mVdis))
  result <- list(
    mean = mean(R2m),
    posterior_mode = posterior.mode(R2m),
    credible_interval = HPDinterval(R2m)
  )
  return(result)
}



###
#
#Function for calculating R2_conditional along with associated CI for Poisson response variables. Works on models with 1 random effect (=phylogeny)
MCMCglmm_R2_conditional_with_CI_poisson <- function(mmF, mm0) { #mmF is model including fixed effects, while mm0 is the intercept only model
  mVdis <- log(1 + 1/exp(mm0$Sol[,1] + 0.5*(rowSums(mmF$VCV))))
  
  mmF_list <- split(mmF$Sol, seq(nrow(mmF$Sol)))
  vmVarF1 <- map(mmF_list, ~.x %*% t(mmF$X)) %>% map_dbl(~var(as.vector(.x)))
  
  R2c <- (vmVarF1 + mmF$VCV[,1]) / (vmVarF1 + mmF$VCV[,1] + mmF$VCV[,2] + mean(mVdis)) #Should this be mean(mVdis)? - I think so
  result <- list(
    mean = mean(R2c),
    posterior_mode = posterior.mode(R2c),
    credible_interval = HPDinterval(R2c)
  )
  return(result)
}


###
#
#Function that converts marginal or conditional R2 list result to a dataframe
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
#Discrete count response variables
#########################################################

#
############
#MF_Caste
############

#Read in all of the individual models
mcmc_list_Int_MF_Caste <- model_reader_function(pair = "MF_Caste", type = "Int")
mcmc_list_Fix_MF_Caste <- model_reader_function(pair = "MF_Caste", type = "Fix")

####
#Marginal
#Calculate R2_marginal values across the 400 trees along with their associated CI
R2_marginal_single_MF_Caste <- MCMCglmm_R2_marginal_with_CI_poisson(mmF = mcmc_list_Fix_MF_Caste[[1]], mm0 = mcmc_list_Int_MF_Caste[[1]]) #0.03000627 #Apply to a single MCMCglmm model

R2_marginal_results_list_MF_Caste <- lapply(1:length(mcmc_list_Fix_MF_Caste), function(i) {
  MCMCglmm_R2_marginal_with_CI_poisson(mmF = mcmc_list_Fix_MF_Caste[[i]], mm0 = mcmc_list_Int_MF_Caste[[i]])
}) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_marginal_distribution_MF_Caste <- convert_results_list_to_dataframe(R2_marginal_results_list_MF_Caste)
quantile(R2_marginal_distribution_MF_Caste$post_mode_R2) #0.0013692827


####
#Conditional
#Calculate R2_conditional values across the 400 trees along with their associated CI
R2_conditional_single_MF_Caste <- MCMCglmm_R2_conditional_with_CI_poisson(mmF = mcmc_list_Fix_MF_Caste[[1]], mm0 = mcmc_list_Int_MF_Caste[[1]]) #Apply to a single MCMCglmm model
R2_conditional_results_list_MF_Caste <- lapply(1:length(mcmc_list_Fix_MF_Caste), function(i) {
  MCMCglmm_R2_conditional_with_CI_poisson(mmF = mcmc_list_Fix_MF_Caste[[i]], mm0 = mcmc_list_Int_MF_Caste[[i]])
}) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_conditional_distribution_MF_Caste <- convert_results_list_to_dataframe(R2_conditional_results_list_MF_Caste)
quantile(R2_conditional_distribution_MF_Caste$post_mode_R2) #0.05752129


#
############
#MFcat_Caste
############

#Read in all of the individual models
mcmc_list_Int_MFcat_Caste <- model_reader_function(pair = "MFcat_Caste", type = "Int")
mcmc_list_Fix_MFcat_Caste <- model_reader_function(pair = "MFcat_Caste", type = "Fix")

####
#Marginal
R2_marginal_results_list_MFcat_Caste <- lapply(1:length(mcmc_list_Fix_MFcat_Caste), function(i) {
  MCMCglmm_R2_marginal_with_CI_poisson(mmF = mcmc_list_Fix_MFcat_Caste[[i]], mm0 = mcmc_list_Int_MFcat_Caste[[i]])
})

#Create R2 marginal distribution data frame
R2_marginal_distribution_MFcat_Caste <- convert_results_list_to_dataframe(R2_marginal_results_list_MFcat_Caste)
quantile(R2_marginal_distribution_MFcat_Caste$post_mode_R2) #0.08466961


####
#Conditional
R2_conditional_results_list_MFcat_Caste <- lapply(1:length(mcmc_list_Fix_MFcat_Caste), function(i) {
  MCMCglmm_R2_conditional_with_CI_poisson(mmF = mcmc_list_Fix_MFcat_Caste[[i]], mm0 = mcmc_list_Int_MFcat_Caste[[i]])
})

#Create R2 marginal distribution data frame
R2_conditional_distribution_MFcat_Caste <- convert_results_list_to_dataframe(R2_conditional_results_list_MFcat_Caste)
quantile(R2_conditional_distribution_MFcat_Caste$post_mode_R2) #0.10912229


#
############
#CS_Caste
############

#Read in all of the individual models
mcmc_list_Int_CS_Caste <- model_reader_function(pair = "CS_Caste", type = "Int")
mcmc_list_Fix_CS_Caste <- model_reader_function(pair = "CS_Caste", type = "Fix")

####
#Marginal
#Calculate R2_marginal values across the 400 trees along with their associated CI
R2_marginal_single_CS_Caste <- MCMCglmm_R2_marginal_with_CI_poisson(mmF = mcmc_list_Fix_CS_Caste[[1]], mm0 = mcmc_list_Int_CS_Caste[[1]]) #Apply to a single MCMCglmm model

R2_marginal_results_list_CS_Caste <- lapply(1:length(mcmc_list_Fix_CS_Caste), function(i) {
  MCMCglmm_R2_marginal_with_CI_poisson(mmF = mcmc_list_Fix_CS_Caste[[i]], mm0 = mcmc_list_Int_CS_Caste[[i]])
}) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_marginal_distribution_CS_Caste <- convert_results_list_to_dataframe(R2_marginal_results_list_CS_Caste)
quantile(R2_marginal_distribution_CS_Caste$post_mode_R2) #0.05978971


####
#Conditional
#Calculate R2_conditional values across the 400 trees along with their associated CI
R2_conditional_single_CS_Caste <- MCMCglmm_R2_conditional_with_CI_poisson(mmF = mcmc_list_Fix_CS_Caste[[1]], mm0 = mcmc_list_Int_CS_Caste[[1]]) #Apply to a single MCMCglmm model
R2_conditional_results_list_CS_Caste <- lapply(1:length(mcmc_list_Fix_CS_Caste), function(i) {
  MCMCglmm_R2_conditional_with_CI_poisson(mmF = mcmc_list_Fix_CS_Caste[[i]], mm0 = mcmc_list_Int_CS_Caste[[i]])
}) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_conditional_distribution_CS_Caste <- convert_results_list_to_dataframe(R2_conditional_results_list_CS_Caste)
quantile(R2_conditional_distribution_CS_Caste$post_mode_R2) #0.06784183



#
############
#PG_Caste
############

#Read in all of the individual models
mcmc_list_Int_PG_Caste <- model_reader_function(pair = "PG_Caste", type = "Int")
mcmc_list_Fix_PG_Caste <- model_reader_function(pair = "PG_Caste", type = "Fix")

####
#Marginal
R2_marginal_results_list_PG_Caste <- lapply(1:length(mcmc_list_Fix_PG_Caste), function(i) {
  MCMCglmm_R2_marginal_with_CI_poisson(mmF = mcmc_list_Fix_PG_Caste[[i]], mm0 = mcmc_list_Int_PG_Caste[[i]])
})

#Create R2 marginal distribution data frame
R2_marginal_distribution_PG_Caste <- convert_results_list_to_dataframe(R2_marginal_results_list_PG_Caste)
quantile(R2_marginal_distribution_PG_Caste$post_mode_R2) #0.000190165


####
#Conditional
R2_conditional_results_list_PG_Caste <- lapply(1:length(mcmc_list_Fix_PG_Caste), function(i) {
  MCMCglmm_R2_conditional_with_CI_poisson(mmF = mcmc_list_Fix_PG_Caste[[i]], mm0 = mcmc_list_Int_PG_Caste[[i]])
})

#Create R2 marginal distribution data frame
R2_conditional_distribution_PG_Caste <- convert_results_list_to_dataframe(R2_conditional_results_list_PG_Caste)
quantile(R2_conditional_distribution_PG_Caste$post_mode_R2) #0.0036445321



#
############
#PGcat_Caste
############

#Read in all of the individual models
mcmc_list_Int_PGcat_Caste <- model_reader_function(pair = "PGcat_Caste", type = "Int")
mcmc_list_Fix_PGcat_Caste <- model_reader_function(pair = "PGcat_Caste", type = "Fix")

####
#Marginal
R2_marginal_results_list_PGcat_Caste <- lapply(1:length(mcmc_list_Fix_PGcat_Caste), function(i) {
  MCMCglmm_R2_marginal_with_CI_poisson(mmF = mcmc_list_Fix_PGcat_Caste[[i]], mm0 = mcmc_list_Int_PGcat_Caste[[i]])
})

#Create R2 marginal distribution data frame
R2_marginal_distribution_PGcat_Caste <- convert_results_list_to_dataframe(R2_marginal_results_list_PGcat_Caste)
quantile(R2_marginal_distribution_PGcat_Caste$post_mode_R2) #0.0017810772


####
#Conditional
R2_conditional_results_list_PGcat_Caste <- lapply(1:length(mcmc_list_Fix_PGcat_Caste), function(i) {
  MCMCglmm_R2_conditional_with_CI_poisson(mmF = mcmc_list_Fix_PGcat_Caste[[i]], mm0 = mcmc_list_Int_PGcat_Caste[[i]])
})

#Create R2 marginal distribution data frame
R2_conditional_distribution_PGcat_Caste <- convert_results_list_to_dataframe(R2_conditional_results_list_PGcat_Caste)
quantile(R2_conditional_distribution_PGcat_Caste$post_mode_R2) #0.020463815



#########################################################
#Continuous response variables
#########################################################

#
############
#MF_siz_var
############

#Read in all of the individual models
mcmc_list_MF_siz_var <- model_reader_function(pair = "MF_siz_var", type = "Fix")

####
#Marginal
R2_marginal_results_list_MF_siz_var <- lapply(mcmc_list_MF_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MF_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_MF_siz_var)
quantile(R2_marginal_distribution_MF_siz_var$post_mode_R2) #0.0013689590

####
#Conditional
R2_conditional_results_list_MF_siz_var <- lapply(mcmc_list_MF_siz_var, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_MF_siz_var <- convert_results_list_to_dataframe(R2_conditional_results_list_MF_siz_var)
quantile(R2_conditional_distribution_MF_siz_var$post_mode_R2) #0.8171965



############
#MFcat_siz_var
############

#Read in all of the individual models
mcmc_list_MFcat_siz_var <- model_reader_function(pair = "MFcat_siz_var", type = "Fix")

####
#Marginal
R2_marginal_results_list_MFcat_siz_var <- lapply(mcmc_list_MFcat_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MFcat_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_MFcat_siz_var)
quantile(R2_marginal_distribution_MFcat_siz_var$post_mode_R2) #0.10363226

####
#Conditional
R2_conditional_results_list_MFcat_siz_var <- lapply(mcmc_list_MFcat_siz_var, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_MFcat_siz_var <- convert_results_list_to_dataframe(R2_conditional_results_list_MFcat_siz_var)
quantile(R2_conditional_distribution_MFcat_siz_var$post_mode_R2) #0.8490185


#
############
#CS_siz_var
############

#Read in all of the individual models
mcmc_list_CS_siz_var <- model_reader_function(pair = "CS_siz_var", type = "Fix")

####
#Marginal
#Calculate R2_marginal values across the 400 trees along with their associated CI
R2_marginal_single_CS_siz_var <- MCMCglmm_R2_marginal_with_CI_gaussian(mcmc_list_CS_siz_var[[1]]) #Apply to a single MCMCglmm model
R2_marginal_results_list_CS_siz_var <- lapply(mcmc_list_CS_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_marginal_distribution_CS_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_CS_siz_var)
quantile(R2_marginal_distribution_CS_siz_var$post_mode_R2) #0.14194722

####
#Conditional
#Calculate R2_conditional values across the 400 trees along with their associated CI
R2_conditional_single_CS_siz_var <- MCMCglmm_R2_conditional_with_CI_gaussian(mcmc_list_CS_siz_var[[1]]) #Apply to a single MCMCglmm model
R2_conditional_results_list_CS_siz_var <- lapply(mcmc_list_CS_siz_var, MCMCglmm_R2_conditional_with_CI_gaussian) #Apply across a list of MCMCglmm models

#Create R2 marginal distribution data frame
R2_conditional_distribution_CS_siz_var <- convert_results_list_to_dataframe(R2_conditional_results_list_CS_siz_var)
quantile(R2_conditional_distribution_CS_siz_var$post_mode_R2) #0.7240469



#
############
#PG_siz_var
############

#Read in all of the individual models
mcmc_list_PG_siz_var <- model_reader_function(pair = "PG_siz_var", type = "Fix")

####
#Marginal
R2_marginal_results_list_PG_siz_var <- lapply(mcmc_list_PG_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_PG_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_PG_siz_var)
quantile(R2_marginal_distribution_PG_siz_var$post_mode_R2) #0.0001601481

####
#Conditional
R2_conditional_results_list_PG_siz_var <- lapply(mcmc_list_PG_siz_var, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_PG_siz_var <- convert_results_list_to_dataframe(R2_conditional_results_list_PG_siz_var)
quantile(R2_conditional_distribution_PG_siz_var$post_mode_R2) #0.7258377


############
#PGcat_siz_var
############

#Read in all of the individual models
mcmc_list_PGcat_siz_var <- model_reader_function(pair = "PGcat_siz_var", type = "Fix")

####
#Marginal
R2_marginal_results_list_PGcat_siz_var <- lapply(mcmc_list_PGcat_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_PGcat_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_PGcat_siz_var)
quantile(R2_marginal_distribution_PGcat_siz_var$post_mode_R2) #0.0017369324

####
#Conditional
R2_conditional_results_list_PGcat_siz_var <- lapply(mcmc_list_PGcat_siz_var, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_PGcat_siz_var <- convert_results_list_to_dataframe(R2_conditional_results_list_PGcat_siz_var)
quantile(R2_conditional_distribution_PGcat_siz_var$post_mode_R2) #0.7279908


#
############
#MF_mono_siz_var
############

#Read in all of the individual models
mcmc_list_MF_mono_siz_var <- model_reader_function(pair = "MF_mono_siz_var", type = "Fix")

####
#Marginal
R2_marginal_results_list_MF_mono_siz_var <- lapply(mcmc_list_MF_mono_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MF_mono_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_MF_mono_siz_var)
quantile(R2_marginal_distribution_MF_mono_siz_var$post_mode_R2) #0.0003314781

####
#Conditional
R2_conditional_results_list_MF_mono_siz_var <- lapply(mcmc_list_MF_mono_siz_var, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_MF_mono_siz_var <- convert_results_list_to_dataframe(R2_conditional_results_list_MF_mono_siz_var)
quantile(R2_conditional_distribution_MF_mono_siz_var$post_mode_R2) #0.08714506



############
#MFcat_mono_siz_var
############

#Read in all of the individual models
mcmc_list_MFcat_mono_siz_var <- model_reader_function(pair = "MFcat_mono_siz_var", type = "Fix")

####
#Marginal
R2_marginal_results_list_MFcat_mono_siz_var <- lapply(mcmc_list_MFcat_mono_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MFcat_mono_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_MFcat_mono_siz_var)
quantile(R2_marginal_distribution_MFcat_mono_siz_var$post_mode_R2) #0.040838872

####
#Conditional
R2_conditional_results_list_MFcat_mono_siz_var <- lapply(mcmc_list_MFcat_mono_siz_var, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_MFcat_mono_siz_var <- convert_results_list_to_dataframe(R2_conditional_results_list_MFcat_mono_siz_var)
quantile(R2_conditional_distribution_MFcat_mono_siz_var$post_mode_R2) #0.16868801


#
############
#CS_mono_siz_var
############

#Read in all of the individual models
mcmc_list_CS_mono_siz_var <- model_reader_function(pair = "CS_mono_siz_var", type = "Fix")

####
#Marginal
#Calculate R2_marginal values across the 400 trees along with their associated CI
R2_marginal_results_list_CS_mono_siz_var <- lapply(mcmc_list_CS_mono_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_CS_mono_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_CS_mono_siz_var)
quantile(R2_marginal_distribution_CS_mono_siz_var$post_mode_R2) #0.07591466

####
#Conditional
R2_conditional_results_list_CS_mono_siz_var <- lapply(mcmc_list_CS_mono_siz_var, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_CS_mono_siz_var <- convert_results_list_to_dataframe(R2_conditional_results_list_CS_mono_siz_var)
quantile(R2_conditional_distribution_CS_mono_siz_var$post_mode_R2) #0.1785738



#
############
#PG_mono_siz_var
############

#Read in all of the individual models
mcmc_list_PG_mono_siz_var <- model_reader_function(pair = "PG_mono_siz_var", type = "Fix")

####
#Marginal
R2_marginal_results_list_PG_mono_siz_var <- lapply(mcmc_list_PG_mono_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_PG_mono_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_PG_mono_siz_var)
quantile(R2_marginal_distribution_PG_mono_siz_var$post_mode_R2) #0.0003373458

####
#Conditional
R2_conditional_results_list_PG_mono_siz_var <- lapply(mcmc_list_PG_mono_siz_var, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_PG_mono_siz_var <- convert_results_list_to_dataframe(R2_conditional_results_list_PG_mono_siz_var)
quantile(R2_conditional_distribution_PG_mono_siz_var$post_mode_R2) #0.09160169


############
#PGcat_mono_siz_var
############

#Read in all of the individual models
mcmc_list_PGcat_mono_siz_var <- model_reader_function(pair = "PGcat_mono_siz_var", type = "Fix")

####
#Marginal
R2_marginal_results_list_PGcat_mono_siz_var <- lapply(mcmc_list_PGcat_mono_siz_var, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_PGcat_mono_siz_var <- convert_results_list_to_dataframe(R2_marginal_results_list_PGcat_mono_siz_var)
quantile(R2_marginal_distribution_PGcat_mono_siz_var$post_mode_R2) #0.0038788658

####
#Conditional
R2_conditional_results_list_PGcat_mono_siz_var <- lapply(mcmc_list_PGcat_mono_siz_var, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_PGcat_mono_siz_var <- convert_results_list_to_dataframe(R2_conditional_results_list_PGcat_mono_siz_var)
quantile(R2_conditional_distribution_PGcat_mono_siz_var$post_mode_R2) #0.1216039



#
############
#MF_PG
############

#Read in all of the individual models
mcmc_list_MF_PG <- model_reader_function(pair = "MF_PG", type = "Fix")

####
#Marginal
R2_marginal_results_list_MF_PG <- lapply(mcmc_list_MF_PG, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MF_PG <- convert_results_list_to_dataframe(R2_marginal_results_list_MF_PG)
quantile(R2_marginal_distribution_MF_PG$post_mode_R2) #0.0003545792

####
#Conditional
R2_conditional_results_list_MF_PG <- lapply(mcmc_list_MF_PG, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_MF_PG <- convert_results_list_to_dataframe(R2_conditional_results_list_MF_PG)
quantile(R2_conditional_distribution_MF_PG$post_mode_R2) #0.8484439


#
############
#MFcat_PG
############

#Read in all of the individual models
mcmc_list_MFcat_PG <- model_reader_function(pair = "MFcat_PG", type = "Fix")

####
#Marginal
R2_marginal_results_list_MFcat_PG <- lapply(mcmc_list_MFcat_PG, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MFcat_PG <- convert_results_list_to_dataframe(R2_marginal_results_list_MFcat_PG)
quantile(R2_marginal_distribution_MFcat_PG$post_mode_R2) #0.0426302698

####
#Conditional
R2_conditional_results_list_MFcat_PG <- lapply(mcmc_list_MFcat_PG, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_MFcat_PG <- convert_results_list_to_dataframe(R2_conditional_results_list_MFcat_PG)
quantile(R2_conditional_distribution_MFcat_PG$post_mode_R2) #0.5576917


#
############
#MF_CS
############

#Read in all of the individual models
mcmc_list_MF_CS <- model_reader_function(pair = "MF_CS", type = "Fix")

####
#Marginal
R2_marginal_results_list_MF_CS <- lapply(mcmc_list_MF_CS, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MF_CS <- convert_results_list_to_dataframe(R2_marginal_results_list_MF_CS)
quantile(R2_marginal_distribution_MF_CS$post_mode_R2) #0.076688517

####
#Conditional
R2_conditional_results_list_MF_CS <- lapply(mcmc_list_MF_CS, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_MF_CS <- convert_results_list_to_dataframe(R2_conditional_results_list_MF_CS)
quantile(R2_conditional_distribution_MF_CS$post_mode_R2) #0.9110338


#
############
#MFcat_CS
############

#Read in all of the individual models
mcmc_list_MFcat_CS <- model_reader_function(pair = "MFcat_CS", type = "Fix")

####
#Marginal
R2_marginal_results_list_MFcat_CS <- lapply(mcmc_list_MFcat_CS, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_MFcat_CS <- convert_results_list_to_dataframe(R2_marginal_results_list_MFcat_CS)
quantile(R2_marginal_distribution_MFcat_CS$post_mode_R2) #0.11838140

####
#Conditional
R2_conditional_results_list_MFcat_CS <- lapply(mcmc_list_MFcat_CS, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_MFcat_CS <- convert_results_list_to_dataframe(R2_conditional_results_list_MFcat_CS)
quantile(R2_conditional_distribution_MFcat_CS$post_mode_R2) #0.9069618


#
############
#CS_PG
############

#Read in all of the individual models
mcmc_list_CS_PG <- model_reader_function(pair = "CS_PG", type = "Fix")

####
#Marginal
R2_marginal_results_list_CS_PG <- lapply(mcmc_list_CS_PG, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_CS_PG <- convert_results_list_to_dataframe(R2_marginal_results_list_CS_PG)
quantile(R2_marginal_distribution_CS_PG$post_mode_R2) #0.00005126519

####
#Conditional
R2_conditional_results_list_CS_PG <- lapply(mcmc_list_CS_PG, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_CS_PG <- convert_results_list_to_dataframe(R2_conditional_results_list_CS_PG)
quantile(R2_conditional_distribution_CS_PG$post_mode_R2) #0.8179485


#
############
#CS_PGcat
############

#Read in all of the individual models
mcmc_list_CS_PGcat <- model_reader_function(pair = "CS_PGcat", type = "Fix")

####
#Marginal
R2_marginal_results_list_CS_PGcat <- lapply(mcmc_list_CS_PGcat, MCMCglmm_R2_marginal_with_CI_gaussian)
R2_marginal_distribution_CS_PGcat <- convert_results_list_to_dataframe(R2_marginal_results_list_CS_PGcat)
quantile(R2_marginal_distribution_CS_PGcat$post_mode_R2) #0.0004636427

####
#Conditional
R2_conditional_results_list_CS_PGcat <- lapply(mcmc_list_CS_PGcat, MCMCglmm_R2_conditional_with_CI_gaussian)
R2_conditional_distribution_CS_PGcat <- convert_results_list_to_dataframe(R2_conditional_results_list_CS_PGcat)
quantile(R2_conditional_distribution_CS_PGcat$post_mode_R2) #0.8121099



########################################################################

#Create table of the R2 results
R2_marginal_vec <- c(median(R2_marginal_distribution_MF_Caste$post_mode_R2),
                     median(R2_marginal_distribution_MFcat_Caste$post_mode_R2),
                     median(R2_marginal_distribution_CS_Caste$post_mode_R2),
                     median(R2_marginal_distribution_PG_Caste$post_mode_R2),
                     median(R2_marginal_distribution_PGcat_Caste$post_mode_R2),
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

R2_conditional_vec <- c(median(R2_conditional_distribution_MF_Caste$post_mode_R2),
                        median(R2_conditional_distribution_MFcat_Caste$post_mode_R2),
                        median(R2_conditional_distribution_CS_Caste$post_mode_R2),
                        median(R2_conditional_distribution_PG_Caste$post_mode_R2),
                        median(R2_conditional_distribution_PGcat_Caste$post_mode_R2),
                        median(R2_conditional_distribution_MF_siz_var$post_mode_R2),
                        median(R2_conditional_distribution_MFcat_siz_var$post_mode_R2),
                        median(R2_conditional_distribution_CS_siz_var$post_mode_R2),
                        median(R2_conditional_distribution_PG_siz_var$post_mode_R2),
                        median(R2_conditional_distribution_PGcat_siz_var$post_mode_R2),
                        median(R2_conditional_distribution_MF_mono_siz_var$post_mode_R2),
                        median(R2_conditional_distribution_MFcat_mono_siz_var$post_mode_R2),
                        median(R2_conditional_distribution_CS_mono_siz_var$post_mode_R2),
                        median(R2_conditional_distribution_PG_mono_siz_var$post_mode_R2),
                        median(R2_conditional_distribution_PGcat_mono_siz_var$post_mode_R2),
                        median(R2_conditional_distribution_MF_CS$post_mode_R2),
                        median(R2_conditional_distribution_MFcat_CS$post_mode_R2),
                        median(R2_conditional_distribution_MF_PG$post_mode_R2),
                        median(R2_conditional_distribution_MFcat_PG$post_mode_R2),
                        median(R2_conditional_distribution_CS_PG$post_mode_R2),
                        median(R2_conditional_distribution_CS_PGcat$post_mode_R2))

analysis_vec <- c("Queen mating frequency - Caste", "Queen mating frequency (categorical) - Caste", "Colony size - Caste", "Number of queens - Caste", "Number of queens (categorical) - Caste", "Queen mating frequency - Variation in worker size", "Queen mating frequency (categorical) - Variation in worker size", "Colony size - Variation in worker size", "Number of queens - Variation in worker size", "Number of queens (categorical) - Variation in worker size", "Queen mating frequency - Variation in worker size (Monomorphic)", "Queen mating frequency (categorical) - Variation in worker size (Monomorphic)", "Colony size - Variation in worker size (Monomorphic)", "Number of queens - Variation in worker size (Monomorphic)", "Number of queens (categorical) - Variation in worker size (Monomorphic)", "Queen mating frequency - Colony size", "Queen mating frequency (categorical) - Colony size", "Queen mating frequency - Number of queens", "Queen mating frequency (categorical) - Number of queens", "Colony size - Number of queens", "Colony size - Number of queens (categorical)")

# Create dataframe
R2_df <- data.frame(Analysis = analysis_vec, R2_marginal = R2_marginal_vec, R2_conditional = R2_conditional_vec)

R2_df$R2_marginal <- format(R2_df$R2_marginal, scientific = FALSE)
R2_df$R2_conditional <- format(R2_df$R2_conditional, scientific = FALSE)
write.csv(R2_df, file = "/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/MCMC_regressions/R2/R2_result.csv", row.names = F)








