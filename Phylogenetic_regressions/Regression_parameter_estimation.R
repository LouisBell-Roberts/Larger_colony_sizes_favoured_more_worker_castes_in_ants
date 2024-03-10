############################### Parameter estimation for BPMM analyses ########################################
###To quantifying phylogenetic uncertainty in parameter estimates, analyses run over a sample of 400 phylogenetic trees
##Louis Bell-Roberts
#01/02/2024


library(tidyverse)
library(ape)
library(phytools)
library(ggplot2)
library(MCMCglmm)
library(arm)
library(lme4)
library(coda)


#########
#Customs functions
#########

#Function for reading in MCMCglmm models and combining model outputs for subsequent use in parameter estimation
combine_and_analyse <- function(variable) {
  # Construct directory path
  directory_path <- file.path("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_regressions/Model_outputs", variable, "1st_chain")
  
  # Get list of all RDS files in the specified directory
  model_files <- list.files(path = directory_path, pattern = "\\.rds", full.names = TRUE)
  
  # Read in all models using lapply()
  mcmc_list <- lapply(model_files, readRDS)
  
  # Extract the Sol objects from each model output and create a list of the mcmc objects
  sol_matrices <- lapply(mcmc_list, function(model) model$Sol)
  
  # Combine the mcmc objects together
  combined_sol <- do.call(rbind, sol_matrices)
  
  # Assign as an mcmc object
  MCMC_combined <- as.mcmc(combined_sol)
  
  #Check length
  print(length(MCMC_combined[,1]))
  
  #Summary outputs
  MCMC_combined_ESS <- effectiveSize(MCMC_combined)
  MCMC_combined_mode <- posterior.mode(MCMC_combined)
  print(MCMC_combined_mode)
  
  #Compute 95% credible interval
  MCMC_combined_HPD <- HPDinterval(MCMC_combined)
  print(MCMC_combined_HPD)
  
  return(list(
    ESS=MCMC_combined_ESS,
    mode=MCMC_combined_mode,
    HPD=MCMC_combined_HPD
  ))
}




####################################################################################################



############################################################
#Pairwise analyses predicting caste number using MCMCglmm
############################################################

##############################################################################
#Analyses predicting number of worker castes

############
#MF_Caste#
############
MF_Caste_result <- combine_and_analyse(variable = "MF_Caste")

############
#MFcat_Caste#
############
MFcat_Caste_result <- combine_and_analyse(variable = "MFcat_Caste")

############
#CS_Caste#
############
CS_Caste_result <- combine_and_analyse(variable = "CS_Caste")

############
#PG_Caste#
############
PG_Caste_result <- combine_and_analyse(variable = "PG_Caste")

############
#PGbinary_Caste#
############
PGbinary_Caste_result <- combine_and_analyse(variable = "PGbinary_Caste")

############
#PGcat_Caste#
#############
PGcat_Caste_result <- combine_and_analyse(variable = "PGcat_Caste")

##############################################################################
#Analyses predicting variation in worker size

############
#MF_siz_var#
############
MF_siz_var_result <- combine_and_analyse(variable = "MF_siz_var")

############
#MFcat_siz_var#
############
MFcat_siz_var_result <- combine_and_analyse(variable = "MFcat_siz_var")

############
#CS_siz_var#
############
CS_siz_var_result <- combine_and_analyse(variable = "CS_siz_var")

############
#PG_siz_var#
############
PG_siz_var_result <- combine_and_analyse(variable = "PG_siz_var")

############
#PGbinary_siz_var#
############
PGbinary_siz_var_result <- combine_and_analyse(variable = "PGbinary_siz_var")

############
#PGcat_siz_var#
############
PGcat_siz_var_result <- combine_and_analyse(variable = "PGcat_siz_var")

##############################################################################
#Analyses predicting variation in worker size in species with a single worker caste

############
#MF_mono_siz_var#
############
MF_mono_siz_var_result <- combine_and_analyse(variable = "MF_mono_siz_var")

############
#MFcat_mono_siz_var#
############
MFcat_mono_siz_var_result <- combine_and_analyse(variable = "MFcat_mono_siz_var")

############
#CS_mono_siz_var#
############
CS_mono_siz_var_result <- combine_and_analyse(variable = "CS_mono_siz_var")

############
#PG_mono_siz_var#
############
PG_mono_siz_var_result <- combine_and_analyse(variable = "PG_mono_siz_var")

############
#PGbinary_mono_siz_var#
############
PGbinary_mono_siz_var_result <- combine_and_analyse(variable = "PGbinary_mono_siz_var")

############
#PGcat_mono_siz_var#
############
PGcat_mono_siz_var_result <- combine_and_analyse(variable = "PGcat_mono_siz_var")

##############################################################################
##Pairwise analysis among predictor variables

############
#MF_CS#
############
MF_CS_result <- combine_and_analyse(variable = "MF_CS")

############
#MFcat_CS#
############
MFcat_CS_result <- combine_and_analyse(variable = "MFcat_CS")

############
#MF_PG#
############
MF_PG_result <- combine_and_analyse(variable = "MF_PG")

############
#MF_PGbinary#
############
MF_PGbinary_result <- combine_and_analyse(variable = "MF_PGbinary")

############
#MF_PGcat#
############
MF_PGcat_result <- combine_and_analyse(variable = "MF_PGcat")

############
#MFcat_PG#
############
MFcat_PG_result <- combine_and_analyse(variable = "MFcat_PG")

############
#CS_PG#
############
CS_PG_result <- combine_and_analyse(variable = "CS_PG")

############
#CS_PGbinary#
############
CS_PGbinary_result <- combine_and_analyse(variable = "CS_PGbinary")

############
#CS_PGcat#
############
CS_PGcat_result <- combine_and_analyse(variable = "CS_PGcat")



############ END ##########################







