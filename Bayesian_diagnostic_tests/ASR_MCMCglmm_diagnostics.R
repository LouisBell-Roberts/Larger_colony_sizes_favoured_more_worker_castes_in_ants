######################## MCMCglmm ASR analyses diagnostic checks over 400 models ######################
##Run diagnostic checks: ESS, Convergence test (Gelman-Rubin scale reduction factor test (SRF)), autocorrelation
#Diagnostic checks need to take place on the nodes as well as the fixed and random effects
#Louis Bell-Roberts
#01/02/2024

library(tidyverse)
library(MCMCglmm)
library(ggplot2)
library(mulTree)
library(gtools)

##ESS was not lower than 300 for any parameter for any model



#################
#Custom functions
#################

#Write a function that reads in all of the individual models
model_reader_function <- function() {
  
  # get list of all RDS files in the working directory
  model_files <- list.files(pattern = "\\.rds$")

  # sort the file names based on the corresponding numbers
  sorted_file_names <- mixedsort(model_files)
  
  # read in all models using lapply()
  mcmc_list <- lapply(sorted_file_names, readRDS)
  #mcmc_list is the list of all of the models read in
  return(mcmc_list)
}



###
#Write a function that tests levels of autocorrelation for the fixed effects for each model - provides the values in the second row of the autocorr.diag function. These should be <0.1

#Function for calculating autocorrelation for the transition category fixed effects
auto_cor_function_fix_trans <- function(mcmc_list){
  
  # Apply the autocorr.diag function to each element of the list
  autocorr_list <- lapply(mcmc_list, function(x) {
    autocorr.diag(x$Sol[, 1:4])
  })
  
  # Use lapply to extract the second row of each matrix
  second_row_list <- lapply(autocorr_list, function(x) {
    x[2,]
  })
  
  # Convert the list to a vector
  second_row_vec <- unlist(second_row_list)
  return(second_row_vec)
}


#Function for calculating autocorrelation for the node fixed effects for the CS_Caste ASR
auto_cor_function_fix_nodes_CS <- function(mcmc_list){
  
  # Apply the autocorr.diag function to each element of the list
  autocorr_list <- lapply(mcmc_list, function(x) {
    autocorr.diag(x$Sol[, 5:438])
  })
  
  # Use lapply to extract the second row of each matrix
  second_row_list <- lapply(autocorr_list, function(x) {
    x[2,]
  })
  
  # Convert the list to a vector
  second_row_vec <- unlist(second_row_list)
  return(second_row_vec)
}


#Function for calculating autocorrelation for the node fixed effects for the MF_Caste ASR
auto_cor_function_fix_nodes_MF <- function(mcmc_list){
  
  # Apply the autocorr.diag function to each element of the list
  autocorr_list <- lapply(mcmc_list, function(x) {
    autocorr.diag(x$Sol[, 5:106]) #columns 5:106 represents the nodes in the $Sol object
  })
  
  # Use lapply to extract the second row of each matrix
  second_row_list <- lapply(autocorr_list, function(x) {
    x[2,]
  })
  
  # Convert the list to a vector
  second_row_vec <- unlist(second_row_list)
  return(second_row_vec)
}





#Write a function that tests levels of autocorrelation for the random effects for each model - provides the values in the second row of the autocorr.diag function. These should be <0.1
auto_cor_function_rand <- function(mcmc_list){
  
  # Apply the autocorr.diag function to each element of the list
  autocorr_list <- lapply(mcmc_list, function(x) {
    autocorr.diag(x$VCV)
  })
  
  # Use lapply to extract the second row of each matrix
  second_row_list <- lapply(autocorr_list, function(x) {
    x[2,]
  })
  
  # Convert the list to a vector
  second_row_vec <- unlist(second_row_list)
  return(second_row_vec)
}



#Two functions: one that calculates the ESS values for the fixed effects for each model across the full sample of models and one that calculates the ESS values for the random effects for each model across the full sample of models

##Fixed - transition categories
ESS_function_fix_trans <- function(mcmc_list) {
  
  #Create an empty list to store values
  result_list <- list()
  
  #Calculate ESS values for each fixed effect and node for each of the 400 models
  for (i in seq_along(mcmc_list)) {
    result_list[[i]] <- effectiveSize(mcmc_list[[i]]$Sol[, 1:4])
  }
  
  #Unlist the values to transform them into a vector
  ESS_vec <- unlist(result_list)
  
  return(ESS_vec)
}


##Fixed - nodes for the CS_Caste ASR
ESS_function_fix_nodes_CS <- function(mcmc_list) {
  
  #Create an empty list to store values
  result_list <- list()
  
  #Calculate ESS values for each fixed effect and node for each of the 400 models
  for (i in seq_along(mcmc_list)) {
    result_list[[i]] <- effectiveSize(mcmc_list[[i]]$Sol[, 5:438])
  }
  
  #Unlist the values to transform them into a vector
  ESS_vec <- unlist(result_list)
  
  return(ESS_vec)
}


##Fixed - nodes for the MF_Caste ASR
ESS_function_fix_nodes_MF <- function(mcmc_list) {
  
  #Create an empty list to store values
  result_list <- list()
  
  #Calculate ESS values for each fixed effect and node for each of the 400 models
  for (i in seq_along(mcmc_list)) {
    result_list[[i]] <- effectiveSize(mcmc_list[[i]]$Sol[, 5:106]) #columns 5:106 represents the nodes in the $Sol object
  }
  
  #Unlist the values to transform them into a vector
  ESS_vec <- unlist(result_list)
  
  return(ESS_vec)
}



##Random
ESS_function_rand <- function(mcmc_list) {
  
  #Create an empty list to store values
  result_list <- list()
  
  #Calculate ESS values for each fixed effect and node for each of the 400 models
  for (i in seq_along(mcmc_list)) {
    result_list[[i]] <- effectiveSize(mcmc_list[[i]]$VCV)
  }
  
  #Unlist the values to transform them into a vector
  ESS_vec <- unlist(result_list)
  
  return(ESS_vec)
}

##############################################################################################################################






############################################
#MF_Caste
############################################

#Read in all of the individual models

#Set working directory
setwd("/Volumes/ADATA SE800/DOL_worker_castes/ASR/MCMCglmm/MF_Caste/1st_run/")
mcmc_list_MF_Caste <- model_reader_function()

#Read in second set of models
setwd("/Volumes/ADATA SE800/DOL_worker_castes/ASR/MCMCglmm/MF_Caste/2nd_run/")
mcmc_list_MF_Caste_2nd_run <- model_reader_function()



#######
#Calculate whether scale-reduction factor (Gelman-Rubin) is <1.1 for all 400 models
#######

# Create a list to store the paired lists
MF_paired_lists_Sol <- list()

# Iterate over the indices of the lists
for (i in 1:length(mcmc_list_MF_Caste)) {
  # Select the "$Sol" element from each object in list_a and list_b
  MF_pair_Sol <- list(mcmc_list_MF_Caste[[i]]$Sol, mcmc_list_MF_Caste_2nd_run[[i]]$Sol)
  
  # Add the pair to the list of paired lists
  MF_paired_lists_Sol[[i]] <- MF_pair_Sol
}

MF_gelman_results_Sol <- lapply(MF_paired_lists_Sol, gelman.diag) #Takes a few minutes to run

quantile(unlist(unlist(MF_gelman_results_Sol)))

#Quantify whether each Sol passes srf test
##Create a single large data frame which holds all of the srf scores, and can then check whether they exceed 1.1

# Initialize an empty data frame
MF_Sol_df <- data.frame()
# Combine matrices into a single data frame
for (i in 1:length(MF_gelman_results_Sol)) {
  # Convert the matrix to a data frame
  MF_Sol_temp_df <- as.data.frame(MF_gelman_results_Sol[[i]]$psrf)
  
  # Add an identifier column to keep track of the original matrix
  MF_Sol_temp_df$matrix_id <- i
  
  # Append the temporary data frame to the combined data frame
  MF_Sol_df <- rbind(MF_Sol_df, MF_Sol_temp_df)
} #Takes 30 seconds to run
#Nearly all Sols <1.1 - should be checking the upper CI not the point estimate


##Repeat process for the VCV objects

MF_paired_lists_VCV <- list()

# Iterate over the indices of the lists
for (i in 1:length(mcmc_list_MF_Caste)) {
  # Select the "$Sol" element from each object in list_a and list_b
  MF_pair_VCV <- list(mcmc_list_MF_Caste[[i]]$VCV, mcmc_list_MF_Caste_2nd_run[[i]]$VCV)
  
  # Add the pair to the list of paired lists
  MF_paired_lists_VCV[[i]] <- MF_pair_VCV
}

MF_gelman_results_VCV <- lapply(MF_paired_lists_VCV, gelman.diag)
quantile(unlist(MF_gelman_results_VCV))

#Quantify whether each VCV passes srf test
##Create a single large data frame which holds all of the srf scores, and can then check whether they exceed 1.1

# Initialize an empty data frame
MF_VCV_df <- data.frame()
# Combine matrices into a single data frame
for (i in 1:length(MF_gelman_results_VCV)) {
  # Convert the matrix to a data frame
  MF_VCV_temp_df <- as.data.frame(MF_gelman_results_VCV[[i]]$psrf)
  
  # Add an identifier column to keep track of the original matrix
  MF_VCV_temp_df$matrix_id <- i
  
  # Append the temporary data frame to the combined data frame
  MF_VCV_df <- rbind(MF_VCV_df, MF_VCV_temp_df)
} 
#All Sols <1.1 - should be checking the upper CI, not the point estimate



########
#Test levels of autocorrelation for the fixed and random effects for each model 
########


#Estimate levels of autocorrelation for both the fixed and random effects for each model
##Supply the list of models to the function

#Fixed effects - transition categories
MF_Caste_auto_cor_est_fix_trans <- auto_cor_function_fix_trans(mcmc_list_MF_Caste) 
quantile(MF_Caste_auto_cor_est_fix_trans) #Check the quantiles for the estimates.
View(as.data.frame(MF_Caste_auto_cor_est_fix_trans))

#Fixed effects - nodes (columns 5:106)
MF_Caste_auto_cor_est_fix_nodes <- auto_cor_function_fix_nodes_MF(mcmc_list_MF_Caste) #Takes a few minutes to run
quantile(MF_Caste_auto_cor_est_fix_nodes) #Check the quantiles for the estimates
View(as.data.frame(MF_Caste_auto_cor_est_fix_nodes))

#Random effects
MF_Caste_auto_cor_est_rand <- auto_cor_function_rand(mcmc_list_MF_Caste)
quantile(MF_Caste_auto_cor_est_rand) #Check the quantiles for the estimates
View(as.data.frame(MF_Caste_auto_cor_est_rand))




#######
#Calculate effective sample size for both the nodes, fixed and random effects for each model
#######

##Fixed effects
#Transition categories
ESS_vec_MF_Caste_fix_trans <- ESS_function_fix_trans(mcmc_list_MF_Caste)
quantile(ESS_vec_MF_Caste_fix_trans)
View(as.data.frame(ESS_vec_MF_Caste_fix_trans)) 

#Nodes (columns 5:106)
ESS_vec_MF_Caste_fix_nodes <- ESS_function_fix_nodes_MF(mcmc_list_MF_Caste)
quantile(ESS_vec_MF_Caste_fix_nodes) 
View(as.data.frame(ESS_vec_MF_Caste_fix_nodes)) 

##Random effects
#Only two columns for the random effects (animal and units)
ESS_vec_MF_Caste_rand <- ESS_function_rand(mcmc_list_MF_Caste)
quantile(ESS_vec_MF_Caste_rand) 
View(as.data.frame(ESS_vec_MF_Caste_rand))






############################################
#CS_Caste
############################################

#Read in the models

#Set working directory
setwd("/Volumes/ADATA SE800/DOL_worker_castes/ASR/MCMCglmm/CS_Caste/1st_run/")
mcmc_list_CS_Caste <- model_reader_function()

#Read in second set of models
setwd("/Volumes/ADATA SE800/DOL_worker_castes/ASR/MCMCglmm/CS_Caste/2nd_run/")
mcmc_list_CS_Caste_2nd_run <- model_reader_function()



#######
#Calculate whether scale-reduction factor (Gelman-Rubin) is <1.1 for all 400 models
#######

#Create a list of paired mcmc objects

# Create a list to store the paired lists
paired_lists_Sol <- list()

# Iterate over the indices of the lists
for (i in 1:length(mcmc_list_CS_Caste)) {
  # Select the "$Sol" element from each object in list_a and list_b
  CS_pair_Sol <- list(mcmc_list_CS_Caste[[i]]$Sol, mcmc_list_CS_Caste_2nd_run[[i]]$Sol)
  
  # Add the pair to the list of paired lists
  paired_lists_Sol[[i]] <- CS_pair_Sol
}

#Run gelman.diag function over the list of paired chain runs
gelman_results_Sol <- lapply(paired_lists_Sol, gelman.diag) #Takes a while to run - max 20 minutes

quantile(unlist(gelman_results_Sol))

#Quantify whether each Sol passes srf test
##Create a single large data frame which holds all of the srf scores, and can then check whether they exceed 1.1

# Initialize an empty data frame
CS_Sol_df <- data.frame()
# Combine matrices into a single data frame
for (i in 1:length(gelman_results_Sol)) {
  # Convert the matrix to a data frame
  CS_Sol_temp_df <- as.data.frame(gelman_results_Sol[[i]]$psrf)
  
  # Add an identifier column to keep track of the original matrix
  CS_Sol_temp_df$matrix_id <- i
  
  # Append the temporary data frame to the combined data frame
  CS_Sol_df <- rbind(CS_Sol_df, CS_Sol_temp_df)
} #Takes 30 seconds to run
#All Sols <1.1 NEED TO LOOK AT UPPER CI


##Repeat process for the VCV objects

#Create a list of paired mcmc objects

CS_paired_lists_VCV <- list()

# Iterate over the indices of the lists
for (i in 1:length(mcmc_list_CS_Caste)) {
  # Select the "$Sol" element from each object in list_a and list_b
  CS_pair_VCV <- list(mcmc_list_CS_Caste[[i]]$VCV, mcmc_list_CS_Caste_2nd_run[[i]]$VCV)
  
  # Add the pair to the list of paired lists
  CS_paired_lists_VCV[[i]] <- CS_pair_VCV
}

#Calculate srf scores for each mcmc object pair in the list
CS_gelman_results_VCV <- lapply(CS_paired_lists_VCV, gelman.diag)
quantile(unlist(CS_gelman_results_VCV))

#Quantify whether each VCV passes srf test
##Create a single large data frame which holds all of the srf scores, and can then check whether they exceed 1.1

# Initialize an empty data frame
CS_VCV_df <- data.frame()
# Combine matrices into a single data frame
for (i in 1:length(CS_gelman_results_VCV)) {
  # Convert the matrix to a data frame
  CS_VCV_temp_df <- as.data.frame(CS_gelman_results_VCV[[i]]$psrf)
  
  # Add an identifier column to keep track of the original matrix
  CS_VCV_temp_df$matrix_id <- i
  
  # Append the temporary data frame to the combined data frame
  CS_VCV_df <- rbind(CS_VCV_df, CS_VCV_temp_df)
} #All VCVs <1.1 NEED TO LOOK AT UPPER CI


########
#Test levels of autocorrelation for the fixed and random effects for each model 
########


#Estimate levels of autocorrelation for both the fixed and random effects for each model
##Supply to the function the list of models

#Fixed effects - transition categories
CS_Caste_auto_cor_est_fix_trans <- auto_cor_function_fix_trans(mcmc_list_CS_Caste) 
quantile(CS_Caste_auto_cor_est_fix_trans) #Check the quantiles for the estimates
View(as.data.frame(CS_Caste_auto_cor_est_fix_trans))

#Fixed effects - nodes
CS_Caste_auto_cor_est_fix_nodes <- auto_cor_function_fix_nodes_CS(mcmc_list_CS_Caste) #Takes a long time to run but does work
quantile(CS_Caste_auto_cor_est_fix_nodes) #Check the quantiles for the estimates. 0.1373388220 is the largest value
View(as.data.frame(CS_Caste_auto_cor_est_fix_nodes))

#Random effects
CS_Caste_auto_cor_est_rand <- auto_cor_function_rand(mcmc_list_CS_Caste)
quantile(CS_Caste_auto_cor_est_rand) #Check the quantiles for the estimates. 0.1318994383 is the largest value
View(as.data.frame(CS_Caste_auto_cor_est_rand))




#######
#Calculate effective sample size for both the nodes, fixed and random effects for each model
#######

##Fixed effects
#Transition categories
#438 is the last node of the phylogeny and there are 4 transition categories
ESS_vec_CS_Caste_fix_trans <- ESS_function_fix_trans(mcmc_list_CS_Caste)
quantile(ESS_vec_CS_Caste_fix_trans) #lowest = 446.2081
View(as.data.frame(ESS_vec_CS_Caste_fix_trans)) #majority of values are >= 1000

#Nodes
ESS_vec_CS_Caste_fix_nodes <- ESS_function_fix_nodes_CS(mcmc_list_CS_Caste)
quantile(ESS_vec_CS_Caste_fix_nodes) #lowest = 229.0684
View(as.data.frame(ESS_vec_CS_Caste_fix_nodes)) #majority of values are >= 1000

##Random effects
#Only two columns for the random effects (animal and units)
ESS_vec_CS_Caste_rand <- ESS_function_rand(mcmc_list_CS_Caste)
quantile(ESS_vec_CS_Caste_rand) #lowest = 543.2235
View(as.data.frame(ESS_vec_CS_Caste_rand)) #majority of values are >= 1000









