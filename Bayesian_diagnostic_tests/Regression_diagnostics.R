######################## BPMM regression analysis diagnostic checks over 400 models ######################
###For all pairwise regression models
##Run diagnostic checks: ESS, Convergence test (Gelman-Rubin scale reduction factor test (SRF)), autocorrelation
#Louis Bell-Roberts
#01/02/2024

#ESS never lower than 400 for any model and typically at 1000


library(tidyverse)
library(MCMCglmm)
library(ggplot2)
library(gtools)

############
#Custom functions
############

###
#Write a function that reads in all of the individual models
model_reader_function <- function(variable, run) {
  # Construct directory path
  directory_path <- file.path("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_regressions/Model_outputs", variable, run)
  
  # Get list of all rda files in the specified directory
  model_files <- list.files(path = directory_path, pattern = ".rds")
  
  # Sort the file names alphanumerically
  sorted_files <- mixedsort(model_files)
  
  # Read in all models
  # Create an empty list to store the loaded data
  mcmc_list <- list()
  
  # Read in each MCMCglmm model
  mcmc_list <- lapply(sorted_files, function(file) {
    readRDS(file.path(directory_path, file))
  })  # mcmc_list is the list of all the models read in
  
  return(mcmc_list)
}



###
#Convergence test function
gelman_rubin_function <- function(mcmc_list_1, mcmc_list_2, component = "Sol") {
  # Create a list to store the paired lists
  paired_lists <- list()
  
  # Iterate over the indices of the lists
  for (i in seq_along(mcmc_list_1)) {
    if (component == "Sol") {
      # Select the "$Sol" element from each object in list_a and list_b
      pair_component <- list(mcmc_list_1[[i]]$Sol, mcmc_list_2[[i]]$Sol)
    } else if (component == "VCV") {
      # Select the "$VCV" element from each object in list_a and list_b
      pair_component <- list(mcmc_list_1[[i]]$VCV, mcmc_list_2[[i]]$VCV)
    } else {
      stop("Invalid component argument. Use 'Sol' or 'VCV'.")
    }
    
    # Add the pair to the list of paired lists
    paired_lists[[i]] <- pair_component
  }
  
  # Run gelman.diag function over the list of paired chain runs
  gelman_results <- lapply(paired_lists, gelman.diag)
  
  # Initialize an empty data frame
  component_df <- data.frame()
  
  # Combine matrices into a single data frame
  for (i in seq_along(gelman_results)) {
    # Convert the matrix to a data frame
    component_temporary_df <- as.data.frame(gelman_results[[i]]$psrf)
    
    # Add an identifier column to keep track of the original matrix
    component_temporary_df$matrix_id <- i
    
    # Append the temporary data frame to the combined data frame
    component_df <- rbind(component_df, component_temporary_df)
  }
  
  # Return a list containing both the component data frame and gelman_results
  return(list(component_df = component_df, gelman_results = gelman_results))
}



###
#Autocorrelation test function
autocorr_analysis_function <- function(mcmc_list, component = "Sol") {
# Apply the autocorr.diag function to each element of the list
autocorr_list <- lapply(mcmc_list, function(x) {
  if (component == "Sol") {
    autocorr.diag(x$Sol)
  } else if (component == "VCV") {
    autocorr.diag(x$VCV)
  } else {
    stop("Invalid component argument. Use 'Sol' or 'VCV'.")
  }
})

# Use lapply to extract the second row of each matrix
second_row_list <- lapply(autocorr_list, function(x) {
  x[2,]
})

# Convert the list to a vector
second_row_vec <- unlist(second_row_list)

# Return the quantiles of the second row vector
return(quantile(second_row_vec))
}


###
#ESS calculation function
calculate_ess_function <- function(mcmc_list, component = "Sol") {
  # Apply the effectiveSize function to each element of the list
  ess_list <- lapply(mcmc_list, function(x) {
    if (component == "Sol") {
      effectiveSize(x$Sol)
    } else if (component == "VCV") {
      effectiveSize(x$VCV)
    } else {
      stop("Invalid component argument. Use 'Sol' or 'VCV'.")
    }
  })
  
  # Convert the list to a vector
  ess_vec <- unlist(ess_list)
  
  # Return the quantiles of the ESS vector
  return(quantile(ess_vec))
}

###
# Function that runs the entire suite of diagnostic checks on a model
perform_analysis <- function(variable_pair) {
  # Read in the models
  mcmc_list_1 <- model_reader_function(variable_pair, "1st_chain")
  mcmc_list_2 <- model_reader_function(variable_pair, "2nd_chain")
  
  # Convergence tests
  sol_conv <- gelman_rubin_function(mcmc_list_1, mcmc_list_2, component = "Sol")
  vcv_conv <- gelman_rubin_function(mcmc_list_1, mcmc_list_2, component = "VCV")
  
  cat("Quantiles of Gelman-Rubin convergence tests for Sol component:\n")
  quantiles_sol_conv <- quantile(unlist(sol_conv$gelman_results))
  print(quantiles_sol_conv)
  
  cat("Quantiles of Gelman-Rubin convergence tests for VCV component:\n")
  quantiles_vcv_conv <- quantile(unlist(vcv_conv$gelman_results))
  print(quantiles_vcv_conv)
  
  cat("Gelman-Rubin component data frames for Sol and VCV components:\n")
  sol_conv_df <- sol_conv$component_df
  vcv_conv_df <- vcv_conv$component_df
  
  # Autocorrelation test
  sol_autocorr <- autocorr_analysis_function(mcmc_list_1, component = "Sol")
  vcv_autocorr <- autocorr_analysis_function(mcmc_list_1, component = "VCV")
  
  cat("Quantiles of autocorrelation tests for Sol component:\n")
  quantiles_sol_autocorr <- quantile(unlist(sol_autocorr))
  print(quantiles_sol_autocorr)
  
  cat("Quantiles of autocorrelation tests for VCV component:\n")
  quantiles_vcv_autocorr <- quantile(unlist(vcv_autocorr))
  print(quantiles_vcv_autocorr)
  
  # ESS test
  sol_ESS <- calculate_ess_function(mcmc_list = mcmc_list_1, component = "Sol")
  vcv_ESS <- calculate_ess_function(mcmc_list = mcmc_list_1, component = "VCV")
  
  cat("Quantiles of ESS tests for Sol component:\n")
  quantiles_sol_ESS <- quantile(unlist(sol_ESS))
  print(quantiles_sol_ESS)
  
  cat("Quantiles of ESS tests for VCV component:\n")
  quantiles_vcv_ESS <- quantile(unlist(vcv_ESS))
  print(quantiles_vcv_ESS)
  
  # Return the results as a list
  return(list(
    quantiles_sol_conv = quantiles_sol_conv,
    quantiles_vcv_conv = quantiles_vcv_conv,
    sol_conv_df = sol_conv_df,
    vcv_conv_df = vcv_conv_df,
    quantiles_sol_autocorr = quantiles_sol_autocorr,
    quantiles_vcv_autocorr = quantiles_vcv_autocorr,
    quantiles_sol_ESS = quantiles_sol_ESS,
    quantiles_vcv_ESS = quantiles_vcv_ESS
  ))
}

############################################################################################################################################################



############################################################
#Pairwise analyses predicting caste number using MCMCglmm
############################################################

##############################################################################
#Analyses predicting number of worker castes

############
#MF_Caste - performed step by step to test that it's working
############
#Read in the models
MF_Caste_mcmc_list_1 <- model_reader_function("MF_Caste","1st_chain")
MF_Caste_mcmc_list_2 <- model_reader_function("MF_Caste","2nd_chain")
#Convergence tests
MF_Caste_sol_conv <- gelman_rubin_function(mcmc_list_1 = MF_Caste_mcmc_list_1, mcmc_list_2 = MF_Caste_mcmc_list_2, component = "Sol")
MF_Caste_vcv_conv <- gelman_rubin_function(mcmc_list_1 = MF_Caste_mcmc_list_1, mcmc_list_2 = MF_Caste_mcmc_list_2, component = "VCV")
quantile(unlist(MF_Caste_sol_conv$gelman_results))
quantile(unlist(MF_Caste_vcv_conv$gelman_results))
View(MF_Caste_sol_conv$component_df)
View(MF_Caste_vcv_conv$component_df)
#Autocorrelation test
MF_Caste_sol_autocorr <- autocorr_analysis_function(MF_Caste_mcmc_list_1, component = "Sol")
MF_Caste_vcv_autocorr <- autocorr_analysis_function(MF_Caste_mcmc_list_1, component = "VCV")
#ESS test
MF_Caste_sol_ESS <- calculate_ess_function(mcmc_list = MF_Caste_mcmc_list_1, component = "Sol")
MF_Caste_vcv_ESS <- calculate_ess_function(mcmc_list = MF_Caste_mcmc_list_1, component = "VCV")
#Remove the models from the environment before moving onto the next analysis as they are large files
rm(MF_Caste_mcmc_list_1)
rm(MF_Caste_mcmc_list_2)

############
#MFcat_Caste#
############
MFcat_Caste_result <- perform_analysis("MFcat_Caste") #This should be sufficient for each pairwise analysis as long as printing the quantile is sufficient

############
#CS_Caste#
############
CS_Caste_result <- perform_analysis("CS_Caste")

############
#PG_Caste#
############
PG_Caste_result <- perform_analysis("PG_Caste")

############
#PGbinary_Caste#
############
PGbinary_Caste_result <- perform_analysis("PGbinary_Caste")

############
#PGcat_Caste#
#############
PGcat_Caste_result <- perform_analysis("PGcat_Caste")


##############################################################################
#Analyses predicting variation in worker size

############
#MF_siz_var#
############
MF_siz_var_result <- perform_analysis("MF_siz_var")

############
#MFcat_siz_var#
############
MFcat_siz_var_result <- perform_analysis("MFcat_siz_var")

############
#CS_siz_var#
############
CS_siz_var_result <- perform_analysis("CS_siz_var")

############
#PG_siz_var#
############
PG_siz_var_result <- perform_analysis("PG_siz_var")

############
#PGbinary_siz_var#
############
PGbinary_siz_var_result <- perform_analysis("PGbinary_siz_var")

############
#PGcat_siz_var#
############
PGcat_siz_var_result <- perform_analysis("PGcat_siz_var")

##############################################################################
#Analyses predicting variation in worker size in species with a single worker caste

############
#MF_mono_siz_var#
############
MF_mono_siz_var_result <- perform_analysis("MF_mono_siz_var")

############
#MFcat_mono_siz_var#
############
MFcat_mono_siz_var_result <- perform_analysis("MFcat_mono_siz_var")

############
#CS_mono_siz_var#
############
CS_mono_siz_var_result <- perform_analysis("CS_mono_siz_var")
View(CS_mono_siz_var_result$vcv_conv_df)

############
#PG_mono_siz_var#
############
PG_mono_siz_var_result <- perform_analysis("PG_mono_siz_var")
View(PG_mono_siz_var_result$vcv_conv_df)

############
#PGbinary_mono_siz_var#
############
PGbinary_mono_siz_var_result <- perform_analysis("PGbinary_mono_siz_var")

############
#PGcat_mono_siz_var#
############
PGcat_mono_siz_var_result <- perform_analysis("PGcat_mono_siz_var")

##############################################################################
##Pairwise analysis among predictor variables

############
#MF_CS#
############
MF_CS_result <- perform_analysis("MF_CS")

############
#MFcat_CS#
############
MFcat_CS_result <- perform_analysis("MFcat_CS")

############
#MF_PG#
############
MF_PG_result <- perform_analysis("MF_PG")

############
#MF_PGbinary#
############
MF_PGbinary_result <- perform_analysis("MF_PGbinary")

############
#MF_PGcat#
############
MF_PGcat_result <- perform_analysis("MF_PGcat")
View(MF_PGcat_result$sol_conv_df)

############
#MFcat_PG#
############
MFcat_PG_result <- perform_analysis("MFcat_PG")

############
#CS_PG#
############
CS_PG_result <- perform_analysis("CS_PG")

############
#CS_PGbinary#
############
CS_PGbinary_result <- perform_analysis("CS_PGbinary")
# View(CS_PGbinary_result$sol_conv_df)

############
#CS_PGcat#
############
CS_PGcat_result <- perform_analysis("CS_PGcat")


