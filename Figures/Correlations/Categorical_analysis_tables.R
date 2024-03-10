############################### ANOVA table for BPMM sensitivity analyses ########################################
###To quantifying phylogenetic uncertainty in parameter estimates, analyses run over a sample of 400 phylogenetic trees
##Louis Bell-Roberts
#01/03/2024


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

#Calculate statistics and produce table for the output summar of linear models with categorical predictor
calculate_stats <- function(result, data_df, predictor_var, response_var, analysis_labels, factor_labels) {
  # Convert result to dataframe
  result_df <- as.data.frame(result$HPD)
  
  # Calculate N_vector
  N_vector <- as.vector(table(data_df[predictor_var]))
  
  #Calculate min and max values for the predictor variable
  min_max_result <- by(data_df[response_var], data_df[predictor_var], FUN = function(x) c(Min = min(x), Max = max(x)))
  min_max_result <- do.call(rbind, min_max_result)
  
  # Bind parameter estimates to the dataframe
  result_df <- cbind(result_df, B = result$mode, N = N_vector, min_max_result[, 1:2],
                     Analysis = analysis_labels, Factor_level = factor_labels) 
  
  # Set row names to NULL
  rownames(result_df) <- NULL
  
  # Select relevant columns and rename them
  stats_final <- result_df %>%
    dplyr::select(Factor_level, B, lower, upper, N, Min, Max, Analysis) %>%
    rename(CI_lower = lower, CI_upper = upper)
  
  return(stats_final)
}


####################################################################################################


#Read in data file
ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Following_review/ant_data.csv")

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

# ant_data$colony.size <- log10(ant_data$colony.size)
# ant_data$queen.mating.frequency <- log10(ant_data$queen.mating.frequency)
# ant_data$queen.number.continuous <- log10(ant_data$queen.number.continuous)
# ant_data$worker.size.variation <- sqrt(ant_data$worker.size.variation)

#Rename 'species' column as 'animal'
ant_data <- ant_data %>% dplyr::rename(animal = species)

#Read in a single phylogenetic tree
ant.trees <- read.tree(file ="/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/15k_NCuniform_stem_mcc.tre")

###Subsets of the variables
##Analysis predicting number of worker castes
data_MFcat_caste <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency.categorical), complete.cases(caste.number))
data_PGcat_caste <- dplyr::filter(ant_data, complete.cases(queen.number.categorical), complete.cases(caste.number))

##Analysis predicting variation in worker size
data_MFcat_siz_var <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency.categorical), complete.cases(worker.size.variation))
data_PGcat_siz_var <- dplyr::filter(ant_data, complete.cases(queen.number.categorical), complete.cases(worker.size.variation))

##Analysis predicting variation in worker size in species with a single worker caste
data_MFcat_mono_siz_var <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency.categorical), complete.cases(worker.size.variation), caste.number <2)
data_PGcat_mono_siz_var <- dplyr::filter(ant_data, complete.cases(queen.number.categorical), complete.cases(worker.size.variation), caste.number <2)

##Pairwise analysis among predictor variables
data_MFcat_CS <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency.categorical), complete.cases(colony.size))
data_MF_PGcat <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(queen.number.categorical))
data_MFcat_PG <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency.categorical), complete.cases(queen.number.continuous))
data_CS_PGcat <- dplyr::filter(ant_data, complete.cases(colony.size), complete.cases(queen.number.categorical))

##Prune multiphylo objects for each of the different sets of predictor variables
#Analyses predicting caste number
PT_data_MFcat_caste <- drop.tip(ant.trees, setdiff(ant.trees$tip.label, data_MFcat_caste$animal))
PT_data_PGcat_caste <- drop.tip(ant.trees, setdiff(ant.trees$tip.label, data_PGcat_caste$animal))

#Analyses predicting variation in worker size
PT_data_MFcat_siz_var <- drop.tip(ant.trees, setdiff(ant.trees$tip.label, data_MFcat_siz_var$animal))
PT_data_PGcat_siz_var <- drop.tip(ant.trees, setdiff(ant.trees$tip.label, data_PGcat_siz_var$animal))

#Analyses predicting variation in worker size in species with a single worker caste
PT_data_MFcat_mono_siz_var <- drop.tip(ant.trees, setdiff(ant.trees$tip.label, data_MFcat_mono_siz_var$animal))
PT_data_PGcat_mono_siz_var <- drop.tip(ant.trees, setdiff(ant.trees$tip.label, data_PGcat_mono_siz_var$animal))

#Pairwise analyses among the predictor variables
PT_data_MFcat_CS <- drop.tip(ant.trees, setdiff(ant.trees$tip.label, data_MFcat_CS$animal))
PT_data_MF_PGcat <- drop.tip(ant.trees, setdiff(ant.trees$tip.label, data_MF_PGcat$animal))
PT_data_MFcat_PG <- drop.tip(ant.trees, setdiff(ant.trees$tip.label, data_MFcat_PG$animal))
PT_data_CS_PGcat <- drop.tip(ant.trees, setdiff(ant.trees$tip.label, data_CS_PGcat$animal))

##Prune database to match the tree
#Filter through dataframe and select only the rows that match the tips of the tree

#Analyses predicting caste number
data_MFcat_caste <- filter(data_MFcat_caste, animal %in% PT_data_MFcat_caste$tip.label)
data_PGcat_caste <- filter(data_PGcat_caste, animal %in% PT_data_PGcat_caste$tip.label)

#Analyses predicting variation in worker size
data_MFcat_siz_var <- filter(data_MFcat_siz_var, animal %in% PT_data_MFcat_siz_var$tip.label)
data_PGcat_siz_var <- filter(data_PGcat_siz_var, animal %in% PT_data_PGcat_siz_var$tip.label)

#Analyses predicting variation in worker size in species with a single worker caste
data_MFcat_mono_siz_var <- filter(data_MFcat_mono_siz_var, animal %in% PT_data_MFcat_mono_siz_var$tip.label)
data_PGcat_mono_siz_var <- filter(data_PGcat_mono_siz_var, animal %in% PT_data_PGcat_mono_siz_var$tip.label)

#Pairwise analyses among the predictor variables
data_MFcat_CS <- filter(data_MFcat_CS, animal %in% PT_data_MFcat_CS$tip.label)
data_MF_PGcat <- filter(data_MF_PGcat, animal %in% PT_data_MF_PGcat$tip.label)
data_MFcat_PG <- filter(data_MFcat_PG, animal %in% PT_data_MFcat_PG$tip.label)
data_CS_PGcat <- filter(data_CS_PGcat, animal %in% PT_data_CS_PGcat$tip.label)


############################################################
#Pairwise analyses predicting caste number using MCMCglmm
############################################################

##############################################################################
#Analyses predicting number of worker castes

############
#MFcat_Caste#
############
MFcat_Caste_result <- combine_and_analyse(variable = "MFcat_Caste")

############
#PGcat_Caste#
#############
PGcat_Caste_result <- combine_and_analyse(variable = "PGcat_Caste")

##############################################################################
#Analyses predicting variation in worker size

############
#MFcat_siz_var#
############
MFcat_siz_var_result <- combine_and_analyse(variable = "MFcat_siz_var")

############
#PGcat_siz_var#
############
PGcat_siz_var_result <- combine_and_analyse(variable = "PGcat_siz_var")

##############################################################################
#Analyses predicting variation in worker size in species with a single worker caste

############
#MFcat_mono_siz_var#
############
MFcat_mono_siz_var_result <- combine_and_analyse(variable = "MFcat_mono_siz_var")

############
#PGcat_mono_siz_var#
############
PGcat_mono_siz_var_result <- combine_and_analyse(variable = "PGcat_mono_siz_var")

##############################################################################
##Pairwise analysis among predictor variables

############
#MFcat_CS#
############
MFcat_CS_result <- combine_and_analyse(variable = "MFcat_CS")

############
#MF_PGcat#
############
MF_PGcat_result <- combine_and_analyse(variable = "MF_PGcat")

############
#MFcat_PG#
############
MFcat_PG_result <- combine_and_analyse(variable = "MFcat_PG")

############
#CS_PGcat#
############
CS_PGcat_result <- combine_and_analyse(variable = "CS_PGcat")



############################################################
#Produces summary statistic tables, including the global intercept as we want to compare the differences between factor levels, rather than between factor levels and the intercept
############################################################

##############################################################################
#Analyses predicting number of worker castes

############
#MFcat_Caste#
############
# test_mod <- readRDS(file = "/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_regressions/Model_outputs/MFcat_Caste/1st_chain/MFcat_Caste_1M_100k_1k_1_1stRun.rds")
# summary(test_mod)

# MFcat_Caste_result_df <- as.data.frame(MFcat_Caste_result$HPD)
# # MFcat_Caste_result_df <- rownames_to_column(MFcat_Caste_result_df, var = "Factor_level") #Assign row names to a column
# # MFcat_Caste_par_est <- c(MFcat_Caste_result$mode[1], MFcat_Caste_result$mode[1] + MFcat_Caste_result$mode[2], MFcat_Caste_result$mode[1] + MFcat_Caste_result$mode[3]) # Transform parameter estimates to their absolute values, rather than relative
# 
# N_vector <- as.vector(table(data_MFcat_caste$queen.mating.frequency.categorical))
# 
# # Use aggregate function to calculate min and max values for caste.number for each factor level of queen.mating.frequency.categorical
# min_max_result <- aggregate(caste.number ~ queen.mating.frequency.categorical, data = data_MFcat_caste, FUN = function(x) c(Min_caste = min(x), Max_caste = max(x)))
# 
# MFcat_Caste_result_df <- cbind(MFcat_Caste_result_df, B = MFcat_Caste_result$mode, N = N_vector, min_max_result$caste.number[, 1:2], Analysis = c("MF_Caste", "", ""), Factor_level = c("Monandry", "Facultative polyandry", "Obligate polyandry")) #Bind parameter estimates to the dataframe
# rownames(MFcat_Caste_result_df) <- NULL
# MFcat_Caste_stats_final <- MFcat_Caste_result_df %>% dplyr::select(Factor_level, B, lower, upper, N, Min_caste, Max_caste, Analysis) %>% rename(CI_lower = lower, CI_upper = upper)


####################################################

MFcat_Caste_stats <- calculate_stats(result = MFcat_Caste_result, data_df = data_MFcat_caste, predictor_var = "queen.mating.frequency.categorical", response_var = "caste.number", analysis_labels = c("Mating_frequency_Caste","",""), factor_labels = c("Monandry", "Facultative polyandry", "Obligate polyandry"))


############
#PGcat_Caste#
############
PGcat_Caste_stats <- calculate_stats(result = PGcat_Caste_result, data_df = data_PGcat_caste, predictor_var = "queen.number.categorical", response_var = "caste.number", analysis_labels = c("Queen_number_Caste","",""), factor_labels = c("Monogyny", "Facultative polygyny", "Obligate polygyny"))


##############################################################################
#Analyses predicting variation in worker size

############
#MFcat_siz_var#
############
MFcat_siz_var_stats <- calculate_stats(result = MFcat_siz_var_result, data_df = data_MFcat_siz_var, predictor_var = "queen.mating.frequency.categorical", response_var = "worker.size.variation", analysis_labels = c("Mating_frequency_Size_variation","",""), factor_labels = c("Monandry", "Facultative polyandry", "Obligate polyandry"))

############
#PGcat_siz_var#
############
PGcat_siz_var_stats <- calculate_stats(result = PGcat_siz_var_result, data_df = data_PGcat_siz_var, predictor_var = "queen.number.categorical", response_var = "worker.size.variation", analysis_labels = c("Queen_number_Size_variation","",""), factor_labels = c("Monogyny", "Facultative polygyny", "Obligate polygyny"))


##############################################################################
#Analyses predicting variation in worker size only in species with a single worker caste

############
#MFcat_mono_siz_var#
############
MFcat_mono_siz_var_stats <- calculate_stats(result = MFcat_mono_siz_var_result, data_df = data_MFcat_mono_siz_var, predictor_var = "queen.mating.frequency.categorical", response_var = "worker.size.variation", analysis_labels = c("Mating_frequency_Size_variation_monomorphic","",""), factor_labels = c("Monandry", "Facultative polyandry", "Obligate polyandry"))

############
#PGcat_mono_siz_var#
############
PGcat_mono_siz_var_stats <- calculate_stats(result = PGcat_mono_siz_var_result, data_df = data_PGcat_mono_siz_var, predictor_var = "queen.number.categorical", response_var = "worker.size.variation", analysis_labels = c("Queen_number_Size_variation_monomorphic","",""), factor_labels = c("Monogyny", "Facultative polygyny", "Obligate polygyny"))

##############################################################################
##Pairwise analysis among predictor variables

############
#MFcat_CS#
############
MFcat_CS_stats <- calculate_stats(result = MFcat_CS_result, data_df = data_MFcat_CS, predictor_var = "queen.mating.frequency.categorical", response_var = "colony.size", analysis_labels = c("Mating_frequency_Colony_size","",""), factor_labels = c("Monandry", "Facultative polyandry", "Obligate polyandry"))

############
#MFcat_PG#
############
MFcat_PG_stats <- calculate_stats(result = MFcat_PG_result, data_df = data_MFcat_PG, predictor_var = "queen.mating.frequency.categorical", response_var = "queen.number.continuous", analysis_labels = c("Mating_frequency_Queen_number","",""), factor_labels = c("Monandry", "Facultative polyandry", "Obligate polyandry"))


############
#CS_PGcat#
############
CS_PGcat_stats <- calculate_stats(result = CS_PGcat_result, data_df = data_CS_PGcat, predictor_var = "queen.number.categorical", response_var = "colony.size", analysis_labels = c("Queen_number_Colony_size","",""), factor_labels = c("Monogyny", "Facultative polygyny", "Obligate polygyny"))



############################################################
#Write tables to csv
############################################################

#Produce individual tables
write.csv(MFcat_Caste_stats, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Correlations/Tables_categorical_predictors/MFcat_Caste_stats.csv", row.names = F)
write.csv(PGcat_Caste_stats, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Correlations/Tables_categorical_predictors/PGcat_Caste_stats.csv", row.names = F)
write.csv(MFcat_siz_var_stats, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Correlations/Tables_categorical_predictors/MFcat_siz_var_stats.csv", row.names = F)
write.csv(PGcat_siz_var_stats, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Correlations/Tables_categorical_predictors/PGcat_siz_var_stats.csv", row.names = F)
write.csv(MFcat_mono_siz_var_stats, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Correlations/Tables_categorical_predictors/MFcat_mono_siz_var_stats.csv", row.names = F)
write.csv(PGcat_mono_siz_var_stats, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Correlations/Tables_categorical_predictors/PGcat_mono_siz_var_stats.csv", row.names = F)
write.csv(MFcat_CS_stats, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Correlations/Tables_categorical_predictors/MFcat_CS_stats.csv", row.names = F)
write.csv(MFcat_PG_stats, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Correlations/Tables_categorical_predictors/MFcat_PG_stats.csv", row.names = F)
write.csv(CS_PGcat_stats, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Correlations/Tables_categorical_predictors/CS_PGcat_stats.csv", row.names = F)

#Combine tables and write to csv file
combined_stats_df <- rbind(MFcat_Caste_stats, PGcat_Caste_stats, MFcat_siz_var_stats, PGcat_siz_var_stats, MFcat_mono_siz_var_stats, PGcat_mono_siz_var_stats, MFcat_CS_stats, MFcat_PG_stats, CS_PGcat_stats)

#Round values to 2 decminal places
combined_stats_df <- combined_stats_df %>%
  mutate_at(vars(B, CI_lower, CI_upper, N, Min, Max), round, 2)

combined_stats_df$Max <- format(combined_stats_df$Max, scientific = FALSE)

write.csv(combined_stats_df, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Correlations/Tables_categorical_predictors/combined_stats_df.csv", row.names = F)



############ END ##########################







