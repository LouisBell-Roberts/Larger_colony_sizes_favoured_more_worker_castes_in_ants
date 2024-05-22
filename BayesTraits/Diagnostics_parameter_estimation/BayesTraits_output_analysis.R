########################Analysing outputs from BayesTraits######################
###MF~Caste
#20/01/2024

library(tidyverse)
library(coda)
setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/BayesTraits/BayesTraitsV4.0.1-OSX/Post_review_analyses/")

##################
#Functions
bayestraits_statistics <- function(ind_1st_file_path, ind_2nd_file_path, dep_1st_file_path, dep_2nd_file_path,
                                   dep_columns = c("Lh", "q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43"), #Select the parameters to analyse for the gelman.diag function for the independent and dependent model runs
                                   ind_columns = c("Lh", "alpha1", "beta1", "alpha2", "beta2")) {
  # Read in independent and dependent models
  BT_dep_1st <- read.csv(dep_1st_file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  BT_dep_2nd <- read.csv(dep_2nd_file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  BT_ind_1st <- read.csv(ind_1st_file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  BT_ind_2nd <- read.csv(ind_2nd_file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Calculate the number of unique trees sampled in the posterior
  dep_tree_count <- length(unique(BT_dep_1st$Tree.No))
  ind_tree_count <- length(unique(BT_ind_1st$Tree.No))
  
  # Convert data frames to mcmc objects
  BT_dep_1st_mcmc <- BT_dep_1st %>% dplyr::select(Lh, q12, q13, q21, q24, q31, q34, q42, q43)
  BT_dep_2nd_mcmc <- BT_dep_2nd %>% dplyr::select(Lh, q12, q13, q21, q24, q31, q34, q42, q43)
  BT_dep_1st_mcmc <- as.mcmc(BT_dep_1st_mcmc)
  BT_dep_2nd_mcmc <- as.mcmc(BT_dep_2nd_mcmc)
  
  BT_ind_1st_mcmc <- BT_ind_1st %>% dplyr::select(Lh, alpha1, beta1, alpha2, beta2)
  BT_ind_2nd_mcmc <- BT_ind_2nd %>% dplyr::select(Lh, alpha1, beta1, alpha2, beta2)
  BT_ind_1st_mcmc <- as.mcmc(BT_ind_1st_mcmc)
  BT_ind_2nd_mcmc <- as.mcmc(BT_ind_2nd_mcmc)
  
  # Calculating effective sample size
  dep_ESS <- effectiveSize(BT_dep_1st_mcmc)
  ind_ESS <- effectiveSize(BT_ind_1st_mcmc)
  
  # Calculating autocorrelation between samples
  dep_autocorr <- autocorr.diag(BT_dep_1st_mcmc)
  ind_autocorr <- autocorr.diag(BT_ind_1st_mcmc)
  
  # Calculate HPD interval
  dep_hpd_interval <- HPDinterval(BT_dep_1st_mcmc)
  ind_hpd_interval <- HPDinterval(BT_ind_1st_mcmc)
  
  # Calculate mean and median
  BT_dep_stats <- BT_dep_1st %>% dplyr::select(Lh, q12, q13, q21, q24, q31, q34, q42, q43)
  BT_ind_stats <- BT_ind_1st %>% dplyr::select(Lh, alpha1, beta1, alpha2, beta2)
  
  dep_mean <- apply(BT_dep_stats, 2, mean)
  dep_median <- apply(BT_dep_stats, 2, median)
  
  ind_mean <- apply(BT_ind_stats, 2, mean)
  ind_median <- apply(BT_ind_stats, 2, median)
  
  # Calculate the percentage of models in the posterior in which a given parameter is estimated as equal to 0 (% Zero)
  dep_percent_zero <- apply(BT_dep_stats, 2, function(x) sum(x == 0) / length(x))
  ind_percent_zero <- apply(BT_ind_stats, 2, function(x) sum(x == 0) / length(x))
  
  # Select BayesTraits parameters to apply the gelman test to - this section allows the removal of parameters that are estimated to all have a value of 0
  BT_dep_1st_gel_mcmc <- BT_dep_1st_mcmc[, dep_columns]
  BT_dep_2nd_gel_mcmc <- BT_dep_2nd_mcmc[, dep_columns]
  
  BT_ind_1st_gel_mcmc <- BT_ind_1st_mcmc[, ind_columns]
  BT_ind_2nd_gel_mcmc <- BT_ind_2nd_mcmc[, ind_columns]
  
  # Test for convergence between chains
  dep_conv_test <- mcmc.list(list(BT_dep_1st_gel_mcmc, BT_dep_2nd_gel_mcmc))  # May need to select columns to analyse if estimates for a particular parameter (e.g. q12) are all 0 values for each iteration of the chain for both chains.
  ind_conv_test <- mcmc.list(list(BT_ind_1st_gel_mcmc, BT_ind_2nd_gel_mcmc))
  
  tryCatch({
    # Gelman-Rubin convergence diagnostic
    dep_gelman_diag <- gelman.diag(dep_conv_test)
    ind_gelman_diag <- gelman.diag(ind_conv_test)
  }, error = function(e) {
    print("Reduce parameters in gelman.diag")
    # Display heads of MCMC objects
    print("BT_dep_1st_mcmc:")
    print(head(BT_dep_1st_gel_mcmc))
    print("BT_dep_2nd_mcmc:")
    print(head(BT_dep_2nd_gel_mcmc))
    print("BT_ind_1st_mcmc:")
    print(head(BT_ind_1st_gel_mcmc))
    print("BT_ind_2nd_mcmc:")
    print(head(BT_ind_2nd_gel_mcmc))
  })
  
  # Store the results in a list
  results_list <- list(
    dep_tree_count = dep_tree_count,
    ind_tree_count = ind_tree_count,
    dep_ESS = dep_ESS,
    ind_ESS = ind_ESS,
    dep_autocorr = dep_autocorr,
    ind_autocorr = ind_autocorr,
    dep_hpd_interval = dep_hpd_interval,
    ind_hpd_interval = ind_hpd_interval,
    dep_mean = dep_mean,
    dep_median = dep_median,
    ind_mean = ind_mean,
    ind_median = ind_median,
    dep_percent_zero = dep_percent_zero,
    ind_percent_zero = ind_percent_zero,
    dep_gelman_diag = dep_gelman_diag,
    ind_gelman_diag = ind_gelman_diag
  )
  
  return(results_list)
}

bayes_factors <- function(ind_stone_path, dep_stone_path) {
  # Read in stepping stones files
  BT_dep_stone <- read.csv(dep_stone_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  BT_ind_stone <- read.csv(ind_stone_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Calculate Bayes Factors for pairs of trees: independent vs dependent
  BT_dep_stone_select <- as.numeric(BT_dep_stone[nrow(BT_dep_stone), ])
  BT_ind_stone_select <- as.numeric(BT_ind_stone[nrow(BT_ind_stone), ])
  
  bayes_factors <- 2 * (BT_dep_stone_select - BT_ind_stone_select)
  
  return(bayes_factors)
}

##########################################################################################
#Repeat analyses with more iterations: MFpolice_Caste (done), MFcat_fac_ob_Caste (done), MFcat_ob_PGcat_fac (done), MFcat_fac_ob_PGcat_fac (done)
#Don't run (not enough variation in traits - 2 species with obligate polygyny): MFcat_ob_PGcat_ob, MFcat_ob_PGcat_fac_ob

##Run again: MFpolice models (done), CS_median_Caste (gamma models - done)

#MF_Caste
##########################################################################################

###############
#MFpolice_Caste
###############
MFpolice_Caste_result <- bayestraits_statistics("MF_Caste/Police_thresh/Exp/Independent/1st_run/ant_data_BayesTraits_MF_police_thresh_Caste_reduced.txt.Log.txt",
                       "MF_Caste/Police_thresh/Exp/Independent/2nd_run/ant_data_BayesTraits_MF_police_thresh_Caste_reduced.txt.Log.txt",
                       "MF_Caste/Police_thresh/Exp/Dependent/1st_run/ant_data_BayesTraits_MF_police_thresh_Caste_reduced.txt.Log.txt",
                       "MF_Caste/Police_thresh/Exp/Dependent/2nd_run/ant_data_BayesTraits_MF_police_thresh_Caste_reduced.txt.Log.txt")
MFpolice_Caste_result

MFpolice_Caste_BFresult <- bayes_factors("MF_Caste/Police_thresh/Exp/Independent/1st_run/ant_data_BayesTraits_MF_police_thresh_Caste.txt.Stones.txt",
              "MF_Caste/Police_thresh/Exp/Dependent/1st_run/ant_data_BayesTraits_MF_police_thresh_Caste.txt.Stones.txt")
MFpolice_Caste_BFresult

#Create csv
MFpolice_Caste_result_df <- data.frame(Rate = names(MFpolice_Caste_result$ind_median), Transition = c("Empty", "Low mating to high mating", "High mating to low mating", "Single caste to multiple castes", "Multiple castes to single caste"), Median = MFpolice_Caste_result$ind_median, Mean = MFpolice_Caste_result$ind_mean, `95% HPD interval` = apply(round(MFpolice_Caste_result$ind_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFpolice_Caste_result$ind_percent_zero, ESS = MFpolice_Caste_result$ind_ESS)
MFpolice_Caste_result_csv <- MFpolice_Caste_result_df[-1, ] %>% mutate_at(c(3,4,6,7), ~ round(., digits = 2))
write.csv(MFpolice_Caste_result_csv, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/BayesTraits/Diagnostics_parameter_estimation/Results_csv/MFpolice_Caste_result.csv")

#Gamma
MFpolice_Caste_result_gamma <- bayestraits_statistics("MF_Caste/Police_thresh/Gamma/Independent/1st_run/ant_data_BayesTraits_MF_police_thresh_Caste_reduced.txt.Log.txt",
                                                "MF_Caste/Police_thresh/Gamma/Independent/2nd_run/ant_data_BayesTraits_MF_police_thresh_Caste_reduced.txt.Log.txt",
                                                "MF_Caste/Police_thresh/Gamma/Dependent/1st_run/ant_data_BayesTraits_MF_police_thresh_Caste_reduced.txt.Log.txt",
                                                "MF_Caste/Police_thresh/Gamma/Dependent/2nd_run/ant_data_BayesTraits_MF_police_thresh_Caste_reduced.txt.Log.txt")
MFpolice_Caste_result_gamma

MFpolice_Caste_BFresult_gamma <- bayes_factors("MF_Caste/Police_thresh/Gamma/Independent/1st_run/ant_data_BayesTraits_MF_police_thresh_Caste.txt.Stones.txt",
                                         "MF_Caste/Police_thresh/Gamma/Dependent/1st_run/ant_data_BayesTraits_MF_police_thresh_Caste.txt.Stones.txt")
MFpolice_Caste_BFresult_gamma

MFpolice_Caste_result$ind_median
MFpolice_Caste_result_gamma$ind_median

#Create csv
MFpolice_Caste_result_gamma_df <- data.frame(Rate = names(MFpolice_Caste_result_gamma$ind_median), Transition = c("Empty", "Low mating to high mating", "High mating to low mating", "Single caste to multiple castes", "Multiple castes to single caste"), Median = MFpolice_Caste_result_gamma$ind_median, Mean = MFpolice_Caste_result_gamma$ind_mean, `95% HPD interval` = apply(round(MFpolice_Caste_result_gamma$ind_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFpolice_Caste_result_gamma$ind_percent_zero, ESS = MFpolice_Caste_result_gamma$ind_ESS)
MFpolice_Caste_result_gamma_csv <- MFpolice_Caste_result_gamma_df[-1, ] %>% mutate_at(c(3,4,6,7), ~ round(., digits = 2))
write.csv(MFpolice_Caste_result_gamma_csv, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/BayesTraits/Diagnostics_parameter_estimation/Results_csv/")

###############
#MFcat_fac_Caste
###############
MFcat_fac_Caste_result <- bayestraits_statistics("MF_Caste/Fac/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_Caste_reduced.txt.Log.txt",
                                                "MF_Caste/Fac/Exp/Independent/2nd_run/ant_data_BayesTraits_MFcat_fac_Caste_reduced.txt.Log.txt",
                                                "MF_Caste/Fac/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_Caste_reduced.txt.Log.txt",
                                                "MF_Caste/Fac/Exp/Dependent/2nd_run/ant_data_BayesTraits_MFcat_fac_Caste_reduced.txt.Log.txt")
MFcat_fac_Caste_result

MFcat_fac_Caste_BFresult <- bayes_factors("MF_Caste/Fac/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_Caste.txt.Stones.txt",
                                         "MF_Caste/Fac/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_Caste.txt.Stones.txt")
MFcat_fac_Caste_BFresult

#Create csv
MFcat_fac_Caste_result_df <- data.frame(Rate = names(MFcat_fac_Caste_result$ind_median), Transition = c("Empty", "Monandrous to facultatively polyandrous", "Facultatively polyandrous to monandrous", "Single caste to multiple castes", "Multiple castes to single caste"), Median = MFcat_fac_Caste_result$ind_median, Mean = MFcat_fac_Caste_result$ind_mean, `95% HPD interval` = apply(round(MFcat_fac_Caste_result$ind_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_fac_Caste_result$ind_percent_zero, ESS = MFcat_fac_Caste_result$ind_ESS)
MFcat_fac_Caste_result_csv <- MFcat_fac_Caste_result_df[-1, ] %>% mutate_at(c(3,4,6,7), ~ round(., digits = 2))
# write.csv(MFcat_fac_Caste_result_csv, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/BayesTraits/BayesTraitsV4.0.1-OSX/Post_review_results_csv/MFcat_fac_Caste_result.csv")


###############
#MFcat_ob_Caste
###############
MFcat_ob_Caste_result <- bayestraits_statistics("MF_Caste/Ob/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_ob_Caste_reduced.txt.Log.txt",
                                                "MF_Caste/Ob/Exp/Independent/2nd_run/ant_data_BayesTraits_MFcat_ob_Caste_reduced.txt.Log.txt",
                                                "MF_Caste/Ob/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_ob_Caste_reduced.txt.Log.txt",
                                                "MF_Caste/Ob/Exp/Dependent/2nd_run/ant_data_BayesTraits_MFcat_ob_Caste_reduced.txt.Log.txt")
MFcat_ob_Caste_result

MFcat_ob_Caste_BFresult <- bayes_factors("MF_Caste/Ob/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_ob_Caste.txt.Stones.txt",
                                         "MF_Caste/Ob/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_ob_Caste.txt.Stones.txt")
MFcat_ob_Caste_BFresult

#Create csv
MFcat_ob_Caste_result_df <- data.frame(Rate = names(MFcat_ob_Caste_result$ind_median), Transition = c("Empty", "Monandrous to obligately polyandrous", "Obligately polyandrous to monandrous", "Single caste to multiple castes", "Multiple castes to single caste"), Median = MFcat_ob_Caste_result$ind_median, Mean = MFcat_ob_Caste_result$ind_mean, `95% HPD interval` = apply(round(MFcat_ob_Caste_result$ind_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_ob_Caste_result$ind_percent_zero, ESS = MFcat_ob_Caste_result$ind_ESS)
MFcat_ob_Caste_result_csv <- MFcat_ob_Caste_result_df[-1, ] %>% mutate_at(c(3,4,6,7), ~ round(., digits = 2))
write.csv(MFcat_ob_Caste_result_csv, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/BayesTraits/Diagnostics_parameter_estimation/Results_csv/MFcat_ob_Caste_result.csv")

###############
#MFcat_fac_ob_Caste
###############
MFcat_fac_ob_Caste_result <- bayestraits_statistics("MF_Caste/Fac_Ob/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_Caste_reduced.txt.Log.txt",
                                                "MF_Caste/Fac_Ob/Exp/Independent/2nd_run/ant_data_BayesTraits_MFcat_fac_ob_Caste_reduced.txt.Log.txt",
                                                "MF_Caste/Fac_Ob/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_Caste_reduced.txt.Log.txt",
                                                "MF_Caste/Fac_Ob/Exp/Dependent/2nd_run/ant_data_BayesTraits_MFcat_fac_ob_Caste_reduced.txt.Log.txt")
MFcat_fac_ob_Caste_result

MFcat_fac_ob_Caste_BFresult <- bayes_factors("MF_Caste/Fac_Ob/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_Caste.txt.Stones.txt",
                                         "MF_Caste/Fac_Ob/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_Caste.txt.Stones.txt")
MFcat_fac_ob_Caste_BFresult

#Create csv
MFcat_fac_ob_Caste_result_df <- data.frame(Rate = names(MFcat_fac_ob_Caste_result$ind_median), Transition = c("Empty", "Monandrous/facultatively polyandrous to obligately polyandrous", "Obligately polyandrous to monandrous/facultatively polyandrous", "Single caste to multiple castes", "Multiple castes to single caste"), Median = MFcat_fac_ob_Caste_result$ind_median, Mean = MFcat_fac_ob_Caste_result$ind_mean, `95% HPD interval` = apply(round(MFcat_fac_ob_Caste_result$ind_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_fac_ob_Caste_result$ind_percent_zero, ESS = MFcat_fac_ob_Caste_result$ind_ESS)
MFcat_fac_ob_Caste_result_csv <- MFcat_fac_ob_Caste_result_df[-1, ] %>% mutate_at(c(3,4,6,7), ~ round(., digits = 2))
# write.csv(MFcat_fac_ob_Caste_result_csv, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/BayesTraits/Diagnostics_parameter_estimation/Results_csv/MFcat_fac_ob_Caste_result.csv")


#CS_Caste
##########################################################################################

###############
#CS_low_Caste
###############
CS_low_Caste_result <- bayestraits_statistics("CS_Caste/Low/Exp/Independent/1st_run/ant_data_BayesTraits_CS_low_thresh_Caste_reduced.txt.Log.txt",
                                                 "CS_Caste/Low/Exp/Independent/2nd_run/ant_data_BayesTraits_CS_low_thresh_Caste_reduced.txt.Log.txt",
                                                 "CS_Caste/Low/Exp/Dependent/1st_run/ant_data_BayesTraits_CS_low_thresh_Caste_reduced.txt.Log.txt",
                                                 "CS_Caste/Low/Exp/Dependent/2nd_run/ant_data_BayesTraits_CS_low_thresh_Caste_reduced.txt.Log.txt")
CS_low_Caste_result

CS_low_Caste_BFresult <- bayes_factors("CS_Caste/Low/Exp/Independent/1st_run/ant_data_BayesTraits_CS_low_thresh_Caste.txt.Stones.txt",
                                          "CS_Caste/Low/Exp/Dependent/1st_run/ant_data_BayesTraits_CS_low_thresh_Caste.txt.Stones.txt")
CS_low_Caste_BFresult

#Create csv
CS_low_Caste_result_df <- data.frame(Rate = names(CS_low_Caste_result$dep_median), `Number of castes` = c("Empty", "Single to multiple", "Single", "Multiple to single", "Multiple", "Single", "Single to multiple", "Multiple", "Multiple to single"), `Colony size` = c("Empty", "Small", "Small to large", "Small", "Small to large", "Large to small", "Large", "Large to small", "Large"), Median = CS_low_Caste_result$dep_median, Mean = CS_low_Caste_result$dep_mean, `95% HPD interval` = apply(round(CS_low_Caste_result$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = CS_low_Caste_result$dep_percent_zero, ESS = CS_low_Caste_result$dep_ESS)
CS_low_Caste_result_csv <- CS_low_Caste_result_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
write.csv(CS_low_Caste_result_csv, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/BayesTraits/Diagnostics_parameter_estimation/Results_csv/CS_low_Caste_result.csv")


###############
#CS_median_Caste
###############
CS_median_Caste_result <- bayestraits_statistics("CS_Caste/Median/Exp/Independent/1st_run/ant_data_BayesTraits_CS_median_thresh_Caste_reduced.txt.Log.txt",
                                              "CS_Caste/Median/Exp/Independent/2nd_run/ant_data_BayesTraits_CS_median_thresh_Caste_reduced.txt.Log.txt",
                                              "CS_Caste/Median/Exp/Dependent/1st_run/ant_data_BayesTraits_CS_median_thresh_Caste_reduced.txt.Log.txt",
                                              "CS_Caste/Median/Exp/Dependent/2nd_run/ant_data_BayesTraits_CS_median_thresh_Caste_reduced.txt.Log.txt", dep_columns = c("q13", "q21", "q24", "q31", "q34", "q42", "q43"))
CS_median_Caste_result

CS_median_Caste_BFresult <- bayes_factors("CS_Caste/Median/Exp/Independent/1st_run/ant_data_BayesTraits_CS_median_thresh_Caste.txt.Stones.txt",
                                       "CS_Caste/Median/Exp/Dependent/1st_run/ant_data_BayesTraits_CS_median_thresh_Caste.txt.Stones.txt")
CS_median_Caste_BFresult

#Create csv
CS_median_Caste_result_df <- data.frame(Rate = names(CS_median_Caste_result$dep_median), `Number of castes` = c("Empty", "Single to multiple", "Single", "Multiple to single", "Multiple", "Single", "Single to multiple", "Multiple", "Multiple to single"), `Colony size` = c("Empty", "Small", "Small to large", "Small", "Small to large", "Large to small", "Large", "Large to small", "Large"), Median = CS_median_Caste_result$dep_median, Mean = CS_median_Caste_result$dep_mean, `95% HPD interval` = apply(round(CS_median_Caste_result$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = CS_median_Caste_result$dep_percent_zero, ESS = CS_median_Caste_result$dep_ESS)
CS_median_Caste_result_csv <- CS_median_Caste_result_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
# write.csv(CS_median_Caste_result_csv, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/BayesTraits/Diagnostics_parameter_estimation/Results_csv/CS_median_Caste_result.csv")

#Gamma
CS_median_Caste_result_gamma <- bayestraits_statistics("CS_Caste/Median/Gamma/Independent/1st_run/ant_data_BayesTraits_CS_median_thresh_Caste_reduced.txt.Log.txt",
                                                 "CS_Caste/Median/Gamma/Independent/2nd_run/ant_data_BayesTraits_CS_median_thresh_Caste_reduced.txt.Log.txt",
                                                 "CS_Caste/Median/Gamma/Dependent/1st_run/ant_data_BayesTraits_CS_median_thresh_Caste_reduced.txt.Log.txt",
                                                 "CS_Caste/Median/Gamma/Dependent/2nd_run/ant_data_BayesTraits_CS_median_thresh_Caste_reduced.txt.Log.txt", dep_columns = c("q13", "q21", "q24", "q31", "q34", "q42", "q43"))
CS_median_Caste_result_gamma

CS_median_Caste_BFresult_gamma <- bayes_factors("CS_Caste/Median/Gamma/Independent/1st_run/ant_data_BayesTraits_CS_median_thresh_Caste.txt.Stones.txt",
                                          "CS_Caste/Median/Gamma/Dependent/1st_run/ant_data_BayesTraits_CS_median_thresh_Caste.txt.Stones.txt")
CS_median_Caste_BFresult_gamma

CS_median_Caste_result$dep_median
CS_median_Caste_result_gamma$dep_median

#Create csv
CS_median_Caste_result_gamma_df <- data.frame(Rate = names(CS_median_Caste_result_gamma$dep_median), `Number of castes` = c("Empty", "Single to multiple", "Single", "Multiple to single", "Multiple", "Single", "Single to multiple", "Multiple", "Multiple to single"), `Colony size` = c("Empty", "Small", "Small to large", "Small", "Small to large", "Large to small", "Large", "Large to small", "Large"), Median = CS_median_Caste_result_gamma$dep_median, Mean = CS_median_Caste_result_gamma$dep_mean, `95% HPD interval` = apply(round(CS_median_Caste_result_gamma$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = CS_median_Caste_result_gamma$dep_percent_zero, ESS = CS_median_Caste_result_gamma$dep_ESS)
CS_median_Caste_result_gamma_csv <- CS_median_Caste_result_gamma_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
# write.csv(CS_median_Caste_result_gamma_csv, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/BayesTraits/Diagnostics_parameter_estimation/Results_csv/CS_median_Caste_result_gamma.csv")

###############
#CS_high_Caste
###############
CS_high_Caste_result <- bayestraits_statistics("CS_Caste/High/Exp/Independent/1st_run/ant_data_BayesTraits_CS_high_thresh_Caste_reduced.txt.Log.txt",
                                              "CS_Caste/High/Exp/Independent/2nd_run/ant_data_BayesTraits_CS_high_thresh_Caste_reduced.txt.Log.txt",
                                              "CS_Caste/High/Exp/Dependent/1st_run/ant_data_BayesTraits_CS_high_thresh_Caste_reduced.txt.Log.txt",
                                              "CS_Caste/High/Exp/Dependent/2nd_run/ant_data_BayesTraits_CS_high_thresh_Caste_reduced.txt.Log.txt",
                                              dep_columns = c("q13", "q21", "q24", "q31", "q34", "q42", "q43"))
CS_high_Caste_result

CS_high_Caste_BFresult <- bayes_factors("CS_Caste/High/Exp/Independent/1st_run/ant_data_BayesTraits_CS_high_thresh_Caste.txt.Stones.txt",
                                       "CS_Caste/High/Exp/Dependent/1st_run/ant_data_BayesTraits_CS_high_thresh_Caste.txt.Stones.txt")
CS_high_Caste_BFresult

#Create csv
CS_high_Caste_result_df <- data.frame(Rate = names(CS_high_Caste_result$dep_median), `Number of castes` = c("Empty", "Single to multiple", "Single", "Multiple to single", "Multiple", "Single", "Single to multiple", "Multiple", "Multiple to single"), `Colony size` = c("Empty", "Small", "Small to large", "Small", "Small to large", "Large to small", "Large", "Large to small", "Large"), Median = CS_high_Caste_result$dep_median, Mean = CS_high_Caste_result$dep_mean, `95% HPD interval` = apply(round(CS_high_Caste_result$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = CS_high_Caste_result$dep_percent_zero, ESS = CS_high_Caste_result$dep_ESS)
CS_high_Caste_result_csv <- CS_high_Caste_result_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
# write.csv(CS_high_Caste_result_csv, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/BayesTraits/Diagnostics_parameter_estimation/Results_csv/CS_high_Caste_result.csv")



#MFcat_PGcat
##########################################################################################

###############
#MFcat_fac_PGcat_fac
###############
MFcat_fac_PGcat_fac_result <- bayestraits_statistics("MF_PG/MF_fac_PG_fac/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_reduced.txt.Log.txt",
                                              "MF_PG/MF_fac_PG_fac/Exp/Independent/2nd_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_reduced.txt.Log.txt",
                                              "MF_PG/MF_fac_PG_fac/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_reduced.txt.Log.txt",
                                              "MF_PG/MF_fac_PG_fac/Exp/Dependent/2nd_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_reduced.txt.Log.txt")
MFcat_fac_PGcat_fac_result

MFcat_fac_PGcat_fac_BFresult <- bayes_factors("MF_PG/MF_fac_PG_fac/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac.txt.Stones.txt",
                                       "MF_PG/MF_fac_PG_fac/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac.txt.Stones.txt")
MFcat_fac_PGcat_fac_BFresult

#Create csv
MFcat_fac_PGcat_fac_result_df <- data.frame(Rate = names(MFcat_fac_PGcat_fac_result$dep_median), `Number of queens` = c("Empty", "Monogynous to facultatively polygynous", "Monogynous", "Facultatively polygynous to monogynous", "Facultatively polygynous", "Monogynous", "Monogynous to facultatively polygynous", "Facultatively polygynous", "Facultatively polygynous to monogynous"), `Queen mating frequency` = c("Empty", "Monandrous", "Monandrous to facultatively polyandrous", "Monandrous", "Monandrous to facultatively polyandrous", "Facultatively polyandrous to monandrous", "Facultatively polyandrous", "Facultatively polyandrous to monandrous", "Facultatively polyandrous"), Median = MFcat_fac_PGcat_fac_result$dep_median, Mean = MFcat_fac_PGcat_fac_result$dep_mean, `95% HPD interval` = apply(round(MFcat_fac_PGcat_fac_result$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_fac_PGcat_fac_result$dep_percent_zero, ESS = MFcat_fac_PGcat_fac_result$dep_ESS)
MFcat_fac_PGcat_fac_result_csv <- MFcat_fac_PGcat_fac_result_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
# write.csv(MFcat_fac_PGcat_fac_result_csv, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/BayesTraits/BayesTraitsV4.0.1-OSX/Post_review_results_csv/MFcat_fac_PGcat_fac_result.csv")

#Gamma
MFcat_fac_PGcat_fac_result_gamma <- bayestraits_statistics("MF_PG/MF_fac_PG_fac/Gamma/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_reduced.txt.Log.txt",
                                                     "MF_PG/MF_fac_PG_fac/Gamma/Independent/2nd_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_reduced.txt.Log.txt",
                                                     "MF_PG/MF_fac_PG_fac/Gamma/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_reduced.txt.Log.txt",
                                                     "MF_PG/MF_fac_PG_fac/Gamma/Dependent/2nd_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_reduced.txt.Log.txt")
MFcat_fac_PGcat_fac_result_gamma

MFcat_fac_PGcat_fac_BFresult_gamma <- bayes_factors("MF_PG/MF_fac_PG_fac/Gamma/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac.txt.Stones.txt",
                                              "MF_PG/MF_fac_PG_fac/Gamma/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac.txt.Stones.txt")
MFcat_fac_PGcat_fac_BFresult_gamma

MFcat_fac_PGcat_fac_result$dep_median
MFcat_fac_PGcat_fac_result_gamma$dep_median

#Create csv
MFcat_fac_PGcat_fac_result_gamma_df <- data.frame(Rate = names(MFcat_fac_PGcat_fac_result_gamma$dep_median), `Number of queens` = c("Empty", "Monogynous to facultatively polygynous", "Monogynous", "Facultatively polygynous to monogynous", "Facultatively polygynous", "Monogynous", "Monogynous to facultatively polygynous", "Facultatively polygynous", "Facultatively polygynous to monogynous"), `Queen mating frequency` = c("Empty", "Monandrous", "Monandrous to facultatively polyandrous", "Monandrous", "Monandrous to facultatively polyandrous", "Facultatively polyandrous to monandrous", "Facultatively polyandrous", "Facultatively polyandrous to monandrous", "Facultatively polyandrous"), Median = MFcat_fac_PGcat_fac_result_gamma$dep_median, Mean = MFcat_fac_PGcat_fac_result_gamma$dep_mean, `95% HPD interval` = apply(round(MFcat_fac_PGcat_fac_result_gamma$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_fac_PGcat_fac_result_gamma$dep_percent_zero, ESS = MFcat_fac_PGcat_fac_result_gamma$dep_ESS)
MFcat_fac_PGcat_fac_result_gamma_csv <- MFcat_fac_PGcat_fac_result_gamma_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
# write.csv(MFcat_fac_PGcat_fac_result_gamma_csv, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/BayesTraits/BayesTraitsV4.0.1-OSX/Post_review_result_gammas_csv/MFcat_fac_PGcat_fac_result_gamma.csv")

###############
#MFcat_ob_PGcat_fac
###############
MFcat_ob_PGcat_fac_result <- bayestraits_statistics("MF_PG/MF_ob_PG_fac/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_ob_PGcat_fac_reduced.txt.Log.txt",
                                                     "MF_PG/MF_ob_PG_fac/Exp/Independent/2nd_run/ant_data_BayesTraits_MFcat_ob_PGcat_fac_reduced.txt.Log.txt",
                                                     "MF_PG/MF_ob_PG_fac/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_ob_PGcat_fac_reduced.txt.Log.txt",
                                                     "MF_PG/MF_ob_PG_fac/Exp/Dependent/2nd_run/ant_data_BayesTraits_MFcat_ob_PGcat_fac_reduced.txt.Log.txt")
MFcat_ob_PGcat_fac_result

MFcat_ob_PGcat_fac_BFresult <- bayes_factors("MF_PG/MF_ob_PG_fac/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_ob_PGcat_fac.txt.Stones.txt",
                                              "MF_PG/MF_ob_PG_fac/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_ob_PGcat_fac.txt.Stones.txt")
MFcat_ob_PGcat_fac_BFresult

#Create csv
MFcat_ob_PGcat_fac_result_df <- data.frame(Rate = names(MFcat_ob_PGcat_fac_result$dep_median), `Number of queens` = c("Empty", "Monogynous to facultatively polygynous", "Monogynous", "Facultatively polygynous to monogynous", "Facultatively polygynous", "Monogynous", "Monogynous to facultatively polygynous", "Facultatively polygynous", "Facultatively polygynous to monogynous"), `Queen mating frequency` = c("Empty", "Monandrous", "Monandrous to obligately polyandrous", "Monandrous", "Monandrous to obligately polyandrous", "Obligately polyandrous to monandrous", "Obligately polyandrous", "Obligately polyandrous to monandrous", "Obligately polyandrous"), Median = MFcat_ob_PGcat_fac_result$dep_median, Mean = MFcat_ob_PGcat_fac_result$dep_mean, `95% HPD interval` = apply(round(MFcat_ob_PGcat_fac_result$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_ob_PGcat_fac_result$dep_percent_zero, ESS = MFcat_ob_PGcat_fac_result$dep_ESS)
MFcat_ob_PGcat_fac_result_csv <- MFcat_ob_PGcat_fac_result_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
# write.csv(MFcat_ob_PGcat_fac_result_csv, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/BayesTraits/BayesTraitsV4.0.1-OSX/Post_review_results_csv/MFcat_ob_PGcat_fac_result_csv.csv")

#Gamma
MFcat_ob_PGcat_fac_result_gamma <- bayestraits_statistics("MF_PG/MF_ob_PG_fac/Gamma/Independent/1st_run/ant_data_BayesTraits_MFcat_ob_PGcat_fac_reduced.txt.Log.txt",
                                                    "MF_PG/MF_ob_PG_fac/Gamma/Independent/2nd_run/ant_data_BayesTraits_MFcat_ob_PGcat_fac_reduced.txt.Log.txt",
                                                    "MF_PG/MF_ob_PG_fac/Gamma/Dependent/1st_run/ant_data_BayesTraits_MFcat_ob_PGcat_fac_reduced.txt.Log.txt",
                                                    "MF_PG/MF_ob_PG_fac/Gamma/Dependent/2nd_run/ant_data_BayesTraits_MFcat_ob_PGcat_fac_reduced.txt.Log.txt")
MFcat_ob_PGcat_fac_result_gamma

MFcat_ob_PGcat_fac_BFresult_gamma <- bayes_factors("MF_PG/MF_ob_PG_fac/Gamma/Independent/1st_run/ant_data_BayesTraits_MFcat_ob_PGcat_fac.txt.Stones.txt",
                                             "MF_PG/MF_ob_PG_fac/Gamma/Dependent/1st_run/ant_data_BayesTraits_MFcat_ob_PGcat_fac.txt.Stones.txt")
MFcat_ob_PGcat_fac_BFresult_gamma

MFcat_ob_PGcat_fac_result$dep_median
MFcat_ob_PGcat_fac_result_gamma$dep_median

#Create csv
MFcat_ob_PGcat_fac_result_gamma_df <- data.frame(Rate = names(MFcat_ob_PGcat_fac_result_gamma$dep_median), `Number of queens` = c("Empty", "Monogynous to facultatively polygynous", "Monogynous", "Facultatively polygynous to monogynous", "Facultatively polygynous", "Monogynous", "Monogynous to facultatively polygynous", "Facultatively polygynous", "Facultatively polygynous to monogynous"), `Queen mating frequency` = c("Empty", "Monandrous", "Monandrous to obligately polyandrous", "Monandrous", "Monandrous to obligately polyandrous", "Obligately polyandrous to monandrous", "Obligately polyandrous", "Obligately polyandrous to monandrous", "Obligately polyandrous"), Median = MFcat_ob_PGcat_fac_result_gamma$dep_median, Mean = MFcat_ob_PGcat_fac_result_gamma$dep_mean, `95% HPD interval` = apply(round(MFcat_ob_PGcat_fac_result_gamma$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_ob_PGcat_fac_result_gamma$dep_percent_zero, ESS = MFcat_ob_PGcat_fac_result_gamma$dep_ESS)
MFcat_ob_PGcat_fac_result_gamma_csv <- MFcat_ob_PGcat_fac_result_gamma_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
# write.csv(MFcat_ob_PGcat_fac_result_gamma_csv, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/BayesTraits/BayesTraitsV4.0.1-OSX/Post_review_result_gammas_csv/MFcat_ob_PGcat_fac_result_gamma_csv.csv")


###############
#MFcat_fac_ob_PGcat_fac
###############
MFcat_fac_ob_PGcat_fac_result <- bayestraits_statistics("MF_PG/MF_fac_ob_PG_fac/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_reduced.txt.Log.txt",
                                                     "MF_PG/MF_fac_ob_PG_fac/Exp/Independent/2nd_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_reduced.txt.Log.txt",
                                                     "MF_PG/MF_fac_ob_PG_fac/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_reduced.txt.Log.txt",
                                                     "MF_PG/MF_fac_ob_PG_fac/Exp/Dependent/2nd_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_reduced.txt.Log.txt")
MFcat_fac_ob_PGcat_fac_result

MFcat_fac_ob_PGcat_fac_BFresult <- bayes_factors("MF_PG/MF_fac_ob_PG_fac/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac.txt.Stones.txt",
                                              "MF_PG/MF_fac_ob_PG_fac/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac.txt.Stones.txt")
MFcat_fac_ob_PGcat_fac_BFresult

MFcat_fac_ob_PGcat_fac_result_df <- data.frame(Rate = names(MFcat_fac_ob_PGcat_fac_result$dep_median), `Number of queens` = c("Empty", "Monogynous to facultatively polygynous", "Monogynous", "Facultatively polygynous to monogynous", "Facultatively polygynous", "Monogynous", "Monogynous to facultatively polygynous", "Facultatively polygynous", "Facultatively polygynous to monogynous"), `Queen mating frequency` = c("Empty", "Monandrous/facultatively polyandrous", "Monandrous/facultatively polyandrous to obligately polyandrous", "Monandrous/facultatively polyandrous", "Monandrous/facultatively polyandrous to obligately polyandrous", "Obligately polyandrous to monandrous/facultatively polyandrous", "Obligately polyandrous", "Obligately polyandrous to monandrous/facultatively polyandrous", "Obligately polyandrous"), Median = MFcat_fac_ob_PGcat_fac_result$dep_median, Mean = MFcat_fac_ob_PGcat_fac_result$dep_mean, `95% HPD interval` = apply(round(MFcat_fac_ob_PGcat_fac_result$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_fac_ob_PGcat_fac_result$dep_percent_zero, ESS = MFcat_fac_ob_PGcat_fac_result$dep_ESS)
MFcat_fac_ob_PGcat_fac_result_csv <- MFcat_fac_ob_PGcat_fac_result_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
# write.csv(MFcat_fac_ob_PGcat_fac_result_csv, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/BayesTraits/BayesTraitsV4.0.1-OSX/Post_review_results_csv/MFcat_fac_ob_PGcat_fac_result_csv.csv")

###############
#MFcat_fac_PGcat_ob
###############
MFcat_fac_PGcat_ob_result <- bayestraits_statistics("MF_PG/MF_fac_PG_ob/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_ob_reduced.txt.Log.txt",
                                                     "MF_PG/MF_fac_PG_ob/Exp/Independent/2nd_run/ant_data_BayesTraits_MFcat_fac_PGcat_ob_reduced.txt.Log.txt",
                                                     "MF_PG/MF_fac_PG_ob/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_ob_reduced.txt.Log.txt",
                                                     "MF_PG/MF_fac_PG_ob/Exp/Dependent/2nd_run/ant_data_BayesTraits_MFcat_fac_PGcat_ob_reduced.txt.Log.txt")
MFcat_fac_PGcat_ob_result

MFcat_fac_PGcat_ob_BFresult <- bayes_factors("MF_PG/MF_fac_PG_ob/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_ob.txt.Stones.txt",
                                              "MF_PG/MF_fac_PG_ob/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_ob.txt.Stones.txt")
MFcat_fac_PGcat_ob_BFresult

#Create csv
MFcat_fac_PGcat_ob_result_df <- data.frame(Rate = names(MFcat_fac_PGcat_ob_result$dep_median), `Number of queens` = c("Empty", "Monogynous to obligately polygynous", "Monogynous", "Obligately polygynous to monogynous", "Obligately polygynous", "Monogynous", "Monogynous to obligately polygynous", "Obligately polygynous", "Obligately polygynous to monogynous"), `Queen mating frequency` = c("Empty", "Monandrous", "Monandrous to facultatively polyandrous", "Monandrous", "Monandrous to facultatively polyandrous", "Facultatively polyandrous to monandrous", "Facultatively polyandrous", "Facultatively polyandrous to monandrous", "Facultatively polyandrous"), Median = MFcat_fac_PGcat_ob_result$dep_median, Mean = MFcat_fac_PGcat_ob_result$dep_mean, `95% HPD interval` = apply(round(MFcat_fac_PGcat_ob_result$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_fac_PGcat_ob_result$dep_percent_zero, ESS = MFcat_fac_PGcat_ob_result$dep_ESS)
MFcat_fac_PGcat_ob_result_csv <- MFcat_fac_PGcat_ob_result_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
# write.csv(MFcat_fac_PGcat_ob_result_csv, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/BayesTraits/BayesTraitsV4.0.1-OSX/Post_review_results_csv/MFcat_fac_PGcat_ob_result.csv")

#Gamma
MFcat_fac_PGcat_ob_result_gamma <- bayestraits_statistics("MF_PG/MF_fac_PG_ob/Gamma/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_ob_reduced.txt.Log.txt",
                                                    "MF_PG/MF_fac_PG_ob/Gamma/Independent/2nd_run/ant_data_BayesTraits_MFcat_fac_PGcat_ob_reduced.txt.Log.txt",
                                                    "MF_PG/MF_fac_PG_ob/Gamma/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_ob_reduced.txt.Log.txt",
                                                    "MF_PG/MF_fac_PG_ob/Gamma/Dependent/2nd_run/ant_data_BayesTraits_MFcat_fac_PGcat_ob_reduced.txt.Log.txt")
MFcat_fac_PGcat_ob_result_gamma

MFcat_fac_PGcat_ob_BFresult_gamma <- bayes_factors("MF_PG/MF_fac_PG_ob/Gamma/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_ob.txt.Stones.txt",
                                             "MF_PG/MF_fac_PG_ob/Gamma/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_ob.txt.Stones.txt")
MFcat_fac_PGcat_ob_BFresult_gamma

MFcat_fac_PGcat_ob_result$dep_median
MFcat_fac_PGcat_ob_result_gamma$dep_median

#Create csv
MFcat_fac_PGcat_ob_result_gamma_df <- data.frame(Rate = names(MFcat_fac_PGcat_ob_result_gamma$dep_median), `Number of queens` = c("Empty", "Monogynous to obligately polygynous", "Monogynous", "Obligately polygynous to monogynous", "Obligately polygynous", "Monogynous", "Monogynous to obligately polygynous", "Obligately polygynous", "Obligately polygynous to monogynous"), `Queen mating frequency` = c("Empty", "Monandrous", "Monandrous to facultatively polyandrous", "Monandrous", "Monandrous to facultatively polyandrous", "Facultatively polyandrous to monandrous", "Facultatively polyandrous", "Facultatively polyandrous to monandrous", "Facultatively polyandrous"), Median = MFcat_fac_PGcat_ob_result_gamma$dep_median, Mean = MFcat_fac_PGcat_ob_result_gamma$dep_mean, `95% HPD interval` = apply(round(MFcat_fac_PGcat_ob_result_gamma$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_fac_PGcat_ob_result_gamma$dep_percent_zero, ESS = MFcat_fac_PGcat_ob_result_gamma$dep_ESS)
MFcat_fac_PGcat_ob_result_gamma_csv <- MFcat_fac_PGcat_ob_result_gamma_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
# write.csv(MFcat_fac_PGcat_ob_result_gamma_csv, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/BayesTraits/BayesTraitsV4.0.1-OSX/Post_review_result_gammas_csv/MFcat_fac_PGcat_ob_result_gamma.csv")

###############
#MFcat_fac_ob_PGcat_ob
###############
MFcat_fac_ob_PGcat_ob_result <- bayestraits_statistics("MF_PG/MF_fac_ob_PG_ob/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_ob_reduced.txt.Log.txt",
                                                    "MF_PG/MF_fac_ob_PG_ob/Exp/Independent/2nd_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_ob_reduced.txt.Log.txt",
                                                    "MF_PG/MF_fac_ob_PG_ob/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_ob_reduced.txt.Log.txt",
                                                    "MF_PG/MF_fac_ob_PG_ob/Exp/Dependent/2nd_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_ob_reduced.txt.Log.txt")
MFcat_fac_ob_PGcat_ob_result

MFcat_fac_ob_PGcat_ob_BFresult <- bayes_factors("MF_PG/MF_fac_ob_PG_ob/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_ob.txt.Stones.txt",
                                             "MF_PG/MF_fac_ob_PG_ob/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_ob.txt.Stones.txt")
MFcat_fac_ob_PGcat_ob_BFresult

#Create csv
MFcat_fac_ob_PGcat_ob_result_df <- data.frame(Rate = names(MFcat_fac_ob_PGcat_ob_result$dep_median), `Number of queens` = c("Empty", "Monogynous to obligately polygynous", "Monogynous", "Obligately polygynous to monogynous", "Obligately polygynous", "Monogynous", "Monogynous to obligately polygynous", "Obligately polygynous", "Obligately polygynous to monogynous"), `Queen mating frequency` = c("Empty", "Monandrous/facultatively polyandrous", "Monandrous/facultatively polyandrous to obligately polyandrous", "Monandrous/facultatively polyandrous", "Monandrous/facultatively polyandrous to obligately polyandrous", "Obligately polyandrous to monandrous/facultatively polyandrous", "Obligately polyandrous", "Obligately polyandrous to monandrous/facultatively polyandrous", "Obligately polyandrous"), Median = MFcat_fac_ob_PGcat_ob_result$dep_median, Mean = MFcat_fac_ob_PGcat_ob_result$dep_mean, `95% HPD interval` = apply(round(MFcat_fac_ob_PGcat_ob_result$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_fac_ob_PGcat_ob_result$dep_percent_zero, ESS = MFcat_fac_ob_PGcat_ob_result$dep_ESS)
MFcat_fac_ob_PGcat_ob_result_csv <- MFcat_fac_ob_PGcat_ob_result_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
write.csv(MFcat_fac_ob_PGcat_ob_result_csv, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/BayesTraits/Diagnostics_parameter_estimation/Results_csv/MFcat_fac_ob_PGcat_ob_result.csv")

###############
#MFcat_fac_PGcat_fac_ob
###############
MFcat_fac_PGcat_fac_ob_result <- bayestraits_statistics("MF_PG/MF_fac_PG_fac_ob/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_ob_reduced.txt.Log.txt",
                                                       "MF_PG/MF_fac_PG_fac_ob/Exp/Independent/2nd_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_ob_reduced.txt.Log.txt",
                                                       "MF_PG/MF_fac_PG_fac_ob/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_ob_reduced.txt.Log.txt",
                                                       "MF_PG/MF_fac_PG_fac_ob/Exp/Dependent/2nd_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_ob_reduced.txt.Log.txt")
MFcat_fac_PGcat_fac_ob_result

MFcat_fac_PGcat_fac_ob_BFresult <- bayes_factors("MF_PG/MF_fac_PG_fac_ob/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_ob.txt.Stones.txt",
                                                "MF_PG/MF_fac_PG_fac_ob/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_PGcat_fac_ob.txt.Stones.txt")
MFcat_fac_PGcat_fac_ob_BFresult

#Create csv
MFcat_fac_PGcat_fac_ob_result_df <- data.frame(Rate = names(MFcat_fac_PGcat_fac_ob_result$ind_median), Transition = c("Empty", "Monandrous to facultatively polyandrous", "Facultatively polyandrous to monandrous", "Monogynous/facultatively polygynous to obligately polygynous", "Obligately polygynous to monogynous/facultatively polygynous"), Median = MFcat_fac_PGcat_fac_ob_result$ind_median, Mean = MFcat_fac_PGcat_fac_ob_result$ind_mean, `95% HPD interval` = apply(round(MFcat_fac_PGcat_fac_ob_result$ind_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_fac_PGcat_fac_ob_result$ind_percent_zero, ESS = MFcat_fac_PGcat_fac_ob_result$ind_ESS)
MFcat_fac_PGcat_fac_ob_result_csv <- MFcat_fac_PGcat_fac_ob_result_df[-1, ] %>% mutate_at(c(3,4,6,7), ~ round(., digits = 2))
# write.csv(MFcat_fac_PGcat_fac_ob_result_csv, file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/BayesTraits/BayesTraitsV4.0.1-OSX/Post_review_results_csv/MFcat_fac_PGcat_fac_ob_result.csv")


###############
#MFcat_fac_ob_PGcat_fac_ob
###############
MFcat_fac_ob_PGcat_fac_ob_result <- bayestraits_statistics("MF_PG/MF_fac_ob_PG_fac_ob/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_ob_reduced.txt.Log.txt",
                                                        "MF_PG/MF_fac_ob_PG_fac_ob/Exp/Independent/2nd_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_ob_reduced.txt.Log.txt",
                                                        "MF_PG/MF_fac_ob_PG_fac_ob/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_ob_reduced.txt.Log.txt",
                                                        "MF_PG/MF_fac_ob_PG_fac_ob/Exp/Dependent/2nd_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_ob_reduced.txt.Log.txt")
MFcat_fac_ob_PGcat_fac_ob_result

MFcat_fac_ob_PGcat_fac_ob_BFresult <- bayes_factors("MF_PG/MF_fac_ob_PG_fac_ob/Exp/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_ob.txt.Stones.txt",
                                                 "MF_PG/MF_fac_ob_PG_fac_ob/Exp/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_ob.txt.Stones.txt")
MFcat_fac_ob_PGcat_fac_ob_BFresult

#Create csv
MFcat_fac_ob_PGcat_fac_ob_result_df <- data.frame(Rate = names(MFcat_fac_ob_PGcat_fac_ob_result$dep_median), `Number of queens` = c("Empty", "Monogynous/facultatively polygynous to obligately polygynous", "Monogynous/facultatively polygynous", "Obligately polygynous to monogynous/facultatively polygynous", "Obligately polygynous", "Monogynous/facultatively polygynous", "Monogynous/facultatively polygynous to obligately polygynous", "Obligately polygynous", "Obligately polygynous to monogynous/facultatively polygynous"), `Queen mating frequency` = c("Empty", "Monandrous/facultatively polyandrous", "Monandrous/facultatively polyandrous to obligately polyandrous", "Monandrous/facultatively polyandrous", "Monandrous/facultatively polyandrous to obligately polyandrous", "Obligately polyandrous to monandrous/facultatively polyandrous", "Obligately polyandrous", "Obligately polyandrous to monandrous/facultatively polyandrous", "Obligately polyandrous"), Median = MFcat_fac_ob_PGcat_fac_ob_result$dep_median, Mean = MFcat_fac_ob_PGcat_fac_ob_result$dep_mean, `95% HPD interval` = apply(round(MFcat_fac_ob_PGcat_fac_ob_result$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_fac_ob_PGcat_fac_ob_result$dep_percent_zero, ESS = MFcat_fac_ob_PGcat_fac_ob_result$dep_ESS)
MFcat_fac_ob_PGcat_fac_ob_result_csv <- MFcat_fac_ob_PGcat_fac_ob_result_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
# write.csv(MFcat_fac_ob_PGcat_fac_ob_result_csv, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/BayesTraits/Diagnostics_parameter_estimation/Results_csv/MFcat_fac_ob_PGcat_fac_ob_result.csv")

#Gamma
MFcat_fac_ob_PGcat_fac_ob_result_gamma <- bayestraits_statistics("MF_PG/MF_fac_ob_PG_fac_ob/Gamma/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_ob_reduced.txt.Log.txt",
                                                           "MF_PG/MF_fac_ob_PG_fac_ob/Gamma/Independent/2nd_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_ob_reduced.txt.Log.txt",
                                                           "MF_PG/MF_fac_ob_PG_fac_ob/Gamma/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_ob_reduced.txt.Log.txt",
                                                           "MF_PG/MF_fac_ob_PG_fac_ob/Gamma/Dependent/2nd_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_ob_reduced.txt.Log.txt")
MFcat_fac_ob_PGcat_fac_ob_result_gamma

MFcat_fac_ob_PGcat_fac_ob_BFresult_gamma <- bayes_factors("MF_PG/MF_fac_ob_PG_fac_ob/Gamma/Independent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_ob.txt.Stones.txt",
                                                    "MF_PG/MF_fac_ob_PG_fac_ob/Gamma/Dependent/1st_run/ant_data_BayesTraits_MFcat_fac_ob_PGcat_fac_ob.txt.Stones.txt")
MFcat_fac_ob_PGcat_fac_ob_BFresult_gamma

#Create csv
MFcat_fac_ob_PGcat_fac_ob_result_gamma_df <- data.frame(Rate = names(MFcat_fac_ob_PGcat_fac_ob_result_gamma$dep_median), `Number of queens` = c("Empty", "Monogynous/facultatively polygynous to obligately polygynous", "Monogynous/facultatively polygynous", "Obligately polygynous to monogynous/facultatively polygynous", "Obligately polygynous", "Monogynous/facultatively polygynous", "Monogynous/facultatively polygynous to obligately polygynous", "Obligately polygynous", "Obligately polygynous to monogynous/facultatively polygynous"), `Queen mating frequency` = c("Empty", "Monandrous/facultatively polyandrous", "Monandrous/facultatively polyandrous to obligately polyandrous", "Monandrous/facultatively polyandrous", "Monandrous/facultatively polyandrous to obligately polyandrous", "Obligately polyandrous to monandrous/facultatively polyandrous", "Obligately polyandrous", "Obligately polyandrous to monandrous/facultatively polyandrous", "Obligately polyandrous"), Median = MFcat_fac_ob_PGcat_fac_ob_result_gamma$dep_median, Mean = MFcat_fac_ob_PGcat_fac_ob_result_gamma$dep_mean, `95% HPD interval` = apply(round(MFcat_fac_ob_PGcat_fac_ob_result_gamma$dep_hpd_interval, digits = 2), 1, function(row) paste0("[", paste(row, collapse = ", "), "]")), `% Zero` = MFcat_fac_ob_PGcat_fac_ob_result_gamma$dep_percent_zero, ESS = MFcat_fac_ob_PGcat_fac_ob_result_gamma$dep_ESS)
MFcat_fac_ob_PGcat_fac_ob_result_gamma_csv <- MFcat_fac_ob_PGcat_fac_ob_result_gamma_df[-1, ] %>% mutate_at(c(4,5,7,8), ~ round(., digits = 2))
write.csv(MFcat_fac_ob_PGcat_fac_ob_result_gamma_csv, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/BayesTraits/Diagnostics_parameter_estimation/Results_csv/MFcat_fac_ob_PGcat_fac_ob_result_gamma.csv")


