###Phylogenetic path analysis analysing variation in worker size.
##Analysing queen mating frequency and queen number as binary variables
#Louis Bell-Roberts
#15/11/2023


library(tidyverse)
library(ape)
library(phylolm)
library(phytools)
library(ggplot2)
library(phylopath)
library(car)

#Read in data file
ant_data <- read.csv("ant_data.csv")

#Set variables so that they're in the correct structure and apply transformations
ant_data[ant_data == ""] <- NA #Replace blank by NA
class(ant_data$colony.size) #numeric
class(ant_data$queen.mating.frequency) #numeric
class(ant_data$queen.number.continuous) #numeric
class(ant_data$worker.size.variation) #numeric
class(ant_data$caste.number) #numeric

ant_data$colony.size <- log10(ant_data$colony.size)
ant_data$queen.mating.frequency <- log10(ant_data$queen.mating.frequency)
ant_data$queen.number.continuous <- log10(ant_data$queen.number.continuous)
ant_data$worker.size.variation <- sqrt(ant_data$worker.size.variation)

#Set rownames as species names for phylolm package
rownames(ant_data) <- ant_data$species

#Read in phylogenetic trees
NCuniform_stem <- read.tree(file = "15k_NCuniform_stem_mcc.tre")
NCuniform_crown <- read.tree(file = "15K_NCuniform_crown_mcc.tre")
FBD_stem <- read.tree(file = "15K_FBD_stem_mcc.tre")
FBD_crown <- read.tree(file = "15K_FBD_crown_mcc.tre")

#Filter data
all_variables <- dplyr::filter(ant_data, complete.cases(worker.size.variation), complete.cases(queen.mating.frequency.categorical), complete.cases(colony.size), complete.cases(queen.number.categorical))

#Transform queen.mating.frequency.categorical and queen.number.categorical into binary variables
all_variables$queen.mating.frequency.categorical <- as.factor(ifelse(all_variables$queen.mating.frequency.categorical > 1,1,0))
all_variables$queen.number.categorical <- as.factor(ifelse(all_variables$queen.number.categorical > 1,1,0))

#Prune tree
NCuniform_stem_pruned <- drop.tip(NCuniform_stem, setdiff(NCuniform_stem$tip.label, all_variables$species))
NCuniform_crown_pruned <- drop.tip(NCuniform_crown, setdiff(NCuniform_crown$tip.label, all_variables$species))
FBD_stem_pruned <- drop.tip(FBD_stem, setdiff(FBD_stem$tip.label, all_variables$species))
FBD_crown_pruned <- drop.tip(FBD_crown, setdiff(FBD_crown$tip.label, all_variables$species))

#Prune database
all_variables <- filter(all_variables, species %in% NCuniform_stem_pruned$tip.label)

#Rename the variables used in the analysis
all_variables <- all_variables %>%
  rename(Mating_frequency=queen.mating.frequency.categorical,
         Colony_size=colony.size,
         Size_variation=worker.size.variation,
         Queen_number=queen.number.categorical)

#Select variables of interest
all_variables <- all_variables %>% dplyr::select(Mating_frequency, Colony_size, Size_variation, Queen_number)


### Define the models for the analysis ###
models <- define_model_set(
  one = c(Size_variation ~ Colony_size, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  two = c(Colony_size ~ Size_variation, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  three = c(Size_variation ~ Colony_size, Size_variation ~ Mating_frequency, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  four = c(Size_variation ~ Mating_frequency, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number)
)

### Main analysis ###
##NCuniform_stem
NCuniform_stem_result <- phylo_path(models, data = all_variables, tree = NCuniform_stem_pruned, method = "logistic_MPLE", model = "lambda", btol = 30)
summary(NCuniform_stem_result)
plot(summary(NCuniform_stem_result))
NCuniform_stem_result_average_model_full <- average(NCuniform_stem_result, avg_method = "full")
plot(NCuniform_stem_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"))
NC_stem_plot <- plot(NCuniform_stem_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"))

##NCuniform_crown
NCuniform_crown_result <- phylo_path(models, data = all_variables, tree = NCuniform_crown_pruned, method = "logistic_MPLE", model = "lambda", btol = 30)
summary(NCuniform_crown_result)
plot(summary(NCuniform_crown_result))
NCuniform_crown_result_average_model_full <- average(NCuniform_crown_result, avg_method = "full")
plot(NCuniform_crown_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"))
NC_crown_plot <- plot(NCuniform_crown_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"))

##FBD_stem
FBD_stem_result <- phylo_path(models, data = all_variables, tree = FBD_stem_pruned, method = "logistic_MPLE", model = "lambda", btol = 30)
summary(FBD_stem_result)
plot(summary(FBD_stem_result))
FBD_stem_result_average_model_full <- average(FBD_stem_result, avg_method = "full")
plot(FBD_stem_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"))
FBD_stem_plot <- plot(FBD_stem_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"))

##FBD_crown
FBD_crown_result <- phylo_path(models, data = all_variables, tree = FBD_crown_pruned, method = "logistic_MPLE", model = "lambda", btol = 30)
summary(FBD_crown_result)
plot(summary(FBD_crown_result))
FBD_crown_result_average_model_full <- average(FBD_crown_result, avg_method = "full")
plot(FBD_crown_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"))
FBD_crown_plot <- plot(FBD_crown_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"))


######################################################################


###Summarise the models
##k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight, p=p-value

#NCuniform_stem
NCuniform_stem_summary <- NCuniform_stem_result %>% summary() %>% as.data.frame() %>% select(2:9) %>% round(digits = 2) %>% mutate(model = row.names(.)) %>% select(model, everything()) %>% rename("CICc difference" = delta_CICc)
#k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight
rownames(NCuniform_stem_summary) <- NULL
NCuniform_stem_summary <- NCuniform_stem_summary %>%
  mutate(phylogeny = "NC uniform stem") %>%
  select(phylogeny, everything())
# write.csv(NCuniform_stem_summary, file = "NCuniform_stem_summary.csv", row.names = FALSE)

#NCuniform_crown
NCuniform_crown_summary <- NCuniform_crown_result %>% summary() %>% as.data.frame() %>% select(2:9) %>% round(digits = 2) %>% mutate(model = row.names(.)) %>% select(model, everything()) %>% rename("CICc difference" = delta_CICc)
#k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight
rownames(NCuniform_crown_summary) <- NULL
NCuniform_crown_summary <- NCuniform_crown_summary %>%
  mutate(phylogeny = "NC uniform crown") %>%
  select(phylogeny, everything())
# write.csv(NCuniform_crown_summary, file = "NCuniform_crown_summary.csv", row.names = FALSE)

#FBD_stem
FBD_stem_summary <- FBD_stem_result %>% summary() %>% as.data.frame() %>% select(2:9) %>% round(digits = 2) %>% mutate(model = row.names(.)) %>% select(model, everything()) %>% rename("CICc difference" = delta_CICc)
#k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight
rownames(FBD_stem_summary) <- NULL
FBD_stem_summary <- FBD_stem_summary %>%
  mutate(phylogeny = "FBD stem") %>%
  select(phylogeny, everything())
# write.csv(FBD_stem_summary, file = "FBD_stem_summary.csv", row.names = FALSE)

#FBD_crown
FBD_crown_summary <- FBD_crown_result %>% summary() %>% as.data.frame() %>% select(2:9) %>% round(digits = 2) %>% mutate(model = row.names(.)) %>% select(model, everything()) %>% rename("CICc difference" = delta_CICc)
#k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight
rownames(FBD_crown_summary) <- NULL
FBD_crown_summary <- FBD_crown_summary %>%
  mutate(phylogeny = "FBD crown") %>%
  select(phylogeny, everything())
# write.csv(FBD_crown_summary, file = "FBD_crown_summary.csv", row.names = FALSE)

##Combine the four data frames
combined_summaries <- rbind(NCuniform_stem_summary, NCuniform_crown_summary, FBD_stem_summary, FBD_crown_summary)
# write.csv(combined_summaries, file = "Path_analysis_summary_siz_var_cat.csv", row.names = FALSE)


######################################################################
#Function for extracting Path coefficient, SE, and 95% confidence interval from the Bootstrapped model results
generate_stats <- function(result, value, MCC) {
  CS_SV <- c(round(result$coef["Colony_size", "Size_variation"], 2), 
             round(result$se["Colony_size", "Size_variation"], 2), 
             round(result$lower["Colony_size", "Size_variation"], 2), 
             round(result$upper["Colony_size", "Size_variation"], 2))
  
  MF_SV <- c(round(result$coef["Mating_frequency", "Size_variation"], 2), 
             round(result$se["Mating_frequency", "Size_variation"], 2), 
             round(result$lower["Mating_frequency", "Size_variation"], 2), 
             round(result$upper["Mating_frequency", "Size_variation"], 2))
  
  SV_CS <- c(round(result$coef["Size_variation", "Colony_size"], 2), 
             round(result$se["Size_variation", "Colony_size"], 2), 
             round(result$lower["Size_variation", "Colony_size"], 2), 
             round(result$upper["Size_variation", "Colony_size"], 2))
  
  CS_MF <- c(round(result$coef["Colony_size", "Mating_frequency"], 2), 
             round(result$se["Colony_size", "Mating_frequency"], 2), 
             round(result$lower["Colony_size", "Mating_frequency"], 2), 
             round(result$upper["Colony_size", "Mating_frequency"], 2))
  
  PG_MF <- c(round(result$coef["Queen_number", "Mating_frequency"], 2), 
             round(result$se["Queen_number", "Mating_frequency"], 2), 
             round(result$lower["Queen_number", "Mating_frequency"], 2), 
             round(result$upper["Queen_number", "Mating_frequency"], 2))
  
  stats_noname <- as.data.frame(t(data.frame(CS_affects_SV = CS_SV, 
                                             MF_affects_SV = MF_SV, 
                                             SV_affects_CS = SV_CS,
                                             CS_affects_MF = CS_MF,
                                             PG_affects_MF = PG_MF)))
  
  stats <- rename(stats_noname, 
                  Path_coefficient = V1, 
                  SE = V2, 
                  Lower_95_CI = V3, 
                  Upper_95_CI = V4)
  
  # Add the column "Model_number" to the first position
  stats <- mutate(stats, Model_number := value, Tree := MCC)
  
  # Add row names to a new column
  stats <- rownames_to_column(stats, var = "Causal_relationship")
  
  return(stats)
}

### Calculate confidence intervals ###
#NCuniform_stem
NCuniform_stem_result_one <- choice(NCuniform_stem_result, "one", boot = 500)
NCuniform_stem_result_two <- choice(NCuniform_stem_result, "two", boot = 500)
NCuniform_stem_result_three <- choice(NCuniform_stem_result, "three", boot = 500)
NCuniform_stem_result_four <- choice(NCuniform_stem_result, "four", boot = 500)
NCuniform_stem_result_average_model_full #CI for the average model

NCuniform_stem_coef_stats_one <- generate_stats(result = NCuniform_stem_result_one, value = "One", MCC = "NCuniform_stem")
NCuniform_stem_coef_stats_two <- generate_stats(result = NCuniform_stem_result_two, value = "Two", MCC = "NCuniform_stem")
NCuniform_stem_coef_stats_three <- generate_stats(result = NCuniform_stem_result_three, value = "Three", MCC = "NCuniform_stem")
NCuniform_stem_coef_stats_four <- generate_stats(result = NCuniform_stem_result_four, value = "Four", MCC = "NCuniform_stem")
NCuniform_stem_coef_stats_all_mod <- rbind(NCuniform_stem_coef_stats_one, NCuniform_stem_coef_stats_two, NCuniform_stem_coef_stats_three, NCuniform_stem_coef_stats_four)

#NCuniform_crown
NCuniform_crown_result_one <- choice(NCuniform_crown_result, "one", boot = 500)
NCuniform_crown_result_two <- choice(NCuniform_crown_result, "two", boot = 500)
NCuniform_crown_result_three <- choice(NCuniform_crown_result, "three", boot = 500)
NCuniform_crown_result_four <- choice(NCuniform_crown_result, "four", boot = 500)
NCuniform_crown_result_average_model_full #CI for the average model

NCuniform_crown_coef_stats_one <- generate_stats(result = NCuniform_crown_result_one, value = "One", MCC = "NCuniform_crown")
NCuniform_crown_coef_stats_two <- generate_stats(result = NCuniform_crown_result_two, value = "Two", MCC = "NCuniform_crown")
NCuniform_crown_coef_stats_three <- generate_stats(result = NCuniform_crown_result_three, value = "Three", MCC = "NCuniform_crown")
NCuniform_crown_coef_stats_four <- generate_stats(result = NCuniform_crown_result_four, value = "Four", MCC = "NCuniform_crown")
NCuniform_crown_coef_stats_all_mod <- rbind(NCuniform_crown_coef_stats_one, NCuniform_crown_coef_stats_two, NCuniform_crown_coef_stats_three, NCuniform_crown_coef_stats_four)

#FBD_stem
FBD_stem_result_one <- choice(FBD_stem_result, "one", boot = 500)
FBD_stem_result_two <- choice(FBD_stem_result, "two", boot = 500)
FBD_stem_result_three <- choice(FBD_stem_result, "three", boot = 500)
FBD_stem_result_four <- choice(FBD_stem_result, "four", boot = 500)
FBD_stem_result_average_model_full #CI for the average model

FBD_stem_coef_stats_one <- generate_stats(result = FBD_stem_result_one, value = "One", MCC = "FBD_stem")
FBD_stem_coef_stats_two <- generate_stats(result = FBD_stem_result_two, value = "Two", MCC = "FBD_stem")
FBD_stem_coef_stats_three <- generate_stats(result = FBD_stem_result_three, value = "Three", MCC = "FBD_stem")
FBD_stem_coef_stats_four <- generate_stats(result = FBD_stem_result_four, value = "Four", MCC = "FBD_stem")
FBD_stem_coef_stats_all_mod <- rbind(FBD_stem_coef_stats_one, FBD_stem_coef_stats_two, FBD_stem_coef_stats_three, FBD_stem_coef_stats_four)

#FBD_crown
FBD_crown_result_one <- choice(FBD_crown_result, "one", boot = 500)
FBD_crown_result_two <- choice(FBD_crown_result, "two", boot = 500)
FBD_crown_result_three <- choice(FBD_crown_result, "three", boot = 500)
FBD_crown_result_four <- choice(FBD_crown_result, "four", boot = 500)
FBD_crown_result_average_model_full #CI for the average model

FBD_crown_coef_stats_one <- generate_stats(result = FBD_crown_result_one, value = "One", MCC = "FBD_crown")
FBD_crown_coef_stats_two <- generate_stats(result = FBD_crown_result_two, value = "Two", MCC = "FBD_crown")
FBD_crown_coef_stats_three <- generate_stats(result = FBD_crown_result_three, value = "Three", MCC = "FBD_crown")
FBD_crown_coef_stats_four <- generate_stats(result = FBD_crown_result_four, value = "Four", MCC = "FBD_crown")
FBD_crown_coef_stats_all_mod <- rbind(FBD_crown_coef_stats_one, FBD_crown_coef_stats_two, FBD_crown_coef_stats_three, FBD_crown_coef_stats_four)

Siz_var_coef_stats <- rbind(NCuniform_stem_coef_stats_all_mod, NCuniform_crown_coef_stats_all_mod, FBD_stem_coef_stats_all_mod, FBD_crown_coef_stats_all_mod)
# write.csv(Siz_var_coef_stats, file = "Path_analysis_siz_var_coef_stats_cat.csv", row.names = F)




