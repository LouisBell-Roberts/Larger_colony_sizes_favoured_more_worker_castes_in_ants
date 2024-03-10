#Phylogenetic path analysis analysing the number of worker castes
#Louis Bell-Roberts
#15/11/2023

library(tidyverse)
library(ape)
library(phylolm)
library(phytools)
library(ggplot2)
library(phylopath)
library(car)
library(grid)
library(gridExtra)

#Read in data file
ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Following_review/ant_data.csv")

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
NCuniform_stem <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/15k_NCuniform_stem_mcc.tre")
NCuniform_crown <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/15K_NCuniform_crown_mcc.tre")
FBD_stem <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/15K_FBD_stem_mcc.tre")
FBD_crown <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/15K_FBD_crown_mcc.tre")

#Filter data
all_variables <- dplyr::filter(ant_data, complete.cases(caste.number), complete.cases(queen.mating.frequency.categorical), complete.cases(colony.size), complete.cases(queen.number.categorical))

#Transform caste.number, queen.mating.frequency.categorical and queen.number.categorical into binary variables
all_variables$CasteBin <- as.factor(ifelse(all_variables$caste.number > 1,1,0))
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
         Caste_number=CasteBin,
         Queen_number=queen.number.categorical)

#Select variables of interest
all_variables <- all_variables %>% dplyr::select(Mating_frequency, Colony_size, Caste_number, Queen_number)


### Define the models for the analysis ###
models <- define_model_set(
  one = c(Caste_number ~ Colony_size, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  two = c(Colony_size ~ Caste_number, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  three = c(Caste_number ~ Colony_size, Caste_number ~ Mating_frequency, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  four = c(Caste_number ~ Mating_frequency, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number)
)


### Main analysis ###
##NCuniform_stem
NCuniform_stem_result <- phylo_path(models, data = all_variables, tree = NCuniform_stem_pruned, method = "logistic_MPLE", model = "lambda", btol = 20)
NCuniform_stem_result$warnings
summary(NCuniform_stem_result)
plot(summary(NCuniform_stem_result))
NCuniform_stem_result_average_model_full <- average(NCuniform_stem_result, avg_method = "full")
plot(NCuniform_stem_result_average_model_full, algorithm = 'mds', curvature = 0.1, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Caste_number = "CN", Queen_number = "QN"))
NC_stem_plot <- plot(NCuniform_stem_result_average_model_full, algorithm = 'mds', curvature = 0.1, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Caste_number = "CN", Queen_number = "QN"))

##NCuniform_crown
NCuniform_crown_result <- phylo_path(models, data = all_variables, tree = NCuniform_crown_pruned, method = "logistic_MPLE", model = "lambda", btol = 20)
NCuniform_crown_result$warnings
summary(NCuniform_crown_result)
plot(summary(NCuniform_crown_result))
NCuniform_crown_result_average_model_full <- average(NCuniform_crown_result, avg_method = "full")
plot(NCuniform_crown_result_average_model_full, algorithm = 'mds', curvature = 0.1, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Caste_number = "CN", Queen_number = "QN"))
NC_crown_plot <- plot(NCuniform_crown_result_average_model_full, algorithm = 'mds', curvature = 0.1, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Caste_number = "CN", Queen_number = "QN"))

##FBD_stem
FBD_stem_result <- phylo_path(models, data = all_variables, tree = FBD_stem_pruned, method = "logistic_MPLE", model = "lambda", btol = 48) #This model produces a warning regardless of the value that btol is set to. Results should be interpreted with caution.
FBD_stem_result$warnings
summary(FBD_stem_result)
plot(summary(FBD_stem_result))
FBD_stem_result_average_model_full <- average(FBD_stem_result, avg_method = "full")
plot(FBD_stem_result_average_model_full, algorithm = 'mds', curvature = 0.1, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Caste_number = "CN", Queen_number = "QN"))
FBD_stem_plot <- plot(FBD_stem_result_average_model_full, algorithm = 'mds', curvature = 0.1, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Caste_number = "CN", Queen_number = "QN"))

##FBD_crown
FBD_crown_result <- phylo_path(models, data = all_variables, tree = FBD_crown_pruned, method = "logistic_MPLE", model = "lambda", btol = 20)
FBD_crown_result$warnings
summary(FBD_crown_result)
plot(summary(FBD_crown_result))
FBD_crown_result_average_model_full <- average(FBD_crown_result, avg_method = "full")
plot(FBD_crown_result_average_model_full, algorithm = 'mds', curvature = 0.1, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Caste_number = "CN", Queen_number = "QN"))
FBD_crown_plot <- plot(FBD_crown_result_average_model_full, algorithm = 'mds', curvature = 0.1, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Caste_number = "CN", Queen_number = "QN"))


######################################################################


#Create 4-panelled plot
# Create a PDF file
pdf("/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Path_analysis/Plots/Multi_panels/Path_analysis_4panel_discrete_caste_cat.pdf", width = 13, height = 12)
# Arrange and label plots
grid.arrange(
  NC_stem_plot, NC_crown_plot, FBD_crown_plot,
  ncol = 2, nrow = 2
)

grid.text("a", x = 0.01, y = 0.97, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("b", x = 0.51, y = 0.97, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("c", x = 0.01, y = 0.475, gp = gpar(fontsize = 18, fontface = "bold"))
# Close the PDF device
dev.off()


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
write.csv(combined_summaries, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Path_analysis/Discrete_castes/Result_CSV/Categorical/Path_analysis_sumary_discrete_castes_cat.csv", row.names = FALSE)


######################################################################
#Function for extracting Path coefficient, SE, and 95% confidence interval from the Bootstrapped model results
generate_stats <- function(result, value, MCC) {
  CS_CN <- c(round(result$coef["Colony_size", "Caste_number"], 2), 
             round(result$se["Colony_size", "Caste_number"], 2), 
             round(result$lower["Colony_size", "Caste_number"], 2), 
             round(result$upper["Colony_size", "Caste_number"], 2))
  
  MF_CN <- c(round(result$coef["Mating_frequency", "Caste_number"], 2), 
             round(result$se["Mating_frequency", "Caste_number"], 2), 
             round(result$lower["Mating_frequency", "Caste_number"], 2), 
             round(result$upper["Mating_frequency", "Caste_number"], 2))
  
  CN_CS <- c(round(result$coef["Caste_number", "Colony_size"], 2), 
             round(result$se["Caste_number", "Colony_size"], 2), 
             round(result$lower["Caste_number", "Colony_size"], 2), 
             round(result$upper["Caste_number", "Colony_size"], 2))
  
  CS_MF <- c(round(result$coef["Colony_size", "Mating_frequency"], 2), 
             round(result$se["Colony_size", "Mating_frequency"], 2), 
             round(result$lower["Colony_size", "Mating_frequency"], 2), 
             round(result$upper["Colony_size", "Mating_frequency"], 2))
  
  PG_MF <- c(round(result$coef["Queen_number", "Mating_frequency"], 2), 
             round(result$se["Queen_number", "Mating_frequency"], 2), 
             round(result$lower["Queen_number", "Mating_frequency"], 2), 
             round(result$upper["Queen_number", "Mating_frequency"], 2))
  
  stats_noname <- as.data.frame(t(data.frame(CS_affects_CN = CS_CN, 
                                             MF_affects_CN = MF_CN, 
                                             CN_affects_CS = CN_CS,
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

FBD_crown_coef_stats_one <- generate_stats(result = FBD_crown_result_one, value = "One", MCC = "FBD_crown")
FBD_crown_coef_stats_two <- generate_stats(result = FBD_crown_result_two, value = "Two", MCC = "FBD_crown")
FBD_crown_coef_stats_three <- generate_stats(result = FBD_crown_result_three, value = "Three", MCC = "FBD_crown")
FBD_crown_coef_stats_four <- generate_stats(result = FBD_crown_result_four, value = "Four", MCC = "FBD_crown")
FBD_crown_coef_stats_all_mod <- rbind(FBD_crown_coef_stats_one, FBD_crown_coef_stats_two, FBD_crown_coef_stats_three, FBD_crown_coef_stats_four)

Discrete_caste_coef_stats <- rbind(NCuniform_stem_coef_stats_all_mod, NCuniform_crown_coef_stats_all_mod, FBD_stem_coef_stats_all_mod, FBD_crown_coef_stats_all_mod)
write.csv(Discrete_caste_coef_stats, file = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Path_analysis/Discrete_castes/Result_CSV/Categorical/Path_analysis_discrete_castes_coef_stats_cat.csv", row.names = F)
