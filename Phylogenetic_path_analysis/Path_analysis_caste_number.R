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
library(readxl)

#Read in data file
ant_data <- as.data.frame(read_xlsx("ant_data.xlsx", col_names = T))

#Set variables so that they're in the correct structure and apply transformations
ant_data[ant_data == ""] <- NA #Replace blank by NA
class(ant_data$colony.size) #numeric
class(ant_data$queen.mating.frequency) #numeric
ant_data$queen.number <- as.factor(ant_data$queen.number) #Assign as a factor
class(ant_data$queen.number) #factor
class(ant_data$worker.size.variation) #numeric
class(ant_data$caste.number) #numeric

ant_data$colony.size <- log10(ant_data$colony.size)
ant_data$queen.mating.frequency <- log10(ant_data$queen.mating.frequency)
ant_data$worker.size.variation <- sqrt(ant_data$worker.size.variation)

#Set rownames as species names for phylolm package
rownames(ant_data) <- ant_data$species

#Read in phylogenetic trees
NCuniform_stem <- read.tree(file = "15k_NCuniform_stem_mcc.tre")
NCuniform_crown <- read.tree(file = "15K_NCuniform_crown_mcc.tre")
FBD_stem <- read.tree(file = "15K_FBD_stem_mcc.tre")
FBD_crown <- read.tree(file = "15K_FBD_crown_mcc.tre")

#Filter data
all_variables_PGbinary <- dplyr::filter(ant_data, complete.cases(caste.number), complete.cases(queen.mating.frequency), complete.cases(colony.size), complete.cases(queen.number))

#Transform caste.number into a binary variable
all_variables_PGbinary$CasteBin <- as.factor(ifelse(all_variables_PGbinary$caste.number > 1,1,0))

#Prune tree
NCuniform_stem_pruned <- drop.tip(NCuniform_stem, setdiff(NCuniform_stem$tip.label, all_variables_PGbinary$species))
NCuniform_crown_pruned <- drop.tip(NCuniform_crown, setdiff(NCuniform_crown$tip.label, all_variables_PGbinary$species))
FBD_stem_pruned <- drop.tip(FBD_stem, setdiff(FBD_stem$tip.label, all_variables_PGbinary$species))
FBD_crown_pruned <- drop.tip(FBD_crown, setdiff(FBD_crown$tip.label, all_variables_PGbinary$species))

#Prune database
all_variables_PGbinary <- filter(all_variables_PGbinary, species %in% NCuniform_stem_pruned$tip.label)

#Rename the variables used in the analysis
all_variables_PGbinary <- all_variables_PGbinary %>% 
  rename(Mating_frequency=queen.mating.frequency,
         Colony_size=colony.size,
         Caste_number=CasteBin,
         Queen_number=queen.number)

#Select variables of interest
all_variables_PGbinary <- all_variables_PGbinary %>% dplyr::select(Mating_frequency, Colony_size, Caste_number, Queen_number)


### Define the models for the analysis ###
models <- define_model_set(
  one = c(Caste_number ~ Colony_size, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  two = c(Colony_size ~ Caste_number, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  three = c(Caste_number ~ Colony_size, Caste_number ~ Mating_frequency, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  four = c(Caste_number ~ Mating_frequency, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number)
)


### Main analysis ###
##NCuniform_stem
NCuniform_stem_result <- phylo_path(models, data = all_variables_PGbinary, tree = NCuniform_stem_pruned, method = "logistic_MPLE", model = "lambda")
summary(NCuniform_stem_result)
plot(summary(NCuniform_stem_result))
NCuniform_stem_result_average_model_full <- average(NCuniform_stem_result, avg_method = "full")
plot(NCuniform_stem_result_average_model_full, algorithm = 'mds', curvature = 0.1, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Caste_number = "CN", Queen_number = "QN"))

##NCuniform_crown
NCuniform_crown_result <- phylo_path(models, data = all_variables_PGbinary, tree = NCuniform_crown_pruned, method = "logistic_MPLE", model = "lambda")
summary(NCuniform_crown_result)
plot(summary(NCuniform_crown_result))
NCuniform_crown_result_average_model_full <- average(NCuniform_crown_result, avg_method = "full")
plot(NCuniform_crown_result_average_model_full, algorithm = 'mds', curvature = 0.1, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Caste_number = "CN", Queen_number = "QN"))

##FBD_stem
FBD_stem_result <- phylo_path(models, data = all_variables_PGbinary, tree = FBD_stem_pruned, method = "logistic_MPLE", model = "lambda")
summary(FBD_stem_result)
plot(summary(FBD_stem_result))
FBD_stem_result_average_model_full <- average(FBD_stem_result, avg_method = "full")
plot(FBD_stem_result_average_model_full, algorithm = 'mds', curvature = 0.1, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Caste_number = "CN", Queen_number = "QN"))

##FBD_crown
FBD_crown_result <- phylo_path(models, data = all_variables_PGbinary, tree = FBD_crown_pruned, method = "logistic_MPLE", model = "lambda")
summary(FBD_crown_result)
plot(summary(FBD_crown_result))
FBD_crown_result_average_model_full <- average(FBD_crown_result, avg_method = "full")
plot(FBD_crown_result_average_model_full, algorithm = 'mds', curvature = 0.1, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Caste_number = "CN", Queen_number = "QN"))


### Calculate confidence intervals ###
#NCuniform_stem
NCuniform_stem_result_one <- choice(NCuniform_stem_result, "one", boot = 500)

#NCuniform_crown
NCuniform_crown_result_one <- choice(NCuniform_crown_result, "one", boot = 500)

#FBD_stem
FBD_stem_result_one <- choice(FBD_stem_result, "one", boot = 500)

#FBD_crown
FBD_crown_result_one <- choice(FBD_crown_result, "one", boot = 500)
FBD_crown_result_two <- choice(FBD_crown_result, "two", boot = 500)

