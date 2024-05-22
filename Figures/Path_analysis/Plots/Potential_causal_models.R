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
all_variables <- dplyr::filter(ant_data, complete.cases(caste.number), complete.cases(queen.mating.frequency), complete.cases(colony.size), complete.cases(queen.number.continuous))

#Transform caste.number into a binary variable
all_variables$CasteBin <- as.factor(ifelse(all_variables$caste.number > 1,1,0))

#Prune tree
NCuniform_stem_pruned <- drop.tip(NCuniform_stem, setdiff(NCuniform_stem$tip.label, all_variables$species))
NCuniform_crown_pruned <- drop.tip(NCuniform_crown, setdiff(NCuniform_crown$tip.label, all_variables$species))
FBD_stem_pruned <- drop.tip(FBD_stem, setdiff(FBD_stem$tip.label, all_variables$species))
FBD_crown_pruned <- drop.tip(FBD_crown, setdiff(FBD_crown$tip.label, all_variables$species))

#Prune database
all_variables <- filter(all_variables, species %in% NCuniform_stem_pruned$tip.label)

#Rename the variables used in the analysis
all_variables <- all_variables %>% 
  rename(Mating_frequency=queen.mating.frequency,
         Colony_size=colony.size,
         Caste_number=CasteBin,
         Queen_number=queen.number.continuous)

#Select variables of interest
all_variables <- all_variables %>% dplyr::select(Mating_frequency, Colony_size, Caste_number, Queen_number)


### Define the models for the analysis ###
models <- define_model_set(
  one = c(Caste_number ~ Colony_size, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  two = c(Colony_size ~ Caste_number, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  three = c(Caste_number ~ Colony_size, Caste_number ~ Mating_frequency, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  four = c(Caste_number ~ Mating_frequency, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number)
)

#Plot models
potential_models <- plot_model_set(models, labels = c(Caste_number = "Division_of_labour", Colony_size = "Colony_size", Mating_frequency = "Mating_frequency", Queen_number = "Queen_number"), text_size = 2.9, box_x = 24.5, box_y = 7)
ggsave("/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Path_analysis/Plots/Multi_panels/Potential_models.jpg", device = "jpeg", plot = potential_models, width = 180, height = 130, units = "mm", dpi = "retina")