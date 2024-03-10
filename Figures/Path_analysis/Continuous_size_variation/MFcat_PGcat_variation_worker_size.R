#Phylogenetic path analysis analysing variation in worker size
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
all_variables_PGbinary <- dplyr::filter(ant_data, complete.cases(worker.size.variation), complete.cases(queen.mating.frequency.categorical), complete.cases(colony.size), complete.cases(queen.number.categorical))

#Transform queen.mating.frequency.categorical and queen.number.categorical into binary variables
all_variables_PGbinary$queen.mating.frequency.categorical <- as.factor(ifelse(all_variables_PGbinary$queen.mating.frequency.categorical > 1,1,0))
all_variables_PGbinary$queen.number.categorical <- as.factor(ifelse(all_variables_PGbinary$queen.number.categorical > 0,1,0))

#Prune tree
NCuniform_stem_pruned <- drop.tip(NCuniform_stem, setdiff(NCuniform_stem$tip.label, all_variables_PGbinary$species))
NCuniform_crown_pruned <- drop.tip(NCuniform_crown, setdiff(NCuniform_crown$tip.label, all_variables_PGbinary$species))
FBD_stem_pruned <- drop.tip(FBD_stem, setdiff(FBD_stem$tip.label, all_variables_PGbinary$species))
FBD_crown_pruned <- drop.tip(FBD_crown, setdiff(FBD_crown$tip.label, all_variables_PGbinary$species))

#Prune database
all_variables_PGbinary <- filter(all_variables_PGbinary, species %in% NCuniform_stem_pruned$tip.label)

#Rename the variables used in the analysis
all_variables_PGbinary <- all_variables_PGbinary %>%
  rename(Mating_frequency=queen.mating.frequency.categorical,
         Colony_size=colony.size,
         Size_variation=worker.size.variation,
         Queen_number=queen.number.categorical)

#Select variables of interest
all_variables_PGbinary <- all_variables_PGbinary %>% dplyr::select(Mating_frequency, Colony_size, Size_variation, Queen_number)


### Define the models for the analysis ###
models <- define_model_set(
  one = c(Size_variation ~ Colony_size, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  two = c(Colony_size ~ Size_variation, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  three = c(Size_variation ~ Colony_size, Size_variation ~ Mating_frequency, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number),
  four = c(Size_variation ~ Mating_frequency, Mating_frequency ~ Colony_size, Mating_frequency ~ Queen_number)
)

### Main analysis ###
##NCuniform_stem
NCuniform_stem_result <- phylo_path(models, data = all_variables_PGbinary, tree = NCuniform_stem_pruned, method = "logistic_MPLE", model = "lambda")
summary(NCuniform_stem_result)
plot(summary(NCuniform_stem_result))
NCuniform_stem_result_average_model_full <- average(NCuniform_stem_result, avg_method = "full")
plot(NCuniform_stem_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"))
ggsave(filename = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Path_analysis/Plots/NC_uniform_stem_cont_siz_var_MFcat_PGcat.pdf", width = 6.5, height = 6)


##NCuniform_crown
NCuniform_crown_result <- phylo_path(models, data = all_variables_PGbinary, tree = NCuniform_crown_pruned, method = "logistic_MPLE", model = "lambda")
summary(NCuniform_crown_result)
plot(summary(NCuniform_crown_result))
NCuniform_crown_result_average_model_full <- average(NCuniform_crown_result, avg_method = "full")
plot(NCuniform_crown_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"), type = "colour")

##FBD_stem
FBD_stem_result <- phylo_path(models, data = all_variables_PGbinary, tree = FBD_stem_pruned, method = "logistic_MPLE", model = "lambda")
summary(FBD_stem_result)
plot(summary(FBD_stem_result))
FBD_stem_result_average_model_full <- average(FBD_stem_result, avg_method = "full")
plot(FBD_stem_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"))
ggsave(filename = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Path_analysis/Plots/FBD_stem_cont_siz_var_MFcat_PGcat.pdf", width = 6.5, height = 6)


##FBD_crown
FBD_crown_result <- phylo_path(models, data = all_variables_PGbinary, tree = FBD_crown_pruned, method = "logistic_MPLE", model = "lambda")
summary(FBD_crown_result)
plot(summary(FBD_crown_result))
FBD_crown_result_average_model_full <- average(FBD_crown_result, avg_method = "full")
plot(FBD_crown_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"))

### Calculate confidence intervals ###
#NCuniform_stem
NCuniform_stem_result_one <- choice(NCuniform_stem_result, "one", boot = 500)
NCuniform_stem_result_two <- choice(NCuniform_stem_result, "two", boot = 500)
NCuniform_stem_result_three <- choice(NCuniform_stem_result, "three", boot = 500)

#NCuniform_crown
NCuniform_crown_result_one <- choice(NCuniform_crown_result, "one", boot = 500)
NCuniform_crown_result_two <- choice(NCuniform_crown_result, "two", boot = 500)
NCuniform_crown_result_three <- choice(NCuniform_crown_result, "three", boot = 500)

#FBD_stem
FBD_stem_result_one <- choice(FBD_stem_result, "one", boot = 500)
FBD_stem_result_two <- choice(FBD_stem_result, "two", boot = 500)
FBD_stem_result_three <- choice(FBD_stem_result, "three", boot = 500)

#FBD_crown
FBD_crown_result_one <- choice(FBD_crown_result, "one", boot = 500)
FBD_crown_result_two <- choice(FBD_crown_result, "two", boot = 500)
FBD_crown_result_three <- choice(FBD_crown_result, "three", boot = 500)

