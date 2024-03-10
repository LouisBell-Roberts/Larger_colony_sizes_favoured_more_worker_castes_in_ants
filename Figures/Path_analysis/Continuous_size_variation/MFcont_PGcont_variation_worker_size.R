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
ant_data$queen.number.binary <- as.factor(ant_data$queen.number.binary)

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
all_variables_PGbinary <- dplyr::filter(ant_data, complete.cases(worker.size.variation), complete.cases(queen.mating.frequency), complete.cases(colony.size), complete.cases(queen.number.continuous))

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
         Size_variation=worker.size.variation,
         Queen_number=queen.number.continuous)

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
ggsave(filename = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Path_analysis/Plots/NC_uniform_stem_cont_siz_var_MFcont_PGcont.pdf", width = 6.5, height = 6)


##NCuniform_crown
NCuniform_crown_result <- phylo_path(models, data = all_variables_PGbinary, tree = NCuniform_crown_pruned, method = "logistic_MPLE", model = "lambda")
summary(NCuniform_crown_result)
plot(summary(NCuniform_crown_result))
NCuniform_crown_result_average_model_full <- average(NCuniform_crown_result, avg_method = "full")
plot(NCuniform_crown_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"), type = "colour")
ggsave(filename = "NC_uniform_crown.pdf", width = 6.5, height = 6)

##FBD_stem
FBD_stem_result <- phylo_path(models, data = all_variables_PGbinary, tree = FBD_stem_pruned, method = "logistic_MPLE", model = "lambda")
summary(FBD_stem_result)
plot(summary(FBD_stem_result))
FBD_stem_result_average_model_full <- average(FBD_stem_result, avg_method = "full")
plot(FBD_stem_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"))
ggsave(filename = "FBD_stem.pdf", width = 6.5, height = 6)

##FBD_crown
FBD_crown_result <- phylo_path(models, data = all_variables_PGbinary, tree = FBD_crown_pruned, method = "logistic_MPLE", model = "lambda")
summary(FBD_crown_result)
plot(summary(FBD_crown_result))
FBD_crown_result_average_model_full <- average(FBD_crown_result, avg_method = "full")
plot(FBD_crown_result_average_model_full, algorithm = 'sugiyama', curvature = 0.04, box_x = 10, box_y = 10, text_size = 5, labels = c(Mating_frequency = "MF", Colony_size = "CS", Size_variation = "SV", Queen_number = "QN"))
ggsave(filename = "FBD_crown.pdf", width = 6.5, height = 6)


###
#Create 4-panel plot
# Create a PDF file
pdf("Path_analysis_4panel_size_var_mono.pdf", width = 13, height = 12)
# Arrange and label plots
grid.arrange(
  NC_stem_plot, NC_crown_plot, FBD_stem_plot, FBD_crown_plot,
  ncol = 2, nrow = 2
)

grid.text("a", x = 0.01, y = 0.97, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("b", x = 0.51, y = 0.97, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("c", x = 0.01, y = 0.475, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("d", x = 0.51, y = 0.475, gp = gpar(fontsize = 18, fontface = "bold"))
# Close the PDF device
dev.off()


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


###################################


###Summarise the models
#NCuniform_stem
NCuniform_stem_summary <- NCuniform_stem_result %>% summary() %>% as.data.frame() %>% select(2:9) %>% round(digits = 2) %>% mutate(model = row.names(.)) %>% select(model, everything()) %>% rename("CICc difference" = delta_CICc) %>% filter(`CICc difference` <= 2)
#k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight
rownames(NCuniform_stem_summary) <- NULL
NCuniform_stem_summary <- NCuniform_stem_summary %>%
  mutate(phylogeny = "NC uniform stem") %>%
  select(phylogeny, everything())
# write.csv(NCuniform_stem_summary, file = "NCuniform_stem_summary.csv", row.names = FALSE)

#NCuniform_crown
NCuniform_crown_summary <- NCuniform_crown_result %>% summary() %>% as.data.frame() %>% select(2:9) %>% round(digits = 2) %>% mutate(model = row.names(.)) %>% select(model, everything()) %>% rename("CICc difference" = delta_CICc) %>% filter(`CICc difference` <= 2)
#k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight
rownames(NCuniform_crown_summary) <- NULL
NCuniform_crown_summary <- NCuniform_crown_summary %>%
  mutate(phylogeny = "NC uniform crown") %>%
  select(phylogeny, everything())
# write.csv(NCuniform_crown_summary, file = "NCuniform_crown_summary.csv", row.names = FALSE)

#FBD_stem
FBD_stem_summary <- FBD_stem_result %>% summary() %>% as.data.frame() %>% select(2:9) %>% round(digits = 2) %>% mutate(model = row.names(.)) %>% select(model, everything()) %>% rename("CICc difference" = delta_CICc) %>% filter(`CICc difference` <= 2)
#k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight
rownames(FBD_stem_summary) <- NULL
FBD_stem_summary <- FBD_stem_summary %>%
  mutate(phylogeny = "FBD stem") %>%
  select(phylogeny, everything())
# write.csv(FBD_stem_summary, file = "FBD_stem_summary.csv", row.names = FALSE)

#FBD_crown
FBD_crown_summary <- FBD_crown_result %>% summary() %>% as.data.frame() %>% select(2:9) %>% round(digits = 2) %>% mutate(model = row.names(.)) %>% select(model, everything()) %>% rename("CICc difference" = delta_CICc) %>% filter(`CICc difference` <= 2)
#k=number of conditional independencies tested; q=number of parameters estimated; l=likelihood; w=CICc weight
rownames(FBD_crown_summary) <- NULL
FBD_crown_summary <- FBD_crown_summary %>%
  mutate(phylogeny = "FBD crown") %>%
  select(phylogeny, everything())
# write.csv(FBD_crown_summary, file = "FBD_crown_summary.csv", row.names = FALSE)

##Combine the four data frames
combined_summaries <- rbind(NCuniform_stem_summary, NCuniform_crown_summary, FBD_stem_summary, FBD_crown_summary)
# write.csv(combined_summaries, file = "Path_analysis_summary.csv", row.names = FALSE)


