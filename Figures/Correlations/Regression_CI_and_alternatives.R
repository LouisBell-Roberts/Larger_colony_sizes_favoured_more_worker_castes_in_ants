#Plotting regression models
#Louis Bell-Roberts
#13/02/2024

#Packages
library(tidyverse)
library(ape)
library(phytools)
library(MCMCglmm)
library(scales)
library(grid)
library(gridExtra)

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

ant_data$colony.size <- log10(ant_data$colony.size)
ant_data$queen.mating.frequency <- log10(ant_data$queen.mating.frequency)
ant_data$queen.number.continuous <- log10(ant_data$queen.number.continuous)
ant_data$worker.size.variation <- sqrt(ant_data$worker.size.variation)

#Rename 'species' column as 'animal'
ant_data <- ant_data %>% dplyr::rename(animal = species)

#Read in sample of 400 phylogenetic trees
ant.trees <- read.tree(file ="/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/Economo_2018_400.tre")

###Subsets of the variables
##Analysis predicting number of worker castes
data_MF_caste <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(caste.number))
data_CS_caste <- dplyr::filter(ant_data, complete.cases(colony.size), complete.cases(caste.number))
data_PG_caste <- dplyr::filter(ant_data, complete.cases(queen.number.continuous), complete.cases(caste.number))

##Analysis predicting variation in worker size
data_MF_siz_var <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(worker.size.variation))
data_CS_siz_var <- dplyr::filter(ant_data, complete.cases(colony.size), complete.cases(worker.size.variation))
data_PG_siz_var <- dplyr::filter(ant_data, complete.cases(queen.number.continuous), complete.cases(worker.size.variation))

##Analysis predicting variation in worker size in species with a single worker caste
data_MF_mono_siz_var <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(worker.size.variation), caste.number <2)
data_CS_mono_siz_var <- dplyr::filter(ant_data, complete.cases(colony.size), complete.cases(worker.size.variation), caste.number <2)
data_PG_mono_siz_var <- dplyr::filter(ant_data, complete.cases(queen.number.continuous), complete.cases(worker.size.variation), caste.number <2)

##Pairwise analysis among predictor variables
data_MF_CS <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(colony.size))
data_MF_PG <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency), complete.cases(queen.number.continuous))
data_CS_PG <- dplyr::filter(ant_data, complete.cases(colony.size), complete.cases(queen.number.continuous))

##Prune multiphylo objects for each of the different sets of predictor variables
#Analyses predicting caste number
PT_data_MF_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_caste$animal))
PT_data_CS_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_caste$animal))
PT_data_PG_caste <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PG_caste$animal))

#Analyses predicting variation in worker size
PT_data_MF_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_siz_var$animal))
PT_data_CS_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_siz_var$animal))
PT_data_PG_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PG_siz_var$animal))

#Analyses predicting variation in worker size in species with a single worker caste
PT_data_MF_mono_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_mono_siz_var$animal))
PT_data_CS_mono_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_mono_siz_var$animal))
PT_data_PG_mono_siz_var <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_PG_mono_siz_var$animal))

#Pairwise analyses among the predictor variables
PT_data_MF_CS <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_CS$animal))
PT_data_MF_PG <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_MF_PG$animal))
PT_data_CS_PG <- drop.tip.multiPhylo(ant.trees, setdiff(ant.trees[[1]]$tip.label, data_CS_PG$animal))

##Prune database to match the tree
#Filter through dataframe and select only the rows that match the tips of the tree

#Analyses predicting caste number
data_MF_caste <- filter(data_MF_caste, animal %in% PT_data_MF_caste[[1]]$tip.label)
data_CS_caste <- filter(data_CS_caste, animal %in% PT_data_CS_caste[[1]]$tip.label)
data_PG_caste <- filter(data_PG_caste, animal %in% PT_data_PG_caste[[1]]$tip.label)

#Analyses predicting variation in worker size
data_MF_siz_var <- filter(data_MF_siz_var, animal %in% PT_data_MF_siz_var[[1]]$tip.label)
data_CS_siz_var <- filter(data_CS_siz_var, animal %in% PT_data_CS_siz_var[[1]]$tip.label)
data_PG_siz_var <- filter(data_PG_siz_var, animal %in% PT_data_PG_siz_var[[1]]$tip.label)

#Analyses predicting variation in worker size in species with a single worker caste
data_MF_mono_siz_var <- filter(data_MF_mono_siz_var, animal %in% PT_data_MF_mono_siz_var[[1]]$tip.label)
data_CS_mono_siz_var <- filter(data_CS_mono_siz_var, animal %in% PT_data_CS_mono_siz_var[[1]]$tip.label)
data_PG_mono_siz_var <- filter(data_PG_mono_siz_var, animal %in% PT_data_PG_mono_siz_var[[1]]$tip.label)

#Pairwise analyses among the predictor variables
data_MF_CS <- filter(data_MF_CS, animal %in% PT_data_MF_CS[[1]]$tip.label)
data_MF_PG <- filter(data_MF_PG, animal %in% PT_data_MF_PG[[1]]$tip.label)
data_CS_PG <- filter(data_CS_PG, animal %in% PT_data_CS_PG[[1]]$tip.label)


##########################################################################################
#Regression plots
##########################################################################################

#############
#MF_CS
#############

##Load MCMCglmm regression model
MF_CS_model <- readRDS("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_regressions/Model_outputs_old/MF_CS/1st_chain/MF_CS_1M_100k_1k_1_1stRun.rds")
summary(MF_CS_model)

#Calculate 95% credible interval around regression line using the predict function
MF_CS_ci95 <- as.data.frame(predict(MF_CS_model, data_MF_CS, interval = "confidence", level = 0.95))

# Joining with data_MF_CS_plotting data
data_MF_CS_plotting <- cbind(data_MF_CS,  MF_CS_ci95)

#Plot
# (MF_CS_scatterplot_ci <- ggplot(data_MF_CS_plotting, aes(x = queen.mating.frequency, y = colony.size)) +
#   geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
#   labs(x = ~log[10]~"Queen mating frequency", y = ~log[10]~"Colony size") +
#   theme_classic() +
#   theme(axis.title.x = element_text(colour = "black"),
#         axis.title.y = element_text(colour = "black"),
#         axis.text=element_text(size=17, colour = "black"),
#         axis.title=element_text(size=17),
#         plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
#    geom_smooth(aes(y = lwr), lty = 2, lwd = 0.5, colour = "black", se = F) +
#    geom_smooth(aes(y = upr), lty = 2, lwd = 0.5, colour = "black", se = F) +
#   geom_line(aes(y = fit), size = 0.6) +
#   scale_y_continuous(expand = c(0, 0.1))
# )

(MF_CS_scatterplot_ci <- ggplot(data_MF_CS_plotting, aes(x = 10^queen.mating.frequency, y = 10^colony.size)) +
    geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
    labs(x = "Queen mating frequency", y = "Colony size") +
    theme_classic() +
    theme(axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text=element_text(size=17, colour = "black"),
          axis.title=element_text(size=17),
          plot.margin = unit(c(1, 1, 1.5, 1), "lines")) + 
  scale_y_log10(breaks = c(10, 1000, 100000, 10000000), labels = comma_format(), expand = c(0, 0.45)) +
  scale_x_log10(labels = comma_format()) +
  geom_smooth(aes(y = 10^lwr), lty = 2, lwd = 0.5, colour = "black", se = F) +
  geom_smooth(aes(y = 10^upr), lty = 2, lwd = 0.5, colour = "black", se = F) +
  geom_line(aes(y = 10^fit), size = 0.6)
)


#############
#CS_PG
#############

##Load MCMCglmm regression model
CS_PG_model <- readRDS("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_regressions/Model_outputs/CS_PG/1st_chain/CS_PG_1M_100k_1k_1_1stRun.rds")
summary(CS_PG_model)

#Calculate 95% credible interval around regression line using the predict function
CS_PG_ci95 <- as.data.frame(predict(CS_PG_model, data_CS_PG, interval = "confidence", level = 0.95))

# Joining with data_PG_CS_plotting data
data_CS_PG_plotting <- cbind(data_CS_PG,  CS_PG_ci95)

#Plot
# (CS_PG_scatterplot_ci <- ggplot(data_CS_PG_plotting, aes(x = queen.number.continuous, y = colony.size)) +
#     geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
#     labs(x = ~log[10]~"Queen number", y = ~log[10]~"Colony size") +
#     theme_classic() +
#     theme(axis.title.x = element_text(colour = "black"),
#           axis.title.y = element_text(colour = "black"),
#           axis.text=element_text(size=17, colour = "black"),
#           axis.title=element_text(size=17),
#           plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
#     geom_smooth(aes(y = lwr), lty = 2, lwd = 0.5, colour = "black", se = F) +
#     geom_smooth(aes(y = upr), lty = 2, lwd = 0.5, colour = "black", se = F) +
#     geom_line(aes(y = fit), size = 0.6) +
#     scale_y_continuous(expand = c(0, 0.1))
# )

(CS_PG_scatterplot_ci <- ggplot(data_CS_PG_plotting, aes(x = 10^queen.number.continuous, y = 10^colony.size)) +
    geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
    labs(x = "Number of queens per colony", y = "Colony size") +
    theme_classic() +
    theme(axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text=element_text(size=17, colour = "black"),
          axis.title=element_text(size=17),
          plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
    scale_y_log10(breaks = c(10, 1000, 100000, 10000000), labels = comma_format(), expand = c(0, 0.5)) +
    scale_x_log10(labels = comma_format()) +
    geom_smooth(aes(y = 10^lwr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_smooth(aes(y = 10^upr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_line(aes(y = 10^fit), size = 0.6))




#############
#MF_siz_var
#############

##Load MCMCglmm regression model
MF_siz_var_model <- readRDS("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_regressions/Model_outputs_old/MF_siz_var/1st_chain/MF_siz_var_1M_100k_1k_1_1stRun.rds")
summary(MF_siz_var_model)

#Calculate 95% credible interval around regression line using the predict function
MF_siz_var_ci95 <- as.data.frame(predict(MF_siz_var_model, interval = "confidence", level = 0.95))

# Joining with data_MF_siz_var data
data_MF_siz_var_plotting <- cbind(data_MF_siz_var,  MF_siz_var_ci95)

#Plot
(MF_siz_var_scatterplot_ci <- ggplot(data_MF_siz_var_plotting, aes(x = 10^queen.mating.frequency, y = worker.size.variation)) +
    geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
    labs(x = "Queen mating frequency", y = expression(sqrt("Variation in worker size"))) +
    theme_classic() +
    theme(axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text=element_text(size=17, colour = "black"),
          axis.title=element_text(size=17),
          plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
    geom_smooth(aes(y = lwr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_smooth(aes(y = upr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_line(aes(y = fit), size = 0.6) +
    scale_y_continuous(expand = c(0, 0.1)) +
    scale_x_log10(labels = comma_format())
)


#############
#CS_siz_var
#############

##Load MCMCglmm regression model
CS_siz_var_model <- readRDS("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_regressions/Model_outputs_old/CS_siz_var/1st_chain/CS_siz_var_1M_100k_1k_1_1stRun.rds")
summary(CS_siz_var_model)

#Calculate 95% credible interval around regression line using the predict function
CS_siz_var_ci95 <- as.data.frame(predict(CS_siz_var_model, interval = "confidence", level = 0.95))

# Joining with data_CS_siz_var data
data_CS_siz_var_plotting <- cbind(data_CS_siz_var,  CS_siz_var_ci95)

#Plot
(CS_siz_var_scatterplot_ci <- ggplot(data_CS_siz_var_plotting, aes(x = 10^colony.size, y = worker.size.variation)) +
    geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
    labs(x = "Colony size", y = expression(sqrt("Variation in worker size"))) +
    theme_classic() +
    theme(axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text=element_text(size=17, colour = "black"),
          axis.title=element_text(size=17),
          plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
    geom_smooth(aes(y = lwr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_smooth(aes(y = upr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_line(aes(y = fit), size = 0.6) +
    scale_y_continuous(expand = c(0, 0.1)) +
    scale_x_log10(breaks = c(10, 1000, 100000, 10000000), labels = comma_format(), expand = c(0, 0.5))
)


#############
#PG_siz_var
#############

##Load MCMCglmm regression model
PG_siz_var_model <- readRDS("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_regressions/Model_outputs/PG_siz_var/1st_chain/PG_siz_var_1M_100k_1k_1_1stRun.rds")
summary(PG_siz_var_model)

#Calculate 95% credible interval around regression line using the predict function
PG_siz_var_ci95 <- as.data.frame(predict(PG_siz_var_model, interval = "confidence", level = 0.95))

# Joining with data_PG_siz_var data
data_PG_siz_var_plotting <- cbind(data_PG_siz_var,  PG_siz_var_ci95)

#Plot
(PG_siz_var_scatterplot_ci <- ggplot(data_PG_siz_var_plotting, aes(x = 10^queen.number.continuous, y = worker.size.variation)) +
    geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
    labs(x = "Number of queens per colony", y = expression(sqrt("Variation in worker size"))) +
    theme_classic() +
    theme(axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text=element_text(size=17, colour = "black"),
          axis.title=element_text(size=17),
          plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
    geom_smooth(aes(y = lwr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_smooth(aes(y = upr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_line(aes(y = fit), size = 0.6) +
    scale_y_continuous(expand = c(0, 0.1)) +
    scale_x_log10(labels = comma_format())
)



#############
#MF_mono_siz_var
#############

##Load MCMCglmm regression model
MF_mono_siz_var_model <- readRDS("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_regressions/Model_outputs/MF_mono_siz_var/1st_chain/MF_mono_siz_var_1M_100k_1k_1_1stRun.rds")
summary(MF_mono_siz_var_model)

#Calculate 95% credible interval around regression line using the predict function
MF_mono_siz_var_ci95 <- as.data.frame(predict(MF_mono_siz_var_model, interval = "confidence", level = 0.95))

# Joining with data_MF_siz_var data
data_MF_mono_siz_var_plotting <- cbind(data_MF_mono_siz_var,  MF_mono_siz_var_ci95)

#Plot
(MF_mono_siz_var_scatterplot_ci <- ggplot(data_MF_mono_siz_var_plotting, aes(x = 10^queen.mating.frequency, y = worker.size.variation)) +
    geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
    labs(x = "Queen mating frequency", y = expression(sqrt("Variation in worker size"))) +
    theme_classic() +
    theme(axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text=element_text(size=17, colour = "black"),
          axis.title=element_text(size=17),
          plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
    geom_smooth(aes(y = lwr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_smooth(aes(y = upr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_line(aes(y = fit), size = 0.6) +
    scale_y_continuous(expand = c(0, 0.1)) +
    scale_x_log10(labels = comma_format())
)


#############
#CS_siz_var_mono
#############

##Load MCMCglmm regression model
CS_mono_siz_var_model <- readRDS("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_regressions/Model_outputs_old/CS_mono_siz_var/1st_chain/CS_mono_siz_var_1M_100k_1k_1_1stRun.rds")
summary(CS_mono_siz_var_model)

#Calculate 95% credible interval around regression line using the predict function
CS_mono_siz_var_ci95 <- as.data.frame(predict(CS_mono_siz_var_model, interval = "confidence", level = 0.95))

# Joining with data_CS_mono_siz_var data
data_CS_mono_siz_var_plotting <- cbind(data_CS_mono_siz_var,  CS_mono_siz_var_ci95)

#Plot
(CS_mono_siz_var_scatterplot_ci <- ggplot(data_CS_mono_siz_var_plotting, aes(x = 10^colony.size, y = worker.size.variation)) +
    geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
    labs(x = "Colony size", y = expression(sqrt("Variation in worker size"))) +
    theme_classic() +
    theme(axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text=element_text(size=17, colour = "black"),
          axis.title=element_text(size=17),
          plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
    geom_smooth(aes(y = lwr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_smooth(aes(y = upr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_line(aes(y = fit), size = 0.6) +
    scale_y_continuous(expand = c(0, 0.1)) +
    scale_x_log10(labels = comma_format())
)


#############
#PG_mono_siz_var
#############

##Load MCMCglmm regression model
PG_mono_siz_var_model <- readRDS("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_regressions/Model_outputs/PG_mono_siz_var/1st_chain/PG_mono_siz_var_1M_100k_1k_1_1stRun.rds")
summary(PG_mono_siz_var_model)

#Calculate 95% credible interval around regression line using the predict function
PG_mono_siz_var_ci95 <- as.data.frame(predict(PG_mono_siz_var_model, interval = "confidence", level = 0.95))

# Joining with data_PG_siz_var data
data_PG_mono_siz_var_plotting <- cbind(data_PG_mono_siz_var,  PG_mono_siz_var_ci95)

#Plot
(PG_mono_siz_var_scatterplot_ci <- ggplot(data_PG_mono_siz_var_plotting, aes(x = 10^queen.number.continuous, y = worker.size.variation)) +
    geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
    labs(x = "Number of queens per colony", y = expression(sqrt("Variation in worker size"))) +
    theme_classic() +
    theme(axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text=element_text(size=17, colour = "black"),
          axis.title=element_text(size=17),
          plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
    geom_smooth(aes(y = lwr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_smooth(aes(y = upr), lty = 2, lwd = 0.5, colour = "black", se = F) +
    geom_line(aes(y = fit), size = 0.6) +
    scale_y_continuous(expand = c(0, 0.1)) +
    scale_x_log10(labels = comma_format())
)




##########################################################################################
#Discrete plots
##########################################################################################

#############
#CS_Caste
#############

#Transform the data so that caste number is a categorical variable
data_CS_caste$caste.number <- as.factor(data_CS_caste$caste.number)
# data_CS_caste$colony.size <- 10^data_CS_caste$colony.size

##Scatter and error bar
# 
# #No transformation applied
# (CS_Caste_scatter <- ggplot(data_CS_caste, aes(x = caste.number, y = colony.size)) +
#     geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
#     stat_summary(fun.data = "mean_se", 
#                  geom = c("errorbar"), 
#                  width = 0.1,
#                  position = position_nudge(x = 0.3)) +
#     stat_summary(fun.y = "mean", 
#                  geom = "point",
#                  size = 1.3,
#                  position = position_nudge(x = 0.3)) +
#     theme_classic() +
#     theme(axis.title.x = element_text(colour = "black"),
#           axis.title.y = element_text(colour = "black"),
#           axis.text=element_text(size=22, colour = "black"),
#           axis.title=element_text(size=22),
#           plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
#     labs(y= expression("Log"[10]*" (Colony size)"), x = "Number of worker castes") +
#     scale_y_continuous(breaks = c(1, 1000000, 5000000, 10000000, 15000000), labels = comma_format(), expand = c(0, 1000000))
# )
# 
# #SE calculated before transformation
# (CS_Caste_scatter <- ggplot(data_CS_caste, aes(x = caste.number, y = colony.size)) +
#     geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
#     stat_summary(fun.data = "mean_se", 
#                  geom = c("errorbar"), 
#                  width = 0.1,
#                  position = position_nudge(x = 0.3)) +
#     stat_summary(fun.y = "mean", 
#                  geom = "point",
#                  size = 1.3,
#                  position = position_nudge(x = 0.3)) +
#     theme_classic() +
#     theme(axis.title.x = element_text(colour = "black"),
#           axis.title.y = element_text(colour = "black"),
#           axis.text=element_text(size=22, colour = "black"),
#           axis.title=element_text(size=22),
#           plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
#     labs(y= expression("Log"[10]*" (Colony size)"), x = "Number of worker castes") +
#     coord_trans(y = "log10") +
#     scale_y_continuous(breaks = c(1, 100, 10000, 500000, 1000000), labels = comma_format(), expand = c(0, 1.5))
# )

#SE calculated after transformation - Standard error should be calculated following log10 transformation
(CS_Caste_scatter <- ggplot(data_CS_caste, aes(x = caste.number, y = 10^colony.size)) +
  geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
  stat_summary(fun.data = "mean_se", 
               geom = c("errorbar"), 
               width = 0.19,
               position = position_nudge(x = 0.3)) +
  stat_summary(fun.y = "mean", 
               geom = "point",
               size = 1.3,
               position = position_nudge(x = 0.3)) +
  theme_classic() +
  theme(axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"),
        axis.text=element_text(size=17, colour = "black"),
        axis.title=element_text(size=17),
        plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
  labs(y= "Colony size", x = "Number of worker castes") +
    scale_y_log10(breaks = c(1, 100, 10000, 1000000), labels = comma_format(), expand = c(0, 1.5)) #expand function is necessary to ensure that the 1 tick on the y axis can be fitted in, otherwise it isn't plotted
  )



#############
#MF_Caste
#############

#Transform the data so that caste number is a categorical variable
data_MF_caste$caste.number <- as.factor(data_MF_caste$caste.number)

##Scatter and error bar
(MF_Caste_scatter <- ggplot(data_MF_caste, aes(x = caste.number, y = 10^queen.mating.frequency)) +
    geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
    stat_summary(fun.data = "mean_se", 
                 geom = c("errorbar"), 
                 width = 0.16,
                 position = position_nudge(x = 0.3)) +
    stat_summary(fun.y = "mean", 
                 geom = "point",
                 size = 1.3,
                 position = position_nudge(x = 0.3)) +
    theme_classic() +
    theme(axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text=element_text(size=17, colour = "black"),
          axis.title=element_text(size=17),
          plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
    labs(y= "Queen mating frequency", x = "Number of worker castes") +
    scale_y_log10(breaks = c(1, 3, 10, 30), labels = comma_format(), expand = c(0, 0.11)) #expand function is necessary to ensure that the 1 tick on the y axis can be fitted in, otherwise it isn't plotted
)



#############
#PG_Caste
#############

#Transform the data so that caste number is a categorical variable
data_PG_caste$caste.number <- as.factor(data_PG_caste$caste.number)

##Scatter and error bar
(PG_Caste_scatter <- ggplot(data_PG_caste, aes(x = caste.number, y = 10^queen.number.continuous)) +
    geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
    stat_summary(fun.data = "mean_se", 
                 geom = c("errorbar"), 
                 width = 0.16,
                 position = position_nudge(x = 0.3)) +
    stat_summary(fun.y = "mean", 
                 geom = "point",
                 size = 1.3,
                 position = position_nudge(x = 0.3)) +
    theme_classic() +
    theme(axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text=element_text(size=17, colour = "black"),
          axis.title=element_text(size=17),
          plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
    labs(y= "Number of queens per colony", x = "Number of worker castes") +
    scale_y_log10(labels = comma_format(), expand = c(0, 0.11)) #expand function is necessary to ensure that the 1 tick on the y axis can be fitted in, otherwise it isn't plotted
)



#############
#MF_PG - make a plot where MF is categorical
#############

#Assign queen mating frequency as a categorical variable
data_MF_PG$queen.mating.frequency.categorical <- as.factor(data_MF_PG$queen.mating.frequency.categorical)

##Scatter and error bar
(MF_PG_scatter <- ggplot(data_MF_PG, aes(x = queen.mating.frequency.categorical, y = 10^queen.number.continuous)) +
    geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
    stat_summary(fun.data = "mean_se", 
                 geom = c("errorbar"), 
                 width = 0.12,
                 position = position_nudge(x = 0.3)) +
    stat_summary(fun.y = "mean", 
                 geom = "point",
                 size = 1.3,
                 position = position_nudge(x = 0.3)) +
    theme_classic() +
    theme(axis.title.x = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black"),
          axis.text=element_text(size=17, colour = "black"),
          axis.title=element_text(size=17),
          plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
    labs(y= "Number of queens per colony", x = "Queen mating frequency") +
    scale_y_log10(labels = comma_format(), expand = c(0, 0.11))
)


# #Assign queen mating frequency as a categorical variable
# data_MF_PG$queen.number.categorical <- as.factor(data_MF_PG$queen.number.categorical)
# 
# ##Scatter and error bar
# (MF_PG_scatter <- ggplot(data_MF_PG, aes(x = queen.number.categorical, y = 10^queen.mating.frequency)) +
#     geom_point(size = 2, alpha = 0.4, color = "#2466A9", position = position_jitter(width = 0.1)) +
#     stat_summary(fun.data = "mean_se", 
#                  geom = c("errorbar"), 
#                  width = 0.1,
#                  position = position_nudge(x = 0.3)) +
#     stat_summary(fun.y = "mean", 
#                  geom = "point",
#                  size = 1.3,
#                  position = position_nudge(x = 0.3)) +
#     theme_classic() +
#     theme(axis.title.x = element_text(colour = "black"),
#           axis.title.y = element_text(colour = "black"),
#           axis.text=element_text(size=22, colour = "black"),
#           axis.title=element_text(size=22),
#           plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
#     labs(y= expression("Log"[10]*" (Queen mating frequency)"), x = "Number of queens per colony") +
#     scale_y_log10(labels = comma_format(), expand = c(0, 0.11)) #expand function is necessary to ensure that the 1 tick on the y axis can be fitted in, otherwise it isn't plotted
# )


##########################################################################################
#Overlaid colony size ~ size variation plot
##########################################################################################

#############
#CS_siz_var

##Load MCMCglmm regression model
CS_siz_var_model <- readRDS("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_regressions/Model_outputs_old/CS_siz_var/1st_chain/CS_siz_var_1M_100k_1k_1_1stRun.rds")
summary(CS_siz_var_model)

#Calculate 95% credible interval around regression line using the predict function
CS_siz_var_ci95 <- as.data.frame(predict(CS_siz_var_model, interval = "confidence", level = 0.95))

# Joining with data_CS_siz_var data
data_CS_siz_var_plotting <- cbind(data_CS_siz_var,  CS_siz_var_ci95)

#Assign a value of 0 to all species that have a single caste. All other species are assigned 1, even if they have NA values
data_CS_siz_var_plotting$caste_colour <- ifelse(is.na(data_CS_siz_var_plotting$caste.number) | data_CS_siz_var_plotting$caste.number > 1, 1, 0)

#############
#CS_siz_var_mono
##Load MCMCglmm regression model
CS_mono_siz_var_model <- readRDS("/Volumes/ADATA SE800/DOL_worker_castes/Phylogenetic_regressions/Model_outputs_old/CS_mono_siz_var/1st_chain/CS_mono_siz_var_1M_100k_1k_1_1stRun.rds")
summary(CS_mono_siz_var_model)

#Calculate 95% credible interval around regression line using the predict function
CS_mono_siz_var_ci95 <- as.data.frame(predict(CS_mono_siz_var_model, interval = "confidence", level = 0.95))

# Joining with data_CS_mono_siz_var data
data_CS_mono_siz_var_plotting <- cbind(data_CS_mono_siz_var,  CS_mono_siz_var_ci95)

#Merge dataframes for all species and species with only a single worker caste
CS_siz_var_merged <- merge(data_CS_siz_var_plotting, data_CS_mono_siz_var_plotting, by = "animal", all.x = T)
CS_siz_var_merged_select <- dplyr::select(CS_siz_var_merged, animal, colony.size.x, worker.size.variation.x, lwr.x, upr.x, fit.x, colony.size.y, worker.size.variation.y, lwr.y, upr.y, fit.y, caste_colour)


#############
#Plot
(CS_mono_poly_siz_var_scatterplot_ci_ <- ggplot(CS_siz_var_merged_select, aes(x = 10^colony.size.x, y = worker.size.variation.x, color = factor(caste_colour))) +
   geom_point(size = 2.1, alpha = 1, position = position_jitter(width = 0)) +
   labs(x = "Colony size", y = expression(sqrt("Variation in worker size"))) +
   theme_classic() +
   theme(axis.title.x = element_text(colour = "black"),
         axis.title.y = element_text(colour = "black"),
         axis.text=element_text(size=16, colour = "black"),
         axis.title=element_text(size=16),
         plot.margin = unit(c(1, 1, 1.5, 1), "lines")) +
   # geom_smooth(aes(y = lwr.x), lty = 2, lwd = 0.5, colour = "black", se = F) + #These lines plot separate confidence intervals
   # geom_smooth(aes(y = upr.x), lty = 2, lwd = 0.5, colour = "black", se = F) +
   geom_smooth(aes(y = fit.x), size = 0.6, colour = "black") +
   scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), expand = c(0, 0.03), limits = c(0, 0.9)) +
   scale_x_log10(breaks = c(100, 10000, 1000000), labels = comma_format(), expand = c(0, 0), limits = c(10, 10000000)) +
   # geom_smooth(aes(y = lwr.y), lty = 3, lwd = 0.5, colour = "black", se = F) + #These lines plot separate confidence intervals
   # geom_smooth(aes(y = upr.y), lty = 3, lwd = 0.5, colour = "black", se = F) +
   geom_smooth(aes(y = fit.y), size = 0.6, colour = "#A6CFFA") +
   scale_color_manual(values = c("#A6CFFA", "#2E4B81")) +
   guides(color = "none")
)

# ggsave(filename = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Correlations/Plots/size_var_combined.pdf", plot = CS_mono_poly_siz_var_scatterplot_ci_, width = 4.6, height = 4.7)

##########################################################################################
#Panel plots
##########################################################################################

#Create 6-panel plot
# Create a PDF file
pdf("/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Regressions_and_alternatives/main_text.pdf", width = 11, height = 14)

# Arrange and label plots
grid.arrange(
  CS_Caste_scatter, MF_Caste_scatter, PG_Caste_scatter, MF_PG_scatter, MF_CS_scatterplot_ci, CS_PG_scatterplot_ci,
  ncol = 2, nrow = 3
)

grid.text("a", x = 0.01, y = 0.98, gp = gpar(fontsize = 21, fontface = "bold"))
grid.text("b", x = 0.51, y = 0.98, gp = gpar(fontsize = 21, fontface = "bold"))
grid.text("c", x = 0.01, y = 0.645, gp = gpar(fontsize = 21, fontface = "bold"))
grid.text("d", x = 0.51, y = 0.645, gp = gpar(fontsize = 21, fontface = "bold"))
grid.text("e", x = 0.01, y = 0.315, gp = gpar(fontsize = 21, fontface = "bold"))
grid.text("f", x = 0.51, y = 0.315, gp = gpar(fontsize = 21, fontface = "bold"))
# Close the PDF device
dev.off()


#############
#Create 4-panel plot
# Create a PDF file
pdf("/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Regressions_and_alternatives/supplementary.pdf", width = 10, height = 13)

# Arrange and label plots
grid.arrange(
  MF_siz_var_scatterplot_ci, PG_siz_var_scatterplot_ci, MF_mono_siz_var_scatterplot_ci, PG_mono_siz_var_scatterplot_ci,
  ncol = 2, nrow = 3
)

grid.text("a", x = 0.01, y = 0.98, gp = gpar(fontsize = 21, fontface = "bold"))
grid.text("b", x = 0.51, y = 0.98, gp = gpar(fontsize = 21, fontface = "bold"))
grid.text("c", x = 0.01, y = 0.645, gp = gpar(fontsize = 21, fontface = "bold"))
grid.text("d", x = 0.51, y = 0.645, gp = gpar(fontsize = 21, fontface = "bold"))
# Close the PDF device
dev.off()



