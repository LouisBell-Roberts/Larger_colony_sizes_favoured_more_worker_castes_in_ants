########Producing figures for the MCMCglmm ASR using reconstructions from corHMM######
###Plotting CS~Caste and MF~Caste results
#Louis Bell-Roberts
#19/02/2024

# packages
library(ape)
library(coda)
library(MCMCglmm)
library(boot)
library(phytools)
library(tidyverse)
library(gridExtra)
library(grid)
library(scales)
library(gtools)

#########
#MF_Caste
#########

#Set working directory
setwd("/Volumes/ADATA SE800/DOL_worker_castes/ASR/MCMCglmm/MF_Caste/1st_run/") 

# get list of all RDS files in the working directory
model_files <- mixedsort(list.files(pattern = "\\.rds$"))

# read in all models using lapply()
mcmc_list <- lapply(model_files, readRDS)

##Look at node estimates for each model: 10^ transformation should be applied to each estimate
# View(as.data.frame(apply((mcmc_list[[1]]$Sol), 2, median)))

#Combine all 400 MCMCglmm Sol model outputs
# extract the Sol objects from each model output in the list
chain_list <- lapply(mcmc_list, function(x) x$Sol)

#Check the number of columns in each model output
num_cols <- lapply(chain_list, ncol)
num_cols==210

#Select only the columns of interest before combining posterior distributions
##Need to do this because the "polymorphic.both" transition is not in every model output. Therefore, the number of columns is not the same in all model outputs
###Columns of interest: "monomorphic.both" and "monomorphic.onlymonomorphic"
chain_list_sub <- lapply(chain_list, function(x) x[, 1:2])

# combine the chains using rbind
MCMC_combined <- do.call(rbind, chain_list_sub)
MCMC_combined <- as.mcmc(MCMC_combined)

#Obtain point estimates - 10^ transformation used as log10 transformation was used prior to the analysis
10^(posterior.mode(MCMC_combined[,1:2])) #CAT2monomorphic.both = 2.462114; CAT2monomorphic.only monomorphic = 1.994259
10^(HPDinterval(MCMC_combined[,1:2])) #CAT2monomorphic.both = 1.1040772, 5.551794; CAT2monomorphic.only monomorphic = 0.9450067, 3.904074

#Plot histograms of the distributions
##Apply 10^ transformation
MCMC_data_frame <- as.data.frame(10^(MCMC_combined))

#Plotting both histograms in the same panel
df_plot_combined <- gather(MCMC_data_frame, key = "column", value = "value") %>%
  mutate(`Transition type` = ifelse(column == "CAT2monomorphic.both", "Single to multiple", "Single to single"))

#Violin plots instead of density plots
df_plot_combined$`Transition type` <- as.factor(df_plot_combined$`Transition type`)

(MF_caste_violin <- ggplot(df_plot_combined, aes(x = fct_rev(`Transition type`), y = value)) +
  geom_violin(fill = "#2466A9", colour = "black") +
  labs(title = "", x = "", y = "Queen mating frequency") +
  theme_classic(base_size = 14) +
  annotate("segment", x = 0.9, xend = 1.1, y = 1.994259, yend = 1.994259,
           color = "black", size = 0.5) +
  annotate("text", x = 0.99, y = 1.994259, label = "1.99",
           color = "black", vjust = -0.8, size = 2.8) +
  annotate("segment", x = 1.9, xend = 2.1, y = 2.462114, yend = 2.462114,
           color = "black", size = 0.5) +
  annotate("text", x = 2, y = 2.462114, label = "2.46",
           color = "black", vjust = -0.8, size = 2.8) +
  scale_y_log10(labels = comma_format()) +
  theme(
    axis.text.x = element_text(size = 10, colour = "black", family = "Helvetica"),
    axis.text.y = element_text(size = 10, colour = "black", family = "Helvetica"),
    axis.title.x = element_text(size = 10, colour = "black", family = "Helvetica"),
    axis.title.y = element_text(size = 10, colour = "black", family = "Helvetica")
  )
)

# ggsave(filename = "/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/ASR_causal/Plots/MF_Caste.pdf", width = 4, height = 4)


########
#CS_Caste
########

#Set working directory
setwd("/Volumes/ADATA SE800/DOL_worker_castes/ASR/MCMCglmm/CS_Caste/1st_run/") 

##Multiple models
# get list of all RDS files in the working directory
model_files <- mixedsort(list.files(pattern = "\\.rds$"))

# read in all models using lapply()
mcmc_list <- lapply(model_files, readRDS)

##Look at node estimates for each model: 10^ transformation should be applied to each estimate
# View(as.data.frame(apply((mcmc_list[[1]]$Sol), 2, median)))

#Combine all 400 MCMCglmm Sol model outputs
# extract the Sol objects from each model output in the list
chain_list <- lapply(mcmc_list, function(x) x$Sol)

#Check the number of columns in each model output
num_cols <- lapply(chain_list, ncol)
num_cols==874

#Select only the columns of interest before combining posterior distributions
##Need to do this because the "polymorphic.both" transition is not in every model output. Therefore, the number of columns is not the same in all model outputs
###Columns of interest: "monomorphic.both" and "monomorphic.onlymonomorphic"
chain_list_sub <- lapply(chain_list, function(x) x[, 1:4])

# combine the chains using rbind
MCMC_combined <- do.call(rbind, chain_list_sub)
MCMC_combined <- as.mcmc(MCMC_combined)

#Obtain point estimates - 10^ transformation used as log10 transformation was used prior to the analysis
10^(posterior.mode(MCMC_combined[,1:4])) #CAT2monomorphic.both = 1025.0506; CAT2monomorphic.only monomorphic = 281.6716
10^(HPDinterval(MCMC_combined[,1:4])) #CAT2monomorphic.only monomorphic = 40.03965, 1590.898; CAT2monomorphic.both = 137.16739, 9064.821

# Compare monomorphic.both (GAIN OF POLYMORPHISM) with monomorphic.only monomorphic (NO CHANGE)
# to test if high colony size makes the origin of worker polymorphism more likely
table(MCMC_combined[,1] > MCMC_combined[,2]) / length(MCMC_combined[,2])	# TRUE = 0.9966. Therefore, colony size is significantly higher in ancestors where transitions to multiple castes occurred

#Plot histograms of the distributions
##Apply 10^ transformation
MCMC_data_frame <- as.data.frame(10^(MCMC_combined))

#Plot
##Plotting both histograms in the same panel
MCMC_data_frame_selected <- MCMC_data_frame %>% dplyr::select(CAT2monomorphic.both, `CAT2monomorphic.only monomorphic`)
df_plot_combined <- gather(MCMC_data_frame_selected, key = "column", value = "value") %>%
  mutate(`Transition type` = ifelse(column == "CAT2monomorphic.both", "Single to multiple", "Single to single"))

#Violin plot instead of density plots
df_plot_combined$`Transition type` <- as.factor(df_plot_combined$`Transition type`)

(CS_caste_violin <- ggplot(df_plot_combined, aes(x = fct_rev(`Transition type`), y = value, fill = `Transition type`)) +
  geom_violin(fill = "#2466A9", color = "black") +
  labs(title = "", x = "", y = "Colony size") +
  theme_classic(base_size = 13) +
  annotate("segment", x = 0.9, xend = 1.1, y = 281.6716, yend = 281.6716,
           color = "black", size = 0.5) +
  annotate("text", x = 1, y = 281.6716, label = "281",
           color = "black", vjust = -0.8, size = 2.8) +
  annotate("segment", x = 1.9, xend = 2.1, y = 1025.0506, yend = 1025.0506,
           color = "black", size = 0.5) +
  annotate("text", x = 1.99, y = 1025.0506, label = "1025",
           color = "black", vjust = -0.8, size = 2.8) +
  scale_y_log10(labels = comma_format()) +
  scale_fill_manual(values = c("white", "white")) +
  guides(fill = "none") +
  theme(
    axis.text.x = element_text(size = 10, colour = "black", family = "Helvetica"),
    axis.text.y = element_text(size = 10, colour = "black", family = "Helvetica"),
    axis.title.x = element_text(size = 10, colour = "black", family = "Helvetica"),
    axis.title.y = element_text(size = 10, colour = "black", family = "Helvetica")
  )
)


##################
#Multi-panel plot
##################

# Create 2-panelled plot
# Create a PDF file
# pdf("/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/ASR_causal/Plots/ASR_2_panel.pdf", width =10, height = 5)
jpeg("/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/ASR_causal/Plots/ASR_2_panel.jpg", width = 180, height = 90, units = "mm", res = 640, quality = 100)

grid.arrange(CS_caste_violin, MF_caste_violin, ncol = 2, nrow = 1)

grid.text("a", x = 0.02, y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("b", x = 0.52, y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
# Close the PDF device
dev.off()

