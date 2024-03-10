#Frequency distributions plots for summary statistics
#Louis Bell-Roberts
#16/02/2024


library(ggplot2)
library(scales)
library(tidyverse)
library(ape)
library(phylolm)
library(phytools)
library(gridExtra)
library(grid)


#Read in data file
ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Following_review/ant_data.csv") #794 species

###########
#Statistics
###########

#Calculate how many different genera in the sample of 794 species and how many different subfamilies

# Split the strings by "_"
genera <- strsplit(ant_data$species, "_")

# Keep only the first element of each split
genera <- sapply(genera, function(x) x[1])

# Get unique elements
genera_unique <- unique(genera)

# Count unique elements
length(genera_unique) #160

#Filter data - different options available
MatingF <- dplyr::filter(ant_data, complete.cases(queen.mating.frequency)) #126
ColonyS <- dplyr::filter(ant_data, complete.cases(colony.size)) #541
CasteN <- dplyr::filter(ant_data, complete.cases(caste.number)) #699
QueenN <- dplyr::filter(ant_data, complete.cases(queen.number.continuous)) #252
Worker_poly <- dplyr::filter(ant_data, complete.cases(worker.size.variation)) #152

#Summary stats - values are for the 794 species
range(MatingF$queen.mating.frequency) #1-25.90
median(MatingF$queen.mating.frequency) #1.1235
range(ColonyS$colony.size) #7-14,750,000
range(CasteN$caste.number) #1-4
range(QueenN$queen.number.continuous) #1-70.48
median(QueenN$queen.number.continuous) #1

quantile(ColonyS$colony.size) #median=300
length(ColonyS$colony.size) #541

length(CasteN$caste.number) #699
length(which(CasteN$caste.number==1))/length(CasteN$caste.number) #0.7625179

length(MatingF$queen.mating.frequency) #126
length(which(MatingF$queen.mating.frequency==1))/length(MatingF$queen.mating.frequency) #0.3095238
length(which(MatingF$queen.mating.frequency > 1 & MatingF$queen.mating.frequency < 2)) / length(MatingF$queen.mating.frequency) #0.3571429
length(which(MatingF$queen.mating.frequency >= 2)) / length(MatingF$queen.mating.frequency) #0.3333333


length(QueenN$queen.number.continuous) #252
length(which(QueenN$queen.number.continuous==1))/length(QueenN$queen.number.continuous) #0.6269841

length(which(QueenN$queen.number.continuous > 1 & QueenN$queen.number.continuous <2))/length(QueenN$queen.number.continuous) #0.1785714
length(which(QueenN$queen.number.continuous >= 2))/length(QueenN$queen.number.continuous) #0.1944444

#Frequency distributions



#MF
(MF_plot <- ggplot(MatingF, aes(queen.mating.frequency)) +
  geom_histogram(binwidth = 1, fill = "white", color = "#CC79A7", size = 1) +
  theme_classic() +
  theme(
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15),
    plot.margin = unit(c(1, 1, 1.5, 1), "lines")
  ) +
  labs(y = "Number of species", x = "Queen mating frequency") + 
  scale_x_continuous(breaks = c(1, 10, 20))
)

#CS
(CS_plot <- ggplot(ColonyS, aes(colony.size)) +
  geom_histogram(binwidth = 0.3, fill = "white", color = "#009E73", size = 1) +
  theme_classic() +
  scale_x_log10(breaks = c(1, 100, 10000, 1000000), labels = comma_format(), expand = c(0, 1.3)) +
  theme(
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15),
    plot.margin = unit(c(1, 1, 1.5, 1), "lines")
  ) +
  labs(y = "Number of species", x = "Colony size")
)

#Worker polymorphism
(WP_plot <- ggplot(Worker_poly, aes(worker.size.variation)) +
  geom_histogram(binwidth = 0.04, fill = "white", color = "#505050", size = 1) +
  theme_classic() +
  theme(
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15),
    plot.margin = unit(c(1, 1, 1.5, 1), "lines")
  ) +
  labs(y = "Number of species", x = "Variation in worker size")
)



#Caste
(Caste_plot <- ggplot(CasteN, aes(caste.number)) +
  geom_histogram(binwidth = 1, color = "#2466A9", fill = "white", size = 0.9) +
  theme_classic() +
  theme(
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15),
    plot.margin = unit(c(1, 1, 1.5, 1), "lines")
  ) +
  scale_y_continuous(breaks = c(0, 200, 400, 600)) +
  labs(y = "Number of species", x = "Number of worker castes") +
  guides(fill = "none")
)

#PG
(PG_plot <- ggplot(QueenN, aes(queen.number.continuous)) +
    geom_histogram(binwidth = 0.2, fill = "white", color = "#FFC300", size = 1) +
    theme_classic() +
    scale_x_log10(breaks = c(1, 3, 10, 30), labels = comma_format(), expand = c(0, 0.1)) +
    # scale_x_log10(labels = comma_format(), expand = c(0, 0.1)) +
    theme(
      axis.title.x = element_text(colour = "black"),
      axis.title.y = element_text(colour = "black"),
      axis.text = element_text(size = 15, colour = "black"),
      axis.title = element_text(size = 15),
      plot.margin = unit(c(1, 1, 1.5, 1), "lines")
    ) +
    labs(y = "Number of species", x = "Number of queens per colony")
)


#############
#Create 6-panel plot
# Create a PDF file
pdf("/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Frequency_summary_stats/Plots/summary_stats.pdf", width = 11, height = 10)

# Arrange and label plots
grid.arrange(Caste_plot, CS_plot, MF_plot, PG_plot, WP_plot, ncol = 2, nrow = 3)

grid.text("a", x = 0.01, y = 0.97, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("b", x = 0.51, y = 0.97, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("c", x = 0.01, y = 0.64, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("d", x = 0.51, y = 0.64, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("e", x = 0.01, y = 0.31, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text("f", x = 0.51, y = 0.31, gp = gpar(fontsize = 18, fontface = "bold"))

# Close the PDF device
dev.off()






