########################Plot transitions of the corHMM ASR for caste number and calculate summary statistics for each transition type######################
#Louis Bell-Roberts
#16/02/2024

#Load packages
library(tidyverse)
library(ape)
library(phytools)
library(ggplot2)
library(geiger)
library(corHMM)
library(phangorn)
library(gtools)

setwd("/Volumes/ADATA SE800/DOL_worker_castes/ASR/corHMM/Caste_only_discrete_1to4/2_rat_SYM/")

#Load custom functions

#Function that takes a matrix or data frame as input and returns a vector of the column indices corresponding to the maximum values in each row - replicates the behaviour of the max.col function
my_max_col <- function(x) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("Input must be a matrix or data frame")
  }
  apply(x, 1, which.max)
}

#The extract_frequency function is designed to extract frequencies from a list of data frames (table_list) based on a specified pattern found in the 'Var1' column of each data frame. Here's a summary of what the function does:
extract_frequency <- function(table_list, pattern) {
  lapply(table_list, function(df) {
    if (pattern %in% df$Var1) {
      subset_df <- df[grep(pattern, df$Var1), ]
      # Extract the 'Freq' column from the subsetted data frame
      freq <- subset_df$Freq
    } else {
      # If 'pattern' is not present, return 0
      freq <- 0
    }
    return(freq)
  })
}

###############################################################

#Read in data file
ant_data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Following_review/ant_data.csv")
ant_data <- ant_data %>% rename(animal = species) %>% filter(complete.cases(caste.number))

#Read in sample of 400 phylogenetic trees
ant_trees <- read.tree(file ="/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Master_cloud_data/Publication/Trees/Economo_2018_400.tre")

#Prune trees to match data
ant_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, ant_data$animal))

#Prune data to match trees
sxData <- filter(ant_data, animal %in% ant_trees_pruned[[1]]$tip.label)


###############
###Read in the 400 corHMM model results in ascending order
###############

# create an empty list to hold the results
ASR_models <- list()

# set folder path
folder_path <- "/Volumes/ADATA SE800/DOL_worker_castes/ASR/corHMM/Caste_only_discrete_1to4/2_rat_SYM/"

# get the file names
file_names <- list.files(folder_path)

# sort the file names based on the corresponding numbers
sorted_file_names <- mixedsort(file_names)

# read in the files and save them to a list
ASR_models <- lapply(sorted_file_names, function(x) readRDS(file.path(folder_path, x)))


###############
###For each tree, estimate the most likely number of castes for each node of the phylogeny
###############

#Convert ASR estimates to list of data frames
ASR_states_list <- lapply(ASR_models, function(df) as.data.frame(df$states))

# Create an empty list to store the results
TransDat_neg_rep_list <- vector(mode = "list", length = length(ASR_states_list))

##############
for (i in seq_along(ASR_states_list)) {
  ASR_states <- ASR_states_list[[i]]
  
  #Combine the probability scores for the two different rate classes for each caste number
  ASR_states_summed <- data.frame(caste_1 = ASR_states$`(1,R1)` + ASR_states$`(1,R2)`, caste_2 = ASR_states$`(2,R1)` + ASR_states$`(2,R2)`, caste_3 = ASR_states$`(3,R1)` + ASR_states$`(3,R2)`, caste_4 = ASR_states$`(4,R1)` + ASR_states$`(4,R2)`)
  
  #Identify which caste number is the highest probability for each node
  ASR_states_summed$max_col <- colnames(ASR_states_summed)[my_max_col(ASR_states_summed)] #originally was using the max.col() function however it was making rounding errors
  
  #Assign node number to each row of the data frame, given that node number corresponds to row number in the ASR_states_summed dataframe
  ASR_states_summed$node_number <- paste0("animal.Node", seq_len(nrow(ASR_states_summed)))
  
  #Create new list of data frames with just the node_state and node_number columns
  ASR_states_summed <- ASR_states_summed %>% dplyr::select(max_col, node_number)
  
  #######
  
  #Get data for each species tip and its caste data
  ##Create 'nodecode' column" which pastes animal and then also the species tip
  sxData$nodecode <- paste("animal", sxData$animal, sep=".")
  # Create new data frame with just the caste.number and nodecode columns
  sxData_minimal <- sxData %>% dplyr::select(caste.number, nodecode)
  #Make the row names equal to nodecode column
  rownames(sxData_minimal) <- sxData_minimal[, "nodecode"]
  
  #Join 'ASR_states_summed' on top of 'sxData_minimal' using rbind
  ##Use lapply to bind sxData_minimal to each element in ASR_states_summed
  ASR_states_summed <- ASR_states_summed %>% dplyr::rename(poly2 = max_col, animal = node_number)
  sxData_minimal <- sxData_minimal %>% dplyr::rename(poly2 = caste.number, animal = nodecode)
  
  sxData_minimal$poly2 <- ifelse(sxData_minimal$poly2 == 1, "caste_1",
                                 ifelse(sxData_minimal$poly2 == 2, "caste_2",
                                        ifelse(sxData_minimal$poly2 == 3, "caste_3",
                                               ifelse(sxData_minimal$poly2 == 4, "caste_4", NA))))
  
  Poly <- rbind(ASR_states_summed, sxData_minimal)
  
  
  #Turn tree edges into a data frame - displays the parent node (V1) and the two offspring nodes (V2). 653 is the root node number as there are 652 tips in the tree
  TransDat <- as.data.frame(ASR_models[[i]]$phy$edge)
  
  #Name each of the nodes in the phylogeny with correct name
  TransDat$V1name <- paste("animal.Node", TransDat$V1-length(ASR_models[[i]]$phy$tip.label), sep = "")
  
  #Do the same but for the descendant nodes
  TransDat$V2name <- paste("animal.Node", TransDat$V2-length(ASR_models[[i]]$phy$tip.label), sep = "")
  
  #Number each of the species in the tree by their tip number
  treesp <- data.frame(tip.label=paste("animal.", ASR_models[[i]]$phy$tip.label,sep=""), no=1:length(ASR_models[[i]]$phy$tip.label))
  
  #Match species names with the descendant nodes
  TransDat_neg_rep <- data.frame(TransDat, ancestors=NA, descendents=treesp$tip.label[match(TransDat$V2, treesp$no)])
  
  # make sure these are of class 'character'
  TransDat_neg_rep$ancestors <- as.character(TransDat_neg_rep$ancestors)
  TransDat_neg_rep$descendents <- as.character(TransDat_neg_rep$descendents)
  
  # ancestors currently NAs, so give these V1name 
  TransDat_neg_rep$ancestors <- as.character(ifelse(!is.na(TransDat_neg_rep$ancestors),TransDat_neg_rep$ancestors,TransDat_neg_rep$V1name))
  
  # replace any descendents (i.e. those that aren't tips) with V2 name
  TransDat_neg_rep$descendents <- as.character(ifelse(!is.na(TransDat_neg_rep$descendents),TransDat_neg_rep$descendents,TransDat_neg_rep$V2name))
  
  
  ################
  #Step that assigns caste number to ancestor and descendant
  TransDat_neg_rep <- data.frame(TransDat_neg_rep, 
                                 ancPoly=Poly$poly2[match(TransDat_neg_rep$ancestors, Poly$animal)],
                                 desPoly=Poly$poly2[match(TransDat_neg_rep$descendents, Poly$animal)])
  
  #create new column with concatenated values from ancPoly and desPoly
  TransDat_neg_rep$trans <- paste0(TransDat_neg_rep$ancPoly, "_to_", TransDat_neg_rep$desPoly)
  TransDat_neg_rep$trans <- as.factor(TransDat_neg_rep$trans)
  
  #Append TransDat_neg_rep to a list
  TransDat_neg_rep_list[[i]] <- TransDat_neg_rep
}
#######


#############
#Calculate summary statistics for the number of transitions of each type over the 400 trees
#############

#For a single dataframe
table(TransDat_neg_rep_list[[1]]$trans)

#For all 400 data frames

# First, apply the table() function to the 'trans' column of each data frame
table_list <- lapply(TransDat_neg_rep_list, function(df) data.frame(table(df$trans)))

##caste_1_to_caste_1
#Calculate the mean frequency of occurrences of 'caste_1_to_caste_1' across the table_list list of dataframes. If the column 'caste_1_to_caste_1' exists in a data frame, extract the frequencies associated with it and compute the mean. If the column doesn't exist in a data frame, assign a frequency of 0 for that data frame.

freq_caste_1_to_caste_1 <- extract_frequency(table_list, "caste_1_to_caste_1")
mean(unlist(freq_caste_1_to_caste_1)) #1003.567
median(unlist(freq_caste_1_to_caste_1)) #1003
range(unlist(freq_caste_1_to_caste_1)) #988 to 1030

##caste_1_to_caste_2
freq_caste_1_to_caste_2 <- extract_frequency(table_list, "caste_1_to_caste_2")
mean(unlist(freq_caste_1_to_caste_2)) #33.13
median(unlist(freq_caste_1_to_caste_2)) #32
range(unlist(freq_caste_1_to_caste_2)) #26 to 49

##caste_1_to_caste_3
freq_caste_1_to_caste_3 <- extract_frequency(table_list, "caste_1_to_caste_3")
mean(unlist(freq_caste_1_to_caste_3)) #5.205
median(unlist(freq_caste_1_to_caste_3)) #5
range(unlist(freq_caste_1_to_caste_3)) #1 to 10

##caste_1_to_caste_4
freq_caste_1_to_caste_4 <- extract_frequency(table_list, "caste_1_to_caste_4")
mean(unlist(freq_caste_1_to_caste_4)) #0.3325
median(unlist(freq_caste_1_to_caste_4)) #0
range(unlist(freq_caste_1_to_caste_4)) #0 to 2

##caste_2_to_caste_2
freq_caste_2_to_caste_2 <- extract_frequency(table_list, "caste_2_to_caste_2")
mean(unlist(freq_caste_2_to_caste_2)) #221.0375
median(unlist(freq_caste_2_to_caste_2)) #222
range(unlist(freq_caste_2_to_caste_2)) #185 to 233

##caste_2_to_caste_3
freq_caste_2_to_caste_3 <- extract_frequency(table_list, "caste_2_to_caste_3")
mean(unlist(freq_caste_2_to_caste_3)) #6.46
median(unlist(freq_caste_2_to_caste_3)) #6
range(unlist(freq_caste_2_to_caste_3)) #2 to 12

##caste_2_to_caste_4
freq_caste_2_to_caste_4 <- extract_frequency(table_list, "caste_2_to_caste_4")
mean(unlist(freq_caste_2_to_caste_4)) #1.24
median(unlist(freq_caste_2_to_caste_4)) #1
range(unlist(freq_caste_2_to_caste_4)) #1 to 4

##caste_3_to_caste_3
freq_caste_3_to_caste_3 <- extract_frequency(table_list, "caste_3_to_caste_3")
mean(unlist(freq_caste_3_to_caste_3)) #11.895
median(unlist(freq_caste_3_to_caste_3)) #12
range(unlist(freq_caste_3_to_caste_3)) #2 to 22

##caste_3_to_caste_4
freq_caste_3_to_caste_4 <- extract_frequency(table_list, "caste_3_to_caste_4")
mean(unlist(freq_caste_3_to_caste_4)) #2.28
median(unlist(freq_caste_3_to_caste_4)) #2
range(unlist(freq_caste_3_to_caste_4)) #1 to 3

##caste_4_to_caste_4
freq_caste_4_to_caste_4 <- extract_frequency(table_list, "caste_4_to_caste_4")
mean(unlist(freq_caste_4_to_caste_4)) #0.4825
median(unlist(freq_caste_4_to_caste_4)) #0
range(unlist(freq_caste_4_to_caste_4)) #0 to 4

##Losses of caste

##caste_2_to_caste_1
freq_caste_2_to_caste_1 <- extract_frequency(table_list, "caste_2_to_caste_1")
mean(unlist(freq_caste_2_to_caste_1)) #14.8875
median(unlist(freq_caste_2_to_caste_1)) #15
range(unlist(freq_caste_2_to_caste_1)) #8 to 22

##caste_3_to_caste_1
freq_caste_3_to_caste_1 <- extract_frequency(table_list, "caste_3_to_caste_1")
mean(unlist(freq_caste_3_to_caste_1)) #0.6525
median(unlist(freq_caste_3_to_caste_1)) #0
range(unlist(freq_caste_3_to_caste_1)) #0 to 5

##caste_3_to_caste_2
freq_caste_3_to_caste_2 <- extract_frequency(table_list, "caste_3_to_caste_2")
mean(unlist(freq_caste_3_to_caste_2)) #0.6425
median(unlist(freq_caste_3_to_caste_2)) #0
range(unlist(freq_caste_3_to_caste_2)) #0 to 4

##caste_4_to_caste_1
freq_caste_4_to_caste_1 <- extract_frequency(table_list, "caste_4_to_caste_1")
mean(unlist(freq_caste_4_to_caste_1)) #0.01
median(unlist(freq_caste_4_to_caste_1)) #0
range(unlist(freq_caste_4_to_caste_1)) #0 to 3

##caste_4_to_caste_2
freq_caste_4_to_caste_2 <- extract_frequency(table_list, "caste_4_to_caste_2")
mean(unlist(freq_caste_4_to_caste_2)) #0.0025
median(unlist(freq_caste_4_to_caste_2)) #0
range(unlist(freq_caste_4_to_caste_2)) #0 to 1

##caste_4_to_caste_3
freq_caste_4_to_caste_3 <- extract_frequency(table_list, "caste_4_to_caste_3")
mean(unlist(freq_caste_4_to_caste_3)) #0.175
median(unlist(freq_caste_4_to_caste_3)) #0
range(unlist(freq_caste_4_to_caste_3)) #0 to 4






#############

#############Plotting ASR with transition nodes###########

#############

#Subset the first dataframe from the TransDat_neg_rep_list list of dataframes for rows with particular transition types
one_to_two <- subset(TransDat_neg_rep_list[[1]], grepl("caste_1_to_caste_2", trans))
one_to_three <- subset(TransDat_neg_rep_list[[1]], grepl("caste_1_to_caste_3", trans))
one_to_four <- subset(TransDat_neg_rep_list[[1]], grepl("caste_1_to_caste_4", trans))

two_to_one <- subset(TransDat_neg_rep_list[[1]], grepl("caste_2_to_caste_1", trans))
two_to_three <- subset(TransDat_neg_rep_list[[1]], grepl("caste_2_to_caste_3", trans))
two_to_four <- subset(TransDat_neg_rep_list[[1]], grepl("caste_2_to_caste_4", trans))

three_to_one <- subset(TransDat_neg_rep_list[[1]], grepl("caste_3_to_caste_1", trans))
three_to_two <- subset(TransDat_neg_rep_list[[1]], grepl("caste_3_to_caste_2", trans))
three_to_four <- subset(TransDat_neg_rep_list[[1]], grepl("caste_3_to_caste_4", trans))

four_to_one <- subset(TransDat_neg_rep_list[[1]], grepl("caste_4_to_caste_1", trans))
four_to_two <- subset(TransDat_neg_rep_list[[1]], grepl("caste_4_to_caste_2", trans))
four_to_three <- subset(TransDat_neg_rep_list[[1]], grepl("caste_4_to_caste_3", trans))

#Load ggtree
library(ggtree)
library(ggtreeExtra)

#Annotating phylogeny to highlight transitional nodes

#Identify transitional nodes where caste number increases
one_to_two_trans <- one_to_two$V2
one_to_three_trans <- one_to_three$V2
one_to_four_trans <- one_to_four$V2
two_to_three_trans <- two_to_three$V2
two_to_four_trans <- two_to_four$V2
three_to_four_trans <- three_to_four$V2

p <- ggtree(ASR_models[[1]]$phy) # + geom_tiplab(size = 0.2)


##########
#Plotting transitions with colour corresponding to the number of castes that the descendant has e.g. if a species ends up with 3 castes (regardless of the number of castes that it started with) its colour is X
##########

###
#Create plot WITH subfamily labels

#Change caste number column name
sxData_rename <- sxData %>% dplyr::rename(`Number of worker castes` = caste.number)
sxData_rename$`Number of worker castes` <- as.factor(sxData_rename$`Number of worker castes`)

#Add node labels to phylogeny - nodes start at 653 and there are 651 of them. Final node is 1303
ASR_models[[1]]$phy$node.label <- c(seq(653,(653+650)))

nicer_tree <- ggtree(ASR_models[[1]]$phy, layout="fan", branch.length = "none", open.angle = 25, size = 0.1) + geom_rootedge(TRUE, size = 0.1) + geom_tiplab(size = 0.9, offset = 2) + geom_nodelab(size = 0.5) #Plot including node and tip labels which can be helpful when plotting subfamily labels

# nicer_tree <- ggtree(ASR_models[[1]]$phy, layout="fan", branch.length = "none", open.angle = 25, size = 0.1) + geom_rootedge(TRUE, size = 0.1)# + geom_tiplab(size = 0.9, offset = 2) + geom_nodelab(size = 0.5)

nicer_tree_rotated <- rotate_tree(nicer_tree, angle=90)

tree_styled <- nicer_tree_rotated +
  geom_fruit(data=sxData_rename, # Data
             geom=geom_tile, # Plots 'tiles' (squares) onto each tip
             width=1.6, #Alters the length of the tiles
             position=position_identityx(hexpand=30), # Adjusts how far away the tiles are from the tre
             mapping=aes(y=animal, fill=`Number of worker castes`), # Analogous to ggplot aes()
             color = "white", # Colour of outline round the tiles (white makes it look like a gap)
             lwd = 0.2, # Width of line between tiles
             linetype = 1, # Default - other numbers make the line dashes, dotted etc.
             axis.params=list( # Add label to the geom_tile - by adding an x-axis
               axis="x",
               text = " ", # Label to plot
               text.size = 3.5, # Size of text
               hjust = 0, # Adjust position of text relative to the geom_tile
               vjust = 0.8,
               text.angle=360,
               line.colour="white"), # Set to white so axis line is not visible
             offset = 4, # Fine-scale adjustment of space between tree & layers
             pwidth=100) + # Width of whole plot
  scale_fill_manual(values=c("#E5F5FF", "#A6CFFA", "#5278B7", "#000005")) +
  theme_tree(legend.position = "none") #This line removes the legend

tree_styled

#Add trasitions from 1 to 2
tree_trans <- tree_styled + ggtree::geom_point2(aes(subset = node %in% one_to_two_trans),
                                                size = 6, colour = '#A6CFFA', shape = 20, alpha = 0.65) #When plotting for multiple nodes, %in% is crucial. we assigning node as == to something. Therefore, to make node == to multiple values, we need the %in% operator.

#Add transitions from 1 to 3
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% one_to_three_trans),
                                               size = 6, colour = '#5278B7', shape = 20, alpha = 0.65)

#Add transitions from 1 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% one_to_four_trans),
                                               size = 6, colour = '#000005', shape = 20, alpha = 0.6)

#Add transitions from 2 to 3
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% two_to_three_trans),
                                               size = 6, colour = '#5278B7', shape = 20, alpha = 0.65)

#Add transitions from 2 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% two_to_four_trans),
                                               size = 6, colour = '#000005', shape = 20, alpha = 0.6)
#Add transitions from 3 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% three_to_four_trans),
                                               size = 6, colour = '#000005', shape = 20, alpha = 0.6)
tree_trans

tree_clad_lab <- tree_trans + geom_cladelab(node=764, label="Pheidole\n(big-headed ants)", align=TRUE, fontsize = 3, angle="auto",
                                            offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=839, label="Fungus-growing \nants (including \nAtta and\nAcromyrmex)", align=TRUE, fontsize = 3, angle="auto",
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=693, label="Carebara", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=982, label="Camponotus &\nCalomyrmex \n(including \ncarpenter ants)", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1063, label="Cataglyphis\n(desert ants)", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1030, label="Formica (wood ants)", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1084, label="Lasius &\nMyrmecocystus", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=923, label="Pogonomyrmex", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=724, label="Temnothorax", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=669, label="Crematogaster", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1206, label="Ponerinae", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1178, label="Dorylinae (army ants)", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')


tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1106, label="Dolichoderinae", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=881, label="Solenopsidini \n(including fire ants)", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=890, label="Stenammini", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=945, label="Myrmicini", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=964, label="Ectatommini", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1011, label="Polyrhachis", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1101, label="Brachymyrmex", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1142, label="Pseudomyrmecini", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1286, label="Amblyoponini", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=690, label="Acanthomyrmex", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=708, label="Tetramorium", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=813, label="Cephalotes \n(turtle ants)", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=819, label="Strumigenys", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1161, label="Myrmeciini", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=879, label="Daceton & \nOrectognathus", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab






################################################
###
#Create plot WITHOUT subfamily labels so that they can be added in PowerPoint

#Change caste number column name
sxData_rename <- sxData %>% dplyr::rename(`Number of worker castes` = caste.number)
sxData_rename$`Number of worker castes` <- as.factor(sxData_rename$`Number of worker castes`)

#Add node labels to phylogeny - nodes start at 653 and there are 651 of them. Final node is 1303
ASR_models[[1]]$phy$node.label <- c(seq(653,(653+650)))

# nicer_tree <- ggtree(ASR_models[[1]]$phy, layout="fan", branch.length = "none", open.angle = 25, size = 0.1) + geom_rootedge(TRUE, size = 0.1) + geom_tiplab(size = 0.9, offset = 2) + geom_nodelab(size = 0.5) #Plot including node and tip labels which can be helpful when plotting subfamily labels

nicer_tree <- ggtree(ASR_models[[1]]$phy, layout="fan", branch.length = "none", open.angle = 25, size = 0.1) + geom_rootedge(TRUE, size = 0.1)# + geom_tiplab(size = 0.9, offset = 2) + geom_nodelab(size = 0.5)

nicer_tree_rotated <- rotate_tree(nicer_tree, angle=90)

tree_styled <- nicer_tree_rotated +
  geom_fruit(data=sxData_rename, # Data
             geom=geom_tile, # Plots 'tiles' (squares) onto each tip
             width=1.6, #Alters the length of the tiles
             position=position_identityx(hexpand=30), # Adjusts how far away the tiles are from the tre
             mapping=aes(y=animal, fill=`Number of worker castes`), # Analogous to ggplot aes()
             color = "white", # Colour of outline round the tiles (white makes it look like a gap)
             lwd = 0.2, # Width of line between tiles
             linetype = 1, # Default - other numbers make the line dashes, dotted etc.
             axis.params=list( # Add label to the geom_tile - by adding an x-axis
               axis="x",
               text = " ", # Label to plot
               text.size = 3.5, # Size of text
               hjust = 0, # Adjust position of text relative to the geom_tile
               vjust = 0.8,
               text.angle=360,
               line.colour="white"), # Set to white so axis line is not visible
             offset = 4, # Fine-scale adjustment of space between tree & layers
             pwidth=100) + # Width of whole plot
  scale_fill_manual(values=c("#F0FAFF", "#A6CFFA", "#4766A8", "black")) +
  theme_tree(legend.position = "none") #This line removes the legend

tree_styled


#Add trasitions from 1 to 2
tree_trans <- tree_styled + ggtree::geom_point2(aes(subset = node %in% one_to_two_trans),
                                                size = 11, colour = '#A6CFFA', shape = 20, alpha = 0.65) #When plotting for multiple nodes, %in% is crucial. we assigning node as == to something. Therefore, to make node == to multiple values, we need the %in% operator.

#Add transitions from 1 to 3
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% one_to_three_trans),
                                               size = 11, colour = '#4766A8', shape = 20, alpha = 0.95)

#Add transitions from 1 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% one_to_four_trans),
                                               size = 11, colour = 'black', shape = 20, alpha = 0.95)

#Add transitions from 2 to 3
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% two_to_three_trans),
                                               size = 11, colour = '#4766A8', shape = 20, alpha = 0.95)

#Add transitions from 2 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% two_to_four_trans),
                                               size = 11, colour = 'black', shape = 20, alpha = 0.95)
#Add transitions from 3 to 4
tree_trans <- tree_trans + ggtree::geom_point2(aes(subset = node %in% three_to_four_trans),
                                               size = 11, colour = 'black', shape = 20, alpha = 0.95)
tree_trans

tree_clad_lab <- tree_trans + geom_cladelab(node=764, label="", align=TRUE, fontsize = 5, angle="auto",
                                            offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=839, label="", align=TRUE, fontsize = 5, angle="auto",
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=693, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=982, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1063, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1030, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1084, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=923, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=724, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=669, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1206, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1178, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')


tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1106, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=881, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=890, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=945, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=964, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1011, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1101, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1142, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1286, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=690, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=708, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=813, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=819, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=1161, label="", align=TRUE, fontsize = 5, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab <- tree_clad_lab + geom_cladelab(node=879, label="", align=TRUE, fontsize = 3, angle='auto',
                                               offset = 2, offset.text = 0.5, textcolor='black', barcolor='black')

tree_clad_lab

ggsave("/Users/louis.bell-roberts/Documents/Github/Testing_the_size_complexity_hypothesis_in_ants/Figures/Caste_transitions/Plots/Phylogeny_caste_discrete_no_labels.pdf", plot = tree_clad_lab, width = 13, height = 13, units = "in", dpi = 300)


####################
# Estimating the ancestral number of worker castes at the root of the phylogeny across 400 trees - run this code after line 200
int_list <- list()

for (i in seq_along(ASR_states_list)) {
  ASR_states <- ASR_states_list[[i]]

  # Combine the probability scores for the two different rate classes for each caste number
  ASR_states_summed <- data.frame(
    caste_1 = ASR_states$`(1,R1)` + ASR_states$`(1,R2)`,
    caste_2 = ASR_states$`(2,R1)` + ASR_states$`(2,R2)`,
    caste_3 = ASR_states$`(3,R1)` + ASR_states$`(3,R2)`,
    caste_4 = ASR_states$`(4,R1)` + ASR_states$`(4,R2)`
  )

  int_list[[i]] <- ASR_states_summed
}

#Create a vector that is all of the probability estimates across the 400 trees that the MRCA for all ants had a single worker caste
values_vector <- sapply(int_list, function(df) df[1, 1])
mean(values_vector) #0.9930264%
