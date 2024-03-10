########ASR with MCMCglmm using ASR from corHMM######
#Louis Bell-Roberts
#29/01/2024

.libPaths(c(.libPaths(), "/drives/4tb/modules/R"))

analysis <- "CS_Caste"


setwd(file.path("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/ASR/Model_outputs/MCMCglmm", paste0(analysis), "1st_run/"))

# packages
library(ape)
library(coda)
library(MCMCglmm)
library(boot)
library(phytools)
library(tidyverse)
library(gridExtra)
library(doParallel)


registerDoParallel(50)


################
#Custom functions
################

#Calculate the most probable state for each node in the corHMM ASR
##Define function to find index of maximum value in a vector using the max_column function
max_column <- function(x) {
  return(names(x)[which.max(x)])
}

# Define a function to assign monomorphic or polymorphic to each node in my corHMM ASR - this function is for when using corHMM ASRs with 2 rate categories
assign_node_state <- function(df) {
  df$node_state <- ifelse(df$MaxColumn %in% c("X.1.R1.", "X.1.R2."), "monomorphic", "polymorphic")
  return(df)
}

# Define function to rename columns as poly2 and animal
rename_cols <- function(df) {
  names(df)[1] <- "poly2"
  names(df)[2] <- "animal"
  return(df)
}


################################################
#Data preparation

#Read in data file
ant_data <- read.csv("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/Data/Trait_database/ant_data.csv")

#Set variables so that they're in the correct structure and apply transformations
class(ant_data$colony.size) #numeric
class(ant_data$queen.mating.frequency) #numeric
class(ant_data$caste.number) #integer

#Transform caste.number into a binary variable
ant_data$CasteBin <- ifelse(ant_data$caste.number > 1,"polymorphic","monomorphic")
ant_data$CasteBin <- as.factor(ant_data$CasteBin)

#Rename 'species' column as 'animal'
ant_data <- ant_data %>% dplyr::rename(animal = species)

#Read in sample of 400 phylogenetic trees
ant_trees <- read.tree(file ="/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/Data/15k_Economo_trees/Economo_2018_400.tre")

#Filter data
sxData <- dplyr::filter(ant_data, caste.number >=1, colony.size >=1)
sxData <- dplyr::select(sxData, animal, colony.size, CasteBin)

#Prune trees
ant_trees_pruned <- drop.tip.multiPhylo(ant_trees, setdiff(ant_trees[[1]]$tip.label, sxData$animal))

#Add node labels as they were not originally included with the phylo object
##Edges and node labels used later to work out where transitions in the number of castes happens
for (i in 1:length(ant_trees_pruned)) {
  ant_trees_pruned[[i]]$node.label <- paste("Node", 1:ant_trees_pruned[[i]]$Nnode, sep = "")
}

#Prune database
sxData <- dplyr::filter(sxData, animal %in% ant_trees_pruned[[1]]$tip.label)

#Check that species names in the tree and the data match
sxData$animal[which((sxData$animal %in% ant_trees_pruned[[1]]$tip.label) == FALSE)]	# 0 means all match
ant_trees_pruned[[1]]$tip.label[which((ant_trees_pruned[[1]]$tip.label %in% sxData$animal) == FALSE)]	# 0 means all match



##########################################################################################
#ANCESTRAL STATE RECONSTRUCTION WITH COLONY SIZE AND CASTE NUMBER
##########################################################################################

###############
#Read in the 400 corHMM model results
##corHMM files must be read into R in a particular order - the particular tree used for each corHMM model must correspond to the tree used for each MCMCglmm model later on
###Create an empty list to hold the results
corHMM_results <- list()

#Set folder path
folder_path <- paste0("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/ASR/Model_outputs/corHMM/",analysis,"/")

#Get the file names
file_names <- list.files(folder_path)

#Extract the numbers from the file names using regular expressions
model_numbers <- as.numeric(gsub(".*_(\\d+)\\.rda", "\\1", file_names))

#Sort the file names based on the corresponding numbers
sorted_file_names <- file_names[order(model_numbers)]

#Read in the files and save them to a list
corHMM_results <- lapply(sorted_file_names, function(x) readRDS(file.path(folder_path, x)))

#Extracting the states component from each model and create a new list of just the states
corHMM_states <- lapply(corHMM_results, function(model) model$states)


###############
#For each tree, estimate the most likely number of castes (monomorphic/polymorphic) for each node of the phylogeny

#First, convert each of the "matrix" "array" objects that are in the corHMM_states list to data frame objects before running the subsequent functions
##Apply function to each matrix in the list
corHMM_states_conv <- lapply(corHMM_states, data.frame)

#Apply max_column custom function to each row of the data frame to find the most likely state at each node in the corHMM ASRs
##Apply function to each data frame in the list
corHMM_states_conv <- lapply(corHMM_states_conv, function(df) {
  df$MaxColumn <- apply(df, 1, max_column)
  return(df)
})

###############
#Classify nodes as monomorphic if in state X.1.R1 or X.1.R2, and polymorphic if in state X.2.R1 or X.2.R2
##Apply the function to each dataframe in the list
corHMM_states_conv <- lapply(corHMM_states_conv, assign_node_state)


###############
#Data preparation
##Create "animal.Node" number column pasted before the node number for each row
for (i in seq_along(corHMM_states_conv)) {
  corHMM_states_conv[[i]]$node_number <- paste0("animal.Node", seq_len(nrow(corHMM_states_conv[[i]])))
}

#Replace "animal.Node1" with "(Intercept)" in each data frame in the list
for (i in seq_along(corHMM_states_conv)) {
  corHMM_states_conv[[i]]$node_number <- gsub("\\banimal\\.Node1\\b", "(Intercept)", corHMM_states_conv[[i]]$node_number)
}

#Create new list of data frames extracting just the node_state and node_number columns of 'corHMM_states_conv'
corHMM_states_conv_minimal <- lapply(corHMM_states_conv, function(df) df[, c("node_state", "node_number")])

#Set row names to match the values in node_number column
corHMM_states_conv_minimal <- lapply(corHMM_states_conv_minimal, function(df) {
  rownames(df) <- df[, "node_number"]
  return(df)
})

#Get data for each species tip and its caste data
##Create 'nodecode' column" for sxData by pasting animal and then also the species tip
sxData$nodecode <- paste("animal", sxData$animal,sep=".")

#Create new data frame extracting just the CasteBin and nodecode columns
sxData_minimal <- sxData[, c("CasteBin", "nodecode")]

#Make the row names equal to nodecode column
rownames(sxData_minimal) <- sxData_minimal[, "nodecode"]

#The column names of corHMM_states_conv_minimal and sxData_minimal must match for rbind to work
##Use lapply to change column names for each data frame in the list of corHMM_states_conv_minimal using the rename_cols custom function
corHMM_states_conv_minimal <- lapply(corHMM_states_conv_minimal, rename_cols)

#Change column names using names function for sxData_minimal
names(sxData_minimal) <- c("poly2", "animal")

#Append the rows from the 'corHMM_states_conv_minimal' data frame to the 'sxData_minimal' data frame using the rbind operation
##Use lapply to bind sxData_minimal to each element in corHMM_states_conv_minimal
Poly <- lapply(corHMM_states_conv_minimal, function(df) {
  rbind(df, sxData_minimal)
})


#############
#Create Transition Dataset
##Use the structure of the phylogeny to figure out how nodes are connected
#Extract edge information for each tree in the multiPhylo object and create a list of dataframes containing information on the relationship between parent and offspring nodes across the tree
TransDat <- lapply(ant_trees_pruned, function(x) as.data.frame(x$edge)) ##Each value appears twice in column V1 as nodes have two descendants

#Subtract tips and create column 'V1name' for each data frame in the TransDat list to identify node names for each row (V1 = parent nodes)
for (i in seq_along(TransDat)) {
  TransDat[[i]]$V1name <- paste("animal.Node", TransDat[[i]]$V1-length(ant_trees_pruned[[1]]$tip.label), sep = "")
}

#Repeat for V2 (offspring nodes/tips)
for (i in seq_along(TransDat)) {
  TransDat[[i]]$V2name <- paste("animal.Node", TransDat[[i]]$V2-length(ant_trees_pruned[[1]]$tip.label), sep = "")
} #'V2name' will include some negative numbers since some of these are tips rather than nodes. Swap these for the correct tip labels later on

#Assign tip number to each species and paste 'animal' in front of it to match with row names in the 'Poly' dataframe created above
treesp <- data.frame(tip.label=paste("animal.", ant_trees_pruned[[1]]$tip.label,sep=""), no=1:length(ant_trees_pruned[[1]]$tip.label)) #While tree topology changes across the sample of 400 trees, the tip number that a species is labelled with always stays the same even if its position on the tree moves

#Replace negative values with tip labels by matching V2 in TransDat with 'no' in treesp - creates the 'descendants' column which contains this information
TransDat_neg_rep <- list()
for (i in seq_along(TransDat)) {
  TransDat_neg_rep[[i]] <- data.frame(TransDat[[i]], ancestors=NA, descendants=treesp$tip.label[match(TransDat[[i]]$V2, treesp$no)])
}

#Assign class 'character' to 'ancestors' and 'descendants'
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(df) {df$ancestors <- as.character(df$ancestors); return(df)})
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(df) {df$descendants <- as.character(df$descendants); return(df)})

#Values for ancestors are currently NAs, so assign V1name 
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(df) {
  df$ancestors <- as.character(ifelse(!is.na(df$ancestors), df$ancestors, df$V1name))
  return(df)
})

#Replace any descendants (i.e. those that aren't tips) with V2 name
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(df) {
  df$descendants <- as.character(ifelse(!is.na(df$descendants), df$descendants, df$V2name))
  return(df)
}) #The names of each ancestor and descendant are now identified. Since they have 'animal.' pasted in front, they will match with row names in the 'Poly' dataframe.


#############
#Combine colony size estimates with transition dataset
##Match up ancestor and descendant colony size predictions based on row names
TransDat_neg_rep <- lapply(seq_along(TransDat_neg_rep), function(i) {
  data.frame(TransDat_neg_rep[[i]], 
             ancPoly=Poly[[i]]$poly2[match(TransDat_neg_rep[[i]]$ancestors,Poly[[i]]$animal)],
             desPoly=Poly[[i]]$poly2[match(TransDat_neg_rep[[i]]$descendants,Poly[[i]]$animal)])
})

#############
#Calculate the number of different types of transitions between caste number
##Make a table of the number of types of descendant each ancestor has
###Possible transition types:
# 1 polymorphic and 1 monomorphic descendant, 2 poly vs. 0 mono, 0 poly vs. 2 mono
obs_list <- lapply(TransDat_neg_rep, function(df) data.frame(table(df$desPoly, df$ancestors)))
obs_list <- lapply(obs_list, function(df){
  df$CAT <- df$Freq
  return(df)
})

#Assign transition category types (CAT) based on the number of castes that each descendant has
for(i in seq_along(obs_list)) {
  obs <- obs_list[[i]]
  
  # Modify the 'CAT' column based on conditions
  obs$CAT[obs$Freq == 2 & obs$Var1 == "monomorphic"] <- "only monomorphic"
  obs$CAT[obs$Freq == 0 & obs$Var1 == "polymorphic"] <- "only monomorphic"
  obs$CAT[obs$Freq == 2 & obs$Var1 == "polymorphic"] <- "only polymorphic"
  obs$CAT[obs$Freq == 0 & obs$Var1 == "monomorphic"] <- "only polymorphic"
  obs$CAT[obs$Freq == 1] <- "both"
  
  # Assign the modified data frame back to the list
  obs_list[[i]] <- obs
}

#############
#Assign Types of Transitions Between caste number to each Node

#Remove NAs in the transitions dataset
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(x) {
  x <- x[!is.na(x$ancPoly),] #This deletes node 1 which doesn't have an ancestor
  x <- x[!is.na(x$desPoly),] #Deletes rows where caste number is uncertain
  return(x)
})

#Assign transition types in TransDat based on the values calculated in obs_list
for (i in seq_along(TransDat_neg_rep)) {
  # match ancestors with corresponding CAT value from obs_list
  TransDat_neg_rep[[i]] <- data.frame(TransDat_neg_rep[[i]], CAT = obs_list[[i]]$CAT[match(TransDat_neg_rep[[i]]$ancestors, obs_list[[i]]$Var2)])
}

#Create a new column that classifies nodes based on ancestor state & transition type. 
for (i in seq_along(TransDat_neg_rep)) {
  # create a new column in each data frame and concatenate ancPoly and CAT columns
  TransDat_neg_rep[[i]]$CAT2 <- paste(TransDat_neg_rep[[i]]$ancPoly, TransDat_neg_rep[[i]]$CAT, sep = ".")
}

#This can produce up to 6 transition type (although, not all 6 will necessarily happen)
table(TransDat_neg_rep[[1]]$CAT2)

#Assign 'double transitions' (e.g. monomorphic ancestor to only polymorphic descendants) as only single transition events - therefore, there will only be 4 transition types
for (i in seq_along(TransDat_neg_rep)) {
  TransDat_neg_rep[[i]]$CAT2[TransDat_neg_rep[[i]]$CAT2 == "polymorphic.only monomorphic"] <- "polymorphic.both"
  TransDat_neg_rep[[i]]$CAT2[TransDat_neg_rep[[i]]$CAT2 == "monomorphic.only polymorphic"] <- "monomorphic.both"
} #To identify the number of transitions in the dataset, divide the total number by 2. This is because each node is represented twice twice as an ancestor has two descendants


#############
#Add estimates for the continuous trait for the extant species to the transition datasets

#Match up tips with colony size data based on nodecode column
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(x) {
  x <- data.frame(x, continuous.trait = sxData$colony.size[match(x$descendants, sxData$nodecode)])
  return(x)
})

#############
#Assign ancestral nodes to the 'animal' column
#When running the model, this will allow us to estimate colony size in 'ancestors'
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(df) {
  df$animal <- gsub("animal.", "", df$ancestors)
  return(df)
})

#Assign transition types as a factor
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(df) {
  df$CAT2 <- as.factor(df$CAT2)
  return(df)
})

#Remove rownames from each data frame in the list
TransDat_neg_rep <- lapply(TransDat_neg_rep, function(x) {
  rownames(x) <- NULL
  return(x)
})

#############
#Model to estimate colony size values across each transition type

#Prior for model with Gaussian response variable
prior1 <- list(R = list(V = 1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 0.002)))

# Parallel models with outputs saved as .rds files. Remove the global intercept and run over 400 trees.
foreach(i = 1:400) %dopar% {
  #1st model
  model1 <- MCMCglmm(log10(continuous.trait) ~ CAT2-1, random = ~animal, family = "gaussian", nodes = "ALL", prior = prior1, pedigree = ant_trees_pruned[[i]], data = TransDat_neg_rep[[i]], nitt = 1100000, burnin = 100000, thin = 1000, pr=TRUE, verbose = F)
  #Run 2nd model
  model2 <- MCMCglmm(log10(continuous.trait) ~ CAT2-1, random = ~animal, family = "gaussian", nodes = "ALL", prior = prior1, pedigree = ant_trees_pruned[[i]], data = TransDat_neg_rep[[i]], nitt = 1100000, burnin = 100000, thin = 1000, pr=TRUE, verbose = F)
  
  # Save the model as an .rds file for 1st and 2nd chain
  saveRDS(model1, file = file.path("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/ASR/Model_outputs/MCMCglmm", analysis, "1st_run", paste0(analysis, "_1M_100k_1k_", i, "_1stRun.rds")))
  saveRDS(model2, file = file.path("/drives/4tb/Louis/Worker_polymorphism/Post_review_analysis/ASR/Model_outputs/MCMCglmm", analysis, "2nd_run", paste0(analysis, "_1M_100k_1k_", i, "_2ndRun.rds")))
}

