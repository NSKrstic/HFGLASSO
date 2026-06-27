#################################################
# Script: HFGL_US_Colleges_Application_DataCleaning.R
# Author: Nikolas Krstic
# Purpose: Apply HFG LASSO to the USA Colleges Dataset (see Cerda and Varoquaux (2020))
# Paper with Dataset:: Cerda, P., & Varoquaux, G. (2020). Encoding high-cardinality string categorical variables. 
## IEEE Transactions on Knowledge and Data Engineering, 34(3), 1164-1176.
#################################################

library(tidyverse)
library(gridExtra)

source("HFG_LASSO_func.R")

#Region Assignments
Southeast = c("AL", "AR", "FL", "GA", "KY", "LA", "MS", "NC", "SC", "TN", "VA", "WV")
FarWest = c("AK", "CA", "HI", "NV", "OR", "WA")
Southwest = c("AZ", "NM", "OK", "TX")
Plains = c("IA", "KS", "MN", "MO", "NE", "ND", "SD")
RockyMountains = c("CO", "ID", "MT", "UT", "WY")
NewEngland = c("CT", "ME", "MA", "NH", "RI", "VT")
MidEast = c("DE", "DC", "MD", "NJ", "NY", "PA")
GreatLakes = c("IL", "IN", "MI", "OH", "WI")

#Supraregion Assignments
West = c("Far West", "Southwest", "Rocky Mountains")
Central = c("Plains", "Great Lakes")
East = c("Southeast", "New England", "Mid East")



#Outcome to use: Percent Pell Grant
College_Data = read_tsv("Data/VueBigData_DataFiles_Colleges.txt")

College_Data = as.data.frame(College_Data)

# Data Cleaning (remove missing outcome observations, unusual categories, etc.)

College_Data = College_Data[which(!College_Data$Region %in% c("US Services Schools", "Outlying Areas (AS, FM,GU, MH, MP, PR, PW, VI)")),]
College_Data = College_Data[which(!is.na(College_Data$`Percent Pell Grant`)),]


### Properly reassign observations to their proper region (some are NA for some strange reason, even though the state is known)
College_Data[which(College_Data$State %in% Southeast),]$Region = "Southeast"
College_Data[which(College_Data$State %in% FarWest),]$Region = "Far West"
College_Data[which(College_Data$State %in% Southwest),]$Region = "Southwest"
College_Data[which(College_Data$State %in% Plains),]$Region = "Plains"
College_Data[which(College_Data$State %in% RockyMountains),]$Region = "Rocky Mountains"
College_Data[which(College_Data$State %in% NewEngland),]$Region = "New England"
College_Data[which(College_Data$State %in% MidEast),]$Region = "Mid East"
College_Data[which(College_Data$State %in% GreatLakes),]$Region = "Great Lakes"


# Make them factors and collect all states in a vector
College_Data$State = as.factor(College_Data$State)
All_States = as.character(unique(College_Data$State))

College_Data$Region = as.factor(College_Data$Region)


########################################################################

# Identify HCP Layers (provide list of layers from highest to lowest)

### Generate a "super-region" variable, merging together regions into "West", "Central" and "East" US

College_Data$SuperRegion = NA
College_Data[which(College_Data$Region %in% c("Southeast", "New England", "Mid East")),]$SuperRegion = "East"
College_Data[which(College_Data$Region %in% c("Far West", "Southwest", "Rocky Mountains")),]$SuperRegion = "West"
College_Data[which(College_Data$Region %in% c("Plains", "Great Lakes")),]$SuperRegion = "Central"

College_Data$SuperRegion = as.factor(College_Data$SuperRegion)
College_Data$USA = as.factor("USA")

#Create HCP Tree (List of the layers of nodes)
HCP_Geography = list(unique(College_Data$USA), unique(College_Data$SuperRegion), unique(College_Data$Region), unique(College_Data$State))

#Create Design Matrix
DesignMat_Geog_HCP = model.matrix(~0+State, data=College_Data)
colnames(DesignMat_Geog_HCP) = gsub("State", "", colnames(DesignMat_Geog_HCP))


# First Step, identify unique combinations of the categories, to form the "Path List" to each leaf node
Path_List = unique(College_Data[,c("USA", "SuperRegion", "Region", "State")])


#####################################################################

# Second Step, generate S Matrix Pairs (Categ_Pair_Inds and Categ_Inds) based on hierarchical tree

#Initialize list to store S matrix layers
S_Matrix_Layers = c()

#Identify the set of nodes that are parent nodes (have at least one descendent), since there will be one S matrix per parent node
Parent_Node_Set = unlist(HCP_Geography[c(1,2,3)])

Categ_Pair_Inds = list()

#Iterate over the parent nodes
for(s in 1:length(Parent_Node_Set)){
  
  Categ_Pairs_Total = c()
  
  Current_Parent_Node = Parent_Node_Set[s]
  
  #Identify the leaf paths in which the parent node can be found
  Desc_Leaf_Paths = apply(Path_List, 1, function(x){if((Current_Parent_Node %in% x)){x}})
  
  if(is.list(Desc_Leaf_Paths)){
    Desc_Leaf_Paths = do.call("rbind", Desc_Leaf_Paths[lengths(Desc_Leaf_Paths) != 0])
  }else{
    Desc_Leaf_Paths = t(Desc_Leaf_Paths)
  }
  
  #Identify the child nodes of this current parent node
  Children_Nodes = unique(apply(Desc_Leaf_Paths, 1, function(x){x[which(x %in% Current_Parent_Node)+1]}))
  
  #Store the leaf nodes for each of the children nodes
  Children_Leaf_Sets = list()
  
  for(j in 1:length(Children_Nodes)){
    
    #Identify the current child node
    Curr_Child_Node = Children_Nodes[j]
    
    #Identify the leaf nodes connected to this current child node
    Child_Leaf_Set = as.vector(unlist(apply(Desc_Leaf_Paths, 1, function(x){if((Curr_Child_Node %in% x)){x[length(x)]}})))
    
    #Store the leaf node set
    Children_Leaf_Sets[[j]] = Child_Leaf_Set
    
  }
  
  #If the number of children nodes is greater than 2, then create an S matrix, otherwise do not
  if(length(Children_Leaf_Sets)>=2){
    
    #Identify all possible pairwise combinations of the child nodes
    Pair_Child_Combs = combn(1:length(Children_Leaf_Sets), 2)
    
    #For each combination, compute the corresponding row of the S matrix for the current parent node
    for(k in 1:ncol(Pair_Child_Combs)){
      
      #Identify the current combination as well as the currently selected pair of child node leaf sets
      Curr_Comb = Pair_Child_Combs[,k]
      Child_1_Set = Children_Leaf_Sets[[Curr_Comb[1]]]
      Child_2_Set = Children_Leaf_Sets[[Curr_Comb[2]]]
      
      #Record all of the fusion pairs in the S Matrix for these two pairs of child nodes (works for second to last layer, and any higher layers)
      
      First_Coords = rep(Child_1_Set, times=length(Child_2_Set))
      Second_Coords = rep(Child_2_Set, each=length(Child_1_Set))
      
      Categ_Pairs = cbind(First_Coords, Second_Coords)
      Categ_Pairs_Total = rbind(Categ_Pairs_Total, Categ_Pairs)
      
    }
    
    Categ_Pair_Inds[[s]] = Categ_Pairs_Total
    S_Matrix_Layers = c(S_Matrix_Layers, which(sapply(HCP_Geography, function(x){Current_Parent_Node %in% x})))
    
  }else{
    
    Categ_Pair_Inds[[s]] = NULL
    
  }
  
}

#Compile the sets of fusion pair indices (both pairs, as well as the categories involved in each S matrix)
Categ_Pair_Inds = Filter(Negate(is.null), Categ_Pair_Inds)
Categ_Inds = lapply(Categ_Pair_Inds, function(x){sort(unique(as.vector(x)))})

Categ_Pair_Inds = lapply(Categ_Pair_Inds, function(x){t(apply(x, 1, function(y){which(colnames(DesignMat_Geog_HCP) %in% y)}))})
Categ_Inds = lapply(Categ_Inds, function(x){which(colnames(DesignMat_Geog_HCP) %in% x)})

#Create the S Matrices (to use if necessary, such as for computing Lambda_Min)
S_Matrices = lapply(Categ_Pair_Inds, function(x){t(apply(x, 1, function(y){S_Row=rep(0, ncol(DesignMat_Geog_HCP)); S_Row[y[1]]=1; S_Row[y[2]]=-1; return(S_Row)}))})


