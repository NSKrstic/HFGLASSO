#################################################
# Script: HCP_creation_SCAL.R
# Author: Nikolas Krstic
# Purpose: Create a hierarchical categorical predictor and its objects (design matrix, S matrices, etc.)
#################################################


source("HFG_LASSO_func.R")

library(MASS)

#################################################
# Generate Hierarchical Categorical Dataset

###Total number of layers of hierarchical categorical predictor
Categ_Num = 4

###Number of levels within each layer of the hierarchical categorical predictors (starting with the least granular)
if(Sim_Scal_Type == "SS"){
  Hier_Categ_Levels = c(1, 2, 5, 20)
}else{
  Hier_Categ_Levels = c(1, 5, 30, Cardinality)
}

print(paste("CURRENT SEED: ", Seed, sep=""))
set.seed(Seed)

#Store the number of categorical levels for each categorical variable
Categ_Counts = c()
#Store the indices of the categorical variables as a list
Cat_Var_Col_Indices = list()


#Create Hierarchical Categorical Predictors

for(i in 1:Categ_Num){
  
  #Record the indices and level count for each layer in hierarchy
  Categ_Count = Hier_Categ_Levels[i]
  Cat_Var_Col_Indices[[i]] = sum(Categ_Counts)+(1:(Categ_Count))
  Categ_Counts = c(Categ_Counts, Categ_Count)
  
  #If it's the final/bottom layer, generate the design matrix
  if(i == Categ_Num){
    
    #Randomly generate the probability of each class of the variable
    Probs = runif(Categ_Count, min=0.1, max=0.9)
    
    #Standardize probabilities
    Probs = Probs/sum(Probs)
    
    #Create the design matrix
    DesignMatrix = t(rmultinom(Samp_Size, 1, prob=Probs))[,1:Categ_Count]
    Valid_DesignMatrix = t(rmultinom(Valid_Samp_Size, 1, prob=Probs))[,1:Categ_Count]
    Test_DesignMatrix = t(rmultinom(Test_Samp_Size, 1, prob=Probs))[,1:Categ_Count]
  }
  
}


#Build random nested tree structure (i.e. assigning a parent node to each lower layer node, starting with the layer right below the root node)

Assignment_List = list()

for(j in 1:(length(Categ_Counts)-1)){
  
  #Identify the indices for each of the two layers (the upper layer and the lower layer)
  Curr_Upp_Layer_Inds = Cat_Var_Col_Indices[[j]]
  Curr_Low_Layer_Inds = Cat_Var_Col_Indices[[j+1]]
  
  #Assign the first few lower layer nodes to one of each of the upper layer nodes, then randomly assign the remaining lower layer nodes
  # (e.g. If there are two parent nodes, then the first two lower layer nodes should be individually assigned to each of these parent nodes)
  # The reason this is done is to ensure a complete tree is made (rather than having disconnected subtrees)
  Curr_Assignment = c(Curr_Upp_Layer_Inds, sample(Curr_Upp_Layer_Inds, length(Curr_Low_Layer_Inds)-length(Curr_Upp_Layer_Inds), replace=TRUE))
  
  #Name each assignment based on the corresponding lower layer node, then save in the list
  names(Curr_Assignment) = paste("Node: ", Curr_Low_Layer_Inds, sep="")
  Assignment_List[[j]] = Curr_Assignment
  
}

#The final assignment list contains the list of parents for each node (the corresponding child node is the name of the vector element)



############################################################################################################
#Identify all tree leave paths (i.e. the sequence of nodes starting from the root to reach each leaf node):

Path_List = list()

#Iterate over each leaf node
for(k in 1:Categ_Counts[length(Categ_Counts)]){
  
  #Identify current leaf node
  Curr_Index = Cat_Var_Col_Indices[[length(Cat_Var_Col_Indices)]][k]
  
  #Initiate path
  Curr_Path = c(Curr_Index)
  
  for(l in length(Assignment_List):1){
    
    #Extract parent of current index
    Parent_Index = as.numeric(Assignment_List[[l]][names(Assignment_List[[l]])==paste("Node: ", Curr_Index, sep="")])
    
    #Add to current path
    Curr_Path = c(Parent_Index, Curr_Path)
    
    #Assign new current index as the parent index just identified
    Curr_Index = Parent_Index
    
  }
  
  #Add the current path to the list
  Path_List[[k]] = Curr_Path
  
}


#####################################################################################################
# Randomly generate the active nodes in the tree

Ancestor_Count = sum(Categ_Counts[1:(length(Categ_Counts)-1)])

#Identify which ancestor nodes in the tree are active
if(Sim_Scal_Type=="SS"){
  
  Active_Nodes_Set = sample(2:Ancestor_Count, size=2)
  
  #Remove active nodes that are simply descendants of other active nodes
  if(any(c(2,3) %in% Active_Nodes_Set)){
    Rem_Nodes = substring(names(Assignment_List[[2]][which(Assignment_List[[2]] %in% Active_Nodes_Set)]), 7)
    Active_Nodes_Set = Active_Nodes_Set[which(!Active_Nodes_Set %in% Rem_Nodes)]
  }
  
}else if(Sim_Scal_Type=="HCPC"){
  
  Active_Nodes_Set = sample(2:Ancestor_Count, size=8)
  
  #Remove active nodes that are simply descendants of other active nodes
  if(any(c(2,3) %in% Active_Nodes_Set)){
    Rem_Nodes = substring(names(Assignment_List[[2]][which(Assignment_List[[2]] %in% Active_Nodes_Set)]), 7)
    Active_Nodes_Set = Active_Nodes_Set[which(!Active_Nodes_Set %in% Rem_Nodes)]
  }
  
}


#Identify which leaves in the tree are inactive because of the active ancestor nodes
Inactive_Leaves_Sets = lapply(Active_Nodes_Set, function(y){unlist(sapply(Path_List, function(x){if(any(y %in% x)){x[length(x)]}}))})

#Identify 
##Curr_Inactive_Nodes_Set = sort(unique(unlist(sapply(Path_List, function(x){if(any(Active_Nodes_Set %in% x)){x}}))))


#################################################################################################################################################
#### Make True Model Matrix


#Identify the active leaves that "survived" the merging above due to active nodes
Active_Leaves = which(!Cat_Var_Col_Indices[[length(Cat_Var_Col_Indices)]] %in% unlist(Inactive_Leaves_Sets))
Ancestor_Node_Count = (sum(Categ_Counts)-Categ_Counts[length(Categ_Counts)])

#Construct the "merged" binary indicators based on the union of the "inactive leaves sets" (with the indices adjusted to start from 1)
if(length(Active_Leaves) != Categ_Counts[length(Categ_Counts)]){
  
  True_Inactive_Leaves_Sets = lapply(Inactive_Leaves_Sets, function(x){x-Ancestor_Node_Count})
  Merged_Binary_Indicators = do.call(cbind, lapply(Inactive_Leaves_Sets, function(x){if(length(x)>1){as.numeric(rowSums(DesignMatrix[,x-Ancestor_Node_Count])>0)}else{as.numeric(DesignMatrix[,x-Ancestor_Node_Count]>0)}}))
  Merged_VALID_Binary_Indicators = do.call(cbind, lapply(Inactive_Leaves_Sets, function(x){if(length(x)>1){as.numeric(rowSums(Valid_DesignMatrix[,x-Ancestor_Node_Count])>0)}else{as.numeric(Valid_DesignMatrix[,x-Ancestor_Node_Count]>0)}}))
  Merged_TEST_Binary_Indicators = do.call(cbind, lapply(Inactive_Leaves_Sets, function(x){if(length(x)>1){as.numeric(rowSums(Test_DesignMatrix[,x-Ancestor_Node_Count])>0)}else{as.numeric(Test_DesignMatrix[,x-Ancestor_Node_Count]>0)}}))
  
  #Construct the "true" design matrix
  Final_True_Design_Matrix = cbind(DesignMatrix[, c(Active_Leaves)], Merged_Binary_Indicators)
  True_Valid_Design_Matrix = cbind(Valid_DesignMatrix[, c(Active_Leaves)], Merged_VALID_Binary_Indicators)
  True_Test_Design_Matrix = cbind(Test_DesignMatrix[, c(Active_Leaves)], Merged_TEST_Binary_Indicators)
  
}else{
  
  Final_True_Design_Matrix = DesignMatrix
  True_Valid_Design_Matrix = Valid_DesignMatrix
  True_Test_Design_Matrix = Test_DesignMatrix
  
}


##################################################################################################################################################################################
# Lastly, construct the S matrices (i.e. the matrices that will be used in the regularization term that dictate the fusion structure of the leaves at each parent (ancestor) node)

#Initialize to store S Matrix Layer Memberships
S_Matrix_Layers = c()

#Identify the set of nodes that are parent nodes (have at least one descendant), since there will be one S matrix per parent (ancestor) node
Parent_Node_Set = unlist(Cat_Var_Col_Indices[1:(length(Cat_Var_Col_Indices)-1)])

#Storage list for the S Matrices (just the fusion pairs for each row, to save on memory)
Categ_Pair_Inds = list()

#Iterate over the parent (ancestor) nodes
for(s in 1:length(Parent_Node_Set)){
  
  #Storage list for current S Matrix
  Categ_Pairs_Total = c()
  
  #Identify the leaf paths in which the parent node can be found
  Desc_Leaf_Paths = lapply(Path_List, function(x){if((s %in% x)){x}})
  Desc_Leaf_Paths = Desc_Leaf_Paths[lengths(Desc_Leaf_Paths) != 0]
  
  #Identify the child nodes of this current parent node
  Children_Nodes = unique(sapply(Desc_Leaf_Paths, function(x){x[which(x %in% s)+1]}))
  
  #Store the leaf nodes for each of the children nodes
  Children_Leaf_Sets = list()
  
  for(j in 1:length(Children_Nodes)){
    
    #Identify the current child node
    Curr_Child_Node = Children_Nodes[j]
    
    #Identify the leaf nodes connected to this current child node
    Child_Leaf_Set = unlist(sapply(Desc_Leaf_Paths, function(x){if((Curr_Child_Node %in% x)){x[length(x)]}}))-(sum(Categ_Counts)-Categ_Counts[length(Categ_Counts)])
    
    #Store the leaf node set
    Children_Leaf_Sets[[j]] = Child_Leaf_Set
    
  }
  
  #If the number of children nodes is greater than 2, then create an S matrix, otherwise do not
  if(length(Children_Leaf_Sets)>=2){
    
    #Identify all possible pairwise combinations of the child nodes
    Pair_Child_Combs = combn(1:length(Children_Leaf_Sets), 2)
    
    #For each combination, compute the corresponding rows of the S matrix for the current parent node
    for(k in 1:ncol(Pair_Child_Combs)){
      
      #Identify the current combination as well as the currently selected pair of child node leaf sets
      Curr_Comb = Pair_Child_Combs[,k]
      Child_1_Set = Children_Leaf_Sets[[Curr_Comb[1]]]
      Child_2_Set = Children_Leaf_Sets[[Curr_Comb[2]]]
      
      #Generate the fusion pairs for the current pair of child node leaf sets
      First_Coords = rep(Child_1_Set, times=length(Child_2_Set))
      Second_Coords = rep(Child_2_Set, each=length(Child_1_Set))
      
      Categ_Pairs = cbind(First_Coords, Second_Coords)
      Categ_Pairs_Total = rbind(Categ_Pairs_Total, Categ_Pairs)
      
    }
    
    #Record the fusion pairs
    Categ_Pair_Inds[[s]] = Categ_Pairs_Total
    #Record the tree layer that each S Matrix corresponds to
    S_Matrix_Layers = c(S_Matrix_Layers, which(sapply(Cat_Var_Col_Indices, function(x){s %in% x})))
    
  }else{
    
    #Assign NULL for an S Matrix that won't be created for the current parent node
    Categ_Pair_Inds[[s]] = NULL
  }
  
}


#Compile the sets of fusion pair indices (both pairs, as well as the categories involved in each S matrix)
Categ_Pair_Inds = Filter(Negate(is.null), Categ_Pair_Inds)
Categ_Inds = lapply(Categ_Pair_Inds, function(x){sort(unique(as.vector(x)))})

