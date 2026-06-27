#################################################
# Script: simulate_data_SCAL.R
# Author: Nikolas Krstic
# Purpose: Simulate a dataset with only one hierarchical categorical predictor (create true model coefficients and response as well)
#################################################


source("Scalability_Simulations/HCP_creation_SCAL.R")

set.seed(Seed)

#Active Coefficients
if(Sim_Scal_Type=="HCPC"){
  
  ## Layer 2 Nodes
  N2 = 4
  N3 = 8
  N4 = 12
  N5 = 16
  N6 = 20
  
  L3_Node_Means = c()
  
  # Setup the Layer 3 node mean values
  for(i in 1:length(Assignment_List[[2]])){
    
    Curr_Inter_Node = Assignment_List[[2]][i]
    
    if(Curr_Inter_Node==2){
      L2_Node = N2
    }else if(Curr_Inter_Node==3){
      L2_Node = N3
    }else if(Curr_Inter_Node==4){
      L2_Node = N4
    }else if(Curr_Inter_Node==5){
      L2_Node = N5
    }else{
      L2_Node = N6
    }
    
    L3_Node_Means = c(L3_Node_Means, rnorm(1, mean=L2_Node, sd = 1))
    
  }
  
  #Active Leaf Coeffs
  True_Beta = c()
  
  if(!length(Active_Leaves)==0){
    for(i in 1:length(Active_Leaves)){
      
      Curr_Active_Leaf = Active_Leaves[i]
      
      Curr_Path = Path_List[[Curr_Active_Leaf]]
      
      Parent_Node = Curr_Path[3]
      
      True_Beta = c(True_Beta, rnorm(1, mean=L3_Node_Means[Parent_Node-6], sd=0.2))
      
    }
  }
  
  if(!length(Active_Nodes_Set)==0){
    for(i in 1:length(Active_Nodes_Set)){
      
      Curr_Active_Node = Active_Nodes_Set[i]
      
      if(Curr_Active_Node %in% 7:36){
        True_Beta = c(True_Beta, L3_Node_Means[Curr_Active_Node-6])
      }else if(Curr_Active_Node==6){
        True_Beta = c(True_Beta, N6)
      }else if(Curr_Active_Node==5){
        True_Beta = c(True_Beta, N5)
      }else if(Curr_Active_Node==4){
        True_Beta = c(True_Beta, N4)
      }else if(Curr_Active_Node==3){
        True_Beta = c(True_Beta, N3)
      }else if(Curr_Active_Node==2){
        True_Beta = c(True_Beta, N2)
      }else if(Curr_Active_Node==1){
        True_Beta = c(True_Beta, runif(1, min=1, max=5)*sample(c(-1,1), 1))
      }
      
    }
  }
  
  
}else if(Sim_Scal_Type=="SS"){
  
  ## Layer 2 Nodes
  N2 = 4
  N3 = 8
  
  L3_Node_Means = c()
  
  # Setup the Layer 3 node mean values
  for(i in 1:length(Assignment_List[[2]])){
    
    Curr_Inter_Node = Assignment_List[[2]][i]
    
    if(Curr_Inter_Node==2){
      L2_Node = N2
    }else{
      L2_Node = N3
    }
    
    L3_Node_Means = c(L3_Node_Means, rnorm(1, mean=L2_Node, sd = 1))
    
  }
  
  #Active Leaf Coeffs
  True_Beta = c()
  
  if(!length(Active_Leaves)==0){
    for(i in 1:length(Active_Leaves)){
      
      Curr_Active_Leaf = Active_Leaves[i]
      
      Curr_Path = Path_List[[Curr_Active_Leaf]]
      
      Parent_Node = Curr_Path[3]
      
      True_Beta = c(True_Beta, rnorm(1, mean=L3_Node_Means[Parent_Node-3], sd=0.2))
      
    }
  }
  
  if(!length(Active_Nodes_Set)==0){
    for(i in 1:length(Active_Nodes_Set)){
      
      Curr_Active_Node = Active_Nodes_Set[i]
      
      if(Curr_Active_Node %in% 4:8){
        True_Beta = c(True_Beta, L3_Node_Means[Curr_Active_Node-3])
      }else if(Curr_Active_Node==3){
        True_Beta = c(True_Beta, N3)
      }else if(Curr_Active_Node==2){
        True_Beta = c(True_Beta, N2)
      }else if(Curr_Active_Node==1){
        True_Beta = c(True_Beta, runif(1, min=1, max=5)*sample(c(-1,1), 1))
      }
      
    }
  }
  
  
  
}

#True_Beta = rnorm(ncol(Final_True_Design_Matrix), mean=5, sd=2)

# Specify Error/Noise Factor using SNR and signal, then obtain errors
Error_Factor = as.numeric(sqrt(var(Final_True_Design_Matrix %*% True_Beta))/sqrt(SNR))
Errors = rnorm(Samp_Size)*Error_Factor

#Generate the Response
Response = Final_True_Design_Matrix %*% True_Beta + Errors

#Generate the Response for the validation and test sets
Validation_Response = True_Valid_Design_Matrix %*% True_Beta + rnorm(Valid_Samp_Size)*Error_Factor
Test_Response = True_Test_Design_Matrix %*% True_Beta + rnorm(Test_Samp_Size)*Error_Factor


##Identify the true beta values we'll expect to see from the analysis (i.e. across all variables of the HCP design matrix)

Beta_Comp = rep(0, Categ_Counts[length(Categ_Counts)])

#First few coefficients of HCP are for the leaf nodes
Beta_Comp[Active_Leaves] = True_Beta[1:length(Active_Leaves)]

#Remaining coefficients of HCP are based on the merged binary indicators.
if(length(Active_Leaves) != Categ_Counts[length(Categ_Counts)]){
  
  for(i in 1:ncol(Merged_Binary_Indicators)){
    
    Beta_Comp[True_Inactive_Leaves_Sets[[i]]] = True_Beta[length(Active_Leaves)+i]
    
  }
}

