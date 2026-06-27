#################################################
# Script: simulate_data.R
# Author: Nikolas Krstic
# Purpose: Simulate a dataset with only one hierarchical categorical predictor (create true model coefficients and response as well)
#################################################


source("Simulations/HCP_creation.R")

set.seed(Seed)

#Active Coefficients
if(Scenario=="A"){
  True_Beta = c(3, 3.5, 4, 4.5, 5, 1, 6, 2, 7)
}else if(Scenario=="B"){
  True_Beta = c(1, 7, 8.5, 2.5, 3, 1.5, 9, 3.5, 4, 7.5, 4.5, 9.5, 10, 5, 6)
}else if(Scenario=="C"){
  True_Beta = c(1.5, 2, 2.5, 3, 6, 4, 5)
}


# Specify Error/Noise Factor using SNR and signal, then obtain errors
Error_Factor = as.numeric(sqrt(var(Final_True_Design_Matrix %*% True_Beta))/sqrt(SNR))
Errors = rnorm(Samp_Size)*Error_Factor

#Generate the Response
Response = Final_True_Design_Matrix %*% True_Beta + Errors

#Generate the Response for the validation/test set
Validation_Response = True_Valid_Design_Matrix %*% True_Beta + rnorm(Valid_Samp_Size)*Error_Factor


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

