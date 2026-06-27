#################################################
# Script: simulation_results_SCAL.R
# Author: Nikolas Krstic
# Purpose: Create the results of the simulation
#################################################

library(ClustAssess)
library(ggplot2)
library(Matrix)
library(tidyr)
library(gridExtra)

#Set the SNR and the Scale Type
Seeds = seq(10, 200, by=10)
SNR = 1
Weight_Scheme = "Layer_Based"
Sim_Scal_Type = "HCPC"

if(Sim_Scal_Type == "SS"){
  Series = seq(1000, 10000, by=1000)
}else if(Sim_Scal_Type == "HCPC"){
  Series = seq(50, 500, 50)
}

### Location to Find Results
Model_Results_Loc = paste("Scalability_Simulations/Results/", Sim_Scal_Type, "_Results/", sep="")

#Record performance measures in matrix
BISM_Matrix = matrix(NA, nrow=length(Seeds), ncol=10)
colnames(BISM_Matrix) = c("ECCS_Min", "ECCS_TBA_Min", "EE_Min", "EE_TBA_Min",
                          "Opt_Min_Lambda_RMSE", "Opt_Min_TBA_RMSE",
                          "Fuse_Check_HFGL_Min", "Fuse_Check_TBA_Min",
                          "HFG_L_Run_Time", "TBA_Run_Time")

for(j in 1:length(Series)){
  
  for(i in 1:length(Seeds)){
    
    Seed = Seeds[i]
    print(Seed)
    
    #Load the Simulation Results
    if(Sim_Scal_Type=="SS"){
      Samp_Size = Series[j]
      load(paste(Model_Results_Loc, "HFGL_Seed_", Seed, "_Samp_Size_", Samp_Size, "_SNR_", SNR, "_", Weight_Scheme, ".Rdata", sep=""))
    }else if(Sim_Scal_Type=="HCPC"){
      Cardinality = Series[j]
      load(paste(Model_Results_Loc, "HFGL_Seed_", Seed, "_Cardinality_", Cardinality, "_SNR_", SNR, "_", Weight_Scheme, ".Rdata", sep=""))
    }
    
    ##### HFG LASSO Regularization Path Data
    
    #Obtain the HFG LASSO coefficients
    Coefficient_Set = t(HFGL_Model_Results$Coefficients[,c(-1,-2)])
    Lambda_Vals = as.numeric(sapply(strsplit(rownames(Coefficient_Set), "_"), function(x){x[[2]]}))
    Coefficient_Set = cbind(Lambda_Vals, Coefficient_Set)
    
    #Obtain HFG LASSO Selected Lambdas
    Opt_Lambda = HFGL_Model_Results$Selected_Lambda[1]
    #HFGL_1SE_Lambda = HFGL_Model_Results$Selected_Lambda[2]
    
    #Identify True Coefficients 
    True_Coefficients = HFGL_Model_Results$Coefficients[,1]
    
    ####### TBA Regularization Path Data
    
    #Obtain the TBA coefficients
    TBA_Coefficient_Set = t(HFGL_Model_Results$TBA_Final_Model_Coefficients)
    TBA_Lambda_Vals = HFGL_Model_Results$TBA_Final_Model_Lambda
    TBA_Coefficient_Set = cbind(TBA_Lambda_Vals, TBA_Coefficient_Set)
    colnames(TBA_Coefficient_Set) = c(colnames(TBA_Coefficient_Set)[1], paste("V", 1:length(True_Coefficients), sep=""))
    
    #TBA Regularization Path
    TBA_Coefficient_Long_Set = gather(as.data.frame(TBA_Coefficient_Set), key="Coefficient", val="Value", -TBA_Lambda_Vals)
    colnames(TBA_Coefficient_Long_Set) = c("Lambda", "Coefficient", "Value")
    TBA_Coefficient_Long_Set$Coefficient = factor(sapply(strsplit(TBA_Coefficient_Long_Set$Coefficient, "V"), function(x){x[[2]]}), levels=as.character(1:length(True_Coefficients)), labels=1:length(True_Coefficients))
    
    #Obtain TBA Selected Lambdas
    TBA_Opt_Lambda = HFGL_Model_Results$TBA_Opt_Lambda
    
    ###################
    ## Compute ECCS/EE/RMSE
    
    #Obtain the coefficients for each method
    HFG_LASSO_Min_Coeffs = HFGL_Model_Results$Coefficients[,paste("Lambda_", HFGL_Model_Results$Selected_Lambda, sep="")]
    TBA_Min_Coeffs = HFGL_Model_Results$TBA_Final_Model_Coefficients[, which(HFGL_Model_Results$TBA_Final_Model_Lambda == HFGL_Model_Results$TBA_Opt_Lambda)]
    
    #Round the coefficients to the fifth decimal place (to assess group fusion and account for imprecision due to tolerance)
    True_Coeff_Aggregs = as.numeric(as.factor(round(True_Coefficients, 5)))
    TBA_Min_Coeff_Aggregs = as.numeric(as.factor(round(TBA_Min_Coeffs, 5)))
    HFG_LASSO_Min_Coeff_Aggregs = as.numeric(as.factor(round(HFG_LASSO_Min_Coeffs, 5)))
    
    #Compute the ECCS values for each method
    ECCS_Min = round(element_sim(True_Coeff_Aggregs, HFG_LASSO_Min_Coeff_Aggregs), 3)
    ECCS_TBA_Min = round(element_sim(True_Coeff_Aggregs, TBA_Min_Coeff_Aggregs), 3)
    
    #Compute the EE values for each method
    EE_Min = norm(True_Coefficients - HFG_LASSO_Min_Coeffs, type="2")^2/length(True_Coefficients)
    EE_TBA_Min = norm(True_Coefficients - TBA_Min_Coeffs, type="2")^2/length(True_Coefficients)
    EE_OLS = norm(True_Coefficients - HFGL_Model_Results$Coefficients[,2], type="2")^2/length(True_Coefficients)
    
    #Obtain the RMSE values for each method
    Opt_Min_Lambda_RMSE = HFGL_Model_Results$Final_Test_Errors[2]
    Opt_Min_TBA_RMSE = HFGL_Model_Results$Final_Test_Errors[4]
    
    ###################
    # Custom Comparison Measure - Complete Group Fusion Index (CGFI)
    
    #Storage for descendant coefficient indices for each ancestor node requiring fusion
    Desc_Coeff_Inds = list()
    
    #Storage for the binary checks of whether a group fusion occurred at each ancestor node
    HFG_Checks = c()
    TBA_Checks = c()
    True_Coeff_Checks = c()
    
    #Storage for Path List
    Path_List = list()
    
    #Tree Structure
    Tree_Struct = HFGL_Model_Results$Tree_Structure
    
    #Number of Ancestor Nodes in Tree (non-leaf nodes)
    Ances_Node_Num = sum(sapply(Tree_Struct, length)[1:2])+1
    
    #Regenerate Path List
    for(k in 1:length(True_Coefficients)){
      
      Curr_Coeff_Ind = k+Ances_Node_Num
      
      Ancestor_Node_1 = Tree_Struct[[3]][which(sapply(strsplit(names(Tree_Struct[[3]]), " "), function(x){as.numeric(x[2])})==Curr_Coeff_Ind)]
      Ancestor_Node_2 = Tree_Struct[[2]][which(sapply(strsplit(names(Tree_Struct[[2]]), " "), function(x){as.numeric(x[2])})==Ancestor_Node_1)]
      Ancestor_Node_3 = 1
      
      Path_List[[k]] = as.numeric(c(Ancestor_Node_3, Ancestor_Node_2, Ancestor_Node_1, Curr_Coeff_Ind))
      
    }
    
    #Ancestor Node Values
    Ances_Node_Vals = 1:Ances_Node_Num
    
    #Conduct the checks to see if the group fusions occurred for each ancestor, for every method (including the truth)
    for(k in 1:length(Ances_Node_Vals)){
      
      Curr_Ances_Node = Ances_Node_Vals[k]
      Desc_Coeff_Inds[[k]] = which(sapply(Path_List, function(x){Curr_Ances_Node %in% x}))
      
      HFG_Checks[k] = length(unique(round(HFG_LASSO_Min_Coeffs[Desc_Coeff_Inds[[k]]], 5)))==1
      
      TBA_Checks[k] = length(unique(round(TBA_Min_Coeffs[Desc_Coeff_Inds[[k]]], 5)))==1
      
      True_Coeff_Checks[k] = length(unique(round(True_Coefficients[Desc_Coeff_Inds[[k]]], 5)))==1
      
    }
    
    #CGFI Measures
    HFG_Fuse_Check = sum(HFG_Checks==True_Coeff_Checks)/length(HFG_Checks)
    TBA_Fuse_Check = sum(TBA_Checks==True_Coeff_Checks)/length(TBA_Checks)
    
    #Run Times
    HFG_LASSO_Run_Time = HFGL_Model_Results$HFG_LASSO_Run_Time
    TBA_Run_Time = HFGL_Model_Results$TBA_Run_Time
    
    #Compile all of the measures together and store them
    BISM_Matrix[i,] = c(ECCS_Min, ECCS_TBA_Min, EE_Min, EE_TBA_Min,
                        Opt_Min_Lambda_RMSE, Opt_Min_TBA_RMSE,
                        HFG_Fuse_Check, TBA_Fuse_Check, HFG_LASSO_Run_Time, TBA_Run_Time)
    
  }
  
  
  
  #Save the matrix to create Tables 1-3
  if(Sim_Scal_Type=="SS"){
    save(BISM_Matrix, file=paste(Model_Results_Loc, "BISM_Matrix_Samp_Size_", Samp_Size, "_SNR_", SNR, "_", Weight_Scheme, ".Rdata", sep=""))
  }else if(Sim_Scal_Type=="HCPC"){
    save(BISM_Matrix, file=paste(Model_Results_Loc, "BISM_Matrix_Cardinality_", Cardinality, "SNR", SNR, "_", Weight_Scheme, ".Rdata", sep=""))
  }
  
  #Mean Measures Across Replications
  colMeans(BISM_Matrix)
  
  #Standard Error of Mean Measures Across Replications
  apply(BISM_Matrix, 2, sd)/sqrt(length(Seeds))
  
  
}
