#################################################
# Script: simulation_results_SPARSE.R
# Author: Nikolas Krstic
# Purpose: Create the results of the simulation to understand how HFG LASSO is performing
#################################################

library(ClustAssess)
library(ggplot2)
library(Matrix)
library(tidyr)
library(gridExtra)

#Set the SNR and Scenario
Seeds = seq(100, 1090, by=10)
SNR = 4
Scenario = "C"
Weight_Scheme = "Layer_Based"

### Location to Find Results
Model_Results_Loc = paste("Sparse_Data_Simulation/Results/Scenario_", Scenario, "/", sep="")

#Record performance measures in matrix
BISM_Matrix = matrix(NA, nrow=length(Seeds), ncol=16)
colnames(BISM_Matrix) = c("ECCS_Min", "ECCS_1SE", "ECCS_TBA_Min", "ECCS_TBA_1SE", "EE_Min", "EE_1SE", "EE_TBA_Min", "EE_TBA_1SE",
                          "Opt_Min_Lambda_RMSE", "Opt_1SE_Lambda_RMSE", "Opt_Min_TBA_RMSE", "Opt_1SE_TBA_RMSE",
                          "Fuse_Check_HFGL_Min", "Fuse_Check_HFGL_1SE", "Fuse_Check_TBA_Min", "Fuse_Check_TBA_1SE")

EE_Mins = c()
EE_OLSs = c()


# Save plots generated in the for loop to this location
pdf(paste(Model_Results_Loc, "Scenario_", Scenario, "_SNR_", SNR, "_", Weight_Scheme, "_RegPaths.pdf", sep=""), width=12, height=8)

for(i in 1:length(Seeds)){
  
  Seed = Seeds[i]
  print(Seed)
  
  #Load the Simulation Results
  load(paste(Model_Results_Loc, "HFGL_CV_Seed_", Seed, "_SNR_", SNR, "_", Weight_Scheme, ".Rdata", sep=""))
  
  ##### HFG LASSO Regularization Path Data
  
  #Obtain the HFG LASSO coefficients
  Coefficient_Set = t(HFGL_CV_Model_Results$Coefficients[,c(-1,-2)])
  Lambda_Vals = as.numeric(sapply(strsplit(rownames(Coefficient_Set), "_"), function(x){x[[2]]}))
  Coefficient_Set = cbind(Lambda_Vals, Coefficient_Set)
  
  #HFG LASSO Regularization Path
  Coefficient_Long_Set = gather(as.data.frame(Coefficient_Set), key="Coefficient", val="Value", -Lambda_Vals)
  colnames(Coefficient_Long_Set) = c("Lambda", "Coefficient", "Value")
  Coefficient_Long_Set$Coefficient = factor(sapply(strsplit(Coefficient_Long_Set$Coefficient, "X"), function(x){x[[2]]}), levels=as.character(1:20), labels=1:20)
  
  ####### Adjust lambda grid to be standardized to sample size
  Coefficient_Long_Set$Lambda = Coefficient_Long_Set$Lambda/nrow(HFGL_CV_Model_Results$Design_Matrix)
  #######
  
  #Obtain HFG LASSO Selected Lambdas
  Opt_Lambda = HFGL_CV_Model_Results$Selected_Lambda[1]
  HFGL_1SE_Lambda = HFGL_CV_Model_Results$Selected_Lambda[2]
  
  #Identify True Coefficients and whether certain coefficients should merge/fuse
  True_Coefficients = HFGL_CV_Model_Results$Coefficients[,1]
  Non_Merge_Check = which(9:28 %in% HFGL_CV_Model_Results$Active_Nodes)
  if(length(Non_Merge_Check)==0){
    Coefficient_Long_Set$`Merged?` = TRUE
  }else{
    Coefficient_Long_Set$`Merged?` = !Coefficient_Long_Set$Coefficient %in% Non_Merge_Check
  }
  
  ####### TBA Regularization Path Data
  
  #Obtain the TBA coefficients
  TBA_Coefficient_Set = t(HFGL_CV_Model_Results$TBA_Final_Model$beta[[1]])
  TBA_Lambda_Vals = HFGL_CV_Model_Results$TBA_Final_Model$lambda
  TBA_Coefficient_Set = cbind(TBA_Lambda_Vals, TBA_Coefficient_Set)
  colnames(TBA_Coefficient_Set) = c(colnames(TBA_Coefficient_Set)[1], paste("V", 1:20, sep=""))
  
  #TBA Regularization Path
  TBA_Coefficient_Long_Set = gather(as.data.frame(TBA_Coefficient_Set), key="Coefficient", val="Value", -TBA_Lambda_Vals)
  colnames(TBA_Coefficient_Long_Set) = c("Lambda", "Coefficient", "Value")
  TBA_Coefficient_Long_Set$Coefficient = factor(sapply(strsplit(TBA_Coefficient_Long_Set$Coefficient, "V"), function(x){x[[2]]}), levels=as.character(1:20), labels=1:20)
  
  #Obtain TBA Selected Lambdas
  TBA_Opt_Lambda = HFGL_CV_Model_Results$TBA_Opt_Lambdas[1]
  TBA_1SE_Lambda = HFGL_CV_Model_Results$TBA_Opt_Lambdas[2]
  
  #Identify whether certain coefficients should merge/fuse
  if(length(Non_Merge_Check)==0){
    TBA_Coefficient_Long_Set$`Merged?` = TRUE
  }else{
    TBA_Coefficient_Long_Set$`Merged?` = !TBA_Coefficient_Long_Set$Coefficient %in% Non_Merge_Check
  }
  
  ###################
  ## Compute ECCS/EE/RMSE
  
  #Obtain the coefficients for each method
  HFG_LASSO_Min_Coeffs = HFGL_CV_Model_Results$Coefficients[,paste("Lambda_", HFGL_CV_Model_Results$Selected_Lambda[1], sep="")]
  TBA_Min_Coeffs = HFGL_CV_Model_Results$TBA_Final_Model$beta[[1]][, which(HFGL_CV_Model_Results$TBA_Final_Model$lambda == HFGL_CV_Model_Results$TBA_Opt_Lambdas[1])]
  HFG_LASSO_1SE_Coeffs = HFGL_CV_Model_Results$Coefficients[,paste("Lambda_", HFGL_CV_Model_Results$Selected_Lambda[2], sep="")]
  TBA_1SE_Coeffs = HFGL_CV_Model_Results$TBA_Final_Model$beta[[1]][, which(HFGL_CV_Model_Results$TBA_Final_Model$lambda == HFGL_CV_Model_Results$TBA_Opt_Lambdas[2])]
  
  #Round the coefficients to the fifth decimal place (to assess group fusion and account for imprecision due to tolerance)
  True_Coeff_Aggregs = as.numeric(as.factor(round(True_Coefficients, 5)))
  TBA_Min_Coeff_Aggregs = as.numeric(as.factor(round(TBA_Min_Coeffs, 5)))
  HFG_LASSO_Min_Coeff_Aggregs = as.numeric(as.factor(round(HFG_LASSO_Min_Coeffs, 5)))
  HFG_LASSO_1SE_Coeff_Aggregs = as.numeric(as.factor(round(HFG_LASSO_1SE_Coeffs, 5)))
  TBA_1SE_Coeff_Aggregs = as.numeric(as.factor(round(TBA_1SE_Coeffs, 5)))
  
  #Compute the ECCS values for each method
  ECCS_Min = round(element_sim(True_Coeff_Aggregs, HFG_LASSO_Min_Coeff_Aggregs), 3)
  ECCS_TBA_Min = round(element_sim(True_Coeff_Aggregs, TBA_Min_Coeff_Aggregs), 3)
  ECCS_1SE = round(element_sim(True_Coeff_Aggregs, HFG_LASSO_1SE_Coeff_Aggregs), 3)
  ECCS_TBA_1SE = round(element_sim(True_Coeff_Aggregs, TBA_1SE_Coeff_Aggregs), 3)
  
  #Compute the EE values for each method
  EE_Min = norm(True_Coefficients - HFG_LASSO_Min_Coeffs, type="2")^2/length(True_Coefficients)
  EE_TBA_Min = norm(True_Coefficients - TBA_Min_Coeffs, type="2")^2/length(True_Coefficients)
  EE_1SE = norm(True_Coefficients - HFG_LASSO_1SE_Coeffs, type="2")^2/length(True_Coefficients)
  EE_TBA_1SE = norm(True_Coefficients - TBA_1SE_Coeffs, type="2")^2/length(True_Coefficients)
  EE_OLS = norm(True_Coefficients - HFGL_CV_Model_Results$Coefficients[,2], type="2")^2/length(True_Coefficients)
  EE_Mins = c(EE_Mins, EE_Min)
  EE_OLSs = c(EE_OLSs, EE_OLS)
  
  #Obtain the RMSE values for each method
  Opt_Min_Lambda_RMSE = HFGL_CV_Model_Results$Final_Test_Errors[2]
  Opt_Min_TBA_RMSE = HFGL_CV_Model_Results$Final_Test_Errors[5]
  Opt_1SE_Lambda_RMSE = HFGL_CV_Model_Results$Final_Test_Errors[3]
  Opt_1SE_TBA_RMSE = HFGL_CV_Model_Results$Final_Test_Errors[6]
  
  ###################
  # Custom Comparison Measure - Complete Group Fusion Index (CGFI)
  
  #Storage for descendant coefficient indices for each ancestor node requiring fusion
  Desc_Coeff_Inds = list()
  
  #Storage for the binary checks of whether a group fusion occurred at each ancestor node
  HFG_Checks = c()
  HFG_1SE_Checks = c()
  TBA_Checks = c()
  TBA_1SE_Checks = c()
  True_Coeff_Checks = c()
  
  #Storage for Path List
  Path_List = list()
  
  #Tree Structure
  Tree_Struct = HFGL_CV_Model_Results$Tree_Structure
  
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
    HFG_1SE_Checks[k] = length(unique(round(HFG_LASSO_1SE_Coeffs[Desc_Coeff_Inds[[k]]], 5)))==1
    
    TBA_Checks[k] = length(unique(round(TBA_Min_Coeffs[Desc_Coeff_Inds[[k]]], 5)))==1
    TBA_1SE_Checks[k] = length(unique(round(TBA_1SE_Coeffs[Desc_Coeff_Inds[[k]]], 5)))==1
    
    True_Coeff_Checks[k] = length(unique(round(True_Coefficients[Desc_Coeff_Inds[[k]]], 5)))==1
    
  }
  
  #CGFI Measures
  HFG_Fuse_Check = sum(HFG_Checks==True_Coeff_Checks)/length(HFG_Checks)
  HFG_1SE_Fuse_Check = sum(HFG_1SE_Checks==True_Coeff_Checks)/length(HFG_1SE_Checks)

  TBA_Fuse_Check = sum(TBA_Checks==True_Coeff_Checks)/length(TBA_Checks)
  TBA_1SE_Fuse_Check = sum(TBA_1SE_Checks==True_Coeff_Checks)/length(TBA_1SE_Checks)
  
  #Compile all of the measures together and store them
  BISM_Matrix[i,] = c(ECCS_Min, ECCS_1SE, ECCS_TBA_Min, ECCS_TBA_1SE, EE_Min, EE_1SE, EE_TBA_Min, EE_TBA_1SE,
                      Opt_Min_Lambda_RMSE, Opt_1SE_Lambda_RMSE, Opt_Min_TBA_RMSE, Opt_1SE_TBA_RMSE,
                      HFG_Fuse_Check, HFG_1SE_Fuse_Check, TBA_Fuse_Check, TBA_1SE_Fuse_Check)
  
  
  #OPTIONAL!!!!! (Used for Figure 3) Reduce the Lambda grid range so that the plot has a more reasonable range of lambda values
  #Coefficient_Long_Set = Coefficient_Long_Set[!(Coefficient_Long_Set$Lambda %in% unique(sort(Coefficient_Long_Set$Lambda))[51:60]),]
  
  #Create regularization path plot for HFG LASSO
  Reg_Path_Plot = ggplot(Coefficient_Long_Set, aes(group=Coefficient, x=Lambda, y=Value))+
    geom_line(aes(colour=`Merged?`))+
    theme_bw()+
    labs(color="Merged?", title="HFG LASSO Regularization Path")+
    theme(axis.text.x=element_text(size=20, face="bold"),
          plot.title=element_text(hjust = 0.5, size=20, face="bold"),
          axis.text.y=element_text(size=20, face="bold"),
          axis.title.x=element_text(size=24),
          axis.title.y=element_text(size=24),
          legend.title = element_text(size=16),
          legend.text = element_text(size=12),
          strip.text.x=element_text(size=24),
          strip.text.y=element_text(size=24))+
    scale_color_manual(breaks=c(FALSE, TRUE), values = c("black", "blue"))+
    geom_vline(xintercept=Opt_Lambda/nrow(HFGL_CV_Model_Results$Design_Matrix), linetype="dotted")+
    geom_vline(xintercept=HFGL_1SE_Lambda/nrow(HFGL_CV_Model_Results$Design_Matrix), linetype="dotted", colour="orange")+
    geom_hline(yintercept=unique(True_Coefficients), linetype="dotted", colour="red", size=0.6)
  
  #OPTIONAL!!!!! (Used for Figure 3) Reduce the Lambda grid range so that the plot has a more reasonable range of lambda values
  #TBA_Coefficient_Long_Set = TBA_Coefficient_Long_Set[!TBA_Coefficient_Long_Set$Lambda %in% TBA_Lambda_Vals[1:10],]
  
  #Create regularization path plot for TBA
  Reg_Path_Plot2 = ggplot(TBA_Coefficient_Long_Set, aes(group=Coefficient, x=Lambda, y=Value))+
    geom_line(aes(colour=`Merged?`))+
    theme_bw()+
    labs(color="Merged?", title="TBA Regularization Path")+
    theme(axis.text.x=element_text(size=20, face="bold"),
          plot.title=element_text(hjust = 0.5, size=20, face="bold"),
          axis.text.y=element_text(size=20, face="bold"),
          axis.title.x=element_text(size=24),
          axis.title.y=element_text(size=24),
          legend.title = element_text(size=16),
          legend.text = element_text(size=12),
          strip.text.x=element_text(size=24),
          strip.text.y=element_text(size=24))+
    scale_color_manual(breaks=c(FALSE, TRUE), values = c("black", "blue"))+
    geom_vline(xintercept=TBA_Opt_Lambda, linetype="dotted")+
    geom_vline(xintercept=TBA_1SE_Lambda, linetype="dotted", colour="orange")+
    geom_hline(yintercept=unique(True_Coefficients), linetype="dotted", colour="red", size=0.6)
  
  grid.arrange(Reg_Path_Plot, Reg_Path_Plot2, ncol=2)
  ##################
  
}


dev.off()


#Save the matrix to create Tables 1-3
save(BISM_Matrix, file=paste(Model_Results_Loc, "BISM_Matrix_SNR", SNR, "_", Weight_Scheme, ".Rdata", sep=""))

#Mean Measures Across Replications
colMeans(BISM_Matrix)

#Standard Error of Mean Measures Across Replications
apply(BISM_Matrix, 2, sd)/sqrt(length(Seeds))



