#################################################
# Script: HFGL_US_Colleges_Application_CV.R
# Author: Nikolas Krstic
# Purpose: Apply HFG LASSO to the US Colleges Dataset (see Cerda and Varoquaux (2020))
#################################################

Rcpp::sourceCpp("HFGL_Fitting_Algorithm_1_Step_3_Functions.cpp")
source("US_Colleges_Application/HFGL_US_Colleges_Application_DataCleaning.R")
source("HFG_LASSO_func.R")
source("US_Colleges_Application/US_Colleges_func.R")

library(rare)
library(Matrix)
library(CVXR)

#############################################################################################################
# CV Version of the Modelling

Seeds = seq(10, 1000, by=10)

for(i in 1:length(Seeds)){

  Seed = Seeds[i]
  
  print(paste("Seed: ", Seed, sep=""))
  
  set.seed(Seed)
  
  #Form Validation Set
  Validation_Set_Indices = sample(1:nrow(College_Data), ceiling(nrow(College_Data)*0.1))
  Validation_Set = College_Data[Validation_Set_Indices,]
  
  #Form Training Set
  College_Training_Set = College_Data[which(!(1:nrow(College_Data) %in% Validation_Set_Indices)),]
  
  #Create new design matrix from training set
  DesignMat_Geog_HCP_Train = model.matrix(~0+State, data=College_Training_Set)
  colnames(DesignMat_Geog_HCP_Train) = gsub("State", "", colnames(DesignMat_Geog_HCP_Train))
  
  #Weights for each of the penalty terms (using layer-based weighting scheme)
  S_Matrices_Weights = sapply(S_Matrix_Layers, function(x){if(x==3){1}else if(x==2){1/8}else if(x==1){1/24}})
  
  ## Arguments for the HFG LASSO function
  HierDesignMatrix = DesignMat_Geog_HCP_Train
  Response = College_Training_Set$`Percent Pell Grant`
  
  Lambda_Min = HFG_LASSO_Lambda_Min(HierDesignMatrix, Response, S_Matrices, S_Matrices_Weights)
  Lambda_Set = exp(seq(0, log(0.3*Lambda_Min), length.out=60))
  TBA_Lambda_Count = 60
  
  #Fit tolerance to use when assessing convergence of the fit algorithm (Overall Coefficient convergence)
  Proximal_Grad_Tol = 1e-8
  #SDMM Setting
  Gamma = 1
  
  
  #############################################################################
  # Set up settings for Tree-Based Aggregation Method (to compare to HFG LASSO)
  
  ## Lambda Set to use (based on how TBA generates the set in rarefit function)
  TBA_Lambdas = max(abs(t(HierDesignMatrix) %*% Response))/nrow(HierDesignMatrix) * exp(seq(0, log(1e-4), len = TBA_Lambda_Count))
  
  ## Tree ID Dataframe (Assign numeric node value to nodes in tree)
  Tree_ID_DF = data.frame(ID = unique(unlist(Path_List)), Val=as.numeric(unique(unlist(Path_List))))
  
  ## Generate "A" matrix based on tree
  A = matrix(0, nrow=ncol(HierDesignMatrix), ncol=nrow(Tree_ID_DF))
  
  ##Assign 1 values to appropriate locations in the "A" matrix (based on tree structure)
  for(j in 1:nrow(A)){
    A[j,Tree_ID_DF[which(Tree_ID_DF[,"ID"] %in% as.vector(unlist(Path_List[j,]))),"Val"]] = 1
  }
  
  #Has to be a sparse matrix
  A = Matrix(A, sparse=TRUE)
  
  Temp_Train_Group_DF = data.frame(Response = Response, State = as.factor(colnames(HierDesignMatrix)[HierDesignMatrix %*% 1:ncol(HierDesignMatrix)]))
  
  ## Cross-Validation Approach
  #Fold Number
  Num_Folds = 5
  #Min Lambda or Min+1SE Lambda
  Lambda_1SE_Bool = TRUE
  
  
  # Make cross-validation folds:
  
  #Shuffle observation indices
  Obs_Indices = sample(1:nrow(HierDesignMatrix), replace=FALSE)
  
  #Set up folds for cross-validation
  Folds = split(Obs_Indices, cut(seq_along(Obs_Indices), Num_Folds, labels=FALSE))
  
  ####################################################
  # Compile all results
  
  #Compute Final Models for Each Method
  ### HFG LASSO
  HFG_LASSO_CV_Results = PredPerf_CV_Computation(Method="HFG LASSO", Lambda_Set, Folds, HierDesignMatrix, Response)
  HFG_LASSO_Thetas = HFG_LASSO_CV_Results[[1]]
  Valid_Errors = HFG_LASSO_CV_Results[[2]]
  Opt_Lambdas = HFG_LASSO_CV_Results[[3]]
  
  ### TBA
  TBA_CV_Results = PredPerf_CV_Computation(Method="TBA", TBA_Lambdas, Folds, HierDesignMatrix, Response)
  TBA_Model_Thetas = TBA_CV_Results[[1]]
  TBA_Valid_Errors = TBA_CV_Results[[2]]
  TBA_Opt_Lambdas = TBA_CV_Results[[3]]
  
  ###OLS Results for Coefficients
  DF = as.data.frame(cbind(Response, HierDesignMatrix))
  
  LM = lm(Response~.-1, data=DF)
  LM_Coeffs = coefficients(LM)
  
  
  ##### Compute Test Errors for Each Method
  
  DesignMat_Geog_HCP_Test = model.matrix(~0+State, data=Validation_Set)
  colnames(DesignMat_Geog_HCP_Test) = gsub("State", "", colnames(DesignMat_Geog_HCP_Test))
  
  Final_Test_Response = Validation_Set$`Percent Pell Grant`
  
  
  HFG_Final_Test_Error_Min = sqrt(sum((DesignMat_Geog_HCP_Test %*% HFG_LASSO_Thetas[,1] - Final_Test_Response)^2)/length(Final_Test_Response))
  HFG_Final_Test_Error_1SE = sqrt(sum((DesignMat_Geog_HCP_Test %*% HFG_LASSO_Thetas[,2] - Final_Test_Response)^2)/length(Final_Test_Response))
  
  TBA_Final_Test_Error_Min = sqrt(sum((DesignMat_Geog_HCP_Test %*% TBA_Model_Thetas[,1] - Final_Test_Response)^2)/length(Final_Test_Response))
  TBA_Final_Test_Error_1SE = sqrt(sum((DesignMat_Geog_HCP_Test %*% TBA_Model_Thetas[,2] - Final_Test_Response)^2)/length(Final_Test_Response))
  
  OLS_Final_Test_Error = sqrt(sum((DesignMat_Geog_HCP_Test %*% LM_Coeffs - Final_Test_Response)^2)/length(Final_Test_Response))
  
  
  Final_Test_Errors = c(OLS_Final_Test_Error, HFG_Final_Test_Error_Min, HFG_Final_Test_Error_1SE,
                        TBA_Final_Test_Error_Min, TBA_Final_Test_Error_1SE)
  
  names(Final_Test_Errors) = c("OLS Test Error", "Min CV HFGL Test Error", "1SE CV HFGL Test Error",
                               "Min CV TBA Test Error", "1SE CV TBA Test Error")
  
  #Compile Results
  #Compile Coefficient Results
  Thetas_Set = cbind(LM_Coeffs=LM_Coeffs, data.frame(HFG_LASSO_Thetas), data.frame(TBA_Model_Thetas))
  names(Thetas_Set) = c("LM_Coeffs", "HFGL_Min_Coeffs", "HFGL_1SE_Coeffs", "TBA_Min_Coeffs", "TBA_1SE_Coeffs")
  
  #Calculate Number of Unique Coefficients to the sixth decimal place in each case
  Unique_Coeff_Nums = apply(Thetas_Set, 2, function(x){length(unique(round(x, 6)))})
  
  US_College_Data_HFGL_CV_Model_Results = list(Coefficients = Thetas_Set, 
                                               Design_Matrix = HierDesignMatrix,
                                               Response = Response, 
                                               Validation_Set = Validation_Set,
                                               HFGL_Selected_Lambdas = Opt_Lambdas,
                                               TBA_Selected_Lambdas = TBA_Opt_Lambdas,
                                               HFGL_Validation_Errors = Valid_Errors,
                                               TBA_Validation_Errors = TBA_Valid_Errors,
                                               Final_Test_Errors = Final_Test_Errors,
                                               Unique_Coeff_Nums = Unique_Coeff_Nums)
  
  save(US_College_Data_HFGL_CV_Model_Results,
       file=paste("./US_Colleges_Application/Results/US_College_Data_CV_HFGL_Results_Holdout_Seed", Seed , "_Layer_Based.Rdata", sep=""))
  
}


###############
### Example Repetition Results:

Seed = 10
load(file=paste("./US_Colleges_Application/Results/US_College_Data_CV_HFGL_Results_Holdout_Seed", Seed , "_Layer_Based.Rdata", sep=""))

US_College_Data_HFGL_CV_Model_Results$Coefficients
US_College_Data_HFGL_CV_Model_Results$Coefficients$HFGL_1SE_Coeffs

##Path_List

#Far West, Rocky Mountains, New England, Mid East, and Plains All Group Fuse for HFG LASSO 1 SE Rule
Group_Result = merge(Path_List, 
      cbind(State=rownames(US_College_Data_HFGL_CV_Model_Results$Coefficients),
            HFG_LASSO_Group=as.numeric(as.factor(round(US_College_Data_HFGL_CV_Model_Results$Coefficients$HFGL_1SE_Coeffs,5))),
            HFG_LASSO_Coeff=as.numeric(round(US_College_Data_HFGL_CV_Model_Results$Coefficients$HFGL_1SE_Coeffs,5)),
            TBA_Group=as.numeric(as.factor(round(US_College_Data_HFGL_CV_Model_Results$Coefficients$TBA_1SE_Coeffs,5))),
            TBA_Coeff=as.numeric(round(US_College_Data_HFGL_CV_Model_Results$Coefficients$TBA_1SE_Coeffs,5))),
      by="State")

Group_Result[order(Group_Result$Region),]


#######################################################
### Creation of Table 4

# Load the method results for each seed, and extract the final test RMSE and the number of unique coefficients for each method
Seeds = seq(10, 1000, by=10)

for(i in 1:length(Seeds)){
  
  Seed = Seeds[i]
  
  load(paste("./US_Colleges_Application/Results/US_College_Data_CV_HFGL_Results_Holdout_Seed", Seed, "_Layer_Based.Rdata", sep=""))
  
  Unique_Coeff_Nums_Check = apply(US_College_Data_HFGL_CV_Model_Results$Coefficients[1:51,], 2, function(x){length(unique(round(x, 6)))})
  
  if(i==1){
    Final_Test_Errors_Set = US_College_Data_HFGL_CV_Model_Results$Final_Test_Errors
    
    Unique_Coeff_Nums = Unique_Coeff_Nums_Check
  }else{
    Final_Test_Errors_Set = rbind(Final_Test_Errors_Set, US_College_Data_HFGL_CV_Model_Results$Final_Test_Errors)
    
    Unique_Coeff_Nums = rbind(Unique_Coeff_Nums, Unique_Coeff_Nums_Check)
  }
  
  
}

#Compute the Final Test Error Means (Average Across Repetitions)
Final_Test_Error_Means = colMeans(Final_Test_Errors_Set)

#Compute Average Percentage Model Sparsity
Unique_Coeff_Nums = paste(round((1-(colMeans(Unique_Coeff_Nums)/51))*100, 3), "\\%", sep="")

#Compute the Final Test Error Standard Errors (Across Repetitions)
Final_Test_Error_SEs = apply(Final_Test_Errors_Set, 2, sd)/sqrt(nrow(Final_Test_Errors_Set))

#Create Table 4
TestErrors_Tab = data.frame(`RMSE Mean Estimate` = as.vector(Final_Test_Error_Means), `RMSE Standard Error` = as.vector(Final_Test_Error_SEs), `Coefficient Sparsity` = as.vector(Unique_Coeff_Nums))

rownames(TestErrors_Tab) = c("OLS", paste("HFG LASSO - MCVE", sep=""), paste("HFG LASSO - 1 SE Rule", sep=""), paste("TBA - MCVE", sep=""), paste("TBA - 1 SE Rule", sep=""))

colnames(TestErrors_Tab) = c("RMSE Mean", "RMSE SE", "HCP Sparsity")









