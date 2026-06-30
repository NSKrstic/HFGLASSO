#################################################
# Script: simulation_study_SPARSE.R
# Author: Nikolas Krstic
# Purpose: Perform "Hierarchical Fused Group LASSO", with cross-validation AND comparison with Tree-Based Aggregation by Yan and Bien (2021)
#################################################

Rcpp::sourceCpp("HFGL_Fitting_Algorithm_1_Step_3_Functions.cpp")
source("HFG_LASSO_func.R")
source("Sparse_Data_Simulation/simulate_data_SPARSE.R")

library(rare)
library(Matrix)
library(CVXR)

#########################################################################
#########################################################################
#########################################################################
#Fit Hierarchical Fused Group LASSO

#############################################################################
# Set up general settings (and those for HFG LASSO)

#Set up the weights for the HFG LASSO penalty (i.e. the weights for each S Matrix)
if(Weight_Scheme=="Layer_Based"){
  S_Matrices_Weights = sapply(S_Matrix_Layers, function(x){if(x==3){1}else if(x==2){1/5}else if(x==1){1/10}})
  S_Matrices_Weights = S_Matrices_Weights/sum(S_Matrices_Weights)
}else if(Weight_Scheme=="LB_SQRT"){
  S_Matrices_Weights = sapply(S_Matrix_Layers, function(x){if(x==3){1}else if(x==2){1/sqrt(5)}else if(x==1){1/sqrt(10)}})
  S_Matrices_Weights = S_Matrices_Weights/sum(S_Matrices_Weights)
}else if(Weight_Scheme=="Col_Dim"){
  S_Matrices_Weights = (1/sapply(S_Matrices, function(x){sum(apply(x, 2, function(x){all(x==0)})==FALSE)}))
  S_Matrices_Weights = S_Matrices_Weights/sum(S_Matrices_Weights)
}else if(Weight_Scheme=="CD_SQRT"){
  S_Matrices_Weights = sqrt(1/sapply(S_Matrices, function(x){sum(apply(x, 2, function(x){all(x==0)})==FALSE)}))
  S_Matrices_Weights = S_Matrices_Weights/sum(S_Matrices_Weights)
}

#Re-assign design matrix to different name (to distinguish from non-HCP design matrix that could be considered here)
HierDesignMatrix = DesignMatrix
#Tolerance to use when assessing convergence of the HFG LASSO fitting algorithm
Proximal_Grad_Tol = 1e-8
#SDMM Hyperparameter Setting
Gamma = 0.5

## Identify Lambda_Max using Theorem 1
Lambda_Max = HFG_LASSO_Lambda_Max(DesignMatrix, Response, S_Matrices, S_Matrices_Weights)

## Set the Lambda Grid for HFG LASSO and the Number of Lambda to Use for TBA
Lambda_Set = exp(seq(0, log(0.3*Lambda_Max), length.out=60))
TBA_Lambda_Count = 60

#Obtain the category frequencies, the category indices for non-zero columns of the S Matrices, and the fusion pairs for each of the S Matrices
Categ_Freqs = colSums(HierDesignMatrix)
Categ_Inds = lapply(S_Matrices, function(x){which(apply(x, 2, function(x){all(x==0)})==FALSE)})
Categ_Pair_Inds = lapply(S_Matrices, function(x){t(apply(x, 1, function(x){which(x!=0)}))})


#############################################################################
# Set up settings for Tree-Based Aggregation Method (to compare to HFG LASSO)

## Lambda Set to use (based on the default method used to generate the lambda grid in the "rarefit" function)
TBA_Lambdas = max(abs(t(HierDesignMatrix) %*% Response))/Samp_Size * exp(seq(0, log(1e-4), len = TBA_Lambda_Count))

## Generate "A" matrix based on tree
A = matrix(0, nrow=ncol(HierDesignMatrix), ncol=sum(Categ_Counts))

for(i in 1:nrow(A)){
  A[i,Path_List[[i]]] = 1
}

A = Matrix(A, sparse=TRUE)


#############################################################################
## Cross-Validation Approach

#Fold Number
Num_Folds = 5

# Make cross-validation folds:

#Shuffle observation indices
Obs_Indices = sample(1:nrow(HierDesignMatrix), replace=FALSE)

#Set up folds for cross-validation
Folds = split(Obs_Indices, cut(seq_along(Obs_Indices), Num_Folds, labels=FALSE))
#Set up fold vector
Fold_Vec = sapply(1:Samp_Size, function(x){as.vector(which(sapply(Folds, function(z){any(x %in% z)})))})

#Validation Error Lists
Test_Set_Validation_Errors_Set = list()
TBA_Val_Errors = list()

for(q in 1:Num_Folds){
  
  print(paste("Current Fold: ", q, sep=""))
  
  #Set up the training and test sets for the current fold
  Curr_Test_Folds_Indices = Folds[[q]]
  Curr_Train_Folds_Indices = as.vector(unlist(Folds[-q]))
  
  #Build Train and Test Predictor Matrices and Responses
  Train_DesignMatrix = HierDesignMatrix[Curr_Train_Folds_Indices,]
  Train_Response = Response[Curr_Train_Folds_Indices,]
  Test_DesignMatrix = HierDesignMatrix[Curr_Test_Folds_Indices,]
  Test_Response = Response[Curr_Test_Folds_Indices,]
  
  # HFG LASSO
  Thetas_Set_List = HFG_LASSO_Fit(Train_DesignMatrix, Train_Response, Categ_Inds, Categ_Pair_Inds, S_Matrices_Weights,
                                  Lambda_Set, Proximal_Grad_Tol, Gamma)
  
  # TBA
  TBA_Model = rarefit(y=Train_Response, X=Train_DesignMatrix, A = A, intercept = FALSE, lambda = TBA_Lambdas, alpha = 1)
  
  #Predictions
  HFG_LASSO_Preds = (Test_DesignMatrix %*% Thetas_Set_List - replicate(length(Lambda_Set), Test_Response))
  TBA_Preds = (Test_DesignMatrix %*% TBA_Model$beta[[1]] - replicate(length(TBA_Lambdas), Test_Response))
  
  #Validation Set Errors
  Test_Set_Validation_Errors_Set[[q]] = apply(HFG_LASSO_Preds, 2, function(x){sqrt(sum(x^2)/nrow(Test_DesignMatrix))})
  TBA_Val_Errors[[q]] = apply(TBA_Preds, 2, function(x){sqrt(sum(x^2)/nrow(Test_DesignMatrix))})
}

# HFG LASSO Validation Summaries
Valid_Errors = do.call(cbind, Test_Set_Validation_Errors_Set)
Mean_Val_Errors = round(rowMeans(Valid_Errors, na.rm=TRUE), 4)
SE_Val_Errors = apply(Valid_Errors, 1, sd)/sqrt(Num_Folds)

# TBA Validation Summaries
TBA_Valid_Errors = do.call(cbind, TBA_Val_Errors)
TBA_Mean_Val_Errors = round(rowMeans(TBA_Valid_Errors, na.rm=TRUE), 4)
TBA_SE_Val_Errors = apply(TBA_Valid_Errors, 1, sd)/sqrt(Num_Folds)

#### HFG LASSO Selected Lambdas

# Lambda for Minimum Validation Error:
Opt_Min_Lambda = Lambda_Set[which.min(Mean_Val_Errors)]

# Lambda for Min+1SE Validation Error
Min_Plus1SE_Check = Mean_Val_Errors < Mean_Val_Errors[which.min(Mean_Val_Errors)] + SE_Val_Errors[which.min(Mean_Val_Errors)]

#Find final index that falls below the threshold
Fus_Loc = which(Min_Plus1SE_Check)
Opt_1SE_Lambda = Lambda_Set[Fus_Loc[length(Fus_Loc)]]

Opt_Lambdas = c(Opt_Min_Lambda, Opt_1SE_Lambda)
names(Opt_Lambdas) = c("Opt Min Lambda", "Opt 1SE Lambda")


#### TBA Selected Lambdas

# Lambda for Minimum Validation Error:
TBA_Opt_Min_Lambda = TBA_Lambdas[which.min(TBA_Mean_Val_Errors)]

# Lambda for Min+1SE Validation Error
Min_Plus1SE_Check = TBA_Mean_Val_Errors < TBA_Mean_Val_Errors[which.min(TBA_Mean_Val_Errors)] + TBA_SE_Val_Errors[which.min(TBA_Mean_Val_Errors)]

#Find final index that falls below the threshold
Fus_Loc = which(Min_Plus1SE_Check)
TBA_Opt_1SE_Lambda = TBA_Lambdas[Fus_Loc[1]]

TBA_Opt_Lambdas = c(TBA_Opt_Min_Lambda, TBA_Opt_1SE_Lambda)
names(TBA_Opt_Lambdas) = c("Opt Min Lambda", "Opt 1SE Lambda")


## Fit the Final Models
Thetas_Set = HFG_LASSO_Fit(HierDesignMatrix, Response, Categ_Inds, Categ_Pair_Inds, S_Matrices_Weights,
                                      Lambda_Set, Proximal_Grad_Tol, Gamma)

TBA_Final_Model = rarefit(y=Response, X=HierDesignMatrix, A = A, intercept = FALSE, lambda = TBA_Lambdas, alpha = 1)


#OLS Results for Coefficients (temporary baseline comparison)
DF = as.data.frame(cbind(Response, DesignMatrix))
names(DF) = c("Response", paste("X", rep(1:ncol(DesignMatrix)), sep=""))

LM = lm(Response~.-1, data=DF)
LM_Coeffs = coefficients(LM)


#### Compute the predictive performance (RMSE) of the final models on the holdout set that was constructed, and store them

OLS_RMSE = sqrt(sum((Valid_DesignMatrix %*% LM_Coeffs - Validation_Response)^2)/nrow(Valid_DesignMatrix))

Opt_Min_Lambda_RMSE = sqrt(sum((Valid_DesignMatrix %*% Thetas_Set[,paste("Lambda_", Opt_Lambdas[1], sep="")] - Validation_Response)^2)/nrow(Valid_DesignMatrix))
Opt_1SE_Lambda_RMSE = sqrt(sum((Valid_DesignMatrix %*% Thetas_Set[,paste("Lambda_", Opt_Lambdas[2], sep="")] - Validation_Response)^2)/nrow(Valid_DesignMatrix))

Oracle_RMSE = sqrt(sum((Valid_DesignMatrix %*% Beta_Comp - Validation_Response)^2)/nrow(Valid_DesignMatrix))

TBA_Min_Lambda_RMSE = sqrt(sum((Valid_DesignMatrix %*% TBA_Final_Model$beta[[1]][,which(TBA_Lambdas == TBA_Opt_Lambdas[1])] - Validation_Response)^2)/nrow(Valid_DesignMatrix))
TBA_1SE_Lambda_RMSE = sqrt(sum((Valid_DesignMatrix %*% TBA_Final_Model$beta[[1]][,which(TBA_Lambdas == TBA_Opt_Lambdas[2])] - Validation_Response)^2)/nrow(Valid_DesignMatrix))


Final_Test_Errors = c(OLS_RMSE, Opt_Min_Lambda_RMSE, Opt_1SE_Lambda_RMSE, Oracle_RMSE, TBA_Min_Lambda_RMSE, TBA_1SE_Lambda_RMSE)

names(Final_Test_Errors) = c("OLS Test Error", "Min CV HFGL Test Error", "1SE CV HFGL Test Error", "Oracle Test Error", "Min CV TBA Test Error", "1SE CV TBA Test Error")


#Compile and Store Results

Thetas_Set = cbind(True_Coeffs = Beta_Comp, LM_Coeffs=LM_Coeffs, data.frame(Thetas_Set))

HFGL_CV_Model_Results = list(Coefficients = Thetas_Set, Design_Matrix = HierDesignMatrix, Response = Response, Selected_Lambda = Opt_Lambdas,
                             Validation_Errors = Valid_Errors, Final_Test_Errors = Final_Test_Errors,
                             Tree_Structure = Assignment_List, S_Matrices = S_Matrices, SNR=SNR,
                             Active_Nodes = unique(c(Active_Nodes_Set, Active_Leaves+sum(Categ_Counts[1:(length(Categ_Counts)-1)]))),
                             TBA_Opt_Lambdas = TBA_Opt_Lambdas, TBA_Final_Model = TBA_Final_Model)


save(HFGL_CV_Model_Results, file=paste(Model_Results_Loc, "HFGL_CV_Seed_", Seed, "_SNR_", SNR, "_", Weight_Scheme, ".Rdata", sep=""))


