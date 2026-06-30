#################################################
# Script: simulation_study_SCAL.R
# Author: Nikolas Krstic
# Purpose: Perform "Hierarchical Fused Group LASSO", with hold-out test set AND comparison with Tree-Based Aggregation by Yan and Bien (2021)
#################################################

Rcpp::sourceCpp("HFGL_Fitting_Algorithm_1_Step_3_Functions.cpp")
source("HFG_LASSO_func.R")
source("Scalability_Simulations/simulate_data_SCAL.R")

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
S_Matrices_Weights = sapply(S_Matrix_Layers, function(x){if(x==3){1}else if(x==2){1/5}else if(x==1){1/10}})
S_Matrices_Weights = S_Matrices_Weights/sum(S_Matrices_Weights)

#Re-assign design matrix to different name (to distinguish from non-HCP design matrix that could be considered here)
HierDesignMatrix = DesignMatrix
#Tolerance to use when assessing convergence of the HFG LASSO fitting algorithm
if(Sim_Scal_Type=="SS"){
  Proximal_Grad_Tol = 1e-8
}else if(Sim_Scal_Type=="HCPC"){
  Proximal_Grad_Tol = 1e-5
}

#SDMM Hyperparameter Setting
if(Sim_Scal_Type=="SS"){
  Gamma = 0.5
}else if(Sim_Scal_Type=="HCPC"){
  Gamma = 0.5*Cardinality/10
}

## Set the Lambda Grid for HFG LASSO and the Number of Lambda to Use for TBA
if(Sim_Scal_Type=="SS"){
  #Create the S Matrices (to use if necessary, such as for computing Lambda_Max)
  S_Matrices = lapply(Categ_Pair_Inds, function(x){t(apply(x, 1, function(y){S_Row=rep(0, ncol(HierDesignMatrix)); S_Row[y[1]]=1; S_Row[y[2]]=-1; return(S_Row)}))})
  
  ## Identify Lambda_Max using Theorem 1
  Lambda_Max = HFG_LASSO_Lambda_Max(DesignMatrix, Response, S_Matrices, S_Matrices_Weights)
  
  Lambda_Set = exp(seq(0, log(0.3*Lambda_Max), length.out=60))
  
}else if(Sim_Scal_Type=="HCPC"){
  Lambda_Set = exp(seq(0, 10, length.out=60))
  
}

TBA_Lambda_Count = 60

#Obtain the category frequencies, the category indices for non-zero columns of the S Matrices, and the fusion pairs for each of the S Matrices
Categ_Freqs = colSums(HierDesignMatrix)
###Categ_Inds = lapply(S_Matrices, function(x){which(apply(x, 2, function(x){all(x==0)})==FALSE)})
###Categ_Pair_Inds = lapply(S_Matrices, function(x){t(apply(x, 1, function(x){which(x!=0)}))})


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
## Hold-Out Set Approach

## Fit the Final Models

### Measure the computation time:

HFG_LASSO_Start = Sys.time()

Thetas_Set = HFG_LASSO_Fit(HierDesignMatrix, Response, Categ_Inds, Categ_Pair_Inds, S_Matrices_Weights,
                           Lambda_Set, Proximal_Grad_Tol, Gamma)

HFG_LASSO_End = Sys.time()
HFG_LASSO_Run_Time = as.numeric(difftime(HFG_LASSO_End, HFG_LASSO_Start, units="secs"))


TBA_Start = Sys.time()

TBA_Final_Model = rarefit(y=Response, X=HierDesignMatrix, A = A, intercept = FALSE, lambda = TBA_Lambdas, alpha = 1)

TBA_End = Sys.time()
TBA_Run_Time = as.numeric(difftime(TBA_End, TBA_Start, units="secs"))


#### Compute the Validation Set Performances

HFG_LASSO_Preds = (Valid_DesignMatrix %*% Thetas_Set - replicate(length(Lambda_Set), Validation_Response, simplify=TRUE))
TBA_Preds = (Valid_DesignMatrix %*% TBA_Final_Model$beta[[1]] - replicate(length(TBA_Lambdas), Validation_Response, simplify=TRUE))
Mean_Val_Errors = apply(HFG_LASSO_Preds, 2, function(x){sqrt(sum(x^2)/nrow(Valid_DesignMatrix))})
TBA_Mean_Val_Errors = apply(TBA_Preds, 2, function(x){sqrt(sum(x^2)/nrow(Valid_DesignMatrix))})


#### HFG LASSO Selected Lambda

# Lambda for Minimum Validation Error:
Opt_Min_Lambda = Lambda_Set[which.min(Mean_Val_Errors)]

#### TBA Selected Lambda

# Lambda for Minimum Validation Error:
TBA_Opt_Min_Lambda = TBA_Lambdas[which.min(TBA_Mean_Val_Errors)]



#OLS Results for Coefficients (temporary baseline comparison)
DF = as.data.frame(cbind(Response, DesignMatrix))
names(DF) = c("Response", paste("X", rep(1:ncol(DesignMatrix)), sep=""))

LM = lm(Response~.-1, data=DF)
LM_Coeffs = coefficients(LM)


#### Compute the predictive performance (RMSE) of the final models on the test holdout set that was also constructed, and store them

OLS_RMSE = sqrt(sum((Test_DesignMatrix %*% LM_Coeffs - Test_Response)^2)/nrow(Test_DesignMatrix))

Opt_Min_Lambda_RMSE = sqrt(sum((Test_DesignMatrix %*% Thetas_Set[,paste("Lambda_", Opt_Min_Lambda, sep="")] - Test_Response)^2)/nrow(Test_DesignMatrix))

Oracle_RMSE = sqrt(sum((Test_DesignMatrix %*% Beta_Comp - Test_Response)^2)/nrow(Test_DesignMatrix))

TBA_Min_Lambda_RMSE = sqrt(sum((Test_DesignMatrix %*% TBA_Final_Model$beta[[1]][,which(TBA_Lambdas == TBA_Opt_Min_Lambda)] - Test_Response)^2)/nrow(Test_DesignMatrix))


Final_Test_Errors = c(OLS_RMSE, Opt_Min_Lambda_RMSE, Oracle_RMSE, TBA_Min_Lambda_RMSE)

names(Final_Test_Errors) = c("OLS Test Error", "Min HFGL Test Error", "Oracle Test Error", "Min TBA Test Error")


#Compile and Store Results

Thetas_Set = cbind(True_Coeffs = Beta_Comp, LM_Coeffs=LM_Coeffs, data.frame(Thetas_Set))

HFGL_Model_Results = list(Coefficients = Thetas_Set, Selected_Lambda = Opt_Min_Lambda,
                             Final_Test_Errors = Final_Test_Errors,
                             Tree_Structure = Assignment_List, SNR=SNR,
                             Active_Nodes = unique(c(Active_Nodes_Set, Active_Leaves+sum(Categ_Counts[1:(length(Categ_Counts)-1)]))),
                             TBA_Opt_Lambda = TBA_Opt_Min_Lambda, TBA_Final_Model_Coefficients = TBA_Final_Model$beta[[1]], 
                             TBA_Final_Model_Lambda = TBA_Final_Model$lambda, HFG_LASSO_Run_Time = HFG_LASSO_Run_Time,
                             TBA_Run_Time = TBA_Run_Time)


if(Sim_Scal_Type=="SS"){
  save(HFGL_Model_Results, file=paste(Model_Results_Loc, "HFGL_Seed_", Seed, "_Samp_Size_", Samp_Size, "_SNR_", SNR, "_", Weight_Scheme, ".Rdata", sep=""))
}else if(Sim_Scal_Type=="HCPC"){
  save(HFGL_Model_Results, file=paste(Model_Results_Loc, "HFGL_Seed_", Seed, "_Cardinality_", Cardinality, "_SNR_", SNR, "_", Weight_Scheme, ".Rdata", sep=""))
}

