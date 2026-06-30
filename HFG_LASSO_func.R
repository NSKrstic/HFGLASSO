#################################################
# Script: HFG_LASSO_func.R
# Author: Nikolas Krstic
# Purpose: Functions to fit Hierarchical Fused Group LASSO
#################################################



################################################
## Function to fit Hierarchical Fused Group LASSO (WITHOUT ADDITIONAL COVARIATES) (using Simultaneous-Direction Method of Multipliers as found in Section 7 of https://arxiv.org/pdf/0912.3522.pdf)
## This version of the function also uses an Rcpp function for Step 3 of Algorithm 1 in the manuscript, please see the function "Thetas_SDMM_Computation" in "HFGL_Fitting_Algorithm_Iteration_StepV3.cpp"

## Note that because the least squares loss in the HFG LASSO objective function is not scaled by the sample size, then this means a set containing large $\lambda$ is required for large sample sizes.
## There is an argument called "Scale_Lambda" that will automatically scale the grid for you, if you prefer.

## Arguments:

## HierDesignMatrix = The binary indicator design matrix of the hierarchical categorical predictor (specifically for the most granular/lowest layer of the hierarchy)
## Response = The vector of responses for each observation
## Categ_Inds = A list containing the indices of the non-zero columns in each of the S Matrices
## Categ_Pair_Inds = A list containing two-column matrices, which contain all of the fusion pair indices for each given S Matrix
## S_Matrices_Weights = The vector of weights for each penalization term
## Lambda_Set = The vector of lambda hyperparameters to fit the model with
## Proximal_Grad_Tol = Tolerance for the overall proximal gradient algorithm
## Gamma = The postive gamma hyperparameter for SDMM (can be set to 1 for simplicity, but empirically values between 0.5 and 5 work well depending on the problem)
## Scale_Lambda = A boolean (TRUE/FALSE) that indicates where to scale Lambda_Set by the sample size (i.e. multiply by sample size). Default is FALSE.
## Print_Lambda = A boolean (TRUE/FALSE) that indicates whether to print the current lambda value being used for solving HFG LASSO. Default is FALSE.


HFG_LASSO_Fit = function(HierDesignMatrix, Response, Categ_Inds, Categ_Pair_Inds, S_Matrices_Weights,
                                  Lambda_Set, Proximal_Grad_Tol, Gamma, Scale_Lambda=FALSE, Print_Lambda=FALSE){
  
  #Scale Lambda Values?
  if(Scale_Lambda){
    Lambda_Set = Lambda_Set*length(Response)
  }
  
  #Store within this list the thetas obtained for each given lambda hyperparameter in the lambda set.
  Thetas_List = list()
  
  #################################
  #Compute inverse of Q matrix for algorithm
  print("Computing Q Matrix for SDMM...")
  
  # Compute X^T*X or S^T*S for all the matrices
  XTX_Mats = list()
  
  for(i in 1:(length(Categ_Inds)+1)){
    
    if(i==1){
      
      Categ_Freqs = colSums(HierDesignMatrix)
      XTX_Mats[[i]] = diag(Categ_Freqs)
      
    }else{
      # Assign the initial entries of the matrix
      Int_XTX_Mat = matrix(data=0, nrow=ncol(HierDesignMatrix), ncol=ncol(HierDesignMatrix))
      Int_XTX_Mat[Categ_Pair_Inds[[i-1]]] = -1
      Int_XTX_Mat[Categ_Pair_Inds[[i-1]][,c(2,1), drop=FALSE]] = -1
      
      # Populate the diagonal of the matrix with the correct values in the correct locations
      Diag_Ents = table(Categ_Pair_Inds[[i-1]])
      Diag_Vals = as.numeric(names(Diag_Ents))
      diag(Int_XTX_Mat)[Diag_Vals] = Diag_Ents
      XTX_Mats[[i]] = Int_XTX_Mat
    }
    
  }
  
  # Conduct the summation and solve for inverse
  Q_Mat = Reduce('+', XTX_Mats)
  Q_Inv_Mat = solve(Q_Mat)
  
  # Remove XTX Matrices and clean-up
  rm(XTX_Mats)
  gc()
  
  ##############################
  ### Additional Auxiliary Items (Used To Optimize Speed For the Rest of the Algorithm, i.e. avoid redundant calculations)
  
  print("Compiling Auxiliary Algorithm Items...")
  
  #### Design Matrix 1 indices and Gamma*Response
  DM1s = which(HierDesignMatrix==1, arr.ind=TRUE)
  DM1s_Inds = lapply(1:ncol(HierDesignMatrix), function(x){DM1s[which(DM1s[,2]==x),1]})
  DM1s_RowByCol_Inds = DM1s[order(DM1s[,1]),2]
  
  Gamma_Response = Gamma*Response
  
  #### Row Numbers of Each Matrix
  S_Mat_Rows = sapply(Categ_Pair_Inds, nrow)
  S_Mat_Rows = c(nrow(HierDesignMatrix), S_Mat_Rows)
  
  ### Dimensions and Sizes
  NCol_HDM = ncol(HierDesignMatrix)
  Length_Categ_Inds = length(Categ_Inds)
  
  #### Set Up 1s and -1s row indices for each S Matrix column (improves computation performance later when computing YZ_Diffs)
  Ones_RowInds = list()
  NegOnes_RowInds = list()
  
  for(q in 2:length(S_Mat_Rows)){
    
    #Identify current set of fusion pairs
    CPS = Categ_Pair_Inds[[q-1]]
    C_Inds = Categ_Inds[[q-1]]
    
    Ones_RowInds[[q-1]] = lapply(C_Inds, function(x){which(CPS[,1]==x)})
    names(Ones_RowInds[[q-1]]) = C_Inds
    NegOnes_RowInds[[q-1]] = lapply(C_Inds, function(x){which(CPS[,2]==x)})
    names(NegOnes_RowInds[[q-1]]) = C_Inds
  }
  
  #### Character version of Categ_Inds
  Char_Categ_Inds = lapply(Categ_Inds, as.character)
  
  
  ############################
  #Create the y and z sets of vectors for Algorithm 7.9 of Combettes and Pesquet (2009)
  ##set.seed(100)
  
  print("Initializing Y Vectors, Z Vectors, and Thetas...")
  
  Y_Vecs = lapply(S_Mat_Rows, function(x){rnorm(x)})
  Z_Vecs = lapply(S_Mat_Rows, function(x){rnorm(x)})
  
  YZ_Diffs = list()
  TEMP_YZ_Diffs = lapply(1:length(Y_Vecs), function(x){Y_Vecs[[x]]-Z_Vecs[[x]]})
  
  for(i in 1:length(S_Mat_Rows)){
    
    if(i==1){
      YZ_Diffs[[i]] = sapply(1:ncol(HierDesignMatrix), function(x){sum((TEMP_YZ_Diffs[[i]])[DM1s_Inds[[x]]])})
      
    }else{
      #Identify current set of fusion pairs
      C_Inds = Categ_Inds[[i-1]]
      Char_C_Inds = Char_Categ_Inds[[i-1]]
      
      #Compute the "y-z" differences by simplifying the above expression
      YZ_Diffs[[i]] = rep(0, ncol(HierDesignMatrix))
      
      YZ_Diffs[[i]][C_Inds] = sapply(Char_C_Inds,
                                     function(x){sum((TEMP_YZ_Diffs[[i]])[Ones_RowInds[[i-1]][[x]]])-
                                         sum((TEMP_YZ_Diffs[[i]])[NegOnes_RowInds[[i-1]][[x]]])})
    }
    
  }
  
  #Compute Initial Thetas
  Thetas = Q_Inv_Mat %*% Reduce('+', YZ_Diffs)
  
  #### Added -1 to function in lapply for C++ indexing
  Categ_Inds_C = lapply(Categ_Inds, function(x){x-1})
  
  #################################################
  print("Commencing Iterative Looping...")
  
  for(p in 1:length(Lambda_Set)){
    
    #Select the current lambda
    Curr_Lambda = Lambda_Set[p]
    if(Print_Lambda){
      print(paste("Current Lambda:", Curr_Lambda, sep=""))
    }
    
    #Auxiliary Item (to improve computational speed)
    Gamma_Lambda_S_Weights = Gamma*Curr_Lambda*(S_Matrices_Weights)
    
    #Primary Function for Iterating During Step 3 of Algorithm 1, Outputs a List of Final Results.
    Final_Results = Thetas_SDMM_Computation(Gamma, Gamma_Lambda_S_Weights, Thetas, Categ_Inds_C,
                                            Char_Categ_Inds, DM1s_RowByCol_Inds, Z_Vecs, Y_Vecs, Gamma_Response,
                                            Categ_Pair_Inds, DM1s_Inds, Ones_RowInds, NegOnes_RowInds, Q_Inv_Mat,
                                            Proximal_Grad_Tol, NCol_HDM, Length_Categ_Inds)
    
    Thetas = Final_Results[[1]]
    Y_Vecs = Final_Results[[2]]
    Z_Vecs = Final_Results[[3]]
    
    Thetas_List[[p]] = Thetas
    
  }
  
  Thetas_Set = do.call(cbind, Thetas_List)
  colnames(Thetas_Set) = paste("Lambda_", Lambda_Set, sep="")
  
  #Output a matrix of coefficient solutions for the different lambda values.
  return(Thetas_Set)
  
}






################################################
## Function to fit Hierarchical Fused Group LASSO - VERSION: Additional Non-Hierarchical predictors in Design matrix
## (using Simultaneous-Direction Method of Multipliers as found in Section 7 of https://arxiv.org/pdf/0912.3522.pdf)
## This version of the function also uses an Rcpp function for Step 3 of Algorithm 1 in the manuscript, please see the function "Thetas_SDMM_Computation" in "HFGL_Fitting_Algorithm_Iteration_StepV3.cpp"


## Note that because the least squares loss in the HFG LASSO objective function is not scaled by the sample size, then this means a set containing large $\lambda$ is required for large sample sizes.
## There is an argument called "Scale_Lambda" that will automatically scale the grid for you, if you prefer.

## Arguments:

## HierDesignMatrix = The complete design matrix (first set of columns are for the hierarchical predictors, second set of columns are for the other covariates)
## HDM_NHCP_Length = The number of columns (starting from the final column) in the complete design matrix that correspond to the other covariates
## Response = The vector of responses for each observation
## Categ_Inds = A list containing the indices of the non-zero columns in the S Matrices
## Categ_Pair_Inds = A list containing two-column matrices, which contain all of the fusion pair indices for each given S Matrix
## S_Matrices_Weights = The vector of weights for each penalization term
## Lambda_Set = The vector of lambda hyperparameters to fit the model with
## Proximal_Grad_Tol = Tolerance for the overall proximal gradient algorithm
## Gamma = The positive gamma hyperparameter for SDMM (likely just set to 1 for simplicity)
## Scale_Lambda = A boolean (TRUE/FALSE) that indicates where to scale Lambda_Set by the sample size (i.e. multiply by sample size). Default is FALSE.
## Print_Lambda = A boolean (TRUE/FALSE) that indicates whether to print the current lambda value being used for solving HFG LASSO. Default is FALSE.


HFG_LASSO_NonHCPVer_Fit = function(HierDesignMatrix, HDM_NHCP_Length, Response, Categ_Inds, Categ_Pair_Inds,
                                    S_Matrices_Weights, Lambda_Set, Proximal_Grad_Tol, Gamma, Scale_Lambda=FALSE, Print_Lambda=FALSE){
  
  #Scale Lambda Values?
  if(Scale_Lambda){
    Lambda_Set = Lambda_Set*length(Response)
  }
  
  #Store within this list the thetas obtained for each given lambda hyperparameter in the lambda set.
  Thetas_List = list()
  
  #Identify the indices for the HCP and the Non-HCP covariates
  HCP_Inds = 1:(ncol(HierDesignMatrix)-HDM_NHCP_Length)
  NonHCP_Inds = ((ncol(HierDesignMatrix)-HDM_NHCP_Length)+1):ncol(HierDesignMatrix)
  
  #Create separate matrices for the HCP and Non-HCP
  HCP_Des_Mat = HierDesignMatrix[,HCP_Inds]
  NHCP_Des_Mat = HierDesignMatrix[,NonHCP_Inds]
  
  #################################
  #Compute inverse of Q matrix for algorithm
  print("Computing Q Matrix for SDMM...")
  
  XTX_Mats = list()
  
  # Compute X^T*X or S^T*S for all of the matrices
  for(i in 1:(length(Categ_Inds)+1)){
    
    if(i==1){
      Categ_Freqs = colSums(HCP_Des_Mat)
      XTX_Mats[[i]] = cbind(rbind(diag(Categ_Freqs), t(NHCP_Des_Mat) %*% HCP_Des_Mat),
                            rbind(t(HCP_Des_Mat) %*% NHCP_Des_Mat, t(NHCP_Des_Mat) %*% NHCP_Des_Mat))
      
    }else{
      # Assign the initial entries of the matrix
      Int_XTX_Mat = matrix(data=0, nrow=ncol(HierDesignMatrix), ncol=ncol(HierDesignMatrix))
      Int_XTX_Mat[Categ_Pair_Inds[[i-1]]] = -1
      Int_XTX_Mat[Categ_Pair_Inds[[i-1]][,c(2,1), drop=FALSE]] = -1
      
      # Populate the diagonal of the matrix with the correct values in the correct locations
      Diag_Ents = table(Categ_Pair_Inds[[i-1]])
      Diag_Vals = as.numeric(names(Diag_Ents))
      diag(Int_XTX_Mat)[Diag_Vals] = Diag_Ents
      XTX_Mats[[i]] = Int_XTX_Mat
    }
    
  }
  
  # Conduct the summation and solve for inverse
  Q_Mat = Reduce('+', XTX_Mats)
  Q_Inv_Mat = solve(Q_Mat)
  
  #Remove XTX Matrices and clean-up
  rm(XTX_Mats)
  gc()
  
  ##############################
  ### Additional Auxiliary Items (Used To Optimize Speed For the Rest of the Algorithm, i.e. avoid redundant calculations)
  
  print("Compiling Auxiliary Algorithm Items...")
  
  #### Design Matrix 1 indices and Gamma*Response
  DM1s = which(HCP_Des_Mat==1, arr.ind=TRUE)
  DM1s_Inds = lapply(1:(ncol(HCP_Des_Mat)), function(x){DM1s[which(DM1s[,2]==x),1]})
  DM1s_RowByCol_Inds = DM1s[order(DM1s[,1]),2]
  
  Gamma_Response = Gamma*Response
  
  #### Row Numbers of Comp_List
  S_Mat_Rows = sapply(Categ_Pair_Inds, nrow)
  S_Mat_Rows = c(nrow(HierDesignMatrix), S_Mat_Rows)
  
  ### Dimensions and Sizes
  NCol_HDM = ncol(HierDesignMatrix)
  Length_Categ_Inds = length(Categ_Inds)
  
  #### Set Up 1s and -1s row indices for each S Matrix column (improves computation performance later when computing YZ_Diffs)
  Ones_RowInds = list()
  NegOnes_RowInds = list()
  
  for(q in 2:length(S_Mat_Rows)){
    
    #Identify current set of fusion pairs
    CPS = Categ_Pair_Inds[[q-1]]
    C_Inds = Categ_Inds[[q-1]]
    
    Ones_RowInds[[q-1]] = lapply(C_Inds, function(x){which(CPS[,1]==x)})
    names(Ones_RowInds[[q-1]]) = C_Inds
    NegOnes_RowInds[[q-1]] = lapply(C_Inds, function(x){which(CPS[,2]==x)})
    names(NegOnes_RowInds[[q-1]]) = C_Inds
  }
  
  #### Character version of Categ_Inds
  Char_Categ_Inds = lapply(Categ_Inds, as.character)
  
  
  ############################
  #Create the y and z sets of vectors for Algorithm 7.9 of Combettes and Pesquet (2009)
  ##set.seed(100)
  
  print("Initializing Y Vectors, Z Vectors, and Thetas...")
  
  #Randomly initialize the Y and Z Vectors
  Y_Vecs = lapply(S_Mat_Rows, function(x){rnorm(x)})
  Z_Vecs = lapply(S_Mat_Rows, function(x){rnorm(x)})
  
  YZ_Diffs = list()
  TEMP_YZ_Diffs = lapply(1:length(Y_Vecs), function(x){Y_Vecs[[x]]-Z_Vecs[[x]]})
  
  for(i in 1:length(S_Mat_Rows)){
    
    if(i==1){
      YZ_Diffs[[i]] = c(sapply(1:(ncol(HierDesignMatrix)-HDM_NHCP_Length), function(x){sum((TEMP_YZ_Diffs[[i]])[DM1s_Inds[[x]]])}), 
                        t(NHCP_Des_Mat) %*% TEMP_YZ_Diffs[[i]])
      
    }else{
      #Identify current set of fusion pairs
      C_Inds = Categ_Inds[[i-1]]
      Char_C_Inds = Char_Categ_Inds[[i-1]]
      
      #Compute the "y-z" differences by simplifying the above expression
      YZ_Diffs[[i]] = rep(0, ncol(HierDesignMatrix))
      
      YZ_Diffs[[i]][C_Inds] = sapply(Char_C_Inds,
                                     function(x){sum((TEMP_YZ_Diffs[[i]])[Ones_RowInds[[i-1]][[x]]])-
                                         sum((TEMP_YZ_Diffs[[i]])[NegOnes_RowInds[[i-1]][[x]]])})
    }
    
  }
  
  #Compute Initial Thetas
  Thetas = Q_Inv_Mat %*% Reduce('+', YZ_Diffs)
  
  #### !!!!!! ADDED -1 to function in lapply for C++ indexing
  Categ_Inds_C = lapply(Categ_Inds, function(x){x-1})
  
  #################################################
  print("Commencing Iterative Looping...")
  
  for(p in 1:length(Lambda_Set)){
    
    #Select the current lambda and reset iteration count
    Curr_Lambda = Lambda_Set[p]
    if(Print_Lambda){
      print(paste("Current Lambda:", Curr_Lambda, sep=""))
    }
    
    #Auxiliary Item (to improve computational speed)
    Gamma_Lambda_S_Weights = Gamma*Curr_Lambda*(S_Matrices_Weights)
    
    #Primary Function for Iterating During Step 3 of Algorithm 1, Outputs a List of Final Results.
    Final_Results = Thetas_SDMM_Computation_NonHCPVer(Gamma, Gamma_Lambda_S_Weights, Thetas, Categ_Inds_C,
                                                      Char_Categ_Inds, DM1s_RowByCol_Inds, Z_Vecs, Y_Vecs, Gamma_Response,
                                                      Categ_Pair_Inds, DM1s_Inds, Ones_RowInds, NegOnes_RowInds, Q_Inv_Mat,
                                                      Proximal_Grad_Tol, NCol_HDM, Length_Categ_Inds, HDM_NHCP_Length,
                                                      NonHCP_Inds, HCP_Des_Mat, NHCP_Des_Mat)
    
    Thetas = Final_Results[[1]]
    Y_Vecs = Final_Results[[2]]
    Z_Vecs = Final_Results[[3]]
    
    Thetas_List[[p]] = Thetas
    
  }
  
  Thetas_Set = do.call(cbind, Thetas_List)
  colnames(Thetas_Set) = paste("Lambda_", Lambda_Set, sep="")
  
  #Output a matrix of coefficient solutions for the different lambda values.
  return(Thetas_Set)
  
}





################################################
## Function to identify Lambda_Max for Hierarchical Fused Group LASSO, using the CVXR R package

## Note that because the least squares loss in the HFG LASSO objective function is not scaled by the sample size,
## then this means a Lambda_Max that is very large is not unusual for large sample sizes, and can be manually scaled down by dividing by the sample size.

## Arguments:

## HierDesignMatrix = The complete design matrix (first set of columns are for the hierarchical predictors, second set of columns are for the other covariates)
## Response = The vector of responses for each observation
## S_Matrices = The S Matrices within the HFG LASSO penalty
## S_Matrices_Weights = The vector of weights for each penalization term
## HDM_NHCP_Length = The number of columns (starting from the final column) in the complete design matrix that correspond to the other covariates (zero by default)

HFG_LASSO_Lambda_Max = function(DesignMatrix, Response, S_Matrices, S_Matrices_Weights, HDM_NHCP_Length=0){
  
  #Compute the "b" vector from Theorem 1
  if(HDM_NHCP_Length==0){
    
    LHS_Vec = t(DesignMatrix) %*% (Response-DesignMatrix %*% rep(mean(Response), ncol(DesignMatrix)))
    
  }else{
    
    NHCP_DesMat = DesignMatrix[,(ncol(DesignMatrix)-HDM_NHCP_Length+1):(ncol(DesignMatrix))]
    LM_DF = data.frame(Response=Response, NHCP_DesMat)
    
    Fitted_Vals = cbind(rep(1, nrow(NHCP_DesMat)), NHCP_DesMat) %*% lm(Response~.,data=LM_DF)$coefficients
    LHS_Vec = t(DesignMatrix) %*% (Response-Fitted_Vals)
    
  }
  
  #Compute the components of the "A" matrix from Theorem 1
  SW_Mats = lapply(1:length(S_Matrices), function(x){t(S_Matrices[[x]])*S_Matrices_Weights[[x]]})
  
  SW_Mat_Prods = list()
  Z_Lambda_Vecs = list()
  
  for(i in 1:length(S_Matrices)){
    
    Z_Lambda_Vecs[[i]] = Variable(nrow(S_Matrices[[i]]))
    SW_Mat_Prods[[i]] = SW_Mats[[i]] %*% Z_Lambda_Vecs[[i]]
    
  }
  
  #Construct the cost/objective function
  Cost_Norms = (lapply(Z_Lambda_Vecs, function(x){cvxr_norm(x, p=2)}))
  Cost = do.call(max_elemwise, Cost_Norms)
  
  #Construct the constraints
  Constraints_All = list((LHS_Vec)-Reduce('+', SW_Mat_Prods)==0)
  
  #Solve for Lambda_max
  Objective = Minimize(Cost)
  Problem = Problem(Objective, Constraints_All)
  Solution = solve(Problem, feastol=1e-10, abstol=1e-10, reltol=1e-10, num_iter=10000, verbose=TRUE, solver="SCS")
  
  
  #Check accuracy of solution
  #Final_Solution_List = lapply(Z_Lambda_Vecs, function(x){Solution$getValue(x)})
  
  #LHS_Vec - (SW_Mats[[1]] %*% Final_Solution_List[[1]] + SW_Mats[[2]] %*% Final_Solution_List[[2]] + SW_Mats[[3]] %*% Final_Solution_List[[3]] +
  #  SW_Mats[[4]] %*% Final_Solution_List[[4]] + SW_Mats[[5]] %*% Final_Solution_List[[5]] + SW_Mats[[6]] %*% Final_Solution_List[[6]] +
  #  SW_Mats[[7]] %*% Final_Solution_List[[7]] + SW_Mats[[8]] %*% Final_Solution_List[[8]])
  
  return(Solution$value)
  
}









