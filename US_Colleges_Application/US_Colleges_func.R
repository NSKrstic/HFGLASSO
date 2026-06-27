#################################################
# Script: US_Colleges_func.R
# Author: Nikolas Krstic
# Purpose: Functions for conducting application of different methods to US Colleges Dataset
#################################################


#######################################################################################
# Function to compute validation errors, optimal lambdas (MCVE and 1SE Rule) and other elements for each method type
## WILL ONLY WORK WHEN OTHER ARGUMENTS FOR SOME OF THE METHODS (e.g., HFG LASSO) ARE AVAILABLE IN THE GLOBAL ENVIRONMENT (e.g., S_Matrices_Weights)

## Arguments:

## Method - A character which indicates which method to pursue (one of "HFG LASSO", "TBA")
## Lambda_Grid - The grid of lambda over which to fit the method
## Folds - A list of fold indices, where each element of the list corresponds to the observations indices of that fold from the dataset
## HierDesignMatrix - The binary indicator design matrix of the HCP (namely from the Training Set)
## Response - A numeric vector of the response

PredPerf_CV_Computation = function(Method, Lambda_Grid, Folds, HierDesignMatrix, Response){
  
  #Conducting Cross-Validation, Store the Validation Errors at Each Fold for the Method
  Test_Set_Validation_Errors_Set = list()
  
  for(q in 1:length(Folds)){
    
    print(paste("Current Fold: ", q, sep=""))
    
    Curr_Test_Folds_Indices = Folds[[q]]
    Curr_Train_Folds_Indices = as.vector(unlist(Folds[-q]))
    
    #Build Train and Test Predictor Matrices and Responses
    Train_DesignMatrix = HierDesignMatrix[Curr_Train_Folds_Indices,]
    Train_Response = Response[Curr_Train_Folds_Indices]
    Test_DesignMatrix = HierDesignMatrix[Curr_Test_Folds_Indices,]
    Test_Response = Response[Curr_Test_Folds_Indices]
    
    #Build a temporary data.frame (for use in some specific methods, based on their expected arguments)
    Train_Group_DF = data.frame(Response = Train_Response, State = as.factor(colnames(Train_DesignMatrix)[Train_DesignMatrix %*% 1:ncol(Train_DesignMatrix)]))
    
    if(Method=="HFG LASSO"){
      
      Thetas_Set = HFG_LASSO_Fit(Train_DesignMatrix, Train_Response, Categ_Inds, Categ_Pair_Inds,
                                      S_Matrices_Weights, Lambda_Grid, Proximal_Grad_Tol, Gamma)
      Final_Coeffs = Thetas_Set
      
    }else if(Method=="TBA"){
      
      TBA_Model = rarefit(y=Train_Response, X=Train_DesignMatrix, A = A, intercept = FALSE, lambda = Lambda_Grid, alpha = 1)
      Final_Coeffs = TBA_Model$beta[[1]]
      
    }
    
    #Predictions And Validation Errors
    Test_Preds = (Test_DesignMatrix %*% Final_Coeffs - replicate(length(Lambda_Grid), Test_Response))
    Test_Set_Validation_Errors_Set[[q]] = apply(Test_Preds, 2, function(x){sqrt(sum(x^2)/nrow(Test_DesignMatrix))})
    
  }
  
  
  ######################################################
  # Identify Optimal Lambdas and then Fit Final Models
  
  College_Group_DF = data.frame(Response = Response, State = as.factor(colnames(HierDesignMatrix)[HierDesignMatrix %*% 1:ncol(HierDesignMatrix)]))
  
  #Compute summaries of validation errors
  Valid_Errors = do.call(cbind, Test_Set_Validation_Errors_Set)
  Mean_Val_Errors = rowMeans(Valid_Errors, na.rm=TRUE)
  SE_Val_Errors = apply(Valid_Errors, 1, sd)/sqrt(length(Folds))
  
  # Lambda for Minimum Validation Error:
  Opt_Min_Lambda = Lambda_Grid[which.min(Mean_Val_Errors)]
  
  
  # Lambda for Min+1SE Validation Error
  
  Min_Plus1SE_Check = Mean_Val_Errors < Mean_Val_Errors[which.min(Mean_Val_Errors)] + SE_Val_Errors[which.min(Mean_Val_Errors)]
  
  #Find final index that falls below the threshold
  Fus_Loc = which(Min_Plus1SE_Check)
  
  #Have to do it backwards for TBA, because TBA fits from largest to smallest lambda 
  # (rarefit DOCUMENTATION DOESN'T SAY THIS AND THE LAMBDA INPUT ORDER IS IGNORED,
  # BUT LAMBDA GRID OUTPUT IS RETURNED AS THE SAME AS ORGINAL INPUT!!!! 
  # CAN LEAD TO CONFUSION SINCE COEFFICIENT SOLUTION CAN BE BACKWARDS WITHOUT INFORMING USER/CORRECTING)
  if(!Method=="TBA"){
    Opt_1SE_Lambda = Lambda_Grid[Fus_Loc[length(Fus_Loc)]]
  }else{
    Opt_1SE_Lambda = Lambda_Grid[Fus_Loc[1]]
  }
  
  #Optimal Lambdas for each lambda selection method
  Opt_Lambdas = c(Opt_Min_Lambda, Opt_1SE_Lambda)
  names(Opt_Lambdas) = c("Opt Min Lambda", "Opt 1SE Lambda")
  
  
  
  if(Method=="HFG LASSO"){
    
    Final_Model_Coeffs = HFG_LASSO_Fit(HierDesignMatrix, Response,
                                       Categ_Inds, Categ_Pair_Inds, S_Matrices_Weights,
                                       Opt_Lambdas, Proximal_Grad_Tol, Gamma)
    
  }else if(Method=="TBA"){
    
    TBA_Final_Model = rarefit(y=Response, X=HierDesignMatrix, A = A, intercept = FALSE, lambda = rev(Opt_Lambdas), alpha = 1)
    Final_Model_Coeffs = TBA_Final_Model$beta[[1]]
    Final_Model_Coeffs = Final_Model_Coeffs[,ncol(Final_Model_Coeffs):1]
    
  }
  
  #Return the method coefficients, validation errors during cross-validation, and the optimal lambdas selected
  return(list(Final_Model_Coeffs, Test_Set_Validation_Errors_Set, Opt_Lambdas))
  
  
}


