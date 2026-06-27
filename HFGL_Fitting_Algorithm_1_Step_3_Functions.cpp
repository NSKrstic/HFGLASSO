#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace Rcpp;


//Function to subset a vector (in Rcpp, which is important because indexing is from 0 instead of 1)

// [[Rcpp::export]]
Rcpp::NumericVector subset_vector(NumericVector x, IntegerVector y){
  
  //Code taken from https://stackoverflow.com/questions/17743746/how-to-extract-multiple-values-from-a-numericvector-in-c
  
  // Length of the index vector
  int n = y.size();
  // Initialize output vector
  NumericVector out(n);
  // Loop through index vector and extract values of x at the given positions
  for(int i = 0; i < n; i++){
    out[i] = x[y[i]-1];
  }
  // Return output
  return out;
}


///////////////////////////////////////////////////////SDMM Functions


// [[Rcpp::export]]
Rcpp::List Thetas_SDMM_Computation(double Gamma, NumericVector Gamma_Lambda_S_Weights, NumericVector Thetas, Rcpp::List Categ_Inds,
                                      Rcpp::List Char_Categ_Inds, IntegerVector DM1s_RowByCol_Inds,
                                      Rcpp::List Z_Vecs, Rcpp::List Y_Vecs, NumericVector Gamma_Response,
                                      Rcpp::List Categ_Pair_Inds, Rcpp::List DM1s_Inds, Rcpp::List Ones_RowInds, Rcpp::List NegOnes_RowInds,
                                      NumericMatrix Q_Inv_Mat, double Proximal_Grad_Tol, int NCol_HDM, int Length_Categ_Inds){
  
  double Outer_Tolerance_Check = 1000;
  
  NumericVector New_Thetas = clone(Thetas);
  
  Rcpp::List New_Y_Vecs = clone(Y_Vecs);
  Rcpp::List New_Z_Vecs = clone(Z_Vecs);
  
  for(int iteration=0; (Outer_Tolerance_Check > Proximal_Grad_Tol) && (iteration<100000); iteration++){
    
    NumericVector Old_Thetas = clone(New_Thetas);
    
    NumericVector New_YZ_Diff_Sum (NCol_HDM);
    
    ////////First Iteration Code
    // Y Vec and Z Vec Code
    
    NumericVector Temp_Z_Vec = New_Z_Vecs[0];
    
    NumericVector Input = subset_vector(New_Thetas, DM1s_RowByCol_Inds) + Temp_Z_Vec;
    
    NumericVector Temp_New_Y_Vec = (Input+(Gamma_Response))/(1+Gamma);
    New_Y_Vecs[0] = Temp_New_Y_Vec;
    
    NumericVector Temp_New_Z_Vec = Input - Temp_New_Y_Vec;
    New_Z_Vecs[0] = Temp_New_Z_Vec;
    
    // YZ Diff Code
    
    NumericVector TYZ_Diff = Temp_New_Y_Vec-Temp_New_Z_Vec;
    
    NumericVector Temp_New_YZ_Diff (NCol_HDM, 0.0);
    
    for(int q=0; q<(NCol_HDM); ++q){
      
      IntegerVector Temp_DM1s_Ind = DM1s_Inds[q];
      
      Temp_New_YZ_Diff[q] = sum(subset_vector(TYZ_Diff, Temp_DM1s_Ind));
      
    }
    
    New_YZ_Diff_Sum = New_YZ_Diff_Sum + Temp_New_YZ_Diff;
    
    
    ////////Further Iteration Code
    
    for(int j=1; j<(Length_Categ_Inds+1); ++j){
      
      // Y Vec and Z Vec Code

      NumericVector Temp_Z_Vec = New_Z_Vecs[j];
      
      IntegerMatrix Categ_Pair_Ind_Mat = Categ_Pair_Inds[j-1];
      IntegerVector FirstCoords = Categ_Pair_Ind_Mat( _ , 0 );
      IntegerVector SecondCoords = Categ_Pair_Ind_Mat( _ , 1 );

      NumericVector Input = (subset_vector(New_Thetas, FirstCoords) - subset_vector(New_Thetas, SecondCoords)) + Temp_Z_Vec;

      NumericVector Temp_New_Y_Vec = max(NumericVector::create(0, (1-(Gamma_Lambda_S_Weights[j-1] / sqrt(sum(pow(Input, 2)))))))*Input;
      New_Y_Vecs[j] = Temp_New_Y_Vec;
        
      NumericVector Temp_New_Z_Vec = Input - Temp_New_Y_Vec;
      New_Z_Vecs[j] = Temp_New_Z_Vec;
      
      
      // YZ Diff Code
      
      NumericVector TYZ_Diff = Temp_New_Y_Vec-Temp_New_Z_Vec;
      
      NumericVector Temp_New_YZ_Diff (NCol_HDM, 0.0);
      
      NumericVector C_Inds = Categ_Inds[j-1];
      CharacterVector Char_C_Inds = Char_Categ_Inds[j-1];
        
      NumericVector Temp (Char_C_Inds.size(), 0.0);
      
      List ORI = Ones_RowInds[j-1];
      List NORI = NegOnes_RowInds[j-1];
        
      for(int q=0; q<Char_C_Inds.size(); q++){
          
        String ChI = Char_C_Inds[q];
        
        IntegerVector Ones_Vec = ORI[ChI];
        IntegerVector NegOnes_Vec = NORI[ChI];
          
        Temp[q] = sum(subset_vector(TYZ_Diff, Ones_Vec))-sum(subset_vector(TYZ_Diff, NegOnes_Vec));
          
      }
        
      Temp_New_YZ_Diff[C_Inds] = Temp;
      
      New_YZ_Diff_Sum = New_YZ_Diff_Sum + Temp_New_YZ_Diff;

    }
    
    
    //New Thetas Code
    
    for(int m=0; m < Q_Inv_Mat.nrow(); ++m){

      New_Thetas[m] = sum(Q_Inv_Mat( m , _ ) * New_YZ_Diff_Sum);
    }
    
    Outer_Tolerance_Check = sqrt(sum(pow(Old_Thetas-New_Thetas, 2)));
    
    
    if((iteration+1) % 1000 == 0){
      Rprintf("  Iteration: %i \n", (iteration+1));
      Rprintf("  Tolerance Check: %.3e \n", Outer_Tolerance_Check);
    }
        
  }
  
  Rcpp::List Final_Results;
  Final_Results["Thetas"] = New_Thetas;
  Final_Results["Y_Vecs"] = New_Y_Vecs;
  Final_Results["Z_Vecs"] = New_Z_Vecs;
  
  return Final_Results;
  
  
}






// [[Rcpp::export]]
Rcpp::List Thetas_SDMM_Computation_NonHCPVer(double Gamma, NumericVector Gamma_Lambda_S_Weights, NumericVector Thetas, Rcpp::List Categ_Inds,
                                   Rcpp::List Char_Categ_Inds, IntegerVector DM1s_RowByCol_Inds,
                                   Rcpp::List Z_Vecs, Rcpp::List Y_Vecs, NumericVector Gamma_Response,
                                   Rcpp::List Categ_Pair_Inds, Rcpp::List DM1s_Inds, Rcpp::List Ones_RowInds, Rcpp::List NegOnes_RowInds,
                                   NumericMatrix Q_Inv_Mat, double Proximal_Grad_Tol, int NCol_HDM, int Length_Categ_Inds, int HDM_NHCP_Length,
                                   IntegerVector NonHCP_Inds, NumericMatrix HCP_Des_Mat, NumericMatrix NHCP_Des_Mat){
  
  double Outer_Tolerance_Check = 1000;
  
  NumericVector New_Thetas = clone(Thetas);
  
  Rcpp::List New_Y_Vecs = clone(Y_Vecs);
  Rcpp::List New_Z_Vecs = clone(Z_Vecs);
  
  for(int iteration=0; (Outer_Tolerance_Check > Proximal_Grad_Tol) && (iteration<10000); iteration++){
    
    NumericVector Old_Thetas = clone(New_Thetas);
    
    NumericVector New_YZ_Diff_Sum (NCol_HDM);
    
    ////////First Iteration Code
    // Y Vec and Z Vec Code
    
    NumericVector Temp_Z_Vec = New_Z_Vecs[0];
    
    // NONHCP Component
    NumericVector NonHCP_Thetas = subset_vector(New_Thetas, NonHCP_Inds);
    NumericVector NHCP_S_Vec (NHCP_Des_Mat.nrow());
    
    for(int k=0; k < NHCP_Des_Mat.nrow(); ++k){
      
      NHCP_S_Vec[k] = sum(NHCP_Des_Mat( k , _ ) * NonHCP_Thetas);
    }
    
    NumericVector Input = subset_vector(New_Thetas, DM1s_RowByCol_Inds) + NHCP_S_Vec + Temp_Z_Vec;
    
    NumericVector Temp_New_Y_Vec = (Input+(Gamma_Response))/(1+Gamma);
    New_Y_Vecs[0] = Temp_New_Y_Vec;
    
    NumericVector Temp_New_Z_Vec = Input - Temp_New_Y_Vec;
    New_Z_Vecs[0] = Temp_New_Z_Vec;
    
    // YZ Diff Code
    
    NumericVector TYZ_Diff = Temp_New_Y_Vec-Temp_New_Z_Vec;
    
    NumericVector Temp_New_YZ_Diff (NCol_HDM, 0.0);
    
    for(int q=0; q<(NCol_HDM); ++q){
      
      if(q<HCP_Des_Mat.ncol()){
      
        IntegerVector Temp_DM1s_Ind = DM1s_Inds[q];
        
        Temp_New_YZ_Diff[q] = sum(subset_vector(TYZ_Diff, Temp_DM1s_Ind));
        
      }else{
        // NONHCP Component
        Temp_New_YZ_Diff[q] = sum(NHCP_Des_Mat( _ , q-HCP_Des_Mat.ncol()) * TYZ_Diff);
        
      }
      
    }
    
    New_YZ_Diff_Sum = New_YZ_Diff_Sum + Temp_New_YZ_Diff;
    
    
    ////////Further Iteration Code
    
    for(int j=1; j<(Length_Categ_Inds+1); ++j){
      
      // Y Vec and Z Vec Code
      
      NumericVector Temp_Z_Vec = New_Z_Vecs[j];
      
      IntegerMatrix Categ_Pair_Ind_Mat = Categ_Pair_Inds[j-1];
      IntegerVector FirstCoords = Categ_Pair_Ind_Mat( _ , 0 );
      IntegerVector SecondCoords = Categ_Pair_Ind_Mat( _ , 1 );
      
      NumericVector Input = (subset_vector(New_Thetas, FirstCoords) - subset_vector(New_Thetas, SecondCoords)) + Temp_Z_Vec;
      
      NumericVector Temp_New_Y_Vec = max(NumericVector::create(0, (1-(Gamma_Lambda_S_Weights[j-1] / sqrt(sum(pow(Input, 2)))))))*Input;
      New_Y_Vecs[j] = Temp_New_Y_Vec;
      
      NumericVector Temp_New_Z_Vec = Input - Temp_New_Y_Vec;
      New_Z_Vecs[j] = Temp_New_Z_Vec;
      
      
      // YZ Diff Code
      
      NumericVector TYZ_Diff = Temp_New_Y_Vec-Temp_New_Z_Vec;
      
      NumericVector Temp_New_YZ_Diff (NCol_HDM, 0.0);
      
      NumericVector C_Inds = Categ_Inds[j-1];
      CharacterVector Char_C_Inds = Char_Categ_Inds[j-1];
      
      NumericVector Temp (Char_C_Inds.size(), 0.0);
      
      List ORI = Ones_RowInds[j-1];
      List NORI = NegOnes_RowInds[j-1];
      
      for(int q=0; q<Char_C_Inds.size(); q++){
        
        String ChI = Char_C_Inds[q];
        
        IntegerVector Ones_Vec = ORI[ChI];
        IntegerVector NegOnes_Vec = NORI[ChI];
        
        Temp[q] = sum(subset_vector(TYZ_Diff, Ones_Vec))-sum(subset_vector(TYZ_Diff, NegOnes_Vec));
        
      }
      
      Temp_New_YZ_Diff[C_Inds] = Temp;
      
      New_YZ_Diff_Sum = New_YZ_Diff_Sum + Temp_New_YZ_Diff;
      
    }
    
    
    // New Thetas Code
    
    for(int m=0; m < Q_Inv_Mat.nrow(); ++m){
      
      New_Thetas[m] = sum(Q_Inv_Mat( m , _ ) * New_YZ_Diff_Sum);
    }
    
    Outer_Tolerance_Check = sqrt(sum(pow(Old_Thetas-New_Thetas, 2)));
    
    
    if((iteration+1) % 1000 == 0){
      Rprintf("  Iteration: %i \n", (iteration+1));
      Rprintf("  Tolerance Check: %f \n", Outer_Tolerance_Check);
    }
    
  }
  
  Rcpp::List Final_Results;
  Final_Results["Thetas"] = New_Thetas;
  Final_Results["Y_Vecs"] = New_Y_Vecs;
  Final_Results["Z_Vecs"] = New_Z_Vecs;
  
  return Final_Results;
  
  
}

