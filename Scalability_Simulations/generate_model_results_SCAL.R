#################################################
# Script: generate_model_results_SCAL.R
# Author: Nikolas Krstic
# Purpose: Overall driver script to vary the sample sizes or HCP cardinality
# and generate model results for multiple different seeds
# RUN THIS TO CONDUCT SIMULATION AND STORE RESULTS, THEN RUN "simulation_results_SCAL.R" TO RECOVER RESULTS
#################################################

source("HFG_LASSO_func.R")

Seeds = seq(10, 200, by=10)


Modelling_Script = "Scalability_Simulations/simulation_study_SCAL.R"

### Specify SNR Level (one of 0.25, 1 or 4)
SNR = 1

### Specify Weighting Scheme (one of "Layer_Based", "LB_SQRT", "Col_Dim" or "CD_SQRT")
Weight_Scheme = "Layer_Based"

### Specify Whether Varying Sample Size ("SS") or Varying HCP Cardinality ("HCPC")?
Sim_Scal_Type = "HCPC"

### Location to Store Simulation Results
Model_Results_Loc = paste("Scalability_Simulations/Results/", Sim_Scal_Type, "_Results/", sep="")


###Sample sizes of dataset (ITERATE ACROSS EACH OF THESE SAMPLE SIZES TO ASSESS PERFORMANCE/COMPUTATION TIMES)
if(Sim_Scal_Type == "SS"){
  Samp_Sizes = seq(1000, 10000, by=1000)
}else if(Sim_Scal_Type == "HCPC"){
  Samp_Sizes = 10000
}

###Sample size of validation dataset
Valid_Samp_Size = 5000

###Sample size of test dataset
Test_Samp_Size = 5000

###Cardinality of the HCP
if(Sim_Scal_Type == "SS"){
  Cardinalities = 20
}else if(Sim_Scal_Type == "HCPC"){
  Cardinalities = seq(50, 500, 50)
}



#################################################################
# Conduct Simulation over the Seeds, Sample Sizes and Cardinalities specified above

for(a in 1:length(Seeds)){
  
  Seed = Seeds[a]
  
  for(b in 1:length(Samp_Sizes)){
    
    for(c in 1:length(Cardinalities)){
      
      Samp_Size = Samp_Sizes[b]
      Cardinality = Cardinalities[c]
      
      source(Modelling_Script)
      
    }
    
  }
  
}




