#################################################
# Script: generate_model_results_SPARSE.R
# Author: Nikolas Krstic
# Purpose: Overall driver script to choose Scenario and generate model results for multiple different seeds.
# This to examine how HFG LASSO performs when some of the categories are sparse (i.e. have few observations)
# RUN THIS TO CONDUCT SIMULATION AND STORE RESULTS, THEN RUN "simulation_results.R" TO RECOVER RESULTS
#################################################

source("HFG_LASSO_func.R")

Seeds = seq(100, 1090, by=10)


Modelling_Script = "Sparse_Data_Simulation/simulation_study_SPARSE.R"

### Specify whether doing Scenario A, B, or C
Scenario = "A"

### Specify SNR Level (one of 0.25, 1 or 4)
SNR = 4

### Specify Weighting Scheme (one of "Layer_Based", "LB_SQRT", "Col_Dim" or "CD_SQRT")
Weight_Scheme = "Layer_Based"

### Location to Store Simulation Results
Model_Results_Loc = paste("Sparse_Data_Simulation/Results/Scenario_", Scenario, "/", sep="")


#################################################################
# Conduct Simulation over Seeds specified above

for(i in 1:length(Seeds)){
  
  Seed = Seeds[i]
  
  source(Modelling_Script)
  
}




