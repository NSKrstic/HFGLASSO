#################################################
# Script: generate_model_results.R
# Author: Nikolas Krstic
# Purpose: Overall driver script to choose Scenario and generate model results for multiple different seeds
# RUN THIS TO CONDUCT SIMULATION AND STORE RESULTS, THEN RUN "simulation_results.R" TO RECOVER RESULTS
#################################################

source("HFG_LASSO_func.R")

Seeds = seq(100, 1090, by=10)


Modelling_Script = "Simulations/simulation_study.R"

### Specify whether doing Scenario A, B, or C
Scenario = "C"

### Specify SNR Level (one of 0.25, 1 or 4)
SNR = 4

### Specify Weighting Scheme (one of "Layer_Based", "LB_SQRT", "Col_Dim" or "CD_SQRT")
Weight_Scheme = "Layer_Based"

### Location to Store Simulation Results
Model_Results_Loc = paste("Simulations/Results/Scenario_", Scenario, "/", sep="")


#################################################################
# Conduct Simulation over Seeds specified above

for(i in 1:length(Seeds)){
  
  Seed = Seeds[i]
  
  source(Modelling_Script)
  
}




