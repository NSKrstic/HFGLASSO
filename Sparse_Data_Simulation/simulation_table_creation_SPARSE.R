#################################################
# Script: simulation_table_creation_SPARSE.R
# Author: Nikolas Krstic
# Purpose: Create the simulation table results from the manuscript
#################################################

#Set the scenario to create the table for
Scenario = "C"
Weight_Scheme = "Layer_Based"
# Location to Find Results
Model_Results_Loc = paste("Sparse_Data_Simulation/Results/Scenario_", Scenario, "/", sep="")


#######################################################

#Load SNR=0.25 results
load(file=paste(Model_Results_Loc, "BISM_Matrix_SNR0.25_", Weight_Scheme, ".Rdata", sep=""))

# Conduct average of the measures across the 100 repetitions for each of the 4 methods and 4 measures, and incorporate the standard errors for each of the averages
RESULTS_TABLE_P1 = as.data.frame(matrix(paste(round(colMeans(BISM_Matrix), 4), " (", round(apply(BISM_Matrix, 2, sd)/sqrt(100), 4) , ")", sep=""), ncol=4, nrow=4, byrow=TRUE))

#Rearrange the rows
RESULTS_TABLE_P1[c(1:4),] = RESULTS_TABLE_P1[c(1,4,2:3),]

#Rearrange the columns
RESULTS_TABLE_P1[,c(1:4)] = RESULTS_TABLE_P1[,c(1,3,2,4)]

rownames(RESULTS_TABLE_P1) = c("ECCS", "CGFI", "EE", "RMSE")
colnames(RESULTS_TABLE_P1) = c("HFG LASSO (Min)", "TBA (Min)", "HFG LASSO (1 SE)", "TBA (1 SE)")

#######################################################

#Load SNR=1 results
load(file=paste(Model_Results_Loc, "BISM_Matrix_SNR1_", Weight_Scheme, ".Rdata", sep=""))

# Conduct average of the measures across the 100 repetitions for each of the 4 methods and 4 measures, and incorporate the standard errors for each of the averages
RESULTS_TABLE_P2 = as.data.frame(matrix(paste(round(colMeans(BISM_Matrix), 4), " (", round(apply(BISM_Matrix, 2, sd)/sqrt(100), 4) , ")", sep=""), ncol=4, nrow=4, byrow=TRUE))

#Rearrange the rows
RESULTS_TABLE_P2[c(1:4),] = RESULTS_TABLE_P2[c(1,4,2:3),]

#Rearrange the columns
RESULTS_TABLE_P2[,c(1:4)] = RESULTS_TABLE_P2[,c(1,3,2,4)]

rownames(RESULTS_TABLE_P2) = c("ECCS", "CGFI", "EE", "RMSE")
colnames(RESULTS_TABLE_P2) = c("HFG LASSO (Min)", "TBA (Min)", "HFG LASSO (1 SE)", "TBA (1 SE)")

###############################

#Load SNR=4 results
load(file=paste(Model_Results_Loc, "BISM_Matrix_SNR4_", Weight_Scheme, ".Rdata", sep=""))

# Conduct average of the measures across the 100 repetitions for each of the 4 methods and 4 measures, and incorporate the standard errors for each of the averages
RESULTS_TABLE_P3 = as.data.frame(matrix(paste(round(colMeans(BISM_Matrix), 4), " (", round(apply(BISM_Matrix, 2, sd)/sqrt(100), 4) , ")", sep=""), ncol=4, nrow=4, byrow=TRUE))

#Rearrange the rows
RESULTS_TABLE_P3[c(1:4),] = RESULTS_TABLE_P3[c(1,4,2:3),]

#Rearrange the columns
RESULTS_TABLE_P3[,c(1:4)] = RESULTS_TABLE_P3[,c(1,3,2,4)]

rownames(RESULTS_TABLE_P3) = c("ECCS", "CGFI", "EE", "RMSE")
colnames(RESULTS_TABLE_P3) = c("HFG LASSO (Min)", "TBA (Min)", "HFG LASSO (1 SE)", "TBA (1 SE)")

############################################



