# HFGLASSO

Code developed to fit the method Hierarchical Fused Group (HFG) LASSO, as seen in the paper "HFG LASSO: a Method for Models with Hierarchical Categorical Predictors"

## Repository Structure

The repository is organized as follows:

- **Data**: Folder to place the US Colleges dataset (see Cerda and Varoquaux (2020)).
- **Scalability_Simulations**: Folder containing the scalability simulation code (and results subfolders) for varying sample size and varying HCP cardinality. See Section 5.3 of the paper.
- **Simulations**: Folder containing the primary simulation code (and results subfolders). See Sections 5.1, 5.2 and Appendix C (for the alternative weighting scheme simulation setting) of the paper.
- **Sparse_Data_Simulation**: Folder containing the simulation code (and results subfolders) for the sparse HCP category simulation. See Appendix D of the paper.
- **US_Colleges_Application**: Folder containing the code to conduct the data application with the US Colleges dataset. See Section 6 of the paper.
- **HFG_LASSO_func.R**: The R script containing the HFG LASSO function code. One function is for fitting HFG LASSO without additional non-HCP covariates, another function is for fitting HFG LASSO with additional non-HCP covariates and the last function allows computation of lambda_min (see Section 4.6 of the paper).
- **HFGL_Fitting_Algorithm_1_Step_3_Functions.cpp**: The C++ script containing function code to conduct Step 3 of Algorithm 1 to fit HFG LASSO (these functions are called by the functions in HFG_LASSO_func.R). We recommend installing the Rcpp R package to be able to compile this code, as it is compiled directly within the simulation and application code.

### Simulation Code

For each of the simulation code folders, the general structure is as follows:

- **Results**: Folder to store the simulation results.
- **generate_model_results.R**: The main driver R script, to run the simulations and store the results.
- **HCP_creation.R**: The R script that builds the HCP, including its design matrices, tree structure, active node selection, S_l matrices, ancestor-to-leaf path list, etc.
- **simulate_data.R**: The R script that builds the other parts of the data, primarily the response data but also the underlying true coefficients of the model.
- **simulation_study.R**: The R script to actually conduct the analyses and apply each of the methods of interest to the simulated data. Stores results in the Results folder.
- **simulation_results.R**: The R script used to compute the peformance measures based on the results and summarize the results in general.
- **simulation_table_creation.R** OR **simulation_figures.R**: R scripts used to create the results tables or figures as presented in the paper. 

### Application Code

For the application code folder, the general structure is as follows:

- **Results**: Folder to store the simulation results.
- **HFGL_US_Colleges_Application_DataCleaning.R**: The R script used to clean the US Colleges dataset and prepare other elements for applying HFG LASSO or TBA (e.g., building the S_l matrices, creating the ancestor-to-leaf node path list, etc.)
- **HFGL_US_Colleges_Application_CV.R**: The R script used to conduct the analyses on the US Colleges dataset.
- **US_Colleges_func.R**: The R script containing some function code to conduct the analyses on the US Colleges dataset more efficiently.

# General Notes

We recommend reviewing the R packages called upon at the start of each R script and installing those corresponding R packages for use.


# References

- Cerda, P., and Varoquaux, G. (2020). Encoding High-Cardinality String Categorical Variables. *IEEE Transactions on Knowledge and Data Engineering*, *34*(3), 1164-1176.


