#################################################
# Script: simulation_figures_SCAL.R
# Author: Nikolas Krstic
# Purpose: Script to produce the figures from the scalability simulations
#################################################

library(tidyverse)
library(cowplot)

Sim_Scal_Type = "HCPC"

if(Sim_Scal_Type=="SS"){
  Series = seq(1000, 10000, by=1000)
  ### Location to Find Results
  Model_Results_Loc = paste("Scalability_Simulations/Results/", Sim_Scal_Type, "_Results/", sep="")
}else if(Sim_Scal_Type=="HCPC"){
  Series = seq(50, 500, 50)
  ### Location to Find Results
  Model_Results_Loc = paste("Scalability_Simulations/Results/", Sim_Scal_Type, "_Results/", sep="")
}

#Begin Loading the Simulation Results

Overall_Simulation_Results = list()

for(i in 1:length(Series)){
  
  if(Sim_Scal_Type=="SS"){
    load(file=paste(Model_Results_Loc, "BISM_Matrix_Samp_Size_", Series[i], "_SNR_", SNR, "_", Weight_Scheme, ".Rdata", sep=""))
    Overall_Simulation_Results[[i]] = cbind("Sample_Size" = Series[i], BISM_Matrix)
  }else if(Sim_Scal_Type=="HCPC"){
    load(file=paste(Model_Results_Loc, "BISM_Matrix_Cardinality_", Series[i], "SNR", SNR, "_", Weight_Scheme, ".Rdata", sep=""))
    Overall_Simulation_Results[[i]] = cbind("Cardinality" = Series[i], BISM_Matrix)
  }
}

Overall_Simulation_Results_DF = as.data.frame(do.call(rbind, Overall_Simulation_Results))

Long_Simulation_Results_DF = pivot_longer(Overall_Simulation_Results_DF, cols=-1, names_to="Measure", values_to="Measure Score")

Long_Simulation_Results_DF$`Measure Type` = rep(c(rep("ECCS", 2), rep("EE",2 ), rep("RMSE", 2), rep("CGFI", 2), rep("Run Time (s)", 2)), times=nrow(Long_Simulation_Results_DF)/10)
Long_Simulation_Results_DF$`Measure` = rep(c("HFG LASSO", "TBA"), times=nrow(Long_Simulation_Results_DF)/2)

#### Plot the Simulation Figure

if(Sim_Scal_Type=="SS"){
  
  pdf(paste(Model_Results_Loc, "Samp_Size_Scalability_Sim_Results.pdf", sep=""))
  
  Long_Simulation_Results_DF_Perf = Long_Simulation_Results_DF[Long_Simulation_Results_DF$`Measure Type`!="Run Time (s)",]
  Long_Simulation_Results_DF_SS_CompTimes = Long_Simulation_Results_DF[Long_Simulation_Results_DF$`Measure Type`=="Run Time (s)",]
  
  ### Take the means and standard errors:
  Long_Simulation_Results_DF_FINAL = Long_Simulation_Results_DF_Perf %>% group_by(Sample_Size, `Measure Type`, Measure) %>% summarize(`Mean Measure`= mean(`Measure Score`), `SE Measure`=sd(`Measure Score`)/n())
  names(Long_Simulation_Results_DF_FINAL)[1] = "Sample Size"
  
  Sim_Plot = ggplot(data=Long_Simulation_Results_DF_FINAL, aes(x=`Sample Size`, y=`Mean Measure`, colour=Measure))+
    geom_line()+
    geom_point()+
    geom_errorbar(aes(ymin=`Mean Measure`-`SE Measure`, ymax=`Mean Measure`+`SE Measure`), width=100)+
    facet_wrap(~`Measure Type`, nrow=2, ncol=2, scales="free")+
    theme_bw()
  
  print(Sim_Plot)
  
  dev.off()
  
  
  
}else if(Sim_Scal_Type=="HCPC"){
  
  pdf(paste(Model_Results_Loc, "Cardinality_Scalability_Sim_Results.pdf", sep=""))
  
  Long_Simulation_Results_DF_Perf = Long_Simulation_Results_DF[Long_Simulation_Results_DF$`Measure Type`!="Run Time (s)",]
  Long_Simulation_Results_DF_HCPC_CompTimes = Long_Simulation_Results_DF[Long_Simulation_Results_DF$`Measure Type`=="Run Time (s)",]
  
  ### Take the means and standard errors:
  Long_Simulation_Results_DF_FINAL = Long_Simulation_Results_DF_Perf %>% group_by(Cardinality, `Measure Type`, Measure) %>% summarize(`Mean Measure`= mean(`Measure Score`), `SE Measure`=sd(`Measure Score`)/n())
  
  Sim_Plot = ggplot(data=Long_Simulation_Results_DF_FINAL, aes(x=Cardinality, y=`Mean Measure`, colour=Measure))+
    geom_line()+
    geom_point()+
    geom_errorbar(aes(ymin=`Mean Measure`-`SE Measure`, ymax=`Mean Measure`+`SE Measure`), width=5)+
    facet_wrap(~`Measure Type`, nrow=3, ncol=2, scales="free")+
    theme_bw()
  
  print(Sim_Plot)
    
  dev.off()
  
  
  
}



############################################################################
############################################################################
############################################################################
# After running both "SS" and "HCPC" versions of the code above, to produce "Long_Simulation_Results_DF_SS_CompTimes" and "Long_Simulation_Results_DF_HCPC_CompTimes", generate the computation time figure

Long_Simulation_Results_DF_SS_CompTimes$`Measure Type` = "Varying Sample Size"
Long_Simulation_Results_DF_HCPC_CompTimes$`Measure Type` = "Varying Cardinality"

### Take the means and standard errors (SS):
Long_Simulation_Results_DF_FINAL = Long_Simulation_Results_DF_SS_CompTimes %>% group_by(Sample_Size, `Measure Type`, Measure) %>% summarize(`Mean Measure`= mean(`Measure Score`), `SE Measure`=sd(`Measure Score`)/n())
names(Long_Simulation_Results_DF_FINAL)[1] = "Sample Size"

Sim_Plot_1 = ggplot(data=Long_Simulation_Results_DF_FINAL, aes(x=`Sample Size`, y=`Mean Measure`, colour=Measure))+
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin=`Mean Measure`-`SE Measure`, ymax=`Mean Measure`+`SE Measure`), width=100)+
  theme_bw()+
  labs(y="Computation Time (s)")+
  theme(legend.position="none")

Legend = get_legend(ggplot(data=Long_Simulation_Results_DF_FINAL, aes(x=`Sample Size`, y=`Mean Measure`, colour=Measure))+
                      geom_line()+
                      geom_point()+
                      geom_errorbar(aes(ymin=`Mean Measure`-`SE Measure`, ymax=`Mean Measure`+`SE Measure`), width=100)+
                      theme_bw()+
                      labs(y="Computation Time (s)"))

### Take the means and standard errors:
Long_Simulation_Results_DF_FINAL = Long_Simulation_Results_DF_HCPC_CompTimes %>% group_by(Cardinality, `Measure Type`, Measure) %>% summarize(`Mean Measure`= mean(`Measure Score`), `SE Measure`=sd(`Measure Score`)/n())

Sim_Plot_2 = ggplot(data=Long_Simulation_Results_DF_FINAL, aes(x=Cardinality, y=`Mean Measure`, colour=Measure))+
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin=`Mean Measure`-`SE Measure`, ymax=`Mean Measure`+`SE Measure`), width=5)+
  theme_bw()+
  labs(y="")+
  theme(legend.position="none")

Figure_6 = plot_grid(Sim_Plot_1, Sim_Plot_2, ncol=2, align="h")

Figure_6 = plot_grid(Figure_6, Legend, rel_widths=c(3,0.6))

pdf("Scalability_Simulations/Results/Average_Computation_Times.pdf", height=4, width=8)

print(Figure_6)

dev.off()




