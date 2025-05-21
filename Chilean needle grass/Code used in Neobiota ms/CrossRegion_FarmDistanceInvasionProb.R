##################################################
###Calculate cross-region invasion probability for individual farms
###Permits modelling scenario of monitoring and management
###for cross-region incursions
###Use same Pareto kernel as used to generate adjacency matrices 
###between farms in the same region
##################################################


library(INA)
library(tidyverse)
library(dplyr)
library(sp)
library(actuar)
memory.limit(size=600000) 

ResultsDir= paste0(main.dir,paste0("/",ClimateScenarios[cs],"/Inputs/CrossRegionDistance/"))
dir.create(ResultsDir,recursive = T)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
for(region1 in 1:(length(Regions)-1))
{
for(region2 in (region1+1):length(Regions))
{
#----load distances-------------------------------------------------------------
dist_mat = readRDS(paste0(ResultsDir,Regions[region1],"_",Regions[region2],"_dist_mat_farm.rds"))
k= 7 #shape parameter
b= 500 #scale parameter
thresh= 1/20 #1 in 20 year spread risk between adjacent farms
dist_mat[dist_mat==0]<-0.01
farm2farm_probs <- thresh*actuar::ppareto(dist_mat, shape=k,  scale=b, lower.tail=F,log.p = F)
Region1InvasionRisk = colSums(farm2farm_probs)
Region2InvasionRisk = rowSums(farm2farm_probs)
saveRDS(farm2farm_probs, paste0(ResultsDir,Regions[region1],"_",Regions[region2],"_InvasionProb.rds"))
saveRDS(t(farm2farm_probs), paste0(ResultsDir,Regions[region2],"_",Regions[region1],"_InvasionProb.rds"))
saveRDS(Region1InvasionRisk, paste0(ResultsDir,Regions[region1],"InvasionRiskFrom_",Regions[region2],".rds"))
saveRDS(Region2InvasionRisk, paste0(ResultsDir,Regions[region2],"InvasionRiskFrom_",Regions[region1],".rds"))
}
}

