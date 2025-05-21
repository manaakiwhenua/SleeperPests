##################################################
###Calculate long distance invasion probability between individual farms
###within and between regions
###Long distance dispersal assumed to be human-mediated
###Use data on cattle movements from NAIT to weigh LDD probability by distance
##################################################

library(tidyverse)
library(dplyr)
library(sp)
memory.limit(size=600000) 

ResultsDir= paste0(main.dir,paste0("/",ClimateScenarios[cs],"/Inputs/CrossRegionDistance/"))
dir.create(ResultsDir,recursive = T)

#######################################################################################
###Loop through all possible region combinations
###Including all region 1 = region2 combinations
###This is required for generating LDD weights linking each farm to all other at-risk farms in NZ
#######################################################################################

Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
for(region1 in 1:length(Regions))
{
for(region2 in (region1):length(Regions))
{
#----load data-------------------------------------------------------------
if(region1 != region2)
	dist_mat = readRDS(paste0(ResultsDir,Regions[region1],"_",Regions[region2],"_dist_mat_farm.rds"))
if(region1 == region2)
	dist_mat = readRDS(paste0(main.dir,"/",ClimateScenarios[cs],"/Inputs/",Regions[region1],"_dist_mat_farm.rds"))
######################################################
###Use NZ cattle movement data to build a
###Long distance dispersal matrix
###Sourced from https://www.mpi.govt.nz/dmsdocument/49114-Analysis-of-New-Zealand-Stock-Movement-Data
###Figure 13a
######################################################
######################################################
LongDistProbPerFarm = 0.05
CattleMovements = read.csv("NZcattleMovements.csv",as.is = T, header = T)
plot(CattleMovements$Distance_m, CattleMovements$Proportion.of.movements)
X = CattleMovements$Distance_m
Y = CattleMovements$Proportion.of.movements
loess30 <- loess(Y ~ X,span=.3)

smooth30 <- predict(loess30) 
#plot(X, Y, pch=19, main='Loess Regression Models')
#lines(smooth30, x=X, col='red')
LDDmatrix = dist_mat
for(i in 1:ncol(dist_mat))
      {
      LDDmatrix[,i] = predict(loess30,newdata = dist_mat[,i])
      LDDmatrix[,i] = LDDmatrix[,i]*LongDistProbPerFarm
      }
if(region1 == region2)
	diag(LDDmatrix) = 0
Region1InvasionRisk = colSums(LDDmatrix)
Region2InvasionRisk = rowSums(LDDmatrix)
saveRDS(LDDmatrix, paste0(ResultsDir,Regions[region1],"_",Regions[region2],"_InvasionProb_LDD.rds"))
if(region1 != region2)
	saveRDS(t(LDDmatrix), paste0(ResultsDir,Regions[region2],"_",Regions[region1],"_InvasionProb_LDD.rds"))
saveRDS(Region1InvasionRisk, paste0(ResultsDir,Regions[region1],"InvasionRiskFrom_",Regions[region2],"_LDD.rds"))
if(region1 != region2)
	saveRDS(Region2InvasionRisk, paste0(ResultsDir,Regions[region2],"InvasionRiskFrom_",Regions[region1],"_LDD.rds"))
}
}

########################################################
###End of loop
########################################################

