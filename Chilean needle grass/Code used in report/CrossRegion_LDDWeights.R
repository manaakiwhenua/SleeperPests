##################################################
###Weight long distance dispersal probabilities so that
###probabilities add to 1 for each source farm nationally
###Essentially shares dispersal events across potential sink farms
###Using long distance dispersal probabilies
##################################################

library(INA)
library(tidyverse)
library(dplyr)
library(sp)
library(actuar)
memory.limit(size=600000) 
options(digits = 12)
ResultsDir= paste0(main.dir,paste0("/",ClimateScenarios[cs],"/Inputs/CrossRegionDistance/"))
dir.create(ResultsDir,recursive = T)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
for(region1 in 1:length(Regions))
{
SelfLDD = readRDS(paste0(ResultsDir,Regions[region1],"_",Regions[region1],"_InvasionProb_LDD.rds"))
SourceFarmLDDsum = colSums(SelfLDD)
for(region2 in 1:length(Regions))
if(region1 != region2)
{
#----load data-------------------------------------------------------------
###Source region comes first in file name
LDDmatrix = readRDS(paste0(ResultsDir,Regions[region1],"_",Regions[region2],"_InvasionProb_LDD.rds"))
SourceFarmLDDsum = SourceFarmLDDsum + colSums(LDDmatrix)
}

for(region2 in 1:length(Regions))
{
LDDmatrix = readRDS(paste0(ResultsDir,Regions[region1],"_",Regions[region2],"_InvasionProb_LDD.rds"))
LDDweight <- sweep(LDDmatrix,2,SourceFarmLDDsum,`/`)
saveRDS(LDDweight,paste0(ResultsDir,Regions[region1],"_",Regions[region2],"_InvasionProb_LDDweight.rds"))
}
###Check weights add to 1 for each farm
#TestWeight = rep(0,length = length(SourceFarmLDDsum))
#for(region2 in 1:length(Regions))
#{
#LDDweight = readRDS(paste0(ResultsDir,Regions[region1],"_",Regions[region2],"_InvasionProb_LDDweight.rds"))
#TestWeight = TestWeight+colSums(LDDweight)
#}
#TestWeight
#mean(TestWeight)
}
