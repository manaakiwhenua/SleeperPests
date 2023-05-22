rm(list = ls())
# set working directory 
setwd("R:\\Projects/SL2116_MPI_SleeperPests/Analysis/TRSM/Release Code and Data")


###Test May 2023 version with code for propagule production, dispersal and population growth 
###timesteps passed as separate function
DoClimateChange=F
DoHumanSpread = T
source("TRSM_INA_OngoingTest.r")

DoClimateChange=T
DoHumanSpread = T
source("TRSM_INA_OngoingTest.r")


###Examine effect of ongoing threat on management effectiveness
DoClimateChange=F
DoHumanSpread = T
source("TRSM_INA_MultDetectionProbsInfoSpreadOngoingExternal.r")
source("TRSM_INA_MultDetectionProbsInfoSpreadOngoingExternalHeatMaps.r")

DoClimateChange=T
DoHumanSpread = T
source("TRSM_INA_MultDetectionProbsInfoSpreadOngoingExternal.r")
source("TRSM_INA_MultDetectionProbsInfoSpreadOngoingExternalHeatMaps.r")


