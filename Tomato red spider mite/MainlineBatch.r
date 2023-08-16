rm(list = ls())
# set working directory 
setwd(r"[C:\Users\MasonN\OneDrive - MWLR\Documents\GitHub\SleeperPests\Tomato red spider mite]")
install.packages("gnorm")

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

###Test August 2023 version with code for parallel processing
###and external information sources
DoClimateChange=F
DoHumanSpread = T
source("TRSM_INA_ExternalInfoParallel.r")

###Test August 2023 version with code for parallel processing
###and reduced environmenta establishment prob and local survival
###between timesteps
DoClimateChange=F
DoHumanSpread = T
source("TRSM_INA_EnvEstabSurvivalParallel.r")


