rm(list = ls())
# set working directory 
setwd(r"[C:\Users\MasonN\OneDrive - MWLR\Documents\GitHub\SleeperPests\Tomato red spider mite]")
#install.packages("gnorm")
CoreFunctionDir <- r"[C:\Users\MasonN\OneDrive - MWLR\Documents\GitHub\SleeperPests\INApest INApestMeta core function code\]"

source(paste0(CoreFunctionDir,"INApestMeta.r"))
source(paste0(CoreFunctionDir,"INApestMetaParallel.r"))
source(paste0(CoreFunctionDir,"INApestMetaMultipleLandUse.r"))
source(paste0(CoreFunctionDir,"INApestMetaParallelMultipleLandUse.r"))

###Test effect of ongoing external invasions (e.g. border incursions)
###Under Current climate
DoClimateChange=F
DoHumanSpread = T
source("TRSM_INA_OngoingTest May2025.r")

###Test effect of ongoing external invasions (e.g. border incursions)
###Under Future climate
DoClimateChange=T
DoHumanSpread = T
source("TRSM_INA_OngoingTest May2025.r")

###These tests use parallel versions of core functions to reduce processing time

###Test multiple land use function and effect of ongoing external communication
###Under Current climate
DoClimateChange=F
DoHumanSpread = T
source("MultipleLandUseExternalInfo May2025.r")

###Test multiple land use function and effect of ongoing external communication
###Under Future climate
DoClimateChange=T
DoHumanSpread = T
source("MultipleLandUseExternalInfo May2025.r")


###Examine effect of ongoing external invasion  on management effectiveness
DoClimateChange=F
DoHumanSpread = T
source("MultDetectionProbsInfoSpreadOngoingExternal_May2025.r")
#source("TRSM_INA_MultDetectionProbsInfoSpreadOngoingExternalHeatMaps.r")



DoClimateChange=T
DoHumanSpread = T
source("MultDetectionProbsInfoSpreadOngoingExternal_May2025.r")
#source("TRSM_INA_MultDetectionProbsInfoSpreadOngoingExternalHeatMaps.r")



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


