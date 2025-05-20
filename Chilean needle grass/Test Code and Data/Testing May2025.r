###########################################################################
###########################################################################
###Define default parameters
###########################################################################
###########################################################################
rm(list = ls())
###Set root directory where testing data and scripts stored
main.dir = r"[C:\Users\MasonN\OneDrive - MWLR\Documents\GitHub\SleeperPests\Chilean needle grass\Test Code and Data]"
dir.create(main.dir,showWarnings = F)
setwd(main.dir)
ClimexVars = c("EI","EI_Niwa_204","EI_Niwa209")
ClimateScenarios = c("Current","Future_2040","Future_2090")
EI_Prob_CurveType = "Logit"
#EI_Prob_CurveType = "SplitLinear"

###Install required packages
#install.packages("actuar")
#install.packages("magick")
#install.packages("sf")

library(tidyverse)
library(dplyr)
library(boot)
library(abind)
library(igraph) 

###########################################################################
###########################################################################
###Read in INApest functions from source files
###########################################################################
###########################################################################
CoreFunctionDir <- r"[C:\Users\MasonN\OneDrive - MWLR\Documents\GitHub\SleeperPests\INApest INApestMeta core function code\]"
source(paste0(CoreFunctionDir,"INApestParallel.R"))
source(paste0(CoreFunctionDir,"INApestParallelINAScene.R"))
source(paste0(CoreFunctionDir,"INApestINAScene.R"))
source(paste0(CoreFunctionDir,"INApest.R"))
source("Scripts/INApestHeatMaps.R")

#############################################################
###Run Marlborough CNG simulations
#############################################################

####test different version of INApest give similar results within different climate scenarios
for(cs in 1:3)
{
###Prepare data for relevant climate scenario (cs)
source("Scripts/dataprep.R")
###Run test source code testing for comparable results
###between INApest function versions
source("Scripts/INApestParallelTestingMay2025.R")
}


###Run test source code testing external info functionality
###for each version of INApest function
cs = 1
source("Scripts/dataprep.R")
Function = "INApest"
source("Scripts/ExternalInfo TestingMay2025.R")
Function = "INApestParallel"
source("Scripts/ExternalInfo TestingMay2025.R")
Function = "INApestINAscene"
source("Scripts/ExternalInfo TestingMay2025.R")
Function = "INApestParallelINAscene"
source("Scripts/ExternalInfo TestingMay2025.R")

###Run test source code testing external invasion functionality
###for each version of INApest function with INAscene functionality done manually
Function = "INApest"
source("Scripts/ExternalInvasion TestingMay2025.R")
Function = "INApestParallel"
source("Scripts/ExternalInvasion TestingMay2025.R")

###Test effect of long distance dispersal
source("Scripts/NoLDDParallel May2025.R")

###Test effect of communication in social networks
source("Scripts/SocialNetworks TestingMay2025.R")

###Test functionality varying management and establishment over time
source("Scripts/ManagementParamInput_TestingMay2025.R")

