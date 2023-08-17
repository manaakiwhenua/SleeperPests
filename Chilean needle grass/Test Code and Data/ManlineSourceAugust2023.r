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

######################################################
######################################################
###Estimate annual eradication probability required to 
###achieve 0.95 prob of achieving eradication at the PATCH scale
###following existing model of seedbank decline
###with management (using glyphosate) to prevent seed production presented in:
###https://www.landcareresearch.co.nz/uploads/public/Discover-Our-Research/Biosecurity/Biocontrol-ecology-of-weeds/3-applications/CNG-review-SFF-project-2010.pdf
###Moderate this by the probability of detecting all known infested patches on a farm
###Use this as management efficay in INAscene function call
######################################################
######################################################
###Seeds per metre square
Nseeds = 1500
###Seed loss from seed bank due to attrition
SeedBankLoss = 0.61
###Seed loss from seed bank due to germination
Germination = 0.064

###Model reduction in seed bank
Ntimesteps = 80
Nseeds = vector(length = Ntimesteps)
N0 = 1500
for(i in 1:Ntimesteps)
{
  Nseeds[i] = N0*(1-(SeedBankLoss+Germination))
  N0 = Nseeds[i] 
}
###Find timestep where seed density <1 seed per hectare
###Consider this effective eradiation 
EradicationYear = length(Nseeds[Nseeds>1/10^4])

###Annual eradication prob to achieve 0.95 errodication
###after 14 (timestep when Nseeds <1 per hectare) timesteps
EradicationProb = 1-(1-0.95)^(1/EradicationYear)

###Moderate eradication prob by probability of identifying all patches on farm
###Value derived from survey data indicating ~50% of farmers are confident in identifying CNG
###And allowing a 50% chance of detecting all patches for those farmers
CompleteDetectionProb = 0.25
EradicationProb = EradicationProb*CompleteDetectionProb 

###########################################################################
###########################################################################
###Read in INApest functions from source files
###########################################################################
###########################################################################

source("Scripts/INApestParallelNoINAscene.R")
source("Scripts/INApestParallel.R")
source("Scripts/INApest.R")
source("Scripts/INApestNoINAscene.R")
source("Scripts/INApestHeatMaps.R")

#############################################################
#############################################################

###Run test source code testing for comparable results
###between INApest function versions
for(cs in 1:3)
source("Scripts/INA__BlindRiver_TestingAugust2023ParallelNoINAScene.R")

###Run test source code testing external info functionality
###for each version of INApest function
cs = 1
Function = "INApest"
source("Scripts/INA__BlindRiver_TestingAugust2023ExternalInfo.R")
Function = "INApestParallel"
source("Scripts/INA__BlindRiver_TestingAugust2023ExternalInfo.R")
Function = "INApestNoINAscene"
source("Scripts/INA__BlindRiver_TestingAugust2023ExternalInfo.R")
Function = "INApestParallelNoINAscene"
source("Scripts/INA__BlindRiver_TestingAugust2023ExternalInfo.R")



###Run test source code
for(cs in 1:3)
source("Scripts/INA__BlindRiver_TestingApril2023v2.R")

###Run test source code implementing parallel processing
###This uses rapid custom function to generate heat maps/gifs 
###for visualisation of results
###Maps not intended for publication
for(cs in 1:3)
  source("Scripts/INA__BlindRiver_TestingApril2023Parallel.R")

###Run test source code implementing simulations
###without calling INAscene
for(cs in 1:3)
  source("Scripts/INA__BlindRiver_TestingApril2023ParallelNoINAscene.R")

###Run test source code illustrating effect of long distance dispersal (LDD)
for(cs in 1:3)
  source("Scripts/INA__BlindRiver_NoLDDParallel.R")

###Run test source code illustrating effect of agency-led communication vs. 
###"Passive" information spread through social networks 
for(cs in 1:3)
  source("Scripts/INA__BlindRiver_SocialNetworksParallel.R")


###Print heat maps illustrating effect long-distance dispersal 
###WARNING: Only runs in R version 4.1.3 on Manaaki Whenua spatial lab system
###Requires various spatial layers and is sensitive to changes in version of 
###some geospatial packages
for(cs in 1:3)
  source("Scripts/NoLDDHeatMaps.R")


###Print heat maps illustrating effect of agency-led communication vs. 
###"Passive" information spread through social networks 
###WARNING: Only runs in R version 4.1.3 on Manaaki Whenua spatial lab system
###Requires various spatial layers and is sensitive to changes in version of 
###some geospatial packages
for(cs in 1:3)
  source("Scripts/SocialNetworkHeatMaps.R")


