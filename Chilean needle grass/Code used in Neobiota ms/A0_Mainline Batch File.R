##################################################################################
###NOTE CAREFULLY: this file refers to scripts used in producing results presented
###in the Neobiota manuscript.
###Due to the size of input data an processing times it is not designed to repeat the analyses
###used to produce results
##################################################################################

##################################################################################
###When SLMACC final report completed will link analyses/source files to relevant sections
##################################################################################

#############################################################################
#############################################################################
###Section 1
###This section sets R library version and defines location of 
###key data sources - Agribase shape file and CLIMEX values under current and future climates
###It also sets the main directory where all outputs are stored
###And defines the type of relationship converting CLIMEX ecoclimatic index values 
###to establishment probability
#############################################################################
#############################################################################
rm(list = ls())
val.path <- paste0( "N:/Projects/BaseData/RPackages/R-4.1.3")
if (dir.exists(val.path))
{
  .libPaths(val.path)
  .libPaths()
} else
{
  cat(val.path, "\r\n")
  stop( "Cannot find the lib path specified")
}
setwd(r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CurrentFutureClimate]")
main.dir = r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CurrentFutureClimate]"
memory.limit(size = 600000)
#dir.create("Inputs/")
ClimexFile = r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\NewClimexData\CNG_NiwaA1B_CS.shp]"
ClimexVars = c("EI","EI_Niwa_204","EI_Niwa209")
ClimateScenarios = c("Current","Future_2040","Future_2090")
AgribaseLocation = r"[N:\Projects\BaseData\NZ\LandUse\Agribase\Agribase_202210\agribase-october-2022.shp]"
EI_Prob_CurveType = "Logit"
#EI_Prob_CurveType = "SplitLinear"



#############################################################################
#############################################################################
###End of section
#############################################################################
#############################################################################


#############################################################################
#############################################################################
###Section 2
###This section generates various inputs required for sumulations
#############################################################################
#############################################################################

##############################################
###Identifies climatically suitable farms for each region
###Converts ecoclimatic index values to establishment probability
###Calculates distance between farms within each region
###Prints maps of ecoclimatic index values and establishment probability
###Outputs for different climate scenarios stored in separate folders
#############################################
for(cs in 1:3)
{
dir.create(paste0(ClimateScenarios[cs],"/Inputs/"),recursive = T)
climexvar = ClimexVars[cs]
source("AgribaseClimexDistMatrix.R")
}

##################################################
###Generates distance matrices between farms in different regions
###Loops through each pair of regions and stores results for each pair separately
###Calculates minimum Euclidean distance between farm boundaries
###so that farms with parcels in multiple regions may have close neighbours in different regions
##################################################
for(cs in 1:3)
source("CrossRegion_FarmDistance.R")

##############################################
###Cross region distances and distance-based invasion prob 
###between at risk farms
###Uses same Pareto kernel as that used to generate adjaceny matrices
###between farms in the same region 
###Very time consuming so best to run separate climate scenarios in parallel
##############################################
for(cs in 1:3)
source("CrossRegion_FarmDistanceInvasionProb.R")

##################################################
###Calculate long distance invasion probability between individual farms
###within and between regions
###Long distance dispersal assumed to be human-mediated
###Use data on cattle movements from NAIT to weight LDD probability by distance
##################################################

for(cs in 1:3)
source("CrossRegion_FarmDistanceInvasionProbLDD.R")

##################################################
###Weight long distance dispersal probabilities so that
###probabilities add to 1 for each source farm nationally
###Essentially shares dispersal events across potential sink farms
###Using long distance dispersal probabilities as weights
##################################################

for(cs in 1:3)
source("CrossRegion_LDDWeights.R")

#############################################################################
#############################################################################
###End of section
#############################################################################
#############################################################################

#############################################################################
#############################################################################
###Section 3
###This section implements simulations for a historic example
###in the Marlborough region.
#############################################################################
#############################################################################

##############################################
###Historic invasion reproduction
###BlindRiver (Marlborough region) as example
###Implements previous version of INApest() core function (basic functionality identical but 
###without some parameter definition options of final version)
###Vary annual detection probability
###Control applied as dispersal probability reduction
###And annual local eradication probability 
###Allow long distance dispersal scaled by NZ cattle movement vs distance
###Permit info spread to nodes within threshold distance of 
###known infestations
###Assign communication rate by authorities and owners of infested nodes
###as information transfer probability within threshold distance
###This implementation also moderates erradication probablity by the probability 
###that all infested patches with in a node will be detected
###Stores outputs for post-hoc analyses
###Prints summary line graphs and heat maps of changing invasion probability through time
##############################################

DetectionProbs = c(0,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75) 
AnnualEradicationProb = 0.04815904

###Run model for multiple detection probability levels under each climate scenario
###DetectionProb = 0 and cs = 1 produce results presented in Figure 2.
for(cs in 1:3)
source("INA_HistoricExample_BlindRiverLDDMatrixInfoSpreadv2.R")

###Produce Figure 2 in manuscript
source("MarlboroughValidation.R")

#############################################################################
#############################################################################
###End of section
#############################################################################
#############################################################################

#############################################################################
#############################################################################
###Section 4
###Explore sensitivity of outcomes to variation in multiple control parameters
###Use Blind River historic invasion as an example
###Produces results presented in Figure 4
###And used in xgboost analyses (presented in Figure 3 and Table 1)
#############################################################################
#############################################################################

###############################
###Historic invasion
###varying multiple management parameters
################################


for(cs in 1:3)
source("INA_BlindRiverLDDMatrixInfoSpreadMultManage.R")

###Summarize resultsfot post hoc analyses in xgboost
for(cs in 1:3)
	source("INA_BlindRiverLDDMatrixInfoSpreadMultManagePostHocProcessing.R")

###Perform xgboost analyses for infested farm years against management and climate variables
###(presented in Figure 3 and Table 1)
source("MultManage.gbm.R")

###Perform xgboost analyses for management success against management and climate variables
###(presented in Figure 3 and Table 1)
source("MultManage.gbmLogistic.R")

#############################################################################
#############################################################################
###End of section
#############################################################################
#############################################################################

########################################################
########################################################
###Section 5: Apply model model validated in Blind River (Marlborough) to existing infestation in Hawkes Bay
########################################################
########################################################

for(cs in 1:3)
source("INA_HistoricExample_HBAYLDDMatrixInfoSpreadv2.R")


###Produce heat maps using model outputs
Region = "HBAY"
for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/HistoricExamplesLDDMatrixInfoSpread_v2/",ClimateScenarios[cs],"/")
source("INA_HistoricExampleHeatMaps.R")
}



###Calculate invasion risk from source region to other regions
###Produces results presented in Figure 5
Region  = "HBAY"
for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/HistoricExamplesLDDMatrixInfoSpread_v2/",ClimateScenarios[cs],"/")
source("INA_HistoricExampleCrossRegionRisk.R")
}

###Produce invasion risk maps from HBAY region to MNWG region
###Produces maps presented in Figure 5
Region  = "HBAY"
ThreatSink = "MNWG"
DetectionProbs = c(0,0.2)
for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/HistoricExamplesLDDMatrixInfoSpread_v2/",ClimateScenarios[cs],"/")
source("INA_HistoricExampleCrossRegionThreatHeatMaps.R")
}

#############################################################################
#############################################################################
###End of section
#############################################################################
#############################################################################



#############################################################################
#############################################################################
###Section 6
###This section implements simulations for spread in the sink region (MNWG) 
###with ongoing cross-region invasions from the source region (HBAY)
#############################################################################
#############################################################################

###############################
###Cross-region invasions
###Implementing model developed for Blind River
###with previous version of INApest() core function (core functionality identical but 
###without some parameter definition options of final version)
###Farm-level incursion risk scaled by proximity to at-risk farms in source region 
###And climatic suitability
###Permit info spread to nodes within threshold distance of 
###known infestations
###Assign communication rate by authorities and owners of infested nodes
###as information transfer probability within threshold distance
###This implementation also moderates eradication probability by the probability 
###that all infested patches within a node will be detected
##############################################

DetectionProbs = c(0,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75) 
AnnualEradicationProb = 0.04815904

############################################################################################
###Use dynamic external invasion threat derived from zero management (detection prob = 0) scenario in source region (HBAY)
######Produces results for "No control at source" scenarios presented in Figures 6 and 7
############################################################################################

SinkRegion = "MNWG"
SourceRegion = "HBAY"

###Simulate spread in sink region given external invasion risk from source region
###Under no management scenario (detection prob = 0) in source region and with 
###different levels of detection prob in sink region
###Produces results for "No control at source" scenarios presented in Figures 6 and 7
for(cs in 1:3)
source("INASourceSinkLDDMatrixInfoSpreadConstantExternal.R")

###Produces maps for "No control at source" scenarios presented in Figure 7
for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/CrossRegionConstantThreat/",ClimateScenarios[cs],"/")
source("INASourceSinkHeatMaps.R")
}

############################################################################################
###Use dynamic external invasion threat derived from moderate management (detection prob = 0.2) scenario in source region (HBAY)
###Produces results for "With control at source" scenarios presented in Figures 6 and 7
###And "With communication" scenarios in Figure 8
############################################################################################

SinkRegion = "MNWG"
SourceRegion = "HBAY"

###Simulate spread in sink region given external invasion risk from source region
###Under no management scenario (detection prob = 0) in source region and with 
###different levels of detection prob in sink region
###Produces results for "No control at source" scenarios presented in Figures 6 and 7
###And "With communication" scenarios in Figure 8
for(cs in 1:3)
  source("INASourceSinkLDDMatrixInfoSpreadConstantExternalv2.R")

###Produces maps for "With control at source" scenarios presented in Figure 7
###And "With communication" scenarios in Figure 8
for(cs in 1:3)
{
  ResultsDir= paste0(main.dir,"/CrossRegionConstantThreatv2/",ClimateScenarios[cs],"/")
  source("INASourceSinkHeatMaps.R")
}


##############################################
###Simulate spread in sink region  (MWNG) with no info spread
###between farms in the sink region
###and with ongoing invasion from source region (HBAY)
##############################################

############################################################################################
###Uses external invasion derived from detection prob = 0.2 management scenario in source region (HBAY)
############################################################################################

SinkRegion = "MNWG"
SourceRegion = "HBAY"

###Simulate spread in sink region given external invasion risk from source region
###Under moderate management scenario (detection prob = 0.2) in source region and with
###no communication between farms in sink region 
###different levels of detection prob in sink region and 
###Produces results for "No communication" scenarios in Figure 8
for(cs in 1:3)
source("INASourceSinkLDDMatrixInfoSpreadConstantExternalv2.R")

###Produces maps for "No communication" scenarios in Figure 8
for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/CrossRegionConstantThreatv2/",ClimateScenarios[cs],"/")
source("INASourceSinkHeatMaps.R")
}

#############################################################################
#############################################################################
###End of section
#############################################################################
#############################################################################






