##################################################################################
###NOTE CAREFULLY: this file contains multiple extremely time-consuming processes
###And should not be called as a source file
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
###Current climate
cs = 1
###Future climate 2040
#cs = 2
###Future climate 2090
#cs = 3
source("CrossRegion_FarmDistance.R")

##############################################
###Cross region distances and distance-based invasion prob 
###between at risk farms
###Uses same Pareto kernel as that used to generate adjaceny matrices
###between farms in the same region 
###Very time consuming so best to run separate climate scenarios in parallel
##############################################
###Current climate
cs = 1
###Future climate 2040
#cs = 2
###Future climate 2090
#cs = 3
source("CrossRegion_FarmDistanceInvasionProb.R")

##################################################
###Calculate long distance invasion probability between individual farms
###within and between regions
###Long distance dispersal assumed to be human-mediated
###Use data on cattle movements from NAIT to weight LDD probability by distance
##################################################

###Current climate
cs = 1
###Future climate 2040
#cs = 2
###Future climate 2090
#cs = 3
source("CrossRegion_FarmDistanceInvasionProbLDD.R")

##################################################
###Weight long distance dispersal probabilities so that
###probabilities add to 1 for each source farm nationally
###Essentially shares dispersal events across potential sink farms
###Using long distance dispersal probabilities as weights
##################################################

###Current climate
#cs = 1
###Future climate 2040
cs = 2
###Future climate 2090
#cs = 3
source("CrossRegion_LDDWeights.R")

#############################################################################
#############################################################################
###End of section
#############################################################################
#############################################################################

#############################################################################
#############################################################################
###Section 3
###This section implements simulations for historic examples
###i.e. where data on known infestations are available
###So far examples are implemented for Marlborough and Hawkes Bay
###Canterbury is the other region were infestations are known. 
#############################################################################
#############################################################################

##############################################
###Historic invasion reproduction
###BlindRiver (Marlborough region) as example
###Include "realistic" control parameters
###Implements previous version of INApest() core function (without some parameter definition options of final version)
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
for(cs in 1:3)
source("INA_HistoricExample_BlindRiverLDDMatrixInfoSpreadv2.R")


###Produce heat maps using model outputs
Region = "MARL"
for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/HistoricExamplesLDDMatrixInfoSpread_v2/",ClimateScenarios[cs],"/")
source("INA_HistoricExampleHeatMaps.R")
}

###Calculate invasion risk from source region to other regions
Region  = "MARL"
for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/HistoricExamplesLDDMatrixInfoSpread_v2/",ClimateScenarios[cs],"/")
source("INA_HistoricExampleCrossRegionRisk.R")
}


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
#############################################################################
#############################################################################

###############################
###Historic invasion
###varying multiple management parameters
################################


#for(cs in 1:3)
cs = 1
source("INA_BlindRiverLDDMatrixInfoSpreadMultManage.R")

###Summarize results
for(cs in 1:3)
	source("INA_BlindRiverLDDMatrixInfoSpreadMultManagePostHocProcessing.R")

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
Region  = "HBAY"
for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/HistoricExamplesLDDMatrixInfoSpread_v2/",ClimateScenarios[cs],"/")
source("INA_HistoricExampleCrossRegionRisk.R")
}

###Produce invasion risk maps from HBAY region to MNWG region
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
###This section implements simulations for cross-region invasions
###i.e. from regions where CNG already occurs to neighbouring regions
###So far examples are implemented for invasions from Marlborough and Hawkes Bay
#############################################################################
#############################################################################

###############################
###Cross-region invasions
###Implementing model developed for Blind River
###with previous version of INApest() core function (without some parameter definition options of final version)
###Farm-level incursion risk scaled by proximity to at-risk farms in source region 
###And climatic suitability
###Use Gisborne as "sink" region to trial code - relatively quick to run
###Permit info spread to nodes within threshold distance of 
###known infestations
###Assign communication rate by authorities and owners of infested nodes
###as information transfer probability within threshold distance
###This implementation also moderates eradication probablity by the probability 
###that all infested patches within a node will be detected
##############################################

DetectionProbs = c(0,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75) 
AnnualEradicationProb = 0.04815904

SinkRegion = "GISB"
SourceRegion = "HBAY"
for(cs in 1:3)
source("INASourceSinkLDDMatrixInfoSpreadv2.R")

for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/CrossRegionInvasionLDDmatrixInfoSpreadv2/",ClimateScenarios[cs],"/")
source("INASourceSinkHeatMaps.R")
}

SinkRegion = "WELL"
SourceRegion = "HBAY"

for(cs in 1:3)
source("INASourceSinkLDDMatrixInfoSpreadv2.R")

for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/CrossRegionInvasionLDDmatrixInfoSpreadv2/",ClimateScenarios[cs],"/")
source("INASourceSinkHeatMaps.R")
}

SinkRegion = "MNWG"
SourceRegion = "HBAY"

###Current climate
cs = 1
###Future climate 2040
#cs = 2
###Future climate 2090
#cs = 3
source("INASourceSinkLDDMatrixInfoSpreadv2.R")

for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/CrossRegionInvasionLDDmatrixInfoSpreadv2/",ClimateScenarios[cs],"/")
source("INASourceSinkHeatMaps.R")
}


SinkRegion = "CANT"
SourceRegion = "MARL"

###Current climate
cs = 1
###Future climate 2040
#cs = 2
###Future climate 2090
#cs = 3
source("INASourceSinkLDDMatrixInfoSpreadv2.R")

for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/CrossRegionInvasionLDDmatrixInfoSpreadv2/",ClimateScenarios[cs],"/")
source("INASourceSinkHeatMaps.R")
}

##############################################
###Try using external invasion threat to initiate invasion 
###and as dynamic ongoing source of invasion
##############################################


DetectionProbs = c(0,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75) 
AnnualEradicationProb = 0.04815904

############################################################################################
###Uses dynamic external invasion threat derived from zero management scenario in source region (HBAY)
############################################################################################

SinkRegion = "GISB"
SourceRegion = "HBAY"
for(cs in 1:3)
source("INASourceSinkLDDMatrixInfoSpreadConstantExternal.R")

for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/CrossRegionConstantThreat/",ClimateScenarios[cs],"/")
source("INASourceSinkHeatMaps.R")
}


SinkRegion = "MNWG"
SourceRegion = "HBAY"
#for(cs in 1:3)
cs = 1
source("INASourceSinkLDDMatrixInfoSpreadConstantExternal.R")

#cs = 2
source("INASourceSinkLDDMatrixInfoSpreadConstantExternal.R")

#cs = 3
source("INASourceSinkLDDMatrixInfoSpreadConstantExternal.R")

for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/CrossRegionConstantThreat/",ClimateScenarios[cs],"/")
source("INASourceSinkHeatMaps.R")
}


##############################################
###Do MWNG with no info spread
###and with ongoing invasion from HBAY
##############################################

SinkRegion = "MNWG"
SourceRegion = "HBAY"
#for(cs in 1:3)
cs = 1
source("INASourceSinkLDDMatrixNoInfoSpreadConstantExternal.R")

#cs = 2
source("INASourceSinkLDDMatrixNoInfoSpreadConstantExternal.R")
#cs = 3
source("INASourceSinkLDDMatrixNoInfoSpreadConstantExternal.R")


for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/CrossRegionConstantThreatNoInfoSpread/",ClimateScenarios[cs],"/")
source("INASourceSinkHeatMaps.R")
}


############################################################################################
###Uses external invasion derived from detection prob = 0.2 management scenario in source region (HBAY)
############################################################################################

SinkRegion = "GISB"
SourceRegion = "HBAY"
for(cs in 1:3)
source("INASourceSinkLDDMatrixInfoSpreadConstantExternalv2.R")

for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/CrossRegionConstantThreatv2/",ClimateScenarios[cs],"/")
source("INASourceSinkHeatMaps.R")
}


SinkRegion = "MNWG"
SourceRegion = "HBAY"
#for(cs in 1:3)
cs = 1
source("INASourceSinkLDDMatrixInfoSpreadConstantExternalv2.R")

#cs = 2
source("INASourceSinkLDDMatrixInfoSpreadConstantExternalv2.R")

#cs = 3
source("INASourceSinkLDDMatrixInfoSpreadConstantExternalv2.R")

for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/CrossRegionConstantThreatv2/",ClimateScenarios[cs],"/")
source("INASourceSinkHeatMaps.R")
}

##############################################
###Do MWNG with no info spread
###and with ongoing invasion from HBAY
##############################################


############################################################################################
###Uses external invasion derived from detection prob = 0.2 management scenario in source region (HBAY)
############################################################################################

SinkRegion = "MNWG"
SourceRegion = "HBAY"
#for(cs in 1:3)
cs = 1
source("INASourceSinkLDDMatrixNoInfoSpreadConstantExternalv2.R")

#cs = 2
source("INASourceSinkLDDMatrixNoInfoSpreadConstantExternalv2.R")
#cs = 3
source("INASourceSinkLDDMatrixNoInfoSpreadConstantExternalv2.R")

for(cs in 1:3)
{
ResultsDir= paste0(main.dir,"/CrossRegionConstantThreatNoInfoSpreadv2/",ClimateScenarios[cs],"/")
source("INASourceSinkHeatMaps.R")
}


#############################################################################
#############################################################################
###End of section
#############################################################################
#############################################################################






