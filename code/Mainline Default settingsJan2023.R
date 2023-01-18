##################################################################################
###NOTE CAREFULLY: this file contains multiple extremely time-consuming processes
###And should not be called as a source file
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
###Prints maps ecoclimatic index values and establishment probability
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
#cs = 1
###Future climate 2040
cs = 2
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
#cs = 1
###Future climate 2040
cs = 2
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
#cs = 1
###Future climate 2040
cs = 2
###Future climate 2090
#cs = 3
source("CrossRegion_FarmDistanceInvasionProbLDD.R")

##################################################
###Weight long distance dispersal probabilities so that
###probabilities add to 1 for each source farm nationally
###Essentially shares dispersal events across potential sink farms
###Using long distance dispersal probabilities
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
###Canterbury is the other region were infestation are known. 
#############################################################################
#############################################################################

##############################################
###Historic invasion reproduction
###BlindRiver (Marlborough region) as example
###Include "realistic" control parameters
###Vary annual detection probability
###Control applied as dispersal probability reduction
###And annual local erradication probability 
###Allow long distance dispersal scaled by NZ cattle movement vs distance
###Stores large outputs for post-hoc analyses
###Prints summary line graphs and heat maps of changing invasion probability through time
##############################################
###Current climate
#cs = 1
###Future climate 2040
#cs = 2
###Future climate 2090
#cs = 3
###Long-distance dispersal independent of distance (equiprobable for all pairs of farms)
#for(cs in 1:3)
#source("INA_HistoricExample_BlindRiverWithCoreFunction.R")

###This version uses a LDD matrix based in NZ cattle movements
for(cs in 1:3)
source("INA_HistoricExample_BlindRiverWithCoreFunctionLDDmatrix.R")

###############################
###Produce heat map jpegs and gifs
##for Blind River historic example
################################
region = 10
for(cs in 1:3)
source("INA_HistoricExampleHeatMaps.R")

###############################
###Produce line graphs of impacts on
##farm area and livestock number
################################
region = 10
for(cs in 1:3)
source("INA_HistoricExampleImpacts.R")

##############################################
##############################################


##############################################
###Historic invasion 
###HBAY as example
###Implements model developed for Blind River (Marlborough region)
###With infestation data as of 2021
###Run for 50 years to reach ~2070
###Roughly match Marlborough timeframe
##############################################
###Current climate
cs = 1
###Future climate 2040
#cs = 2
###Future climate 2090
#cs = 3

#for(cs in 1:3)
source("INA_HistoricExample_HBAYWithCoreFunctionLDDMatrix.R")


###############################
###Produce heat map jpegs and gifs
##for Hawkes Bay historic example
################################
#for(cs in 1:3)
source("INA_HistoricExampleHeatMapsHBAY.R")

###############################
###Produce line graphs of impacts on
##farm area and livestock number
################################
region = 7
for(cs in 1:3)
source("INA_HistoricExampleImpacts.R")

##############################################
##############################################

##############################################
###Historic invasion reproduction
###BlindRiver (Marlborough region) as example
###Include "realistic" control parameters
###Permit info spread to nodes within threshold distance of 
###known infestations
###Assign communication rate by authorities and owners of infested nodes
###as information transfer probability within threshold distance
###This implementation also moderates erradication probablity by the probability 
###that all infested patches with in a node will be detected
##############################################

for(cs in 1:3)
source("INA_HistoricExample_BlindRiverLDDMatrixInfoSpread.R")






#############################################################################
#############################################################################
###End of section
#############################################################################
#############################################################################


#############################################################################
#############################################################################
###Section 4
###This section implements simulations for cross-region invasions
###i.e. from regions where CNG already occurs to neighbouring regions
###So far examples are implemented for invasions from Marlborough and Hawkes Bay
#############################################################################
#############################################################################

###############################
###Cross-region invasions
###Implementing model developed for Blind River
###with a standardised core function
###Farm-level incursion risk scaled by proximity to farms in source region 
###And climatic suitability
###Use Gisborne as "sink" region to trial code - relatively quick to run
################################

###Current climate
cs = 1
###Future climate 2040
#cs = 2
###Future climate 2090
#cs = 3

###This version uses a LDD matrix based in NZ cattle movements
SinkRegion = "GISB"
SourceRegion = "HBAY"
#for(cs in 1:3)
source("INASourceSinkWithCoreFunctionLDDmatrix.R")
source("INASourceSinkHeatMaps.R")

for(cs in 1:3)
source("INASourceSinkImpacts.R")

###Other options for invasion from Hawkes Bay are 
###Greater Wellington and Manawatu-Wanganui regions
###MNWG in particular is more time-consuming

SinkRegion = "WELL"
SourceRegion = "HBAY"
cs = 1
#cs = 2
#cs = 3
source("INASourceSinkWithCoreFunctionLDDMatrix.R")
source("INASourceSinkHeatMaps.R")

SinkRegion = "MNWG"
SourceRegion = "HBAY"
cs = 1
#cs = 2
#cs = 3
source("INASourceSinkWithCoreFunctionLDDMatrix.R")
source("INASourceSinkHeatMaps.R")


###Examine invasion to Canterbury from Marlborough
SinkRegion = "CANT"
SourceRegion = "MARL"
cs = 1
#cs = 2
#cs = 3
source("INASourceSinkWithCoreFunctionLDDMatrix.R")
source("INASourceSinkHeatMaps.R")

###############################
###Produce heat map jpegs and gifs
################################

SinkRegion = "GISB"
SourceRegion = "HBAY"
for(cs in 1:3)
source("INASourceSinkHeatMaps.R")

SinkRegion = "MNWG"
SourceRegion = "HBAY"
for(cs in 1:3)
source("INASourceSinkHeatMaps.R")

SinkRegion = "WELL"
SourceRegion = "HBAY"
for(cs in 1:3)
source("INASourceSinkHeatMaps.R")

SinkRegion = "CANT"
SourceRegion = "MARL"
for(cs in 1:3)
source("INASourceSinkHeatMaps.R")


#############################################################################
#############################################################################
###End of section
#############################################################################
#############################################################################



#############################################################################
#############################################################################
###Section 5
###Explore sensitivity of outcomes to variation in multiple control parameters
###Use Blind River hostoric invasion as an example
#############################################################################
#############################################################################

###############################
###Historic invasions
###varying multiple management parameters
################################

for(cs in 1:3)
source("INA_HistoricExample_MultipleManagement.R")

#############################################################################
#############################################################################
###End of section
#############################################################################
#############################################################################


