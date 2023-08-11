##################################################
###Print invasion heat maps from spread simulated spread from
###infested farms in Blind River
###Initial invasion approximates infestation documented here:
###https://www.tandfonline.com/doi/abs/10.1080/0028825X.1989.10414122
###WARNING: Only runs in R version 4.1.3 on Manaaki Whenua spatial lab system
###Requires various spatial layers and is sensitive to changes in version of 
###some geospatial packages
##################################################
#install.packages("colorblindcheck")
#install.packages("terra")
#install.packages("sp")

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
library(INA)
library(tidyverse)
library(dplyr)
###########################################################################
###########################################################################
###Read in invasion heat map function
###########################################################################
###########################################################################

source("Scripts/INA_Function_Plot_InfestationHeatMaps_V03.R")
Region = "MARL"
#ResultsDir= paste0(main.dir,"/HistoricExamplesCoreFunctionLDDMatrix/",ClimateScenarios[cs],"/")
#dir.create(ResultsDir,recursive = T)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
ResultsDir= paste0(main.dir,"/Results/NoLDD/",ClimateScenarios[cs],"/")
RegionResultsDir= paste0(ResultsDir,Region,"/")
###Read in farm polygons
RegionFarms = sf::st_read(paste0(main.dir,"/",ClimateScenarios[cs],"/",Region,"SheepBeefAtRisk",ClimateScenarios[cs], ".wgs84.shp"))

####################################################
###Input paramaters for INA
####################################################
###Set number of realisations
Nperm = 100
###Set simulation duration
Ntimesteps = 80

###Set probability of detection
DetectionProbs = c(0,0.025) 


##################################################
###Loop through different detection probabilities
##################################################
for(detprob in 1:length(DetectionProbs))
{
DetectionProb = DetectionProbs[detprob]
###Standard settings
###Declare name for storing results
ModelName = paste0("NoLDD","DetProb_",DetectionProb)
OutputDir = paste0(RegionResultsDir,ModelName,"/")  
TitleStem = ModelName
InvasionFileNameStem =paste0(OutputDir,ModelName)
InvasionFileName <- paste0(InvasionFileNameStem,"InvasionProb.rds")
 
InvasionProbResults = readRDS(InvasionFileName)
                              
  ImageFileNameStem <- ModelName
  
  ImageDir = paste0(RegionResultsDir,"\\", "InvasionHeatMaps")
  cols.invasion = fun.map.colors.tealtored(breaks = fun.map.breaks(), dir = ImageDir)
  #..Call----
  InfestationHeatMaps(
      Nyears = Ntimesteps,
      InvasionProbResults = InvasionProbResults,
      ImageDir = ImageDir,
      FarmPolygons = RegionFarms,
      TitleStem = TitleStem,
      FileNameStem = ImageFileNameStem,
      region = Region,
      Coastline = fun.map.coastline()
      , Hillshade = fun.map.hillshade()
      , breaks = fun.map.breaks()
      , cols.invasion = cols.invasion
      , spf = 0.75
      , outlines = F
    )
ModelName = paste0("WithLDD","DetProb_",DetectionProb)
OutputDir = paste0(RegionResultsDir,ModelName,"/")
  TitleStem = ModelName
  InvasionFileNameStem =paste0(OutputDir,ModelName)
  InvasionFileName <- paste0(InvasionFileNameStem,"InvasionProb.rds")
 
  InvasionProbResults = readRDS(InvasionFileName)
                              
  ImageFileNameStem <- ModelName
  
  ImageDir = paste0(RegionResultsDir,"\\", "InvasionHeatMaps")
  cols.invasion = fun.map.colors.tealtored(breaks = fun.map.breaks(), dir = ImageDir)
  #..Call----
  InfestationHeatMaps(
      Nyears = Ntimesteps,
      InvasionProbResults = InvasionProbResults,
      ImageDir = ImageDir,
      FarmPolygons = RegionFarms,
      TitleStem = TitleStem,
      FileNameStem = ImageFileNameStem,
      region = Region,
      Coastline = fun.map.coastline()
      , Hillshade = fun.map.hillshade()
      , breaks = fun.map.breaks()
      , cols.invasion = cols.invasion
      , spf = 0.75
      , outlines = F
    )
}




