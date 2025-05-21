##################################################
###Print invasion heat maps from spread simulated spread from
###infested farms in Blind River
###Initial invasion approximates infestation documented here:
###https://www.tandfonline.com/doi/abs/10.1080/0028825X.1989.10414122
##################################################

library(INA)
library(tidyverse)
library(dplyr)
memory.limit(size=600000) 

###########################################################################
###########################################################################
###Read in invasion heat map function
###########################################################################
###########################################################################

source("INA_Function_Plot_InfestationHeatMaps_V03.R")


######################################################
######################################################
###Estimate annual eradication probability required to 
###achieve 0.95 prob of achieving eradication
###following existing model of seedbank decline
###with management (using glyphosate) to prevent seed production
###use this as management efficay in INAscene function call
######################################################
######################################################
###Seeds per metre square
Nseeds = 1500
###Seed loss from seed bank due to attrition
SeedBankLoss = 0.61
###Seed loss from seed bank due to germination
Germination = 0.064

###Model reduction in seed bank
Nyears = 50
Nseeds = vector(length = Nyears)
N0 = 1500
for(i in 1:Nyears)
 {
 Nseeds[i] = N0*(1-(SeedBankLoss+Germination))
 N0 = Nseeds[i] 
 }
###Find year where seed density <1 seed per hectare
###Consider this effective eradiation 
EradicationYear = length(Nseeds[Nseeds>1/10^4])
 
###Annual eradication prob to achieve 0.95 errodication
###after 14 (year when Nseeds <1 per hectare) years
#AnnualEradicationProb = 1-(1-0.95)^(1/EradicationYear)

#############################################################
#############################################################


#ResultsDir= paste0(main.dir,"/HistoricExamplesCoreFunctionLDDMatrix/",ClimateScenarios[cs],"/")
#dir.create(ResultsDir,recursive = T)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
RegionResultsDir= paste0(ResultsDir,Region,"/")
####################################################
###Input paramaters for INA
####################################################
###Set number of realisations
Nperm = 100
###Set simulation duration
Nyears = 80

###Probability of management adoption when infestation detected
ManageProb = 0.9

###Set spread reduction under management
###Assume 1000-fold reduction in spread prob under management
SpreadReduction = 0.999

###Management efficacy defined as annual eradication probability
ManEfficacy = AnnualEradicationProb

###Set probability of detection
###50% farmers confident of CNG ID 
###Default assumes those with ID knowledge have 50% chance of detection annually
###Declare vector to explore other probs
#DetectionProbs = c(0,0.05,0.1,0.15,0.2,0.25,0.3) 

CrossRegionThreatDir = paste0(RegionResultsDir,"CrossRegionThreat/") 
###Loop through different detection probabilities
###Read in farm polygons

ThreatSource = Region
#for(sink in 1:length(Regions))
#if(Regions[sink] != Region)
#{
#ThreatSink = Regions[sink]
RegionFarms = sf::st_read(paste0(ClimateScenarios[cs],"/Inputs/",ThreatSink,"SheepBeefAtRisk",ClimateScenarios[cs], ".wgs84.shp"))
for(detprob in 1:length(DetectionProbs))
{
DetectionProb = DetectionProbs[detprob]
 
  TitleStem = Title = paste0("Source ", Region," Detection Prob. ",round(mean(DetectionProb),digits = 3))
   OutputFileNameStem = paste0(CrossRegionThreatDir,"DetProb_",DetectionProb,"_EradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
  
   InvasionFileName = paste0(OutputFileNameStem,ThreatSource,"_",ThreatSink,"FarmRisk.rds")
   InvasionProbResults = readRDS(InvasionFileName)
   plot(1:Nyears,colSums(InvasionProbResults))
   hist(InvasionProbResults)  
   max(InvasionProbResults)
   min(InvasionProbResults)                         
  ImageFileNameStem <- paste0("DetProb_",DetectionProb)
  
  ImageDir = paste0(CrossRegionThreatDir,"\\", "InvasionHeatMaps","\\",ThreatSink)
  dir.create(ImageDir, showWarnings = F)
  cols.invasion = fun.map.colors.tealtored(breaks = c(0.00001,0.0001,0.001,0.01,0.05,0.1,1), dir = ImageDir)
  #..Call----
  InfestationHeatMaps(
      Nyears = Nyears,
      InvasionProbResults = InvasionProbResults,
      ImageDir = ImageDir,
      FarmPolygons = RegionFarms,
      TitleStem = TitleStem,
      FileNameStem = ImageFileNameStem,
      region = ThreatSink,
      , Coastline = fun.map.coastline()
      , Hillshade = fun.map.hillshade()
      , breaks =  c(0.00001,0.0001,0.001,0.01,0.05,0.1,1)
      , cols.invasion = cols.invasion
      , spf = 0.75
      , outlines = F
    )
}




