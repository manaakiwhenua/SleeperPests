##################################################
###Print invasion heat maps
###for cross-region invasions
##################################################

library(INA)
library(tidyverse)
library(dplyr)
library(actuar)
memory.limit(size=600000) 

###########################################################################
###########################################################################
###Read in function from source file
###This simulates invasion from multiple randomised start points
###Start point probabilities can be weighted by a risk vector
###For instance, spatial or ownership links of farms in a "sink" region
###to farms in a "source region"
###########################################################################

source("INA_Function_Plot_InfestationHeatMaps.R")


######################################################
######################################################
###Estimate annual erradication probability required to 
###achieve 0.95 prob of achieving erradication
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
###Consider this effective erradiation 
ErradicationYear = length(Nseeds[Nseeds>1/10^4])
 
###Annual erradication prob to achieve 0.95 errodication
###after 14 (year when Nseeds <1 per hectare) years
AnnualErradicationProb = 1-(1-0.95)^(1/ErradicationYear)

#############################################################
#############################################################

ResultsDir= paste0(main.dir,"/CrossRegionInvasionCoreFunctionLDDmatrix/",ClimateScenarios[cs],"/")
#dir.create(ResultsDir,recursive = T)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
RegionResultsDir= paste0(ResultsDir,SinkRegion,"_From_",SourceRegion,"/")

###Read in farm polygons
RegionFarms = sf::st_read(paste0(ClimateScenarios[cs],"/Inputs/",SinkRegion,"SheepBeefAtRisk",ClimateScenarios[cs], ".wgs84.shp"))

####################################################
###Input paramaters for INA
####################################################
###Set number of initial configurations
Nconfig = 10
###Set number of realisations
Nperm = 10
###Set simulation duration
Nyears = 80

###Probability of management adoption when infestation detected
ManageProb = 0.9

###Set spread reduction under management
###Assume 1000-fold reduction in spread prob under management
SpreadReduction = 0.999

###Management efficacy defined as annual erradication probability
ManEfficacy = AnnualErradicationProb

###Set probability of detection
###50% farmers confident of CNG ID 
###Default assumes those with ID knowledge have 50% chance of detection annually
###i.e. detection prob = 0.25
###Declare vector to explore other detection probabilities
DetectionProbs = c(0,0.05,0.1,0.15,0.2,0.25,0.3) 



###Loop through different detection probabilities
###Can loop through other levels of management variables too
###Function incorporates management variables in filenames of outputs
for(detprob in 1:length(DetectionProbs))
{
  DetectionProb = DetectionProbs[detprob]
  TitleStem = paste0("Detection Prob. ",DetectionProb," Erradication Prob. ",round(AnnualErradicationProb,digits =3),
                     "\nSpread Reduction ",SpreadReduction)
  InvasionFileNameStem = paste0(RegionResultsDir, "\\"
                                 , fun.basename.DES( d = DetectionProb, e = AnnualErradicationProb, s = SpreadReduction, short = F, norm = T ))

  InvasionFileName <- paste0(InvasionFileNameStem,"_InvasionProb.rds")
  file.exists(InvasionFileName)
  
  InvasionProbResults = readRDS(InvasionFileName)
 
                             
  ImageFileNameStem <- fun.basename.DES( d = DetectionProb, e = AnnualErradicationProb, s = SpreadReduction, short = T,norm = F )
  
  ImageDir = paste0(RegionResultsDir,"\\", "ImpactFigures")
  dir.create(ImageDir)
  #..Call----
  InfestationHeatMaps(
      Nyears = ncol(InvasionProbResults),
      InvasionProbResults = InvasionProbResults,
      ImageDir = ImageDir,
      FarmPolygons = RegionFarms,
      TitleStem = TitleStem,
      FileNameStem = ImageFileNameStem
      , Coastline = fun.map.coastline()
      , Hillshade = fun.map.hillshade()
      , breaks = fun.map.breaks()
      , col = fun.map.colors.tealtored(breaks = fun.map.breaks(), dir = ImageDir)
      , spf = 0.75
      , outlines = F
    )
}

