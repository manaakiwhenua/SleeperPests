##################################################
###Simulate spread from infested farms in Blind River
###Initial invasion approximates infestation documented here:
###https://www.tandfonline.com/doi/abs/10.1080/0028825X.1989.10414122
###Apply management as 1000-fold reduction in dispersal 
###And local eradication following 
###https://www.landcareresearch.co.nz/uploads/public/Discover-Our-Research/Biosecurity/Biocontrol-ecology-of-weeds/3-applications/CNG-review-SFF-project-2010.pdf
###Calls core function and provides binary vector of initial infestations
###Records farms infested, farms under management, number and proportion of extant infestations detected 
###and invasive threat to farms in neighbouring regions
###This version implements information spread to
###farms within a threshold distance of known infestations
###Represents management scenario where authorities communicate with 
###neighbours of farms with known infestations
##################################################

library(INA)
library(tidyverse)
library(dplyr)
memory.limit(size=600000) 

###########################################################################
###########################################################################
###Read in the function from source file
###########################################################################
###########################################################################

source("INApest.R")
source("CrossRegionInvasionThreat.R")


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
Nyears = 80
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
AnnualEradicationProb = 1-(1-0.95)^(1/EradicationYear)

###Moderate eradication prob by probability of identifying all patches on farm
###Value derived from survey data indicating ~50% of farmers are confident in identifying CNG
###And allowing a 50% chance of detecting all patches for those farmers
CompleteDetectionProb = 0.25
AnnualEradicationProb = AnnualEradicationProb*CompleteDetectionProb 

#############################################################
#############################################################


#############################################################
###Define directories for storing results
#############################################################

#ResultsDir= paste0(main.dir,"/HistoricExamplesLDDMatrixInfoSpread_v2/",ClimateScenarios[cs],"/")
dir.create(ResultsDir,recursive = T)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
RegionResultsDir= paste0(ResultsDir,Region,"/")
dir.create(RegionResultsDir)
CrossRegionThreatDir = paste0(RegionResultsDir,"CrossRegionThreat/") 
dir.create(CrossRegionThreatDir)

####################################################
###Input paramaters for INA
####################################################
###Set number of realisations
Nperm = 100
###Set simulation duration
Nyears = 80

LongDistProbPerFarm = 0.05
###Probability of management adoption when infestation detected
ManageProb = 0.9

###Set spread reduction under management
###Assume 1000-fold reduction in spread prob under management
SpreadReduction = 0.999

###Set probability of detection
###50% farmers confident of CNG ID 
###Default assumes those with ID knowledge have 50% chance of detection annually
###Declare vector to explore other probs
#DetectionProbs = c(0,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75) 

##################################################
###Call cross-region invasion function
###Looping through different detection probabilities
##################################################
for(detprob in 1:length(DetectionProbs))
{
DetectionProb = DetectionProbs[detprob]
##################################################
###Use INApest outputs to calculate 
###Invasion threat to other regions
##################################################
LongDistProbPerFarm = 0.05
ThreatSource = Region
for(sink in 1:length(Regions))
if(Regions[sink] != Region)
{
ThreatSink = Regions[sink]
InputDir = paste0(main.dir,"/",ClimateScenarios[cs],"/Inputs/")
OutputDir = CrossRegionThreatDir
dir.create(OutputDir)
CrossRegionInvasionDir = paste0(main.dir,paste0("/",ClimateScenarios[cs],"/Inputs/CrossRegionDistance/"))
InvasionFileNameStem = paste0(RegionResultsDir,"DetProb_",DetectionProb,"_EradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
OutputFileNameStem = paste0(OutputDir,"DetProb_",DetectionProb,"_EradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
SourceInvasionProb = readRDS(paste0(InvasionFileNameStem,"InvasionProb.rds"))

CrossRegionInvasionThreat(
ThreatSource = ThreatSource, ###Source region
ThreatSink = ThreatSink, ###Sink region
LDDrate = LongDistProbPerFarm, ###annual long distance dispersal events per source farm
SourceInvasionProb = SourceInvasionProb, ###Annual farm-level invasion prob from simulations
CrossRegionInvasionDir=CrossRegionInvasionDir , ###Directory where cross-region invasion probs stored
EIDir = InputDir, ##Directory where Ecoclimatic Index data stored
OutputFileNameStem = OutputFileNameStem
)
}
##################################################
##################################################

}
##################################################
###End of detection prob loop
##################################################


c1 <- rbind( 10, 10, 10, 10)
c2 <- rbind( 1, 1, 2, 2)
df <- as.data.frame( cbind( c1, c2)) 
df %>% dplyr::group_by(V2 ) %>% summarise(total = sum( V1)) 
df %>% dplyr::group_by(V2 ) %>% summarise(total = sum( V1, na.rm = T))