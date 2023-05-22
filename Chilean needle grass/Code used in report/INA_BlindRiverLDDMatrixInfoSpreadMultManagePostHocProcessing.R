##################################################
###Post processing of results from 
###"INA_BlindRiverLDDMatrixInfoSpreadMultManage.R"
###Source code
##################################################


library(tidyverse)
library(dplyr)
memory.limit(size=600000) 

######################################################
######################################################
###Estimate annual Eradication probability required to 
###achieve 0.95 prob of achieving Eradication
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
Nyears = 80
Nseeds = vector(length = Nyears)
N0 = 1500
for(i in 1:Nyears)
 {
 Nseeds[i] = N0*(1-(SeedBankLoss+Germination))
 N0 = Nseeds[i] 
 }
###Find year where seed density <1 seed per hectare
###Consider this effective erradiation 
EradicationYear = length(Nseeds[Nseeds>1/10^4])
 
###Annual Eradication prob to achieve 0.95 errodication
###after 14 (year when Nseeds <1 per hectare) years
AnnualEradicationProb = 1-(1-0.95)^(1/EradicationYear)

###Moderate Eradication prob by probability of identifying all patches on farm
CompleteDetectionProb = 0.25
AnnualEradicationProb = AnnualEradicationProb*CompleteDetectionProb 
#############################################################
#############################################################

##############################################################################
###Set region results directory and read in farm agribase data
##############################################################################

ResultsDir= paste0(main.dir,"/HistoricExamplesLDDMatrixInfoSpread/",ClimateScenarios[cs],"/")
dir.create(ResultsDir,recursive = T)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
region = 10
RegionResultsDir= paste0(ResultsDir,Regions[region],"/")
dir.create(RegionResultsDir)
CrossRegionThreatDir = paste0(RegionResultsDir,"CrossRegionThreat/") 
dir.create(CrossRegionThreatDir)

###Read in farm polygons
RegionFarms = sf::st_read(paste0(ClimateScenarios[cs],"/Inputs/",Regions[region],"SheepBeefAtRisk",ClimateScenarios[cs], ".wgs84.shp"))

########################################################
###Calculate distance-based dispersal probability
###and set long distance dispersal rate per farm
#############################################################

###Read in distance matrix
dist_mat = readRDS(paste0(ClimateScenarios[cs],"/Inputs/",Regions[region],"_dist_mat_farm.rds"))



###Store names for nodes in network
FarmNames = row.names(dist_mat)



####################################################################
###Read in historic invasion info
####################################################################
InputDir = r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CurrentFutureClimate\Inputs\]"
HistoricInfestation = read.csv(paste0(InputDir,"BlindRiverInfestedFarms.txt"),header = T, as.is =T)


####################################################
###Input paramaters for INA
####################################################
###Set number of realisations
Nperm = 30
###Set simulation duration
Nyears = 80

###Probability of management adoption when infestation detected
ManageProb = 0.9

###Management efficacy defined as annual Eradication probability
ManEfficacy = AnnualEradicationProb

###Set prob of communication with neighbours of detected infestations
CommProbs = c(0,0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.75)


###Set Eradication prob and spread reduction under management
SpreadReductions = c(0.05,0.1,0.3,0.5,0.9,0.99,0.999)
EradicationProbs = c(0,0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.75)



###Set probability of detection
###50% farmers confident of CNG ID 
###Default assumes those with ID knowledge have 50% chance of detection annually
###Declare vector to explore other probs
DetectionProbs = c(0,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75) 


##################################################
###Loop through different parameter combinations
##################################################
INApestSummary = matrix(ncol = 7,nrow = length(CommProbs)*length(SpreadReductions)*length(EradicationProbs)*length(DetectionProbs))
colnames(INApestSummary) = c("CommunicationProb","SpreadReduction","EradicationProb","DetectionProb","Invasion","Management","Detection")
counter = 0
for(commprob in 1:length(CommProbs))
{
CommunicationProb = CommProbs[commprob] #Proportion of nodes within threshold distance contacted each year
CommProbDir = paste0(RegionResultsDir,"CommProb_",CommunicationProb,"/")
#CrossRegionThreatDir = paste0(CommProbDir,"CrossRegionThreat/") 
#dir.create(CrossRegionThreatDir)
for(spreadreduction in 1:length(SpreadReductions))
for(erradprob in 1:length(EradicationProbs))
for(detprob in 1:length(DetectionProbs))
{
counter = counter+1
###Set spread reduction under management
SpreadReduction = SpreadReductions[spreadreduction]
###Management efficacy defined as annual Eradication probability
AnnualEradicationProb = EradicationProbs[erradprob]
DetectionProb = DetectionProbs[detprob]

FileNameStem = paste0(CommProbDir,"DetProb_",DetectionProb,"_erradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
ManagingResults = readRDS(paste0(FileNameStem,"InfoLargeOut.rds"))
InvasionResults = readRDS(paste0(FileNameStem,"InvasionLargeOut.rds"))
DetectedResults = readRDS(paste0(FileNameStem,"DetectedLargeOut.rds"))
InvasionProb = readRDS(paste0(FileNameStem,"InvasionProb.rds"))
Dim = as.array(dim(InvasionResults))
Invasion = sum(InvasionResults)/(Dim[1]*Dim[2]*Dim[3])
Managing = sum(ManagingResults)/(Dim[1]*Dim[2]*Dim[3])
Detected = sum(DetectedResults)/(Dim[1]*Dim[2]*Dim[3])
INApestSummary[counter,] = c(CommunicationProb,SpreadReduction,AnnualEradicationProb,DetectionProb,Invasion,Managing,Detected)
}
}
head(INApestSummary)
Filename = paste0(RegionResultsDir,"MultiManagementSummary_ClimateScenario_",cs,".csv")
write.table(INApestSummary,Filename,sep = ",",row.names = F)