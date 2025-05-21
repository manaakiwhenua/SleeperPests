##################################################
###Simulate spread in sink region with constant invasion threat from the source region
###This version uses results from a zero management scenario (detection prob = 0)
###In the source region
##################################################


library(INA)
library(tidyverse)
library(dplyr)
library(actuar)
memory.limit(size=600000) 

###########################################################################
###########################################################################
###Read in function from source file
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

ResultsDir= paste0(main.dir,"/CrossRegionConstantThreatv2/",ClimateScenarios[cs],"/")
dir.create(ResultsDir,recursive = T,showWarnings = F)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
RegionResultsDir= paste0(ResultsDir,SinkRegion,"_From_",SourceRegion,"/")
dir.create(RegionResultsDir,showWarnings = F)
CrossRegionThreatDir = paste0(RegionResultsDir,"CrossRegionThreat/",showWarnings = F) 
dir.create(CrossRegionThreatDir)

########################################################
###Read in simple XY points for INA
###And Ecoclimatic index values for each farm
########################################################

data<-read_csv(paste0(ClimateScenarios[cs],"/Inputs/",SinkRegion,"_Simple_points_for_INA.csv"))

#Points for plotting in INA
geocoords<-matrix(c(data$X, data$Y), byrow=F, ncol=2)

########################################################
###Calculate distance-based dispersal probability
#############################################################

###Read in distance matrix
dist_mat = readRDS(paste0(ClimateScenarios[cs],"/Inputs/",SinkRegion,"_dist_mat_farm.rds"))

###Use custom distance kernel to obtain probabilities
k= 7 #shape parameter
b= 500 #scale parameter
thresh= 1/20 #1 in 20 year spread risk between adjacent farms
adj = thresh*actuar::ppareto(dist_mat, shape=k,  scale=b, lower.tail=F,log.p = F)

###Store names for nodes in network
FarmNames = row.names(adj)

########################################################
###Define a maximum threshold distance for informaton spread
###from known infestations
###and communication rate (by authorities and owners of infested farms) to farms within 
###the threshold distance
###Set socioeconomic adjacency matrix (SEAM) so that nodes within threshold distance
###are assigned a fixed probability of info transfer equal to the communication rate
###set all other SEAM values to zero
#############################################################

InfoMaxDist = 150 #Info trnasfer threshold distance in metres
                  #Can be set to maximum annual (non-human aided) dispersal distance of pest
CommunicationRate = 0.25 #Proportion of nodes within threshold distance contacted each year
SEAM = dist_mat
SEAM = ifelse(dist_mat < InfoMaxDist,CommunicationRate,0)
diag(SEAM) = 1

######################################################
######################################################
###Use NZ cattle movement data to build a
###Long distance dispersal matrix
###Sourced from https://www.mpi.govt.nz/dmsdocument/49114-Analysis-of-New-Zealand-Stock-Movement-Data
###Figure 13a
######################################################
######################################################

###Set mean annual number of incoming long distance events per farm year
###This is the number of events is expected in a fully stocked matrix
###i.e all farms infested. Number of events will increase as more farms are infested
LongDistProbPerFarm = 0.05
LDDweightDir = paste0(main.dir,paste0("/",ClimateScenarios[cs],"/Inputs/CrossRegionDistance/"))
LDDweight = readRDS(paste0(LDDweightDir,SinkRegion,"_",SinkRegion,"_InvasionProb_LDDweight.rds"))
LDDmatrix = LDDweight*LongDistProbPerFarm
LDDprob = ppois(0,LDDmatrix,lower.tail = F)



##############################################################
###Assign climate based establishment prob based on ecoclimatic index
###Default function is logit
##############################################################
prob_est<-as.vector(data$Probability_Estab)
if(EI_Prob_CurveType == "SplitLinear")
	prob_est<-as.vector(data$Probability_Estab_SplitLinear)

##############################################
###Read in farm-level cross-region invasion risk
###based on simulation results from source region
###Load for all years
##############################################
SourceResultsDir = r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CurrentFutureClimate\HistoricExamplesLDDMatrixInfoSpread_v2\]"
SourceThreatDir = paste0(SourceResultsDir,ClimateScenarios[cs],"/",SourceRegion,"/CrossRegionThreat/")

CrossRegionThreatFile = paste0(SourceThreatDir,"DetProb_0.2_EradProb_0.05_SpreadReduction_0.999_",SourceRegion,"_",SinkRegion,"FarmRisk.rds") 
ExternalFarmRisk = readRDS(CrossRegionThreatFile)

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

###Set probability of detection
###50% farmers confident of CNG ID 
###Default assumes those with ID knowledge have 50% chance of detection annually
###i.e. detection prob = 0.25
###Declare vector to explore other detection probabilities
#DetectionProbs = c(0,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75) 



##################################################
###Call INApest function 
###Looping through different detection probabilities
###Function incorporates management variables in filenames of outputs
##################################################

for(detprob in 1:length(DetectionProbs))
{
DetectionProb = DetectionProbs[detprob]
###Declare name for storing results
ModelName = paste0("DetProb_",DetectionProb,"_eradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
INApest(
ModelName = ModelName,          #Name for storing results to disk 
Nperm = Nperm,                  #Number of permutations per paramteter combination
Nyears = Nyears,                 #Simulation duration
DetectionProb = DetectionProb,          #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
ManageProb = ManageProb,             #Annual Probability of node adopting management upon detection
AnnualEradicationProb = AnnualEradicationProb, #Annual probability of eradication (must be between 0 and 1)
SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
InitBioP = NA,			#Proportion of nodes infested at start of simulations
InvasionRisk = ExternalFarmRisk,           #Vector of probabilities for weighting random assignment of initial invasion occurrences
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #distance-based disperal probability
SEAM = SEAM,			#Option to provide socioeconomic adjacency matrix
LDDprob = LDDprob,              #Option to provide long distance dispersal probability matrix
				#e.g. could be weighted by human visitation law or data on stock movements
OngoingExternal = T,    #Allow ongoing invasion from external sources
geocoords = geocoords,              #XY points for INAscene
OutputDir = RegionResultsDir			#Directory for storing results	
)

##################################################
###Use INApest outputs to calculate 
###Invasion threat to other regions
##################################################
if(1 == 2)
{
ThreatSource = SinkRegion
for(sink in 1:length(Regions))
if(Regions[sink] != ThreatSource)
{
ThreatSink = Regions[sink]
InputDir = paste0(main.dir,"/",ClimateScenarios[cs],"/Inputs/")
OutputDir = CrossRegionThreatDir
dir.create(OutputDir)
CrossRegionInvasionDir = paste0(main.dir,paste0("/",ClimateScenarios[cs],"/Inputs/CrossRegionDistance/"))
InvasionFileNameStem = paste0(RegionResultsDir,"DetProb_",DetectionProb,"_EradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
OutputFileNameStem = paste0(OutputDir,"DP_",DetectionProb,"_EP_",round(AnnualEradicationProb,digits =2),
		"_SR_",SpreadReduction,"_")
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
}
##################################################
##################################################

}

##################################################
###End of detection prob loop
##################################################