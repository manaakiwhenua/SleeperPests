##################################################
###Simulate spread from infested farms HBAY as of 2021
###Apply management as 1000-fold reduction in dispersal 
###And local eradication following 
###https://www.landcareresearch.co.nz/uploads/public/Discover-Our-Research/Biosecurity/Biocontrol-ecology-of-weeds/3-applications/CNG-review-SFF-project-2010.pdf
###Calls a core function for fixed initial infestations
###Records farms infested, farms under management 
###and invasive threat to farms in neighbouring regions
###KEY UPDATE: this version implements information spread to
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
###Read in fixed start function invasion simulation function
###And cross region invasion risk function
###########################################################################
###########################################################################

source("INApest.R")
source("CrossRegionInvasionThreat.R")


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
CompleteDetectionProb = 0.25
AnnualEradicationProb = AnnualEradicationProb*CompleteDetectionProb 
#############################################################
#############################################################


ResultsDir= paste0(main.dir,"/HistoricExamplesLDDMatrixInfoSpread_v2/",ClimateScenarios[cs],"/")
dir.create(ResultsDir,recursive = T)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
region = 7
RegionResultsDir= paste0(ResultsDir,Regions[region],"/")
dir.create(RegionResultsDir)
CrossRegionThreatDir = paste0(RegionResultsDir,"CrossRegionThreat/") 
dir.create(CrossRegionThreatDir)

########################################################
###Read in simple XY points for INA
###Ecoclimatic index values and climate based extablishment for each farm
########################################################

data<-read_csv(paste0(ClimateScenarios[cs],"/Inputs/",Regions[region],"_Simple_points_for_INA.csv"))

#Points for plotting in INA
geocoords<-matrix(c(data$X, data$Y), byrow=F, ncol=2)

########################################################
###Calculate distance-based dispersal probability
#############################################################

###Read in distance matrix
dist_mat = readRDS(paste0(ClimateScenarios[cs],"/Inputs/",Regions[region],"_dist_mat_farm.rds"))

###Use custom distance kernel to obtain probabilities
library(actuar)
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
LDDweight = readRDS(paste0(LDDweightDir,Regions[region],"_",Regions[region],"_InvasionProb_LDDweight.rds"))
LDDmatrix = LDDweight*LongDistProbPerFarm
LDDprob = ppois(0,LDDmatrix,lower.tail = F)

##############################################################
###Assign climate based establishment prob based on ecoclimatic index
###Default function is logit
##############################################################
prob_est<-as.vector(data$Probability_Estab)
if(EI_Prob_CurveType == "SplitLinear")
	prob_est<-as.vector(data$Probability_Estab_SplitLinear)

####################################################################
###Read in historic invasion info
####################################################################
InputDir = r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CurrentFutureClimate\Inputs\]"
HistoricInfestation = read.delim(paste0(InputDir,"HBRC_CNG_points.txt"),header = T, as.is =T)
InfestedFarms = unique(HistoricInfestation$FarmID)


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
DetectionProbs = c(0,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75) 


###Declare initial infestation vector
###Start with infested farms from historic data
InitBio = rep(0,times = nrow(adj))
InitBio[FarmNames %in% InfestedFarms] = 1 
##################################################
###Call INApest function 
###Looping through different detection probabilities
##################################################
for(detprob in 1:length(DetectionProbs))
{
DetectionProb = DetectionProbs[detprob]
###Declare name for storing results
ModelName = paste0("DetProb_",DetectionProb,"_EradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
INApest(
ModelName = ModelName,		#Name for writing outputs to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Nyears = Nyears,                 #Simulation duration
DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be betwee 0 and 1)
ManageProb = ManageProb,         #Annual Probability of node adopting management upon detection
AnnualEradicationProb = AnnualEradicationProb, #Annual probability of eradication (must be betwee 0 and 1)
SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be betwee 0 and 1)
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic sdjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Option to provide long-distance (e.g. human-mediated) dispersal probability matrix
				#could be weighted by human visitation law or data on stock movements between farms
geocoords = geocoords,              #XY points for INAscene
OutputDir = RegionResultsDir			#Directory for storing results	
)


##################################################
###Use INApest outputs to calculate 
###Invasion threat to other regions
##################################################

ThreatSource = Regions[region]
for(sink in 1:length(Regions))
if(sink != region)
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