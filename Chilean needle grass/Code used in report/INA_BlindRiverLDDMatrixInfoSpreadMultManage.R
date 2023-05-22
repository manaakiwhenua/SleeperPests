##################################################
###Simulate spread from infested farms in Blind River
###Initial invasion approximates infestation documented here:
###https://www.tandfonline.com/doi/abs/10.1080/0028825X.1989.10414122
###Explores multiple combinations of the following management variabels:
###Detection probability
###Eradication probability
###Spread reduction
###Communication probability (i.e. annual probability of contacting farms within threshold distance of known extant infestations)
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

source("INApestFixedStart.R")
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
###Consider this effective erradiation 
EradicationYear = length(Nseeds[Nseeds>1/10^4])
 
###Annual Eradication prob to achieve 0.95 errodication
###after 14 (year when Nseeds <1 per hectare) years
AnnualEradicationProb = 1-(1-0.95)^(1/EradicationYear)

###Moderate eradication prob by probability of identifying all patches on farm
###Value derived from survey data indicating ~50% of farmers are confident in identifying CNG
###And allowing a 50% chance of detecting all patches for those farmers
CompleteDetectionProb = 0.25
AnnualEradicationProb = AnnualEradicationProb*CompleteDetectionProb 

#############################################################
#############################################################


ResultsDir= paste0(main.dir,"/HistoricExamplesLDDMatrixInfoSpread/",ClimateScenarios[cs],"/")
dir.create(ResultsDir,recursive = T)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
region = 10
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
###and set long distance dispersal rate per farm
#############################################################

###Set mean annual number of incoming long distance events per farm year
###This ist the number of events is expected in a fully stocked matrix
###i.e all farms infested. Number of events will increase as more farms are infested
LongDistProbPerFarm = 0.05

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


######################################################
######################################################
###Use NZ cattle movement data to build a
###Long distance dispersal matrix
###Sourced from https://www.mpi.govt.nz/dmsdocument/49114-Analysis-of-New-Zealand-Stock-Movement-Data
###Figure 13a
######################################################
######################################################

LDDweightDir = paste0(main.dir,paste0("/",ClimateScenarios[cs],"/Inputs/CrossRegionDistance/"))
LDDweight = readRDS(paste0(LDDweightDir,Regions[region],"_",Regions[region],"_InvasionProb_LDDweight.rds"))

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

###Set spread reduction under management
###Assume 1000-fold reduction in spread prob under management
SpreadReduction = 0.999

###Management efficacy defined as annual Eradication probability
ManEfficacy = AnnualEradicationProb

###Set Eradication prob of cummunication with neighbours of detected infestations
CommProbs = c(0,0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.75)


###Set Eradication prob and spread reduction under management
SpreadReductions = c(0.05,0.1,0.3,0.5,0.9,0.99,0.999)
EradicationProbs = c(0,0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.75)



###Set probability of detection
###50% farmers confident of CNG ID 
###Default assumes those with ID knowledge have 50% chance of detection annually
###Declare vector to explore other probs
DetectionProbs = c(0,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75) 


###Declare initial infestation vector
###Start with infested farms from historic data
InitBio = rep(0,times = nrow(adj))
InitBio[FarmNames %in% HistoricInfestation$farm_id] = 1 
##################################################
###Call INApest function 
###Looping through different detection probabilities
##################################################
for(commprob in 1:length(CommProbs))
{
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
CommunicationProb = CommProbs[commprob] #Proportion of nodes within threshold distance contacted each year
SEAM = dist_mat
SEAM = ifelse(dist_mat < InfoMaxDist,CommunicationProb,0)
diag(SEAM) = 1
CommProbDir = paste0(RegionResultsDir,"CommProb_",CommunicationProb,"/")
dir.create(CommProbDir)

CrossRegionThreatDir = paste0(CommProbDir,"CrossRegionThreat/") 
dir.create(CrossRegionThreatDir)

for(spreadreduction in 1:length(SpreadReductions))
for(erradprob in 1:length(EradicationProbs))
for(detprob in 1:length(DetectionProbs))
{
###Set spread reduction under management
SpreadReduction = SpreadReductions[spreadreduction]
###Management efficacy defined as annual Eradication probability
AnnualEradicationProb = EradicationProbs[erradprob]
ManEfficacy = AnnualEradicationProb


DetectionProb = DetectionProbs[detprob]
INApestFixedStart(
Nperm = Nperm,                  #Number of permutations per paramteter combination
Nyears = Nyears,                 #Simulation duration
DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be betwee 0 and 1)
ManageProb = ManageProb,         #Annual Probability of node adopting management upon detection
AnnualEradicationProb = AnnualEradicationProb, #Annual probability of Eradication (must be betwee 0 and 1)
SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be betwee 0 and 1)
LDDrate = LongDistProbPerFarm,            #Mean annual per-farm distant-independent dispersal rate. Value reflects rate for fully invaded network and zero
                        #environmental limitation on establishment
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
BPAM = adj,                   #Biophysical sdjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic sdjaceny matrix - info transfer probability between farms
LDDmatrix = LDDweight*LongDistProbPerFarm,         #Option to provide long distance dispersal matrix instead of distance-independent dispesal rate
				#e.g. could be weighted by human visitation law or data on stock movements between farms
geocoords = geocoords,              #XY points for INAscene
OutputDir = CommProbDir			#Directory for storing results	
)


##################################################
###Use INApest outputs to calculate 
###Invasion threat to other regions
##################################################

ThreatSource = Regions[region]
for(sink in 1:length(Regions))
if(sink != region && 1==2)
{
ThreatSink = Regions[sink]
InputDir = paste0(main.dir,"/",ClimateScenarios[cs],"/Inputs/")
OutputDir = CrossRegionThreatDir
dir.create(OutputDir)
CrossRegionInvasionDir = paste0(main.dir,paste0("/",ClimateScenarios[cs],"/Inputs/CrossRegionDistance/"))
InvasionFileNameStem = paste0(RegionResultsDir,"DetProb_",DetectionProb,"_ErradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
OutputFileNameStem = paste0(OutputDir,"DetProb_",DetectionProb,"_ErradProb_",round(AnnualEradicationProb,digits =2),
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
}
##################################################
###End of detection prob loop
##################################################