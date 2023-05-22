##################################################
###Simulate spread from infested farms in Blind River
###Initial invasion approximates infestation documented here:
###https://www.tandfonline.com/doi/abs/10.1080/0028825X.1989.10414122
###Testing new functionality:
###1) User-defined stadard deviation for management variables
###2) EnvEstab prob provided as a matrix rather than vector 
###   Applicable for changes thruugh time in factors affecting establishment e.g. climate change timeseries)
###3) Local extinction probability, provided as single number, vector (nodes) or matrix (nodes x years)
##################################################

library(INA)
library(tidyverse)
library(dplyr)
library(boot)
memory.limit(size=600000) 






###########################################################################
###########################################################################
###Read in the function from source file
###########################################################################
###########################################################################

source("Scripts/INApest.R")

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

ResultsDir= paste0(main.dir,"/INA_HistoricTesting/",ClimateScenarios[cs],"/")
dir.create(ResultsDir,recursive = T)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
region = 10
RegionResultsDir= paste0(ResultsDir,Regions[region],"/")
dir.create(RegionResultsDir)
CrossRegionThreatDir = paste0(RegionResultsDir,"CrossRegionThreat/") 
dir.create(CrossRegionThreatDir)

########################################################
###Read in simple XY points for INA and
###Ecoclimatic index values and climate based extablishment probability for each farm
########################################################

data<-read_csv(paste0(ClimateScenarios[cs],"/",Regions[region],"_Simple_points_for_INA.csv"))

#Points for plotting in INA
geocoords<-matrix(c(data$X, data$Y), byrow=F, ncol=2)

########################################################
###Calculate distance-based dispersal probability between farms
#############################################################


###Read in distance matrix
dist_mat = readRDS(paste0(ClimateScenarios[cs],"/",Regions[region],"_dist_mat_farm.rds"))

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

InfoMaxDist = 150 #Info transfer threshold distance in metres
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
###Produced using following source files:
###"CrossRegion_FarmDistanceInvasionProbLDD.R"
###"CrossRegion_LDDWeights.R"
######################################################
######################################################

###Set mean annual number of incoming long distance events per farm year
###This is the number of events is expected in a fully stocked matrix
###i.e all farms infested. Number of events will increase as more farms are infested
LongDistProbPerFarm = 0.05
LDDweightDir = paste0(main.dir,paste0("/",ClimateScenarios[cs],"/Inputs/CrossRegionDistance/"))
LDDweight = readRDS(paste0(ClimateScenarios[cs],"/",Regions[region],"_",Regions[region],"_InvasionProb_LDDweight.rds"))
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
###Uses Farm ID from Agribase 2022 as property identifier
###Farms selected approximate area infested as presented in:
###https://www.tandfonline.com/doi/abs/10.1080/0028825X.1989.10414122
####################################################################
InputDir = r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CurrentFutureClimate\Inputs\]"
HistoricInfestation = read.csv("BlindRiverInfestedFarms.txt",header = T, as.is =T)


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

###Set probability of detection
###50% farmers confident of CNG ID 
###Default assumes those with ID knowledge have 50% chance of detection annually
###Declare vector to explore other probs
DetectionProbs = c(0,0.05,0.2) 


###Declare initial infestation vector
###Start with infested farms from historic data
InitBio = rep(0,times = nrow(adj))
InitBio[FarmNames %in% HistoricInfestation$farm_id] = 1 

#######################################################
###Generate hypothetical climate time series with increasing 
###Prob est. Try 1% annual increase
###Would be possible to explore extreme years where 
###conditions particularly favourable for establishment
#######################################################
prob_est_matrix = matrix(nrow = nrow(adj), ncol = Nyears)
AnnulIncrease = 0.01
prob_est_matrix[,1] = prob_est
for(year in 2:Nyears)
{
prob_est_matrix[,year] = prob_est_matrix[,(year-1)]*(1+AnnulIncrease*(year-1))
prob_est_matrix[,year] = ifelse(prob_est_matrix[,year]<1,prob_est_matrix[,year],1)
}
#plot(1:Nyears,colMeans(prob_est_matrix)) 

#######################################################
###Generate hypothetical climate time series with decreasing
###local extinction. Try 5% annual decrease
###Would be possible to explore extreme years where 
###conditions drive particularly high or low local extinction
#######################################################
ExtinctionMatrix = matrix(nrow = nrow(adj), ncol = Nyears)
AnnulDecrease = 0.01
ExtinctionMatrix[,1] = runif(nrow(adj),0,0.2)
for(year in 2:Nyears)
{
ExtinctionMatrix[,year] = ExtinctionMatrix[,(year-1)]*(1-AnnulDecrease*(year-1))
}
#plot(1:Nyears,colMeans(ExtinctionMatrix)) 



##################################################
###Call INApest function 
###Looping through different detection probabilities
##################################################
for(detprob in 1:length(DetectionProbs))
{
DetectionProb = DetectionProbs[detprob]
###Standard settings
###Declare name for storing results
ModelName = paste0("Standard","DetProb_",DetectionProb,"_eradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Nyears = Nyears,                 #Simulation duration
DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManageProb,         
ManageSD = NULL,
AnnualEradicationProb = AnnualEradicationProb, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = RegionResultsDir	#Directory for storing results to disk	
)


###Add in SD for management variables
###Declare name for storing results
ModelName = paste0("ManSD_User","DetProb_",DetectionProb,"_eradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Nyears = Nyears,                 #Simulation duration
DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = runif(nrow(adj), 0.01,0.05), #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManageProb,         
ManageSD = runif(nrow(adj), 0.01,0.05),
AnnualEradicationProb = AnnualEradicationProb, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = runif(1, 0.01,0.05),
SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = runif(nrow(adj), 0.01,0.05),
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = RegionResultsDir	#Directory for storing results to disk	
)

###Provide EnvEstabProb as matrix
ModelName = paste0("ProbEstMatrix","DetProb_",DetectionProb,"_eradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Nyears = Nyears,                 #Simulation duration
DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManageProb,         
ManageSD = NULL,
AnnualEradicationProb = AnnualEradicationProb, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est_matrix,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = RegionResultsDir	#Directory for storing results to disk	
)


ModelName = paste0("LocalExtinctionFixed","DetProb_",DetectionProb,"_eradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Nyears = Nyears,                 #Simulation duration
DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManageProb,         
ManageSD = NULL,
AnnualEradicationProb = AnnualEradicationProb, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
LocalExtinction = ExtinctionMatrix[,1],     #Annual local population extinction probability. Set to zero for no local extinction. Can be vector (nodes) or matrix (nodes x years)
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = RegionResultsDir	#Directory for storing results to disk	
)

ModelName = paste0("LocalExtinctionDeclining","DetProb_",DetectionProb,"_eradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Nyears = Nyears,                 #Simulation duration
DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManageProb,         
ManageSD = NULL,
AnnualEradicationProb = AnnualEradicationProb, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
LocalExtinction = ExtinctionMatrix,     #Annual local population extinction probability. Set to zero for no local extinction. Can be vector (nodes) or matrix (nodes x years)
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = RegionResultsDir	#Directory for storing results to disk	
)


##################################################
##################################################

}
##################################################
###End of detection prob loop
##################################################