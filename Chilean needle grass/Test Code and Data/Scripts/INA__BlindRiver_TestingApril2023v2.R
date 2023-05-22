##################################################
###Simulate spread from infested farms in Blind River
###Initial invasion approximates infestation documented here:
###https://www.tandfonline.com/doi/abs/10.1080/0028825X.1989.10414122
###Testing new functionality:
###1) Detection, spread reduction and management probability provided as matrix (nodes x timesteps)
###2) Erradication probability provided as vector (timesteps) 
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
Ntimesteps = 80
Nseeds = vector(length = Ntimesteps)
N0 = 1500
for(i in 1:Ntimesteps)
 {
 Nseeds[i] = N0*(1-(SeedBankLoss+Germination))
 N0 = Nseeds[i] 
 }
###Find timestep where seed density <1 seed per hectare
###Consider this effective eradiation 
EradicationTimestep = length(Nseeds[Nseeds>1/10^4])
 
###Annual eradication prob to achieve 0.95 errodication
###after 14 (timestep when Nseeds <1 per hectare) timesteps
AnnualEradicationProb = 1-(1-0.95)^(1/EradicationTimestep)

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

ResultsDir= paste0(main.dir,"/ResultsApril2023v2/",ClimateScenarios[cs],"/")
dir.create(ResultsDir,recursive = T,showWarnings = F)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
region = 10
RegionResultsDir= paste0(ResultsDir,Regions[region],"/")
dir.create(RegionResultsDir,showWarnings = F)

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
thresh= 1/20 #1 in 20 timestep spread risk between adjacent farms
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
CommunicationRate = 0.25 #Proportion of nodes within threshold distance contacted each timestep
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

###Set mean annual number of incoming long distance events per farm timestep
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
Ntimesteps = 80

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

######################################################
###Generate hypothetical vector (nodes) for detection, eradication, management prob and spread reduction
###Weighted by initial risk 
######################################################
StandardDetection = 0.1
StandardManage = 0.1
StandardSpreadReduction = 0.1
StandardEradication = AnnualEradicationProb

DispProb =1-(1-adj)*(1-LDDprob)
BPAM = sweep(DispProb,2,prob_est,`*`)
diag(BPAM) = 1

Risk =  sweep(BPAM,1,InitBio,`*`)
FarmRisk = colSums(Risk)
FarmRisk = FarmRisk^(1/3.1)
FarmRisk = FarmRisk/sum(FarmRisk)*length(FarmRisk)
DetectionVector = StandardDetection*FarmRisk
hist(DetectionVector)
mean(DetectionVector)

FarmRisk = colSums(Risk)
FarmRisk = FarmRisk^(1/3.1)
FarmRisk = FarmRisk/sum(FarmRisk)*length(FarmRisk)
hist(FarmRisk)
ManageVector = StandardManage*FarmRisk
hist(ManageVector)
max(ManageVector)
mean(ManageVector)

FarmRisk = colSums(Risk)
FarmRisk = FarmRisk^(1/3.1)
FarmRisk = FarmRisk/sum(FarmRisk)*length(FarmRisk)
hist(FarmRisk)
SpreadReductionVector = StandardSpreadReduction*FarmRisk
hist(SpreadReductionVector)
max(SpreadReductionVector)
mean(SpreadReductionVector)

FarmRisk = colSums(Risk)
FarmRisk = FarmRisk^(1/2)
FarmRisk = FarmRisk/sum(FarmRisk)*length(FarmRisk)
hist(FarmRisk)
EradicationVector = StandardEradication*FarmRisk
hist(EradicationVector)
max(EradicationVector)
mean(EradicationVector)

#######################################################
###Generate hypothetical climate time series with increasing 
###Prob est. Try 1% annual increase
###Would be possible to explore extreme timesteps where 
###conditions particularly favourable for establishment
#######################################################
prob_est_matrix = matrix(nrow = nrow(adj), ncol = Ntimesteps)
AnnulIncrease = 0.01
prob_est_matrix[,1] = prob_est
for(timestep in 2:Ntimesteps)
{
prob_est_matrix[,timestep] = prob_est_matrix[,(timestep-1)]*(1+AnnulIncrease*(timestep-1))
prob_est_matrix[,timestep] = ifelse(prob_est_matrix[,timestep]<1,prob_est_matrix[,timestep],1)
}
#plot(1:Ntimesteps,colMeans(prob_est_matrix)) 

#######################################################
###Generate hypothetical detection and management probability time series
###Try 1% annual increase
###May reflect awareness raising scenario
#######################################################
DetectionMatrix = matrix(nrow = nrow(adj), ncol = Ntimesteps)
AnnualIncrease = 0.01
DetectionMatrix[,1] = rnorm(nrow(adj),0.025,0.025/10)
for(timestep in 2:Ntimesteps)
{
DetectionMatrix[,timestep] = DetectionMatrix[,(timestep-1)]+AnnualIncrease
}
plot(1:Ntimesteps,colMeans(DetectionMatrix)) 

ManagementMatrix = matrix(nrow = nrow(adj), ncol = Ntimesteps)
AnnualIncrease = 0.01
ManagementMatrix[,1] = rnorm(nrow(adj),0.1,0.1/10)
for(timestep in 2:Ntimesteps)
{
ManagementMatrix[,timestep] = ManagementMatrix[,(timestep-1)]+AnnualIncrease
}
plot(1:Ntimesteps,colMeans(ManagementMatrix)) 

EradicationMatrix = matrix(nrow = nrow(adj), ncol = Ntimesteps)
AnnualIncrease = 0.01
EradicationMatrix[,1] = rnorm(nrow(adj),StandardEradication,StandardEradication/10)
for(timestep in 2:Ntimesteps)
{
EradicationMatrix[,timestep] = EradicationMatrix[,(timestep-1)]+AnnualIncrease
}
plot(1:Ntimesteps,colMeans(EradicationMatrix))

#######################################################
###Generate hypothetical spread reduction timeseries
###Try logistic growth curve with R = 0.5
###May reflect rapid increase in hygeine of movement restricton
#######################################################

N0 = 10
Nt1 = N0
K = 100
R = 0.5
SR = Nt1/100
for(timesteps in 2:80)
{
Nt2 = Nt1+R*Nt1*((K-Nt1)/K)
SR = c(SR,Nt2/100)
Nt1 = Nt2
}
plot(1:Ntimesteps,SR)
SpreadReductionMatrix = matrix(nrow = nrow(adj), ncol = Ntimesteps)
AnnualIncrease = 0.01
SpreadReductionMatrix[,1] = SR[1]
for(timestep in 2:Ntimesteps)
{
SpreadReductionMatrix[,timestep] = SR[timestep]
}
plot(1:Ntimesteps,colMeans(SpreadReductionMatrix)) 


##################################################
###Call INApest function 
###Comparing different trials with a standad model
##################################################

###Standard settings
###Declare name for storing results
ModelName = paste0("Standard")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
AnnualEradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = OutputDir	#Directory for storing results to disk	
)

###Vector (nodes) of annual eradication prob weighted by inital risk
###Declare name for storing results
ModelName = paste0("EradicationVector")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
AnnualEradicationProb = EradicationVector, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = OutputDir	#Directory for storing results to disk	
)

###Vector (nodes) of detection prob weighted by inital risk
ModelName = paste0("DetectionVector")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = DetectionVector,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
AnnualEradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = OutputDir	#Directory for storing results to disk	
)

###Vector (nodes) of manage prob weighted by inital risk
ModelName = paste0("ManageVector")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManageVector,         
ManageSD = NULL,
AnnualEradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = OutputDir	#Directory for storing results to disk	
)

###Vector (nodes) of spread reduction weighted by inital risk
ModelName = paste0("SpreadReductionVector")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
AnnualEradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = SpreadReductionVector,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = OutputDir	#Directory for storing results to disk	
)


###Vector (nodes) of detection, management, eradication and spread reduction weighted by inital risk
###Prioritises all spatially varying managemet variables by risk
ModelName = paste0("AllManagementPrioritised")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = DetectionVector,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManageVector,         
ManageSD = NULL,
AnnualEradicationProb = EradicationVector, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = SpreadReductionVector,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = OutputDir	#Directory for storing results to disk	
)


###Matrix (nodes x timesteps) of annual eradication prob randomly allocated to nodes, 
###increasing over timeme
###Declare name for storing results
ModelName = paste0("EradicationMatrix")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
AnnualEradicationProb = EradicationMatrix, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = OutputDir	#Directory for storing results to disk	
)



###Matrix (nodes x timesteps) of annually increasing detection probs 
ModelName = paste0("DetectionMatrix")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = DetectionMatrix,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
AnnualEradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = OutputDir	#Directory for storing results to disk	
)

###Matrix (nodes x timesteps) of annually increasing management probs
###Declare name for storing results
ModelName = paste0("ManagementMatrix")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManagementMatrix,         
ManageSD = NULL,
AnnualEradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = OutputDir	#Directory for storing results to disk	
)

###Spread reduction matrix
###Declare name for storing results
ModelName = paste0("SpreadReductionMatrix")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
AnnualEradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = SpreadReductionMatrix,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = OutputDir	#Directory for storing results to disk	
)

###All management variables improving over time
###Declare name for storing results
ModelName = paste0("AllManagementIncreasing")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApest(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = DetectionMatrix,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManagementMatrix,         
ManageSD = NULL,
AnnualEradicationProb = EradicationMatrix, #Annual probability of eradication (must be between 0 and 1)
AnnualEradicationSD = NULL,
SpreadReduction = SpreadReductionMatrix,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = OutputDir	#Directory for storing results to disk	
)


