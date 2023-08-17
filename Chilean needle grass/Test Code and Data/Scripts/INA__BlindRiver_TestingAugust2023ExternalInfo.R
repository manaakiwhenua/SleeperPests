##################################################
###Simulate spread from infested farms in Blind River
###Initial invasion approximates infestation documented here:
###https://www.tandfonline.com/doi/abs/10.1080/0028825X.1989.10414122
###Testing new functionality:
### 1) Communication of info from external sources
##################################################



#############################################################
###Define directories for storing results
#############################################################

ResultsDir= paste0(main.dir,"/ResultsAugust2023ExternalInfo/",ClimateScenarios[cs],"/")
dir.create(ResultsDir,recursive = T,showWarnings = F)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
region = 10
RegionResultsDir= paste0(ResultsDir,Regions[region],"/")
dir.create(RegionResultsDir,showWarnings = F)

###Read in farm shape file for plotting heat maps
FarmShapeFile <- paste0(main.dir,"/",ClimateScenarios[cs],"/","MARLSheepBeefAtRisk",ClimateScenarios[cs],".wgs84.shp")
FarmShape <-sf::st_read(FarmShapeFile)

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
Nperm = 10
###Set simulation duration
Ntimesteps = 80

###Probability of management adoption when infestation detected
ManageProb = 0.9

###Set spread reduction under management
###Assume 1000-fold reduction in spread prob under management
SpreadReduction = 0.999

###Set probability of detection
DetectionProbs = c(0,0.025) 


###Declare initial infestation vector
###Start with infested farms from historic data
InitBio = rep(0,times = nrow(adj))
InitBio[FarmNames %in% HistoricInfestation$farm_id] = 1 




##################################################
###Call different versions of INApest function
###With zero management
###And different external info settings
###Check outputs to ensure different versions of INApest
###Give same results
##################################################

InitInfoP = 5/nrow(adj)
set.seed(42)
InitialInfo = rep(0,times = nrow(adj))
Sample <- sample(1:nrow(adj),5)
InitialInfo[Sample] = 1

ExternalInfoProb = rep(5/nrow(adj),times = nrow(adj))

ExternalInfoProbMatrix = matrix(ncol = Ntimesteps, nrow  = nrow(adj))
for(i in 1:20)
  ExternalInfoProbMatrix[,i] = ExternalInfoProb
for(i in 21:Ntimesteps)
  ExternalInfoProbMatrix[,i] = ExternalInfoProb*10

DetectionProb <- 0
if(Function == "INApest")
{
######################################  
###Initial info as binary vector
######################################
  
###Declare name for storing results
set.seed(42)
ModelName = paste0("INApestInitialInfo","DetProb_",DetectionProb,"_")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
ModelName = ""
INApest(
  ModelName = ModelName,      #Name for storing results to disk
  Nperm = Nperm,                  #Number of permutations per paramteter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
  DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
  ManageProb = ManageProb,         
  ManageSD = NULL,
  EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
  EradicationSD = NULL,
  SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
  SpreadReductionSD = NULL,
  InitialInvasion = InitBio,        #Nodes infested at start of simulations
  InitialInfo = InitialInfo,
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  geocoords = geocoords,              #XY points for INAscene
  OutputDir = OutputDir	#Directory for storing results to disk	
)

######################################  
###Initial info as proportion of nodes
######################################

ModelName = paste0("INApestInitInfoP","DetProb_",DetectionProb,"_")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
ModelName = ""
INApest(
  ModelName = ModelName,      #Name for storing results to disk
  Nperm = Nperm,                  #Number of permutations per paramteter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
  DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
  ManageProb = ManageProb,         
  ManageSD = NULL,
  EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
  EradicationSD = NULL,
  SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
  SpreadReductionSD = NULL,
  InitialInvasion = InitBio,        #Nodes infested at start of simulations
  InitInfoP = InitInfoP,
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  geocoords = geocoords,              #XY points for INAscene
  OutputDir = OutputDir	#Directory for storing results to disk	
)

################################
###Initial info as vector of probabilities
################################
ModelName = paste0("INApestExternalInfoProb","DetProb_",DetectionProb,"_")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
ModelName = ""
INApest(
  ModelName = ModelName,      #Name for storing results to disk
  Nperm = Nperm,                  #Number of permutations per paramteter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
  DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
  ManageProb = ManageProb,         
  ManageSD = NULL,
  EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
  EradicationSD = NULL,
  SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
  SpreadReductionSD = NULL,
  InitialInvasion = InitBio,        #Nodes infested at start of simulations
  ExternalInfoProb = ExternalInfoProb,
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  geocoords = geocoords,              #XY points for INAscene
  OutputDir = OutputDir	#Directory for storing results to disk	
)

################################
###Initial info as vector of probabilities
###and proportion of nodes
################################
ModelName = paste0("INApestExternalInfoProbInitInfoP","DetProb_",DetectionProb,"_")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
ModelName = ""
INApest(
  ModelName = ModelName,      #Name for storing results to disk
  Nperm = Nperm,                  #Number of permutations per paramteter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
  DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
  ManageProb = ManageProb,         
  ManageSD = NULL,
  EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
  EradicationSD = NULL,
  SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
  SpreadReductionSD = NULL,
  InitialInvasion = InitBio,        #Nodes infested at start of simulations
  InitInfoP = InitInfoP,
  ExternalInfoProb = ExternalInfoProb,
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  geocoords = geocoords,              #XY points for INAscene
  OutputDir = OutputDir	#Directory for storing results to disk	
)

################################
###Initial info as vector of probabilities
###with ongoing external communication
################################
ModelName = paste0("INApestOngoingExternalInfo","DetProb_",DetectionProb,"_")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
ModelName = ""
INApest(
  ModelName = ModelName,      #Name for storing results to disk
  Nperm = Nperm,                  #Number of permutations per paramteter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
  DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
  ManageProb = ManageProb,         
  ManageSD = NULL,
  EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
  EradicationSD = NULL,
  SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
  SpreadReductionSD = NULL,
  InitialInvasion = InitBio,        #Nodes infested at start of simulations
  InitInfoP = InitInfoP,
  ExternalInfoProb = ExternalInfoProb,
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  geocoords = geocoords,              #XY points for INAscene
  OngoingExternalInfo = T,
  OutputDir = OutputDir	#Directory for storing results to disk	
)

################################
###Initial info as matrix of probabilities
###with ongoing external communication
################################
ModelName = paste0("INApestOngoingExternalInfoMatrix","DetProb_",DetectionProb,"_")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
ModelName = ""
INApest(
  ModelName = ModelName,      #Name for storing results to disk
  Nperm = Nperm,                  #Number of permutations per paramteter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
  DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
  ManageProb = ManageProb,         
  ManageSD = NULL,
  EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
  EradicationSD = NULL,
  SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
  SpreadReductionSD = NULL,
  InitialInvasion = InitBio,        #Nodes infested at start of simulations
  InitInfoP = InitInfoP,
  ExternalInfoProb = ExternalInfoProbMatrix,
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  geocoords = geocoords,              #XY points for INAscene
  OngoingExternalInfo = T,
  OutputDir = OutputDir	#Directory for storing results to disk	
)
}

DetectionProb <- 0
if(Function == "INApestParallel")
{
  
  
  ######################################  
  ###Initial info as binary vector
  ######################################
  
  ###Declare name for storing results
  set.seed(42)
  ModelName = paste0("INApestParallelInitialInfo","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  Start <- Sys.time()
  
  INApestParallel(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitialInfo = InitialInfo,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    geocoords = geocoords,              #XY points for INAscene
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ######################################  
  ###Initial info as proportion of nodes
  ######################################
  
  ModelName = paste0("INApestParallelInitInfoP","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  Start <- Sys.time()
  INApestParallel(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitInfoP = InitInfoP,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    geocoords = geocoords,              #XY points for INAscene
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial info as vector of probabilities
  ################################
  ModelName = paste0("INApestParallelExternalInfoProb","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestParallel(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    ExternalInfoProb = ExternalInfoProb,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    geocoords = geocoords,              #XY points for INAscene
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial info as vector of probabilities
  ###and proportion of nodes
  ################################
  ModelName = paste0("INApestParallelExternalInfoProbInitInfoP","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestParallel(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitInfoP = InitInfoP,
    ExternalInfoProb = ExternalInfoProb,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    geocoords = geocoords,              #XY points for INAscene
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial info as vector of probabilities
  ###with ongoing external communication
  ################################
  ModelName = paste0("INApestParallelOngoingExternalInfo","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestParallel(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitInfoP = InitInfoP,
    ExternalInfoProb = ExternalInfoProb,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    geocoords = geocoords,              #XY points for INAscene
    OngoingExternalInfo = T,
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial info as matrix of probabilities
  ###with ongoing external communication
  ################################
  ModelName = paste0("INApestParallelOngoingExternalInfoMatrix","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestParallel(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitInfoP = InitInfoP,
    ExternalInfoProb = ExternalInfoProbMatrix,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    geocoords = geocoords,              #XY points for INAscene
    OngoingExternalInfo = T,
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
}



DetectionProb <- 0
if(Function == "INApestNoINAscene")
{
  ######################################  
  ###Initial info as binary vector
  ######################################
  
  ###Declare name for storing results
  set.seed(42)
  ModelName = paste0("INApestNoINAsceneInitialInfo","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestNoINAscene(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitialInfo = InitialInfo,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ######################################  
  ###Initial info as proportion of nodes
  ######################################
  
  ModelName = paste0("INApestNoINAsceneInitInfoP","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestNoINAscene(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitInfoP = InitInfoP,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial info as vector of probabilities
  ################################
  ModelName = paste0("INApestNoINAsceneExternalInfoProb","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestNoINAscene(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    ExternalInfoProb = ExternalInfoProb,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial info as vector of probabilities
  ###and proportion of nodes
  ################################
  ModelName = paste0("INApestNoINAsceneExternalInfoProbInitInfoP","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestNoINAscene(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitInfoP = InitInfoP,
    ExternalInfoProb = ExternalInfoProb,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial info as vector of probabilities
  ###with ongoing external communication
  ################################
  ModelName = paste0("INApestNoINAsceneOngoingExternalInfo","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestNoINAscene(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitInfoP = InitInfoP,
    ExternalInfoProb = ExternalInfoProb,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    
    OngoingExternalInfo = T,
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial info as matrix of probabilities
  ###with ongoing external communication
  ################################
  ModelName = paste0("INApestNoINAsceneOngoingExternalInfoMatrix","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestNoINAscene(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitInfoP = InitInfoP,
    ExternalInfoProb = ExternalInfoProbMatrix,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    
    OngoingExternalInfo = T,
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
}

DetectionProb <- 0
if(Function == "INApestParallelNoINAscene")
{
  
  
  ######################################  
  ###Initial info as binary vector
  ######################################
  
  ###Declare name for storing results
  set.seed(42)
  ModelName = paste0("INApestParallelNoINAsceneInitialInfo","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  Start <- Sys.time()
  
  INApestParallelNoINAscene(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitialInfo = InitialInfo,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ######################################  
  ###Initial info as proportion of nodes
  ######################################
  
  ModelName = paste0("INApestParallelNoINAsceneInitInfoP","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  Start <- Sys.time()
  INApestParallelNoINAscene(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitInfoP = InitInfoP,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial info as vector of probabilities
  ################################
  ModelName = paste0("INApestParallelNoINAsceneExternalInfoProb","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestParallelNoINAscene(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    ExternalInfoProb = ExternalInfoProb,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial info as vector of probabilities
  ###and proportion of nodes
  ################################
  ModelName = paste0("INApestParallelNoINAsceneExternalInfoProbInitInfoP","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestParallelNoINAscene(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitInfoP = InitInfoP,
    ExternalInfoProb = ExternalInfoProb,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial info as vector of probabilities
  ###with ongoing external communication
  ################################
  ModelName = paste0("INApestParallelNoINAsceneOngoingExternalInfo","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestParallelNoINAscene(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitInfoP = InitInfoP,
    ExternalInfoProb = ExternalInfoProb,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    
    OngoingExternalInfo = T,
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial info as matrix of probabilities
  ###with ongoing external communication
  ################################
  ModelName = paste0("INApestParallelNoINAsceneOngoingExternalInfoMatrix","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
  INApestParallelNoINAscene(
    ModelName = ModelName,      #Name for storing results to disk
    Nperm = Nperm,                  #Number of permutations per paramteter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
    DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
    ManageProb = ManageProb,         
    ManageSD = NULL,
    EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
    EradicationSD = NULL,
    SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
    SpreadReductionSD = NULL,
    InitialInvasion = InitBio,        #Nodes infested at start of simulations
    InitInfoP = InitInfoP,
    ExternalInfoProb = ExternalInfoProbMatrix,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    
    OngoingExternalInfo = T,
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
}












