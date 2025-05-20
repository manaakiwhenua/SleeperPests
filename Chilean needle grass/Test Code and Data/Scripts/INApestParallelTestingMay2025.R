##################################################
###Simulate spread from infested farms in Blind River
###Initial invasion approximates infestation documented here:
###https://www.tandfonline.com/doi/abs/10.1080/0028825X.1989.10414122
###Testing new functionality:
### 1) Parallel processing of permutation loop
### 2) Manually implementing steps formerly performed by INAscene
###Benchmarking tests below indicate ~3 fold speed increase from parallel processing
###with 8 cores for Marlborough
###Removal of INAscene call doesn't significantly change processing time
##################################################





####################################################
###Input paramaters for INApest
####################################################
###Set number of realisations
Nperm = 100
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
###And detection prob
###Check outputs to ensure different versions of INApest
###Give same results
##################################################

for(detprob in 1:length(DetectionProbs))
{
DetectionProb <- DetectionProbs[detprob]
###INApest
###Declare name for storing results
set.seed(42)
ModelName = paste0("INApest","DetProb_",DetectionProb,"_")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
Start <- Sys.time()
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
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  OutputDir = OutputDir	#Directory for storing results to disk	
)
End <- Sys.time()
SerialTime <- End-Start

###INApestParallel
###Declare name for storing results
set.seed(42)
ModelName = paste0("INApestParallel","DetProb_",DetectionProb,"_")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
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
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  OutputDir = OutputDir	#Directory for storing results to disk	
)
End <- Sys.time()
ParallelTime <- End-Start


###INApestINAscene
###Declare name for storing results
set.seed(42)
ModelName = paste0("INApestINAscene","DetProb_",DetectionProb,"_")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
Start <- Sys.time()
INApestINAscene(
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
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix
  geocoords = geocoords,              #XY points for INAscene
  OutputDir = OutputDir	#Directory for storing results to disk	
)
End <- Sys.time()
SerialTimeINAscene <- End-Start

###INApestParallelINAscene
###Declare name for storing results
set.seed(42)
ModelName = paste0("INApestParallelINAscene","DetProb_",DetectionProb,"_")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
Start <- Sys.time()
INApestParallelINAscene(
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
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
geocoords = geocoords,              #XY points for INAscene
OutputDir = OutputDir	#Directory for storing results to disk	
)
End <- Sys.time()
ParallelTimeINAscene <- End-Start
}
