

####################################################
###Input paramaters for INA
####################################################
###Set number of realisations
Nperm = 50
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
    OngoingExternalInfo = T,
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
}



DetectionProb <- 0
if(Function == "INApestINAscene")
{
  ######################################  
  ###Initial info as binary vector
  ######################################
  
  ###Declare name for storing results
  set.seed(42)
  ModelName = paste0("INApestINAsceneInitialInfo","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
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
  
  ModelName = paste0("INApestINAsceneInitInfoP","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
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
  ModelName = paste0("INApestINAsceneExternalInfoProb","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
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
  ModelName = paste0("INApestINAsceneExternalInfoProbInitInfoP","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
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
  ModelName = paste0("INApestINAsceneOngoingExternalInfo","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
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
  ModelName = paste0("INApestINAsceneOngoingExternalInfoMatrix","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
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
if(Function == "INApestParallelINAscene")
{
  
  
  ######################################  
  ###Initial info as binary vector
  ######################################
  
  ###Declare name for storing results
  set.seed(42)
  ModelName = paste0("INApestParallelINAsceneInitialInfo","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
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
  
  ModelName = paste0("INApestParallelINAsceneInitInfoP","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
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
  ModelName = paste0("INApestParallelINAsceneExternalInfoProb","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
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
  ModelName = paste0("INApestParallelINAsceneExternalInfoProbInitInfoP","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
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
  ModelName = paste0("INApestParallelINAsceneOngoingExternalInfo","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
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
  ModelName = paste0("INApestParallelINAsceneOngoingExternalInfoMatrix","DetProb_",DetectionProb,"_")
  OutputDir = paste0(RegionResultsDir,ModelName,"/")
  dir.create(OutputDir,showWarnings = F)
  ModelName = ""
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












