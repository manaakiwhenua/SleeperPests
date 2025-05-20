

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

###Declare initial infestation vector
###Start with infested farms from historic data
InitBio = rep(0,times = nrow(adj))
InitBio[FarmNames %in% HistoricInfestation$farm_id] = 1 




##################################################
###Call different versions of INApest function
###With zero management
###And different external Inv settings
###Check outputs to ensure different versions of INApest
###Give same results
##################################################

InitInvP = 5/nrow(adj)
set.seed(42)
InitialInv = rep(0,times = nrow(adj))
Sample <- sample(1:nrow(adj),5)
InitialInv[Sample] = 1

ExternalInvProb = rep(5/nrow(adj),times = nrow(adj))

ExternalInvProbMatrix = matrix(ncol = Ntimesteps, nrow  = nrow(adj))
for(i in 1:20)
  ExternalInvProbMatrix[,i] = ExternalInvProb
for(i in 21:Ntimesteps)
  ExternalInvProbMatrix[,i] = ExternalInvProb*10

DetectionProb <- 0
if(Function == "INApest")
{
######################################  
###Initial Inv as binary vector
######################################
  
###Declare name for storing results
set.seed(42)
ModelName = paste0("INApestInitialInv","DetProb_",DetectionProb,"_")
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
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = 0,                  ##Socioeconomic adjaceny matrix - Inv transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  OutputDir = OutputDir	#Directory for storing results to disk	
)

######################################  
###Initial Inv as proportion of nodes
######################################

ModelName = paste0("INApestInitInvP","DetProb_",DetectionProb,"_")
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
  InitialInvasion = NA,        #Nodes infested at start of simulations
  InitBioP = sum(InitBio)/length(InitBio),
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = 0,                  ##Socioeconomic adjaceny matrix - Inv transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  OutputDir = OutputDir	#Directory for storing results to disk	
)

################################
###Initial Inv as vector of probabilities
################################
ModelName = paste0("INApestExternalInvProb","DetProb_",DetectionProb,"_")
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
  InvasionRisk = ExternalInvProb,
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = 0,                  ##Socioeconomic adjaceny matrix - Inv transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  OutputDir = OutputDir	#Directory for storing results to disk	
)

################################
###Initial Inv as vector of probabilities
###and proportion of nodes
################################
ModelName = paste0("INApestExternalInvProbInitInvP","DetProb_",DetectionProb,"_")
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
  InitBioP = InitInvP,
  InvasionRisk = ExternalInvProb,
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = 0,                  ##Socioeconomic adjaceny matrix - Inv transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  OutputDir = OutputDir	#Directory for storing results to disk	
)

################################
###Initial Inv as vector of probabilities
###with ongoing external invasion
################################
ModelName = paste0("INApestOngoingExternalInv","DetProb_",DetectionProb,"_")
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
  InvasionRisk = ExternalInvProb,
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = 0,                  ##Socioeconomic adjaceny matrix - Inv transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  OngoingExternalInvasion = T,
  OutputDir = OutputDir	#Directory for storing results to disk	
)

################################
###Initial Inv as matrix of probabilities
###with ongoing external communication
################################
ModelName = paste0("INApestOngoingExternalInvMatrix","DetProb_",DetectionProb,"_")
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
  InvasionRisk =  ExternalInvProbMatrix,
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = 0,                  ##Socioeconomic adjaceny matrix - Inv transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  OngoingExternalInvasion = T,
  OutputDir = OutputDir	#Directory for storing results to disk	
)
}

DetectionProb <- 0
if(Function == "INApestParallel")
{
  ######################################  
  ###Initial Inv as binary vector
  ######################################
  
  ###Declare name for storing results
  set.seed(42)
  ModelName = paste0("INApestParallelInitialInv","DetProb_",DetectionProb,"_")
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
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - Inv transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ######################################  
  ###Initial Inv as proportion of nodes
  ######################################
  
  ModelName = paste0("INApestParallelInitInvP","DetProb_",DetectionProb,"_")
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
    InitialInvasion = NA,        #Nodes infested at start of simulations
    InitBioP = sum(InitBio)/length(InitBio),
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - Inv transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial Inv as vector of probabilities
  ################################
  ModelName = paste0("INApestParallelExternalInvProb","DetProb_",DetectionProb,"_")
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
    InvasionRisk = ExternalInvProb,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - Inv transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial Inv as vector of probabilities
  ###and proportion of nodes
  ################################
  ModelName = paste0("INApestParallelExternalInvProbInitInvP","DetProb_",DetectionProb,"_")
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
    InitBioP = InitInvP,
    InvasionRisk = ExternalInvProb,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - Inv transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial Inv as vector of probabilities
  ###with ongoing external invasion
  ################################
  ModelName = paste0("INApestParallelOngoingExternalInv","DetProb_",DetectionProb,"_")
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
    InvasionRisk = ExternalInvProb,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - Inv transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    OngoingExternalInvasion = T,
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
  
  ################################
  ###Initial Inv as matrix of probabilities
  ###with ongoing external communication
  ################################
  ModelName = paste0("INApestParallelOngoingExternalInvMatrix","DetProb_",DetectionProb,"_")
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
    InvasionRisk =  ExternalInvProbMatrix,
    EnvEstabProb = prob_est,           #Environmentally determined establishment probability
    SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
    SEAM = 0,                  ##Socioeconomic adjaceny matrix - Inv transfer probability between farms
    LDDprob = LDDprob,         #Long distance dispersal probability matrix 
    OngoingExternalInvasion = T,
    OutputDir = OutputDir	#Directory for storing results to disk	
  )
}
