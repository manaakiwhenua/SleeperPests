#####################################################################################
###Test external invasion functionality
#####################################################################################
Nperm = 30
Ntimesteps = 50
DetectionProb = 0.1

###Make this high for obvious effect
AnnualIncursionRate = 0.1

ExternalInvasionDir <- paste0(ResultsDir,"ExternalInvasion/")
dir.create(ExternalInvasionDir)

################################
###Detection prob 0.1 no external invasion
################################
InitialPop = rep(0,times = nrow(d))
Sample <- sample(1:nrow(d),10)
InitialPop[Sample] = 10
ModelName = paste0("CurrentClim_DetProb_",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClim_DetProb_",DetectionProb)
OutputDir = paste0(ExternalInvasionDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
  ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = 0.99,           #Annual mortality probability under management
  SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = InitialPop,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  #InitialInfo = sample(1:nrow(d),10),
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = K,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
Time <- End-Start

################################
###Initial invasion as binary vector
################################
InitialPop = rep(0,times = nrow(d))
Sample <- sample(1:nrow(d),10)
InitialPop[Sample] = 1
ModelName = paste0("CurrentClim_DetProb_InitialInvasion",DetectionProb)
if(DoClimateChange == TRUE)
ModelName = paste0("FutureClim_DetProb_InitialInvasion",DetectionProb)
OutputDir = paste0(ExternalInvasionDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
ModelName = ModelName,
Nperm = Nperm,                  #Number of permutations per parameter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = DetectionProb,  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
MortalityProb = 0.99,           #Annual mortality probability under management
SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
InitialPopulation = InitialPop,        #Population size at start if simulations
InitBioP = NA,		#Proportion of nodes infested at start of simulations
InvasionRisk = NA,  #Vector of probabilities of external invasion risk
InitialInfo = NA, #
EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
K = K,		       #Population carrying capacity
PropaguleProduction = alpha*d$R0, #Propagules produced per individual
PropaguleEstablishment = estab, #Propagules establishment rate
IncursionStartPop=10,      #option to set population size for new incursions
SDDprob = nd,                   #Natural disperal probability between each pair of nodes
SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
			      #e.g. could be weighted by law of human visitation or data on stock movements
LDDrate = H_vectors,         #Proportion of available propagules entering LDD
OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
OngoingExternalInfo = F,   ##Option to include ongoing invasion from external sources
OutputDir = OutputDir,		      #Directory for storing results
DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
Time <- End-Start

################################
###Initial invasion as proportion of nodes
################################
InitInvP = 10/nrow(d)
ModelName = paste0("CurrentClim_DetProb_InitInvP",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClim_DetProb_InitInvP",DetectionProb)
OutputDir = paste0(ExternalInvasionDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
  ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = 0.99,           #Annual mortality probability under management
  SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = NA,        #Population size at start if simulations
  InitBioP = InitInvP,		#Proportion of nodes infested at start of simulations
  InvasionRisk = NA,  #Vector of probabilities of external invasion risk
  InitialInfo = NA, #
  InitInfoP = NA, #
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = K,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
  OngoingExternalInfo = F,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
Time <- End-Start

################################
###Initial invasion as vector of probabilities
################################
ModelName = paste0("CurrentClim_DetProb_ExternalInvProb",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClim_DetProb_ExternalInvProb",DetectionProb)
OutputDir = paste0(ExternalInvasionDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
  ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = 0.99,           #Annual mortality probability under management
  SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = NA,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  InitialInfo = NA, #
  InitInfoP = NA, #
  ExternalInfoProb = NA, #
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = K,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
  OngoingExternalInfo = F,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
Time <- End-Start


################################
###Initial invasion as vector of probabilities
###and proportion of nodes
################################
InitInvP = 10/nrow(d)
ModelName = paste0("CurrentClim_DetProb_ExternalInvProbInitInvP",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClim_DetProb_ExternalInvProbInitInvP",DetectionProb)
OutputDir = paste0(ExternalInvasionDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
  ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = 0.99,           #Annual mortality probability under management
  SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = NA,        #Population size at start if simulations
  InitBioP = InitInvP,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  InitInfoP = NA, #
  ExternalInfoProb = NA, #
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = K,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
  OngoingExternalInfo = F,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
Time <- End-Start

################################
###Initial invasion as vector of probabilities
###with ongoing external invasion
################################
ModelName = paste0("CurrentClim_DetProb_OngoingExternalInv",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClim_DetProb_OngoingExternalInv",DetectionProb)
OutputDir = paste0(ExternalInvasionDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
  ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = 0.99,           #Annual mortality probability under management
  SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = NA,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  ExternalInfoProb = NA, #
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = K,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = T,   ##Option to include ongoing invasion from external sources
  OngoingExternalInfo = F,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
Time <- End-Start

################################
###Initial invasion as matrix of probabilities
###with ongoing external invasion
################################
ExternalInvProbMatrix = matrix(ncol = Ntimesteps, nrow  = nrow(d))
for(i in 1:20)
  ExternalInvProbMatrix[,i] = PropTotalHumans*AnnualIncursionRate
for(i in 21:Ntimesteps)
  ExternalInvProbMatrix[,i] = PropTotalHumans*AnnualIncursionRate*10
ModelName = paste0("CurrentClim_DetProb_ExternalInvMatrix",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClim_DetProb_ExternalInvMatrix",DetectionProb)
OutputDir = paste0(ExternalInvasionDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
  ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = 0.99,           #Annual mortality probability under management
  SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = NA,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = ExternalInvProbMatrix,  #Vector of probabilities of external invasion risk
  ExternalInfoProb = NA, #
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = K,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = T,   ##Option to include ongoing invasion from external sources
  OngoingExternalInfo = F,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
Time <- End-Start