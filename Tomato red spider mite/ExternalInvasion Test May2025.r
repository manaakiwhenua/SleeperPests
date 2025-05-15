#####################################################################################
###Test initial invasion and ongoing external invasion functionality
#####################################################################################
Nperm = 30
Ntimesteps = 50
DetectionProb = 0.1

###Make this high for obvious effect
AnnualIncursionRate = 0.1

InvasionDir <- paste0(ResultsDir,"InvasionInputTest/")
dir.create(InvasionDir)

################################
###Detection prob 0.1 no external invasion
################################
InitialPop = rep(0,times = nrow(d))
Sample <- sample(1:nrow(d),10)
InitialPop[Sample] = 10
ModelName = paste0("CurrentClimInitInvPopVector",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClimInitInvPopVector",DetectionProb)
OutputDir = paste0(InvasionDir,ModelName,"/")
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
InitialPopBinary <- InitialPop
InitialPopBinary[InitialPop>0] = 1
ModelName = paste0("CurrentClimInitInvBinary",DetectionProb)
if(DoClimateChange == TRUE)
ModelName = paste0("FutureClimInitInvBinary",DetectionProb)
OutputDir = paste0(InvasionDir,ModelName,"/")
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
InitialPopulation = InitialPopBinary,        #Population size at start if simulations
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
ModelName = paste0("CurrentClimInitInvPnodes",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClimInitInvPnodes",DetectionProb)
OutputDir = paste0(InvasionDir,ModelName,"/")
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
ModelName = paste0("CurrentClimInitInvProbVector",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClimInitInvProbVector",DetectionProb)
OutputDir = paste0(InvasionDir,ModelName,"/")
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
ModelName = paste0("CurrentClimInitInvProbVectorxPnodes",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClimInitInvProbVectorxPnodes",DetectionProb)
OutputDir = paste0(InvasionDir,ModelName,"/")
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
ModelName = paste0("CurrentClimOngoingExternalInvVector",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClimOngoingExternalInvVector",DetectionProb)
OutputDir = paste0(InvasionDir,ModelName,"/")
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
ModelName = paste0("CurrentClimExternalInvMatrix",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClimExternalInvMatrix",DetectionProb)
OutputDir = paste0(InvasionDir,ModelName,"/")
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