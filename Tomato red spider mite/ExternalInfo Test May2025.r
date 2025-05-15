#####################################################################################
###Test initial info and ongoing external communication functionality. 
#####################################################################################
Nperm = 30
Ntimesteps = 50
DetectionProb = 0.1
AnnualIncursionRate = 0.1
InfoDir <- paste0(ResultsDir,"InfoInputTest/")
dir.create(InfoDir)

################################
###Detection prob 0.1 no initial info
################################
ModelName = paste0("CurrentClimNoInitInfo",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClimNoInitInfo",DetectionProb)
OutputDir = paste0(InfoDir,ModelName,"/")
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
  InitialPopulation = d$N0,        #Population size at start if simulations
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
###Initial info as binary vector
################################
InitialInfo = rep(0,times = nrow(d))
Sample <- sample(1:nrow(d),10)
InitialInfo[Sample] = 1
ModelName = paste0("CurrentClimInitInfoBinaryVector",DetectionProb)
if(DoClimateChange == TRUE)
ModelName = paste0("FutureClimInitInfoBinaryVector",DetectionProb)
OutputDir = paste0(InfoDir,ModelName,"/")
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
InitialPopulation = d$N0,        #Population size at start if simulations
InitBioP = NA,		#Proportion of nodes infested at start of simulations
InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
InitialInfo = InitialInfo, #
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
###Initial info as proportion of nodes
################################
InitInfoP = 10/nrow(d)
ModelName = paste0("CurrentClimInitInfoPNodes",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClimInitInfoPnodes",DetectionProb)
OutputDir = paste0(InfoDir,ModelName,"/")
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
  InitialPopulation = d$N0,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  InitialInfo = NA, #
  InitInfoP = InitInfoP, #
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
###Initial info as vector of probabilities
################################
ExternalInfoProb = rep(10/nrow(d),times = nrow(d))
ModelName = paste0("CurrentClimInitInfoProbVector",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClimInitInfoProbVector",DetectionProb)
OutputDir = paste0(InfoDir,ModelName,"/")
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
  InitialPopulation = d$N0,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  InitialInfo = NA, #
  InitInfoP = NA, #
  ExternalInfoProb = ExternalInfoProb, #
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
###Initial info as vector of probabilities
###and proportion of nodes
################################
ModelName = paste0("CurrentClimInitInfoProbVectorPnodes",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClimInitInfoProbVectorPnodes",DetectionProb)
OutputDir = paste0(InfoDir,ModelName,"/")
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
  InitialPopulation = d$N0,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  InitInfoP = 10/nrow(d), #
  ExternalInfoProb = ExternalInfoProb, #
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
###Initial info as vector of probabilities
###with ongoing external communication
################################
ModelName = paste0("CurrentClimOngoingExternalInfoVector",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClimOngoingExternalInfoVector",DetectionProb)
OutputDir = paste0(InfoDir,ModelName,"/")
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
  InitialPopulation = d$N0,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  ExternalInfoProb = ExternalInfoProb, #
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
  OngoingExternalInfo = T,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
Time <- End-Start

################################
###Initial info as matrix of probabilities
###with ongoing external communication
################################
ExternalInfoProbMatrix = matrix(ncol = Ntimesteps, nrow  = nrow(d))
for(i in 1:20)
  ExternalInfoProbMatrix[,i] = ExternalInfoProb
for(i in 21:Ntimesteps)
  ExternalInfoProbMatrix[,i] = ExternalInfoProb*10
ModelName = paste0("CurrentClimOngoingExternalInfoMatrix",DetectionProb)
if(DoClimateChange == TRUE)
  ModelName = paste0("FutureClimOngoingExternalInfoMatrix",DetectionProb)
OutputDir = paste0(InfoDir,ModelName,"/")
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
  InitialPopulation = d$N0,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  ExternalInfoProb = ExternalInfoProbMatrix, #
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
  OngoingExternalInfo = T,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
Time <- End-Start