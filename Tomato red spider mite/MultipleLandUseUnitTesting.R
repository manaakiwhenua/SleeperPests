Nperm = 30
Ntimesteps = 50

############################################################################
###Provide K as matrix of nodes x land uses
###Just use artificial example to test new functionality
############################################################################
Planduses = c(0.5,0.5)
Klanduse = matrix(nrow = length(K),ncol = 2)
Klanduse[,1] = floor(K*Planduses[1])
Klanduse[,2] = ceiling(K*Planduses[2])
N0landuse  = cbind(floor(d$N0*Planduses[1]),ceiling(d$N0*Planduses[2]))
#Klanduse[118,1] = 10
#N0landuse[118,] = 10
###########################################################
###Serial function
###########################################################

ModelName = "MultipleLandUseZeroMan"
OutputDir = paste0(ResultsDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
set.seed(42)
INApestMetaMultipleLandUse(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  Nlanduses = 2,
  DetectionProb = c(0,0),  #Annual per-individual detection Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  ManageProb = c(0.5,0.99),             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = c(0.5,0.99),           #Annual mortality probability under management
  SpreadReduction = c(0,0.9),        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = N0landuse,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  #InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  #InitialInfo = sample(1:nrow(d),10),
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = Klanduse,		       #Population carrying capacity
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
ParallelTime <- End-Start

ModelName = "MultipleLandUseDet_0.5"
OutputDir = paste0(ResultsDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
set.seed(42)
INApestMetaMultipleLandUse(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  Nlanduses = 2,
  DetectionProb = c(0.5,0.5),  #Annual per-individual detection Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  ManageProb = c(0.5,0.99),             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = c(0.5,0.99),           #Annual mortality probability under management
  SpreadReduction = c(0,0.9),        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = N0landuse,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  #InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  #InitialInfo = sample(1:nrow(d),10),
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = Klanduse,		       #Population carrying capacity
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
ParallelTime <- End-Start

###########################################################
###Parallel function
###########################################################

ModelName = "MultipleLandUseZeroManParallel"
OutputDir = paste0(ResultsDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
set.seed(42)
INApestMetaParallelMultipleLandUse(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  Nlanduses = 2,
  DetectionProb = c(0,0),  #Annual per-individual detection Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  ManageProb = c(0.5,0.99),             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = c(0.5,0.99),           #Annual mortality probability under management
  SpreadReduction = c(0,0.9),        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = N0landuse,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  #InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  #InitialInfo = sample(1:nrow(d),10),
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = Klanduse,		       #Population carrying capacity
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
ParallelTime <- End-Start

ModelName = "MultipleLandUseDet_0.5_Parallel"
OutputDir = paste0(ResultsDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
set.seed(42)
INApestMetaParallelMultipleLandUse(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  Nlanduses = 2,
  DetectionProb = c(0.5,0.5),  #Annual per-individual detection Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  ManageProb = c(0.5,0.99),             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = c(0.5,0.99),           #Annual mortality probability under management
  SpreadReduction = c(0,0.9),        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = N0landuse,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  #InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  #InitialInfo = sample(1:nrow(d),10),
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = Klanduse,		       #Population carrying capacity
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
ParallelTime <- End-Start