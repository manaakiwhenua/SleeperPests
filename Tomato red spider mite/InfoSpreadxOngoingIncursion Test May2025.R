############################################################
############################################################
###Test ongoing incursion and info spread functionality
############################################################
############################################################

Nperm = 30
Ntimesteps = 50
DetectionProbs <- c(0.05,0.1,0.2,0.3,0.5)

###Make this high for obvious effect
AnnualIncursionRate = 0.5

InfoSpreadDir <-  paste0(ResultsDir,"InfoSpread_x_OngoingIncursion/")
dir.create(InfoSpreadDir)
for(detprob in 1:length(DetectionProbs))
{
################################################
###With OngoingExternalInvasion == F
################################################  
  
###Don't provide SEAM for scenario of no info spread  
  ModelName = paste0("NoInfoSpreadStandardCurrentClim_DetProb_",DetectionProbs[detprob])
  if(DoClimateChange == TRUE)
    ModelName = paste0("NoInfoSpreadStandardFutureClim_DetProb_",DetectionProbs[detprob])
  OutDir <- paste0(InfoSpreadDir,ModelName,"/")
  dir.create(OutDir)
  INApestMetaParallel(
    ModelName = ModelName,
    Nperm = Nperm,                  #Number of permutations per parameter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProbs[detprob],  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
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
    OutputDir = OutDir,		      #Directory for storing results
    DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
  )
  
  ###Provide SEAM to permit info spread
  ModelName = paste0("InfoSpreadStandardCurrentClim_DetProb_",DetectionProbs[detprob])
  if(DoClimateChange == TRUE)
    ModelName = paste0("InfoSpreadStandardFutureClim_DetProb_",DetectionProbs[detprob])
  OutDir <- paste0(InfoSpreadDir,ModelName,"/")
  dir.create(OutDir)
  INApestMetaParallel(
    ModelName = ModelName,
    Nperm = Nperm,                  #Number of permutations per parameter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProbs[detprob],  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
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
    SEAM = SEAM,			#Option to provide socioeconomic adjacency matrix for information spread
    LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
    #e.g. could be weighted by law of human visitation or data on stock movements
    LDDrate = H_vectors,         #Proportion of available propagules entering LDD
    OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
    OutputDir = OutDir,		      #Directory for storing results
    DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
  )

################################################
###With OngoingExternalInvasion == T
################################################  
  ###Don't provide SEAM for scenario of no info spread     
  ModelName = paste0("NoInfoSpreadOngoingCurrentClim_DetProb_",DetectionProbs[detprob])
  if(DoClimateChange == TRUE)
    ModelName = paste0("NoInfoSpreadOngoingFutureClim_DetProb_",DetectionProbs[detprob])
  OutDir <- paste0(InfoSpreadDir,ModelName,"/")
  dir.create(OutDir)
  INApestMetaParallel(
    ModelName = ModelName,
    Nperm = Nperm,                  #Number of permutations per parameter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProbs[detprob],  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
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
    OngoingExternalInvasion = T,   ##Option to include ongoing invasion from external sources
    OutputDir = OutDir,		      #Directory for storing results
    DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
  )
  
  ###Provide SEAM to permit info spread
  ModelName = paste0("InfoSpreadOngoingCurrentClim_DetProb_",DetectionProbs[detprob])
  if(DoClimateChange == TRUE)
    ModelName = paste0("InfoSpreadOngoingFutureClim_DetProb_",DetectionProbs[detprob])
  OutDir <- paste0(InfoSpreadDir,ModelName,"/")
  dir.create(OutDir)
  INApestMetaParallel(
    ModelName = ModelName,
    Nperm = Nperm,                  #Number of permutations per parameter combination
    Ntimesteps = Ntimesteps,                 #Simulation duration
    DetectionProb = DetectionProbs[detprob],  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
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
    SEAM = SEAM,			#Option to provide socioeconomic adjacency matrix for information spread
    LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
    #e.g. could be weighted by law of human visitation or data on stock movements
    LDDrate = H_vectors,         #Proportion of available propagules entering LDD
    OngoingExternalInvasion = T,   ##Option to include ongoing invasion from external sources
    OutputDir = OutDir,		      #Directory for storing results
    DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
  )
}
