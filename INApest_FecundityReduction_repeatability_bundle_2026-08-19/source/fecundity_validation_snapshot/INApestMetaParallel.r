###########################################################################
###########################################################################
###Declares a function overlaying management on a metapopulation spread model 
###Key inputs are: 
###1) matrix of natural dispersal probibilities between each pair of sites (i.e. nodes of the network). Matrix can be non-symmetrical (i.e. can have source and sink nodes) 
###2) Envionmentally-determined per-capita propagule production
###3) Envionmentally-determined carrying capacity (K)
###4) Management parameters
### a) Annual detection probability
### b) Annual management adoption probability subsequent to detection
### c) Annual mortality probability under management
### d) Spread reduction when management adopted
###5) A binary vector of intially infested nodes OR
###   A vector of invasion risk probabilities and proportion of nodes initially invaded
###Key outputs are:
###3-dimensional arrays of invasion, management and detection status for each node in each timestep of each permuation
###2-dimensional array of invasion probability (i.e. proportion of permutations pest present) for each node in each timestep
###Line graphs summarising number of total population (as a proportion of carrying capacity) nodes infested, where infestations
###are detected, the proportion of infestations detected and nodes under management against time 
###########################################################################
###########################################################################

#######################################################################
###This version implements parallel processing with a PSOCK cluster and parallel::parLapply for the permutation loop.
###See file "ParallelSetup.r" for notes on steps for setting up parallel processing
#######################################################################



INApestMetaParallel = function(
ModelName, #Name for storing results to file 
Nperm,                  #Number of permutations per parameter combination
Ntimesteps,                 #Simulation duration timesteps can be any length of time
DetectionProb,          #Per-individual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be single number or vector (nodes)
ManageProb,             #Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
ManageSD = NULL, #Option to provide standard deviation for management probability. Can be single number or vector (nodes)
MortalityProb,           #Mortality probability under management
MortalitySD = NULL, #Option to provide standard deviation for mortality probability. Can be single number or vector (nodes)
FecundityReduction = 0, #Proportional reduction in per-capita fecundity under management: scalar, vector (nodes), or matrix (nodes x timesteps)
SpreadReduction,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
SpreadReductionSD = NULL, #Option to provide standard deviation for spread reduction. Can be single number or vector (nodes)
InitialPopulation = NA,        #Vector of population sizes at start of simulations
InitBioP = NA,		#Proportion of nodes infested at start of simulations
InvasionRisk = NA,           #Vector or matrix (nodes x timesteps) of probabilities of invasion from external sources
InitialInfo = NA,        #Vector or of nodes with information at start of simulations
InitInfoP = NA,		#Proportion of nodes with information at start of simulations
ExternalInfoProb = NA,           #Vector of probabilities of communication from external sources
InfoRetentionProb = 1,       #Probability that existing information is retained between timesteps. Can be single number, vector (nodes) or matrix (nodes x timesteps)
InfoPersistenceSteps = NA,    #Number of timesteps information persists after last known local presence. Can be single number, vector (nodes) or matrix (nodes x timesteps); NA uses InfoRetentionProb
EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
K,		       #Population carrying capacity - vector (nodes)
PropaguleProduction, #Propagules produced per individual, can be single value, vector (nodes) or matrix (nodes x years)
PropaguleEstablishment, #Propagules establishment probability. The likelihood of a dispersing propagule encountering a single
                        #host plant or establishment site within a node. Can be a ratio of search radius or patch size to node area
IncursionStartPop=NA,      #option to set population size for new incursions
SDDprob,                   #Natural dispersal probability matrix, or 3D array (nodes x nodes x timesteps)
SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
LDDprob = NA,         #Option to provide long distance (human-mediated) dispersal matrix or 3D array (nodes x nodes x timesteps) instead of distance-independent dispesal rate
			      #e.g. could be weighted by law of human visitation or data on stock movements
LDDrate = 0,         #Proportion of available propagules entering LDD
OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
OngoingExternalInfo = F,   ##Option to include ongoing communication from external sources
OutputDir = NA,		      #Directory for storing results
DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
{
###POTENTIAL ADDITIONS
###1) Make detection prob a function of population size. Could be based on individual detection prob so that DetectionProb = 1-(1-DPindividual)^N)
###   DPindividual could vary between nodes


###Allow SDD and LDD connectivity to vary through time
if(length(dim(SDDprob)) == 3 && (dim(SDDprob)[1] != dim(SDDprob)[2] || dim(SDDprob)[3] != Ntimesteps))
  stop("SDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")
if(length(dim(LDDprob)) == 3 && (dim(LDDprob)[1] != nrow(SDDprob) || dim(LDDprob)[2] != nrow(SDDprob) || dim(LDDprob)[3] != Ntimesteps))
  stop("LDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")
###Allow management-induced fecundity reduction to vary by node and through time
if(is.matrix(FecundityReduction)) {
  if(nrow(FecundityReduction) != nrow(SDDprob) || ncol(FecundityReduction) != Ntimesteps)
    stop("FecundityReduction matrix must have dimensions nodes x Ntimesteps")
} else if(!(length(FecundityReduction) == 1 || length(FecundityReduction) == nrow(SDDprob))) {
  stop("FecundityReduction must be a single value, vector of length nodes, or matrix nodes x Ntimesteps")
}
if(any(!is.finite(FecundityReduction)) || any(FecundityReduction < 0) || any(FecundityReduction > 1))
  stop("FecundityReduction values must be between 0 and 1")

###Allow information retention to vary by node and through time
if(is.matrix(InfoRetentionProb) == T && (nrow(InfoRetentionProb) != nrow(SDDprob) || ncol(InfoRetentionProb) != Ntimesteps))
  stop("InfoRetentionProb matrix must have dimensions nodes x Ntimesteps")
if(is.matrix(InfoRetentionProb) == F && !(length(InfoRetentionProb) == 1 || length(InfoRetentionProb) == nrow(SDDprob)))
  stop("InfoRetentionProb must be a single value, vector of length nodes, or matrix nodes x Ntimesteps")
if(any(is.na(InfoRetentionProb)) || any(InfoRetentionProb < 0) || any(InfoRetentionProb > 1))
  stop("InfoRetentionProb values must be between 0 and 1")

###Allow programmed information persistence after last known local presence to vary by node and through time
if(is.matrix(InfoPersistenceSteps) == T && (nrow(InfoPersistenceSteps) != nrow(SDDprob) || ncol(InfoPersistenceSteps) != Ntimesteps))
  stop("InfoPersistenceSteps matrix must have dimensions nodes x Ntimesteps")
if(is.matrix(InfoPersistenceSteps) == F && !(length(InfoPersistenceSteps) == 1 || length(InfoPersistenceSteps) == nrow(SDDprob)))
  stop("InfoPersistenceSteps must be a single value, vector of length nodes, or matrix nodes x Ntimesteps")
if(any(!is.na(InfoPersistenceSteps) & (!is.finite(InfoPersistenceSteps) | InfoPersistenceSteps < 0 | InfoPersistenceSteps != floor(InfoPersistenceSteps))))
  stop("InfoPersistenceSteps values must be non-negative whole numbers or NA")
UseInfoPersistence = any(!is.na(InfoPersistenceSteps))
if(UseInfoPersistence == T && any(InfoRetentionProb < 1))
  warning("Both InfoPersistenceSteps and InfoRetentionProb specify information loss. Programmed stopping takes priority where InfoPersistenceSteps is not NA; InfoRetentionProb is only used where InfoPersistenceSteps is NA.",call. = F)

# pre-evaluate some variables for efficiency
if(is.matrix(K) == FALSE)
{
K_is_0 <- K<=0
inv_K <- 1 / sum(K)
NodeK = K
}

if(is.matrix(PropaguleProduction) == FALSE)
 NodePropaguleProduction = PropaguleProduction

if(is.matrix(PropaguleEstablishment) == FALSE)
  NodePropaguleEstablishment = PropaguleEstablishment

if(is.matrix(EnvEstabProb) == F)
  NodeEnvEstabProb <- EnvEstabProb

if(is.matrix(Survival) == F)
  NodeSurvival <- Survival


###Declare matrix for information spread simulations
if(is.matrix(SEAM) == T)
     {
     diag(SEAM) = 0
     RandSEAM <- matrix(NA,nrow = nrow(SDDprob),ncol=nrow(SDDprob))
     }

###Assign standard deviation value to management in no value provided
if(is.null(ManageSD) == T)
	ManageSD = mean(ManageProb)/10
if(is.null(SpreadReductionSD) == T)
	SpreadReductionSD = (1-mean(SpreadReduction))/10
if(is.null(DetectionSD) == T)
	DetectionSD = mean(DetectionProb)/10
if(is.null(MortalitySD) == T)
	MortalitySD = mean(MortalityProb)/10


###########################################################
###Start of simulation
###########################################################
    
detected_cores <- parallel::detectCores()
if (is.na(detected_cores)) detected_cores <- 2L
n_cores <- max(1L, min(Nperm, detected_cores - 1L))

###Run one stochastic realisation. Function arguments and local helpers are
###captured in this closure, avoiding fragile manual worker export lists.
PermutationWorker <- function(i_perm)
  {
  ###Max integer for propagule dispersal using rmultinom
  MaxInteger <- .Machine$integer.max  

  ###Set initial dispersal connectivity
  NodeSDDprob = SDDprob
  if(length(dim(SDDprob)) == 3)
    NodeSDDprob = SDDprob[,,1]
  NodeLDDprob = LDDprob
  if(length(dim(LDDprob)) == 3)
    NodeLDDprob = LDDprob[,,1]
  NodeFecundityReduction <- if(is.matrix(FecundityReduction)) FecundityReduction[,1] else FecundityReduction
  if(length(NodeFecundityReduction) == 1)
    NodeFecundityReduction <- rep(NodeFecundityReduction, nrow(SDDprob))
  
  InvasionResultsLoop <- array(dim = c(nrow(SDDprob),Ntimesteps))
  PopulationResultsLoop <- InvasionResultsLoop
  ManagingResultsLoop <- InvasionResultsLoop
  DetectedResultsLoop <- InvasionResultsLoop

  
  ###If carrying capacity provided as matrix assign values from first timestep for population initialisation
  if(is.matrix(K) == TRUE)
    {
    K_is_0 <- K[,1]<=0
    inv_K <- 1 / sum(K[,1])
    NodeK = K[,1] 
    }      
  
  
###Assign initial infestations according either to "InitialInvasion" binary vector OR
###"InvasionRisk" probabilities and/or initial proportion of nodes infested ("InitBioP") OR
###just "InitBioP" if neither "InitialInvasion" or "InvasionRisk" supplied by user
InitBio = rep(0,times = nrow(SDDprob))

if(length(InitialPopulation) == nrow(SDDprob))
  InitBio = InitialPopulation

if(length(InitialPopulation) != nrow(SDDprob))
{
risk = NULL
if(is.matrix(InvasionRisk) == T && nrow(InvasionRisk) == nrow(SDDprob))
  risk = InvasionRisk[,1]
if(is.matrix(InvasionRisk) == F && length(InvasionRisk) == nrow(SDDprob))
  risk = InvasionRisk

if(is.na(InitBioP) == F)
  Infested = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitBioP),prob = risk)
if(is.na(InitBioP) == T && is.null(risk) == F)
  {
  Infested = rbinom(1:nrow(SDDprob),size = 1,prob = risk)
  Infested = which(Infested == 1)
  }
if(is.na(InitBioP) == T && is.null(risk) == T)
  Infested = integer(0)

if(is.na(IncursionStartPop) == T)
  InitBio[Infested] = 1
if(is.na(IncursionStartPop) == F)
  InitBio[Infested] = IncursionStartPop
}

###Ensure initial population not greater than carrying capacity
InitBio[InitBio > NodeK] = NodeK[InitBio > NodeK] 

# initialise the population
N <- InitBio
if(sum(N) == 0 && OngoingExternalInvasion == F)
  warning("No initial populations and no future external invasions")

###Select nodes with information at start of simulation  according either to "InitialInfo" binary vector OR
###"ExternalInfoProb" probabilities and/or initial proportion of nodes with information ("InitInfoP") OR
###just "InitInfoP" if neither "InitialInfo" or "ExternalInfoProb" supplied by user.
###If no initial info variables provided, no nodes have info at start of simulations
InitInfo = rep(0,times = nrow(SDDprob))
if(length(InitialInfo) == nrow(SDDprob) || (is.na(InitInfoP) == F && InitInfoP>0) || is.na(sum(ExternalInfoProb)) == F )
{
if(length(InitialInfo) != nrow(SDDprob))
  {
  if(length(ExternalInfoProb) == nrow(SDDprob))
    {
    if(is.na(InitInfoP) == F)
      Info = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitInfoP),prob = ExternalInfoProb)
    if(is.na(InitInfoP) == T)
      {
      Info = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb) 
      Info = which(Info == 1)
      } 
    }
  if(length(ExternalInfoProb) != nrow(SDDprob))
    {
    if(is.matrix(ExternalInfoProb) == F)
      Info = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitInfoP))
    if(is.matrix(ExternalInfoProb) == T)
      {
      Info = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb[,1])
      Info = which(Info == 1)
      }
    }
  InitInfo[Info] = 1
  
  }
if(length(InitialInfo) == nrow(SDDprob))
  InitInfo = InitialInfo  
}

###Randomly assign annual detection probability, based on mean and sd
###If DetectionProb given as single value or vector (nodes)
if(is.matrix(DetectionProb)==FALSE &&(length(DetectionProb) == 1 ||length(DetectionProb) == nrow(SDDprob) ))
      {
      NodeDetectionProb = rnorm(DetectionProb,DetectionSD,n = nrow(SDDprob))
      NodeDetectionProb[NodeDetectionProb<0] = 0
      NodeDetectionProb[NodeDetectionProb>1] = 1
      }

###If DetectionProb given as matrix (nodes x timesteps) use values for first timestep to get initial detections
if(is.matrix(DetectionProb)==TRUE && nrow(DetectionProb) == nrow(SDDprob) && ncol(DetectionProb) == Ntimesteps)
      {
      NodeDetectionProb = rnorm(DetectionProb[,1],DetectionSD,n = nrow(SDDprob))
      NodeDetectionProb[NodeDetectionProb<0] = 0
      NodeDetectionProb[NodeDetectionProb>1] = 1
      }


###Randomly assign probability of mangement adoption upon detection of infestation
###If ManageProb given as single value or vector (nodes)
if(is.matrix(ManageProb)==FALSE &&(length(ManageProb) == 1 ||length(ManageProb) == nrow(SDDprob) ))
      {
      NodeManageProb = rnorm(ManageProb,ManageSD,n = nrow(SDDprob))
      NodeManageProb[NodeManageProb<0] = 0
      NodeManageProb[NodeManageProb>1] = 1
      }

###Randomly assign spread reduction factor when management adopted
###If SpreadReduction given as single value or vector (nodes)
if(is.matrix(SpreadReduction)==FALSE &&(length(SpreadReduction) == 1 ||length(SpreadReduction) == nrow(SDDprob) ))
      {
      NodeSpreadReduction = rnorm(SpreadReduction,ManageSD,n = nrow(SDDprob))
      NodeSpreadReduction[NodeSpreadReduction<0] = 0
      NodeSpreadReduction[NodeSpreadReduction>1] = 1
      }

###Randomly assign mortality probability when management applied
###If MortalityProb given as single value or vector (nodes)
if(is.matrix(MortalityProb)==FALSE &&(length(MortalityProb) == 1 ||length(MortalityProb) == nrow(SDDprob) ))
      {
      NodeMortalityProb = rnorm(MortalityProb,MortalitySD,n = nrow(SDDprob))
      NodeMortalityProb[NodeMortalityProb<0] = 0
      NodeMortalityProb[NodeMortalityProb>1] = 1
      }

###Populate invasion status vector ahead of timestep loop
Invaded = ifelse(InitBio>0,1,0) 

###Probability of info at start of simulation depends on
###Presence of pest and detection probability
###Select nodes that have detected infestation 
InitDetection = rbinom(1:nrow(SDDprob),size = 1,prob = 1-(1-NodeDetectionProb)^InitBio)
InitInfo[InitInfo == 0] = InitDetection[InitInfo == 0]
###Populate information status vector ahead of timestep loop
HaveInfo = InitInfo

###Track the most recent timestep with known local presence
LastKnownPresence = rep(NA,nrow(SDDprob))
if(UseInfoPersistence == T)
  {
  InitialKnownPresence = which(InitDetection == 1)
  if(length(InitialKnownPresence) > 0)
    LastKnownPresence[InitialKnownPresence] = 0
  }

    
  # run simulation
for(timestep in 1:Ntimesteps) 
  { 
 
  ###Allow for variation in dispersal connectivity through time
  if(length(dim(SDDprob)) == 3)
    NodeSDDprob = SDDprob[,,timestep]
  if(length(dim(LDDprob)) == 3)
    NodeLDDprob = LDDprob[,,timestep]
  if(is.matrix(FecundityReduction))
    NodeFecundityReduction <- FecundityReduction[,timestep]

  ###Allow for variation in establishment through time
  ###e.g.  climate change predictions
  ###Note: could be done outside loop, but would take heaps of memory to store 
  if(is.matrix(EnvEstabProb) == T)
    NodeEnvEstabProb <- EnvEstabProb[,timestep]
   
  if(is.matrix(Survival) == T)
    NodeSurvival <- Survival[,timestep]
    
    
  ###If carrying capacity provided as matrix assign values for relevant timestep
  if(is.matrix(K) == TRUE)
    {
    K_is_0 <- K[,timestep]<=0
    inv_K <- 1 / sum(K[,timestep])
    NodeK = K[,timestep] 
    }  

  ###If propagule production provided as matrix assign values for relevant timestep
  if(is.matrix(PropaguleProduction) == TRUE)
    NodePropaguleProduction = PropaguleProduction[,timestep] 
  
  if(is.matrix(PropaguleEstablishment) == TRUE)
    NodePropaguleEstablishment = PropaguleEstablishment[,timestep]
      
  ###Randomly assign annual detection probability, based on mean and sd
  ###If DetectionProb given as matrix (nodes x timesteps)
  if(is.matrix(DetectionProb)==TRUE && nrow(DetectionProb) == nrow(SDDprob) && ncol(DetectionProb) == Ntimesteps)
   	{	
   	NodeDetectionProb = rnorm(DetectionProb[,timestep],DetectionSD,n = nrow(SDDprob))
   	NodeDetectionProb[NodeDetectionProb<0] = 0
   	NodeDetectionProb[NodeDetectionProb>1] = 1
   	}

  ###Randomly assign probability of mangement adoption upon detection of infestation
  ###If ManageProb given as matrix (nodes x timesteps)
  if(is.matrix(ManageProb)==TRUE && nrow(ManageProb) == nrow(SDDprob) && ncol(ManageProb) == Ntimesteps)
   	{	
   	NodeManageProb = rnorm(ManageProb[,timestep],ManageSD,n = nrow(SDDprob))
   	NodeManageProb[NodeManageProb<0] = 0
   	NodeManageProb[NodeManageProb>1] = 1
   	}

  ###Randomly assign spread reduction factor when management adopted
  ###If SpreadReduction given as matrix (nodes x timesteps)
  if(is.matrix(SpreadReduction)==TRUE && nrow(SpreadReduction) == nrow(SDDprob) && ncol(SpreadReduction) == Ntimesteps)
   	{	
   	NodeSpreadReduction = rnorm(SpreadReduction[,timestep],ManageSD,n = nrow(SDDprob))
   	NodeSpreadReduction[NodeSpreadReduction<0] = 0
   	NodeSpreadReduction[NodeSpreadReduction>1] = 1
   	}
  
  ###Randomly assign annual mortality probability when management applied
  ###If MortalityProb given as matrix (nodes x timesteps)
  if(is.matrix(MortalityProb)==TRUE && nrow(MortalityProb) == nrow(SDDprob) && ncol(MortalityProb) == Ntimesteps)
      {
      NodeMortalityProb = rnorm(MortalityProb[,timestep],MortalitySD,n = nrow(SDDprob))
      NodeMortalityProb[NodeMortalityProb<0] = 0
      NodeMortalityProb[NodeMortalityProb>1] = 1
      }


  ###Assign management status to nodes   
  ###Management is only applied to nodes which have information
  ###i.e. where pest has been detected or following communication of information
  ###from neighbouring infested farms 
  Managing = rbinom(1:nrow(SDDprob),size = 1,prob = NodeManageProb*HaveInfo)
  
  ###Identify nodes with known extant infestations 
  Detected = Invaded*HaveInfo
  
  ###Adjust starting population for natural and managed mortality
  N0 = rbinom(nrow(SDDprob),N,NodeSurvival*(1-NodeMortalityProb*Managing))

  ###Track known local presence from actual management mortality
  ###Condition on realised total deaths so the existing population draw is unchanged
  if(UseInfoPersistence == T)
    {
    ManagementMortality = NodeMortalityProb*Managing
    TotalMortalityProb = 1-NodeSurvival*(1-ManagementMortality)
    ConditionalManagementMortality = rep(0,nrow(SDDprob))
    ManagementMortalityNodes = which((N-N0) > 0 & ManagementMortality > 0 & TotalMortalityProb > 0)
    if(length(ManagementMortalityNodes) > 0)
      ConditionalManagementMortality[ManagementMortalityNodes] = (NodeSurvival*ManagementMortality)[ManagementMortalityNodes]/TotalMortalityProb[ManagementMortalityNodes]
    ManagementKillProb = 1-(1-ConditionalManagementMortality)^(N-N0)
    ManagementKillProb[ManagementKillProb < 0] = 0
    ManagementKillProb[ManagementKillProb > 1] = 1
    CertainManagementKillNodes = which(ManagementKillProb >= 1)
    if(length(CertainManagementKillNodes) > 0)
      LastKnownPresence[CertainManagementKillNodes] = timestep
    PotentialManagementKillNodes = which(ManagementKillProb > 0 & ManagementKillProb < 1)
    if(length(PotentialManagementKillNodes) > 0)
      {
      ManagementKilled = rbinom(n = length(PotentialManagementKillNodes),size = 1,prob = ManagementKillProb[PotentialManagementKillNodes])
      KnownPresence = PotentialManagementKillNodes[ManagementKilled == 1]
      if(length(KnownPresence) > 0)
        LastKnownPresence[KnownPresence] = timestep
      }
    }
  if(sum(N0)<=0 )
    N = N0
  Pin <-0
  Qin <- 0  
    # natural dispersal 
  if(sum(N0)>0 ) 
  {
    
  EffectiveReproductivePopulation <- N0 * (1-NodeFecundityReduction*Managing)
  Propagules <- rpois(nrow(SDDprob), NodePropaguleProduction * EffectiveReproductivePopulation)# propagules are produced after management effects on fecundity
  
  ###natural dispersal
  ###adjust propagule spread for environmentally determined establishment probability
  ###in receiving nodes
  Pout <- Propagules*(1-LDDrate)
  if(sum(Pout)>0 && sum(Pout)< MaxInteger) 
    Pin <- t(rmultinom(1, size=sum(Pout*rowSums(NodeSDDprob)), prob=Pout %*% NodeSDDprob))  # propagules are dispersed
  if(sum(Pout) >= MaxInteger) 
    Pin <- colSums(sweep(NodeSDDprob,1,Pout,`*`))
  ###human-mediated spread
  ###adjust propagule spread for environmentally determined establishment probability
  ###in receiving nodes
  if (is.matrix(NodeLDDprob)==T) 
    {
    Qout  = Propagules*LDDrate *(1-NodeSpreadReduction*Managing)       
    if(sum(Qout)>0 && sum(Qout) < MaxInteger)  
      Qin <- t(rmultinom(1, size=sum(Qout*rowSums(NodeLDDprob)), prob=Qout %*% NodeLDDprob))    # propagules are dispersed
    if(sum(Qout) >= MaxInteger) 
      Qin <- colSums(sweep(NodeLDDprob,1,Qout,`*`))
    }
  
  
  #Propagule establishment probability can be viewed as a spatial process indicating the likelihood
  #that a dispersing propagule will encounter and individual site or host.
  #The simplest approach is to express the search area as a proportion of the node area
  #For a weed this could be the area of individual sites (patches) relative to node area
  #Where K is the number of available patches within the node
  #This also incorporates density dependence since propagule success depends on availability of uninfested host plants 
  #or unoccupied patches
  ###Note can produce NAs if NodeK-N0 very large and using old version of R
  N <- ifelse(K_is_0, 0, N0 + rbinom(nrow(SDDprob), NodeK-N0, 1 - exp(-NodePropaguleEstablishment*NodeEnvEstabProb*(Pin+Qin))))

  } 
 ###Apply programmed stopping after last known local presence
NodeInfoPersistenceSteps = InfoPersistenceSteps
if(is.matrix(InfoPersistenceSteps) == T)
  NodeInfoPersistenceSteps = InfoPersistenceSteps[,timestep]
if(length(NodeInfoPersistenceSteps) == 1)
  NodeInfoPersistenceSteps = rep(NodeInfoPersistenceSteps,nrow(SDDprob))
ProgrammedInfoNodes = which(HaveInfo == 1 & !is.na(NodeInfoPersistenceSteps))
if(length(ProgrammedInfoNodes) > 0)
  {
  TimeSinceKnownPresence = timestep-LastKnownPresence
  InfoStopNodes = ProgrammedInfoNodes[is.na(LastKnownPresence[ProgrammedInfoNodes]) | TimeSinceKnownPresence[ProgrammedInfoNodes] >= NodeInfoPersistenceSteps[ProgrammedInfoNodes]]
  if(length(InfoStopNodes) > 0)
    HaveInfo[InfoStopNodes] = 0
  }

###Allow information to decay after management and spread where no programmed stop is supplied
NodeInfoRetentionProb = InfoRetentionProb
if(is.matrix(InfoRetentionProb) == T)
  NodeInfoRetentionProb = InfoRetentionProb[,timestep]
if(length(NodeInfoRetentionProb) == 1)
  NodeInfoRetentionProb = rep(NodeInfoRetentionProb,nrow(SDDprob))
InfoDecayNodes = which(HaveInfo == 1 & is.na(NodeInfoPersistenceSteps) & NodeInfoRetentionProb < 1)
if(length(InfoDecayNodes) > 0)
  HaveInfo[InfoDecayNodes] = rbinom(n = length(InfoDecayNodes),size = 1,prob = NodeInfoRetentionProb[InfoDecayNodes])

###Update info vector for any info spread (if SEAM supplied)
 ###Only zero values updated here so information can refresh nodes that lost information
 if(is.matrix(SEAM) == T)
  {
  RandSEAM[] <- rbinom(n=nrow(SDDprob)^2, size=1, prob = SEAM*Detected)
  InfoTransferred = ifelse(colSums(RandSEAM)>0,1,0)
  HaveInfo[HaveInfo == 0] = InfoTransferred[HaveInfo == 0]
  }
 
 
 ###Add invasion resulting from colonisation from external sources
 if(OngoingExternalInvasion == T)
  {
  if(is.matrix(InvasionRisk) == F)
   ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk)
  if(is.matrix(InvasionRisk) == T)
   ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk[,timestep])
  Invaded[Invaded == 0] = ExternalInvasion[Invaded==0]
  if(is.na(IncursionStartPop) == T) 
	N = N+ExternalInvasion
  if(is.na(IncursionStartPop) == F) 
	N = N+ExternalInvasion*IncursionStartPop
  N[N > NodeK] = NodeK[N > NodeK] 
  }
 
  ###Add nodes with information resulting from external sources
  if(OngoingExternalInfo == T)
    {
    if(is.matrix(ExternalInfoProb) == F)
      ExternalInfo = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb)
    if(is.matrix(ExternalInfoProb) == T)
      ExternalInfo = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb[,timestep])
    HaveInfo[HaveInfo == 0] = ExternalInfo[HaveInfo==0]
     }
 ###Update infestation vector
 Invaded = ifelse(N>0,1,0)
 
 ###Record nodes adopting management
 ManagingResultsLoop[,timestep] = Managing
  
 ###Record infested nodes
 InvasionResultsLoop[,timestep] = Invaded

 ###Record populations
 PopulationResultsLoop[,timestep] = N

 ###Select new nodes where infestation detected
 NewHaveInfo =  rbinom(1:nrow(SDDprob),size = 1,prob = 1-(1-NodeDetectionProb)^N)
 
 ###Record newly detected infestations as known local presence
 if(UseInfoPersistence == T)
   {
   KnownPresence = which(NewHaveInfo == 1)
   if(length(KnownPresence) > 0)
     LastKnownPresence[KnownPresence] = timestep
   }
 
 ###Add newly detected infestations to info vector
 ###Only zero values updated here so information can refresh nodes that lost information
 HaveInfo[HaveInfo==0] = NewHaveInfo[HaveInfo==0]  
 
 ###Record detection status
 DetectedResultsLoop[,timestep] = HaveInfo*Invaded 
 }
 simplify2array(list(InvasionResultsLoop, PopulationResultsLoop, ManagingResultsLoop, DetectedResultsLoop), higher = TRUE)
}

###Use a common PSOCK/parLapply architecture across parallel INApest variants.
###Static scheduling and L'Ecuyer-CMRG worker streams support reproducible
###parallel simulations when the caller fixes the R seed.
if(n_cores == 1L)
  {
  PermutationResults <- lapply(seq_len(Nperm), PermutationWorker)
  } else {
  cluster <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(if(inherits(cluster, "cluster")) parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterSetRNGStream(cluster)
  PermutationResults <- parallel::parLapply(cluster, seq_len(Nperm), PermutationWorker)
  parallel::stopCluster(cluster)
  cluster <- NULL
  }
PermOut <- simplify2array(PermutationResults, higher = TRUE)
###########################################################
###End of Simulation
###########################################################

###Extract results for invasion, managing and dectection status into separate 3d arrays
###Ensure a permutation dimension is retained when Nperm = 1
if(length(dim(PermOut)) == 3)
  dim(PermOut) = c(dim(PermOut),1)
InvasionResults <- array(PermOut[,,1,,drop=FALSE],dim = c(nrow(SDDprob),Ntimesteps,Nperm))
PopulationResults <- array(PermOut[,,2,,drop=FALSE],dim = c(nrow(SDDprob),Ntimesteps,Nperm))
ManagingResults <- array(PermOut[,,3,,drop=FALSE],dim = c(nrow(SDDprob),Ntimesteps,Nperm))
DetectedResults <- array(PermOut[,,4,,drop=FALSE],dim = c(nrow(SDDprob),Ntimesteps,Nperm))

###########################################################
###Save results for post-hoc analyses
###########################################################
###ModelName used to generate filenames
###Use standard format for ease of reading results to produce heat maps 
###and conduct post-hoc stats comparing managment settings/scenarios 
if(is.na(OutputDir) == T)
	OutputDir = ""
FileNameStem = paste0(OutputDir,ModelName)

###These are 3D arrays with dimensions (Nodes,Timesteps,Realisations)
saveRDS(ManagingResults, paste0(FileNameStem,"InfoLargeOut.rds"))
saveRDS(PopulationResults, paste0(FileNameStem,"PopulationLargeOut.rds"))
saveRDS(InvasionResults, paste0(FileNameStem,"InvasionLargeOut.rds"))
saveRDS(DetectedResults, paste0(FileNameStem,"DetectedLargeOut.rds"))

##########################################################
###Store annual node-level invasion probs for heat maps
###and estimation of invasion threat to other regions
##########################################################

InvasionProb = matrix(ncol = Ntimesteps, nrow = nrow(SDDprob))
for(timestep in 1:Ntimesteps)
{
TimestepData = matrix(InvasionResults[,timestep,,drop=FALSE],nrow=nrow(SDDprob),ncol=Nperm)
InvasionProb[,timestep] = rowSums(TimestepData)/Nperm
}
saveRDS(InvasionProb, paste0(FileNameStem,"InvasionProb.rds"))
if(DoPlots == T)
{
###########################################################
###Produce summary figs when processing completed
###########################################################

Title = ModelName


###Change in total population with time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

PopulationSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(PopulationSummary) = c("Realisation",   "Timestep",  "NodesInfested")

if(is.matrix(K) == TRUE)
    inv_K <- 1 / colSums(K)

for(perm in 1:Nperm)
{
PopulationData = PopulationResults[,,perm]
dim(PopulationData)
NodesInfested = colSums(PopulationData)*inv_K
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesInfested)
PopulationSummary = rbind(PopulationSummary,Results)
}


Filename = paste0(FileNameStem,"PopulationRaw.png")
png(Filename)
plot(PopulationSummary$Timestep,PopulationSummary$NodesInfested,ylim = c(0,1),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Total population (proportion of K)", main = Title)

for(perm in 1:Nperm)
{
Sub = PopulationSummary[PopulationSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesInfested,col  = perm)
}
dev.off()

Quantiles = as.data.frame(aggregate(PopulationSummary$NodesInfested, by = list(PopulationSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"PopulationSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,1), xlab = "Time since incursion detected (timesteps)",
ylab = "Total population (proportion of K)", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()

###Change in number of nodes infested with time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

InvasionSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(InvasionSummary) = c("Realisation",   "Timestep",  "NodesInfested")

for(perm in 1:Nperm)
{
InvasionData = InvasionResults[,,perm]
dim(InvasionData)
NodesInfested = colSums(InvasionData)
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesInfested)
InvasionSummary = rbind(InvasionSummary,Results)
}


Filename = paste0(FileNameStem,"InvasionRaw.png")
png(Filename)
plot(InvasionSummary$Timestep,InvasionSummary$NodesInfested,ylim = c(0,max(InvasionSummary$NodesInfested)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Number of nodes infested", main = Title)

for(perm in 1:Nperm)
{
Sub = InvasionSummary[InvasionSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesInfested,col  = perm)
}
dev.off()

Quantiles = as.data.frame(aggregate(InvasionSummary$NodesInfested, by = list(InvasionSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"InvasionSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Number of nodes infested", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


###Change in number of farms managing through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided

ManagingSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(ManagingSummary) = c("Realisation",   "Timestep",  "NodesManaging")

for(perm in 1:Nperm)
{
ManagingData = ManagingResults[,,perm]
dim(ManagingData)
NodesManaging = colSums(ManagingData)
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesManaging)
ManagingSummary = rbind(ManagingSummary,Results)
}



 
Filename = paste0(FileNameStem,"ManagingRaw.png")
png(Filename)
plot(ManagingSummary$Timestep,ManagingSummary$NodesManaging,ylim = c(0,max(ManagingSummary$NodesManaging)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes under management", main = Title)

for(perm in 1:Nperm)
{
Sub = ManagingSummary[ManagingSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesManaging,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(ManagingSummary$NodesManaging, by = list(ManagingSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"ManagingSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes under management", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


###Change in number of known extant infestations through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

DetectedSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(DetectedSummary) = c("Realisation",   "Timestep",  "NodesDetected")

for(perm in 1:Nperm)
{
DetectedData = DetectedResults[,,perm]
dim(DetectedData)
NodesDetected = colSums(DetectedData)
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesDetected)
DetectedSummary = rbind(DetectedSummary,Results)
}


Filename = paste0(FileNameStem,"DetectedRaw.png")
png(Filename)
plot(DetectedSummary$Timestep,DetectedSummary$NodesDetected,ylim = c(0,max(DetectedSummary$NodesDetected)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes pest detected", main = Title)

for(perm in 1:Nperm)
{
Sub = DetectedSummary[DetectedSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesDetected,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(DetectedSummary$NodesDetected, by = list(DetectedSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"DetectedSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes pest detected", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


###Change in proportion of extant infestations detected through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided
 
DetectedProportionSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(DetectedProportionSummary) = c("Realisation",   "Timestep",  "DetectedProportion")

for(perm in 1:Nperm)
{
InvasionData = InvasionResults[,,perm]
DetectedData = DetectedResults[,,perm]
NodesDetected = colSums(DetectedData)
NodesInvaded = colSums(InvasionData)
DetectedProportion = NodesDetected/NodesInvaded
DetectedProportion[is.na(DetectedProportion)==T] = 1
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,DetectedProportion)
DetectedProportionSummary = rbind(DetectedProportionSummary,Results)
}


Filename = paste0(FileNameStem,"DetectedProportionRaw.png")
png(Filename)
plot(DetectedProportionSummary$Timestep,DetectedProportionSummary$DetectedProportion,ylim = c(0,max(DetectedProportionSummary$DetectedProportion)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Proportion of infested nodes detected", main = Title)

for(perm in 1:Nperm)
{
Sub = DetectedProportionSummary[DetectedProportionSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$DetectedProportion,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(DetectedProportionSummary$DetectedProportion, by = list(DetectedProportionSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"DetectedProportionSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Proportion infested nodes detected", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()
}
}


################################################################
################################################################
###End of function
################################################################
################################################################
