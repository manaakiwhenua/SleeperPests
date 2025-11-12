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
local.dynamics.transition.matrix = function(nodetransition = NodeTransition,weights = Weights,sddprob = SDDprob, nodeenvestabprob = NodeEnvEstabProb,n0=N0,
                                            lddprob = LDDprob, lddrate = LDDrate,k_is_0 = K_is_0, nodeK = NodeK,nodepropaguleestablishment = NodePropaguleEstablishment,
                                            nodespreadreduction = NodeSpreadReduction,managing = Managing, MaxInteger = MaxInteger)
{
  
  n_pops <- nrow(n0)
  
  
  S <- ncol(n0)
  w <- if (is.null(weights)) rep(1, S) else weights
  
  
  # create new population matrix
  # to be updated
  n <- t(n0)
  # --- STEP 1: Propagule production (vectorized) ---
  if (is.list(nodetransition)) {
    # One matrix per population
    fec_means <- vapply(seq_len(n_pops), function(p) sum(nodetransition[[p]][1, ] * n[, p]), numeric(1))
  } else {
    # Single matrix for all populations
    fec_means <- as.numeric(t(nodetransition[1, ]) %*% n)
  }
  
  # Poisson draws, only where mean > 0
  propagules <- ifelse(fec_means > 0, rpois(n_pops, fec_means), 0) 
  
  # Pre-extract all transition matrices
  nodetransition_list <- if (is.list(nodetransition)) nodetransition else replicate(n_pops, nodetransition, simplify = FALSE)
  
  
  # Loop over stages (must remain sequential)
  for (s in S:2) {
    
    # Stage-specific vectors across populations
    N_prev <- n[s - 1, ]
    trans_prob <- vapply(nodetransition_list, function(Ap) Ap[s, s - 1], numeric(1))      # s-1 -> s
    surv_prob  <- vapply(nodetransition_list, function(Ap) Ap[s, s], numeric(1))          # stasis
    
    # Stage transitions (s-1 -> s)
    n_trans_candidates <- rbinom(n_pops, size = N_prev, prob = trans_prob)
    
    total_pop       <- colSums(n * w)
    slots_available <- pmax(0, floor((nodeK - total_pop) / w[s]))
    n_trans_limits  <- pmin(slots_available, n_trans_candidates)
    
    # Establishment probability per population
    est_prob <- 1 - exp(-nodepropaguleestablishment* n_trans_candidates)
    
    # Vectorized transition
    n_trans_actual <- rbinom(n_pops, size = n_trans_limits, prob = est_prob)
    
    
    # Add transitioned individuals to stage pops
    n[s, ] <- n[s, ] +  n_trans_actual
    
    # Remove transitioned individuals from previous stage
    n[s - 1, ] <- N_prev - n_trans_actual
    
    # Apply survival to all individuals currently in stage s
    # This includes both existing and newly transitioned individuals
    n[s, ] <- rbinom(n_pops, size = n[s, ], prob = surv_prob)
  }
  
  # --- STEP 3: Propagule dispersal (integer allocation) ---
  Pin <- numeric(n_pops)
  Qin <- numeric(n_pops)
  if (sum(propagules) > 0) {
    ###self-mediated spread
    total_p <- sum((propagules*(1-lddrate))* rowSums(sddprob))
    sddprob_matrix <- (propagules*(1-lddrate)) %*% sddprob
    if(floor(total_p) < MaxInteger)
      Pin <- as.numeric(t(rmultinom(1, size = floor(total_p), prob = sddprob_matrix)))
    else
      Pin <- floor(colSums(sweep(sddprob,1,floor(total_p) ,`*`)))
    ###human-mediated spread
    if (is.matrix(lddprob)){
      total_q <- sum((propagules*lddrate)* rowSums(lddprob))  
      lddprob_matrix <- (propagules*lddrate) %*% lddprob  
      if(floor(total_q) < MaxInteger)
        Qin <- as.numeric(t(rmultinom(1, size = floor(total_q*(1-nodespreadreduction*managing)), prob = lddprob_matrix)))  
      else
        Qin <- floor(colSums(sweep(lddprob,1,floor(total_q) ,`*`)))
    } 
  }
  
  
  # --- STEP 4: Recruitment after dispersal ---
  
  # Pre-calculate available slots per population
  slots <- pmax(0, floor((nodeK - colSums(n * w)) / w[1]))
  
  # Limit potential recruits to available slots
  max_recruits <- pmin(slots, Pin+Qin)
  
  # Establishment probability per population
  est_prob <- 1 - exp(-nodepropaguleestablishment * nodeenvestabprob * (Pin+Qin))
  
  # Vectorized recruitment (R automatically handles size = 0 or prob = 0)
  recruits <- rbinom(n_pops, size = max_recruits, prob = est_prob)
  
  # Add recruits to first stage
  n[1, ] <- n[1, ] + recruits
  
  
  #return updated population matrix
  
  return(t(n))
}


INApestMetaTansitionMatrix = function(
ModelName, #Name for storing results to file 
Nperm,                  #Number of permutations per parameter combination
Ntimesteps,                 #Simulation duration timesteps can be any length of time
Nstages,               #Number of stages in transition matrix
Weights,               #Weight for converting stage populations to total populations
Transition,            #Transition matrix (N stages x N stages), list of matrices (length = N nodes)
                       #or 4D array (Nstages x N stages x N nodes x N timesteps)
LocalDynamics = local.dynamics.transition.matrix, #Local population growth,dispersal and management function
DetectionProb,          #Vector of Per-individual detection probability for each stage, or matrix of probabilities per stage per node (e.g. farm) 
                        #or 3D array of probabilities per stage per node per year (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability. Can be single number or vector (nodes)
ManageProb,             #Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
ManageSD = NULL, #Option to provide standard deviation for management probability. Can be single number or vector (nodes)
MortalityProb,           ##Vector of Per-individual mortality probability for each stage, or matrix of probabilities per stage per node (e.g. farm) 
                         #or 3D array of probabilities per stage per node per year (must be between 0 and 1)
MortalitySD = NULL, #Option to provide standard deviation for mortality probability. Can be single number or vector (nodes)
SpreadReduction,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
SpreadReductionSD = NULL, #Option to provide standard deviation for spread reduction probability can be single number or vector (nodes)
InitialPopulation = NA,        #matrix (nodes x stages) of population sizes at start of simulations
InitBioP = NA,		#Proportion of nodes infested at start of simulations
InvasionRisk = NA,           #Vector of probabilities of invasion from external sources
InitialInfo = NA,        #Vector or of nodes with information at start of simulations
InitInfoP = NA,		#Proportion of nodes with information at start of simulations
ExternalInfoProb = 0.001,           #Vector of probabilities of communication from external sources
EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
K,		       #Population carrying capacity - vector (nodes)
PropaguleEstablishment, #Propagules establishment probability
IncursionStartPop=NA,      #option to set population size for new incursions
SDDprob,                   #Natural disperal probability between each pair of nodes
SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
LDDprob = NA,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
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
###2) Allow provision of natural mortality rate to permit extinction of local populations (may happen in climates where R0 is very low?)
  
###Max integer for propagule dispersal using rmultinom
MaxInteger <- .Machine$integer.max  
  
# pre-evaluate some variables for efficiency
if(is.matrix(K) == FALSE)
{
K_is_0 <- K<=0
inv_K <- 1 / sum(K)
NodeK = K
}

###If carrying capacity provided as matrix assign values from first timestep for population initialisation
if(is.matrix(K) == TRUE)
  {
  K_is_0 <- K[,1]<=0
  inv_K <- 1 / sum(K[,1])
  NodeK = K[,1] 
}    
  
if(is.matrix(PropaguleEstablishment) == FALSE)
  NodePropaguleEstablishment = PropaguleEstablishment

if(is.matrix(EnvEstabProb) == F)
  NodeEnvEstabProb <- EnvEstabProb

###Declare array tracking population size
###of individual nodes in each timestep of each realisation
PopulationResults = array(dim = c(nrow(SDDprob),Ntimesteps,Nperm))

###Declare array tracking population size
###of individual nodes in each timestep of each realisation
PopulationStageResults = array(dim = c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))

###Declare array tracking invasion status
###of individual nodes in each timestep of each realisation
InvasionResults = array(dim = c(nrow(SDDprob),Ntimesteps,Nperm))


###Declare array tracking detection status 
###of individual nodes in each timestep of each realisation
DetectedResults = InvasionResults

###Declare array for tracking management adoption status 
###of individual nodes in each timestep of each realisation
###This is a measure of potential disruption to farm businesses
###or ongoing management burden (surveillance and removal)
###for publicly-owned lands
ManagingResults = InvasionResults


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

n_nodes <- nrow(SDDprob)


# --- Detection SD ---
if (is.null(DetectionSD)) {
  if (is.matrix(DetectionProb)) {
    # Calculate mean per stage (column)
    stage_means <- colMeans(DetectionProb, na.rm = TRUE)
    # Replicate per node
    DetectionSD <- matrix(stage_means / 10, nrow = n_nodes, ncol = ncol(DetectionProb), byrow = TRUE)
  } else {
    # If scalar or vector, repeat across all nodes
    DetectionSD <- matrix(mean(DetectionProb, na.rm = TRUE) / 10,
                          nrow = n_nodes, ncol = Nstages)
  }
}

# --- Mortality SD ---
if (is.null(MortalitySD)) {
  if (is.matrix(MortalityProb)) {
    # Calculate mean per stage (column)
    stage_means <- colMeans(MortalityProb, na.rm = TRUE)
    # Replicate per node
    MortalitySD <- matrix(stage_means / 10, nrow = n_nodes, ncol = ncol(MortalityProb), byrow = TRUE)
  } else {
    # If scalar or vector, repeat across all nodes
    MortalitySD <- matrix(mean(MortalityProb, na.rm = TRUE) / 10,
                          nrow = n_nodes, ncol = Nstages)
  }
}




###########################################################
###Start of simulation
###########################################################
    
for (perm in 1:Nperm) 
{ 
###Assign initial infestations according either to "InitialInvasion" binary vector OR
###"InvasionRisk" probabilities and/or initial proportion of nodes infested ("InitBioP") OR
###just "InitBioP" if neither "InitialInvasion" or "InvasionRisk" supplied by user
InitBio = matrix(nrow = nrow(SDDprob),ncol = Nstages)
if(nrow(InitialPopulation) != nrow(SDDprob))
{
 if(length(InvasionRisk) == nrow(SDDprob))
        {
 	if(is.na(InitBioP) == F)
	  Infested = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitBioP),prob = InvasionRisk)
	if(is.na(InitBioP) == T)
          {
	  Infested = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk)
          Infested = which(Infested == 1)
	  } 
	}
 if(length(InvasionRisk) != nrow(SDDprob))
        {
 	if(is.matrix(InvasionRisk) == F)
	 Infested = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitBioP))
        if(is.matrix(InvasionRisk) == T)
          {
          Infested = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk[,1])
          Infested = which(Infested == 1)
	  }
	}
if(is.na(IncursionStartPop) == T) 
	InitBio[Infested,1] = 1
if(is.na(IncursionStartPop) == F) 
	InitBio[Infested,1] = IncursionStartPop
}

if(nrow(InitialPopulation) == nrow(SDDprob))
	InitBio = InitialPopulation

# --- Ensure initial population (weighted) does not exceed carrying capacity ---

if (exists("Weights") && length(Weights) == Nstages) {
  # Calculate weighted population per node
  weighted_pop <- InitBio %*% Weights  # (n_nodes x 1)
  
  # Identify nodes exceeding K
  overcap <- which(weighted_pop > NodeK)
  
  if (length(overcap) > 0) {
    # Scale down all stage values proportionally
    scale_factor <- NodeK[overcap] / weighted_pop[overcap]
    InitBio[overcap, ] <- InitBio[overcap, , drop = FALSE] * scale_factor
  }
  
} else {
  # Fallback: simple elementwise comparison (if Weights missing)
  for (i in seq_len(nrow(InitBio))) {
    if (any(InitBio[i, ] > NodeK[i])) {
      # Clip all stages proportionally to match K
      total_pop <- sum(InitBio[i, ])
      if (total_pop > 0) {
        InitBio[i, ] <- InitBio[i, ] * (NodeK[i] / total_pop)
      }
    }
  }
}

# Ensure all stage populations are integers
InitBio <- floor(InitBio)

# initialise the population
N <- InitBio
if(sum(N) == 0 && OngoingExternalInvasion == F)
  warning("No initial populations and no future external invasions")

###Select nodes with information at start of simulation  according either to "InitialInfo" binary vector OR
###"ExternalInfoProb" probabilities and/or initial proportion of nodes with information ("InitInfoP") OR
###just "InitInfoP" if neither "InitialInfo" or "ExternalInfoProb" supplied by user.
###If no initial info variables provided, no nodes have info at start of simulations
# Robust InitInfo setup
n_nodes <- nrow(SDDprob)  # number of nodes
InitInfo <- rep(0, n_nodes)  # default: no nodes initially informed

# Only proceed if any initial info or probabilities are provided
if(!all(is.na(InitialInfo)) || !is.na(InitInfoP) || !all(is.na(ExternalInfoProb))) {
  
  # If InitialInfo is valid length, just use it
  if(length(InitialInfo) == n_nodes) {
    InitInfo <- InitialInfo
    
  } else {
    # Ensure InitInfoP is numeric and in [0,1]
    if(is.na(InitInfoP)) InitInfoP <- 0
    
    # Ensure ExternalInfoProb is numeric of length n_nodes
    if(is.na(sum(ExternalInfoProb))) ExternalInfoProb <- rep(1, n_nodes)
    if(length(ExternalInfoProb) != n_nodes) ExternalInfoProb <- rep(1, n_nodes)
    
    # Determine number of nodes to sample
    n_sample <- ceiling(n_nodes * InitInfoP)
    
    # Only sample if n_sample > 0
    if(n_sample > 0) {
      Info <- sample(1:n_nodes, size = n_sample, prob = ExternalInfoProb)
      InitInfo[Info] <- 1
    }
  }
}





# --- Determine dimensions ---
n_nodes <- nrow(SDDprob)

###Randomly assign detection probabilities
# --- Case 1: Single value or vector (stages only) ---
if (!is.array(DetectionProb) && (length(DetectionProb) == 1 || length(DetectionProb) == Nstages)) {
  
  NodeDetectionProb <- array(NA, dim = c(n_nodes, Nstages, Ntimesteps))
  
  for (s in 1:Nstages) {
    for (t in 1:Ntimesteps) {
      NodeDetectionProb[, s, t] <- pmax(0, pmin(1,
                                                rnorm(n_nodes,
                                                      mean = if (length(DetectionProb) == 1) DetectionProb else DetectionProb[s],
                                                      sd   = if (is.matrix(DetectionSD)) DetectionSD[, s]   # <--- per node & stage
                                                      else if (length(DetectionSD) == 1) DetectionSD
                                                      else DetectionSD[s]
                                                )
      ))
    }
  }
}


# --- Case 2: Matrix (nodes × stages) ---
if (is.matrix(DetectionProb) && all(dim(DetectionProb) == c(n_nodes, Nstages))) {
  
  NodeDetectionProb <- array(NA, dim = c(n_nodes, Nstages, Ntimesteps))
  
  for (t in 1:Ntimesteps) {
    NodeDetectionProb[, , t] <- pmax(0, pmin(1,
                                             matrix(
                                               rnorm(n_nodes * Nstages,
                                                     mean = as.vector(DetectionProb),
                                                     sd   = as.vector(DetectionSD)),
                                               nrow = n_nodes, ncol = Nstages)
    ))
  }
}

# --- Case 3: 3D array (nodes × stages × timesteps) ---
if (length(dim(DetectionProb)) == 3 && all(dim(DetectionProb)[1:2] == c(n_nodes, Nstages))) {
  
  NodeDetectionProb <- array(NA, dim = dim(DetectionProb))
  
  for (t in 1:Ntimesteps) {
    NodeDetectionProb[, , t] <- pmax(0, pmin(1,
                                             matrix(
                                               rnorm(n_nodes * Nstages,
                                                     mean = as.vector(DetectionProb[, , t]),
                                                     sd   = as.vector(DetectionSD)),
                                               nrow = n_nodes, ncol = Nstages)
    ))
  }
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
      NodeSpreadReduction = rnorm(SpreadReduction,SpreadReductionSD,n = nrow(SDDprob))
      NodeSpreadReduction[NodeSpreadReduction<0] = 0
      NodeSpreadReduction[NodeSpreadReduction>1] = 1
      }

###Randomly assign mortality probability when management applied
# --- Determine dimensions ---
n_nodes <- nrow(SDDprob)

# --- Case 1: Single value or vector (stages only) ---
if (!is.array(MortalityProb) && (length(MortalityProb) == 1 || length(MortalityProb) == Nstages)) {
  
  NodeMortalityProb <- array(NA, dim = c(n_nodes, Nstages, Ntimesteps))
  
  for (s in 1:Nstages) {
    for (t in 1:Ntimesteps) {
      NodeMortalityProb[, s, t] <- pmax(0, pmin(1,
                                                rnorm(n_nodes,
                                                      mean = if (length(MortalityProb) == 1) MortalityProb else MortalityProb[s],
                                                      sd   = if (is.matrix(MortalitySD)) MortalitySD[, s]   # <--- per node & stage
                                                      else if (length(MortalitySD) == 1) MortalitySD
                                                      else MortalitySD[s]
                                                )
      ))
    }
  }
}


# --- Case 2: Matrix (nodes × stages) ---
if (is.matrix(MortalityProb) && all(dim(MortalityProb) == c(n_nodes, Nstages))) {
  
  NodeMortalityProb <- array(NA, dim = c(n_nodes, Nstages, Ntimesteps))
  
  for (t in 1:Ntimesteps) {
    NodeMortalityProb[, , t] <- pmax(0, pmin(1,
                                             matrix(
                                               rnorm(n_nodes * Nstages,
                                                     mean = as.vector(MortalityProb),
                                                     sd   = as.vector(MortalitySD)),
                                               nrow = n_nodes, ncol = Nstages)
    ))
  }
}

# --- Case 3: 3D array (nodes × stages × timesteps) ---
if (length(dim(MortalityProb)) == 3 && all(dim(MortalityProb)[1:2] == c(n_nodes, Nstages))) {
  
  NodeMortalityProb <- array(NA, dim = dim(MortalityProb))
  
  for (t in 1:Ntimesteps) {
    NodeMortalityProb[, , t] <- pmax(0, pmin(1,
                                             matrix(
                                               rnorm(n_nodes * Nstages,
                                                     mean = as.vector(MortalityProb[, , t]),
                                                     sd   = as.vector(MortalitySD)),
                                               nrow = n_nodes, ncol = Nstages)
    ))
  }
}






###Populate invasion status vector ahead of timestep loop
Invaded <- ifelse(rowSums(InitBio) > 0, 1, 0)

###Probability of info at start of simulation depends on
###Presence of pest and detection probability
###Select nodes that have detected infestation 




# Calculate probability of detection per node
prob_detect <- 1 - apply((1 - NodeDetectionProb[,,1])^InitBio, 1, prod)

# Draw initial detection (0/1) per node
InitDetection <- rbinom(n = nrow(InitBio), size = 1, prob = prob_detect)


###Add detections to nodes which already have info (e.g. pre-emptive control and hygiene measures)
InitInfo[InitInfo == 0] = InitDetection[InitInfo == 0]

###Populate information status vector ahead of timestep loop
HaveInfo = InitInfo


  # run simulation
  for (timestep in 1:Ntimesteps) 
    { 
    ###Print progress
    cat("\r", "Realisation ", perm, "Timestep ", timestep, "...")
 
    ###Allow for variation in recruit establishment through time
    ###e.g.  climate change predictions
    ###Note: could be done outside loop, but would take heaps of memory to store 
    if(is.matrix(EnvEstabProb) == T)
      NodeEnvEstabProb <- EnvEstabProb[,timestep]
    
    if (is.array(Transition) && length(dim(Transition)) == 3) {
      NodeTransition <- Transition[,,timestep]
      } else if (is.list(Transition)) {
        NodeTransition <- Transition
      } else {
      NodeTransition <- Transition
      }
    
    
    ###If carrying capacity provided as matrix assign values for relevant timestep
    if(is.matrix(K) == TRUE)
      {
      K_is_0 <- K[,timestep]<=0
      inv_K <- 1 / sum(K[,timestep])
      NodeK = K[,timestep] 
      }  
    
    if(is.matrix(PropaguleEstablishment) == TRUE)
      NodePropaguleEstablishment = PropaguleEstablishment[,timestep]

      
  
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
  
  ###Assign management status to nodes   
  ###Management is only applied to nodes which have information
  ###i.e. where pest has been detected or following communication of information
  ###from neighbouring infested farms 
  Managing = rbinom(1:nrow(SDDprob),size = 1,prob = NodeManageProb*HaveInfo)
  
  ###Management is only applied to nodes which have information
  ###i.e. where pest has been detected or following communication of information
  ###from neighbouring infested farms 
  Managing = Managing*HaveInfo
  
  ###Identify nodes with known extant infestations 
  Detected = Invaded*HaveInfo
  
  ###Apply management
  # N: current population matrix (nodes x stages)
  # Managing: vector of management adoption per node (0/1) for current timestep
  # NodeMortalityProb: array (nodes x stages x timesteps)
  # timestep: current timestep index
  
  # Apply management-driven mortality
  N0 <- matrix(NA, nrow = nrow(N), ncol = ncol(N))
  for (s in 1:Nstages) {
    N0[, s] <- rbinom(n = nrow(N), size = N[, s], prob = 1 - (NodeMortalityProb[, s,t] * Managing))
  }
  
  # Update population matrix
  if (sum(N0) <= 0) {
    N <- N0
  } else if (sum(N0) > 0) {
    N <- N0
  }
  
  if(sum(N0)>0 ) 
  {
      
  N <- LocalDynamics(nodetransition = NodeTransition,weights = Weights,sddprob = SDDprob, nodeenvestabprob = NodeEnvEstabProb,n0=N0,
                     lddprob = LDDprob, lddrate = LDDrate,k_is_0 = K_is_0, nodeK = NodeK,nodepropaguleestablishment = NodePropaguleEstablishment,
                     nodespreadreduction = NodeSpreadReduction,managing = Managing,MaxInteger=MaxInteger)
  } 
 ###Update info vector for any info spread (if SEAM supplied)
 ###Note once nodes obtain info they always have info (only zero values updated)
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
	  N[,1] = N[,1]+ExternalInvasion
  if(is.na(IncursionStartPop) == F) 
	  N[,1] = N[,1]+ExternalInvasion*IncursionStartPop
  if (exists("Weights") && length(Weights) == Nstages) {
    
    
    # Calculate weighted population per node
    weighted_pop <- N %*% Weights  # (n_nodes x 1)
    # Identify nodes exceeding K
    overcap <- which(weighted_pop > NodeK)
    
    if (length(overcap) > 0) {
      # Scale down all stage values proportionally
      scale_factor <- NodeK[overcap] / weighted_pop[overcap]
      N[overcap, ] <- N[overcap, , drop = FALSE] * scale_factor
    }
    
  } else {
    # Fallback: simple elementwise comparison (if weights missing)
    for (i in seq_len(nrow(N))) {
      if (any(N[i, ] > NodeK[i])) {
        # Clip all stages proportionally to match K
        total_pop <- sum(N[i, ])
        if (total_pop > 0) {
          N[i, ] <- N[i, ] * (NodeK[i] / total_pop)
        }
      }
    }
  }
  
  # Ensure all stage populations are integers
  N <- floor(N)
  
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
 Invaded = ifelse(rowSums(N)>0,1,0)
 
 ###Record nodes adopting management
 ManagingResults[,timestep,perm] = Managing
  
 ###Record infested nodes
 InvasionResults[,timestep,perm] = Invaded

 ###Record populations
 PopulationResults[,timestep,perm] = N %*% Weights
 
 ###Record stage populations
 PopulationStageResults[,,timestep,perm] = N

 ###Select new nodes where infestation detected
 # Probability of detection per node per stage
 DetectionProbPerStage <- 1 - (1 - NodeDetectionProb[,,t])^N  # element-wise OK: both nodes x stages
 
 # Combine stages to get probability of detecting at least one stage
 ProbDetectNode <- 1 - apply(1 - DetectionProbPerStage, 1, prod)  # multiply across stages
 
 # Sample new detections per node
 NewHaveInfo <- rbinom(n = nrow(SDDprob), size = 1, prob = ProbDetectNode)
 
 
 ###Add newly detected infestations to info vector
 ###Note once nodes obtain info they always have info (only zero values updated)
 HaveInfo[HaveInfo==0] = NewHaveInfo[HaveInfo==0]  
 
 ###Record detection status
 DetectedResults[,timestep,perm] = HaveInfo*Invaded 
 }
}
###########################################################
###End of Simulation
###########################################################

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
saveRDS(PopulationStageResults, paste0(FileNameStem,"PopulationStageLargeOut.rds"))
saveRDS(InvasionResults, paste0(FileNameStem,"InvasionLargeOut.rds"))
saveRDS(DetectedResults, paste0(FileNameStem,"DetectedLargeOut.rds"))

##########################################################
###Store annual node-level invasion probs for heat maps
###and estimation of invasion threat to other regions
##########################################################

InvasionProb = matrix(ncol = Ntimesteps, nrow = nrow(SDDprob))
for(timestep in 1:Ntimesteps)
{
TimestepData = InvasionResults[,timestep,]
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
