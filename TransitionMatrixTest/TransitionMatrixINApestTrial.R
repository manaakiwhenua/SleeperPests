setwd(r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\TransitionMatrix]")
source(r"[C:\Users\MasonN\OneDrive - MWLR\Documents\GitHub\SleeperPests\INApest INApestMeta core function code/INApestMetaTransitionMatrixParallel.R]")
source(r"[C:\Users\MasonN\OneDrive - MWLR\Documents\GitHub\SleeperPests\INApest INApestMeta core function code/INApestMetaTransitionMatrix.R]")
# -----------------------------
# Basic example setup
# -----------------------------
generate_transition_matrix_dd <- function(n_subadult = 3,
                                          surv_subadult = 0.8,
                                          surv_adult = 0.9,
                                          fecundity = 2,
                                          surv_to_adulthood = 0.5) {
  n_stages <- n_subadult + 1
  
  # Special case: 0 subadult stages
  if(n_subadult == 0) {
    A <- matrix(0, nrow = 2, ncol = 1)
    A[1, 1] <- fecundity     # top row = fecundity
    A[2, 1] <- surv_to_adulthood    # bottom row = adult survival
    return(A)
  }
  
  # Normal case
  A <- matrix(0, nrow = n_stages, ncol = n_stages)
  
  # Adult fecundity
  A[1, n_stages] <- fecundity
  
  # Subadult transitions
  g <- surv_to_adulthood^(1 / n_subadult)
  
  for(i in 1:n_subadult) {
    A[i, i] <- surv_subadult       # survival/stasis
    A[i+1, i] <- g                 # transition to next stage
  }
  
  # Adult survival
  A[n_stages, n_stages] <- surv_adult
  
  return(A)
}

####Generate stage weights
generate_stage_weights <- function(Nstages, AdultSize = 1) {
  # The first stage is smallest, last stage = AdultSize
  # Smallest stage = AdultSize / 2^(Nstages - 1)
  smallest <- AdultSize / 2^(Nstages - 1)
  weights <- smallest * 2^((0:(Nstages - 1)))
  return(weights)
}


n_nodes <- 3
n_stages <- 3
Nperm <- 30
Ntimesteps <- 100
K <- c(1000,10000,100000)  # Carrying capacity per node

weights = generate_stage_weights(Nstages = 3, AdultSize = 10)

# Initial population (nodes x stages)
InitialPopulation <- matrix(c(
  100, 10, 20,
  0, 0, 0,
  0, 0, 0
), nrow = n_stages, byrow = TRUE)
# Transition matrix (single matrix for all nodes)
A_input <- matrix(c(
  0.7, 0.2, 0.0,
  0.1, 0.6, 0.2,
  0.0, 0.1, 0.8
), nrow = n_nodes, byrow = TRUE)

# Detection and management probabilities (can be vectors of length stages)
DetectionProb <- rep(0, n_stages)
ManageProb <- rep(0, n_nodes)

# Mortality probability (vector length stages)
MortalityProb <- rep(0.1, n_stages)

# Carrying capacity for each node
K <- c(1000,10000,100000)

# Example SDDprob (dispersal matrix)
SDDprob <- matrix(1/n_nodes, nrow = n_nodes, ncol = n_nodes)

# Run
results <- INApestMetaTansitionMatrix(
  ModelName = "Subadult_2_",
  Nperm = Nperm,
  Ntimesteps = Ntimesteps,
  Nstages = n_stages,
  Weights = weights,
  Transition = A_input,
  LocalDynamics = local.dynamics.transition.matrix,
  DetectionProb = DetectionProb,
  ManageProb = ManageProb,
  MortalityProb = MortalityProb,
  SpreadReduction = 0,
  InitialPopulation = InitialPopulation,
  K = K,
  PropaguleEstablishment = 50/K,
  SDDprob = SDDprob,
  DoPlots = TRUE
)

# Run
INApestMetaTansitionMatrixParallel(
  ModelName = "Subadult_2_parallel",
  Nperm = Nperm,
  Ntimesteps = Ntimesteps,
  Nstages = n_stages,
  Weights = weights,
  Transition = A_input,
  LocalDynamics = local.dynamics.transition.matrix,
  DetectionProb = DetectionProb,
  ManageProb = ManageProb,
  MortalityProb = MortalityProb,
  SpreadReduction = 0,
  InitialPopulation = InitialPopulation,
  K = K,
  PropaguleEstablishment = 50/K,
  SDDprob = SDDprob,
  DoPlots = TRUE
)

###Try with only 1 subadult stage
A1 <- generate_transition_matrix_dd(n_subadult = 1, surv_subadult = 0.8,
                                    surv_adult = 0.9, fecundity = 2,
                                    surv_to_adulthood = 0.5)
weights = generate_stage_weights(Nstages = 2, AdultSize = 10)
DetectionProb <- rep(0, 2)
MortalityProb <- rep(0.1, 2)

# Run
INApestMetaTansitionMatrix(
  ModelName = "Subadult_1",
  Nperm = Nperm,
  Ntimesteps = Ntimesteps,
  Nstages = 2,
  Weights = weights,
  Transition = A1,
  LocalDynamics = local.dynamics.transition.matrix,
  DetectionProb = DetectionProb,
  ManageProb = 0,
  MortalityProb = MortalityProb,
  SpreadReduction = 0,
  InitialPopulation = InitialPopulation[,1:2],
  K = K,
  PropaguleEstablishment = 50/K,
  SDDprob = SDDprob,
  DoPlots = TRUE
)

INApestMetaTansitionMatrixParallel(
  ModelName = "Subadult_1_parallel",
  Nperm = Nperm,
  Ntimesteps = Ntimesteps,
  Nstages = 2,
  Weights = c(0.1,10),
  Transition = A1,
  LocalDynamics = local.dynamics.transition.matrix,
  DetectionProb = DetectionProb,
  ManageProb = ManageProb,
  MortalityProb = MortalityProb,
  SpreadReduction = 0,
  InitialPopulation = InitialPopulation[,1:2],
  K = K,
  PropaguleEstablishment = 50/K,
  SDDprob = SDDprob,
  DoPlots = TRUE
)

###Try with only 9 subadult stages
A9 <- generate_transition_matrix_dd(n_subadult = 9, surv_subadult = 0.8,
                                    surv_adult = 0.9, fecundity = 2,
                                    surv_to_adulthood = 0.5)
InitialPopulation <- matrix(c(
  100, 10, 20,
  0, 0, 0,
  0, 0, 0
), nrow = n_stages, byrow = TRUE)
AddYears <- matrix(0,nrow = 3,ncol = 7)
InitialPopulation <- cbind(AddYears,InitialPopulation)

DetectionProb <- rep(0, 10)
MortalityProb <- rep(0.1, 10)
weights = generate_stage_weights(Nstages = 10, AdultSize = 10)
INApestMetaTansitionMatrix(
  ModelName = "Subadult_9",
  Nperm = Nperm,
  Ntimesteps = Ntimesteps,
  Nstages = 10,
  Weights = weights,
  Transition = A9,
  LocalDynamics = local.dynamics.transition.matrix,
  DetectionProb = DetectionProb,
  ManageProb = ManageProb,
  MortalityProb = MortalityProb,
  SpreadReduction = 0,
  InitialPopulation = InitialPopulation,
  K = K,
  PropaguleEstablishment = 50/K,
  SDDprob = SDDprob,
  DoPlots = TRUE
)

INApestMetaTansitionMatrixParallel(
  ModelName = "Subadult_9_parallel",
  Nperm = Nperm,
  Ntimesteps = Ntimesteps,
  Nstages = 10,
  Weights = weights,
  Transition = A9,
  LocalDynamics = local.dynamics.transition.matrix,
  DetectionProb = DetectionProb,
  ManageProb = ManageProb,
  MortalityProb = MortalityProb,
  SpreadReduction = 0,
  InitialPopulation = InitialPopulation,
  K = K,
  PropaguleEstablishment = 50/K,
  SDDprob = SDDprob,
  DoPlots = TRUE
)

