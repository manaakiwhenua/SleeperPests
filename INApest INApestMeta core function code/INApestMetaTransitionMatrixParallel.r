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
# Local stage dynamics with passive- or active-propagule dispersal.
#
# Backward-compatible options added for animal applications:
#   * BlockedTransitionMortality controls mortality among individuals that
#     attempt a stage transition but fail to obtain a target-stage slot.
#   * nodepropaguleestablishment = 1 lets active dispersers search every
#     available slot after they select a destination cell.
#   * DispersalDensityFactor > 0 biases movement toward cells that were less
#     occupied when propagules were produced, before mortality or transitions.
local.dynamics.transition.matrix <- function(
    nodetransition = NodeTransition,
    weights = Weights,
    sddprob = SDDprob,
    nodeenvestabprob = NodeEnvEstabProb,
    n0 = N0,
    lddprob = LDDprob,
    lddrate = LDDrate,
    nodeK = NodeK,
    node.seedbankK = NodeSeedbankK,
    nodepropaguleestablishment = NodePropaguleEstablishment,
    nodespreadreduction = NodeSpreadReduction,
    managing = Managing,
    MaxInteger = MaxInteger,
    ApplyFootprintToTransitions = FALSE,
    BlockedTransitionMortality = 0,
    DispersalDensityFactor = 0
) {
  n_pops <- nrow(n0)
  S <- ncol(n0)
  n <- t(n0)
  
  # Stage weights may be supplied in the original form as one value per
  # stage, or as an n_pops x S matrix for node-specific stage weights.
  # A vector is expanded to an identical row for every population, exactly
  # preserving previous scale-invariance and density-dispersal applications.
  if (is.null(weights)) {
    w <- matrix(1, nrow = n_pops, ncol = S)
  } else if (is.matrix(weights)) {
    if (!identical(dim(weights), c(n_pops, S))) {
      stop("A weights matrix must have n_pops rows and S columns")
    }
    w <- weights
  } else {
    if (length(weights) != S) {
      stop("A weights vector must contain one value per stage")
    }
    w <- matrix(
      rep(weights, each = n_pops),
      nrow = n_pops,
      ncol = S
    )
  }
  
  if (any(!is.finite(w)) || any(w <= 0)) {
    stop("All stage weights must be finite and greater than zero")
  }
  
  # Recycle node-level inputs once rather than repeatedly inside loops.
  nodeK <- rep_len(nodeK, n_pops)
  node.seedbankK <- rep_len(node.seedbankK, n_pops)
  nodeenvestabprob <- rep_len(nodeenvestabprob, n_pops)
  nodepropaguleestablishment <- rep_len(nodepropaguleestablishment, n_pops)
  
  # Probability that a progression candidate dies when it fails to secure
  # a target-stage slot. A scalar applies to every source stage and node; a
  # vector of length S - 1 supplies one probability for each source stage;
  # and an n_pops x (S - 1) matrix permits node-by-stage values. The default
  # zero exactly preserves the former treatment of blocking as stasis.
  if (is.matrix(BlockedTransitionMortality)) {
    if (!identical(dim(BlockedTransitionMortality), c(n_pops, S - 1L))) {
      stop(
        "BlockedTransitionMortality must be scalar, length S - 1, or an n_pops x (S - 1) matrix"
      )
    }
    blocked_transition_mortality <- BlockedTransitionMortality
  } else if (length(BlockedTransitionMortality) == 1L) {
    blocked_transition_mortality <- matrix(
      BlockedTransitionMortality,
      nrow = n_pops,
      ncol = S - 1L
    )
  } else if (length(BlockedTransitionMortality) == S - 1L) {
    blocked_transition_mortality <- matrix(
      rep(BlockedTransitionMortality, each = n_pops),
      nrow = n_pops,
      ncol = S - 1L
    )
  } else {
    stop(
      "BlockedTransitionMortality must be scalar, length S - 1, or an n_pops x (S - 1) matrix"
    )
  }
  
  if (any(!is.finite(blocked_transition_mortality)) ||
      any(blocked_transition_mortality < 0 |
          blocked_transition_mortality > 1)) {
    stop("BlockedTransitionMortality probabilities must be between 0 and 1")
  }
  
  # ------------------------------------------------------------
  # STEP 1: Propagule production and transition extraction
  # ------------------------------------------------------------
  
  transition_is_list <- is.list(nodetransition)
  
  if (transition_is_list) {
    if (length(nodetransition) != n_pops) {
      stop("A nodetransition list must contain one matrix per population")
    }
    
    fecundity <- matrix(
      vapply(
        nodetransition,
        function(A) A[1, -1],
        numeric(S - 1)
      ),
      nrow = S - 1L,
      ncol = n_pops
    )
    
    transition_probabilities <- matrix(
      vapply(
        nodetransition,
        function(A) A[cbind(2:S, 1:(S - 1))],
        numeric(S - 1)
      ),
      nrow = S - 1L,
      ncol = n_pops
    )
    
    stasis_probabilities <- matrix(
      vapply(
        nodetransition,
        function(A) diag(A)[1:(S - 1)],
        numeric(S - 1)
      ),
      nrow = S - 1L,
      ncol = n_pops
    )
    
    terminal_survival_prob <- vapply(
      nodetransition,
      function(A) A[S, S],
      numeric(1)
    )
    
    fec_means <- colSums(fecundity * n[-1, , drop = FALSE])
    mother_counts <- colSums(
      n[-1, , drop = FALSE] * (fecundity > 0)
    )
  } else {
    fecundity <- nodetransition[1, -1]
    transition_probabilities <-
      nodetransition[cbind(2:S, 1:(S - 1))]
    stasis_probabilities <- diag(nodetransition)[1:(S - 1)]
    terminal_survival_prob <- nodetransition[S, S]
    
    fec_means <- as.numeric(
      fecundity %*% n[-1, , drop = FALSE]
    )
    
    reproductive_stages <- which(fecundity > 0) + 1L
    mother_counts <- if (length(reproductive_stages)) {
      colSums(n[reproductive_stages, , drop = FALSE])
    } else {
      numeric(n_pops)
    }
  }
  
  if (any(!is.finite(fec_means)) || any(fec_means < 0)) {
    stop("Propagule-production means must be finite and non-negative")
  }
  
  propagules <- rpois(n_pops, fec_means)
  
  # ------------------------------------------------------------
  # STEP 1b: Optional density-dependent animal dispersal
  # ------------------------------------------------------------
  #
  # Reweight dispersal before mortality and stage transitions. Propagules
  # produced in STEP 1 therefore respond to the same pre-transition census
  # that produced them. For the hornet model, gynes can consequently respond
  # to main nests that are present during gyne production, even though those
  # annual nests subsequently die during STEP 2.
  #
  # NA or 0 leaves the original dispersal matrix unchanged. For any positive
  # value, DispersalDensityFactor is the exponent alpha. The source-cell
  # diagonal is treated like every other destination. Probability absent from
  # an original row remains an outside-landscape option, permitting export.
  
  if (length(DispersalDensityFactor) != 1L) {
    stop("DispersalDensityFactor must be a single value")
  }
  
  density_dependent_dispersal <-
    !is.na(DispersalDensityFactor) &&
    DispersalDensityFactor != 0
  
  if (density_dependent_dispersal) {
    if (!is.finite(DispersalDensityFactor) ||
        DispersalDensityFactor <= 0) {
      stop(
        "DispersalDensityFactor must be finite and positive, or 0/NA to disable it"
      )
    }
    
    if (any(propagules > 0)) {
      # Weighted occupancy of stages 2:S in the pre-transition census.
      pre_transition_population <- rowSums(
        n0[, 2:S, drop = FALSE] * w[, 2:S, drop = FALSE]
      )
      
      # Cells with K <= 0 cannot attract dispersers. Calculating the ratio
      # only for positive capacities also avoids division by zero.
      relative_occupancy <- rep(1, n_pops)
      positive_capacity <- nodeK > 0
      relative_occupancy[positive_capacity] <- pmin(
        1,
        pmax(
          0,
          pre_transition_population[positive_capacity] /
            nodeK[positive_capacity]
        )
      )
      
      # Zero occupancy gives attractiveness one; occupancy at or above K
      # gives attractiveness zero.
      cell_attractiveness <-
        pmax(0, 1 - relative_occupancy) ^ DispersalDensityFactor
      
      # The outside landscape retains constant attractiveness one.
      base_export_probability <- pmax(0, 1 - rowSums(sddprob))
      
      # R stores matrices column-wise. Repeating each destination value
      # n_pops times scales columns without sweep().
      cell_choice_weights <-
        sddprob * rep(cell_attractiveness, each = n_pops)
      
      total_choice_weights <-
        rowSums(cell_choice_weights) + base_export_probability
      
      # A row-length vector is recycled down every column, normalising the
      # destination and export choices without sweep() or matrix subsetting.
      row_multiplier <- numeric(n_pops)
      valid_rows <- total_choice_weights > 0
      row_multiplier[valid_rows] <- 1 / total_choice_weights[valid_rows]
      
      # If all reachable cells are full and no original export probability
      # exists, the adjusted row remains zero and propagules find no modelled
      # destination.
      sddprob <- cell_choice_weights * row_multiplier
    }
  }
  
  # ------------------------------------------------------------
  # STEP 2: Stage transitions
  # ------------------------------------------------------------
  
  if (any(!is.finite(terminal_survival_prob)) ||
      any(terminal_survival_prob < 0 | terminal_survival_prob > 1)) {
    stop("Terminal-stage survival probabilities must be between 0 and 1")
  }
  
  n[S, ] <- rbinom(
    n = n_pops,
    size = n[S, ],
    prob = terminal_survival_prob
  )
  
  transition_footprint <- if (ApplyFootprintToTransitions) {
    pmin(1, pmax(0, nodepropaguleestablishment))
  } else {
    1
  }
  
  # Weighted population in stages already processed above the target stage.
  capacity_above <- numeric(n_pops)
  
  # The descending order prevents multiple stage transitions per timestep.
  for (s in S:2) {
    N_prev <- n[s - 1, ]
    
    trans_prob <- if (transition_is_list) {
      transition_probabilities[s - 1, ]
    } else {
      transition_probabilities[s - 1]
    }
    
    stay_prob <- if (transition_is_list) {
      stasis_probabilities[s - 1, ]
    } else {
      stasis_probabilities[s - 1]
    }
    
    surv_total_prob <- pmin(trans_prob + stay_prob, 1)
    
    N_surv <- rbinom(
      n = n_pops,
      size = N_prev,
      prob = surv_total_prob
    )
    
    prog_cond_prob <- ifelse(
      surv_total_prob > 0,
      trans_prob / surv_total_prob,
      0
    )
    
    n_trans_candidates <- rbinom(
      n = n_pops,
      size = N_surv,
      prob = prog_cond_prob
    )
    
    # STEP 2c: target stage and all higher stages consume capacity.
    # stage_weight is node-specific when weights was supplied as a matrix.
    stage_weight <- w[, s]
    total_pop <- capacity_above + n[s, ] * stage_weight
    slots_available <- pmax(
      0,
      floor((nodeK - total_pop) / stage_weight)
    )
    
    max_slots <- nodeK / stage_weight
    free_fraction <- ifelse(
      max_slots > 0,
      slots_available / max_slots,
      0
    )
    
    candidate_prob <- pmin(
      1,
      pmax(0, transition_footprint * free_fraction)
    )
    
    n_trans_actual <- rbinom(
      n = n_pops,
      size = n_trans_candidates,
      prob = candidate_prob
    )
    
    n_trans_actual <- pmin(n_trans_actual, slots_available)
    
    # STEP 2e: resolve stasis and mortality of blocked candidates.
    # Individuals allocated to stasis are unaffected. Only candidates that
    # attempted progression but did not obtain a target-stage slot are exposed
    # to BlockedTransitionMortality for their source stage (s - 1).
    n_stay <- N_surv - n_trans_candidates
    n_blocked <- n_trans_candidates - n_trans_actual
    blocked_mortality_prob <- blocked_transition_mortality[, s - 1L]
    n_blocked_survivors <- n_blocked
    
    # Avoid unnecessary random draws at the exact 0 and 1 limits. In
    # particular, BlockedTransitionMortality = 0 preserves the previous RNG
    # sequence as well as the previous biological behaviour.
    certain_death <- n_blocked > 0 & blocked_mortality_prob == 1
    n_blocked_survivors[certain_death] <- 0
    
    stochastic_blocking <-
      n_blocked > 0 &
      blocked_mortality_prob > 0 &
      blocked_mortality_prob < 1
    
    if (any(stochastic_blocking)) {
      n_blocked_survivors[stochastic_blocking] <- rbinom(
        n = sum(stochastic_blocking),
        size = n_blocked[stochastic_blocking],
        prob = 1 - blocked_mortality_prob[stochastic_blocking]
      )
    }
    
    n[s, ] <- n[s, ] + n_trans_actual
    n[s - 1, ] <- n_stay + n_blocked_survivors
    
    # This now equals the weighted population in stages s:S.
    capacity_above <- total_pop + n_trans_actual * stage_weight
  }
  
  # ------------------------------------------------------------
  # STEP 3: Propagule dispersal
  # ------------------------------------------------------------
  
  Pin <- numeric(n_pops)
  Qin <- numeric(n_pops)
  
  if (any(propagules > 0)) {
    # Self-mediated spread.
    sdd_sources <- propagules * (1 - lddrate)
    sdd_destination_weights <- as.numeric(sdd_sources %*% sddprob)
    total_p <- sum(sdd_destination_weights)
    
    if (total_p > 0) {
      Pin <- if (floor(total_p) < MaxInteger) {
        as.numeric(
          rmultinom(
            n = 1,
            size = floor(total_p),
            prob = sdd_destination_weights
          )
        )
      } else {
        # Algebraically equivalent to colSums(sweep(...)), without
        # allocating another n_pops x n_pops matrix.
        as.numeric(floor(sdd_sources) %*% sddprob)
      }
    }
    
    # Human-mediated spread. Management is applied at each source before
    # calculating both the destination weights and the multinomial size.
    if (is.matrix(lddprob)) {
      spread_reduction <- pmin(
        1,
        pmax(
          0,
          rep_len(nodespreadreduction, n_pops) *
            rep_len(managing, n_pops)
        )
      )
      
      ldd_sources <- propagules * lddrate * (1 - spread_reduction)
      ldd_destination_weights <- as.numeric(ldd_sources %*% lddprob)
      total_q <- sum(ldd_destination_weights)
      
      if (total_q > 0) {
        Qin <- if (floor(total_q) < MaxInteger) {
          as.numeric(
            rmultinom(
              n = 1,
              size = floor(total_q),
              prob = ldd_destination_weights
            )
          )
        } else {
          as.numeric(floor(ldd_sources) %*% lddprob)
        }
      }
    }
  }
  
  # ------------------------------------------------------------
  # STEP 4: Recruitment into seedbank
  # ------------------------------------------------------------
  
  seedbank_slots <- pmax(0, floor(node.seedbankK - n[1, ]))
  Pin <- pmax(0, floor(Pin))
  Qin <- pmax(0, floor(Qin))
  env_prob <- pmin(1, pmax(0, nodeenvestabprob))
  
  # A value below one retains the passive-propagule maternal-footprint model.
  # Setting every value to one removes mother dependence: any naturally
  # dispersed propagule that reaches a cell can search all its free slots.
  unrestricted_natural_search <- all(
    is.finite(nodepropaguleestablishment) &
      nodepropaguleestablishment >= 1
  )
  
  if (unrestricted_natural_search) {
    accessible_fraction <- as.numeric(Pin > 0)
    coverage_pressure <- accessible_fraction
  } else {
    footprint_area <- nodepropaguleestablishment * nodeK
    coverage_pressure <- numeric(n_pops)
    
    if (any(Pin > 0) && any(mother_counts > 0)) {
      coverage_numerator <- as.numeric(
        (mother_counts * footprint_area) %*% sddprob
      )
      
      positive_capacity <- nodeK > 0
      coverage_pressure[positive_capacity] <-
        coverage_numerator[positive_capacity] /
        nodeK[positive_capacity]
    }
    
    accessible_fraction <- pmin(
      1,
      pmax(0, -expm1(-coverage_pressure))
    )
  }
  
  natural_slots <- floor(seedbank_slots * accessible_fraction)
  
  # Retain rare, non-zero colonisation links without creating slots when
  # no naturally dispersed propagule arrived.
  natural_slots <- ifelse(
    Pin > 0 & accessible_fraction > 0,
    pmax(1, natural_slots),
    0
  )
  
  natural_slots <- pmin(natural_slots, seedbank_slots)
  other_slots <- pmax(0, seedbank_slots - natural_slots)
  
  lambda_P <- ifelse(
    natural_slots > 0,
    env_prob * Pin / pmax(natural_slots, 1),
    0
  )
  
  lambda_Q <- ifelse(
    seedbank_slots > 0,
    env_prob * Qin / pmax(seedbank_slots, 1),
    0
  )
  
  # Qin contributes inside and outside the natural footprint; Pin only inside.
  p_natural_slots <- -expm1(-(lambda_P + lambda_Q))
  p_other_slots <- -expm1(-lambda_Q)
  
  recruits <-
    rbinom(n_pops, natural_slots, p_natural_slots) +
    rbinom(n_pops, other_slots, p_other_slots)
  
  recruits <- pmin(recruits, Pin + Qin)
  n[1, ] <- n[1, ] + recruits
  
  t(n)
}



# local.dynamics.transition.matrix.skip <- function(
#     nodetransition = NodeTransition,
#     weights = Weights,
#     sddprob = SDDprob,
#     nodeenvestabprob = NodeEnvEstabProb,
#     n0 = N0,
#     lddprob = LDDprob,
#     lddrate = LDDrate,
#     k_is_0 = K_is_0,
#     nodeK = NodeK,
#     node.seedbankK = NodeSeedbankK,
#     nodepropaguleestablishment = NodePropaguleEstablishment,
#     nodespreadreduction = NodeSpreadReduction,
#     managing = Managing,
#     MaxInteger = MaxInteger
# ) {
#   
#   n_pops <- nrow(n0)
#   S <- ncol(n0)
#   w <- if (is.null(weights)) rep(1, S) else weights
#   
#   # --- Create new population matrix ---
#   n <- t(n0)
#   
#   # --- STEP 1: Propagule production (skip seedbank stasis A[1,1]) ---
#   if (is.list(nodetransition)) {
#     fec_means <- vapply(seq_len(n_pops), function(p) {
#       sum(nodetransition[[p]][1, -1] * n[-1, p])
#     }, numeric(1))
#   } else {
#     fec_means <- as.numeric(nodetransition[1, -1] %*% n[-1, , drop = FALSE])
#   }
#   propagules <- ifelse(fec_means > 0, rpois(n_pops, fec_means), 0) 
#   
#   # Pre-extract transition matrices
#   nodetransition_list <- if (is.list(nodetransition)) nodetransition else replicate(n_pops, nodetransition, simplify = FALSE)
#   
#   # --- STEP 0: Terminal stage survival ---
#   surv_prob_S <- vapply(nodetransition_list, function(Ap) Ap[S, S], numeric(1))
#   n[S, ] <- rbinom(n_pops, size = n[S, ], prob = surv_prob_S)
#   
#   # --- STEP 1: Stage transitions with skipping ---
#   total_biomass <- colSums(n * w)
#   
#   for (j in (S-1):1) {
#     N_prev <- n[j, ]
#     trans_mat <- vapply(nodetransition_list, function(Ap) Ap[j:S, j], numeric(S - j + 1))
#     surv_prob <- pmin(colSums(trans_mat), 1)
#     N_surv <- rbinom(n_pops, size = N_prev, prob = surv_prob)
#     
#     dest_prob <- trans_mat
#     pos <- surv_prob > 0
#     dest_prob[, pos] <- sweep(trans_mat[, pos, drop = FALSE], 2, surv_prob[pos], "/")
#     dest_prob[, !pos] <- 0
#     
#     
#     # Multinomial allocation
#     moves <- matrix(0, nrow = S - j + 1, ncol = n_pops)
#     nonzero <- which(N_surv > 0)
#     if (length(nonzero) > 0) {
#       for (p in nonzero) moves[, p] <- rmultinom(1, size = N_surv[p], prob = dest_prob[, p])
#     }
#     
#     # Capacity-aware constraint (sequential stages)
#     free_capacity <- pmax(0, nodeK - total_biomass)
#     stages <- S - j + 1
#     slots <- matrix(0, nrow = stages, ncol = n_pops)
#     for (k in seq_len(stages)) {
#       slots[k, ] <- floor(free_capacity / w[j + k - 1])
#     }
#     
#     # Limit movers by slots
#     moves <- pmin(moves, slots)
#     
#     # Update populations
#     n[j, ] <- N_prev - colSums(moves[-1, , drop = FALSE])
#     for (k in 2:stages) n[j + k - 1, ] <- n[j + k - 1, ] + moves[k, ]
#     
#     # Update biomass ledger
#     delta_biomass <- colSums(moves[-1, , drop = FALSE] * w[(j + 1):(j + nrow(moves)-1)])
#     total_biomass <- total_biomass + delta_biomass
#   }
#   
#   # --- STEP 2: Propagule dispersal ---
#   Pin <- numeric(n_pops)
#   Qin <- numeric(n_pops)
#   if (sum(propagules) > 0) {
#     # Self-mediated dispersal (SDD)
#     total_p <- sum((propagules * (1 - lddrate)) * rowSums(sddprob))
#     sddprob_matrix <- (propagules * (1 - lddrate)) %*% sddprob
#     if (floor(total_p) < MaxInteger) {
#       Pin <- as.numeric(t(rmultinom(1, size = floor(total_p), prob = sddprob_matrix)))
#     } else {
#       Pin <- colSums(sweep(sddprob, 1, floor(propagules * (1 - lddrate)), `*`))
#     }
#     
#     # Human-mediated dispersal (LDD)
#     if (is.matrix(lddprob)) {
#       total_q <- sum((propagules * lddrate) * rowSums(lddprob))
#       lddprob_matrix <- (propagules * lddrate) %*% lddprob
#       if (floor(total_q) < MaxInteger) {
#         Qin <- as.numeric(t(rmultinom(1, size = floor(total_q * (1 - nodespreadreduction * managing)), prob = lddprob_matrix)))
#       } else {
#         Qin <- colSums(sweep(lddprob, 1, floor(propagules * lddrate * (1 - nodespreadreduction * managing)), `*`))
#       }
#     }
#   }
#   
#   # --- STEP 3: Recruitment after dispersal ---
#   slots <- pmax(0, floor((nodeK - colSums(n * w)) / w[1]))
#   max_recruits <- pmin(slots, Pin + Qin)
#   est_prob <- 1 - exp(-nodepropaguleestablishment * nodeenvestabprob * (Pin + Qin))
#   recruits <- rbinom(n_pops, size = max_recruits, prob = est_prob)
#   n[1, ] <- n[1, ] + recruits
#   
#   return(t(n))
# }

INApestMetaTransitionMatrixParallel = function(
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
InfoRetentionProb = 1,       #Probability that existing information is retained between timesteps. Can be single number, vector (nodes) or matrix (nodes x timesteps)
InfoPersistenceSteps = NA,    #Number of timesteps information persists after last known local presence. Can be single number, vector (nodes) or matrix (nodes x timesteps); NA uses InfoRetentionProb
EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
K,		       #Population carrying capacity - vector (nodes)
SeedbankK,        # Seedbank carrying capacity - vector (nodes) or matrix (nodes x timesteps)
PropaguleEstablishment, #Propagules establishment probability
IncursionStartPop=NA,      #option to set population size for new incursions
SDDprob,                   #Natural dispersal probability matrix, or 3D array (nodes x nodes x timesteps)
SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
LDDprob = NA,         #Option to provide long distance (human-mediated) dispersal matrix or 3D array (nodes x nodes x timesteps) instead of distance-independent dispesal rate
			      #e.g. could be weighted by law of human visitation or data on stock movements
LDDrate = 0,         #Proportion of available propagules entering LDD
DispersalDensityFactor = 0,
BlockedTransitionMortality = 0,
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

###Allow SDD and LDD connectivity to vary through time
if(length(dim(SDDprob)) == 3 && (dim(SDDprob)[1] != dim(SDDprob)[2] || dim(SDDprob)[3] != Ntimesteps))
  stop("SDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")
if(length(dim(LDDprob)) == 3 && (dim(LDDprob)[1] != nrow(SDDprob) || dim(LDDprob)[2] != nrow(SDDprob) || dim(LDDprob)[3] != Ntimesteps))
  stop("LDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")

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

###If carrying capacity provided as matrix assign values from first timestep for population initialisation
if(is.matrix(K) == TRUE)
  {
  K_is_0 <- K[,1]<=0
  inv_K <- 1 / sum(K[,1])
  NodeK = K[,1] 
}    
  
  # --- Seedbank carrying capacity ---
if (is.matrix(SeedbankK) == FALSE)
  {
    NodeSeedbankK <- SeedbankK
  }
  
if (is.matrix(SeedbankK) == TRUE)
  {
    NodeSeedbankK <- SeedbankK[, 1]
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
    
## ---- Parallel permutation loop (replace the serial for-loop) ----
library(parallel)

# number of worker cores to use (adjust if needed)
n_cores <- max(1, detectCores() - 1)
cl <- makeCluster(n_cores)
# ensure cluster is stopped on function exit / error
#on.exit(stopCluster(cl), add = TRUE)

## Export variables and functions needed by workers.
## Add/remove names here to match what your worker block references.
vars_to_export <- c(
  "Ntimesteps", "Nstages", "InitialPopulation", "InitBioP", "InvasionRisk",
  "IncursionStartPop", "SDDprob", "EnvEstabProb", "Transition", "LocalDynamics",
  "K", "SeedbankK","PropaguleEstablishment", "SEAM", "ExternalInfoProb", "InfoRetentionProb", "InfoPersistenceSteps", "UseInfoPersistence", "OngoingExternalInvasion",
  "OngoingExternalInfo", "DetectionProb", "DetectionSD", "ManageProb", "ManageSD",
  "MortalityProb", "MortalitySD", "SpreadReduction", "SpreadReductionSD",
  "Weights", "NodeK", "NodeSeedbankK","LDDprob", "LDDrate", "NodePropaguleEstablishment", 
  "DispersalDensityFactor","BlockedTransitionMortality"
  )

clusterExport(cl, vars_to_export, envir = environment())

## If LocalDynamics is a function in the calling env, also export it (already in list).
## If your function depends on packages, load them on each worker:
clusterEvalQ(cl, {
  # load packages used inside worker, e.g. none here but if you use e.g. 'Matrix' add library(Matrix)
  NULL
})

# Run permutations in parallel
perm_results <- parLapply(cl, seq_len(Nperm), function(i_perm) {
  # worker environment
  set.seed(12345 + i_perm)  # reproducible per-perm RNG
  
  ###Max integer for propagule dispersal using rmultinom
  MaxInteger <- .Machine$integer.max  
  # Local containers for this permutation
  n_nodes <- nrow(SDDprob)
  PopulationResults_local      <- matrix(0, nrow = n_nodes, ncol = Ntimesteps)
  PopulationStageResults_local <- array(0, dim = c(n_nodes, Nstages, Ntimesteps))
  InvasionResults_local        <- matrix(0, nrow = n_nodes, ncol = Ntimesteps)
  DetectedResults_local        <- matrix(0, nrow = n_nodes, ncol = Ntimesteps)
  ManagingResults_local        <- matrix(0, nrow = n_nodes, ncol = Ntimesteps)
  
  ## ---------- BEGIN exactly your block (minimal edits) ----------
  # Assign initial infestations
  n_nodes <- nrow(SDDprob)
  InitBio <- matrix(0, n_nodes, Nstages)
  
  if (is.matrix(InitialPopulation) &&
      nrow(InitialPopulation) == n_nodes &&
      ncol(InitialPopulation) == Nstages) {
    
    InitBio <- InitialPopulation
    
  } else {
    
    risk <- if (is.matrix(InvasionRisk) && nrow(InvasionRisk) == n_nodes)
      InvasionRisk[,1] else if (length(InvasionRisk) == n_nodes) InvasionRisk else NULL
    
    if (!is.na(InitBioP)) {
      Infested <- sample.int(n_nodes, ceiling(n_nodes * InitBioP), prob = risk)
    } else if (!is.null(risk)) {
      Infested <- which(rbinom(n_nodes, 1, risk) == 1)
    } else {
      Infested <- integer(0)
    }
    
    InitBio[Infested,1] <- if (is.na(IncursionStartPop)) 1 else IncursionStartPop
  }
  
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
  
  InitBio <- floor(InitBio)   # integers
  N <- InitBio
  if (sum(N) == 0 && OngoingExternalInvasion == FALSE) {
    warning("No initial populations and no future external invasions")
  }
  
  # Robust InitInfo setup
  n_nodes <- nrow(SDDprob)
  InitInfo <- rep(0, n_nodes)
  if (!all(is.na(InitialInfo)) || !is.na(InitInfoP) || !all(is.na(ExternalInfoProb))) {
    if (length(InitialInfo) == n_nodes) {
      InitInfo <- InitialInfo
    } else {
      if (is.na(InitInfoP)) InitInfoP <- 0
      if (is.na(sum(ExternalInfoProb))) ExternalInfoProb <- rep(1, n_nodes)
      if (length(ExternalInfoProb) != n_nodes) ExternalInfoProb <- rep(1, n_nodes)
      n_sample <- ceiling(n_nodes * InitInfoP)
      if (n_sample > 0) {
        Info <- sample(1:n_nodes, size = n_sample, prob = ExternalInfoProb)
        InitInfo[Info] <- 1
      }
    }
  }
  
  # --- Determine dimensions and pre-sample detection/mortality arrays (as in your code) ---
  # Build NodeDetectionProb array (n_nodes x Nstages x Ntimesteps)
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
  } else if (is.matrix(DetectionProb) && all(dim(DetectionProb) == c(n_nodes, Nstages))) {
    NodeDetectionProb <- array(NA, dim = c(n_nodes, Nstages, Ntimesteps))
    for (t in 1:Ntimesteps) {
      NodeDetectionProb[, , t] <- pmax(0, pmin(1,
                                               matrix(
                                                 rnorm(n_nodes * Nstages, mean = as.vector(DetectionProb), sd = as.vector(DetectionSD)),
                                                 nrow = n_nodes, ncol = Nstages
                                               )
      ))
    }
  } else if (length(dim(DetectionProb)) == 3 && all(dim(DetectionProb)[1:2] == c(n_nodes, Nstages))) {
    NodeDetectionProb <- array(NA, dim = dim(DetectionProb))
    for (t in 1:Ntimesteps) {
      NodeDetectionProb[, , t] <- pmax(0, pmin(1,
                                               matrix(
                                                 rnorm(n_nodes * Nstages, mean = as.vector(DetectionProb[, , t]), sd = as.vector(DetectionSD)),
                                                 nrow = n_nodes, ncol = Nstages
                                               )
      ))
    }
  } else {
    stop("Unsupported DetectionProb shape in worker")
  }
  
  # NodeMortalityProb: analogous logic (kept concise - follow same pattern)
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
  } else if (is.matrix(MortalityProb) && all(dim(MortalityProb) == c(n_nodes, Nstages))) {
    NodeMortalityProb <- array(NA, dim = c(n_nodes, Nstages, Ntimesteps))
    for (t in 1:Ntimesteps) {
      NodeMortalityProb[, , t] <- pmax(0, pmin(1,
                                               matrix(
                                                 rnorm(n_nodes * Nstages, mean = as.vector(MortalityProb), sd = as.vector(MortalitySD)),
                                                 nrow = n_nodes, ncol = Nstages
                                               )
      ))
    }
  } else if (length(dim(MortalityProb)) == 3 && all(dim(MortalityProb)[1:2] == c(n_nodes, Nstages))) {
    NodeMortalityProb <- array(NA, dim = dim(MortalityProb))
    for (t in 1:Ntimesteps) {
      NodeMortalityProb[, , t] <- pmax(0, pmin(1,
                                               matrix(
                                                 rnorm(n_nodes * Nstages, mean = as.vector(MortalityProb[, , t]), sd = as.vector(MortalitySD)),
                                                 nrow = n_nodes, ncol = Nstages
                                               )
      ))
    }
  } else {
    stop("Unsupported MortalityProb shape in worker")
  }
  
  # Populate initial invaded and detection
  Invaded <- ifelse(rowSums(InitBio) > 0, 1, 0)
  
  # Calculate probability of detection per node at timestep 1 and draw InitDetection
  prob_detect <- 1 - apply((1 - NodeDetectionProb[, , 1])^InitBio, 1, prod)
  InitDetection <- rbinom(n = nrow(InitBio), size = 1, prob = prob_detect)
  InitInfo[InitInfo == 0] <- InitDetection[InitInfo == 0]
  HaveInfo <- InitInfo
  
  # Track the most recent timestep with known local presence
  LastKnownPresence <- rep(NA, nrow(SDDprob))
  if (UseInfoPersistence == TRUE) {
    InitialKnownPresence <- which(InitDetection == 1)
    if (length(InitialKnownPresence) > 0) LastKnownPresence[InitialKnownPresence] <- 0
  }
  
  # ---------- Now run timesteps (use i_perm in prints) ----------
  for (timestep in 1:Ntimesteps) {
    # minimal progress print (safe inside worker)
    cat(sprintf("Worker %d: timestep %d\n", i_perm, timestep))
    
    # update dispersal connectivity if supplied as nodes x nodes x timesteps arrays
    NodeSDDprob <- SDDprob
    if (length(dim(SDDprob)) == 3) NodeSDDprob <- SDDprob[, , timestep]
    NodeLDDprob <- LDDprob
    if (length(dim(LDDprob)) == 3) NodeLDDprob <- LDDprob[, , timestep]

    # update NodeEnvEstabProb, NodeTransition, NodeK, NodePropaguleEstablishment if matrices/arrays
    if (is.matrix(EnvEstabProb) == TRUE) NodeEnvEstabProb <- EnvEstabProb[, timestep]
    if (is.array(Transition) && length(dim(Transition)) == 3) {
      NodeTransition <- Transition[, , timestep]
    } else if (is.list(Transition)) {
      NodeTransition <- Transition
    } else {
      NodeTransition <- Transition
    }
    if (is.matrix(K) == TRUE) NodeK <- K[, timestep]
    
    ###If seedbank carrying capacity provided as matrix assign values for relevant timestep
    if (is.matrix(SeedbankK) == TRUE) NodeSeedbankK <- SeedbankK[, timestep]
    
    
    # ManageProb / SpreadReduction per timestep if provided as matrices
    if (is.matrix(ManageProb) == TRUE && nrow(ManageProb) == nrow(SDDprob) && ncol(ManageProb) == Ntimesteps) {
      NodeManageProb <- pmin(1, pmax(0, rnorm(n = nrow(SDDprob), mean = ManageProb[, timestep], sd = ManageSD)))
    } else {
      NodeManageProb <- pmin(1, pmax(0, rnorm(n = nrow(SDDprob), mean = ManageProb, sd = ManageSD)))
    }
    if (is.matrix(SpreadReduction) == TRUE && nrow(SpreadReduction) == nrow(SDDprob) && ncol(SpreadReduction) == Ntimesteps) {
      NodeSpreadReduction <- pmin(1, pmax(0, rnorm(n = nrow(SDDprob), mean = SpreadReduction[, timestep], sd = SpreadReductionSD)))
    } else {
      NodeSpreadReduction <- pmin(1, pmax(0, rnorm(n = nrow(SDDprob), mean = SpreadReduction, sd = SpreadReductionSD)))
    }
    
    # Assign management (only where HaveInfo)
    Managing <- rbinom(n = nrow(SDDprob), size = 1, prob = NodeManageProb * HaveInfo)
    Managing <- Managing * HaveInfo
    Detected <- Invaded * HaveInfo
    
    # Apply management mortality by stage (use NodeMortalityProb[,,timestep])
    N0 <- matrix(NA, nrow = nrow(N), ncol = ncol(N))
    for (s in 1:Nstages) {
      N0[, s] <- rbinom(n = nrow(N), size = N[, s], prob = 1 - (NodeMortalityProb[, s, timestep] * Managing))
    }
    
    # Track known local presence from actual management mortality
    if (UseInfoPersistence == TRUE) {
      KnownPresence <- which(rowSums(N - N0) > 0)
      if (length(KnownPresence) > 0) LastKnownPresence[KnownPresence] <- timestep
    }
    N <- N0
    if (sum(N0) > 0) {
      N <- LocalDynamics(nodetransition = NodeTransition, weights = Weights, sddprob = NodeSDDprob,
                         nodeenvestabprob = NodeEnvEstabProb, n0 = N0, lddprob = NodeLDDprob,
                         lddrate = LDDrate,  nodeK = NodeK, node.seedbankK = NodeSeedbankK,
                         nodepropaguleestablishment = NodePropaguleEstablishment,
                         nodespreadreduction = NodeSpreadReduction, managing = Managing,
                         MaxInteger=MaxInteger,BlockedTransitionMortality = BlockedTransitionMortality,
                         DispersalDensityFactor=DispersalDensityFactor)
    }
    
    # Apply programmed stopping after last known local presence
    NodeInfoPersistenceSteps <- InfoPersistenceSteps
    if (is.matrix(InfoPersistenceSteps) == TRUE) NodeInfoPersistenceSteps <- InfoPersistenceSteps[, timestep]
    if (length(NodeInfoPersistenceSteps) == 1) NodeInfoPersistenceSteps <- rep(NodeInfoPersistenceSteps, nrow(SDDprob))
    ProgrammedInfoNodes <- which(HaveInfo == 1 & !is.na(NodeInfoPersistenceSteps))
    if (length(ProgrammedInfoNodes) > 0) {
      TimeSinceKnownPresence <- timestep - LastKnownPresence
      InfoStopNodes <- ProgrammedInfoNodes[is.na(LastKnownPresence[ProgrammedInfoNodes]) | TimeSinceKnownPresence[ProgrammedInfoNodes] >= NodeInfoPersistenceSteps[ProgrammedInfoNodes]]
      if (length(InfoStopNodes) > 0) HaveInfo[InfoStopNodes] <- 0
    }
    
    # Allow information to decay after management and spread where no programmed stop is supplied
    NodeInfoRetentionProb <- InfoRetentionProb
    if (is.matrix(InfoRetentionProb) == TRUE) NodeInfoRetentionProb <- InfoRetentionProb[, timestep]
    if (length(NodeInfoRetentionProb) == 1) NodeInfoRetentionProb <- rep(NodeInfoRetentionProb, nrow(SDDprob))
    InfoDecayNodes <- which(HaveInfo == 1 & is.na(NodeInfoPersistenceSteps) & NodeInfoRetentionProb < 1)
    if (length(InfoDecayNodes) > 0) {
      HaveInfo[InfoDecayNodes] <- rbinom(n = length(InfoDecayNodes), size = 1, prob = NodeInfoRetentionProb[InfoDecayNodes])
    }
    
    # Info spread via SEAM
    if (is.matrix(SEAM) == TRUE) {
      RandSEAM[] <- rbinom(n = nrow(SDDprob)^2, size = 1, prob = SEAM * Detected)
      InfoTransferred <- ifelse(colSums(RandSEAM) > 0, 1, 0)
      HaveInfo[HaveInfo == 0] <- InfoTransferred[HaveInfo == 0]
    }
    
    # External invasion
    if (OngoingExternalInvasion == TRUE) {
      if (is.matrix(InvasionRisk) == FALSE) {
        ExternalInvasion <- rbinom(1:nrow(SDDprob), size = 1, prob = InvasionRisk)
      } else {
        ExternalInvasion <- rbinom(1:nrow(SDDprob), size = 1, prob = InvasionRisk[, timestep])
      }
      Invaded[Invaded == 0] <- ExternalInvasion[Invaded == 0]
      if (is.na(IncursionStartPop) == TRUE) N[, 1] <- N[, 1] + ExternalInvasion
      else N[, 1] <- N[, 1] + ExternalInvasion * IncursionStartPop
      
      
    }
    N <- floor(N)
    # External info
    if (OngoingExternalInfo == TRUE) {
      if (is.matrix(ExternalInfoProb) == FALSE) ExternalInfo <- rbinom(1:nrow(SDDprob), size = 1, prob = ExternalInfoProb)
      else ExternalInfo <- rbinom(1:nrow(SDDprob), size = 1, prob = ExternalInfoProb[, timestep])
      HaveInfo[HaveInfo == 0] <- ExternalInfo[HaveInfo == 0]
    }
    
    # Update invaded vector
    Invaded <- ifelse(rowSums(N) > 0, 1, 0)
    
    # Record results for this timestep into local containers
    ManagingResults_local[, timestep] <- Managing
    InvasionResults_local[, timestep] <- Invaded
    weighted_population <- if (is.matrix(Weights)) {
      rowSums(N[, 2:Nstages, drop = FALSE] *
                Weights[, 2:Nstages, drop = FALSE])
    } else {
      as.numeric(N[, 2:Nstages, drop = FALSE] %*% Weights[2:Nstages])
    }
    PopulationResults_local[, timestep] <- weighted_population
    PopulationStageResults_local[, , timestep] <- N
    
    # detection: compute per-stage detection probs at current N and NodeDetectionProb[,,timestep]
    DetectionProbPerStage <- 1 - (1 - NodeDetectionProb[, , timestep])^N  # nodes x stages
    ProbDetectNode <- 1 - apply(1 - DetectionProbPerStage, 1, prod)      # nodes
    NewHaveInfo <- rbinom(n = nrow(SDDprob), size = 1, prob = ProbDetectNode)
    # Record newly detected infestations as known local presence
    if (UseInfoPersistence == TRUE) {
      KnownPresence <- which(NewHaveInfo == 1)
      if (length(KnownPresence) > 0) LastKnownPresence[KnownPresence] <- timestep
    }
    HaveInfo[HaveInfo == 0] <- NewHaveInfo[HaveInfo == 0]
    DetectedResults_local[, timestep] <- HaveInfo * Invaded
  } # end timesteps
  
  ## ---------- END your block ----------
  
  # Return local results for this permutation
  list(
    PopulationResults = PopulationResults_local,
    PopulationStageResults = PopulationStageResults_local,
    InvasionResults = InvasionResults_local,
    DetectedResults = DetectedResults_local,
    ManagingResults = ManagingResults_local
  )
}) # end parLapply

# combine the returned per-perm results into master arrays
for (i in seq_along(perm_results)) {
  PopulationResults[,, i]       <- perm_results[[i]]$PopulationResults
  PopulationStageResults[,,, i] <- perm_results[[i]]$PopulationStageResults
  InvasionResults[,, i]         <- perm_results[[i]]$InvasionResults
  DetectedResults[,, i]         <- perm_results[[i]]$DetectedResults
  ManagingResults[,, i]         <- perm_results[[i]]$ManagingResults
}

# stopCluster will be called by on.exit when function ends or errors
# --- Shut down parallel cluster ---
if (exists("cl")) {
  parallel::stopCluster(cl)
  rm(cl)
  cat("Parallel cluster stopped and removed from memory.\n")
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
