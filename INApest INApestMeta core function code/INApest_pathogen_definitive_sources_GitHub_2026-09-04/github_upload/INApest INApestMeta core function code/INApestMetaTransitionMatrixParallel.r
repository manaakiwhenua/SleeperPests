###############################################################################
### INApestMetaTransitionMatrixParallel -- parallel node x life-stage engine
###
### Parallel implementation of the Transition Matrix architecture. Independent
### permutations use the same stage-structured host dynamics, transition movement,
### information, management and optional pathogen processes as the serial engine.
###
### Parallelisation changes execution only; biological model semantics are the
### same as INApestMetaTransitionMatrix.
###############################################################################

# Default node x stage survival, progression, reproduction and dispersal process.
local.dynamics.transition.matrix <- function(
    nodetransition = NodeTransition,
    weights = Weights,
    sddprob = SDDprob,
    nodeenvestabprob = NodeEnvEstabProb,
    n0 = N0,
    lddprob = LDDprob,
    lddrate = LDDrate,
    transition_sddprob = NULL,
    transition_lddprob = NULL,
    transition_lddrate = 0,
    nodeK = NodeK,
    node.seedbankK = NodeSeedbankK,
    nodepropaguleestablishment = NodePropaguleEstablishment,
    nodespreadreduction = NodeSpreadReduction,
    nodefecundityreduction = NodeFecundityReduction,
    managing = Managing,
    MaxInteger = MaxInteger,
    ApplyFootprintToTransitions = FALSE,
    BlockedTransitionMortality = 0,
    DispersalDensityFactor = 0
) {
  n_pops <- nrow(n0)
  S <- ncol(n0)
  n <- t(n0)

  if(is.null(dim(nodefecundityreduction))) {
    if(length(nodefecundityreduction) > 1 && n_pops == S && length(nodefecundityreduction) == n_pops)
      stop("nodefecundityreduction vector is ambiguous because n_pops equals S; supply an n_pops x S matrix")
    if(length(nodefecundityreduction) == 1) {
      nodefecundityreduction <- matrix(nodefecundityreduction, nrow = n_pops, ncol = S)
    } else if(length(nodefecundityreduction) == n_pops) {
      nodefecundityreduction <- matrix(rep(nodefecundityreduction, S), nrow = n_pops, ncol = S)
    } else if(length(nodefecundityreduction) == S) {
      nodefecundityreduction <- matrix(rep(nodefecundityreduction, each = n_pops), nrow = n_pops, ncol = S)
    } else {
      stop("nodefecundityreduction must be scalar, length n_pops, length S, or an n_pops x S matrix")
    }
  } else if(!identical(dim(nodefecundityreduction), c(n_pops, S))) {
    stop("nodefecundityreduction matrix must have dimensions n_pops x S")
  }
  if(any(!is.finite(nodefecundityreduction)) || any(nodefecundityreduction < 0) || any(nodefecundityreduction > 1))
    stop("nodefecundityreduction values must be between 0 and 1")
  EffectiveFecundityMultiplier <- 1 - nodefecundityreduction * managing
  
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
  # Optional spatial movement during stage progression
  # ------------------------------------------------------------
  # A bare n_pops x n_pops matrix applies to every progression transition.
  # A list of length S - 1 permits transition-specific matrices; NULL list
  # elements retain the historical local-only transition for that source stage.
  # The top-level INApestMetaTransitionMatrix function resolves any time-varying
  # node x node x timestep arrays before calling LocalDynamics.
  .normalise_transition_movement <- function(x, name) {
    if (is.null(x)) {
      return(rep(list(NULL), S - 1L))
    }

    if (is.list(x)) {
      if (length(x) != S - 1L) {
        stop(name, " list must have length S - 1 (one element per progression transition)")
      }
      out <- x
    } else {
      out <- rep(list(x), S - 1L)
    }

    for (k in seq_len(S - 1L)) {
      P <- out[[k]]
      if (is.null(P)) next
      if (!is.matrix(P) || !identical(dim(P), c(n_pops, n_pops))) {
        stop(name, " entries must be NULL or n_pops x n_pops matrices")
      }
      if (any(!is.finite(P)) || any(P < 0)) {
        stop(name, " matrices must contain finite non-negative probabilities")
      }
      rs <- rowSums(P)
      if (any(rs > 1 + 1e-10)) {
        stop(name, " matrix row sums must not exceed 1; unused row mass represents movement outside the modelled landscape")
      }
      # Remove tiny floating-point excess so the residual outside probability
      # used by rmultinom is always non-negative.
      if (any(rs > 1)) {
        P <- P / pmax(1, rs)
        out[[k]] <- P
      }
    }
    out
  }

  transition_sdd <- .normalise_transition_movement(
    transition_sddprob,
    "transition_sddprob"
  )
  transition_ldd <- .normalise_transition_movement(
    transition_lddprob,
    "transition_lddprob"
  )

  if (!(length(transition_lddrate) %in% c(1L, S - 1L))) {
    stop("transition_lddrate must be scalar or length S - 1")
  }
  transition_lddrate <- rep_len(transition_lddrate, S - 1L)
  if (any(!is.finite(transition_lddrate)) ||
      any(transition_lddrate < 0 | transition_lddrate > 1)) {
    stop("transition_lddrate values must be between 0 and 1")
  }

  transition_movement_active <- vapply(
    seq_len(S - 1L),
    function(k) !is.null(transition_sdd[[k]]) || !is.null(transition_ldd[[k]]),
    logical(1)
  )

  # Allocate transitioning individuals from each source node to internal
  # destinations plus an implicit outside-landscape destination. As for the
  # existing transition-model dispersal semantics, missing row mass is export.
  .transition_flow_matrix <- function(source_counts, P) {
    flows <- matrix(0, nrow = n_pops, ncol = n_pops)
    exported <- numeric(n_pops)

    active_sources <- which(source_counts > 0)
    if (length(active_sources) == 0L) {
      return(list(flows = flows, exported = exported))
    }

    for (i in active_sources) {
      n_i <- floor(source_counts[i])
      if (n_i <= 0) next
      p_internal <- pmax(0, P[i, ])
      p_outside <- max(0, 1 - sum(p_internal))
      probs <- c(p_internal, p_outside)
      prob_total <- sum(probs)

      if (prob_total <= 0) next
      probs <- probs / prob_total

      if (n_i < MaxInteger) {
        allocation <- as.numeric(rmultinom(1, size = n_i, prob = probs))
      } else {
        # Preserve integer-valued stage abundances when counts exceed the
        # rmultinom integer limit. Rounding remainder is assigned to export.
        allocation <- floor(n_i * probs)
        remainder <- n_i - sum(allocation)
        allocation[length(allocation)] <- allocation[length(allocation)] + remainder
      }

      flows[i, ] <- allocation[seq_len(n_pops)]
      exported[i] <- allocation[n_pops + 1L]
    }

    list(flows = flows, exported = exported)
  }

  # Uniformly select accepted transitioners from source-specific arrivals at a
  # destination without allowing any source to contribute more than arrived.
  .take_without_replacement <- function(counts, n_take) {
    counts <- pmax(0, floor(counts))
    n_take <- min(floor(n_take), sum(counts))
    out <- numeric(length(counts))
    if (n_take <= 0 || sum(counts) <= 0) return(out)
    if (n_take >= sum(counts)) return(counts)

    if (sum(counts) >= MaxInteger) {
      expected <- n_take * counts / sum(counts)
      out <- floor(expected)
      remainder <- n_take - sum(out)
      if (remainder > 0) {
        spare <- counts - out
        ord <- order(expected - out, decreasing = TRUE)
        for (idx in ord) {
          if (remainder <= 0) break
          add <- min(spare[idx], remainder)
          if (add > 0) {
            out[idx] <- out[idx] + add
            remainder <- remainder - add
          }
        }
      }
      return(out)
    }

    remaining_take <- n_take
    remaining_total <- sum(counts)
    last_positive <- max(which(counts > 0))

    for (i in seq_along(counts)) {
      if (remaining_take <= 0) break
      if (counts[i] <= 0) {
        remaining_total <- remaining_total - counts[i]
        next
      }
      if (i == last_positive) {
        out[i] <- remaining_take
        break
      }
      out[i] <- rhyper(
        1,
        m = counts[i],
        n = remaining_total - counts[i],
        k = remaining_take
      )
      remaining_take <- remaining_take - out[i]
      remaining_total <- remaining_total - counts[i]
    }
    out
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
    
    EffectiveFecundity <- fecundity * t(EffectiveFecundityMultiplier[, 2:S, drop = FALSE])
    fec_means <- colSums(EffectiveFecundity * n[-1, , drop = FALSE])
    mother_counts <- colSums(
      n[-1, , drop = FALSE] * (fecundity > 0)
    )
  } else {
    fecundity <- nodetransition[1, -1]
    transition_probabilities <-
      nodetransition[cbind(2:S, 1:(S - 1))]
    stasis_probabilities <- diag(nodetransition)[1:(S - 1)]
    terminal_survival_prob <- nodetransition[S, S]
    
    EffectiveFecundity <- matrix(rep(fecundity, n_pops), nrow = S - 1L, ncol = n_pops) *
      t(EffectiveFecundityMultiplier[, 2:S, drop = FALSE])
    fec_means <- colSums(EffectiveFecundity * n[-1, , drop = FALSE])
    
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
    
    # A transition with no transition-specific movement matrices follows the
    # historical local-only path exactly, including its random-number sequence.
    if (!transition_movement_active[s - 1L]) {
      n_trans_actual <- rbinom(
        n = n_pops,
        size = n_trans_candidates,
        prob = candidate_prob
      )
      n_trans_actual <- pmin(n_trans_actual, slots_available)
      n_blocked <- n_trans_candidates - n_trans_actual
      transition_additions <- n_trans_actual
    } else {
      P_sdd <- transition_sdd[[s - 1L]]
      P_ldd <- transition_ldd[[s - 1L]]

      # If both movement modes are supplied, transition_lddrate gives the
      # fraction entering LDD. If only one matrix is supplied, every moving
      # transition candidate uses that matrix and the rate is immaterial.
      if (!is.null(P_sdd) && !is.null(P_ldd)) {
        r_ldd <- transition_lddrate[s - 1L]
        if (r_ldd <= 0) {
          ldd_candidates <- numeric(n_pops)
          sdd_candidates <- n_trans_candidates
        } else if (r_ldd >= 1) {
          ldd_candidates <- n_trans_candidates
          sdd_candidates <- numeric(n_pops)
        } else {
          ldd_candidates <- rbinom(
            n = n_pops,
            size = n_trans_candidates,
            prob = r_ldd
          )
          sdd_candidates <- n_trans_candidates - ldd_candidates
        }
      } else if (!is.null(P_sdd)) {
        sdd_candidates <- n_trans_candidates
        ldd_candidates <- numeric(n_pops)
      } else {
        sdd_candidates <- numeric(n_pops)
        ldd_candidates <- n_trans_candidates
      }

      sdd_flow <- if (!is.null(P_sdd)) {
        .transition_flow_matrix(sdd_candidates, P_sdd)
      } else {
        list(flows = matrix(0, n_pops, n_pops), exported = numeric(n_pops))
      }
      ldd_flow <- if (!is.null(P_ldd)) {
        .transition_flow_matrix(ldd_candidates, P_ldd)
      } else {
        list(flows = matrix(0, n_pops, n_pops), exported = numeric(n_pops))
      }

      internal_flows <- sdd_flow$flows + ldd_flow$flows
      exported_by_source <- sdd_flow$exported + ldd_flow$exported
      arrivals <- colSums(internal_flows)

      # Capacity and optional transition footprint are evaluated at the
      # destination of the moving transition. Candidates exported beyond the
      # modelled landscape are successful departures and do not compete for
      # internal slots.
      accepted_by_destination <- rbinom(
        n = n_pops,
        size = arrivals,
        prob = candidate_prob
      )
      accepted_by_destination <- pmin(
        accepted_by_destination,
        slots_available,
        arrivals
      )

      accepted_by_source <- numeric(n_pops)
      destination_nodes <- which(arrivals > 0 & accepted_by_destination > 0)
      for (j in destination_nodes) {
        accepted_sources_j <- .take_without_replacement(
          internal_flows[, j],
          accepted_by_destination[j]
        )
        accepted_by_source <- accepted_by_source + accepted_sources_j
      }

      internal_by_source <- rowSums(internal_flows)
      n_blocked <- pmax(0, internal_by_source - accepted_by_source)
      transition_additions <- accepted_by_destination

      # Sanity check: every progression candidate is either accepted into an
      # internal destination, blocked internally, or exported from the model.
      if (any(abs(
        n_trans_candidates -
          (accepted_by_source + n_blocked + exported_by_source)
      ) > 1e-8)) {
        stop("Internal error while allocating movement during stage transition")
      }
    }
    
    # STEP 2e: resolve stasis and mortality of blocked candidates.
    # Individuals allocated to stasis are unaffected. Only candidates that
    # attempted progression, remained inside the modelled landscape, and did
    # not obtain a target-stage slot are exposed to BlockedTransitionMortality.
    # Blocked survivors return to their source node and source stage. Movement
    # outside the modelled landscape is treated as export, not blocking.
    n_stay <- N_surv - n_trans_candidates
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
    
    n[s, ] <- n[s, ] + transition_additions
    n[s - 1, ] <- n_stay + n_blocked_survivors
    
    # This now equals the weighted population in stages s:S at each
    # destination node. For local transitions transition_additions is the
    # historical n_trans_actual vector; for moving transitions it is the
    # accepted destination total.
    capacity_above <- total_pop + transition_additions * stage_weight
  }
  
  # ------------------------------------------------------------
  # STEP 3: Propagule dispersal
  # ------------------------------------------------------------
  # Route realised integer propagules, rather than converting an expected
  # fractional number of movers into rmultinom(size). This is the same
  # correction validated for INApestMeta: first assign each realised
  # propagule to SDD/LDD, then allocate it from its source row to an internal
  # destination or the residual outside-landscape category.
  Pin <- numeric(n_pops)
  Qin <- numeric(n_pops)

  # Route realised reproductive propagules through a transition-matrix dispersal surface.
  .route_integer_propagules_tm <- function(source_counts, P, label) {
    arrivals <- numeric(n_pops)
    if(!is.matrix(P)) return(arrivals)
    if(!identical(dim(P),c(n_pops,n_pops))) stop(label," must be an n_pops x n_pops matrix")
    active_sources <- which(source_counts > 0)
    if(!length(active_sources)) return(arrivals)
    for(i in active_sources) {
      n_i <- floor(source_counts[i])
      if(n_i <= 0) next
      p_internal <- pmax(0,as.numeric(P[i,]))
      rs <- sum(p_internal)
      if(!is.finite(rs) || rs > 1 + 1e-10) stop(label," source-row probabilities must be finite, non-negative, and sum to at most 1")
      if(rs > 1) p_internal <- p_internal / rs
      probs <- c(p_internal,max(0,1-sum(p_internal)))
      if(sum(probs) <= 0) next
      if(n_i < MaxInteger) {
        allocation <- as.numeric(rmultinom(1,size=n_i,prob=probs))
      } else {
        # Retain the previous large-count fallback without fractional rmultinom sizes.
        allocation <- floor(n_i * probs)
        allocation[length(allocation)] <- allocation[length(allocation)] + (n_i-sum(allocation))
      }
      arrivals <- arrivals + allocation[seq_len(n_pops)]
    }
    arrivals
  }

  if(any(propagules > 0)) {
    # Integer SDD/LDD branching gives the intended Poisson-thinning
    # marginals when propagule production itself is Poisson.
    if(lddrate <= 0) {
      ldd_sources <- numeric(n_pops)
      sdd_sources <- propagules
    } else if(lddrate >= 1) {
      ldd_sources <- propagules
      sdd_sources <- numeric(n_pops)
    } else {
      ldd_sources <- rbinom(n_pops,size=propagules,prob=lddrate)
      sdd_sources <- propagules - ldd_sources
    }

    Pin <- .route_integer_propagules_tm(sdd_sources,sddprob,"sddprob")

    if(is.matrix(lddprob) && any(ldd_sources > 0)) {
      spread_reduction <- pmin(1,pmax(0,rep_len(nodespreadreduction,n_pops) * rep_len(managing,n_pops)))
      keep_prob <- 1-spread_reduction
      keep_sources <- ldd_sources
      stochastic_keep <- ldd_sources > 0 & keep_prob > 0 & keep_prob < 1
      keep_sources[keep_prob <= 0] <- 0
      if(any(stochastic_keep))
        keep_sources[stochastic_keep] <- rbinom(sum(stochastic_keep),size=ldd_sources[stochastic_keep],prob=keep_prob[stochastic_keep])
      Qin <- .route_integer_propagules_tm(keep_sources,lddprob,"lddprob")
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
# nodetransition = NodeTransition,
# weights = Weights,
# sddprob = SDDprob,
# nodeenvestabprob = NodeEnvEstabProb,
# n0 = N0,
# lddprob = LDDprob,
# lddrate = LDDrate,
# k_is_0 = K_is_0,
# nodeK = NodeK,
# node.seedbankK = NodeSeedbankK,
# nodepropaguleestablishment = NodePropaguleEstablishment,
# nodespreadreduction = NodeSpreadReduction,
# managing = Managing,
# MaxInteger = MaxInteger
# ) {
#   
# n_pops <- nrow(n0)
# S <- ncol(n0)
# w <- if (is.null(weights)) rep(1, S) else weights
#   
# # --- Create new population matrix ---
# n <- t(n0)
#   
# # --- STEP 1: Propagule production (skip seedbank stasis A[1,1]) ---
# if (is.list(nodetransition)) {
# fec_means <- vapply(seq_len(n_pops), function(p) {
# sum(nodetransition[[p]][1, -1] * n[-1, p])
# }, numeric(1))
# } else {
# fec_means <- as.numeric(nodetransition[1, -1] %*% n[-1, , drop = FALSE])
# }
# propagules <- ifelse(fec_means > 0, rpois(n_pops, fec_means), 0)
#   
# # Pre-extract transition matrices
# nodetransition_list <- if (is.list(nodetransition)) nodetransition else replicate(n_pops, nodetransition, simplify = FALSE)
#   
# # --- STEP 0: Terminal stage survival ---
# surv_prob_S <- vapply(nodetransition_list, function(Ap) Ap[S, S], numeric(1))
# n[S, ] <- rbinom(n_pops, size = n[S, ], prob = surv_prob_S)
#   
# # --- STEP 1: Stage transitions with skipping ---
# total_biomass <- colSums(n * w)
#   
# for (j in (S-1):1) {
# N_prev <- n[j, ]
# trans_mat <- vapply(nodetransition_list, function(Ap) Ap[j:S, j], numeric(S - j + 1))
# surv_prob <- pmin(colSums(trans_mat), 1)
# N_surv <- rbinom(n_pops, size = N_prev, prob = surv_prob)
#     
# dest_prob <- trans_mat
# pos <- surv_prob > 0
# dest_prob[, pos] <- sweep(trans_mat[, pos, drop = FALSE], 2, surv_prob[pos], "/")
# dest_prob[, !pos] <- 0
#     
#     
# # Multinomial allocation
# moves <- matrix(0, nrow = S - j + 1, ncol = n_pops)
# nonzero <- which(N_surv > 0)
# if (length(nonzero) > 0) {
# for (p in nonzero) moves[, p] <- rmultinom(1, size = N_surv[p], prob = dest_prob[, p])
# }
#     
# # Capacity-aware constraint (sequential stages)
# free_capacity <- pmax(0, nodeK - total_biomass)
# stages <- S - j + 1
# slots <- matrix(0, nrow = stages, ncol = n_pops)
# for (k in seq_len(stages)) {
# slots[k, ] <- floor(free_capacity / w[j + k - 1])
# }
#     
# # Limit movers by slots
# moves <- pmin(moves, slots)
#     
# # Update populations
# n[j, ] <- N_prev - colSums(moves[-1, , drop = FALSE])
# for (k in 2:stages) n[j + k - 1, ] <- n[j + k - 1, ] + moves[k, ]
#     
# # Update biomass ledger
# delta_biomass <- colSums(moves[-1, , drop = FALSE] * w[(j + 1):(j + nrow(moves)-1)])
# total_biomass <- total_biomass + delta_biomass
# }
#   
# # --- STEP 2: Propagule dispersal ---
# Pin <- numeric(n_pops)
# Qin <- numeric(n_pops)
# if (sum(propagules) > 0) {
# # Self-mediated dispersal (SDD)
# total_p <- sum((propagules * (1 - lddrate)) * rowSums(sddprob))
# sddprob_matrix <- (propagules * (1 - lddrate)) %*% sddprob
# if (floor(total_p) < MaxInteger) {
# Pin <- as.numeric(t(rmultinom(1, size = floor(total_p), prob = sddprob_matrix)))
# } else {
# Pin <- colSums(sweep(sddprob, 1, floor(propagules * (1 - lddrate)), `*`))
# }
#     
# # Human-mediated dispersal (LDD)
# if (is.matrix(lddprob)) {
# total_q <- sum((propagules * lddrate) * rowSums(lddprob))
# lddprob_matrix <- (propagules * lddrate) %*% lddprob
# if (floor(total_q) < MaxInteger) {
# Qin <- as.numeric(t(rmultinom(1, size = floor(total_q * (1 - nodespreadreduction * managing)), prob = lddprob_matrix)))
# } else {
# Qin <- colSums(sweep(lddprob, 1, floor(propagules * lddrate * (1 - nodespreadreduction * managing)), `*`))
# }
# }
# }
#   
# # --- STEP 3: Recruitment after dispersal ---
# slots <- pmax(0, floor((nodeK - colSums(n * w)) / w[1]))
# max_recruits <- pmin(slots, Pin + Qin)
# est_prob <- 1 - exp(-nodepropaguleestablishment * nodeenvestabprob * (Pin + Qin))
# recruits <- rbinom(n_pops, size = max_recruits, prob = est_prob)
# n[1, ] <- n[1, ] + recruits
#   
# return(t(n))
# }


# -----------------------------------------------------------------------------
# Optional custom LocalDynamics arguments
# -----------------------------------------------------------------------------
# Ordinary LocalDynamicsArgs entries are passed unchanged. Wrap a vector,
# matrix, array, or list in INApestLocalDynamicsTimeArg() when its final
# dimension/list position is indexed by simulation timestep. This avoids
# guessing whether an arbitrary custom matrix is static or time-varying.
INApestLocalDynamicsTimeArg <- function(x) {
  structure(list(values = x), class = "INApestLocalDynamicsTimeArg")
}

# Resolve custom LocalDynamics arguments for the current timestep.
.resolve_INApest_LocalDynamicsArgs <- function(LocalDynamicsArgs, timestep, Ntimesteps) {
  if (is.null(LocalDynamicsArgs)) LocalDynamicsArgs <- list()
  if (!is.list(LocalDynamicsArgs))
    stop("LocalDynamicsArgs must be a named list")
  if (!length(LocalDynamicsArgs)) return(list())
  if (is.null(names(LocalDynamicsArgs)) || any(!nzchar(names(LocalDynamicsArgs))))
    stop("Every LocalDynamicsArgs entry must have a non-empty name")
  if (anyDuplicated(names(LocalDynamicsArgs)))
    stop("LocalDynamicsArgs names must be unique")

  resolve_one <- function(x, name) {
    if (inherits(x, "INApestLocalDynamicsTimeArg")) {
      values <- x$values
      if (is.list(values)) {
        if (length(values) != Ntimesteps)
          stop(name, " wrapped list must have length Ntimesteps")
        return(values[[timestep]])
      }
      d <- dim(values)
      if (is.null(d)) {
        if (length(values) != Ntimesteps)
          stop(name, " wrapped vector must have length Ntimesteps")
        return(values[[timestep]])
      }
      if (tail(d, 1L) != Ntimesteps)
        stop(name, " wrapped array/matrix must have Ntimesteps in its final dimension")
      index <- lapply(d, seq_len)
      index[[length(d)]] <- timestep
      return(do.call(`[`, c(list(values), index, list(drop = TRUE))))
    }

    # Resolver functions are evaluated by the parent model. The returned
    # current value, not timestep itself, is passed to LocalDynamics.
    if (is.function(x)) {
      fm <- names(formals(x))
      call_args <- list(timestep = timestep, Ntimesteps = Ntimesteps)
      if (!is.null(fm) && !("..." %in% fm))
        call_args <- call_args[intersect(names(call_args), fm)]
      return(do.call(x, call_args))
    }

    x
  }

  out <- Map(resolve_one, LocalDynamicsArgs, names(LocalDynamicsArgs))
  names(out) <- names(LocalDynamicsArgs)
  out
}

INApestMetaTransitionMatrixParallel = function(
ModelName,                                      # Model and output name
Nperm,                                          # Number of stochastic simulation runs
Ntimesteps,                                     # Timesteps in each simulation
Nstages,                                        # Number of demographic/life stages
Weights,                                        # Stage weights used for population/capacity totals
Transition,                                     # Stage-transition matrix or node/time-varying equivalent
LocalDynamics = local.dynamics.transition.matrix, # Local population growth, movement and management function
LocalDynamicsArgs = list(),                     # Named extra arguments passed to LocalDynamics
DetectionProb,                                  # Background-surveillance detection probability
DetectionSD = NULL,                             # Variation in background detection probability
ManageProb,                                     # Management probability when information is available
ManageSD = NULL,                                # Variation in management probability
MortalityProb,                                  # Management-driven host mortality probability
MortalitySD = NULL,                             # Variation in management mortality probability
FecundityReduction = 0,                         # Proportional reduction in reproduction under management
SpreadReduction,                                # Proportional reduction in dispersal under management
SpreadReductionSD = NULL,                       # Variation in spread reduction
InitialPopulation = NA,                         # Starting host abundance by node x stage
InitBioP = NA,                                  # Proportion of nodes initially invaded
InvasionRisk = NA,                              # External invasion probability or node weighting
InitialInfo = NA,                               # Starting information state
InitInfoP = NA,                                 # Proportion of nodes initially with information
ExternalInfoProb = 0.0,                         # Information arriving from outside the modelled system
InfoRetentionProb = 1,                          # Probability existing information persists one timestep
InfoPersistenceSteps = NA,                      # Programmed timesteps information persists after evidence
EnvEstabProb = 1,                               # Environmental establishment probability
K,                                              # Host carrying capacity by node
SeedbankK,                                      # Stage-1 / seedbank carrying capacity by node
PropaguleEstablishment,                         # Establishment probability for arriving propagules
IncursionStartPop=NA,                           # Host abundance assigned to new external incursions
SDDprob,                                        # Short-distance source-to-target dispersal
SEAM = 0,                                       # Information-transfer adjacency between nodes
LDDprob = NA,                                   # Long-distance source-to-target dispersal
LDDrate = 0,                                    # Fraction of propagules entering long-distance dispersal
TransitionSDDprob = NULL,                       # Short-distance movement during stage progression
TransitionLDDprob = NULL,                       # Long-distance movement during stage progression
TransitionLDDrate = 0,                          # LDD fraction for moving stage transitions
DispersalDensityFactor = 0,                     # Strength of density-dependent dispersal redistribution
BlockedTransitionMortality = 0,                 # Mortality when stage progression is blocked
OngoingExternalInvasion = F,                    # Allow new host incursions after initialisation
OngoingExternalInfo = F,                        # Allow new external information after initialisation
OutputDir = NA,                                 # Directory for saved outputs
DoPlots = TRUE,                                 # Legacy plotting option; plotting is post-processing
Pathogen = NULL,                                # Optional pathogen process specification
InitialPathogenState = NULL,                    # Starting node x stage x pathogen-state counts
Cores = NULL,                                   # Number of worker processes
Seed = NULL,                                    # Random seed for reproducible simulations
ExternalPathogenStateProb = NULL,               # Pathogen-state distribution for external host arrivals
InfoTriggeredDetectionProb = 0,                 # Detection probability where information already exists
InfoTriggeredDetectionSD = NULL,                # Variation in information-triggered detection
ReturnResults = FALSE                           # Return the in-memory result object
)
{

# ---------------------------------------------------------------------------
# Set up and validate the node x life-stage simulation.
# ---------------------------------------------------------------------------
if(!is.function(LocalDynamics))
  stop("LocalDynamics must be a function")
if(!is.null(Pathogen) && identical(LocalDynamics, local.dynamics.transition.matrix)) {
  if(!exists("local.dynamics.transition.matrix.pathogen", mode="function"))
    stop("Source INApestPathogenTransitionMatrix.R before using Pathogen")
  LocalDynamics <- local.dynamics.transition.matrix.pathogen
}
if (is.null(LocalDynamicsArgs)) LocalDynamicsArgs <- list()
if (!is.list(LocalDynamicsArgs))
  stop("LocalDynamicsArgs must be a named list")
if (length(LocalDynamicsArgs) &&
    (is.null(names(LocalDynamicsArgs)) || any(!nzchar(names(LocalDynamicsArgs)))))
  stop("Every LocalDynamicsArgs entry must have a non-empty name")
if (anyDuplicated(names(LocalDynamicsArgs)))
  stop("LocalDynamicsArgs names must be unique")
# Force the argument before any parallel worker closure is created. This keeps
# the selected default or user-supplied function as an explicit model input.
force(LocalDynamics)
force(LocalDynamicsArgs)
UserLocalDynamicsArgs <- LocalDynamicsArgs
# collect lexical bindings needed by custom LocalDynamics on PSOCK workers.
.CollectLocalDynamicsPSOCKBindings <- function(fun) {
  collected <- list(); seen <- character(0)
  find_user_binding <- function(nm,env) {
    ee <- env
    while(!identical(ee,emptyenv())) {
      enm <- environmentName(ee)
      if(isNamespace(ee) || identical(ee,baseenv()) || startsWith(enm,"package:")) return(list(found=FALSE,value=NULL))
      if(exists(nm,envir=ee,inherits=FALSE)) return(list(found=TRUE,value=get(nm,envir=ee,inherits=FALSE)))
      ee <- parent.env(ee)
    }
    list(found=FALSE,value=NULL)
  }
  collect_fun <- function(f) {
    if(!is.function(f)) return(invisible(NULL))
    fg <- codetools::findGlobals(f,merge=FALSE)
    refs <- unique(c(fg$variables,fg$functions))
    for(nm in refs) {
      if(nm %in% seen) next
      hit <- find_user_binding(nm,environment(f))
      if(!isTRUE(hit$found)) next
      val <- hit$value
      seen <<- c(seen,nm); collected[[nm]] <<- val
      if(is.function(val)) collect_fun(val)
    }
    invisible(NULL)
  }
  collect_any <- function(x) {
    if(is.function(x)) return(collect_fun(x))
    if(is.list(x)) for(v in x) collect_any(v)
    invisible(NULL)
  }
  collect_any(fun); collect_any(UserLocalDynamicsArgs); collected
}
LocalDynamicsPSOCKBindings <- .CollectLocalDynamicsPSOCKBindings(LocalDynamics)
force(LocalDynamicsPSOCKBindings)
# Validate the host-engine contract before starting stochastic simulation.
n_nodes <- nrow(SDDprob)
if(!is.numeric(Nperm) || length(Nperm)!=1L || !is.finite(Nperm) || Nperm<1 || Nperm!=floor(Nperm)) stop("Nperm must be a positive integer")
if(!is.numeric(Ntimesteps) || length(Ntimesteps)!=1L || !is.finite(Ntimesteps) || Ntimesteps<1 || Ntimesteps!=floor(Ntimesteps)) stop("Ntimesteps must be a positive integer")
if(!is.numeric(Nstages) || length(Nstages)!=1L || !is.finite(Nstages) || Nstages<2 || Nstages!=floor(Nstages)) stop("Nstages must be an integer of at least 2")
if(is.null(dim(SDDprob))) stop("SDDprob must be a nodes x nodes matrix or nodes x nodes x Ntimesteps array")
n_nodes <- dim(SDDprob)[1L]
# Validate fixed or time-varying node connectivity.
.ValidateConnectivityTM <- function(x,name,allow_disabled_scalar=FALSE) {
  if(allow_disabled_scalar && length(x)==1L && (is.na(x) || identical(as.numeric(x),0))) return(invisible(TRUE))
  d <- dim(x)
  if(is.null(d)) stop(name," must be a nodes x nodes matrix or nodes x nodes x Ntimesteps array")
  if(length(d)==2L) {
    if(!all(d==c(n_nodes,n_nodes))) stop(name," matrix must have dimensions nodes x nodes")
    mats <- list(x)
  } else if(length(d)==3L) {
    if(!all(d==c(n_nodes,n_nodes,Ntimesteps))) stop(name," 3D array must have dimensions nodes x nodes x Ntimesteps")
    mats <- lapply(seq_len(Ntimesteps),function(tt) matrix(x[,,tt],nrow=n_nodes,ncol=n_nodes))
  } else stop(name," must be a matrix or 3D array")
  for(M in mats) {
    if(any(!is.finite(M)) || any(M<0)) stop(name," must contain finite non-negative probabilities")
    if(any(rowSums(M)>1+1e-10)) stop(name," source-row sums must not exceed 1")
  }
  invisible(TRUE)
}
.ValidateConnectivityTM(SDDprob,"SDDprob",FALSE)
.ValidateConnectivityTM(LDDprob,"LDDprob",TRUE)
# Return the connectivity matrix used in the current timestep.
.SliceConnectivityTM <- function(x,timestep) {
  if(length(x)==1L && is.na(x)) return(x)
  if(length(dim(x))==3L) return(matrix(x[,,timestep],nrow=n_nodes,ncol=n_nodes))
  x
}
# Validate supported demographic transition-matrix forms.
.ValidateTransitionTM <- function(x) {
  if(is.matrix(x)) {
    if(!all(dim(x)==c(Nstages,Nstages))) stop("Transition matrix must have dimensions Nstages x Nstages")
    if(any(!is.finite(x)) || any(x<0)) stop("Transition must contain finite non-negative values")
    return(invisible(TRUE))
  }
  if(is.list(x)) {
    if(length(x)!=n_nodes) stop("Transition list must contain one Nstages x Nstages matrix per node")
    for(ii in seq_len(n_nodes)) {
      if(!is.matrix(x[[ii]]) || !all(dim(x[[ii]])==c(Nstages,Nstages))) stop("Every Transition list entry must be an Nstages x Nstages matrix")
      if(any(!is.finite(x[[ii]])) || any(x[[ii]]<0)) stop("Transition must contain finite non-negative values")
    }
    return(invisible(TRUE))
  }
  d <- dim(x)
  if(length(d)==3L) {
    if(!all(d==c(Nstages,Nstages,Ntimesteps))) stop("3D Transition array must have dimensions Nstages x Nstages x Ntimesteps")
  } else if(length(d)==4L) {
    if(!all(d==c(Nstages,Nstages,n_nodes,Ntimesteps))) stop("4D Transition array must have dimensions Nstages x Nstages x nodes x Ntimesteps")
  } else stop("Transition must be a matrix, a per-node list, a Nstages x Nstages x Ntimesteps array, or a Nstages x Nstages x nodes x Ntimesteps array")
  if(any(!is.finite(x)) || any(x<0)) stop("Transition must contain finite non-negative values")
  invisible(TRUE)
}
# Resolve demographic transition matrices for the current timestep.
.ResolveTransitionTM <- function(x,timestep) {
  if(is.matrix(x) || is.list(x)) return(x)
  d <- dim(x)
  if(length(d)==3L) return(matrix(x[,,timestep],nrow=Nstages,ncol=Nstages))
  lapply(seq_len(n_nodes),function(ii) matrix(x[,,ii,timestep],nrow=Nstages,ncol=Nstages))
}
.ValidateTransitionTM(Transition)
if(is.matrix(Weights)) {
  if(!all(dim(Weights)==c(n_nodes,Nstages))) stop("Weights matrix must have dimensions nodes x Nstages")
} else if(length(Weights)!=Nstages) stop("Weights must be length Nstages or a nodes x Nstages matrix")
if(any(!is.finite(Weights)) || any(Weights<=0)) stop("Weights must contain finite values greater than zero")
# Validate node-level scalar/vector/time-varying inputs.
.ValidateNodeTimeTM <- function(x,name,unit_interval=FALSE,allow_na_scalar=FALSE,whole=FALSE) {
  if(allow_na_scalar && length(x)==1L && is.na(x)) return(invisible(TRUE))
  d <- dim(x); ok <- FALSE
  if(is.null(d)) ok <- length(x) %in% c(1L,n_nodes)
  if(!is.null(d) && length(d)==2L) ok <- all(d==c(n_nodes,Ntimesteps))
  if(!ok) stop(name," must be scalar, length nodes, or nodes x Ntimesteps")
  if(any(!is.finite(x))) stop(name," must contain finite values")
  if(any(x<0)) stop(name," must contain non-negative values")
  if(unit_interval && any(x>1)) stop(name," values must be between 0 and 1")
  if(whole && any(x!=floor(x))) stop(name," must contain whole numbers")
  invisible(TRUE)
}
.ValidateNodeTimeTM(K,"K",FALSE,FALSE,FALSE)
.ValidateNodeTimeTM(SeedbankK,"SeedbankK",FALSE,FALSE,FALSE)
.ValidateNodeTimeTM(EnvEstabProb,"EnvEstabProb",TRUE,FALSE,FALSE)
.ValidateNodeTimeTM(PropaguleEstablishment,"PropaguleEstablishment",FALSE,FALSE,FALSE)
.ValidateNodeTimeTM(ManageProb,"ManageProb",TRUE,FALSE,FALSE)
.ValidateNodeTimeTM(SpreadReduction,"SpreadReduction",TRUE,FALSE,FALSE)
.ValidateNodeTimeTM(ExternalInfoProb,"ExternalInfoProb",TRUE,TRUE,FALSE)
.ValidateNodeTimeTM(InvasionRisk,"InvasionRisk",TRUE,TRUE,FALSE)
if(!is.na(InitBioP) && (!is.numeric(InitBioP) || length(InitBioP)!=1L || !is.finite(InitBioP) || InitBioP<0 || InitBioP>1)) stop("InitBioP must be NA or one probability between 0 and 1")
if(!(length(IncursionStartPop)==1L && is.na(IncursionStartPop)) && (!is.numeric(IncursionStartPop) || length(IncursionStartPop)!=1L || !is.finite(IncursionStartPop) || IncursionStartPop<0 || IncursionStartPop!=floor(IncursionStartPop))) stop("IncursionStartPop must be NA or one non-negative whole number")
# Validate stage-specific probability inputs.
.ValidateStageProbTM <- function(x,name) {
  d <- dim(x); ok <- FALSE
  if(is.null(d)) ok <- length(x) %in% c(1L,Nstages)
  if(!is.null(d) && length(d)==2L) ok <- all(d==c(n_nodes,Nstages))
  if(!is.null(d) && length(d)==3L) ok <- all(d==c(n_nodes,Nstages,Ntimesteps))
  if(!ok) stop(name," must be scalar, length Nstages, nodes x Nstages, or nodes x Nstages x Ntimesteps")
  if(any(!is.finite(x)) || any(x<0 | x>1)) stop(name," values must be probabilities between 0 and 1")
  invisible(TRUE)
}
.ValidateStageProbTM(DetectionProb,"DetectionProb")
.ValidateStageProbTM(MortalityProb,"MortalityProb")
if(!(length(InitialPopulation)==1L && is.na(InitialPopulation))) {
  if(!is.matrix(InitialPopulation) || !all(dim(InitialPopulation)==c(n_nodes,Nstages))) stop("InitialPopulation must be NA or a nodes x Nstages matrix")
  if(any(!is.finite(InitialPopulation)) || any(InitialPopulation<0) || any(InitialPopulation!=floor(InitialPopulation))) stop("InitialPopulation must contain non-negative whole numbers")
}
if(isTRUE(OngoingExternalInfo) && length(ExternalInfoProb)==1L && is.na(ExternalInfoProb)) stop("OngoingExternalInfo requires finite ExternalInfoProb")
if(isTRUE(OngoingExternalInvasion) && length(InvasionRisk)==1L && is.na(InvasionRisk)) stop("OngoingExternalInvasion requires finite InvasionRisk")
if(!is.numeric(LDDrate) || length(LDDrate)!=1L || !is.finite(LDDrate) || LDDrate<0 || LDDrate>1) stop("LDDrate must be one probability between 0 and 1")
if(!is.numeric(TransitionLDDrate) || any(!is.finite(TransitionLDDrate)) || any(TransitionLDDrate<0 | TransitionLDDrate>1)) stop("TransitionLDDrate values must be probabilities between 0 and 1")
if(!(length(DispersalDensityFactor)==1L && (is.na(DispersalDensityFactor) || (is.finite(DispersalDensityFactor) && DispersalDensityFactor>=0)))) stop("DispersalDensityFactor must be NA, zero, or one finite positive value")
# Host-engine contract validation is complete.


# Optional pathogen-state composition of external stage-1 host immigrants.
# This does not change the timing, number or demographic-stage placement of
# external incursions; it only records the pathogen state of hosts that arrive.
ExternalPathogenStateProbResolved <- NULL
if(!is.null(ExternalPathogenStateProb)) {
  if(is.null(Pathogen))
    stop("ExternalPathogenStateProb requires Pathogen to be supplied")
  if(!is.numeric(ExternalPathogenStateProb) || is.null(names(ExternalPathogenStateProb)) ||
     any(!nzchar(names(ExternalPathogenStateProb))) || anyDuplicated(names(ExternalPathogenStateProb)))
    stop("ExternalPathogenStateProb must be a uniquely named numeric vector of pathogen-state probabilities")
  if(any(!names(ExternalPathogenStateProb) %in% Pathogen$States))
    stop("ExternalPathogenStateProb names must be pathogen states: ", paste(Pathogen$States, collapse = ", "))
  if(any(!is.finite(ExternalPathogenStateProb)) || any(ExternalPathogenStateProb < 0))
    stop("ExternalPathogenStateProb values must be finite and non-negative")
  ExternalPathogenStateProbResolved <- setNames(rep(0, length(Pathogen$States)), Pathogen$States)
  ExternalPathogenStateProbResolved[names(ExternalPathogenStateProb)] <- ExternalPathogenStateProb
  if(abs(sum(ExternalPathogenStateProbResolved) - 1) > 1e-10)
    stop("ExternalPathogenStateProb must sum to 1")
}
force(ExternalPathogenStateProbResolved)

# Connectivity was validated with the other host-engine inputs above.

# Optional movement while individuals progress between life-history stages.
# A bare matrix/3D array applies to every progression transition. A list of
# length Nstages-1 allows each source-stage transition to have its own matrix,
# time-varying array, or NULL (local-only transition).
n_nodes <- nrow(SDDprob)
# Validate stage-transition movement inputs.
ValidateTransitionMovement <- function(x, name) {
  validate_one <- function(z, label) {
    if(is.null(z)) return(invisible(NULL))
    dz <- dim(z)
    if(is.null(dz))
      stop(label, " must be a nodes x nodes matrix or nodes x nodes x Ntimesteps array")
    if(length(dz) == 2L) {
      if(!all(dz == c(n_nodes, n_nodes)))
        stop(label, " matrix must have dimensions nodes x nodes")
    } else if(length(dz) == 3L) {
      if(!all(dz == c(n_nodes, n_nodes, Ntimesteps)))
        stop(label, " 3D array must have dimensions nodes x nodes x Ntimesteps")
    } else {
      stop(label, " must be a matrix or 3D array")
    }
    if(any(!is.finite(z)) || any(z < 0))
      stop(label, " must contain finite non-negative probabilities")
    if(length(dz) == 2L) {
      if(any(rowSums(z) > 1 + 1e-10))
        stop(label, " row sums must not exceed 1")
    } else {
      for(tt in seq_len(Ntimesteps))
        if(any(rowSums(matrix(z[,,tt],nrow=n_nodes,ncol=n_nodes)) > 1 + 1e-10))
          stop(label, " row sums must not exceed 1 in any timestep")
    }
    invisible(NULL)
  }

  if(is.null(x)) return(invisible(NULL))
  if(is.list(x)) {
    if(length(x) != Nstages - 1L)
      stop(name, " list must have length Nstages - 1")
    for(k in seq_len(Nstages - 1L))
      validate_one(x[[k]], paste0(name, "[[", k, "]]"))
  } else {
    validate_one(x, name)
  }
  invisible(NULL)
}

# Resolve stage-transition movement for the current timestep.
ResolveTransitionMovement <- function(x, timestep) {
  resolve_one <- function(z) {
    if(is.null(z)) return(NULL)
    if(length(dim(z)) == 3L) return(matrix(z[,,timestep],nrow=n_nodes,ncol=n_nodes))
    z
  }
  if(is.null(x)) return(NULL)
  if(is.list(x)) return(lapply(x, resolve_one))
  resolve_one(x)
}

# Check whether stage-transition movement is active.
HasTransitionMovement <- function(x) {
  if(is.null(x)) return(FALSE)
  if(is.list(x)) return(any(!vapply(x, is.null, logical(1))))
  TRUE
}

ValidateTransitionMovement(TransitionSDDprob, "TransitionSDDprob")
ValidateTransitionMovement(TransitionLDDprob, "TransitionLDDprob")
if(!(length(TransitionLDDrate) %in% c(1L, Nstages - 1L)))
  stop("TransitionLDDrate must be scalar or length Nstages - 1")
if(any(!is.finite(TransitionLDDrate)) || any(TransitionLDDrate < 0) || any(TransitionLDDrate > 1))
  stop("TransitionLDDrate values must be between 0 and 1")
TransitionMovementConfigured <- HasTransitionMovement(TransitionSDDprob) || HasTransitionMovement(TransitionLDDprob)

# Allow information retention to vary by node and through time
if(is.matrix(InfoRetentionProb) == T && (nrow(InfoRetentionProb) != nrow(SDDprob) || ncol(InfoRetentionProb) != Ntimesteps))
  stop("InfoRetentionProb matrix must have dimensions nodes x Ntimesteps")
if(is.matrix(InfoRetentionProb) == F && !(length(InfoRetentionProb) == 1 || length(InfoRetentionProb) == nrow(SDDprob)))
  stop("InfoRetentionProb must be a single value, vector of length nodes, or matrix nodes x Ntimesteps")
if(any(is.na(InfoRetentionProb)) || any(InfoRetentionProb < 0) || any(InfoRetentionProb > 1))
  stop("InfoRetentionProb values must be between 0 and 1")

# Allow programmed information persistence after last known local presence to vary by node and through time
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

# Use first-timestep carrying capacity when K varies through time.
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

# Allocate host-abundance histories.
PopulationResults = array(dim = c(nrow(SDDprob),Ntimesteps,Nperm))

# Allocate host-abundance histories.
PopulationStageResults = array(dim = c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))

PathogenStageResults <- NULL
PathogenDeathResults <- NULL
NewInfectionResults <- NULL
PathogenDetectedResults <- NULL
if(!is.null(Pathogen)) {
  if(!inherits(Pathogen, "INApestPathogen")) stop("Pathogen must be created by INApestPathogen()")
  PathogenStageResults <- array(0, dim=c(nrow(SDDprob),Nstages,length(Pathogen$States),Ntimesteps,Nperm), dimnames=list(NULL,NULL,Pathogen$States,NULL,NULL))
  PathogenDeathResults <- array(0, dim=c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))
  NewInfectionResults <- array(0, dim=c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))
  PathogenDetectedResults <- array(0L, dim=c(nrow(SDDprob),Ntimesteps,Nperm))
}

# Allocate host/pest presence histories.
InvasionResults = array(dim = c(nrow(SDDprob),Ntimesteps,Nperm))


# Allocate known-presence and surveillance histories.
DetectedResults = InvasionResults
BackgroundDetectedResults = array(0L,dim=c(nrow(SDDprob),Ntimesteps,Nperm))
InfoTriggeredDetectedResults = array(0L,dim=c(nrow(SDDprob),Ntimesteps,Nperm))
BackgroundDetectionProbabilityResults = array(0,dim=c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))
InfoTriggeredDetectionProbabilityResults = array(0,dim=c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))
InformationStateBeforeSurveillanceResults = array(0L,dim=c(nrow(SDDprob),Ntimesteps,Nperm))
HaveInfoResults = array(0L,dim=c(nrow(SDDprob),Ntimesteps,Nperm))

# Allocate management histories.
ManagingResults = InvasionResults


# Validate the socioeconomic information-transfer network when supplied.
if(is.matrix(SEAM) == T)
     {
     diag(SEAM) = 0
     RandSEAM <- matrix(NA,nrow = nrow(SDDprob),ncol=nrow(SDDprob))
     }

# Validate management-induced fecundity reduction
n_nodes <- nrow(SDDprob)
FecundityReductionDims <- dim(FecundityReduction)
if(is.null(FecundityReductionDims)) {
  valid_lengths <- unique(c(1L, n_nodes, Nstages))
  if(!(length(FecundityReduction) %in% valid_lengths))
    stop("FecundityReduction must be scalar, length nodes, length Nstages, nodes x Ntimesteps, or nodes x Nstages x Ntimesteps")
  if(length(FecundityReduction) > 1 && n_nodes == Nstages && length(FecundityReduction) == n_nodes)
    stop("FecundityReduction vector is ambiguous because number of nodes equals Nstages; use a nodes x Ntimesteps matrix for node-specific values or a nodes x Nstages x Ntimesteps array for stage-specific values")
} else if(length(FecundityReductionDims) == 2) {
  if(!all(FecundityReductionDims == c(n_nodes, Ntimesteps)))
    stop("FecundityReduction matrix must have dimensions nodes x Ntimesteps")
} else if(length(FecundityReductionDims) == 3) {
  if(!all(FecundityReductionDims == c(n_nodes, Nstages, Ntimesteps)))
    stop("FecundityReduction array must have dimensions nodes x Nstages x Ntimesteps")
} else {
  stop("FecundityReduction must be scalar, length nodes, length Nstages, nodes x Ntimesteps, or nodes x Nstages x Ntimesteps")
}
if(any(!is.finite(FecundityReduction)) || any(FecundityReduction < 0) || any(FecundityReduction > 1))
  stop("FecundityReduction values must be between 0 and 1")

# Resolve fecundity reduction to node x stage values.
ResolveFecundityReductionTM <- function(timestep) {
  if(is.null(dim(FecundityReduction))) {
    if(length(FecundityReduction) == 1)
      return(matrix(FecundityReduction, nrow = n_nodes, ncol = Nstages))
    if(length(FecundityReduction) == n_nodes)
      return(matrix(rep(FecundityReduction, Nstages), nrow = n_nodes, ncol = Nstages))
    return(matrix(rep(FecundityReduction, each = n_nodes), nrow = n_nodes, ncol = Nstages))
  }
  if(length(dim(FecundityReduction)) == 2)
    return(matrix(rep(FecundityReduction[,timestep], Nstages), nrow = n_nodes, ncol = Nstages))
  FecundityReduction[,,timestep]
}

# Assign standard deviation value to management in no value provided
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

# --- Information-triggered detection SD and shape validation ---
.InfoTriggeredTMShapeOK <- function(x) {
  d <- dim(x)
  if(is.null(d)) return(length(x) == 1 || length(x) == Nstages)
  if(length(d) == 2) return(all(d == c(n_nodes,Nstages)))
  if(length(d) == 3) return(all(d == c(n_nodes,Nstages,Ntimesteps)))
  FALSE
}
if(!.InfoTriggeredTMShapeOK(InfoTriggeredDetectionProb))
  stop("InfoTriggeredDetectionProb must be scalar, length Nstages, nodes x stages, or nodes x stages x Ntimesteps")
if(any(!is.finite(InfoTriggeredDetectionProb)) || any(InfoTriggeredDetectionProb < 0) || any(InfoTriggeredDetectionProb > 1))
  stop("InfoTriggeredDetectionProb values must be between 0 and 1")
if (is.null(InfoTriggeredDetectionSD)) {
  if (is.matrix(InfoTriggeredDetectionProb)) {
    stage_means <- colMeans(InfoTriggeredDetectionProb, na.rm = TRUE)
    InfoTriggeredDetectionSD <- matrix(stage_means / 10, nrow=n_nodes, ncol=Nstages, byrow=TRUE)
  } else {
    InfoTriggeredDetectionSD <- matrix(mean(InfoTriggeredDetectionProb, na.rm=TRUE)/10, nrow=n_nodes, ncol=Nstages)
  }
}
# Check supported shapes for targeted-detection uncertainty.
.InfoTriggeredTMSDShapeOK <- function(x) {
  d <- dim(x)
  if(is.null(d)) return(length(x) == 1 || length(x) == Nstages)
  if(length(d) == 2) return(all(d == c(n_nodes,Nstages)))
  FALSE
}
if(!.InfoTriggeredTMSDShapeOK(InfoTriggeredDetectionSD))
  stop("InfoTriggeredDetectionSD must be scalar, length Nstages, or nodes x stages; time variation belongs in InfoTriggeredDetectionProb")
if(any(!is.finite(InfoTriggeredDetectionSD)) || any(InfoTriggeredDetectionSD < 0))
  stop("InfoTriggeredDetectionSD must contain finite non-negative values")
UseInfoTriggeredSurveillance <- any(InfoTriggeredDetectionProb != 0) || any(InfoTriggeredDetectionSD != 0)

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
### Start of simulation
###########################################################
    
# ---- Parallel permutation loop ----
detected_cores <- parallel::detectCores()
if (is.na(detected_cores)) detected_cores <- 2L
if (is.null(Cores)) {
  n_cores <- max(1L, min(Nperm, detected_cores - 1L))
} else {
  if(!is.numeric(Cores) || length(Cores)!=1L || !is.finite(Cores) || Cores < 1 || Cores != floor(Cores)) stop("Cores must be a positive integer or NULL")
  n_cores <- max(1L, min(Nperm, as.integer(Cores)))
}
if(!is.null(Seed)) {
  if(!is.numeric(Seed) || length(Seed)!=1L || !is.finite(Seed)) stop("Seed must be one finite number or NULL")
  set.seed(as.integer(Seed))
}

# Run one stochastic realisation. Function arguments and LocalDynamics are
# captured in this closure, so new parameters do not require a clusterExport list.
# Capture the LocalDynamicsArgs resolver in this call environment so Windows
# PSOCK workers do not depend on a helper that exists only in the master session.
LocalDynamicsArgsResolver <- .resolve_INApest_LocalDynamicsArgs
force(LocalDynamicsArgsResolver)

# ---------------------------------------------------------------------------
# Run one complete stochastic stage-structured history on a worker.
# ---------------------------------------------------------------------------
PermutationWorker <- function(i_perm) {
  
  # Max integer for propagule dispersal using rmultinom
  MaxInteger <- .Machine$integer.max  
  # Local containers for this permutation
  n_nodes <- nrow(SDDprob)
  PopulationResults_local      <- matrix(0, nrow = n_nodes, ncol = Ntimesteps)
  PopulationStageResults_local <- array(0, dim = c(n_nodes, Nstages, Ntimesteps))
  PathogenStageResults_local <- if(is.null(Pathogen)) NULL else array(0, dim=c(n_nodes,Nstages,length(Pathogen$States),Ntimesteps), dimnames=list(NULL,NULL,Pathogen$States,NULL))
  PathogenDeathResults_local <- if(is.null(Pathogen)) NULL else array(0, dim=c(n_nodes,Nstages,Ntimesteps))
  NewInfectionResults_local <- if(is.null(Pathogen)) NULL else array(0, dim=c(n_nodes,Nstages,Ntimesteps))
  PathogenDetectedResults_local <- if(is.null(Pathogen)) NULL else matrix(0L,nrow=n_nodes,ncol=Ntimesteps)
  InvasionResults_local        <- matrix(0, nrow = n_nodes, ncol = Ntimesteps)
  DetectedResults_local        <- matrix(0, nrow = n_nodes, ncol = Ntimesteps)
  BackgroundDetectedResults_local <- matrix(0L,nrow=n_nodes,ncol=Ntimesteps)
  InfoTriggeredDetectedResults_local <- matrix(0L,nrow=n_nodes,ncol=Ntimesteps)
  BackgroundDetectionProbabilityResults_local <- array(0,dim=c(n_nodes,Nstages,Ntimesteps))
  InfoTriggeredDetectionProbabilityResults_local <- array(0,dim=c(n_nodes,Nstages,Ntimesteps))
  InformationStateBeforeSurveillanceResults_local <- matrix(0L,nrow=n_nodes,ncol=Ntimesteps)
  HaveInfoResults_local <- matrix(0L,nrow=n_nodes,ncol=Ntimesteps)
  ManagingResults_local        <- matrix(0, nrow = n_nodes, ncol = Ntimesteps)
  
  # ---------- BEGIN exactly your block (minimal edits) ----------
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
  PathogenStageState <- NULL
  if(!is.null(Pathogen)) {
    if(!exists("INApestPathogenStageState", mode="function")) stop("Source INApestPathogenTransitionMatrix.R before using Pathogen")
    PathogenStageState <- INApestPathogenStageState(N, Pathogen, InitialPathogenState, Ntimesteps)
  }
  if (sum(N) == 0 && OngoingExternalInvasion == FALSE) {
    warning("No initial populations and no future external invasions")
  }
  
  # Robust InitInfo setup
n_nodes <- nrow(SDDprob)
HasExplicitInitialInfo <- !(length(InitialInfo)==1L && is.na(InitialInfo))
if(HasExplicitInitialInfo) {
  if(!is.numeric(InitialInfo) || length(InitialInfo)!=n_nodes || any(!is.finite(InitialInfo)) || any(!InitialInfo %in% c(0,1)))
    stop("InitialInfo must be scalar NA or a binary vector of length nodes")
  InitInfo <- as.integer(InitialInfo)
} else {
  InitInfo <- integer(n_nodes)
  ExternalInfoMissing <- length(ExternalInfoProb)==1L && is.na(ExternalInfoProb)
  ext0 <- if(ExternalInfoMissing) { if(!is.na(InitInfoP)) rep(1,n_nodes) else rep(0,n_nodes) } else if(is.matrix(ExternalInfoProb)) ExternalInfoProb[,1] else rep_len(ExternalInfoProb,n_nodes)
  if(any(!is.finite(ext0)) || any(ext0<0 | ext0>1)) stop("ExternalInfoProb values must be probabilities between 0 and 1")
  if(!is.na(InitInfoP)) {
    if(!is.numeric(InitInfoP) || length(InitInfoP)!=1L || !is.finite(InitInfoP) || InitInfoP<0 || InitInfoP>1) stop("InitInfoP must be NA or one probability between 0 and 1")
    n_sample <- min(n_nodes,ceiling(n_nodes*InitInfoP))
    if(n_sample>0L) {
      weights <- ext0
      if(sum(weights)<=0) weights <- rep(1,n_nodes)
      Info <- sample.int(n_nodes,size=n_sample,prob=weights)
      InitInfo[Info] <- 1L
    }
  } else if(any(ext0>0)) {
    InitInfo <- rbinom(n_nodes,size=1,prob=ext0)
  }
}


# --- Determine dimensions and pre-sample detection/mortality arrays (as in your code) ---
  # Build NodeDetectionProb array (n_nodes x Nstages x Ntimesteps)
  # Randomly assign detection probabilities
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
# Randomly assign information-triggered detection probabilities. This is an
# additive observation pathway; when disabled no additional RNG is consumed.
if(UseInfoTriggeredSurveillance) {
  if (!is.array(InfoTriggeredDetectionProb) && (length(InfoTriggeredDetectionProb) == 1 || length(InfoTriggeredDetectionProb) == Nstages)) {
    NodeInfoTriggeredDetectionProb <- array(NA, dim=c(n_nodes,Nstages,Ntimesteps))
    for (s in 1:Nstages) for (t in 1:Ntimesteps)
      NodeInfoTriggeredDetectionProb[,s,t] <- pmax(0,pmin(1,rnorm(n_nodes,
        mean=if(length(InfoTriggeredDetectionProb)==1) InfoTriggeredDetectionProb else InfoTriggeredDetectionProb[s],
        sd=if(is.matrix(InfoTriggeredDetectionSD)) InfoTriggeredDetectionSD[,s] else if(length(InfoTriggeredDetectionSD)==1) InfoTriggeredDetectionSD else InfoTriggeredDetectionSD[s])))
  } else if (is.matrix(InfoTriggeredDetectionProb) && all(dim(InfoTriggeredDetectionProb)==c(n_nodes,Nstages))) {
    NodeInfoTriggeredDetectionProb <- array(NA,dim=c(n_nodes,Nstages,Ntimesteps))
    for(t in 1:Ntimesteps) NodeInfoTriggeredDetectionProb[,,t] <- pmax(0,pmin(1,matrix(rnorm(n_nodes*Nstages,
      mean=as.vector(InfoTriggeredDetectionProb),sd=as.vector(InfoTriggeredDetectionSD)),nrow=n_nodes,ncol=Nstages)))
  } else if (length(dim(InfoTriggeredDetectionProb))==3 && all(dim(InfoTriggeredDetectionProb)==c(n_nodes,Nstages,Ntimesteps))) {
    NodeInfoTriggeredDetectionProb <- array(NA,dim=dim(InfoTriggeredDetectionProb))
    for(t in 1:Ntimesteps) NodeInfoTriggeredDetectionProb[,,t] <- pmax(0,pmin(1,matrix(rnorm(n_nodes*Nstages,
      mean=as.vector(InfoTriggeredDetectionProb[,,t]),sd=as.vector(InfoTriggeredDetectionSD)),nrow=n_nodes,ncol=Nstages)))
  }
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
  prob_detect <- 1 - apply((1 - matrix(NodeDetectionProb[,,1],nrow=n_nodes,ncol=Nstages))^InitBio, 1, prod)
  InitDetection <- rbinom(n = nrow(InitBio), size = 1, prob = prob_detect)
  InitInfo[InitInfo == 0] <- InitDetection[InitInfo == 0]
  HaveInfo <- InitInfo
  
  # Track the most recent timestep with known local presence
  LastKnownPresence <- rep(NA, nrow(SDDprob))
  if (UseInfoPersistence == TRUE) {
    InitialKnownPresence <- which(InitDetection == 1)
    if (length(InitialKnownPresence) > 0) LastKnownPresence[InitialKnownPresence] <- 0
  }
  
  # static management/spread uncertainty is one parameter draw per realisation, matching serial.
  if(!is.matrix(ManageProb)) NodeManageProb <- pmin(1,pmax(0,rnorm(n=n_nodes,mean=ManageProb,sd=ManageSD)))
  if(!is.matrix(SpreadReduction)) NodeSpreadReduction <- pmin(1,pmax(0,rnorm(n=n_nodes,mean=SpreadReduction,sd=SpreadReductionSD)))

  # ---------- Now run timesteps (use i_perm in prints) ----------

  # ---------------------------------------------------------------------------
  # Advance stage dynamics, pathogen state, information and response through time.
  # ---------------------------------------------------------------------------
  for (timestep in 1:Ntimesteps) {
    # minimal progress print (safe inside worker)
    cat(sprintf("Worker %d: timestep %d\n", i_perm, timestep))
    
    # Resolve short- and long-distance connectivity for the current timestep.
    NodeSDDprob <- .SliceConnectivityTM(SDDprob,timestep)
    NodeLDDprob <- .SliceConnectivityTM(LDDprob,timestep)
    NodeTransitionSDDprob <- ResolveTransitionMovement(TransitionSDDprob,timestep)
    NodeTransitionLDDprob <- ResolveTransitionMovement(TransitionLDDprob,timestep)
    NodeFecundityReduction <- ResolveFecundityReductionTM(timestep)

    # update NodeEnvEstabProb, NodeTransition, NodeK, NodePropaguleEstablishment if matrices/arrays
    if (is.matrix(EnvEstabProb) == TRUE) NodeEnvEstabProb <- EnvEstabProb[, timestep]
    NodeTransition <- .ResolveTransitionTM(Transition,timestep)

    if (is.matrix(K) == TRUE) NodeK <- K[, timestep]
    
    # If seedbank carrying capacity provided as matrix assign values for relevant timestep
    if (is.matrix(SeedbankK) == TRUE) NodeSeedbankK <- SeedbankK[, timestep]
    if (is.matrix(PropaguleEstablishment) == TRUE) NodePropaguleEstablishment <- PropaguleEstablishment[, timestep]
    
    
    # ManageProb / SpreadReduction per timestep when explicitly time-varying.
    if (is.matrix(ManageProb) == TRUE && nrow(ManageProb) == nrow(SDDprob) && ncol(ManageProb) == Ntimesteps)
      NodeManageProb <- pmin(1,pmax(0,rnorm(n=nrow(SDDprob),mean=ManageProb[,timestep],sd=ManageSD)))
    if (is.matrix(SpreadReduction) == TRUE && nrow(SpreadReduction) == nrow(SDDprob) && ncol(SpreadReduction) == Ntimesteps)
      NodeSpreadReduction <- pmin(1,pmax(0,rnorm(n=nrow(SDDprob),mean=SpreadReduction[,timestep],sd=SpreadReductionSD)))
    
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
    if(!is.null(Pathogen)) PathogenStageState <- .iptm_reconcile(PathogenStageState, N0)
    N <- N0
    if (sum(N0) > 0) {
      LocalDynamicsArgs <- list(nodetransition = NodeTransition, weights = Weights, sddprob = NodeSDDprob,
                                nodeenvestabprob = NodeEnvEstabProb, n0 = N0, lddprob = NodeLDDprob,
                                lddrate = LDDrate, nodeK = NodeK, node.seedbankK = NodeSeedbankK,
                                nodepropaguleestablishment = NodePropaguleEstablishment,
                                nodespreadreduction = NodeSpreadReduction, managing = Managing,
                                MaxInteger = MaxInteger, BlockedTransitionMortality = BlockedTransitionMortality,
                                DispersalDensityFactor = DispersalDensityFactor)
      LocalDynamicsFormals <- names(formals(LocalDynamics))
      LocalDynamicsAcceptsFecundityReduction <- "nodefecundityreduction" %in% LocalDynamicsFormals || "..." %in% LocalDynamicsFormals
      if(LocalDynamicsAcceptsFecundityReduction)
        LocalDynamicsArgs$nodefecundityreduction <- NodeFecundityReduction
      else if(any(NodeFecundityReduction * Managing > 0))
        stop("Custom LocalDynamics must accept a 'nodefecundityreduction' argument (or ...) when FecundityReduction is active")

      if(!is.null(Pathogen)) {
        if(!all(c("pathogen_state","Pathogen","timestep","Ntimesteps") %in% LocalDynamicsFormals) && !("..." %in% LocalDynamicsFormals))
          stop("Pathogen-aware LocalDynamics must accept pathogen_state, Pathogen, timestep and Ntimesteps (or ...)")
        LocalDynamicsArgs$pathogen_state <- PathogenStageState
        LocalDynamicsArgs$Pathogen <- Pathogen
        LocalDynamicsArgs$timestep <- timestep
        LocalDynamicsArgs$Ntimesteps <- Ntimesteps
      }
      LocalDynamicsAcceptsTransitionMovement <-
        "..." %in% LocalDynamicsFormals ||
        all(c("transition_sddprob", "transition_lddprob", "transition_lddrate") %in% LocalDynamicsFormals)
      if(LocalDynamicsAcceptsTransitionMovement) {
        LocalDynamicsArgs$transition_sddprob <- NodeTransitionSDDprob
        LocalDynamicsArgs$transition_lddprob <- NodeTransitionLDDprob
        LocalDynamicsArgs$transition_lddrate <- TransitionLDDrate
      } else if(TransitionMovementConfigured) {
        stop("Custom LocalDynamics must accept 'transition_sddprob', 'transition_lddprob' and 'transition_lddrate' arguments (or ...) when transition movement is active")
      }
      ResolvedLocalDynamicsArgs <- LocalDynamicsArgsResolver(
        UserLocalDynamicsArgs, timestep = timestep, Ntimesteps = Ntimesteps
      )
      if (length(ResolvedLocalDynamicsArgs)) {
        duplicate_args <- intersect(names(ResolvedLocalDynamicsArgs), names(LocalDynamicsArgs))
        if (length(duplicate_args))
          stop("LocalDynamicsArgs may not override INApest core LocalDynamics argument(s): ",
               paste(duplicate_args, collapse = ", "))
        LocalDynamicsFormals <- names(formals(LocalDynamics))
        unknown_args <- setdiff(names(ResolvedLocalDynamicsArgs), LocalDynamicsFormals)
        if (!("..." %in% LocalDynamicsFormals) && length(unknown_args))
          stop("Custom LocalDynamics does not accept LocalDynamicsArgs entry/entries: ",
               paste(unknown_args, collapse = ", "))
        LocalDynamicsArgs <- c(LocalDynamicsArgs, ResolvedLocalDynamicsArgs)
      }
      LocalDynamicsResult <- do.call(LocalDynamics, LocalDynamicsArgs)
      if(!is.null(Pathogen)) {
        if(!is.list(LocalDynamicsResult) || is.null(LocalDynamicsResult$N) || is.null(LocalDynamicsResult$PathogenState)) stop("Pathogen-aware LocalDynamics must return list(N=..., PathogenState=...)")
        N <- LocalDynamicsResult$N; PathogenStageState <- LocalDynamicsResult$PathogenState
        PathogenDeathResults_local[,,timestep] <- if(is.null(LocalDynamicsResult$PathogenDeaths)) 0 else LocalDynamicsResult$PathogenDeaths
        NewInfectionResults_local[,,timestep] <- if(is.null(LocalDynamicsResult$NewInfections)) 0 else LocalDynamicsResult$NewInfections
      } else N <- LocalDynamicsResult
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
    
    # External invasion. Host timing and stage-1 incursion logic are unchanged.
    NBeforeExternalInvasion <- N
    if (OngoingExternalInvasion == TRUE) {
      if (is.matrix(InvasionRisk) == FALSE) {
        ExternalInvasion <- rbinom(1:nrow(SDDprob), size = 1, prob = InvasionRisk)
      } else {
        ExternalInvasion <- rbinom(1:nrow(SDDprob), size = 1, prob = InvasionRisk[, timestep])
      }
      Invaded[Invaded == 0] <- ExternalInvasion[Invaded == 0]
      if (is.na(IncursionStartPop) == TRUE) N[, 1] <- N[, 1] + ExternalInvasion
      else N[, 1] <- N[, 1] + ExternalInvasion * IncursionStartPop
      if(!is.null(Pathogen)) {
        if(is.null(ExternalPathogenStateProbResolved)) {
          PathogenStageState <- .iptm_reconcile(PathogenStageState, N)
        } else {
          NBeforeStage1 <- as.integer(NBeforeExternalInvasion[,1])
          NAfterStage1 <- as.integer(N[,1])
          ExternalAccepted <- pmax(0L, NAfterStage1-NBeforeStage1)
          for(ii in which(ExternalAccepted > 0L)) {
            ExternalByState <- as.integer(rmultinom(1L, size=ExternalAccepted[ii], prob=ExternalPathogenStateProbResolved))
            PathogenStageState[ii,1,Pathogen$States] <-
              PathogenStageState[ii,1,Pathogen$States] + ExternalByState
          }
          storage.mode(PathogenStageState) <- "integer"
          NTarget <- floor(N)
          if(any(apply(PathogenStageState,c(1,2),sum) != NTarget))
            stop("External host pathogen-state assignment violated stage-wise pathogen totals = N")
        }
      }
      
      
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
    if(!is.null(Pathogen)) {
      PathogenStageResults_local[,,,timestep] <- PathogenStageState
      pdet <- Pathogen$Engine$Resolve(Pathogen$DetectionProb, timestep, list(n_nodes=nrow(SDDprob),Ntimesteps=Ntimesteps), "DetectionProb")
      I_node <- apply(PathogenStageState[, , "I", drop=FALSE], 1, sum)
      PathogenDetectedNow <- rbinom(nrow(SDDprob),1,1-(1-pdet)^I_node)
      PathogenDetectedResults_local[,timestep] <- PathogenDetectedNow
      if(isTRUE(Pathogen$DetectionTriggersInfo)) {
        if(UseInfoPersistence == T) LastKnownPresence[PathogenDetectedNow == 1] <- timestep
        HaveInfo[HaveInfo == 0] <- PathogenDetectedNow[HaveInfo == 0]
      }
    }
    
    InfoBeforeSurveillance <- as.integer(HaveInfo != 0)
    InformationStateBeforeSurveillanceResults_local[,timestep] <- InfoBeforeSurveillance
    BackgroundDetectionProbabilityResults_local[,,timestep] <- matrix(NodeDetectionProb[,,timestep],nrow=n_nodes,ncol=Nstages)
    DetectionProbPerStage <- 1 - (1 - matrix(NodeDetectionProb[,,timestep],nrow=n_nodes,ncol=Nstages))^N
    ProbDetectNode <- 1 - apply(1 - DetectionProbPerStage, 1, prod)
    BackgroundDetection <- rbinom(n=nrow(SDDprob),size=1,prob=ProbDetectNode)
    InfoTriggeredDetection <- integer(nrow(SDDprob))
    if(UseInfoTriggeredSurveillance) {
      InfoTriggeredDetectionProbabilityResults_local[,,timestep] <- matrix(NodeInfoTriggeredDetectionProb[,,timestep],nrow=n_nodes,ncol=Nstages)
      InfoDetectionProbPerStage <- 1 - (1 - matrix(NodeInfoTriggeredDetectionProb[,,timestep],nrow=n_nodes,ncol=Nstages))^N
      ProbInfoDetectNode <- 1 - apply(1 - InfoDetectionProbPerStage,1,prod)
      InfoTriggeredDetection <- rbinom(n=nrow(SDDprob),size=1,prob=ProbInfoDetectNode*InfoBeforeSurveillance)
    }
    BackgroundDetectedResults_local[,timestep] <- BackgroundDetection
    InfoTriggeredDetectedResults_local[,timestep] <- InfoTriggeredDetection
    HostDetectionEvidence <- pmax(BackgroundDetection,InfoTriggeredDetection)
    if (UseInfoPersistence == TRUE) {
      KnownPresence <- which(HostDetectionEvidence == 1)
      if (length(KnownPresence) > 0) LastKnownPresence[KnownPresence] <- timestep
    }
    HaveInfo[HaveInfo == 0] <- HostDetectionEvidence[HaveInfo == 0]
    HaveInfoResults_local[,timestep] <- HaveInfo
    DetectedResults_local[, timestep] <- HaveInfo * Invaded
  } # end timesteps
  
  # ---------- END your block ----------
  
  # Return local results for this permutation
  list(
    PopulationResults = PopulationResults_local,
    PopulationStageResults = PopulationStageResults_local,
    InvasionResults = InvasionResults_local,
    DetectedResults = DetectedResults_local,
    BackgroundDetectedResults = BackgroundDetectedResults_local,
    InfoTriggeredDetectedResults = InfoTriggeredDetectedResults_local,
    BackgroundDetectionProbabilityResults = BackgroundDetectionProbabilityResults_local,
    InfoTriggeredDetectionProbabilityResults = InfoTriggeredDetectionProbabilityResults_local,
    InformationStateBeforeSurveillanceResults = InformationStateBeforeSurveillanceResults_local,
    HaveInfoResults = HaveInfoResults_local,
    ManagingResults = ManagingResults_local,
    PathogenStageResults = PathogenStageResults_local,
    PathogenDeathResults = PathogenDeathResults_local,
    NewInfectionResults = NewInfectionResults_local,
    PathogenDetected = PathogenDetectedResults_local
  )
 }

if (n_cores == 1L) {
  perm_results <- lapply(seq_len(Nperm), PermutationWorker)
} else {
  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(if (inherits(cl, "cluster")) parallel::stopCluster(cl), add = TRUE)
  if (is.null(Seed)) parallel::clusterSetRNGStream(cl) else parallel::clusterSetRNGStream(cl, iseed=as.integer(Seed))
  if(!is.null(Pathogen)) {
    helper_names <- c(".iptm_take", ".iptm_reconcile", "INApestPathogenStageState",
                      ".iptm_resolve", ".iptm_node_contact", ".iptm_pathogen_step")
    helper_env <- environment(sys.function())
    missing_helpers <- helper_names[!vapply(helper_names, function(nm)
      exists(nm, envir=helper_env, mode="function", inherits=TRUE), logical(1))]
    if(length(missing_helpers)) stop("Missing pathogen transition helper(s): ", paste(missing_helpers, collapse=", "), ". Source INApestPathogenTransitionMatrix.R first.")
    for(nm in helper_names) assign(nm, get(nm, envir=helper_env, mode="function", inherits=TRUE), envir=environment())
    parallel::clusterExport(cl, helper_names, envir=environment())
  }
  # materialise custom LocalDynamics lexical bindings on PSOCK workers.
  if(length(LocalDynamicsPSOCKBindings)) parallel::clusterCall(cl,function(bindings){list2env(bindings,envir=.GlobalEnv);NULL},LocalDynamicsPSOCKBindings)

  # ---------------------------------------------------------------------------
  # Run independent permutations across PSOCK workers.
  # ---------------------------------------------------------------------------
  perm_results <- parallel::parLapply(cl, seq_len(Nperm), PermutationWorker)
  parallel::stopCluster(cl)
  cl <- NULL
}

# combine the returned per-perm results into master arrays
for (i in seq_along(perm_results)) {
  PopulationResults[,, i]       <- perm_results[[i]]$PopulationResults
  PopulationStageResults[,,, i] <- perm_results[[i]]$PopulationStageResults
  InvasionResults[,, i]         <- perm_results[[i]]$InvasionResults
  DetectedResults[,, i]         <- perm_results[[i]]$DetectedResults
  BackgroundDetectedResults[,,i] <- perm_results[[i]]$BackgroundDetectedResults
  InfoTriggeredDetectedResults[,,i] <- perm_results[[i]]$InfoTriggeredDetectedResults
  InformationStateBeforeSurveillanceResults[,,i] <- perm_results[[i]]$InformationStateBeforeSurveillanceResults
  HaveInfoResults[,,i] <- perm_results[[i]]$HaveInfoResults
  BackgroundDetectionProbabilityResults[,,,i] <- perm_results[[i]]$BackgroundDetectionProbabilityResults
  InfoTriggeredDetectionProbabilityResults[,,,i] <- perm_results[[i]]$InfoTriggeredDetectionProbabilityResults
  ManagingResults[,, i]         <- perm_results[[i]]$ManagingResults
  if(!is.null(Pathogen)) {
    PathogenStageResults[,,,,i] <- perm_results[[i]]$PathogenStageResults
    PathogenDeathResults[,,,i] <- perm_results[[i]]$PathogenDeathResults
    NewInfectionResults[,,,i] <- perm_results[[i]]$NewInfectionResults
    PathogenDetectedResults[,,i] <- perm_results[[i]]$PathogenDetected
  }
}


###########################################################
### End of Simulation
###########################################################

###########################################################
### Save results for post-hoc analyses
###########################################################
### ModelName used to generate filenames
# Use standard format for ease of reading results to produce heat maps
# Support post-processing comparisons among management scenarios.
if(is.na(OutputDir) == T)
	OutputDir = ""
FileNameStem = paste0(OutputDir,ModelName)

# These are 3D arrays with dimensions (Nodes,Timesteps,Realisations)
saveRDS(ManagingResults, paste0(FileNameStem,"InfoLargeOut.rds"))
saveRDS(PopulationResults, paste0(FileNameStem,"PopulationLargeOut.rds"))
saveRDS(PopulationStageResults, paste0(FileNameStem,"PopulationStageLargeOut.rds"))
if(!is.null(Pathogen)) {
  saveRDS(PathogenStageResults, paste0(FileNameStem,"PathogenStageLargeOut.rds"))
  saveRDS(PathogenDeathResults, paste0(FileNameStem,"PathogenDeathLargeOut.rds"))
  saveRDS(NewInfectionResults, paste0(FileNameStem,"NewInfectionLargeOut.rds"))
  saveRDS(PathogenDetectedResults, paste0(FileNameStem,"PathogenDetectedLargeOut.rds"))
}
saveRDS(InvasionResults, paste0(FileNameStem,"InvasionLargeOut.rds"))
saveRDS(DetectedResults, paste0(FileNameStem,"DetectedLargeOut.rds"))
saveRDS(BackgroundDetectedResults, paste0(FileNameStem,"BackgroundDetectedLargeOut.rds"))
saveRDS(InfoTriggeredDetectedResults, paste0(FileNameStem,"InfoTriggeredDetectedLargeOut.rds"))
saveRDS(InformationStateBeforeSurveillanceResults, paste0(FileNameStem,"InformationStateBeforeSurveillanceLargeOut.rds"))
saveRDS(HaveInfoResults, paste0(FileNameStem,"HaveInfoLargeOut.rds"))
saveRDS(BackgroundDetectionProbabilityResults, paste0(FileNameStem,"BackgroundDetectionProbabilityLargeOut.rds"))
saveRDS(InfoTriggeredDetectionProbabilityResults, paste0(FileNameStem,"InfoTriggeredDetectionProbabilityLargeOut.rds"))

##########################################################
### Store node-level invasion probabilities for each timestep.
### and estimation of invasion threat to other regions
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
### Produce summary figs when processing completed
###########################################################

Title = ModelName


# Change in total population with time
# Plots of raw values for each realisation and summaries (median and 95% CI) provided

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

# Change in number of nodes infested with time
# Plots of raw values for each realisation and summaries (median and 95% CI) provided

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


# Summarise the number of nodes under management through time.
# Plots of raw values for each realisation and summaries (median and 95% CI) provided

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


# Summarise occupied nodes known to the response system through time.
# Plots of raw values for each realisation and summaries (median and 95% CI) provided

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


# Summarise the proportion of occupied nodes known through time.
# Plots of raw values for each realisation and summaries (median and 95% CI) provided
 
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
if(ReturnResults) {
  ResultObject <- list(ModelName=ModelName,PopulationResults=PopulationResults,PopulationStageResults=PopulationStageResults,
    InvasionResults=InvasionResults,ManagingResults=ManagingResults,DetectedResults=DetectedResults,
    BackgroundDetectedResults=BackgroundDetectedResults,InfoTriggeredDetectedResults=InfoTriggeredDetectedResults,
    BackgroundDetectionProbabilityResults=BackgroundDetectionProbabilityResults,InfoTriggeredDetectionProbabilityResults=InfoTriggeredDetectionProbabilityResults,
    InformationStateBeforeSurveillanceResults=InformationStateBeforeSurveillanceResults,HaveInfoResults=HaveInfoResults,
    InvasionProb=InvasionProb)
  if(!is.null(Pathogen)) { ResultObject$PathogenStageResults<-PathogenStageResults; ResultObject$PathogenDeathResults<-PathogenDeathResults; ResultObject$NewInfectionResults<-NewInfectionResults; ResultObject$PathogenDetectedResults<-PathogenDetectedResults }
  class(ResultObject)<-c("INApestMetaTransitionMatrixParallel","list")
  return(invisible(ResultObject))
}
invisible(NULL)
}


################################################################
################################################################
### End of function
################################################################
################################################################
