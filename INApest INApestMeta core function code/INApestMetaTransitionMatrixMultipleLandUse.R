###############################################################################
### INApest MLU x Transition Matrix definitive H+P+B engine
### Repository source-hygiene release v0.5.4 -- 2026-10-01
###
### This file is the exact sequential composition of the validated MLUTM
### development lineage v0.1.1 -> v0.2 -> v0.3 -> v0.4.2. No biological
### statements were reordered: each historical layer appears below in the same
### source order used by the frozen validation bundles.
###############################################################################


###############################################################################
### BEGIN HISTORICAL LAYER: INApestMetaTransitionMatrixMultipleLandUse_core_v0.1.1.R
###############################################################################

###############################################################################
### INApestMetaTransitionMatrixMultipleLandUse host biology core v0.1.1
### Date: 2026-10-01
###
### Purpose
### -------
### First additive implementation checkpoint for a combined
### node x land-use x demographic-stage INApest host state.
###
### This checkpoint is deliberately HOST ONLY and SERIAL ONLY.
### Pathogen, biocontrol, RK continuous local dynamics, analytical wrappers,
### proof of absence/freedom, observation-history conditioning and PSOCK are
### reserved for later gated checkpoints.
###
### Canonical parent source pins (GitHub main, verified 2026-10-01):
###   INApestMetaMultipleLandUse.r  blob 7f1e47949bf48c12c24ee32bcb3f0553b39fce18
###   INApestMetaTransitionMatrix.r blob 265a880dad3102c8e41675be9f298281ddee60fa
###
### Design principle
### ----------------
### For >=2 demographic stages, the existing definitive Transition Matrix local
### biology is reused on a flattened node x land-use cell axis. This carries
### forward the current TM contracts for stage transitions, transition movement,
### blocked-transition mortality, density-dependent dispersal, integer
### propagule routing and seedbank recruitment. Land use supplies cell-specific
### capacities, management modifiers and target-land-use recruitment weights.
###
### Exact boundary reductions are explicit:
###   * Nlanduses = 1 -> definitive Transition Matrix local biology.
###   * Nstages    = 1 -> definitive Multiple Land Use local biology.
###
### No parent source is modified by this file.
###############################################################################

.inatmlu_stop_missing <- function(name) {
  stop(name, " is required for this reduction/core path", call. = FALSE)
}

.inatmlu_require_parent_functions <- function() {
  if (!exists("local.dynamicsLU", mode = "function"))
    stop("Source INApestMetaMultipleLandUse.r before the MLU x TM core", call. = FALSE)
  if (!exists("local.dynamics.transition.matrix", mode = "function"))
    stop("Source INApestMetaTransitionMatrix.r before the MLU x TM core", call. = FALSE)
  invisible(TRUE)
}

.inatmlu_validate_array_state <- function(n0) {
  d <- dim(n0)
  if (length(d) != 3L)
    stop("n0 must be a nodes x land-uses x stages array", call. = FALSE)
  if (any(d < 1L))
    stop("Every n0 dimension must be at least one", call. = FALSE)
  if (any(!is.finite(n0)) || any(n0 < 0) || any(abs(n0 - round(n0)) > 1e-8))
    stop("n0 must contain finite non-negative integer-valued host counts", call. = FALSE)
  as.integer(d)
}

.inatmlu_surface <- function(x, n_nodes, n_landuses, name, allow_vector = TRUE) {
  if (is.null(x)) stop(name, " cannot be NULL", call. = FALSE)
  d <- dim(x)
  if (!is.null(d)) {
    if (length(d) != 2L || !identical(as.integer(d), c(n_nodes, n_landuses)))
      stop(name, " matrix must have dimensions nodes x land-uses", call. = FALSE)
    out <- matrix(as.numeric(x), nrow = n_nodes, ncol = n_landuses)
  } else {
    x <- as.numeric(x)
    if (length(x) == 1L) {
      out <- matrix(x, nrow = n_nodes, ncol = n_landuses)
    } else if (isTRUE(allow_vector) && length(x) == n_nodes && length(x) != n_landuses) {
      out <- matrix(rep(x, n_landuses), nrow = n_nodes, ncol = n_landuses)
    } else if (isTRUE(allow_vector) && length(x) == n_landuses && length(x) != n_nodes) {
      out <- matrix(rep(x, each = n_nodes), nrow = n_nodes, ncol = n_landuses)
    } else if (isTRUE(allow_vector) && length(x) == n_nodes && n_nodes == n_landuses) {
      stop(name, " vector is ambiguous because nodes equals land-uses; supply a nodes x land-uses matrix", call. = FALSE)
    } else {
      stop(name, " must be scalar, length nodes, length land-uses, or nodes x land-uses", call. = FALSE)
    }
  }
  if (any(!is.finite(out))) stop(name, " must be finite", call. = FALSE)
  out
}

.inatmlu_probability_surface <- function(x, n_nodes, n_landuses, name) {
  out <- .inatmlu_surface(x, n_nodes, n_landuses, name)
  if (any(out < 0 | out > 1)) stop(name, " values must be between 0 and 1", call. = FALSE)
  out
}

.inatmlu_flat_index <- function(node, landuse, n_landuses) {
  (node - 1L) * n_landuses + landuse
}

.inatmlu_flatten_state <- function(x) {
  d <- dim(x)
  n_nodes <- d[1L]; n_landuses <- d[2L]; n_stages <- d[3L]
  out <- matrix(0, nrow = n_nodes * n_landuses, ncol = n_stages)
  for (i in seq_len(n_nodes)) {
    for (l in seq_len(n_landuses)) {
      out[.inatmlu_flat_index(i, l, n_landuses), ] <- x[i, l, ]
    }
  }
  out
}

.inatmlu_unflatten_state <- function(x, n_nodes, n_landuses, n_stages) {
  if (!is.matrix(x) || !identical(dim(x), c(n_nodes * n_landuses, n_stages)))
    stop("Internal error: flattened state has unexpected dimensions", call. = FALSE)
  out <- array(0, dim = c(n_nodes, n_landuses, n_stages),
               dimnames = list(node = seq_len(n_nodes),
                               landuse = seq_len(n_landuses),
                               stage = seq_len(n_stages)))
  for (i in seq_len(n_nodes)) {
    for (l in seq_len(n_landuses)) {
      out[i, l, ] <- x[.inatmlu_flat_index(i, l, n_landuses), ]
    }
  }
  out
}

.inatmlu_flatten_surface <- function(x) {
  n_nodes <- nrow(x); n_landuses <- ncol(x)
  out <- numeric(n_nodes * n_landuses)
  for (i in seq_len(n_nodes)) {
    for (l in seq_len(n_landuses)) {
      out[.inatmlu_flat_index(i, l, n_landuses)] <- x[i, l]
    }
  }
  out
}

.inatmlu_weights <- function(weights, n_nodes, n_landuses, n_stages) {
  n_cells <- n_nodes * n_landuses
  if (is.null(weights)) return(rep(1, n_stages))
  d <- dim(weights)
  if (is.null(d)) {
    w <- as.numeric(weights)
    if (length(w) != n_stages)
      stop("weights vector must have length stages", call. = FALSE)
    if (any(!is.finite(w)) || any(w <= 0))
      stop("weights must be finite and > 0", call. = FALSE)
    return(w)
  }
  if (length(d) == 2L && identical(as.integer(d), c(n_nodes, n_stages))) {
    # Node-specific weights are repeated across land-use classes.
    out <- matrix(0, nrow = n_cells, ncol = n_stages)
    for (i in seq_len(n_nodes)) for (l in seq_len(n_landuses))
      out[.inatmlu_flat_index(i, l, n_landuses), ] <- weights[i, ]
  } else if (length(d) == 3L && identical(as.integer(d), c(n_nodes, n_landuses, n_stages))) {
    out <- .inatmlu_flatten_state(weights)
  } else {
    stop("weights must be length stages, nodes x stages, or nodes x land-uses x stages", call. = FALSE)
  }
  if (any(!is.finite(out)) || any(out <= 0)) stop("weights must be finite and > 0", call. = FALSE)
  out
}

.inatmlu_fecundity_reduction <- function(x, n_nodes, n_landuses, n_stages) {
  n_cells <- n_nodes * n_landuses
  d <- dim(x)
  if (is.null(d)) {
    z <- as.numeric(x)
    if (length(z) == 1L) return(z)
    if (length(z) == n_stages && !(n_stages %in% c(n_nodes, n_landuses))) {
      out <- matrix(rep(z, each = n_cells), nrow = n_cells, ncol = n_stages)
      return(out)
    }
    # Other vectors are intentionally rejected to avoid silent ambiguity.
    stop("nodefecundityreduction must be scalar, nodes x land-uses, or nodes x land-uses x stages", call. = FALSE)
  }
  if (length(d) == 2L && identical(as.integer(d), c(n_nodes, n_landuses))) {
    surf <- .inatmlu_probability_surface(x, n_nodes, n_landuses, "nodefecundityreduction")
    v <- .inatmlu_flatten_surface(surf)
    return(matrix(rep(v, n_stages), nrow = n_cells, ncol = n_stages))
  }
  if (length(d) == 3L && identical(as.integer(d), c(n_nodes, n_landuses, n_stages))) {
    out <- .inatmlu_flatten_state(x)
    if (any(!is.finite(out)) || any(out < 0 | out > 1))
      stop("nodefecundityreduction values must be between 0 and 1", call. = FALSE)
    return(out)
  }
  stop("nodefecundityreduction must be scalar, nodes x land-uses, or nodes x land-uses x stages", call. = FALSE)
}

.inatmlu_transition <- function(x, n_nodes, n_landuses, n_stages) {
  n_cells <- n_nodes * n_landuses
  if (is.matrix(x)) {
    if (!identical(dim(x), c(n_stages, n_stages)))
      stop("nodetransition matrix must be stages x stages", call. = FALSE)
    if (any(!is.finite(x)) || any(x < 0))
      stop("nodetransition entries must be finite and non-negative", call. = FALSE)
    return(x)
  }
  if (is.list(x)) {
    if (length(x) != n_cells)
      stop("nodetransition list must contain one stages x stages matrix per node x land-use cell", call. = FALSE)
    for (k in seq_len(n_cells)) {
      if (!is.matrix(x[[k]]) || !identical(dim(x[[k]]), c(n_stages, n_stages)) ||
          any(!is.finite(x[[k]])) || any(x[[k]] < 0))
        stop("Every nodetransition list entry must be a finite non-negative stages x stages matrix", call. = FALSE)
    }
    return(x)
  }
  d <- dim(x)
  if (length(d) == 3L && identical(as.integer(d), c(n_stages, n_stages, n_landuses))) {
    out <- vector("list", n_cells)
    for (i in seq_len(n_nodes)) for (l in seq_len(n_landuses))
      out[[.inatmlu_flat_index(i, l, n_landuses)]] <- matrix(x[, , l], n_stages, n_stages)
    return(out)
  }
  if (length(d) == 4L && identical(as.integer(d), c(n_stages, n_stages, n_nodes, n_landuses))) {
    out <- vector("list", n_cells)
    for (i in seq_len(n_nodes)) for (l in seq_len(n_landuses))
      out[[.inatmlu_flat_index(i, l, n_landuses)]] <- matrix(x[, , i, l], n_stages, n_stages)
    return(out)
  }
  stop("nodetransition must be a stages x stages matrix, a cell list, stages x stages x land-uses, or stages x stages x nodes x land-uses", call. = FALSE)
}

.inatmlu_validate_node_connectivity <- function(P, n_nodes, name, allow_disabled = FALSE) {
  if (allow_disabled && length(P) == 1L && (is.na(P) || identical(as.numeric(P), 0))) return(invisible(TRUE))
  if (!is.matrix(P) || !identical(dim(P), c(n_nodes, n_nodes)))
    stop(name, " must be a nodes x nodes matrix", call. = FALSE)
  if (any(!is.finite(P)) || any(P < 0) || any(rowSums(P) > 1 + 1e-10))
    stop(name, " must contain finite non-negative probabilities with source-row sums <= 1", call. = FALSE)
  invisible(TRUE)
}

.inatmlu_recruitment_weights <- function(x, seedbankK, n_nodes, n_landuses) {
  if (is.null(x)) {
    w <- seedbankK
    w[w < 0] <- 0
  } else {
    w <- .inatmlu_surface(x, n_nodes, n_landuses, "LandUseRecruitmentWeights")
    if (any(w < 0)) stop("LandUseRecruitmentWeights must be non-negative", call. = FALSE)
  }
  rs <- rowSums(w)
  out <- matrix(0, nrow = n_nodes, ncol = n_landuses)
  pos <- rs > 0
  if (any(pos)) out[pos, ] <- w[pos, , drop = FALSE] / rs[pos]
  out
}

# Expand node-to-node reproductive dispersal to node x land-use cells.
# Every source land-use uses its source node's spatial row. At a target node,
# LandUseRecruitmentWeights allocate represented arrival probability among
# land-use classes. A zero target row means propagules routed to that node have
# no represented receiving land-use class and therefore become residual loss.
.inatmlu_expand_reproductive_dispersal <- function(Pnode, target_weights) {
  n_nodes <- nrow(Pnode); n_landuses <- ncol(target_weights)
  n_cells <- n_nodes * n_landuses
  out <- matrix(0, nrow = n_cells, ncol = n_cells)
  for (i in seq_len(n_nodes)) {
    for (ls in seq_len(n_landuses)) {
      src <- .inatmlu_flat_index(i, ls, n_landuses)
      for (j in seq_len(n_nodes)) {
        if (Pnode[i, j] <= 0) next
        for (lt in seq_len(n_landuses)) {
          if (target_weights[j, lt] <= 0) next
          dst <- .inatmlu_flat_index(j, lt, n_landuses)
          out[src, dst] <- Pnode[i, j] * target_weights[j, lt]
        }
      }
    }
  }
  out
}

.inatmlu_landuse_mixing <- function(x, n_landuses) {
  if (is.null(x)) return(diag(n_landuses))
  if (!is.matrix(x) || !identical(dim(x), c(n_landuses, n_landuses)))
    stop("TransitionLandUseMixing must be land-uses x land-uses", call. = FALSE)
  if (any(!is.finite(x)) || any(x < 0) || any(rowSums(x) > 1 + 1e-10))
    stop("TransitionLandUseMixing must be finite, non-negative, with row sums <= 1", call. = FALSE)
  x
}

# Expand stage-transition movement. By default TransitionLandUseMixing is the
# identity, so moving individuals retain their land-use class while changing
# node. A supplied mixing matrix can allow simultaneous land-use change.
.inatmlu_expand_transition_movement_one <- function(Pnode, LU) {
  n_nodes <- nrow(Pnode); n_landuses <- nrow(LU)
  n_cells <- n_nodes * n_landuses
  out <- matrix(0, nrow = n_cells, ncol = n_cells)
  for (i in seq_len(n_nodes)) for (ls in seq_len(n_landuses)) {
    src <- .inatmlu_flat_index(i, ls, n_landuses)
    for (j in seq_len(n_nodes)) {
      if (Pnode[i, j] <= 0) next
      for (lt in seq_len(n_landuses)) {
        if (LU[ls, lt] <= 0) next
        dst <- .inatmlu_flat_index(j, lt, n_landuses)
        out[src, dst] <- Pnode[i, j] * LU[ls, lt]
      }
    }
  }
  out
}

.inatmlu_expand_transition_movement <- function(x, LU, n_nodes, n_stages, name) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) {
    if (length(x) != n_stages - 1L)
      stop(name, " list must have length stages - 1", call. = FALSE)
    return(lapply(x, function(P) {
      if (is.null(P)) return(NULL)
      .inatmlu_validate_node_connectivity(P, n_nodes, name)
      .inatmlu_expand_transition_movement_one(P, LU)
    }))
  }
  .inatmlu_validate_node_connectivity(x, n_nodes, name)
  .inatmlu_expand_transition_movement_one(x, LU)
}

.inatmlu_blocked_mortality <- function(x, n_nodes, n_landuses, n_stages) {
  if (length(x) == 1L || (is.null(dim(x)) && length(x) == n_stages - 1L)) return(x)
  d <- dim(x)
  if (length(d) == 3L && identical(as.integer(d), c(n_nodes, n_landuses, n_stages - 1L))) {
    n_cells <- n_nodes * n_landuses
    out <- matrix(0, nrow = n_cells, ncol = n_stages - 1L)
    for (s in seq_len(n_stages - 1L)) out[, s] <- .inatmlu_flatten_surface(x[, , s])
    return(out)
  }
  stop("BlockedTransitionMortality must be scalar, length stages - 1, or nodes x land-uses x (stages - 1)", call. = FALSE)
}

#' One-step host biology for the MLU x Transition Matrix architecture.
#'
#' State is nodes x land-uses x demographic stages. This first checkpoint does
#' not implement surveillance/information state and does not call pathogen,
#' biocontrol or RK components.
INApestMLUTMHostStep <- function(
    n0,
    nodetransition,
    weights = NULL,
    sddprob,
    nodeenvestabprob = 1,
    lddprob = NA,
    lddrate = 0,
    nodeK,
    node.seedbankK = nodeK,
    nodepropaguleestablishment = 1,
    nodespreadreduction = 0,
    nodefecundityreduction = 0,
    managing = 0,
    transition_sddprob = NULL,
    transition_lddprob = NULL,
    transition_lddrate = 0,
    ApplyFootprintToTransitions = FALSE,
    BlockedTransitionMortality = 0,
    DispersalDensityFactor = 0,
    MaxInteger = .Machine$integer.max,
    LandUseRecruitmentWeights = NULL,
    TransitionLandUseMixing = NULL,
    MLUPropaguleProduction = NULL) {

  .inatmlu_require_parent_functions()
  dims <- .inatmlu_validate_array_state(n0)
  n_nodes <- dims[1L]; n_landuses <- dims[2L]; n_stages <- dims[3L]

  .inatmlu_validate_node_connectivity(sddprob, n_nodes, "sddprob")
  .inatmlu_validate_node_connectivity(lddprob, n_nodes, "lddprob", allow_disabled = TRUE)

  # -------------------------------------------------------------------------
  # Exact MLU reduction: one demographic stage.
  # -------------------------------------------------------------------------
  if (n_stages == 1L) {
    if (is.null(MLUPropaguleProduction)) .inatmlu_stop_missing("MLUPropaguleProduction")
    K_lu <- .inatmlu_surface(nodeK, n_nodes, n_landuses, "nodeK")
    manage_lu <- .inatmlu_probability_surface(managing, n_nodes, n_landuses, "managing")
    spread_lu <- .inatmlu_probability_surface(nodespreadreduction, n_nodes, n_landuses, "nodespreadreduction")
    fec_lu <- .inatmlu_probability_surface(nodefecundityreduction, n_nodes, n_landuses, "nodefecundityreduction")
    n_lu <- matrix(n0[, , 1L], nrow = n_nodes, ncol = n_landuses)
    out <- local.dynamicsLU(
      sddprob = sddprob,
      nodepropaguleproduction = MLUPropaguleProduction,
      nodeenvestabprob = nodeenvestabprob,
      n = n_lu,
      lddprob = lddprob,
      lddrate = lddrate,
      k_is_0 = rowSums(K_lu) <= 0,
      nodeK = K_lu,
      nodepropaguleestablishment = nodepropaguleestablishment,
      nodespreadreduction = spread_lu,
      nodefecundityreduction = fec_lu,
      managing = manage_lu
    )
    return(array(out, dim = c(n_nodes, n_landuses, 1L),
                 dimnames = list(node = seq_len(n_nodes), landuse = seq_len(n_landuses), stage = 1L)))
  }

  # -------------------------------------------------------------------------
  # Exact Transition Matrix reduction: one land-use class.
  # -------------------------------------------------------------------------
  if (n_landuses == 1L) {
    n_tm <- matrix(n0[, 1L, ], nrow = n_nodes, ncol = n_stages)
    K_tm <- as.numeric(.inatmlu_surface(nodeK, n_nodes, 1L, "nodeK")[, 1L])
    seed_tm <- as.numeric(.inatmlu_surface(node.seedbankK, n_nodes, 1L, "node.seedbankK")[, 1L])
    env_tm <- as.numeric(.inatmlu_surface(nodeenvestabprob, n_nodes, 1L, "nodeenvestabprob")[, 1L])
    estab_tm <- as.numeric(.inatmlu_surface(nodepropaguleestablishment, n_nodes, 1L, "nodepropaguleestablishment")[, 1L])
    spread_tm <- as.numeric(.inatmlu_probability_surface(nodespreadreduction, n_nodes, 1L, "nodespreadreduction")[, 1L])
    manage_tm <- as.numeric(.inatmlu_probability_surface(managing, n_nodes, 1L, "managing")[, 1L])

    # Preserve scalar fecundity-reduction input exactly where possible because
    # the parent TM function treats scalar input as its canonical legacy path.
    fec_tm <- nodefecundityreduction
    if (!is.null(dim(nodefecundityreduction))) {
      dfr <- dim(nodefecundityreduction)
      if (length(dfr) == 2L && identical(as.integer(dfr), c(n_nodes, 1L)))
        fec_tm <- as.numeric(nodefecundityreduction[, 1L])
      else if (length(dfr) == 3L && identical(as.integer(dfr), c(n_nodes, 1L, n_stages)))
        fec_tm <- matrix(nodefecundityreduction[, 1L, ], nrow = n_nodes, ncol = n_stages)
    }

    out <- local.dynamics.transition.matrix(
      nodetransition = nodetransition,
      weights = weights,
      sddprob = sddprob,
      nodeenvestabprob = env_tm,
      n0 = n_tm,
      lddprob = lddprob,
      lddrate = lddrate,
      transition_sddprob = transition_sddprob,
      transition_lddprob = transition_lddprob,
      transition_lddrate = transition_lddrate,
      nodeK = K_tm,
      node.seedbankK = seed_tm,
      nodepropaguleestablishment = estab_tm,
      nodespreadreduction = spread_tm,
      nodefecundityreduction = fec_tm,
      managing = manage_tm,
      MaxInteger = MaxInteger,
      ApplyFootprintToTransitions = ApplyFootprintToTransitions,
      BlockedTransitionMortality = BlockedTransitionMortality,
      DispersalDensityFactor = DispersalDensityFactor
    )
    return(array(out, dim = c(n_nodes, 1L, n_stages),
                 dimnames = list(node = seq_len(n_nodes), landuse = 1L, stage = seq_len(n_stages))))
  }

  # -------------------------------------------------------------------------
  # Combined node x land-use x stage path.
  # -------------------------------------------------------------------------
  n_cells <- n_nodes * n_landuses
  n_flat <- .inatmlu_flatten_state(n0)
  K_lu <- .inatmlu_surface(nodeK, n_nodes, n_landuses, "nodeK")
  seed_lu <- .inatmlu_surface(node.seedbankK, n_nodes, n_landuses, "node.seedbankK")
  if (any(K_lu < 0) || any(seed_lu < 0)) stop("nodeK and node.seedbankK must be non-negative", call. = FALSE)

  env_lu <- .inatmlu_probability_surface(nodeenvestabprob, n_nodes, n_landuses, "nodeenvestabprob")
  estab_lu <- .inatmlu_probability_surface(nodepropaguleestablishment, n_nodes, n_landuses, "nodepropaguleestablishment")
  spread_lu <- .inatmlu_probability_surface(nodespreadreduction, n_nodes, n_landuses, "nodespreadreduction")
  manage_lu <- .inatmlu_probability_surface(managing, n_nodes, n_landuses, "managing")
  fec_red <- .inatmlu_fecundity_reduction(nodefecundityreduction, n_nodes, n_landuses, n_stages)
  w_flat <- .inatmlu_weights(weights, n_nodes, n_landuses, n_stages)
  A_flat <- .inatmlu_transition(nodetransition, n_nodes, n_landuses, n_stages)

  target_w <- .inatmlu_recruitment_weights(LandUseRecruitmentWeights, seed_lu, n_nodes, n_landuses)
  SDD_flat <- .inatmlu_expand_reproductive_dispersal(sddprob, target_w)
  LDD_flat <- if (length(lddprob) == 1L && (is.na(lddprob) || identical(as.numeric(lddprob), 0))) {
    NA
  } else {
    .inatmlu_expand_reproductive_dispersal(lddprob, target_w)
  }

  LU_mix <- .inatmlu_landuse_mixing(TransitionLandUseMixing, n_landuses)
  TSDD_flat <- .inatmlu_expand_transition_movement(transition_sddprob, LU_mix, n_nodes, n_stages, "transition_sddprob")
  TLDD_flat <- .inatmlu_expand_transition_movement(transition_lddprob, LU_mix, n_nodes, n_stages, "transition_lddprob")
  blocked_flat <- .inatmlu_blocked_mortality(BlockedTransitionMortality, n_nodes, n_landuses, n_stages)

  out_flat <- local.dynamics.transition.matrix(
    nodetransition = A_flat,
    weights = w_flat,
    sddprob = SDD_flat,
    nodeenvestabprob = .inatmlu_flatten_surface(env_lu),
    n0 = n_flat,
    lddprob = LDD_flat,
    lddrate = lddrate,
    transition_sddprob = TSDD_flat,
    transition_lddprob = TLDD_flat,
    transition_lddrate = transition_lddrate,
    nodeK = .inatmlu_flatten_surface(K_lu),
    node.seedbankK = .inatmlu_flatten_surface(seed_lu),
    nodepropaguleestablishment = .inatmlu_flatten_surface(estab_lu),
    nodespreadreduction = .inatmlu_flatten_surface(spread_lu),
    nodefecundityreduction = fec_red,
    managing = .inatmlu_flatten_surface(manage_lu),
    MaxInteger = MaxInteger,
    ApplyFootprintToTransitions = ApplyFootprintToTransitions,
    BlockedTransitionMortality = blocked_flat,
    DispersalDensityFactor = DispersalDensityFactor
  )

  out <- .inatmlu_unflatten_state(out_flat, n_nodes, n_landuses, n_stages)
  if (any(!is.finite(out)) || any(out < 0) || any(abs(out - round(out)) > 1e-8))
    stop("Internal error: combined host step returned invalid counts", call. = FALSE)
  out
}

# Descriptive alias matching the eventual architecture name.
local.dynamics.transition.matrix.multiple.land.use <- INApestMLUTMHostStep


### END HISTORICAL LAYER: INApestMetaTransitionMatrixMultipleLandUse_core_v0.1.1.R


###############################################################################
### BEGIN HISTORICAL LAYER: INApestMetaTransitionMatrixMultipleLandUse_response_v0.2.R
###############################################################################

###############################################################################
### INApest MLU x Transition Matrix serial response layer v0.2
### Date: 2026-10-01
###
### Parent checkpoint
### -----------------
### INApest_MLUTM_HostCore_v0.1.1_FROZEN_20261001 (native R 4.4.1: 10/10 PASS)
###
### Scope of this checkpoint
### ------------------------
### Adds a SERIAL one-timestep response/observation wrapper around the frozen
### node x land-use x demographic-stage host biology core. This layer owns:
###   * node x land-use management adoption conditional on node information;
###   * stage-specific management mortality;
###   * programmed information persistence and probabilistic information decay;
###   * SEAM transfer from currently known occupied nodes;
###   * background host surveillance;
###   * information-triggered host surveillance; and
###   * end-of-step information updating.
###
### Deliberately NOT active here:
###   pathogen, biocontrol, RK continuous local dynamics, external invasion,
###   external information, parameter SD draws, analytical, PoA/PoF,
###   ObservationHistory and PSOCK.
###
### Event order follows the current MLU/TM parent engines:
###   1. current information activates management;
###   2. management mortality is applied;
###   3. frozen v0.1.1 host biology advances;
###   4. programmed information stopping / retention are applied;
###   5. SEAM transfers existing response information;
###   6. surveillance observes the post-biology host state;
###   7. detections refresh/add information for the NEXT response step.
###
### This means a detection at the end of timestep t cannot trigger management
### in timestep t. It can trigger management in timestep t+1.
###############################################################################

.inatmlu_response_require <- function() {
  if (!exists("INApestMLUTMHostStep", mode = "function"))
    stop("Source the frozen MLUTM host core before the response layer", call. = FALSE)
  invisible(TRUE)
}

.inatmlu_binary_info <- function(x, n_nodes, name) {
  x <- as.numeric(x)
  if (length(x) == 1L) x <- rep(x, n_nodes)
  if (length(x) != n_nodes || any(!is.finite(x)) || any(!x %in% c(0, 1)))
    stop(name, " must be binary scalar or vector of length nodes", call. = FALSE)
  as.integer(x)
}

.inatmlu_node_prob <- function(x, n_nodes, name) {
  x <- as.numeric(x)
  if (length(x) == 1L) x <- rep(x, n_nodes)
  if (length(x) != n_nodes || any(!is.finite(x)) || any(x < 0 | x > 1))
    stop(name, " must be a probability scalar or vector of length nodes", call. = FALSE)
  x
}

.inatmlu_persistence <- function(x, n_nodes) {
  x <- as.numeric(x)
  if (length(x) == 1L) x <- rep(x, n_nodes)
  if (length(x) != n_nodes)
    stop("info_persistence_steps must be scalar or vector of length nodes", call. = FALSE)
  bad <- !is.na(x) & (!is.finite(x) | x < 0 | x != floor(x))
  if (any(bad))
    stop("info_persistence_steps values must be non-negative whole numbers or NA", call. = FALSE)
  x
}

.inatmlu_response_surface <- function(x, n_nodes, n_landuses, name) {
  d <- dim(x)
  if (is.null(d)) {
    z <- as.numeric(x)
    if (length(z) == 1L) out <- matrix(z, n_nodes, n_landuses)
    else if (length(z) == n_landuses && n_landuses != n_nodes)
      out <- matrix(rep(z, each = n_nodes), n_nodes, n_landuses)
    else if (length(z) == n_nodes && n_nodes != n_landuses)
      out <- matrix(rep(z, n_landuses), n_nodes, n_landuses)
    else if (length(z) == n_nodes && n_nodes == n_landuses)
      stop(name, " vector is ambiguous because nodes equals land uses; supply a matrix", call. = FALSE)
    else stop(name, " must be scalar, length nodes, length land uses, or nodes x land uses", call. = FALSE)
  } else {
    if (length(d) != 2L || !identical(as.integer(d), c(n_nodes, n_landuses)))
      stop(name, " matrix must be nodes x land uses", call. = FALSE)
    out <- matrix(as.numeric(x), n_nodes, n_landuses)
  }
  if (any(!is.finite(out)) || any(out < 0 | out > 1))
    stop(name, " values must be finite probabilities in [0,1]", call. = FALSE)
  out
}

.inatmlu_response_cube <- function(x, n_nodes, n_landuses, n_stages, name) {
  d <- dim(x)
  if (is.null(d)) {
    z <- as.numeric(x)
    if (length(z) == 1L) {
      out <- array(z, c(n_nodes, n_landuses, n_stages))
    } else {
      matches <- c(nodes = length(z) == n_nodes,
                   landuses = length(z) == n_landuses,
                   stages = length(z) == n_stages)
      if (sum(matches) != 1L)
        stop(name, " vector is unsupported or ambiguous; use an explicit matrix/array", call. = FALSE)
      if (matches["nodes"])
        out <- array(rep(z, times = n_landuses * n_stages), c(n_nodes, n_landuses, n_stages))
      if (matches["landuses"])
        out <- array(rep(rep(z, each = n_nodes), times = n_stages), c(n_nodes, n_landuses, n_stages))
      if (matches["stages"])
        out <- array(rep(z, each = n_nodes * n_landuses), c(n_nodes, n_landuses, n_stages))
    }
  } else if (length(d) == 2L) {
    if (identical(as.integer(d), c(n_nodes, n_landuses))) {
      surf <- matrix(as.numeric(x), n_nodes, n_landuses)
      out <- array(rep(as.numeric(surf), times = n_stages), c(n_nodes, n_landuses, n_stages))
    } else if (identical(as.integer(d), c(n_nodes, n_stages))) {
      m <- matrix(as.numeric(x), n_nodes, n_stages)
      out <- array(0, c(n_nodes, n_landuses, n_stages))
      for (l in seq_len(n_landuses)) out[, l, ] <- m
    } else {
      stop(name, " matrix must be nodes x land uses or nodes x stages", call. = FALSE)
    }
  } else if (length(d) == 3L && identical(as.integer(d), c(n_nodes, n_landuses, n_stages))) {
    out <- array(as.numeric(x), c(n_nodes, n_landuses, n_stages))
  } else {
    stop(name, " must be scalar, an unambiguous vector, nodes x land uses, nodes x stages, or nodes x land uses x stages", call. = FALSE)
  }
  if (any(!is.finite(out)) || any(out < 0 | out > 1))
    stop(name, " values must be finite probabilities in [0,1]", call. = FALSE)
  out
}

INApestMLUTMHostDetectionProbability <- function(n, detection_prob) {
  dims <- .inatmlu_validate_array_state(n)
  n_nodes <- dims[1L]; n_landuses <- dims[2L]; n_stages <- dims[3L]
  p <- .inatmlu_response_cube(detection_prob, n_nodes, n_landuses, n_stages, "detection_prob")
  p_not <- (1 - p)^n
  out <- numeric(n_nodes)
  for (i in seq_len(n_nodes)) out[i] <- 1 - prod(p_not[i, , ])
  pmin(1, pmax(0, out))
}

INApestMLUTMResponseStep <- function(
    n,
    have_info,
    timestep = 1L,
    last_known_presence = NULL,
    manage_prob = 0,
    mortality_prob = 0,
    detection_prob = 0,
    info_triggered_detection_prob = 0,
    info_retention_prob = 1,
    info_persistence_steps = NA,
    seam = 0,
    biology_args = list()) {

  .inatmlu_response_require()
  dims <- .inatmlu_validate_array_state(n)
  n_nodes <- dims[1L]; n_landuses <- dims[2L]; n_stages <- dims[3L]
  if (length(timestep) != 1L || !is.finite(timestep) || timestep < 1 || timestep != floor(timestep))
    stop("timestep must be a positive integer", call. = FALSE)
  have_info <- .inatmlu_binary_info(have_info, n_nodes, "have_info")
  if (is.null(last_known_presence)) last_known_presence <- rep(NA_real_, n_nodes)
  last_known_presence <- as.numeric(last_known_presence)
  if (length(last_known_presence) != n_nodes || any(!is.na(last_known_presence) & !is.finite(last_known_presence)))
    stop("last_known_presence must be NULL or length nodes", call. = FALSE)
  retention <- .inatmlu_node_prob(info_retention_prob, n_nodes, "info_retention_prob")
  persistence <- .inatmlu_persistence(info_persistence_steps, n_nodes)
  use_persistence <- any(!is.na(persistence))

  if (!is.list(biology_args)) stop("biology_args must be a named list", call. = FALSE)
  if (length(biology_args) && (is.null(names(biology_args)) || any(!nzchar(names(biology_args)))))
    stop("biology_args must be a named list", call. = FALSE)
  if (anyDuplicated(names(biology_args))) stop("biology_args names must be unique", call. = FALSE)
  reserved <- intersect(names(biology_args), c("n0", "managing"))
  if (length(reserved)) stop("biology_args may not override n0 or managing", call. = FALSE)

  manage_surface <- .inatmlu_response_surface(manage_prob, n_nodes, n_landuses, "manage_prob")
  mortality_cube <- .inatmlu_response_cube(mortality_prob, n_nodes, n_landuses, n_stages, "mortality_prob")
  detection_cube <- .inatmlu_response_cube(detection_prob, n_nodes, n_landuses, n_stages, "detection_prob")
  info_detection_cube <- .inatmlu_response_cube(info_triggered_detection_prob, n_nodes, n_landuses, n_stages, "info_triggered_detection_prob")

  node_invaded_start <- as.integer(vapply(seq_len(n_nodes), function(i) sum(n[i, , ]) > 0, logical(1)))
  known_occupied_start <- node_invaded_start * have_info

  # Management adoption is node x land use, gated by node-level information.
  manage_p <- as.numeric(manage_surface) * rep(have_info, times = n_landuses)
  managing <- matrix(rbinom(n_nodes * n_landuses, size = 1L, prob = manage_p),
                     nrow = n_nodes, ncol = n_landuses)

  # Stage-specific management mortality. Flattening order is node -> land use ->
  # stage, so the one-land-use boundary consumes RNG in the same stage-major
  # order as the definitive Transition Matrix management loop.
  managing_cube <- array(rep(as.numeric(managing), times = n_stages),
                         dim = c(n_nodes, n_landuses, n_stages))
  survive_p <- 1 - mortality_cube * managing_cube
  n_after_management <- array(
    rbinom(length(n), size = as.integer(n), prob = as.numeric(survive_p)),
    dim = dim(n), dimnames = dimnames(n)
  )
  management_deaths <- n - n_after_management

  if (use_persistence) {
    killed_by_node <- vapply(seq_len(n_nodes), function(i) sum(management_deaths[i, , ]) > 0, logical(1))
    last_known_presence[killed_by_node] <- timestep
  }

  # Frozen biological core. Current management state is also supplied to the
  # biology step because it controls fecundity/spread reductions there.
  call_args <- c(list(n0 = n_after_management, managing = managing), biology_args)
  n_after_biology <- do.call(INApestMLUTMHostStep, call_args)

  # Programmed information stopping has priority where supplied.
  programmed <- which(have_info == 1L & !is.na(persistence))
  if (length(programmed)) {
    elapsed <- timestep - last_known_presence
    stop_nodes <- programmed[is.na(last_known_presence[programmed]) |
                             elapsed[programmed] >= persistence[programmed]]
    if (length(stop_nodes)) have_info[stop_nodes] <- 0L
  }

  # Probabilistic decay applies only to nodes without a programmed clock.
  decay <- which(have_info == 1L & is.na(persistence) & retention < 1)
  if (length(decay))
    have_info[decay] <- rbinom(length(decay), size = 1L, prob = retention[decay])

  # SEAM transfer uses the start-of-step known-occupied state, matching the
  # parent engines' stored Detected = Invaded * HaveInfo event ordering.
  seam_transferred <- integer(n_nodes)
  if (!(length(seam) == 1L && identical(as.numeric(seam), 0))) {
    if (!is.matrix(seam) || !identical(dim(seam), c(n_nodes, n_nodes)))
      stop("seam must be 0 or a nodes x nodes matrix", call. = FALSE)
    if (any(!is.finite(seam)) || any(seam < 0 | seam > 1))
      stop("seam entries must be probabilities in [0,1]", call. = FALSE)
    S <- seam
    diag(S) <- 0
    p_transfer <- S * known_occupied_start
    draws <- matrix(rbinom(n_nodes * n_nodes, 1L, prob = as.numeric(p_transfer)),
                    nrow = n_nodes, ncol = n_nodes)
    seam_transferred <- as.integer(colSums(draws) > 0)
    have_info[have_info == 0L] <- seam_transferred[have_info == 0L]
  }

  info_before_surveillance <- as.integer(have_info != 0L)

  p_background <- INApestMLUTMHostDetectionProbability(n_after_biology, detection_cube)
  background_detected <- rbinom(n_nodes, size = 1L, prob = p_background)

  p_info <- INApestMLUTMHostDetectionProbability(n_after_biology, info_detection_cube)
  if (any(info_detection_cube > 0)) {
    info_triggered_detected <- rbinom(n_nodes, size = 1L,
                                      prob = p_info * info_before_surveillance)
  } else {
    info_triggered_detected <- integer(n_nodes)
  }
  host_detection_evidence <- pmax(background_detected, info_triggered_detected)

  if (use_persistence) last_known_presence[host_detection_evidence == 1L] <- timestep
  have_info[have_info == 0L] <- host_detection_evidence[have_info == 0L]

  invaded_by_landuse <- matrix(0L, n_nodes, n_landuses)
  for (i in seq_len(n_nodes)) for (l in seq_len(n_landuses))
    invaded_by_landuse[i, l] <- as.integer(sum(n_after_biology[i, l, ]) > 0)
  node_invaded_end <- as.integer(rowSums(invaded_by_landuse) > 0)
  known_present_by_landuse <- invaded_by_landuse * have_info

  list(
    N = n_after_biology,
    NAfterManagement = n_after_management,
    ManagementDeaths = management_deaths,
    Managing = managing,
    HaveInfo = have_info,
    LastKnownPresence = last_known_presence,
    InformationStateBeforeSurveillance = info_before_surveillance,
    BackgroundDetectionProbability = p_background,
    InfoTriggeredDetectionProbability = p_info,
    BackgroundDetected = background_detected,
    InfoTriggeredDetected = info_triggered_detected,
    HostDetectionEvidence = host_detection_evidence,
    SEAMTransferred = seam_transferred,
    InvadedByLandUse = invaded_by_landuse,
    Invaded = node_invaded_end,
    KnownPresentByLandUse = known_present_by_landuse
  )
}


### END HISTORICAL LAYER: INApestMetaTransitionMatrixMultipleLandUse_response_v0.2.R


###############################################################################
### BEGIN HISTORICAL LAYER: INApestMetaTransitionMatrixMultipleLandUse_biocontrol_v0.3.R
###############################################################################

###############################################################################
### INApest MLU x Transition Matrix biocontrol adapter v0.3
### Date: 2026-10-01
###
### Frozen parents
### --------------
### * INApest_MLUTM_HostCore_v0.1.1_FROZEN_20261001 (10/10 native PASS)
### * INApest_MLUTM_Response_v0.2_FROZEN_20261001 (14/14 native PASS)
### * definitive central INApestBiocontrol.R, Git blob
###   afdc55fa132b3074e91f70225aa36b74a8b9f549
###
### Scope
### -----
### Adds serial H+B coupling to the combined node x land-use x pest-stage host.
### The biocontrol population deliberately retains the validated central state:
###
###     Q[node, agent-stage]
###
### rather than silently introducing a new land-use-structured agent state.
### Host attack components are the cross-product land-use x pest-stage cells.
### Agent TargetStage is therefore expanded across all land uses. This preserves
### exact reduction to current Transition-Matrix + biocontrol at one land use and
### current MLU + biocontrol at one pest stage.
###
### Event order in the active H+B response step:
###   1. information-gated management adoption;
###   2. management mortality;
###   3. frozen v0.1.1 host biology;
###   4. programmed information stopping / retention;
###   5. SEAM information transfer;
###   6. biocontrol release -> attack -> agent transition -> recruit -> movement;
###   7. host surveillance on the post-biocontrol host state;
###   8. detections refresh information for the next response step.
###
### Pathogen, RK, analytical, PoA/PoF and PSOCK remain outside this checkpoint.
###############################################################################

.inatmlu_bc_require <- function() {
  needed <- c(
    "INApestMLUTMResponseStep", "INApestMLUTMHostDetectionProbability",
    "INApestBiocontrol", "INApestBiocontrolAgent", "INApestBiocontrolContext"
  )
  miss <- needed[!vapply(needed, exists, logical(1), mode = "function")]
  if (length(miss))
    stop("Required source function(s) missing: ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

.inatmlu_bc_stage_names <- function(n) {
  S <- dim(n)[3L]
  dn <- dimnames(n)
  if (!is.null(dn) && length(dn) >= 3L && !is.null(dn[[3L]]) &&
      length(dn[[3L]]) == S && all(nzchar(dn[[3L]]))) return(as.character(dn[[3L]]))
  as.character(seq_len(S))
}

.inatmlu_bc_component_labels <- function(n_landuses, stage_names) {
  unlist(lapply(seq_len(n_landuses), function(l)
    paste0("landuse_", l, ":stage_", stage_names)), use.names = FALSE)
}

.inatmlu_bc_host_to_components <- function(n) {
  d <- dim(n); nn <- d[1L]; L <- d[2L]; S <- d[3L]
  stage_names <- .inatmlu_bc_stage_names(n)
  out <- matrix(0L, nn, L * S,
                dimnames = list(NULL, .inatmlu_bc_component_labels(L, stage_names)))
  cc <- 0L
  for (l in seq_len(L)) for (s in seq_len(S)) {
    cc <- cc + 1L
    out[, cc] <- as.integer(n[, l, s])
  }
  out
}

.inatmlu_bc_components_to_host <- function(x, n_nodes, n_landuses, n_stages, template = NULL) {
  if (!is.matrix(x) || !identical(dim(x), c(n_nodes, n_landuses * n_stages)))
    stop("Internal error: biocontrol host component matrix has unexpected dimensions", call. = FALSE)
  out <- array(0L, c(n_nodes, n_landuses, n_stages), dimnames = dimnames(template))
  cc <- 0L
  for (l in seq_len(n_landuses)) for (s in seq_len(n_stages)) {
    cc <- cc + 1L
    out[, l, s] <- as.integer(x[, cc])
  }
  out
}

.inatmlu_bc_base_target_stage <- function(a, stage_names, n_stages) {
  ts <- a$TargetStage
  if (is.null(ts)) {
    if (n_stages != 1L)
      stop("Every biocontrol agent must specify TargetStage when the host has more than one demographic stage", call. = FALSE)
    return(1L)
  }
  if (is.character(ts)) ts <- match(ts, stage_names)
  ts <- as.integer(ts)
  if (!length(ts) || anyNA(ts) || any(ts < 1L | ts > n_stages))
    stop("Biocontrol TargetStage contains an invalid host stage", call. = FALSE)
  unique(ts)
}

# Transform only the host-target indexing. Agent state, releases, transition,
# attack rates and node movement remain byte/semantic descendants of the central
# validated companion contract.
INApestMLUTMPrepareBiocontrol <- function(Biocontrol, n, Ntimesteps = 1L) {
  .inatmlu_bc_require()
  if (!inherits(Biocontrol, "INApestBiocontrol"))
    stop("Biocontrol must be created by INApestBiocontrol()", call. = FALSE)
  d <- .inatmlu_validate_array_state(n)
  nn <- d[1L]; L <- d[2L]; S <- d[3L]
  Ntimesteps <- as.integer(Ntimesteps)
  if (length(Ntimesteps) != 1L || is.na(Ntimesteps) || Ntimesteps < 1L)
    stop("Ntimesteps must be a positive integer", call. = FALSE)
  stage_names <- .inatmlu_bc_stage_names(n)
  comp_labels <- .inatmlu_bc_component_labels(L, stage_names)

  agents <- Biocontrol$Agents
  transformed <- lapply(agents, function(a) {
    base_ts <- .inatmlu_bc_base_target_stage(a, stage_names, S)
    aa <- a
    aa$TargetStage <- as.integer(unlist(lapply(seq_len(L), function(l)
      (l - 1L) * S + base_ts), use.names = FALSE))
    class(aa) <- class(a)
    aa
  })
  names(transformed) <- names(agents)

  wrapped <- INApestBiocontrol(transformed, TimestepLength = Biocontrol$TimestepLength)
  context <- INApestBiocontrolContext(
    Architecture = "transition", n_nodes = nn, Ntimesteps = Ntimesteps,
    host_stages = comp_labels
  )
  wrapped$Engine$Validate(context)
  structure(list(
    Biocontrol = wrapped, Context = context,
    OriginalBiocontrol = Biocontrol,
    n_nodes = nn, n_landuses = L, n_stages = S,
    stage_names = stage_names, component_labels = comp_labels
  ), class = "INApestMLUTMPreparedBiocontrol")
}

INApestMLUTMBiocontrolInitial <- function(Biocontrol, n, Ntimesteps = 1L, Prepared = NULL) {
  if (is.null(Prepared)) Prepared <- INApestMLUTMPrepareBiocontrol(Biocontrol, n, Ntimesteps)
  target <- .inatmlu_bc_host_to_components(n)
  Prepared$Biocontrol$Engine$Initial(target, Prepared$Context)
}

.inatmlu_bc_reshape_attacks <- function(attacks, n_nodes, n_landuses, n_stages, template) {
  out <- lapply(attacks, function(x) {
    .inatmlu_bc_components_to_host(x, n_nodes, n_landuses, n_stages, template)
  })
  names(out) <- names(attacks)
  out
}

INApestMLUTMBiocontrolStep <- function(
    n, Biocontrol, BiocontrolState = NULL, timestep = 1L,
    Ntimesteps = max(1L, as.integer(timestep)), Prepared = NULL) {
  .inatmlu_bc_require()
  d <- .inatmlu_validate_array_state(n)
  nn <- d[1L]; L <- d[2L]; S <- d[3L]
  if (is.null(Prepared)) Prepared <- INApestMLUTMPrepareBiocontrol(Biocontrol, n, Ntimesteps)
  if (!inherits(Prepared, "INApestMLUTMPreparedBiocontrol"))
    stop("Prepared must come from INApestMLUTMPrepareBiocontrol()", call. = FALSE)
  if (!identical(c(Prepared$n_nodes, Prepared$n_landuses, Prepared$n_stages), c(nn, L, S)))
    stop("Prepared biocontrol dimensions do not match host state", call. = FALSE)
  target <- .inatmlu_bc_host_to_components(n)
  if (is.null(BiocontrolState))
    BiocontrolState <- Prepared$Biocontrol$Engine$Initial(target, Prepared$Context)
  z <- Prepared$Biocontrol$Engine$Step(
    target, BiocontrolState, as.integer(timestep), Prepared$Context
  )
  host <- .inatmlu_bc_components_to_host(z$Target, nn, L, S, n)
  impact <- z$Impact
  impact$AttacksByAgentMLUTM <- .inatmlu_bc_reshape_attacks(
    z$Impact$AttacksByAgent, nn, L, S, n
  )
  list(N = host, State = z$State, Impact = impact, Prepared = Prepared)
}

# Full one-step H+B response wrapper. Biocontrol=NULL delegates directly to the
# frozen v0.2 function, preserving exact behaviour and random-number use.
INApestMLUTMBiocontrolResponseStep <- function(
    n,
    have_info,
    timestep = 1L,
    last_known_presence = NULL,
    manage_prob = 0,
    mortality_prob = 0,
    detection_prob = 0,
    info_triggered_detection_prob = 0,
    info_retention_prob = 1,
    info_persistence_steps = NA,
    seam = 0,
    biology_args = list(),
    Biocontrol = NULL,
    BiocontrolState = NULL,
    Ntimesteps = max(1L, as.integer(timestep)),
    BiocontrolPrepared = NULL) {

  if (!exists("INApestMLUTMResponseStep", mode = "function"))
    stop("Source the frozen v0.2 response layer before the biocontrol layer", call. = FALSE)
  if (is.null(Biocontrol)) {
    return(INApestMLUTMResponseStep(
      n = n, have_info = have_info, timestep = timestep,
      last_known_presence = last_known_presence,
      manage_prob = manage_prob, mortality_prob = mortality_prob,
      detection_prob = detection_prob,
      info_triggered_detection_prob = info_triggered_detection_prob,
      info_retention_prob = info_retention_prob,
      info_persistence_steps = info_persistence_steps,
      seam = seam, biology_args = biology_args
    ))
  }

  .inatmlu_bc_require()
  dims <- .inatmlu_validate_array_state(n)
  n_nodes <- dims[1L]; n_landuses <- dims[2L]; n_stages <- dims[3L]
  if (length(timestep) != 1L || !is.finite(timestep) || timestep < 1 || timestep != floor(timestep))
    stop("timestep must be a positive integer", call. = FALSE)
  have_info <- .inatmlu_binary_info(have_info, n_nodes, "have_info")
  if (is.null(last_known_presence)) last_known_presence <- rep(NA_real_, n_nodes)
  last_known_presence <- as.numeric(last_known_presence)
  if (length(last_known_presence) != n_nodes || any(!is.na(last_known_presence) & !is.finite(last_known_presence)))
    stop("last_known_presence must be NULL or length nodes", call. = FALSE)
  retention <- .inatmlu_node_prob(info_retention_prob, n_nodes, "info_retention_prob")
  persistence <- .inatmlu_persistence(info_persistence_steps, n_nodes)
  use_persistence <- any(!is.na(persistence))

  if (!is.list(biology_args)) stop("biology_args must be a named list", call. = FALSE)
  if (length(biology_args) && (is.null(names(biology_args)) || any(!nzchar(names(biology_args)))))
    stop("biology_args must be a named list", call. = FALSE)
  if (anyDuplicated(names(biology_args))) stop("biology_args names must be unique", call. = FALSE)
  reserved <- intersect(names(biology_args), c("n0", "managing"))
  if (length(reserved)) stop("biology_args may not override n0 or managing", call. = FALSE)

  manage_surface <- .inatmlu_response_surface(manage_prob, n_nodes, n_landuses, "manage_prob")
  mortality_cube <- .inatmlu_response_cube(mortality_prob, n_nodes, n_landuses, n_stages, "mortality_prob")
  detection_cube <- .inatmlu_response_cube(detection_prob, n_nodes, n_landuses, n_stages, "detection_prob")
  info_detection_cube <- .inatmlu_response_cube(info_triggered_detection_prob, n_nodes, n_landuses, n_stages, "info_triggered_detection_prob")

  node_invaded_start <- as.integer(vapply(seq_len(n_nodes), function(i) sum(n[i, , ]) > 0, logical(1)))
  known_occupied_start <- node_invaded_start * have_info

  manage_p <- as.numeric(manage_surface) * rep(have_info, times = n_landuses)
  managing <- matrix(rbinom(n_nodes * n_landuses, size = 1L, prob = manage_p),
                     nrow = n_nodes, ncol = n_landuses)
  managing_cube <- array(rep(as.numeric(managing), times = n_stages),
                         dim = c(n_nodes, n_landuses, n_stages))
  survive_p <- 1 - mortality_cube * managing_cube
  n_after_management <- array(
    rbinom(length(n), size = as.integer(n), prob = as.numeric(survive_p)),
    dim = dim(n), dimnames = dimnames(n)
  )
  management_deaths <- n - n_after_management

  if (use_persistence) {
    killed_by_node <- vapply(seq_len(n_nodes), function(i) sum(management_deaths[i, , ]) > 0, logical(1))
    last_known_presence[killed_by_node] <- timestep
  }

  call_args <- c(list(n0 = n_after_management, managing = managing), biology_args)
  n_after_biology <- do.call(INApestMLUTMHostStep, call_args)

  programmed <- which(have_info == 1L & !is.na(persistence))
  if (length(programmed)) {
    elapsed <- timestep - last_known_presence
    stop_nodes <- programmed[is.na(last_known_presence[programmed]) |
                             elapsed[programmed] >= persistence[programmed]]
    if (length(stop_nodes)) have_info[stop_nodes] <- 0L
  }
  decay <- which(have_info == 1L & is.na(persistence) & retention < 1)
  if (length(decay))
    have_info[decay] <- rbinom(length(decay), size = 1L, prob = retention[decay])

  seam_transferred <- integer(n_nodes)
  if (!(length(seam) == 1L && identical(as.numeric(seam), 0))) {
    if (!is.matrix(seam) || !identical(dim(seam), c(n_nodes, n_nodes)))
      stop("seam must be 0 or a nodes x nodes matrix", call. = FALSE)
    if (any(!is.finite(seam)) || any(seam < 0 | seam > 1))
      stop("seam entries must be probabilities in [0,1]", call. = FALSE)
    Sx <- seam; diag(Sx) <- 0
    p_transfer <- Sx * known_occupied_start
    draws <- matrix(rbinom(n_nodes * n_nodes, 1L, prob = as.numeric(p_transfer)),
                    nrow = n_nodes, ncol = n_nodes)
    seam_transferred <- as.integer(colSums(draws) > 0)
    have_info[have_info == 0L] <- seam_transferred[have_info == 0L]
  }
  info_before_surveillance <- as.integer(have_info != 0L)

  bc <- INApestMLUTMBiocontrolStep(
    n_after_biology, Biocontrol = Biocontrol, BiocontrolState = BiocontrolState,
    timestep = timestep, Ntimesteps = Ntimesteps, Prepared = BiocontrolPrepared
  )
  n_after_biocontrol <- bc$N

  p_background <- INApestMLUTMHostDetectionProbability(n_after_biocontrol, detection_cube)
  background_detected <- rbinom(n_nodes, size = 1L, prob = p_background)
  p_info <- INApestMLUTMHostDetectionProbability(n_after_biocontrol, info_detection_cube)
  if (any(info_detection_cube > 0)) {
    info_triggered_detected <- rbinom(n_nodes, size = 1L,
                                      prob = p_info * info_before_surveillance)
  } else {
    info_triggered_detected <- integer(n_nodes)
  }
  host_detection_evidence <- pmax(background_detected, info_triggered_detected)
  if (use_persistence) last_known_presence[host_detection_evidence == 1L] <- timestep
  have_info[have_info == 0L] <- host_detection_evidence[have_info == 0L]

  invaded_by_landuse <- matrix(0L, n_nodes, n_landuses)
  for (i in seq_len(n_nodes)) for (l in seq_len(n_landuses))
    invaded_by_landuse[i, l] <- as.integer(sum(n_after_biocontrol[i, l, ]) > 0)
  node_invaded_end <- as.integer(rowSums(invaded_by_landuse) > 0)
  known_present_by_landuse <- invaded_by_landuse * have_info

  list(
    N = n_after_biocontrol,
    NAfterManagement = n_after_management,
    NAfterBiology = n_after_biology,
    ManagementDeaths = management_deaths,
    Managing = managing,
    HaveInfo = have_info,
    LastKnownPresence = last_known_presence,
    InformationStateBeforeSurveillance = info_before_surveillance,
    BackgroundDetectionProbability = p_background,
    InfoTriggeredDetectionProbability = p_info,
    BackgroundDetected = background_detected,
    InfoTriggeredDetected = info_triggered_detected,
    HostDetectionEvidence = host_detection_evidence,
    SEAMTransferred = seam_transferred,
    InvadedByLandUse = invaded_by_landuse,
    Invaded = node_invaded_end,
    KnownPresentByLandUse = known_present_by_landuse,
    BiocontrolState = bc$State,
    BiocontrolImpact = bc$Impact,
    BiocontrolPrepared = bc$Prepared
  )
}


### END HISTORICAL LAYER: INApestMetaTransitionMatrixMultipleLandUse_biocontrol_v0.3.R


###############################################################################
### BEGIN HISTORICAL LAYER: INApestMetaTransitionMatrixMultipleLandUse_pathogen_biocontrol_v0.4.R
###############################################################################

###############################################################################
### INApest MLU x Transition Matrix pathogen + biocontrol adapter v0.4
### Date: 2026-10-01
###
### Frozen parents
### --------------
### * MLUTM host core v0.1.1 (native 10/10 PASS)
### * MLUTM response v0.2 (native 14/14 PASS)
### * MLUTM biocontrol v0.3 (native 14/14 PASS)
### * definitive central INApestBiocontrol.R
### * definitive INApestPathogen.R
### * definitive INApestPathogenTransitionMatrix.R
###
### Scope
### -----
### Adds serial H+P+B coupling to the combined host state
###
###   H[node, land-use, host-stage]
###   P[node, land-use, host-stage, pathogen-state]
###   Q[node, agent-stage]
###
### For the genuine cross-product path (>=2 land uses and >=2 host stages),
### each node x land-use host cell becomes a pseudo-node for the existing
### stage-aware pathogen helper. Pathogen state therefore follows demographic
### survival, stage progression and transition-associated movement rather than
### being reassigned after the host step.
###
### Biocontrol retains the frozen v0.3 node x agent-stage state. After attack,
### pathogen compartments are reconciled without replacement to the surviving
### hosts in each node x land-use x host-stage cell.
###
### RK, analytical, PoA/PoF and PSOCK remain outside this checkpoint.
###############################################################################

.inatmlu_p_require <- function() {
  needed <- c(
    "INApestMLUTMHostStep", "INApestMLUTMBiocontrolResponseStep",
    "INApestMLUTMBiocontrolStep", "INApestPathogenStageState",
    "local.dynamics.transition.matrix.pathogen", ".iptm_reconcile"
  )
  miss <- needed[!vapply(needed, exists, logical(1), mode = "function")]
  if (length(miss))
    stop("Required source function(s) missing: ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

.inatmlu_p_validate_pathogen <- function(Pathogen) {
  if (!inherits(Pathogen, "INApestPathogen"))
    stop("Pathogen must be created by INApestPathogen()", call. = FALSE)
  if (is.null(Pathogen$States) || !("S" %in% Pathogen$States) || !("I" %in% Pathogen$States))
    stop("v0.4 requires a compartmental INApestPathogen with S and I states", call. = FALSE)
  invisible(TRUE)
}

.inatmlu_p_validate_state <- function(PathogenState, n, Pathogen) {
  .inatmlu_p_validate_pathogen(Pathogen)
  d <- .inatmlu_validate_array_state(n)
  expected <- c(d, length(Pathogen$States))
  if (length(dim(PathogenState)) != 4L || !identical(as.integer(dim(PathogenState)), expected))
    stop("PathogenState must be nodes x land-uses x host-stages x pathogen-states", call. = FALSE)
  if (any(!is.finite(PathogenState)) || any(PathogenState < 0) ||
      any(PathogenState != floor(PathogenState)))
    stop("PathogenState must contain finite non-negative integer counts", call. = FALSE)
  if (any(apply(PathogenState, c(1, 2, 3), sum) != n))
    stop("PathogenState must sum to host N within every node x land-use x host-stage cell", call. = FALSE)
  dn <- dimnames(PathogenState)
  if (is.null(dn)) dn <- vector("list", 4L)
  dn[[4L]] <- Pathogen$States
  dimnames(PathogenState) <- dn
  PathogenState
}

# MLUTM host core uses node-major / land-use-fast pseudo-node ordering.
.inatmlu_p_flatten_tm <- function(PathogenState) {
  d <- dim(PathogenState); nn <- d[1L]; L <- d[2L]; S <- d[3L]; P <- d[4L]
  out <- array(0L, c(nn * L, S, P),
               dimnames = list(NULL, dimnames(PathogenState)[[3L]], dimnames(PathogenState)[[4L]]))
  for (i in seq_len(nn)) for (l in seq_len(L))
    out[.inatmlu_flat_index(i, l, L), , ] <- PathogenState[i, l, , ]
  out
}

.inatmlu_p_unflatten_tm <- function(x, nn, L, S, states, template = NULL) {
  if (length(dim(x)) != 3L || !identical(as.integer(dim(x)), c(nn * L, S, length(states))))
    stop("Internal error: flattened pathogen state has unexpected dimensions", call. = FALSE)
  dn <- if (is.null(template)) list(NULL, NULL, NULL, states) else dimnames(template)
  if (is.null(dn)) dn <- vector("list", 4L)
  dn[[4L]] <- states
  out <- array(0L, c(nn, L, S, length(states)), dimnames = dn)
  for (i in seq_len(nn)) for (l in seq_len(L))
    out[i, l, , ] <- x[.inatmlu_flat_index(i, l, L), , ]
  out
}

# Existing MLU pathogen engine uses R matrix flattening: node varies fastest.
.inatmlu_p_flatten_mlu_one_stage <- function(PathogenState) {
  d <- dim(PathogenState); nn <- d[1L]; L <- d[2L]; P <- d[4L]
  out <- matrix(0L, nn * L, P, dimnames = list(NULL, dimnames(PathogenState)[[4L]]))
  for (q in seq_len(P)) out[, q] <- as.integer(PathogenState[, , 1L, q])
  out
}

.inatmlu_p_unflatten_mlu_one_stage <- function(x, nn, L, states, template = NULL) {
  if (!is.matrix(x) || !identical(dim(x), c(nn * L, length(states))))
    stop("Internal error: MLU pathogen state has unexpected dimensions", call. = FALSE)
  dn <- if (is.null(template)) list(NULL, NULL, NULL, states) else dimnames(template)
  if (is.null(dn)) dn <- vector("list", 4L)
  dn[[4L]] <- states
  out <- array(0L, c(nn, L, 1L, length(states)), dimnames = dn)
  for (q in seq_along(states)) out[, , 1L, q] <- matrix(x[, q], nn, L)
  out
}

.inatmlu_p_unflatten_metric <- function(x, nn, L, S) {
  if (is.null(x)) return(NULL)
  if (is.null(dim(x))) {
    if (length(x) != nn * L * S) stop("Internal pathogen metric length mismatch", call. = FALSE)
    x <- matrix(x, nrow = nn * L, ncol = S)
  }
  if (!identical(as.integer(dim(x)), c(nn * L, S)))
    stop("Internal pathogen metric dimensions mismatch", call. = FALSE)
  out <- array(0L, c(nn, L, S))
  for (i in seq_len(nn)) for (l in seq_len(L))
    out[i, l, ] <- as.integer(x[.inatmlu_flat_index(i, l, L), ])
  out
}

.inatmlu_p_landuse_mixing <- function(x, L) {
  if (is.null(x)) return(diag(L))
  if (!is.matrix(x) || !identical(dim(x), c(L, L)) ||
      any(!is.finite(x)) || any(x < 0))
    stop("PathogenLandUseMixing must be a finite non-negative land-uses x land-uses matrix", call. = FALSE)
  x
}

# Static v0.4 parameter expansion for the genuine cross-product path.
.inatmlu_p_expand_parameter <- function(x, nn, L, S, name) {
  if (is.function(x))
    stop(name, " resolver functions are reserved for a later dynamic-input gate; use a static value in v0.4", call. = FALSE)
  nc <- nn * L
  d <- dim(x)
  if (is.null(d)) {
    z <- as.numeric(x)
    if (length(z) == 1L) return(z)
    labels <- c(stage = S, node = nn, landuse = L, cell = nc)
    hit <- names(labels)[labels == length(z)]
    if (length(hit) != 1L)
      stop(name, " vector length is unsupported or ambiguous; supply an explicit array/matrix", call. = FALSE)
    if (hit == "stage") return(z)
    out <- matrix(0, nc, S)
    if (hit == "node") {
      for (i in seq_len(nn)) for (l in seq_len(L)) out[.inatmlu_flat_index(i,l,L), ] <- z[i]
    } else if (hit == "landuse") {
      for (i in seq_len(nn)) for (l in seq_len(L)) out[.inatmlu_flat_index(i,l,L), ] <- z[l]
    } else {
      out[] <- rep(z, S)
    }
    return(out)
  }
  if (length(d) == 2L && identical(as.integer(d), c(nc, S))) return(matrix(as.numeric(x), nc, S))
  if (length(d) == 2L && identical(as.integer(d), c(nn, S))) {
    out <- matrix(0, nc, S)
    for (i in seq_len(nn)) for (l in seq_len(L)) out[.inatmlu_flat_index(i,l,L), ] <- x[i, ]
    return(out)
  }
  if (length(d) == 2L && identical(as.integer(d), c(nn, L))) {
    z <- .inatmlu_flatten_surface(matrix(as.numeric(x), nn, L))
    return(matrix(rep(z, S), nrow = nc, ncol = S))
  }
  if (length(d) == 3L && identical(as.integer(d), c(nn, L, S)))
    return(.inatmlu_flatten_state(array(as.numeric(x), c(nn, L, S))))
  stop(name, " must be scalar, an unambiguous static vector, nodes x stages, nodes x land-uses, nodes x land-uses x stages, or cells x stages", call. = FALSE)
}

.inatmlu_p_prepare_combined <- function(Pathogen, nn, L, S, PathogenLandUseMixing = NULL) {
  .inatmlu_p_validate_pathogen(Pathogen)
  out <- Pathogen
  for (nm in c("Beta", "RecoveryProb", "ProgressionProb", "PathogenMortalityProb",
               "ImmunityLossProb", "IntroductionProb", "IntroductionNumber", "DensityScale")) {
    out[[nm]] <- .inatmlu_p_expand_parameter(Pathogen[[nm]], nn, L, S, nm)
  }
  C <- Pathogen$ContactMatrix
  if (is.function(C))
    stop("Pathogen ContactMatrix resolver functions are reserved for a later dynamic-input gate in the crossed path", call. = FALSE)
  if (!is.null(C)) {
    if (!is.matrix(C) || any(!is.finite(C)) || any(C < 0))
      stop("Pathogen ContactMatrix must be a finite non-negative matrix", call. = FALSE)
    nc <- nn * L
    if (identical(dim(C), c(nn, nn))) {
      LU <- .inatmlu_p_landuse_mixing(PathogenLandUseMixing, L)
      # MLUTM pseudo-node order is node-major, land-use-fast.
      out$ContactMatrix <- kronecker(C, LU)
    } else if (identical(dim(C), c(nc, nc))) {
      out$ContactMatrix <- C
    } else {
      stop("For crossed MLUTM, Pathogen ContactMatrix must be nodes x nodes or MLUTM-cells x MLUTM-cells", call. = FALSE)
    }
  } else {
    out$ContactMatrix <- NULL
  }
  class(out) <- class(Pathogen)
  out
}

INApestMLUTMPathogenInitial <- function(n, Pathogen, InitialState = NULL, Ntimesteps = 1L) {
  .inatmlu_p_require(); .inatmlu_p_validate_pathogen(Pathogen)
  d <- .inatmlu_validate_array_state(n); nn <- d[1L]; L <- d[2L]; S <- d[3L]
  states <- Pathogen$States
  if (!is.null(InitialState)) return(.inatmlu_p_validate_state(InitialState, n, Pathogen))

  if (L == 1L) {
    z <- INApestPathogenStageState(matrix(n[,1L,], nn, S), Pathogen, NULL, Ntimesteps)
    return(.inatmlu_p_unflatten_tm(z, nn, 1L, S, states))
  }
  if (S == 1L) {
    ctx <- list(n_nodes = nn, n_landuses = L, Ntimesteps = as.integer(Ntimesteps))
    Pathogen$Engine$Validate(ctx)
    z <- Pathogen$Engine$Initial(matrix(n[,,1L], nn, L), ctx)
    return(.inatmlu_p_unflatten_mlu_one_stage(z, nn, L, states))
  }

  simple_zero <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == 0
  if (!all(vapply(list(Pathogen$InitialInfected, Pathogen$InitialExposed, Pathogen$InitialRecovered),
                  simple_zero, logical(1))))
    stop("For crossed MLUTM v0.4, non-zero/structured initial infection must be supplied explicitly as InitialState", call. = FALSE)
  out <- array(0L, c(nn, L, S, length(states)), dimnames = list(NULL,NULL,NULL,states))
  out[,,,"S"] <- n
  out
}

INApestMLUTMPathogenReconcile <- function(PathogenState, n, Pathogen) {
  .inatmlu_p_require()
  current_n <- array(apply(PathogenState,c(1,2,3),sum),dim=dim(PathogenState)[1:3])
  state <- .inatmlu_p_validate_state(PathogenState, current_n, Pathogen)
  d <- .inatmlu_validate_array_state(n); nn <- d[1L]; L <- d[2L]; S <- d[3L]
  # Preserve the definitive MLU random-thinning order at the one-stage boundary.
  if (S == 1L) {
    ctx <- list(n_nodes=nn,n_landuses=L,Ntimesteps=1L)
    pf <- .inatmlu_p_flatten_mlu_one_stage(state)
    pf <- Pathogen$Engine$Reconcile(pf,matrix(n[,,1L],nn,L),ctx)
    return(.inatmlu_p_unflatten_mlu_one_stage(pf,nn,L,Pathogen$States,state))
  }
  flat <- .inatmlu_p_flatten_tm(state)
  nflat <- .inatmlu_flatten_state(n)
  flat <- .iptm_reconcile(flat, nflat)
  .inatmlu_p_unflatten_tm(flat, nn, L, S, Pathogen$States, state)
}

# One pathogen-aware host biology step. Host management mortality belongs to the
# response wrapper and must be reconciled before entering this function.
INApestMLUTMPathogenHostStep <- function(
    n0, PathogenState, Pathogen, timestep = 1L, Ntimesteps = max(1L, as.integer(timestep)),
    StageMixing = NULL, PathogenLandUseMixing = NULL,
    nodetransition, weights = NULL, sddprob, nodeenvestabprob = 1,
    lddprob = NA, lddrate = 0, nodeK, node.seedbankK = nodeK,
    nodepropaguleestablishment = 1, nodespreadreduction = 0,
    nodefecundityreduction = 0, managing = 0,
    transition_sddprob = NULL, transition_lddprob = NULL, transition_lddrate = 0,
    ApplyFootprintToTransitions = FALSE, BlockedTransitionMortality = 0,
    DispersalDensityFactor = 0, MaxInteger = .Machine$integer.max,
    LandUseRecruitmentWeights = NULL, TransitionLandUseMixing = NULL,
    MLUPropaguleProduction = NULL) {

  .inatmlu_p_require(); .inatmlu_p_validate_pathogen(Pathogen)
  dims <- .inatmlu_validate_array_state(n0)
  nn <- dims[1L]; L <- dims[2L]; S <- dims[3L]
  PathogenState <- .inatmlu_p_validate_state(PathogenState, n0, Pathogen)

  # Exact one-stage MLU pathogen boundary.
  if (S == 1L) {
    host <- INApestMLUTMHostStep(
      n0=n0, nodetransition=nodetransition, weights=weights, sddprob=sddprob,
      nodeenvestabprob=nodeenvestabprob, lddprob=lddprob, lddrate=lddrate,
      nodeK=nodeK, node.seedbankK=node.seedbankK,
      nodepropaguleestablishment=nodepropaguleestablishment,
      nodespreadreduction=nodespreadreduction, nodefecundityreduction=nodefecundityreduction,
      managing=managing, transition_sddprob=transition_sddprob,
      transition_lddprob=transition_lddprob, transition_lddrate=transition_lddrate,
      ApplyFootprintToTransitions=ApplyFootprintToTransitions,
      BlockedTransitionMortality=BlockedTransitionMortality,
      DispersalDensityFactor=DispersalDensityFactor, MaxInteger=MaxInteger,
      LandUseRecruitmentWeights=LandUseRecruitmentWeights,
      TransitionLandUseMixing=TransitionLandUseMixing,
      MLUPropaguleProduction=MLUPropaguleProduction
    )
    ctx <- list(n_nodes=nn, n_landuses=L, Ntimesteps=as.integer(Ntimesteps))
    Pathogen$Engine$Validate(ctx)
    pf <- .inatmlu_p_flatten_mlu_one_stage(PathogenState)
    hmat <- matrix(host[,,1L], nn, L)
    pf <- Pathogen$Engine$Reconcile(pf, hmat, ctx)
    ps <- Pathogen$Engine$Step(pf, hmat, as.integer(timestep), ctx)
    hout <- array(ps$N, c(nn,L,1L), dimnames=dimnames(host))
    pout <- .inatmlu_p_unflatten_mlu_one_stage(ps$State, nn, L, Pathogen$States, PathogenState)
    deaths <- if (is.null(ps$Deaths)) array(0L,c(nn,L,1L)) else array(ps$Deaths,c(nn,L,1L))
    newinf <- if (is.null(ps$NewInfections)) array(0L,c(nn,L,1L)) else array(ps$NewInfections,c(nn,L,1L))
    return(list(N=hout, PathogenState=pout, PathogenDeaths=deaths, NewInfections=newinf,
                PathogenIntroduced=if(is.null(ps$Introduced)) NULL else array(ps$Introduced,c(nn,L,1L))))
  }

  # Exact one-land-use Transition-Matrix pathogen boundary.
  if (L == 1L) {
    n_tm <- matrix(n0[,1L,], nn, S)
    p_tm <- array(PathogenState[,1L,,], c(nn,S,length(Pathogen$States)),
                  dimnames=list(NULL,dimnames(PathogenState)[[3L]],Pathogen$States))
    K <- as.numeric(.inatmlu_surface(nodeK,nn,1L,"nodeK")[,1L])
    seed <- as.numeric(.inatmlu_surface(node.seedbankK,nn,1L,"node.seedbankK")[,1L])
    env <- as.numeric(.inatmlu_surface(nodeenvestabprob,nn,1L,"nodeenvestabprob")[,1L])
    estab <- as.numeric(.inatmlu_surface(nodepropaguleestablishment,nn,1L,"nodepropaguleestablishment")[,1L])
    spread <- as.numeric(.inatmlu_probability_surface(nodespreadreduction,nn,1L,"nodespreadreduction")[,1L])
    manage <- as.numeric(.inatmlu_probability_surface(managing,nn,1L,"managing")[,1L])
    fec <- nodefecundityreduction
    if (!is.null(dim(fec))) {
      if (identical(dim(fec),c(nn,1L))) fec <- as.numeric(fec[,1L])
      else if (identical(dim(fec),c(nn,1L,S))) fec <- matrix(fec[,1L,],nn,S)
    }
    z <- local.dynamics.transition.matrix.pathogen(
      nodetransition=nodetransition, weights=weights, sddprob=sddprob,
      nodeenvestabprob=env, n0=n_tm, lddprob=lddprob, lddrate=lddrate,
      transition_sddprob=transition_sddprob, transition_lddprob=transition_lddprob,
      transition_lddrate=transition_lddrate, nodeK=K, node.seedbankK=seed,
      nodepropaguleestablishment=estab, nodespreadreduction=spread,
      nodefecundityreduction=fec, managing=manage, MaxInteger=MaxInteger,
      BlockedTransitionMortality=BlockedTransitionMortality,
      DispersalDensityFactor=DispersalDensityFactor,
      pathogen_state=p_tm, Pathogen=Pathogen, timestep=as.integer(timestep),
      Ntimesteps=as.integer(Ntimesteps), StageMixing=StageMixing
    )
    return(list(
      N=array(z$N,c(nn,1L,S),dimnames=dimnames(n0)),
      PathogenState=.inatmlu_p_unflatten_tm(z$PathogenState,nn,1L,S,Pathogen$States,PathogenState),
      PathogenDeaths=array(z$PathogenDeaths,c(nn,1L,S)),
      NewInfections=array(z$NewInfections,c(nn,1L,S)),
      PathogenIntroduced=if(is.null(z$PathogenIntroduced)) NULL else array(z$PathogenIntroduced,c(nn,1L,S))
    ))
  }

  # Genuine crossed path: flatten node x land-use to pseudo-nodes.
  .inatmlu_validate_node_connectivity(sddprob, nn, "sddprob")
  .inatmlu_validate_node_connectivity(lddprob, nn, "lddprob", allow_disabled=TRUE)
  nc <- nn * L
  nflat <- .inatmlu_flatten_state(n0)
  pflat <- .inatmlu_p_flatten_tm(PathogenState)
  Klu <- .inatmlu_surface(nodeK,nn,L,"nodeK")
  seedlu <- .inatmlu_surface(node.seedbankK,nn,L,"node.seedbankK")
  envlu <- .inatmlu_probability_surface(nodeenvestabprob,nn,L,"nodeenvestabprob")
  establu <- .inatmlu_probability_surface(nodepropaguleestablishment,nn,L,"nodepropaguleestablishment")
  spreadlu <- .inatmlu_probability_surface(nodespreadreduction,nn,L,"nodespreadreduction")
  managelu <- .inatmlu_probability_surface(managing,nn,L,"managing")
  fec <- .inatmlu_fecundity_reduction(nodefecundityreduction,nn,L,S)
  w <- .inatmlu_weights(weights,nn,L,S)
  A <- .inatmlu_transition(nodetransition,nn,L,S)
  rw <- .inatmlu_recruitment_weights(LandUseRecruitmentWeights,seedlu,nn,L)
  sdd <- .inatmlu_expand_reproductive_dispersal(sddprob,rw)
  ldd <- if(length(lddprob)==1L && (is.na(lddprob) || identical(as.numeric(lddprob),0))) NA else .inatmlu_expand_reproductive_dispersal(lddprob,rw)
  lumix <- .inatmlu_landuse_mixing(TransitionLandUseMixing,L)
  tsdd <- .inatmlu_expand_transition_movement(transition_sddprob,lumix,nn,S,"transition_sddprob")
  tldd <- .inatmlu_expand_transition_movement(transition_lddprob,lumix,nn,S,"transition_lddprob")
  blocked <- .inatmlu_blocked_mortality(BlockedTransitionMortality,nn,L,S)
  P2 <- .inatmlu_p_prepare_combined(Pathogen,nn,L,S,PathogenLandUseMixing)

  z <- local.dynamics.transition.matrix.pathogen(
    nodetransition=A, weights=w, sddprob=sdd,
    nodeenvestabprob=.inatmlu_flatten_surface(envlu), n0=nflat,
    lddprob=ldd, lddrate=lddrate,
    transition_sddprob=tsdd, transition_lddprob=tldd,
    transition_lddrate=transition_lddrate,
    nodeK=.inatmlu_flatten_surface(Klu), node.seedbankK=.inatmlu_flatten_surface(seedlu),
    nodepropaguleestablishment=.inatmlu_flatten_surface(establu),
    nodespreadreduction=.inatmlu_flatten_surface(spreadlu),
    nodefecundityreduction=fec, managing=.inatmlu_flatten_surface(managelu),
    MaxInteger=MaxInteger, BlockedTransitionMortality=blocked,
    DispersalDensityFactor=DispersalDensityFactor,
    pathogen_state=pflat, Pathogen=P2, timestep=as.integer(timestep),
    Ntimesteps=as.integer(Ntimesteps), StageMixing=StageMixing
  )
  host <- .inatmlu_unflatten_state(z$N,nn,L,S)
  dimnames(host) <- dimnames(n0)
  pout <- .inatmlu_p_unflatten_tm(z$PathogenState,nn,L,S,Pathogen$States,PathogenState)
  list(N=host, PathogenState=pout,
       PathogenDeaths=.inatmlu_p_unflatten_metric(z$PathogenDeaths,nn,L,S),
       NewInfections=.inatmlu_p_unflatten_metric(z$NewInfections,nn,L,S),
       PathogenIntroduced=.inatmlu_p_unflatten_metric(z$PathogenIntroduced,nn,L,S))
}

.inatmlu_p_detection_cube <- function(x, nn, L, S) {
  .inatmlu_response_cube(x, nn, L, S, "Pathogen DetectionProb")
}

INApestMLUTMPathogenDetectionProbability <- function(PathogenState, DetectionProb) {
  d <- dim(PathogenState); nn <- d[1L]; L <- d[2L]; S <- d[3L]
  st <- dimnames(PathogenState)[[4L]]
  if (is.null(st) || !("I" %in% st)) stop("PathogenState must name an I compartment", call. = FALSE)
  p <- .inatmlu_p_detection_cube(DetectionProb,nn,L,S)
  I <- PathogenState[,,,"I",drop=TRUE]
  if (length(dim(I)) != 3L) I <- array(I,c(nn,L,S))
  out <- numeric(nn)
  for(i in seq_len(nn)) out[i] <- 1 - prod((1-p[i,,])^I[i,,])
  pmin(1,pmax(0,out))
}

# Full serial H+P+B response step. Pathogen=NULL delegates to frozen v0.3 exactly.
INApestMLUTMPathogenBiocontrolResponseStep <- function(
    n, have_info, timestep=1L, last_known_presence=NULL,
    manage_prob=0, mortality_prob=0,
    detection_prob=0, info_triggered_detection_prob=0,
    info_retention_prob=1, info_persistence_steps=NA, seam=0,
    biology_args=list(),
    Pathogen=NULL, PathogenState=NULL, StageMixing=NULL,
    PathogenLandUseMixing=NULL,
    Biocontrol=NULL, BiocontrolState=NULL,
    Ntimesteps=max(1L,as.integer(timestep)), BiocontrolPrepared=NULL) {

  if (is.null(Pathogen)) {
    return(INApestMLUTMBiocontrolResponseStep(
      n=n, have_info=have_info, timestep=timestep, last_known_presence=last_known_presence,
      manage_prob=manage_prob, mortality_prob=mortality_prob,
      detection_prob=detection_prob, info_triggered_detection_prob=info_triggered_detection_prob,
      info_retention_prob=info_retention_prob, info_persistence_steps=info_persistence_steps,
      seam=seam, biology_args=biology_args,
      Biocontrol=Biocontrol, BiocontrolState=BiocontrolState,
      Ntimesteps=Ntimesteps, BiocontrolPrepared=BiocontrolPrepared
    ))
  }

  .inatmlu_p_require(); .inatmlu_p_validate_pathogen(Pathogen)
  dims <- .inatmlu_validate_array_state(n); nn<-dims[1L]; L<-dims[2L]; S<-dims[3L]
  if (length(timestep)!=1L || !is.finite(timestep) || timestep<1 || timestep!=floor(timestep))
    stop("timestep must be a positive integer", call. = FALSE)
  have_info <- .inatmlu_binary_info(have_info,nn,"have_info")
  if(is.null(last_known_presence)) last_known_presence <- rep(NA_real_,nn)
  last_known_presence <- as.numeric(last_known_presence)
  if(length(last_known_presence)!=nn || any(!is.na(last_known_presence)&!is.finite(last_known_presence)))
    stop("last_known_presence must be NULL or length nodes",call.=FALSE)
  retention <- .inatmlu_node_prob(info_retention_prob,nn,"info_retention_prob")
  persistence <- .inatmlu_persistence(info_persistence_steps,nn)
  use_persistence <- any(!is.na(persistence))

  if(is.null(PathogenState)) PathogenState <- INApestMLUTMPathogenInitial(n,Pathogen,NULL,Ntimesteps)
  PathogenState <- .inatmlu_p_validate_state(PathogenState,n,Pathogen)
  if(!is.list(biology_args)) stop("biology_args must be a named list",call.=FALSE)
  if(length(biology_args) && (is.null(names(biology_args)) || any(!nzchar(names(biology_args))))) stop("biology_args must be a named list",call.=FALSE)
  if(anyDuplicated(names(biology_args))) stop("biology_args names must be unique",call.=FALSE)
  reserved <- intersect(names(biology_args),c("n0","managing","PathogenState","Pathogen","timestep","Ntimesteps","StageMixing","PathogenLandUseMixing"))
  if(length(reserved)) stop("biology_args may not override pathogen/response core arguments: ",paste(reserved,collapse=", "),call.=FALSE)

  manage_surface <- .inatmlu_response_surface(manage_prob,nn,L,"manage_prob")
  mortality_cube <- .inatmlu_response_cube(mortality_prob,nn,L,S,"mortality_prob")
  detection_cube <- .inatmlu_response_cube(detection_prob,nn,L,S,"detection_prob")
  info_detection_cube <- .inatmlu_response_cube(info_triggered_detection_prob,nn,L,S,"info_triggered_detection_prob")
  node_invaded_start <- as.integer(vapply(seq_len(nn),function(i) sum(n[i,,])>0,logical(1)))
  known_occupied_start <- node_invaded_start * have_info

  manage_p <- as.numeric(manage_surface) * rep(have_info,times=L)
  managing <- matrix(rbinom(nn*L,1L,manage_p),nn,L)
  managing_cube <- array(rep(as.numeric(managing),times=S),c(nn,L,S))
  survive_p <- 1 - mortality_cube * managing_cube
  n_after_management <- array(rbinom(length(n),size=as.integer(n),prob=as.numeric(survive_p)),dim=dim(n),dimnames=dimnames(n))
  management_deaths <- n - n_after_management
  if(use_persistence) {
    killed <- vapply(seq_len(nn),function(i) sum(management_deaths[i,,])>0,logical(1))
    last_known_presence[killed] <- timestep
  }
  p_after_management <- INApestMLUTMPathogenReconcile(PathogenState,n_after_management,Pathogen)

  call_args <- c(list(n0=n_after_management,PathogenState=p_after_management,Pathogen=Pathogen,
                      timestep=timestep,Ntimesteps=Ntimesteps,StageMixing=StageMixing,
                      PathogenLandUseMixing=PathogenLandUseMixing,managing=managing),biology_args)
  bio <- do.call(INApestMLUTMPathogenHostStep,call_args)
  n_after_biology <- bio$N
  p_after_biology <- bio$PathogenState

  programmed <- which(have_info==1L & !is.na(persistence))
  if(length(programmed)) {
    elapsed <- timestep-last_known_presence
    stop_nodes <- programmed[is.na(last_known_presence[programmed]) | elapsed[programmed]>=persistence[programmed]]
    if(length(stop_nodes)) have_info[stop_nodes] <- 0L
  }
  decay <- which(have_info==1L & is.na(persistence) & retention<1)
  if(length(decay)) have_info[decay] <- rbinom(length(decay),1L,retention[decay])

  seam_transferred <- integer(nn)
  if(!(length(seam)==1L && identical(as.numeric(seam),0))) {
    if(!is.matrix(seam) || !identical(dim(seam),c(nn,nn))) stop("seam must be 0 or a nodes x nodes matrix",call.=FALSE)
    if(any(!is.finite(seam)) || any(seam<0|seam>1)) stop("seam entries must be probabilities in [0,1]",call.=FALSE)
    Sx<-seam; diag(Sx)<-0; draws<-matrix(rbinom(nn*nn,1L,as.numeric(Sx*known_occupied_start)),nn,nn)
    seam_transferred<-as.integer(colSums(draws)>0); have_info[have_info==0L]<-seam_transferred[have_info==0L]
  }

  if(is.null(Biocontrol)) {
    n_after_biocontrol <- n_after_biology
    p_after_biocontrol <- p_after_biology
    bc_state <- NULL; bc_impact <- NULL; bc_prepared <- NULL
  } else {
    bc <- INApestMLUTMBiocontrolStep(n_after_biology,Biocontrol=Biocontrol,BiocontrolState=BiocontrolState,
                                     timestep=timestep,Ntimesteps=Ntimesteps,Prepared=BiocontrolPrepared)
    n_after_biocontrol <- bc$N
    p_after_biocontrol <- INApestMLUTMPathogenReconcile(p_after_biology,n_after_biocontrol,Pathogen)
    bc_state <- bc$State; bc_impact <- bc$Impact; bc_prepared <- bc$Prepared
  }

  p_path_det <- INApestMLUTMPathogenDetectionProbability(p_after_biocontrol,Pathogen$DetectionProb)
  pathogen_detected <- rbinom(nn,1L,p_path_det)
  if(isTRUE(Pathogen$DetectionTriggersInfo)) {
    if(use_persistence) last_known_presence[pathogen_detected==1L] <- timestep
    have_info[have_info==0L] <- pathogen_detected[have_info==0L]
  }

  # Parent H+P engines place pathogen detection before host surveillance, so a
  # pathogen detection can activate info-triggered host surveillance immediately.
  info_before_surveillance <- as.integer(have_info!=0L)
  p_background <- INApestMLUTMHostDetectionProbability(n_after_biocontrol,detection_cube)
  background_detected <- rbinom(nn,1L,p_background)
  p_info <- INApestMLUTMHostDetectionProbability(n_after_biocontrol,info_detection_cube)
  if(any(info_detection_cube>0)) info_triggered_detected <- rbinom(nn,1L,p_info*info_before_surveillance) else info_triggered_detected <- integer(nn)
  host_detection_evidence <- pmax(background_detected,info_triggered_detected)
  if(use_persistence) last_known_presence[host_detection_evidence==1L] <- timestep
  have_info[have_info==0L] <- host_detection_evidence[have_info==0L]

  invaded_by_landuse <- matrix(0L,nn,L)
  for(i in seq_len(nn)) for(l in seq_len(L)) invaded_by_landuse[i,l] <- as.integer(sum(n_after_biocontrol[i,l,])>0)
  node_invaded_end <- as.integer(rowSums(invaded_by_landuse)>0)

  list(
    N=n_after_biocontrol, NAfterManagement=n_after_management, NAfterBiology=n_after_biology,
    ManagementDeaths=management_deaths, Managing=managing,
    PathogenState=p_after_biocontrol, PathogenStateAfterBiology=p_after_biology,
    PathogenDeaths=bio$PathogenDeaths, NewInfections=bio$NewInfections,
    PathogenDetectionProbability=p_path_det, PathogenDetected=pathogen_detected,
    HaveInfo=have_info, LastKnownPresence=last_known_presence,
    InformationStateBeforeSurveillance=info_before_surveillance,
    BackgroundDetectionProbability=p_background, InfoTriggeredDetectionProbability=p_info,
    BackgroundDetected=background_detected, InfoTriggeredDetected=info_triggered_detected,
    HostDetectionEvidence=host_detection_evidence, SEAMTransferred=seam_transferred,
    InvadedByLandUse=invaded_by_landuse, Invaded=node_invaded_end,
    KnownPresentByLandUse=invaded_by_landuse*have_info,
    BiocontrolState=bc_state, BiocontrolImpact=bc_impact, BiocontrolPrepared=bc_prepared
  )
}


### END HISTORICAL LAYER: INApestMetaTransitionMatrixMultipleLandUse_pathogen_biocontrol_v0.4.R
