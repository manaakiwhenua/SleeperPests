###############################################################################
### INApest vertebrate host analytical upgrade: specialist point mechanisms
### 27 August 2026
###
### Mathematical scope
### ------------------
### * Exact current-geometry marginal exposure for the default point-device
###   control model.
### * Scalable point host dynamics are contracted to analysis-cell x stage x
###   persistent-state types. Existing animals retain persistent state through
###   stage movement; default matrix births enter an NA/default state stratum,
###   matching the stochastic engine (custom Birth is responsible for child
###   attributes).
### * Routine vertebrate control precedes response management. Response
###   management is integrated as the exact unmanaged/managed first-moment
###   mixture, not by multiplying marginal survival and fecundity reductions.
### * Arbitrary custom Birth and Interaction hooks do not determine their own
###   expectations. They are represented only when an explicit analytical mean
###   callback is supplied; otherwise the analytical function stops rather than
###   inventing an expectation.
### * Nonlinear hooks use a deterministic expectation map F(x) and numerical
###   Jacobian rho(dF/dx). This is a mean-field closure when the stochastic hook
###   depends on pair availability, local density, mate presence, or other
###   realised configuration properties.
### * Group-correlated branching is provided separately with the social group as
###   the branching unit; individual-lineage independence is not assumed.
###############################################################################

.ivps_clip01 <- function(x) {
  y <- x
  y[y < 0] <- 0
  y[y > 1] <- 1
  y
}

.ivps_call <- function(fun, args, name = "vertebrate point hook") {
  if (!is.function(fun)) stop(name, " must be a function")
  fm <- names(formals(fun))
  if (is.null(fm) || "..." %in% fm) return(do.call(fun, args))
  do.call(fun, args[intersect(names(args), fm)])
}

.ivps_module <- function(Vertebrate, name) {
  if (is.null(Vertebrate) || is.null(Vertebrate[[name]])) return(NULL)
  Vertebrate[[name]]
}

.ivps_point_frame <- function(points) {
  p <- as.data.frame(points, stringsAsFactors = FALSE)
  if (!all(c("x", "y") %in% names(p))) stop("Points must contain x and y")
  if (!"id" %in% names(p)) p$id <- seq_len(nrow(p))
  if (!"stage" %in% names(p)) p$stage <- 1L
  p$stage <- as.integer(p$stage)
  for (nm in c("parent_id", "birth_timestep", "last_known_timestep"))
    if (!nm %in% names(p)) p[[nm]] <- NA_integer_
  for (nm in c("have_info", "detected", "managing"))
    if (!nm %in% names(p)) p[[nm]] <- FALSE
  p
}

.ivps_home_range <- function(module, points, timestep, counts = NULL,
                             types = NULL, phase = "analytical") {
  if (is.null(module) || !nrow(points)) return(NULL)
  afun <- if (is.function(module)) attr(module, "INApestAnalyticalHomeRange") else NULL
  ans <- if (is.function(afun)) {
    .ivps_call(afun, list(points = points, state = points, counts = counts,
                          types = types, timestep = timestep, perm = 1L,
                          context = list(TypeCounts = counts, TypeTable = types,
                                         phase = phase)),
               "Vertebrate$HomeRange analytical callback")
  } else if (is.function(module)) {
    .ivps_call(module, list(points = points, state = points, timestep = timestep,
                            perm = 1L,
                            context = list(TypeCounts = counts, TypeTable = types,
                                           phase = phase)),
               "Vertebrate$HomeRange")
  } else module
  if (is.data.frame(ans)) {
    if (!all(c("id", "sigma") %in% names(ans)))
      stop("Point HomeRange data.frame output must contain id and sigma")
    sigma <- ans$sigma[match(points$id, ans$id)]
  } else {
    sigma <- as.numeric(ans)
    if (length(sigma) == 1L) sigma <- rep(sigma, nrow(points))
    if (length(sigma) != nrow(points))
      stop("Point HomeRange must resolve to scalar or one sigma per analytical type")
  }
  if (any(is.na(sigma)) || any(!is.finite(sigma)) || any(sigma < 0))
    stop("Home-range sigma values must be finite and non-negative")
  as.numeric(sigma)
}

.ivps_active_devices <- function(devices, timestep) {
  if (is.null(devices)) return(NULL)
  if (!is.data.frame(devices)) stop("Vertebrate$Control$Devices must be a data.frame")
  d <- as.data.frame(devices, stringsAsFactors = FALSE)
  if (!nrow(d)) return(d)
  sn <- if ("start" %in% names(d)) "start" else if ("active_from" %in% names(d)) "active_from" else NULL
  en <- if ("end" %in% names(d)) "end" else if ("active_to" %in% names(d)) "active_to" else NULL
  keep <- rep(TRUE, nrow(d))
  if (!is.null(sn)) keep <- keep & (is.na(d[[sn]]) | d[[sn]] <= timestep)
  if (!is.null(en)) keep <- keep & (is.na(d[[en]]) | d[[en]] >= timestep)
  d[keep, , drop = FALSE]
}

.ivps_empty_effects <- function(points) data.frame(
  id = points$id,
  kill_prob = rep(0, nrow(points)),
  detect_prob = rep(0, nrow(points)),
  fecundity_reduction = rep(0, nrow(points)),
  stringsAsFactors = FALSE)

.ivps_normalise_effects <- function(ans, points, name) {
  cost <- 0
  if (is.null(ans)) return(list(effects = .ivps_empty_effects(points), cost = 0))
  if (is.list(ans) && !is.data.frame(ans) && !is.null(ans$Effects)) {
    cost <- if (is.null(ans$Cost)) 0 else as.numeric(ans$Cost)[1L]
    ans <- ans$Effects
  }
  if (is.numeric(ans) && !is.data.frame(ans)) {
    z <- as.numeric(ans)
    if (length(z) == 1L) z <- rep(z, nrow(points))
    if (length(z) != nrow(points)) stop(name, " numeric output must be scalar or one kill probability per type")
    ans <- data.frame(id = points$id, kill_prob = z)
  }
  if (!is.data.frame(ans)) stop(name, " must return numeric, data.frame, or list(Effects=..., Cost=...)")
  if (!"id" %in% names(ans)) {
    if (nrow(ans) != nrow(points)) stop(name, " output without id must have one row per type")
    ans$id <- points$id
  }
  idx <- match(points$id, ans$id)
  out <- .ivps_empty_effects(points)
  for (nm in c("kill_prob", "detect_prob", "fecundity_reduction")) {
    if (nm %in% names(ans)) {
      v <- as.numeric(ans[[nm]][idx]); v[is.na(v)] <- 0
      out[[nm]] <- v
    }
    if (any(!is.finite(out[[nm]])) || any(out[[nm]] < 0 | out[[nm]] > 1))
      stop(name, " ", nm, " values must be in [0,1]")
  }
  if (!is.finite(cost) || cost < 0) stop(name, " Cost must be finite and non-negative")
  list(effects = out, cost = cost)
}

.ivps_combine_effects <- function(a, b) {
  out <- a
  out$kill_prob <- 1 - (1 - a$kill_prob) * (1 - b$kill_prob)
  out$detect_prob <- 1 - (1 - a$detect_prob) * (1 - b$detect_prob)
  out$fecundity_reduction <- 1 - (1 - a$fecundity_reduction) * (1 - b$fecundity_reduction)
  out
}

.ivps_device_effects <- function(points, devices, home_range) {
  out <- .ivps_empty_effects(points)
  if (is.null(devices) || !nrow(devices) || !nrow(points))
    return(list(effects = out, cost = 0))
  if (!all(c("x", "y") %in% names(devices)))
    stop("Default point-device control requires x and y columns")
  d <- devices
  defaults <- list(g0 = 1, effort = 1, kill = 1, detect = 0,
                   fecundity_reduction = 0, cost = 0)
  for (nm in names(defaults)) if (!nm %in% names(d)) d[[nm]] <- defaults[[nm]]
  for (nm in c("g0", "kill", "detect", "fecundity_reduction")) {
    d[[nm]] <- as.numeric(d[[nm]])
    if (any(!is.finite(d[[nm]])) || any(d[[nm]] < 0 | d[[nm]] > 1))
      stop("Device ", nm, " must be in [0,1]")
  }
  for (nm in c("effort", "cost")) {
    d[[nm]] <- as.numeric(d[[nm]])
    if (any(!is.finite(d[[nm]])) || any(d[[nm]] < 0))
      stop("Device ", nm, " must be finite and non-negative")
  }
  for (j in seq_len(nrow(d))) {
    sigma <- if ("sigma" %in% names(d) && is.finite(as.numeric(d$sigma[j]))) {
      rep(as.numeric(d$sigma[j]), nrow(points))
    } else {
      if (is.null(home_range)) stop("Default point-device control requires HomeRange or device sigma")
      home_range
    }
    d2 <- (points$x - as.numeric(d$x[j]))^2 + (points$y - as.numeric(d$y[j]))^2
    base <- numeric(nrow(points)); pos <- sigma > 0
    base[pos] <- d$g0[j] * exp(-d2[pos] / (2 * sigma[pos]^2))
    base[!pos & d2 == 0] <- d$g0[j]
    # Exact repeated-encounter semantics of the point engine.
    encounter <- 1 - (1 - .ivps_clip01(base))^d$effort[j]
    out$kill_prob <- 1 - (1 - out$kill_prob) * (1 - encounter * d$kill[j])
    out$detect_prob <- 1 - (1 - out$detect_prob) * (1 - encounter * d$detect[j])
    out$fecundity_reduction <- 1 - (1 - out$fecundity_reduction) *
      (1 - encounter * d$fecundity_reduction[j])
  }
  list(effects = out, cost = sum(d$cost * d$effort))
}

.ivps_area_effects <- function(points, area, timestep, counts, types) {
  if (is.null(area)) return(list(effects = .ivps_empty_effects(points), cost = 0))
  context <- list(TypeCounts = counts, TypeTable = types, phase = "pre_response_control")
  if (is.function(area)) {
    ans <- .ivps_call(area, list(points = points, state = points, timestep = timestep,
                                 perm = 1L, context = context),
                      "Vertebrate$Control$Area")
    return(.ivps_normalise_effects(ans, points, "Vertebrate$Control$Area"))
  }
  if (is.numeric(area)) return(.ivps_normalise_effects(area, points, "Vertebrate$Control$Area"))
  if (!is.list(area)) stop("Vertebrate$Control$Area must be numeric, function or list")
  out <- .ivps_empty_effects(points)
  resolve <- function(x, label) {
    if (is.null(x)) return(rep(0, nrow(points)))
    if (is.function(x)) x <- .ivps_call(x, list(points = points, state = points,
      timestep = timestep, perm = 1L, context = context), label)
    x <- as.numeric(x)
    if (length(x) == 1L) x <- rep(x, nrow(points))
    if (length(x) != nrow(points) || any(!is.finite(x)) || any(x < 0 | x > 1))
      stop(label, " must resolve to scalar or one [0,1] value per type")
    x
  }
  out$kill_prob <- resolve(area$Kill, "Control$Area$Kill")
  out$detect_prob <- resolve(area$Detect, "Control$Area$Detect")
  out$fecundity_reduction <- resolve(area$FecundityReduction, "Control$Area$FecundityReduction")
  cost <- if (is.null(area$Cost)) 0 else {
    z <- if (is.function(area$Cost)) .ivps_call(area$Cost, list(points = points,
      state = points, timestep = timestep, perm = 1L, context = context),
      "Control$Area$Cost") else area$Cost
    sum(as.numeric(z))
  }
  if (!is.finite(cost) || cost < 0) stop("Control$Area$Cost must be finite and non-negative")
  list(effects = out, cost = cost)
}

.ivps_control <- function(module, points, home_range, timestep, counts = NULL,
                          types = NULL) {
  if (is.null(module) || !nrow(points))
    return(list(effects = .ivps_empty_effects(points), cost = 0))
  if (!is.list(module)) stop("Vertebrate$Control must be NULL or a list")
  context <- list(TypeCounts = counts, TypeTable = types, phase = "pre_response_control")
  devices <- .ivps_active_devices(module$Devices, timestep)
  dev <- if (is.function(module$Model)) {
    ans <- .ivps_call(module$Model, list(points = points, state = points,
      devices = devices, home_range = home_range, timestep = timestep, perm = 1L,
      context = context), "Vertebrate$Control$Model")
    .ivps_normalise_effects(ans, points, "Vertebrate$Control$Model")
  } else .ivps_device_effects(points, devices, home_range)
  area <- .ivps_area_effects(points, module$Area, timestep, counts, types)
  list(effects = .ivps_combine_effects(dev$effects, area$effects),
       cost = dev$cost + area$cost)
}

# Public exact current-geometry exposure helper. No population averaging or grid
# contraction occurs here; each supplied activity centre is evaluated directly.
INApestVertebratePointControlExposure <- function(Points, Vertebrate,
                                                  timestep = 1L) {
  p <- .ivps_point_frame(Points)
  hr <- .ivps_home_range(.ivps_module(Vertebrate, "HomeRange"), p, timestep,
                         counts = rep(1, nrow(p)), types = p,
                         phase = "pre_response_control")
  z <- .ivps_control(.ivps_module(Vertebrate, "Control"), p, hr, timestep,
                     counts = rep(1, nrow(p)), types = p)
  out <- z$effects
  out$survival_prob <- 1 - out$kill_prob
  attr(out, "control_cost") <- z$cost
  out
}

.ivps_key <- function(df, cols) {
  if (!length(cols)) return(rep(".all", nrow(df)))
  zz <- lapply(cols, function(nm) {
    z <- df[[nm]]
    ifelse(is.na(z), "<NA>", paste0(typeof(z), ":", as.character(z)))
  })
  do.call(paste, c(zz, sep = "\r"))
}

.ivps_state_levels <- function(InitialPoints, StateColumns, StateLevels = NULL,
                               include_na = FALSE) {
  if (!length(StateColumns)) return(data.frame(.state_key = ".all", .stratum = 1L, stringsAsFactors = FALSE))
  miss <- setdiff(StateColumns, names(InitialPoints))
  if (length(miss)) stop("StateColumns missing from InitialPoints: ", paste(miss, collapse = ", "))
  x <- InitialPoints[, StateColumns, drop = FALSE]
  if (!is.null(StateLevels)) {
    sl <- as.data.frame(StateLevels, stringsAsFactors = FALSE)
    miss2 <- setdiff(StateColumns, names(sl))
    if (length(miss2)) stop("StateLevels must contain every StateColumns field")
    x <- rbind(x, sl[, StateColumns, drop = FALSE])
  }
  if (include_na) {
    na_row <- x[rep(1L, 1L), , drop = FALSE]
    for (nm in StateColumns) na_row[[nm]] <- NA
    x <- rbind(x, na_row)
  }
  key <- .ivps_key(x, StateColumns)
  keep <- !duplicated(key)
  out <- x[keep, , drop = FALSE]
  out$.state_key <- key[keep]
  out$.stratum <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}

.ivps_type_table <- function(InitialPoints, Nstages, analysis_grid,
                             StateColumns, StateLevels = NULL,
                             include_na = FALSE) {
  base <- .ina_pt_transition_representatives(InitialPoints, analysis_grid, Nstages)
  states <- .ivps_state_levels(InitialPoints, StateColumns, StateLevels, include_na)
  blocks <- lapply(seq_len(nrow(states)), function(k) {
    z <- base
    z$.base_type <- seq_len(nrow(base))
    z$.stratum <- k
    if (length(StateColumns)) {
      for (nm in StateColumns) z[[nm]] <- states[[nm]][k]
    }
    z
  })
  types <- do.call(rbind, blocks); rownames(types) <- NULL
  types$.type <- seq_len(nrow(types))
  types$id <- types$.type
  for (nm in c("parent_id", "birth_timestep", "last_known_timestep")) types[[nm]] <- NA_integer_
  for (nm in c("have_info", "detected", "managing")) types[[nm]] <- FALSE
  list(types = types, states = states, base = base)
}

.ivps_initial_counts <- function(InitialPoints, type_obj, analysis_grid,
                                 Nstages, StateColumns) {
  p <- .ivps_point_frame(InitialPoints)
  stage <- p$stage
  if (any(stage < 1L | stage > Nstages)) stop("InitialPoints$stage must be 1..Nstages")
  nc <- if (is.null(analysis_grid)) 1L else analysis_grid$nrow * analysis_grid$ncol
  cell <- if (is.null(analysis_grid)) rep(1L, nrow(p)) else .ina_pt_xy_to_cell(p$x, p$y, analysis_grid)
  sk <- .ivps_key(p, StateColumns)
  sm <- match(sk, type_obj$states$.state_key)
  bt <- (stage - 1L) * nc + cell
  idx <- (sm - 1L) * nrow(type_obj$base) + bt
  tabulate(idx[!is.na(idx)], nbins = nrow(type_obj$types))
}

.ivps_expand_prob <- function(x, types, timestep, Ntimesteps, name,
                              counts = NULL, default = 0) {
  if (is.null(x)) return(rep(default, nrow(types)))
  if (is.function(x)) {
    z <- .ivps_call(x, list(points = types, state = types, timestep = timestep,
      perm = 1L, counts = counts,
      context = list(TypeCounts = counts, TypeTable = types)), name)
  } else if (length(x) == 1L) z <- rep(as.numeric(x), nrow(types))
  else if (is.numeric(x) && length(x) == Ntimesteps) z <- rep(as.numeric(x[timestep]), nrow(types))
  else if (is.numeric(x) && length(x) == nrow(types)) z <- as.numeric(x)
  else stop(name, " must resolve to scalar, Ntimesteps vector, analytical-type vector, or function")
  if (length(z) == 1L) z <- rep(z, nrow(types))
  if (length(z) != nrow(types) || any(!is.finite(z))) stop(name, " must return one finite value per analytical type")
  .ivps_clip01(as.numeric(z))
}

.ivps_transition_export <- function(st, representatives, managed = FALSE) {
  A <- st$A; out <- numeric(nrow(representatives))
  for (i in seq_len(nrow(representatives))) {
    s <- representatives$stage[i]
    if (s < nrow(A)) {
      pr <- st$progress[[s]]; rr <- match(i, pr$rows)
      q <- if (managed) st$branch[[i]]$q1 else st$branch[[i]]$q0
      out[i] <- q * A[s + 1L, s] * pr$outside[rr]
    }
  }
  out
}

.ivps_base_step <- function(timestep, Ntimesteps, representatives, analysis_grid,
    Transition, Nstages, SDDkernel, LDDkernel, LDDrate,
    PropaguleEstablishment, EnvEstabProb,
    TransitionKernels, TransitionHabitatSearch, ApplyHabitatToTransitions,
    TransitionEstablishment, BlockedTransitionMortality,
    HabitatSuitability, HabitatSearchRadius, HabitatSearchCandidates,
    MortalityProb, MortalitySpatial, FecundityReduction,
    FecundityReductionSpatial, SpreadReduction, SpreadReductionSpatial,
    SpreadReductionAppliesTo, KernelSamples, PointSeed, managed) {

  .ina_pt_transition_operator_step(
    timestep = timestep, Ntimesteps = Ntimesteps,
    representatives = representatives, analysis_grid = analysis_grid,
    Transition = Transition, Nstages = Nstages,
    SDDkernel = SDDkernel, LDDkernel = LDDkernel, LDDrate = LDDrate,
    PropaguleEstablishment = PropaguleEstablishment, EnvEstabProb = EnvEstabProb,
    TransitionKernels = TransitionKernels,
    TransitionHabitatSearch = TransitionHabitatSearch,
    ApplyHabitatToTransitions = ApplyHabitatToTransitions,
    TransitionEstablishment = TransitionEstablishment,
    BlockedTransitionMortality = BlockedTransitionMortality,
    HabitatSuitability = HabitatSuitability,
    HabitatSearchRadius = HabitatSearchRadius,
    HabitatSearchCandidates = HabitatSearchCandidates,
    DetectionProb = 0, DetectionSpatial = NULL,
    ManageProb = if (managed) 1 else 0, ManageSpatial = NULL,
    MortalityProb = if (managed) MortalityProb else 0,
    MortalitySpatial = if (managed) MortalitySpatial else NULL,
    FecundityReduction = if (managed) FecundityReduction else 0,
    FecundityReductionSpatial = if (managed) FecundityReductionSpatial else NULL,
    SpreadReduction = if (managed) SpreadReduction else 0,
    SpreadReductionSpatial = if (managed) SpreadReductionSpatial else NULL,
    SpreadReductionAppliesTo = SpreadReductionAppliesTo,
    InfoRetentionProb = 1, InfoRadius = 0, InfoTransferProb = 0,
    InfoKernel = NULL, KernelSamples = KernelSamples,
    PointSeed = PointSeed, mode = if (managed) "all_informed" else "none",
    OutsideEstablishmentProb = 1)
}

.ivps_default_state_stratum <- function(type_obj, StateColumns, DefaultOffspringState) {
  if (!length(StateColumns)) return(1L)
  if (is.null(DefaultOffspringState)) {
    k <- .ivps_key(type_obj$states, StateColumns)
    target <- paste(rep("<NA>", length(StateColumns)), collapse = "\r")
    # Key encodes typeof for non-NA but literal <NA> for NA.
    j <- match(target, k)
    if (is.na(j)) stop("Internal error: NA offspring state stratum was not created")
    return(j)
  }
  d <- as.data.frame(DefaultOffspringState, stringsAsFactors = FALSE)
  if (nrow(d) != 1L || !all(StateColumns %in% names(d)))
    stop("DefaultOffspringState must be a one-row data.frame containing StateColumns")
  key <- .ivps_key(d[, StateColumns, drop = FALSE], StateColumns)
  j <- match(key, type_obj$states$.state_key)
  if (is.na(j)) stop("DefaultOffspringState is not represented by InitialPoints/StateLevels")
  j
}

.ivps_custom_birth_fun <- function(BirthModule) {
  if (is.null(BirthModule)) return(NULL)
  f <- attr(BirthModule, "INApestAnalyticalMean")
  if (is.null(f) && is.list(BirthModule) && is.function(BirthModule$AnalyticalMean)) f <- BirthModule$AnalyticalMean
  if (!is.function(f))
    stop("Custom point Vertebrate$Birth requires an explicit analytical expectation callback in attr(Birth, 'INApestAnalyticalMean'). The stochastic Birth hook returns realised offspring, so its expectation cannot be inferred safely.")
  f
}

.ivps_interaction_fun <- function(InteractionModule) {
  if (is.null(InteractionModule)) return(NULL)
  if (is.list(InteractionModule) && is.function(InteractionModule$AnalyticalMap))
    return(InteractionModule$AnalyticalMap)
  if (is.function(InteractionModule)) {
    f <- attr(InteractionModule, "INApestAnalyticalMap")
    if (is.function(f)) return(f)
  }
  if (is.list(InteractionModule) && !is.null(InteractionModule$AnalyticalMatrix)) {
    M <- as.matrix(InteractionModule$AnalyticalMatrix)
    return(function(counts, ...) as.numeric(M %*% counts))
  }
  stop("Point Vertebrate$Interaction requires AnalyticalMap or AnalyticalMatrix for analytical state dynamics; arbitrary realised Contact/Update hooks do not identify an expectation.")
}

.ivps_step_builder <- function(timestep, counts, type_obj, analysis_grid,
    Ntimesteps, Nstages, Transition, SDDkernel, LDDkernel, LDDrate,
    PropaguleEstablishment, EnvEstabProb, TransitionKernels,
    TransitionHabitatSearch, ApplyHabitatToTransitions, TransitionEstablishment,
    BlockedTransitionMortality, HabitatSuitability, HabitatSearchRadius,
    HabitatSearchCandidates, MortalityProb, MortalitySpatial,
    ManagementExposure, FecundityReduction, FecundityReductionSpatial,
    SpreadReduction, SpreadReductionSpatial, SpreadReductionAppliesTo,
    OutsideEstablishmentProb, Vertebrate, StateColumns, DefaultOffspringState,
    KernelSamples, PointSeed) {

  types <- type_obj$types; states <- type_obj$states; base <- type_obj$base
  nb <- nrow(base); ns <- nrow(states); K <- nrow(types)
  HomeRangeModule <- .ivps_module(Vertebrate, "HomeRange")
  ControlModule <- .ivps_module(Vertebrate, "Control")
  BirthModule <- .ivps_module(Vertebrate, "Birth")
  InteractionModule <- .ivps_module(Vertebrate, "Interaction")
  birth_fun <- .ivps_custom_birth_fun(BirthModule)
  interaction_fun <- .ivps_interaction_fun(InteractionModule)

  # Routine control is evaluated at analysis-type representative geometry.
  hr <- .ivps_home_range(HomeRangeModule, types, timestep, counts, types,
                         phase = "pre_response_control")
  ce <- .ivps_control(ControlModule, types, hr, timestep, counts, types)
  qctrl <- 1 - ce$effects$kill_prob
  cfr <- ce$effects$fecundity_reduction
  a <- .ivps_expand_prob(ManagementExposure, types, timestep, Ntimesteps,
                         "ManagementExposure", counts, 0)

  Parent <- matrix(0, K, K)
  Recruit <- matrix(0, K, K)
  export_transition <- numeric(K)
  export_recruit <- numeric(K)

  default_child_stratum <- .ivps_default_state_stratum(type_obj, StateColumns,
                                                        DefaultOffspringState)

  for (h in seq_len(ns)) {
    rows <- (h - 1L) * nb + seq_len(nb)
    reps <- types[rows, , drop = FALSE]
    st0 <- .ivps_base_step(timestep, Ntimesteps, reps, analysis_grid,
      Transition, Nstages, SDDkernel, LDDkernel, LDDrate,
      PropaguleEstablishment, EnvEstabProb, TransitionKernels,
      TransitionHabitatSearch, ApplyHabitatToTransitions, TransitionEstablishment,
      BlockedTransitionMortality, HabitatSuitability, HabitatSearchRadius,
      HabitatSearchCandidates, MortalityProb, MortalitySpatial,
      FecundityReduction, FecundityReductionSpatial, SpreadReduction,
      SpreadReductionSpatial, SpreadReductionAppliesTo, KernelSamples,
      PointSeed + h * 1000003L, FALSE)
    st1 <- .ivps_base_step(timestep, Ntimesteps, reps, analysis_grid,
      Transition, Nstages, SDDkernel, LDDkernel, LDDrate,
      PropaguleEstablishment, EnvEstabProb, TransitionKernels,
      TransitionHabitatSearch, ApplyHabitatToTransitions, TransitionEstablishment,
      BlockedTransitionMortality, HabitatSuitability, HabitatSearchRadius,
      HabitatSearchCandidates, MortalityProb, MortalitySpatial,
      FecundityReduction, FecundityReductionSpatial, SpreadReduction,
      SpreadReductionSpatial, SpreadReductionAppliesTo, KernelSamples,
      PointSeed + h * 1000003L, TRUE)

    # Existing animals retain state stratum. Routine control kills before both
    # response management and stage transition. Response management is an exact
    # unmanaged/managed mixture at the first-moment level.
    Pmix <- sweep(st0$Parent0, 2, 1 - a[rows], `*`) +
      sweep(st1$ParentH, 2, a[rows], `*`)
    Pmix <- sweep(Pmix, 2, qctrl[rows], `*`)
    Parent[rows, rows] <- Parent[rows, rows] + Pmix

    tr0 <- .ivps_transition_export(st0, reps, FALSE)
    tr1 <- .ivps_transition_export(st1, reps, TRUE)
    export_transition[rows] <- qctrl[rows] * ((1 - a[rows]) * tr0 + a[rows] * tr1)

    if (is.null(birth_fun)) {
      # Default stochastic births have no custom point-state attributes, so
      # persistent extra fields become NA in the new recruit. Map all default
      # recruits to the designated NA/default state stratum.
      Rmix <- sweep(st0$Recruit0, 2, (1 - a[rows]), `*`) +
        sweep(st1$RecruitH, 2, a[rows], `*`)
      Rmix <- sweep(Rmix, 2, qctrl[rows] * (1 - cfr[rows]), `*`)
      target <- (default_child_stratum - 1L) * nb + seq_len(nb)
      Recruit[target, rows] <- Recruit[target, rows] + Rmix

      rr0 <- pmax(0, st0$export - tr0)
      rr1 <- pmax(0, st1$export - tr1)
      pout <- if (length(OutsideEstablishmentProb) == 1L)
        rep(.ivps_clip01(as.numeric(OutsideEstablishmentProb)), nb)
      else {
        z <- as.numeric(OutsideEstablishmentProb)
        if (length(z) == K) z[rows] else if (length(z) == nb) z else
          stop("OutsideEstablishmentProb must be scalar, base-type, or full analytical-type length")
      }
      # st$export was built with outside establishment = 1. Apply the user's
      # outside establishment probability here after separating transition and
      # reproductive export.
      export_recruit[rows] <- qctrl[rows] * (1 - cfr[rows]) *
        ((1 - a[rows]) * rr0 + a[rows] * rr1) * pout
      export_transition[rows] <- export_transition[rows] * pout
    }
  }

  x0 <- as.numeric(counts)
  parent_next <- as.numeric(Parent %*% x0)
  expected_export <- sum(x0 * (export_transition + export_recruit))

  if (is.null(birth_fun)) {
    recruit_next <- as.numeric(Recruit %*% x0)
  } else {
    # Custom point Birth returns realised offspring in the stochastic engine;
    # analytical expectation is supplied explicitly. The callback returns an
    # absolute expected established-recruit vector and optional expected export.
    # It sees the expected post-routine-control/post-response-mortality census.
    # For nonlinear mating/pair processes this is a mean-field closure.
    # Survivor first moment by analytical type is read from the column sums of
    # the Parent carrier before stage redistribution; compute it directly from
    # routine control and response mortality probabilities for the context.
    # Response mortality may be stage/spatial/function-valued, so obtain the
    # conditional managed survival from a pair of unit-parent carrier sums.
    surv0 <- numeric(K); surv1 <- numeric(K)
    for (h in seq_len(ns)) {
      rows <- (h - 1L) * nb + seq_len(nb)
      reps <- types[rows, , drop = FALSE]
      st0 <- .ivps_base_step(timestep, Ntimesteps, reps, analysis_grid,
        Transition, Nstages, SDDkernel, LDDkernel, LDDrate,
        PropaguleEstablishment, EnvEstabProb, TransitionKernels,
        TransitionHabitatSearch, ApplyHabitatToTransitions, TransitionEstablishment,
        BlockedTransitionMortality, HabitatSuitability, HabitatSearchRadius,
        HabitatSearchCandidates, MortalityProb, MortalitySpatial,
        FecundityReduction, FecundityReductionSpatial, SpreadReduction,
        SpreadReductionSpatial, SpreadReductionAppliesTo, KernelSamples,
        PointSeed + h * 1000003L, FALSE)
      st1 <- .ivps_base_step(timestep, Ntimesteps, reps, analysis_grid,
        Transition, Nstages, SDDkernel, LDDkernel, LDDrate,
        PropaguleEstablishment, EnvEstabProb, TransitionKernels,
        TransitionHabitatSearch, ApplyHabitatToTransitions, TransitionEstablishment,
        BlockedTransitionMortality, HabitatSuitability, HabitatSearchRadius,
        HabitatSearchCandidates, MortalityProb, MortalitySpatial,
        FecundityReduction, FecundityReductionSpatial, SpreadReduction,
        SpreadReductionSpatial, SpreadReductionAppliesTo, KernelSamples,
        PointSeed + h * 1000003L, TRUE)
      # Existing carrier column sum excludes reproductive offspring and gives
      # survival/stage-transition retention inside the analysis domain. For
      # Birth timing we need mortality survival before stage transition, so use
      # branch q0/q1 rather than Parent column sums.
      surv0[rows] <- vapply(st0$branch, `[[`, numeric(1), "q0")
      surv1[rows] <- vapply(st1$branch, `[[`, numeric(1), "q1")
    }
    survivor_counts <- x0 * qctrl * ((1 - a) * surv0 + a * surv1)
    managed_survivor <- numeric(K)
    den <- (1 - a) * surv0 + a * surv1
    ok <- den > 0
    managed_survivor[ok] <- a[ok] * surv1[ok] / den[ok]
    response_fec <- .ina_pt_expected_parameter(
      FecundityReduction, FecundityReductionSpatial, types, timestep, 1L,
      Ntimesteps, "FecundityReduction")
    response_fec <- .ivps_clip01(as.numeric(response_fec))
    fec_mult_survivor <- (1 - response_fec * managed_survivor) * (1 - cfr)
    birth_home_range <- .ivps_home_range(
      HomeRangeModule, types, timestep, survivor_counts, types, phase = "birth")
    context <- list(
      TypeTable = types, TypeCountsBefore = x0,
      SurvivorCounts = survivor_counts,
      ManagementExposure = a,
      ManagementProbabilityAmongSurvivors = managed_survivor,
      ResponseFecundityReduction = response_fec,
      ControlFecundityReduction = cfr,
      FecundityMultiplier = fec_mult_survivor,
      fecundity_multiplier = fec_mult_survivor,
      ControlSurvival = qctrl,
      HomeRange = birth_home_range,
      home_range = birth_home_range,
      timestep = timestep,
      Approximation = "plug-in expectation on contracted point types")
    bz <- .ivps_call(birth_fun, list(
      counts = survivor_counts, types = types, timestep = timestep,
      context = context), "Vertebrate$Birth analytical expectation")
    if (is.numeric(bz)) bz <- list(recruits = bz, export = 0)
    if (!is.list(bz) || is.null(bz$recruits))
      stop("Birth analytical callback must return numeric recruits or list(recruits=..., export=...)")
    recruit_next <- as.numeric(bz$recruits)
    if (length(recruit_next) != K || any(!is.finite(recruit_next)) || any(recruit_next < 0))
      stop("Birth analytical recruits must be a non-negative vector with one value per analytical type")
    bx <- if (is.null(bz$export)) 0 else as.numeric(bz$export)[1L]
    if (!is.finite(bx) || bx < 0) stop("Birth analytical export must be one non-negative finite value")
    expected_export <- sum(x0 * export_transition) + bx
  }

  y <- parent_next + recruit_next
  if (!is.null(interaction_fun)) {
    interaction_home_range <- .ivps_home_range(
      HomeRangeModule, types, timestep, y, types, phase = "interaction")
    iz <- .ivps_call(interaction_fun, list(counts = y, types = types,
      timestep = timestep, home_range = interaction_home_range,
      context = list(TypeTable = types, TypeCountsBeforeInteraction = y,
                     HomeRange = interaction_home_range,
                     home_range = interaction_home_range, phase = "interaction")),
      "Vertebrate$Interaction analytical map")
    y <- as.numeric(iz)
    if (length(y) != K || any(!is.finite(y)) || any(y < 0))
      stop("Interaction analytical map must return one non-negative finite count per analytical type")
  }

  list(next_state = y, pre_interaction_state = parent_next + recruit_next,
       expected_recruits = recruit_next, export = expected_export, Parent = Parent,
       Recruit = if (is.null(birth_fun)) Recruit else NULL,
       export_transition_by_source = export_transition,
       export_recruit_by_source = export_recruit,
       control = ce, control_detection = ce$effects$detect_prob,
       type_table = types)
}

INApestVertebrateAnalyticalPointSpecialist <- function(
    Ntimesteps = 10,
    Nstages,
    Transition,
    InitialPoints,
    SDDkernel,
    LDDkernel = NULL,
    LDDrate = 0,
    PropaguleEstablishment = 1,
    EnvEstabProb = 1,
    TransitionKernels = NULL,
    TransitionHabitatSearch = FALSE,
    ApplyHabitatToTransitions = FALSE,
    TransitionEstablishment = 1,
    BlockedTransitionMortality = 0,
    LocalK = Inf,
    KRadius = 0,
    Weights = rep(1, Nstages),
    HabitatSuitability = NULL,
    HabitatSearchRadius = 0,
    HabitatSearchCandidates = 128,
    MortalityProb = 0,
    MortalitySpatial = NULL,
    ManagementExposure = 0,
    FecundityReduction = 0,
    FecundityReductionSpatial = NULL,
    SpreadReduction = 0,
    SpreadReductionSpatial = NULL,
    SpreadReductionAppliesTo = c("LDD", "all"),
    OutsideEstablishmentProb = 1,
    Vertebrate = NULL,
    StateColumns = NULL,
    StateLevels = NULL,
    DefaultOffspringState = NULL,
    PointAnalysisGrid = NULL,
    KernelSamples = 2000L,
    PointSeed = 1L,
    JacobianReference = c("zero", "initial", "trajectory"),
    JacobianEpsilon = NULL,
    AlleeScales = NULL,
    ReturnOperators = FALSE) {

  .ina_pt_require_point_helpers(TRUE)
  SpreadReductionAppliesTo <- match.arg(SpreadReductionAppliesTo)
  JacobianReference <- match.arg(JacobianReference)
  Ntimesteps <- as.integer(Ntimesteps); Nstages <- as.integer(Nstages)
  if (Ntimesteps < 1L || Nstages < 2L) stop("Ntimesteps >= 1 and Nstages >= 2 required")
  p <- .ivps_point_frame(InitialPoints)
  analysis_grid <- .ina_pt_analysis_grid(PointAnalysisGrid, HabitatSuitability,
                                          MortalitySpatial = MortalitySpatial,
                                          FecundityReductionSpatial = FecundityReductionSpatial,
                                          SpreadReductionSpatial = SpreadReductionSpatial)
  spatial_specialist <- !is.null(.ivps_module(Vertebrate, "Control")) ||
    !is.null(.ivps_module(Vertebrate, "HomeRange")) || length(StateColumns)
  if (spatial_specialist && is.null(analysis_grid) && nrow(p) > 1L) {
    # Homogeneous one-cell analysis remains valid for total population only if
    # exact activity-centre geometry is not needed. Explicit devices do need it.
    ctrl <- .ivps_module(Vertebrate, "Control")
    if (!is.null(ctrl) && !is.null(ctrl$Devices))
      stop("Explicit point devices require PointAnalysisGrid for future-population analytical dynamics. Use INApestVertebratePointControlExposure() for exact current-geometry exposure.")
  }
  include_na <- length(StateColumns) && is.null(DefaultOffspringState)
  to <- .ivps_type_table(p, Nstages, analysis_grid, StateColumns, StateLevels,
                         include_na = include_na)
  x0 <- .ivps_initial_counts(p, to, analysis_grid, Nstages, StateColumns)
  K <- length(x0)

  step_fun <- function(x, tt) {
    .ivps_step_builder(tt, x, to, analysis_grid, Ntimesteps, Nstages,
      Transition, SDDkernel, LDDkernel, LDDrate,
      PropaguleEstablishment, EnvEstabProb, TransitionKernels,
      TransitionHabitatSearch, ApplyHabitatToTransitions, TransitionEstablishment,
      BlockedTransitionMortality, HabitatSuitability, HabitatSearchRadius,
      HabitatSearchCandidates, MortalityProb, MortalitySpatial,
      ManagementExposure, FecundityReduction, FecundityReductionSpatial,
      SpreadReduction, SpreadReductionSpatial, SpreadReductionAppliesTo,
      OutsideEstablishmentProb, Vertebrate, StateColumns, DefaultOffspringState,
      KernelSamples, PointSeed)$next_state
  }

  traj <- matrix(0, Ntimesteps + 1L, K); traj[1L, ] <- x0
  exports <- numeric(Ntimesteps); control_cost <- numeric(Ntimesteps)
  control_detection <- vector("list", Ntimesteps)
  linear_ops <- vector("list", Ntimesteps)
  for (tt in seq_len(Ntimesteps)) {
    st <- .ivps_step_builder(tt, traj[tt, ], to, analysis_grid, Ntimesteps, Nstages,
      Transition, SDDkernel, LDDkernel, LDDrate,
      PropaguleEstablishment, EnvEstabProb, TransitionKernels,
      TransitionHabitatSearch, ApplyHabitatToTransitions, TransitionEstablishment,
      BlockedTransitionMortality, HabitatSuitability, HabitatSearchRadius,
      HabitatSearchCandidates, MortalityProb, MortalitySpatial,
      ManagementExposure, FecundityReduction, FecundityReductionSpatial,
      SpreadReduction, SpreadReductionSpatial, SpreadReductionAppliesTo,
      OutsideEstablishmentProb, Vertebrate, StateColumns, DefaultOffspringState,
      KernelSamples, PointSeed)
    traj[tt + 1L, ] <- st$next_state
    exports[tt] <- st$export
    control_cost[tt] <- st$control$cost
    control_detection[[tt]] <- st$control_detection
    if (!is.null(st$Recruit)) linear_ops[[tt]] <- st$Parent + st$Recruit
  }

  # Numerical Jacobians make the nonlinear definition explicit. For the fully
  # linear default-birth/no-interaction case, this also supplies a strong
  # identity check against the assembled first-moment operator.
  refs <- vector("list", Ntimesteps)
  for (tt in seq_len(Ntimesteps)) {
    refs[[tt]] <- switch(JacobianReference,
      zero = rep(0, K), initial = x0, trajectory = traj[tt, ])
  }
  J <- lapply(seq_len(Ntimesteps), function(tt)
    .ivspec_num_jacobian(function(z) step_fun(pmax(0, z), tt), refs[[tt]],
                          eps = JacobianEpsilon))
  cyc <- diag(K)
  for (tt in seq_len(Ntimesteps)) cyc <- J[[tt]] %*% cyc
  lambda <- max(Mod(eigen(cyc, only.values = TRUE)$values))^(1 / Ntimesteps)

  finite_mult <- if (sum(x0) > 0) sum(step_fun(x0, 1L)) / sum(x0) else NA_real_
  escape <- 1 - exp(-sum(pmax(0, exports)))

  allee <- NULL
  if (!is.null(AlleeScales)) {
    sc <- as.numeric(AlleeScales)
    allee <- data.frame(scale = sc, initial_total = NA_real_,
                        one_step_multiplier = NA_real_)
    for (i in seq_along(sc)) {
      z <- pmax(0, x0 * sc[i]); allee$initial_total[i] <- sum(z)
      allee$one_step_multiplier[i] <- if (sum(z) > 0) sum(step_fun(z, 1L)) / sum(z) else NA_real_
    }
  }

  nonlinear <- character(0)
  if (is.function(.ivps_module(Vertebrate, "HomeRange"))) nonlinear <- c(nonlinear, "HomeRange function")
  if (!is.null(.ivps_module(Vertebrate, "Birth"))) nonlinear <- c(nonlinear, "custom Birth expectation")
  if (!is.null(.ivps_module(Vertebrate, "Interaction"))) nonlinear <- c(nonlinear, "Interaction analytical map")
  if (is.function(ManagementExposure)) nonlinear <- c(nonlinear, "state-dependent management exposure")
  nonlinear <- unique(nonlinear)
  local_crowding_active <- !(length(LocalK) == 1L && is.numeric(LocalK) && is.infinite(LocalK)) || KRadius > 0

  diagnostics <- c(
    "Point specialist event order is routine vertebrate control -> response-management mixture -> birth -> stage survival/progression/movement -> analytical Interaction map.",
    "Default point-device exposure uses the exact fixed-device repeated-encounter formula at analytical type representative coordinates; refine PointAnalysisGrid when device exposure changes rapidly within cells.",
    if (length(StateColumns)) paste0("Persistent point state is stratified by: ", paste(StateColumns, collapse = ", "), ". Existing animals retain stratum through stage movement; default matrix offspring enter the configured default/NA stratum, matching the stochastic engine's lack of inherited extra fields without a custom Birth hook.") else "No extra persistent-state stratification was requested.",
    if (length(nonlinear)) paste0("Nonlinear/closure mechanisms: ", paste(nonlinear, collapse = ", "), ". GrowthRate is therefore a local Jacobian multiplier, not a global constant growth rate.") else "No explicit nonlinear specialist callback was detected; the contracted first-moment map is linear at fixed analytical geometry.",
    "Custom stochastic Birth and Contact/Update hooks are never reverse-engineered. Explicit analytical expectation callbacks are required because realised-hook code does not uniquely determine an expectation.",
    if (local_crowding_active) "Finite LocalK/KRadius uses random sequential spatial acceptance of recruits in the stochastic engine. That process is not silently replaced by an ordinary carrying-capacity term here; the rare/low-density growth operator remains valid, while finite-abundance trajectories require stochastic simulation once local crowding is material." else "Local recruit crowding is inactive (LocalK=Inf).",
    "EscapeProbability is a cumulative expected-export Poisson hazard for the specialist nonlinear map. Social/group dependence invalidates individual-lineage branching independence; use INApestVertebrateGroupBranching when groups are the natural branching unit."
  )

  result <- list(
    Variant = "animal_point_specialist",
    GrowthRate = lambda,
    LocalJacobianMultiplier = lambda,
    JacobianReference = JacobianReference,
    EndPopulation = sum(traj[Ntimesteps + 1L, ]),
    EndPopulationByType = traj[Ntimesteps + 1L, ],
    EscapeProbability = escape,
    EscapeHazardApproximation = escape,
    ExpectedTrajectory = rowSums(traj),
    ExpectedStateTrajectory = traj,
    TypeTable = to$types,
    StateLevels = to$states,
    OneStepFiniteAbundanceMultiplier = finite_mult,
    Jacobians = J,
    ExpectedExportByTimestep = exports,
    RoutineControlCost = control_cost,
    RoutineControlDetectionProbability = control_detection,
    NonlinearMechanisms = nonlinear,
    AlleeProfile = allee,
    LocalCrowdingAnalyticalStatus = if (local_crowding_active) "low-density only; finite crowding stochastic" else "inactive",
    Diagnostics = diagnostics,
    PointApproximation = list(
      SpatialMode = if (is.null(analysis_grid)) "homogeneous" else "analysis-grid contraction",
      AnalysisGrid = analysis_grid,
      KernelSamples = if (is.null(analysis_grid)) NA_integer_ else as.integer(KernelSamples)))
  if (ReturnOperators) result$LinearOperators <- linear_ops
  class(result) <- c("INApestVertebrateAnalyticalPointSpecialist", "list")
  result
}

###############################################################################
### Group-correlated branching
###############################################################################

# Treat a social group as the branching unit. MeanGroupOffspring[j,i] is the
# expected number of descendant groups of type j produced by one group of type
# i per timestep. This deliberately retains within-group dependence rather than
# pretending the animals inside one group are independent lineages.
INApestVertebrateGroupBranching <- function(
    MeanGroupOffspring,
    InitialGroups,
    EscapeProb = 0,
    Ntimesteps = 10L,
    GroupPGF = NULL,
    ExtinctionGenerations = 500L,
    tolerance = 1e-12) {

  M <- as.matrix(MeanGroupOffspring)
  if (nrow(M) != ncol(M) || any(!is.finite(M)) || any(M < 0))
    stop("MeanGroupOffspring must be a square non-negative finite matrix")
  K <- nrow(M); init <- as.numeric(InitialGroups)
  if (length(init) != K || any(!is.finite(init)) || any(init < 0) || any(init != floor(init)))
    stop("InitialGroups must be one non-negative whole-number count per group type")
  e <- as.numeric(EscapeProb); if (length(e) == 1L) e <- rep(e, K)
  if (length(e) != K || any(!is.finite(e)) || any(e < 0 | e > 1))
    stop("EscapeProb must be scalar or one [0,1] probability per group type")
  Ntimesteps <- as.integer(Ntimesteps)
  if (Ntimesteps < 1L) stop("Ntimesteps must be >= 1")

  pgf <- if (is.null(GroupPGF)) {
    function(z) exp(colSums(M * (z - 1)))
  } else {
    if (!is.function(GroupPGF)) stop("GroupPGF must be NULL or a function(z)")
    function(z) {
      q <- as.numeric(GroupPGF(z))
      if (length(q) != K || any(!is.finite(q)) || any(q < 0 | q > 1))
        stop("GroupPGF(z) must return K probabilities")
      q
    }
  }

  # Finite-horizon extinction (escape ignored).
  q <- rep(0, K)
  for (tt in seq_len(Ntimesteps)) q <- pgf(q)
  extinction_h <- exp(sum(init * log(pmax(q, .Machine$double.xmin))))

  # Eventual extinction minimal fixed point.
  qe <- rep(0, K)
  for (gg in seq_len(as.integer(ExtinctionGenerations))) {
    old <- qe; qe <- pgf(qe)
    if (max(abs(qe - old)) < tolerance) break
  }
  extinction_eventual <- exp(sum(init * log(pmax(qe, .Machine$double.xmin))))

  # No-escape recursion. A group escape event is correlated across all animals
  # represented by that group and therefore happens once at the group level.
  g <- rep(1, K)
  for (tt in seq_len(Ntimesteps)) g <- (1 - e) * pgf(g)
  noescape <- exp(sum(init * log(pmax(g, .Machine$double.xmin))))

  list(
    GrowthRate = max(Mod(eigen(M, only.values = TRUE)$values)),
    ExtinctionProbabilityByHorizon = extinction_h,
    EventualExtinctionProbability = extinction_eventual,
    EscapeProbabilityByHorizon = 1 - noescape,
    TypeExtinctionProbability = qe,
    TypeNoEscapeProbability = g,
    BranchingUnit = "social group",
    Method = if (is.null(GroupPGF))
      "multitype Poisson group branching" else "user-supplied multitype group PGF",
    Diagnostics = c(
      "Groups, not animals, are independent branching lineages in this calculation.",
      "Use group types that retain the composition/behaviour features responsible for correlated movement, reproduction or survival."))
}

###############################################################################
### User-facing wrappers / dispatcher refinement
###############################################################################

INApestVertebrateAnalyticalPoint_pre_specialist <- INApestVertebrateAnalyticalPoint
INApestVertebrateAnalyticalPoint <- function(
    Ntimesteps,
    Nstages,
    Transition,
    InitialStageCounts = NULL,
    ReproductiveMovement = NULL,
    TransitionMovement = NULL,
    PropaguleEstablishment = 1,
    EnvEstabProb = 1,
    MortalityProb = 0,
    ManagementExposure = 0,
    FecundityReduction = 0,
    OutsideEstablishmentProb = 1,
    Vertebrate = NULL,
    InitialPoints = NULL,
    SDDkernel = NULL,
    LDDkernel = NULL,
    LDDrate = 0,
    TransitionKernels = NULL,
    TransitionHabitatSearch = FALSE,
    ApplyHabitatToTransitions = FALSE,
    TransitionEstablishment = 1,
    BlockedTransitionMortality = 0,
    LocalK = Inf,
    KRadius = 0,
    Weights = rep(1, Nstages),
    HabitatSuitability = NULL,
    HabitatSearchRadius = 0,
    HabitatSearchCandidates = 128,
    MortalitySpatial = NULL,
    FecundityReductionSpatial = NULL,
    SpreadReduction = 0,
    SpreadReductionSpatial = NULL,
    SpreadReductionAppliesTo = c("LDD", "all"),
    StateColumns = NULL,
    StateLevels = NULL,
    DefaultOffspringState = NULL,
    PointAnalysisGrid = NULL,
    KernelSamples = 2000L,
    PointSeed = 1L,
    JacobianReference = "zero",
    JacobianEpsilon = NULL,
    AlleeScales = NULL,
    ReturnOperators = FALSE) {

  specialist <- !is.null(InitialPoints) ||
    (!is.null(Vertebrate) && any(vapply(c("Birth", "HomeRange", "Control", "Interaction"),
      function(nm) !is.null(Vertebrate[[nm]]), logical(1)))) || length(StateColumns)
  if (!specialist) {
    return(INApestVertebrateAnalyticalPoint_pre_specialist(
      Ntimesteps = Ntimesteps, Nstages = Nstages, Transition = Transition,
      InitialStageCounts = InitialStageCounts,
      ReproductiveMovement = ReproductiveMovement,
      TransitionMovement = TransitionMovement,
      PropaguleEstablishment = PropaguleEstablishment,
      EnvEstabProb = EnvEstabProb, MortalityProb = MortalityProb,
      ManagementExposure = ManagementExposure,
      FecundityReduction = FecundityReduction,
      OutsideEstablishmentProb = OutsideEstablishmentProb,
      Vertebrate = Vertebrate))
  }
  if (is.null(InitialPoints)) stop("Specialist point analytical mode requires InitialPoints")
  if (is.null(SDDkernel)) stop("Specialist point analytical mode requires SDDkernel")
  INApestVertebrateAnalyticalPointSpecialist(
    Ntimesteps = Ntimesteps, Nstages = Nstages, Transition = Transition,
    InitialPoints = InitialPoints, SDDkernel = SDDkernel,
    LDDkernel = LDDkernel, LDDrate = LDDrate,
    PropaguleEstablishment = PropaguleEstablishment,
    EnvEstabProb = EnvEstabProb, TransitionKernels = TransitionKernels,
    TransitionHabitatSearch = TransitionHabitatSearch,
    ApplyHabitatToTransitions = ApplyHabitatToTransitions,
    TransitionEstablishment = TransitionEstablishment,
    BlockedTransitionMortality = BlockedTransitionMortality,
    LocalK = LocalK, KRadius = KRadius, Weights = Weights,
    HabitatSuitability = HabitatSuitability,
    HabitatSearchRadius = HabitatSearchRadius,
    HabitatSearchCandidates = HabitatSearchCandidates,
    MortalityProb = MortalityProb, MortalitySpatial = MortalitySpatial,
    ManagementExposure = ManagementExposure,
    FecundityReduction = FecundityReduction,
    FecundityReductionSpatial = FecundityReductionSpatial,
    SpreadReduction = SpreadReduction,
    SpreadReductionSpatial = SpreadReductionSpatial,
    SpreadReductionAppliesTo = SpreadReductionAppliesTo,
    OutsideEstablishmentProb = OutsideEstablishmentProb,
    Vertebrate = Vertebrate, StateColumns = StateColumns,
    StateLevels = StateLevels, DefaultOffspringState = DefaultOffspringState,
    PointAnalysisGrid = PointAnalysisGrid, KernelSamples = KernelSamples,
    PointSeed = PointSeed, JacobianReference = JacobianReference,
    JacobianEpsilon = JacobianEpsilon, AlleeScales = AlleeScales,
    ReturnOperators = ReturnOperators)
}

INApestAnalytical_pre_vertebrate_point_specialist <- INApestAnalytical
INApestAnalytical <- function(...) {
  args <- list(...)
  Model <- if (!is.null(args$Model)) as.character(args$Model)[1L] else "INApest"
  if (!identical(Model, "INApestVertebratePoint"))
    return(do.call(INApestAnalytical_pre_vertebrate_point_specialist, args))

  InitialPoints <- args$InitialPoints
  if (!is.null(InitialPoints)) {
    S <- args$Nstages
    if (is.null(S)) S <- max(as.integer(InitialPoints$stage %||% 1L))
    ans <- INApestVertebrateAnalyticalPoint(
      Ntimesteps = args$Ntimesteps %||% 10L,
      Nstages = S,
      Transition = args$Transition,
      InitialPoints = InitialPoints,
      SDDkernel = args$SDDkernel,
      LDDkernel = args$LDDkernel,
      LDDrate = args$LDDrate %||% 0,
      PropaguleEstablishment = args$PropaguleEstablishment %||% 1,
      EnvEstabProb = args$EnvEstabProb %||% 1,
      TransitionKernels = args$TransitionKernels,
      TransitionHabitatSearch = args$TransitionHabitatSearch %||% FALSE,
      ApplyHabitatToTransitions = args$ApplyHabitatToTransitions %||% FALSE,
      TransitionEstablishment = args$TransitionEstablishment %||% 1,
      BlockedTransitionMortality = args$BlockedTransitionMortality %||% 0,
      LocalK = args$LocalK %||% Inf, KRadius = args$KRadius %||% 0,
      Weights = args$Weights %||% rep(1, S),
      HabitatSuitability = args$HabitatSuitability,
      HabitatSearchRadius = args$HabitatSearchRadius %||% 0,
      HabitatSearchCandidates = args$HabitatSearchCandidates %||% 128L,
      MortalityProb = args$MortalityProb %||% 0,
      MortalitySpatial = args$MortalitySpatial,
      ManagementExposure = args$ManagementExposure %||% args$ManageProb %||% 0,
      FecundityReduction = args$FecundityReduction %||% 0,
      FecundityReductionSpatial = args$FecundityReductionSpatial,
      SpreadReduction = args$SpreadReduction %||% 0,
      SpreadReductionSpatial = args$SpreadReductionSpatial,
      SpreadReductionAppliesTo = args$SpreadReductionAppliesTo %||% "LDD",
      OutsideEstablishmentProb = args$OutsideEstablishmentProb %||% 1,
      Vertebrate = args$Vertebrate,
      StateColumns = args$StateColumns,
      StateLevels = args$StateLevels,
      DefaultOffspringState = args$DefaultOffspringState,
      PointAnalysisGrid = args$PointAnalysisGrid,
      KernelSamples = args$KernelSamples %||% 2000L,
      PointSeed = args$PointSeed %||% 1L,
      JacobianReference = args$JacobianReference %||% "zero",
      JacobianEpsilon = args$JacobianEpsilon,
      AlleeScales = args$AlleeScales,
      ReturnOperators = args$ReturnOperators %||% FALSE)
  } else {
    # Preserve the previous contracted-matrix API exactly.
    keep <- setdiff(names(args), c("Model"))
    ans <- do.call(INApestVertebrateAnalyticalPoint,
                   args[intersect(keep, names(formals(INApestVertebrateAnalyticalPoint)))])
  }
  ans$Model <- Model
  ans$HeadlineEstimands <- list(GrowthRate = ans$GrowthRate,
                                EndPopulation = ans$EndPopulation,
                                EscapeProbability = ans$EscapeProbability)
  ans
}
