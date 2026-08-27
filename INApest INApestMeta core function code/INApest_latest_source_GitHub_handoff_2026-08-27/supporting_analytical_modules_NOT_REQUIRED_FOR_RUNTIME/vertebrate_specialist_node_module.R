###############################################################################
### INApest vertebrate host analytical upgrade: specialist node mechanisms
### 27 August 2026
###
### This module is designed to be sourced after INApestAnalytical.R.  It adds
### a simulator-order deterministic expectation map for the specialist
### Vertebrate = list(Birth, HomeRange, Control, Interaction) architecture,
### plus local Jacobians, nonlinear finite-horizon growth diagnostics and an
### optional mate-limited birth helper / Allee profile.
###
### Analytical contract
### -------------------
### * Static linear mechanisms: exact first-moment map under the low-density
###   demographic/dispersal assumptions already used by the vertebrate
###   analytical companion.
### * State-dependent HomeRange/Control/Birth/Interaction or density-dependent
###   dispersal: plug-in deterministic mean map F(E[N]); generally not equal to
###   E[F(N)] and therefore labelled a mean-field closure.
### * Growth in nonlinear models is local: rho(J_F(x*)) for a one-step static
###   map, or rho(J_T ... J_1)^(1/T) for a time-varying trajectory.
### * Mate limitation can make the birth derivative at zero vanish.  Therefore
###   the low-density Jacobian and finite-abundance growth profile are both
###   reported; no single lambda is allowed to hide a strong Allee effect.
###############################################################################

.ivspec_clip01 <- function(x) {
  # Preserve matrix/array dimensions. pmin/pmax with a scalar first argument
  # can drop dimensions, which is unsafe for node x stage effects.
  y <- x
  y[y < 0] <- 0
  y[y > 1] <- 1
  y
}

.ivspec_call_hook <- function(fun, args, name = "vertebrate hook") {
  if (!is.function(fun)) stop(name, " must be a function")
  fm <- names(formals(fun))
  if (is.null(fm) || "..." %in% fm) return(do.call(fun, args))
  do.call(fun, args[intersect(names(args), fm)])
}

.ivspec_module <- function(Vertebrate, name) {
  if (is.null(Vertebrate) || is.null(Vertebrate[[name]])) return(NULL)
  Vertebrate[[name]]
}

.ivspec_node_time <- function(x, t, T, n, name, default = NULL, allow_inf = FALSE) {
  if (is.null(x)) {
    if (!is.null(default)) return(rep(default, n))
    stop(name, " may not be NULL")
  }
  if (is.function(x)) x <- .ivspec_call_hook(x, list(timestep = t), name)
  if (is.matrix(x)) {
    if (!all(dim(x) == c(n, T))) stop(name, " matrix must be nodes x Ntimesteps")
    x <- x[, t]
  }
  x <- as.numeric(x)
  if (length(x) == 1L) x <- rep(x, n)
  else if (length(x) == T && length(x) != n) x <- rep(x[t], n)
  else if (length(x) != n) stop(name, " must resolve to scalar, node vector, timestep vector, nodes x time matrix, or function")
  if (any(is.na(x)) || any(is.nan(x)) || (!allow_inf && any(!is.finite(x))))
    stop(name, if (allow_inf) " must not contain NA/NaN" else " must be finite")
  x
}

.ivspec_weights <- function(Weights, n, S) {
  if (is.null(Weights)) return(matrix(1, n, S))
  if (is.matrix(Weights)) {
    if (!all(dim(Weights) == c(n, S))) stop("Weights matrix must be nodes x stages")
    return(Weights)
  }
  z <- as.numeric(Weights)
  if (length(z) == 1L) return(matrix(z, n, S))
  if (length(z) == S) return(matrix(rep(z, each = n), n, S))
  if (length(z) == n && n != S) return(matrix(rep(z, S), n, S))
  stop("Weights must be scalar, length stages, length nodes (when unambiguous), or nodes x stages")
}

.ivspec_home_range_node <- function(module, population, timestep, context = list()) {
  n <- nrow(population); S <- ncol(population)
  if (is.null(module)) return(NULL)
  ans <- if (is.function(module)) {
    .ivspec_call_hook(module, list(
      population = population, state = population, timestep = timestep,
      perm = 1L, context = context
    ), "Vertebrate$HomeRange")
  } else module
  if (is.matrix(ans)) {
    if (!identical(dim(ans), c(n, S))) stop("Node HomeRange matrix must have dimensions nodes x stages")
    sigma <- ans
  } else {
    z <- as.numeric(ans)
    if (length(z) == 1L) sigma <- matrix(z, n, S)
    else if (length(z) == S && length(z) != n) sigma <- matrix(rep(z, each = n), n, S)
    else if (length(z) == n && length(z) != S) sigma <- matrix(rep(z, S), n, S)
    else if (length(z) == n && length(z) == S)
      stop("Node HomeRange vector is ambiguous because nodes equals stages; supply a nodes x stages matrix")
    else stop("Node HomeRange must be scalar, length nodes, length stages, or nodes x stages")
  }
  if (any(!is.finite(sigma)) || any(sigma < 0)) stop("Home-range sigma must be finite and non-negative")
  sigma
}

.ivspec_active_devices <- function(devices, timestep) {
  if (is.null(devices)) return(NULL)
  if (!is.data.frame(devices)) stop("Vertebrate$Control$Devices must be a data.frame")
  d <- as.data.frame(devices, stringsAsFactors = FALSE)
  if (!nrow(d)) return(d)
  start_name <- if ("start" %in% names(d)) "start" else if ("active_from" %in% names(d)) "active_from" else NULL
  end_name <- if ("end" %in% names(d)) "end" else if ("active_to" %in% names(d)) "active_to" else NULL
  keep <- rep(TRUE, nrow(d))
  if (!is.null(start_name)) keep <- keep & (is.na(d[[start_name]]) | d[[start_name]] <= timestep)
  if (!is.null(end_name)) keep <- keep & (is.na(d[[end_name]]) | d[[end_name]] >= timestep)
  d[keep, , drop = FALSE]
}

.ivspec_empty_effects <- function(n, S) list(
  kill_prob = matrix(0, n, S),
  detect_prob = matrix(0, n, S),
  fecundity_reduction = matrix(0, n, S),
  cost = 0
)

.ivspec_prob_matrix <- function(x, n, S, name, default = 0) {
  if (is.null(x)) return(matrix(default, n, S))
  if (is.matrix(x)) {
    if (!all(dim(x) == c(n, S))) stop(name, " must be nodes x stages")
    z <- x
  } else {
    v <- as.numeric(x)
    if (length(v) == 1L) z <- matrix(v, n, S)
    else if (length(v) == S && length(v) != n) z <- matrix(rep(v, each = n), n, S)
    else if (length(v) == n && length(v) != S) z <- matrix(rep(v, S), n, S)
    else if (length(v) == n && length(v) == S) stop(name, " vector is ambiguous because nodes equals stages; supply a matrix")
    else stop(name, " must resolve to scalar, node vector, stage vector, or nodes x stages matrix")
  }
  if (any(!is.finite(z)) || any(z < 0 | z > 1)) stop(name, " must resolve to probabilities in [0,1]")
  z
}

.ivspec_normalise_node_effects <- function(ans, population, name) {
  n <- nrow(population); S <- ncol(population)
  out <- .ivspec_empty_effects(n, S)
  if (is.null(ans)) return(out)
  if (is.numeric(ans) || is.matrix(ans)) {
    out$kill_prob <- .ivspec_prob_matrix(ans, n, S, paste0(name, " kill probability"))
    return(out)
  }
  if (!is.list(ans)) stop(name, " must return numeric/matrix kill probability or a list")
  if (!is.null(ans$Effects)) ans <- c(ans$Effects, list(cost = ans$Cost))
  out$kill_prob <- .ivspec_prob_matrix(ans$kill_prob %||% ans$Kill, n, S, paste0(name, "$kill_prob"), 0)
  out$detect_prob <- .ivspec_prob_matrix(ans$detect_prob %||% ans$Detect, n, S, paste0(name, "$detect_prob"), 0)
  out$fecundity_reduction <- .ivspec_prob_matrix(ans$fecundity_reduction %||% ans$FecundityReduction, n, S, paste0(name, "$fecundity_reduction"), 0)
  cc <- ans$cost %||% ans$Cost
  out$cost <- if (is.null(cc)) 0 else sum(as.numeric(cc))
  if (!is.finite(out$cost) || out$cost < 0) stop(name, " cost must be finite and non-negative")
  out
}

# base-R null coalescing helper local to this module
`%||%` <- function(a, b) if (!is.null(a)) a else b

.ivspec_combine_effects <- function(a, b) {
  list(
    kill_prob = 1 - (1 - a$kill_prob) * (1 - b$kill_prob),
    detect_prob = 1 - (1 - a$detect_prob) * (1 - b$detect_prob),
    fecundity_reduction = 1 - (1 - a$fecundity_reduction) * (1 - b$fecundity_reduction),
    cost = a$cost + b$cost
  )
}

.ivspec_node_device_default <- function(population, devices, home_range) {
  n <- nrow(population); S <- ncol(population)
  out <- .ivspec_empty_effects(n, S)
  if (is.null(devices) || !nrow(devices)) return(out)
  d <- as.data.frame(devices, stringsAsFactors = FALSE)
  if (!"node" %in% names(d)) stop("Default node-device control requires a node column")
  defaults <- list(density = 1, g0 = 1, effort = 1, kill = 1, detect = 0,
                   fecundity_reduction = 0, cost = 0)
  for (nm in names(defaults)) if (!nm %in% names(d)) d[[nm]] <- defaults[[nm]]
  if (any(is.na(d$node)) || any(d$node < 1 | d$node > n)) stop("Device node values must be 1..number of nodes")
  for (nm in c("g0", "kill", "detect", "fecundity_reduction")) {
    z <- as.numeric(d[[nm]])
    if (any(!is.finite(z)) || any(z < 0 | z > 1)) stop("Device ", nm, " values must be in [0,1]")
    d[[nm]] <- z
  }
  for (nm in c("density", "effort", "cost")) {
    z <- as.numeric(d[[nm]])
    if (any(!is.finite(z)) || any(z < 0)) stop("Device ", nm, " values must be finite and non-negative")
    d[[nm]] <- z
  }
  for (j in seq_len(nrow(d))) {
    node <- as.integer(d$node[j])
    sigma <- if ("sigma" %in% names(d) && is.finite(as.numeric(d$sigma[j]))) {
      rep(as.numeric(d$sigma[j]), S)
    } else {
      if (is.null(home_range)) stop("Default node-device control requires Vertebrate$HomeRange or device sigma")
      home_range[node, ]
    }
    encounter <- .ivspec_clip01(1 - exp(-2 * pi * d$g0[j] * sigma^2 * d$density[j] * d$effort[j]))
    out$kill_prob[node, ] <- 1 - (1 - out$kill_prob[node, ]) * (1 - encounter * d$kill[j])
    out$detect_prob[node, ] <- 1 - (1 - out$detect_prob[node, ]) * (1 - encounter * d$detect[j])
    out$fecundity_reduction[node, ] <- 1 - (1 - out$fecundity_reduction[node, ]) *
      (1 - encounter * d$fecundity_reduction[j])
  }
  out$cost <- sum(d$cost * d$effort)
  out
}

.ivspec_node_area_default <- function(population, area, timestep, context) {
  n <- nrow(population); S <- ncol(population)
  if (is.null(area)) return(.ivspec_empty_effects(n, S))
  if (is.function(area)) {
    ans <- .ivspec_call_hook(area, list(population = population, state = population,
      timestep = timestep, perm = 1L, context = context), "Vertebrate$Control$Area")
    return(.ivspec_normalise_node_effects(ans, population, "Vertebrate$Control$Area"))
  }
  if (is.numeric(area) || is.matrix(area))
    return(.ivspec_normalise_node_effects(area, population, "Vertebrate$Control$Area"))
  if (!is.list(area)) stop("Vertebrate$Control$Area must be numeric/matrix, function or list")
  resolve <- function(x, nm) {
    if (is.function(x)) x <- .ivspec_call_hook(x, list(population = population, state = population,
      timestep = timestep, perm = 1L, context = context), nm)
    .ivspec_prob_matrix(x, n, S, nm, 0)
  }
  out <- .ivspec_empty_effects(n, S)
  out$kill_prob <- resolve(area$Kill, "Control$Area$Kill")
  out$detect_prob <- resolve(area$Detect, "Control$Area$Detect")
  out$fecundity_reduction <- resolve(area$FecundityReduction, "Control$Area$FecundityReduction")
  if (!is.null(area$Cost)) {
    z <- if (is.function(area$Cost)) .ivspec_call_hook(area$Cost, list(population = population,
      state = population, timestep = timestep, perm = 1L, context = context), "Control$Area$Cost") else area$Cost
    out$cost <- sum(as.numeric(z))
    if (!is.finite(out$cost) || out$cost < 0) stop("Control$Area$Cost must be finite and non-negative")
  }
  out
}

.ivspec_node_control <- function(module, population, home_range, timestep, context = list()) {
  n <- nrow(population); S <- ncol(population)
  if (is.null(module)) return(.ivspec_empty_effects(n, S))
  if (!is.list(module)) stop("Vertebrate$Control must be NULL or a list")
  devices <- .ivspec_active_devices(module$Devices, timestep)
  dev <- if (is.function(module$Model)) {
    ans <- .ivspec_call_hook(module$Model, list(population = population, state = population,
      devices = devices, home_range = home_range, timestep = timestep, perm = 1L,
      context = context), "Vertebrate$Control$Model")
    .ivspec_normalise_node_effects(ans, population, "Vertebrate$Control$Model")
  } else .ivspec_node_device_default(population, devices, home_range)
  area <- .ivspec_node_area_default(population, module$Area, timestep, context)
  .ivspec_combine_effects(dev, area)
}

.ivspec_density_adjust_sdd <- function(SDD, population, K, Weights, alpha) {
  if (is.na(alpha) || alpha == 0) return(SDD)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0)
    stop("DispersalDensityFactor must be finite positive, or 0/NA")
  n <- nrow(SDD); S <- ncol(population)
  if (nrow(population) != n) stop("population/SDD node mismatch")
  # Mirrors the stochastic transition-matrix engine: reproductive dispersal
  # sees the same pre-transition census that produced the offspring.
  W <- .ivspec_weights(Weights, n, S)
  adult <- if (S >= 2L) rowSums(population[, 2:S, drop = FALSE] * W[, 2:S, drop = FALSE]) else rep(0, n)
  rel <- rep(1, n)
  pos <- K > 0
  rel[pos] <- pmin(1, pmax(0, adult[pos] / K[pos]))
  attract <- pmax(0, 1 - rel)^alpha
  export <- pmax(0, 1 - rowSums(SDD))
  choice <- SDD * rep(attract, each = n)
  den <- rowSums(choice) + export
  mult <- numeric(n); ok <- den > 0; mult[ok] <- 1 / den[ok]
  choice * mult
}

.ivspec_birth_movement <- function(SDDprob, LDDprob, LDDrate, population, K,
                                    Weights, DispersalDensityFactor, t, T) {
  n <- nrow(population)
  sdd <- .ivanal_resolve_time(SDDprob, t, T, "SDDprob")
  if (!all(dim(sdd) == c(n, n))) stop("SDDprob must resolve to nodes x nodes")
  if (!is.na(DispersalDensityFactor) && DispersalDensityFactor != 0)
    sdd <- .ivspec_density_adjust_sdd(sdd, population, K, Weights, DispersalDensityFactor)
  ldd <- .ivanal_resolve_time(LDDprob, t, T, "LDDprob", allow_null = TRUE)
  rr <- if (length(LDDrate) == T) LDDrate[t] else LDDrate
  .ivanal_mix_movement(sdd, ldd, rr, n)
}

.ivspec_call_birth <- function(module, population, transition, timestep, context) {
  if (is.null(module)) return(NULL)
  ans <- .ivspec_call_hook(module, list(population = population, state = population,
    transition = transition, timestep = timestep, perm = 1L, context = context),
    "Vertebrate$Birth")
  n <- nrow(population)
  if (is.null(ans)) return(list(mean = rep(0, n), mother_counts = NULL))
  if (is.numeric(ans) && !is.list(ans)) ans <- list(mean = ans)
  if (!is.list(ans) || is.null(ans$mean)) stop("Node Vertebrate$Birth must return numeric births or list(mean=..., mother_counts=...)")
  mean <- as.numeric(ans$mean)
  if (length(mean) == 1L) mean <- rep(mean, n)
  if (length(mean) != n || any(!is.finite(mean)) || any(mean < 0))
    stop("Node Vertebrate$Birth mean must resolve to one finite non-negative value per node")
  mothers <- ans$mother_counts
  if (!is.null(mothers)) {
    mothers <- as.numeric(mothers)
    if (length(mothers) == 1L) mothers <- rep(mothers, n)
    if (length(mothers) != n || any(!is.finite(mothers)) || any(mothers < 0))
      stop("Node Vertebrate$Birth mother_counts must resolve to one finite non-negative value per node")
  }
  list(mean = mean, mother_counts = mothers)
}

.ivspec_interaction_raw <- function(module, population, home_range, timestep, context = list(),
                                    floor_output = FALSE) {
  if (is.null(module)) return(population)
  # Optional exact linear analytical metadata.  Extra list elements are ignored
  # by the stochastic hook, so users can attach AnalyticalMatrix next to Update.
  if (is.list(module) && !is.null(module$AnalyticalMatrix)) {
    R <- module$AnalyticalMatrix
    if (is.function(R)) R <- .ivspec_call_hook(R, list(population = population,
      state = population, home_range = home_range, timestep = timestep,
      perm = 1L, context = context), "Interaction$AnalyticalMatrix")
    R <- as.matrix(R)
    K <- length(population)
    if (!all(dim(R) == c(K, K))) stop("Interaction$AnalyticalMatrix must be (nodes*stages) squared")
    z <- as.numeric(R %*% as.vector(t(population)))
    ans <- matrix(z, nrow(population), ncol(population), byrow = TRUE)
  } else {
    fun <- if (is.function(module)) module else if (is.list(module) && is.function(module$Update)) module$Update else NULL
    if (is.null(fun)) return(population)
    ans <- .ivspec_call_hook(fun, list(population = population, state = population,
      home_range = home_range, timestep = timestep, perm = 1L, context = context),
      "Vertebrate$Interaction")
    ans <- as.matrix(ans)
    if (!identical(dim(ans), dim(population))) stop("Vertebrate$Interaction must return nodes x stages")
  }
  if (any(!is.finite(ans)) || any(ans < 0)) stop("Vertebrate$Interaction returned non-finite or negative state")
  if (floor_output) ans <- floor(ans)
  ans
}

.ivspec_transition_existing <- function(N0, t, T, Nstages, Transition,
                                        TransitionSDDprob, TransitionLDDprob,
                                        TransitionLDDrate, OutsideEstablishmentProb) {
  n <- nrow(N0); S <- Nstages
  oe <- as.numeric(OutsideEstablishmentProb)
  if (length(oe) == 1L) oe <- rep(oe, n)
  if (length(oe) != n || any(!is.finite(oe)))
    stop("OutsideEstablishmentProb must resolve to scalar or node vector")
  oe <- .ivspec_clip01(oe)
  inside <- matrix(0, n, S)
  export <- 0
  trans <- .ivspec_transition_list(Transition, t, n, S, T)
  for (node in seq_len(n)) {
    A <- trans[[node]]
    if (any(!is.finite(A)) || any(A < 0)) stop("Transition entries must be finite and non-negative")
    for (stage in seq_len(S)) {
      x <- N0[node, stage]
      if (x <= 0) next
      if (stage == S) {
        inside[node, S] <- inside[node, S] + x * .ivspec_clip01(A[S, S])
      } else {
        p_stasis <- .ivspec_clip01(A[stage, stage])
        p_prog <- .ivspec_clip01(A[stage + 1L, stage])
        if (p_stasis + p_prog > 1 + 1e-10) {
          z <- p_stasis + p_prog
          p_stasis <- p_stasis / z
          p_prog <- p_prog / z
        }
        inside[node, stage] <- inside[node, stage] + x * p_stasis
        if (p_prog > 0) {
          Ps <- .ivanal_movement_for_stage(TransitionSDDprob, stage, t, T, n, S, "TransitionSDDprob")
          Pl <- .ivanal_movement_for_stage(TransitionLDDprob, stage, t, T, n, S, "TransitionLDDprob")
          trr <- if (length(TransitionLDDrate) == S - 1L) TransitionLDDrate[stage] else TransitionLDDrate
          Pm <- .ivanal_mix_movement(Ps, Pl, trr, n)
          if (is.null(Pm)) {
            inside[node, stage + 1L] <- inside[node, stage + 1L] + x * p_prog
          } else {
            pin <- pmax(0, Pm[node, ])
            if (sum(pin) > 1 + 1e-10) stop("Transition movement row sums may not exceed 1")
            inside[, stage + 1L] <- inside[, stage + 1L] + x * p_prog * pin
            export <- export + x * p_prog * pmax(0, 1 - sum(pin)) * oe[node]
          }
        }
      }
    }
  }
  list(inside = inside, export = export)
}

.ivspec_default_birth_mean <- function(N0, transition, response_fec_mult,
                                       control_fec_mult, n, S) {
  b <- numeric(n)
  if (S < 2L) return(b)
  for (i in seq_len(n)) {
    A <- transition[[i]]
    for (s in 2:S) {
      fec <- pmax(0, A[1, s])
      b[i] <- b[i] + N0[i, s] * fec * response_fec_mult[i, s] * control_fec_mult[i, s]
    }
  }
  b
}

.ivspec_disperse_birth_mean <- function(birth_mean, birth_move, estab,
                                        OutsideEstablishmentProb, n, S) {
  out <- matrix(0, n, S); export <- 0
  if (is.null(birth_move)) birth_move <- diag(n)
  oe <- as.numeric(OutsideEstablishmentProb)
  if (length(oe) == 1L) oe <- rep(oe, n)
  if (length(oe) != n) stop("OutsideEstablishmentProb must be scalar or node vector in specialist node solution")
  for (i in seq_len(n)) {
    if (birth_mean[i] <= 0) next
    pin <- pmax(0, birth_move[i, ])
    if (sum(pin) > 1 + 1e-10) stop("Birth movement row sums may not exceed 1")
    out[, 1L] <- out[, 1L] + birth_mean[i] * pin * estab
    export <- export + birth_mean[i] * pmax(0, 1 - sum(pin)) * .ivspec_clip01(oe[i])
  }
  list(inside = out, export = export)
}

.ivspec_transition_list <- function(Transition, t, n, S, T) {
  lapply(seq_len(n), function(i) .ivanal_transition_for_node(Transition, i, t, n, S, T))
}

.ivspec_transition_for_hook <- function(Transition, t, n, S, T) {
  trans <- .ivspec_transition_list(Transition, t, n, S, T)
  # Match the stochastic NodeTransition contract: homogeneous/time-sliced
  # transition input is a matrix; node-specific input is a list of matrices.
  if (is.list(Transition)) return(trans)
  trans[[1L]]
}

INApestVertebrateNodeMeanStep <- function(
    population,
    timestep = 1L,
    Ntimesteps = 1L,
    Nstages = ncol(population),
    Transition,
    SDDprob,
    LDDprob = NULL,
    LDDrate = 0,
    TransitionSDDprob = NULL,
    TransitionLDDprob = NULL,
    TransitionLDDrate = 0,
    PropaguleEstablishment = 1,
    EnvEstabProb = 1,
    MortalityProb = 0,
    ManagementExposure = 0,
    FecundityReduction = 0,
    OutsideEstablishmentProb = 1,
    Vertebrate = NULL,
    K = Inf,
    Weights = rep(1, Nstages),
    DispersalDensityFactor = 0,
    InteractionFloor = FALSE,
    ReturnComponents = FALSE) {

  N <- as.matrix(population)
  n <- nrow(N); S <- Nstages; t <- timestep; T <- Ntimesteps
  if (ncol(N) != S || any(!is.finite(N)) || any(N < 0)) stop("population must be finite non-negative nodes x stages")
  Knode <- .ivspec_node_time(K, t, T, n, "K", default = Inf, allow_inf = TRUE)
  W <- .ivspec_weights(Weights, n, S)
  trans <- .ivspec_transition_list(Transition, t, n, S, T)
  transition_hook <- .ivspec_transition_for_hook(Transition, t, n, S, T)
  HomeRangeModule <- .ivspec_module(Vertebrate, "HomeRange")
  ControlModule <- .ivspec_module(Vertebrate, "Control")
  BirthModule <- .ivspec_module(Vertebrate, "Birth")
  InteractionModule <- .ivspec_module(Vertebrate, "Interaction")

  # 1) Routine vertebrate control is independent of HaveInfo and precedes the
  #    ordinary response-management mortality in the stochastic engine.
  control_context <- list(transition = transition_hook, K = Knode, Nstages = S,
                          Weights = W, phase = "pre_response_control")
  hr_control <- .ivspec_home_range_node(HomeRangeModule, N, t, control_context)
  ce <- .ivspec_node_control(ControlModule, N, hr_control, t, control_context)
  N_control <- N * (1 - ce$kill_prob)

  # 2) Information-triggered INApest response mortality. ManagementExposure is
  #    the externally supplied probability/exposure envelope for this analytical
  #    host calculation. Shared information/adoption correlations are not hidden.
  response_mort <- .ivanal_expand_node_stage(MortalityProb, n, S, t, T, "MortalityProb", 0)
  manage <- .ivanal_expand_node_stage(ManagementExposure, n, S, t, T, "ManagementExposure", 0)
  response_mort <- .ivspec_clip01(response_mort)
  manage <- .ivspec_clip01(manage)
  # A node is either managed or not.  If management can kill animals before
  # reproduction, the managed fraction among survivors is not the original
  # adoption probability.  Conditioning on survival avoids the spurious
  # (1-a*m)(1-a*f) cross term and gives the exact first-moment mixture
  # (1-a) + a*(1-m)*(1-f) for linear fecundity.
  response_surv <- (1 - manage) + manage * (1 - response_mort)
  N0 <- N_control * response_surv
  manage_survivor <- matrix(0, n, S)
  nz_surv <- response_surv > 0
  manage_survivor[nz_surv] <- manage[nz_surv] * (1 - response_mort[nz_surv]) / response_surv[nz_surv]

  # 3) Birth hook is evaluated on the expected survivor census.  The context
  # reports the survivor-conditioned management probability because the
  # stochastic hook only sees animals that survived response mortality.
  response_fec <- .ivanal_expand_node_stage(FecundityReduction, n, S, t, T, "FecundityReduction", 0)
  response_fec <- .ivspec_clip01(response_fec)
  response_fec_mult <- 1 - response_fec * manage_survivor
  control_fec_mult <- 1 - ce$fecundity_reduction
  birth_context <- list(
    fecundity_multiplier = response_fec_mult * control_fec_mult,
    response_fecundity_reduction = response_fec,
    control_fecundity_reduction = ce$fecundity_reduction,
    managing = manage_survivor,
    management_probability_original = manage,
    management_probability_among_survivors = manage_survivor,
    home_range = .ivspec_home_range_node(HomeRangeModule, N0, t,
      list(transition = transition_hook, K = Knode, phase = "birth")),
    K = Knode, Nstages = S, Weights = W, phase = "birth")
  bspec <- .ivspec_call_birth(BirthModule, N0, transition_hook, t, birth_context)
  birth_mean <- if (is.null(bspec)) {
    .ivspec_default_birth_mean(N0, trans, response_fec_mult, control_fec_mult, n, S)
  } else bspec$mean

  # 4) Reproductive dispersal sees the same pre-transition N0 census. This also
  #    preserves the source engine's density-dependent dispersal timing.
  birth_move <- .ivspec_birth_movement(SDDprob, LDDprob, LDDrate, N0, Knode,
                                       W, DispersalDensityFactor, t, T)
  env <- .ivspec_node_time(EnvEstabProb, t, T, n, "EnvEstabProb")
  pe <- .ivspec_node_time(PropaguleEstablishment, t, T, n, "PropaguleEstablishment")
  outside_estab <- .ivspec_node_time(OutsideEstablishmentProb, t, T, n,
                                     "OutsideEstablishmentProb")
  estab <- .ivspec_clip01(env * pe)
  births <- .ivspec_disperse_birth_mean(birth_mean, birth_move, estab,
                                        outside_estab, n, S)

  # 5) Existing-animal stage survival/progression and optional transition
  #    movement are categorical and are applied once to N0.
  existing <- .ivspec_transition_existing(N0, t, T, S, Transition,
    TransitionSDDprob, TransitionLDDprob, TransitionLDDrate,
    outside_estab)
  NpreInteraction <- existing$inside + births$inside

  # 6) In the stochastic node engine Interaction is after local dynamics and
  #    external invasion. External invasion is outside the current intrinsic
  #    host estimands, so the analytical host map applies Interaction last.
  interaction_context <- list(transition = transition_hook, K = Knode, Nstages = S,
                              Weights = W, phase = "interaction")
  hr_interaction <- .ivspec_home_range_node(HomeRangeModule, NpreInteraction, t,
                                             interaction_context)
  Nnext <- .ivspec_interaction_raw(InteractionModule, NpreInteraction,
    hr_interaction, t, interaction_context, floor_output = InteractionFloor)

  # Expected node-level detection from routine control. In the simulator this is
  # registered after current-timestep information transfer, so it affects
  # response management no earlier than the next timestep.
  p_control_detect_node <- 1 - apply((1 - ce$detect_prob)^N, 1, prod)

  if (!ReturnComponents) return(Nnext)
  list(
    population = Nnext,
    routine_control = ce,
    expected_after_routine_control = N_control,
    expected_after_response_mortality = N0,
    expected_births_by_source_node = birth_mean,
    expected_births_inside = births$inside,
    expected_birth_export = births$export,
    expected_existing_inside = existing$inside,
    expected_existing_export = existing$export,
    expected_control_detection_probability = .ivspec_clip01(p_control_detect_node),
    expected_export = births$export + existing$export,
    birth_mother_counts = if (is.null(bspec)) NULL else bspec$mother_counts
  )
}

.ivspec_num_jacobian <- function(fun, x, eps = NULL, nonnegative = TRUE) {
  x <- as.numeric(x); p <- length(x)
  f0 <- as.numeric(fun(x)); m <- length(f0)
  J <- matrix(NA_real_, m, p)
  if (is.null(eps)) eps <- sqrt(.Machine$double.eps) * pmax(1, abs(x))
  if (length(eps) == 1L) eps <- rep(eps, p)
  if (length(eps) != p || any(!is.finite(eps)) || any(eps <= 0)) stop("Jacobian epsilon must be positive")
  for (j in seq_len(p)) {
    h <- eps[j]
    xp <- x; xp[j] <- xp[j] + h
    if (!nonnegative || x[j] > h) {
      xm <- x; xm[j] <- xm[j] - h
      J[, j] <- (as.numeric(fun(xp)) - as.numeric(fun(xm))) / (2 * h)
    } else {
      J[, j] <- (as.numeric(fun(xp)) - f0) / h
    }
  }
  J
}

.ivspec_static_inputs <- function(...) {
  xs <- list(...)
  !any(vapply(xs, function(x) {
    if (is.function(x)) return(TRUE)
    d <- dim(x)
    !is.null(d) && length(d) >= 3L
  }, logical(1)))
}

.ivspec_nonlinear_reasons <- function(Vertebrate, DispersalDensityFactor) {
  z <- character()
  if (!is.null(Vertebrate)) {
    if (is.function(Vertebrate$HomeRange)) z <- c(z, "state-dependent HomeRange hook")
    if (is.function(Vertebrate$Birth)) z <- c(z, "custom Birth hook")
    if (!is.null(Vertebrate$Control)) {
      C <- Vertebrate$Control
      if (is.function(C$Model)) z <- c(z, "custom Control$Model")
      if (is.function(C$Area) || (is.list(C$Area) && any(vapply(C$Area, is.function, logical(1)))))
        z <- c(z, "state-dependent Control$Area")
      if (is.function(Vertebrate$HomeRange) && !is.null(C$Devices))
        z <- c(z, "home-range-dependent device exposure")
    }
    if (!is.null(Vertebrate$Interaction)) {
      I <- Vertebrate$Interaction
      if (!(is.list(I) && !is.null(I$AnalyticalMatrix))) z <- c(z, "nonlinear/custom Interaction hook")
    }
  }
  if (!is.na(DispersalDensityFactor) && DispersalDensityFactor != 0)
    z <- c(z, "density-dependent reproductive dispersal")
  unique(z)
}

INApestVertebrateAnalyticalNodeSpecialist <- function(
    Ntimesteps,
    Nstages,
    Transition,
    InitialPopulation,
    SDDprob,
    LDDprob = NULL,
    LDDrate = 0,
    TransitionSDDprob = NULL,
    TransitionLDDprob = NULL,
    TransitionLDDrate = 0,
    PropaguleEstablishment = 1,
    EnvEstabProb = 1,
    MortalityProb = 0,
    ManagementExposure = 0,
    FecundityReduction = 0,
    OutsideEstablishmentProb = 1,
    Vertebrate = NULL,
    K = Inf,
    Weights = rep(1, Nstages),
    DispersalDensityFactor = 0,
    JacobianReference = c("zero", "initial", "final"),
    JacobianEpsilon = NULL,
    InteractionFloorTrajectory = FALSE,
    AlleeScales = NULL) {

  JacobianReference <- match.arg(JacobianReference)
  X0 <- as.matrix(InitialPopulation); n <- nrow(X0); S <- Nstages
  if (ncol(X0) != S) stop("InitialPopulation must be nodes x Nstages")
  p <- n * S
  stepfun <- function(x, tt, floor_interaction = InteractionFloorTrajectory, components = FALSE) {
    X <- matrix(x, n, S, byrow = TRUE)
    INApestVertebrateNodeMeanStep(
      population = X, timestep = tt, Ntimesteps = Ntimesteps, Nstages = S,
      Transition = Transition, SDDprob = SDDprob, LDDprob = LDDprob,
      LDDrate = LDDrate, TransitionSDDprob = TransitionSDDprob,
      TransitionLDDprob = TransitionLDDprob, TransitionLDDrate = TransitionLDDrate,
      PropaguleEstablishment = PropaguleEstablishment, EnvEstabProb = EnvEstabProb,
      MortalityProb = MortalityProb, ManagementExposure = ManagementExposure,
      FecundityReduction = FecundityReduction,
      OutsideEstablishmentProb = OutsideEstablishmentProb,
      Vertebrate = Vertebrate, K = K, Weights = Weights,
      DispersalDensityFactor = DispersalDensityFactor,
      InteractionFloor = floor_interaction, ReturnComponents = components)
  }

  traj <- matrix(NA_real_, Ntimesteps + 1L, p)
  traj[1, ] <- as.vector(t(X0)); exports <- numeric(Ntimesteps)
  control_detect <- matrix(NA_real_, Ntimesteps, n)
  control_cost <- numeric(Ntimesteps)
  for (tt in seq_len(Ntimesteps)) {
    comp <- stepfun(traj[tt, ], tt, components = TRUE)
    traj[tt + 1L, ] <- as.vector(t(comp$population))
    exports[tt] <- comp$expected_export
    control_detect[tt, ] <- comp$expected_control_detection_probability
    control_cost[tt] <- comp$routine_control$cost
  }

  ref <- switch(JacobianReference,
    zero = rep(0, p),
    initial = traj[1, ],
    final = traj[Ntimesteps + 1L, ])

  # For Jacobians use the smooth pre-floor Interaction map: floor() is a
  # stochastic implementation detail with derivative zero almost everywhere.
  Jsteps <- vector("list", Ntimesteps)
  jstate <- ref
  for (tt in seq_len(Ntimesteps)) {
    f <- function(x) as.vector(t(stepfun(x, tt, floor_interaction = FALSE)))
    Jsteps[[tt]] <- .ivspec_num_jacobian(f, jstate, JacobianEpsilon, TRUE)
    # A final-reference Jacobian is local to final state. For initial/zero we
    # retain the same reference through time; users can inspect each J_t.
  }
  if (Ntimesteps == 1L) {
    local_growth <- .ivanal_spectral_radius(Jsteps[[1L]])
    cycle_growth <- local_growth
  } else {
    P <- diag(p)
    for (tt in seq_len(Ntimesteps)) P <- Jsteps[[tt]] %*% P
    cycle_growth <- .ivanal_spectral_radius(P)
    local_growth <- cycle_growth^(1 / Ntimesteps)
  }

  nonlinear <- .ivspec_nonlinear_reasons(Vertebrate, DispersalDensityFactor)
  finite_mult <- if (sum(traj[1, ]) > 0) sum(traj[2, ]) / sum(traj[1, ]) else NA_real_

  # Mean-field first-passage export hazard. For a truly linear individual
  # branching model the older PGF solution remains preferable. With nonlinear
  # social/mating/home-range mechanisms this is deliberately labelled a hazard
  # approximation rather than a branching escape probability.
  escape_hazard <- .ivspec_clip01(1 - exp(-sum(pmax(0, exports))))

  allee <- NULL
  if (!is.null(AlleeScales)) {
    scales <- sort(unique(as.numeric(AlleeScales)))
    scales <- scales[is.finite(scales) & scales >= 0]
    if (!length(scales)) stop("AlleeScales must contain non-negative finite values")
    comp0 <- as.vector(t(X0))
    if (sum(comp0) <= 0) stop("AlleeScales requires a positive InitialPopulation composition")
    base_comp <- comp0 / sum(comp0)
    rows <- lapply(scales, function(a) {
      x <- a * base_comp
      y <- as.vector(t(stepfun(x, 1L, floor_interaction = FALSE)))
      data.frame(
        abundance = a,
        one_step_total = sum(y),
        one_step_multiplier = if (a > 0) sum(y) / a else NA_real_,
        net_change = sum(y) - a,
        stringsAsFactors = FALSE)
    })
    allee <- do.call(rbind, rows)
    # First sign-change interpolation for multiplier-1 / net change.
    d <- allee$net_change
    cross <- which(d[-length(d)] * d[-1L] <= 0 & allee$abundance[-length(d)] > 0)
    threshold <- NA_real_
    if (length(cross)) {
      k <- cross[1L]
      x1 <- allee$abundance[k]; x2 <- allee$abundance[k + 1L]
      y1 <- d[k]; y2 <- d[k + 1L]
      if (is.finite(y1) && is.finite(y2) && y2 != y1)
        threshold <- x1 - y1 * (x2 - x1) / (y2 - y1)
    }
    attr(allee, "ApproximateThreshold") <- threshold
  }

  diagnostics <- c(
    "Specialist vertebrate event order: routine control -> information-triggered response mortality -> birth/local demographic dynamics -> aggregate Interaction.",
    "Routine-control detections are reported as current-timestep detection probabilities but are not allowed to alter current response management; the simulator registers them for subsequent timesteps.",
    "Custom Birth is evaluated on the post-routine-control, post-response-mortality census and receives the same combined response/control fecundity multiplier as the stochastic hook.",
    "Reproductive density-dependent dispersal, when enabled, uses the pre-transition census and preserves original outside-row mass as an export option.",
    "Finite-horizon trajectory is a plug-in expectation map. It is an exact first moment only when the mechanism is linear/state-independent and low-density recruitment assumptions hold.",
    "Jacobian growth uses a smooth Interaction hook output (before simulator floor) because floor() has no useful infinitesimal derivative.",
    "For nonlinear mating/social/home-range mechanisms, local Jacobian growth and finite-abundance growth are both needed; lambda at zero alone is not a sufficient persistence diagnostic.",
    if (length(nonlinear)) paste0("Nonlinear/closure mechanisms detected: ", paste(nonlinear, collapse = "; "), ".") else
      "No obvious state-dependent specialist hook was detected; first-moment dynamics are linear under the supplied analytical inputs.",
    "EscapeHazardApproximation is based on cumulative expected successful export and is not an independent-lineage branching probability when social/mating/group dependence is active."
  )

  result <- list(
    Variant = "animal_node_specialist",
    GrowthRate = local_growth,
    LocalJacobianMultiplier = local_growth,
    JacobianReference = JacobianReference,
    EndPopulation = sum(traj[Ntimesteps + 1L, ]),
    EndPopulationByStage = colSums(matrix(traj[Ntimesteps + 1L, ], n, S, byrow = TRUE)),
    EscapeProbability = escape_hazard,
    EscapeHazardApproximation = escape_hazard,
    ExpectedTrajectory = rowSums(traj),
    ExpectedStateTrajectory = traj,
    OneStepFiniteAbundanceMultiplier = finite_mult,
    Jacobians = Jsteps,
    ExpectedExportByTimestep = exports,
    RoutineControlDetectionProbability = control_detect,
    RoutineControlCost = control_cost,
    NonlinearMechanisms = nonlinear,
    AlleeProfile = allee,
    Diagnostics = diagnostics)
  class(result) <- c("INApestVertebrateAnalyticalSpecialist", "list")
  result
}

# Helper that returns a stochastic-compatible node Birth function with explicit
# mate limitation. It can be passed directly as Vertebrate$Birth.
INApestVertebrateMateLimitedBirth <- function(
    FemaleStages,
    MaleStages,
    Fecundity,
    MateEncounter = c("poisson", "half_saturation", "presence"),
    MateEncounterRate = 1,
    HalfSaturation = 1,
    MaleWeights = NULL) {

  MateEncounter <- match.arg(MateEncounter)
  FemaleStages <- as.integer(FemaleStages); MaleStages <- as.integer(MaleStages)
  if (!length(FemaleStages) || !length(MaleStages) || any(FemaleStages < 1L) || any(MaleStages < 1L))
    stop("FemaleStages and MaleStages must be positive stage indices")
  fec <- as.numeric(Fecundity)
  if (length(fec) == 1L) fec <- rep(fec, length(FemaleStages))
  if (length(fec) != length(FemaleStages) || any(!is.finite(fec)) || any(fec < 0))
    stop("Fecundity must be scalar or one non-negative value per FemaleStage")
  mw <- if (is.null(MaleWeights)) rep(1, length(MaleStages)) else as.numeric(MaleWeights)
  if (length(mw) != length(MaleStages) || any(!is.finite(mw)) || any(mw < 0))
    stop("MaleWeights must be non-negative and match MaleStages")
  fun <- function(population, context = list(), ...) {
    X <- as.matrix(population); n <- nrow(X); S <- ncol(X)
    if (max(c(FemaleStages, MaleStages)) > S) stop("Mate-limited birth stage index exceeds population columns")
    mult <- context$fecundity_multiplier
    if (is.null(mult)) mult <- matrix(1, n, S)
    if (!all(dim(mult) == c(n, S))) stop("context$fecundity_multiplier must be nodes x stages")
    female_birth_base <- rowSums(sweep(X[, FemaleStages, drop = FALSE] * mult[, FemaleStages, drop = FALSE],
                                       2L, fec, `*`))
    males <- as.numeric(X[, MaleStages, drop = FALSE] %*% mw)
    g <- switch(MateEncounter,
      poisson = 1 - exp(-MateEncounterRate * males),
      half_saturation = males / (males + HalfSaturation),
      presence = as.numeric(males > 0))
    mean <- pmax(0, female_birth_base * g)
    mothers <- rowSums(X[, FemaleStages, drop = FALSE])
    list(mean = mean, mother_counts = mothers)
  }
  attr(fun, "INApestMateLimitedBirth") <- list(
    FemaleStages = FemaleStages, MaleStages = MaleStages,
    Fecundity = fec, MateEncounter = MateEncounter,
    MateEncounterRate = MateEncounterRate, HalfSaturation = HalfSaturation,
    MaleWeights = mw)
  class(fun) <- c("INApestVertebrateMateLimitedBirth", class(fun))
  fun
}

# Update existing animal analytical function without disturbing the old result
# for Vertebrate=NULL. Specialist hooks now trigger the simulator-order map.
INApestVertebrateAnalyticalAnimal_pre_specialist <- INApestVertebrateAnalyticalAnimal
INApestVertebrateAnalyticalAnimal <- function(
    Ntimesteps,
    Nstages,
    Transition,
    InitialPopulation,
    SDDprob,
    LDDprob = NULL,
    LDDrate = 0,
    TransitionSDDprob = NULL,
    TransitionLDDprob = NULL,
    TransitionLDDrate = 0,
    PropaguleEstablishment = 1,
    EnvEstabProb = 1,
    MortalityProb = 0,
    ManagementExposure = 0,
    FecundityReduction = 0,
    OutsideEstablishmentProb = 1,
    Vertebrate = NULL,
    K = Inf,
    Weights = rep(1, Nstages),
    DispersalDensityFactor = 0,
    JacobianReference = "zero",
    JacobianEpsilon = NULL,
    InteractionFloorTrajectory = FALSE,
    AlleeScales = NULL) {

  specialist <- !is.null(Vertebrate) && any(vapply(c("Birth", "HomeRange", "Control", "Interaction"),
    function(nm) !is.null(Vertebrate[[nm]]), logical(1)))
  specialist <- specialist || (!is.na(DispersalDensityFactor) && DispersalDensityFactor != 0)
  if (!specialist) return(INApestVertebrateAnalyticalAnimal_pre_specialist(
    Ntimesteps = Ntimesteps, Nstages = Nstages, Transition = Transition,
    InitialPopulation = InitialPopulation, SDDprob = SDDprob, LDDprob = LDDprob,
    LDDrate = LDDrate, TransitionSDDprob = TransitionSDDprob,
    TransitionLDDprob = TransitionLDDprob, TransitionLDDrate = TransitionLDDrate,
    PropaguleEstablishment = PropaguleEstablishment, EnvEstabProb = EnvEstabProb,
    MortalityProb = MortalityProb, ManagementExposure = ManagementExposure,
    FecundityReduction = FecundityReduction,
    OutsideEstablishmentProb = OutsideEstablishmentProb, Vertebrate = Vertebrate))

  INApestVertebrateAnalyticalNodeSpecialist(
    Ntimesteps = Ntimesteps, Nstages = Nstages, Transition = Transition,
    InitialPopulation = InitialPopulation, SDDprob = SDDprob, LDDprob = LDDprob,
    LDDrate = LDDrate, TransitionSDDprob = TransitionSDDprob,
    TransitionLDDprob = TransitionLDDprob, TransitionLDDrate = TransitionLDDrate,
    PropaguleEstablishment = PropaguleEstablishment, EnvEstabProb = EnvEstabProb,
    MortalityProb = MortalityProb, ManagementExposure = ManagementExposure,
    FecundityReduction = FecundityReduction,
    OutsideEstablishmentProb = OutsideEstablishmentProb, Vertebrate = Vertebrate,
    K = K, Weights = Weights, DispersalDensityFactor = DispersalDensityFactor,
    JacobianReference = JacobianReference, JacobianEpsilon = JacobianEpsilon,
    InteractionFloorTrajectory = InteractionFloorTrajectory, AlleeScales = AlleeScales)
}

# Final dispatcher refinement for specialist vertebrate-node arguments.
INApestAnalytical_pre_vertebrate_specialist <- INApestAnalytical
INApestAnalytical <- function(...) {
  args <- list(...)
  Model <- if (!is.null(args$Model)) as.character(args$Model)[1L] else "INApest"
  if (!identical(Model, "INApestVertebrateNode"))
    return(do.call(INApestAnalytical_pre_vertebrate_specialist, args))

  InitialState <- args$InitialState
  if (is.null(InitialState)) stop("INApestVertebrateNode requires InitialState")
  Transition <- args$Transition
  if (is.null(Transition)) stop("INApestVertebrateNode requires Transition")
  SDDprob <- args$SDDprob
  if (is.null(SDDprob)) stop("INApestVertebrateNode requires SDDprob")
  S <- args$Nstages %||% ncol(as.matrix(InitialState))
  LDD <- args$LDDprob
  if (!is.null(LDD) && length(LDD) == 1L && isTRUE(LDD == 0)) LDD <- NULL
  ans <- INApestVertebrateAnalyticalAnimal(
    Ntimesteps = args$Ntimesteps %||% 10L,
    Nstages = S,
    Transition = Transition,
    InitialPopulation = as.matrix(InitialState),
    SDDprob = SDDprob,
    LDDprob = LDD,
    LDDrate = args$LDDrate %||% 0,
    TransitionSDDprob = args$TransitionSDDprob,
    TransitionLDDprob = args$TransitionLDDprob,
    TransitionLDDrate = args$TransitionLDDrate %||% 0,
    PropaguleEstablishment = args$PropaguleEstablishment %||% 1,
    EnvEstabProb = args$EnvEstabProb %||% 1,
    MortalityProb = args$MortalityProb %||% 0,
    ManagementExposure = args$ManagementExposure %||% args$ManageProb %||% 0,
    FecundityReduction = args$FecundityReduction %||% 0,
    OutsideEstablishmentProb = args$OutsideEstablishmentProb %||% 1,
    Vertebrate = args$Vertebrate,
    K = args$K %||% Inf,
    Weights = args$Weights %||% rep(1, S),
    DispersalDensityFactor = args$DispersalDensityFactor %||% 0,
    JacobianReference = args$JacobianReference %||% "zero",
    JacobianEpsilon = args$JacobianEpsilon,
    InteractionFloorTrajectory = args$InteractionFloorTrajectory %||% FALSE,
    AlleeScales = args$AlleeScales)
  ans$Model <- Model
  ans$HeadlineEstimands <- list(
    GrowthRate = ans$GrowthRate,
    EndPopulation = ans$EndPopulation,
    EscapeProbability = ans$EscapeProbability)
  ans
}
