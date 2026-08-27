###############################################################################
### Vertebrate + pathogen analytical extension
###
### Rare-pathogen and finite-background calculations that compose the validated
### vertebrate host carrier with the validated INApest SIS/SIR/SEIR disease
### process. Disease remains an Interaction process; no Vertebrate$Disease API
### is introduced.
###############################################################################

.ivp_spec <- function(Pathogen) {
  if (inherits(Pathogen, "INApestPointPathogenInteraction")) return(Pathogen$Pathogen)
  if (inherits(Pathogen, "INApestPathogen")) return(Pathogen)
  stop("Pathogen must be an INApestPathogen or INApestPointPathogenInteraction object")
}

.ivp_rho <- function(M) {
  M <- as.matrix(M)
  if (!nrow(M)) return(0)
  max(Mod(eigen(M, only.values = TRUE)$values))
}

.ivp_class <- function(x, tol = 1e-12) {
  if (x > 1 + tol) "growing" else if (x < 1 - tol) "declining" else "threshold"
}

.ivp_ngm <- function(G, T0) {
  G <- as.matrix(G); T0 <- as.matrix(T0)
  if (!identical(dim(G), dim(T0))) stop("G and T0 dimensions differ")
  F <- G - T0
  # Numerical construction can produce epsilon-scale negative entries.
  F[abs(F) < 1e-13] <- 0
  rt <- .ivp_rho(T0)
  if (rt >= 1 - 1e-12) {
    K <- matrix(NA_real_, nrow(G), ncol(G)); r0 <- Inf
  } else {
    K <- F %*% solve(diag(nrow(G)) - T0)
    r0 <- .ivp_rho(K)
  }
  list(TransitionWithoutNewInfection = T0,
       NewInfectionOperator = F,
       NextGenerationOperator = K,
       R0 = r0,
       RhoTransition = rt)
}

.ivp_ordered_growth <- function(ops) {
  if (!is.list(ops) || !length(ops)) stop("ops must be a non-empty list")
  d <- dim(ops[[1L]])
  P <- diag(d[1L])
  for (G in ops) {
    if (!all(dim(G) == d)) stop("all operators must have the same dimensions")
    P <- G %*% P
  }
  rho <- .ivp_rho(P)
  list(CycleOperator = P, CycleMultiplier = rho,
       PerTimestepMultiplier = rho^(1 / length(ops)),
       Classification = .ivp_class(rho^(1 / length(ops))))
}

.ivp_poisson_branching <- function(Operators, InitialActive,
    EscapeMean = NULL, Eventual = TRUE, tolerance = 1e-12, maxiter = 100000L) {
  ops <- lapply(Operators, as.matrix)
  K <- nrow(ops[[1L]])
  init <- as.numeric(InitialActive)
  if (length(init) != K || any(!is.finite(init)) || any(init < 0))
    stop("InitialActive must be one non-negative finite value per active type")
  if (is.null(EscapeMean)) EscapeMean <- lapply(ops, function(M) rep(0, K))
  if (!is.list(EscapeMean)) EscapeMean <- rep(list(as.numeric(EscapeMean)), length(ops))
  if (length(EscapeMean) != length(ops)) stop("EscapeMean length must match Operators")
  EscapeMean <- lapply(EscapeMean, function(e) {
    e <- as.numeric(e); if (length(e) == 1L) e <- rep(e, K)
    if (length(e) != K || any(!is.finite(e)) || any(e < 0)) stop("invalid EscapeMean")
    e
  })
  pgf <- function(M, z) exp(colSums(M * (z - 1)))

  # Extinction by the finite horizon, with chronological PGFs nested backward.
  q <- rep(0, K)
  for (tt in rev(seq_along(ops))) q <- pgf(ops[[tt]], q)
  pext_h <- exp(sum(init * log(pmax(q, .Machine$double.xmin))))

  # Probability of no escape by horizon: terminal future is safe (g=1), then
  # include a Poisson active-export event at each source type and descendant PGF.
  g <- rep(1, K)
  for (tt in rev(seq_along(ops)))
    g <- exp(-EscapeMean[[tt]] + colSums(ops[[tt]] * (g - 1)))
  p_noescape <- exp(sum(init * log(pmax(g, .Machine$double.xmin))))

  pevent <- NA_real_; qevent <- rep(NA_real_, K)
  if (isTRUE(Eventual) && length(ops) == 1L) {
    qevent <- rep(0, K); M <- ops[[1L]]
    for (ii in seq_len(maxiter)) {
      qn <- pgf(M, qevent)
      if (max(abs(qn - qevent)) < tolerance) { qevent <- qn; break }
      qevent <- qn
    }
    pevent <- exp(sum(init * log(pmax(qevent, .Machine$double.xmin))))
  }
  list(ExtinctionProbabilityByHorizon = pext_h,
       EventualExtinctionProbability = pevent,
       TypeEventualExtinctionProbability = qevent,
       EscapeProbabilityByHorizon = 1 - p_noescape,
       TypeNoEscapeProbabilityByHorizon = g,
       Method = "multitype Poisson branching matched to the rare-pathogen mean operator")
}

.ivp_vertical_matrix <- function(R, p, K, name = "VerticalTransmissionProb") {
  if (is.null(R)) return(NULL)
  R <- as.matrix(R)
  if (!all(dim(R) == c(K, K))) stop("recruit operator has wrong dimensions")
  if (is.function(p)) stop(name, " resolver functions require an explicit VerticalRecruitOperator")
  if (is.matrix(p)) {
    if (!all(dim(p) == c(K, K))) stop(name, " matrix must be analytical types x analytical types")
    V <- R * p
  } else {
    z <- as.numeric(p)
    if (length(z) == 1L) V <- R * z
    else if (length(z) == K) V <- sweep(R, 2, z, `*`)
    else stop(name, " must be scalar, source-type vector, or target x source matrix")
  }
  if (any(!is.finite(V)) || any(V < 0) || any(V - R > 1e-12))
    stop(name, " must imply probabilities in [0,1]")
  V
}

.ivp_active_carrier <- function(Model, H, V = NULL,
    VerticalSourceState = "I", VerticalTargetState = NULL) {
  Model <- match.arg(Model, c("SIS", "SIR", "SEIR")); H <- as.matrix(H); K <- nrow(H)
  if (ncol(H) != K) stop("H must be square")
  if (is.null(VerticalTargetState)) VerticalTargetState <- if (Model == "SEIR") "E" else "I"
  if (Model %in% c("SIS", "SIR")) {
    if (!VerticalSourceState %in% "I" || !VerticalTargetState %in% "I")
      stop("SIS/SIR vertical transmission can only map I mothers to I offspring")
    C <- H
    if (!is.null(V)) C <- C + V
    return(list(PreDiseaseCarrier = C, NoNewInfectionCarrier = H,
                ActiveStates = "I", VerticalTargetState = VerticalTargetState))
  }
  states <- c("E", "I")
  if (!VerticalSourceState %in% states || !VerticalTargetState %in% states)
    stop("SEIR VerticalSourceState/VerticalTargetState must be E or I")
  Z <- matrix(0, K, K)
  C <- rbind(cbind(H, Z), cbind(Z, H))
  H0 <- C
  if (!is.null(V)) {
    rs <- if (VerticalTargetState == "E") seq_len(K) else K + seq_len(K)
    cs <- if (VerticalSourceState == "E") seq_len(K) else K + seq_len(K)
    C[rs, cs] <- C[rs, cs] + V
  }
  list(PreDiseaseCarrier = C, NoNewInfectionCarrier = H0,
       ActiveStates = states, VerticalTargetState = VerticalTargetState)
}

.ivp_disease_step_from_B <- function(Model, B, RecoveryProb, PathogenMortalityProb,
    ProgressionProb = 1) {
  Model <- match.arg(Model, c("SIS", "SIR", "SEIR")); B <- as.matrix(B); K <- nrow(B)
  if (ncol(B) != K) stop("B must be square")
  rec <- as.numeric(RecoveryProb); mort <- as.numeric(PathogenMortalityProb)
  prog <- as.numeric(ProgressionProb)
  if (length(rec) == 1L) rec <- rep(rec, K); if (length(mort) == 1L) mort <- rep(mort, K)
  if (length(prog) == 1L) prog <- rep(prog, K)
  if (length(rec) != K || length(mort) != K || length(prog) != K ||
      any(!is.finite(c(rec, mort, prog))) || any(rec < 0 | mort < 0 | prog < 0) ||
      any(rec + mort > 1 + 1e-12) || any(prog > 1)) stop("invalid pathogen transition probabilities")
  DI <- diag(1 - rec - mort, K)
  if (Model %in% c("SIS", "SIR")) {
    D <- DI + B; D0 <- DI
  } else {
    P <- diag(prog, K); E0 <- diag(1 - prog, K); Z <- matrix(0, K, K)
    D <- rbind(cbind(E0, B), cbind(P, DI))
    D0 <- rbind(cbind(E0, Z), cbind(P, DI))
  }
  list(DiseaseOperator = D, DiseaseOperatorNoHorizontalInfection = D0,
       TransmissionBlock = B)
}

.ivp_point_active_initial <- function(InitialPoints, type_obj, analysis_grid,
    Nstages, StateColumns, interaction, InitialInfectedByType = NULL,
    InitialExposedByType = NULL) {
  K <- nrow(type_obj$types); model <- interaction$Pathogen$Model
  count_subset <- function(z) {
    if (!nrow(z)) return(rep(0, K))
    .ivps_initial_counts(z, type_obj, analysis_grid, Nstages, StateColumns)
  }
  sf <- interaction$PathogenStateField
  if (!is.null(InitialInfectedByType)) I <- as.numeric(InitialInfectedByType)
  else if (sf %in% names(InitialPoints)) I <- count_subset(InitialPoints[as.character(InitialPoints[[sf]]) == "I", , drop = FALSE])
  else I <- rep(0, K)
  if (length(I) != K) stop("InitialInfectedByType must have one value per analytical type")
  if (model == "SEIR") {
    if (!is.null(InitialExposedByType)) E <- as.numeric(InitialExposedByType)
    else if (sf %in% names(InitialPoints)) E <- count_subset(InitialPoints[as.character(InitialPoints[[sf]]) == "E", , drop = FALSE])
    else E <- rep(0, K)
    if (length(E) != K) stop("InitialExposedByType must have one value per analytical type")
    c(E, I)
  } else I
}

.ivp_point_disease_operator <- function(Model, Types, HostPreInteraction,
    PathogenInteraction, timestep, H, R = NULL, VerticalTransmissionProb = 0,
    VerticalRecruitOperator = NULL, VerticalSourceState = "I",
    VerticalTargetState = NULL) {
  ps <- .ivp_spec(PathogenInteraction); K <- nrow(Types)
  cr <- .ina_point_interaction_contact(PathogenInteraction, "ContactRadius", Inf)
  ck <- .ina_point_interaction_contact(PathogenInteraction, "ContactKernel", NULL)
  cp <- .ina_point_interaction_contact(PathogenInteraction, "ContactProb", 1)
  Q <- INApestPointPathogenEdgeMatrix(Types, PathogenInteraction,
    ContactRadius = cr, ContactKernel = ck, ContactProb = cp,
    timestep = timestep, perm = 1L, IncludeDiagonal = TRUE)
  # Infinitesimal disease-free Jacobian: a source I is an infinitesimal
  # reclassification, so the susceptible background remains N_j to first order.
  B <- sweep(Q, 1, as.numeric(HostPreInteraction), `*`)
  rec <- .ina_point_resolve(ps$RecoveryProb, Types, timestep, 1L, "RecoveryProb")
  mort <- .ina_point_resolve(ps$PathogenMortalityProb, Types, timestep, 1L, "PathogenMortalityProb")
  prog <- .ina_point_resolve(ps$ProgressionProb, Types, timestep, 1L, "ProgressionProb")
  V <- if (!is.null(VerticalRecruitOperator)) {
    z <- if (is.function(VerticalRecruitOperator)) VerticalRecruitOperator(timestep = timestep, types = Types, Recruit = R) else VerticalRecruitOperator
    as.matrix(z)
  } else if (any(as.numeric(VerticalTransmissionProb) != 0)) {
    if (is.null(R)) stop("Vertical transmission with custom Birth requires VerticalRecruitOperator")
    .ivp_vertical_matrix(R, VerticalTransmissionProb, K)
  } else matrix(0, K, K)
  if (!all(dim(V) == c(K, K)) || any(!is.finite(V)) || any(V < 0)) stop("invalid VerticalRecruitOperator")
  car <- .ivp_active_carrier(Model, H, V, VerticalSourceState, VerticalTargetState)
  dis <- .ivp_disease_step_from_B(Model, B, rec, mort, prog)
  G <- dis$DiseaseOperator %*% car$PreDiseaseCarrier
  T0 <- dis$DiseaseOperatorNoHorizontalInfection %*% car$NoNewInfectionCarrier
  ng <- .ivp_ngm(G, T0)
  # A finite one-carrier lineage sees one fewer susceptible host in its own
  # analytical type. This is distinct from the infinitesimal Jacobian B above.
  Ssus <- matrix(as.numeric(HostPreInteraction), nrow=K, ncol=K)
  diag(Ssus) <- pmax(0, as.numeric(HostPreInteraction)-1)
  B1 <- Q * Ssus
  dis1 <- .ivp_disease_step_from_B(Model, B1, rec, mort, prog)
  if (max(abs(V)) < 1e-14) {
    if (Model %in% c("SIS","SIR")) G1 <- dis1$DiseaseOperator %*% H
    else {Z<-matrix(0,K,K);H2<-rbind(cbind(H,Z),cbind(Z,H));G1<-dis1$DiseaseOperator%*%H2}
    one_note <- "Finite one-carrier mean removes the carrier itself from its own susceptible target type."
  } else {
    # Vertical offspring can create several pre-Interaction active carriers from
    # one mother; their horizontal descendants are not independent. Keep the
    # exact mean G but do not label the resulting offspring law finite-binomial.
    G1 <- G
    one_note <- "With vertical transmission, one maternal lineage can contain multiple pre-Interaction active carriers; OneCarrierMeanOperator is therefore mean-matched rather than a finite-binomial offspring law."
  }
  list(Operator = G, OneCarrierMeanOperator=G1, OneCarrierMultiplier=.ivp_rho(G1),
       NoNewInfectionOperator = T0, NextGeneration = ng,
       TransmissionBlock = B, OneCarrierTransmissionBlock=B1, OneSourceTransmissionProbability = Q,
       ExistingHostCarrier = H, VerticalRecruitCarrier = V,
       PreDiseaseActiveCarrier = car$PreDiseaseCarrier,
       DiseaseOperator = dis$DiseaseOperator, OneCarrierDiagnostic=one_note,
       Lambda = .ivp_rho(G), R0 = ng$R0)
}

INApestVertebratePathogenAnalyticalPoint <- function(
    Ntimesteps = 10L, Nstages, Transition, InitialPoints, Pathogen,
    SDDkernel, LDDkernel = NULL, LDDrate = 0,
    PropaguleEstablishment = 1, EnvEstabProb = 1,
    TransitionKernels = NULL, TransitionHabitatSearch = FALSE,
    ApplyHabitatToTransitions = FALSE, TransitionEstablishment = 1,
    BlockedTransitionMortality = 0, HabitatSuitability = NULL,
    HabitatSearchRadius = 0, HabitatSearchCandidates = 128,
    MortalityProb = 0, MortalitySpatial = NULL, ManagementExposure = 0,
    FecundityReduction = 0, FecundityReductionSpatial = NULL,
    SpreadReduction = 0, SpreadReductionSpatial = NULL,
    SpreadReductionAppliesTo = c("LDD", "all"),
    OutsideEstablishmentProb = 1, Vertebrate = NULL,
    StateColumns = NULL, StateLevels = NULL, DefaultOffspringState = NULL,
    PointAnalysisGrid = NULL, KernelSamples = 2000L, PointSeed = 1L,
    VerticalTransmissionProb = 0, VerticalRecruitOperator = NULL,
    VerticalSourceState = "I", VerticalTargetState = NULL,
    InitialInfectedByType = NULL, InitialExposedByType = NULL,
    ReturnOperators = FALSE) {

  if (!inherits(Pathogen, "INApestPointPathogenInteraction"))
    stop("Vertebrate Point pathogen analysis requires INApestPointPathogenInteraction so contact semantics match the stochastic Interaction hook")
  Ntimesteps <- as.integer(Ntimesteps); Nstages <- as.integer(Nstages)
  if (Ntimesteps < 1L || Nstages < 2L) stop("Ntimesteps >= 1 and Nstages >= 2 required")
  ps <- .ivp_spec(Pathogen); Model <- ps$Model
  if (!Model %in% c("SIS", "SIR", "SEIR")) stop("SIS/SIR/SEIR only")
  if (!is.null(Vertebrate) && !is.null(Vertebrate$Interaction))
    stop("Pathogen occupies Vertebrate$Interaction in the stochastic Point engine. Supply social/contact heterogeneity through persistent StateColumns/Pathogen ContactProb, or use a custom joint analytical map; do not supply a second Vertebrate$Interaction here.")
  SpreadReductionAppliesTo <- match.arg(SpreadReductionAppliesTo)
  if (Pathogen$PathogenStateField %in% StateColumns)
    stop("Do not include the pathogen-state field in StateColumns: pathogen state is an analytical compartment layered on top of persistent vertebrate host types.")
  p <- .ivps_point_frame(InitialPoints)
  analysis_grid <- .ina_pt_analysis_grid(PointAnalysisGrid, HabitatSuitability,
    MortalitySpatial = MortalitySpatial, FecundityReductionSpatial = FecundityReductionSpatial,
    SpreadReductionSpatial = SpreadReductionSpatial)
  include_na <- length(StateColumns) && is.null(DefaultOffspringState)
  to <- .ivps_type_table(p, Nstages, analysis_grid, StateColumns, StateLevels, include_na = include_na)
  x <- .ivps_initial_counts(p, to, analysis_grid, Nstages, StateColumns)
  K <- length(x); z0 <- .ivp_point_active_initial(p, to, analysis_grid, Nstages,
    StateColumns, Pathogen, InitialInfectedByType, InitialExposedByType)
  adim <- length(z0)
  host_traj <- matrix(0, Ntimesteps + 1L, K); host_traj[1L, ] <- x
  active_traj <- matrix(0, Ntimesteps + 1L, adim); active_traj[1L, ] <- z0
  ops <- vector("list", Ntimesteps); ngs <- vector("list", Ntimesteps)
  host_pre <- vector("list", Ntimesteps); escape_mean <- vector("list", Ntimesteps)
  vert_host <- Vertebrate
  if (is.null(vert_host)) vert_host <- list()
  vert_host$Interaction <- NULL

  for (tt in seq_len(Ntimesteps)) {
    st <- .ivps_step_builder(tt, x, to, analysis_grid, Ntimesteps, Nstages,
      Transition, SDDkernel, LDDkernel, LDDrate,
      PropaguleEstablishment, EnvEstabProb, TransitionKernels,
      TransitionHabitatSearch, ApplyHabitatToTransitions, TransitionEstablishment,
      BlockedTransitionMortality, HabitatSuitability, HabitatSearchRadius,
      HabitatSearchCandidates, MortalityProb, MortalitySpatial,
      ManagementExposure, FecundityReduction, FecundityReductionSpatial,
      SpreadReduction, SpreadReductionSpatial, SpreadReductionAppliesTo,
      OutsideEstablishmentProb, vert_host, StateColumns, DefaultOffspringState,
      KernelSamples, PointSeed)
    pre <- as.numeric(st$pre_interaction_state)
    host_pre[[tt]] <- pre
    op <- .ivp_point_disease_operator(Model, st$type_table, pre, Pathogen, tt,
      st$Parent, st$Recruit, VerticalTransmissionProb, VerticalRecruitOperator,
      VerticalSourceState, VerticalTargetState)
    ops[[tt]] <- op; ngs[[tt]] <- op$NextGeneration

    # Active host export occurs before disease Interaction. Vertical infected
    # recruit export is included when the default recruit source decomposition is available.
    eH <- as.numeric(st$export_transition_by_source)
    if (length(eH) != K) eH <- rep(0, K)
    eVsrc <- rep(0, K)
    if (any(as.numeric(VerticalTransmissionProb) != 0) && !is.null(st$Recruit)) {
      er <- as.numeric(st$export_recruit_by_source)
      vp <- as.numeric(VerticalTransmissionProb)
      if (length(vp) == 1L) eVsrc <- er * vp
      else if (length(vp) == K) eVsrc <- er * vp
    }
    if (Model %in% c("SIS", "SIR")) escape_mean[[tt]] <- eH + eVsrc
    else {
      e <- c(eH, eH)
      srcblock <- if (VerticalSourceState == "E") seq_len(K) else K + seq_len(K)
      e[srcblock] <- e[srcblock] + eVsrc
      escape_mean[[tt]] <- e
    }

    active_traj[tt + 1L, ] <- as.numeric(op$Operator %*% active_traj[tt, ])
    x <- as.numeric(st$next_state); host_traj[tt + 1L, ] <- x
  }
  ord <- .ivp_ordered_growth(lapply(ops, `[[`, "Operator"))
  br <- .ivp_poisson_branching(lapply(ops, `[[`, "OneCarrierMeanOperator"), z0, escape_mean,
                               Eventual = length(ops) == 1L)
  out <- list(
    Model = "INApestVertebratePoint", PathogenModel = Model,
    GrowthRate = ord$PerTimestepMultiplier,
    PathogenGrowthRate = ord$PerTimestepMultiplier,
    Classification = ord$Classification,
    StepLambda = vapply(ops, `[[`, numeric(1), "Lambda"),
    StepOneCarrierMultiplier = vapply(ops, `[[`, numeric(1), "OneCarrierMultiplier"),
    StepR0 = vapply(ops, function(o) o$R0, numeric(1)),
    HostEndPopulationDiseaseFree = sum(host_traj[Ntimesteps + 1L, ]),
    EndActivePathogen = sum(active_traj[Ntimesteps + 1L, ]),
    HostTrajectoryDiseaseFree = rowSums(host_traj),
    HostStateTrajectoryDiseaseFree = host_traj,
    ActivePathogenTrajectory = rowSums(active_traj),
    ActivePathogenStateTrajectory = active_traj,
    HostPopulationAtDiseaseStep = host_pre,
    EscapeProbability = br$EscapeProbabilityByHorizon,
    ExtinctionProbabilityByHorizon = br$ExtinctionProbabilityByHorizon,
    EventualExtinctionProbability = br$EventualExtinctionProbability,
    Branching = br,
    TypeTable = to$types,
    VerticalTransmission = list(Probability = VerticalTransmissionProb,
      SourceState = VerticalSourceState,
      TargetState = if (is.null(VerticalTargetState)) if (Model == "SEIR") "E" else "I" else VerticalTargetState),
    Diagnostics = c(
      "Host demography/control is evaluated on the disease-free vertebrate trajectory; active pathogen states are then rare reclassifications carried through the exact existing-host first-moment operator.",
      "Routine vertebrate control and information-triggered response mortality occur before birth/stage transport; pathogen transmission/progression/recovery occurs in the Interaction slot after host transport.",
      "Vertically infected recruits enter the pre-Interaction active carrier. Consequently vertical I recruits can transmit/recover/die in their birth timestep, and vertical E recruits can progress in that timestep; horizontally new infections cannot.",
      "R0 counts all new infected lineages, including vertical transmission and same-timestep horizontal descendants of vertically infectious recruits, by defining F = G - T0 rather than treating vertical transmission as persistence.",
      "For Point models, the infinitesimal invasion multiplier and finite one-carrier branching multiplier are reported separately: a real carrier removes itself from its own susceptible type, which matters in small vertebrate populations.",
      "Dynamic pathogen-detection -> shared HaveInfo -> later response management is not closed under independent lineage branching; ManagementExposure should be used as a none-informed/all-informed or otherwise externally specified host-management envelope.",
      "Finite prevalence can feed back on density-dependent HomeRange, mating, control exposure and social state. The rare-pathogen operator does not silently approximate that nonlinear feedback; use stochastic simulation for epidemic-burden trajectories.",
      "Individual-lineage Poisson branching is invalid when social groups create correlated transmission/movement. Use INApestVertebratePathogenGroupBranching or a custom joint group PGF in that case."))
  if (ReturnOperators) out$StepOperators <- ops
  class(out) <- c("INApestVertebratePathogenAnalyticalPoint", "list")
  out
}

.ivp_node_resolve <- function(x, timestep, Ntimesteps, n, S, name, prob = FALSE) {
  if (is.function(x)) {
    fm <- names(formals(x)); a <- list(timestep = timestep, Ntimesteps = Ntimesteps,
      n_nodes = n, n_stages = S)
    if (!is.null(fm) && !"..." %in% fm) a <- a[intersect(names(a), fm)]
    x <- do.call(x, a)
  }
  .ina_tm_ns_matrix(x, n, S, name, prob)
}

.ivp_node_contact <- function(x, timestep, Ntimesteps, n) {
  if (is.null(x)) return(diag(n))
  if (is.function(x)) {
    fm <- names(formals(x)); a <- list(timestep = timestep, Ntimesteps = Ntimesteps, n_nodes = n)
    if (!is.null(fm) && !"..." %in% fm) a <- a[intersect(names(a), fm)]
    x <- do.call(x, a)
  }
  z <- as.matrix(x)
  if (!all(dim(z) == c(n, n)) || any(!is.finite(z)) || any(z < 0)) stop("ContactMatrix must resolve to non-negative nodes x nodes")
  z
}

.ivp_node_parts <- function(population, timestep, Ntimesteps, Nstages,
    Transition, SDDprob, LDDprob, LDDrate, TransitionSDDprob,
    TransitionLDDprob, TransitionLDDrate, PropaguleEstablishment,
    EnvEstabProb, MortalityProb, ManagementExposure, FecundityReduction,
    OutsideEstablishmentProb, Vertebrate, K, Weights, DispersalDensityFactor) {
  if (!is.null(Vertebrate) && !is.null(Vertebrate$Interaction))
    stop("Separate Pathogen analysis occupies the node Interaction slot; supply no second Vertebrate$Interaction unless you provide a custom joint analytical map")
  N <- as.matrix(population); n <- nrow(N); S <- Nstages; t <- timestep; T <- Ntimesteps
  comp <- INApestVertebrateNodeMeanStep(N, t, T, S, Transition, SDDprob, LDDprob, LDDrate,
    TransitionSDDprob, TransitionLDDprob, TransitionLDDrate,
    PropaguleEstablishment, EnvEstabProb, MortalityProb, ManagementExposure,
    FecundityReduction, OutsideEstablishmentProb, Vertebrate, K, Weights,
    DispersalDensityFactor, ReturnComponents = TRUE)
  Knode <- .ivspec_node_time(K, t, T, n, "K", default = Inf, allow_inf = TRUE)
  W <- .ivspec_weights(Weights, n, S)
  transition_hook <- .ivspec_transition_for_hook(Transition, t, n, S, T)
  trans <- .ivspec_transition_list(Transition, t, n, S, T)
  HomeRangeModule <- .ivspec_module(Vertebrate, "HomeRange")
  ControlModule <- .ivspec_module(Vertebrate, "Control")
  BirthModule <- .ivspec_module(Vertebrate, "Birth")
  hr <- .ivspec_home_range_node(HomeRangeModule, N, t,
    list(transition = transition_hook, K = Knode, Nstages = S, Weights = W, phase = "pre_response_control"))
  ce <- .ivspec_node_control(ControlModule, N, hr, t,
    list(transition = transition_hook, K = Knode, Nstages = S, Weights = W, phase = "pre_response_control"))
  rm <- .ivanal_expand_node_stage(MortalityProb, n, S, t, T, "MortalityProb", 0)
  a <- .ivanal_expand_node_stage(ManagementExposure, n, S, t, T, "ManagementExposure", 0)
  rm <- .ivspec_clip01(rm); a <- .ivspec_clip01(a)
  rs <- (1 - a) + a * (1 - rm)
  q <- (1 - ce$kill_prob) * rs
  U <- n * S; idx <- function(i,s) (i - 1L) * S + s
  H <- matrix(0, U, U); eH <- numeric(U)
  for (i in seq_len(n)) for (s in seq_len(S)) {
    src <- idx(i,s); unit <- matrix(0, n, S); unit[i,s] <- q[i,s]
    ex <- .ivspec_transition_existing(unit, t, T, S, Transition,
      TransitionSDDprob, TransitionLDDprob, TransitionLDDrate,
      OutsideEstablishmentProb)
    H[,src] <- as.numeric(t(ex$inside)); eH[src] <- ex$export
  }
  R <- NULL; eR <- numeric(U)
  if (is.null(BirthModule)) {
    R <- matrix(0, U, U)
    af <- matrix(0, n, S); ok <- rs > 0
    af[ok] <- a[ok] * (1 - rm[ok]) / rs[ok]
    rf <- .ivanal_expand_node_stage(FecundityReduction, n, S, t, T, "FecundityReduction", 0)
    fm <- 1 - .ivspec_clip01(rf) * af
    bm <- .ivspec_birth_movement(SDDprob, LDDprob, LDDrate,
      N * q, Knode, W, DispersalDensityFactor, t, T)
    if (is.null(bm)) bm <- diag(n)
    env <- .ivspec_node_time(EnvEstabProb, t, T, n, "EnvEstabProb")
    pe <- .ivspec_node_time(PropaguleEstablishment, t, T, n, "PropaguleEstablishment")
    estab <- .ivspec_clip01(env * pe)
    oe <- .ivspec_node_time(OutsideEstablishmentProb, t, T, n, "OutsideEstablishmentProb")
    for (i in seq_len(n)) for (s in 2:S) {
      src <- idx(i,s); fec <- pmax(0, trans[[i]][1,s])
      b <- q[i,s] * fec * fm[i,s] * (1 - ce$fecundity_reduction[i,s])
      if (b <= 0) next
      pin <- pmax(0, bm[i,])
      for (j in seq_len(n)) R[idx(j,1L),src] <- R[idx(j,1L),src] + b * pin[j] * estab[j]
      eR[src] <- b * pmax(0, 1 - sum(pin)) * .ivspec_clip01(oe[i])
    }
  }
  list(next_host = comp$population,
       pre_interaction_host = comp$expected_existing_inside + comp$expected_births_inside,
       Parent = H, Recruit = R, export_transition_by_source = eH,
       export_recruit_by_source = eR, control = ce)
}

.ivp_node_disease_operator <- function(Model, HostPre, Pathogen, timestep, Ntimesteps,
    H, R, StageMixing = NULL, VerticalTransmissionProb = 0,
    VerticalRecruitOperator = NULL, VerticalSourceState = "I", VerticalTargetState = NULL) {
  ps <- .ivp_spec(Pathogen); Nmat <- as.matrix(HostPre); n <- nrow(Nmat); S <- ncol(Nmat); K <- n*S
  beta <- .ivp_node_resolve(ps$Beta, timestep, Ntimesteps, n, S, "Beta")
  rec <- .ivp_node_resolve(ps$RecoveryProb, timestep, Ntimesteps, n, S, "RecoveryProb", TRUE)
  mort <- .ivp_node_resolve(ps$PathogenMortalityProb, timestep, Ntimesteps, n, S, "PathogenMortalityProb", TRUE)
  prog <- .ivp_node_resolve(ps$ProgressionProb, timestep, Ntimesteps, n, S, "ProgressionProb", TRUE)
  ds <- .ivp_node_resolve(ps$DensityScale, timestep, Ntimesteps, n, S, "DensityScale")
  Cn <- .ivp_node_contact(ps$ContactMatrix, timestep, Ntimesteps, n)
  dl <- INApestTransitionPathogenDiseaseLinearOperator(Model, Nmat, beta, rec, prog, mort,
    Cn, StageMixing, ps$Transmission, ds)
  V <- if (!is.null(VerticalRecruitOperator)) {
    z <- if (is.function(VerticalRecruitOperator)) VerticalRecruitOperator(timestep=timestep, HostPre=HostPre, Recruit=R) else VerticalRecruitOperator
    as.matrix(z)
  } else if (any(as.numeric(VerticalTransmissionProb) != 0)) {
    if (is.null(R)) stop("Vertical transmission with custom node Birth requires VerticalRecruitOperator")
    .ivp_vertical_matrix(R, VerticalTransmissionProb, K)
  } else matrix(0, K, K)
  car <- .ivp_active_carrier(Model, H, V, VerticalSourceState, VerticalTargetState)
  # Rebuild D0 from B=0 using exactly the resolved probabilities.
  B <- dl$TransmissionBlock
  dis <- .ivp_disease_step_from_B(Model, B, as.numeric(t(rec)), as.numeric(t(mort)), as.numeric(t(prog)))
  G <- dis$DiseaseOperator %*% car$PreDiseaseCarrier
  T0 <- dis$DiseaseOperatorNoHorizontalInfection %*% car$NoNewInfectionCarrier
  ng <- .ivp_ngm(G, T0)
  if(max(abs(V))<1e-14){
    N<-as.numeric(t(Nmat));beta_v<-as.numeric(t(beta));rec_v<-as.numeric(t(rec));mort_v<-as.numeric(t(mort));prog_v<-as.numeric(t(prog));ds_v<-as.numeric(t(ds))
    Sm<-if(is.null(StageMixing))matrix(1,S,S)else as.matrix(StageMixing);C<-kronecker(Cn,Sm)
    pinf_one<-function(v){
      out<-numeric(K);if(ps$Transmission=="frequency"){den<-as.numeric(crossprod(N,C));for(w in seq_len(K))if(den[w]>0&&beta_v[w]>0&&C[v,w]>0)out[w]<--expm1(-beta_v[w]*C[v,w]/den[w])}
      else for(w in seq_len(K))if(beta_v[w]>0&&C[v,w]>0)out[w]<--expm1(-beta_v[w]*C[v,w]/ds_v[w]);pmin(1,pmax(0,out))
    }
    if(Model%in%c("SIS","SIR")){
      G1<-matrix(0,K,K);for(u in seq_len(K))for(v in which(H[,u]>0)){pp<-pinf_one(v);ss<-pmax(0,N-as.numeric(seq_len(K)==v));G1[,u]<-G1[,u]+H[v,u]*(ss*pp);G1[v,u]<-G1[v,u]+H[v,u]*(1-rec_v[v]-mort_v[v])}
    }else{
      G1<-matrix(0,2*K,2*K);for(u in seq_len(K))for(v in which(H[,u]>0)){G1[v,u]<-G1[v,u]+H[v,u]*(1-prog_v[v]);G1[K+v,u]<-G1[K+v,u]+H[v,u]*prog_v[v];pp<-pinf_one(v);ss<-pmax(0,N-as.numeric(seq_len(K)==v));G1[seq_len(K),K+u]<-G1[seq_len(K),K+u]+H[v,u]*(ss*pp);G1[K+v,K+u]<-G1[K+v,K+u]+H[v,u]*(1-rec_v[v]-mort_v[v])}
    }
    one_note<-"Finite one-carrier mean uses the stochastic abundance engine's exact exponential infection probability and removes the carrier from its own susceptible unit."
  } else {G1<-G;one_note<-"With vertical transmission the one-mother lineage can contain multiple pre-Interaction active carriers, so the one-carrier mean is mean-matched to G rather than a finite-binomial offspring law."}
  list(Operator=G, OneCarrierMeanOperator=G1, OneCarrierMultiplier=.ivp_rho(G1),
       NoNewInfectionOperator=T0, NextGeneration=ng,
       TransmissionBlock=B, ExistingHostCarrier=H, VerticalRecruitCarrier=V,
       PreDiseaseActiveCarrier=car$PreDiseaseCarrier, DiseaseOperator=dis$DiseaseOperator,
       OneCarrierDiagnostic=one_note, Lambda=.ivp_rho(G), R0=ng$R0)
}

INApestVertebratePathogenAnalyticalNode <- function(
    Ntimesteps = 10L, Nstages, Transition, InitialPopulation, Pathogen,
    SDDprob, LDDprob = NULL, LDDrate = 0,
    TransitionSDDprob = NULL, TransitionLDDprob = NULL, TransitionLDDrate = 0,
    PropaguleEstablishment = 1, EnvEstabProb = 1, MortalityProb = 0,
    ManagementExposure = 0, FecundityReduction = 0,
    OutsideEstablishmentProb = 1, Vertebrate = NULL, K = Inf,
    Weights = rep(1, Nstages), DispersalDensityFactor = 0,
    StageMixing = NULL, VerticalTransmissionProb = 0,
    VerticalRecruitOperator = NULL, VerticalSourceState = "I",
    VerticalTargetState = NULL, InitialInfected = NULL, InitialExposed = NULL,
    ReturnOperators = FALSE) {
  Ntimesteps <- as.integer(Ntimesteps); Nstages <- as.integer(Nstages)
  if (Ntimesteps < 1L || Nstages < 2L) stop("Ntimesteps >= 1 and Nstages >= 2 required")
  ps <- .ivp_spec(Pathogen); Model <- ps$Model
  if (!Model %in% c("SIS","SIR","SEIR")) stop("SIS/SIR/SEIR only")
  x <- as.matrix(InitialPopulation); n <- nrow(x); S <- Nstages; U <- n*S
  if (ncol(x) != S) stop("InitialPopulation must be nodes x Nstages")
  I0 <- if (is.null(InitialInfected)) matrix(0,n,S) else .ina_tm_ns_matrix(InitialInfected,n,S,"InitialInfected")
  if (any(I0 > x + 1e-12)) stop("InitialInfected exceeds host population")
  if (Model == "SEIR") {
    E0 <- if (is.null(InitialExposed)) matrix(0,n,S) else .ina_tm_ns_matrix(InitialExposed,n,S,"InitialExposed")
    if (any(E0 + I0 > x + 1e-12)) stop("Initial active pathogen counts exceed host population")
    z0 <- c(as.numeric(t(E0)), as.numeric(t(I0)))
  } else z0 <- as.numeric(t(I0))
  adim <- length(z0)
  host_traj <- matrix(0,Ntimesteps+1L,U); host_traj[1L,] <- as.numeric(t(x))
  active_traj <- matrix(0,Ntimesteps+1L,adim); active_traj[1L,] <- z0
  ops <- vector("list",Ntimesteps); prehost <- vector("list",Ntimesteps); escape_mean <- vector("list",Ntimesteps)
  for (tt in seq_len(Ntimesteps)) {
    st <- .ivp_node_parts(x,tt,Ntimesteps,S,Transition,SDDprob,LDDprob,LDDrate,
      TransitionSDDprob,TransitionLDDprob,TransitionLDDrate,PropaguleEstablishment,
      EnvEstabProb,MortalityProb,ManagementExposure,FecundityReduction,
      OutsideEstablishmentProb,Vertebrate,K,Weights,DispersalDensityFactor)
    op <- .ivp_node_disease_operator(Model,st$pre_interaction_host,Pathogen,tt,Ntimesteps,
      st$Parent,st$Recruit,StageMixing,VerticalTransmissionProb,VerticalRecruitOperator,
      VerticalSourceState,VerticalTargetState)
    ops[[tt]] <- op; prehost[[tt]] <- st$pre_interaction_host
    eH <- st$export_transition_by_source; eV <- rep(0,U)
    if (any(as.numeric(VerticalTransmissionProb) != 0) && !is.null(st$Recruit)) {
      vp <- as.numeric(VerticalTransmissionProb)
      if (length(vp)==1L) eV <- st$export_recruit_by_source*vp
      else if (length(vp)==U) eV <- st$export_recruit_by_source*vp
    }
    if (Model %in% c("SIS","SIR")) escape_mean[[tt]] <- eH+eV
    else {
      ee <- c(eH,eH); src <- if(VerticalSourceState=="E")seq_len(U) else U+seq_len(U)
      ee[src] <- ee[src]+eV; escape_mean[[tt]] <- ee
    }
    active_traj[tt+1L,] <- as.numeric(op$Operator %*% active_traj[tt,])
    x <- as.matrix(st$next_host); host_traj[tt+1L,] <- as.numeric(t(x))
  }
  ord <- .ivp_ordered_growth(lapply(ops,`[[`,"Operator"))
  br <- .ivp_poisson_branching(lapply(ops,`[[`,"OneCarrierMeanOperator"),z0,escape_mean,Eventual=length(ops)==1L)
  out <- list(Model="INApestVertebrateNode",PathogenModel=Model,
    GrowthRate=ord$PerTimestepMultiplier,PathogenGrowthRate=ord$PerTimestepMultiplier,
    Classification=ord$Classification,StepLambda=vapply(ops,`[[`,numeric(1),"Lambda"),
    StepOneCarrierMultiplier=vapply(ops,`[[`,numeric(1),"OneCarrierMultiplier"),
    StepR0=vapply(ops,function(o)o$R0,numeric(1)),
    HostEndPopulationDiseaseFree=sum(host_traj[Ntimesteps+1L,]),
    EndActivePathogen=sum(active_traj[Ntimesteps+1L,]),
    HostTrajectoryDiseaseFree=rowSums(host_traj),HostStateTrajectoryDiseaseFree=host_traj,
    ActivePathogenTrajectory=rowSums(active_traj),ActivePathogenStateTrajectory=active_traj,
    HostPopulationAtDiseaseStep=prehost,EscapeProbability=br$EscapeProbabilityByHorizon,
    ExtinctionProbabilityByHorizon=br$ExtinctionProbabilityByHorizon,
    EventualExtinctionProbability=br$EventualExtinctionProbability,Branching=br,
    VerticalTransmission=list(Probability=VerticalTransmissionProb,SourceState=VerticalSourceState,
      TargetState=if(is.null(VerticalTargetState))if(Model=="SEIR")"E"else"I"else VerticalTargetState),
    Diagnostics=c(
      "Node vertebrate disease analysis composes the validated routine-control/response-management/stage-transport carrier before the generic pathogen Interaction step.",
      "Demographic births determine the susceptible host background but do not create infected offspring unless vertical transmission is explicitly supplied.",
      "Vertical transmission is represented in the pre-Interaction carrier, so its within-birth-timestep progression/recovery/transmission timing differs from horizontal new infection and matches the vertebrate Interaction architecture.",
      "Growth and R0 are rare-pathogen quantities evaluated along the disease-free host trajectory. Disease-mortality feedback on density-dependent HomeRange, mating, social state and future demography is intentionally left to stochastic simulation at material prevalence.",
      "The infinitesimal lambda and finite one-carrier multiplier are both returned. They can differ materially in small vertebrate populations because the carrier itself is not susceptible and the stochastic infection probability is exponential rather than its infinitesimal linearisation.",
      "Pathogen-triggered shared information and future response management break independent-lineage closure; ManagementExposure is an externally specified management envelope."))
  if(ReturnOperators)out$StepOperators<-ops
  class(out)<-c("INApestVertebratePathogenAnalyticalNode","list");out
}

INApestVertebratePathogenGroupBranching <- function(
    MeanInfectedGroupOffspring, InitialInfectedGroups, EscapeProb = 0,
    Ntimesteps = 10L, GroupPGF = NULL, ExtinctionGenerations = 500L,
    tolerance = 1e-12) {
  z <- INApestVertebrateGroupBranching(MeanInfectedGroupOffspring,
    InitialInfectedGroups, EscapeProb, Ntimesteps, GroupPGF,
    ExtinctionGenerations, tolerance)
  z$PathogenInterpretation <- paste(
    "Each lineage is an infected social group rather than an animal.",
    "MeanInfectedGroupOffspring should integrate the within-group epidemic/contact process and count newly infected descendant groups.",
    "Use a supplied GroupPGF when group offspring are not adequately Poisson.")
  z$R0Group <- z$GrowthRate
  class(z) <- c("INApestVertebratePathogenGroupBranching","list")
  z
}

###############################################################################
### Unified dispatcher refinement
###############################################################################
INApestAnalytical_pre_vertebrate_pathogen_specialist <- INApestAnalytical
INApestAnalytical <- function(...) {
  args <- list(...)
  Model <- if (!is.null(args$Model)) as.character(args$Model)[1L] else "INApest"
  Pathogen <- args$Pathogen
  if (is.null(Pathogen) || !Model %in% c("INApestVertebrateNode","INApestVertebratePoint"))
    return(do.call(INApestAnalytical_pre_vertebrate_pathogen_specialist,args))
  if (Model == "INApestVertebrateNode") {
    keep <- names(formals(INApestVertebratePathogenAnalyticalNode))
    aa <- args[intersect(names(args),keep)]
    if (is.null(aa$InitialPopulation)) aa$InitialPopulation <- args$InitialState
    if (is.null(aa$ManagementExposure)) aa$ManagementExposure <- args$ManageProb %||% 0
    ans <- do.call(INApestVertebratePathogenAnalyticalNode,aa)
  } else {
    keep <- names(formals(INApestVertebratePathogenAnalyticalPoint))
    aa <- args[intersect(names(args),keep)]
    if (is.null(aa$ManagementExposure)) aa$ManagementExposure <- args$ManageProb %||% 0
    ans <- do.call(INApestVertebratePathogenAnalyticalPoint,aa)
  }
  ans$HeadlineEstimands <- list(PathogenGrowthRate=ans$PathogenGrowthRate,
    HostEndPopulationDiseaseFree=ans$HostEndPopulationDiseaseFree,
    EndActivePathogen=ans$EndActivePathogen,
    PathogenEscapeProbability=ans$EscapeProbability)
  ans
}

###############################################################################
### Finite-prevalence deterministic mean closure
###############################################################################

.ivp_point_initial_full <- function(InitialPoints, type_obj, analysis_grid, Nstages,
    StateColumns, PathogenInteraction, InitialInfectedByType=NULL,
    InitialExposedByType=NULL, InitialRecoveredByType=NULL) {
  ps <- .ivp_spec(PathogenInteraction); states <- ps$States; K <- nrow(type_obj$types)
  total <- .ivps_initial_counts(InitialPoints,type_obj,analysis_grid,Nstages,StateColumns)
  sf <- PathogenInteraction$PathogenStateField
  count_state <- function(st) {
    if (!(sf %in% names(InitialPoints))) return(rep(0,K))
    z <- InitialPoints[as.character(InitialPoints[[sf]])==st,,drop=FALSE]
    if (!nrow(z)) return(rep(0,K))
    .ivps_initial_counts(z,type_obj,analysis_grid,Nstages,StateColumns)
  }
  I <- if(is.null(InitialInfectedByType))count_state("I") else as.numeric(InitialInfectedByType)
  E <- if("E"%in%states)if(is.null(InitialExposedByType))count_state("E") else as.numeric(InitialExposedByType) else rep(0,K)
  R <- if("R"%in%states)if(is.null(InitialRecoveredByType))count_state("R") else as.numeric(InitialRecoveredByType) else rep(0,K)
  if(any(vapply(list(I,E,R),length,integer(1))!=K)||any(I+E+R>total+1e-10))stop("Initial pathogen state counts are incompatible with host type counts")
  out<-list(S=pmax(0,total-I-E-R),I=I);if("E"%in%states)out$E<-E;if("R"%in%states)out$R<-R
  out[states]
}

.ivp_point_horizontal_mean <- function(S,I,Q) {
  K<-length(S);p<-numeric(K)
  for(j in seq_len(K)){
    q<-pmin(1,pmax(0,Q[j,]));ii<-pmax(0,I)
    if(any(q>=1 & ii>0))p[j]<-1 else p[j]<-1-exp(sum(ii*log1p(-q)))
  }
  S*p
}

.ivp_point_joint_disease_step <- function(comp, Types, PathogenInteraction, timestep, H,
    RecruitTotal, Rop=NULL, VerticalTransmissionProb=0, VerticalRecruitOperator=NULL,
    VerticalSourceState="I", VerticalTargetState=NULL) {
  ps<-.ivp_spec(PathogenInteraction);model<-ps$Model;states<-ps$States;K<-nrow(Types)
  pre<-lapply(states,function(st)as.numeric(H%*%comp[[st]]));names(pre)<-states
  V<-if(!is.null(VerticalRecruitOperator)){
    z<-if(is.function(VerticalRecruitOperator))VerticalRecruitOperator(timestep=timestep,types=Types,Recruit=Rop)else VerticalRecruitOperator
    as.matrix(z)
  } else if(any(as.numeric(VerticalTransmissionProb)!=0)){
    if(is.null(Rop))stop("Finite vertical transmission with custom Birth requires VerticalRecruitOperator")
    .ivp_vertical_matrix(Rop,VerticalTransmissionProb,K)
  } else matrix(0,K,K)
  vt<-if(is.null(VerticalTargetState))if(model=="SEIR")"E"else"I"else VerticalTargetState
  vs<-VerticalSourceState
  vb<-as.numeric(V%*%comp[[vs]])
  susceptible_births<-as.numeric(RecruitTotal)-vb
  if(any(susceptible_births < -1e-8))stop("Vertical infected recruits exceed total expected recruits")
  pre$S<-pre$S+pmax(0,susceptible_births);pre[[vt]]<-pre[[vt]]+vb
  cr<-.ina_point_interaction_contact(PathogenInteraction,"ContactRadius",Inf)
  ck<-.ina_point_interaction_contact(PathogenInteraction,"ContactKernel",NULL)
  cp<-.ina_point_interaction_contact(PathogenInteraction,"ContactProb",1)
  Q<-INApestPointPathogenEdgeMatrix(Types,PathogenInteraction,ContactRadius=cr,ContactKernel=ck,ContactProb=cp,timestep=timestep,perm=1L,IncludeDiagonal=TRUE)
  new<-.ivp_point_horizontal_mean(pre$S,pre$I,Q)
  ip<-.ina_point_resolve(ps$IntroductionProb,Types,timestep,1L,"IntroductionProb");ip<-pmin(1,pmax(0,ip))
  intro<-(pre$S-new)*ip
  rec<-.ina_point_resolve(ps$RecoveryProb,Types,timestep,1L,"RecoveryProb")
  mort<-.ina_point_resolve(ps$PathogenMortalityProb,Types,timestep,1L,"PathogenMortalityProb")
  prog<-.ina_point_resolve(ps$ProgressionProb,Types,timestep,1L,"ProgressionProb")
  wan<-.ina_point_resolve(ps$ImmunityLossProb,Types,timestep,1L,"ImmunityLossProb")
  stayI<-pre$I*(1-rec-mort);recover<-pre$I*rec;deaths<-pre$I*mort
  if(model=="SIS"){
    out<-list(S=pre$S-new-intro+recover,I=stayI+new+intro)
  } else if(model=="SIR"){
    lose<-pre$R*wan
    out<-list(S=pre$S-new-intro+lose,I=stayI+new+intro,R=pre$R-lose+recover)
  } else {
    progress<-pre$E*prog;lose<-pre$R*wan
    out<-list(S=pre$S-new-intro+lose,E=pre$E-progress+new+intro,I=stayI+progress,R=pre$R-lose+recover)
  }
  list(State=out,PreInteraction=pre,NewInfections=new,Introduced=intro,PathogenDeaths=deaths,Q=Q,VerticalBirths=vb)
}

INApestVertebratePathogenMeanPoint <- function(
    Ntimesteps=10L,Nstages,Transition,InitialPoints,Pathogen,SDDkernel,
    LDDkernel=NULL,LDDrate=0,PropaguleEstablishment=1,EnvEstabProb=1,
    TransitionKernels=NULL,TransitionHabitatSearch=FALSE,ApplyHabitatToTransitions=FALSE,
    TransitionEstablishment=1,BlockedTransitionMortality=0,HabitatSuitability=NULL,
    HabitatSearchRadius=0,HabitatSearchCandidates=128,MortalityProb=0,MortalitySpatial=NULL,
    ManagementExposure=0,FecundityReduction=0,FecundityReductionSpatial=NULL,
    SpreadReduction=0,SpreadReductionSpatial=NULL,SpreadReductionAppliesTo=c("LDD","all"),
    OutsideEstablishmentProb=1,Vertebrate=NULL,StateColumns=NULL,StateLevels=NULL,
    DefaultOffspringState=NULL,PointAnalysisGrid=NULL,KernelSamples=2000L,PointSeed=1L,
    VerticalTransmissionProb=0,VerticalRecruitOperator=NULL,VerticalSourceState="I",
    VerticalTargetState=NULL,InitialInfectedByType=NULL,InitialExposedByType=NULL,
    InitialRecoveredByType=NULL) {
  if(!inherits(Pathogen,"INApestPointPathogenInteraction"))stop("Point mean trajectory requires INApestPointPathogenInteraction")
  if(!is.null(Vertebrate)&&!is.null(Vertebrate$Interaction))stop("Finite joint mean currently requires pathogen to be the sole Interaction update")
  Ntimesteps<-as.integer(Ntimesteps);Nstages<-as.integer(Nstages);SpreadReductionAppliesTo<-match.arg(SpreadReductionAppliesTo)
  if(Pathogen$PathogenStateField %in% StateColumns)stop("Do not include the pathogen-state field in StateColumns")
  p<-.ivps_point_frame(InitialPoints);ag<-.ina_pt_analysis_grid(PointAnalysisGrid,HabitatSuitability,MortalitySpatial=MortalitySpatial,FecundityReductionSpatial=FecundityReductionSpatial,SpreadReductionSpatial=SpreadReductionSpatial)
  to<-.ivps_type_table(p,Nstages,ag,StateColumns,StateLevels,include_na=length(StateColumns)&&is.null(DefaultOffspringState))
  comp<-.ivp_point_initial_full(p,to,ag,Nstages,StateColumns,Pathogen,InitialInfectedByType,InitialExposedByType,InitialRecoveredByType)
  states<-names(comp);K<-nrow(to$types);hist<-lapply(states,function(st)matrix(0,Ntimesteps+1L,K));names(hist)<-states
  for(st in states)hist[[st]][1L,]<-comp[[st]]
  death_hist<-matrix(0,Ntimesteps,K);new_hist<-matrix(0,Ntimesteps,K);intro_hist<-matrix(0,Ntimesteps,K);export_hazard<-numeric(Ntimesteps)
  vh<-Vertebrate;if(is.null(vh))vh<-list();vh$Interaction<-NULL
  for(tt in seq_len(Ntimesteps)){
    total<-Reduce(`+`,comp)
    st<-.ivps_step_builder(tt,total,to,ag,Ntimesteps,Nstages,Transition,SDDkernel,LDDkernel,LDDrate,
      PropaguleEstablishment,EnvEstabProb,TransitionKernels,TransitionHabitatSearch,ApplyHabitatToTransitions,
      TransitionEstablishment,BlockedTransitionMortality,HabitatSuitability,HabitatSearchRadius,HabitatSearchCandidates,
      MortalityProb,MortalitySpatial,ManagementExposure,FecundityReduction,FecundityReductionSpatial,
      SpreadReduction,SpreadReductionSpatial,SpreadReductionAppliesTo,OutsideEstablishmentProb,vh,StateColumns,
      DefaultOffspringState,KernelSamples,PointSeed)
    dj<-.ivp_point_joint_disease_step(comp,st$type_table,Pathogen,tt,st$Parent,st$expected_recruits,st$Recruit,
      VerticalTransmissionProb,VerticalRecruitOperator,VerticalSourceState,VerticalTargetState)
    # Expected active export before Interaction; use exact source decomposition for existing hosts.
    active0<-comp$I+if("E"%in%states)comp$E else 0
    export_hazard[tt]<-sum(active0*st$export_transition_by_source)
    comp<-dj$State;for(s in states)hist[[s]][tt+1L,]<-comp[[s]]
    death_hist[tt,]<-dj$PathogenDeaths;new_hist[tt,]<-dj$NewInfections;intro_hist[tt,]<-dj$Introduced
  }
  total_hist<-Reduce(`+`,hist)
  list(Model="INApestVertebratePoint",PathogenModel=.ivp_spec(Pathogen)$Model,
    StateTrajectory=hist,HostTrajectory=rowSums(total_hist),HostStateTrajectory=total_hist,
    ActivePathogenTrajectory=if("E"%in%states)rowSums(hist$E+hist$I)else rowSums(hist$I),
    EndHostPopulation=sum(total_hist[Ntimesteps+1L,]),
    EndActivePathogen=if("E"%in%states)sum(hist$E[Ntimesteps+1L,]+hist$I[Ntimesteps+1L,])else sum(hist$I[Ntimesteps+1L,]),
    ExpectedPathogenDeathsByTimestep=rowSums(death_hist),ExpectedNewInfectionsByTimestep=rowSums(new_hist),
    ExpectedIntroductionsByTimestep=rowSums(intro_hist),EscapeHazardApproximation=1-exp(-sum(export_hazard)),
    TypeTable=to$types,Method="finite-prevalence deterministic first-moment / independent-contact closure",
    Diagnostics=c("Unlike the rare-pathogen operator, this trajectory feeds pathogen mortality back into the total host census before the next vertebrate demographic/control step.",
      "Host demographic/control rates are assumed to depend on total host type abundance and persistent host attributes, not directly on pathogen state. Disease-state-dependent control, movement or fecundity requires a state-specific custom analytical map.",
      "Point horizontal infection uses the exact independent source-target one-source infection probabilities, extended to fractional expected infectious counts by the standard product-PGF closure.",
      "Finite LocalK/KRadius and arbitrary social Interaction remain stochastic when material; they are not silently converted to ordinary density dependence."))
}

.ivp_node_initial_full<-function(InitialPopulation,Model,InitialInfected=NULL,InitialExposed=NULL,InitialRecovered=NULL){
  N<-as.matrix(InitialPopulation);n<-nrow(N);S<-ncol(N);I<-if(is.null(InitialInfected))matrix(0,n,S)else .ina_tm_ns_matrix(InitialInfected,n,S,"InitialInfected")
  E<-if(Model=="SEIR")if(is.null(InitialExposed))matrix(0,n,S)else .ina_tm_ns_matrix(InitialExposed,n,S,"InitialExposed")else matrix(0,n,S)
  R<-if(Model!="SIS")if(is.null(InitialRecovered))matrix(0,n,S)else .ina_tm_ns_matrix(InitialRecovered,n,S,"InitialRecovered")else matrix(0,n,S)
  if(any(I+E+R>N+1e-10))stop("Initial pathogen states exceed InitialPopulation")
  out<-list(S=N-I-E-R,I=I);if(Model=="SEIR")out$E<-E;if(Model!="SIS")out$R<-R;out[c("S",if(Model=="SEIR")"E", "I",if(Model!="SIS")"R")]
}

.ivp_node_joint_disease_step<-function(comp,HostPre,Pathogen,timestep,Ntimesteps,H,RecruitTotal,Rop=NULL,StageMixing=NULL,
    VerticalTransmissionProb=0,VerticalRecruitOperator=NULL,VerticalSourceState="I",VerticalTargetState=NULL){
  ps<-.ivp_spec(Pathogen);model<-ps$Model;states<-ps$States;n<-nrow(HostPre);S<-ncol(HostPre);K<-n*S
  vf<-function(M)as.numeric(t(M));mf<-function(v)matrix(v,n,S,byrow=TRUE)
  prev<-lapply(states,function(st)as.numeric(H%*%vf(comp[[st]])));names(prev)<-states
  V<-if(!is.null(VerticalRecruitOperator)){z<-if(is.function(VerticalRecruitOperator))VerticalRecruitOperator(timestep=timestep,HostPre=HostPre,Recruit=Rop)else VerticalRecruitOperator;as.matrix(z)}else if(any(as.numeric(VerticalTransmissionProb)!=0)){if(is.null(Rop))stop("Vertical transmission with custom Birth requires VerticalRecruitOperator");.ivp_vertical_matrix(Rop,VerticalTransmissionProb,K)}else matrix(0,K,K)
  vt<-if(is.null(VerticalTargetState))if(model=="SEIR")"E"else"I"else VerticalTargetState;vs<-VerticalSourceState
  vb<-as.numeric(V%*%vf(comp[[vs]]));rb<-vf(RecruitTotal);sb<-rb-vb;if(any(sb< -1e-8))stop("Vertical births exceed total recruits")
  prev$S<-prev$S+pmax(0,sb);prev[[vt]]<-prev[[vt]]+vb
  live<-Reduce(`+`,prev);I<-prev$I
  beta<-vf(.ivp_node_resolve(ps$Beta,timestep,Ntimesteps,n,S,"Beta"));rec<-vf(.ivp_node_resolve(ps$RecoveryProb,timestep,Ntimesteps,n,S,"RecoveryProb",TRUE));mort<-vf(.ivp_node_resolve(ps$PathogenMortalityProb,timestep,Ntimesteps,n,S,"PathogenMortalityProb",TRUE));prog<-vf(.ivp_node_resolve(ps$ProgressionProb,timestep,Ntimesteps,n,S,"ProgressionProb",TRUE));wan<-vf(.ivp_node_resolve(ps$ImmunityLossProb,timestep,Ntimesteps,n,S,"ImmunityLossProb",TRUE));ds<-vf(.ivp_node_resolve(ps$DensityScale,timestep,Ntimesteps,n,S,"DensityScale"))
  Cn<-.ivp_node_contact(ps$ContactMatrix,timestep,Ntimesteps,n);Sm<-if(is.null(StageMixing))matrix(1,S,S)else as.matrix(StageMixing);C<-kronecker(Cn,Sm)
  press<-as.numeric(crossprod(I,C));if(ps$Transmission=="frequency"){den<-as.numeric(crossprod(live,C));foi<-beta*ifelse(den>0,press/den,0)}else foi<-beta*press/ds
  pinf<-pmin(1,pmax(0,-expm1(-pmax(0,foi))));new<-prev$S*pinf
  ip<-vf(.ivp_node_resolve(ps$IntroductionProb,timestep,Ntimesteps,n,S,"IntroductionProb",TRUE));intro_n<-vf(.ivp_node_resolve(ps$IntroductionNumber,timestep,Ntimesteps,n,S,"IntroductionNumber"));intro<-ip*pmin(pmax(0,prev$S-new),intro_n)
  stay<-prev$I*(1-rec-mort);recover<-prev$I*rec;deaths<-prev$I*mort
  if(model=="SIS")out<-list(S=prev$S-new-intro+recover,I=stay+new+intro)
  else if(model=="SIR"){lose<-prev$R*wan;out<-list(S=prev$S-new-intro+lose,I=stay+new+intro,R=prev$R-lose+recover)}
  else {pr<-prev$E*prog;lose<-prev$R*wan;out<-list(S=prev$S-new-intro+lose,E=prev$E-pr+new+intro,I=stay+pr,R=prev$R-lose+recover)}
  list(State=lapply(out,mf),NewInfections=mf(new),Introduced=mf(intro),PathogenDeaths=mf(deaths),VerticalBirths=mf(vb))
}

INApestVertebratePathogenMeanNode<-function(Ntimesteps=10L,Nstages,Transition,InitialPopulation,Pathogen,SDDprob,LDDprob=NULL,LDDrate=0,
    TransitionSDDprob=NULL,TransitionLDDprob=NULL,TransitionLDDrate=0,PropaguleEstablishment=1,EnvEstabProb=1,MortalityProb=0,ManagementExposure=0,FecundityReduction=0,OutsideEstablishmentProb=1,Vertebrate=NULL,K=Inf,Weights=rep(1,Nstages),DispersalDensityFactor=0,StageMixing=NULL,VerticalTransmissionProb=0,VerticalRecruitOperator=NULL,VerticalSourceState="I",VerticalTargetState=NULL,InitialInfected=NULL,InitialExposed=NULL,InitialRecovered=NULL){
  ps<-.ivp_spec(Pathogen);model<-ps$Model;if(!model%in%c("SIS","SIR","SEIR"))stop("SIS/SIR/SEIR only");Ntimesteps<-as.integer(Ntimesteps);Nstages<-as.integer(Nstages)
  comp<-.ivp_node_initial_full(InitialPopulation,model,InitialInfected,InitialExposed,InitialRecovered);states<-names(comp);n<-nrow(InitialPopulation);S<-Nstages;U<-n*S
  hist<-lapply(states,function(st)matrix(0,Ntimesteps+1L,U));names(hist)<-states;for(st in states)hist[[st]][1L,]<-as.numeric(t(comp[[st]]))
  deaths<-numeric(Ntimesteps);new<-numeric(Ntimesteps);intro<-numeric(Ntimesteps);exhaz<-numeric(Ntimesteps)
  for(tt in seq_len(Ntimesteps)){
    total<-Reduce(`+`,comp);st<-.ivp_node_parts(total,tt,Ntimesteps,S,Transition,SDDprob,LDDprob,LDDrate,TransitionSDDprob,TransitionLDDprob,TransitionLDDrate,PropaguleEstablishment,EnvEstabProb,MortalityProb,ManagementExposure,FecundityReduction,OutsideEstablishmentProb,Vertebrate,K,Weights,DispersalDensityFactor)
    recruit_total<-st$pre_interaction_host-matrix(as.numeric(st$Parent%*%as.numeric(t(total))),n,S,byrow=TRUE)
    dj<-.ivp_node_joint_disease_step(comp,st$pre_interaction_host,Pathogen,tt,Ntimesteps,st$Parent,recruit_total,st$Recruit,StageMixing,VerticalTransmissionProb,VerticalRecruitOperator,VerticalSourceState,VerticalTargetState)
    active0<-comp$I+if("E"%in%states)comp$E else 0;exhaz[tt]<-sum(as.numeric(t(active0))*st$export_transition_by_source)
    comp<-dj$State;for(ss in states)hist[[ss]][tt+1L,]<-as.numeric(t(comp[[ss]]));deaths[tt]<-sum(dj$PathogenDeaths);new[tt]<-sum(dj$NewInfections);intro[tt]<-sum(dj$Introduced)
  }
  totalhist<-Reduce(`+`,hist);list(Model="INApestVertebrateNode",PathogenModel=model,StateTrajectory=hist,HostTrajectory=rowSums(totalhist),HostStateTrajectory=totalhist,ActivePathogenTrajectory=if("E"%in%states)rowSums(hist$E+hist$I)else rowSums(hist$I),EndHostPopulation=sum(totalhist[Ntimesteps+1L,]),EndActivePathogen=if("E"%in%states)sum(hist$E[Ntimesteps+1L,]+hist$I[Ntimesteps+1L,])else sum(hist$I[Ntimesteps+1L,]),ExpectedPathogenDeathsByTimestep=deaths,ExpectedNewInfectionsByTimestep=new,ExpectedIntroductionsByTimestep=intro,EscapeHazardApproximation=1-exp(-sum(exhaz)),Method="finite-prevalence deterministic first-moment / mass-action closure",Diagnostics=c("Pathogen mortality feeds back into the next timestep's vertebrate host census and therefore into density-dependent host mechanisms.","Host control, HomeRange, birth and movement are assumed pathogen-state independent conditional on total host type abundance. State-specific disease effects beyond PathogenMortalityProb require a custom joint map.","Arbitrary social Interaction and shared pathogen-detection information feedback remain stochastic/non-closed at finite prevalence."))
}

# Wrap rare functions so the default user result also carries the finite-mean
# trajectory, while retaining the rare invasion/branching operator as the
# primary threshold calculation.
INApestVertebratePathogenAnalyticalPoint_rare <- INApestVertebratePathogenAnalyticalPoint
INApestVertebratePathogenAnalyticalPoint <- function(..., FinitePrevalence=TRUE, InitialRecoveredByType=NULL){
  aa<-list(...);r<-do.call(INApestVertebratePathogenAnalyticalPoint_rare,aa)
  if(isTRUE(FinitePrevalence)){
    fm<-names(formals(INApestVertebratePathogenMeanPoint));ma<-aa[intersect(names(aa),fm)];ma$InitialRecoveredByType<-InitialRecoveredByType
    r$FinitePrevalenceMean<-do.call(INApestVertebratePathogenMeanPoint,ma)
    r$EndHostPopulation<-r$FinitePrevalenceMean$EndHostPopulation
    r$EndActivePathogenFiniteMean<-r$FinitePrevalenceMean$EndActivePathogen
  }
  r
}

INApestVertebratePathogenAnalyticalNode_rare <- INApestVertebratePathogenAnalyticalNode
INApestVertebratePathogenAnalyticalNode <- function(..., FinitePrevalence=TRUE, InitialRecovered=NULL){
  aa<-list(...);r<-do.call(INApestVertebratePathogenAnalyticalNode_rare,aa)
  if(isTRUE(FinitePrevalence)){
    fm<-names(formals(INApestVertebratePathogenMeanNode));ma<-aa[intersect(names(aa),fm)];ma$InitialRecovered<-InitialRecovered
    r$FinitePrevalenceMean<-do.call(INApestVertebratePathogenMeanNode,ma)
    r$EndHostPopulation<-r$FinitePrevalenceMean$EndHostPopulation
    r$EndActivePathogenFiniteMean<-r$FinitePrevalenceMean$EndActivePathogen
  }
  r
}

# Final dispatcher, aware of the finite-prevalence wrapper arguments.
INApestAnalytical_pre_vertebrate_pathogen_finite <- INApestAnalytical
INApestAnalytical <- function(...) {
  args<-list(...);Model<-if(!is.null(args$Model))as.character(args$Model)[1L]else"INApest";Pathogen<-args$Pathogen
  if(is.null(Pathogen)||!Model%in%c("INApestVertebrateNode","INApestVertebratePoint"))return(do.call(INApestAnalytical_pre_vertebrate_pathogen_finite,args))
  if(Model=="INApestVertebrateNode"){
    aa<-args[intersect(names(args),c(names(formals(INApestVertebratePathogenAnalyticalNode_rare)),"FinitePrevalence","InitialRecovered"))]
    if(is.null(aa$InitialPopulation))aa$InitialPopulation<-args$InitialState;if(is.null(aa$ManagementExposure))aa$ManagementExposure<-args$ManageProb%||%0
    ans<-do.call(INApestVertebratePathogenAnalyticalNode,aa)
  }else{
    aa<-args[intersect(names(args),c(names(formals(INApestVertebratePathogenAnalyticalPoint_rare)),"FinitePrevalence","InitialRecoveredByType"))]
    if(is.null(aa$ManagementExposure))aa$ManagementExposure<-args$ManageProb%||%0
    ans<-do.call(INApestVertebratePathogenAnalyticalPoint,aa)
  }
  ans$HeadlineEstimands<-list(PathogenGrowthRate=ans$PathogenGrowthRate,HostEndPopulation=ans$EndHostPopulation%||%ans$HostEndPopulationDiseaseFree,EndActivePathogen=ans$EndActivePathogenFiniteMean%||%ans$EndActivePathogen,PathogenEscapeProbability=ans$EscapeProbability)
  ans
}
