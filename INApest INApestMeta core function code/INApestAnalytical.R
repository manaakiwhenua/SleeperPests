###############################################################################
### Analytical / semi-analytical approximations for binary INApest
### Mirrors INApest process ordering: management -> survival/eradication ->
### spread -> information update/detection.
###############################################################################

.inapest_recycle <- function(x, n, name) {
  if (length(x) == 1L) return(rep(as.numeric(x), n))
  if (length(x) != n) stop(name, " must have length 1 or number of nodes")
  as.numeric(x)
}

inapest_combine_dispersal <- function(SDDprob, LDDprob = 0) {
  SDDprob <- as.matrix(SDDprob)
  n <- nrow(SDDprob)
  if (ncol(SDDprob) != n) stop("SDDprob must be square")
  if (length(LDDprob) == 1L) {
    if (LDDprob == 0) return(SDDprob)
    stop("Scalar LDDprob is only supported as 0")
  }
  LDDprob <- as.matrix(LDDprob)
  if (!all(dim(LDDprob) == c(n, n))) stop("LDDprob must match SDDprob")
  1 - (1 - SDDprob) * (1 - LDDprob)
}

inapest_edge_prob <- function(SDDprob, LDDprob = 0, EnvEstabProb = 1) {
  D <- inapest_combine_dispersal(SDDprob, LDDprob)
  n <- nrow(D)
  env <- .inapest_recycle(EnvEstabProb, n, "EnvEstabProb")
  P <- sweep(D, 2, env, `*`)
  diag(P) <- 0 # self-dispersal cannot restore a source that died before spread
  P
}

inapest_exogenous_operator <- function(SDDprob, LDDprob = 0,
                                        EnvEstabProb = 1, Survival = 1,
                                        ManageProb = 0,
                                        EradicationProb = 0,
                                        SpreadReduction = 0) {
  P <- inapest_edge_prob(SDDprob, LDDprob, EnvEstabProb)
  n <- nrow(P)
  s <- .inapest_recycle(Survival, n, "Survival")
  a <- .inapest_recycle(ManageProb, n, "ManageProb")
  e <- .inapest_recycle(EradicationProb, n, "EradicationProb")
  r <- .inapest_recycle(SpreadReduction, n, "SpreadReduction")

  # Expected persistence of an occupied informed node.
  u <- s * (1 - a * e)
  # Expected source weight for a successful off-diagonal colonisation.
  # Same Bernoulli management draw controls eradication and spread reduction.
  v <- s * ((1 - a) + a * (1 - e) * (1 - r))

  G <- t(P) %*% diag(v, nrow = n)
  diag(G) <- u
  dimnames(G) <- list(NULL, NULL)
  attr(G, "self_persistence") <- u
  attr(G, "spread_source_weight") <- v
  G
}

inapest_spectral_radius <- function(G) {
  max(Mod(eigen(G, only.values = TRUE)$values))
}

inapest_finite_mean <- function(G, initial, timesteps) {
  x <- as.numeric(initial)
  out <- matrix(NA_real_, nrow = length(x), ncol = timesteps)
  for (tt in seq_len(timesteps)) {
    x <- as.numeric(G %*% x)
    out[, tt] <- x
  }
  out
}

inapest_meanfield_exogenous <- function(SDDprob, LDDprob = 0,
                                         EnvEstabProb = 1, Survival = 1,
                                         ManageProb = 0,
                                         EradicationProb = 0,
                                         SpreadReduction = 0,
                                         initial, timesteps) {
  P <- inapest_edge_prob(SDDprob, LDDprob, EnvEstabProb)
  n <- nrow(P)
  s <- .inapest_recycle(Survival, n, "Survival")
  a <- .inapest_recycle(ManageProb, n, "ManageProb")
  e <- .inapest_recycle(EradicationProb, n, "EradicationProb")
  r <- .inapest_recycle(SpreadReduction, n, "SpreadReduction")
  u <- s * (1 - a * e)
  v <- s * ((1 - a) + a * (1 - e) * (1 - r))

  q <- as.numeric(initial)
  out <- matrix(NA_real_, nrow = n, ncol = timesteps)
  for (tt in seq_len(timesteps)) {
    self <- q * u
    edge <- sweep(P, 1, q * v, `*`)
    diag(edge) <- 0
    no_col <- apply(1 - edge, 2, prod)
    new_col <- 1 - no_col
    q <- self + (1 - self) * new_col
    q <- pmin(1, pmax(0, q))
    out[, tt] <- q
  }
  out
}

inapest_exogenous_extinction <- function(SDDprob, LDDprob = 0,
                                          EnvEstabProb = 1, Survival = 1,
                                          ManageProb = 0,
                                          EradicationProb = 0,
                                          SpreadReduction = 0,
                                          generations = 50,
                                          tolerance = 1e-12) {
  P <- inapest_edge_prob(SDDprob, LDDprob, EnvEstabProb)
  n <- nrow(P)
  s <- .inapest_recycle(Survival, n, "Survival")
  a <- .inapest_recycle(ManageProb, n, "ManageProb")
  e <- .inapest_recycle(EradicationProb, n, "EradicationProb")
  r <- .inapest_recycle(SpreadReduction, n, "SpreadReduction")

  q <- rep(0, n)
  history <- matrix(NA_real_, nrow = n, ncol = generations)
  for (tt in seq_len(generations)) {
    q_old <- q
    q_new <- numeric(n)
    for (i in seq_len(n)) {
      off <- seq_len(n) != i
      nochild0 <- prod(1 - P[i, off] + P[i, off] * q[off])
      p1 <- P[i, off] * (1 - r[i])
      nochild1 <- prod(1 - p1 + p1 * q[off])
      f0 <- (1 - s[i]) + s[i] * q[i] * nochild0
      sm <- s[i] * (1 - e[i])
      f1 <- (1 - sm) + sm * q[i] * nochild1
      q_new[i] <- (1 - a[i]) * f0 + a[i] * f1
    }
    q <- q_new
    history[, tt] <- q
    if (max(abs(q - q_old)) < tolerance && tt < generations) {
      if (tt < generations) history[, (tt + 1):generations] <- q
      break
    }
  }
  list(extinction = q, history = history)
}

inapest_detection_operator <- function(SDDprob, LDDprob = 0,
                                        EnvEstabProb = 1, Survival = 1,
                                        DetectionProb = 0,
                                        ManageProb = 0,
                                        EradicationProb = 0,
                                        SpreadReduction = 0,
                                        SEAM = NULL,
                                        InfoRetentionProb = 1) {
  P <- inapest_edge_prob(SDDprob, LDDprob, EnvEstabProb)
  n <- nrow(P)
  s <- .inapest_recycle(Survival, n, "Survival")
  d <- .inapest_recycle(DetectionProb, n, "DetectionProb")
  a <- .inapest_recycle(ManageProb, n, "ManageProb")
  e <- .inapest_recycle(EradicationProb, n, "EradicationProb")
  r <- .inapest_recycle(SpreadReduction, n, "SpreadReduction")
  ir <- .inapest_recycle(InfoRetentionProb, n, "InfoRetentionProb")

  if (is.null(SEAM)) {
    C <- matrix(0, n, n)
  } else {
    C <- as.matrix(SEAM)
    if (!all(dim(C) == c(n, n))) stop("SEAM must match SDDprob")
    diag(C) <- 0
  }

  u <- s * (1 - a * e)
  v <- s * ((1 - a) + a * (1 - e) * (1 - r))
  G <- matrix(0, 2 * n, 2 * n)
  U <- seq_len(n)
  H <- n + seq_len(n)

  for (i in seq_len(n)) {
    # Uninformed source: no management, self and colonists can be detected at end.
    G[U[i], U[i]] <- G[U[i], U[i]] + s[i] * (1 - d[i])
    G[H[i], U[i]] <- G[H[i], U[i]] + s[i] * d[i]
    for (j in seq_len(n)) if (j != i && P[i, j] > 0) {
      w <- s[i] * P[i, j]
      G[U[j], U[i]] <- G[U[j], U[i]] + w * (1 - d[j])
      G[H[j], U[i]] <- G[H[j], U[i]] + w * d[j]
    }

    # Informed source: management applies probabilistically.
    # At self, information may decay after spread, then detection can refresh it.
    hself <- ir[i] + (1 - ir[i]) * d[i]
    G[U[i], H[i]] <- G[U[i], H[i]] + u[i] * (1 - hself)
    G[H[i], H[i]] <- G[H[i], H[i]] + u[i] * hself
    for (j in seq_len(n)) if (j != i && P[i, j] > 0) {
      w <- v[i] * P[i, j]
      # Direct-child SEAM approximation: if this H source colonises j, the
      # same source can also communicate information to j before detection.
      hchild <- 1 - (1 - d[j]) * (1 - C[i, j])
      G[U[j], H[i]] <- G[U[j], H[i]] + w * (1 - hchild)
      G[H[j], H[i]] <- G[H[j], H[i]] + w * hchild
    }
  }
  attr(G, "note") <- "2N low-density operator; SEAM captures direct source-to-child information only, not information-only preconditioning of empty nodes"
  G
}

inapest_detection_extinction <- function(SDDprob, LDDprob = 0,
                                          EnvEstabProb = 1, Survival = 1,
                                          DetectionProb = 0,
                                          ManageProb = 0,
                                          EradicationProb = 0,
                                          SpreadReduction = 0,
                                          SEAM = NULL,
                                          InfoRetentionProb = 1,
                                          generations = 50) {
  P <- inapest_edge_prob(SDDprob, LDDprob, EnvEstabProb)
  n <- nrow(P)
  s <- .inapest_recycle(Survival, n, "Survival")
  d <- .inapest_recycle(DetectionProb, n, "DetectionProb")
  a <- .inapest_recycle(ManageProb, n, "ManageProb")
  e <- .inapest_recycle(EradicationProb, n, "EradicationProb")
  r <- .inapest_recycle(SpreadReduction, n, "SpreadReduction")
  ir <- .inapest_recycle(InfoRetentionProb, n, "InfoRetentionProb")
  if (is.null(SEAM)) C <- matrix(0, n, n) else {
    C <- as.matrix(SEAM); diag(C) <- 0
  }

  qU <- rep(0, n); qH <- rep(0, n)
  hU <- matrix(NA_real_, n, generations); hH <- hU
  for (tt in seq_len(generations)) {
    nU <- numeric(n); nH <- numeric(n)
    for (i in seq_len(n)) {
      off <- seq_len(n) != i
      childU <- (1 - d) * qU + d * qH
      prodU <- prod(1 - P[i, off] + P[i, off] * childU[off])
      selfU <- (1 - d[i]) * qU[i] + d[i] * qH[i]
      nU[i] <- (1 - s[i]) + s[i] * selfU * prodU

      hself <- ir[i] + (1 - ir[i]) * d[i]
      selfH <- (1 - hself) * qU[i] + hself * qH[i]
      hchild <- 1 - (1 - d) * (1 - C[i, ])
      childH <- (1 - hchild) * qU + hchild * qH
      prodH0 <- prod(1 - P[i, off] + P[i, off] * childH[off])
      Pm <- P[i, ] * (1 - r[i])
      prodH1 <- prod(1 - Pm[off] + Pm[off] * childH[off])
      f0 <- (1 - s[i]) + s[i] * selfH * prodH0
      sm <- s[i] * (1 - e[i])
      f1 <- (1 - sm) + sm * selfH * prodH1
      nH[i] <- (1 - a[i]) * f0 + a[i] * f1
    }
    qU <- nU; qH <- nH
    hU[, tt] <- qU; hH[, tt] <- qH
  }
  list(U = qU, H = qH, historyU = hU, historyH = hH,
       note = "SEAM extinction uses direct source-to-child approximation only")
}

inapest_temporal_cycle <- function(operators) {
  if (!is.list(operators) || length(operators) < 1) stop("operators must be a non-empty list")
  n <- nrow(operators[[1]])
  B <- diag(n)
  for (G in operators) {
    if (!all(dim(G) == c(n, n))) stop("all operators must have same dimensions")
    B <- G %*% B
  }
  rho_cycle <- inapest_spectral_radius(B)
  list(CycleOperator = B,
       CycleMultiplier = rho_cycle,
       GeometricPerTimestepMultiplier = rho_cycle^(1 / length(operators)))
}

# Expected export intensity requires explicit source -> outside probabilities.
# They cannot in general be reconstructed as 1-rowSums(SDD/LDD), because
# INApest treats each destination edge as an independent Bernoulli opportunity.
inapest_expected_exports <- function(state, ExportProb, source_weight = 1) {
  n <- length(state)
  w <- .inapest_recycle(source_weight, n, "source_weight")
  X <- as.matrix(ExportProb)
  if (nrow(X) != n) stop("ExportProb rows must equal number of internal nodes")
  sum(sweep(X, 1, state * w, `*`))
}

# Nonlinear four-state mean-field approximation for information-limited INApest.
# States at the start of each timestep are U = invaded/uninformed,
# H = invaded/informed, I = uninvaded/informed, E = neither.
# This captures information-only SEAM preconditioning, which is second-order
# around the pest-free state and therefore absent from the 2N linear operator.
inapest_meanfield_information <- function(SDDprob, LDDprob = 0,
                                           EnvEstabProb = 1, Survival = 1,
                                           DetectionProb = 0,
                                           ManageProb = 0,
                                           EradicationProb = 0,
                                           SpreadReduction = 0,
                                           SEAM = NULL,
                                           InfoRetentionProb = 1,
                                           initialU, initialH, initialI = NULL,
                                           timesteps) {
  P <- inapest_edge_prob(SDDprob, LDDprob, EnvEstabProb)
  n <- nrow(P)
  s <- .inapest_recycle(Survival, n, "Survival")
  d <- .inapest_recycle(DetectionProb, n, "DetectionProb")
  a <- .inapest_recycle(ManageProb, n, "ManageProb")
  e <- .inapest_recycle(EradicationProb, n, "EradicationProb")
  r <- .inapest_recycle(SpreadReduction, n, "SpreadReduction")
  ir <- .inapest_recycle(InfoRetentionProb, n, "InfoRetentionProb")
  if (is.null(initialI)) initialI <- rep(0, n)
  U <- as.numeric(initialU); H <- as.numeric(initialH); I <- as.numeric(initialI)
  if (any(U + H + I > 1 + 1e-12)) stop("Initial state probabilities exceed 1")
  if (is.null(SEAM)) C <- matrix(0, n, n) else {
    C <- as.matrix(SEAM)
    if (!all(dim(C) == c(n, n))) stop("SEAM must match SDDprob")
    diag(C) <- 0
  }
  u <- s * (1 - a * e)
  v <- s * ((1 - a) + a * (1 - e) * (1 - r))

  out <- data.frame(timestep=seq_len(timesteps), invaded=NA_real_, informed=NA_real_, managing=NA_real_)
  stateU <- matrix(NA_real_, n, timesteps); stateH <- stateU; stateI <- stateU
  for (tt in seq_len(timesteps)) {
    expected_managing <- sum(a * (H + I))
    SU <- U * s
    SH <- H * u
    pest_surv <- SU + SH

    # Marginal probability that each source generates a successful pest edge.
    source_edge_weight <- U * s + H * v
    edge <- sweep(P, 1, source_edge_weight, `*`)
    diag(edge) <- 0
    new_col <- 1 - apply(1 - edge, 2, prod)
    pest <- pest_surv + (1 - pest_surv) * new_col

    # Information already present is retained after management/spread.
    info_base <- (H + I) * ir

    # SEAM sources are informed pest populations that survive local mortality.
    info_edge <- sweep(C, 1, SH, `*`)
    diag(info_edge) <- 0
    seam_receive <- 1 - apply(1 - info_edge, 2, prod)
    info_before_detection <- info_base + (1 - info_base) * seam_receive

    # Approximate pest/info independence at a destination before direct detection.
    hprob_given_pest <- info_before_detection + (1 - info_before_detection) * d
    Hnew <- pest * hprob_given_pest
    Unew <- pest - Hnew
    Inew <- (1 - pest) * info_before_detection

    U <- pmax(0, pmin(1, Unew)); H <- pmax(0, pmin(1, Hnew)); I <- pmax(0, pmin(1, Inew))
    stateU[,tt] <- U; stateH[,tt] <- H; stateI[,tt] <- I
    out$invaded[tt] <- sum(U + H)
    out$informed[tt] <- sum(H + I)
    out$managing[tt] <- expected_managing
  }
  list(summary=out, U=stateU, H=stateH, I=stateI,
       note="Four-state nonlinear mean-field; captures information-only SEAM preconditioning but approximates pest-information correlations")
}

# Branching approximation to probability of at least one successful export to
# explicitly represented outside destinations. ExportProb must be source-node x
# outside-destination probabilities derived from the same movement kernel or an
# extended landscape; it cannot generally be reconstructed from 1-rowSums(D).
inapest_escape_branching <- function(SDDprob, ExportProb, LDDprob = 0,
                                      EnvEstabProb = 1, Survival = 1,
                                      ManageProb = 0,
                                      EradicationProb = 0,
                                      SpreadReduction = 0,
                                      timesteps = 10) {
  P <- inapest_edge_prob(SDDprob, LDDprob, EnvEstabProb)
  n <- nrow(P)
  X <- as.matrix(ExportProb)
  if (nrow(X) != n) stop("ExportProb rows must equal number of internal nodes")
  if (any(!is.finite(X)) || any(X < 0 | X > 1)) stop("ExportProb must be in [0,1]")
  s <- .inapest_recycle(Survival, n, "Survival")
  a <- .inapest_recycle(ManageProb, n, "ManageProb")
  e <- .inapest_recycle(EradicationProb, n, "EradicationProb")
  r <- .inapest_recycle(SpreadReduction, n, "SpreadReduction")

  # g_i(t) = probability that a lineage starting at source i produces no
  # outside colonisation during the next t timesteps (branching approximation).
  g <- rep(1, n)
  history <- matrix(NA_real_, n, timesteps)
  for (tt in seq_len(timesteps)) {
    old <- g; new <- numeric(n)
    for (i in seq_len(n)) {
      off <- seq_len(n) != i
      internal0 <- prod(1 - P[i,off] + P[i,off] * old[off])
      outside0 <- prod(1 - X[i,])
      f0 <- (1 - s[i]) + s[i] * old[i] * internal0 * outside0

      sm <- s[i] * (1 - e[i])
      Pint <- P[i,off] * (1 - r[i])
      Xman <- X[i,] * (1 - r[i])
      internal1 <- prod(1 - Pint + Pint * old[off])
      outside1 <- prod(1 - Xman)
      f1 <- (1 - sm) + sm * old[i] * internal1 * outside1
      new[i] <- (1 - a[i]) * f0 + a[i] * f1
    }
    g <- new; history[,tt] <- g
  }
  list(no_escape=g, escape=1-g, history_escape=1-history)
}

# First-moment export trajectory and Poisson at-least-one approximation.
inapest_export_first_moment <- function(SDDprob, ExportProb, LDDprob = 0,
                                         EnvEstabProb = 1, Survival = 1,
                                         ManageProb = 0,
                                         EradicationProb = 0,
                                         SpreadReduction = 0,
                                         initial, timesteps) {
  G <- inapest_exogenous_operator(SDDprob,LDDprob,EnvEstabProb,Survival,
                                  ManageProb,EradicationProb,SpreadReduction)
  X <- as.matrix(ExportProb); n <- nrow(G)
  if (nrow(X) != n) stop("ExportProb rows must equal number of internal nodes")
  s <- .inapest_recycle(Survival,n,"Survival")
  a <- .inapest_recycle(ManageProb,n,"ManageProb")
  e <- .inapest_recycle(EradicationProb,n,"EradicationProb")
  r <- .inapest_recycle(SpreadReduction,n,"SpreadReduction")
  v <- s * ((1-a) + a*(1-e)*(1-r))
  x <- as.numeric(initial)
  mu <- numeric(timesteps)
  for(tt in seq_len(timesteps)) {
    mu[tt] <- sum(x * v * rowSums(X))
    x <- as.numeric(G %*% x)
  }
  data.frame(timestep=seq_len(timesteps), expected_exports=mu,
             cumulative_expected_exports=cumsum(mu),
             poisson_escape_probability=1-exp(-cumsum(mu)))
}
###############################################################################
### Analytical / semi-analytical approximations for abundance/stage INApest
### families: INApestMeta, INApestMetaTransitionMatrix, and multiple-land-use.
###
### The primary quantity is the low-density one-timestep mean offspring/operator
### matrix G.  rho(G) < 1 indicates decline when rare; rho(G) > 1 indicates
### growth when rare.  Nonlinear mean-field recursions are also supplied where
### the simulator's density dependence can be represented compactly.
###############################################################################

.ina_recycle <- function(x, n, name) {
  if (length(x) == 1L) return(rep(as.numeric(x), n))
  if (length(x) != n) stop(name, " must have length 1 or the required dimension")
  as.numeric(x)
}

.ina_mat <- function(x, n, name, zero_ok = TRUE) {
  if (length(x) == 1L) {
    if ((zero_ok && is.na(x)) || x == 0) return(matrix(0, n, n))
    stop(name, " scalar form is only supported as 0/NA")
  }
  x <- as.matrix(x)
  if (!all(dim(x) == c(n, n))) stop(name, " must be nodes x nodes")
  x
}

ina_spectral_radius <- function(G) max(Mod(eigen(G, only.values = TRUE)$values))

ina_temporal_cycle <- function(operators) {
  if (!is.list(operators) || !length(operators)) stop("operators must be a non-empty list")
  n <- nrow(operators[[1]])
  B <- diag(n)
  for (G in operators) {
    if (!all(dim(G) == c(n, n))) stop("all operators must have the same dimensions")
    B <- G %*% B
  }
  rc <- ina_spectral_radius(B)
  list(CycleOperator = B, CycleMultiplier = rc,
       GeometricPerTimestepMultiplier = rc^(1 / length(operators)))
}

ina_finite_mean <- function(G, initial, timesteps) {
  x <- as.numeric(initial)
  out <- matrix(NA_real_, nrow = length(x), ncol = timesteps)
  for (tt in seq_len(timesteps)) {
    x <- as.numeric(G %*% x)
    out[, tt] <- x
  }
  out
}

# =============================================================================
# INApestMeta (single population per node)
# =============================================================================

meta_components <- function(SDDprob, LDDprob = 0, LDDrate = 0,
                            EnvEstabProb = 1, Survival = 1, K,
                            PropaguleProduction,
                            PropaguleEstablishment,
                            ManageProb = 0, MortalityProb = 0,
                            SpreadReduction = 0) {
  SDD <- as.matrix(SDDprob)
  n <- nrow(SDD)
  if (ncol(SDD) != n) stop("SDDprob must be square")
  LDD <- .ina_mat(LDDprob, n, "LDDprob")
  r <- as.numeric(LDDrate)
  if (length(r) != 1L || !is.finite(r) || r < 0 || r > 1) stop("LDDrate must be in [0,1]")

  env <- .ina_recycle(EnvEstabProb, n, "EnvEstabProb")
  surv <- .ina_recycle(Survival, n, "Survival")
  cap <- .ina_recycle(K, n, "K")
  prod <- .ina_recycle(PropaguleProduction, n, "PropaguleProduction")
  pest <- .ina_recycle(PropaguleEstablishment, n, "PropaguleEstablishment")
  a <- .ina_recycle(ManageProb, n, "ManageProb")
  m <- .ina_recycle(MortalityProb, n, "MortalityProb")
  g <- .ina_recycle(SpreadReduction, n, "SpreadReduction")

  # Management is a shared node-level Bernoulli draw. Natural/Self-mediated
  # spread is not reduced by SpreadReduction in INApestMeta; LDD is.
  q0 <- surv
  q1 <- surv * (1 - m)
  qbar <- (1 - a) * q0 + a * q1

  K0 <- (1 - r) * SDD + r * LDD
  K1 <- (1 - r) * SDD + r * sweep(LDD, 1, 1 - g, `*`)

  # Expected arrival kernel per original individual (source x destination).
  Arr <- sweep(K0, 1, (1 - a) * q0 * prod, `*`) +
         sweep(K1, 1, a * q1 * prod, `*`)

  # If A~Poisson(lambda), E[1-exp(-alpha*A)] =
  # 1-exp(-lambda*(1-exp(-alpha))).  Thus c=1-exp(-alpha) is the exact
  # low-density slope under the Poisson-arrival approximation.
  alpha <- pest * env
  c_est <- 1 - exp(-alpha)

  list(n = n, SDD = SDD, LDD = LDD, LDDrate = r,
       env = env, survival = surv, K = cap, production = prod,
       prop_est = pest, adoption = a, mortality = m, spread_reduction = g,
       q0 = q0, q1 = q1, qbar = qbar, kernel0 = K0, kernel1 = K1,
       arrivals = Arr, alpha = alpha, recruit_slope = c_est)
}

meta_operator <- function(SDDprob, LDDprob = 0, LDDrate = 0,
                          EnvEstabProb = 1, Survival = 1, K,
                          PropaguleProduction, PropaguleEstablishment,
                          ManageProb = 0, MortalityProb = 0,
                          SpreadReduction = 0) {
  z <- meta_components(SDDprob, LDDprob, LDDrate, EnvEstabProb, Survival, K,
                       PropaguleProduction, PropaguleEstablishment,
                       ManageProb, MortalityProb, SpreadReduction)
  dest_gain <- z$K * z$recruit_slope
  G <- sweep(t(z$arrivals), 1, dest_gain, `*`)
  diag(G) <- diag(G) + z$qbar
  attr(G, "components") <- z
  G
}

meta_meanfield <- function(SDDprob, LDDprob = 0, LDDrate = 0,
                           EnvEstabProb = 1, Survival = 1, K,
                           PropaguleProduction, PropaguleEstablishment,
                           ManageProb = 0, MortalityProb = 0,
                           SpreadReduction = 0,
                           initial, timesteps) {
  z <- meta_components(SDDprob, LDDprob, LDDrate, EnvEstabProb, Survival, K,
                       PropaguleProduction, PropaguleEstablishment,
                       ManageProb, MortalityProb, SpreadReduction)
  x <- as.numeric(initial)
  out <- matrix(NA_real_, z$n, timesteps)
  for (tt in seq_len(timesteps)) {
    n0 <- z$qbar * x
    lambda <- as.numeric(x %*% z$arrivals)
    p_rec <- 1 - exp(-z$recruit_slope * lambda)
    x <- n0 + pmax(0, z$K - n0) * p_rec
    x <- pmin(z$K, pmax(0, x))
    out[, tt] <- x
  }
  out
}

meta_extinction <- function(SDDprob, LDDprob = 0, LDDrate = 0,
                            EnvEstabProb = 1, Survival = 1, K,
                            PropaguleProduction, PropaguleEstablishment,
                            ManageProb = 0, MortalityProb = 0,
                            SpreadReduction = 0,
                            generations = 100, tolerance = 1e-12) {
  z <- meta_components(SDDprob, LDDprob, LDDrate, EnvEstabProb, Survival, K,
                       PropaguleProduction, PropaguleEstablishment,
                       ManageProb, MortalityProb, SpreadReduction)
  n <- z$n
  gain <- z$K * z$recruit_slope
  mu0 <- sweep(z$kernel0, 1, z$production, `*`)
  mu0 <- sweep(mu0, 2, gain, `*`)
  mu1 <- sweep(z$kernel1, 1, z$production, `*`)
  mu1 <- sweep(mu1, 2, gain, `*`)

  q <- rep(0, n)
  hist <- matrix(NA_real_, n, generations)
  for (tt in seq_len(generations)) {
    qo <- q
    qn <- numeric(n)
    for (i in seq_len(n)) {
      R0 <- exp(sum(mu0[i, ] * (q - 1)))
      R1 <- exp(sum(mu1[i, ] * (q - 1)))
      f0 <- (1 - z$q0[i]) + z$q0[i] * q[i] * R0
      f1 <- (1 - z$q1[i]) + z$q1[i] * q[i] * R1
      qn[i] <- (1 - z$adoption[i]) * f0 + z$adoption[i] * f1
    }
    q <- pmin(1, pmax(0, qn)); hist[, tt] <- q
    if (max(abs(q - qo)) < tolerance) {
      if (tt < generations) hist[, (tt + 1):generations] <- q
      break
    }
  }
  list(extinction = q, history = hist,
       note = "Multitype branching approximation; recruits approximated as independent Poisson offspring with simulator-matched low-density mean")
}

meta_detection_operator <- function(SDDprob, LDDprob = 0, LDDrate = 0,
                                    EnvEstabProb = 1, Survival = 1, K,
                                    PropaguleProduction, PropaguleEstablishment,
                                    DetectionProb = 0, ManageProb = 0,
                                    MortalityProb = 0, SpreadReduction = 0,
                                    SEAM = NULL, InfoRetentionProb = 1) {
  n <- nrow(as.matrix(SDDprob))
  d <- .ina_recycle(DetectionProb, n, "DetectionProb")
  ir <- .ina_recycle(InfoRetentionProb, n, "InfoRetentionProb")
  C <- if (is.null(SEAM) || length(SEAM) == 1L) matrix(0, n, n) else as.matrix(SEAM)
  if (!all(dim(C) == c(n, n))) stop("SEAM must be nodes x nodes")
  diag(C) <- 0

  # U: no management; H: expected management.
  G0 <- meta_operator(SDDprob, LDDprob, LDDrate, EnvEstabProb, Survival, K,
                      PropaguleProduction, PropaguleEstablishment,
                      ManageProb = 0, MortalityProb = 0, SpreadReduction = 0)
  GH <- meta_operator(SDDprob, LDDprob, LDDrate, EnvEstabProb, Survival, K,
                      PropaguleProduction, PropaguleEstablishment,
                      ManageProb, MortalityProb, SpreadReduction)
  G <- matrix(0, 2*n, 2*n)
  U <- seq_len(n); H <- n + seq_len(n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      w <- G0[j, i]
      if (w != 0) {
        hj <- d[j]
        G[U[j], U[i]] <- G[U[j], U[i]] + w * (1 - hj)
        G[H[j], U[i]] <- G[H[j], U[i]] + w * hj
      }
      w <- GH[j, i]
      if (w != 0) {
        if (j == i) hj <- ir[j] + (1 - ir[j]) * d[j]
        else hj <- 1 - (1 - d[j]) * (1 - C[i, j])
        G[U[j], H[i]] <- G[U[j], H[i]] + w * (1 - hj)
        G[H[j], H[i]] <- G[H[j], H[i]] + w * hj
      }
    }
  }
  attr(G, "note") <- "2N individual-type low-density approximation; node-level detection/SEAM correlations and information-only preconditioning are omitted"
  G
}

# Expected raw propagules sent to explicitly represented outside destinations.
meta_expected_exports <- function(state, ExportSDDprob, ExportLDDprob = 0,
                                  LDDrate = 0, Survival = 1,
                                  PropaguleProduction,
                                  ManageProb = 0, MortalityProb = 0,
                                  SpreadReduction = 0) {
  x <- as.numeric(state); n <- length(x)
  XS <- as.matrix(ExportSDDprob); if (nrow(XS) != n) stop("ExportSDDprob rows must equal nodes")
  XL <- if (length(ExportLDDprob) == 1L && ExportLDDprob == 0) matrix(0, n, ncol(XS)) else as.matrix(ExportLDDprob)
  if (nrow(XL) != n) stop("ExportLDDprob rows must equal nodes")
  s <- .ina_recycle(Survival,n,"Survival"); p <- .ina_recycle(PropaguleProduction,n,"PropaguleProduction")
  a <- .ina_recycle(ManageProb,n,"ManageProb"); m <- .ina_recycle(MortalityProb,n,"MortalityProb")
  g <- .ina_recycle(SpreadReduction,n,"SpreadReduction"); r <- LDDrate
  q0 <- s; q1 <- s*(1-m)
  e0 <- rowSums((1-r)*XS + r*XL)
  e1 <- rowSums((1-r)*XS + r*sweep(XL,1,1-g,`*`))
  sum(x * p * ((1-a)*q0*e0 + a*q1*e1))
}

# =============================================================================
# INApestMetaTransitionMatrix
# =============================================================================

.transition_list <- function(Transition, n, S) {
  if (is.list(Transition)) {
    if (length(Transition) != n) stop("Transition list must contain one matrix per node")
    A <- lapply(Transition, as.matrix)
  } else {
    Tm <- as.matrix(Transition)
    if (!all(dim(Tm) == c(S,S))) stop("Transition must be S x S")
    A <- replicate(n, Tm, simplify = FALSE)
  }
  if (any(vapply(A, function(x) !all(dim(x)==c(S,S)), logical(1)))) stop("All transition matrices must be S x S")
  A
}

transition_zero_density_sdd <- function(SDDprob, DispersalDensityFactor = 0) {
  S <- as.matrix(SDDprob)
  if (is.na(DispersalDensityFactor) || DispersalDensityFactor == 0) return(S)
  outside <- pmax(0, 1 - rowSums(S))
  tot <- rowSums(S) + outside
  mult <- ifelse(tot > 0, 1/tot, 0)
  sweep(S, 1, mult, `*`)
}

transition_components <- function(Transition, Nstages, SDDprob,
                                  LDDprob = 0, LDDrate = 0,
                                  EnvEstabProb = 1,
                                  PropaguleEstablishment = 1,
                                  ManageProb = 0, MortalityProb = 0,
                                  SpreadReduction = 0,
                                  DispersalDensityFactor = 0,
                                  K = 1, SeedbankK = 1) {
  SDD0 <- transition_zero_density_sdd(SDDprob, DispersalDensityFactor)
  n <- nrow(SDD0); S <- Nstages
  LDD <- .ina_mat(LDDprob,n,"LDDprob")
  A <- .transition_list(Transition,n,S)
  env <- .ina_recycle(EnvEstabProb,n,"EnvEstabProb")
  pe <- .ina_recycle(PropaguleEstablishment,n,"PropaguleEstablishment")
  a <- .ina_recycle(ManageProb,n,"ManageProb")
  g <- .ina_recycle(SpreadReduction,n,"SpreadReduction")
  cap <- .ina_recycle(K,n,"K"); sb <- .ina_recycle(SeedbankK,n,"SeedbankK")

  if (length(MortalityProb) == 1L) M <- matrix(MortalityProb,n,S)
  else if (length(MortalityProb) == S) M <- matrix(rep(MortalityProb, each=n),n,S)
  else if (is.matrix(MortalityProb) && all(dim(MortalityProb)==c(n,S))) M <- MortalityProb
  else stop("MortalityProb must be scalar, length Nstages, or nodes x Nstages")

  # One-arrival recruitment success into stage 1. Under the simulator's
  # Poisson hazard and recruits<=arrivals cap, a single accessible propagule
  # recruits with probability 1-exp(-EnvEstabProb).  Positive footprint is
  # required for SDD; LDD can search the whole seedbank.
  h <- (1 - exp(-env)) * as.numeric(sb > 0)
  sdd_on <- as.numeric(pe > 0 & cap > 0)
  r <- LDDrate

  list(n=n,S=S,A=A,SDD=SDD0,LDD=LDD,env=env,prop_est=pe,adoption=a,
       mortality=M,spread_reduction=g,K=cap,SeedbankK=sb,recruit_success=h,
       sdd_enabled=sdd_on,LDDrate=r)
}

transition_operator <- function(Transition, Nstages, SDDprob,
                                LDDprob = 0, LDDrate = 0,
                                EnvEstabProb = 1,
                                PropaguleEstablishment = 1,
                                ManageProb = 0, MortalityProb = 0,
                                SpreadReduction = 0,
                                DispersalDensityFactor = 0,
                                K = 1, SeedbankK = 1) {
  z <- transition_components(Transition,Nstages,SDDprob,LDDprob,LDDrate,
                             EnvEstabProb,PropaguleEstablishment,ManageProb,
                             MortalityProb,SpreadReduction,DispersalDensityFactor,
                             K,SeedbankK)
  n<-z$n; S<-z$S; G<-matrix(0,n*S,n*S)
  idx <- function(i,s) (i-1L)*S+s
  for(i in seq_len(n)) {
    Ai<-z$A[[i]]; ai<-z$adoption[i]; gi<-z$spread_reduction[i]
    for(k in seq_len(S)) {
      qbar <- 1 - ai*z$mortality[i,k]
      src<-idx(i,k)
      # Simulator uses only stasis, adjacent progression, and terminal stasis.
      if(k < S) {
        G[idx(i,k),src] <- G[idx(i,k),src] + qbar*Ai[k,k]
        G[idx(i,k+1L),src] <- G[idx(i,k+1L),src] + qbar*Ai[k+1L,k]
      } else {
        G[idx(i,S),src] <- G[idx(i,S),src] + qbar*Ai[S,S]
      }
      if(k >= 2L && Ai[1,k] > 0) {
        f<-Ai[1,k]
        q0<-1; q1<-1-z$mortality[i,k]
        for(j in seq_len(n)) {
          nat <- (1-z$LDDrate)*z$sdd_enabled[i]*z$SDD[i,j]
          hum0 <- z$LDDrate*z$LDD[i,j]
          hum1 <- z$LDDrate*(1-gi)*z$LDD[i,j]
          w <- f*z$recruit_success[j]*((1-ai)*q0*(nat+hum0)+ai*q1*(nat+hum1))
          G[idx(j,1L),src] <- G[idx(j,1L),src] + w
        }
      }
    }
  }
  attr(G,"components")<-z
  G
}

transition_intrinsic_managed_lambda <- function(Transition, MortalityProb = 0,
                                                ManageProb = 1) {
  A <- as.matrix(Transition); S <- nrow(A)
  m <- if(length(MortalityProb)==1) rep(MortalityProb,S) else as.numeric(MortalityProb)
  if(length(m)!=S) stop("MortalityProb must be scalar or one value per stage")
  surv <- 1 - ManageProb*m
  M <- A %*% diag(surv)
  list(ManagedTransition=M, Lambda=max(Mod(eigen(M,only.values=TRUE)$values)))
}

transition_extinction <- function(Transition, Nstages, SDDprob,
                                  LDDprob = 0, LDDrate = 0,
                                  EnvEstabProb = 1,
                                  PropaguleEstablishment = 1,
                                  ManageProb = 0, MortalityProb = 0,
                                  SpreadReduction = 0,
                                  DispersalDensityFactor = 0,
                                  K = 1, SeedbankK = 1,
                                  generations = 100) {
  z <- transition_components(Transition,Nstages,SDDprob,LDDprob,LDDrate,
                             EnvEstabProb,PropaguleEstablishment,ManageProb,
                             MortalityProb,SpreadReduction,DispersalDensityFactor,
                             K,SeedbankK)
  n<-z$n; S<-z$S; nt<-n*S; idx<-function(i,s)(i-1L)*S+s
  q<-rep(0,nt); hist<-matrix(NA_real_,nt,generations)
  for(tt in seq_len(generations)) {
    qn<-numeric(nt)
    for(i in seq_len(n)) {
      Ai<-z$A[[i]]; ai<-z$adoption[i]; gi<-z$spread_reduction[i]
      for(k in seq_len(S)) {
        src<-idx(i,k); f<-if(k>=2) Ai[1,k] else 0
        calcM<-function(M) {
          surv<-1-z$mortality[i,k]*M
          if(k<S) local <- (1-Ai[k,k]-Ai[k+1,k]) + Ai[k,k]*q[idx(i,k)] + Ai[k+1,k]*q[idx(i,k+1L)]
          else local <- (1-Ai[S,S]) + Ai[S,S]*q[idx(i,S)]
          muq<-0
          if(f>0) for(j in seq_len(n)) {
            nat<-(1-z$LDDrate)*z$sdd_enabled[i]*z$SDD[i,j]
            hum<-z$LDDrate*(1-gi*M)*z$LDD[i,j]
            mu<-f*(nat+hum)*z$recruit_success[j]
            muq<-muq+mu*(q[idx(j,1L)]-1)
          }
          (1-surv)+surv*local*exp(muq)
        }
        qn[src]<-(1-ai)*calcM(0)+ai*calcM(1)
      }
    }
    q<-pmin(1,pmax(0,qn)); hist[,tt]<-q
  }
  list(extinction=q,history=hist,
       note="Branching approximation; exact local stasis/progression structure and Poisson fecundity, approximate independent successful recruits")
}

transition_detection_operator <- function(Transition, Nstages, SDDprob,
                                           LDDprob = 0, LDDrate = 0,
                                           EnvEstabProb = 1,
                                           PropaguleEstablishment = 1,
                                           DetectionProb = 0,
                                           ManageProb = 0, MortalityProb = 0,
                                           SpreadReduction = 0,
                                           SEAM = NULL, InfoRetentionProb = 1,
                                           DispersalDensityFactor = 0,
                                           K = 1, SeedbankK = 1) {
  base <- transition_operator(Transition,Nstages,SDDprob,LDDprob,LDDrate,
                              EnvEstabProb,PropaguleEstablishment,0,0,0,
                              DispersalDensityFactor,K,SeedbankK)
  managed <- transition_operator(Transition,Nstages,SDDprob,LDDprob,LDDrate,
                                 EnvEstabProb,PropaguleEstablishment,ManageProb,
                                 MortalityProb,SpreadReduction,
                                 DispersalDensityFactor,K,SeedbankK)
  z<-attr(managed,"components"); n<-z$n; S<-z$S; nt<-n*S
  if(length(DetectionProb)==1) D<-matrix(DetectionProb,n,S)
  else if(length(DetectionProb)==S) D<-matrix(rep(DetectionProb,each=n),n,S)
  else if(is.matrix(DetectionProb)&&all(dim(DetectionProb)==c(n,S))) D<-DetectionProb
  else stop("DetectionProb must be scalar, length Nstages, or nodes x Nstages")
  ir<-.ina_recycle(InfoRetentionProb,n,"InfoRetentionProb")
  C<-if(is.null(SEAM)||length(SEAM)==1L) matrix(0,n,n) else as.matrix(SEAM); diag(C)<-0
  idx<-function(i,s)(i-1L)*S+s
  G<-matrix(0,2*nt,2*nt); U<-seq_len(nt); H<-nt+seq_len(nt)
  for(i in seq_len(n)) for(k in seq_len(S)) {
    src<-idx(i,k)
    for(j in seq_len(n)) for(s in seq_len(S)) {
      dst<-idx(j,s)
      w<-base[dst,src]
      if(w!=0){h<-D[j,s];G[U[dst],U[src]]<-G[U[dst],U[src]]+w*(1-h);G[H[dst],U[src]]<-G[H[dst],U[src]]+w*h}
      w<-managed[dst,src]
      if(w!=0){
        h<-if(j==i) ir[j]+(1-ir[j])*D[j,s] else 1-(1-D[j,s])*(1-C[i,j])
        G[U[dst],H[src]]<-G[U[dst],H[src]]+w*(1-h);G[H[dst],H[src]]<-G[H[dst],H[src]]+w*h
      }
    }
  }
  attr(G,"note")<-"2NS low-density individual-type approximation; shared node information creates correlations not represented here"
  G
}

# =============================================================================
# INApestMetaMultipleLandUse
# =============================================================================

.mlu_matrix <- function(x,n,L,name) {
  if(length(x)==1L) return(matrix(x,n,L))
  if(length(x)==L) return(matrix(rep(x,each=n),n,L))
  if(is.matrix(x)&&all(dim(x)==c(n,L))) return(x)
  stop(name," must be scalar, length Nlanduses, or nodes x Nlanduses")
}

mlu_components <- function(SDDprob,LDDprob=0,LDDrate=0,EnvEstabProb=1,
                           Survival=1,K,PropaguleProduction,
                           PropaguleEstablishment,
                           ManageProb=0,MortalityProb=0,SpreadReduction=0,
                           current_code=TRUE) {
  SDD<-as.matrix(SDDprob); n<-nrow(SDD); LDD<-.ina_mat(LDDprob,n,"LDDprob")
  K<-as.matrix(K); if(nrow(K)!=n) stop("K must have one row per node"); L<-ncol(K)
  s<-.ina_recycle(Survival,n,"Survival"); p<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction")
  env<-.ina_recycle(EnvEstabProb,n,"EnvEstabProb"); pe<-.ina_recycle(PropaguleEstablishment,n,"PropaguleEstablishment")
  A<-.mlu_matrix(ManageProb,n,L,"ManageProb"); M<-.mlu_matrix(MortalityProb,n,L,"MortalityProb")
  Gm<-.mlu_matrix(SpreadReduction,n,L,"SpreadReduction")
  qbar<-matrix(s,n,L)*(1-A*M)
  r<-LDDrate; alpha<-env*pe; c_est<-1-exp(-alpha)

  # Emission[ i, l, j ] = expected arrival count at destination j per one
  # original individual in source node i / land use l.
  Em<-array(0,c(n,L,n))
  for(i in seq_len(n)) for(l in seq_len(L)) {
    natcoef<-qbar[i,l]*(1-r)
    if(current_code) {
      # Match current local.dynamicsLU: Qout is calculated as
      # Propagules*r*sum_l(1-g_il*M_il) before the later Pn-weighted Qout is
      # computed (but not used). Empty land-use management draws therefore
      # affect LDD and the factor can scale with Nlanduses.
      selfterm<-s[i]*((1-A[i,l])+A[i,l]*(1-M[i,l])*(1-Gm[i,l]))
      other<-0
      if(L>1) for(h in setdiff(seq_len(L),l)) other<-other+qbar[i,l]*(1-A[i,h]*Gm[i,h])
      lddcoef<-selfterm+other
    } else {
      # Intended composition-weighted one-individual limit from the unused Pn line.
      lddcoef<-s[i]*((1-A[i,l])+A[i,l]*(1-M[i,l])*(1-Gm[i,l]))
    }
    Em[i,l,]<-p[i]*(natcoef*SDD[i,] + r*lddcoef*LDD[i,])
  }
  list(n=n,L=L,SDD=SDD,LDD=LDD,K=K,Ktot=rowSums(K),survival=s,production=p,
       env=env,prop_est=pe,adoption=A,mortality=M,spread_reduction=Gm,
       qbar=qbar,c_est=c_est,emission=Em,current_code=current_code)
}

mlu_operator <- function(SDDprob,LDDprob=0,LDDrate=0,EnvEstabProb=1,
                         Survival=1,K,PropaguleProduction,PropaguleEstablishment,
                         ManageProb=0,MortalityProb=0,SpreadReduction=0,
                         current_code=TRUE) {
  z<-mlu_components(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                    PropaguleProduction,PropaguleEstablishment,ManageProb,
                    MortalityProb,SpreadReduction,current_code)
  n<-z$n;L<-z$L;idx<-function(i,l)(i-1L)*L+l;G<-matrix(0,n*L,n*L)
  for(i in seq_len(n)) for(l in seq_len(L)) {
    src<-idx(i,l);G[src,src]<-G[src,src]+z$qbar[i,l]
    for(j in seq_len(n)) for(h in seq_len(L))
      G[idx(j,h),src]<-G[idx(j,h),src]+z$K[j,h]*z$c_est[j]*z$emission[i,l,j]
  }
  attr(G,"components")<-z;G
}

mlu_meanfield <- function(SDDprob,LDDprob=0,LDDrate=0,EnvEstabProb=1,
                          Survival=1,K,PropaguleProduction,PropaguleEstablishment,
                          ManageProb=0,MortalityProb=0,SpreadReduction=0,
                          current_code=TRUE,initial,timesteps) {
  z<-mlu_components(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                    PropaguleProduction,PropaguleEstablishment,ManageProb,
                    MortalityProb,SpreadReduction,current_code)
  X<-as.matrix(initial);if(!all(dim(X)==c(z$n,z$L)))stop("initial must be nodes x landuses")
  out<-array(NA_real_,c(z$n,z$L,timesteps))
  for(tt in seq_len(timesteps)) {
    N0<-z$qbar*X
    lam<-numeric(z$n)
    for(i in seq_len(z$n)) for(l in seq_len(z$L)) lam<-lam+X[i,l]*z$emission[i,l,]
    p_rec<-1-exp(-z$c_est*lam)
    free<-pmax(0,z$K-N0);free_tot<-rowSums(free)
    total_rec<-free_tot*p_rec
    share<-free/free_tot;share[!is.finite(share)]<-0
    X<-N0+share*total_rec
    X<-pmin(z$K,pmax(0,X));out[,,tt]<-X
  }
  out
}

mlu_detection_operator <- function(SDDprob,LDDprob=0,LDDrate=0,EnvEstabProb=1,
                                   Survival=1,K,PropaguleProduction,PropaguleEstablishment,
                                   DetectionProb=0,ManageProb=0,MortalityProb=0,
                                   SpreadReduction=0,SEAM=NULL,InfoRetentionProb=1,
                                   current_code=TRUE) {
  Kmat<-as.matrix(K);n<-nrow(Kmat);L<-ncol(Kmat);nt<-n*L
  D<-.mlu_matrix(DetectionProb,n,L,"DetectionProb")
  ir<-.ina_recycle(InfoRetentionProb,n,"InfoRetentionProb")
  C<-if(is.null(SEAM)||length(SEAM)==1L)matrix(0,n,n)else as.matrix(SEAM);diag(C)<-0
  G0<-mlu_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                   PropaguleProduction,PropaguleEstablishment,0,0,0,current_code)
  GH<-mlu_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                   PropaguleProduction,PropaguleEstablishment,ManageProb,
                   MortalityProb,SpreadReduction,current_code)
  idx<-function(i,l)(i-1L)*L+l;G<-matrix(0,2*nt,2*nt);U<-seq_len(nt);H<-nt+seq_len(nt)
  for(i in seq_len(n))for(l in seq_len(L)){src<-idx(i,l);for(j in seq_len(n))for(h in seq_len(L)){
    dst<-idx(j,h);w<-G0[dst,src];if(w!=0){hh<-D[j,h];G[U[dst],U[src]]<-G[U[dst],U[src]]+w*(1-hh);G[H[dst],U[src]]<-G[H[dst],U[src]]+w*hh}
    w<-GH[dst,src];if(w!=0){hh<-if(j==i)ir[j]+(1-ir[j])*D[j,h]else 1-(1-D[j,h])*(1-C[i,j]);G[U[dst],H[src]]<-G[U[dst],H[src]]+w*(1-hh);G[H[dst],H[src]]<-G[H[dst],H[src]]+w*hh}
  }}
  attr(G,"note")<-"2NL low-density individual-type approximation; node-shared information correlations omitted; current_code=TRUE matches existing LDD Qout implementation"
  G
}

# --- Simulator-faithful rare-lineage operator for INApestMeta -----------------
# The simulator supplies a generally fractional value to rmultinom(size=...),
# which R truncates to an integer.  This creates a genuine small-population
# nonlinearity when dispersal row sums or LDDrate are fractional.  The function
# below calculates the expected one-timestep offspring from ONE original
# individual exactly with respect to Poisson propagule production, multinomial
# allocation, and the binomial recruitment mean.  It is therefore the preferred
# growth-when-rare operator for the current INApestMeta implementation.

.meta_poisson_support <- function(lambda, tail = 1e-12) {
  if (lambda <= 0) return(list(k=0, p=1))
  kmax <- max(20L, as.integer(qpois(1-tail, lambda)))
  k <- 0:kmax; p <- dpois(k,lambda)
  list(k=k,p=p/sum(p))
}

meta_single_parent_operator <- function(SDDprob, LDDprob = 0, LDDrate = 0,
                                        EnvEstabProb = 1, Survival = 1, K,
                                        PropaguleProduction,
                                        PropaguleEstablishment,
                                        ManageProb = 0, MortalityProb = 0,
                                        SpreadReduction = 0) {
  SDD<-as.matrix(SDDprob);n<-nrow(SDD);LDD<-.ina_mat(LDDprob,n,"LDDprob")
  env<-.ina_recycle(EnvEstabProb,n,"EnvEstabProb");s<-.ina_recycle(Survival,n,"Survival")
  cap<-.ina_recycle(K,n,"K");prod<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction")
  pe<-.ina_recycle(PropaguleEstablishment,n,"PropaguleEstablishment")
  a<-.ina_recycle(ManageProb,n,"ManageProb");m<-.ina_recycle(MortalityProb,n,"MortalityProb")
  g<-.ina_recycle(SpreadReduction,n,"SpreadReduction");r<-LDDrate; alpha<-pe*env
  rs<-rowSums(SDD); rl<-rowSums(LDD); G<-matrix(0,n,n)
  for(i in seq_len(n)) {
    ps <- if(rs[i]>0) SDD[i,]/rs[i] else rep(0,n)
    pl <- if(rl[i]>0) LDD[i,]/rl[i] else rep(0,n)
    supp <- .meta_poisson_support(prod[i]); k<-supp$k; pk<-supp$p
    for(M in 0:1) {
      pm <- if(M==0) 1-a[i] else a[i]
      if(pm==0) next
      surv <- s[i]*(1-m[i]*M)
      if(surv==0) next
      ms <- floor(k * ((1-r)*rs[i]))
      ml <- floor(k * (r*(1-g[i]*M)*rl[i]))
      free <- cap; free[i] <- pmax(0,free[i]-1)
      for(j in seq_len(n)) {
        bS <- 1-ps[j]*(1-exp(-alpha[j]))
        bL <- 1-pl[j]*(1-exp(-alpha[j]))
        nohaz <- sum(pk * (bS^ms) * (bL^ml))
        recruits <- free[j] * (1-nohaz)
        G[j,i] <- G[j,i] + pm*surv*(as.numeric(j==i)+recruits)
      }
    }
  }
  attr(G,"note") <- "One-parent mean operator matching current rmultinom integer truncation and binomial recruitment mean; preferred for growth when rare. Multiple-source cooperation through floor(sum(...)) is necessarily omitted."
  G
}

meta_detection_single_parent_operator <- function(SDDprob, LDDprob = 0, LDDrate = 0,
                                                   EnvEstabProb = 1, Survival = 1, K,
                                                   PropaguleProduction, PropaguleEstablishment,
                                                   DetectionProb = 0, ManageProb = 0,
                                                   MortalityProb = 0, SpreadReduction = 0,
                                                   SEAM = NULL, InfoRetentionProb = 1) {
  n<-nrow(as.matrix(SDDprob));d<-.ina_recycle(DetectionProb,n,"DetectionProb");ir<-.ina_recycle(InfoRetentionProb,n,"InfoRetentionProb")
  C<-if(is.null(SEAM)||length(SEAM)==1L)matrix(0,n,n)else as.matrix(SEAM);diag(C)<-0
  G0<-meta_single_parent_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,PropaguleProduction,PropaguleEstablishment,0,0,0)
  GH<-meta_single_parent_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,PropaguleProduction,PropaguleEstablishment,ManageProb,MortalityProb,SpreadReduction)
  G<-matrix(0,2*n,2*n);U<-seq_len(n);H<-n+seq_len(n)
  for(i in seq_len(n))for(j in seq_len(n)){
    w<-G0[j,i];if(w!=0){h<-d[j];G[U[j],U[i]]<-G[U[j],U[i]]+w*(1-h);G[H[j],U[i]]<-G[H[j],U[i]]+w*h}
    w<-GH[j,i];if(w!=0){h<-if(j==i)ir[j]+(1-ir[j])*d[j]else 1-(1-d[j])*(1-C[i,j]);G[U[j],H[i]]<-G[U[j],H[i]]+w*(1-h);G[H[j],H[i]]<-G[H[j],H[i]]+w*h}
  }
  attr(G,"note")<-"2N one-parent operator with simulator integer truncation; information is individual-typed approximation"
  G
}

meta_single_parent_recruit_means <- function(SDDprob, LDDprob = 0, LDDrate = 0,
                                             EnvEstabProb = 1, K,
                                             PropaguleProduction,
                                             PropaguleEstablishment,
                                             SpreadReduction = 0) {
  SDD<-as.matrix(SDDprob);n<-nrow(SDD);LDD<-.ina_mat(LDDprob,n,"LDDprob")
  env<-.ina_recycle(EnvEstabProb,n,"EnvEstabProb");cap<-.ina_recycle(K,n,"K")
  prod<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction");pe<-.ina_recycle(PropaguleEstablishment,n,"PropaguleEstablishment")
  g<-.ina_recycle(SpreadReduction,n,"SpreadReduction");r<-LDDrate;alpha<-pe*env;rs<-rowSums(SDD);rl<-rowSums(LDD)
  out<-list(matrix(0,n,n),matrix(0,n,n))
  for(i in seq_len(n)){
    ps<-if(rs[i]>0)SDD[i,]/rs[i]else rep(0,n);pl<-if(rl[i]>0)LDD[i,]/rl[i]else rep(0,n)
    supp<-.meta_poisson_support(prod[i]);k<-supp$k;pk<-supp$p
    for(M in 0:1){ms<-floor(k*((1-r)*rs[i]));ml<-floor(k*(r*(1-g[i]*M)*rl[i]));free<-cap;free[i]<-pmax(0,free[i]-1)
      for(j in seq_len(n)){bS<-1-ps[j]*(1-exp(-alpha[j]));bL<-1-pl[j]*(1-exp(-alpha[j]));nohaz<-sum(pk*(bS^ms)*(bL^ml));out[[M+1L]][i,j]<-free[j]*(1-nohaz)}
    }
  }
  out
}

meta_single_parent_extinction <- function(SDDprob,LDDprob=0,LDDrate=0,
                                          EnvEstabProb=1,Survival=1,K,
                                          PropaguleProduction,PropaguleEstablishment,
                                          ManageProb=0,MortalityProb=0,SpreadReduction=0,
                                          generations=100,tolerance=1e-12){
  n<-nrow(as.matrix(SDDprob));s<-.ina_recycle(Survival,n,"Survival");a<-.ina_recycle(ManageProb,n,"ManageProb");m<-.ina_recycle(MortalityProb,n,"MortalityProb")
  mus<-meta_single_parent_recruit_means(SDDprob,LDDprob,LDDrate,EnvEstabProb,K,PropaguleProduction,PropaguleEstablishment,SpreadReduction)
  q<-rep(0,n);hist<-matrix(NA_real_,n,generations)
  for(tt in seq_len(generations)){qo<-q;qn<-numeric(n);for(i in seq_len(n)){
    q0<-s[i];q1<-s[i]*(1-m[i]);R0<-exp(sum(mus[[1]][i,]*(q-1)));R1<-exp(sum(mus[[2]][i,]*(q-1)))
    f0<-(1-q0)+q0*q[i]*R0;f1<-(1-q1)+q1*q[i]*R1;qn[i]<-(1-a[i])*f0+a[i]*f1}
    q<-pmin(1,pmax(0,qn));hist[,tt]<-q;if(max(abs(q-qo))<tolerance){if(tt<generations)hist[,(tt+1):generations]<-q;break}}
  list(extinction=q,history=hist,note="Branching PGF with simulator-faithful one-parent recruit means; recruit correlations approximated as independent Poisson")
}

# --- Simulator-faithful rare-lineage operator for Multiple Land Use ----------
.mlu_management_states <- function(a) {
  L<-length(a); if(L>12) stop("Exact current-code MLU one-parent operator currently supports up to 12 land uses")
  M<-as.matrix(expand.grid(rep(list(0:1),L))); if(L==1)M<-matrix(M,ncol=1)
  p<-apply(M,1,function(z)prod(ifelse(z==1,a,1-a)))
  list(M=M,p=p)
}

mlu_single_parent_operator <- function(SDDprob,LDDprob=0,LDDrate=0,EnvEstabProb=1,
                                       Survival=1,K,PropaguleProduction,
                                       PropaguleEstablishment,
                                       ManageProb=0,MortalityProb=0,
                                       SpreadReduction=0,current_code=FALSE) {
  SDD<-as.matrix(SDDprob);n<-nrow(SDD);LDD<-.ina_mat(LDDprob,n,"LDDprob");K<-as.matrix(K);L<-ncol(K)
  s<-.ina_recycle(Survival,n,"Survival");prod<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction")
  env<-.ina_recycle(EnvEstabProb,n,"EnvEstabProb");pe<-.ina_recycle(PropaguleEstablishment,n,"PropaguleEstablishment")
  A<-.mlu_matrix(ManageProb,n,L,"ManageProb");Mort<-.mlu_matrix(MortalityProb,n,L,"MortalityProb");Gr<-.mlu_matrix(SpreadReduction,n,L,"SpreadReduction")
  alpha<-env*pe;r<-LDDrate;rs<-rowSums(SDD);rl<-rowSums(LDD);idx<-function(i,l)(i-1L)*L+l;G<-matrix(0,n*L,n*L)
  for(i in seq_len(n)) {
    ps<-if(rs[i]>0)SDD[i,]/rs[i]else rep(0,n);pl<-if(rl[i]>0)LDD[i,]/rl[i]else rep(0,n);supp<-.meta_poisson_support(prod[i]);k<-supp$k;pk<-supp$p
    st<-.mlu_management_states(A[i,])
    for(l in seq_len(L)) {
      src<-idx(i,l)
      for(z in seq_len(nrow(st$M))) {
        Mv<-st$M[z,];pm<-st$p[z];if(pm==0)next
        surv<-s[i]*(1-Mort[i,l]*Mv[l]);if(surv==0)next
        sf<-if(current_code)sum(1-Gr[i,]*Mv)else(1-Gr[i,l]*Mv[l])
        ms<-floor(k*((1-r)*rs[i]));ml<-floor(k*(r*sf*rl[i]))
        for(j in seq_len(n)) {
          bS<-1-ps[j]*(1-exp(-alpha[j]));bL<-1-pl[j]*(1-exp(-alpha[j]));nohaz<-sum(pk*(bS^ms)*(bL^ml))
          free<-K[j,];if(j==i)free[l]<-pmax(0,free[l]-1)
          for(h in seq_len(L)) G[idx(j,h),src]<-G[idx(j,h),src]+pm*surv*free[h]*(1-nohaz)
        }
        G[src,src]<-G[src,src]+pm*surv
      }
    }
  }
  attr(G,"note")<-if(current_code)"Legacy one-parent MLU operator matching the former unweighted LDD Qout" else "One-parent MLU operator using the population-share-weighted intended LDD Qout"
  G
}

mlu_detection_single_parent_operator <- function(SDDprob,LDDprob=0,LDDrate=0,EnvEstabProb=1,
                                                  Survival=1,K,PropaguleProduction,PropaguleEstablishment,
                                                  DetectionProb=0,ManageProb=0,MortalityProb=0,
                                                  SpreadReduction=0,SEAM=NULL,InfoRetentionProb=1,
                                                  current_code=FALSE) {
  Kmat<-as.matrix(K);n<-nrow(Kmat);L<-ncol(Kmat);nt<-n*L;D<-.mlu_matrix(DetectionProb,n,L,"DetectionProb")
  ir<-.ina_recycle(InfoRetentionProb,n,"InfoRetentionProb");C<-if(is.null(SEAM)||length(SEAM)==1L)matrix(0,n,n)else as.matrix(SEAM);diag(C)<-0
  G0<-mlu_single_parent_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,PropaguleProduction,PropaguleEstablishment,0,0,0,current_code)
  GH<-mlu_single_parent_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,PropaguleProduction,PropaguleEstablishment,ManageProb,MortalityProb,SpreadReduction,current_code)
  idx<-function(i,l)(i-1L)*L+l;G<-matrix(0,2*nt,2*nt);U<-seq_len(nt);H<-nt+seq_len(nt)
  for(i in seq_len(n))for(l in seq_len(L)){src<-idx(i,l);for(j in seq_len(n))for(h in seq_len(L)){dst<-idx(j,h)
    w<-G0[dst,src];if(w!=0){hh<-D[j,h];G[U[dst],U[src]]<-G[U[dst],U[src]]+w*(1-hh);G[H[dst],U[src]]<-G[H[dst],U[src]]+w*hh}
    w<-GH[dst,src];if(w!=0){hh<-if(j==i)ir[j]+(1-ir[j])*D[j,h]else 1-(1-D[j,h])*(1-C[i,j]);G[U[dst],H[src]]<-G[U[dst],H[src]]+w*(1-hh);G[H[dst],H[src]]<-G[H[dst],H[src]]+w*hh}
  }}
  G
}

mlu_single_parent_extinction <- function(..., generations=100,current_code=FALSE){
  dots<-list(...);G<-do.call(mlu_single_parent_operator,c(dots,list(current_code=current_code)))
  nt<-nrow(G); self<-diag(G); rec<-G;diag(rec)<-0
  # Include same-type recruitment in the conditional offspring mean too.
  for(i in seq_len(nt)) rec[i,i]<-pmax(0,G[i,i]-pmin(1,self[i]))
  # Approximate the probability that the original individual survives by the
  # expected self-survivor term, then condition Poisson recruit means on it.
  qsurv<-pmin(1,diag(G)); q<-rep(0,nt);hist<-matrix(NA_real_,nt,generations)
  for(tt in seq_len(generations)){qn<-numeric(nt);for(i in seq_len(nt)){
    mu<-if(qsurv[i]>0)G[,i]/qsurv[i]else rep(0,nt);mu[i]<-pmax(0,mu[i]-1)
    qn[i]<-(1-qsurv[i])+qsurv[i]*q[i]*exp(sum(mu*(q-1)))};q<-pmin(1,pmax(0,qn));hist[,tt]<-q}
  list(extinction=q,history=hist,note="Mean-matched branching approximation based on simulator-faithful MLU one-parent operator")
}

# =============================================================================
# Generic first-moment escape approximation
# =============================================================================

# Given an internal next-generation operator G and expected numbers of raw
# propagules exported per source type/timestep, approximate the probability of
# >=1 successful escape by treating successful exports over time as Poisson.
# This is a first-moment approximation; it is not a full branching no-escape
# recursion.  export_vector is indexed by the same source types as state0.
ina_escape_first_moment <- function(G, state0, export_vector, steps,
                                    OutsideEstablishmentProb = 1) {
  G <- as.matrix(G); x <- as.numeric(state0); e <- as.numeric(export_vector)
  if (nrow(G) != ncol(G) || length(x) != nrow(G) || length(e) != nrow(G))
    stop("G, state0 and export_vector have incompatible dimensions")
  if (length(OutsideEstablishmentProb) == 1L) {
    p_out <- rep(OutsideEstablishmentProb, length(e))
  } else {
    p_out <- as.numeric(OutsideEstablishmentProb)
    if (length(p_out) != length(e)) stop("OutsideEstablishmentProb must be scalar or match export_vector")
  }
  if (any(p_out < 0 | p_out > 1, na.rm=TRUE)) stop("OutsideEstablishmentProb must be in [0,1]")
  by_step <- numeric(steps)
  for (t in seq_len(steps)) {
    by_step[t] <- sum(x * e * p_out)
    x <- as.numeric(G %*% x)
  }
  cum <- cumsum(by_step)
  list(ExpectedSuccessfulEscapesByStep=by_step,
       ExpectedSuccessfulEscapesCumulative=cum,
       EscapeProbabilityByStep=1-exp(-cum),
       EscapeProbability=if(length(cum)) tail(1-exp(-cum),1) else 0)
}

.ina_residual_export <- function(InternalProb, ExternalProb=NULL, name="dispersal") {
  I <- as.matrix(InternalProb); n <- nrow(I)
  if (!is.null(ExternalProb)) {
    E <- as.matrix(ExternalProb)
    if (nrow(E) != n) stop("External probability matrix rows must match internal nodes")
    return(rowSums(E))
  }
  rs <- rowSums(I)
  if (any(rs > 1 + 1e-10))
    warning(name, " rows exceed 1; residual 1-rowSum cannot be interpreted as export. Clamping residual export at zero.")
  pmax(0, 1-rs)
}

meta_export_vector <- function(SDDprob, LDDprob=0, LDDrate=0,
                               Survival=1, PropaguleProduction,
                               ManageProb=0, MortalityProb=0,
                               SpreadReduction=0,
                               ExportSDDprob=NULL, ExportLDDprob=NULL) {
  SDD<-as.matrix(SDDprob); n<-nrow(SDD); LDD<-.ina_mat(LDDprob,n,"LDDprob")
  os<-.ina_residual_export(SDD,ExportSDDprob,"SDDprob")
  ol<-.ina_residual_export(LDD,ExportLDDprob,"LDDprob")
  s<-.ina_recycle(Survival,n,"Survival"); p<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction")
  a<-.ina_recycle(ManageProb,n,"ManageProb"); m<-.ina_recycle(MortalityProb,n,"MortalityProb")
  g<-.ina_recycle(SpreadReduction,n,"SpreadReduction"); r<-LDDrate
  p*((1-a)*s*((1-r)*os+r*ol) + a*s*(1-m)*((1-r)*os+r*(1-g)*ol))
}

meta_detection_export_vector <- function(SDDprob,LDDprob=0,LDDrate=0,
                                         Survival=1,PropaguleProduction,
                                         ManageProb=0,MortalityProb=0,
                                         SpreadReduction=0,
                                         ExportSDDprob=NULL,ExportLDDprob=NULL) {
  e0<-meta_export_vector(SDDprob,LDDprob,LDDrate,Survival,PropaguleProduction,
                         0,0,0,ExportSDDprob,ExportLDDprob)
  e1<-meta_export_vector(SDDprob,LDDprob,LDDrate,Survival,PropaguleProduction,
                         ManageProb,MortalityProb,SpreadReduction,
                         ExportSDDprob,ExportLDDprob)
  c(e0,e1)
}

transition_export_vector <- function(Transition,Nstages,SDDprob,LDDprob=0,LDDrate=0,
                                     ManageProb=0,MortalityProb=0,SpreadReduction=0,
                                     DispersalDensityFactor=0,
                                     ExportSDDprob=NULL,ExportLDDprob=NULL) {
  SDD<-transition_zero_density_sdd(SDDprob,DispersalDensityFactor); n<-nrow(SDD); S<-Nstages
  LDD<-.ina_mat(LDDprob,n,"LDDprob"); A<-.transition_list(Transition,n,S)
  a<-.ina_recycle(ManageProb,n,"ManageProb"); g<-.ina_recycle(SpreadReduction,n,"SpreadReduction")
  if(length(MortalityProb)==1L) M<-matrix(MortalityProb,n,S)
  else if(length(MortalityProb)==S) M<-matrix(rep(MortalityProb,each=n),n,S)
  else if(is.matrix(MortalityProb)&&all(dim(MortalityProb)==c(n,S))) M<-MortalityProb
  else stop("MortalityProb must be scalar, length Nstages, or nodes x Nstages")
  os<-.ina_residual_export(SDD,ExportSDDprob,"SDDprob"); ol<-.ina_residual_export(LDD,ExportLDDprob,"LDDprob")
  out<-numeric(n*S); idx<-function(i,k)(i-1L)*S+k; r<-LDDrate
  for(i in seq_len(n)) for(k in 2:S) {
    f<-A[[i]][1,k]; if(f<=0) next
    out[idx(i,k)]<-f*((1-a[i])*((1-r)*os[i]+r*ol[i]) +
                         a[i]*(1-M[i,k])*((1-r)*os[i]+r*(1-g[i])*ol[i]))
  }
  out
}

transition_detection_export_vector <- function(Transition,Nstages,SDDprob,LDDprob=0,LDDrate=0,
                                               ManageProb=0,MortalityProb=0,SpreadReduction=0,
                                               DispersalDensityFactor=0,
                                               ExportSDDprob=NULL,ExportLDDprob=NULL) {
  e0<-transition_export_vector(Transition,Nstages,SDDprob,LDDprob,LDDrate,0,0,0,
                               DispersalDensityFactor,ExportSDDprob,ExportLDDprob)
  e1<-transition_export_vector(Transition,Nstages,SDDprob,LDDprob,LDDrate,ManageProb,MortalityProb,SpreadReduction,
                               DispersalDensityFactor,ExportSDDprob,ExportLDDprob)
  c(e0,e1)
}

mlu_export_vector <- function(SDDprob,LDDprob=0,LDDrate=0,Survival=1,
                              K,PropaguleProduction,ManageProb=0,MortalityProb=0,
                              SpreadReduction=0,current_code=FALSE,
                              ExportSDDprob=NULL,ExportLDDprob=NULL) {
  SDD<-as.matrix(SDDprob);n<-nrow(SDD);LDD<-.ina_mat(LDDprob,n,"LDDprob");K<-as.matrix(K);L<-ncol(K)
  s<-.ina_recycle(Survival,n,"Survival");prod<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction")
  A<-.mlu_matrix(ManageProb,n,L,"ManageProb");Mort<-.mlu_matrix(MortalityProb,n,L,"MortalityProb");Gr<-.mlu_matrix(SpreadReduction,n,L,"SpreadReduction")
  os<-.ina_residual_export(SDD,ExportSDDprob,"SDDprob");ol<-.ina_residual_export(LDD,ExportLDDprob,"LDDprob")
  idx<-function(i,l)(i-1L)*L+l;out<-numeric(n*L);r<-LDDrate
  for(i in seq_len(n)) {
    st<-.mlu_management_states(A[i,])
    for(l in seq_len(L)) {
      zsum<-0
      for(z in seq_len(nrow(st$M))) {
        Mv<-st$M[z,];pm<-st$p[z];surv<-s[i]*(1-Mort[i,l]*Mv[l])
        sf<-if(current_code)sum(1-Gr[i,]*Mv)else(1-Gr[i,l]*Mv[l])
        zsum<-zsum+pm*surv*prod[i]*((1-r)*os[i]+r*sf*ol[i])
      }
      out[idx(i,l)]<-zsum
    }
  }
  out
}

mlu_detection_export_vector <- function(SDDprob,LDDprob=0,LDDrate=0,Survival=1,
                                        K,PropaguleProduction,ManageProb=0,MortalityProb=0,
                                        SpreadReduction=0,current_code=FALSE,
                                        ExportSDDprob=NULL,ExportLDDprob=NULL) {
  e0<-mlu_export_vector(SDDprob,LDDprob,LDDrate,Survival,K,PropaguleProduction,
                        0,0,0,current_code,ExportSDDprob,ExportLDDprob)
  e1<-mlu_export_vector(SDDprob,LDDprob,LDDrate,Survival,K,PropaguleProduction,
                        ManageProb,MortalityProb,SpreadReduction,current_code,
                        ExportSDDprob,ExportLDDprob)
  c(e0,e1)
}

###############################################################################
### Unified user-facing analytical screening interface
###############################################################################

.ina_any_nonzero <- function(x) {
  if (is.null(x) || !length(x)) return(FALSE)
  any(is.finite(as.numeric(x)) & as.numeric(x) != 0, na.rm = TRUE)
}

.ina_normalize_ldd <- function(x) {
  if (length(x) == 1L && is.na(x)) 0 else x
}

.ina_slice_connection <- function(x, t, n, T, name) {
  x <- .ina_normalize_ldd(x)
  if (length(dim(x)) == 3L) {
    if (!identical(dim(x), c(n, n, T)))
      stop(name, " 3D array must have dimensions nodes x nodes x Ntimesteps")
    return(x[, , t, drop = TRUE])
  }
  if (length(x) == 1L) return(x)
  x <- as.matrix(x)
  if (!all(dim(x) == c(n, n))) stop(name, " must be nodes x nodes")
  x
}

.ina_slice_export <- function(x, t, n, T, name) {
  if (is.null(x)) return(NULL)
  dx <- dim(x)
  if (length(dx) == 3L) {
    if (dx[1] != n || dx[3] != T)
      stop(name, " 3D array must have dimensions nodes x outside-destinations x Ntimesteps")
    return(x[, , t, drop = TRUE])
  }
  x <- as.matrix(x)
  if (nrow(x) != n) stop(name, " rows must equal nodes")
  x
}

.ina_slice_node <- function(x, t, n, T, name) {
  if (is.null(x)) return(x)
  if (is.matrix(x) && all(dim(x) == c(n, T))) return(x[, t])
  if (length(x) == 1L || length(x) == n) return(x)
  stop(name, " must be scalar, length nodes, or nodes x Ntimesteps")
}

.ina_slice_stage <- function(x, t, n, S, T, name) {
  if (is.null(x)) return(x)
  dx <- dim(x)
  if (length(dx) == 3L) {
    if (!identical(dx, c(n, S, T)))
      stop(name, " 3D array must have dimensions nodes x stages x Ntimesteps")
    return(x[, , t, drop = TRUE])
  }
  if (length(x) == 1L || length(x) == S) return(x)
  if (is.matrix(x) && all(dim(x) == c(n, S))) return(x)
  stop(name, " must be scalar, length stages, nodes x stages, or nodes x stages x Ntimesteps")
}

.ina_slice_mlu <- function(x, t, n, L, T, name) {
  if (is.null(x)) return(x)
  dx <- dim(x)
  if (length(dx) == 3L) {
    if (!identical(dx, c(n, L, T)))
      stop(name, " 3D array must have dimensions nodes x land uses x Ntimesteps")
    return(x[, , t, drop = TRUE])
  }
  if (length(x) == 1L || length(x) == L) return(x)
  if (is.matrix(x) && all(dim(x) == c(n, L))) return(x)
  stop(name, " must be scalar, length land uses, nodes x land uses, or nodes x land uses x Ntimesteps")
}

.ina_slice_K_mlu <- function(K, t, n, L, T) {
  if (length(dim(K)) == 3L) {
    if (!identical(dim(K), c(n, L, T)))
      stop("K 3D array must have dimensions nodes x land uses x Ntimesteps")
    return(K[, , t, drop = TRUE])
  }
  K <- as.matrix(K)
  if (!all(dim(K) == c(n, L))) stop("K must be nodes x land uses")
  K
}

.ina_is_temporal <- function(x, n, T, kind = c("node", "connection", "stage", "mlu"), extra = NULL) {
  kind <- match.arg(kind)
  if (is.null(x)) return(FALSE)
  dx <- dim(x)
  if (kind == "connection") return(length(dx) == 3L)
  if (kind == "node") return(is.matrix(x) && all(dx == c(n, T)))
  if (kind == "stage") return(length(dx) == 3L)
  if (kind == "mlu") return(length(dx) == 3L)
  FALSE
}

.ina_apply_operators <- function(operators, initial) {
  x <- as.numeric(initial)
  out <- matrix(NA_real_, nrow = length(x), ncol = length(operators))
  for (t in seq_along(operators)) {
    x <- as.numeric(operators[[t]] %*% x)
    out[, t] <- x
  }
  out
}

.ina_cycle_growth <- function(operators) {
  if (length(operators) == 1L) {
    r <- ina_spectral_radius(operators[[1]])
    return(list(CycleMultiplier = r, EquivalentPerTimestepMultiplier = r,
                CycleOperator = operators[[1]]))
  }
  z <- ina_temporal_cycle(operators)
  z$EquivalentPerTimestepMultiplier <- z$GeometricPerTimestepMultiplier
  z
}

.ina_classify_growth <- function(r, tolerance = 1e-6) {
  if (!is.finite(r)) return(NA_character_)
  if (r < 1 - tolerance) return("decline when rare")
  if (r > 1 + tolerance) return("growth when rare")
  "approximately replacement"
}

.ina_overall_extinction <- function(q, initial) {
  q <- pmin(1, pmax(0, as.numeric(q)))
  x <- as.numeric(initial)
  if (length(q) != length(x)) stop("Extinction vector and initial state differ in length")
  if (any(x < 0 | !is.finite(x))) stop("InitialState must be finite and non-negative")
  # Independent-lineage branching approximation. Fractional x is allowed as a
  # continuous approximation to an expected initial state, with a diagnostic.
  if (any(q == 0 & x > 0)) return(0)
  exp(sum(x * log(pmax(q, .Machine$double.xmin))))
}

.ina_overall_no_escape <- function(g, initial) {
  g <- pmin(1, pmax(0, as.numeric(g)))
  x <- as.numeric(initial)
  if (length(g) != length(x)) stop("No-escape vector and initial state differ in length")
  if (any(g == 0 & x > 0)) return(0)
  exp(sum(x * log(pmax(g, .Machine$double.xmin))))
}

# Generic multitype Poisson branching fallback. This is deliberately distinct
# from the older first-moment Poisson approximation: the recursion propagates
# the full lineage no-event/extinction probability through descendants.
ina_poisson_branching_extinction <- function(G, generations = 100, tolerance = 1e-12) {
  G <- as.matrix(G); nt <- nrow(G)
  if (ncol(G) != nt) stop("G must be square")
  q <- rep(0, nt); hist <- matrix(NA_real_, nt, generations)
  for (tt in seq_len(generations)) {
    old <- q
    q <- exp(colSums(sweep(G, 1, q - 1, `*`)))
    q <- pmin(1, pmax(0, q)); hist[, tt] <- q
    if (max(abs(q - old)) < tolerance) {
      if (tt < generations) hist[, (tt + 1):generations] <- q
      break
    }
  }
  list(extinction = q, history = hist,
       note = "Mean-matched multitype Poisson branching fallback")
}

ina_poisson_branching_extinction_horizon <- function(operators) {
  if (!length(operators)) stop("operators must be non-empty")
  q <- rep(0, nrow(operators[[1]]))
  for (tt in rev(seq_along(operators)))
    q <- exp(colSums(sweep(operators[[tt]], 1, q - 1, `*`)))
  pmin(1, pmax(0, q))
}

ina_poisson_branching_no_escape <- function(G, successful_export_mean,
                                             generations = 10) {
  G <- as.matrix(G); e <- as.numeric(successful_export_mean); nt <- nrow(G)
  if (ncol(G) != nt || length(e) != nt) stop("G and successful_export_mean are incompatible")
  g <- rep(1, nt); hist <- matrix(NA_real_, nt, generations)
  for (tt in seq_len(generations)) {
    g <- exp(-e + colSums(sweep(G, 1, g - 1, `*`)))
    g <- pmin(1, pmax(0, g)); hist[, tt] <- g
  }
  list(no_escape = g, escape = 1 - g, history_escape = 1 - hist,
       note = "Mean-matched multitype branching no-escape recursion; successful external events are Poisson conditional on source type")
}

ina_poisson_branching_no_escape_horizon <- function(operators, successful_export_means) {
  if (length(operators) != length(successful_export_means))
    stop("operators and successful_export_means must have the same length")
  g <- rep(1, nrow(operators[[1]]))
  for (tt in rev(seq_along(operators))) {
    e <- as.numeric(successful_export_means[[tt]])
    g <- exp(-e + colSums(sweep(operators[[tt]], 1, g - 1, `*`)))
    g <- pmin(1, pmax(0, g))
  }
  g
}

ina_escape_first_moment_temporal <- function(operators, state0, export_vectors) {
  if (length(operators) != length(export_vectors))
    stop("operators and export_vectors must have the same length")
  x <- as.numeric(state0); by_step <- numeric(length(operators))
  for (tt in seq_along(operators)) {
    by_step[tt] <- sum(x * as.numeric(export_vectors[[tt]]))
    x <- as.numeric(operators[[tt]] %*% x)
  }
  cum <- cumsum(by_step)
  list(ExpectedSuccessfulEscapesByStep = by_step,
       ExpectedSuccessfulEscapesCumulative = cum,
       EscapeProbabilityByStep = 1 - exp(-cum),
       EscapeProbability = if (length(cum)) tail(1 - exp(-cum), 1) else 0)
}

.inapest_export_vector <- function(ExportProb, Survival = 1, ManageProb = 0,
                                    EradicationProb = 0, SpreadReduction = 0,
                                    dynamic_information = FALSE) {
  X <- as.matrix(ExportProb); n <- nrow(X)
  s <- .inapest_recycle(Survival, n, "Survival")
  a <- .inapest_recycle(ManageProb, n, "ManageProb")
  e <- .inapest_recycle(EradicationProb, n, "EradicationProb")
  r <- .inapest_recycle(SpreadReduction, n, "SpreadReduction")
  e0 <- s * rowSums(X)
  e1 <- s * ((1 - a) + a * (1 - e) * (1 - r)) * rowSums(X)
  if (dynamic_information) c(e0, e1) else e1
}

# Binary-INApest branching no-escape recursion with dynamic information.
inapest_detection_escape_branching <- function(SDDprob, ExportProb, LDDprob = 0,
                                                EnvEstabProb = 1, Survival = 1,
                                                DetectionProb = 0, ManageProb = 0,
                                                EradicationProb = 0,
                                                SpreadReduction = 0,
                                                SEAM = NULL,
                                                InfoRetentionProb = 1,
                                                timesteps = 10) {
  P <- inapest_edge_prob(SDDprob, LDDprob, EnvEstabProb)
  n <- nrow(P); X <- as.matrix(ExportProb)
  if (nrow(X) != n) stop("ExportProb rows must equal nodes")
  s <- .inapest_recycle(Survival, n, "Survival")
  d <- .inapest_recycle(DetectionProb, n, "DetectionProb")
  a <- .inapest_recycle(ManageProb, n, "ManageProb")
  e <- .inapest_recycle(EradicationProb, n, "EradicationProb")
  r <- .inapest_recycle(SpreadReduction, n, "SpreadReduction")
  ir <- .inapest_recycle(InfoRetentionProb, n, "InfoRetentionProb")
  C <- if (is.null(SEAM) || length(SEAM) == 1L) matrix(0, n, n) else as.matrix(SEAM)
  diag(C) <- 0
  gU <- rep(1, n); gH <- rep(1, n)
  hU <- matrix(NA_real_, n, timesteps); hH <- hU
  for (tt in seq_len(timesteps)) {
    nU <- numeric(n); nH <- numeric(n)
    for (i in seq_len(n)) {
      off <- seq_len(n) != i
      childU <- (1 - d) * gU + d * gH
      internal0 <- prod(1 - P[i, off] + P[i, off] * childU[off])
      outside0 <- prod(1 - X[i, ])
      selfU <- (1 - d[i]) * gU[i] + d[i] * gH[i]
      nU[i] <- (1 - s[i]) + s[i] * selfU * internal0 * outside0

      hself <- ir[i] + (1 - ir[i]) * d[i]
      selfH <- (1 - hself) * gU[i] + hself * gH[i]
      hchild <- 1 - (1 - d) * (1 - C[i, ])
      childH <- (1 - hchild) * gU + hchild * gH
      internal_unmanaged <- prod(1 - P[i, off] + P[i, off] * childH[off])
      outside_unmanaged <- prod(1 - X[i, ])
      Pm <- P[i, ] * (1 - r[i]); Xm <- X[i, ] * (1 - r[i])
      internal_managed <- prod(1 - Pm[off] + Pm[off] * childH[off])
      outside_managed <- prod(1 - Xm)
      f0 <- (1 - s[i]) + s[i] * selfH * internal_unmanaged * outside_unmanaged
      sm <- s[i] * (1 - e[i])
      f1 <- (1 - sm) + sm * selfH * internal_managed * outside_managed
      nH[i] <- (1 - a[i]) * f0 + a[i] * f1
    }
    gU <- pmin(1, pmax(0, nU)); gH <- pmin(1, pmax(0, nH))
    hU[, tt] <- gU; hH[, tt] <- gH
  }
  list(no_escape_U = gU, no_escape_H = gH,
       escape_U = 1 - gU, escape_H = 1 - gH,
       history_escape_U = 1 - hU, history_escape_H = 1 - hH,
       note = "Binary edge-based branching recursion with direct-child SEAM approximation")
}

.meta_export_fraction <- function(InternalProb, ExternalProb, allow_residual, name) {
  I <- as.matrix(InternalProb); n <- nrow(I)
  if (!is.null(ExternalProb)) {
    E <- as.matrix(ExternalProb)
    if (nrow(E) != n) stop(name, " external matrix rows must match internal nodes")
    return(rowSums(E))
  }
  if (!allow_residual) return(NULL)
  pmax(0, 1 - rowSums(I))
}

.meta_no_export_probability <- function(PropaguleProduction,
                                        sdd_internal_mass, os,
                                        ldd_internal_mass, ol, LDDrate,
                                        SpreadReduction, OutsideEstablishmentProb,
                                        managed = FALSE) {
  n <- length(PropaguleProduction); out <- numeric(n)
  for (i in seq_len(n)) {
    supp <- .meta_poisson_support(PropaguleProduction[i]); k <- supp$k; pk <- supp$p

    # Reconstruct each cropped multinomial over internal + external destinations.
    # Integer disperser counts are calculated from the combined kernel mass,
    # then the no-escape probability follows from the external share of that
    # multinomial.  Separately flooring k * outside_mass would incorrectly make
    # small but real export pathways disappear at low propagule production.
    sdd_total <- sdd_internal_mass[i] + os[i]
    sdd_size <- floor(k * (1 - LDDrate) * sdd_total)
    sdd_out_share <- if (sdd_total > 0) os[i] / sdd_total else 0
    no_sdd <- (1 - sdd_out_share * OutsideEstablishmentProb[i])^sdd_size

    ldd_total <- ldd_internal_mass[i] + ol[i]
    ldd_size <- floor(k * LDDrate * (1 - SpreadReduction[i] * managed) * ldd_total)
    ldd_out_share <- if (ldd_total > 0) ol[i] / ldd_total else 0
    no_ldd <- (1 - ldd_out_share * OutsideEstablishmentProb[i])^ldd_size

    out[i] <- sum(pk * no_sdd * no_ldd)
  }
  out
}

# INApestMeta branching no-escape approximation using simulator-faithful
# one-parent recruit means plus an integrated one-parent no-export probability.
meta_single_parent_escape_branching <- function(SDDprob, LDDprob = 0, LDDrate = 0,
                                                EnvEstabProb = 1, Survival = 1, K,
                                                PropaguleProduction,
                                                PropaguleEstablishment,
                                                ManageProb = 0, MortalityProb = 0,
                                                SpreadReduction = 0,
                                                ExportSDDprob = NULL,
                                                ExportLDDprob = NULL,
                                                OutsideEstablishmentProb = 1,
                                                AssumeResidualExport = FALSE,
                                                timesteps = 10) {
  SDD <- as.matrix(SDDprob); n <- nrow(SDD); LDD <- .ina_mat(LDDprob, n, "LDDprob")
  os <- .meta_export_fraction(SDD, ExportSDDprob, AssumeResidualExport, "SDD")
  ol <- .meta_export_fraction(LDD, ExportLDDprob, AssumeResidualExport, "LDD")
  if (is.null(os) && is.null(ol)) stop("Explicit export matrices are required unless AssumeResidualExport=TRUE")
  if (is.null(os)) os <- rep(0, n); if (is.null(ol)) ol <- rep(0, n)
  s <- .ina_recycle(Survival, n, "Survival")
  a <- .ina_recycle(ManageProb, n, "ManageProb")
  m <- .ina_recycle(MortalityProb, n, "MortalityProb")
  g <- .ina_recycle(SpreadReduction, n, "SpreadReduction")
  prod <- .ina_recycle(PropaguleProduction, n, "PropaguleProduction")
  pout <- .ina_recycle(OutsideEstablishmentProb, n, "OutsideEstablishmentProb")
  mus <- meta_single_parent_recruit_means(SDD, LDD, LDDrate, EnvEstabProb, K,
                                           PropaguleProduction, PropaguleEstablishment,
                                           SpreadReduction)
  no0 <- .meta_no_export_probability(prod, rowSums(SDD), os, rowSums(LDD), ol,
                                         LDDrate, g, pout, FALSE)
  no1 <- .meta_no_export_probability(prod, rowSums(SDD), os, rowSums(LDD), ol,
                                         LDDrate, g, pout, TRUE)
  h <- rep(1, n); hist <- matrix(NA_real_, n, timesteps)
  for (tt in seq_len(timesteps)) {
    hn <- numeric(n)
    for (i in seq_len(n)) {
      q0 <- s[i]; q1 <- s[i] * (1 - m[i])
      R0 <- exp(sum(mus[[1]][i, ] * (h - 1)))
      R1 <- exp(sum(mus[[2]][i, ] * (h - 1)))
      f0 <- (1 - q0) + q0 * h[i] * R0 * no0[i]
      f1 <- (1 - q1) + q1 * h[i] * R1 * no1[i]
      hn[i] <- (1 - a[i]) * f0 + a[i] * f1
    }
    h <- pmin(1, pmax(0, hn)); hist[, tt] <- h
  }
  list(no_escape = h, escape = 1 - h, history_escape = 1 - hist,
       note = "Branching PGF using simulator-faithful one-parent recruitment and integrated integer export counts; recruit/export correlation remains approximated")
}

transition_escape_branching <- function(Transition, Nstages, SDDprob,
                                        LDDprob = 0, LDDrate = 0,
                                        EnvEstabProb = 1,
                                        PropaguleEstablishment = 1,
                                        ManageProb = 0, MortalityProb = 0,
                                        SpreadReduction = 0,
                                        DispersalDensityFactor = 0,
                                        K = 1, SeedbankK = 1,
                                        ExportSDDprob = NULL,
                                        ExportLDDprob = NULL,
                                        OutsideEstablishmentProb = 1,
                                        AssumeResidualExport = TRUE,
                                        timesteps = 10) {
  z <- transition_components(Transition, Nstages, SDDprob, LDDprob, LDDrate,
                             EnvEstabProb, PropaguleEstablishment, ManageProb,
                             MortalityProb, SpreadReduction, DispersalDensityFactor,
                             K, SeedbankK)
  n <- z$n; S <- z$S; nt <- n * S; idx <- function(i, s) (i - 1L) * S + s
  os <- .meta_export_fraction(z$SDD, ExportSDDprob, AssumeResidualExport, "SDD")
  ol <- .meta_export_fraction(z$LDD, ExportLDDprob, AssumeResidualExport, "LDD")
  if (is.null(os)) os <- rep(0, n); if (is.null(ol)) ol <- rep(0, n)
  pout <- .ina_recycle(OutsideEstablishmentProb, n, "OutsideEstablishmentProb")
  h <- rep(1, nt); hist <- matrix(NA_real_, nt, timesteps)
  for (tt in seq_len(timesteps)) {
    hn <- numeric(nt)
    for (i in seq_len(n)) {
      Ai <- z$A[[i]]; ai <- z$adoption[i]; gi <- z$spread_reduction[i]
      for (k in seq_len(S)) {
        src <- idx(i, k); fec <- if (k >= 2L) Ai[1, k] else 0
        calcM <- function(M) {
          surv <- 1 - z$mortality[i, k] * M
          if (k < S) {
            local <- (1 - Ai[k, k] - Ai[k + 1L, k]) +
              Ai[k, k] * h[idx(i, k)] + Ai[k + 1L, k] * h[idx(i, k + 1L)]
          } else {
            local <- (1 - Ai[S, S]) + Ai[S, S] * h[idx(i, S)]
          }
          muq <- 0
          if (fec > 0) for (j in seq_len(n)) {
            nat <- (1 - z$LDDrate) * z$sdd_enabled[i] * z$SDD[i, j]
            hum <- z$LDDrate * (1 - gi * M) * z$LDD[i, j]
            mu <- fec * (nat + hum) * z$recruit_success[j]
            muq <- muq + mu * (h[idx(j, 1L)] - 1)
          }
          outprob <- (1 - z$LDDrate) * z$sdd_enabled[i] * os[i] +
                     z$LDDrate * (1 - gi * M) * ol[i]
          noexp <- exp(-fec * outprob * pout[i])
          (1 - surv) + surv * local * exp(muq) * noexp
        }
        hn[src] <- (1 - ai) * calcM(0) + ai * calcM(1)
      }
    }
    h <- pmin(1, pmax(0, hn)); hist[, tt] <- h
  }
  list(no_escape = h, escape = 1 - h, history_escape = 1 - hist,
       note = "Stage x node branching PGF; Poisson fecundity is thinned into internal recruitment and successful external export")
}

.ina_initial_detection_binary <- function(x, d) pmin(1, pmax(0, d))

.ina_initial_detection_meta <- function(x, d) {
  x <- pmax(0, as.numeric(x)); d <- pmin(1, pmax(0, d))
  1 - (1 - d)^x
}

.ina_initial_detection_transition <- function(X, D) {
  1 - apply((1 - D)^X, 1, prod)
}

.ina_initial_detection_mlu <- function(X, D) {
  LU <- 1 - (1 - D)^X
  1 - apply(1 - LU, 1, prod)
}

.ina_split_information <- function(base_state, node_info, n, block) {
  info_by_type <- rep(node_info, each = block)
  c(base_state * (1 - info_by_type), base_state * info_by_type)
}

.ina_programmed_layout <- function(base_types, max_age) {
  max_age <- max(1L, as.integer(max_age))
  blocks <- max_age + 3L
  list(
    base_types = base_types,
    max_age = max_age,
    blocks = blocks,
    U = seq_len(base_types),
    X = base_types + seq_len(base_types),
    H = lapply(seq_len(max_age), function(a) (1L + a) * base_types + seq_len(base_types)),
    Overflow = (max_age + 2L) * base_types + seq_len(base_types),
    size = blocks * base_types
  )
}

.ina_programmed_global_max_age <- function(InfoPersistenceSteps) {
  x <- as.numeric(InfoPersistenceSteps)
  finite <- x[!is.na(x)]
  if (!length(finite)) return(NULL)
  if (any(!is.finite(finite) | finite < 0 | finite != floor(finite)))
    stop("InfoPersistenceSteps values must be non-negative whole numbers or NA")
  max(1L, as.integer(max(finite)))
}

.ina_programmed_initial_state <- function(base_state, InitialInfo,
                                          InitialDetectionProb,
                                          type_node, layout,
                                          binary_known_presence = FALSE) {
  B <- length(base_state)
  n <- max(type_node)
  p0 <- .ina_recycle(InitialInfo, n, "InitialInfo")
  pd <- .ina_recycle(InitialDetectionProb, n, "InitialDetectionProb")
  p0 <- pmin(1, pmax(0, p0)); pd <- pmin(1, pmax(0, pd))
  out <- numeric(layout$size)
  for (b in seq_len(B)) {
    i <- type_node[b]
    if (binary_known_presence) {
      ph <- p0[i] + (1 - p0[i]) * pd[i]
      pu <- 1 - ph
      px <- 0
    } else {
      ph <- pd[i]
      px <- p0[i] * (1 - pd[i])
      pu <- (1 - p0[i]) * (1 - pd[i])
    }
    out[layout$U[b]] <- base_state[b] * pu
    out[layout$X[b]] <- base_state[b] * px
    out[layout$H[[1L]][b]] <- base_state[b] * ph
  }
  out
}

.ina_programmed_expand_vector <- function(v, layout) {
  B <- layout$base_types
  if (length(v) != 2L * B)
    stop("Programmed-information vector expansion requires a 2B uninformed/informed vector")
  out <- numeric(layout$size)
  out[layout$U] <- v[seq_len(B)]
  informed <- v[B + seq_len(B)]
  out[layout$X] <- informed
  for (h in layout$H) out[h] <- informed
  out[layout$Overflow] <- informed
  out
}

# Expand the ordinary uninformed/informed low-density operator into explicit
# time-since-last-local-evidence states.  The state at the start of a timestep
# is U (uninformed), X (informed without a currently valid local-evidence
# clock), H1...Hmax (1...Hmax timesteps since local evidence), or Overflow
# (>Hmax timesteps since local evidence).  Management is applied to every
# informed state before the programmed stop is evaluated, matching the
# simulation ordering.  A local detection resets the next-timestep state to H1.
# For binary INApest, an informed extant infestation is itself known local
# presence, so surviving same-node lineages reset to H1 every timestep.
.ina_programmed_information_operator <- function(G_unmanaged, G_managed,
                                                  type_node,
                                                  DetectionProbByType,
                                                  SEAM = NULL,
                                                  InfoRetentionProb = 1,
                                                  InfoPersistenceSteps = NA,
                                                  layout,
                                                  binary_known_presence = FALSE) {
  G0 <- as.matrix(G_unmanaged); GH <- as.matrix(G_managed)
  if (!all(dim(G0) == dim(GH)) || nrow(G0) != ncol(G0))
    stop("G_unmanaged and G_managed must be square matrices of equal dimension")
  B <- nrow(G0)
  if (length(type_node) != B) stop("type_node must identify the node for every base type")
  n <- max(type_node)
  d <- as.numeric(DetectionProbByType)
  if (length(d) == 1L) d <- rep(d, B)
  if (length(d) != B) stop("DetectionProbByType must be scalar or length equal to base types")
  d <- pmin(1, pmax(0, d))
  ir <- .ina_recycle(InfoRetentionProb, n, "InfoRetentionProb")
  K <- .ina_recycle(InfoPersistenceSteps, n, "InfoPersistenceSteps")
  finite <- K[!is.na(K)]
  if (length(finite) && any(!is.finite(finite) | finite < 0 | finite != floor(finite)))
    stop("InfoPersistenceSteps values must be non-negative whole numbers or NA")
  if (is.null(SEAM) || length(SEAM) == 1L) C <- matrix(0, n, n) else C <- as.matrix(SEAM)
  if (!all(dim(C) == c(n, n))) stop("SEAM must be nodes x nodes")
  diag(C) <- 0

  G <- matrix(0, layout$size, layout$size)
  Hblock <- function(age, b) layout$H[[age]][b]

  add_new_destination <- function(dst, src_col, weight, informed_source, src_node) {
    if (weight == 0) return(invisible(NULL))
    jnode <- type_node[dst]
    pdet <- d[dst]
    if (informed_source && jnode != src_node) {
      pseam <- C[src_node, jnode]
      G[Hblock(1L, dst), src_col] <<- G[Hblock(1L, dst), src_col] + weight * pdet
      G[layout$X[dst], src_col] <<- G[layout$X[dst], src_col] + weight * (1 - pdet) * pseam
      G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + weight * (1 - pdet) * (1 - pseam)
    } else {
      G[Hblock(1L, dst), src_col] <<- G[Hblock(1L, dst), src_col] + weight * pdet
      G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + weight * (1 - pdet)
    }
    invisible(NULL)
  }

  add_same_node_informed <- function(dst, src_col, weight, source_kind, source_age, node) {
    if (weight == 0) return(invisible(NULL))
    if (binary_known_presence) {
      G[Hblock(1L, dst), src_col] <<- G[Hblock(1L, dst), src_col] + weight
      return(invisible(NULL))
    }
    pdet <- d[dst]
    G[Hblock(1L, dst), src_col] <<- G[Hblock(1L, dst), src_col] + weight * pdet
    rem <- weight * (1 - pdet)
    if (rem == 0) return(invisible(NULL))

    if (identical(source_kind, "X")) {
      if (is.na(K[node])) {
        G[layout$X[dst], src_col] <<- G[layout$X[dst], src_col] + rem * ir[node]
        G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + rem * (1 - ir[node])
      } else {
        # No valid local-evidence clock: programmed stopping occurs after this
        # timestep's management/spread unless local evidence is generated.
        G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + rem
      }
      return(invisible(NULL))
    }

    if (identical(source_kind, "Overflow")) {
      if (is.na(K[node])) {
        G[layout$Overflow[dst], src_col] <<- G[layout$Overflow[dst], src_col] + rem * ir[node]
        G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + rem * (1 - ir[node])
      } else {
        G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + rem
      }
      return(invisible(NULL))
    }

    age <- source_age
    if (is.na(K[node])) {
      next_idx <- if (age < layout$max_age) Hblock(age + 1L, dst) else layout$Overflow[dst]
      G[next_idx, src_col] <<- G[next_idx, src_col] + rem * ir[node]
      G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + rem * (1 - ir[node])
    } else if (age >= K[node]) {
      # The node was informed for management in this timestep, then the
      # programmed window expires before the next timestep.
      G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + rem
    } else {
      next_idx <- if (age < layout$max_age) Hblock(age + 1L, dst) else layout$Overflow[dst]
      G[next_idx, src_col] <<- G[next_idx, src_col] + rem
    }
    invisible(NULL)
  }

  for (src in seq_len(B)) {
    inode <- type_node[src]
    # Uninformed source: no management. Detection creates genuine local evidence.
    scol <- layout$U[src]
    nz <- which(G0[, src] != 0)
    for (dst in nz) add_new_destination(dst, scol, G0[dst, src], FALSE, inode)

    # X source: informed for this timestep but without a valid local-evidence clock.
    scol <- layout$X[src]
    nz <- which(GH[, src] != 0)
    for (dst in nz) {
      if (type_node[dst] == inode) add_same_node_informed(dst, scol, GH[dst, src], "X", NA_integer_, inode)
      else add_new_destination(dst, scol, GH[dst, src], TRUE, inode)
    }

    # Explicit local-evidence ages.
    for (age in seq_len(layout$max_age)) {
      scol <- layout$H[[age]][src]
      for (dst in nz) {
        if (type_node[dst] == inode) add_same_node_informed(dst, scol, GH[dst, src], "H", age, inode)
        else add_new_destination(dst, scol, GH[dst, src], TRUE, inode)
      }
    }

    # Age older than every finite persistence value supplied anywhere.
    scol <- layout$Overflow[src]
    for (dst in nz) {
      if (type_node[dst] == inode) add_same_node_informed(dst, scol, GH[dst, src], "Overflow", NA_integer_, inode)
      else add_new_destination(dst, scol, GH[dst, src], TRUE, inode)
    }
  }
  attr(G, "note") <- paste0(
    "Finite-state time-since-local-evidence operator. Programmed stopping is exact for the represented information clock; ",
    "same-node shared evidence from other individuals and information-only preconditioning remain low-density approximations."
  )
  attr(G, "programmed_information_layout") <- layout
  G
}

.ina_info_mode <- function(mode, ManageProb, InitialInfo, DetectionProb, SEAM,
                           InfoRetentionProb, InfoPersistenceSteps) {
  mode <- match.arg(mode, c("auto", "dynamic", "all_informed", "none"))
  if (mode != "auto") return(mode)
  if (!.ina_any_nonzero(ManageProb)) return("none")
  if (all(as.numeric(InitialInfo) >= 1, na.rm = TRUE) &&
      !.ina_any_nonzero(DetectionProb) && !.ina_any_nonzero(SEAM) &&
      all(as.numeric(InfoRetentionProb) >= 1, na.rm = TRUE) &&
      all(is.na(InfoPersistenceSteps))) return("all_informed")
  "dynamic"
}

.ina_state_total <- function(state_matrix, dynamic) {
  if (!dynamic) return(colSums(state_matrix))
  nt <- nrow(state_matrix) / 2L
  colSums(state_matrix[seq_len(nt), , drop = FALSE] +
          state_matrix[nt + seq_len(nt), , drop = FALSE])
}

.ina_export_available <- function(Model, ExportProb, ExportSDDprob, ExportLDDprob,
                                  AssumeResidualExport) {
  if (Model == "INApest") return(!is.null(ExportProb))
  !is.null(ExportSDDprob) || !is.null(ExportLDDprob) || isTRUE(AssumeResidualExport)
}

#' Analytical screening companion for INApest model families
#'
#' Returns low-density growth, an initial-condition-specific expected trajectory,
#' branching extinction, and (when export is identifiable) branching escape risk.
#' Finite InfoPersistenceSteps values are represented with explicit time-since-local-evidence states;
#' InfoRetentionProb continues to apply where InfoPersistenceSteps is NA.
#' This is a screening approximation, not a replacement for stochastic simulation.
INApestAnalytical <- function(
    Model = c("INApest", "INApestMeta", "INApestMetaTransitionMatrix",
              "INApestMetaMultipleLandUse"),
    Ntimesteps = 10,
    InitialState,
    InitialInfo = 0,
    InformationMode = c("auto", "dynamic", "all_informed", "none"),
    ApplyInitialDetection = TRUE,
    SDDprob,
    LDDprob = 0,
    LDDrate = 0,
    EnvEstabProb = 1,
    Survival = 1,
    K = NULL,
    PropaguleProduction = NULL,
    PropaguleEstablishment = 1,
    Transition = NULL,
    Nstages = NULL,
    SeedbankK = NULL,
    DetectionProb = 0,
    ManageProb = 0,
    EradicationProb = 0,
    MortalityProb = 0,
    SpreadReduction = 0,
    SEAM = NULL,
    InfoRetentionProb = 1,
    InfoPersistenceSteps = NA,
    DispersalDensityFactor = 0,
    ExportProb = NULL,
    ExportSDDprob = NULL,
    ExportLDDprob = NULL,
    OutsideEstablishmentProb = 1,
    AssumeResidualExport = FALSE,
    ExtinctionGenerations = 100,
    ReturnOperators = FALSE) {

  Model <- match.arg(Model)
  InformationMode <- match.arg(InformationMode)
  if (length(Ntimesteps) != 1L || Ntimesteps < 1 || Ntimesteps != as.integer(Ntimesteps))
    stop("Ntimesteps must be a positive integer")
  Ntimesteps <- as.integer(Ntimesteps)
  LDDprob <- .ina_normalize_ldd(LDDprob)
  SDD0 <- if (length(dim(SDDprob)) == 3L) SDDprob[, , 1] else as.matrix(SDDprob)
  n <- nrow(SDD0)
  if (ncol(SDD0) != n) stop("SDDprob must be square")
  diagnostics <- character(0)
  # Validate the programmed-persistence parameterisation up front.  A finite
  # value activates explicit time-since-local-evidence state classes; NA leaves
  # that node/timestep on the stochastic InfoRetentionProb pathway.
  .ina_slice_node(InfoPersistenceSteps, 1, n, Ntimesteps, "InfoPersistenceSteps")
  persistence_max_age <- .ina_programmed_global_max_age(InfoPersistenceSteps)
  persistence_requested <- !is.null(persistence_max_age)
  if (persistence_requested) {
    diagnostics <- c(diagnostics,
      "InfoPersistenceSteps is represented with explicit time-since-last-local-evidence states. Programmed stopping takes priority over InfoRetentionProb wherever the persistence value is finite.")
    if (any(as.numeric(InfoRetentionProb) < 1, na.rm = TRUE))
      diagnostics <- c(diagnostics,
        "Both InfoPersistenceSteps and InfoRetentionProb are supplied: programmed stopping has priority for finite InfoPersistenceSteps; stochastic retention applies only where InfoPersistenceSteps is NA.")
  }
  if (any(as.numeric(InitialState) %% 1 != 0, na.rm = TRUE)) {
    diagnostics <- c(diagnostics,
      "InitialState contains fractional values. Branching event probabilities treat these as continuous lineage weights and are therefore approximate.")
  }

  mode <- .ina_info_mode(InformationMode, ManageProb, InitialInfo,
                         DetectionProb, SEAM, InfoRetentionProb,
                         InfoPersistenceSteps)
  dynamic <- identical(mode, "dynamic")
  programmed_dynamic <- dynamic && persistence_requested
  information_state_method <- if (programmed_dynamic)
    "explicit time-since-local-evidence states" else if (dynamic)
    "uninformed/informed memoryless states" else mode
  programmed_layout <- NULL

  operators <- vector("list", Ntimesteps)
  export_vectors <- vector("list", Ntimesteps)
  intrinsic <- NULL
  type_count <- NULL
  base_initial <- NULL

  if (Model == "INApest") {
    base_initial <- as.numeric(InitialState)
    if (length(base_initial) != n) stop("InitialState must have length nodes for INApest")
    d1 <- .ina_slice_node(DetectionProb, 1, n, Ntimesteps, "DetectionProb")
    pdet0 <- if (ApplyInitialDetection && dynamic)
      .ina_initial_detection_binary(base_initial, .inapest_recycle(d1, n, "DetectionProb")) else rep(0, n)
    if (programmed_dynamic) {
      programmed_layout <- .ina_programmed_layout(n, persistence_max_age)
      state0 <- .ina_programmed_initial_state(base_initial, InitialInfo, pdet0,
                                               seq_len(n), programmed_layout,
                                               binary_known_presence = TRUE)
    } else {
      info0 <- .inapest_recycle(InitialInfo, n, "InitialInfo")
      if (ApplyInitialDetection && dynamic) info0 <- info0 + (1 - info0) * pdet0
      info0 <- pmin(1, pmax(0, info0))
      state0 <- if (dynamic) c(base_initial * (1 - info0), base_initial * info0) else base_initial
    }
    type_count <- length(state0)

    for (tt in seq_len(Ntimesteps)) {
      SDDt <- .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob")
      LDDt <- .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob")
      Et <- .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb")
      St <- .ina_slice_node(Survival, tt, n, Ntimesteps, "Survival")
      Dt <- .ina_slice_node(DetectionProb, tt, n, Ntimesteps, "DetectionProb")
      At <- .ina_slice_node(ManageProb, tt, n, Ntimesteps, "ManageProb")
      Kt <- .ina_slice_node(EradicationProb, tt, n, Ntimesteps, "EradicationProb")
      Rt <- .ina_slice_node(SpreadReduction, tt, n, Ntimesteps, "SpreadReduction")
      IRt <- .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb")
      IPt <- .ina_slice_node(InfoPersistenceSteps, tt, n, Ntimesteps, "InfoPersistenceSteps")
      if (programmed_dynamic) {
        G0 <- inapest_exogenous_operator(SDDt, LDDt, Et, St, 0, 0, 0)
        GH <- inapest_exogenous_operator(SDDt, LDDt, Et, St, At, Kt, Rt)
        operators[[tt]] <- .ina_programmed_information_operator(
          G0, GH, seq_len(n), Dt, SEAM, IRt, IPt, programmed_layout,
          binary_known_presence = TRUE)
      } else if (mode == "dynamic")
        operators[[tt]] <- inapest_detection_operator(SDDt, LDDt, Et, St, Dt, At, Kt, Rt, SEAM, IRt)
      else if (mode == "all_informed")
        operators[[tt]] <- inapest_exogenous_operator(SDDt, LDDt, Et, St, At, Kt, Rt)
      else
        operators[[tt]] <- inapest_exogenous_operator(SDDt, LDDt, Et, St, 0, 0, 0)

      if (!is.null(ExportProb)) {
        X <- .ina_slice_export(ExportProb, tt, n, Ntimesteps, "ExportProb")
        pout <- .inapest_recycle(OutsideEstablishmentProb, n, "OutsideEstablishmentProb")
        Xeff <- sweep(X, 1, pout, `*`)
        if (mode == "dynamic") {
          e <- .inapest_export_vector(Xeff, St, At, Kt, Rt, TRUE)
          export_vectors[[tt]] <- if (programmed_dynamic) .ina_programmed_expand_vector(e, programmed_layout) else e
        } else if (mode == "all_informed") export_vectors[[tt]] <- .inapest_export_vector(Xeff, St, At, Kt, Rt, FALSE)
        else export_vectors[[tt]] <- .inapest_export_vector(Xeff, St, 0, 0, 0, FALSE)
      }
    }
  }

  if (Model == "INApestMeta") {
    if (is.null(K) || is.null(PropaguleProduction)) stop("K and PropaguleProduction are required for INApestMeta")
    base_initial <- as.numeric(InitialState)
    if (length(base_initial) != n) stop("InitialState must have length nodes for INApestMeta")
    d1 <- .ina_slice_node(DetectionProb, 1, n, Ntimesteps, "DetectionProb")
    pdet0 <- if (ApplyInitialDetection && dynamic)
      .ina_initial_detection_meta(base_initial, .ina_recycle(d1, n, "DetectionProb")) else rep(0, n)
    if (programmed_dynamic) {
      programmed_layout <- .ina_programmed_layout(n, persistence_max_age)
      state0 <- .ina_programmed_initial_state(base_initial, InitialInfo, pdet0,
                                               seq_len(n), programmed_layout,
                                               binary_known_presence = FALSE)
    } else {
      info0 <- .ina_recycle(InitialInfo, n, "InitialInfo")
      if (ApplyInitialDetection && dynamic) info0 <- info0 + (1 - info0) * pdet0
      info0 <- pmin(1, pmax(0, info0))
      state0 <- if (dynamic) c(base_initial * (1 - info0), base_initial * info0) else base_initial
    }
    type_count <- length(state0)

    for (tt in seq_len(Ntimesteps)) {
      SDDt <- .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob")
      LDDt <- .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob")
      Et <- .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb")
      St <- .ina_slice_node(Survival, tt, n, Ntimesteps, "Survival")
      Kt <- .ina_slice_node(K, tt, n, Ntimesteps, "K")
      Pt <- .ina_slice_node(PropaguleProduction, tt, n, Ntimesteps, "PropaguleProduction")
      PEt <- .ina_slice_node(PropaguleEstablishment, tt, n, Ntimesteps, "PropaguleEstablishment")
      Dt <- .ina_slice_node(DetectionProb, tt, n, Ntimesteps, "DetectionProb")
      At <- .ina_slice_node(ManageProb, tt, n, Ntimesteps, "ManageProb")
      Mt <- .ina_slice_node(MortalityProb, tt, n, Ntimesteps, "MortalityProb")
      Rt <- .ina_slice_node(SpreadReduction, tt, n, Ntimesteps, "SpreadReduction")
      IRt <- .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb")
      IPt <- .ina_slice_node(InfoPersistenceSteps, tt, n, Ntimesteps, "InfoPersistenceSteps")
      if (programmed_dynamic) {
        G0 <- meta_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, 0, 0, 0)
        GH <- meta_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, At, Mt, Rt)
        operators[[tt]] <- .ina_programmed_information_operator(
          G0, GH, seq_len(n), Dt, SEAM, IRt, IPt, programmed_layout,
          binary_known_presence = FALSE)
      } else if (mode == "dynamic")
        operators[[tt]] <- meta_detection_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, Dt, At, Mt, Rt, SEAM, IRt)
      else if (mode == "all_informed")
        operators[[tt]] <- meta_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, At, Mt, Rt)
      else
        operators[[tt]] <- meta_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, 0, 0, 0)

      if (.ina_export_available(Model, NULL, ExportSDDprob, ExportLDDprob, AssumeResidualExport)) {
        ESDDt <- .ina_slice_export(ExportSDDprob, tt, n, Ntimesteps, "ExportSDDprob")
        ELDDt <- .ina_slice_export(ExportLDDprob, tt, n, Ntimesteps, "ExportLDDprob")
        if (mode == "dynamic") {
          e <- meta_detection_export_vector(SDDt, LDDt, LDDrate, St, Pt, At, Mt, Rt, ESDDt, ELDDt)
          if (programmed_dynamic) e <- .ina_programmed_expand_vector(e, programmed_layout)
        } else if (mode == "all_informed") e <- meta_export_vector(SDDt, LDDt, LDDrate, St, Pt, At, Mt, Rt, ESDDt, ELDDt)
        else e <- meta_export_vector(SDDt, LDDt, LDDrate, St, Pt, 0, 0, 0,
                                     ExportSDDprob = ESDDt, ExportLDDprob = ELDDt)
        pout_node <- .ina_recycle(OutsideEstablishmentProb, n, "OutsideEstablishmentProb")
        if (programmed_dynamic) {
          pout2 <- rep(pout_node, 2)
          pout <- .ina_programmed_expand_vector(pout2, programmed_layout)
        } else if (dynamic) pout <- rep(pout_node, 2) else pout <- pout_node
        export_vectors[[tt]] <- e * pout
      }
    }
  }

  if (Model == "INApestMetaTransitionMatrix") {
    if (is.null(Transition)) stop("Transition is required for INApestMetaTransitionMatrix")
    if (is.null(Nstages)) Nstages <- if (is.list(Transition)) nrow(as.matrix(Transition[[1]])) else nrow(as.matrix(Transition))
    S <- as.integer(Nstages)
    if (is.null(K)) K <- 1
    if (is.null(SeedbankK)) SeedbankK <- K
    X0 <- if (is.matrix(InitialState)) as.matrix(InitialState) else matrix(as.numeric(InitialState), nrow = n, ncol = S, byrow = TRUE)
    if (!all(dim(X0) == c(n, S))) stop("InitialState must be nodes x stages for INApestMetaTransitionMatrix")
    base_initial <- as.vector(t(X0))
    D1 <- .ina_slice_stage(DetectionProb, 1, n, S, Ntimesteps, "DetectionProb")
    if (length(D1) == 1L) D1m <- matrix(D1, n, S) else if (length(D1) == S) D1m <- matrix(rep(D1, each = n), n, S) else D1m <- as.matrix(D1)
    pdet0 <- if (ApplyInitialDetection && dynamic) .ina_initial_detection_transition(X0, D1m) else rep(0, n)
    type_node <- rep(seq_len(n), each = S)
    if (programmed_dynamic) {
      programmed_layout <- .ina_programmed_layout(n * S, persistence_max_age)
      state0 <- .ina_programmed_initial_state(base_initial, InitialInfo, pdet0,
                                               type_node, programmed_layout,
                                               binary_known_presence = FALSE)
    } else {
      info0 <- .ina_recycle(InitialInfo, n, "InitialInfo")
      if (ApplyInitialDetection && dynamic) info0 <- info0 + (1 - info0) * pdet0
      info0 <- pmin(1, pmax(0, info0))
      state0 <- if (dynamic) .ina_split_information(base_initial, info0, n, S) else base_initial
    }
    type_count <- length(state0)

    if (!dynamic && length(operators) && length(dim(SDDprob)) != 3L && !is.list(Transition)) {
      # Filled below after mortality slices are available; retained as node-level vector.
      intrinsic <- numeric(n)
    }

    for (tt in seq_len(Ntimesteps)) {
      SDDt <- .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob")
      LDDt <- .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob")
      Et <- .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb")
      PEt <- .ina_slice_node(PropaguleEstablishment, tt, n, Ntimesteps, "PropaguleEstablishment")
      Dt <- .ina_slice_stage(DetectionProb, tt, n, S, Ntimesteps, "DetectionProb")
      At <- .ina_slice_node(ManageProb, tt, n, Ntimesteps, "ManageProb")
      Mt <- .ina_slice_stage(MortalityProb, tt, n, S, Ntimesteps, "MortalityProb")
      Rt <- .ina_slice_node(SpreadReduction, tt, n, Ntimesteps, "SpreadReduction")
      IRt <- .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb")
      IPt <- .ina_slice_node(InfoPersistenceSteps, tt, n, Ntimesteps, "InfoPersistenceSteps")
      Kt <- .ina_slice_node(K, tt, n, Ntimesteps, "K")
      SBt <- .ina_slice_node(SeedbankK, tt, n, Ntimesteps, "SeedbankK")
      if (programmed_dynamic) {
        G0 <- transition_operator(Transition, S, SDDt, LDDt, LDDrate, Et, PEt, 0, 0, 0, DispersalDensityFactor, Kt, SBt)
        GH <- transition_operator(Transition, S, SDDt, LDDt, LDDrate, Et, PEt, At, Mt, Rt, DispersalDensityFactor, Kt, SBt)
        Dm <- if (length(Dt) == 1L) matrix(Dt, n, S) else if (length(Dt) == S) matrix(rep(Dt, each = n), n, S) else as.matrix(Dt)
        operators[[tt]] <- .ina_programmed_information_operator(
          G0, GH, type_node, as.vector(t(Dm)), SEAM, IRt, IPt, programmed_layout,
          binary_known_presence = FALSE)
      } else if (mode == "dynamic")
        operators[[tt]] <- transition_detection_operator(Transition, S, SDDt, LDDt, LDDrate, Et, PEt, Dt, At, Mt, Rt, SEAM, IRt, DispersalDensityFactor, Kt, SBt)
      else if (mode == "all_informed")
        operators[[tt]] <- transition_operator(Transition, S, SDDt, LDDt, LDDrate, Et, PEt, At, Mt, Rt, DispersalDensityFactor, Kt, SBt)
      else
        operators[[tt]] <- transition_operator(Transition, S, SDDt, LDDt, LDDrate, Et, PEt, 0, 0, 0, DispersalDensityFactor, Kt, SBt)

      if (.ina_export_available(Model, NULL, ExportSDDprob, ExportLDDprob, TRUE)) {
        ESDDt <- .ina_slice_export(ExportSDDprob, tt, n, Ntimesteps, "ExportSDDprob")
        ELDDt <- .ina_slice_export(ExportLDDprob, tt, n, Ntimesteps, "ExportLDDprob")
        if (mode == "dynamic") {
          e <- transition_detection_export_vector(Transition, S, SDDt, LDDt, LDDrate, At, Mt, Rt,
                                                    DispersalDensityFactor, ESDDt, ELDDt)
          if (programmed_dynamic) e <- .ina_programmed_expand_vector(e, programmed_layout)
        } else if (mode == "all_informed") e <- transition_export_vector(Transition, S, SDDt, LDDt, LDDrate, At, Mt, Rt,
                                                                           DispersalDensityFactor, ESDDt, ELDDt)
        else e <- transition_export_vector(Transition, S, SDDt, LDDt, LDDrate, 0, 0, 0,
                                            DispersalDensityFactor, ESDDt, ELDDt)
        pout_base <- rep(.ina_recycle(OutsideEstablishmentProb, n, "OutsideEstablishmentProb"), each = S)
        if (programmed_dynamic) {
          pout <- .ina_programmed_expand_vector(c(pout_base, pout_base), programmed_layout)
        } else if (dynamic) pout <- c(pout_base, pout_base) else pout <- pout_base
        export_vectors[[tt]] <- e * pout
      }
    }

    if (!is.list(Transition) && length(dim(Transition)) <= 2L &&
        !.ina_is_temporal(MortalityProb, n, Ntimesteps, "stage") &&
        !.ina_is_temporal(ManageProb, n, Ntimesteps, "node")) {
      A <- as.matrix(Transition)
      M0 <- .ina_slice_stage(MortalityProb, 1, n, S, Ntimesteps, "MortalityProb")
      if (length(M0) == 1L) Mmat <- matrix(M0, n, S) else if (length(M0) == S) Mmat <- matrix(rep(M0, each = n), n, S) else Mmat <- as.matrix(M0)
      Avec <- .ina_recycle(.ina_slice_node(ManageProb, 1, n, Ntimesteps, "ManageProb"), n, "ManageProb")
      if (mode == "none") Avec[] <- 0
      intrinsic <- vapply(seq_len(n), function(i) {
        max(Mod(eigen(A %*% diag(1 - Avec[i] * Mmat[i, ]), only.values = TRUE)$values))
      }, numeric(1))
    }
  }

  if (Model == "INApestMetaMultipleLandUse") {
    if (is.null(K) || is.null(PropaguleProduction)) stop("K and PropaguleProduction are required for INApestMetaMultipleLandUse")
    K0 <- if (length(dim(K)) == 3L) K[, , 1] else as.matrix(K); L <- ncol(K0)
    X0 <- if (is.matrix(InitialState)) as.matrix(InitialState) else matrix(as.numeric(InitialState), nrow = n, ncol = L, byrow = TRUE)
    if (!all(dim(X0) == c(n, L))) stop("InitialState must be nodes x land uses")
    base_initial <- as.vector(t(X0))
    D1 <- .ina_slice_mlu(DetectionProb, 1, n, L, Ntimesteps, "DetectionProb")
    if (length(D1) == 1L) D1m <- matrix(D1, n, L) else if (length(D1) == L) D1m <- matrix(rep(D1, each = n), n, L) else D1m <- as.matrix(D1)
    pdet0 <- if (ApplyInitialDetection && dynamic) .ina_initial_detection_mlu(X0, D1m) else rep(0, n)
    type_node <- rep(seq_len(n), each = L)
    if (programmed_dynamic) {
      programmed_layout <- .ina_programmed_layout(n * L, persistence_max_age)
      state0 <- .ina_programmed_initial_state(base_initial, InitialInfo, pdet0,
                                               type_node, programmed_layout,
                                               binary_known_presence = FALSE)
    } else {
      info0 <- .ina_recycle(InitialInfo, n, "InitialInfo")
      if (ApplyInitialDetection && dynamic) info0 <- info0 + (1 - info0) * pdet0
      info0 <- pmin(1, pmax(0, info0))
      state0 <- if (dynamic) .ina_split_information(base_initial, info0, n, L) else base_initial
    }
    type_count <- length(state0)

    for (tt in seq_len(Ntimesteps)) {
      SDDt <- .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob")
      LDDt <- .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob")
      Et <- .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb")
      St <- .ina_slice_node(Survival, tt, n, Ntimesteps, "Survival")
      Kt <- .ina_slice_K_mlu(K, tt, n, L, Ntimesteps)
      Pt <- .ina_slice_node(PropaguleProduction, tt, n, Ntimesteps, "PropaguleProduction")
      PEt <- .ina_slice_node(PropaguleEstablishment, tt, n, Ntimesteps, "PropaguleEstablishment")
      Dt <- .ina_slice_mlu(DetectionProb, tt, n, L, Ntimesteps, "DetectionProb")
      At <- .ina_slice_mlu(ManageProb, tt, n, L, Ntimesteps, "ManageProb")
      Mt <- .ina_slice_mlu(MortalityProb, tt, n, L, Ntimesteps, "MortalityProb")
      Rt <- .ina_slice_mlu(SpreadReduction, tt, n, L, Ntimesteps, "SpreadReduction")
      IRt <- .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb")
      IPt <- .ina_slice_node(InfoPersistenceSteps, tt, n, Ntimesteps, "InfoPersistenceSteps")
      if (programmed_dynamic) {
        G0 <- mlu_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, 0, 0, 0, FALSE)
        GH <- mlu_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, At, Mt, Rt, FALSE)
        Dm <- if (length(Dt) == 1L) matrix(Dt, n, L) else if (length(Dt) == L) matrix(rep(Dt, each = n), n, L) else as.matrix(Dt)
        operators[[tt]] <- .ina_programmed_information_operator(
          G0, GH, type_node, as.vector(t(Dm)), SEAM, IRt, IPt, programmed_layout,
          binary_known_presence = FALSE)
      } else if (mode == "dynamic")
        operators[[tt]] <- mlu_detection_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, Dt, At, Mt, Rt, SEAM, IRt, FALSE)
      else if (mode == "all_informed")
        operators[[tt]] <- mlu_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, At, Mt, Rt, FALSE)
      else
        operators[[tt]] <- mlu_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, 0, 0, 0, FALSE)

      if (.ina_export_available(Model, NULL, ExportSDDprob, ExportLDDprob, AssumeResidualExport)) {
        ESDDt <- .ina_slice_export(ExportSDDprob, tt, n, Ntimesteps, "ExportSDDprob")
        ELDDt <- .ina_slice_export(ExportLDDprob, tt, n, Ntimesteps, "ExportLDDprob")
        if (mode == "dynamic") {
          e <- mlu_detection_export_vector(SDDt, LDDt, LDDrate, St, Kt, Pt, At, Mt, Rt, FALSE,
                                            ESDDt, ELDDt)
          if (programmed_dynamic) e <- .ina_programmed_expand_vector(e, programmed_layout)
        } else if (mode == "all_informed") e <- mlu_export_vector(SDDt, LDDt, LDDrate, St, Kt, Pt, At, Mt, Rt, FALSE,
                                                                    ESDDt, ELDDt)
        else e <- mlu_export_vector(SDDt, LDDt, LDDrate, St, Kt, Pt, 0, 0, 0, FALSE,
                                    ESDDt, ELDDt)
        pout_base <- rep(.ina_recycle(OutsideEstablishmentProb, n, "OutsideEstablishmentProb"), each = L)
        if (programmed_dynamic) {
          pout <- .ina_programmed_expand_vector(c(pout_base, pout_base), programmed_layout)
        } else if (dynamic) pout <- rep(pout_base, 2) else pout <- pout_base
        export_vectors[[tt]] <- e * pout
      }
    }
  }

  growth <- .ina_cycle_growth(operators)
  Rstep <- growth$EquivalentPerTimestepMultiplier
  trajectory_state <- .ina_apply_operators(operators, state0)
  totals <- if (programmed_dynamic) colSums(trajectory_state) else .ina_state_total(trajectory_state, dynamic)
  trajectory <- data.frame(timestep = seq_len(Ntimesteps), expected_state_total = totals)

  # Extinction: use the most simulator-specific PGF available for static simple
  # cases; otherwise use a transparent mean-matched multitype branching fallback.
  static_ops <- all(vapply(operators[-1], function(x) isTRUE(all.equal(x, operators[[1]], tolerance = 0)), logical(1)))
  if (Ntimesteps == 1L) static_ops <- TRUE
  extinction_method <- NULL; qh <- NULL; qe <- NULL
  if (static_ops && !dynamic && Model == "INApest") {
    pars <- list(SDDprob = .ina_slice_connection(SDDprob, 1, n, Ntimesteps, "SDDprob"),
                 LDDprob = .ina_slice_connection(LDDprob, 1, n, Ntimesteps, "LDDprob"),
                 EnvEstabProb = .ina_slice_node(EnvEstabProb, 1, n, Ntimesteps, "EnvEstabProb"),
                 Survival = .ina_slice_node(Survival, 1, n, Ntimesteps, "Survival"),
                 ManageProb = if (mode == "all_informed") .ina_slice_node(ManageProb, 1, n, Ntimesteps, "ManageProb") else 0,
                 EradicationProb = if (mode == "all_informed") .ina_slice_node(EradicationProb, 1, n, Ntimesteps, "EradicationProb") else 0,
                 SpreadReduction = if (mode == "all_informed") .ina_slice_node(SpreadReduction, 1, n, Ntimesteps, "SpreadReduction") else 0,
                 generations = max(ExtinctionGenerations, Ntimesteps))
    ex <- do.call(inapest_exogenous_extinction, pars); qh <- ex$history[, Ntimesteps]; qe <- ex$extinction
    extinction_method <- "binary edge-based multitype branching PGF"
  } else if (static_ops && dynamic && Model == "INApest" && !programmed_dynamic) {
    ex <- inapest_detection_extinction(.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),
                                       .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),
                                       .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                       .ina_slice_node(Survival,1,n,Ntimesteps,"Survival"),
                                       .ina_slice_node(DetectionProb,1,n,Ntimesteps,"DetectionProb"),
                                       .ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb"),
                                       .ina_slice_node(EradicationProb,1,n,Ntimesteps,"EradicationProb"),
                                       .ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction"),
                                       SEAM, .ina_slice_node(InfoRetentionProb,1,n,Ntimesteps,"InfoRetentionProb"),
                                       generations = max(ExtinctionGenerations,Ntimesteps))
    qh <- c(ex$historyU[, Ntimesteps], ex$historyH[, Ntimesteps]); qe <- c(ex$U, ex$H)
    extinction_method <- "binary informed/uninformed branching PGF"
  } else if (static_ops && !dynamic && Model == "INApestMeta") {
    ex <- meta_single_parent_extinction(.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),
                                        .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),LDDrate,
                                        .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                        .ina_slice_node(Survival,1,n,Ntimesteps,"Survival"),
                                        .ina_slice_node(K,1,n,Ntimesteps,"K"),
                                        .ina_slice_node(PropaguleProduction,1,n,Ntimesteps,"PropaguleProduction"),
                                        .ina_slice_node(PropaguleEstablishment,1,n,Ntimesteps,"PropaguleEstablishment"),
                                        if(mode=="all_informed").ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb") else 0,
                                        if(mode=="all_informed").ina_slice_node(MortalityProb,1,n,Ntimesteps,"MortalityProb") else 0,
                                        if(mode=="all_informed").ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction") else 0,
                                        generations=max(ExtinctionGenerations,Ntimesteps))
    qh <- ex$history[, Ntimesteps]; qe <- ex$extinction; extinction_method <- "simulator-faithful one-parent Meta branching PGF"
  } else if (static_ops && !dynamic && Model == "INApestMetaTransitionMatrix") {
    ex <- transition_extinction(Transition,Nstages,.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),
                                .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),LDDrate,
                                .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                .ina_slice_node(PropaguleEstablishment,1,n,Ntimesteps,"PropaguleEstablishment"),
                                if(mode=="all_informed").ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb") else 0,
                                if(mode=="all_informed").ina_slice_stage(MortalityProb,1,n,Nstages,Ntimesteps,"MortalityProb") else 0,
                                if(mode=="all_informed").ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction") else 0,
                                DispersalDensityFactor,
                                .ina_slice_node(K,1,n,Ntimesteps,"K"),
                                .ina_slice_node(SeedbankK,1,n,Ntimesteps,"SeedbankK"),
                                generations=max(ExtinctionGenerations,Ntimesteps))
    qh <- ex$history[, Ntimesteps]; qe <- ex$extinction; extinction_method <- "stage x node branching PGF"
  } else if (static_ops && !dynamic && Model == "INApestMetaMultipleLandUse") {
    ex <- mlu_single_parent_extinction(SDDprob=.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),
                                       LDDprob=.ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),LDDrate=LDDrate,
                                       EnvEstabProb=.ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                       Survival=.ina_slice_node(Survival,1,n,Ntimesteps,"Survival"),
                                       K=.ina_slice_K_mlu(K,1,n,L,Ntimesteps),
                                       PropaguleProduction=.ina_slice_node(PropaguleProduction,1,n,Ntimesteps,"PropaguleProduction"),
                                       PropaguleEstablishment=.ina_slice_node(PropaguleEstablishment,1,n,Ntimesteps,"PropaguleEstablishment"),
                                       ManageProb=if(mode=="all_informed").ina_slice_mlu(ManageProb,1,n,L,Ntimesteps,"ManageProb") else 0,
                                       MortalityProb=if(mode=="all_informed").ina_slice_mlu(MortalityProb,1,n,L,Ntimesteps,"MortalityProb") else 0,
                                       SpreadReduction=if(mode=="all_informed").ina_slice_mlu(SpreadReduction,1,n,L,Ntimesteps,"SpreadReduction") else 0,
                                       generations=max(ExtinctionGenerations,Ntimesteps),current_code=FALSE)
    qh <- ex$history[, Ntimesteps]; qe <- ex$extinction; extinction_method <- "mean-matched MLU branching approximation"
  } else {
    qh <- ina_poisson_branching_extinction_horizon(operators)
    if (static_ops) qe <- ina_poisson_branching_extinction(operators[[1]], ExtinctionGenerations)$extinction
    extinction_method <- if (programmed_dynamic) {
      if (static_ops) "age-structured mean-matched multitype branching" else "time-inhomogeneous age-structured mean-matched multitype branching"
    } else if (static_ops) "mean-matched multitype Poisson branching fallback" else "time-inhomogeneous mean-matched multitype Poisson branching"
  }
  extinction_horizon <- .ina_overall_extinction(qh, state0)
  extinction_eventual <- if (!is.null(qe)) .ina_overall_extinction(qe, state0) else NA_real_

  escape <- NULL
  if (.ina_export_available(Model, ExportProb, ExportSDDprob, ExportLDDprob,
                            if(Model=="INApestMetaTransitionMatrix") TRUE else AssumeResidualExport)) {
    first <- ina_escape_first_moment_temporal(operators, state0, export_vectors)
    noesc <- NULL; escape_method <- NULL
    if (static_ops && !dynamic && Model == "INApest" && !is.null(ExportProb)) {
      X <- .ina_slice_export(ExportProb,1,n,Ntimesteps,"ExportProb"); pout <- .inapest_recycle(OutsideEstablishmentProb,n,"OutsideEstablishmentProb")
      Xeff <- sweep(X,1,pout,`*`)
      br <- inapest_escape_branching(.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),Xeff,
                                     .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),
                                     .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                     .ina_slice_node(Survival,1,n,Ntimesteps,"Survival"),
                                     if(mode=="all_informed").ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb") else 0,
                                     if(mode=="all_informed").ina_slice_node(EradicationProb,1,n,Ntimesteps,"EradicationProb") else 0,
                                     if(mode=="all_informed").ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction") else 0,
                                     Ntimesteps)
      noesc <- br$no_escape; escape_method <- "binary edge-based branching no-escape recursion"
    } else if (static_ops && dynamic && Model == "INApest" && !programmed_dynamic && !is.null(ExportProb)) {
      X <- .ina_slice_export(ExportProb,1,n,Ntimesteps,"ExportProb"); pout <- .inapest_recycle(OutsideEstablishmentProb,n,"OutsideEstablishmentProb")
      Xeff <- sweep(X,1,pout,`*`)
      br <- inapest_detection_escape_branching(.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),Xeff,
                                                .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),
                                                .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                                .ina_slice_node(Survival,1,n,Ntimesteps,"Survival"),
                                                .ina_slice_node(DetectionProb,1,n,Ntimesteps,"DetectionProb"),
                                                .ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb"),
                                                .ina_slice_node(EradicationProb,1,n,Ntimesteps,"EradicationProb"),
                                                .ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction"),SEAM,
                                                .ina_slice_node(InfoRetentionProb,1,n,Ntimesteps,"InfoRetentionProb"),Ntimesteps)
      noesc <- c(br$no_escape_U,br$no_escape_H); escape_method <- "binary informed/uninformed branching no-escape recursion"
    } else if (static_ops && !dynamic && Model == "INApestMeta") {
      br <- meta_single_parent_escape_branching(.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),
                                                .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),LDDrate,
                                                .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                                .ina_slice_node(Survival,1,n,Ntimesteps,"Survival"),
                                                .ina_slice_node(K,1,n,Ntimesteps,"K"),
                                                .ina_slice_node(PropaguleProduction,1,n,Ntimesteps,"PropaguleProduction"),
                                                .ina_slice_node(PropaguleEstablishment,1,n,Ntimesteps,"PropaguleEstablishment"),
                                                if(mode=="all_informed").ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb") else 0,
                                                if(mode=="all_informed").ina_slice_node(MortalityProb,1,n,Ntimesteps,"MortalityProb") else 0,
                                                if(mode=="all_informed").ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction") else 0,
                                                .ina_slice_export(ExportSDDprob,1,n,Ntimesteps,"ExportSDDprob"),
                                                .ina_slice_export(ExportLDDprob,1,n,Ntimesteps,"ExportLDDprob"),
                                                OutsideEstablishmentProb,AssumeResidualExport,Ntimesteps)
      noesc <- br$no_escape; escape_method <- "Meta one-parent branching no-escape PGF"
    } else if (static_ops && !dynamic && Model == "INApestMetaTransitionMatrix") {
      br <- transition_escape_branching(Transition,Nstages,.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),
                                         .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),LDDrate,
                                         .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                         .ina_slice_node(PropaguleEstablishment,1,n,Ntimesteps,"PropaguleEstablishment"),
                                         if(mode=="all_informed").ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb") else 0,
                                         if(mode=="all_informed").ina_slice_stage(MortalityProb,1,n,Nstages,Ntimesteps,"MortalityProb") else 0,
                                         if(mode=="all_informed").ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction") else 0,
                                         DispersalDensityFactor,.ina_slice_node(K,1,n,Ntimesteps,"K"),
                                         .ina_slice_node(SeedbankK,1,n,Ntimesteps,"SeedbankK"),
                                         .ina_slice_export(ExportSDDprob,1,n,Ntimesteps,"ExportSDDprob"),
                                         .ina_slice_export(ExportLDDprob,1,n,Ntimesteps,"ExportLDDprob"),
                                         OutsideEstablishmentProb,TRUE,Ntimesteps)
      noesc <- br$no_escape; escape_method <- "stage x node branching no-escape PGF"
    } else {
      noesc <- ina_poisson_branching_no_escape_horizon(operators, export_vectors)
      escape_method <- if (programmed_dynamic) {
        if (static_ops) "age-structured mean-matched branching no-escape recursion" else "time-inhomogeneous age-structured mean-matched branching no-escape recursion"
      } else if (static_ops) "mean-matched multitype branching no-escape recursion" else "time-inhomogeneous mean-matched multitype branching no-escape recursion"
    }
    branch_prob <- 1 - .ina_overall_no_escape(noesc,state0)
    escape <- list(Method = escape_method,
                   BranchingProbabilityByHorizon = branch_prob,
                   FirstMomentPoissonProbabilityByHorizon = first$EscapeProbability,
                   ExpectedSuccessfulEscapesByStep = first$ExpectedSuccessfulEscapesByStep,
                   ExpectedSuccessfulEscapesCumulative = first$ExpectedSuccessfulEscapesCumulative)
  } else {
    diagnostics <- c(diagnostics,
      if (Model == "INApest")
        "Escape not calculated: binary INApest requires explicit source x outside-destination ExportProb; outside probability cannot be recovered from the internal edge matrix."
      else
        "Escape not calculated: provide ExportSDDprob/ExportLDDprob, or set AssumeResidualExport=TRUE only when missing row mass genuinely represents outside dispersal.")
  }

  if (dynamic && Model != "INApest") diagnostics <- c(diagnostics,
    "For information-limited Meta/Transition/MLU models, branching extinction/escape may use a mean-matched multitype fallback because node-level shared information induces correlations not represented by independent individual types.")
  if (programmed_dynamic && Model != "INApest") diagnostics <- c(diagnostics,
    "The programmed information clock itself is represented explicitly. At low density, however, a management kill that removes the focal individual cannot also leave a surviving focal lineage; information retained by that kill can precondition later immigrants and is therefore a shared-node, higher-order effect best captured by simulation or a nonlinear node-level extension.")
  if (.ina_any_nonzero(DispersalDensityFactor)) diagnostics <- c(diagnostics,
    "DispersalDensityFactor is linearized at the pest-free state. It can strongly alter finite-density trajectories even when the rare-state multiplier is unchanged.")
  diagnostics <- unique(diagnostics)

  result <- list(
    Model = Model,
    InformationMode = mode,
    InformationStateMethod = information_state_method,
    Growth = list(
      EquivalentPerTimestepMultiplier = Rstep,
      CycleMultiplier = growth$CycleMultiplier,
      Classification = .ina_classify_growth(Rstep),
      IntrinsicLocalLambda = intrinsic,
      IntrinsicLocalLambdaInterpretation = if (Model == "INApestMetaTransitionMatrix" && !is.null(intrinsic)) {
        if (mode == "all_informed") "managed local transition multiplier" else if (mode == "none") "unmanaged local transition multiplier" else "conditional-on-information managed local transition multiplier; landscape information limitation is represented separately in the augmented information-state operator"
      } else NULL
    ),
    Trajectory = trajectory,
    Extinction = list(
      Method = extinction_method,
      ProbabilityByHorizon = extinction_horizon,
      BranchingFadeoutProbability = extinction_eventual
    ),
    Escape = escape,
    Diagnostics = diagnostics,
    ApproximationScope = c(
      "Low-density / rare-invasion screening unless otherwise stated",
      "Expected-parameter treatment of management/detection SD rather than annual random-effect integration",
      "Finite carrying capacity, lineage collisions and density dependence become increasingly important away from rarity",
      if (programmed_dynamic) "Programmed stopping uses explicit time-since-local-evidence states; shared-node evidence created by other individuals and information-only preconditioning remain approximations away from rarity" else "Memoryless information retention is represented with uninformed/informed states"
    )
  )
  if (ReturnOperators) result$Operators <- operators
  class(result) <- "INApestAnalyticalResult"
  result
}

print.INApestAnalyticalResult <- function(x, ...) {
  cat("INApest analytical screening result\n")
  cat("  Model:", x$Model, "\n")
  cat("  Information mode:", x$InformationMode, "\n")
  if (!is.null(x$InformationStateMethod)) cat("  Information state method:", x$InformationStateMethod, "\n")
  cat("  Growth multiplier per timestep:", format(round(x$Growth$EquivalentPerTimestepMultiplier, 4), nsmall = 4), "\n")
  cat("  Classification:", x$Growth$Classification, "\n")
  if (!is.null(x$Growth$IntrinsicLocalLambda))
    cat("  Intrinsic local lambda range:", paste(round(range(x$Growth$IntrinsicLocalLambda), 4), collapse = " to "), "\n")
  cat("  Expected state total at horizon:", round(tail(x$Trajectory$expected_state_total, 1), 4), "\n")
  cat("  Branching extinction by horizon:", round(x$Extinction$ProbabilityByHorizon, 4), "\n")
  if (!is.null(x$Escape)) {
    cat("  Branching escape by horizon:", round(x$Escape$BranchingProbabilityByHorizon, 4), "\n")
    cat("  First-moment escape comparator:", round(x$Escape$FirstMomentPoissonProbabilityByHorizon, 4), "\n")
  }
  if (length(x$Diagnostics)) {
    cat("  Diagnostics:", length(x$Diagnostics), "(see $Diagnostics)\n")
  }
  invisible(x)
}

summary.INApestAnalyticalResult <- function(object, ...) {
  data.frame(
    Model = object$Model,
    InformationMode = object$InformationMode,
    GrowthMultiplier = object$Growth$EquivalentPerTimestepMultiplier,
    GrowthClassification = object$Growth$Classification,
    ExpectedTotalAtHorizon = tail(object$Trajectory$expected_state_total, 1),
    ExtinctionByHorizon = object$Extinction$ProbabilityByHorizon,
    EscapeByHorizon = if (is.null(object$Escape)) NA_real_ else object$Escape$BranchingProbabilityByHorizon,
    FirstMomentEscape = if (is.null(object$Escape)) NA_real_ else object$Escape$FirstMomentPoissonProbabilityByHorizon,
    stringsAsFactors = FALSE
  )
}
###############################################################################
### Round 2 analytical extension: management-induced fecundity reduction
###
### This module is designed to be sourced after the 14-Aug-2026 analytical
### baseline.  The baseline functions are retained under *_pre_fecundity names.
### Every replacement below delegates to the original function when
### FecundityReduction is identically zero, providing an explicit regression
### path rather than merely relying on algebraic equivalence.
###############################################################################

.ina_fr_nonzero <- function(x) {
  if (is.null(x) || !length(x)) return(FALSE)
  .ina_fr_validate(x)
  any(as.numeric(x) != 0)
}

.ina_fr_validate <- function(x, name = "FecundityReduction") {
  if (is.function(x) || !is.numeric(x))
    stop(name, " must be numeric")
  v <- as.numeric(x)
  if (!length(v) || any(!is.finite(v)) || any(v < 0) || any(v > 1))
    stop(name, " values must be finite and between 0 and 1")
  invisible(TRUE)
}

.ina_transition_fecundity_matrix <- function(x, n, S,
                                              name = "FecundityReduction") {
  .ina_fr_validate(x, name)
  dx <- dim(x)
  if (is.null(dx)) {
    if (length(x) == 1L) return(matrix(as.numeric(x), n, S))
    if (n == S && length(x) == n)
      stop(name, " vector is ambiguous because number of nodes equals Nstages; use a dimensioned form")
    if (length(x) == n)
      return(matrix(rep(as.numeric(x), S), nrow = n, ncol = S))
    if (length(x) == S)
      return(matrix(rep(as.numeric(x), each = n), nrow = n, ncol = S))
    stop(name, " must be scalar, length nodes, length Nstages, or nodes x Nstages")
  }
  if (!identical(dx, c(n, S))) stop(name, " matrix must have dimensions nodes x Nstages")
  matrix(as.numeric(x), nrow = n, ncol = S)
}

.ina_fecundity_transition_is_temporal <- function(x, n, S, T) {
  dx <- dim(x)
  if (length(dx) == 3L) return(TRUE)
  if (length(dx) == 2L && length(dx) == 2L && all(dx == c(n, T))) return(TRUE)
  FALSE
}

.ina_slice_fecundity_transition <- function(x, t, n, S, T,
                                             name = "FecundityReduction") {
  .ina_fr_validate(x, name)
  dx <- dim(x)
  if (length(dx) == 3L) {
    if (length(dx) != 3L || !all(dx == c(n, S, T)))
      stop(name, " 3D array must have dimensions nodes x stages x Ntimesteps")
    return(x[, , t, drop = TRUE])
  }
  if (length(dx) == 2L) {
    # Match the public INApestMetaTransitionMatrix contract: a 2-D
    # FecundityReduction matrix is nodes x Ntimesteps.  A static nodes x
    # stages matrix is deliberately not accepted here because that shape is
    # not a simulator input form; use a length-Nstages vector for a static
    # stage effect or a nodes x stages x Ntimesteps array for both.
    if (all(dx == c(n, T)))
      return(matrix(rep(as.numeric(x[, t]), S), nrow = n, ncol = S))
    stop(name, " 2D matrix must have dimensions nodes x Ntimesteps")
  }
  .ina_transition_fecundity_matrix(x, n, S, name)
}

###############################################################################
### Preserve baseline implementations
###############################################################################
.meta_components_pre_fecundity <- meta_components
.meta_operator_pre_fecundity <- meta_operator
.meta_meanfield_pre_fecundity <- meta_meanfield
.meta_extinction_pre_fecundity <- meta_extinction
.meta_detection_operator_pre_fecundity <- meta_detection_operator
.meta_expected_exports_pre_fecundity <- meta_expected_exports
.meta_single_parent_operator_pre_fecundity <- meta_single_parent_operator
.meta_detection_single_parent_operator_pre_fecundity <- meta_detection_single_parent_operator
.meta_single_parent_recruit_means_pre_fecundity <- meta_single_parent_recruit_means
.meta_single_parent_extinction_pre_fecundity <- meta_single_parent_extinction
.meta_export_vector_pre_fecundity <- meta_export_vector
.meta_detection_export_vector_pre_fecundity <- meta_detection_export_vector
.meta_single_parent_escape_branching_pre_fecundity <- meta_single_parent_escape_branching

.transition_components_pre_fecundity <- transition_components
.transition_operator_pre_fecundity <- transition_operator
.transition_intrinsic_managed_lambda_pre_fecundity <- transition_intrinsic_managed_lambda
.transition_extinction_pre_fecundity <- transition_extinction
.transition_detection_operator_pre_fecundity <- transition_detection_operator
.transition_export_vector_pre_fecundity <- transition_export_vector
.transition_detection_export_vector_pre_fecundity <- transition_detection_export_vector
.transition_escape_branching_pre_fecundity <- transition_escape_branching

.mlu_components_pre_fecundity <- mlu_components
.mlu_operator_pre_fecundity <- mlu_operator
.mlu_meanfield_pre_fecundity <- mlu_meanfield
.mlu_detection_operator_pre_fecundity <- mlu_detection_operator
.mlu_single_parent_operator_pre_fecundity <- mlu_single_parent_operator
.mlu_detection_single_parent_operator_pre_fecundity <- mlu_detection_single_parent_operator
.mlu_single_parent_extinction_pre_fecundity <- mlu_single_parent_extinction
.mlu_export_vector_pre_fecundity <- mlu_export_vector
.mlu_detection_export_vector_pre_fecundity <- mlu_detection_export_vector

###############################################################################
### Ordinary INApestMeta
###############################################################################
meta_components <- function(SDDprob, LDDprob = 0, LDDrate = 0,
                            EnvEstabProb = 1, Survival = 1, K,
                            PropaguleProduction,
                            PropaguleEstablishment,
                            ManageProb = 0, MortalityProb = 0,
                            SpreadReduction = 0,
                            FecundityReduction = 0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.meta_components_pre_fecundity(
      SDDprob, LDDprob, LDDrate, EnvEstabProb, Survival, K,
      PropaguleProduction, PropaguleEstablishment, ManageProb,
      MortalityProb, SpreadReduction))
  SDD <- as.matrix(SDDprob); n <- nrow(SDD)
  if (ncol(SDD) != n) stop("SDDprob must be square")
  LDD <- .ina_mat(LDDprob, n, "LDDprob")
  r <- as.numeric(LDDrate)
  if (length(r) != 1L || !is.finite(r) || r < 0 || r > 1) stop("LDDrate must be in [0,1]")
  env <- .ina_recycle(EnvEstabProb,n,"EnvEstabProb")
  surv <- .ina_recycle(Survival,n,"Survival")
  cap <- .ina_recycle(K,n,"K")
  prod <- .ina_recycle(PropaguleProduction,n,"PropaguleProduction")
  pest <- .ina_recycle(PropaguleEstablishment,n,"PropaguleEstablishment")
  a <- .ina_recycle(ManageProb,n,"ManageProb")
  m <- .ina_recycle(MortalityProb,n,"MortalityProb")
  g <- .ina_recycle(SpreadReduction,n,"SpreadReduction")
  f <- .ina_recycle(FecundityReduction,n,"FecundityReduction"); .ina_fr_validate(f)
  q0 <- surv; q1 <- surv * (1 - m); qbar <- (1-a)*q0 + a*q1
  K0 <- (1-r)*SDD + r*LDD
  K1 <- (1-r)*SDD + r*sweep(LDD,1,1-g,`*`)
  Arr <- sweep(K0,1,(1-a)*q0*prod,`*`) +
         sweep(K1,1,a*q1*prod*(1-f),`*`)
  alpha <- pest * env; c_est <- 1-exp(-alpha)
  list(n=n,SDD=SDD,LDD=LDD,LDDrate=r,env=env,survival=surv,K=cap,
       production=prod,prop_est=pest,adoption=a,mortality=m,
       spread_reduction=g,fecundity_reduction=f,q0=q0,q1=q1,qbar=qbar,
       kernel0=K0,kernel1=K1,arrivals=Arr,alpha=alpha,recruit_slope=c_est)
}

meta_operator <- function(SDDprob, LDDprob = 0, LDDrate = 0,
                          EnvEstabProb = 1, Survival = 1, K,
                          PropaguleProduction, PropaguleEstablishment,
                          ManageProb = 0, MortalityProb = 0,
                          SpreadReduction = 0, FecundityReduction = 0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.meta_operator_pre_fecundity(SDDprob,LDDprob,LDDrate,EnvEstabProb,
      Survival,K,PropaguleProduction,PropaguleEstablishment,ManageProb,
      MortalityProb,SpreadReduction))
  z <- meta_components(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                       PropaguleProduction,PropaguleEstablishment,ManageProb,
                       MortalityProb,SpreadReduction,FecundityReduction)
  dest_gain <- z$K*z$recruit_slope
  G <- sweep(t(z$arrivals),1,dest_gain,`*`); diag(G) <- diag(G)+z$qbar
  attr(G,"components") <- z; G
}

meta_meanfield <- function(SDDprob, LDDprob = 0, LDDrate = 0,
                           EnvEstabProb = 1, Survival = 1, K,
                           PropaguleProduction, PropaguleEstablishment,
                           ManageProb = 0, MortalityProb = 0,
                           SpreadReduction = 0, initial, timesteps,
                           FecundityReduction = 0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.meta_meanfield_pre_fecundity(SDDprob,LDDprob,LDDrate,EnvEstabProb,
      Survival,K,PropaguleProduction,PropaguleEstablishment,ManageProb,
      MortalityProb,SpreadReduction,initial,timesteps))
  z <- meta_components(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                       PropaguleProduction,PropaguleEstablishment,ManageProb,
                       MortalityProb,SpreadReduction,FecundityReduction)
  x <- as.numeric(initial); out <- matrix(NA_real_,z$n,timesteps)
  for (tt in seq_len(timesteps)) {
    n0 <- z$qbar*x; lambda <- as.numeric(x %*% z$arrivals)
    p_rec <- 1-exp(-z$recruit_slope*lambda)
    x <- n0 + pmax(0,z$K-n0)*p_rec
    x <- pmin(z$K,pmax(0,x)); out[,tt] <- x
  }
  out
}

meta_extinction <- function(SDDprob,LDDprob=0,LDDrate=0,
                            EnvEstabProb=1,Survival=1,K,
                            PropaguleProduction,PropaguleEstablishment,
                            ManageProb=0,MortalityProb=0,SpreadReduction=0,
                            generations=100,tolerance=1e-12,
                            FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.meta_extinction_pre_fecundity(SDDprob,LDDprob,LDDrate,EnvEstabProb,
      Survival,K,PropaguleProduction,PropaguleEstablishment,ManageProb,
      MortalityProb,SpreadReduction,generations,tolerance))
  z <- meta_components(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                       PropaguleProduction,PropaguleEstablishment,ManageProb,
                       MortalityProb,SpreadReduction,FecundityReduction)
  n<-z$n; gain<-z$K*z$recruit_slope
  mu0<-sweep(sweep(z$kernel0,1,z$production,`*`),2,gain,`*`)
  mu1<-sweep(sweep(z$kernel1,1,z$production*(1-z$fecundity_reduction),`*`),2,gain,`*`)
  q<-rep(0,n); hist<-matrix(NA_real_,n,generations)
  for(tt in seq_len(generations)){
    qo<-q; qn<-numeric(n)
    for(i in seq_len(n)){
      R0<-exp(sum(mu0[i,]*(q-1))); R1<-exp(sum(mu1[i,]*(q-1)))
      f0<-(1-z$q0[i])+z$q0[i]*q[i]*R0
      f1<-(1-z$q1[i])+z$q1[i]*q[i]*R1
      qn[i]<-(1-z$adoption[i])*f0+z$adoption[i]*f1
    }
    q<-pmin(1,pmax(0,qn));hist[,tt]<-q
    if(max(abs(q-qo))<tolerance){if(tt<generations)hist[,(tt+1):generations]<-q;break}
  }
  list(extinction=q,history=hist,
       note="Multitype branching approximation with management-induced fecundity reduction in the managed offspring mean")
}

meta_detection_operator <- function(SDDprob,LDDprob=0,LDDrate=0,
                                    EnvEstabProb=1,Survival=1,K,
                                    PropaguleProduction,PropaguleEstablishment,
                                    DetectionProb=0,ManageProb=0,MortalityProb=0,
                                    SpreadReduction=0,SEAM=NULL,InfoRetentionProb=1,
                                    FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.meta_detection_operator_pre_fecundity(SDDprob,LDDprob,LDDrate,
      EnvEstabProb,Survival,K,PropaguleProduction,PropaguleEstablishment,
      DetectionProb,ManageProb,MortalityProb,SpreadReduction,SEAM,InfoRetentionProb))
  n<-nrow(as.matrix(SDDprob));d<-.ina_recycle(DetectionProb,n,"DetectionProb")
  ir<-.ina_recycle(InfoRetentionProb,n,"InfoRetentionProb")
  C<-if(is.null(SEAM)||length(SEAM)==1L)matrix(0,n,n)else as.matrix(SEAM)
  if(!all(dim(C)==c(n,n)))stop("SEAM must be nodes x nodes");diag(C)<-0
  G0<-meta_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                    PropaguleProduction,PropaguleEstablishment,0,0,0,0)
  GH<-meta_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                    PropaguleProduction,PropaguleEstablishment,ManageProb,
                    MortalityProb,SpreadReduction,FecundityReduction)
  G<-matrix(0,2*n,2*n);U<-seq_len(n);H<-n+seq_len(n)
  for(i in seq_len(n))for(j in seq_len(n)){
    w<-G0[j,i];if(w!=0){h<-d[j];G[U[j],U[i]]<-G[U[j],U[i]]+w*(1-h);G[H[j],U[i]]<-G[H[j],U[i]]+w*h}
    w<-GH[j,i];if(w!=0){h<-if(j==i)ir[j]+(1-ir[j])*d[j]else 1-(1-d[j])*(1-C[i,j]);G[U[j],H[i]]<-G[U[j],H[i]]+w*(1-h);G[H[j],H[i]]<-G[H[j],H[i]]+w*h}
  }
  attr(G,"note")<-"2N individual-type low-density approximation with managed fecundity reduction"
  G
}

meta_expected_exports <- function(state,ExportSDDprob,ExportLDDprob=0,
                                  LDDrate=0,Survival=1,PropaguleProduction,
                                  ManageProb=0,MortalityProb=0,
                                  SpreadReduction=0,FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.meta_expected_exports_pre_fecundity(state,ExportSDDprob,ExportLDDprob,
      LDDrate,Survival,PropaguleProduction,ManageProb,MortalityProb,SpreadReduction))
  x<-as.numeric(state);n<-length(x);XS<-as.matrix(ExportSDDprob)
  if(nrow(XS)!=n)stop("ExportSDDprob rows must equal nodes")
  XL<-if(length(ExportLDDprob)==1L&&ExportLDDprob==0)matrix(0,n,ncol(XS))else as.matrix(ExportLDDprob)
  if(nrow(XL)!=n)stop("ExportLDDprob rows must equal nodes")
  s<-.ina_recycle(Survival,n,"Survival");p<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction")
  a<-.ina_recycle(ManageProb,n,"ManageProb");m<-.ina_recycle(MortalityProb,n,"MortalityProb")
  g<-.ina_recycle(SpreadReduction,n,"SpreadReduction");f<-.ina_recycle(FecundityReduction,n,"FecundityReduction")
  r<-LDDrate;q0<-s;q1<-s*(1-m);e0<-rowSums((1-r)*XS+r*XL)
  e1<-rowSums((1-r)*XS+r*sweep(XL,1,1-g,`*`))
  sum(x*p*((1-a)*q0*e0+a*q1*(1-f)*e1))
}

meta_single_parent_operator <- function(SDDprob,LDDprob=0,LDDrate=0,
                                        EnvEstabProb=1,Survival=1,K,
                                        PropaguleProduction,PropaguleEstablishment,
                                        ManageProb=0,MortalityProb=0,
                                        SpreadReduction=0,FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.meta_single_parent_operator_pre_fecundity(SDDprob,LDDprob,LDDrate,
      EnvEstabProb,Survival,K,PropaguleProduction,PropaguleEstablishment,
      ManageProb,MortalityProb,SpreadReduction))
  SDD<-as.matrix(SDDprob);n<-nrow(SDD);LDD<-.ina_mat(LDDprob,n,"LDDprob")
  env<-.ina_recycle(EnvEstabProb,n,"EnvEstabProb");s<-.ina_recycle(Survival,n,"Survival")
  cap<-.ina_recycle(K,n,"K");prod<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction")
  pe<-.ina_recycle(PropaguleEstablishment,n,"PropaguleEstablishment")
  a<-.ina_recycle(ManageProb,n,"ManageProb");m<-.ina_recycle(MortalityProb,n,"MortalityProb")
  g<-.ina_recycle(SpreadReduction,n,"SpreadReduction");f<-.ina_recycle(FecundityReduction,n,"FecundityReduction");.ina_fr_validate(f)
  r<-LDDrate;alpha<-pe*env;rs<-rowSums(SDD);rl<-rowSums(LDD);G<-matrix(0,n,n)
  for(i in seq_len(n)){
    ps<-if(rs[i]>0)SDD[i,]/rs[i]else rep(0,n);pl<-if(rl[i]>0)LDD[i,]/rl[i]else rep(0,n)
    for(M in 0:1){
      pm<-if(M==0)1-a[i]else a[i];if(pm==0)next
      surv<-s[i]*(1-m[i]*M);if(surv==0)next
      lambda<-prod[i]*(1-f[i]*M);supp<-.meta_poisson_support(lambda);k<-supp$k;pk<-supp$p
      ms<-floor(k*((1-r)*rs[i]));ml<-floor(k*(r*(1-g[i]*M)*rl[i]))
      free<-cap;free[i]<-pmax(0,free[i]-1)
      for(j in seq_len(n)){
        bS<-1-ps[j]*(1-exp(-alpha[j]));bL<-1-pl[j]*(1-exp(-alpha[j]))
        nohaz<-sum(pk*(bS^ms)*(bL^ml));recruits<-free[j]*(1-nohaz)
        G[j,i]<-G[j,i]+pm*surv*(as.numeric(j==i)+recruits)
      }
    }
  }
  attr(G,"note")<-"One-parent mean operator matching current integer truncation; managed Poisson propagule mean includes FecundityReduction"
  G
}

meta_detection_single_parent_operator <- function(SDDprob,LDDprob=0,LDDrate=0,
                                                   EnvEstabProb=1,Survival=1,K,
                                                   PropaguleProduction,PropaguleEstablishment,
                                                   DetectionProb=0,ManageProb=0,
                                                   MortalityProb=0,SpreadReduction=0,
                                                   SEAM=NULL,InfoRetentionProb=1,
                                                   FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.meta_detection_single_parent_operator_pre_fecundity(SDDprob,LDDprob,
      LDDrate,EnvEstabProb,Survival,K,PropaguleProduction,PropaguleEstablishment,
      DetectionProb,ManageProb,MortalityProb,SpreadReduction,SEAM,InfoRetentionProb))
  n<-nrow(as.matrix(SDDprob));d<-.ina_recycle(DetectionProb,n,"DetectionProb");ir<-.ina_recycle(InfoRetentionProb,n,"InfoRetentionProb")
  C<-if(is.null(SEAM)||length(SEAM)==1L)matrix(0,n,n)else as.matrix(SEAM);diag(C)<-0
  G0<-meta_single_parent_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                                  PropaguleProduction,PropaguleEstablishment,0,0,0,0)
  GH<-meta_single_parent_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                                  PropaguleProduction,PropaguleEstablishment,
                                  ManageProb,MortalityProb,SpreadReduction,FecundityReduction)
  G<-matrix(0,2*n,2*n);U<-seq_len(n);H<-n+seq_len(n)
  for(i in seq_len(n))for(j in seq_len(n)){
    w<-G0[j,i];if(w!=0){h<-d[j];G[U[j],U[i]]<-G[U[j],U[i]]+w*(1-h);G[H[j],U[i]]<-G[H[j],U[i]]+w*h}
    w<-GH[j,i];if(w!=0){h<-if(j==i)ir[j]+(1-ir[j])*d[j]else 1-(1-d[j])*(1-C[i,j]);G[U[j],H[i]]<-G[U[j],H[i]]+w*(1-h);G[H[j],H[i]]<-G[H[j],H[i]]+w*h}
  }
  attr(G,"note")<-"2N one-parent operator with simulator integer truncation and managed fecundity reduction"
  G
}

meta_single_parent_recruit_means <- function(SDDprob,LDDprob=0,LDDrate=0,
                                             EnvEstabProb=1,K,
                                             PropaguleProduction,
                                             PropaguleEstablishment,
                                             SpreadReduction=0,
                                             FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.meta_single_parent_recruit_means_pre_fecundity(SDDprob,LDDprob,
      LDDrate,EnvEstabProb,K,PropaguleProduction,PropaguleEstablishment,
      SpreadReduction))
  SDD<-as.matrix(SDDprob);n<-nrow(SDD);LDD<-.ina_mat(LDDprob,n,"LDDprob")
  env<-.ina_recycle(EnvEstabProb,n,"EnvEstabProb");cap<-.ina_recycle(K,n,"K")
  prod<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction");pe<-.ina_recycle(PropaguleEstablishment,n,"PropaguleEstablishment")
  g<-.ina_recycle(SpreadReduction,n,"SpreadReduction");f<-.ina_recycle(FecundityReduction,n,"FecundityReduction");.ina_fr_validate(f)
  r<-LDDrate;alpha<-pe*env;rs<-rowSums(SDD);rl<-rowSums(LDD);out<-list(matrix(0,n,n),matrix(0,n,n))
  for(i in seq_len(n)){
    ps<-if(rs[i]>0)SDD[i,]/rs[i]else rep(0,n);pl<-if(rl[i]>0)LDD[i,]/rl[i]else rep(0,n)
    for(M in 0:1){supp<-.meta_poisson_support(prod[i]*(1-f[i]*M));k<-supp$k;pk<-supp$p
      ms<-floor(k*((1-r)*rs[i]));ml<-floor(k*(r*(1-g[i]*M)*rl[i]));free<-cap;free[i]<-pmax(0,free[i]-1)
      for(j in seq_len(n)){bS<-1-ps[j]*(1-exp(-alpha[j]));bL<-1-pl[j]*(1-exp(-alpha[j]));nohaz<-sum(pk*(bS^ms)*(bL^ml));out[[M+1L]][i,j]<-free[j]*(1-nohaz)}
    }
  }
  out
}

meta_single_parent_extinction <- function(SDDprob,LDDprob=0,LDDrate=0,
                                          EnvEstabProb=1,Survival=1,K,
                                          PropaguleProduction,PropaguleEstablishment,
                                          ManageProb=0,MortalityProb=0,
                                          SpreadReduction=0,generations=100,
                                          tolerance=1e-12,FecundityReduction=0){
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.meta_single_parent_extinction_pre_fecundity(SDDprob,LDDprob,LDDrate,
      EnvEstabProb,Survival,K,PropaguleProduction,PropaguleEstablishment,
      ManageProb,MortalityProb,SpreadReduction,generations,tolerance))
  n<-nrow(as.matrix(SDDprob));s<-.ina_recycle(Survival,n,"Survival");a<-.ina_recycle(ManageProb,n,"ManageProb");m<-.ina_recycle(MortalityProb,n,"MortalityProb")
  mus<-meta_single_parent_recruit_means(SDDprob,LDDprob,LDDrate,EnvEstabProb,K,
    PropaguleProduction,PropaguleEstablishment,SpreadReduction,FecundityReduction)
  q<-rep(0,n);hist<-matrix(NA_real_,n,generations)
  for(tt in seq_len(generations)){qo<-q;qn<-numeric(n);for(i in seq_len(n)){
    q0<-s[i];q1<-s[i]*(1-m[i]);R0<-exp(sum(mus[[1]][i,]*(q-1)));R1<-exp(sum(mus[[2]][i,]*(q-1)))
    f0<-(1-q0)+q0*q[i]*R0;f1<-(1-q1)+q1*q[i]*R1;qn[i]<-(1-a[i])*f0+a[i]*f1}
    q<-pmin(1,pmax(0,qn));hist[,tt]<-q;if(max(abs(q-qo))<tolerance){if(tt<generations)hist[,(tt+1):generations]<-q;break}}
  list(extinction=q,history=hist,note="Branching PGF with simulator-faithful one-parent recruit means and managed fecundity reduction")
}

meta_export_vector <- function(SDDprob,LDDprob=0,LDDrate=0,Survival=1,
                               PropaguleProduction,ManageProb=0,
                               MortalityProb=0,SpreadReduction=0,
                               ExportSDDprob=NULL,ExportLDDprob=NULL,
                               FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.meta_export_vector_pre_fecundity(SDDprob,LDDprob,LDDrate,Survival,
      PropaguleProduction,ManageProb,MortalityProb,SpreadReduction,
      ExportSDDprob,ExportLDDprob))
  SDD<-as.matrix(SDDprob);n<-nrow(SDD);LDD<-.ina_mat(LDDprob,n,"LDDprob")
  os<-.ina_residual_export(SDD,ExportSDDprob,"SDDprob");ol<-.ina_residual_export(LDD,ExportLDDprob,"LDDprob")
  s<-.ina_recycle(Survival,n,"Survival");p<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction")
  a<-.ina_recycle(ManageProb,n,"ManageProb");m<-.ina_recycle(MortalityProb,n,"MortalityProb")
  g<-.ina_recycle(SpreadReduction,n,"SpreadReduction");f<-.ina_recycle(FecundityReduction,n,"FecundityReduction");.ina_fr_validate(f);r<-LDDrate
  p*((1-a)*s*((1-r)*os+r*ol)+a*s*(1-m)*(1-f)*((1-r)*os+r*(1-g)*ol))
}

meta_detection_export_vector <- function(SDDprob,LDDprob=0,LDDrate=0,
                                         Survival=1,PropaguleProduction,
                                         ManageProb=0,MortalityProb=0,
                                         SpreadReduction=0,ExportSDDprob=NULL,
                                         ExportLDDprob=NULL,FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.meta_detection_export_vector_pre_fecundity(SDDprob,LDDprob,LDDrate,
      Survival,PropaguleProduction,ManageProb,MortalityProb,SpreadReduction,
      ExportSDDprob,ExportLDDprob))
  e0<-meta_export_vector(SDDprob,LDDprob,LDDrate,Survival,PropaguleProduction,
                         0,0,0,ExportSDDprob,ExportLDDprob,0)
  e1<-meta_export_vector(SDDprob,LDDprob,LDDrate,Survival,PropaguleProduction,
                         ManageProb,MortalityProb,SpreadReduction,
                         ExportSDDprob,ExportLDDprob,FecundityReduction)
  c(e0,e1)
}

meta_single_parent_escape_branching <- function(SDDprob,LDDprob=0,LDDrate=0,
                                                EnvEstabProb=1,Survival=1,K,
                                                PropaguleProduction,
                                                PropaguleEstablishment,
                                                ManageProb=0,MortalityProb=0,
                                                SpreadReduction=0,
                                                ExportSDDprob=NULL,
                                                ExportLDDprob=NULL,
                                                OutsideEstablishmentProb=1,
                                                AssumeResidualExport=FALSE,
                                                timesteps=10,
                                                FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.meta_single_parent_escape_branching_pre_fecundity(SDDprob,LDDprob,
      LDDrate,EnvEstabProb,Survival,K,PropaguleProduction,
      PropaguleEstablishment,ManageProb,MortalityProb,SpreadReduction,
      ExportSDDprob,ExportLDDprob,OutsideEstablishmentProb,
      AssumeResidualExport,timesteps))
  SDD<-as.matrix(SDDprob);n<-nrow(SDD);LDD<-.ina_mat(LDDprob,n,"LDDprob")
  os<-.meta_export_fraction(SDD,ExportSDDprob,AssumeResidualExport,"SDD")
  ol<-.meta_export_fraction(LDD,ExportLDDprob,AssumeResidualExport,"LDD")
  if(is.null(os)&&is.null(ol))stop("Explicit export matrices are required unless AssumeResidualExport=TRUE")
  if(is.null(os))os<-rep(0,n);if(is.null(ol))ol<-rep(0,n)
  s<-.ina_recycle(Survival,n,"Survival");a<-.ina_recycle(ManageProb,n,"ManageProb")
  m<-.ina_recycle(MortalityProb,n,"MortalityProb");g<-.ina_recycle(SpreadReduction,n,"SpreadReduction")
  prod<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction");f<-.ina_recycle(FecundityReduction,n,"FecundityReduction")
  pout<-.ina_recycle(OutsideEstablishmentProb,n,"OutsideEstablishmentProb")
  mus<-meta_single_parent_recruit_means(SDD,LDD,LDDrate,EnvEstabProb,K,
        PropaguleProduction,PropaguleEstablishment,SpreadReduction,FecundityReduction)
  no0<-.meta_no_export_probability(prod,rowSums(SDD),os,rowSums(LDD),ol,
         LDDrate,g,pout,FALSE)
  no1<-.meta_no_export_probability(prod*(1-f),rowSums(SDD),os,rowSums(LDD),ol,
         LDDrate,g,pout,TRUE)
  h<-rep(1,n);hist<-matrix(NA_real_,n,timesteps)
  for(tt in seq_len(timesteps)){hn<-numeric(n);for(i in seq_len(n)){
    q0<-s[i];q1<-s[i]*(1-m[i]);R0<-exp(sum(mus[[1]][i,]*(h-1)));R1<-exp(sum(mus[[2]][i,]*(h-1)))
    f0<-(1-q0)+q0*h[i]*R0*no0[i];f1<-(1-q1)+q1*h[i]*R1*no1[i]
    hn[i]<-(1-a[i])*f0+a[i]*f1};h<-pmin(1,pmax(0,hn));hist[,tt]<-h}
  list(no_escape=h,escape=1-h,history_escape=1-hist,
       note="Branching PGF with simulator-faithful recruitment/export and managed fecundity reduction")
}

###############################################################################
### INApestMetaTransitionMatrix
###############################################################################
transition_components <- function(Transition,Nstages,SDDprob,LDDprob=0,
                                  LDDrate=0,EnvEstabProb=1,
                                  PropaguleEstablishment=1,ManageProb=0,
                                  MortalityProb=0,SpreadReduction=0,
                                  DispersalDensityFactor=0,K=1,SeedbankK=1,
                                  FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.transition_components_pre_fecundity(Transition,Nstages,SDDprob,
      LDDprob,LDDrate,EnvEstabProb,PropaguleEstablishment,ManageProb,
      MortalityProb,SpreadReduction,DispersalDensityFactor,K,SeedbankK))
  z <- .transition_components_pre_fecundity(Transition,Nstages,SDDprob,
      LDDprob,LDDrate,EnvEstabProb,PropaguleEstablishment,ManageProb,
      MortalityProb,SpreadReduction,DispersalDensityFactor,K,SeedbankK)
  z$fecundity_reduction <- .ina_transition_fecundity_matrix(
    FecundityReduction,z$n,z$S,"FecundityReduction")
  z
}

transition_operator <- function(Transition,Nstages,SDDprob,LDDprob=0,LDDrate=0,
                                EnvEstabProb=1,PropaguleEstablishment=1,
                                ManageProb=0,MortalityProb=0,SpreadReduction=0,
                                DispersalDensityFactor=0,K=1,SeedbankK=1,
                                FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.transition_operator_pre_fecundity(Transition,Nstages,SDDprob,
      LDDprob,LDDrate,EnvEstabProb,PropaguleEstablishment,ManageProb,
      MortalityProb,SpreadReduction,DispersalDensityFactor,K,SeedbankK))
  z<-transition_components(Transition,Nstages,SDDprob,LDDprob,LDDrate,
       EnvEstabProb,PropaguleEstablishment,ManageProb,MortalityProb,
       SpreadReduction,DispersalDensityFactor,K,SeedbankK,FecundityReduction)
  n<-z$n;S<-z$S;G<-matrix(0,n*S,n*S);idx<-function(i,s)(i-1L)*S+s
  for(i in seq_len(n)){
    Ai<-z$A[[i]];ai<-z$adoption[i];gi<-z$spread_reduction[i]
    for(k in seq_len(S)){
      qbar<-1-ai*z$mortality[i,k];src<-idx(i,k)
      if(k<S){G[idx(i,k),src]<-G[idx(i,k),src]+qbar*Ai[k,k];G[idx(i,k+1L),src]<-G[idx(i,k+1L),src]+qbar*Ai[k+1L,k]}
      else G[idx(i,S),src]<-G[idx(i,S),src]+qbar*Ai[S,S]
      if(k>=2L&&Ai[1,k]>0){fec<-Ai[1,k];q0<-1;q1<-1-z$mortality[i,k];fr<-z$fecundity_reduction[i,k]
        for(j in seq_len(n)){
          nat<-(1-z$LDDrate)*z$sdd_enabled[i]*z$SDD[i,j]
          hum0<-z$LDDrate*z$LDD[i,j];hum1<-z$LDDrate*(1-gi)*z$LDD[i,j]
          w<-fec*z$recruit_success[j]*((1-ai)*q0*(nat+hum0)+ai*q1*(1-fr)*(nat+hum1))
          G[idx(j,1L),src]<-G[idx(j,1L),src]+w
        }
      }
    }
  }
  attr(G,"components")<-z;G
}

transition_intrinsic_managed_lambda <- function(Transition,MortalityProb=0,
                                                ManageProb=1,
                                                FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.transition_intrinsic_managed_lambda_pre_fecundity(
      Transition,MortalityProb,ManageProb))
  A<-as.matrix(Transition);S<-nrow(A)
  m<-if(length(MortalityProb)==1)rep(MortalityProb,S)else as.numeric(MortalityProb)
  f<-if(length(FecundityReduction)==1)rep(FecundityReduction,S)else as.numeric(FecundityReduction)
  if(length(m)!=S||length(f)!=S)stop("MortalityProb and FecundityReduction must be scalar or one value per stage")
  .ina_fr_validate(f)
  survbar<-(1-ManageProb)+ManageProb*(1-m)
  M<-matrix(0,S,S)
  for(k in seq_len(S)){
    if(k<S){M[k,k]<-A[k,k]*survbar[k];M[k+1L,k]<-A[k+1L,k]*survbar[k]}
    else M[S,S]<-A[S,S]*survbar[S]
    if(k>=2L)M[1L,k]<-A[1L,k]*((1-ManageProb)+ManageProb*(1-m[k])*(1-f[k]))
  }
  list(ManagedTransition=M,Lambda=max(Mod(eigen(M,only.values=TRUE)$values)))
}

transition_extinction <- function(Transition,Nstages,SDDprob,LDDprob=0,
                                  LDDrate=0,EnvEstabProb=1,
                                  PropaguleEstablishment=1,ManageProb=0,
                                  MortalityProb=0,SpreadReduction=0,
                                  DispersalDensityFactor=0,K=1,SeedbankK=1,
                                  generations=100,FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.transition_extinction_pre_fecundity(Transition,Nstages,SDDprob,
      LDDprob,LDDrate,EnvEstabProb,PropaguleEstablishment,ManageProb,
      MortalityProb,SpreadReduction,DispersalDensityFactor,K,SeedbankK,
      generations))
  z<-transition_components(Transition,Nstages,SDDprob,LDDprob,LDDrate,
       EnvEstabProb,PropaguleEstablishment,ManageProb,MortalityProb,
       SpreadReduction,DispersalDensityFactor,K,SeedbankK,FecundityReduction)
  n<-z$n;S<-z$S;nt<-n*S;idx<-function(i,s)(i-1L)*S+s
  q<-rep(0,nt);hist<-matrix(NA_real_,nt,generations)
  for(tt in seq_len(generations)){qn<-numeric(nt);for(i in seq_len(n)){
    Ai<-z$A[[i]];ai<-z$adoption[i];gi<-z$spread_reduction[i]
    for(k in seq_len(S)){src<-idx(i,k);fec<-if(k>=2)Ai[1,k]else 0
      calcM<-function(M){surv<-1-z$mortality[i,k]*M
        if(k<S)local<-(1-Ai[k,k]-Ai[k+1,k])+Ai[k,k]*q[idx(i,k)]+Ai[k+1,k]*q[idx(i,k+1L)]
        else local<-(1-Ai[S,S])+Ai[S,S]*q[idx(i,S)]
        muq<-0
        if(fec>0)for(j in seq_len(n)){
          nat<-(1-z$LDDrate)*z$sdd_enabled[i]*z$SDD[i,j]
          hum<-z$LDDrate*(1-gi*M)*z$LDD[i,j]
          mu<-fec*(1-z$fecundity_reduction[i,k]*M)*(nat+hum)*z$recruit_success[j]
          muq<-muq+mu*(q[idx(j,1L)]-1)}
        (1-surv)+surv*local*exp(muq)}
      qn[src]<-(1-ai)*calcM(0)+ai*calcM(1)}}
    q<-pmin(1,pmax(0,qn));hist[,tt]<-q}
  list(extinction=q,history=hist,
       note="Stage x node branching PGF with management mortality and fecundity reduction represented separately")
}

transition_detection_operator <- function(Transition,Nstages,SDDprob,LDDprob=0,
                                           LDDrate=0,EnvEstabProb=1,
                                           PropaguleEstablishment=1,
                                           DetectionProb=0,ManageProb=0,
                                           MortalityProb=0,SpreadReduction=0,
                                           SEAM=NULL,InfoRetentionProb=1,
                                           DispersalDensityFactor=0,K=1,
                                           SeedbankK=1,FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.transition_detection_operator_pre_fecundity(Transition,Nstages,
      SDDprob,LDDprob,LDDrate,EnvEstabProb,PropaguleEstablishment,
      DetectionProb,ManageProb,MortalityProb,SpreadReduction,SEAM,
      InfoRetentionProb,DispersalDensityFactor,K,SeedbankK))
  base<-transition_operator(Transition,Nstages,SDDprob,LDDprob,LDDrate,
        EnvEstabProb,PropaguleEstablishment,0,0,0,DispersalDensityFactor,K,SeedbankK,0)
  managed<-transition_operator(Transition,Nstages,SDDprob,LDDprob,LDDrate,
        EnvEstabProb,PropaguleEstablishment,ManageProb,MortalityProb,
        SpreadReduction,DispersalDensityFactor,K,SeedbankK,FecundityReduction)
  z<-attr(managed,"components");n<-z$n;S<-z$S;nt<-n*S
  if(length(DetectionProb)==1)D<-matrix(DetectionProb,n,S)
  else if(length(DetectionProb)==S)D<-matrix(rep(DetectionProb,each=n),n,S)
  else if(is.matrix(DetectionProb)&&all(dim(DetectionProb)==c(n,S)))D<-DetectionProb
  else stop("DetectionProb must be scalar, length Nstages, or nodes x Nstages")
  ir<-.ina_recycle(InfoRetentionProb,n,"InfoRetentionProb")
  C<-if(is.null(SEAM)||length(SEAM)==1L)matrix(0,n,n)else as.matrix(SEAM);diag(C)<-0
  idx<-function(i,s)(i-1L)*S+s;G<-matrix(0,2*nt,2*nt);U<-seq_len(nt);H<-nt+seq_len(nt)
  for(i in seq_len(n))for(k in seq_len(S)){src<-idx(i,k);for(j in seq_len(n))for(s in seq_len(S)){
    dst<-idx(j,s);w<-base[dst,src];if(w!=0){h<-D[j,s];G[U[dst],U[src]]<-G[U[dst],U[src]]+w*(1-h);G[H[dst],U[src]]<-G[H[dst],U[src]]+w*h}
    w<-managed[dst,src];if(w!=0){h<-if(j==i)ir[j]+(1-ir[j])*D[j,s]else 1-(1-D[j,s])*(1-C[i,j]);G[U[dst],H[src]]<-G[U[dst],H[src]]+w*(1-h);G[H[dst],H[src]]<-G[H[dst],H[src]]+w*h}}}
  attr(G,"note")<-"2NS low-density approximation with managed fecundity reduction"
  G
}

transition_export_vector <- function(Transition,Nstages,SDDprob,LDDprob=0,
                                     LDDrate=0,ManageProb=0,MortalityProb=0,
                                     SpreadReduction=0,DispersalDensityFactor=0,
                                     ExportSDDprob=NULL,ExportLDDprob=NULL,
                                     FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.transition_export_vector_pre_fecundity(Transition,Nstages,SDDprob,
      LDDprob,LDDrate,ManageProb,MortalityProb,SpreadReduction,
      DispersalDensityFactor,ExportSDDprob,ExportLDDprob))
  SDD<-transition_zero_density_sdd(SDDprob,DispersalDensityFactor);n<-nrow(SDD);S<-Nstages
  LDD<-.ina_mat(LDDprob,n,"LDDprob");A<-.transition_list(Transition,n,S)
  a<-.ina_recycle(ManageProb,n,"ManageProb");g<-.ina_recycle(SpreadReduction,n,"SpreadReduction")
  if(length(MortalityProb)==1L)M<-matrix(MortalityProb,n,S)
  else if(length(MortalityProb)==S)M<-matrix(rep(MortalityProb,each=n),n,S)
  else if(is.matrix(MortalityProb)&&all(dim(MortalityProb)==c(n,S)))M<-MortalityProb
  else stop("MortalityProb must be scalar, length Nstages, or nodes x Nstages")
  F<-.ina_transition_fecundity_matrix(FecundityReduction,n,S)
  os<-.ina_residual_export(SDD,ExportSDDprob,"SDDprob");ol<-.ina_residual_export(LDD,ExportLDDprob,"LDDprob")
  out<-numeric(n*S);idx<-function(i,k)(i-1L)*S+k;r<-LDDrate
  for(i in seq_len(n))for(k in 2:S){fec<-A[[i]][1,k];if(fec<=0)next
    out[idx(i,k)]<-fec*((1-a[i])*((1-r)*os[i]+r*ol[i])+
      a[i]*(1-M[i,k])*(1-F[i,k])*((1-r)*os[i]+r*(1-g[i])*ol[i]))}
  out
}

transition_detection_export_vector <- function(Transition,Nstages,SDDprob,
                                               LDDprob=0,LDDrate=0,ManageProb=0,
                                               MortalityProb=0,SpreadReduction=0,
                                               DispersalDensityFactor=0,
                                               ExportSDDprob=NULL,ExportLDDprob=NULL,
                                               FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.transition_detection_export_vector_pre_fecundity(Transition,Nstages,
      SDDprob,LDDprob,LDDrate,ManageProb,MortalityProb,SpreadReduction,
      DispersalDensityFactor,ExportSDDprob,ExportLDDprob))
  e0<-transition_export_vector(Transition,Nstages,SDDprob,LDDprob,LDDrate,0,0,0,
       DispersalDensityFactor,ExportSDDprob,ExportLDDprob,0)
  e1<-transition_export_vector(Transition,Nstages,SDDprob,LDDprob,LDDrate,
       ManageProb,MortalityProb,SpreadReduction,DispersalDensityFactor,
       ExportSDDprob,ExportLDDprob,FecundityReduction)
  c(e0,e1)
}

transition_escape_branching <- function(Transition,Nstages,SDDprob,LDDprob=0,
                                        LDDrate=0,EnvEstabProb=1,
                                        PropaguleEstablishment=1,ManageProb=0,
                                        MortalityProb=0,SpreadReduction=0,
                                        DispersalDensityFactor=0,K=1,SeedbankK=1,
                                        ExportSDDprob=NULL,ExportLDDprob=NULL,
                                        OutsideEstablishmentProb=1,
                                        AssumeResidualExport=TRUE,timesteps=10,
                                        FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.transition_escape_branching_pre_fecundity(Transition,Nstages,SDDprob,
      LDDprob,LDDrate,EnvEstabProb,PropaguleEstablishment,ManageProb,
      MortalityProb,SpreadReduction,DispersalDensityFactor,K,SeedbankK,
      ExportSDDprob,ExportLDDprob,OutsideEstablishmentProb,
      AssumeResidualExport,timesteps))
  z<-transition_components(Transition,Nstages,SDDprob,LDDprob,LDDrate,
       EnvEstabProb,PropaguleEstablishment,ManageProb,MortalityProb,
       SpreadReduction,DispersalDensityFactor,K,SeedbankK,FecundityReduction)
  n<-z$n;S<-z$S;nt<-n*S;idx<-function(i,s)(i-1L)*S+s
  os<-.meta_export_fraction(z$SDD,ExportSDDprob,AssumeResidualExport,"SDD")
  ol<-.meta_export_fraction(z$LDD,ExportLDDprob,AssumeResidualExport,"LDD")
  if(is.null(os))os<-rep(0,n);if(is.null(ol))ol<-rep(0,n)
  pout<-.ina_recycle(OutsideEstablishmentProb,n,"OutsideEstablishmentProb")
  h<-rep(1,nt);hist<-matrix(NA_real_,nt,timesteps)
  for(tt in seq_len(timesteps)){hn<-numeric(nt);for(i in seq_len(n)){
    Ai<-z$A[[i]];ai<-z$adoption[i];gi<-z$spread_reduction[i]
    for(k in seq_len(S)){src<-idx(i,k);fec<-if(k>=2L)Ai[1,k]else 0
      calcM<-function(M){surv<-1-z$mortality[i,k]*M
        if(k<S)local<-(1-Ai[k,k]-Ai[k+1L,k])+Ai[k,k]*h[idx(i,k)]+Ai[k+1L,k]*h[idx(i,k+1L)]
        else local<-(1-Ai[S,S])+Ai[S,S]*h[idx(i,S)]
        feff<-fec*(1-z$fecundity_reduction[i,k]*M);muq<-0
        if(feff>0)for(j in seq_len(n)){
          nat<-(1-z$LDDrate)*z$sdd_enabled[i]*z$SDD[i,j];hum<-z$LDDrate*(1-gi*M)*z$LDD[i,j]
          mu<-feff*(nat+hum)*z$recruit_success[j];muq<-muq+mu*(h[idx(j,1L)]-1)}
        outprob<-(1-z$LDDrate)*z$sdd_enabled[i]*os[i]+z$LDDrate*(1-gi*M)*ol[i]
        noexp<-exp(-feff*outprob*pout[i]);(1-surv)+surv*local*exp(muq)*noexp}
      hn[src]<-(1-ai)*calcM(0)+ai*calcM(1)}}
    h<-pmin(1,pmax(0,hn));hist[,tt]<-h}
  list(no_escape=h,escape=1-h,history_escape=1-hist,
       note="Stage x node branching no-escape PGF with managed fecundity reduction")
}

###############################################################################
### INApestMetaMultipleLandUse
###############################################################################
mlu_components <- function(SDDprob,LDDprob=0,LDDrate=0,EnvEstabProb=1,
                           Survival=1,K,PropaguleProduction,
                           PropaguleEstablishment,ManageProb=0,
                           MortalityProb=0,SpreadReduction=0,
                           current_code=TRUE,FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.mlu_components_pre_fecundity(SDDprob,LDDprob,LDDrate,EnvEstabProb,
      Survival,K,PropaguleProduction,PropaguleEstablishment,ManageProb,
      MortalityProb,SpreadReduction,current_code))
  SDD<-as.matrix(SDDprob);n<-nrow(SDD);LDD<-.ina_mat(LDDprob,n,"LDDprob")
  K<-as.matrix(K);if(nrow(K)!=n)stop("K must have one row per node");L<-ncol(K)
  s<-.ina_recycle(Survival,n,"Survival");p<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction")
  env<-.ina_recycle(EnvEstabProb,n,"EnvEstabProb");pe<-.ina_recycle(PropaguleEstablishment,n,"PropaguleEstablishment")
  A<-.mlu_matrix(ManageProb,n,L,"ManageProb");M<-.mlu_matrix(MortalityProb,n,L,"MortalityProb")
  Gm<-.mlu_matrix(SpreadReduction,n,L,"SpreadReduction");F<-.mlu_matrix(FecundityReduction,n,L,"FecundityReduction");.ina_fr_validate(F)
  qbar<-matrix(s,n,L)*(1-A*M);r<-LDDrate;alpha<-env*pe;c_est<-1-exp(-alpha)
  Em<-array(0,c(n,L,n))
  for(i in seq_len(n))for(l in seq_len(L)){
    natcoef<-s[i]*((1-A[i,l])+A[i,l]*(1-M[i,l])*(1-F[i,l]))*(1-r)
    if(current_code){
      # Legacy current_code branch retained for historical diagnostics.  The
      # production multiplier is source-LU specific; the old unweighted LDD
      # management sum is preserved only when explicitly requested.
      lddcoef<-s[i]*((1-A[i,l])+A[i,l]*(1-M[i,l])*(1-F[i,l]))
      sf<-sum(1-Gm[i,]*A[i,]); lddcoef<-lddcoef*sf
    } else {
      lddcoef<-s[i]*((1-A[i,l])+A[i,l]*(1-M[i,l])*(1-F[i,l])*(1-Gm[i,l]))
    }
    Em[i,l,]<-p[i]*(natcoef*SDD[i,]+r*lddcoef*LDD[i,])
  }
  list(n=n,L=L,SDD=SDD,LDD=LDD,K=K,Ktot=rowSums(K),survival=s,production=p,
       env=env,prop_est=pe,adoption=A,mortality=M,spread_reduction=Gm,
       fecundity_reduction=F,qbar=qbar,c_est=c_est,emission=Em,
       current_code=current_code)
}

mlu_operator <- function(SDDprob,LDDprob=0,LDDrate=0,EnvEstabProb=1,
                         Survival=1,K,PropaguleProduction,PropaguleEstablishment,
                         ManageProb=0,MortalityProb=0,SpreadReduction=0,
                         current_code=TRUE,FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.mlu_operator_pre_fecundity(SDDprob,LDDprob,LDDrate,EnvEstabProb,
      Survival,K,PropaguleProduction,PropaguleEstablishment,ManageProb,
      MortalityProb,SpreadReduction,current_code))
  z<-mlu_components(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                    PropaguleProduction,PropaguleEstablishment,ManageProb,
                    MortalityProb,SpreadReduction,current_code,FecundityReduction)
  n<-z$n;L<-z$L;idx<-function(i,l)(i-1L)*L+l;G<-matrix(0,n*L,n*L)
  for(i in seq_len(n))for(l in seq_len(L)){
    src<-idx(i,l);G[src,src]<-G[src,src]+z$qbar[i,l]
    for(j in seq_len(n))for(h in seq_len(L))G[idx(j,h),src]<-G[idx(j,h),src]+z$K[j,h]*z$c_est[j]*z$emission[i,l,j]}
  attr(G,"components")<-z;G
}

mlu_meanfield <- function(SDDprob,LDDprob=0,LDDrate=0,EnvEstabProb=1,
                          Survival=1,K,PropaguleProduction,PropaguleEstablishment,
                          ManageProb=0,MortalityProb=0,SpreadReduction=0,
                          current_code=TRUE,initial,timesteps,
                          FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.mlu_meanfield_pre_fecundity(SDDprob,LDDprob,LDDrate,EnvEstabProb,
      Survival,K,PropaguleProduction,PropaguleEstablishment,ManageProb,
      MortalityProb,SpreadReduction,current_code,initial,timesteps))
  z<-mlu_components(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                    PropaguleProduction,PropaguleEstablishment,ManageProb,
                    MortalityProb,SpreadReduction,current_code,FecundityReduction)
  X<-as.matrix(initial);if(!all(dim(X)==c(z$n,z$L)))stop("initial must be nodes x landuses")
  out<-array(NA_real_,c(z$n,z$L,timesteps))
  for(tt in seq_len(timesteps)){N0<-z$qbar*X;lam<-numeric(z$n)
    for(i in seq_len(z$n))for(l in seq_len(z$L))lam<-lam+X[i,l]*z$emission[i,l,]
    p_rec<-1-exp(-z$c_est*lam);free<-pmax(0,z$K-N0);free_tot<-rowSums(free)
    total_rec<-free_tot*p_rec;share<-free/free_tot;share[!is.finite(share)]<-0
    X<-N0+share*total_rec;X<-pmin(z$K,pmax(0,X));out[,,tt]<-X}
  out
}

mlu_detection_operator <- function(SDDprob,LDDprob=0,LDDrate=0,
                                   EnvEstabProb=1,Survival=1,K,
                                   PropaguleProduction,PropaguleEstablishment,
                                   DetectionProb=0,ManageProb=0,MortalityProb=0,
                                   SpreadReduction=0,SEAM=NULL,
                                   InfoRetentionProb=1,current_code=TRUE,
                                   FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.mlu_detection_operator_pre_fecundity(SDDprob,LDDprob,LDDrate,
      EnvEstabProb,Survival,K,PropaguleProduction,PropaguleEstablishment,
      DetectionProb,ManageProb,MortalityProb,SpreadReduction,SEAM,
      InfoRetentionProb,current_code))
  Kmat<-as.matrix(K);n<-nrow(Kmat);L<-ncol(Kmat);nt<-n*L;D<-.mlu_matrix(DetectionProb,n,L,"DetectionProb")
  ir<-.ina_recycle(InfoRetentionProb,n,"InfoRetentionProb");C<-if(is.null(SEAM)||length(SEAM)==1L)matrix(0,n,n)else as.matrix(SEAM);diag(C)<-0
  G0<-mlu_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                   PropaguleProduction,PropaguleEstablishment,0,0,0,current_code,0)
  GH<-mlu_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
                   PropaguleProduction,PropaguleEstablishment,ManageProb,
                   MortalityProb,SpreadReduction,current_code,FecundityReduction)
  idx<-function(i,l)(i-1L)*L+l;G<-matrix(0,2*nt,2*nt);U<-seq_len(nt);H<-nt+seq_len(nt)
  for(i in seq_len(n))for(l in seq_len(L)){src<-idx(i,l);for(j in seq_len(n))for(h in seq_len(L)){
    dst<-idx(j,h);w<-G0[dst,src];if(w!=0){hh<-D[j,h];G[U[dst],U[src]]<-G[U[dst],U[src]]+w*(1-hh);G[H[dst],U[src]]<-G[H[dst],U[src]]+w*hh}
    w<-GH[dst,src];if(w!=0){hh<-if(j==i)ir[j]+(1-ir[j])*D[j,h]else 1-(1-D[j,h])*(1-C[i,j]);G[U[dst],H[src]]<-G[U[dst],H[src]]+w*(1-hh);G[H[dst],H[src]]<-G[H[dst],H[src]]+w*hh}}}
  attr(G,"note")<-"2NL low-density approximation with land-use fecundity reduction";G
}

mlu_single_parent_operator <- function(SDDprob,LDDprob=0,LDDrate=0,
                                       EnvEstabProb=1,Survival=1,K,
                                       PropaguleProduction,PropaguleEstablishment,
                                       ManageProb=0,MortalityProb=0,
                                       SpreadReduction=0,current_code=FALSE,
                                       FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.mlu_single_parent_operator_pre_fecundity(SDDprob,LDDprob,LDDrate,
      EnvEstabProb,Survival,K,PropaguleProduction,PropaguleEstablishment,
      ManageProb,MortalityProb,SpreadReduction,current_code))
  SDD<-as.matrix(SDDprob);n<-nrow(SDD);LDD<-.ina_mat(LDDprob,n,"LDDprob");K<-as.matrix(K);L<-ncol(K)
  s<-.ina_recycle(Survival,n,"Survival");prod<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction")
  env<-.ina_recycle(EnvEstabProb,n,"EnvEstabProb");pe<-.ina_recycle(PropaguleEstablishment,n,"PropaguleEstablishment")
  A<-.mlu_matrix(ManageProb,n,L,"ManageProb");Mort<-.mlu_matrix(MortalityProb,n,L,"MortalityProb");Gr<-.mlu_matrix(SpreadReduction,n,L,"SpreadReduction");F<-.mlu_matrix(FecundityReduction,n,L,"FecundityReduction");.ina_fr_validate(F)
  alpha<-env*pe;r<-LDDrate;rs<-rowSums(SDD);rl<-rowSums(LDD);idx<-function(i,l)(i-1L)*L+l;G<-matrix(0,n*L,n*L)
  for(i in seq_len(n)){
    ps<-if(rs[i]>0)SDD[i,]/rs[i]else rep(0,n);pl<-if(rl[i]>0)LDD[i,]/rl[i]else rep(0,n);st<-.mlu_management_states(A[i,])
    for(l in seq_len(L)){src<-idx(i,l)
      for(z in seq_len(nrow(st$M))){Mv<-st$M[z,];pm<-st$p[z];if(pm==0)next
        surv<-s[i]*(1-Mort[i,l]*Mv[l]);if(surv==0)next
        lambda<-prod[i]*(1-F[i,l]*Mv[l]);supp<-.meta_poisson_support(lambda);k<-supp$k;pk<-supp$p
        # Latest population-share-weighted implementation: a one-parent lineage
        # is entirely in land use l, so only that source land use controls LDD.
        sf<-if(current_code)sum(1-Gr[i,]*Mv)else(1-Gr[i,l]*Mv[l])
        ms<-floor(k*((1-r)*rs[i]));ml<-floor(k*(r*sf*rl[i]))
        for(j in seq_len(n)){bS<-1-ps[j]*(1-exp(-alpha[j]));bL<-1-pl[j]*(1-exp(-alpha[j]));nohaz<-sum(pk*(bS^ms)*(bL^ml));free<-K[j,];if(j==i)free[l]<-pmax(0,free[l]-1)
          for(h in seq_len(L))G[idx(j,h),src]<-G[idx(j,h),src]+pm*surv*free[h]*(1-nohaz)}
        G[src,src]<-G[src,src]+pm*surv}}
  }
  attr(G,"note")<-if(current_code)"Legacy current_code one-parent MLU operator with fecundity reduction" else "One-parent MLU operator matching fecundity-adjusted population-share LDD weighting"
  G
}

mlu_detection_single_parent_operator <- function(SDDprob,LDDprob=0,LDDrate=0,
                                                  EnvEstabProb=1,Survival=1,K,
                                                  PropaguleProduction,
                                                  PropaguleEstablishment,
                                                  DetectionProb=0,ManageProb=0,
                                                  MortalityProb=0,
                                                  SpreadReduction=0,SEAM=NULL,
                                                  InfoRetentionProb=1,
                                                  current_code=FALSE,
                                                  FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.mlu_detection_single_parent_operator_pre_fecundity(SDDprob,LDDprob,
      LDDrate,EnvEstabProb,Survival,K,PropaguleProduction,
      PropaguleEstablishment,DetectionProb,ManageProb,MortalityProb,
      SpreadReduction,SEAM,InfoRetentionProb,current_code))
  Kmat<-as.matrix(K);n<-nrow(Kmat);L<-ncol(Kmat);nt<-n*L;D<-.mlu_matrix(DetectionProb,n,L,"DetectionProb")
  ir<-.ina_recycle(InfoRetentionProb,n,"InfoRetentionProb");C<-if(is.null(SEAM)||length(SEAM)==1L)matrix(0,n,n)else as.matrix(SEAM);diag(C)<-0
  G0<-mlu_single_parent_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
      PropaguleProduction,PropaguleEstablishment,0,0,0,current_code,0)
  GH<-mlu_single_parent_operator(SDDprob,LDDprob,LDDrate,EnvEstabProb,Survival,K,
      PropaguleProduction,PropaguleEstablishment,ManageProb,MortalityProb,
      SpreadReduction,current_code,FecundityReduction)
  idx<-function(i,l)(i-1L)*L+l;G<-matrix(0,2*nt,2*nt);U<-seq_len(nt);H<-nt+seq_len(nt)
  for(i in seq_len(n))for(l in seq_len(L)){src<-idx(i,l);for(j in seq_len(n))for(h in seq_len(L)){dst<-idx(j,h)
    w<-G0[dst,src];if(w!=0){hh<-D[j,h];G[U[dst],U[src]]<-G[U[dst],U[src]]+w*(1-hh);G[H[dst],U[src]]<-G[H[dst],U[src]]+w*hh}
    w<-GH[dst,src];if(w!=0){hh<-if(j==i)ir[j]+(1-ir[j])*D[j,h]else 1-(1-D[j,h])*(1-C[i,j]);G[U[dst],H[src]]<-G[U[dst],H[src]]+w*(1-hh);G[H[dst],H[src]]<-G[H[dst],H[src]]+w*hh}}}
  G
}

mlu_single_parent_extinction <- function(...,generations=100,current_code=FALSE,
                                         FecundityReduction=0){
  dots<-list(...)
  if (!.ina_fr_nonzero(FecundityReduction))
    return(do.call(.mlu_single_parent_extinction_pre_fecundity,
                   c(dots,list(generations=generations,current_code=current_code))))
  G<-do.call(mlu_single_parent_operator,
             c(dots,list(current_code=current_code,
                         FecundityReduction=FecundityReduction)))
  nt<-nrow(G);self<-diag(G);rec<-G;diag(rec)<-0
  for(i in seq_len(nt))rec[i,i]<-pmax(0,G[i,i]-pmin(1,self[i]))
  qsurv<-pmin(1,diag(G));q<-rep(0,nt);hist<-matrix(NA_real_,nt,generations)
  for(tt in seq_len(generations)){qn<-numeric(nt);for(i in seq_len(nt)){
    mu<-if(qsurv[i]>0)G[,i]/qsurv[i]else rep(0,nt);mu[i]<-pmax(0,mu[i]-1)
    qn[i]<-(1-qsurv[i])+qsurv[i]*q[i]*exp(sum(mu*(q-1)))};q<-pmin(1,pmax(0,qn));hist[,tt]<-q}
  list(extinction=q,history=hist,note="Mean-matched branching approximation from the fecundity-aware MLU one-parent operator")
}

mlu_export_vector <- function(SDDprob,LDDprob=0,LDDrate=0,Survival=1,K,
                              PropaguleProduction,ManageProb=0,MortalityProb=0,
                              SpreadReduction=0,current_code=FALSE,
                              ExportSDDprob=NULL,ExportLDDprob=NULL,
                              FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.mlu_export_vector_pre_fecundity(SDDprob,LDDprob,LDDrate,Survival,K,
      PropaguleProduction,ManageProb,MortalityProb,SpreadReduction,current_code,
      ExportSDDprob,ExportLDDprob))
  SDD<-as.matrix(SDDprob);n<-nrow(SDD);LDD<-.ina_mat(LDDprob,n,"LDDprob");K<-as.matrix(K);L<-ncol(K)
  s<-.ina_recycle(Survival,n,"Survival");prod<-.ina_recycle(PropaguleProduction,n,"PropaguleProduction")
  A<-.mlu_matrix(ManageProb,n,L,"ManageProb");Mort<-.mlu_matrix(MortalityProb,n,L,"MortalityProb");Gr<-.mlu_matrix(SpreadReduction,n,L,"SpreadReduction");F<-.mlu_matrix(FecundityReduction,n,L,"FecundityReduction");.ina_fr_validate(F)
  os<-.ina_residual_export(SDD,ExportSDDprob,"SDDprob");ol<-.ina_residual_export(LDD,ExportLDDprob,"LDDprob")
  idx<-function(i,l)(i-1L)*L+l;out<-numeric(n*L);r<-LDDrate
  for(i in seq_len(n)){st<-.mlu_management_states(A[i,]);for(l in seq_len(L)){zsum<-0
    for(z in seq_len(nrow(st$M))){Mv<-st$M[z,];pm<-st$p[z];surv<-s[i]*(1-Mort[i,l]*Mv[l]);lambda<-prod[i]*(1-F[i,l]*Mv[l]);sf<-if(current_code)sum(1-Gr[i,]*Mv)else(1-Gr[i,l]*Mv[l])
      zsum<-zsum+pm*surv*lambda*((1-r)*os[i]+r*sf*ol[i])};out[idx(i,l)]<-zsum}}
  out
}

mlu_detection_export_vector <- function(SDDprob,LDDprob=0,LDDrate=0,Survival=1,
                                        K,PropaguleProduction,ManageProb=0,
                                        MortalityProb=0,SpreadReduction=0,
                                        current_code=FALSE,ExportSDDprob=NULL,
                                        ExportLDDprob=NULL,FecundityReduction=0) {
  if (!.ina_fr_nonzero(FecundityReduction))
    return(.mlu_detection_export_vector_pre_fecundity(SDDprob,LDDprob,LDDrate,
      Survival,K,PropaguleProduction,ManageProb,MortalityProb,SpreadReduction,
      current_code,ExportSDDprob,ExportLDDprob))
  e0<-mlu_export_vector(SDDprob,LDDprob,LDDrate,Survival,K,PropaguleProduction,
        0,0,0,current_code,ExportSDDprob,ExportLDDprob,0)
  e1<-mlu_export_vector(SDDprob,LDDprob,LDDrate,Survival,K,PropaguleProduction,
        ManageProb,MortalityProb,SpreadReduction,current_code,ExportSDDprob,
        ExportLDDprob,FecundityReduction)
  c(e0,e1)
}

###############################################################################
### INApest analytical extension: point-based model families
###
### Integrated point-model analytical extension for the unified INApestAnalytical.R.
### Point-grid/kernel helper functions are supplied by INApestMetaPoint.R and,
### for the stage-structured point model, INApestPointTransitionMatrix.R.
###
### Design principles
###  * preserve the existing rare-lineage/operator/branching interpretation;
###  * use exact count-level branching where continuous-space geometry cancels;
###  * when spatial heterogeneity/boundaries matter, contract continuous kernels
###    onto an INApestSpatialGrid analysis grid by reproducible Monte Carlo;
###  * never silently turn LocalK/KRadius into ordinary node carrying capacity.
###############################################################################

.ina_pt_clip01 <- function(x) pmin(1, pmax(0, x))

.ina_pt_require_point_helpers <- function(transition = FALSE) {
  required <- c(".ipp_resolve", ".ipp_spatial_value", ".ipp_habitat_search",
                ".ipp_draw_displacement", ".ipp_validate_probability_schedule")
  if (transition) required <- c(required, ".ipptm_get_transition",
                                ".ipptm_resolve_stage_schedule",
                                ".ipptm_get_kernel",
                                ".ipptm_get_transition_kernel")
  missing <- required[!vapply(required, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(missing)) {
    stop("Point analytical support requires the point simulation source to be sourced first. Missing helper(s): ",
         paste(missing, collapse = ", "), ".")
  }
  invisible(TRUE)
}

.ina_pt_is_native_grid <- function(x) inherits(x, "INApestSpatialGrid")
.ina_pt_is_schedule <- function(x) inherits(x, "INApestSpatialSchedule")

.ina_pt_first_native_grid <- function(x) {
  if (is.null(x)) return(NULL)
  if (.ina_pt_is_native_grid(x)) return(x)
  if (.ina_pt_is_schedule(x)) {
    if (!length(x$grids)) return(NULL)
    for (g in x$grids) if (.ina_pt_is_native_grid(g)) return(g)
  }
  NULL
}

.ina_pt_analysis_grid <- function(PointAnalysisGrid = NULL,
                                  HabitatSuitability = NULL,
                                  DetectionSpatial = NULL,
                                  ManageSpatial = NULL,
                                  MortalitySpatial = NULL,
                                  FecundityReductionSpatial = NULL,
                                  SpreadReductionSpatial = NULL) {
  if (!is.null(PointAnalysisGrid)) {
    if (!.ina_pt_is_native_grid(PointAnalysisGrid))
      stop("PointAnalysisGrid must be an INApestSpatialGrid.")
    if (PointAnalysisGrid$nlayer != 1L)
      stop("PointAnalysisGrid is geometry only and must have one layer.")
    return(PointAnalysisGrid)
  }
  candidates <- list(HabitatSuitability, DetectionSpatial, ManageSpatial,
                     MortalitySpatial, FecundityReductionSpatial,
                     SpreadReductionSpatial)
  for (z in candidates) {
    g <- .ina_pt_first_native_grid(z)
    if (!is.null(g)) {
      g$values <- matrix(1, nrow = g$nrow, ncol = g$ncol)
      g$nlayer <- 1L
      class(g) <- unique(c("INApestSpatialGrid", class(g)))
      return(g)
    }
  }
  NULL
}

.ina_pt_grid_centres <- function(g) {
  x <- g$xmin + (seq_len(g$ncol) - 0.5) * g$xres
  y <- if (identical(g$row_origin, "ymax")) {
    g$ymax - (seq_len(g$nrow) - 0.5) * g$yres
  } else {
    g$ymin + (seq_len(g$nrow) - 0.5) * g$yres
  }
  zz <- expand.grid(row = seq_len(g$nrow), col = seq_len(g$ncol))
  zz$cell <- (zz$row - 1L) * g$ncol + zz$col
  zz$x <- x[zz$col]
  zz$y <- y[zz$row]
  zz[order(zz$cell), c("cell", "row", "col", "x", "y")]
}

.ina_pt_xy_to_cell <- function(x, y, g) {
  out <- rep(NA_integer_, length(x))
  inside <- is.finite(x) & is.finite(y) &
    x >= g$xmin & x <= g$xmax & y >= g$ymin & y <= g$ymax
  if (!any(inside)) return(out)
  col <- floor((x[inside] - g$xmin) / g$xres) + 1L
  col <- pmin(g$ncol, pmax(1L, col))
  row <- if (identical(g$row_origin, "ymax")) {
    floor((g$ymax - y[inside]) / g$yres) + 1L
  } else {
    floor((y[inside] - g$ymin) / g$yres) + 1L
  }
  row <- pmin(g$nrow, pmax(1L, row))
  out[inside] <- (row - 1L) * g$ncol + col
  out
}

.ina_pt_with_seed <- function(seed, expr) {
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_exists) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_exists) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

.ina_pt_rho <- function(G) {
  if (!length(G)) return(NA_real_)
  max(Mod(eigen(G, only.values = TRUE)$values))
}

.ina_pt_cycle_growth <- function(operators) {
  n <- nrow(operators[[1L]])
  B <- diag(n)
  for (G in operators) B <- G %*% B
  rc <- .ina_pt_rho(B)
  list(CycleOperator = B,
       CycleMultiplier = rc,
       EquivalentPerTimestepMultiplier = rc^(1 / length(operators)))
}

.ina_pt_apply_operators <- function(operators, initial) {
  x <- as.numeric(initial)
  out <- matrix(NA_real_, nrow = length(x), ncol = length(operators))
  for (tt in seq_along(operators)) {
    x <- as.numeric(operators[[tt]] %*% x)
    out[, tt] <- x
  }
  out
}

.ina_pt_overall_event <- function(q, initial) {
  q <- .ina_pt_clip01(q)
  x <- pmax(0, as.numeric(initial))
  exp(sum(x * log(pmax(q, .Machine$double.xmin))))
}

.ina_pt_poisson_extinction_horizon <- function(operators) {
  q <- rep(0, nrow(operators[[1L]]))
  for (tt in rev(seq_along(operators))) {
    G <- operators[[tt]]
    q <- exp(colSums(G * (q - 1)))
    q <- .ina_pt_clip01(q)
  }
  q
}

.ina_pt_poisson_extinction_static <- function(G, generations = 100,
                                               tolerance = 1e-12) {
  q <- rep(0, nrow(G))
  for (tt in seq_len(generations)) {
    qo <- q
    q <- exp(colSums(G * (q - 1)))
    q <- .ina_pt_clip01(q)
    if (max(abs(q - qo)) < tolerance) break
  }
  q
}

.ina_pt_poisson_noescape_horizon <- function(operators, export_vectors) {
  g <- rep(1, nrow(operators[[1L]]))
  for (tt in rev(seq_along(operators))) {
    G <- operators[[tt]]
    e <- export_vectors[[tt]]
    if (is.null(e)) e <- rep(0, nrow(G))
    g <- exp(-pmax(0, e) + colSums(G * (g - 1)))
    g <- .ina_pt_clip01(g)
  }
  g
}

.ina_pt_first_moment_escape <- function(operators, initial, export_vectors) {
  x <- as.numeric(initial)
  mu <- numeric(length(operators))
  for (tt in seq_along(operators)) {
    e <- export_vectors[[tt]]
    if (is.null(e)) e <- rep(0, length(x))
    mu[tt] <- sum(x * e)
    x <- as.numeric(operators[[tt]] %*% x)
  }
  data.frame(
    timestep = seq_along(operators),
    expected_successful_escapes = mu,
    cumulative_expected_successful_escapes = cumsum(mu),
    poisson_escape_probability = 1 - exp(-cumsum(mu))
  )
}

.ina_pt_classify <- function(rho, tol = 1e-6) {
  if (!is.finite(rho)) return(NA_character_)
  if (rho < 1 - tol) "decline when rare"
  else if (rho > 1 + tol) "growth when rare"
  else "approximately replacement when rare"
}

.ina_pt_any_nonzero <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.function(x)) return(TRUE)
  any(as.numeric(x) != 0, na.rm = TRUE)
}

.ina_pt_persistence_requested <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.function(x)) return(TRUE)
  any(!is.na(as.numeric(x)))
}

.ina_pt_retention_loss_possible <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.function(x)) return(TRUE)
  any(as.numeric(x) < 1, na.rm = TRUE)
}

.ina_pt_sd_nonzero <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.function(x)) return(TRUE)
  any(as.numeric(x) != 0, na.rm = TRUE)
}

# The point simulators create stochastic parameter variation when several *SD
# arguments are left NULL. The analytical companion deliberately works at the
# nominal mean parameter values, matching the existing INApestAnalytical
# convention, so flag implicit simulator SDs as well as explicitly supplied SDs.
.ina_pt_implicit_sd_names <- function(DetectionProb, DetectionSD,
                                      ManageProb, ManageSD,
                                      MortalityProb, MortalitySD,
                                      FecundityReduction, FecundityReductionSD,
                                      SpreadReduction, SpreadReductionSD) {
  out <- character(0)
  numeric_nonzero <- function(x) !is.function(x) && .ina_pt_any_nonzero(x)
  if (is.null(DetectionSD) && numeric_nonzero(DetectionProb))
    out <- c(out, "DetectionSD")
  if (is.null(ManageSD) && numeric_nonzero(ManageProb))
    out <- c(out, "ManageSD")
  if (is.null(MortalitySD) && numeric_nonzero(MortalityProb))
    out <- c(out, "MortalitySD")
  if (is.null(FecundityReductionSD) && !is.function(FecundityReduction) &&
      .ina_pt_any_nonzero(FecundityReduction))
    out <- c(out, "FecundityReductionSD")
  # Current point simulators use (1 - SpreadReduction) / 10 when the SD is
  # omitted, so this default is stochastic unless the mean is exactly 1.
  if (is.null(SpreadReductionSD) && !is.function(SpreadReduction) &&
      any((1 - as.numeric(SpreadReduction)) != 0, na.rm = TRUE))
    out <- c(out, "SpreadReductionSD")
  unique(out)
}

.ina_pt_info_mode <- function(InformationMode, ManageProb, InitialInfo,
                              DetectionProb, InfoRadius, InfoTransferProb,
                              InfoKernel, InfoRetentionProb,
                              InfoPersistenceSteps) {
  if (InformationMode != "auto") return(InformationMode)
  info_process <- .ina_pt_any_nonzero(DetectionProb) ||
    .ina_pt_any_nonzero(InitialInfo) || InfoRadius > 0 ||
    .ina_pt_any_nonzero(InfoTransferProb) || !is.null(InfoKernel) ||
    .ina_pt_persistence_requested(InfoPersistenceSteps) ||
    .ina_pt_retention_loss_possible(InfoRetentionProb)
  if (.ina_pt_any_nonzero(ManageProb) && info_process) "dynamic"
  else if (.ina_pt_any_nonzero(ManageProb) && .ina_pt_any_nonzero(InitialInfo)) "dynamic"
  else "none"
}

.ina_pt_expected_parameter <- function(mean_value, SpatialSurface, points,
                                       timestep, perm, Ntimesteps, name) {
  mu <- .ipp_resolve(mean_value, points, timestep, perm, Ntimesteps, name)
  spatial <- .ipp_spatial_value(points$x, points$y, SpatialSurface, timestep,
                                paste0(name, "Spatial"))
  .ina_pt_clip01(mu * spatial)
}

.ina_pt_point_values <- function(x, points, timestep, Ntimesteps, name) {
  .ipp_resolve(x, points, timestep, 1L, Ntimesteps, name)
}

.ina_pt_transition_values <- function(x, points, timestep, Ntimesteps,
                                      Nstages, name) {
  .ipptm_resolve_stage_schedule(x, points, timestep, 1L,
                                Ntimesteps, Nstages, name)
}

.ina_pt_transition_prob_spatial_mean <- function(x, SpatialSurface, points,
                                                 timestep, Ntimesteps,
                                                 Nstages, name) {
  mu <- .ina_pt_transition_values(x, points, timestep, Ntimesteps, Nstages, name)
  spatial <- .ipp_spatial_value(points$x, points$y, SpatialSurface, timestep,
                                paste0(name, "Spatial"))
  .ina_pt_clip01(mu * spatial)
}

.ina_pt_kernel_contract <- function(kernel, representatives, analysis_grid,
                                    HabitatSuitability, HabitatSearchRadius,
                                    HabitatSearchCandidates,
                                    PropaguleEstablishment, EnvEstabProb,
                                    timestep, Ntimesteps,
                                    KernelSamples, seed,
                                    destination_stage = NULL,
                                    apply_habitat_search = TRUE,
                                    apply_establishment = TRUE,
                                    TransitionEstablishment = 1) {
  nt <- nrow(representatives)
  nc <- analysis_grid$nrow * analysis_grid$ncol
  P <- matrix(0, nt, nc)
  outside <- numeric(nt)
  mean_success <- numeric(nt)

  .ina_pt_with_seed(seed, {
    for (i in seq_len(nt)) {
      parent <- representatives[i, , drop = FALSE]
      parents <- parent[rep(1L, KernelSamples), , drop = FALSE]
      if (is.null(kernel)) {
        dx <- dy <- rep(0, KernelSamples)
      } else {
        z <- .ipp_draw_displacement(kernel, parents, timestep, 1L, "PointKernel")
        dx <- z$dx; dy <- z$dy
      }
      px <- parents$x + dx
      py <- parents$y + dy
      if (apply_habitat_search && !is.null(HabitatSuitability) && HabitatSearchRadius > 0) {
        dest <- .ipp_habitat_search(px, py, HabitatSuitability,
                                    HabitatSearchRadius, timestep,
                                    HabitatSearchCandidates)
      } else {
        hv <- .ipp_spatial_value(px, py, HabitatSuitability, timestep,
                                 "HabitatSuitability")
        dest <- data.frame(x = px, y = py, habitat = hv,
                           habitat_nudged = FALSE)
      }
      candidates <- data.frame(x = dest$x, y = dest$y)
      if (!is.null(destination_stage)) candidates$stage <- destination_stage
      if ("stage" %in% names(parents) && !"stage" %in% names(candidates))
        candidates$stage <- parents$stage
      if (apply_establishment) {
        pe <- .ipp_resolve(PropaguleEstablishment, candidates, timestep, 1L,
                           Ntimesteps, "PropaguleEstablishment")
        ee <- .ipp_resolve(EnvEstabProb, candidates, timestep, 1L,
                           Ntimesteps, "EnvEstabProb")
        te <- .ipp_resolve(TransitionEstablishment, candidates, timestep, 1L,
                           Ntimesteps, "TransitionEstablishment")
        w <- .ina_pt_clip01(pe * ee * te * dest$habitat)
      } else {
        w <- rep(1, KernelSamples)
      }
      cell <- .ina_pt_xy_to_cell(dest$x, dest$y, analysis_grid)
      inside <- !is.na(cell)
      if (any(inside)) {
        sums <- rowsum(w[inside], group = cell[inside], reorder = FALSE)
        P[i, as.integer(rownames(sums))] <- as.numeric(sums[, 1L]) / KernelSamples
      }
      outside[i] <- sum(w[!inside]) / KernelSamples
      mean_success[i] <- mean(w)
    }
  })

  list(internal = P, outside = outside, mean_success = mean_success,
       samples = KernelSamples)
}

.ina_pt_direct_info_matrix <- function(representatives, InfoRadius,
                                       InfoTransferProb, InfoKernel,
                                       timestep, Ntimesteps) {
  nt <- nrow(representatives)
  C <- matrix(0, nt, nt)
  if (nt == 0L || (is.null(InfoKernel) && InfoRadius <= 0)) return(C)
  for (i in seq_len(nt)) {
    d <- sqrt((representatives$x - representatives$x[i])^2 +
              (representatives$y - representatives$y[i])^2)
    if (!is.null(InfoKernel)) {
      src <- data.frame(info_id = NA_integer_, source_point_id = i,
                        x = representatives$x[i], y = representatives$y[i],
                        created_timestep = NA_integer_,
                        last_known_timestep = NA_integer_, active = TRUE)
      for (j in seq_len(nt)) {
        p <- InfoKernel(distance = d[j], source = src,
                        target = representatives[j, , drop = FALSE],
                        timestep = timestep, perm = 1L)
        C[i, j] <- .ina_pt_clip01(as.numeric(p)[1L])
      }
    } else {
      p0 <- .ipp_resolve(InfoTransferProb, representatives, timestep, 1L,
                         Ntimesteps, "InfoTransferProb")
      C[i, ] <- ifelse(d <= InfoRadius, p0, 0)
    }
  }
  C
}

.ina_pt_dynamic_operator <- function(Parent0, ParentH, Recruit0, RecruitH,
                                     representatives, DetectionProb,
                                     DetectionSpatial, InfoRetentionProb,
                                     InfoRadius, InfoTransferProb, InfoKernel,
                                     timestep, Ntimesteps) {
  nt <- nrow(Parent0)
  mats <- list(Parent0, ParentH, Recruit0, RecruitH)
  if (any(vapply(mats, function(x) !all(dim(x) == c(nt, nt)), logical(1))))
    stop("Point dynamic-information matrices must all be nt x nt.")

  D <- .ina_pt_expected_parameter(DetectionProb, DetectionSpatial,
                                  representatives, timestep, 1L,
                                  Ntimesteps, "DetectionProb")
  IR <- .ina_pt_point_values(InfoRetentionProb, representatives, timestep,
                             Ntimesteps, "InfoRetentionProb")
  IR <- .ina_pt_clip01(IR)
  C <- .ina_pt_direct_info_matrix(representatives, InfoRadius,
                                  InfoTransferProb, InfoKernel,
                                  timestep, Ntimesteps)
  G <- matrix(0, 2L * nt, 2L * nt)
  U <- seq_len(nt); H <- nt + seq_len(nt)

  for (i in seq_len(nt)) for (j in seq_len(nt)) {
    # Uninformed biological parent and its recruits can become informed only
    # through direct detection by the end of the timestep.
    w <- Parent0[j, i]
    if (w != 0) {
      h <- D[j]
      G[U[j], U[i]] <- G[U[j], U[i]] + w * (1 - h)
      G[H[j], U[i]] <- G[H[j], U[i]] + w * h
    }
    w <- Recruit0[j, i]
    if (w != 0) {
      h <- D[j]
      G[U[j], U[i]] <- G[U[j], U[i]] + w * (1 - h)
      G[H[j], U[i]] <- G[H[j], U[i]] + w * h
    }

    # An informed surviving parent remains the same biological individual, so
    # information retention follows the parent even when analytical types are
    # spatially aggregated. Detection can refresh information after retention.
    w <- ParentH[j, i]
    if (w != 0) {
      h <- IR[j] + (1 - IR[j]) * D[j]
      G[U[j], H[i]] <- G[U[j], H[i]] + w * (1 - h)
      G[H[j], H[i]] <- G[H[j], H[i]] + w * h
    }

    # Recruits are new biological individuals and do not inherit the parent's
    # information merely because they occupy the same analytical cell/type.
    # They can receive direct point-to-point information from the informed
    # source (SEAM analogue) and/or be detected at the end of the timestep.
    w <- RecruitH[j, i]
    if (w != 0) {
      h <- 1 - (1 - D[j]) * (1 - C[i, j])
      G[U[j], H[i]] <- G[U[j], H[i]] + w * (1 - h)
      G[H[j], H[i]] <- G[H[j], H[i]] + w * h
    }
  }
  attr(G, "note") <- paste(
    "Point-grid informed/uninformed low-density operator; parent information",
    "is retained on the same biological individual, while recruits receive",
    "information only through direct transfer/detection; persistent information-only",
    "sites and multi-lineage information correlations are omitted."
  )
  G
}


.ina_pt_persistence_values <- function(InfoPersistenceSteps, representatives,
                                       timestep, Ntimesteps) {
  z <- .ipp_resolve(InfoPersistenceSteps, representatives, timestep, 1L,
                    Ntimesteps, "InfoPersistenceSteps")
  bad <- !is.na(z) & ((!is.finite(z) & !is.infinite(z)) | z < 0 |
                      (is.finite(z) & z != floor(z)))
  if (any(bad))
    stop("InfoPersistenceSteps must resolve to non-negative whole numbers, Inf, or NA.")
  as.numeric(z)
}

.ina_pt_persistence_profile <- function(InfoPersistenceSteps, representatives,
                                        Ntimesteps) {
  values <- lapply(seq_len(Ntimesteps), function(tt)
    .ina_pt_persistence_values(InfoPersistenceSteps, representatives, tt,
                               Ntimesteps))
  allv <- unlist(values, use.names = FALSE)
  requested <- any(!is.na(allv))
  finite <- allv[!is.na(allv) & is.finite(allv)]
  max_age <- if (requested) max(1L, if (length(finite)) as.integer(max(finite)) else 1L) else NULL
  list(values = values, requested = requested, max_age = max_age)
}

# Explicit point-information clock. Parent matrices contain the same biological
# individual after survival/progression; recruit matrices contain new
# individuals. This distinction lets the information clock follow a moving
# parent without incorrectly treating stage movement as information transfer.
.ina_pt_programmed_operator <- function(Parent0, ParentH, Recruit0, RecruitH,
                                        representatives, DetectionProbByType,
                                        InfoRetentionProbByType,
                                        InfoPersistenceStepsByType,
                                        InfoTransferMatrix, layout) {
  B <- nrow(Parent0)
  mats <- list(Parent0, ParentH, Recruit0, RecruitH)
  if (any(vapply(mats, function(x) !all(dim(x) == c(B, B)), logical(1))))
    stop("Point programmed-information matrices must all be B x B.")
  d <- .ina_pt_clip01(rep_len(as.numeric(DetectionProbByType), B))
  ir <- .ina_pt_clip01(rep_len(as.numeric(InfoRetentionProbByType), B))
  K <- rep_len(as.numeric(InfoPersistenceStepsByType), B)
  C <- as.matrix(InfoTransferMatrix)
  if (!all(dim(C) == c(B, B))) stop("InfoTransferMatrix must be B x B.")
  if (layout$base_types != B) stop("Programmed-information layout has the wrong base type count.")

  G <- matrix(0, layout$size, layout$size)
  Hblock <- function(age, b) layout$H[[age]][b]

  add_uninformed <- function(dst, src_col, weight) {
    if (weight == 0) return(invisible(NULL))
    G[Hblock(1L, dst), src_col] <<- G[Hblock(1L, dst), src_col] + weight * d[dst]
    G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + weight * (1 - d[dst])
    invisible(NULL)
  }

  add_parent_informed <- function(dst, src_col, weight, kind, age = NA_integer_) {
    if (weight == 0) return(invisible(NULL))
    G[Hblock(1L, dst), src_col] <<- G[Hblock(1L, dst), src_col] + weight * d[dst]
    rem <- weight * (1 - d[dst])
    if (rem == 0) return(invisible(NULL))

    kval <- K[dst]
    if (identical(kind, "X")) {
      if (is.na(kval)) {
        G[layout$X[dst], src_col] <<- G[layout$X[dst], src_col] + rem * ir[dst]
        G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + rem * (1 - ir[dst])
      } else if (is.infinite(kval)) {
        G[layout$X[dst], src_col] <<- G[layout$X[dst], src_col] + rem
      } else {
        # The source is informed but has no valid local-evidence clock.
        G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + rem
      }
      return(invisible(NULL))
    }

    if (identical(kind, "Overflow")) {
      if (is.na(kval)) {
        G[layout$Overflow[dst], src_col] <<- G[layout$Overflow[dst], src_col] + rem * ir[dst]
        G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + rem * (1 - ir[dst])
      } else if (is.infinite(kval)) {
        G[layout$Overflow[dst], src_col] <<- G[layout$Overflow[dst], src_col] + rem
      } else {
        G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + rem
      }
      return(invisible(NULL))
    }

    if (is.na(kval)) {
      next_idx <- if (age < layout$max_age) Hblock(age + 1L, dst) else layout$Overflow[dst]
      G[next_idx, src_col] <<- G[next_idx, src_col] + rem * ir[dst]
      G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + rem * (1 - ir[dst])
    } else if (is.infinite(kval)) {
      next_idx <- if (age < layout$max_age) Hblock(age + 1L, dst) else layout$Overflow[dst]
      G[next_idx, src_col] <<- G[next_idx, src_col] + rem
    } else if (age >= kval) {
      G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + rem
    } else {
      next_idx <- if (age < layout$max_age) Hblock(age + 1L, dst) else layout$Overflow[dst]
      G[next_idx, src_col] <<- G[next_idx, src_col] + rem
    }
    invisible(NULL)
  }

  add_recruit_from_informed <- function(dst, src, src_col, weight) {
    if (weight == 0) return(invisible(NULL))
    pdet <- d[dst]
    pseam <- .ina_pt_clip01(C[src, dst])
    G[Hblock(1L, dst), src_col] <<- G[Hblock(1L, dst), src_col] + weight * pdet
    G[layout$X[dst], src_col] <<- G[layout$X[dst], src_col] + weight * (1 - pdet) * pseam
    G[layout$U[dst], src_col] <<- G[layout$U[dst], src_col] + weight * (1 - pdet) * (1 - pseam)
    invisible(NULL)
  }

  for (src in seq_len(B)) {
    # Uninformed biological parent and its recruits: only direct detection can
    # create information by the next timestep.
    scol <- layout$U[src]
    for (dst in which(Parent0[, src] != 0)) add_uninformed(dst, scol, Parent0[dst, src])
    for (dst in which(Recruit0[, src] != 0)) add_uninformed(dst, scol, Recruit0[dst, src])

    # Informed source without a valid evidence clock.
    scol <- layout$X[src]
    for (dst in which(ParentH[, src] != 0)) add_parent_informed(dst, scol, ParentH[dst, src], "X")
    for (dst in which(RecruitH[, src] != 0)) add_recruit_from_informed(dst, src, scol, RecruitH[dst, src])

    # Explicit local-evidence ages.
    for (age in seq_len(layout$max_age)) {
      scol <- layout$H[[age]][src]
      for (dst in which(ParentH[, src] != 0)) add_parent_informed(dst, scol, ParentH[dst, src], "H", age)
      for (dst in which(RecruitH[, src] != 0)) add_recruit_from_informed(dst, src, scol, RecruitH[dst, src])
    }

    scol <- layout$Overflow[src]
    for (dst in which(ParentH[, src] != 0)) add_parent_informed(dst, scol, ParentH[dst, src], "Overflow")
    for (dst in which(RecruitH[, src] != 0)) add_recruit_from_informed(dst, src, scol, RecruitH[dst, src])
  }
  attr(G, "note") <- paste(
    "Explicit point time-since-local-evidence operator; the information clock",
    "follows the same biological parent through survival or stage movement.",
    "Persistent information-only sites after parent removal remain a shared-lineage approximation.")
  attr(G, "programmed_information_layout") <- layout
  G
}

.ina_pt_initial_dynamic_state <- function(base_initial, representatives,
                                          InitialInfo, DetectionProb,
                                          DetectionSpatial, ApplyInitialDetection,
                                          Ntimesteps) {
  nt <- length(base_initial)
  if (is.null(InitialInfo)) info <- rep(0, nt)
  else if (length(InitialInfo) == 1L) info <- rep(as.numeric(InitialInfo), nt)
  else if (length(InitialInfo) == nt) info <- as.numeric(InitialInfo)
  else stop("For point analytical types, InitialInfo must be scalar or one value per analytical type.")
  info <- .ina_pt_clip01(info)
  if (ApplyInitialDetection) {
    D <- .ina_pt_expected_parameter(DetectionProb, DetectionSpatial,
                                    representatives, 1L, 1L, Ntimesteps,
                                    "DetectionProb")
    info <- info + (1 - info) * D
  }
  c(base_initial * (1 - info), base_initial * info)
}

###############################################################################
### INApestMetaPoint analytical implementation
###############################################################################


.ina_pt_initial_info_raw <- function(InitialPoints, InitialInfo) {
  n <- nrow(InitialPoints)
  if (is.null(InitialInfo)) {
    if ("have_info" %in% names(InitialPoints)) {
      z <- as.numeric(as.logical(InitialPoints$have_info))
      z[is.na(z)] <- 0
      return(z)
    }
    return(rep(0, n))
  }
  if (is.logical(InitialInfo) && length(InitialInfo) == n) {
    z <- as.numeric(InitialInfo); z[is.na(z)] <- 0; return(z)
  }
  if (is.numeric(InitialInfo)) {
    # Match the point simulator: numeric InitialInfo denotes row indices.
    return(as.numeric(seq_len(n) %in% as.integer(InitialInfo)))
  }
  stop("InitialInfo must be NULL, a logical vector per initial point, or row indices.")
}

.ina_pt_metapoint_representatives <- function(InitialPoints, analysis_grid) {
  if (is.null(analysis_grid)) {
    return(data.frame(
      id = 1L,
      x = if (nrow(InitialPoints)) mean(InitialPoints$x) else 0,
      y = if (nrow(InitialPoints)) mean(InitialPoints$y) else 0,
      stage = if ("stage" %in% names(InitialPoints) && nrow(InitialPoints))
        as.character(InitialPoints$stage[1L]) else "default",
      cell = 1L,
      stringsAsFactors = FALSE
    ))
  }
  cc <- .ina_pt_grid_centres(analysis_grid)
  stage_value <- if ("stage" %in% names(InitialPoints) && nrow(InitialPoints))
    as.character(InitialPoints$stage[1L]) else "default"
  data.frame(id = cc$cell, x = cc$x, y = cc$y,
             stage = stage_value, cell = cc$cell,
             stringsAsFactors = FALSE)
}

.ina_pt_metapoint_initial <- function(InitialPoints, analysis_grid) {
  if (is.null(analysis_grid)) return(c(nrow(InitialPoints)))
  cell <- .ina_pt_xy_to_cell(InitialPoints$x, InitialPoints$y, analysis_grid)
  tabulate(cell[!is.na(cell)], nbins = analysis_grid$nrow * analysis_grid$ncol)
}

.ina_pt_metapoint_operator_step <- function(
    timestep, Ntimesteps, representatives, analysis_grid,
    Survival, PropaguleProduction, PropaguleEstablishment, EnvEstabProb,
    SDDkernel, LDDkernel, LDDrate,
    HabitatSuitability, HabitatSearchRadius, HabitatSearchCandidates,
    DetectionProb, DetectionSpatial,
    ManageProb, ManageSpatial,
    MortalityProb, MortalitySpatial,
    FecundityReduction, FecundityReductionSpatial,
    SpreadReduction, SpreadReductionSpatial, SpreadReductionAppliesTo,
    InfoRetentionProb, InfoRadius, InfoTransferProb, InfoKernel,
    KernelSamples, PointSeed, mode, OutsideEstablishmentProb) {

  nt <- nrow(representatives)
  s <- .ina_pt_clip01(.ina_pt_point_values(Survival, representatives, timestep,
                                            Ntimesteps, "Survival"))
  p <- pmax(0, .ina_pt_point_values(PropaguleProduction, representatives,
                                    timestep, Ntimesteps,
                                    "PropaguleProduction"))
  a <- .ina_pt_expected_parameter(ManageProb, ManageSpatial, representatives,
                                  timestep, 1L, Ntimesteps, "ManageProb")
  m <- .ina_pt_expected_parameter(MortalityProb, MortalitySpatial, representatives,
                                  timestep, 1L, Ntimesteps, "MortalityProb")
  f <- .ina_pt_expected_parameter(FecundityReduction, FecundityReductionSpatial,
                                  representatives, timestep, 1L, Ntimesteps,
                                  "FecundityReduction")
  g <- .ina_pt_expected_parameter(SpreadReduction, SpreadReductionSpatial,
                                  representatives, timestep, 1L, Ntimesteps,
                                  "SpreadReduction")
  r <- as.numeric(LDDrate)
  if (length(r) != 1L || !is.finite(r) || r < 0 || r > 1)
    stop("INApestMetaPoint analytical support requires scalar LDDrate in [0,1], matching the simulator.")

  if (is.null(analysis_grid)) {
    pe <- .ina_pt_point_values(PropaguleEstablishment, representatives, timestep,
                               Ntimesteps, "PropaguleEstablishment")
    ee <- .ina_pt_point_values(EnvEstabProb, representatives, timestep,
                               Ntimesteps, "EnvEstabProb")
    estab <- .ina_pt_clip01(pe * ee)
    PS <- diag(estab, nt)
    PL <- diag(estab, nt)
    outS <- outL <- rep(0, nt)
  } else {
    sdd <- .ina_pt_kernel_contract(
      SDDkernel, representatives, analysis_grid, HabitatSuitability,
      HabitatSearchRadius, HabitatSearchCandidates,
      PropaguleEstablishment, EnvEstabProb,
      timestep, Ntimesteps, KernelSamples,
      seed = PointSeed + 100000L * timestep + 101L,
      apply_habitat_search = TRUE, apply_establishment = TRUE)
    PS <- sdd$internal; outS <- sdd$outside
    if (r > 0) {
      if (is.null(LDDkernel)) stop("LDDkernel is required when LDDrate > 0.")
      ldd <- .ina_pt_kernel_contract(
        LDDkernel, representatives, analysis_grid, HabitatSuitability,
        HabitatSearchRadius, HabitatSearchCandidates,
        PropaguleEstablishment, EnvEstabProb,
        timestep, Ntimesteps, KernelSamples,
        seed = PointSeed + 100000L * timestep + 202L,
        apply_habitat_search = TRUE, apply_establishment = TRUE)
      PL <- ldd$internal; outL <- ldd$outside
    } else {
      PL <- matrix(0, nt, ncol(PS)); outL <- rep(0, nt)
    }
  }

  P0 <- (1 - r) * PS + r * PL
  out0 <- (1 - r) * outS + r * outL
  if (SpreadReductionAppliesTo == "LDD") {
    P1 <- (1 - r) * PS + r * sweep(PL, 1, 1 - g, `*`)
    out1 <- (1 - r) * outS + r * (1 - g) * outL
  } else {
    P1 <- sweep(P0, 1, 1 - g, `*`)
    out1 <- (1 - g) * out0
  }

  q0 <- s
  q1 <- s * (1 - m)
  lam0 <- p
  lam1 <- p * (1 - f)

  Parent0 <- diag(q0, nrow = nt, ncol = nt)
  Parent1 <- diag(q1, nrow = nt, ncol = nt)
  Recruit0 <- t(sweep(P0, 1, q0 * lam0, `*`))
  Recruit1 <- t(sweep(P1, 1, q1 * lam1, `*`))
  G0 <- Parent0 + Recruit0
  G1 <- Parent1 + Recruit1
  ParentH <- sweep(Parent0, 2, 1 - a, `*`) + sweep(Parent1, 2, a, `*`)
  RecruitH <- sweep(Recruit0, 2, 1 - a, `*`) + sweep(Recruit1, 2, a, `*`)
  GH <- ParentH + RecruitH

  pout <- if (length(OutsideEstablishmentProb) == 1L) {
    rep(as.numeric(OutsideEstablishmentProb), nt)
  } else if (length(OutsideEstablishmentProb) == nt) {
    as.numeric(OutsideEstablishmentProb)
  } else stop("OutsideEstablishmentProb must be scalar or one value per point analytical type.")
  pout <- .ina_pt_clip01(pout)
  e0 <- q0 * lam0 * out0 * pout
  e1 <- q1 * lam1 * out1 * pout
  eH <- (1 - a) * e0 + a * e1

  if (mode == "dynamic") {
    G <- .ina_pt_dynamic_operator(Parent0, ParentH, Recruit0, RecruitH,
                                  representatives, DetectionProb,
                                  DetectionSpatial, InfoRetentionProb,
                                  InfoRadius, InfoTransferProb, InfoKernel,
                                  timestep, Ntimesteps)
    e <- c(e0, eH)
  } else if (mode == "all_informed") {
    G <- GH; e <- eH
  } else {
    G <- G0; e <- e0
  }

  list(G = G, export = e, G0 = G0, GH = GH,
       Parent0 = Parent0, ParentH = ParentH,
       Recruit0 = Recruit0, RecruitH = RecruitH,
       branch = list(a = a, q0 = q0, q1 = q1,
                     mu0 = rowSums(P0) * lam0,
                     mu1 = rowSums(P1) * lam1,
                     e0 = e0, e1 = e1))
}

.ina_pt_metapoint_exact_extinction <- function(branch_steps, mode,
                                                generations = 100) {
  if (!(mode %in% c("none", "all_informed"))) return(NULL)
  if (length(branch_steps[[1L]]$q0) != 1L) return(NULL)
  stepfun <- function(br, q) {
    a <- if (mode == "all_informed") br$a else 0
    f0 <- (1 - br$q0) + br$q0 * q * exp(br$mu0 * (q - 1))
    f1 <- (1 - br$q1) + br$q1 * q * exp(br$mu1 * (q - 1))
    .ina_pt_clip01((1 - a) * f0 + a * f1)
  }
  qh <- 0
  for (tt in rev(seq_along(branch_steps))) qh <- stepfun(branch_steps[[tt]], qh)
  static <- all(vapply(branch_steps[-1L], function(z)
    isTRUE(all.equal(z, branch_steps[[1L]], tolerance = 0)), logical(1)))
  if (length(branch_steps) == 1L) static <- TRUE
  qe <- NA_real_
  if (static) {
    q <- 0
    for (gg in seq_len(generations)) {
      qo <- q; q <- stepfun(branch_steps[[1L]], q)
      if (abs(q - qo) < 1e-12) break
    }
    qe <- q
  }
  list(horizon = qh, eventual = qe,
       method = "exact single-type survival + Poisson-recruit branching PGF")
}

INApestMetaPointAnalytical <- function(
    Ntimesteps = 10,
    InitialPoints,
    InitialInfo = NULL,
    InformationMode = c("auto", "dynamic", "all_informed", "none"),
    ApplyInitialDetection = TRUE,
    Survival = 1,
    PropaguleProduction,
    PropaguleEstablishment = 1,
    EnvEstabProb = 1,
    SDDkernel,
    LDDkernel = NULL,
    LDDrate = 0,
    HabitatSuitability = NULL,
    HabitatSearchRadius = 0,
    HabitatSearchCandidates = 128,
    LocalK = Inf,
    KRadius = 0,
    DetectionProb = 0,
    DetectionSD = NULL,
    DetectionSpatial = NULL,
    ManageProb = 0,
    ManageSD = NULL,
    ManageSpatial = NULL,
    MortalityProb = 0,
    MortalitySD = NULL,
    MortalitySpatial = NULL,
    FecundityReduction = 0,
    FecundityReductionSD = NULL,
    FecundityReductionSpatial = NULL,
    SpreadReduction = 0,
    SpreadReductionSD = NULL,
    SpreadReductionSpatial = NULL,
    SpreadReductionAppliesTo = c("LDD", "all"),
    InfoRadius = 0,
    InfoTransferProb = 0,
    InfoKernel = NULL,
    InfoRetentionProb = 1,
    InfoPersistenceSteps = NA,
    ExternalInfoProb = 0,
    OngoingExternalInfo = FALSE,
    ExternalIncursionGenerator = NULL,
    OngoingExternalInvasion = FALSE,
    PointAnalysisGrid = NULL,
    MaxAnalysisCells = 400L,
    KernelSamples = 5000L,
    PointSeed = 1L,
    OutsideEstablishmentProb = 1,
    ExtinctionGenerations = 100,
    ReturnOperators = FALSE) {

  .ina_pt_require_point_helpers(FALSE)
  InformationMode <- match.arg(InformationMode)
  SpreadReductionAppliesTo <- match.arg(SpreadReductionAppliesTo)
  if (!is.data.frame(InitialPoints) || !all(c("x", "y") %in% names(InitialPoints)))
    stop("InitialPoints must be a data.frame with x and y.")
  if (any(!is.finite(InitialPoints$x)) || any(!is.finite(InitialPoints$y)))
    stop("InitialPoints x and y must be finite.")
  if (!is.numeric(Ntimesteps) || length(Ntimesteps) != 1L || !is.finite(Ntimesteps) ||
      Ntimesteps < 1L || Ntimesteps != floor(Ntimesteps))
    stop("Ntimesteps must be a positive integer.")
  Ntimesteps <- as.integer(Ntimesteps)
  if (!is.numeric(KernelSamples) || length(KernelSamples) != 1L ||
      !is.finite(KernelSamples) || KernelSamples < 100L ||
      KernelSamples != floor(KernelSamples))
    stop("KernelSamples must be an integer >= 100.")
  KernelSamples <- as.integer(KernelSamples)
  if (!is.function(SDDkernel))
    stop("SDDkernel must be a dispersal-kernel function, matching INApestMetaPoint.")
  if (!is.null(LDDkernel) && !is.function(LDDkernel))
    stop("LDDkernel must be NULL or a dispersal-kernel function, matching INApestMetaPoint.")
  if (length(LDDrate) != 1L || !is.finite(LDDrate) || LDDrate < 0 || LDDrate > 1)
    stop("LDDrate must be a scalar in [0,1] for INApestMetaPoint.")
  if (LDDrate > 0 && is.null(LDDkernel))
    stop("LDDkernel is required when LDDrate > 0, matching INApestMetaPoint.")
  .ipp_validate_probability_schedule(FecundityReduction, Ntimesteps,
                                     "FecundityReduction")
  .ipp_validate_probability_schedule(FecundityReductionSD, Ntimesteps,
                                     "FecundityReductionSD", allow_null = TRUE)
  if (isTRUE(OngoingExternalInvasion))
    stop("OngoingExternalInvasion is additive immigration and is not represented by the current rare-lineage point analytical solution; use the stochastic point simulator for that process.")
  if (isTRUE(OngoingExternalInfo) && .ina_pt_any_nonzero(ExternalInfoProb))
    stop("OngoingExternalInfo is not yet represented in the point information-state operator; use the stochastic point simulator when ongoing external information is active.")

  analysis_grid <- .ina_pt_analysis_grid(
    PointAnalysisGrid, HabitatSuitability, DetectionSpatial, ManageSpatial,
    MortalitySpatial, FecundityReductionSpatial, SpreadReductionSpatial)
  if (!is.null(analysis_grid)) {
    ncells <- analysis_grid$nrow * analysis_grid$ncol
    if (is.null(PointAnalysisGrid) && ncells > as.integer(MaxAnalysisCells))
      stop("The automatically inherited point grid has ", ncells,
           " cells. Supply a coarser PointAnalysisGrid (<= MaxAnalysisCells) ",
           "rather than unintentionally contracting kernels over the full fine habitat grid.")
  }

  if (is.null(analysis_grid) && any(vapply(list(HabitatSuitability, DetectionSpatial,
                                                   ManageSpatial, MortalitySpatial,
                                                   FecundityReductionSpatial, SpreadReductionSpatial),
                                              function(z) !is.null(z), logical(1))))
    stop("Spatial habitat/management surfaces require PointAnalysisGrid unless a native INApestSpatialGrid can be inherited automatically.")

  if (is.null(analysis_grid) &&
      (is.function(HabitatSuitability) || is.function(Survival) ||
       is.function(PropaguleProduction) || is.function(PropaguleEstablishment) ||
       is.function(EnvEstabProb) || is.function(DetectionProb) ||
       is.function(ManageProb) || is.function(MortalityProb) ||
       is.function(SpreadReduction) || is.function(InfoTransferProb)))
    stop("Location-dependent/function-valued point parameters require PointAnalysisGrid for analytical contraction.")

  stage_levels <- if ("stage" %in% names(InitialPoints)) unique(as.character(InitialPoints$stage)) else "default"
  custom_point_functions <- any(vapply(list(Survival, PropaguleProduction, PropaguleEstablishment,
                                             EnvEstabProb, DetectionProb, ManageProb, MortalityProb,
                                             FecundityReduction, SpreadReduction, InfoTransferProb),
                                        is.function, logical(1)))
  if (length(stage_levels) > 1L && custom_point_functions)
    stop("INApestMetaPoint analytical contraction currently requires one persistent stage label when custom point-valued functions are used; multiple labels would need stage-stratified analytical types.")

  representatives <- .ina_pt_metapoint_representatives(InitialPoints, analysis_grid)
  base_initial <- .ina_pt_metapoint_initial(InitialPoints, analysis_grid)
  initial_outside <- nrow(InitialPoints) - sum(base_initial)

  rawinfo <- .ina_pt_initial_info_raw(InitialPoints, InitialInfo)
  initial_info_by_type <- if (is.null(analysis_grid)) {
    if (nrow(InitialPoints)) mean(rawinfo) else 0
  } else {
    cell <- .ina_pt_xy_to_cell(InitialPoints$x, InitialPoints$y, analysis_grid)
    z <- numeric(nrow(representatives)); n <- numeric(length(z))
    for (i in which(!is.na(cell))) {
      z[cell[i]] <- z[cell[i]] + rawinfo[i]; n[cell[i]] <- n[cell[i]] + 1
    }
    z[n > 0] <- z[n > 0] / n[n > 0]
    z
  }

  mode <- .ina_pt_info_mode(InformationMode, ManageProb, initial_info_by_type,
                            DetectionProb, InfoRadius, InfoTransferProb,
                            InfoKernel, InfoRetentionProb, InfoPersistenceSteps)
  persistence_profile <- if (mode == "dynamic")
    .ina_pt_persistence_profile(InfoPersistenceSteps, representatives, Ntimesteps)
  else list(requested = FALSE, values = NULL, max_age = NULL)
  programmed_dynamic <- mode == "dynamic" && isTRUE(persistence_profile$requested)
  programmed_layout <- if (programmed_dynamic)
    .ina_programmed_layout(nrow(representatives), persistence_profile$max_age) else NULL
  if (is.null(analysis_grid) && mode == "dynamic" &&
      (InfoRadius > 0 || !is.null(InfoKernel)))
    stop("Distance-based information transfer in INApestMetaPoint requires PointAnalysisGrid so source-target distances are represented rather than collapsed to zero distance.")

  operators <- vector("list", Ntimesteps)
  exports <- vector("list", Ntimesteps)
  branches <- vector("list", Ntimesteps)
  point_step_data <- vector("list", Ntimesteps)
  for (tt in seq_len(Ntimesteps)) {
    st <- .ina_pt_metapoint_operator_step(
      tt, Ntimesteps, representatives, analysis_grid,
      Survival, PropaguleProduction, PropaguleEstablishment, EnvEstabProb,
      SDDkernel, LDDkernel, LDDrate,
      HabitatSuitability, HabitatSearchRadius, HabitatSearchCandidates,
      DetectionProb, DetectionSpatial, ManageProb, ManageSpatial,
      MortalityProb, MortalitySpatial,
      FecundityReduction, FecundityReductionSpatial,
      SpreadReduction, SpreadReductionSpatial, SpreadReductionAppliesTo,
      InfoRetentionProb, InfoRadius, InfoTransferProb, InfoKernel,
      KernelSamples, PointSeed, mode, OutsideEstablishmentProb)
    if (programmed_dynamic) {
      Dtt <- .ina_pt_expected_parameter(DetectionProb, DetectionSpatial,
                                        representatives, tt, 1L, Ntimesteps,
                                        "DetectionProb")
      IRtt <- .ina_pt_point_values(InfoRetentionProb, representatives, tt,
                                   Ntimesteps, "InfoRetentionProb")
      Ctt <- .ina_pt_direct_info_matrix(representatives, InfoRadius,
                                        InfoTransferProb, InfoKernel,
                                        tt, Ntimesteps)
      operators[[tt]] <- .ina_pt_programmed_operator(
        st$Parent0, st$ParentH, st$Recruit0, st$RecruitH, representatives,
        Dtt, IRtt, persistence_profile$values[[tt]], Ctt, programmed_layout)
      exports[[tt]] <- .ina_programmed_expand_vector(st$export, programmed_layout)
    } else {
      operators[[tt]] <- st$G
      exports[[tt]] <- st$export
    }
    branches[[tt]] <- st$branch
    point_step_data[[tt]] <- st
  }

  if (programmed_dynamic) {
    D1 <- .ina_pt_expected_parameter(DetectionProb, DetectionSpatial,
                                     representatives, 1L, 1L, Ntimesteps,
                                     "DetectionProb")
    if (!ApplyInitialDetection) D1[] <- 0
    state0 <- .ina_programmed_initial_state(
      base_initial, initial_info_by_type, D1, seq_along(base_initial),
      programmed_layout, binary_known_presence = FALSE)
  } else if (mode == "dynamic") {
    state0 <- .ina_pt_initial_dynamic_state(
      base_initial, representatives, initial_info_by_type,
      DetectionProb, DetectionSpatial, ApplyInitialDetection, Ntimesteps)
  } else state0 <- base_initial

  growth <- .ina_pt_cycle_growth(operators)
  states <- .ina_pt_apply_operators(operators, state0)
  if (programmed_dynamic) totals <- colSums(states)
  else if (mode == "dynamic") {
    nt <- length(base_initial)
    totals <- colSums(states[seq_len(nt), , drop = FALSE] +
                      states[nt + seq_len(nt), , drop = FALSE])
  } else totals <- colSums(states)

  # Model-specific point branching PGF. Continuous movement has already been
  # contracted to analytical spatial types above; conditional on the biological
  # parent branch, Poisson thinning gives the recruit distribution among those
  # types directly. Thus spatial or time-varying point extinction does not need
  # a generic mean-matched Poisson fallback.
  pgf_steps <- lapply(seq_len(Ntimesteps), function(tt)
    .ina_pt_metapoint_branch_step(
      point_step_data[[tt]], representatives, DetectionProb, DetectionSpatial,
      InfoRetentionProb, InfoRadius, InfoTransferProb, InfoKernel,
      tt, Ntimesteps))
  initial_type <- .ina_pt_initial_type_metapoint(InitialPoints, analysis_grid)
  if (programmed_dynamic) {
    qh <- .ina_programmed_branch_horizon(
      pgf_steps, persistence_profile$values, programmed_layout)
    ext_h <- .ina_programmed_initial_extinction_point(
      qh, InitialPoints, InitialInfo, initial_type, DetectionProb, DetectionSpatial,
      Ntimesteps, programmed_layout, transition = FALSE)
    static <- .ina_branch_static(pgf_steps) &&
      (length(persistence_profile$values) == 1L ||
       all(vapply(persistence_profile$values[-1L], function(z)
         isTRUE(all.equal(z, persistence_profile$values[[1L]], tolerance = 0)),
         logical(1))))
    qe <- if (static) .ina_programmed_branch_eventual(
      pgf_steps[[1L]], persistence_profile$values[[1L]], programmed_layout,
      ExtinctionGenerations) else NULL
    ext_e <- if (is.null(qe)) NA_real_ else .ina_programmed_initial_extinction_point(
      qe, InitialPoints, InitialInfo, initial_type, DetectionProb, DetectionSpatial,
      Ntimesteps, programmed_layout, transition = FALSE)
    ext_method <- if (static)
      "point-grid programmed-information branching PGF" else
      "time-inhomogeneous point-grid programmed-information branching PGF"
  } else {
    qh <- .ina_branch_horizon(pgf_steps, mode)
    ext_h <- .ina_pt_initial_event_probability(
      qh, InitialPoints, InitialInfo, initial_type, mode, DetectionProb,
      DetectionSpatial, Ntimesteps, transition = FALSE)
    static <- .ina_branch_static(pgf_steps)
    qe <- if (static) .ina_branch_eventual(pgf_steps[[1L]], mode,
                                            ExtinctionGenerations) else NULL
    ext_e <- if (is.null(qe)) NA_real_ else .ina_pt_initial_event_probability(
      qe, InitialPoints, InitialInfo, initial_type, mode, DetectionProb,
      DetectionSpatial, Ntimesteps, transition = FALSE)
    ext_method <- if (mode == "dynamic") {
      if (static) "point-grid informed/uninformed branching PGF" else
        "time-inhomogeneous point-grid informed/uninformed branching PGF"
    } else {
      if (static) "point-grid survival + Poisson-recruit branching PGF" else
        "time-inhomogeneous point-grid survival + Poisson-recruit branching PGF"
    }
  }

  escape <- NULL
  if (!is.null(analysis_grid)) {
    first <- .ina_pt_first_moment_escape(operators, state0, exports)
    ng <- .ina_pt_poisson_noescape_horizon(operators, exports)
    escape <- list(
      Method = "point-grid mean-matched branching no-escape recursion",
      BranchingProbabilityByHorizon = 1 - .ina_pt_overall_event(ng, state0),
      FirstMomentPoissonProbabilityByHorizon = tail(first$poisson_escape_probability, 1L),
      ExpectedSuccessfulEscapesByStep = first$expected_successful_escapes,
      ExpectedSuccessfulEscapesCumulative = first$cumulative_expected_successful_escapes)
  }

  diagnostics <- character(0)
  if (!is.null(analysis_grid)) diagnostics <- c(diagnostics,
    paste0("Continuous dispersal kernels were contracted to ", analysis_grid$nrow,
           " x ", analysis_grid$ncol, " analysis cells using ", KernelSamples,
           " Monte Carlo draws per source type and kernel per timestep."))
  else diagnostics <- c(diagnostics,
    "Spatially homogeneous point solution: dispersal distances do not affect total rare-population growth or extinction when habitat, management and boundaries are homogeneous.")
  if (!is.infinite(LocalK) || KRadius > 0) diagnostics <- c(diagnostics,
    "Finite LocalK/KRadius is a local density interaction and is not inserted as ordinary node carrying capacity in the rare-lineage operator; stochastic simulation is required once local crowding matters.")
  if (programmed_dynamic) diagnostics <- c(diagnostics,
    "InfoPersistenceSteps is represented with explicit time-since-local-evidence states. Inf is treated as a programmed window that does not expire; InfoRetentionProb applies only where InfoPersistenceSteps is NA.")
  if (mode == "dynamic") diagnostics <- c(diagnostics,
    "Point information transfer is represented as direct source-to-descendant transfer at analytical grid resolution; persistent information-only sites and shared multi-lineage information effects remain higher-order approximations.")
  explicit_sd <- any(vapply(list(DetectionSD, ManageSD, MortalitySD,
                                   FecundityReductionSD, SpreadReductionSD),
                              .ina_pt_sd_nonzero, logical(1)))
  implicit_sd <- .ina_pt_implicit_sd_names(
    DetectionProb, DetectionSD, ManageProb, ManageSD,
    MortalityProb, MortalitySD, FecundityReduction, FecundityReductionSD,
    SpreadReduction, SpreadReductionSD)
  if (explicit_sd || length(implicit_sd)) diagnostics <- c(diagnostics,
    "As in the existing analytical companion, parameter SDs are not integrated; nominal mean parameter values are used.")
  if (length(implicit_sd)) diagnostics <- c(diagnostics,
    paste0("For direct comparison with INApestMetaPoint stochastic runs, note that NULL ",
           paste(implicit_sd, collapse = ", "),
           " activates simulator default variation. Set the relevant SD argument(s) explicitly to 0 when comparing to the nominal-mean analytical solution."))
  if (initial_outside > 0) diagnostics <- c(diagnostics,
    paste(initial_outside, "initial point(s) fall outside PointAnalysisGrid and are not included in the internal analytical state."))

  result <- list(
    Model = "INApestMetaPoint",
    InformationMode = mode,
    InformationStateMethod = if (programmed_dynamic)
      "explicit time-since-local-evidence point-grid states" else if (mode == "dynamic")
      "uninformed/informed point-grid states" else mode,
    Growth = list(
      EquivalentPerTimestepMultiplier = growth$EquivalentPerTimestepMultiplier,
      CycleMultiplier = growth$CycleMultiplier,
      Classification = .ina_pt_classify(growth$EquivalentPerTimestepMultiplier),
      IntrinsicLocalLambda = NULL,
      IntrinsicLocalLambdaInterpretation = NULL),
    Trajectory = data.frame(timestep = seq_len(Ntimesteps),
                            expected_state_total = totals),
    Extinction = list(Method = ext_method,
                      ProbabilityByHorizon = ext_h,
                      BranchingFadeoutProbability = ext_e),
    Escape = escape,
    Diagnostics = unique(diagnostics),
    ApproximationScope = c(
      "Low-density / rare-invasion screening unless otherwise stated",
      "Continuous-space geometry is exact only in the homogeneous count solution; spatial solutions use kernel contraction to an analysis grid",
      "Finite local capacity, lineage collisions and persistent information-only sites become increasingly important away from rarity",
      "Kernel-contraction Monte Carlo error can be reduced by increasing KernelSamples"),
    PointApproximation = list(
      SpatialMode = if (is.null(analysis_grid)) "homogeneous" else "analysis-grid contraction",
      KernelSamples = if (is.null(analysis_grid)) NA_integer_ else KernelSamples,
      AnalysisGrid = analysis_grid)
  )
  if (ReturnOperators) result$Operators <- operators
  class(result) <- "INApestAnalyticalResult"
  result
}

###############################################################################
### INApestPointTransitionMatrix analytical implementation
###############################################################################

.ina_pt_transition_representatives <- function(InitialPoints, analysis_grid,
                                               Nstages) {
  if (is.null(analysis_grid)) {
    xy <- if (nrow(InitialPoints)) c(mean(InitialPoints$x), mean(InitialPoints$y)) else c(0, 0)
    return(data.frame(
      id = seq_len(Nstages),
      x = rep(xy[1L], Nstages), y = rep(xy[2L], Nstages),
      stage = seq_len(Nstages), cell = 1L,
      stringsAsFactors = FALSE))
  }
  cc <- .ina_pt_grid_centres(analysis_grid)
  nc <- nrow(cc)
  out <- do.call(rbind, lapply(seq_len(Nstages), function(s)
    data.frame(id = (s - 1L) * nc + seq_len(nc),
               x = cc$x, y = cc$y, stage = s, cell = cc$cell,
               stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out
}

.ina_pt_transition_initial <- function(InitialPoints, analysis_grid, Nstages) {
  if ("stage" %in% names(InitialPoints)) {
    raw_stage <- InitialPoints$stage
    if (!is.numeric(raw_stage) || any(!is.finite(raw_stage)) ||
        any(raw_stage != floor(raw_stage)) || any(raw_stage < 1 | raw_stage > Nstages))
      stop("InitialPoints$stage must contain integers from 1 to Nstages.")
    stage <- as.integer(raw_stage)
  } else {
    stage <- rep(1L, nrow(InitialPoints))
  }
  if (is.null(analysis_grid)) return(tabulate(stage, nbins = Nstages))
  nc <- analysis_grid$nrow * analysis_grid$ncol
  cell <- .ina_pt_xy_to_cell(InitialPoints$x, InitialPoints$y, analysis_grid)
  idx <- (stage - 1L) * nc + cell
  tabulate(idx[!is.na(idx)], nbins = nc * Nstages)
}

.ina_pt_dynamic_operator_values <- function(G0, GH, representatives,
                                            D, InfoRetentionProb,
                                            InfoRadius, InfoTransferProb,
                                            InfoKernel, timestep,
                                            Ntimesteps) {
  nt <- nrow(G0)
  D <- .ina_pt_clip01(as.numeric(D))
  IR <- .ina_pt_point_values(InfoRetentionProb, representatives, timestep,
                             Ntimesteps, "InfoRetentionProb")
  IR <- .ina_pt_clip01(IR)
  C <- .ina_pt_direct_info_matrix(representatives, InfoRadius,
                                  InfoTransferProb, InfoKernel,
                                  timestep, Ntimesteps)
  G <- matrix(0, 2L * nt, 2L * nt)
  U <- seq_len(nt); H <- nt + seq_len(nt)
  for (i in seq_len(nt)) for (j in seq_len(nt)) {
    w <- G0[j, i]
    if (w != 0) {
      h <- D[j]
      G[U[j], U[i]] <- G[U[j], U[i]] + w * (1 - h)
      G[H[j], U[i]] <- G[H[j], U[i]] + w * h
    }
    w <- GH[j, i]
    if (w != 0) {
      h <- if (i == j) IR[j] + (1 - IR[j]) * D[j]
           else 1 - (1 - D[j]) * (1 - C[i, j])
      G[U[j], H[i]] <- G[U[j], H[i]] + w * (1 - h)
      G[H[j], H[i]] <- G[H[j], H[i]] + w * h
    }
  }
  attr(G, "note") <- paste(
    "Point-transition grid informed/uninformed low-density operator;",
    "shared persistent information sites and exact movement-information correlations are omitted."
  )
  G
}


.ina_pt_transition_dynamic_operator <- function(Parent0, ParentH, Recruit0, RecruitH,
                                                representatives, D,
                                                InfoRetentionProb, InfoRadius,
                                                InfoTransferProb, InfoKernel,
                                                timestep, Ntimesteps) {
  nt <- nrow(Parent0)
  D <- .ina_pt_clip01(as.numeric(D))
  IR <- .ina_pt_point_values(InfoRetentionProb, representatives, timestep,
                             Ntimesteps, "InfoRetentionProb")
  IR <- .ina_pt_clip01(IR)
  C <- .ina_pt_direct_info_matrix(representatives, InfoRadius,
                                  InfoTransferProb, InfoKernel,
                                  timestep, Ntimesteps)
  G <- matrix(0, 2L * nt, 2L * nt)
  U <- seq_len(nt); H <- nt + seq_len(nt)
  for (i in seq_len(nt)) for (j in seq_len(nt)) {
    # The persisting/progressing parent is the same biological individual even
    # when its stage or coordinates change. Information retention therefore
    # follows the parent rather than being treated as a source-target transfer.
    w <- Parent0[j, i]
    if (w != 0) {
      h <- D[j]
      G[U[j], U[i]] <- G[U[j], U[i]] + w * (1 - h)
      G[H[j], U[i]] <- G[H[j], U[i]] + w * h
    }
    w <- ParentH[j, i]
    if (w != 0) {
      h <- IR[j] + (1 - IR[j]) * D[j]
      G[U[j], H[i]] <- G[U[j], H[i]] + w * (1 - h)
      G[H[j], H[i]] <- G[H[j], H[i]] + w * h
    }

    # Recruits are new biological individuals. Direct information transfer is
    # approximated from the source analytical type; when the parent itself
    # moves during transition, the exact post-transition source-child distance
    # is a stochastic correlation and is reported as a limitation.
    w <- Recruit0[j, i]
    if (w != 0) {
      h <- D[j]
      G[U[j], U[i]] <- G[U[j], U[i]] + w * (1 - h)
      G[H[j], U[i]] <- G[H[j], U[i]] + w * h
    }
    w <- RecruitH[j, i]
    if (w != 0) {
      h <- 1 - (1 - D[j]) * (1 - C[i, j])
      G[U[j], H[i]] <- G[U[j], H[i]] + w * (1 - h)
      G[H[j], H[i]] <- G[H[j], H[i]] + w * h
    }
  }
  attr(G, "note") <- paste(
    "Stage x point-grid informed/uninformed operator; parent information follows stage movement;",
    "offspring information transfer uses source-type geometry and therefore approximates",
    "the correlation with a parent's realised stage-transition movement."
  )
  G
}

.ina_pt_transition_repro_contract <- function(kernels, representatives,
                                              Nstages, analysis_grid,
                                              HabitatSuitability,
                                              HabitatSearchRadius,
                                              HabitatSearchCandidates,
                                              PropaguleEstablishment,
                                              EnvEstabProb,
                                              timestep, Ntimesteps,
                                              KernelSamples, seed_base,
                                              active_stages = seq_len(Nstages)) {
  nt <- nrow(representatives)
  nc <- analysis_grid$nrow * analysis_grid$ncol
  P <- matrix(0, nt, nt)
  outside <- numeric(nt)
  stage1_types <- seq_len(nc)
  for (s in active_stages) {
    rows <- which(representatives$stage == s)
    if (!length(rows)) next
    k <- .ipptm_get_kernel(kernels, s, Nstages, "PointReproductiveKernel")
    z <- .ina_pt_kernel_contract(
      k, representatives[rows, , drop = FALSE], analysis_grid,
      HabitatSuitability, HabitatSearchRadius, HabitatSearchCandidates,
      PropaguleEstablishment, EnvEstabProb,
      timestep, Ntimesteps, KernelSamples,
      seed = seed_base + s * 997L,
      destination_stage = 1L,
      apply_habitat_search = TRUE, apply_establishment = TRUE)
    P[rows, stage1_types] <- z$internal
    outside[rows] <- z$outside
  }
  list(internal = P, outside = outside)
}

.ina_pt_transition_progress_contract <- function(TransitionKernels,
                                                 representatives,
                                                 source_stage,
                                                 Nstages, analysis_grid,
                                                 HabitatSuitability,
                                                 HabitatSearchRadius,
                                                 HabitatSearchCandidates,
                                                 TransitionHabitatSearch,
                                                 ApplyHabitatToTransitions,
                                                 TransitionEstablishment,
                                                 EnvEstabProb,
                                                 timestep, Ntimesteps,
                                                 KernelSamples, seed) {
  rows <- which(representatives$stage == source_stage)
  nc <- analysis_grid$nrow * analysis_grid$ncol
  nt <- nrow(representatives)
  P <- matrix(0, length(rows), nt)
  if (!length(rows)) return(list(rows = rows, internal = P,
                                 outside = numeric(0), fail = numeric(0)))
  k <- .ipptm_get_transition_kernel(TransitionKernels, source_stage, Nstages)
  z <- .ina_pt_kernel_contract(
    k, representatives[rows, , drop = FALSE], analysis_grid,
    HabitatSuitability, HabitatSearchRadius, HabitatSearchCandidates,
    PropaguleEstablishment = 1,
    EnvEstabProb = EnvEstabProb,
    timestep = timestep, Ntimesteps = Ntimesteps,
    KernelSamples = KernelSamples, seed = seed,
    destination_stage = source_stage + 1L,
    apply_habitat_search = isTRUE(TransitionHabitatSearch),
    apply_establishment = isTRUE(ApplyHabitatToTransitions),
    TransitionEstablishment = TransitionEstablishment)
  dest_types <- source_stage * nc + seq_len(nc)
  P[, dest_types] <- z$internal
  success <- rowSums(z$internal) + z$outside
  list(rows = rows, internal = P, outside = z$outside,
       fail = pmax(0, 1 - success))
}

.ina_pt_transition_operator_step <- function(
    timestep, Ntimesteps, representatives, analysis_grid,
    Transition, Nstages,
    SDDkernel, LDDkernel, LDDrate,
    PropaguleEstablishment, EnvEstabProb,
    TransitionKernels, TransitionHabitatSearch,
    ApplyHabitatToTransitions, TransitionEstablishment,
    BlockedTransitionMortality,
    HabitatSuitability, HabitatSearchRadius, HabitatSearchCandidates,
    DetectionProb, DetectionSpatial,
    ManageProb, ManageSpatial,
    MortalityProb, MortalitySpatial,
    FecundityReduction, FecundityReductionSpatial,
    SpreadReduction, SpreadReductionSpatial, SpreadReductionAppliesTo,
    InfoRetentionProb, InfoRadius, InfoTransferProb, InfoKernel,
    KernelSamples, PointSeed, mode, OutsideEstablishmentProb) {

  A <- .ipptm_get_transition(Transition, timestep, 1L, Nstages)
  nt <- nrow(representatives)
  nc <- if (is.null(analysis_grid)) 1L else analysis_grid$nrow * analysis_grid$ncol
  a <- .ina_pt_expected_parameter(ManageProb, ManageSpatial, representatives,
                                  timestep, 1L, Ntimesteps, "ManageProb")
  m <- .ina_pt_transition_prob_spatial_mean(MortalityProb, MortalitySpatial,
                                             representatives, timestep,
                                             Ntimesteps, Nstages,
                                             "MortalityProb")
  f_red <- .ina_pt_expected_parameter(FecundityReduction,
                                      FecundityReductionSpatial,
                                      representatives, timestep, 1L,
                                      Ntimesteps, "FecundityReduction")
  g <- .ina_pt_expected_parameter(SpreadReduction, SpreadReductionSpatial,
                                  representatives, timestep, 1L, Ntimesteps,
                                  "SpreadReduction")
  D <- .ina_pt_transition_prob_spatial_mean(DetectionProb, DetectionSpatial,
                                             representatives, timestep,
                                             Ntimesteps, Nstages,
                                             "DetectionProb")
  r <- if (length(LDDrate) == 1L) as.numeric(LDDrate) else as.numeric(LDDrate[timestep])
  if (!is.finite(r) || r < 0 || r > 1) stop("LDDrate must resolve to [0,1].")
  btm <- if (length(BlockedTransitionMortality) == 1L)
    rep(BlockedTransitionMortality, Nstages - 1L) else BlockedTransitionMortality
  if (length(btm) != Nstages - 1L) stop("BlockedTransitionMortality must be scalar or length Nstages-1.")
  btm <- .ina_pt_clip01(as.numeric(btm))

  if (is.null(analysis_grid)) {
    pe <- .ina_pt_point_values(PropaguleEstablishment, representatives,
                               timestep, Ntimesteps, "PropaguleEstablishment")
    ee <- .ina_pt_point_values(EnvEstabProb, representatives,
                               timestep, Ntimesteps, "EnvEstabProb")
    # Recruits always enter stage 1. Destination establishment parameters are
    # evaluated on a representative stage-1 point in homogeneous mode.
    rep1 <- representatives[1L, , drop = FALSE]; rep1$stage <- 1L
    pes <- .ipp_resolve(PropaguleEstablishment, rep1, timestep, 1L,
                        Ntimesteps, "PropaguleEstablishment")[1L]
    ees <- .ipp_resolve(EnvEstabProb, rep1, timestep, 1L,
                        Ntimesteps, "EnvEstabProb")[1L]
    estab <- .ina_pt_clip01(pes * ees)
    PS <- PL <- matrix(0, nt, nt)
    PS[, 1L] <- estab; PL[, 1L] <- estab
    outS <- outL <- rep(0, nt)

    progress <- vector("list", Nstages - 1L)
    for (s in seq_len(Nstages - 1L)) {
      rows <- which(representatives$stage == s)
      success <- rep(1, length(rows))
      if (ApplyHabitatToTransitions) {
        cand <- representatives[rows, c("x", "y", "stage"), drop = FALSE]
        cand$stage <- s + 1L
        pte <- .ipp_resolve(TransitionEstablishment, cand, timestep, 1L,
                            Ntimesteps, "TransitionEstablishment")
        penv <- .ipp_resolve(EnvEstabProb, cand, timestep, 1L,
                             Ntimesteps, "EnvEstabProb")
        success <- .ina_pt_clip01(pte * penv)
      }
      P <- matrix(0, length(rows), nt)
      P[cbind(seq_along(rows), rows + 1L)] <- success
      progress[[s]] <- list(rows = rows, internal = P,
                            outside = rep(0, length(rows)),
                            fail = 1 - success)
    }
  } else {
    reproductive_stages <- which(seq_len(Nstages) >= 2L & A[1L, ] > 0)
    sdd <- .ina_pt_transition_repro_contract(
      SDDkernel, representatives, Nstages, analysis_grid,
      HabitatSuitability, HabitatSearchRadius, HabitatSearchCandidates,
      PropaguleEstablishment, EnvEstabProb,
      timestep, Ntimesteps, KernelSamples,
      PointSeed + 100000L * timestep + 3000L,
      active_stages = reproductive_stages)
    PS <- sdd$internal; outS <- sdd$outside
    if (r > 0) {
      if (is.null(LDDkernel)) stop("LDDkernel is required when LDDrate > 0.")
      ldd <- .ina_pt_transition_repro_contract(
        LDDkernel, representatives, Nstages, analysis_grid,
        HabitatSuitability, HabitatSearchRadius, HabitatSearchCandidates,
        PropaguleEstablishment, EnvEstabProb,
        timestep, Ntimesteps, KernelSamples,
        PointSeed + 100000L * timestep + 6000L,
        active_stages = reproductive_stages)
      PL <- ldd$internal; outL <- ldd$outside
    } else {
      PL <- matrix(0, nt, nt); outL <- rep(0, nt)
    }
    progress <- vector("list", Nstages - 1L)
    for (s in seq_len(Nstages - 1L)) {
      progress[[s]] <- .ina_pt_transition_progress_contract(
        TransitionKernels, representatives, s, Nstages, analysis_grid,
        HabitatSuitability, HabitatSearchRadius, HabitatSearchCandidates,
        TransitionHabitatSearch, ApplyHabitatToTransitions,
        TransitionEstablishment, EnvEstabProb,
        timestep, Ntimesteps, KernelSamples,
        PointSeed + 100000L * timestep + 9000L + s * 997L)
    }
  }

  P0 <- (1 - r) * PS + r * PL
  out0 <- (1 - r) * outS + r * outL
  if (SpreadReductionAppliesTo == "LDD") {
    P1 <- (1 - r) * PS + r * sweep(PL, 1, 1 - g, `*`)
    out1 <- (1 - r) * outS + r * (1 - g) * outL
  } else {
    P1 <- sweep(P0, 1, 1 - g, `*`)
    out1 <- (1 - g) * out0
  }

  Parent0 <- Parent1 <- matrix(0, nt, nt)
  Recruit0 <- Recruit1 <- matrix(0, nt, nt)
  transition_out0 <- transition_out1 <- numeric(nt)
  branch <- vector("list", nt)

  for (i in seq_len(nt)) {
    k <- representatives$stage[i]
    fec <- if (k >= 2L) A[1L, k] else 0
    q0 <- 1
    q1 <- 1 - m[i]

    # Reproductive recruits are produced after management mortality but before
    # the stage transition. Fecundity reduction therefore acts only on fecundity.
    if (fec > 0) {
      Recruit0[, i] <- Recruit0[, i] + q0 * fec * P0[i, ]
      Recruit1[, i] <- Recruit1[, i] + q1 * fec * (1 - f_red[i]) * P1[i, ]
    }

    pself <- pnext <- pdeath <- 0
    if (k < Nstages) {
      stay <- A[k, k]
      prog <- A[k + 1L, k]
      pr <- progress[[k]]
      rr <- match(i, pr$rows)
      fail <- pr$fail[rr]
      # Failed progression remains in the source stage/location unless the
      # optional blocked-transition mortality event removes the individual.
      pself <- stay + prog * fail * (1 - btm[k])
      Parent0[i, i] <- Parent0[i, i] + q0 * pself
      Parent1[i, i] <- Parent1[i, i] + q1 * pself
      if (prog > 0) {
        Parent0[, i] <- Parent0[, i] + q0 * prog * pr$internal[rr, ]
        Parent1[, i] <- Parent1[, i] + q1 * prog * pr$internal[rr, ]
        transition_out0[i] <- q0 * prog * pr$outside[rr]
        transition_out1[i] <- q1 * prog * pr$outside[rr]
      }
      pnext <- prog * sum(pr$internal[rr, ])
      pdeath <- pmax(0, 1 - pself - pnext - prog * pr$outside[rr])
    } else {
      pself <- A[Nstages, Nstages]
      pdeath <- 1 - pself
      Parent0[i, i] <- Parent0[i, i] + q0 * pself
      Parent1[i, i] <- Parent1[i, i] + q1 * pself
    }

    branch[[i]] <- list(stage = k, a = a[i], q0 = q0, q1 = q1,
                        fec = fec, f = f_red[i],
                        mu0 = fec * rowSums(P0)[i],
                        mu1 = fec * (1 - f_red[i]) * rowSums(P1)[i],
                        pself = pself, pnext = pnext, pdeath = pdeath)
  }

  G0 <- Parent0 + Recruit0
  G1 <- Parent1 + Recruit1
  ParentH <- sweep(Parent0, 2, 1 - a, `*`) + sweep(Parent1, 2, a, `*`)
  RecruitH <- sweep(Recruit0, 2, 1 - a, `*`) + sweep(Recruit1, 2, a, `*`)
  GH <- ParentH + RecruitH
  pout <- if (length(OutsideEstablishmentProb) == 1L)
    rep(as.numeric(OutsideEstablishmentProb), nt)
  else if (length(OutsideEstablishmentProb) == nt) as.numeric(OutsideEstablishmentProb)
  else stop("OutsideEstablishmentProb must be scalar or one value per point analytical type.")
  pout <- .ina_pt_clip01(pout)

  fecvec <- ifelse(representatives$stage >= 2L,
                   A[1L, representatives$stage], 0)
  e0 <- fecvec * out0 * pout + transition_out0
  q1_export <- 1 - m
  e1 <- q1_export * fecvec * (1 - f_red) * out1 * pout + transition_out1
  eH <- (1 - a) * e0 + a * e1

  if (mode == "dynamic") {
    G <- .ina_pt_transition_dynamic_operator(
      Parent0, ParentH, Recruit0, RecruitH, representatives, D,
      InfoRetentionProb, InfoRadius, InfoTransferProb, InfoKernel,
      timestep, Ntimesteps)
    e <- c(e0, eH)
  } else if (mode == "all_informed") {
    G <- GH; e <- eH
  } else {
    G <- G0; e <- e0
  }

  # Intrinsic local demographic multiplier for each spatial cell. This strips
  # out dispersal geometry/establishment but retains expected management
  # mortality and fecundity reduction for an informed source when relevant.
  ncells <- nc
  intrinsic <- numeric(ncells)
  for (cell in seq_len(ncells)) {
    Mloc <- matrix(0, Nstages, Nstages)
    for (k in seq_len(Nstages)) {
      ii <- (k - 1L) * ncells + cell
      ai <- if (mode == "none") 0 else a[ii]
      survbar <- (1 - ai) + ai * (1 - m[ii])
      if (k < Nstages) {
        Mloc[k, k] <- A[k, k] * survbar
        Mloc[k + 1L, k] <- A[k + 1L, k] * survbar
      } else Mloc[k, k] <- A[k, k] * survbar
      if (k >= 2L) {
        fecbar <- (1 - ai) + ai * (1 - m[ii]) * (1 - f_red[ii])
        Mloc[1L, k] <- A[1L, k] * fecbar
      }
    }
    intrinsic[cell] <- .ina_pt_rho(Mloc)
  }

  list(G = G, export = e, G0 = G0, GH = GH,
       Parent0 = Parent0, ParentH = ParentH,
       Recruit0 = Recruit0, RecruitH = RecruitH,
       branch = branch, intrinsic = intrinsic, A = A,
       progress = progress)
}

.ina_pt_transition_exact_extinction <- function(step_data, mode, Nstages,
                                                 generations = 100) {
  if (!(mode %in% c("none", "all_informed"))) return(NULL)
  # Exact closed recursion is only used for homogeneous spatial structure.
  if (length(step_data[[1L]]$branch) != Nstages) return(NULL)

  stepfun <- function(st, qnext) {
    A <- st$A
    qsrc <- numeric(Nstages)
    for (k in seq_len(Nstages)) {
      br <- st$branch[[k]]
      a <- if (mode == "all_informed") br$a else 0
      calc <- function(M) {
        qkill <- if (M == 0L) br$q0 else br$q1
        mu <- if (M == 0L) br$mu0 else br$mu1
        if (k < Nstages) {
          # In homogeneous mode pnext is the successful progression probability
          # and pself contains stasis plus nonfatal blocked progression.
          local <- br$pdeath + br$pself * qnext[k] + br$pnext * qnext[k + 1L]
        } else local <- br$pdeath + br$pself * qnext[k]
        (1 - qkill) + qkill * local * exp(mu * (qnext[1L] - 1))
      }
      qsrc[k] <- (1 - a) * calc(0L) + a * calc(1L)
    }
    .ina_pt_clip01(qsrc)
  }

  qh <- rep(0, Nstages)
  for (tt in rev(seq_along(step_data))) qh <- stepfun(step_data[[tt]], qh)
  static <- length(step_data) == 1L || all(vapply(step_data[-1L], function(z)
    isTRUE(all.equal(z$A, step_data[[1L]]$A, tolerance = 0)) &&
      isTRUE(all.equal(z$branch, step_data[[1L]]$branch, tolerance = 0)),
    logical(1)))
  qe <- NULL
  if (static) {
    q <- rep(0, Nstages)
    for (gg in seq_len(generations)) {
      qo <- q; q <- stepfun(step_data[[1L]], q)
      if (max(abs(q - qo)) < 1e-12) break
    }
    qe <- q
  }
  list(horizon = qh, eventual = qe,
       method = "exact homogeneous stage-structured survival/transition + Poisson-recruit branching PGF")
}

INApestPointTransitionMatrixAnalytical <- function(
    Ntimesteps = 10,
    Nstages,
    Weights = rep(1, Nstages),
    Transition,
    InitialPoints,
    InitialInfo = NULL,
    InformationMode = c("auto", "dynamic", "all_informed", "none"),
    ApplyInitialDetection = TRUE,
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
    HabitatSuitability = NULL,
    HabitatSearchRadius = 0,
    HabitatSearchCandidates = 128,
    LocalK = Inf,
    KRadius = 0,
    DetectionProb = 0,
    DetectionSD = NULL,
    DetectionSpatial = NULL,
    ManageProb = 0,
    ManageSD = NULL,
    ManageSpatial = NULL,
    MortalityProb = 0,
    MortalitySD = NULL,
    MortalitySpatial = NULL,
    FecundityReduction = 0,
    FecundityReductionSD = NULL,
    FecundityReductionSpatial = NULL,
    SpreadReduction = 0,
    SpreadReductionSD = NULL,
    SpreadReductionSpatial = NULL,
    SpreadReductionAppliesTo = c("LDD", "all"),
    InfoRadius = 0,
    InfoTransferProb = 0,
    InfoKernel = NULL,
    InfoRetentionProb = 1,
    InfoPersistenceSteps = NA,
    ExternalInfoProb = 0,
    OngoingExternalInfo = FALSE,
    ExternalIncursionGenerator = NULL,
    OngoingExternalInvasion = FALSE,
    PointAnalysisGrid = NULL,
    MaxAnalysisCells = 400L,
    KernelSamples = 5000L,
    PointSeed = 1L,
    OutsideEstablishmentProb = 1,
    ExtinctionGenerations = 100,
    ReturnOperators = FALSE) {

  .ina_pt_require_point_helpers(TRUE)
  InformationMode <- match.arg(InformationMode)
  SpreadReductionAppliesTo <- match.arg(SpreadReductionAppliesTo)
  if (!is.numeric(Nstages) || length(Nstages) != 1L || !is.finite(Nstages) ||
      Nstages < 2L || Nstages != floor(Nstages))
    stop("Nstages must be an integer >= 2.")
  Nstages <- as.integer(Nstages)
  if (length(Weights) != Nstages || any(!is.finite(Weights)) || any(Weights <= 0))
    stop("Weights must contain Nstages finite positive values, matching INApestPointTransitionMatrix.")
  Weights <- as.numeric(Weights)
  if (!is.data.frame(InitialPoints) || !all(c("x", "y") %in% names(InitialPoints)))
    stop("InitialPoints must be a data.frame with x and y.")
  if (any(!is.finite(InitialPoints$x)) || any(!is.finite(InitialPoints$y)))
    stop("InitialPoints x and y must be finite.")
  if (!is.numeric(Ntimesteps) || length(Ntimesteps) != 1L || !is.finite(Ntimesteps) ||
      Ntimesteps < 1L || Ntimesteps != floor(Ntimesteps))
    stop("Ntimesteps must be a positive integer.")
  Ntimesteps <- as.integer(Ntimesteps)
  if (!is.numeric(KernelSamples) || length(KernelSamples) != 1L ||
      !is.finite(KernelSamples) || KernelSamples < 100L ||
      KernelSamples != floor(KernelSamples))
    stop("KernelSamples must be an integer >= 100.")
  KernelSamples <- as.integer(KernelSamples)
  if (!(is.function(SDDkernel) || (is.list(SDDkernel) && length(SDDkernel) == Nstages)))
    stop("SDDkernel must be a function or a list of length Nstages, matching INApestPointTransitionMatrix.")
  if (!is.null(LDDkernel) &&
      !(is.function(LDDkernel) || (is.list(LDDkernel) && length(LDDkernel) == Nstages)))
    stop("LDDkernel must be NULL, a function, or a list of length Nstages, matching INApestPointTransitionMatrix.")
  if (!is.null(TransitionKernels)) {
    if (!is.list(TransitionKernels) || length(TransitionKernels) != Nstages - 1L)
      stop("TransitionKernels must be NULL or a list of length Nstages - 1.")
    bad_tk <- vapply(TransitionKernels, function(k) !is.null(k) && !is.function(k), logical(1))
    if (any(bad_tk)) stop("Every non-NULL TransitionKernels entry must be a dispersal-kernel function.")
  }
  .ipp_validate_probability_schedule(LDDrate, Ntimesteps, "LDDrate")
  .ipp_validate_probability_schedule(FecundityReduction, Ntimesteps,
                                     "FecundityReduction")
  .ipp_validate_probability_schedule(FecundityReductionSD, Ntimesteps,
                                     "FecundityReductionSD", allow_null = TRUE)
  if (any(as.numeric(LDDrate) > 0) && is.null(LDDkernel))
    stop("LDDkernel is required when LDDrate > 0, matching INApestPointTransitionMatrix.")
  if (isTRUE(OngoingExternalInvasion))
    stop("OngoingExternalInvasion is additive immigration and is not represented by the current rare-lineage transition-point analytical solution; use the stochastic point simulator for that process.")
  if (isTRUE(OngoingExternalInfo) && .ina_pt_any_nonzero(ExternalInfoProb))
    stop("OngoingExternalInfo is not yet represented in the transition-point information-state operator; use the stochastic point simulator when ongoing external information is active.")

  analysis_grid <- .ina_pt_analysis_grid(
    PointAnalysisGrid, HabitatSuitability, DetectionSpatial, ManageSpatial,
    MortalitySpatial, FecundityReductionSpatial, SpreadReductionSpatial)
  if (!is.null(analysis_grid)) {
    ncells <- analysis_grid$nrow * analysis_grid$ncol
    if (is.null(PointAnalysisGrid) && ncells > as.integer(MaxAnalysisCells))
      stop("The automatically inherited point grid has ", ncells,
           " cells. Supply a coarser PointAnalysisGrid (<= MaxAnalysisCells) ",
           "rather than unintentionally contracting kernels over the full fine habitat grid.")
  }

  if (is.null(analysis_grid) && any(vapply(list(HabitatSuitability, DetectionSpatial,
                                                   ManageSpatial, MortalitySpatial,
                                                   FecundityReductionSpatial, SpreadReductionSpatial),
                                              function(z) !is.null(z), logical(1))))
    stop("Spatial habitat/management surfaces require PointAnalysisGrid unless a native INApestSpatialGrid can be inherited automatically.")

  if (is.null(analysis_grid) &&
      (!is.null(TransitionKernels) || is.function(HabitatSuitability) ||
       is.function(DetectionProb) || is.function(ManageProb) ||
       is.function(MortalityProb) || is.function(FecundityReduction) ||
       is.function(SpreadReduction))) {
    # Transition kernels themselves do not affect total homogeneous abundance
    # unless there is a spatial boundary/habitat effect. Permit them if habitat
    # is homogeneous; custom point-valued vital rates still require a grid.
    function_rate <- is.function(HabitatSuitability) || is.function(DetectionProb) ||
      is.function(ManageProb) || is.function(MortalityProb) ||
      is.function(FecundityReduction) || is.function(SpreadReduction)
    if (function_rate)
      stop("Location-dependent/function-valued transition-point parameters require PointAnalysisGrid.")
  }

  representatives <- .ina_pt_transition_representatives(InitialPoints,
                                                          analysis_grid,
                                                          Nstages)
  base_initial <- .ina_pt_transition_initial(InitialPoints, analysis_grid, Nstages)
  initial_stage <- if ("stage" %in% names(InitialPoints)) as.integer(InitialPoints$stage)
                   else rep(1L, nrow(InitialPoints))
  initial_outside <- if (is.null(analysis_grid)) 0L else {
    cell0 <- .ina_pt_xy_to_cell(InitialPoints$x, InitialPoints$y, analysis_grid)
    sum(is.na(cell0))
  }

  raw_info <- .ina_pt_initial_info_raw(InitialPoints, InitialInfo)

  info_type <- numeric(length(base_initial)); info_n <- numeric(length(base_initial))
  if (is.null(analysis_grid)) {
    for (i in seq_len(nrow(InitialPoints))) {
      idx <- initial_stage[i]
      info_type[idx] <- info_type[idx] + raw_info[i]
      info_n[idx] <- info_n[idx] + 1
    }
  } else {
    nc <- analysis_grid$nrow * analysis_grid$ncol
    cell <- .ina_pt_xy_to_cell(InitialPoints$x, InitialPoints$y, analysis_grid)
    for (i in which(!is.na(cell))) {
      idx <- (initial_stage[i] - 1L) * nc + cell[i]
      info_type[idx] <- info_type[idx] + raw_info[i]
      info_n[idx] <- info_n[idx] + 1
    }
  }
  info_type[info_n > 0] <- info_type[info_n > 0] / info_n[info_n > 0]

  mode <- .ina_pt_info_mode(InformationMode, ManageProb, info_type,
                            DetectionProb, InfoRadius, InfoTransferProb,
                            InfoKernel, InfoRetentionProb, InfoPersistenceSteps)
  persistence_profile <- if (mode == "dynamic")
    .ina_pt_persistence_profile(InfoPersistenceSteps, representatives, Ntimesteps)
  else list(requested = FALSE, values = NULL, max_age = NULL)
  programmed_dynamic <- mode == "dynamic" && isTRUE(persistence_profile$requested)
  programmed_layout <- if (programmed_dynamic)
    .ina_programmed_layout(length(base_initial), persistence_profile$max_age) else NULL
  if (is.null(analysis_grid) && mode == "dynamic" &&
      (InfoRadius > 0 || !is.null(InfoKernel)))
    stop("Distance-based information transfer in INApestPointTransitionMatrix requires PointAnalysisGrid. This is especially important when TransitionKernels move the informed parent before information transfer.")

  operators <- vector("list", Ntimesteps)
  exports <- vector("list", Ntimesteps)
  step_data <- vector("list", Ntimesteps)
  intrinsic <- vector("list", Ntimesteps)
  for (tt in seq_len(Ntimesteps)) {
    st <- .ina_pt_transition_operator_step(
      tt, Ntimesteps, representatives, analysis_grid,
      Transition, Nstages, SDDkernel, LDDkernel, LDDrate,
      PropaguleEstablishment, EnvEstabProb,
      TransitionKernels, TransitionHabitatSearch,
      ApplyHabitatToTransitions, TransitionEstablishment,
      BlockedTransitionMortality,
      HabitatSuitability, HabitatSearchRadius, HabitatSearchCandidates,
      DetectionProb, DetectionSpatial, ManageProb, ManageSpatial,
      MortalityProb, MortalitySpatial,
      FecundityReduction, FecundityReductionSpatial,
      SpreadReduction, SpreadReductionSpatial, SpreadReductionAppliesTo,
      InfoRetentionProb, InfoRadius, InfoTransferProb, InfoKernel,
      KernelSamples, PointSeed, mode, OutsideEstablishmentProb)
    if (programmed_dynamic) {
      Dtt <- .ina_pt_transition_prob_spatial_mean(
        DetectionProb, DetectionSpatial, representatives, tt, Ntimesteps,
        Nstages, "DetectionProb")
      IRtt <- .ina_pt_point_values(InfoRetentionProb, representatives, tt,
                                   Ntimesteps, "InfoRetentionProb")
      Ctt <- .ina_pt_direct_info_matrix(representatives, InfoRadius,
                                        InfoTransferProb, InfoKernel,
                                        tt, Ntimesteps)
      operators[[tt]] <- .ina_pt_programmed_operator(
        st$Parent0, st$ParentH, st$Recruit0, st$RecruitH, representatives,
        Dtt, IRtt, persistence_profile$values[[tt]], Ctt, programmed_layout)
      exports[[tt]] <- .ina_programmed_expand_vector(st$export, programmed_layout)
    } else {
      operators[[tt]] <- st$G
      exports[[tt]] <- st$export
    }
    step_data[[tt]] <- st
    intrinsic[[tt]] <- st$intrinsic
  }

  if (programmed_dynamic) {
    D1 <- .ina_pt_transition_prob_spatial_mean(
      DetectionProb, DetectionSpatial, representatives, 1L,
      Ntimesteps, Nstages, "DetectionProb")
    if (!ApplyInitialDetection) D1[] <- 0
    state0 <- .ina_programmed_initial_state(
      base_initial, .ina_pt_clip01(info_type), D1, seq_along(base_initial),
      programmed_layout, binary_known_presence = FALSE)
  } else if (mode == "dynamic") {
    # Initial detection is stage-specific in the transition-point simulator.
    D1 <- .ina_pt_transition_prob_spatial_mean(
      DetectionProb, DetectionSpatial, representatives, 1L,
      Ntimesteps, Nstages, "DetectionProb")
    info0 <- .ina_pt_clip01(info_type)
    if (ApplyInitialDetection) info0 <- info0 + (1 - info0) * D1
    state0 <- c(base_initial * (1 - info0), base_initial * info0)
  } else state0 <- base_initial

  growth <- .ina_pt_cycle_growth(operators)
  states <- .ina_pt_apply_operators(operators, state0)
  if (programmed_dynamic) totals <- colSums(states)
  else if (mode == "dynamic") {
    nt <- length(base_initial)
    totals <- colSums(states[seq_len(nt), , drop = FALSE] +
                      states[nt + seq_len(nt), , drop = FALSE])
  } else totals <- colSums(states)

  # Model-specific stage x point-grid branching PGF. The original biological
  # parent transition and the Poisson recruit distribution are retained
  # separately; spatial kernel contraction supplies the transition/recruit
  # probabilities among geographic analysis cells.
  pgf_steps <- lapply(seq_len(Ntimesteps), function(tt)
    .ina_pt_transition_branch_step_from_operator(
      step_data[[tt]], representatives, DetectionProb, DetectionSpatial,
      InfoRetentionProb, InfoRadius, InfoTransferProb, InfoKernel,
      tt, Ntimesteps))
  initial_type <- .ina_pt_initial_type_transition(InitialPoints, analysis_grid,
                                                  Nstages)
  if (programmed_dynamic) {
    qh <- .ina_programmed_branch_horizon(
      pgf_steps, persistence_profile$values, programmed_layout)
    ext_h <- .ina_programmed_initial_extinction_point(
      qh, InitialPoints, InitialInfo, initial_type, DetectionProb, DetectionSpatial,
      Ntimesteps, programmed_layout, transition = TRUE, Nstages = Nstages)
    static <- .ina_branch_static(pgf_steps) &&
      (length(persistence_profile$values) == 1L ||
       all(vapply(persistence_profile$values[-1L], function(z)
         isTRUE(all.equal(z, persistence_profile$values[[1L]], tolerance = 0)),
         logical(1))))
    qe <- if (static) .ina_programmed_branch_eventual(
      pgf_steps[[1L]], persistence_profile$values[[1L]], programmed_layout,
      ExtinctionGenerations) else NULL
    ext_e <- if (is.null(qe)) NA_real_ else .ina_programmed_initial_extinction_point(
      qe, InitialPoints, InitialInfo, initial_type, DetectionProb, DetectionSpatial,
      Ntimesteps, programmed_layout, transition = TRUE, Nstages = Nstages)
    ext_method <- if (static)
      "stage x point-grid programmed-information branching PGF" else
      "time-inhomogeneous stage x point-grid programmed-information branching PGF"
  } else {
    qh <- .ina_branch_horizon(pgf_steps, mode)
    ext_h <- .ina_pt_initial_event_probability(
      qh, InitialPoints, InitialInfo, initial_type, mode, DetectionProb,
      DetectionSpatial, Ntimesteps, transition = TRUE, Nstages = Nstages)
    static <- .ina_branch_static(pgf_steps)
    qe <- if (static) .ina_branch_eventual(pgf_steps[[1L]], mode,
                                            ExtinctionGenerations) else NULL
    ext_e <- if (is.null(qe)) NA_real_ else .ina_pt_initial_event_probability(
      qe, InitialPoints, InitialInfo, initial_type, mode, DetectionProb,
      DetectionSpatial, Ntimesteps, transition = TRUE, Nstages = Nstages)
    ext_method <- if (mode == "dynamic") {
      if (static) "stage x point-grid informed/uninformed branching PGF" else
        "time-inhomogeneous stage x point-grid informed/uninformed branching PGF"
    } else {
      if (static) "stage x point-grid survival/transition + Poisson-recruit branching PGF" else
        "time-inhomogeneous stage x point-grid survival/transition + Poisson-recruit branching PGF"
    }
  }

  escape <- NULL
  if (!is.null(analysis_grid)) {
    first <- .ina_pt_first_moment_escape(operators, state0, exports)
    ng <- .ina_pt_poisson_noescape_horizon(operators, exports)
    escape <- list(
      Method = paste(
        "stage x point-grid mean-matched branching no-escape recursion;",
        "includes reproductive export and successful stage-transition movement outside the analysis grid"),
      BranchingProbabilityByHorizon = 1 - .ina_pt_overall_event(ng, state0),
      FirstMomentPoissonProbabilityByHorizon = tail(first$poisson_escape_probability, 1L),
      ExpectedSuccessfulEscapesByStep = first$expected_successful_escapes,
      ExpectedSuccessfulEscapesCumulative = first$cumulative_expected_successful_escapes)
  }

  diagnostics <- character(0)
  if (!is.null(analysis_grid)) diagnostics <- c(diagnostics,
    paste0("Reproductive and stage-transition kernels were contracted to ",
           analysis_grid$nrow, " x ", analysis_grid$ncol,
           " analysis cells using ", KernelSamples,
           " Monte Carlo draws per source type/kernel/timestep."))
  else diagnostics <- c(diagnostics,
    "Spatially homogeneous transition-point solution: kernel distances do not change total abundance/extinction unless movement interacts with habitat, a boundary, or blocked transition success.")
  if (!is.infinite(LocalK) || KRadius > 0) diagnostics <- c(diagnostics,
    "Finite LocalK/KRadius is not represented in the rare-lineage operator; it becomes a local-density interaction once points accumulate.")
  if (programmed_dynamic) diagnostics <- c(diagnostics,
    "InfoPersistenceSteps is represented with explicit time-since-local-evidence states that follow the biological parent through stage movement. Inf is treated as a programmed window that does not expire.")
  if (mode == "dynamic") diagnostics <- c(diagnostics,
    "Information transfer after stage movement is approximated at analytical-grid type resolution; the exact correlation between realised movement, persistent information sites and later information transfer remains stochastic.")
  if (!is.null(analysis_grid)) diagnostics <- c(diagnostics,
    "For spatial transition-point extinction, lineages that move outside PointAnalysisGrid are outside the internal state. Use Escape alongside Extinction and choose an analysis grid covering the biological domain of interest.")
  explicit_sd <- any(vapply(list(DetectionSD, ManageSD, MortalitySD,
                                   FecundityReductionSD, SpreadReductionSD),
                              .ina_pt_sd_nonzero, logical(1)))
  implicit_sd <- .ina_pt_implicit_sd_names(
    DetectionProb, DetectionSD, ManageProb, ManageSD,
    MortalityProb, MortalitySD, FecundityReduction, FecundityReductionSD,
    SpreadReduction, SpreadReductionSD)
  if (explicit_sd || length(implicit_sd)) diagnostics <- c(diagnostics,
    "Parameter SDs are not integrated in the analytical solution; nominal mean parameter values are used.")
  if (length(implicit_sd)) diagnostics <- c(diagnostics,
    paste0("For direct comparison with INApestPointTransitionMatrix stochastic runs, note that NULL ",
           paste(implicit_sd, collapse = ", "),
           " activates simulator default variation. Set the relevant SD argument(s) explicitly to 0 when comparing to the nominal-mean analytical solution."))
  if (any(Weights != 1)) diagnostics <- c(diagnostics,
    "Stage Weights affect weighted output summaries and finite LocalK accounting in the simulator, but not the rare-lineage count operator reported here; expected_state_total is an individual/colony count, not a weighted total.")
  if (initial_outside > 0) diagnostics <- c(diagnostics,
    paste(initial_outside, "initial point(s) fall outside PointAnalysisGrid and are omitted from the internal analytical state."))

  intrinsic_out <- if (length(intrinsic) == 1L) {
    as.numeric(intrinsic[[1L]])
  } else {
    z <- do.call(cbind, intrinsic)
    colnames(z) <- paste0("timestep", seq_along(intrinsic))
    z
  }
  result <- list(
    Model = "INApestPointTransitionMatrix",
    InformationMode = mode,
    InformationStateMethod = if (programmed_dynamic)
      "explicit time-since-local-evidence stage x point-grid states" else if (mode == "dynamic")
      "uninformed/informed stage x point-grid states" else mode,
    Growth = list(
      EquivalentPerTimestepMultiplier = growth$EquivalentPerTimestepMultiplier,
      CycleMultiplier = growth$CycleMultiplier,
      Classification = .ina_pt_classify(growth$EquivalentPerTimestepMultiplier),
      IntrinsicLocalLambda = intrinsic_out,
      IntrinsicLocalLambdaInterpretation = paste(
        "Local stage-demographic multipliers ignoring dispersal geometry;",
        "management mortality and fecundity reduction are retained according to the selected information mode.",
        "For multi-timestep analyses the result is cells x timesteps.")),
    Trajectory = data.frame(timestep = seq_len(Ntimesteps),
                            expected_state_total = totals),
    Extinction = list(Method = ext_method,
                      ProbabilityByHorizon = ext_h,
                      BranchingFadeoutProbability = ext_e),
    Escape = escape,
    Diagnostics = unique(diagnostics),
    ApproximationScope = c(
      "Low-density / rare-invasion screening unless otherwise stated",
      "Homogeneous stage branching is exact for the represented survival/transition and Poisson reproductive processes",
      "Spatial point solutions use Monte Carlo contraction of continuous kernels to an analysis grid",
      "Finite local capacity and persistent shared information become increasingly important away from rarity"),
    PointApproximation = list(
      SpatialMode = if (is.null(analysis_grid)) "homogeneous" else "analysis-grid contraction",
      KernelSamples = if (is.null(analysis_grid)) NA_integer_ else KernelSamples,
      AnalysisGrid = analysis_grid)
  )
  if (ReturnOperators) result$Operators <- operators
  class(result) <- "INApestAnalyticalResult"
  result
}

###############################################################################
### Development dispatcher for the two point families.
### The final merge should add these Model values directly to INApestAnalytical().
###############################################################################

INApestPointAnalytical <- function(
    Model = c("INApestMetaPoint", "INApestPointTransitionMatrix"), ...) {
  Model <- match.arg(Model)
  if (Model == "INApestMetaPoint") INApestMetaPointAnalytical(...)
  else INApestPointTransitionMatrixAnalytical(...)
}

INApestAnalytical_round2_core <- function(
    Model = c("INApest", "INApestMeta", "INApestMetaTransitionMatrix",
              "INApestMetaMultipleLandUse"),
    Ntimesteps = 10,
    InitialState,
    InitialInfo = 0,
    InformationMode = c("auto", "dynamic", "all_informed", "none"),
    ApplyInitialDetection = TRUE,
    SDDprob,
    LDDprob = 0,
    LDDrate = 0,
    EnvEstabProb = 1,
    Survival = 1,
    K = NULL,
    PropaguleProduction = NULL,
    PropaguleEstablishment = 1,
    Transition = NULL,
    Nstages = NULL,
    SeedbankK = NULL,
    DetectionProb = 0,
    ManageProb = 0,
    EradicationProb = 0,
    MortalityProb = 0,
    SpreadReduction = 0,
    SEAM = NULL,
    InfoRetentionProb = 1,
    InfoPersistenceSteps = NA,
    DispersalDensityFactor = 0,
    ExportProb = NULL,
    ExportSDDprob = NULL,
    ExportLDDprob = NULL,
    OutsideEstablishmentProb = 1,
    AssumeResidualExport = FALSE,
    ExtinctionGenerations = 100,
    ReturnOperators = FALSE,
    FecundityReduction = 0) {

  Model <- match.arg(Model)
  InformationMode <- match.arg(InformationMode)
  if (length(Ntimesteps) != 1L || Ntimesteps < 1 || Ntimesteps != as.integer(Ntimesteps))
    stop("Ntimesteps must be a positive integer")
  Ntimesteps <- as.integer(Ntimesteps)
  LDDprob <- .ina_normalize_ldd(LDDprob)
  SDD0 <- if (length(dim(SDDprob)) == 3L) SDDprob[, , 1] else as.matrix(SDDprob)
  n <- nrow(SDD0)
  if (ncol(SDD0) != n) stop("SDDprob must be square")
  diagnostics <- character(0)
  # Validate the programmed-persistence parameterisation up front.  A finite
  # value activates explicit time-since-local-evidence state classes; NA leaves
  # that node/timestep on the stochastic InfoRetentionProb pathway.
  .ina_slice_node(InfoPersistenceSteps, 1, n, Ntimesteps, "InfoPersistenceSteps")
  persistence_max_age <- .ina_programmed_global_max_age(InfoPersistenceSteps)
  persistence_requested <- !is.null(persistence_max_age)
  if (persistence_requested) {
    diagnostics <- c(diagnostics,
      "InfoPersistenceSteps is represented with explicit time-since-last-local-evidence states. Programmed stopping takes priority over InfoRetentionProb wherever the persistence value is finite.")
    if (any(as.numeric(InfoRetentionProb) < 1, na.rm = TRUE))
      diagnostics <- c(diagnostics,
        "Both InfoPersistenceSteps and InfoRetentionProb are supplied: programmed stopping has priority for finite InfoPersistenceSteps; stochastic retention applies only where InfoPersistenceSteps is NA.")
  }
  if (any(as.numeric(InitialState) %% 1 != 0, na.rm = TRUE)) {
    diagnostics <- c(diagnostics,
      "InitialState contains fractional values. Branching event probabilities treat these as continuous lineage weights and are therefore approximate.")
  }

  mode <- .ina_info_mode(InformationMode, ManageProb, InitialInfo,
                         DetectionProb, SEAM, InfoRetentionProb,
                         InfoPersistenceSteps)
  dynamic <- identical(mode, "dynamic")
  programmed_dynamic <- dynamic && persistence_requested
  information_state_method <- if (programmed_dynamic)
    "explicit time-since-local-evidence states" else if (dynamic)
    "uninformed/informed memoryless states" else mode
  programmed_layout <- NULL

  operators <- vector("list", Ntimesteps)
  export_vectors <- vector("list", Ntimesteps)
  intrinsic <- NULL
  type_count <- NULL
  base_initial <- NULL

  if (Model == "INApest") {
    base_initial <- as.numeric(InitialState)
    if (length(base_initial) != n) stop("InitialState must have length nodes for INApest")
    d1 <- .ina_slice_node(DetectionProb, 1, n, Ntimesteps, "DetectionProb")
    pdet0 <- if (ApplyInitialDetection && dynamic)
      .ina_initial_detection_binary(base_initial, .inapest_recycle(d1, n, "DetectionProb")) else rep(0, n)
    if (programmed_dynamic) {
      programmed_layout <- .ina_programmed_layout(n, persistence_max_age)
      state0 <- .ina_programmed_initial_state(base_initial, InitialInfo, pdet0,
                                               seq_len(n), programmed_layout,
                                               binary_known_presence = TRUE)
    } else {
      info0 <- .inapest_recycle(InitialInfo, n, "InitialInfo")
      if (ApplyInitialDetection && dynamic) info0 <- info0 + (1 - info0) * pdet0
      info0 <- pmin(1, pmax(0, info0))
      state0 <- if (dynamic) c(base_initial * (1 - info0), base_initial * info0) else base_initial
    }
    type_count <- length(state0)

    for (tt in seq_len(Ntimesteps)) {
      SDDt <- .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob")
      LDDt <- .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob")
      Et <- .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb")
      St <- .ina_slice_node(Survival, tt, n, Ntimesteps, "Survival")
      Dt <- .ina_slice_node(DetectionProb, tt, n, Ntimesteps, "DetectionProb")
      At <- .ina_slice_node(ManageProb, tt, n, Ntimesteps, "ManageProb")
      Kt <- .ina_slice_node(EradicationProb, tt, n, Ntimesteps, "EradicationProb")
      Rt <- .ina_slice_node(SpreadReduction, tt, n, Ntimesteps, "SpreadReduction")
      IRt <- .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb")
      IPt <- .ina_slice_node(InfoPersistenceSteps, tt, n, Ntimesteps, "InfoPersistenceSteps")
      if (programmed_dynamic) {
        G0 <- inapest_exogenous_operator(SDDt, LDDt, Et, St, 0, 0, 0)
        GH <- inapest_exogenous_operator(SDDt, LDDt, Et, St, At, Kt, Rt)
        operators[[tt]] <- .ina_programmed_information_operator(
          G0, GH, seq_len(n), Dt, SEAM, IRt, IPt, programmed_layout,
          binary_known_presence = TRUE)
      } else if (mode == "dynamic")
        operators[[tt]] <- inapest_detection_operator(SDDt, LDDt, Et, St, Dt, At, Kt, Rt, SEAM, IRt)
      else if (mode == "all_informed")
        operators[[tt]] <- inapest_exogenous_operator(SDDt, LDDt, Et, St, At, Kt, Rt)
      else
        operators[[tt]] <- inapest_exogenous_operator(SDDt, LDDt, Et, St, 0, 0, 0)

      if (!is.null(ExportProb)) {
        X <- .ina_slice_export(ExportProb, tt, n, Ntimesteps, "ExportProb")
        pout <- .inapest_recycle(OutsideEstablishmentProb, n, "OutsideEstablishmentProb")
        Xeff <- sweep(X, 1, pout, `*`)
        if (mode == "dynamic") {
          e <- .inapest_export_vector(Xeff, St, At, Kt, Rt, TRUE)
          export_vectors[[tt]] <- if (programmed_dynamic) .ina_programmed_expand_vector(e, programmed_layout) else e
        } else if (mode == "all_informed") export_vectors[[tt]] <- .inapest_export_vector(Xeff, St, At, Kt, Rt, FALSE)
        else export_vectors[[tt]] <- .inapest_export_vector(Xeff, St, 0, 0, 0, FALSE)
      }
    }
  }

  if (Model == "INApestMeta") {
    if (is.null(K) || is.null(PropaguleProduction)) stop("K and PropaguleProduction are required for INApestMeta")
    base_initial <- as.numeric(InitialState)
    if (length(base_initial) != n) stop("InitialState must have length nodes for INApestMeta")
    d1 <- .ina_slice_node(DetectionProb, 1, n, Ntimesteps, "DetectionProb")
    pdet0 <- if (ApplyInitialDetection && dynamic)
      .ina_initial_detection_meta(base_initial, .ina_recycle(d1, n, "DetectionProb")) else rep(0, n)
    if (programmed_dynamic) {
      programmed_layout <- .ina_programmed_layout(n, persistence_max_age)
      state0 <- .ina_programmed_initial_state(base_initial, InitialInfo, pdet0,
                                               seq_len(n), programmed_layout,
                                               binary_known_presence = FALSE)
    } else {
      info0 <- .ina_recycle(InitialInfo, n, "InitialInfo")
      if (ApplyInitialDetection && dynamic) info0 <- info0 + (1 - info0) * pdet0
      info0 <- pmin(1, pmax(0, info0))
      state0 <- if (dynamic) c(base_initial * (1 - info0), base_initial * info0) else base_initial
    }
    type_count <- length(state0)

    for (tt in seq_len(Ntimesteps)) {
      SDDt <- .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob")
      LDDt <- .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob")
      Et <- .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb")
      St <- .ina_slice_node(Survival, tt, n, Ntimesteps, "Survival")
      Kt <- .ina_slice_node(K, tt, n, Ntimesteps, "K")
      Pt <- .ina_slice_node(PropaguleProduction, tt, n, Ntimesteps, "PropaguleProduction")
      PEt <- .ina_slice_node(PropaguleEstablishment, tt, n, Ntimesteps, "PropaguleEstablishment")
      Dt <- .ina_slice_node(DetectionProb, tt, n, Ntimesteps, "DetectionProb")
      At <- .ina_slice_node(ManageProb, tt, n, Ntimesteps, "ManageProb")
      Mt <- .ina_slice_node(MortalityProb, tt, n, Ntimesteps, "MortalityProb")
      Rt <- .ina_slice_node(SpreadReduction, tt, n, Ntimesteps, "SpreadReduction")
      Ft <- .ina_slice_node(FecundityReduction, tt, n, Ntimesteps, "FecundityReduction")
      IRt <- .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb")
      IPt <- .ina_slice_node(InfoPersistenceSteps, tt, n, Ntimesteps, "InfoPersistenceSteps")
      if (programmed_dynamic) {
        G0 <- meta_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, 0, 0, 0, FecundityReduction = 0)
        GH <- meta_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, At, Mt, Rt, FecundityReduction = Ft)
        operators[[tt]] <- .ina_programmed_information_operator(
          G0, GH, seq_len(n), Dt, SEAM, IRt, IPt, programmed_layout,
          binary_known_presence = FALSE)
      } else if (mode == "dynamic")
        operators[[tt]] <- meta_detection_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, Dt, At, Mt, Rt, SEAM, IRt, FecundityReduction = Ft)
      else if (mode == "all_informed")
        operators[[tt]] <- meta_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, At, Mt, Rt, FecundityReduction = Ft)
      else
        operators[[tt]] <- meta_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, 0, 0, 0, FecundityReduction = 0)

      if (.ina_export_available(Model, NULL, ExportSDDprob, ExportLDDprob, AssumeResidualExport)) {
        ESDDt <- .ina_slice_export(ExportSDDprob, tt, n, Ntimesteps, "ExportSDDprob")
        ELDDt <- .ina_slice_export(ExportLDDprob, tt, n, Ntimesteps, "ExportLDDprob")
        if (mode == "dynamic") {
          e <- meta_detection_export_vector(SDDt, LDDt, LDDrate, St, Pt, At, Mt, Rt, ESDDt, ELDDt, FecundityReduction = Ft)
          if (programmed_dynamic) e <- .ina_programmed_expand_vector(e, programmed_layout)
        } else if (mode == "all_informed") e <- meta_export_vector(SDDt, LDDt, LDDrate, St, Pt, At, Mt, Rt, ESDDt, ELDDt, FecundityReduction = Ft)
        else e <- meta_export_vector(SDDt, LDDt, LDDrate, St, Pt, 0, 0, 0,
                                     ExportSDDprob = ESDDt, ExportLDDprob = ELDDt,
                                     FecundityReduction = 0)
        pout_node <- .ina_recycle(OutsideEstablishmentProb, n, "OutsideEstablishmentProb")
        if (programmed_dynamic) {
          pout2 <- rep(pout_node, 2)
          pout <- .ina_programmed_expand_vector(pout2, programmed_layout)
        } else if (dynamic) pout <- rep(pout_node, 2) else pout <- pout_node
        export_vectors[[tt]] <- e * pout
      }
    }
  }

  if (Model == "INApestMetaTransitionMatrix") {
    if (is.null(Transition)) stop("Transition is required for INApestMetaTransitionMatrix")
    if (is.null(Nstages)) Nstages <- if (is.list(Transition)) nrow(as.matrix(Transition[[1]])) else nrow(as.matrix(Transition))
    S <- as.integer(Nstages)
    if (is.null(K)) K <- 1
    if (is.null(SeedbankK)) SeedbankK <- K
    X0 <- if (is.matrix(InitialState)) as.matrix(InitialState) else matrix(as.numeric(InitialState), nrow = n, ncol = S, byrow = TRUE)
    if (!all(dim(X0) == c(n, S))) stop("InitialState must be nodes x stages for INApestMetaTransitionMatrix")
    base_initial <- as.vector(t(X0))
    D1 <- .ina_slice_stage(DetectionProb, 1, n, S, Ntimesteps, "DetectionProb")
    if (length(D1) == 1L) D1m <- matrix(D1, n, S) else if (length(D1) == S) D1m <- matrix(rep(D1, each = n), n, S) else D1m <- as.matrix(D1)
    pdet0 <- if (ApplyInitialDetection && dynamic) .ina_initial_detection_transition(X0, D1m) else rep(0, n)
    type_node <- rep(seq_len(n), each = S)
    if (programmed_dynamic) {
      programmed_layout <- .ina_programmed_layout(n * S, persistence_max_age)
      state0 <- .ina_programmed_initial_state(base_initial, InitialInfo, pdet0,
                                               type_node, programmed_layout,
                                               binary_known_presence = FALSE)
    } else {
      info0 <- .ina_recycle(InitialInfo, n, "InitialInfo")
      if (ApplyInitialDetection && dynamic) info0 <- info0 + (1 - info0) * pdet0
      info0 <- pmin(1, pmax(0, info0))
      state0 <- if (dynamic) .ina_split_information(base_initial, info0, n, S) else base_initial
    }
    type_count <- length(state0)

    if (!dynamic && length(operators) && length(dim(SDDprob)) != 3L && !is.list(Transition) &&
        !.ina_fecundity_transition_is_temporal(FecundityReduction, n, S, Ntimesteps)) {
      # Filled below after mortality/fecundity slices are available; retained as node-level vector.
      intrinsic <- numeric(n)
    }

    for (tt in seq_len(Ntimesteps)) {
      SDDt <- .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob")
      LDDt <- .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob")
      Et <- .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb")
      PEt <- .ina_slice_node(PropaguleEstablishment, tt, n, Ntimesteps, "PropaguleEstablishment")
      Dt <- .ina_slice_stage(DetectionProb, tt, n, S, Ntimesteps, "DetectionProb")
      At <- .ina_slice_node(ManageProb, tt, n, Ntimesteps, "ManageProb")
      Mt <- .ina_slice_stage(MortalityProb, tt, n, S, Ntimesteps, "MortalityProb")
      Rt <- .ina_slice_node(SpreadReduction, tt, n, Ntimesteps, "SpreadReduction")
      Ft <- .ina_slice_fecundity_transition(FecundityReduction, tt, n, S, Ntimesteps, "FecundityReduction")
      IRt <- .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb")
      IPt <- .ina_slice_node(InfoPersistenceSteps, tt, n, Ntimesteps, "InfoPersistenceSteps")
      Kt <- .ina_slice_node(K, tt, n, Ntimesteps, "K")
      SBt <- .ina_slice_node(SeedbankK, tt, n, Ntimesteps, "SeedbankK")
      if (programmed_dynamic) {
        G0 <- transition_operator(Transition, S, SDDt, LDDt, LDDrate, Et, PEt, 0, 0, 0, DispersalDensityFactor, Kt, SBt, FecundityReduction = 0)
        GH <- transition_operator(Transition, S, SDDt, LDDt, LDDrate, Et, PEt, At, Mt, Rt, DispersalDensityFactor, Kt, SBt, FecundityReduction = Ft)
        Dm <- if (length(Dt) == 1L) matrix(Dt, n, S) else if (length(Dt) == S) matrix(rep(Dt, each = n), n, S) else as.matrix(Dt)
        operators[[tt]] <- .ina_programmed_information_operator(
          G0, GH, type_node, as.vector(t(Dm)), SEAM, IRt, IPt, programmed_layout,
          binary_known_presence = FALSE)
      } else if (mode == "dynamic")
        operators[[tt]] <- transition_detection_operator(Transition, S, SDDt, LDDt, LDDrate, Et, PEt, Dt, At, Mt, Rt, SEAM, IRt, DispersalDensityFactor, Kt, SBt, FecundityReduction = Ft)
      else if (mode == "all_informed")
        operators[[tt]] <- transition_operator(Transition, S, SDDt, LDDt, LDDrate, Et, PEt, At, Mt, Rt, DispersalDensityFactor, Kt, SBt, FecundityReduction = Ft)
      else
        operators[[tt]] <- transition_operator(Transition, S, SDDt, LDDt, LDDrate, Et, PEt, 0, 0, 0, DispersalDensityFactor, Kt, SBt, FecundityReduction = 0)

      if (.ina_export_available(Model, NULL, ExportSDDprob, ExportLDDprob, TRUE)) {
        ESDDt <- .ina_slice_export(ExportSDDprob, tt, n, Ntimesteps, "ExportSDDprob")
        ELDDt <- .ina_slice_export(ExportLDDprob, tt, n, Ntimesteps, "ExportLDDprob")
        if (mode == "dynamic") {
          e <- transition_detection_export_vector(Transition, S, SDDt, LDDt, LDDrate, At, Mt, Rt,
                                                    DispersalDensityFactor, ESDDt, ELDDt,
                                                    FecundityReduction = Ft)
          if (programmed_dynamic) e <- .ina_programmed_expand_vector(e, programmed_layout)
        } else if (mode == "all_informed") e <- transition_export_vector(Transition, S, SDDt, LDDt, LDDrate, At, Mt, Rt,
                                                                           DispersalDensityFactor, ESDDt, ELDDt,
                                                                           FecundityReduction = Ft)
        else e <- transition_export_vector(Transition, S, SDDt, LDDt, LDDrate, 0, 0, 0,
                                            DispersalDensityFactor, ESDDt, ELDDt,
                                            FecundityReduction = 0)
        pout_base <- rep(.ina_recycle(OutsideEstablishmentProb, n, "OutsideEstablishmentProb"), each = S)
        if (programmed_dynamic) {
          pout <- .ina_programmed_expand_vector(c(pout_base, pout_base), programmed_layout)
        } else if (dynamic) pout <- c(pout_base, pout_base) else pout <- pout_base
        export_vectors[[tt]] <- e * pout
      }
    }

    if (!is.list(Transition) && length(dim(Transition)) <= 2L &&
        !.ina_is_temporal(MortalityProb, n, Ntimesteps, "stage") &&
        !.ina_is_temporal(ManageProb, n, Ntimesteps, "node") &&
        !.ina_fecundity_transition_is_temporal(FecundityReduction, n, S, Ntimesteps)) {
      A <- as.matrix(Transition)
      M0 <- .ina_slice_stage(MortalityProb, 1, n, S, Ntimesteps, "MortalityProb")
      if (length(M0) == 1L) Mmat <- matrix(M0, n, S) else if (length(M0) == S) Mmat <- matrix(rep(M0, each = n), n, S) else Mmat <- as.matrix(M0)
      Fmat <- .ina_slice_fecundity_transition(FecundityReduction, 1, n, S, Ntimesteps, "FecundityReduction")
      Avec <- .ina_recycle(.ina_slice_node(ManageProb, 1, n, Ntimesteps, "ManageProb"), n, "ManageProb")
      if (mode == "none") Avec[] <- 0
      intrinsic <- vapply(seq_len(n), function(i) {
        Mloc <- matrix(0, S, S)
        for (k in seq_len(S)) {
          survbar <- (1 - Avec[i]) + Avec[i] * (1 - Mmat[i, k])
          if (k < S) {
            Mloc[k, k] <- A[k, k] * survbar
            Mloc[k + 1L, k] <- A[k + 1L, k] * survbar
          } else Mloc[S, S] <- A[S, S] * survbar
          if (k >= 2L)
            Mloc[1L, k] <- A[1L, k] * ((1 - Avec[i]) +
              Avec[i] * (1 - Mmat[i, k]) * (1 - Fmat[i, k]))
        }
        max(Mod(eigen(Mloc, only.values = TRUE)$values))
      }, numeric(1))
    }
  }

  if (Model == "INApestMetaMultipleLandUse") {
    if (is.null(K) || is.null(PropaguleProduction)) stop("K and PropaguleProduction are required for INApestMetaMultipleLandUse")
    K0 <- if (length(dim(K)) == 3L) K[, , 1] else as.matrix(K); L <- ncol(K0)
    X0 <- if (is.matrix(InitialState)) as.matrix(InitialState) else matrix(as.numeric(InitialState), nrow = n, ncol = L, byrow = TRUE)
    if (!all(dim(X0) == c(n, L))) stop("InitialState must be nodes x land uses")
    base_initial <- as.vector(t(X0))
    D1 <- .ina_slice_mlu(DetectionProb, 1, n, L, Ntimesteps, "DetectionProb")
    if (length(D1) == 1L) D1m <- matrix(D1, n, L) else if (length(D1) == L) D1m <- matrix(rep(D1, each = n), n, L) else D1m <- as.matrix(D1)
    pdet0 <- if (ApplyInitialDetection && dynamic) .ina_initial_detection_mlu(X0, D1m) else rep(0, n)
    type_node <- rep(seq_len(n), each = L)
    if (programmed_dynamic) {
      programmed_layout <- .ina_programmed_layout(n * L, persistence_max_age)
      state0 <- .ina_programmed_initial_state(base_initial, InitialInfo, pdet0,
                                               type_node, programmed_layout,
                                               binary_known_presence = FALSE)
    } else {
      info0 <- .ina_recycle(InitialInfo, n, "InitialInfo")
      if (ApplyInitialDetection && dynamic) info0 <- info0 + (1 - info0) * pdet0
      info0 <- pmin(1, pmax(0, info0))
      state0 <- if (dynamic) .ina_split_information(base_initial, info0, n, L) else base_initial
    }
    type_count <- length(state0)

    for (tt in seq_len(Ntimesteps)) {
      SDDt <- .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob")
      LDDt <- .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob")
      Et <- .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb")
      St <- .ina_slice_node(Survival, tt, n, Ntimesteps, "Survival")
      Kt <- .ina_slice_K_mlu(K, tt, n, L, Ntimesteps)
      Pt <- .ina_slice_node(PropaguleProduction, tt, n, Ntimesteps, "PropaguleProduction")
      PEt <- .ina_slice_node(PropaguleEstablishment, tt, n, Ntimesteps, "PropaguleEstablishment")
      Dt <- .ina_slice_mlu(DetectionProb, tt, n, L, Ntimesteps, "DetectionProb")
      At <- .ina_slice_mlu(ManageProb, tt, n, L, Ntimesteps, "ManageProb")
      Mt <- .ina_slice_mlu(MortalityProb, tt, n, L, Ntimesteps, "MortalityProb")
      Rt <- .ina_slice_mlu(SpreadReduction, tt, n, L, Ntimesteps, "SpreadReduction")
      Ft <- .ina_slice_mlu(FecundityReduction, tt, n, L, Ntimesteps, "FecundityReduction")
      IRt <- .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb")
      IPt <- .ina_slice_node(InfoPersistenceSteps, tt, n, Ntimesteps, "InfoPersistenceSteps")
      if (programmed_dynamic) {
        G0 <- mlu_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, 0, 0, 0, FALSE, FecundityReduction = 0)
        GH <- mlu_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, At, Mt, Rt, FALSE, FecundityReduction = Ft)
        Dm <- if (length(Dt) == 1L) matrix(Dt, n, L) else if (length(Dt) == L) matrix(rep(Dt, each = n), n, L) else as.matrix(Dt)
        operators[[tt]] <- .ina_programmed_information_operator(
          G0, GH, type_node, as.vector(t(Dm)), SEAM, IRt, IPt, programmed_layout,
          binary_known_presence = FALSE)
      } else if (mode == "dynamic")
        operators[[tt]] <- mlu_detection_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, Dt, At, Mt, Rt, SEAM, IRt, FALSE, FecundityReduction = Ft)
      else if (mode == "all_informed")
        operators[[tt]] <- mlu_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, At, Mt, Rt, FALSE, FecundityReduction = Ft)
      else
        operators[[tt]] <- mlu_single_parent_operator(SDDt, LDDt, LDDrate, Et, St, Kt, Pt, PEt, 0, 0, 0, FALSE, FecundityReduction = 0)

      if (.ina_export_available(Model, NULL, ExportSDDprob, ExportLDDprob, AssumeResidualExport)) {
        ESDDt <- .ina_slice_export(ExportSDDprob, tt, n, Ntimesteps, "ExportSDDprob")
        ELDDt <- .ina_slice_export(ExportLDDprob, tt, n, Ntimesteps, "ExportLDDprob")
        if (mode == "dynamic") {
          e <- mlu_detection_export_vector(SDDt, LDDt, LDDrate, St, Kt, Pt, At, Mt, Rt, FALSE,
                                            ESDDt, ELDDt, FecundityReduction = Ft)
          if (programmed_dynamic) e <- .ina_programmed_expand_vector(e, programmed_layout)
        } else if (mode == "all_informed") e <- mlu_export_vector(SDDt, LDDt, LDDrate, St, Kt, Pt, At, Mt, Rt, FALSE,
                                                                    ESDDt, ELDDt, FecundityReduction = Ft)
        else e <- mlu_export_vector(SDDt, LDDt, LDDrate, St, Kt, Pt, 0, 0, 0, FALSE,
                                    ESDDt, ELDDt, FecundityReduction = 0)
        pout_base <- rep(.ina_recycle(OutsideEstablishmentProb, n, "OutsideEstablishmentProb"), each = L)
        if (programmed_dynamic) {
          pout <- .ina_programmed_expand_vector(c(pout_base, pout_base), programmed_layout)
        } else if (dynamic) pout <- rep(pout_base, 2) else pout <- pout_base
        export_vectors[[tt]] <- e * pout
      }
    }
  }

  growth <- .ina_cycle_growth(operators)
  Rstep <- growth$EquivalentPerTimestepMultiplier
  trajectory_state <- .ina_apply_operators(operators, state0)
  totals <- if (programmed_dynamic) colSums(trajectory_state) else .ina_state_total(trajectory_state, dynamic)
  trajectory <- data.frame(timestep = seq_len(Ntimesteps), expected_state_total = totals)

  # Extinction: use the most simulator-specific PGF available for static simple
  # cases; otherwise use a transparent mean-matched multitype branching fallback.
  static_ops <- all(vapply(operators[-1], function(x) isTRUE(all.equal(x, operators[[1]], tolerance = 0)), logical(1)))
  if (Ntimesteps == 1L) static_ops <- TRUE
  extinction_method <- NULL; qh <- NULL; qe <- NULL
  if (static_ops && !dynamic && Model == "INApest") {
    pars <- list(SDDprob = .ina_slice_connection(SDDprob, 1, n, Ntimesteps, "SDDprob"),
                 LDDprob = .ina_slice_connection(LDDprob, 1, n, Ntimesteps, "LDDprob"),
                 EnvEstabProb = .ina_slice_node(EnvEstabProb, 1, n, Ntimesteps, "EnvEstabProb"),
                 Survival = .ina_slice_node(Survival, 1, n, Ntimesteps, "Survival"),
                 ManageProb = if (mode == "all_informed") .ina_slice_node(ManageProb, 1, n, Ntimesteps, "ManageProb") else 0,
                 EradicationProb = if (mode == "all_informed") .ina_slice_node(EradicationProb, 1, n, Ntimesteps, "EradicationProb") else 0,
                 SpreadReduction = if (mode == "all_informed") .ina_slice_node(SpreadReduction, 1, n, Ntimesteps, "SpreadReduction") else 0,
                 generations = max(ExtinctionGenerations, Ntimesteps))
    ex <- do.call(inapest_exogenous_extinction, pars); qh <- ex$history[, Ntimesteps]; qe <- ex$extinction
    extinction_method <- "binary edge-based multitype branching PGF"
  } else if (static_ops && dynamic && Model == "INApest" && !programmed_dynamic) {
    ex <- inapest_detection_extinction(.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),
                                       .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),
                                       .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                       .ina_slice_node(Survival,1,n,Ntimesteps,"Survival"),
                                       .ina_slice_node(DetectionProb,1,n,Ntimesteps,"DetectionProb"),
                                       .ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb"),
                                       .ina_slice_node(EradicationProb,1,n,Ntimesteps,"EradicationProb"),
                                       .ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction"),
                                       SEAM, .ina_slice_node(InfoRetentionProb,1,n,Ntimesteps,"InfoRetentionProb"),
                                       generations = max(ExtinctionGenerations,Ntimesteps))
    qh <- c(ex$historyU[, Ntimesteps], ex$historyH[, Ntimesteps]); qe <- c(ex$U, ex$H)
    extinction_method <- "binary informed/uninformed branching PGF"
  } else if (static_ops && !dynamic && Model == "INApestMeta") {
    ex <- meta_single_parent_extinction(.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),
                                        .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),LDDrate,
                                        .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                        .ina_slice_node(Survival,1,n,Ntimesteps,"Survival"),
                                        .ina_slice_node(K,1,n,Ntimesteps,"K"),
                                        .ina_slice_node(PropaguleProduction,1,n,Ntimesteps,"PropaguleProduction"),
                                        .ina_slice_node(PropaguleEstablishment,1,n,Ntimesteps,"PropaguleEstablishment"),
                                        if(mode=="all_informed").ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb") else 0,
                                        if(mode=="all_informed").ina_slice_node(MortalityProb,1,n,Ntimesteps,"MortalityProb") else 0,
                                        if(mode=="all_informed").ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction") else 0,
                                        generations=max(ExtinctionGenerations,Ntimesteps),
                                        FecundityReduction=if(mode=="all_informed").ina_slice_node(FecundityReduction,1,n,Ntimesteps,"FecundityReduction") else 0)
    qh <- ex$history[, Ntimesteps]; qe <- ex$extinction; extinction_method <- "simulator-faithful one-parent Meta branching PGF"
  } else if (static_ops && !dynamic && Model == "INApestMetaTransitionMatrix") {
    ex <- transition_extinction(Transition,Nstages,.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),
                                .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),LDDrate,
                                .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                .ina_slice_node(PropaguleEstablishment,1,n,Ntimesteps,"PropaguleEstablishment"),
                                if(mode=="all_informed").ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb") else 0,
                                if(mode=="all_informed").ina_slice_stage(MortalityProb,1,n,Nstages,Ntimesteps,"MortalityProb") else 0,
                                if(mode=="all_informed").ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction") else 0,
                                DispersalDensityFactor,
                                .ina_slice_node(K,1,n,Ntimesteps,"K"),
                                .ina_slice_node(SeedbankK,1,n,Ntimesteps,"SeedbankK"),
                                generations=max(ExtinctionGenerations,Ntimesteps),
                                FecundityReduction=if(mode=="all_informed").ina_slice_fecundity_transition(FecundityReduction,1,n,Nstages,Ntimesteps,"FecundityReduction") else 0)
    qh <- ex$history[, Ntimesteps]; qe <- ex$extinction; extinction_method <- "stage x node branching PGF"
  } else if (static_ops && !dynamic && Model == "INApestMetaMultipleLandUse") {
    ex <- mlu_single_parent_extinction(SDDprob=.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),
                                       LDDprob=.ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),LDDrate=LDDrate,
                                       EnvEstabProb=.ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                       Survival=.ina_slice_node(Survival,1,n,Ntimesteps,"Survival"),
                                       K=.ina_slice_K_mlu(K,1,n,L,Ntimesteps),
                                       PropaguleProduction=.ina_slice_node(PropaguleProduction,1,n,Ntimesteps,"PropaguleProduction"),
                                       PropaguleEstablishment=.ina_slice_node(PropaguleEstablishment,1,n,Ntimesteps,"PropaguleEstablishment"),
                                       ManageProb=if(mode=="all_informed").ina_slice_mlu(ManageProb,1,n,L,Ntimesteps,"ManageProb") else 0,
                                       MortalityProb=if(mode=="all_informed").ina_slice_mlu(MortalityProb,1,n,L,Ntimesteps,"MortalityProb") else 0,
                                       SpreadReduction=if(mode=="all_informed").ina_slice_mlu(SpreadReduction,1,n,L,Ntimesteps,"SpreadReduction") else 0,
                                       generations=max(ExtinctionGenerations,Ntimesteps),current_code=FALSE,
                                       FecundityReduction=if(mode=="all_informed").ina_slice_mlu(FecundityReduction,1,n,L,Ntimesteps,"FecundityReduction") else 0)
    qh <- ex$history[, Ntimesteps]; qe <- ex$extinction; extinction_method <- "mean-matched MLU branching approximation"
  } else {
    qh <- ina_poisson_branching_extinction_horizon(operators)
    if (static_ops) qe <- ina_poisson_branching_extinction(operators[[1]], ExtinctionGenerations)$extinction
    extinction_method <- if (programmed_dynamic) {
      if (static_ops) "age-structured mean-matched multitype branching" else "time-inhomogeneous age-structured mean-matched multitype branching"
    } else if (static_ops) "mean-matched multitype Poisson branching fallback" else "time-inhomogeneous mean-matched multitype Poisson branching"
  }
  extinction_horizon <- .ina_overall_extinction(qh, state0)
  extinction_eventual <- if (!is.null(qe)) .ina_overall_extinction(qe, state0) else NA_real_

  escape <- NULL
  if (.ina_export_available(Model, ExportProb, ExportSDDprob, ExportLDDprob,
                            if(Model=="INApestMetaTransitionMatrix") TRUE else AssumeResidualExport)) {
    first <- ina_escape_first_moment_temporal(operators, state0, export_vectors)
    noesc <- NULL; escape_method <- NULL
    if (static_ops && !dynamic && Model == "INApest" && !is.null(ExportProb)) {
      X <- .ina_slice_export(ExportProb,1,n,Ntimesteps,"ExportProb"); pout <- .inapest_recycle(OutsideEstablishmentProb,n,"OutsideEstablishmentProb")
      Xeff <- sweep(X,1,pout,`*`)
      br <- inapest_escape_branching(.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),Xeff,
                                     .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),
                                     .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                     .ina_slice_node(Survival,1,n,Ntimesteps,"Survival"),
                                     if(mode=="all_informed").ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb") else 0,
                                     if(mode=="all_informed").ina_slice_node(EradicationProb,1,n,Ntimesteps,"EradicationProb") else 0,
                                     if(mode=="all_informed").ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction") else 0,
                                     Ntimesteps)
      noesc <- br$no_escape; escape_method <- "binary edge-based branching no-escape recursion"
    } else if (static_ops && dynamic && Model == "INApest" && !programmed_dynamic && !is.null(ExportProb)) {
      X <- .ina_slice_export(ExportProb,1,n,Ntimesteps,"ExportProb"); pout <- .inapest_recycle(OutsideEstablishmentProb,n,"OutsideEstablishmentProb")
      Xeff <- sweep(X,1,pout,`*`)
      br <- inapest_detection_escape_branching(.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),Xeff,
                                                .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),
                                                .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                                .ina_slice_node(Survival,1,n,Ntimesteps,"Survival"),
                                                .ina_slice_node(DetectionProb,1,n,Ntimesteps,"DetectionProb"),
                                                .ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb"),
                                                .ina_slice_node(EradicationProb,1,n,Ntimesteps,"EradicationProb"),
                                                .ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction"),SEAM,
                                                .ina_slice_node(InfoRetentionProb,1,n,Ntimesteps,"InfoRetentionProb"),Ntimesteps)
      noesc <- c(br$no_escape_U,br$no_escape_H); escape_method <- "binary informed/uninformed branching no-escape recursion"
    } else if (static_ops && !dynamic && Model == "INApestMeta") {
      br <- meta_single_parent_escape_branching(.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),
                                                .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),LDDrate,
                                                .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                                .ina_slice_node(Survival,1,n,Ntimesteps,"Survival"),
                                                .ina_slice_node(K,1,n,Ntimesteps,"K"),
                                                .ina_slice_node(PropaguleProduction,1,n,Ntimesteps,"PropaguleProduction"),
                                                .ina_slice_node(PropaguleEstablishment,1,n,Ntimesteps,"PropaguleEstablishment"),
                                                if(mode=="all_informed").ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb") else 0,
                                                if(mode=="all_informed").ina_slice_node(MortalityProb,1,n,Ntimesteps,"MortalityProb") else 0,
                                                if(mode=="all_informed").ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction") else 0,
                                                .ina_slice_export(ExportSDDprob,1,n,Ntimesteps,"ExportSDDprob"),
                                                .ina_slice_export(ExportLDDprob,1,n,Ntimesteps,"ExportLDDprob"),
                                                OutsideEstablishmentProb,AssumeResidualExport,Ntimesteps,
                                                FecundityReduction=if(mode=="all_informed").ina_slice_node(FecundityReduction,1,n,Ntimesteps,"FecundityReduction") else 0)
      noesc <- br$no_escape; escape_method <- "Meta one-parent branching no-escape PGF"
    } else if (static_ops && !dynamic && Model == "INApestMetaTransitionMatrix") {
      br <- transition_escape_branching(Transition,Nstages,.ina_slice_connection(SDDprob,1,n,Ntimesteps,"SDDprob"),
                                         .ina_slice_connection(LDDprob,1,n,Ntimesteps,"LDDprob"),LDDrate,
                                         .ina_slice_node(EnvEstabProb,1,n,Ntimesteps,"EnvEstabProb"),
                                         .ina_slice_node(PropaguleEstablishment,1,n,Ntimesteps,"PropaguleEstablishment"),
                                         if(mode=="all_informed").ina_slice_node(ManageProb,1,n,Ntimesteps,"ManageProb") else 0,
                                         if(mode=="all_informed").ina_slice_stage(MortalityProb,1,n,Nstages,Ntimesteps,"MortalityProb") else 0,
                                         if(mode=="all_informed").ina_slice_node(SpreadReduction,1,n,Ntimesteps,"SpreadReduction") else 0,
                                         DispersalDensityFactor,.ina_slice_node(K,1,n,Ntimesteps,"K"),
                                         .ina_slice_node(SeedbankK,1,n,Ntimesteps,"SeedbankK"),
                                         .ina_slice_export(ExportSDDprob,1,n,Ntimesteps,"ExportSDDprob"),
                                         .ina_slice_export(ExportLDDprob,1,n,Ntimesteps,"ExportLDDprob"),
                                         OutsideEstablishmentProb,TRUE,Ntimesteps,
                                         FecundityReduction=if(mode=="all_informed").ina_slice_fecundity_transition(FecundityReduction,1,n,Nstages,Ntimesteps,"FecundityReduction") else 0)
      noesc <- br$no_escape; escape_method <- "stage x node branching no-escape PGF"
    } else {
      noesc <- ina_poisson_branching_no_escape_horizon(operators, export_vectors)
      escape_method <- if (programmed_dynamic) {
        if (static_ops) "age-structured mean-matched branching no-escape recursion" else "time-inhomogeneous age-structured mean-matched branching no-escape recursion"
      } else if (static_ops) "mean-matched multitype branching no-escape recursion" else "time-inhomogeneous mean-matched multitype branching no-escape recursion"
    }
    branch_prob <- 1 - .ina_overall_no_escape(noesc,state0)
    escape <- list(Method = escape_method,
                   BranchingProbabilityByHorizon = branch_prob,
                   FirstMomentPoissonProbabilityByHorizon = first$EscapeProbability,
                   ExpectedSuccessfulEscapesByStep = first$ExpectedSuccessfulEscapesByStep,
                   ExpectedSuccessfulEscapesCumulative = first$ExpectedSuccessfulEscapesCumulative)
  } else {
    diagnostics <- c(diagnostics,
      if (Model == "INApest")
        "Escape not calculated: binary INApest requires explicit source x outside-destination ExportProb; outside probability cannot be recovered from the internal edge matrix."
      else
        "Escape not calculated: provide ExportSDDprob/ExportLDDprob, or set AssumeResidualExport=TRUE only when missing row mass genuinely represents outside dispersal.")
  }

  if (dynamic && Model != "INApest") diagnostics <- c(diagnostics,
    "For information-limited Meta/Transition/MLU models, branching extinction/escape may use a mean-matched multitype fallback because node-level shared information induces correlations not represented by independent individual types.")
  if (programmed_dynamic && Model != "INApest") diagnostics <- c(diagnostics,
    "The programmed information clock itself is represented explicitly. At low density, however, a management kill that removes the focal individual cannot also leave a surviving focal lineage; information retained by that kill can precondition later immigrants and is therefore a shared-node, higher-order effect best captured by simulation or a nonlinear node-level extension.")
  if (.ina_any_nonzero(FecundityReduction) && Model != "INApest") diagnostics <- c(diagnostics,
    "FecundityReduction is applied only to reproductive output under management, after management mortality, and is kept separate from SpreadReduction and local survival.")
  if (.ina_any_nonzero(DispersalDensityFactor)) diagnostics <- c(diagnostics,
    "DispersalDensityFactor is linearized at the pest-free state. It can strongly alter finite-density trajectories even when the rare-state multiplier is unchanged.")
  diagnostics <- unique(diagnostics)

  result <- list(
    Model = Model,
    InformationMode = mode,
    InformationStateMethod = information_state_method,
    Growth = list(
      EquivalentPerTimestepMultiplier = Rstep,
      CycleMultiplier = growth$CycleMultiplier,
      Classification = .ina_classify_growth(Rstep),
      IntrinsicLocalLambda = intrinsic,
      IntrinsicLocalLambdaInterpretation = if (Model == "INApestMetaTransitionMatrix" && !is.null(intrinsic)) {
        if (mode == "all_informed") "managed local transition multiplier" else if (mode == "none") "unmanaged local transition multiplier" else "conditional-on-information managed local transition multiplier; landscape information limitation is represented separately in the augmented information-state operator"
      } else NULL
    ),
    Trajectory = trajectory,
    Extinction = list(
      Method = extinction_method,
      ProbabilityByHorizon = extinction_horizon,
      BranchingFadeoutProbability = extinction_eventual
    ),
    Escape = escape,
    Diagnostics = diagnostics,
    ApproximationScope = c(
      "Low-density / rare-invasion screening unless otherwise stated",
      "Expected-parameter treatment of management/detection SD rather than annual random-effect integration",
      "Finite carrying capacity, lineage collisions and density dependence become increasingly important away from rarity",
      if (programmed_dynamic) "Programmed stopping uses explicit time-since-local-evidence states; shared-node evidence created by other individuals and information-only preconditioning remain approximations away from rarity" else "Memoryless information retention is represented with uninformed/informed states"
    )
  )
  if (ReturnOperators) result$Operators <- operators
  class(result) <- "INApestAnalyticalResult"
  result
}
###############################################################################
### Round 2 unified user-facing dispatcher
###
### The original INApestAnalytical() is retained as INApestAnalytical_legacy.
### Calls with FecundityReduction == 0 for the pre-existing model families are
### sent directly through that legacy function. Point-model calls are routed to
### the point analytical engines added in Round 2.
###############################################################################

INApestAnalytical_legacy <- INApestAnalytical

INApestAnalytical <- function(
    Model = c("INApest", "INApestMeta", "INApestMetaTransitionMatrix",
              "INApestMetaMultipleLandUse", "INApestMetaPoint",
              "INApestPointTransitionMatrix"),
    Ntimesteps = 10,
    InitialState,
    InitialInfo = 0,
    InformationMode = c("auto", "dynamic", "all_informed", "none"),
    ApplyInitialDetection = TRUE,
    SDDprob,
    LDDprob = 0,
    LDDrate = 0,
    EnvEstabProb = 1,
    Survival = 1,
    K = NULL,
    PropaguleProduction = NULL,
    PropaguleEstablishment = 1,
    Transition = NULL,
    Nstages = NULL,
    SeedbankK = NULL,
    DetectionProb = 0,
    ManageProb = 0,
    EradicationProb = 0,
    MortalityProb = 0,
    SpreadReduction = 0,
    SEAM = NULL,
    InfoRetentionProb = 1,
    InfoPersistenceSteps = NA,
    DispersalDensityFactor = 0,
    ExportProb = NULL,
    ExportSDDprob = NULL,
    ExportLDDprob = NULL,
    OutsideEstablishmentProb = 1,
    AssumeResidualExport = FALSE,
    ExtinctionGenerations = 100,
    ReturnOperators = FALSE,
    FecundityReduction = 0,
    ...) {

  Model <- match.arg(Model)
  dots <- list(...)

  if (Model %in% c("INApestMetaPoint", "INApestPointTransitionMatrix")) {
    if (!"InitialPoints" %in% names(dots))
      stop("Point analytical models require InitialPoints supplied in ...")
    if (!"SDDkernel" %in% names(dots))
      stop("Point analytical models require SDDkernel supplied in ...")

    if (Model == "INApestMetaPoint") {
      args <- c(list(
        Ntimesteps = Ntimesteps,
        InitialInfo = InitialInfo,
        InformationMode = InformationMode,
        ApplyInitialDetection = ApplyInitialDetection,
        Survival = Survival,
        PropaguleProduction = PropaguleProduction,
        PropaguleEstablishment = PropaguleEstablishment,
        EnvEstabProb = EnvEstabProb,
        LDDrate = LDDrate,
        DetectionProb = DetectionProb,
        ManageProb = ManageProb,
        MortalityProb = MortalityProb,
        FecundityReduction = FecundityReduction,
        SpreadReduction = SpreadReduction,
        InfoRetentionProb = InfoRetentionProb,
        InfoPersistenceSteps = InfoPersistenceSteps,
        OutsideEstablishmentProb = OutsideEstablishmentProb,
        ExtinctionGenerations = ExtinctionGenerations,
        ReturnOperators = ReturnOperators
      ), dots)
      return(do.call(INApestMetaPointAnalytical, args))
    }

    args <- c(list(
      Ntimesteps = Ntimesteps,
      Nstages = Nstages,
      Transition = Transition,
      InitialInfo = InitialInfo,
      InformationMode = InformationMode,
      ApplyInitialDetection = ApplyInitialDetection,
      PropaguleEstablishment = PropaguleEstablishment,
      EnvEstabProb = EnvEstabProb,
      LDDrate = LDDrate,
      DetectionProb = DetectionProb,
      ManageProb = ManageProb,
      MortalityProb = MortalityProb,
      FecundityReduction = FecundityReduction,
      SpreadReduction = SpreadReduction,
      InfoRetentionProb = InfoRetentionProb,
      InfoPersistenceSteps = InfoPersistenceSteps,
      OutsideEstablishmentProb = OutsideEstablishmentProb,
      ExtinctionGenerations = ExtinctionGenerations,
      ReturnOperators = ReturnOperators
    ), dots)
    return(do.call(INApestPointTransitionMatrixAnalytical, args))
  }

  if (length(dots))
    stop("Arguments in ... are reserved for the point-model analytical extensions.")

  legacy_args <- list(
    Model = Model,
    Ntimesteps = Ntimesteps,
    InitialState = InitialState,
    InitialInfo = InitialInfo,
    InformationMode = InformationMode,
    ApplyInitialDetection = ApplyInitialDetection,
    SDDprob = SDDprob,
    LDDprob = LDDprob,
    LDDrate = LDDrate,
    EnvEstabProb = EnvEstabProb,
    Survival = Survival,
    K = K,
    PropaguleProduction = PropaguleProduction,
    PropaguleEstablishment = PropaguleEstablishment,
    Transition = Transition,
    Nstages = Nstages,
    SeedbankK = SeedbankK,
    DetectionProb = DetectionProb,
    ManageProb = ManageProb,
    EradicationProb = EradicationProb,
    MortalityProb = MortalityProb,
    SpreadReduction = SpreadReduction,
    SEAM = SEAM,
    InfoRetentionProb = InfoRetentionProb,
    InfoPersistenceSteps = InfoPersistenceSteps,
    DispersalDensityFactor = DispersalDensityFactor,
    ExportProb = ExportProb,
    ExportSDDprob = ExportSDDprob,
    ExportLDDprob = ExportLDDprob,
    OutsideEstablishmentProb = OutsideEstablishmentProb,
    AssumeResidualExport = AssumeResidualExport,
    ExtinctionGenerations = ExtinctionGenerations,
    ReturnOperators = ReturnOperators
  )

  # Preserve the pre-existing analytical path byte-for-byte at the function
  # level when fecundity management is absent. The Round-2 helper wrappers also
  # delegate to their saved pre-fecundity versions in this case.
  if (!.ina_fr_nonzero(FecundityReduction))
    return(do.call(INApestAnalytical_legacy, legacy_args))

  do.call(INApestAnalytical_round2_core,
          c(legacy_args, list(FecundityReduction = FecundityReduction)))
}

###############################################################################
### Dynamic model-specific branching PGF extension (development, 19-Aug-2026)
###
### Finite-horizon extinction under time-varying habitat, connectivity and
### management is calculated by backward composition of timestep-specific
### offspring PGFs.  For Meta/transition/point-type processes, parent survival
### and the shared management draw are retained explicitly rather than inferred
### from the mean operator.  Dynamic information uses separate parent/recruit
### information transforms.  Shared information held by pest-free nodes remains
### outside the independent-lineage branching state and is diagnosed explicitly.
###############################################################################

.ina_branch_clip <- function(x) pmin(1, pmax(0, as.numeric(x)))

.ina_branch_step_eval <- function(step, qnext,
                                  mode = c("none", "all_informed", "dynamic")) {
  mode <- match.arg(mode)
  B <- length(step$q0)
  reqm <- c("Parent0", "Parent1", "Recruit0", "Recruit1")
  for (nm in reqm) {
    M <- as.matrix(step[[nm]])
    if (!all(dim(M) == c(B, B))) stop(nm, " must be B x B")
  }
  q0 <- .ina_branch_clip(step$q0)
  q1 <- .ina_branch_clip(step$q1)
  a <- .ina_branch_clip(step$adoption)
  if (length(a) == 1L) a <- rep(a, B)
  if (length(a) != B) stop("adoption must be scalar or length B")

  branch_value <- function(i, M, parent_transform, recruit_transform) {
    qs <- if (M == 0L) q0[i] else q1[i]
    if (qs <= 0) return(1)
    Pm <- if (M == 0L) step$Parent0[, i] else step$Parent1[, i]
    Rm <- if (M == 0L) step$Recruit0[, i] else step$Recruit1[, i]
    # Parent0/1 and Recruit0/1 are unconditional on the management-survival
    # gate.  Conditional on surviving that gate, the biological parent may
    # still die/leave the represented state during local transition/movement.
    psum <- sum(Pm)
    local <- 1 - psum / qs + sum((Pm / qs) * parent_transform)
    mu <- Rm / qs
    val <- (1 - qs) + qs * local * exp(sum(mu * (recruit_transform - 1)))
    .ina_branch_clip(val)
  }

  if (mode != "dynamic") {
    q <- .ina_branch_clip(qnext)
    if (length(q) != B) stop("qnext has wrong length for base branching step")
    out <- numeric(B)
    for (i in seq_len(B)) {
      f0 <- branch_value(i, 0L, q, q)
      if (mode == "all_informed" && a[i] > 0) {
        f1 <- branch_value(i, 1L, q, q)
        out[i] <- (1 - a[i]) * f0 + a[i] * f1
      } else out[i] <- f0
    }
    return(.ina_branch_clip(out))
  }

  qnext <- .ina_branch_clip(qnext)
  if (length(qnext) != 2L * B) stop("qnext has wrong length for dynamic branching step")
  qU <- qnext[seq_len(B)]
  qH <- qnext[B + seq_len(B)]
  D <- .ina_branch_clip(step$Detection)
  IR <- .ina_branch_clip(step$Retention)
  if (length(D) == 1L) D <- rep(D, B)
  if (length(IR) == 1L) IR <- rep(IR, B)
  if (length(D) != B || length(IR) != B)
    stop("Detection and Retention must be scalar or length B")
  C <- as.matrix(step$Transfer)
  if (!all(dim(C) == c(B, B))) stop("Transfer must be B x B")
  C[] <- pmin(1, pmax(0, C))

  # Uninformed biological sources: the surviving parent and newborn recruits
  # acquire information only by their own end-of-timestep detection.
  childU <- (1 - D) * qU + D * qH
  outU <- numeric(B)
  for (i in seq_len(B))
    outU[i] <- branch_value(i, 0L, childU, childU)

  # Informed sources: information retention belongs to the same biological
  # parent. New recruits do not inherit that information automatically; they
  # can be informed by direct source-to-recruit transfer and/or detection.
  outH <- numeric(B)
  parentH <- (1 - (IR + (1 - IR) * D)) * qU +
             (IR + (1 - IR) * D) * qH
  for (i in seq_len(B)) {
    hchild <- 1 - (1 - D) * (1 - C[i, ])
    recruitH <- (1 - hchild) * qU + hchild * qH
    f0 <- branch_value(i, 0L, parentH, recruitH)
    f1 <- branch_value(i, 1L, parentH, recruitH)
    outH[i] <- (1 - a[i]) * f0 + a[i] * f1
  }
  .ina_branch_clip(c(outU, outH))
}

.ina_branch_horizon <- function(steps, mode) {
  if (!length(steps)) stop("steps must be non-empty")
  B <- length(steps[[1L]]$q0)
  q <- rep(0, if (mode == "dynamic") 2L * B else B)
  for (tt in rev(seq_along(steps)))
    q <- .ina_branch_step_eval(steps[[tt]], q, mode)
  q
}

.ina_branch_static <- function(steps) {
  if (length(steps) <= 1L) return(TRUE)
  all(vapply(steps[-1L], function(z)
    isTRUE(all.equal(z, steps[[1L]], tolerance = 0, check.attributes = FALSE)),
    logical(1)))
}

.ina_branch_eventual <- function(step, mode, generations = 100,
                                 tolerance = 1e-12) {
  B <- length(step$q0)
  q <- rep(0, if (mode == "dynamic") 2L * B else B)
  for (gg in seq_len(generations)) {
    qo <- q
    q <- .ina_branch_step_eval(step, q, mode)
    if (max(abs(q - qo)) < tolerance) break
  }
  q
}

.ina_node_transfer_types <- function(SEAM, type_node) {
  B <- length(type_node); n <- max(type_node)
  if (is.null(SEAM) || length(SEAM) == 1L) Cn <- matrix(0, n, n)
  else {
    Cn <- as.matrix(SEAM)
    if (!all(dim(Cn) == c(n, n))) stop("SEAM must be nodes x nodes")
  }
  diag(Cn) <- 0
  Cn[type_node, type_node, drop = FALSE]
}

.ina_meta_branch_step <- function(SDDprob, LDDprob, LDDrate,
                                  EnvEstabProb, Survival, K,
                                  PropaguleProduction, PropaguleEstablishment,
                                  DetectionProb, ManageProb, MortalityProb,
                                  SpreadReduction, SEAM, InfoRetentionProb,
                                  FecundityReduction = 0) {
  n <- nrow(as.matrix(SDDprob))
  s <- .ina_recycle(Survival, n, "Survival")
  a <- .ina_recycle(ManageProb, n, "ManageProb")
  m <- .ina_recycle(MortalityProb, n, "MortalityProb")
  d <- .ina_recycle(DetectionProb, n, "DetectionProb")
  ir <- .ina_recycle(InfoRetentionProb, n, "InfoRetentionProb")
  mus <- meta_single_parent_recruit_means(
    SDDprob, LDDprob, LDDrate, EnvEstabProb, K,
    PropaguleProduction, PropaguleEstablishment, SpreadReduction,
    FecundityReduction)
  q0 <- s
  q1 <- s * (1 - m)
  Parent0 <- diag(q0, nrow = n)
  Parent1 <- diag(q1, nrow = n)
  Recruit0 <- sweep(t(mus[[1L]]), 2, q0, `*`)
  Recruit1 <- sweep(t(mus[[2L]]), 2, q1, `*`)
  C <- .ina_node_transfer_types(SEAM, seq_len(n))
  list(q0 = q0, q1 = q1, adoption = a,
       Parent0 = Parent0, Parent1 = Parent1,
       Recruit0 = Recruit0, Recruit1 = Recruit1,
       Detection = d, Retention = ir, Transfer = C,
       approximation = paste(
         "single-parent simulator-faithful recruitment means with",
         "independent-Poisson recruit branching"))
}

.ina_transition_branch_step <- function(Transition, Nstages, SDDprob,
                                        LDDprob, LDDrate, EnvEstabProb,
                                        PropaguleEstablishment,
                                        DetectionProb, ManageProb,
                                        MortalityProb, SpreadReduction,
                                        SEAM, InfoRetentionProb,
                                        DispersalDensityFactor, K, SeedbankK,
                                        FecundityReduction = 0) {
  z <- transition_components(
    Transition, Nstages, SDDprob, LDDprob, LDDrate,
    EnvEstabProb, PropaguleEstablishment, ManageProb, MortalityProb,
    SpreadReduction, DispersalDensityFactor, K, SeedbankK,
    FecundityReduction)
  n <- z$n; S <- z$S; B <- n * S
  if (is.null(z$fecundity_reduction)) z$fecundity_reduction <- matrix(0, n, S)
  idx <- function(i, s) (i - 1L) * S + s
  q0 <- rep(1, B); q1 <- numeric(B); a <- numeric(B)
  Parent0 <- Parent1 <- matrix(0, B, B)
  Recruit0 <- Recruit1 <- matrix(0, B, B)
  for (i in seq_len(n)) {
    Ai <- z$A[[i]]
    for (k in seq_len(S)) {
      src <- idx(i, k)
      q1[src] <- 1 - z$mortality[i, k]
      a[src] <- z$adoption[i]
      if (k < S) {
        Parent0[idx(i, k), src] <- Ai[k, k]
        Parent0[idx(i, k + 1L), src] <- Ai[k + 1L, k]
        Parent1[idx(i, k), src] <- q1[src] * Ai[k, k]
        Parent1[idx(i, k + 1L), src] <- q1[src] * Ai[k + 1L, k]
      } else {
        Parent0[idx(i, S), src] <- Ai[S, S]
        Parent1[idx(i, S), src] <- q1[src] * Ai[S, S]
      }
      if (k >= 2L && Ai[1L, k] > 0) {
        fec <- Ai[1L, k]
        fr <- z$fecundity_reduction[i, k]
        for (j in seq_len(n)) {
          nat <- (1 - z$LDDrate) * z$sdd_enabled[i] * z$SDD[i, j]
          hum0 <- z$LDDrate * z$LDD[i, j]
          hum1 <- z$LDDrate * (1 - z$spread_reduction[i]) * z$LDD[i, j]
          mu0 <- fec * (nat + hum0) * z$recruit_success[j]
          mu1 <- fec * (1 - fr) * (nat + hum1) * z$recruit_success[j]
          Recruit0[idx(j, 1L), src] <- q0[src] * mu0
          Recruit1[idx(j, 1L), src] <- q1[src] * mu1
        }
      }
    }
  }
  if (length(DetectionProb) == 1L) Dm <- matrix(DetectionProb, n, S)
  else if (length(DetectionProb) == S) Dm <- matrix(rep(DetectionProb, each = n), n, S)
  else if (is.matrix(DetectionProb) && all(dim(DetectionProb) == c(n, S))) Dm <- DetectionProb
  else stop("DetectionProb must be scalar, length Nstages, or nodes x Nstages")
  D <- as.vector(t(Dm))
  irn <- .ina_recycle(InfoRetentionProb, n, "InfoRetentionProb")
  IR <- rep(irn, each = S)
  type_node <- rep(seq_len(n), each = S)
  C <- .ina_node_transfer_types(SEAM, type_node)
  list(q0 = q0, q1 = q1, adoption = a,
       Parent0 = Parent0, Parent1 = Parent1,
       Recruit0 = Recruit0, Recruit1 = Recruit1,
       Detection = D, Retention = IR, Transfer = C,
       approximation = paste(
         "exact local stasis/progression conditional on management survival,",
         "with Poisson successful recruits"))
}

.ina_binary_step_eval <- function(step, qnext,
                                  mode = c("none", "all_informed", "dynamic")) {
  mode <- match.arg(mode)
  P <- step$P; n <- nrow(P)
  s <- step$s; a <- step$a; e <- step$e; r <- step$r
  if (mode != "dynamic") {
    q <- .ina_branch_clip(qnext)
    out <- numeric(n)
    for (i in seq_len(n)) {
      off <- seq_len(n) != i
      prod0 <- prod(1 - P[i, off] + P[i, off] * q[off])
      f0 <- (1 - s[i]) + s[i] * q[i] * prod0
      if (mode == "all_informed") {
        sm <- s[i] * (1 - e[i])
        Pm <- P[i, off] * (1 - r[i])
        prod1 <- prod(1 - Pm + Pm * q[off])
        f1 <- (1 - sm) + sm * q[i] * prod1
        out[i] <- (1 - a[i]) * f0 + a[i] * f1
      } else out[i] <- f0
    }
    return(.ina_branch_clip(out))
  }
  q <- .ina_branch_clip(qnext)
  if (length(q) != 2L * n) stop("dynamic binary qnext must have length 2N")
  qU <- q[seq_len(n)]; qH <- q[n + seq_len(n)]
  d <- step$d; ir <- step$ir; C <- step$C
  outU <- outH <- numeric(n)
  childU <- (1 - d) * qU + d * qH
  for (i in seq_len(n)) {
    off <- seq_len(n) != i
    prodU <- prod(1 - P[i, off] + P[i, off] * childU[off])
    selfU <- (1 - d[i]) * qU[i] + d[i] * qH[i]
    outU[i] <- (1 - s[i]) + s[i] * selfU * prodU

    hp <- ir[i] + (1 - ir[i]) * d[i]
    selfH <- (1 - hp) * qU[i] + hp * qH[i]
    hc <- 1 - (1 - d) * (1 - C[i, ])
    childH <- (1 - hc) * qU + hc * qH
    prod0 <- prod(1 - P[i, off] + P[i, off] * childH[off])
    f0 <- (1 - s[i]) + s[i] * selfH * prod0
    sm <- s[i] * (1 - e[i])
    Pm <- P[i, off] * (1 - r[i])
    prod1 <- prod(1 - Pm + Pm * childH[off])
    f1 <- (1 - sm) + sm * selfH * prod1
    outH[i] <- (1 - a[i]) * f0 + a[i] * f1
  }
  .ina_branch_clip(c(outU, outH))
}

.ina_binary_horizon <- function(steps, mode) {
  n <- nrow(steps[[1L]]$P)
  q <- rep(0, if (mode == "dynamic") 2L * n else n)
  for (tt in rev(seq_along(steps))) q <- .ina_binary_step_eval(steps[[tt]], q, mode)
  q
}

.ina_binary_eventual <- function(step, mode, generations = 100,
                                 tolerance = 1e-12) {
  n <- nrow(step$P)
  q <- rep(0, if (mode == "dynamic") 2L * n else n)
  for (gg in seq_len(generations)) {
    qo <- q; q <- .ina_binary_step_eval(step, q, mode)
    if (max(abs(q - qo)) < tolerance) break
  }
  q
}

.ina_dynamic_pgf_initial <- function(Model, mode, Ntimesteps, InitialState,
                                     InitialInfo, ApplyInitialDetection,
                                     DetectionProb, n, Nstages = NULL) {
  if (Model == "INApest") {
    base <- as.numeric(InitialState)
    if (mode != "dynamic") return(base)
    d1 <- .ina_slice_node(DetectionProb, 1, n, Ntimesteps, "DetectionProb")
    d1 <- .inapest_recycle(d1, n, "DetectionProb")
    info <- .inapest_recycle(InitialInfo, n, "InitialInfo")
    if (ApplyInitialDetection) info <- info + (1 - info) * d1
    info <- .ina_branch_clip(info)
    return(c(base * (1 - info), base * info))
  }
  if (Model == "INApestMeta") {
    base <- as.numeric(InitialState)
    if (mode != "dynamic") return(base)
    d1 <- .ina_slice_node(DetectionProb, 1, n, Ntimesteps, "DetectionProb")
    d1 <- .ina_recycle(d1, n, "DetectionProb")
    pdet <- if (ApplyInitialDetection) .ina_initial_detection_meta(base, d1) else rep(0, n)
    info <- .ina_recycle(InitialInfo, n, "InitialInfo")
    if (ApplyInitialDetection) info <- info + (1 - info) * pdet
    info <- .ina_branch_clip(info)
    return(c(base * (1 - info), base * info))
  }
  if (Model == "INApestMetaTransitionMatrix") {
    S <- as.integer(Nstages)
    X0 <- if (is.matrix(InitialState)) as.matrix(InitialState) else
      matrix(as.numeric(InitialState), nrow = n, ncol = S, byrow = TRUE)
    base <- as.vector(t(X0))
    if (mode != "dynamic") return(base)
    D1 <- .ina_slice_stage(DetectionProb, 1, n, S, Ntimesteps, "DetectionProb")
    if (length(D1) == 1L) Dm <- matrix(D1, n, S)
    else if (length(D1) == S) Dm <- matrix(rep(D1, each = n), n, S)
    else Dm <- as.matrix(D1)
    pdet <- if (ApplyInitialDetection) .ina_initial_detection_transition(X0, Dm) else rep(0, n)
    info <- .ina_recycle(InitialInfo, n, "InitialInfo")
    if (ApplyInitialDetection) info <- info + (1 - info) * pdet
    info <- .ina_branch_clip(info)
    return(.ina_split_information(base, info, n, S))
  }
  stop("Unsupported model in .ina_dynamic_pgf_initial")
}

.ina_dynamic_pgf_extinction <- function(
    Model, Ntimesteps, InitialState, InitialInfo, InformationMode,
    ApplyInitialDetection, SDDprob, LDDprob, LDDrate, EnvEstabProb,
    Survival, K, PropaguleProduction, PropaguleEstablishment,
    Transition, Nstages, SeedbankK, DetectionProb, ManageProb,
    EradicationProb, MortalityProb, SpreadReduction, SEAM,
    InfoRetentionProb, InfoPersistenceSteps, DispersalDensityFactor,
    FecundityReduction, ExtinctionGenerations) {

  if (any(!is.na(as.numeric(InfoPersistenceSteps)))) return(NULL)
  SDD0 <- if (length(dim(SDDprob)) == 3L) SDDprob[, , 1L] else as.matrix(SDDprob)
  n <- nrow(SDD0)
  mode <- .ina_info_mode(InformationMode, ManageProb, InitialInfo,
                         DetectionProb, SEAM, InfoRetentionProb,
                         InfoPersistenceSteps)
  if (!(mode %in% c("none", "all_informed", "dynamic"))) return(NULL)
  steps <- vector("list", Ntimesteps)

  if (Model == "INApest") {
    for (tt in seq_len(Ntimesteps)) {
      P <- inapest_edge_prob(
        .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob"),
        .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob"),
        .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb"))
      d <- .inapest_recycle(.ina_slice_node(DetectionProb, tt, n, Ntimesteps, "DetectionProb"), n, "DetectionProb")
      ir <- .inapest_recycle(.ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb"), n, "InfoRetentionProb")
      C <- if (is.null(SEAM) || length(SEAM) == 1L) matrix(0, n, n) else as.matrix(SEAM)
      if (!all(dim(C) == c(n, n))) stop("SEAM must be nodes x nodes"); diag(C) <- 0
      steps[[tt]] <- list(
        P = P,
        s = .inapest_recycle(.ina_slice_node(Survival, tt, n, Ntimesteps, "Survival"), n, "Survival"),
        d = d,
        a = .inapest_recycle(.ina_slice_node(ManageProb, tt, n, Ntimesteps, "ManageProb"), n, "ManageProb"),
        e = .inapest_recycle(.ina_slice_node(EradicationProb, tt, n, Ntimesteps, "EradicationProb"), n, "EradicationProb"),
        r = .inapest_recycle(.ina_slice_node(SpreadReduction, tt, n, Ntimesteps, "SpreadReduction"), n, "SpreadReduction"),
        ir = ir, C = C)
    }
    qh <- .ina_binary_horizon(steps, mode)
    static <- .ina_branch_static(steps)
    qe <- if (static) .ina_binary_eventual(steps[[1L]], mode, ExtinctionGenerations) else NULL
    method <- if (mode == "dynamic")
      "time-inhomogeneous binary informed/uninformed edge-based branching PGF"
    else "time-inhomogeneous binary edge-based branching PGF"
  } else if (Model == "INApestMeta") {
    if (is.null(K) || is.null(PropaguleProduction)) return(NULL)
    for (tt in seq_len(Ntimesteps)) {
      steps[[tt]] <- .ina_meta_branch_step(
        .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob"),
        .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob"),
        LDDrate,
        .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb"),
        .ina_slice_node(Survival, tt, n, Ntimesteps, "Survival"),
        .ina_slice_node(K, tt, n, Ntimesteps, "K"),
        .ina_slice_node(PropaguleProduction, tt, n, Ntimesteps, "PropaguleProduction"),
        .ina_slice_node(PropaguleEstablishment, tt, n, Ntimesteps, "PropaguleEstablishment"),
        .ina_slice_node(DetectionProb, tt, n, Ntimesteps, "DetectionProb"),
        .ina_slice_node(ManageProb, tt, n, Ntimesteps, "ManageProb"),
        .ina_slice_node(MortalityProb, tt, n, Ntimesteps, "MortalityProb"),
        .ina_slice_node(SpreadReduction, tt, n, Ntimesteps, "SpreadReduction"),
        SEAM,
        .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb"),
        .ina_slice_node(FecundityReduction, tt, n, Ntimesteps, "FecundityReduction"))
    }
    qh <- .ina_branch_horizon(steps, mode)
    static <- .ina_branch_static(steps)
    qe <- if (static) .ina_branch_eventual(steps[[1L]], mode, ExtinctionGenerations) else NULL
    method <- if (mode == "dynamic")
      "time-inhomogeneous information-aware one-parent Meta branching PGF"
    else "time-inhomogeneous simulator-faithful one-parent Meta branching PGF"
  } else if (Model == "INApestMetaTransitionMatrix") {
    if (is.null(Transition)) return(NULL)
    if (is.null(Nstages)) Nstages <- if (is.list(Transition))
      nrow(as.matrix(Transition[[1L]])) else nrow(as.matrix(Transition))
    S <- as.integer(Nstages)
    if (is.null(K)) K <- 1
    if (is.null(SeedbankK)) SeedbankK <- K
    for (tt in seq_len(Ntimesteps)) {
      Ft <- .ina_slice_fecundity_transition(FecundityReduction, tt, n, S,
                                            Ntimesteps, "FecundityReduction")
      steps[[tt]] <- .ina_transition_branch_step(
        Transition, S,
        .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob"),
        .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob"),
        LDDrate,
        .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb"),
        .ina_slice_node(PropaguleEstablishment, tt, n, Ntimesteps, "PropaguleEstablishment"),
        .ina_slice_stage(DetectionProb, tt, n, S, Ntimesteps, "DetectionProb"),
        .ina_slice_node(ManageProb, tt, n, Ntimesteps, "ManageProb"),
        .ina_slice_stage(MortalityProb, tt, n, S, Ntimesteps, "MortalityProb"),
        .ina_slice_node(SpreadReduction, tt, n, Ntimesteps, "SpreadReduction"),
        SEAM,
        .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb"),
        DispersalDensityFactor,
        .ina_slice_node(K, tt, n, Ntimesteps, "K"),
        .ina_slice_node(SeedbankK, tt, n, Ntimesteps, "SeedbankK"),
        Ft)
    }
    qh <- .ina_branch_horizon(steps, mode)
    static <- .ina_branch_static(steps)
    qe <- if (static) .ina_branch_eventual(steps[[1L]], mode, ExtinctionGenerations) else NULL
    method <- if (mode == "dynamic")
      "time-inhomogeneous information-aware stage x node branching PGF"
    else "time-inhomogeneous stage x node branching PGF"
  } else return(NULL)

  state0 <- .ina_dynamic_pgf_initial(
    Model, mode, Ntimesteps, InitialState, InitialInfo, ApplyInitialDetection,
    DetectionProb, n, Nstages)
  list(
    mode = mode,
    state0 = state0,
    qh = qh,
    qe = qe,
    static = static,
    ProbabilityByHorizon = .ina_overall_extinction(qh, state0),
    BranchingFadeoutProbability = if (is.null(qe)) NA_real_ else .ina_overall_extinction(qe, state0),
    Method = method,
    note = if (mode == "dynamic" && Model != "INApest")
      paste(
        "Parent survival and management correlations are retained; newborn",
        "recruits receive information through detection/direct transfer rather",
        "than inheriting parental information. Persistent information held by",
        "pest-free nodes and cross-lineage information remain shared-node effects.")
    else NULL)
}

###############################################################################
### User-facing post-processor: retain all existing calculations, replacing
### only extinction when the stronger model-specific time-dependent PGF applies.
###############################################################################

INApestAnalytical_pre_dynamic_pgf <- INApestAnalytical

INApestAnalytical <- function(
    Model = c("INApest", "INApestMeta", "INApestMetaTransitionMatrix",
              "INApestMetaMultipleLandUse", "INApestMetaPoint",
              "INApestPointTransitionMatrix"),
    Ntimesteps = 10,
    InitialState,
    InitialInfo = 0,
    InformationMode = c("auto", "dynamic", "all_informed", "none"),
    ApplyInitialDetection = TRUE,
    SDDprob,
    LDDprob = 0,
    LDDrate = 0,
    EnvEstabProb = 1,
    Survival = 1,
    K = NULL,
    PropaguleProduction = NULL,
    PropaguleEstablishment = 1,
    Transition = NULL,
    Nstages = NULL,
    SeedbankK = NULL,
    DetectionProb = 0,
    ManageProb = 0,
    EradicationProb = 0,
    MortalityProb = 0,
    SpreadReduction = 0,
    SEAM = NULL,
    InfoRetentionProb = 1,
    InfoPersistenceSteps = NA,
    DispersalDensityFactor = 0,
    ExportProb = NULL,
    ExportSDDprob = NULL,
    ExportLDDprob = NULL,
    OutsideEstablishmentProb = 1,
    AssumeResidualExport = FALSE,
    ExtinctionGenerations = 100,
    ReturnOperators = FALSE,
    FecundityReduction = 0,
    ...) {

  Model <- match.arg(Model)
  dots <- list(...)
  if (Model == "INApestMetaPoint") {
    return(do.call(INApestMetaPointAnalytical, c(list(
      Ntimesteps = Ntimesteps, InitialInfo = InitialInfo,
      InformationMode = InformationMode,
      ApplyInitialDetection = ApplyInitialDetection, Survival = Survival,
      PropaguleProduction = PropaguleProduction,
      PropaguleEstablishment = PropaguleEstablishment,
      EnvEstabProb = EnvEstabProb, LDDrate = LDDrate,
      DetectionProb = DetectionProb, ManageProb = ManageProb,
      MortalityProb = MortalityProb, FecundityReduction = FecundityReduction,
      SpreadReduction = SpreadReduction, InfoRetentionProb = InfoRetentionProb,
      InfoPersistenceSteps = InfoPersistenceSteps,
      OutsideEstablishmentProb = OutsideEstablishmentProb,
      ExtinctionGenerations = ExtinctionGenerations,
      ReturnOperators = ReturnOperators), dots)))
  }
  if (Model == "INApestPointTransitionMatrix") {
    return(do.call(INApestPointTransitionMatrixAnalytical, c(list(
      Ntimesteps = Ntimesteps, Nstages = Nstages, Transition = Transition,
      InitialInfo = InitialInfo, InformationMode = InformationMode,
      ApplyInitialDetection = ApplyInitialDetection,
      PropaguleEstablishment = PropaguleEstablishment,
      EnvEstabProb = EnvEstabProb, LDDrate = LDDrate,
      DetectionProb = DetectionProb, ManageProb = ManageProb,
      MortalityProb = MortalityProb, FecundityReduction = FecundityReduction,
      SpreadReduction = SpreadReduction, InfoRetentionProb = InfoRetentionProb,
      InfoPersistenceSteps = InfoPersistenceSteps,
      OutsideEstablishmentProb = OutsideEstablishmentProb,
      ExtinctionGenerations = ExtinctionGenerations,
      ReturnOperators = ReturnOperators), dots)))
  }
  args <- list(
    Model = Model, Ntimesteps = Ntimesteps, InitialState = InitialState,
    InitialInfo = InitialInfo, InformationMode = InformationMode,
    ApplyInitialDetection = ApplyInitialDetection, SDDprob = SDDprob,
    LDDprob = LDDprob, LDDrate = LDDrate, EnvEstabProb = EnvEstabProb,
    Survival = Survival, K = K, PropaguleProduction = PropaguleProduction,
    PropaguleEstablishment = PropaguleEstablishment, Transition = Transition,
    Nstages = Nstages, SeedbankK = SeedbankK, DetectionProb = DetectionProb,
    ManageProb = ManageProb, EradicationProb = EradicationProb,
    MortalityProb = MortalityProb, SpreadReduction = SpreadReduction,
    SEAM = SEAM, InfoRetentionProb = InfoRetentionProb,
    InfoPersistenceSteps = InfoPersistenceSteps,
    DispersalDensityFactor = DispersalDensityFactor, ExportProb = ExportProb,
    ExportSDDprob = ExportSDDprob, ExportLDDprob = ExportLDDprob,
    OutsideEstablishmentProb = OutsideEstablishmentProb,
    AssumeResidualExport = AssumeResidualExport,
    ExtinctionGenerations = ExtinctionGenerations,
    ReturnOperators = ReturnOperators, FecundityReduction = FecundityReduction)
  res <- do.call(INApestAnalytical_pre_dynamic_pgf, c(args, dots))

  if (Model %in% c("INApest", "INApestMeta", "INApestMetaTransitionMatrix")) {
    improved <- .ina_dynamic_pgf_extinction(
      Model, Ntimesteps, InitialState, InitialInfo, InformationMode,
      ApplyInitialDetection, SDDprob, LDDprob, LDDrate, EnvEstabProb,
      Survival, K, PropaguleProduction, PropaguleEstablishment,
      Transition, Nstages, SeedbankK, DetectionProb, ManageProb,
      EradicationProb, MortalityProb, SpreadReduction, SEAM,
      InfoRetentionProb, InfoPersistenceSteps, DispersalDensityFactor,
      FecundityReduction, ExtinctionGenerations)
    # Preserve the existing specialised static solutions byte-for-byte at the
    # result level. Replace extinction only where the new PGF adds capability:
    # time-varying biology/management, dynamic information in Meta/transition,
    # or programmed-information Meta/transition.
    programmed_here <- any(!is.na(as.numeric(InfoPersistenceSteps)))
    use_improved <- !is.null(improved) && (
      !isTRUE(improved$static) ||
      (identical(improved$mode, "dynamic") &&
       Model %in% c("INApestMeta", "INApestMetaTransitionMatrix")) ||
      (programmed_here && Model %in% c("INApestMeta", "INApestMetaTransitionMatrix"))
    )
    if (use_improved) {
      old_info_fallback_diag <- "For information-limited Meta/Transition/MLU models, branching extinction/escape may use a mean-matched multitype fallback because node-level shared information induces correlations not represented by independent individual types."
      res$Diagnostics <- res$Diagnostics[res$Diagnostics != old_info_fallback_diag]
      res$Extinction$Method <- improved$Method
      res$Extinction$ProbabilityByHorizon <- improved$ProbabilityByHorizon
      res$Extinction$BranchingFadeoutProbability <- improved$BranchingFadeoutProbability
      if (!is.null(improved$note)) res$Diagnostics <- unique(c(res$Diagnostics, improved$note))
      res$Diagnostics <- unique(c(res$Diagnostics,
        if (improved$static)
          "Extinction uses the model-specific branching PGF; static conditions also permit a fixed-point fade-out probability."
        else
          "Finite-horizon extinction uses backward composition of timestep-specific model PGFs, so time-varying habitat/connectivity/management no longer requires the generic mean-matched Poisson fallback."))
    }
  }
  res
}

###############################################################################
### Dynamic-information initial-condition mixture correction
###
### Initial information/detection is a node-level event.  For a node containing
### x starting individuals, the extinction probability is therefore a mixture
### of q_U^x and q_H^x, not the geometric mean obtained by splitting x into
### fractional U/H lineage counts.  The fractional split remains appropriate
### for the expected-state trajectory, but not for a distributional event.
###############################################################################

.ina_dynamic_initial_info_probability <- function(Model, Ntimesteps,
                                                   InitialState, InitialInfo,
                                                   ApplyInitialDetection,
                                                   DetectionProb, n,
                                                   Nstages = NULL) {
  p0 <- .ina_recycle(InitialInfo, n, "InitialInfo")
  p0 <- pmin(1, pmax(0, p0))
  if (!ApplyInitialDetection) return(p0)
  if (Model == "INApest") {
    d1 <- .ina_slice_node(DetectionProb, 1, n, Ntimesteps, "DetectionProb")
    pd <- .inapest_recycle(d1, n, "DetectionProb")
  } else if (Model == "INApestMeta") {
    x <- as.numeric(InitialState)
    d1 <- .ina_slice_node(DetectionProb, 1, n, Ntimesteps, "DetectionProb")
    pd <- .ina_initial_detection_meta(x, .ina_recycle(d1, n, "DetectionProb"))
  } else if (Model == "INApestMetaTransitionMatrix") {
    S <- as.integer(Nstages)
    X <- if (is.matrix(InitialState)) as.matrix(InitialState) else
      matrix(as.numeric(InitialState), nrow = n, ncol = S, byrow = TRUE)
    D1 <- .ina_slice_stage(DetectionProb, 1, n, S, Ntimesteps, "DetectionProb")
    if (length(D1) == 1L) Dm <- matrix(D1, n, S)
    else if (length(D1) == S) Dm <- matrix(rep(D1, each = n), n, S)
    else Dm <- as.matrix(D1)
    pd <- .ina_initial_detection_transition(X, Dm)
  } else stop("Unsupported model for dynamic initial-information mixture")
  pmin(1, pmax(0, p0 + (1 - p0) * pd))
}

.ina_dynamic_initial_extinction <- function(q, Model, Ntimesteps,
                                             InitialState, InitialInfo,
                                             ApplyInitialDetection,
                                             DetectionProb, n,
                                             Nstages = NULL) {
  q <- pmin(1, pmax(0, as.numeric(q)))
  pH <- .ina_dynamic_initial_info_probability(
    Model, Ntimesteps, InitialState, InitialInfo, ApplyInitialDetection,
    DetectionProb, n, Nstages)
  if (Model %in% c("INApest", "INApestMeta")) {
    x <- as.numeric(InitialState)
    if (length(x) != n) stop("InitialState length mismatch")
    qU <- q[seq_len(n)]; qH <- q[n + seq_len(n)]
    nodeU <- exp(x * log(pmax(qU, .Machine$double.xmin)))
    nodeH <- exp(x * log(pmax(qH, .Machine$double.xmin)))
    nodeU[qU == 0 & x > 0] <- 0
    nodeH[qH == 0 & x > 0] <- 0
    return(prod((1 - pH) * nodeU + pH * nodeH))
  }
  if (Model == "INApestMetaTransitionMatrix") {
    S <- as.integer(Nstages); B <- n * S
    X <- if (is.matrix(InitialState)) as.matrix(InitialState) else
      matrix(as.numeric(InitialState), nrow = n, ncol = S, byrow = TRUE)
    qU <- q[seq_len(B)]; qH <- q[B + seq_len(B)]
    out <- 1
    for (i in seq_len(n)) {
      ids <- (i - 1L) * S + seq_len(S); xi <- X[i, ]
      pu <- if (any(qU[ids] == 0 & xi > 0)) 0 else
        exp(sum(xi * log(pmax(qU[ids], .Machine$double.xmin))))
      ph <- if (any(qH[ids] == 0 & xi > 0)) 0 else
        exp(sum(xi * log(pmax(qH[ids], .Machine$double.xmin))))
      out <- out * ((1 - pH[i]) * pu + pH[i] * ph)
    }
    return(out)
  }
  stop("Unsupported model for dynamic initial extinction")
}

.ina_dynamic_pgf_extinction_pre_initial_mix <- .ina_dynamic_pgf_extinction
.ina_dynamic_pgf_extinction <- function(
    Model, Ntimesteps, InitialState, InitialInfo, InformationMode,
    ApplyInitialDetection, SDDprob, LDDprob, LDDrate, EnvEstabProb,
    Survival, K, PropaguleProduction, PropaguleEstablishment,
    Transition, Nstages, SeedbankK, DetectionProb, ManageProb,
    EradicationProb, MortalityProb, SpreadReduction, SEAM,
    InfoRetentionProb, InfoPersistenceSteps, DispersalDensityFactor,
    FecundityReduction, ExtinctionGenerations) {
  z <- .ina_dynamic_pgf_extinction_pre_initial_mix(
    Model, Ntimesteps, InitialState, InitialInfo, InformationMode,
    ApplyInitialDetection, SDDprob, LDDprob, LDDrate, EnvEstabProb,
    Survival, K, PropaguleProduction, PropaguleEstablishment,
    Transition, Nstages, SeedbankK, DetectionProb, ManageProb,
    EradicationProb, MortalityProb, SpreadReduction, SEAM,
    InfoRetentionProb, InfoPersistenceSteps, DispersalDensityFactor,
    FecundityReduction, ExtinctionGenerations)
  if (is.null(z) || z$mode != "dynamic") return(z)
  SDD0 <- if (length(dim(SDDprob)) == 3L) SDDprob[, , 1L] else as.matrix(SDDprob)
  n <- nrow(SDD0)
  z$ProbabilityByHorizon <- .ina_dynamic_initial_extinction(
    z$qh, Model, Ntimesteps, InitialState, InitialInfo,
    ApplyInitialDetection, DetectionProb, n, Nstages)
  if (!is.null(z$qe)) z$BranchingFadeoutProbability <- .ina_dynamic_initial_extinction(
    z$qe, Model, Ntimesteps, InitialState, InitialInfo,
    ApplyInitialDetection, DetectionProb, n, Nstages)
  z$note <- paste(z$note,
    "Initial node information/detection is mixed as a shared node-level event rather than as fractional independent U/H lineages.")
  z
}

###############################################################################
### Point-family model-specific branching PGFs
###
### Reproduction in both point families is Poisson.  Once continuous movement
### has been contracted to analytical spatial types, Poisson thinning gives a
### multitype Poisson recruit distribution directly.  The previous generic
### mean-matched point-grid extinction fallback is therefore unnecessary for
### memoryless information states: retain the one-parent survival/transition
### branch and the Poisson recruit PGF explicitly.
###############################################################################

.ina_pt_recover_branch_matrix <- function(M0, MH, adoption) {
  M0 <- as.matrix(M0); MH <- as.matrix(MH); a <- as.numeric(adoption)
  B <- ncol(M0); M1 <- matrix(0, nrow(M0), B)
  for (i in seq_len(B)) {
    if (a[i] > 1e-14) M1[, i] <- (MH[, i] - (1 - a[i]) * M0[, i]) / a[i]
    else M1[, i] <- M0[, i]
  }
  M1[abs(M1) < 1e-14] <- 0
  if (any(M1 < -1e-9)) stop("Could not recover managed point branching matrix")
  M1[] <- pmax(0, M1)
  M1
}

.ina_pt_metapoint_branch_step <- function(st, representatives,
                                           DetectionProb, DetectionSpatial,
                                           InfoRetentionProb, InfoRadius,
                                           InfoTransferProb, InfoKernel,
                                           timestep, Ntimesteps) {
  a <- as.numeric(st$branch$a)
  Parent1 <- .ina_pt_recover_branch_matrix(st$Parent0, st$ParentH, a)
  Recruit1 <- .ina_pt_recover_branch_matrix(st$Recruit0, st$RecruitH, a)
  D <- .ina_pt_expected_parameter(DetectionProb, DetectionSpatial,
                                  representatives, timestep, 1L, Ntimesteps,
                                  "DetectionProb")
  IR <- .ina_pt_point_values(InfoRetentionProb, representatives, timestep,
                             Ntimesteps, "InfoRetentionProb")
  C <- .ina_pt_direct_info_matrix(representatives, InfoRadius,
                                  InfoTransferProb, InfoKernel,
                                  timestep, Ntimesteps)
  list(q0 = as.numeric(st$branch$q0), q1 = as.numeric(st$branch$q1),
       adoption = a, Parent0 = st$Parent0, Parent1 = Parent1,
       Recruit0 = st$Recruit0, Recruit1 = Recruit1,
       Detection = D, Retention = IR, Transfer = C)
}

.ina_pt_transition_branch_step_from_operator <- function(st, representatives,
                                                          DetectionProb,
                                                          DetectionSpatial,
                                                          InfoRetentionProb,
                                                          InfoRadius,
                                                          InfoTransferProb,
                                                          InfoKernel,
                                                          timestep,
                                                          Ntimesteps) {
  a <- vapply(st$branch, function(z) z$a, numeric(1))
  q0 <- vapply(st$branch, function(z) z$q0, numeric(1))
  q1 <- vapply(st$branch, function(z) z$q1, numeric(1))
  Parent1 <- .ina_pt_recover_branch_matrix(st$Parent0, st$ParentH, a)
  Recruit1 <- .ina_pt_recover_branch_matrix(st$Recruit0, st$RecruitH, a)
  D <- .ina_pt_transition_prob_spatial_mean(
    DetectionProb, DetectionSpatial, representatives, timestep,
    Ntimesteps, max(representatives$stage), "DetectionProb")
  IR <- .ina_pt_point_values(InfoRetentionProb, representatives, timestep,
                             Ntimesteps, "InfoRetentionProb")
  C <- .ina_pt_direct_info_matrix(representatives, InfoRadius,
                                  InfoTransferProb, InfoKernel,
                                  timestep, Ntimesteps)
  list(q0 = q0, q1 = q1, adoption = a,
       Parent0 = st$Parent0, Parent1 = Parent1,
       Recruit0 = st$Recruit0, Recruit1 = Recruit1,
       Detection = D, Retention = IR, Transfer = C)
}

.ina_pt_initial_type_metapoint <- function(InitialPoints, analysis_grid) {
  if (is.null(analysis_grid)) return(rep(1L, nrow(InitialPoints)))
  .ina_pt_xy_to_cell(InitialPoints$x, InitialPoints$y, analysis_grid)
}

.ina_pt_initial_type_transition <- function(InitialPoints, analysis_grid,
                                             Nstages) {
  stage <- if ("stage" %in% names(InitialPoints)) as.integer(InitialPoints$stage)
           else rep(1L, nrow(InitialPoints))
  if (is.null(analysis_grid)) return(stage)
  nc <- analysis_grid$nrow * analysis_grid$ncol
  cell <- .ina_pt_xy_to_cell(InitialPoints$x, InitialPoints$y, analysis_grid)
  ifelse(is.na(cell), NA_integer_, (stage - 1L) * nc + cell)
}

.ina_pt_initial_event_probability <- function(q, InitialPoints, InitialInfo,
                                               initial_type, mode,
                                               DetectionProb,
                                               DetectionSpatial,
                                               Ntimesteps,
                                               transition = FALSE,
                                               Nstages = NULL) {
  raw <- .ina_pt_initial_info_raw(InitialPoints, InitialInfo)
  keep <- which(!is.na(initial_type))
  if (!length(keep)) return(1)
  if (mode != "dynamic") {
    qq <- q[initial_type[keep]]
    return(prod(qq))
  }
  B <- length(q) / 2L; qU <- q[seq_len(B)]; qH <- q[B + seq_len(B)]
  pts <- InitialPoints
  if (!"id" %in% names(pts)) pts$id <- seq_len(nrow(pts))
  if (!"stage" %in% names(pts)) pts$stage <- if (transition) 1L else "default"
  if (transition) {
    D <- .ina_pt_transition_prob_spatial_mean(DetectionProb, DetectionSpatial,
                                               pts, 1L, Ntimesteps, Nstages,
                                               "DetectionProb")
  } else {
    D <- .ina_pt_expected_parameter(DetectionProb, DetectionSpatial, pts,
                                    1L, 1L, Ntimesteps, "DetectionProb")
  }
  D <- .ina_pt_clip01(D)
  ans <- 1
  for (ii in keep) {
    pH <- if (raw[ii] > 0) 1 else D[ii]
    b <- initial_type[ii]
    ans <- ans * ((1 - pH) * qU[b] + pH * qH[b])
  }
  ans
}

###############################################################################
### Programmed-information branching PGF
###
### Extends the one-parent PGF to U, X, H1...Hmax and Overflow information
### states.  This preserves parent survival/transition and Poisson recruit
### distributions instead of applying a mean-matched Poisson PGF to the
### age-structured mean operator.  Persistent pest-free information sites and
### evidence generated by other lineages remain outside independent-lineage
### branching and are retained as an explicit limitation.
###############################################################################

.ina_branch_value_core <- function(step, i, managed,
                                   parent_transform, recruit_transform) {
  qs <- if (!managed) step$q0[i] else step$q1[i]
  qs <- .ina_branch_clip(qs)
  if (qs <= 0) return(1)
  Pm <- if (!managed) step$Parent0[, i] else step$Parent1[, i]
  Rm <- if (!managed) step$Recruit0[, i] else step$Recruit1[, i]
  psum <- sum(Pm)
  local <- 1 - psum / qs + sum((Pm / qs) * parent_transform)
  mu <- Rm / qs
  .ina_branch_clip((1 - qs) + qs * local * exp(sum(mu * (recruit_transform - 1))))
}

.ina_prog_uninformed_transform <- function(q, layout, d) {
  B <- layout$base_types
  vapply(seq_len(B), function(j)
    (1 - d[j]) * q[layout$U[j]] + d[j] * q[layout$H[[1L]][j]], numeric(1))
}

.ina_prog_parent_transform <- function(q, layout, d, ir, K,
                                       kind, age = NA_integer_) {
  B <- layout$base_types
  out <- numeric(B)
  for (j in seq_len(B)) {
    h1 <- q[layout$H[[1L]][j]]; u <- q[layout$U[j]]
    rem <- 0
    if (identical(kind, "X")) {
      if (is.na(K[j])) rem <- ir[j] * q[layout$X[j]] + (1 - ir[j]) * u
      else if (is.infinite(K[j])) rem <- q[layout$X[j]]
      else rem <- u
    } else if (identical(kind, "Overflow")) {
      if (is.na(K[j])) rem <- ir[j] * q[layout$Overflow[j]] + (1 - ir[j]) * u
      else if (is.infinite(K[j])) rem <- q[layout$Overflow[j]]
      else rem <- u
    } else {
      nextq <- if (age < layout$max_age) q[layout$H[[age + 1L]][j]] else q[layout$Overflow[j]]
      if (is.na(K[j])) rem <- ir[j] * nextq + (1 - ir[j]) * u
      else if (is.infinite(K[j])) rem <- nextq
      else if (age >= K[j]) rem <- u
      else rem <- nextq
    }
    out[j] <- d[j] * h1 + (1 - d[j]) * rem
  }
  out
}

.ina_prog_recruit_informed_transform <- function(q, layout, d, C, src) {
  B <- layout$base_types
  out <- numeric(B)
  for (j in seq_len(B)) {
    pdet <- d[j]; ptx <- .ina_branch_clip(C[src, j])
    out[j] <- pdet * q[layout$H[[1L]][j]] +
      (1 - pdet) * ptx * q[layout$X[j]] +
      (1 - pdet) * (1 - ptx) * q[layout$U[j]]
  }
  out
}

.ina_programmed_branch_step_eval <- function(step, qnext, layout,
                                              persistence_by_type) {
  B <- layout$base_types
  q <- .ina_branch_clip(qnext)
  if (length(q) != layout$size) stop("programmed qnext has wrong length")
  d <- .ina_branch_clip(rep_len(as.numeric(step$Detection), B))
  ir <- .ina_branch_clip(rep_len(as.numeric(step$Retention), B))
  K <- rep_len(as.numeric(persistence_by_type), B)
  C <- as.matrix(step$Transfer)
  if (!all(dim(C) == c(B, B))) stop("Transfer must be B x B")
  a <- .ina_branch_clip(rep_len(as.numeric(step$adoption), B))
  out <- numeric(layout$size)

  Utr <- .ina_prog_uninformed_transform(q, layout, d)
  for (i in seq_len(B))
    out[layout$U[i]] <- .ina_branch_value_core(step, i, FALSE, Utr, Utr)

  eval_informed <- function(i, kind, age = NA_integer_) {
    Ptr <- .ina_prog_parent_transform(q, layout, d, ir, K, kind, age)
    Rtr <- .ina_prog_recruit_informed_transform(q, layout, d, C, i)
    f0 <- .ina_branch_value_core(step, i, FALSE, Ptr, Rtr)
    f1 <- .ina_branch_value_core(step, i, TRUE, Ptr, Rtr)
    (1 - a[i]) * f0 + a[i] * f1
  }

  for (i in seq_len(B)) out[layout$X[i]] <- eval_informed(i, "X")
  for (age in seq_len(layout$max_age))
    for (i in seq_len(B)) out[layout$H[[age]][i]] <- eval_informed(i, "H", age)
  for (i in seq_len(B)) out[layout$Overflow[i]] <- eval_informed(i, "Overflow")
  .ina_branch_clip(out)
}

.ina_programmed_branch_horizon <- function(steps, persistence_steps, layout) {
  if (length(steps) != length(persistence_steps))
    stop("steps and persistence_steps must have equal length")
  q <- rep(0, layout$size)
  for (tt in rev(seq_along(steps)))
    q <- .ina_programmed_branch_step_eval(steps[[tt]], q, layout,
                                          persistence_steps[[tt]])
  q
}

.ina_programmed_branch_eventual <- function(step, persistence_step, layout,
                                             generations = 100,
                                             tolerance = 1e-12) {
  q <- rep(0, layout$size)
  for (gg in seq_len(generations)) {
    qo <- q
    q <- .ina_programmed_branch_step_eval(step, q, layout, persistence_step)
    if (max(abs(q - qo)) < tolerance) break
  }
  q
}

.ina_programmed_initial_extinction_meta <- function(q, InitialState, InitialInfo,
                                                     ApplyInitialDetection,
                                                     DetectionProb, n,
                                                     Ntimesteps, layout,
                                                     type_node, Nstages = NULL) {
  p0 <- .ina_recycle(InitialInfo, n, "InitialInfo"); p0 <- .ina_branch_clip(p0)
  if (is.null(Nstages)) {
    xmat <- matrix(as.numeric(InitialState), nrow = n, ncol = 1L)
    d1 <- .ina_recycle(.ina_slice_node(DetectionProb, 1L, n, Ntimesteps,
                                       "DetectionProb"), n, "DetectionProb")
    pd <- if (ApplyInitialDetection) .ina_initial_detection_meta(xmat[,1], d1) else rep(0,n)
    S <- 1L
  } else {
    S <- as.integer(Nstages)
    xmat <- if (is.matrix(InitialState)) as.matrix(InitialState) else
      matrix(as.numeric(InitialState), nrow=n, ncol=S, byrow=TRUE)
    D1 <- .ina_slice_stage(DetectionProb,1L,n,S,Ntimesteps,"DetectionProb")
    if (length(D1)==1L) Dm<-matrix(D1,n,S)
    else if(length(D1)==S) Dm<-matrix(rep(D1,each=n),n,S)
    else Dm<-as.matrix(D1)
    pd <- if (ApplyInitialDetection) .ina_initial_detection_transition(xmat,Dm) else rep(0,n)
  }
  ans <- 1
  for (i in seq_len(n)) {
    ids <- which(type_node == i); xi <- if (S==1L) xmat[i,1] else xmat[i,]
    pu_prob <- (1-p0[i])*(1-pd[i])
    px_prob <- p0[i]*(1-pd[i])
    ph_prob <- pd[i]
    powerprod <- function(indices, vals) {
      qq <- q[indices]
      if (any(qq==0 & vals>0)) return(0)
      exp(sum(vals*log(pmax(qq,.Machine$double.xmin))))
    }
    pu <- powerprod(layout$U[ids], xi)
    px <- powerprod(layout$X[ids], xi)
    ph <- powerprod(layout$H[[1L]][ids], xi)
    ans <- ans * (pu_prob*pu + px_prob*px + ph_prob*ph)
  }
  ans
}

.ina_programmed_initial_extinction_point <- function(q, InitialPoints,
                                                      InitialInfo, initial_type,
                                                      DetectionProb,
                                                      DetectionSpatial,
                                                      Ntimesteps, layout,
                                                      transition=FALSE,
                                                      Nstages=NULL) {
  raw <- .ina_pt_initial_info_raw(InitialPoints, InitialInfo)
  keep <- which(!is.na(initial_type)); if (!length(keep)) return(1)
  pts <- InitialPoints
  if (!"id" %in% names(pts)) pts$id <- seq_len(nrow(pts))
  if (!"stage" %in% names(pts)) pts$stage <- if (transition) 1L else "default"
  D <- if (transition)
    .ina_pt_transition_prob_spatial_mean(DetectionProb,DetectionSpatial,pts,1L,
                                         Ntimesteps,Nstages,"DetectionProb")
  else .ina_pt_expected_parameter(DetectionProb,DetectionSpatial,pts,1L,1L,
                                   Ntimesteps,"DetectionProb")
  D <- .ina_branch_clip(D)
  ans <- 1
  for (ii in keep) {
    b <- initial_type[ii]
    if (raw[ii] > 0) {
      ans <- ans * ((1-D[ii])*q[layout$X[b]] + D[ii]*q[layout$H[[1L]][b]])
    } else {
      ans <- ans * ((1-D[ii])*q[layout$U[b]] + D[ii]*q[layout$H[[1L]][b]])
    }
  }
  ans
}

###############################################################################
### Programmed-information model-specific PGF integration for Meta/transition
###
### Finite InfoPersistenceSteps no longer forces extinction onto the generic
### age-structured mean-matched Poisson fallback for these two abundance/stage
### families.  The biological parent branch and Poisson recruits are preserved
### while the information clock is carried as U/X/H-age/Overflow types.
###############################################################################

.ina_dynamic_pgf_extinction_pre_programmed <- .ina_dynamic_pgf_extinction

.ina_dynamic_pgf_extinction <- function(
    Model, Ntimesteps, InitialState, InitialInfo, InformationMode,
    ApplyInitialDetection, SDDprob, LDDprob, LDDrate, EnvEstabProb,
    Survival, K, PropaguleProduction, PropaguleEstablishment,
    Transition, Nstages, SeedbankK, DetectionProb, ManageProb,
    EradicationProb, MortalityProb, SpreadReduction, SEAM,
    InfoRetentionProb, InfoPersistenceSteps, DispersalDensityFactor,
    FecundityReduction, ExtinctionGenerations) {

  persistence_requested <- any(!is.na(as.numeric(InfoPersistenceSteps)))
  if (!persistence_requested ||
      !(Model %in% c("INApestMeta", "INApestMetaTransitionMatrix"))) {
    return(.ina_dynamic_pgf_extinction_pre_programmed(
      Model, Ntimesteps, InitialState, InitialInfo, InformationMode,
      ApplyInitialDetection, SDDprob, LDDprob, LDDrate, EnvEstabProb,
      Survival, K, PropaguleProduction, PropaguleEstablishment,
      Transition, Nstages, SeedbankK, DetectionProb, ManageProb,
      EradicationProb, MortalityProb, SpreadReduction, SEAM,
      InfoRetentionProb, InfoPersistenceSteps, DispersalDensityFactor,
      FecundityReduction, ExtinctionGenerations))
  }

  SDD0 <- if (length(dim(SDDprob)) == 3L) SDDprob[, , 1L] else as.matrix(SDDprob)
  n <- nrow(SDD0)
  mode <- .ina_info_mode(InformationMode, ManageProb, InitialInfo,
                         DetectionProb, SEAM, InfoRetentionProb,
                         InfoPersistenceSteps)
  if (mode != "dynamic") {
    # A finite persistence rule has no biological effect when information is not
    # dynamically represented, so retain the ordinary model-specific PGF.
    no_persistence <- rep(NA_real_, n)
    return(.ina_dynamic_pgf_extinction_pre_programmed(
      Model, Ntimesteps, InitialState, InitialInfo, InformationMode,
      ApplyInitialDetection, SDDprob, LDDprob, LDDrate, EnvEstabProb,
      Survival, K, PropaguleProduction, PropaguleEstablishment,
      Transition, Nstages, SeedbankK, DetectionProb, ManageProb,
      EradicationProb, MortalityProb, SpreadReduction, SEAM,
      InfoRetentionProb, no_persistence, DispersalDensityFactor,
      FecundityReduction, ExtinctionGenerations))
  }

  max_age <- .ina_programmed_global_max_age(InfoPersistenceSteps)
  steps <- vector("list", Ntimesteps)
  persistence <- vector("list", Ntimesteps)

  if (Model == "INApestMeta") {
    if (is.null(K) || is.null(PropaguleProduction)) return(NULL)
    B <- n
    type_node <- seq_len(n)
    layout <- .ina_programmed_layout(B, max_age)
    for (tt in seq_len(Ntimesteps)) {
      steps[[tt]] <- .ina_meta_branch_step(
        .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob"),
        .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob"),
        LDDrate,
        .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb"),
        .ina_slice_node(Survival, tt, n, Ntimesteps, "Survival"),
        .ina_slice_node(K, tt, n, Ntimesteps, "K"),
        .ina_slice_node(PropaguleProduction, tt, n, Ntimesteps, "PropaguleProduction"),
        .ina_slice_node(PropaguleEstablishment, tt, n, Ntimesteps, "PropaguleEstablishment"),
        .ina_slice_node(DetectionProb, tt, n, Ntimesteps, "DetectionProb"),
        .ina_slice_node(ManageProb, tt, n, Ntimesteps, "ManageProb"),
        .ina_slice_node(MortalityProb, tt, n, Ntimesteps, "MortalityProb"),
        .ina_slice_node(SpreadReduction, tt, n, Ntimesteps, "SpreadReduction"),
        SEAM,
        .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb"),
        .ina_slice_node(FecundityReduction, tt, n, Ntimesteps, "FecundityReduction"))
      persistence[[tt]] <- .ina_recycle(
        .ina_slice_node(InfoPersistenceSteps, tt, n, Ntimesteps,
                        "InfoPersistenceSteps"), n, "InfoPersistenceSteps")
    }
    qh <- .ina_programmed_branch_horizon(steps, persistence, layout)
    ext_h <- .ina_programmed_initial_extinction_meta(
      qh, InitialState, InitialInfo, ApplyInitialDetection, DetectionProb,
      n, Ntimesteps, layout, type_node)
  } else {
    if (is.null(Transition)) return(NULL)
    if (is.null(Nstages)) Nstages <- if (is.list(Transition))
      nrow(as.matrix(Transition[[1L]])) else nrow(as.matrix(Transition))
    S <- as.integer(Nstages)
    if (is.null(K)) K <- 1
    if (is.null(SeedbankK)) SeedbankK <- K
    B <- n * S
    type_node <- rep(seq_len(n), each = S)
    layout <- .ina_programmed_layout(B, max_age)
    for (tt in seq_len(Ntimesteps)) {
      Ft <- .ina_slice_fecundity_transition(FecundityReduction, tt, n, S,
                                            Ntimesteps, "FecundityReduction")
      steps[[tt]] <- .ina_transition_branch_step(
        Transition, S,
        .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob"),
        .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob"),
        LDDrate,
        .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb"),
        .ina_slice_node(PropaguleEstablishment, tt, n, Ntimesteps,
                        "PropaguleEstablishment"),
        .ina_slice_stage(DetectionProb, tt, n, S, Ntimesteps, "DetectionProb"),
        .ina_slice_node(ManageProb, tt, n, Ntimesteps, "ManageProb"),
        .ina_slice_stage(MortalityProb, tt, n, S, Ntimesteps, "MortalityProb"),
        .ina_slice_node(SpreadReduction, tt, n, Ntimesteps, "SpreadReduction"),
        SEAM,
        .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb"),
        DispersalDensityFactor,
        .ina_slice_node(K, tt, n, Ntimesteps, "K"),
        .ina_slice_node(SeedbankK, tt, n, Ntimesteps, "SeedbankK"),
        Ft)
      kp <- .ina_recycle(
        .ina_slice_node(InfoPersistenceSteps, tt, n, Ntimesteps,
                        "InfoPersistenceSteps"), n, "InfoPersistenceSteps")
      persistence[[tt]] <- rep(kp, each = S)
    }
    qh <- .ina_programmed_branch_horizon(steps, persistence, layout)
    ext_h <- .ina_programmed_initial_extinction_meta(
      qh, InitialState, InitialInfo, ApplyInitialDetection, DetectionProb,
      n, Ntimesteps, layout, type_node, Nstages = S)
  }

  static <- .ina_branch_static(steps) &&
    (length(persistence) == 1L ||
     all(vapply(persistence[-1L], function(z)
       isTRUE(all.equal(z, persistence[[1L]], tolerance = 0)), logical(1))))
  qe <- if (static) .ina_programmed_branch_eventual(
    steps[[1L]], persistence[[1L]], layout, ExtinctionGenerations) else NULL
  ext_e <- if (is.null(qe)) NA_real_ else
    .ina_programmed_initial_extinction_meta(
      qe, InitialState, InitialInfo, ApplyInitialDetection, DetectionProb,
      n, Ntimesteps, layout, type_node,
      Nstages = if (Model == "INApestMetaTransitionMatrix") S else NULL)

  list(
    mode = mode,
    state0 = NULL,
    qh = qh,
    qe = qe,
    static = static,
    ProbabilityByHorizon = ext_h,
    BranchingFadeoutProbability = ext_e,
    Method = if (Model == "INApestMeta") {
      if (static) "programmed-information one-parent Meta branching PGF" else
        "time-inhomogeneous programmed-information one-parent Meta branching PGF"
    } else {
      if (static) "programmed-information stage x node branching PGF" else
        "time-inhomogeneous programmed-information stage x node branching PGF"
    },
    note = paste(
      "The programmed information clock is carried inside the model-specific branching PGF rather than applied to a mean offspring matrix.",
      "Parent survival/transition and Poisson recruits remain distributional.",
      "Information persisting at pest-free nodes and cross-lineage shared-node information remain outside independent-lineage branching."))
}

###############################################################################
### Current-transition one-parent recruit means
###
### The current production transition simulator floors the total SDD/LDD
### multinomial size after Poisson propagule production.  The older analytical
### transition operator treated this flow as continuously divisible.  At rarity
### that can materially overstate recruitment.  These helpers integrate the
### exact one-parent Poisson support, marginal multinomial counts and the
### current seedbank recruitment rule, while retaining the branching
### approximation that recruit lineages are independent conditional on their
### simulator-faithful mean counts.
###############################################################################

.ina_expected_capped_two_binomials <- function(n1, p1, n2, p2, cap) {
  n1 <- max(0L, as.integer(floor(n1))); n2 <- max(0L, as.integer(floor(n2)))
  cap <- max(0L, as.integer(floor(cap)))
  if (cap <= 0L || (n1 <= 0L && n2 <= 0L)) return(0)
  p1 <- min(1,max(0,p1)); p2 <- min(1,max(0,p2))
  maxz <- min(cap - 1L, n1 + n2)
  if (maxz < 0L) return(cap)
  probs <- numeric(maxz + 1L)
  for (z in 0:maxz) {
    lo <- max(0L, z - n2); hi <- min(n1, z)
    if (lo <= hi) {
      xs <- lo:hi
      probs[z + 1L] <- sum(dbinom(xs, n1, p1) * dbinom(z - xs, n2, p2))
    }
  }
  # E[min(Z,cap)] = cap - sum_{z<cap}(cap-z) P(Z=z)
  z <- 0:maxz
  max(0, cap - sum((cap - z) * probs))
}

.ina_transition_recruit_mean_given_arrivals <- function(pin, qin, j,
                                                         seedbank_slots,
                                                         env_prob,
                                                         accessible_fraction,
                                                         unrestricted) {
  pin <- as.integer(pin); qin <- as.integer(qin)
  slots <- max(0L, as.integer(floor(seedbank_slots)))
  arrivals <- pin + qin
  if (slots <= 0L || arrivals <= 0L) return(0)
  if (unrestricted) {
    natural_slots <- if (pin > 0L) slots else 0L
  } else {
    natural_slots <- floor(slots * accessible_fraction[j])
    if (pin > 0L && accessible_fraction[j] > 0)
      natural_slots <- max(1L, natural_slots)
    natural_slots <- min(natural_slots, slots)
  }
  other_slots <- slots - natural_slots
  lambdaP <- if (natural_slots > 0L) env_prob[j] * pin / natural_slots else 0
  lambdaQ <- env_prob[j] * qin / slots
  pnat <- -expm1(-(lambdaP + lambdaQ))
  pother <- -expm1(-lambdaQ)
  .ina_expected_capped_two_binomials(natural_slots, pnat,
                                      other_slots, pother,
                                      arrivals)
}

.ina_transition_single_parent_recruit_means_current <- function(
    Transition, Nstages, SDDprob, LDDprob = 0, LDDrate = 0,
    EnvEstabProb = 1, PropaguleEstablishment = 1,
    SpreadReduction = 0, K = 1, SeedbankK = 1,
    FecundityReduction = 0) {
  SDD <- as.matrix(SDDprob); n <- nrow(SDD); S <- as.integer(Nstages)
  if (ncol(SDD) != n) stop("SDDprob must be square")
  LDD <- .ina_mat(LDDprob, n, "LDDprob")
  A <- .transition_list(Transition, n, S)
  env <- pmin(1,pmax(0,.ina_recycle(EnvEstabProb,n,"EnvEstabProb")))
  pe <- .ina_recycle(PropaguleEstablishment,n,"PropaguleEstablishment")
  cap <- .ina_recycle(K,n,"K"); sb <- .ina_recycle(SeedbankK,n,"SeedbankK")
  g <- .ina_recycle(SpreadReduction,n,"SpreadReduction")
  if (is.matrix(FecundityReduction) && all(dim(FecundityReduction)==c(n,S))) F <- FecundityReduction
  else if (length(FecundityReduction)==1L) F <- matrix(FecundityReduction,n,S)
  else if (length(FecundityReduction)==S) F <- matrix(rep(FecundityReduction,each=n),n,S)
  else if (length(FecundityReduction)==n && n != S) F <- matrix(rep(FecundityReduction,S),n,S)
  else stop("FecundityReduction must resolve to scalar, stage vector, node vector, or nodes x stages")
  if(any(F<0|F>1|!is.finite(F))) stop("FecundityReduction must be in [0,1]")
  r <- as.numeric(LDDrate); if(length(r)!=1L||r<0||r>1) stop("LDDrate must be in [0,1]")
  rs <- rowSums(SDD); rl <- rowSums(LDD)
  unrestricted <- all(is.finite(pe) & pe >= 1)
  B <- n*S; idx <- function(i,s)(i-1L)*S+s
  out <- list(matrix(0,B,B), matrix(0,B,B)) # unmanaged / managed source branches
  for (i in seq_len(n)) for (k in 2:S) {
    fec0 <- A[[i]][1L,k]
    if (!is.finite(fec0) || fec0 <= 0) next
    src <- idx(i,k)
    # One surviving reproductive parent defines the maternal natural-dispersal footprint.
    if (unrestricted) {
      acc <- rep(1,n)
    } else {
      coverage <- numeric(n)
      positive <- cap > 0
      coverage[positive] <- pe[i] * cap[i] * SDD[i,positive] / cap[positive]
      acc <- pmin(1,pmax(0,-expm1(-coverage)))
    }
    ps <- if(rs[i]>0) SDD[i,]/rs[i] else rep(0,n)
    pl <- if(rl[i]>0) LDD[i,]/rl[i] else rep(0,n)
    for (M in 0:1) {
      lambda <- fec0 * (1 - F[i,k]*M)
      supp <- .meta_poisson_support(lambda)
      meanj <- numeric(n)
      for (zz in seq_along(supp$k)) {
        P <- supp$k[zz]; pk <- supp$p[zz]
        ms <- floor(P * (1-r) * rs[i] + 1e-12)
        ml <- floor(P * r * (1-g[i]*M) * rl[i] + 1e-12)
        for (j in seq_len(n)) {
          # Only destination marginals are needed for the recruit mean.
          bs <- if(ms>0) dbinom(0:ms,ms,ps[j]) else 1
          bl <- if(ml>0) dbinom(0:ml,ml,pl[j]) else 1
          es <- 0
          for (u in 0:ms) for (v in 0:ml) {
            pr <- bs[u+1L] * bl[v+1L]
            if (pr==0) next
            es <- es + pr * .ina_transition_recruit_mean_given_arrivals(
              u,v,j,sb[j],env,acc,unrestricted)
          }
          meanj[j] <- meanj[j] + pk*es
        }
      }
      out[[M+1L]][idx(seq_len(n),1L),src] <- meanj
    }
  }
  out
}

# Override only the transition branching-step constructor.  Existing static
# all-informed/unmanaged user-facing paths remain preserved by the wrapper;
# dynamic/time-varying extinction gains the simulator-faithful low-count means.
.ina_transition_branch_step <- function(Transition, Nstages, SDDprob,
                                        LDDprob, LDDrate, EnvEstabProb,
                                        PropaguleEstablishment,
                                        DetectionProb, ManageProb,
                                        MortalityProb, SpreadReduction,
                                        SEAM, InfoRetentionProb,
                                        DispersalDensityFactor, K, SeedbankK,
                                        FecundityReduction = 0) {
  z <- transition_components(
    Transition, Nstages, SDDprob, LDDprob, LDDrate,
    EnvEstabProb, PropaguleEstablishment, ManageProb, MortalityProb,
    SpreadReduction, DispersalDensityFactor, K, SeedbankK,
    FecundityReduction)
  n <- z$n; S <- z$S; B <- n*S
  if (is.null(z$fecundity_reduction)) z$fecundity_reduction <- matrix(0,n,S)
  idx <- function(i,s)(i-1L)*S+s
  q0 <- rep(1,B); q1 <- numeric(B); a <- numeric(B)
  Parent0 <- Parent1 <- matrix(0,B,B)
  for(i in seq_len(n)) {
    Ai <- z$A[[i]]
    for(k in seq_len(S)) {
      src <- idx(i,k); q1[src] <- 1-z$mortality[i,k]; a[src] <- z$adoption[i]
      if(k<S) {
        Parent0[idx(i,k),src] <- Ai[k,k]
        Parent0[idx(i,k+1L),src] <- Ai[k+1L,k]
        Parent1[idx(i,k),src] <- q1[src]*Ai[k,k]
        Parent1[idx(i,k+1L),src] <- q1[src]*Ai[k+1L,k]
      } else {
        Parent0[idx(i,S),src] <- Ai[S,S]
        Parent1[idx(i,S),src] <- q1[src]*Ai[S,S]
      }
    }
  }
  if (!is.na(DispersalDensityFactor) && DispersalDensityFactor != 0) {
    # The current one-parent density-dependent transition kernel also depends on
    # stage weights, which are not yet exposed in INApestAnalytical(). Retain
    # the previous zero-density recruitment construction for this special case.
    old <- transition_components(Transition,Nstages,SDDprob,LDDprob,LDDrate,
      EnvEstabProb,PropaguleEstablishment,ManageProb,MortalityProb,
      SpreadReduction,DispersalDensityFactor,K,SeedbankK,FecundityReduction)
    Recruit0 <- Recruit1 <- matrix(0,B,B)
    for(i in seq_len(n)) for(k in 2:S) if(old$A[[i]][1,k]>0) {
      src<-idx(i,k);fec<-old$A[[i]][1,k];fr<-old$fecundity_reduction[i,k]
      for(j in seq_len(n)) {
        nat<-(1-old$LDDrate)*old$sdd_enabled[i]*old$SDD[i,j]
        hum0<-old$LDDrate*old$LDD[i,j]
        hum1<-old$LDDrate*(1-old$spread_reduction[i])*old$LDD[i,j]
        Recruit0[idx(j,1L),src]<-fec*(nat+hum0)*old$recruit_success[j]
        Recruit1[idx(j,1L),src]<-q1[src]*fec*(1-fr)*(nat+hum1)*old$recruit_success[j]
      }
    }
  } else {
    rm <- .ina_transition_single_parent_recruit_means_current(
      Transition,S,SDDprob,LDDprob,LDDrate,EnvEstabProb,
      PropaguleEstablishment,SpreadReduction,K,SeedbankK,z$fecundity_reduction)
    Recruit0 <- rm[[1L]]
    Recruit1 <- sweep(rm[[2L]],2,q1,`*`)
  }
  if(length(DetectionProb)==1L) Dm<-matrix(DetectionProb,n,S)
  else if(length(DetectionProb)==S) Dm<-matrix(rep(DetectionProb,each=n),n,S)
  else if(is.matrix(DetectionProb)&&all(dim(DetectionProb)==c(n,S))) Dm<-DetectionProb
  else stop("DetectionProb must be scalar, length Nstages, or nodes x Nstages")
  D<-as.vector(t(Dm)); irn<-.ina_recycle(InfoRetentionProb,n,"InfoRetentionProb")
  IR<-rep(irn,each=S); type_node<-rep(seq_len(n),each=S); C<-.ina_node_transfer_types(SEAM,type_node)
  list(q0=q0,q1=q1,adoption=a,Parent0=Parent0,Parent1=Parent1,
       Recruit0=Recruit0,Recruit1=Recruit1,Detection=D,Retention=IR,Transfer=C,
       approximation=paste("exact one-parent management/stage transition and simulator-faithful",
                           "Poisson/integer disperser recruitment means; recruit lineages mean-matched independent"))
}

###############################################################################
### Simulator-faithful transition mean operators for the cases where the
### stronger branching machinery adds capability.
###
### The same Parent/Recruit decomposition used by the transition PGF is used
### here so that reported rare-state growth and expected trajectories are
### internally consistent with the low-count integer disperser calculation.
###############################################################################

.ina_branch_mean_operator <- function(step,
                                      mode = c("none", "all_informed", "dynamic")) {
  mode <- match.arg(mode)
  B <- length(step$q0)
  mats <- list(step$Parent0, step$Parent1, step$Recruit0, step$Recruit1)
  if (any(vapply(mats, function(x) !all(dim(x) == c(B, B)), logical(1))))
    stop("Branch-step Parent/Recruit matrices must all be B x B")
  a <- .ina_branch_clip(rep_len(as.numeric(step$adoption), B))
  G0 <- step$Parent0 + step$Recruit0
  ParentH <- sweep(step$Parent0, 2L, 1 - a, `*`) +
             sweep(step$Parent1, 2L, a, `*`)
  RecruitH <- sweep(step$Recruit0, 2L, 1 - a, `*`) +
              sweep(step$Recruit1, 2L, a, `*`)
  GH <- ParentH + RecruitH
  if (mode == "none") return(G0)
  if (mode == "all_informed") return(GH)

  D <- .ina_branch_clip(rep_len(as.numeric(step$Detection), B))
  IR <- .ina_branch_clip(rep_len(as.numeric(step$Retention), B))
  C <- as.matrix(step$Transfer)
  if (!all(dim(C) == c(B, B))) stop("Branch-step Transfer matrix must be B x B")
  C[] <- pmin(1, pmax(0, C))

  G <- matrix(0, 2L * B, 2L * B)
  U <- seq_len(B); H <- B + seq_len(B)
  for (i in seq_len(B)) for (j in seq_len(B)) {
    # Uninformed biological parent/recruits: information can be generated only
    # by their own end-of-timestep detection.
    w <- step$Parent0[j, i]
    if (w != 0) {
      G[U[j], U[i]] <- G[U[j], U[i]] + w * (1 - D[j])
      G[H[j], U[i]] <- G[H[j], U[i]] + w * D[j]
    }
    w <- step$Recruit0[j, i]
    if (w != 0) {
      G[U[j], U[i]] <- G[U[j], U[i]] + w * (1 - D[j])
      G[H[j], U[i]] <- G[H[j], U[i]] + w * D[j]
    }

    # An informed surviving/progressing parent remains the same biological
    # individual, so information retention follows that parent. Recruits are
    # new individuals and can be informed only by transfer/detection.
    w <- ParentH[j, i]
    if (w != 0) {
      hp <- IR[j] + (1 - IR[j]) * D[j]
      G[U[j], H[i]] <- G[U[j], H[i]] + w * (1 - hp)
      G[H[j], H[i]] <- G[H[j], H[i]] + w * hp
    }
    w <- RecruitH[j, i]
    if (w != 0) {
      hr <- 1 - (1 - D[j]) * (1 - C[i, j])
      G[U[j], H[i]] <- G[U[j], H[i]] + w * (1 - hr)
      G[H[j], H[i]] <- G[H[j], H[i]] + w * hr
    }
  }
  attr(G, "note") <- paste(
    "Information-aware one-parent mean operator. Parent information follows",
    "the biological parent; recruits receive information only by detection or",
    "direct transfer. Information retained at pest-free nodes remains outside",
    "the independent-lineage representation.")
  G
}

.ina_transition_improved_mean <- function(
    Ntimesteps, InitialState, InitialInfo, InformationMode,
    ApplyInitialDetection, SDDprob, LDDprob, LDDrate, EnvEstabProb,
    K, PropaguleEstablishment, Transition, Nstages, SeedbankK,
    DetectionProb, ManageProb, MortalityProb, SpreadReduction, SEAM,
    InfoRetentionProb, InfoPersistenceSteps, DispersalDensityFactor,
    FecundityReduction) {

  SDD0 <- if (length(dim(SDDprob)) == 3L) SDDprob[, , 1L] else as.matrix(SDDprob)
  n <- nrow(SDD0)
  if (is.null(Transition)) return(NULL)
  if (is.null(Nstages)) Nstages <- if (is.list(Transition))
    nrow(as.matrix(Transition[[1L]])) else nrow(as.matrix(Transition))
  S <- as.integer(Nstages)
  if (is.null(K)) K <- 1
  if (is.null(SeedbankK)) SeedbankK <- K
  mode <- .ina_info_mode(InformationMode, ManageProb, InitialInfo,
                         DetectionProb, SEAM, InfoRetentionProb,
                         InfoPersistenceSteps)
  programmed <- identical(mode, "dynamic") && any(!is.na(as.numeric(InfoPersistenceSteps)))
  steps <- vector("list", Ntimesteps)
  for (tt in seq_len(Ntimesteps)) {
    Ft <- .ina_slice_fecundity_transition(FecundityReduction, tt, n, S,
                                          Ntimesteps, "FecundityReduction")
    steps[[tt]] <- .ina_transition_branch_step(
      Transition, S,
      .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob"),
      .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob"),
      LDDrate,
      .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb"),
      .ina_slice_node(PropaguleEstablishment, tt, n, Ntimesteps,
                      "PropaguleEstablishment"),
      .ina_slice_stage(DetectionProb, tt, n, S, Ntimesteps, "DetectionProb"),
      .ina_slice_node(ManageProb, tt, n, Ntimesteps, "ManageProb"),
      .ina_slice_stage(MortalityProb, tt, n, S, Ntimesteps, "MortalityProb"),
      .ina_slice_node(SpreadReduction, tt, n, Ntimesteps, "SpreadReduction"),
      SEAM,
      .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb"),
      DispersalDensityFactor,
      .ina_slice_node(K, tt, n, Ntimesteps, "K"),
      .ina_slice_node(SeedbankK, tt, n, Ntimesteps, "SeedbankK"),
      Ft)
  }

  X0 <- if (is.matrix(InitialState)) as.matrix(InitialState) else
    matrix(as.numeric(InitialState), nrow = n, ncol = S, byrow = TRUE)
  if (!all(dim(X0) == c(n, S)))
    stop("InitialState must be nodes x stages for INApestMetaTransitionMatrix")
  base <- as.vector(t(X0)); type_node <- rep(seq_len(n), each = S)

  if (!programmed) {
    ops <- lapply(steps, .ina_branch_mean_operator, mode = mode)
    state0 <- .ina_dynamic_pgf_initial(
      "INApestMetaTransitionMatrix", mode, Ntimesteps, InitialState,
      InitialInfo, ApplyInitialDetection, DetectionProb, n, S)
    return(list(mode = mode, programmed = FALSE, operators = ops,
                state0 = state0, steps = steps))
  }

  max_age <- .ina_programmed_global_max_age(InfoPersistenceSteps)
  layout <- .ina_programmed_layout(n * S, max_age)
  ops <- vector("list", Ntimesteps)
  for (tt in seq_len(Ntimesteps)) {
    step <- steps[[tt]]; B <- n * S
    a <- .ina_branch_clip(rep_len(as.numeric(step$adoption), B))
    ParentH <- sweep(step$Parent0, 2L, 1 - a, `*`) +
               sweep(step$Parent1, 2L, a, `*`)
    RecruitH <- sweep(step$Recruit0, 2L, 1 - a, `*`) +
                sweep(step$Recruit1, 2L, a, `*`)
    kp <- .ina_recycle(
      .ina_slice_node(InfoPersistenceSteps, tt, n, Ntimesteps,
                      "InfoPersistenceSteps"), n, "InfoPersistenceSteps")
    ops[[tt]] <- .ina_pt_programmed_operator(
      step$Parent0, ParentH, step$Recruit0, RecruitH,
      representatives = NULL,
      DetectionProbByType = step$Detection,
      InfoRetentionProbByType = step$Retention,
      InfoPersistenceStepsByType = rep(kp, each = S),
      InfoTransferMatrix = step$Transfer,
      layout = layout)
  }
  D1 <- .ina_slice_stage(DetectionProb, 1L, n, S, Ntimesteps, "DetectionProb")
  if (length(D1) == 1L) Dm <- matrix(D1, n, S)
  else if (length(D1) == S) Dm <- matrix(rep(D1, each = n), n, S)
  else Dm <- as.matrix(D1)
  pdet0 <- if (ApplyInitialDetection)
    .ina_initial_detection_transition(X0, Dm) else rep(0, n)
  state0 <- .ina_programmed_initial_state(
    base, InitialInfo, pdet0, type_node, layout,
    binary_known_presence = FALSE)
  list(mode = mode, programmed = TRUE, operators = ops, state0 = state0,
       steps = steps, layout = layout)
}

# Final transition-output harmoniser.  Static unmanaged/all-informed results are
# deliberately retained from the established analytical solution.  The new
# simulator-faithful mean operator is used where it adds capability: dynamic
# information/programmed stopping or genuinely time-varying transition inputs.
INApestAnalytical_pre_transition_mean <- INApestAnalytical

INApestAnalytical <- function(...) {
  args <- list(...)
  res <- do.call(INApestAnalytical_pre_transition_mean, args)
  Model <- if (!is.null(args$Model)) args$Model else "INApest"
  if (length(Model) > 1L) Model <- Model[1L]
  if (!identical(Model, "INApestMetaTransitionMatrix")) return(res)

  # Resolve arguments from the public defaults without altering the established
  # interface.  formals() gives defaults for omitted optional arguments.
  fm <- formals(INApestAnalytical_pre_transition_mean)
  getarg <- function(nm) {
    if (!is.null(args[[nm]])) return(args[[nm]])
    z <- fm[[nm]]
    if (is.null(z)) return(NULL)
    eval(z, envir = parent.frame())
  }
  Ntimesteps <- getarg("Ntimesteps")
  InitialState <- getarg("InitialState")
  InitialInfo <- getarg("InitialInfo")
  InformationMode <- getarg("InformationMode")
  ApplyInitialDetection <- getarg("ApplyInitialDetection")
  SDDprob <- getarg("SDDprob"); LDDprob <- getarg("LDDprob")
  LDDrate <- getarg("LDDrate"); EnvEstabProb <- getarg("EnvEstabProb")
  K <- getarg("K"); PropaguleEstablishment <- getarg("PropaguleEstablishment")
  Transition <- getarg("Transition"); Nstages <- getarg("Nstages")
  SeedbankK <- getarg("SeedbankK"); DetectionProb <- getarg("DetectionProb")
  ManageProb <- getarg("ManageProb"); MortalityProb <- getarg("MortalityProb")
  SpreadReduction <- getarg("SpreadReduction"); SEAM <- getarg("SEAM")
  InfoRetentionProb <- getarg("InfoRetentionProb")
  InfoPersistenceSteps <- getarg("InfoPersistenceSteps")
  DispersalDensityFactor <- getarg("DispersalDensityFactor")
  FecundityReduction <- getarg("FecundityReduction")
  ReturnOperators <- isTRUE(getarg("ReturnOperators"))

  improved <- .ina_transition_improved_mean(
    Ntimesteps, InitialState, InitialInfo, InformationMode,
    ApplyInitialDetection, SDDprob, LDDprob, LDDrate, EnvEstabProb,
    K, PropaguleEstablishment, Transition, Nstages, SeedbankK,
    DetectionProb, ManageProb, MortalityProb, SpreadReduction, SEAM,
    InfoRetentionProb, InfoPersistenceSteps, DispersalDensityFactor,
    FecundityReduction)
  if (is.null(improved)) return(res)

  # Use improved mean outputs only when their corresponding PGF is adding new
  # capability. This preserves historical static none/all-informed outputs.
  timevary <- !.ina_branch_static(improved$steps)
  use_improved <- timevary || identical(improved$mode, "dynamic") || improved$programmed
  if (!use_improved) return(res)

  growth <- .ina_cycle_growth(improved$operators)
  traj_state <- .ina_apply_operators(improved$operators, improved$state0)
  res$Growth$EquivalentPerTimestepMultiplier <- growth$EquivalentPerTimestepMultiplier
  res$Growth$CycleMultiplier <- growth$CycleMultiplier
  res$Growth$Classification <- .ina_classify_growth(growth$EquivalentPerTimestepMultiplier)
  res$Trajectory <- data.frame(
    timestep = seq_len(Ntimesteps),
    expected_state_total = colSums(traj_state))
  if (ReturnOperators) res$Operators <- improved$operators
  res$Diagnostics <- unique(c(res$Diagnostics,
    "Transition growth/trajectory use the simulator-faithful one-parent low-count dispersal operator in dynamic-information or time-varying cases; this is consistent with the transition branching PGF used for extinction.",
    if (!is.na(DispersalDensityFactor) && DispersalDensityFactor != 0)
      "For transition DispersalDensityFactor != 0, the one-parent rare-state recruitment mean retains the previous zero-density approximation because stage weights are not yet exposed through INApestAnalytical()." else NULL))
  res
}

###############################################################################
### Node-shared information refinement for dynamic Meta/transition branching
###
### In the node-based simulators, HaveInfo belongs to the node, not to an
### individual. A retained information state, SEAM transfer, or any detection
### therefore informs all extant descendants at that node. Conditional on the
### one-parent biological branch and the independent-Poisson recruit
### approximation already used by the PGF, this shared detection/information
### event can be integrated analytically. Cross-lineage information retained at
### an empty node remains outside an independent-lineage branching process.
###############################################################################

.ina_branch_step_eval_pre_nodeshared <- .ina_branch_step_eval

.ina_branch_step_eval_nodeshared_dynamic <- function(step, qnext) {
  B <- length(step$q0)
  type_node <- as.integer(step$SharedInfoNode)
  if (length(type_node) != B || any(!is.finite(type_node)) || any(type_node < 1L))
    stop("SharedInfoNode must identify a positive node index for every base type")
  nodes <- sort(unique(type_node))
  qnext <- .ina_branch_clip(qnext)
  if (length(qnext) != 2L * B) stop("qnext has wrong length for dynamic branching step")
  qU <- qnext[seq_len(B)]; qH <- qnext[B + seq_len(B)]
  q0 <- .ina_branch_clip(step$q0); q1 <- .ina_branch_clip(step$q1)
  a <- .ina_branch_clip(rep_len(as.numeric(step$adoption), B))
  D <- .ina_branch_clip(rep_len(as.numeric(step$Detection), B))
  IR <- .ina_branch_clip(rep_len(as.numeric(step$Retention), B))
  C <- as.matrix(step$Transfer)
  if (!all(dim(C) == c(B, B))) stop("Transfer must be B x B")
  C[] <- pmin(1, pmax(0, C))

  node_product <- function(i, M, informed_source) {
    qs <- if (M == 0L) q0[i] else q1[i]
    if (qs <= 0) return(1)
    Pm <- if (M == 0L) step$Parent0[, i] else step$Parent1[, i]
    Rm <- if (M == 0L) step$Recruit0[, i] else step$Recruit1[, i]
    # The mortality/survival gate is outside the local parent transition and
    # recruitment process. Condition on passing that gate here.
    ppar <- pmax(0, as.numeric(Pm) / qs)
    mu <- pmax(0, as.numeric(Rm) / qs)
    # The current conventional Meta/transition implementations keep the
    # surviving biological parent in its source node. If a future analytical
    # path allows parent relocation, fall back to the previous individual-type
    # evaluator until the mutually-exclusive parent-location event is handled.
    srcnode <- type_node[i]
    parent_other <- which(ppar > 1e-14 & type_node != srcnode)
    if (length(parent_other)) return(NA_real_)
    pabs <- pmax(0, 1 - sum(ppar))

    simple_node_gf <- function(types, x) {
      rec <- exp(sum(mu[types] * (x - 1)))
      if (type_node[i] == type_node[types[1L]]) {
        par <- pabs + sum(ppar[types] * x)
        # Any parent probabilities outside this node would have triggered the
        # fallback above. Parent absence includes local transition death.
        rec * par
      } else rec
    }

    ans <- 1
    for (nd in nodes) {
      types <- which(type_node == nd)
      GH <- simple_node_gf(types, qH[types])
      noU <- simple_node_gf(types, (1 - D[types]) * qU[types])
      noH <- simple_node_gf(types, (1 - D[types]) * qH[types])
      # If the node is not already informed before detection, all descendants
      # remain U only when none is detected; otherwise the whole occupied node
      # becomes H. This formula also evaluates to 1 for an empty destination.
      GU <- noU + GH - noH
      if (!informed_source) {
        pinfo <- 0
      } else if (nd == srcnode) {
        # Existing source-node information is retained as one shared event.
        pinfo <- IR[i]
      } else {
        # One SEAM edge event informs the destination node as a whole. The
        # type-level transfer matrix is constant across types sharing a node.
        pinfo <- C[i, types[1L]]
      }
      ans <- ans * (pinfo * GH + (1 - pinfo) * GU)
    }
    .ina_branch_clip((1 - qs) + qs * ans)
  }

  outU <- outH <- numeric(B)
  for (i in seq_len(B)) {
    u <- node_product(i, 0L, FALSE)
    if (is.na(u)) return(NULL)
    outU[i] <- u
    h0 <- node_product(i, 0L, TRUE)
    h1 <- node_product(i, 1L, TRUE)
    if (is.na(h0) || is.na(h1)) return(NULL)
    outH[i] <- (1 - a[i]) * h0 + a[i] * h1
  }
  .ina_branch_clip(c(outU, outH))
}

.ina_branch_step_eval <- function(step, qnext,
                                  mode = c("none", "all_informed", "dynamic")) {
  mode <- match.arg(mode)
  if (mode == "dynamic" && !is.null(step$SharedInfoNode)) {
    z <- .ina_branch_step_eval_nodeshared_dynamic(step, qnext)
    if (!is.null(z)) return(z)
  }
  .ina_branch_step_eval_pre_nodeshared(step, qnext, mode)
}

# Mark only node-based Meta and conventional transition steps as sharing a
# node-level information state. Point models retain individual/site information
# semantics and therefore deliberately do not set SharedInfoNode.
.ina_meta_branch_step_pre_nodeshared <- .ina_meta_branch_step
.ina_meta_branch_step <- function(...) {
  z <- .ina_meta_branch_step_pre_nodeshared(...)
  z$SharedInfoNode <- seq_along(z$q0)
  z$approximation <- paste(z$approximation,
    "; node-level detection/retention/SEAM is integrated as a shared event within each lineage")
  z
}

.ina_transition_branch_step_pre_nodeshared <- .ina_transition_branch_step
.ina_transition_branch_step <- function(...) {
  z <- .ina_transition_branch_step_pre_nodeshared(...)
  # Transition base types are node-major, stage-minor throughout the analytical
  # source. Recover S from the repeated detection/mortality structure by using
  # the explicit Nstages argument supplied positionally/named to the wrapped
  # function.
  dots <- list(...)
  S <- if (!is.null(dots$Nstages)) dots$Nstages else if (length(dots) >= 2L) dots[[2L]] else NULL
  if (is.null(S) || length(S) != 1L || !is.finite(S) || S < 1) return(z)
  S <- as.integer(S)
  B <- length(z$q0); n <- B %/% S
  if (n * S == B) z$SharedInfoNode <- rep(seq_len(n), each = S)
  z$approximation <- paste(z$approximation,
    "; node-level detection/retention/SEAM is integrated as a shared event within each lineage")
  z
}
# Integrated helpers for the simulator-faithful MLU one-parent branching PGF.
# These build on the analytical helpers defined earlier in this same source file.

.ina_mlu_multinom_factor_expect <- function(mS, mL, pS, pL, node_factor) {
  mS <- as.integer(mS); mL <- as.integer(mL)
  pS <- as.numeric(pS); pL <- as.numeric(pL); n <- length(pS)
  if (length(pL) != n) stop("pS/pL length mismatch")
  # Drop destinations that receive no SDD/LDD mass and contribute a neutral
  # factor at zero arrivals. This is exact and makes sparse landscapes much
  # cheaper without changing the offspring distribution.
  original_index <- seq_len(n)
  f0 <- vapply(original_index, function(j) node_factor(j, 0L), numeric(1))
  keep <- pS > 0 | pL > 0 | abs(f0 - 1) > 1e-14
  if (!all(keep)) {
    original_index <- original_index[keep]
    pS <- pS[keep]; pL <- pL[keep]; n <- length(original_index)
  }
  if (n < 1L) return(1)
  if (mS > 0L && abs(sum(pS) - 1) > 1e-9) stop("pS must sum to 1")
  if (mL > 0L && abs(sum(pL) - 1) > 1e-9) stop("pL must sum to 1")
  # Precompute node factors for all possible total arrivals.
  Ftab <- lapply(seq_len(n), function(j)
    vapply(0:(mS + mL), function(a) node_factor(original_index[j], a), numeric(1)))
  tailS <- rev(cumsum(rev(pS))); tailL <- rev(cumsum(rev(pL)))
  memo <- new.env(hash = TRUE, parent = emptyenv())
  rec <- function(j, rS, rL) {
    key <- paste(j, rS, rL, sep = ":")
    if (exists(key, memo, inherits = FALSE)) return(get(key, memo, inherits = FALSE))
    if (j == n) {
      val <- Ftab[[j]][rS + rL + 1L]
      assign(key, val, memo); return(val)
    }
    ps <- if (rS == 0L) 0 else if (tailS[j] > 0) pS[j] / tailS[j] else 0
    pl <- if (rL == 0L) 0 else if (tailL[j] > 0) pL[j] / tailL[j] else 0
    va <- 0:rS; vb <- 0:rL
    pa <- if (rS == 0L) 1 else dbinom(va, rS, ps)
    pb <- if (rL == 0L) 1 else dbinom(vb, rL, pl)
    ans <- 0
    for (ia in seq_along(va)) for (ib in seq_along(vb)) {
      w <- pa[ia] * pb[ib]
      if (w == 0) next
      ans <- ans + w * Ftab[[j]][va[ia] + vb[ib] + 1L] *
        rec(j + 1L, rS - va[ia], rL - vb[ib])
    }
    assign(key, ans, memo); ans
  }
  rec(1L, mS, mL)
}

.ina_mlu_branch_step <- function(SDDprob, LDDprob = 0, LDDrate = 0,
                                 EnvEstabProb = 1, Survival = 1, K,
                                 PropaguleProduction, PropaguleEstablishment,
                                 DetectionProb = 0, ManageProb = 0,
                                 MortalityProb = 0, SpreadReduction = 0,
                                 SEAM = NULL, InfoRetentionProb = 1,
                                 FecundityReduction = 0) {
  SDD <- as.matrix(SDDprob); n <- nrow(SDD)
  if (ncol(SDD) != n) stop("SDDprob must be square")
  LDD <- .ina_mat(LDDprob, n, "LDDprob")
  K <- as.matrix(K); if (nrow(K) != n) stop("K rows must equal nodes")
  L <- ncol(K); B <- n * L
  if (any(!is.finite(K)) || any(K < 0) || any(abs(K - round(K)) > 1e-10))
    stop("Exact MLU branching PGF currently requires non-negative integer K")
  r <- as.numeric(LDDrate)
  if (length(r) != 1L || !is.finite(r) || r < 0 || r > 1) stop("LDDrate must be scalar in [0,1]")
  env <- .ina_recycle(EnvEstabProb, n, "EnvEstabProb")
  surv <- .ina_recycle(Survival, n, "Survival")
  prod <- .ina_recycle(PropaguleProduction, n, "PropaguleProduction")
  pe <- .ina_recycle(PropaguleEstablishment, n, "PropaguleEstablishment")
  A <- .mlu_matrix(ManageProb, n, L, "ManageProb")
  Mort <- .mlu_matrix(MortalityProb, n, L, "MortalityProb")
  Gr <- .mlu_matrix(SpreadReduction, n, L, "SpreadReduction")
  Fec <- .mlu_matrix(FecundityReduction, n, L, "FecundityReduction")
  D <- .mlu_matrix(DetectionProb, n, L, "DetectionProb")
  ir <- .ina_recycle(InfoRetentionProb, n, "InfoRetentionProb")
  C <- if (is.null(SEAM) || length(SEAM) == 1L) matrix(0, n, n) else as.matrix(SEAM)
  if (!all(dim(C) == c(n, n))) stop("SEAM must be nodes x nodes")
  diag(C) <- 0
  if (any(c(A, Mort, Gr, Fec, D, ir, C) < 0) || any(c(A, Mort, Gr, Fec, D, ir, C) > 1))
    stop("MLU probabilities/effects must be in [0,1]")
  type_node <- rep(seq_len(n), each = L); type_lu <- rep(seq_len(L), times = n)
  q0 <- surv[type_node]
  q1 <- q0 * (1 - Mort[cbind(type_node, type_lu)])
  list(EvaluatorType = "MLUExact", n = n, L = L, B = B,
       SDD = SDD, LDD = LDD, LDDrate = r, K = K,
       alpha = env * pe, Survival = surv, Production = prod,
       Adoption = A, Mortality = Mort, Spread = Gr, Fecundity = Fec,
       Detection = as.vector(t(D)), DetectionMatrix = D,
       Retention = ir[type_node], RetentionNode = ir,
       Transfer = .ina_node_transfer_types(C, type_node), SEAMNode = C,
       type_node = type_node, type_lu = type_lu,
       q0 = q0, q1 = q1, adoption = A[cbind(type_node, type_lu)],
       approximation = paste(
         "simulator-faithful one-parent MLU PGF with exact integer SDD/LDD",
         "multinomial allocation and exact finite-slot land-use recruitment"))
}

.ina_mlu_count_factor <- function(step, src, zmat, j, arrivals) {
  i <- step$type_node[src]; l <- step$type_lu[src]
  free <- step$K[j, ]
  if (j == i) free[l] <- max(0, free[l] - 1)
  p <- if (arrivals <= 0 || step$alpha[j] <= 0) 0 else 1 - exp(-step$alpha[j] * arrivals)
  rec <- prod((1 - p + p * zmat[j, ])^free)
  if (j == i) rec <- rec * zmat[j, l]
  rec
}

.ina_mlu_counts_pgf <- function(step, src, managed, node_factor) {
  i <- step$type_node[src]; l <- step$type_lu[src]
  rs <- rowSums(step$SDD); rl <- rowSums(step$LDD)
  lambda <- step$Production[i] * (1 - step$Fecundity[i, l] * managed)
  supp <- .meta_poisson_support(lambda)
  ms <- floor(supp$k * (1 - step$LDDrate) * rs[i])
  ml <- floor(supp$k * step$LDDrate * (1 - step$Spread[i, l] * managed) * rl[i])
  keys <- paste(ms, ml, sep = ":"); out <- 0
  pS <- if (rs[i] > 0) step$SDD[i, ] / rs[i] else c(1, rep(0, step$n - 1L))
  pL <- if (rl[i] > 0) step$LDD[i, ] / rl[i] else c(1, rep(0, step$n - 1L))
  for (key in unique(keys)) {
    ids <- which(keys == key)
    ab <- as.integer(strsplit(key, ":", fixed = TRUE)[[1L]])
    v <- .ina_mlu_multinom_factor_expect(ab[1L], ab[2L], pS, pL, node_factor)
    out <- out + sum(supp$p[ids]) * v
  }
  .ina_branch_clip(out)
}

.ina_mlu_branch_eval_one <- function(step, src, managed, informed_source,
                                     qU = NULL, qH = NULL, qbase = NULL) {
  i <- step$type_node[src]; l <- step$type_lu[src]
  surv <- step$Survival[i] * (1 - step$Mortality[i, l] * managed)
  if (surv <= 0) return(1)
  n <- step$n; L <- step$L
  if (!is.null(qbase)) {
    Z <- matrix(qbase, n, L, byrow = TRUE)
    nf <- function(j, a) .ina_mlu_count_factor(step, src, Z, j, a)
  } else {
    U <- matrix(qU, n, L, byrow = TRUE)
    H <- matrix(qH, n, L, byrow = TRUE)
    D <- step$DetectionMatrix
    noUZ <- (1 - D) * U; noHZ <- (1 - D) * H
    nf <- function(j, a) {
      GH <- .ina_mlu_count_factor(step, src, H, j, a)
      noU <- .ina_mlu_count_factor(step, src, noUZ, j, a)
      noH <- .ina_mlu_count_factor(step, src, noHZ, j, a)
      GU <- noU + GH - noH
      pinfo <- if (!informed_source) 0 else if (j == i)
        step$RetentionNode[i] else step$SEAMNode[i, j]
      .ina_branch_clip(pinfo * GH + (1 - pinfo) * GU)
    }
  }
  (1 - surv) + surv * .ina_mlu_counts_pgf(step, src, managed, nf)
}

.ina_mlu_branch_step_eval <- function(step, qnext,
                                      mode = c("none", "all_informed", "dynamic")) {
  mode <- match.arg(mode); B <- step$B
  if (mode != "dynamic") {
    q <- .ina_branch_clip(qnext); if (length(q) != B) stop("MLU q length mismatch")
    out <- numeric(B)
    for (src in seq_len(B)) {
      i <- step$type_node[src]; l <- step$type_lu[src]
      f0 <- .ina_mlu_branch_eval_one(step, src, 0L, FALSE, qbase = q)
      if (mode == "all_informed") {
        f1 <- .ina_mlu_branch_eval_one(step, src, 1L, TRUE, qbase = q)
        a <- step$Adoption[i, l]; out[src] <- (1 - a) * f0 + a * f1
      } else out[src] <- f0
    }
    return(.ina_branch_clip(out))
  }
  q <- .ina_branch_clip(qnext); if (length(q) != 2L * B) stop("MLU dynamic q length mismatch")
  qU <- q[seq_len(B)]; qH <- q[B + seq_len(B)]
  outU <- outH <- numeric(B)
  for (src in seq_len(B)) {
    i <- step$type_node[src]; l <- step$type_lu[src]
    outU[src] <- .ina_mlu_branch_eval_one(step, src, 0L, FALSE, qU, qH)
    f0 <- .ina_mlu_branch_eval_one(step, src, 0L, TRUE, qU, qH)
    f1 <- .ina_mlu_branch_eval_one(step, src, 1L, TRUE, qU, qH)
    a <- step$Adoption[i, l]; outH[src] <- (1 - a) * f0 + a * f1
  }
  .ina_branch_clip(c(outU, outH))
}

.ina_mlu_initial_extinction <- function(q, mode, InitialState, InitialInfo,
                                        ApplyInitialDetection, DetectionProb,
                                        n, L, Ntimesteps) {
  X <- if (is.matrix(InitialState)) as.matrix(InitialState) else
    matrix(as.numeric(InitialState), nrow = n, ncol = L, byrow = TRUE)
  if (!all(dim(X) == c(n, L))) stop("InitialState must be nodes x land uses")
  B <- n * L
  powerprod <- function(qq, xx) {
    if (any(qq == 0 & xx > 0)) return(0)
    exp(sum(xx * log(pmax(qq, .Machine$double.xmin))))
  }
  if (mode != "dynamic") return(powerprod(.ina_branch_clip(q), as.vector(t(X))))
  D1 <- .ina_slice_mlu(DetectionProb, 1L, n, L, Ntimesteps, "DetectionProb")
  Dm <- if (length(D1) == 1L) matrix(D1, n, L) else if (length(D1) == L)
    matrix(rep(D1, each = n), n, L) else as.matrix(D1)
  pd <- if (ApplyInitialDetection) .ina_initial_detection_mlu(X, Dm) else rep(0, n)
  p0 <- .ina_recycle(InitialInfo, n, "InitialInfo"); pH <- .ina_branch_clip(p0 + (1 - p0) * pd)
  q <- .ina_branch_clip(q); qU <- q[seq_len(B)]; qH <- q[B + seq_len(B)]
  out <- 1
  for (i in seq_len(n)) {
    ids <- (i - 1L) * L + seq_len(L); xi <- X[i, ]
    pu <- powerprod(qU[ids], xi); ph <- powerprod(qH[ids], xi)
    out <- out * ((1 - pH[i]) * pu + pH[i] * ph)
  }
  out
}

# Exact one-parent MLU programmed-information evaluator. Information is shared
# at node level within the focal lineage; information retained at pest-free
# nodes remains outside independent-lineage branching.
.ina_mlu_prog_node_factor <- function(step, src, q, layout, j, arrivals,
                                      pre_state = c("U", "X", "Hnext", "Overflow"),
                                      pre_prob = 0, hnext_index = NULL) {
  pre_state <- match.arg(pre_state)
  B <- step$B; n <- step$n; L <- step$L
  ids <- (j - 1L) * L + seq_len(L)
  D <- step$DetectionMatrix[j, ]
  zU <- q[layout$U[ids]]
  zH1 <- q[layout$H[[1L]][ids]]
  zNoU <- (1 - D) * zU
  zNoH1 <- (1 - D) * zH1
  noU <- .ina_mlu_count_factor(step, src, matrix(replace(rep(1, n*L), ids, zNoU), n, L, byrow=TRUE), j, arrivals)
  GH1 <- .ina_mlu_count_factor(step, src, matrix(replace(rep(1, n*L), ids, zH1), n, L, byrow=TRUE), j, arrivals)
  noH1 <- .ina_mlu_count_factor(step, src, matrix(replace(rep(1, n*L), ids, zNoH1), n, L, byrow=TRUE), j, arrivals)
  if (pre_prob <= 0 || pre_state == "U") return(.ina_branch_clip(noU + GH1 - noH1))
  zPre <- switch(pre_state,
    X = q[layout$X[ids]],
    Hnext = q[hnext_index[ids]],
    Overflow = q[layout$Overflow[ids]])
  noPre <- .ina_mlu_count_factor(step, src,
    matrix(replace(rep(1, n*L), ids, (1-D)*zPre), n,L,byrow=TRUE), j, arrivals)
  .ina_branch_clip(pre_prob * noPre + (1-pre_prob) * noU + GH1 - noH1)
}

.ina_mlu_programmed_branch_step_eval <- function(step, qnext, layout,
                                                  persistence_by_node) {
  q <- .ina_branch_clip(qnext); if (length(q) != layout$size) stop("MLU programmed q length mismatch")
  n <- step$n; L <- step$L; B <- step$B
  Kp <- .ina_recycle(persistence_by_node, n, "InfoPersistenceSteps")

  # Evaluate the biological one-parent branch with a node-factor generator.
  branch_with_nodes <- function(src, managed, node_factor_builder) {
    i <- step$type_node[src]; l <- step$type_lu[src]
    surv <- step$Survival[i] * (1 - step$Mortality[i,l] * managed)
    if (surv <= 0) return(1)
    nf <- function(j,a) node_factor_builder(j,a)
    (1-surv) + surv * .ina_mlu_counts_pgf(step,src,managed,nf)
  }

  # Uninformed source: no management or transfer; detection creates H1.
  out <- numeric(layout$size)
  for (src in seq_len(B)) {
    nodefun <- function(j,a) .ina_mlu_prog_node_factor(
      step,src,q,layout,j,a,pre_state="U",pre_prob=0)
    out[layout$U[src]] <- branch_with_nodes(src,0L,nodefun)
  }

  eval_informed <- function(src, kind, age = NA_integer_) {
    i <- step$type_node[src]; l <- step$type_lu[src]
    # Source-node state if no new local detection occurs at the end of step.
    source_state <- "U"; source_prob <- 0; hnext <- NULL
    if (identical(kind,"X")) {
      if (is.na(Kp[i])) { source_state <- "X"; source_prob <- step$RetentionNode[i] }
    } else if (identical(kind,"Overflow")) {
      if (is.na(Kp[i])) { source_state <- "Overflow"; source_prob <- step$RetentionNode[i] }
    } else {
      next_ids <- if (age < layout$max_age) layout$H[[age+1L]] else layout$Overflow
      if (is.na(Kp[i])) {
        source_state <- if (age < layout$max_age) "Hnext" else "Overflow"
        source_prob <- step$RetentionNode[i]; hnext <- next_ids
      } else if (age < Kp[i]) {
        source_state <- if (age < layout$max_age) "Hnext" else "Overflow"
        source_prob <- 1; hnext <- next_ids
      }
    }
    nodefun <- function(j,a) {
      if (j == i) {
        .ina_mlu_prog_node_factor(step,src,q,layout,j,a,
          pre_state=source_state,pre_prob=source_prob,hnext_index=hnext)
      } else {
        # A SEAM refresh gives information without a local-evidence clock (X).
        .ina_mlu_prog_node_factor(step,src,q,layout,j,a,
          pre_state="X",pre_prob=step$SEAMNode[i,j])
      }
    }
    f0 <- branch_with_nodes(src,0L,nodefun)
    f1 <- branch_with_nodes(src,1L,nodefun)
    a <- step$Adoption[i,l]
    (1-a)*f0+a*f1
  }

  for (src in seq_len(B)) out[layout$X[src]] <- eval_informed(src,"X")
  for (age in seq_len(layout$max_age))
    for (src in seq_len(B)) out[layout$H[[age]][src]] <- eval_informed(src,"H",age)
  for (src in seq_len(B)) out[layout$Overflow[src]] <- eval_informed(src,"Overflow")
  .ina_branch_clip(out)
}

.ina_mlu_programmed_initial_extinction <- function(q, InitialState, InitialInfo,
                                                    ApplyInitialDetection,
                                                    DetectionProb, n, L,
                                                    Ntimesteps, layout) {
  X <- if (is.matrix(InitialState)) as.matrix(InitialState) else
    matrix(as.numeric(InitialState),nrow=n,ncol=L,byrow=TRUE)
  D1 <- .ina_slice_mlu(DetectionProb,1L,n,L,Ntimesteps,"DetectionProb")
  Dm <- if(length(D1)==1L)matrix(D1,n,L) else if(length(D1)==L)
    matrix(rep(D1,each=n),n,L) else as.matrix(D1)
  pd <- if(ApplyInitialDetection).ina_initial_detection_mlu(X,Dm) else rep(0,n)
  p0 <- .ina_branch_clip(.ina_recycle(InitialInfo,n,"InitialInfo"))
  powerprod <- function(ids,xx){qq<-q[ids];if(any(qq==0 & xx>0))return(0);exp(sum(xx*log(pmax(qq,.Machine$double.xmin))))}
  out<-1
  for(i in seq_len(n)){
    ids<-(i-1L)*L+seq_len(L);xi<-X[i,]
    pu<-(1-p0[i])*(1-pd[i]);px<-p0[i]*(1-pd[i]);ph<-pd[i]
    out<-out*(pu*powerprod(layout$U[ids],xi)+px*powerprod(layout$X[ids],xi)+ph*powerprod(layout$H[[1L]][ids],xi))
  }
  out
}

###############################################################################
### Exact one-parent MLU PGF integration
###############################################################################

.ina_branch_step_eval_pre_mlu_exact <- .ina_branch_step_eval
.ina_branch_step_eval <- function(step, qnext,
                                  mode = c("none", "all_informed", "dynamic")) {
  mode <- match.arg(mode)
  if (identical(step$EvaluatorType, "MLUExact"))
    return(.ina_mlu_branch_step_eval(step, qnext, mode))
  .ina_branch_step_eval_pre_mlu_exact(step, qnext, mode)
}

.ina_mlu_exact_tractability <- function(Ntimesteps, SDDprob, LDDprob,
                                                LDDrate, PropaguleProduction,
                                                n, L) {
  B <- n * L
  scale <- sqrt(12 / max(12, B * as.integer(Ntimesteps)))
  worst <- NULL
  for (tt in seq_len(as.integer(Ntimesteps))) {
    SDDt <- .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob")
    LDDt <- .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob")
    Pt <- .ina_recycle(.ina_slice_node(PropaguleProduction, tt, n, Ntimesteps,
                                       "PropaguleProduction"), n,
                       "PropaguleProduction")
    rs <- rowSums(as.matrix(SDDt)); rl <- rowSums(.ina_mat(LDDt, n, "LDDprob"))
    for (i in seq_len(n)) {
      supp <- .meta_poisson_support(Pt[i])
      kmax <- max(supp$k)
      mS <- floor(kmax * (1 - LDDrate) * rs[i])
      mL <- floor(kmax * LDDrate * rl[i])
      active <- (as.matrix(SDDt)[i, ] > 0) |
                (.ina_mat(LDDt, n, "LDDprob")[i, ] > 0)
      # The source node must remain represented even if it receives no
      # dispersal mass because the surviving parent contributes there.
      active[i] <- TRUE
      d <- sum(active)
      base_limit <- if (d <= 1L) Inf else if (d == 2L) 100 else if (d == 3L) 35 else
                    if (d == 4L) 22 else if (d <= 6L) 15 else if (d <= 10L) 10 else
                    if (d <= 15L) 7 else 4
      limit <- max(2, floor(base_limit * scale))
      count <- mS + mL
      if (is.null(worst) || count / limit > worst$ratio)
        worst <- list(ratio = count / limit, timestep = tt, source = i,
                      destinations = d, max_dispersers = count, limit = limit,
                      lambda = Pt[i])
    }
  }
  if (!is.null(worst) && is.finite(worst$ratio) && worst$ratio > 1) {
    return(list(ok = FALSE, reason = paste0(
      "exact MLU PGF integration is computationally large for this landscape ",
      "(timestep ", worst$timestep, ", source node ", worst$source,
      ": ", worst$destinations, " active spatial destinations and up to ",
      worst$max_dispersers, " low-count dispersers in the retained Poisson support; ",
      "automatic exact-integration limit ", worst$limit, ")")))
  }
  list(ok = TRUE)
}

.ina_mlu_dynamic_pgf_extinction <- function(
    Ntimesteps, InitialState, InitialInfo, InformationMode,
    ApplyInitialDetection, SDDprob, LDDprob, LDDrate, EnvEstabProb,
    Survival, K, PropaguleProduction, PropaguleEstablishment,
    DetectionProb, ManageProb, MortalityProb, SpreadReduction, SEAM,
    InfoRetentionProb, InfoPersistenceSteps, FecundityReduction,
    ExtinctionGenerations) {
  Ntimesteps <- as.integer(Ntimesteps)
  if (is.null(K) || is.null(PropaguleProduction)) return(NULL)
  SDD0 <- if (length(dim(SDDprob)) == 3L) SDDprob[,,1L] else as.matrix(SDDprob)
  n <- nrow(SDD0)
  K0 <- if (length(dim(K)) == 3L) K[,,1L] else as.matrix(K)
  L <- ncol(K0); B <- n * L
  mode <- .ina_info_mode(InformationMode, ManageProb, InitialInfo,
                         DetectionProb, SEAM, InfoRetentionProb,
                         InfoPersistenceSteps)
  if (!(mode %in% c("none", "all_informed", "dynamic"))) return(NULL)
  tract <- .ina_mlu_exact_tractability(Ntimesteps, SDDprob, LDDprob, LDDrate,
                                       PropaguleProduction, n, L)
  if (!isTRUE(tract$ok))
    return(list(Unavailable = TRUE, Reason = tract$reason))
  steps <- vector("list", Ntimesteps)
  for (tt in seq_len(Ntimesteps)) {
    steps[[tt]] <- tryCatch(.ina_mlu_branch_step(
      .ina_slice_connection(SDDprob, tt, n, Ntimesteps, "SDDprob"),
      .ina_slice_connection(LDDprob, tt, n, Ntimesteps, "LDDprob"),
      LDDrate,
      .ina_slice_node(EnvEstabProb, tt, n, Ntimesteps, "EnvEstabProb"),
      .ina_slice_node(Survival, tt, n, Ntimesteps, "Survival"),
      .ina_slice_K_mlu(K, tt, n, L, Ntimesteps),
      .ina_slice_node(PropaguleProduction, tt, n, Ntimesteps, "PropaguleProduction"),
      .ina_slice_node(PropaguleEstablishment, tt, n, Ntimesteps, "PropaguleEstablishment"),
      .ina_slice_mlu(DetectionProb, tt, n, L, Ntimesteps, "DetectionProb"),
      .ina_slice_mlu(ManageProb, tt, n, L, Ntimesteps, "ManageProb"),
      .ina_slice_mlu(MortalityProb, tt, n, L, Ntimesteps, "MortalityProb"),
      .ina_slice_mlu(SpreadReduction, tt, n, L, Ntimesteps, "SpreadReduction"),
      SEAM,
      .ina_slice_node(InfoRetentionProb, tt, n, Ntimesteps, "InfoRetentionProb"),
      .ina_slice_mlu(FecundityReduction, tt, n, L, Ntimesteps, "FecundityReduction")),
      error = function(e) structure(list(error = conditionMessage(e)), class = "MLUPGFError"))
    if (inherits(steps[[tt]], "MLUPGFError"))
      return(list(Unavailable = TRUE, Reason = steps[[tt]]$error))
  }

  persistence_requested <- any(!is.na(as.numeric(InfoPersistenceSteps)))
  if (persistence_requested && mode == "dynamic") {
    max_age <- .ina_programmed_global_max_age(InfoPersistenceSteps)
    layout <- .ina_programmed_layout(B, max_age)
    persistence <- lapply(seq_len(Ntimesteps), function(tt)
      .ina_recycle(.ina_slice_node(InfoPersistenceSteps, tt, n, Ntimesteps,
                                   "InfoPersistenceSteps"), n,
                   "InfoPersistenceSteps"))
    qh <- rep(0, layout$size)
    for (tt in rev(seq_len(Ntimesteps)))
      qh <- .ina_mlu_programmed_branch_step_eval(
        steps[[tt]], qh, layout, persistence[[tt]])
    static <- .ina_branch_static(steps) &&
      (length(persistence) <= 1L || all(vapply(persistence[-1L], function(x)
        isTRUE(all.equal(x, persistence[[1L]], tolerance = 0)), logical(1))))
    qe <- NULL
    if (static) {
      qe <- rep(0, layout$size)
      for (gg in seq_len(ExtinctionGenerations)) {
        old <- qe
        qe <- .ina_mlu_programmed_branch_step_eval(
          steps[[1L]], qe, layout, persistence[[1L]])
        if (max(abs(qe - old)) < 1e-12) break
      }
    }
    ext_h <- .ina_mlu_programmed_initial_extinction(
      qh, InitialState, InitialInfo, ApplyInitialDetection, DetectionProb,
      n, L, Ntimesteps, layout)
    ext_e <- if (is.null(qe)) NA_real_ else .ina_mlu_programmed_initial_extinction(
      qe, InitialState, InitialInfo, ApplyInitialDetection, DetectionProb,
      n, L, Ntimesteps, layout)
    return(list(mode = mode, qh = qh, qe = qe, static = static,
      ProbabilityByHorizon = ext_h,
      BranchingFadeoutProbability = ext_e,
      Method = "simulator-faithful MLU programmed-information one-parent branching PGF",
      note = paste(
        "MLU extinction retains integer SDD/LDD multinomial allocation, finite",
        "land-use recruitment slots and node-shared detection/SEAM within each",
        "rare lineage. Shared management/evidence events across multiple",
        "coexisting lineages, including information retained by pest-free nodes,",
        "remain higher-order finite-population effects.")))
  }

  qh <- .ina_branch_horizon(steps, mode)
  static <- .ina_branch_static(steps)
  qe <- if (static) .ina_branch_eventual(steps[[1L]], mode,
                                          ExtinctionGenerations) else NULL
  ext_h <- .ina_mlu_initial_extinction(
    qh, mode, InitialState, InitialInfo, ApplyInitialDetection,
    DetectionProb, n, L, Ntimesteps)
  ext_e <- if (is.null(qe)) NA_real_ else .ina_mlu_initial_extinction(
    qe, mode, InitialState, InitialInfo, ApplyInitialDetection,
    DetectionProb, n, L, Ntimesteps)
  list(mode = mode, qh = qh, qe = qe, static = static,
       ProbabilityByHorizon = ext_h,
       BranchingFadeoutProbability = ext_e,
       Method = if (mode == "dynamic")
         "time-inhomogeneous node-information-aware simulator-faithful MLU one-parent branching PGF"
       else "time-inhomogeneous simulator-faithful MLU one-parent branching PGF",
       note = if (mode == "dynamic") paste(
         "MLU node-level detection, information retention and direct SEAM are",
         "integrated as shared events within each rare lineage. Information",
         "retained by pest-free nodes and management/information shared across",
         "coexisting independent lineages remain higher-order effects.") else
         "MLU finite-horizon extinction uses the simulator-faithful one-parent offspring distribution rather than a mean-matched Poisson fallback.")
}

INApestAnalytical_pre_mlu_exact_pgf <- INApestAnalytical
INApestAnalytical <- function(...) {
  args <- list(...)
  res <- do.call(INApestAnalytical_pre_mlu_exact_pgf, args)
  Model <- if (!is.null(args$Model)) args$Model else "INApest"
  if (length(Model) > 1L) Model <- Model[1L]
  if (!identical(Model, "INApestMetaMultipleLandUse")) return(res)

  fm <- formals(INApestAnalytical_pre_transition_mean)
  getarg <- function(nm) {
    if (!is.null(args[[nm]])) return(args[[nm]])
    z <- fm[[nm]]
    if (is.null(z)) return(NULL)
    eval(z, envir = parent.frame())
  }
  improved <- .ina_mlu_dynamic_pgf_extinction(
    getarg("Ntimesteps"), getarg("InitialState"), getarg("InitialInfo"),
    getarg("InformationMode"), getarg("ApplyInitialDetection"),
    getarg("SDDprob"), getarg("LDDprob"), getarg("LDDrate"),
    getarg("EnvEstabProb"), getarg("Survival"), getarg("K"),
    getarg("PropaguleProduction"), getarg("PropaguleEstablishment"),
    getarg("DetectionProb"), getarg("ManageProb"), getarg("MortalityProb"),
    getarg("SpreadReduction"), getarg("SEAM"), getarg("InfoRetentionProb"),
    getarg("InfoPersistenceSteps"), getarg("FecundityReduction"),
    getarg("ExtinctionGenerations"))
  if (is.null(improved)) return(res)
  if (isTRUE(improved$Unavailable)) {
    res$Diagnostics <- unique(c(res$Diagnostics,
      paste0("Simulator-faithful MLU branching PGF unavailable: ", improved$Reason,
             ". Existing MLU extinction fallback retained.")))
    return(res)
  }
  programmed_here <- any(!is.na(as.numeric(getarg("InfoPersistenceSteps"))))
  fecundity_here <- .ina_fr_nonzero(getarg("FecundityReduction"))
  use_improved <- identical(improved$mode, "dynamic") || programmed_here ||
                  !isTRUE(improved$static) || fecundity_here
  # Preserve the established static, non-fecundity MLU result exactly. The new
  # PGF takes over only where it adds capability requested in this extension.
  if (!use_improved) return(res)
  old_info_fallback_diag <- "For information-limited Meta/Transition/MLU models, branching extinction/escape may use a mean-matched multitype fallback because node-level shared information induces correlations not represented by independent individual types."
  res$Diagnostics <- res$Diagnostics[res$Diagnostics != old_info_fallback_diag]
  res$Extinction$Method <- improved$Method
  res$Extinction$ProbabilityByHorizon <- improved$ProbabilityByHorizon
  res$Extinction$BranchingFadeoutProbability <- improved$BranchingFadeoutProbability
  res$Diagnostics <- unique(c(res$Diagnostics, improved$note,
    if (!improved$static)
      "Finite-horizon MLU extinction is composed backward through timestep-specific PGFs, so changing habitat/connectivity/management does not require the generic mean-matched extinction fallback." else NULL))
  res
}


###############################################################################
### Vertebrate analytical extensions: animal populations and generic state process
### Integrated into the unified analytical source, 24 August 2026
###############################################################################

# INApest analytical extension for vertebrate models
# Prototype analytical companion, 21 August 2026
#
# Scope deliberately mirrors the existing INApest analytical philosophy:
#   1. GrowthRate        - low-density landscape multiplier
#   2. EndPopulation     - finite-horizon expected population
#   3. EscapeProbability - finite-horizon probability of >=1 successful escape
#
# Two analytical variants are provided:
#   * animal   : demographic stage x node branching process
#   * pathogen : generic host-state process, linearised about a supplied host state
#
# This file does NOT silently approximate arbitrary Vertebrate Birth/HomeRange/
# Control/Interaction hooks. Such modules can supply analytical equivalents via
# explicit matrices/functions; otherwise Diagnostics records the limitation.

.ivanal_clip01 <- function(x) pmin(pmax(x, 0), 1)

.ivanal_spectral_radius <- function(A) {
  if (!length(A)) return(0)
  max(Mod(eigen(A, only.values = TRUE)$values))
}

.ivanal_resolve_time <- function(x, t, T, name, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(NULL)
    stop(name, " may not be NULL")
  }
  if (is.function(x)) return(x(timestep = t))
  d <- dim(x)
  if (is.null(d) || length(d) <= 2L) return(x)
  if (length(d) == 3L) {
    if (d[3] != T) stop(name, " third dimension must equal Ntimesteps")
    return(x[, , t, drop = FALSE][, , 1])
  }
  if (length(d) == 4L) {
    if (d[4] != T) stop(name, " fourth dimension must equal Ntimesteps")
    return(x[, , , t, drop = FALSE][, , , 1])
  }
  stop(name, " has unsupported dimensions")
}

.ivanal_expand_node_stage <- function(x, n_nodes, n_stages, t = 1L,
                                      T = 1L, name = deparse(substitute(x)),
                                      default = NULL) {
  if (is.null(x)) {
    if (!is.null(default)) return(matrix(default, n_nodes, n_stages))
    stop(name, " may not be NULL")
  }
  if (is.function(x)) x <- x(timestep = t)
  d <- dim(x)
  if (length(d) == 3L) {
    if (!all(d[1:2] == c(n_nodes, n_stages)) || d[3] != T)
      stop(name, " 3D array must be nodes x stages x Ntimesteps")
    return(x[, , t, drop = FALSE][, , 1])
  }
  if (is.matrix(x)) {
    if (!all(dim(x) == c(n_nodes, n_stages)))
      stop(name, " matrix must be nodes x stages")
    return(x)
  }
  if (length(x) == 1L) return(matrix(x, n_nodes, n_stages))
  if (length(x) == n_nodes && length(x) == n_stages)
    stop(name, " vector is ambiguous because n_nodes equals n_stages; supply a matrix")
  if (length(x) == n_nodes)
    return(matrix(rep(x, n_stages), nrow = n_nodes, ncol = n_stages))
  if (length(x) == n_stages)
    return(matrix(rep(x, each = n_nodes), nrow = n_nodes, ncol = n_stages))
  stop(name, " must be scalar, node vector, stage vector, nodes x stages matrix, or nodes x stages x time array")
}

.ivanal_transition_for_node <- function(Transition, node, t, n_nodes, n_stages, T) {
  if (is.function(Transition)) {
    A <- Transition(node = node, timestep = t)
  } else if (is.list(Transition)) {
    if (length(Transition) != n_nodes) stop("Transition list must contain one matrix per node")
    A <- Transition[[node]]
  } else {
    d <- dim(Transition)
    if (length(d) == 2L) {
      A <- Transition
    } else if (length(d) == 3L && d[3] == T) {
      A <- Transition[, , t, drop = FALSE][, , 1]
    } else if (length(d) == 4L && d[3] == n_nodes && d[4] == T) {
      A <- Transition[, , node, t, drop = FALSE][, , 1, 1]
    } else {
      stop("Transition must be a stage matrix, node list, stage x stage x time array, or stage x stage x node x time array")
    }
  }
  if (!all(dim(A) == c(n_stages, n_stages))) stop("Resolved Transition matrix has wrong dimensions")
  A
}

.ivanal_movement_for_stage <- function(x, stage, t, T, n_nodes, n_stages, name) {
  if (is.null(x)) return(NULL)
  z <- if (is.list(x)) {
    if (length(x) != n_stages - 1L) stop(name, " list must have length Nstages - 1")
    x[[stage]]
  } else x
  if (is.null(z)) return(NULL)
  if (is.function(z)) z <- z(stage = stage, timestep = t)
  d <- dim(z)
  if (length(d) == 2L) {
    if (!all(d == c(n_nodes, n_nodes))) stop(name, " matrix must be nodes x nodes")
    return(z)
  }
  if (length(d) == 3L && d[3] == T) {
    if (!all(d[1:2] == c(n_nodes, n_nodes))) stop(name, " array must be nodes x nodes x Ntimesteps")
    return(z[, , t, drop = FALSE][, , 1])
  }
  stop(name, " entries must be node matrices or node x node x time arrays")
}

.ivanal_mix_movement <- function(Psdd, Pldd, lddrate, n_nodes) {
  if (is.null(Psdd) && is.null(Pldd)) return(NULL)
  if (is.null(Psdd)) Psdd <- matrix(0, n_nodes, n_nodes)
  if (is.null(Pldd)) Pldd <- matrix(0, n_nodes, n_nodes)
  r <- .ivanal_clip01(lddrate)
  (1 - r) * Psdd + r * Pldd
}

# Build the low-density mean operator and one-parent branching ingredients for
# one timestep. Type ordering is node-major within stage:
#   type(node, stage) = (node - 1) * Nstages + stage.
.ivanal_animal_node_step <- function(
    t, T, Nstages, Transition, n_nodes,
    SDDprob, LDDprob = NULL, LDDrate = 0,
    TransitionSDDprob = NULL, TransitionLDDprob = NULL,
    TransitionLDDrate = 0,
    PropaguleEstablishment = 1, EnvEstabProb = 1,
    MortalityProb = 0, ManagementExposure = 0,
    FecundityReduction = 0,
    OutsideEstablishmentProb = 1) {

  S <- Nstages
  K <- n_nodes * S
  idx <- function(node, stage) (node - 1L) * S + stage

  sdd <- .ivanal_resolve_time(SDDprob, t, T, "SDDprob")
  if (!all(dim(sdd) == c(n_nodes, n_nodes))) stop("SDDprob must resolve to nodes x nodes")
  ldd <- .ivanal_resolve_time(LDDprob, t, T, "LDDprob", allow_null = TRUE)
  if (!is.null(ldd) && !all(dim(ldd) == c(n_nodes, n_nodes))) stop("LDDprob must resolve to nodes x nodes")
  rr <- if (length(LDDrate) == T) LDDrate[t] else LDDrate
  birth_move <- .ivanal_mix_movement(sdd, ldd, rr, n_nodes)

  env <- if (is.matrix(EnvEstabProb)) {
    if (nrow(EnvEstabProb) != n_nodes || ncol(EnvEstabProb) != T) stop("EnvEstabProb matrix must be nodes x Ntimesteps")
    EnvEstabProb[, t]
  } else rep_len(EnvEstabProb, n_nodes)
  pest <- if (is.matrix(PropaguleEstablishment)) {
    if (nrow(PropaguleEstablishment) != n_nodes || ncol(PropaguleEstablishment) != T) stop("PropaguleEstablishment matrix must be nodes x Ntimesteps")
    PropaguleEstablishment[, t]
  } else rep_len(PropaguleEstablishment, n_nodes)
  estab <- .ivanal_clip01(env * pest)

  mort <- .ivanal_expand_node_stage(MortalityProb, n_nodes, S, t, T, "MortalityProb", 0)
  manage <- .ivanal_expand_node_stage(ManagementExposure, n_nodes, S, t, T, "ManagementExposure", 0)
  fec_red <- .ivanal_expand_node_stage(FecundityReduction, n_nodes, S, t, T, "FecundityReduction", 0)
  source_survival <- 1 - .ivanal_clip01(mort) * .ivanal_clip01(manage)
  fec_mult <- 1 - .ivanal_clip01(fec_red) * .ivanal_clip01(manage)

  M <- matrix(0, K, K)            # expected inside descendants next timestep
  trans_prob <- vector("list", K) # categorical existing-animal outcomes
  birth_mu <- vector("list", K)   # Poisson inside births by destination type
  birth_export_mu <- numeric(K)    # successful outside births
  trans_export_prob <- numeric(K)  # existing-animal successful outside move

  for (node in seq_len(n_nodes)) {
    A <- .ivanal_transition_for_node(Transition, node, t, n_nodes, S, T)
    if (any(!is.finite(A)) || any(A < 0)) stop("Transition entries must be finite and non-negative")

    for (stage in seq_len(S)) {
      j <- idx(node, stage)
      surv_mgmt <- source_survival[node, stage]
      tp <- numeric(K)
      bm <- numeric(K)
      p_export_existing <- 0

      # Existing animal: terminal stage survival, or stasis/progression.
      if (stage == S) {
        p_stay <- .ivanal_clip01(A[S, S]) * surv_mgmt
        tp[idx(node, S)] <- tp[idx(node, S)] + p_stay
      } else {
        p_stasis <- .ivanal_clip01(A[stage, stage]) * surv_mgmt
        p_prog <- .ivanal_clip01(A[stage + 1L, stage]) * surv_mgmt
        if (p_stasis + p_prog > surv_mgmt + 1e-10) {
          scale <- surv_mgmt / (p_stasis + p_prog)
          p_stasis <- p_stasis * scale
          p_prog <- p_prog * scale
        }
        tp[idx(node, stage)] <- tp[idx(node, stage)] + p_stasis

        if (p_prog > 0) {
          Ps <- .ivanal_movement_for_stage(TransitionSDDprob, stage, t, T, n_nodes, S, "TransitionSDDprob")
          Pl <- .ivanal_movement_for_stage(TransitionLDDprob, stage, t, T, n_nodes, S, "TransitionLDDprob")
          trr <- if (length(TransitionLDDrate) == S - 1L) TransitionLDDrate[stage] else TransitionLDDrate
          Pm <- .ivanal_mix_movement(Ps, Pl, trr, n_nodes)
          if (is.null(Pm)) {
            tp[idx(node, stage + 1L)] <- tp[idx(node, stage + 1L)] + p_prog
          } else {
            pin <- pmax(0, Pm[node, ])
            if (sum(pin) > 1 + 1e-10) stop("Transition movement row sums may not exceed 1")
            for (dest in seq_len(n_nodes))
              tp[idx(dest, stage + 1L)] <- tp[idx(dest, stage + 1L)] + p_prog * pin[dest]
            p_export_existing <- p_prog * pmax(0, 1 - sum(pin)) * .ivanal_clip01(OutsideEstablishmentProb)
          }
        }
      }

      # Births: simulator convention A[1, reproductive stage]. Births are
      # Poisson and routed through SDD/LDD before establishment.
      fec <- if (stage >= 2L) pmax(0, A[1, stage]) else 0
      fec <- fec * fec_mult[node, stage]
      if (fec > 0) {
        pin_move <- pmax(0, birth_move[node, ])
        if (sum(pin_move) > 1 + 1e-10) stop("Birth movement row sums may not exceed 1")
        for (dest in seq_len(n_nodes)) {
          mu <- fec * pin_move[dest] * estab[dest]
          bm[idx(dest, 1L)] <- bm[idx(dest, 1L)] + mu
        }
        birth_export_mu[j] <- fec * pmax(0, 1 - sum(pin_move)) * .ivanal_clip01(OutsideEstablishmentProb)
      }

      trans_prob[[j]] <- tp
      birth_mu[[j]] <- bm
      trans_export_prob[j] <- p_export_existing
      M[, j] <- tp + bm
    }
  }

  list(MeanOperator = M,
       TransitionProb = trans_prob,
       BirthMean = birth_mu,
       BirthExportMean = birth_export_mu,
       TransitionExportProb = trans_export_prob)
}

.ivanal_noescape_backward <- function(steps) {
  T <- length(steps)
  K <- nrow(steps[[1]]$MeanOperator)
  q <- rep(1, K) # no escape after horizon = 1
  q_by_time <- vector("list", T + 1L)
  q_by_time[[T + 1L]] <- q

  for (tt in T:1L) {
    st <- steps[[tt]]
    q0 <- numeric(K)
    for (j in seq_len(K)) {
      tp <- st$TransitionProb[[j]]
      # Existing individual has one categorical next-state outcome or dies.
      # Export is a failure of the no-escape event.
      p_inside <- sum(tp)
      p_export <- st$TransitionExportProb[j]
      p_death_or_other <- pmax(0, 1 - p_inside - p_export)
      existing_factor <- p_death_or_other + sum(tp * q)

      # Poisson births split independently among inside types and successful
      # outside export. PGF evaluated at q inside and zero outside.
      bm <- st$BirthMean[[j]]
      birth_factor <- exp(sum(bm * (q - 1)) - st$BirthExportMean[j])
      q0[j] <- .ivanal_clip01(existing_factor * birth_factor)
    }
    q <- q0
    q_by_time[[tt]] <- q
  }
  list(q0 = q, q_by_time = q_by_time)
}

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
    Vertebrate = NULL) {

  if (!is.matrix(InitialPopulation)) stop("InitialPopulation must be nodes x stages")
  n_nodes <- nrow(InitialPopulation)
  if (ncol(InitialPopulation) != Nstages) stop("InitialPopulation must have Nstages columns")
  if (any(InitialPopulation < 0) || any(!is.finite(InitialPopulation))) stop("InitialPopulation must be finite and non-negative")
  if (Ntimesteps < 1L) stop("Ntimesteps must be >= 1")

  diagnostics <- character()
  if (!is.null(Vertebrate)) {
    special <- intersect(names(Vertebrate), c("Birth", "HomeRange", "Control", "Interaction"))
    if (length(special)) {
      diagnostics <- c(diagnostics,
        paste0("Vertebrate specialist modules present: ", paste(special, collapse = ", "),
               ". Prototype analytical solution uses native transition/movement parameters only unless equivalent effects are supplied explicitly through Transition, movement, ManagementExposure or fecundity inputs."))
    }
  }

  steps <- lapply(seq_len(Ntimesteps), function(t)
    .ivanal_animal_node_step(
      t = t, T = Ntimesteps, Nstages = Nstages, Transition = Transition,
      n_nodes = n_nodes, SDDprob = SDDprob, LDDprob = LDDprob,
      LDDrate = LDDrate, TransitionSDDprob = TransitionSDDprob,
      TransitionLDDprob = TransitionLDDprob,
      TransitionLDDrate = TransitionLDDrate,
      PropaguleEstablishment = PropaguleEstablishment,
      EnvEstabProb = EnvEstabProb, MortalityProb = MortalityProb,
      ManagementExposure = ManagementExposure,
      FecundityReduction = FecundityReduction,
      OutsideEstablishmentProb = OutsideEstablishmentProb))

  static <- all(vapply(steps[-1L], function(z)
    isTRUE(all.equal(z$MeanOperator, steps[[1]]$MeanOperator, tolerance = 1e-12)), logical(1)))

  if (static) {
    growth <- .ivanal_spectral_radius(steps[[1]]$MeanOperator)
    growth_definition <- "dominant eigenvalue of static low-density node x stage mean operator"
  } else {
    cycle <- diag(n_nodes * Nstages)
    for (t in seq_len(Ntimesteps)) cycle <- steps[[t]]$MeanOperator %*% cycle
    growth <- .ivanal_spectral_radius(cycle)^(1 / Ntimesteps)
    growth_definition <- "per-timestep geometric multiplier from the dominant eigenvalue of the finite-horizon operator product"
  }

  x <- as.vector(t(InitialPopulation))
  trajectory <- matrix(NA_real_, nrow = Ntimesteps + 1L, ncol = length(x))
  trajectory[1, ] <- x
  for (t in seq_len(Ntimesteps)) {
    x <- as.numeric(steps[[t]]$MeanOperator %*% x)
    trajectory[t + 1L, ] <- x
  }

  noescape <- .ivanal_noescape_backward(steps)
  x0 <- as.vector(t(InitialPopulation))
  # Independent initial lineages under the branching approximation.
  log_noescape <- sum(x0 * log(pmax(noescape$q0, .Machine$double.xmin)))
  escape <- .ivanal_clip01(1 - exp(log_noescape))

  stage_end <- colSums(matrix(x, nrow = n_nodes, ncol = Nstages, byrow = TRUE))
  names(stage_end) <- paste0("Stage", seq_len(Nstages))

  list(
    Variant = "animal",
    GrowthRate = growth,
    EndPopulation = sum(x),
    EndPopulationByStage = stage_end,
    EscapeProbability = escape,
    ExpectedTrajectory = rowSums(trajectory),
    MeanOperators = lapply(steps, `[[`, "MeanOperator"),
    Diagnostics = c(
      paste0("GrowthRate: ", growth_definition, "."),
      "EndPopulation: finite-horizon low-density expectation; carrying-capacity and density-dependent interactions are omitted.",
      "EscapeProbability: multitype branching no-escape recursion; births are Poisson and existing-animal stage/movement outcomes are categorical.",
      "Escape requires omitted movement row mass to represent movement outside the analysis domain; OutsideEstablishmentProb converts outside movement to successful escape.",
      diagnostics)
  )
}

# -----------------------------------------------------------------------------
# Generic transmissible-state analytical variant
# -----------------------------------------------------------------------------
# This branch is intentionally generic ecological state-process mathematics.
# The active lineage consists of a latent state E (optional) and a transmissible
# state I. S supplies hosts/resources for new E/I production; R and D are sinks.
# The operator is the exact first-order (rare-state) linearisation of the
# stochastic binomial infection rule used by INApestVertebrateDiseaseExtended.
# It preserves E->I progression, competing I stay/recover/die outcomes,
# state-specific movement, event order, time variation, and optional external
# introductions. It does not linearise host demographic feedbacks.

.ivanal_state_param <- function(x, t, T, n, name, nonnegative = FALSE) {
  if (is.function(x)) x <- x(timestep = t, n_nodes = n)
  if (is.matrix(x)) {
    if (!all(dim(x) == c(n, T))) stop(name, " matrix must be nodes x Ntimesteps")
    x <- x[, t]
  }
  x <- as.numeric(x)
  if (length(x) == 1L) x <- rep(x, n)
  else if (length(x) == T && length(x) != n) x <- rep(x[t], n)
  else if (length(x) != n) stop(name, " must resolve to scalar, node vector, timestep vector, nodes x time matrix, or function")
  if (any(!is.finite(x))) stop(name, " must be finite")
  if (nonnegative) pmax(0, x) else .ivanal_clip01(x)
}

.ivanal_state_move_matrix <- function(P, move_prob, t, T, n, name) {
  if (is.null(P)) return(list(Inside = diag(n), Outside = rep(0, n)))
  if (is.function(P)) P <- P(timestep = t)
  d <- dim(P)
  if (length(d) == 3L) {
    if (!all(d == c(n, n, T))) stop(name, " array must be nodes x nodes x Ntimesteps")
    P <- P[, , t, drop = FALSE][, , 1]
  }
  if (!is.matrix(P) || !all(dim(P) == c(n, n))) stop(name, " must resolve to nodes x nodes")
  if (any(!is.finite(P)) || any(P < 0)) stop(name, " entries must be finite and non-negative")
  rs <- rowSums(P)
  if (any(rs > 1 + 1e-10)) stop(name, " row sums may not exceed 1")
  mp <- .ivanal_state_param(move_prob, t, T, n, paste0(name, " move probability"))
  Q <- matrix(0, n, n)
  outside <- numeric(n)
  for (i in seq_len(n)) {
    Q[i, ] <- mp[i] * P[i, ]
    Q[i, i] <- Q[i, i] + (1 - mp[i])
    outside[i] <- mp[i] * pmax(0, 1 - rs[i])
  }
  list(Inside = Q, Outside = outside)
}

.ivanal_state_infection_derivative <- function(beta, susceptible, live, transmission, density_scale) {
  if (transmission == "frequency") {
    # d/dI S{1-exp[-beta I/live]} at I=0
    beta * ifelse(live > 0, susceptible / live, 0)
  } else {
    beta * susceptible / density_scale
  }
}

.ivanal_state_step <- function(t, T, HostPopulation, SusceptiblePopulation,
                               Beta, ProgressionProb, RecoveryProb,
                               StateMortalityProb, Transmission, DensityScale,
                               HasLatentState,
                               MovementMatrixE, MovementMatrixI,
                               MoveProbE, MoveProbI,
                               EventOrder,
                               OutsideEstablishmentProb) {
  n <- length(if (is.matrix(HostPopulation)) HostPopulation[,1] else HostPopulation)
  live <- .ivanal_state_param(HostPopulation, t, T, n, "HostPopulation", nonnegative = TRUE)
  susc <- .ivanal_state_param(SusceptiblePopulation, t, T, n, "SusceptiblePopulation", nonnegative = TRUE)
  if (any(susc > live + 1e-10)) stop("SusceptiblePopulation may not exceed HostPopulation in the rare-state baseline")
  beta <- .ivanal_state_param(Beta, t, T, n, "Beta", nonnegative = TRUE)
  prog <- .ivanal_state_param(ProgressionProb, t, T, n, "ProgressionProb")
  rec <- .ivanal_state_param(RecoveryProb, t, T, n, "RecoveryProb")
  mort <- .ivanal_state_param(StateMortalityProb, t, T, n, "StateMortalityProb")
  if (any(rec + mort > 1 + 1e-10)) stop("RecoveryProb + StateMortalityProb may not exceed 1")
  ds <- .ivanal_state_param(DensityScale, t, T, n, "DensityScale", nonnegative = TRUE)
  if (Transmission == "density" && any(ds <= 0)) stop("DensityScale must be > 0 for density-dependent transmission")
  infect <- .ivanal_state_infection_derivative(beta, susc, live, Transmission, ds)
  stayI <- pmax(0, 1 - rec - mort)

  ME <- .ivanal_state_move_matrix(MovementMatrixE, MoveProbE, t, T, n, "MovementMatrixE")
  MI <- .ivanal_state_move_matrix(MovementMatrixI, MoveProbI, t, T, n, "MovementMatrixI")
  oe <- .ivanal_state_param(OutsideEstablishmentProb, t, T, n, "OutsideEstablishmentProb")

  if (!HasLatentState) {
    # One active type per node: I. New infections enter I immediately.
    if (EventOrder == "move_transmit_progress") {
      A <- diag(stayI + infect, n) %*% t(MI$Inside)
      # An I that leaves before local events is an escape; no local descendants
      # are produced from that outside individual inside the analysis domain.
      p_export <- MI$Outside * oe
    } else {
      A <- t(MI$Inside) %*% diag(stayI + infect, n)
      # Existing/new I are moved after local events. The exact branching PGF
      # below treats their movement individually; mean outside descendants are
      # shown here only as a diagnostic.
      p_export <- MI$Outside * oe
    }
    return(list(MeanOperator=A, InfectDerivative=infect, Progression=prog,
                StayI=stayI, MoveE=ME, MoveI=MI, OutsideEstab=oe,
                HasLatentState=FALSE, EventOrder=EventOrder,
                MeanOutsideDiagnostic=(stayI + infect) * p_export))
  }

  # Type order: E nodes 1..n, then I nodes 1..n.
  if (EventOrder == "move_transmit_progress") {
    A_EE <- diag(1 - prog, n) %*% t(ME$Inside)
    A_EI <- diag(infect, n) %*% t(MI$Inside)
    A_IE <- diag(prog, n) %*% t(ME$Inside)
    A_II <- diag(stayI, n) %*% t(MI$Inside)
  } else {
    A_EE <- t(ME$Inside) %*% diag(1 - prog, n)
    A_EI <- t(ME$Inside) %*% diag(infect, n)
    A_IE <- t(MI$Inside) %*% diag(prog, n)
    A_II <- t(MI$Inside) %*% diag(stayI, n)
  }
  A <- rbind(cbind(A_EE, A_EI), cbind(A_IE, A_II))
  list(MeanOperator=A, InfectDerivative=infect, Progression=prog,
       StayI=stayI, MoveE=ME, MoveI=MI, OutsideEstab=oe,
       HasLatentState=TRUE, EventOrder=EventOrder)
}

# Exact rare-lineage no-escape PGF for the state process. Existing E/I animals
# have categorical state/movement outcomes; newly affected hosts are binomial in
# the stochastic model but converge to Poisson offspring under the first-order
# rare-state branching limit. This is the same approximation level as the mean
# operator rather than an additional mean-matched shortcut.
.ivanal_state_noescape_backward <- function(steps) {
  T <- length(steps); n <- length(steps[[1]]$StayI)
  latent <- isTRUE(steps[[1]]$HasLatentState)
  K <- if (latent) 2L*n else n
  q <- rep(1, K)
  q_by_time <- vector("list", T + 1L); q_by_time[[T+1L]] <- q

  for (tt in T:1L) {
    st <- steps[[tt]]; q0 <- numeric(K)
    if (!latent) {
      for (i in seq_len(n)) {
        s <- st$StayI[i]; mu <- st$InfectDerivative[i]
        QI <- st$MoveI$Inside[i, ]; outI <- st$MoveI$Outside[i] * st$OutsideEstab[i]
        if (st$EventOrder == "move_transmit_progress") {
          # Movement first: outside movement itself escapes. Conditional on an
          # internal destination j, persistence is Bernoulli and new I are
          # Poisson at j.
          no_move_escape <- pmax(0, 1 - sum(QI) - st$MoveI$Outside[i])
          val <- no_move_escape + st$MoveI$Outside[i] * (1 - st$OutsideEstab[i])
          for (j in seq_len(n)) {
            existing <- (1 - s) + s*q[j]
            births <- exp(mu * (q[j] - 1))
            val <- val + QI[j] * existing * births
          }
          q0[i] <- .ivanal_clip01(val)
        } else {
          # Local first: persistent parent (Bernoulli) and Poisson new I then
          # move independently with the I movement kernel.
          zmove <- sum(QI * q) + st$MoveI$Outside[i] * (1 - st$OutsideEstab[i]) + pmax(0, 1 - sum(QI) - st$MoveI$Outside[i])
          existing <- (1 - s) + s*zmove
          births <- exp(mu * (zmove - 1))
          q0[i] <- .ivanal_clip01(existing * births)
        }
      }
    } else {
      qE <- q[seq_len(n)]; qI <- q[n + seq_len(n)]
      for (i in seq_len(n)) {
        # One E parent: no new infections directly; categorical stay-E/progress-I.
        QE <- st$MoveE$Inside[i, ]; outE <- st$MoveE$Outside[i]
        if (st$EventOrder == "move_transmit_progress") {
          valE <- outE * (1 - st$OutsideEstab[i]) + pmax(0, 1 - sum(QE) - outE)
          for (j in seq_len(n))
            valE <- valE + QE[j] * ((1-st$Progression[j])*qE[j] + st$Progression[j]*qI[j])
        } else {
          zE <- sum(QE*qE) + outE*(1-st$OutsideEstab[i]) + pmax(0,1-sum(QE)-outE)
          zI <- sum(st$MoveI$Inside[i,]*qI) + st$MoveI$Outside[i]*(1-st$OutsideEstab[i]) + pmax(0,1-sum(st$MoveI$Inside[i,])-st$MoveI$Outside[i])
          valE <- (1-st$Progression[i])*zE + st$Progression[i]*zI
        }
        q0[i] <- .ivanal_clip01(valE)

        # One I parent: persistence is Bernoulli; new E are Poisson.
        QI <- st$MoveI$Inside[i, ]; outI <- st$MoveI$Outside[i]
        if (st$EventOrder == "move_transmit_progress") {
          valI <- outI*(1-st$OutsideEstab[i]) + pmax(0,1-sum(QI)-outI)
          for (j in seq_len(n)) {
            existing <- (1-st$StayI[j]) + st$StayI[j]*qI[j]
            births <- exp(st$InfectDerivative[j]*(qE[j]-1))
            valI <- valI + QI[j]*existing*births
          }
        } else {
          zI <- sum(QI*qI) + outI*(1-st$OutsideEstab[i]) + pmax(0,1-sum(QI)-outI)
          QEnew <- st$MoveE$Inside[i, ]; outEnew <- st$MoveE$Outside[i]
          zE <- sum(QEnew*qE) + outEnew*(1-st$OutsideEstab[i]) + pmax(0,1-sum(QEnew)-outEnew)
          existing <- (1-st$StayI[i]) + st$StayI[i]*zI
          births <- exp(st$InfectDerivative[i]*(zE-1))
          valI <- existing*births
        }
        q0[n+i] <- .ivanal_clip01(valI)
      }
    }
    q <- q0; q_by_time[[tt]] <- q
  }
  list(q0=q, q_by_time=q_by_time)
}

.ivanal_external_intro <- function(ExternalIntroductionProb, ExternalIntroductionCount,
                                   ExternalIntroductionState, t, T, n, latent) {
  if (is.null(ExternalIntroductionProb)) return(NULL)
  p <- .ivanal_state_param(ExternalIntroductionProb, t, T, n, "ExternalIntroductionProb")
  count <- .ivanal_state_param(ExternalIntroductionCount, t, T, n, "ExternalIntroductionCount", nonnegative=TRUE)
  if (any(abs(count - round(count)) > 1e-10)) stop("ExternalIntroductionCount must resolve to non-negative integer counts")
  count <- round(count)
  state <- toupper(ExternalIntroductionState)
  if (!(state %in% c("E","I"))) stop("ExternalIntroductionState must be E or I")
  if (state == "E" && !latent) state <- "I"
  list(prob=p, count=count, state=state)
}

INApestVertebrateAnalyticalPathogen <- function(
    Ntimesteps,
    HostPopulation,
    InitialAffected = NULL,
    InitialState = NULL,
    SusceptiblePopulation = HostPopulation,
    Beta,
    ProgressionProb = 1,
    RecoveryProb = 0,
    StateMortalityProb = 0,
    Transmission = c("frequency", "density"),
    DensityScale = 1,
    HasLatentState = TRUE,
    MovementMatrix = NULL,
    MovementMatrixE = MovementMatrix,
    MovementMatrixI = MovementMatrix,
    MoveProbAffected = 0,
    MoveProbE = MoveProbAffected,
    MoveProbI = MoveProbAffected,
    EventOrder = c("move_transmit_progress", "transmit_progress_move"),
    OutsideEstablishmentProb = 1,
    ExternalIntroductionProb = NULL,
    ExternalIntroductionCount = 1,
    ExternalIntroductionState = if (HasLatentState) "E" else "I") {

  Transmission <- match.arg(Transmission); EventOrder <- match.arg(EventOrder)
  if (Ntimesteps < 1L) stop("Ntimesteps must be >= 1")
  H0 <- if (is.matrix(HostPopulation)) HostPopulation[,1] else as.numeric(HostPopulation)
  n <- length(H0)
  if (!n || any(!is.finite(H0)) || any(H0 < 0)) stop("HostPopulation must define non-negative host abundance by node")

  latent <- isTRUE(HasLatentState)
  if (!is.null(InitialState)) {
    X <- as.matrix(InitialState)
    need <- if (latent) 2L else 1L
    if (!all(dim(X) == c(n, need))) stop("InitialState must be nodes x active states (E,I when latent; I otherwise)")
    if (latent) { E0 <- X[,1]; I0 <- X[,2] } else { E0 <- numeric(n); I0 <- X[,1] }
  } else {
    if (is.null(InitialAffected)) stop("Supply InitialAffected or InitialState")
    z <- as.numeric(InitialAffected); if (length(z)!=n) stop("InitialAffected must match HostPopulation")
    E0 <- if (latent) z else numeric(n); I0 <- if (latent) numeric(n) else z
  }
  if (any(c(E0,I0) < 0) || any(!is.finite(c(E0,I0)))) stop("Initial active-state abundance must be finite and non-negative")

  steps <- lapply(seq_len(Ntimesteps), function(t) .ivanal_state_step(
    t=t, T=Ntimesteps, HostPopulation=HostPopulation, SusceptiblePopulation=SusceptiblePopulation,
    Beta=Beta, ProgressionProb=ProgressionProb, RecoveryProb=RecoveryProb,
    StateMortalityProb=StateMortalityProb, Transmission=Transmission, DensityScale=DensityScale,
    HasLatentState=latent, MovementMatrixE=MovementMatrixE, MovementMatrixI=MovementMatrixI,
    MoveProbE=MoveProbE, MoveProbI=MoveProbI, EventOrder=EventOrder,
    OutsideEstablishmentProb=OutsideEstablishmentProb))

  static <- if (length(steps)==1L) TRUE else all(vapply(steps[-1L], function(z)
    isTRUE(all.equal(z$MeanOperator, steps[[1]]$MeanOperator, tolerance=1e-12)), logical(1)))
  if (static) {
    growth <- .ivanal_spectral_radius(steps[[1]]$MeanOperator)
    growth_def <- "dominant eigenvalue of the static rare-state operator"
  } else {
    cycle <- diag(nrow(steps[[1]]$MeanOperator))
    for (t in seq_len(Ntimesteps)) cycle <- steps[[t]]$MeanOperator %*% cycle
    growth <- .ivanal_spectral_radius(cycle)^(1/Ntimesteps)
    growth_def <- "per-timestep geometric multiplier of the finite-horizon rare-state operator product"
  }

  x <- if (latent) c(E0,I0) else I0
  traj <- matrix(NA_real_, Ntimesteps+1L, length(x)); traj[1,] <- x
  intro_records <- vector("list", Ntimesteps)
  for (t in seq_len(Ntimesteps)) {
    intro <- .ivanal_external_intro(ExternalIntroductionProb, ExternalIntroductionCount,
                                    ExternalIntroductionState, t, Ntimesteps, n, latent)
    # Parent simulator adds external arrivals before the Interaction in each step.
    if (!is.null(intro)) {
      add <- intro$prob * intro$count
      if (latent && intro$state=="E") x[seq_len(n)] <- x[seq_len(n)] + add
      else if (latent) x[n+seq_len(n)] <- x[n+seq_len(n)] + add
      else x <- x + add
      intro_records[[t]] <- intro
    }
    x <- as.numeric(steps[[t]]$MeanOperator %*% x)
    traj[t+1L,] <- x
  }

  # Baseline no-escape probabilities for a lineage present immediately before
  # each interaction timestep. Compose independent initial lineages and the
  # Bernoulli external-introduction process exactly at the introduction-event level.
  neb <- .ivanal_state_noescape_backward(steps)
  x0 <- if (latent) c(E0,I0) else I0
  log_no <- sum(x0 * log(pmax(neb$q0, .Machine$double.xmin)))
  if (!is.null(ExternalIntroductionProb)) {
    for (t in seq_len(Ntimesteps)) {
      intro <- .ivanal_external_intro(ExternalIntroductionProb, ExternalIntroductionCount,
                                      ExternalIntroductionState, t, Ntimesteps, n, latent)
      q_t <- neb$q_by_time[[t]]
      ids <- if (latent && intro$state=="E") seq_len(n) else if (latent) n+seq_len(n) else seq_len(n)
      # One Bernoulli event per node introduces 'count' independent lineages.
      fac <- (1-intro$prob) + intro$prob * (pmax(q_t[ids],0)^intro$count)
      log_no <- log_no + sum(log(pmax(fac, .Machine$double.xmin)))
    }
  }
  escape <- .ivanal_clip01(1-exp(log_no))

  if (latent) {
    endE <- sum(x[seq_len(n)]); endI <- sum(x[n+seq_len(n)])
    bystate <- c(E=endE, I=endI)
    active_traj <- rowSums(traj)
  } else {
    bystate <- c(I=sum(x)); active_traj <- rowSums(traj)
  }

  list(
    Variant="pathogen",
    GrowthRate=growth,
    EndPopulation=sum(x),
    EndPopulationByState=bystate,
    EscapeProbability=escape,
    ExpectedTrajectory=active_traj,
    ExpectedStateTrajectory=traj,
    MeanOperators=lapply(steps, `[[`, "MeanOperator"),
    Diagnostics=c(
      paste0("GrowthRate: ", growth_def, "; external introductions are excluded from this intrinsic multiplier."),
      paste0("EndPopulation: expected active-state abundance at timestep ", Ntimesteps,
             if (latent) " (E + I)." else " (I)."),
      "The operator is the first-order rare-state linearisation of the stochastic transmission rule around HostPopulation/SusceptiblePopulation.",
      paste0("Event order reproduced as ", EventOrder, "; E and I can have different movement probabilities/matrices."),
      "EscapeProbability: finite-horizon multitype branching no-escape recursion over active-state lineages; external Bernoulli introductions are composed separately when supplied.",
      "Recovered/dead states are sinks and therefore are not included as branching types or in EndPopulation.",
      "Host demographic feedback, susceptible depletion at finite prevalence, contact saturation, and arbitrary Interaction/HomeRange hooks remain stochastic/nonlinear effects rather than being hidden in this low-density operator."
    )
  )
}

INApestAnalyticalVertebrate <- function(Variant = c("animal", "pathogen"), ...) {
  Variant <- match.arg(Variant)
  if (Variant == "animal") return(INApestVertebrateAnalyticalAnimal(...))
  INApestVertebrateAnalyticalPathogen(...)
}

# Convenience wrapper for INApestVertebratePoint.
#
# The continuous-space problem is deliberately separated into two layers:
#   (a) contract continuous kernels onto an analysis grid using the existing
#       INApest point analytical machinery; then
#   (b) solve the resulting cell x stage animal branching process here.
#
# This wrapper accepts those contracted source->destination matrices directly.
# With a homogeneous landscape and no explicit boundary, a one-cell contraction
# is exact for total growth/end population because location does not affect the
# count process. Escape then requires an explicit contraction whose missing row
# mass represents movement outside the analysis domain.
INApestVertebrateAnalyticalPoint <- function(
    Ntimesteps,
    Nstages,
    Transition,
    InitialStageCounts,
    ReproductiveMovement,
    TransitionMovement = NULL,
    PropaguleEstablishment = 1,
    EnvEstabProb = 1,
    MortalityProb = 0,
    ManagementExposure = 0,
    FecundityReduction = 0,
    OutsideEstablishmentProb = 1,
    Vertebrate = NULL) {

  if (!is.matrix(ReproductiveMovement) || nrow(ReproductiveMovement) != ncol(ReproductiveMovement))
    stop("ReproductiveMovement must be a contracted analysis-cell x analysis-cell matrix")
  n_cells <- nrow(ReproductiveMovement)

  if (is.vector(InitialStageCounts) && !is.list(InitialStageCounts)) {
    if (length(InitialStageCounts) != Nstages)
      stop("InitialStageCounts vector must have Nstages entries")
    init <- matrix(0, nrow = n_cells, ncol = Nstages)
    init[1, ] <- InitialStageCounts
  } else {
    init <- as.matrix(InitialStageCounts)
    if (!all(dim(init) == c(n_cells, Nstages)))
      stop("InitialStageCounts matrix must be analysis cells x Nstages")
  }

  transition_sdd <- NULL
  if (!is.null(TransitionMovement)) {
    transition_sdd <- TransitionMovement
  }

  ans <- INApestVertebrateAnalyticalAnimal(
    Ntimesteps = Ntimesteps,
    Nstages = Nstages,
    Transition = Transition,
    InitialPopulation = init,
    SDDprob = ReproductiveMovement,
    LDDprob = NULL,
    LDDrate = 0,
    TransitionSDDprob = transition_sdd,
    TransitionLDDprob = NULL,
    TransitionLDDrate = 0,
    PropaguleEstablishment = PropaguleEstablishment,
    EnvEstabProb = EnvEstabProb,
    MortalityProb = MortalityProb,
    ManagementExposure = ManagementExposure,
    FecundityReduction = FecundityReduction,
    OutsideEstablishmentProb = OutsideEstablishmentProb,
    Vertebrate = Vertebrate)

  ans$Variant <- "animal_point"
  ans$Diagnostics <- c(
    "Point model: continuous-space kernels must first be contracted onto a spatial analysis grid; the resulting cell-to-cell probabilities are the ReproductiveMovement/TransitionMovement inputs.",
    ans$Diagnostics)
  ans
}

###############################################################################
### Unified vertebrate-aware user-facing dispatcher
###############################################################################

# Preserve the complete pre-vertebrate analytical dispatcher, including all
# point, fecundity, dynamic-information and strengthened finite-horizon PGF
# refinements. Existing model calls are delegated to this object unchanged.
INApestAnalytical_prevertebrate <- INApestAnalytical

INApestAnalytical <- function(...) {
  args <- list(...)
  Model <- if (!is.null(args$Model)) as.character(args$Model)[1L] else "INApest"

  legacy_models <- c("INApest", "INApestMeta", "INApestMetaTransitionMatrix",
                     "INApestMetaMultipleLandUse", "INApestMetaPoint",
                     "INApestPointTransitionMatrix")
  vertebrate_models <- c("INApestVertebrateNode", "INApestVertebratePoint",
                         "INApestVertebrateStateProcess")
  all_models <- c(legacy_models, vertebrate_models)
  if (!Model %in% all_models)
    stop("Model must be one of: ", paste(all_models, collapse = ", "))

  if (Model %in% legacy_models)
    return(do.call(INApestAnalytical_prevertebrate, args))

  getarg <- function(name, default = NULL, required = FALSE) {
    if (name %in% names(args)) return(args[[name]])
    if (required) stop(Model, " requires ", name)
    default
  }

  if (Model == "INApestVertebrateNode") {
    InitialState <- getarg("InitialState", required = TRUE)
    Transition <- getarg("Transition", required = TRUE)
    SDDprob <- getarg("SDDprob", required = TRUE)
    Nstages <- getarg("Nstages", NULL)
    if (is.null(Nstages)) {
      if (!is.matrix(InitialState))
        stop("INApestVertebrateNode requires Nstages when InitialState is not a matrix")
      Nstages <- ncol(InitialState)
    }
    LDDprob <- getarg("LDDprob", 0)
    if (is.null(LDDprob) || (length(LDDprob) == 1L && isTRUE(LDDprob == 0))) LDDprob <- NULL
    ManagementExposure <- if ("ManagementExposure" %in% names(args))
      args$ManagementExposure else getarg("ManageProb", 0)

    ans <- INApestVertebrateAnalyticalAnimal(
      Ntimesteps = getarg("Ntimesteps", 10),
      Nstages = Nstages,
      Transition = Transition,
      InitialPopulation = InitialState,
      SDDprob = SDDprob,
      LDDprob = LDDprob,
      LDDrate = getarg("LDDrate", 0),
      TransitionSDDprob = getarg("TransitionSDDprob", NULL),
      TransitionLDDprob = getarg("TransitionLDDprob", NULL),
      TransitionLDDrate = getarg("TransitionLDDrate", 0),
      PropaguleEstablishment = getarg("PropaguleEstablishment", 1),
      EnvEstabProb = getarg("EnvEstabProb", 1),
      MortalityProb = getarg("MortalityProb", 0),
      ManagementExposure = ManagementExposure,
      FecundityReduction = getarg("FecundityReduction", 0),
      OutsideEstablishmentProb = getarg("OutsideEstablishmentProb", 1),
      Vertebrate = getarg("Vertebrate", NULL))

    ans$Model <- Model
    ans$HeadlineEstimands <- list(
      GrowthRate = ans$GrowthRate,
      EndPopulation = ans$EndPopulation,
      EscapeProbability = ans$EscapeProbability)
    return(ans)
  }

  if (Model == "INApestVertebratePoint") {
    InitialStageCounts <- if ("InitialStageCounts" %in% names(args))
      args$InitialStageCounts else getarg("InitialState", required = TRUE)
    ReproductiveMovement <- if ("ReproductiveMovement" %in% names(args))
      args$ReproductiveMovement else getarg("SDDprob", required = TRUE)
    Transition <- getarg("Transition", required = TRUE)
    Nstages <- getarg("Nstages", NULL)
    if (is.null(Nstages))
      Nstages <- if (is.matrix(InitialStageCounts)) ncol(InitialStageCounts) else length(InitialStageCounts)
    ManagementExposure <- if ("ManagementExposure" %in% names(args))
      args$ManagementExposure else getarg("ManageProb", 0)

    ans <- INApestVertebrateAnalyticalPoint(
      Ntimesteps = getarg("Ntimesteps", 10),
      Nstages = Nstages,
      Transition = Transition,
      InitialStageCounts = InitialStageCounts,
      ReproductiveMovement = ReproductiveMovement,
      TransitionMovement = getarg("TransitionMovement", NULL),
      PropaguleEstablishment = getarg("PropaguleEstablishment", 1),
      EnvEstabProb = getarg("EnvEstabProb", 1),
      MortalityProb = getarg("MortalityProb", 0),
      ManagementExposure = ManagementExposure,
      FecundityReduction = getarg("FecundityReduction", 0),
      OutsideEstablishmentProb = getarg("OutsideEstablishmentProb", 1),
      Vertebrate = getarg("Vertebrate", NULL))

    ans$Model <- Model
    ans$HeadlineEstimands <- list(
      GrowthRate = ans$GrowthRate,
      EndPopulation = ans$EndPopulation,
      EscapeProbability = ans$EscapeProbability)
    return(ans)
  }

  # Generic affected/transmissible state process. Keep this user-facing branch
  # ecological and model-generic; the internal solver name is retained for
  # backwards development provenance only.
  state_names <- c("Ntimesteps", "HostPopulation", "InitialAffected", "InitialState",
                   "SusceptiblePopulation", "Beta", "ProgressionProb", "RecoveryProb",
                   "StateMortalityProb", "Transmission", "DensityScale", "HasLatentState",
                   "MovementMatrix", "MovementMatrixE", "MovementMatrixI",
                   "MoveProbAffected", "MoveProbE", "MoveProbI", "EventOrder",
                   "OutsideEstablishmentProb", "ExternalIntroductionProb",
                   "ExternalIntroductionCount", "ExternalIntroductionState")
  state_args <- args[intersect(names(args), state_names)]
  state_args$Ntimesteps <- getarg("Ntimesteps", 10)
  state_args$HostPopulation <- getarg("HostPopulation", required = TRUE)
  state_args$Beta <- getarg("Beta", required = TRUE)
  state_args$OutsideEstablishmentProb <- getarg("OutsideEstablishmentProb", 1)
  state_args$Model <- NULL

  ans <- do.call(INApestVertebrateAnalyticalPathogen, state_args)
  ans$Variant <- "state_process"
  ans$Model <- Model
  ans$HeadlineEstimands <- list(
    GrowthRate = ans$GrowthRate,
    EndPopulation = ans$EndPopulation,
    EscapeProbability = ans$EscapeProbability)
  ans
}
###############################################################################
### Binary INApest pathogen analytical extension
###
### Purpose
###   * exact finite-state propagation for small Binary INApest systems;
###   * intrinsic rare-pathogen growth operator for scalable screening;
###   * finite-horizon branching approximations for extinction and containment
###     escape when exact state propagation is deliberately disabled/too large;
###   * one shared HaveInfo state with user-selectable local evidence source:
###       InformationAcquisition = "host", "pathogen", or "both".
###
### The exact engine mirrors INApest event ordering. It treats management,
### information persistence/programmed stopping, host dispersal, pathogen
### transmission, clearance, introduction, pathogen-associated host extinction,
### and host/pathogen detection as discrete Bernoulli events. Parameter SDs are
### intentionally outside this exact conditional calculation; stochastic INApest
### validation should use SD = 0 when comparing against it.
###############################################################################

.ina_bp_clip01 <- function(x) pmin(1, pmax(0, as.numeric(x)))

.ina_bp_resolve_node <- function(x, timestep, n, Ntimesteps, name,
                                 allow_na = FALSE, default = NULL) {
  if (is.null(x)) {
    if (is.null(default)) stop(name, " may not be NULL")
    x <- default
  }
  if (is.function(x)) {
    fm <- names(formals(x))
    a <- list(timestep = timestep, n_nodes = n, Ntimesteps = Ntimesteps)
    if (!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a), fm)]
    x <- do.call(x, a)
  }
  d <- dim(x)
  if (!is.null(d)) {
    if (length(d) == 2L && all(d == c(n, Ntimesteps))) x <- x[, timestep]
    else stop(name, " must be scalar, length nodes, nodes x Ntimesteps, or resolver function")
  } else {
    x <- as.numeric(x)
    if (length(x) == 1L) x <- rep(x, n)
    else if (length(x) == n) x <- x
    else if (length(x) == Ntimesteps && Ntimesteps != n) x <- rep(x[timestep], n)
    else stop(name, " has unsupported or ambiguous shape")
  }
  x <- as.numeric(x)
  if (length(x) == 1L) x <- rep(x, n)
  if (length(x) != n) stop(name, " resolved to wrong length")
  if (!allow_na && any(is.na(x))) stop(name, " may not contain NA")
  x
}

.ina_bp_resolve_connection <- function(x, timestep, n, Ntimesteps, name,
                                       scalar_zero_ok = TRUE) {
  if (is.function(x)) {
    fm <- names(formals(x))
    a <- list(timestep = timestep, n_nodes = n, Ntimesteps = Ntimesteps)
    if (!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a), fm)]
    x <- do.call(x, a)
  }
  if (length(x) == 1L && is.numeric(x)) {
    if (scalar_zero_ok && identical(as.numeric(x), 0)) return(matrix(0, n, n))
    return(matrix(as.numeric(x), n, n))
  }
  d <- dim(x)
  if (length(d) == 2L && all(d == c(n, n))) return(as.matrix(x))
  if (length(d) == 3L && all(d == c(n, n, Ntimesteps)))
    return(as.matrix(x[, , timestep, drop = TRUE]))
  stop(name, " must be scalar, nodes x nodes, nodes x nodes x Ntimesteps, or resolver function")
}

.ina_bp_contact <- function(Pathogen, timestep, n, Ntimesteps) {
  x <- Pathogen$ContactMatrix
  if (is.null(x)) return(diag(n))
  C <- .ina_bp_resolve_connection(x, timestep, n, Ntimesteps, "Pathogen ContactMatrix", FALSE)
  if (any(!is.finite(C)) || any(C < 0)) stop("Pathogen ContactMatrix entries must be finite and non-negative")
  C
}

.ina_bp_information_acquisition <- function(InformationAcquisition, Pathogen) {
  if (is.null(InformationAcquisition)) {
    return(if (!is.null(Pathogen) && isTRUE(Pathogen$DetectionTriggersInfo)) "both" else "host")
  }
  mode <- match.arg(as.character(InformationAcquisition)[1L], c("host", "pathogen", "both"))
  if (mode %in% c("pathogen", "both") && is.null(Pathogen))
    stop("InformationAcquisition = '", mode, "' requires a Binary Pathogen specification")
  mode
}

.ina_bp_state <- function(H, P, Info, Age, Escaped = 0L, M = NULL, D = NULL, Psource = NULL) {
  H <- as.integer(H != 0)
  P <- as.integer(P != 0) * H
  Info <- as.integer(Info != 0)
  Age <- as.integer(Age)
  Age[Info == 0L] <- -1L
  list(H = H, P = P, Info = Info, Age = Age,
       Escaped = as.integer(Escaped != 0), M = M, D = D, Psource = Psource)
}

.ina_bp_state_key <- function(s) {
  paste(c(s$H, -2L, s$P, -3L, s$Info, -4L, s$Age, -5L,
          s$Escaped,
          if (is.null(s$M)) -6L else c(-6L, s$M),
          if (is.null(s$D)) -7L else c(-7L, s$D),
          if (is.null(s$Psource)) -8L else c(-8L, s$Psource)), collapse = ",")
}

.ina_bp_dist_new <- function() new.env(hash = TRUE, parent = emptyenv())

.ina_bp_dist_add <- function(env, state, prob) {
  if (!is.finite(prob) || prob <= 0) return(invisible(NULL))
  state <- .ina_bp_state(state$H, state$P, state$Info, state$Age,
                         state$Escaped, state$M, state$D, state$Psource)
  key <- .ina_bp_state_key(state)
  if (exists(key, envir = env, inherits = FALSE)) {
    z <- get(key, envir = env, inherits = FALSE)
    z$prob <- z$prob + prob
    assign(key, z, envir = env)
  } else assign(key, list(state = state, prob = prob), envir = env)
  invisible(NULL)
}

.ina_bp_dist_items <- function(env) {
  nm <- ls(env, all.names = TRUE)
  if (!length(nm)) return(list())
  unname(mget(nm, envir = env, inherits = FALSE))
}

.ina_bp_dist_mass <- function(env) {
  it <- .ina_bp_dist_items(env)
  if (!length(it)) return(0)
  sum(vapply(it, function(z) z$prob, numeric(1)))
}

.ina_bp_set <- function(s, field, i, value) {
  x <- s[[field]]
  x[i] <- value
  s[[field]] <- x
  s
}

.ina_bp_branch_nodes <- function(dist, nodes, probability, update) {
  if (!length(nodes)) return(dist)
  cur <- dist
  for (i in nodes) {
    nxt <- .ina_bp_dist_new()
    for (z in .ina_bp_dist_items(cur)) {
      s <- z$state; w <- z$prob
      p <- .ina_bp_clip01(probability(s, i))[1L]
      if (p < 1) .ina_bp_dist_add(nxt, s, w * (1 - p))
      if (p > 0) .ina_bp_dist_add(nxt, update(s, i), w * p)
    }
    cur <- nxt
  }
  cur
}

.ina_bp_apply_deterministic <- function(dist, fun) {
  out <- .ina_bp_dist_new()
  for (z in .ina_bp_dist_items(dist)) .ina_bp_dist_add(out, fun(z$state), z$prob)
  out
}

.ina_bp_use_persistence <- function(InfoPersistenceSteps, n, Ntimesteps) {
  for (tt in seq_len(Ntimesteps)) {
    k <- .ina_bp_resolve_node(InfoPersistenceSteps, tt, n, Ntimesteps,
                              "InfoPersistenceSteps", allow_na = TRUE, default = NA)
    if (any(!is.na(k))) return(TRUE)
  }
  FALSE
}

.ina_bp_initial_distribution <- function(InitialState, InitialInfo,
                                         InitialPathogen, Pathogen,
                                         DetectionProb, ApplyInitialDetection,
                                         InformationAcquisition,
                                         UseInfoPersistence, n, Ntimesteps,
                                         OutsideNodes = integer(0)) {
  Hprob <- .ina_bp_clip01(.ina_bp_resolve_node(InitialState, 1L, n, Ntimesteps, "InitialState"))
  Iprob <- .ina_bp_clip01(.ina_bp_resolve_node(InitialInfo, 1L, n, Ntimesteps, "InitialInfo"))
  Pprob <- if (is.null(InitialPathogen)) {
    .ina_bp_clip01(.ina_bp_resolve_node(Pathogen$InitialPresent, 1L, n, Ntimesteps, "InitialPresent"))
  } else .ina_bp_clip01(.ina_bp_resolve_node(InitialPathogen, 1L, n, Ntimesteps, "InitialPathogen"))

  d <- .ina_bp_dist_new()
  .ina_bp_dist_add(d, .ina_bp_state(rep(0L, n), rep(0L, n), rep(0L, n), rep(-1L, n)), 1)

  # Initial host occupancy, explicit InitialInfo, then pathogen occupancy.
  d <- .ina_bp_branch_nodes(d, seq_len(n),
    function(s, i) Hprob[i],
    function(s, i) .ina_bp_set(s, "H", i, 1L))
  d <- .ina_bp_branch_nodes(d, seq_len(n),
    function(s, i) Iprob[i],
    function(s, i) {
      s <- .ina_bp_set(s, "Info", i, 1L)
      if (UseInfoPersistence && s$H[i] == 1L) s <- .ina_bp_set(s, "Age", i, 0L)
      s
    })
  d <- .ina_bp_branch_nodes(d, seq_len(n),
    function(s, i) if (s$H[i] == 1L) Pprob[i] else 0,
    function(s, i) .ina_bp_set(s, "P", i, 1L))

  host_acq <- InformationAcquisition %in% c("host", "both")
  path_acq <- InformationAcquisition %in% c("pathogen", "both")

  if (isTRUE(ApplyInitialDetection) && host_acq) {
    hd <- .ina_bp_clip01(.ina_bp_resolve_node(DetectionProb, 1L, n, Ntimesteps, "DetectionProb"))
    d <- .ina_bp_branch_nodes(d, seq_len(n),
      function(s, i) if (s$H[i] == 1L) hd[i] else 0,
      function(s, i) {
        s <- .ina_bp_set(s, "Info", i, 1L)
        if (UseInfoPersistence) s <- .ina_bp_set(s, "Age", i, 0L)
        s
      })
  }
  if (isTRUE(ApplyInitialDetection) && path_acq) {
    pd <- .ina_bp_clip01(.ina_bp_resolve_node(Pathogen$DetectionProb, 1L, n, Ntimesteps, "Pathogen DetectionProb"))
    d <- .ina_bp_branch_nodes(d, seq_len(n),
      function(s, i) if (s$P[i] == 1L) pd[i] else 0,
      function(s, i) {
        s <- .ina_bp_set(s, "Info", i, 1L)
        if (UseInfoPersistence) s <- .ina_bp_set(s, "Age", i, 0L)
        s
      })
  }

  if (length(OutsideNodes)) {
    d <- .ina_bp_apply_deterministic(d, function(s) {
      if (any(s$P[OutsideNodes] == 1L)) s$Escaped <- 1L
      s
    })
  }
  d
}

.ina_bp_step_parameters <- function(timestep, n, Ntimesteps,
                                    Pathogen,
                                    DetectionProb, ManageProb, EradicationProb,
                                    SpreadReduction, InfoRetentionProb,
                                    InfoPersistenceSteps, EnvEstabProb, Survival,
                                    SDDprob, LDDprob, SEAM,
                                    OngoingExternalInvasion, InvasionRisk,
                                    OngoingExternalInfo, ExternalInfoProb,
                                    ForceNoPathogenIntroduction = FALSE) {
  SDD <- .ina_bp_resolve_connection(SDDprob, timestep, n, Ntimesteps, "SDDprob")
  LDD <- .ina_bp_resolve_connection(LDDprob, timestep, n, Ntimesteps, "LDDprob")
  D <- 1 - (1 - SDD) * (1 - LDD)
  env <- .ina_bp_clip01(.ina_bp_resolve_node(EnvEstabProb, timestep, n, Ntimesteps, "EnvEstabProb"))
  host_disp <- sweep(D, 2L, env, `*`)
  host_disp <- matrix(.ina_bp_clip01(host_disp), n, n)

  C <- .ina_bp_contact(Pathogen, timestep, n, Ntimesteps)
  tr <- .ina_bp_clip01(.ina_bp_resolve_node(Pathogen$TransmissionProb, timestep, n, Ntimesteps, "TransmissionProb"))
  clear <- .ina_bp_clip01(.ina_bp_resolve_node(Pathogen$ClearanceProb, timestep, n, Ntimesteps, "ClearanceProb"))
  intro <- if (ForceNoPathogenIntroduction) rep(0, n) else
    .ina_bp_clip01(.ina_bp_resolve_node(Pathogen$IntroductionProb, timestep, n, Ntimesteps, "IntroductionProb"))
  pext <- .ina_bp_clip01(.ina_bp_resolve_node(Pathogen$PathogenHostExtinctionProb, timestep, n, Ntimesteps, "PathogenHostExtinctionProb"))
  pdet <- .ina_bp_clip01(.ina_bp_resolve_node(Pathogen$DetectionProb, timestep, n, Ntimesteps, "Pathogen DetectionProb"))

  Seam <- if (is.null(SEAM) || (length(SEAM) == 1L && isTRUE(SEAM == 0))) matrix(0, n, n) else as.matrix(SEAM)
  if (!all(dim(Seam) == c(n, n))) stop("SEAM must be nodes x nodes")
  Seam <- matrix(.ina_bp_clip01(Seam), n, n)

  ext_host <- if (isTRUE(OngoingExternalInvasion))
    .ina_bp_clip01(.ina_bp_resolve_node(InvasionRisk, timestep, n, Ntimesteps, "InvasionRisk")) else rep(0, n)
  ext_info <- if (isTRUE(OngoingExternalInfo))
    .ina_bp_clip01(.ina_bp_resolve_node(ExternalInfoProb, timestep, n, Ntimesteps, "ExternalInfoProb")) else rep(0, n)

  list(
    host_det = .ina_bp_clip01(.ina_bp_resolve_node(DetectionProb, timestep, n, Ntimesteps, "DetectionProb")),
    manage = .ina_bp_clip01(.ina_bp_resolve_node(ManageProb, timestep, n, Ntimesteps, "ManageProb")),
    erad = .ina_bp_clip01(.ina_bp_resolve_node(EradicationProb, timestep, n, Ntimesteps, "EradicationProb")),
    spread_red = .ina_bp_clip01(.ina_bp_resolve_node(SpreadReduction, timestep, n, Ntimesteps, "SpreadReduction")),
    info_ret = .ina_bp_clip01(.ina_bp_resolve_node(InfoRetentionProb, timestep, n, Ntimesteps, "InfoRetentionProb")),
    Kinfo = .ina_bp_resolve_node(InfoPersistenceSteps, timestep, n, Ntimesteps,
                                 "InfoPersistenceSteps", allow_na = TRUE, default = NA),
    survival = .ina_bp_clip01(.ina_bp_resolve_node(Survival, timestep, n, Ntimesteps, "Survival")),
    host_disp = host_disp, seam = Seam, ext_host = ext_host, ext_info = ext_info,
    contact = C, trans = tr, clear = clear, intro = intro, path_ext = pext,
    path_det = pdet)
}

.ina_bp_step_exact <- function(dist, par, InformationAcquisition,
                               UseInfoPersistence, OutsideNodes = integer(0)) {
  it <- .ina_bp_dist_items(dist)
  if (!length(it)) stop("Exact Binary pathogen distribution is empty")
  n <- length(it[[1L]]$state$H)
  host_acq <- InformationAcquisition %in% c("host", "both")
  path_acq <- InformationAcquisition %in% c("pathogen", "both")

  # Age the current local-evidence clock. Binary host-mode semantics refresh an
  # informed extant host to age 0 before management; pathogen-only mode does not.
  dist <- .ina_bp_apply_deterministic(dist, function(s) {
    if (UseInfoPersistence) {
      for (i in seq_len(n)) if (s$Info[i] == 1L) {
        if (host_acq && s$H[i] == 1L) s$Age[i] <- 0L
        else if (s$Age[i] >= 0L) s$Age[i] <- s$Age[i] + 1L
      }
    }
    s
  })

  # Management draws are retained temporarily because the same draw controls
  # eradication and spread reduction.
  dist <- .ina_bp_apply_deterministic(dist, function(s) { s$M <- integer(n); s })
  dist <- .ina_bp_branch_nodes(dist, seq_len(n),
    function(s, i) par$manage[i] * s$Info[i],
    function(s, i) .ina_bp_set(s, "M", i, 1L))

  # Host survival / management eradication. Pathogen cannot persist without host.
  dist <- .ina_bp_branch_nodes(dist, seq_len(n),
    function(s, i) if (s$H[i] == 1L) 1 - par$survival[i] * (1 - par$erad[i] * s$M[i]) else 0,
    function(s, i) {
      s <- .ina_bp_set(s, "H", i, 0L)
      .ina_bp_set(s, "P", i, 0L)
    })

  # Detected snapshot used for SEAM is taken after survival but before host spread
  # and before programmed stopping, exactly as in Binary INApest.
  dist <- .ina_bp_apply_deterministic(dist, function(s) {
    s$D <- as.integer(s$H == 1L & s$Info == 1L)
    s
  })

  # Host dispersal: conditional on surviving sources and their shared management
  # draws, target-level establishment events are independent across target columns.
  dist <- .ina_bp_branch_nodes(dist, seq_len(n),
    function(s, j) {
      if (s$H[j] == 1L) return(0)
      pe <- par$host_disp[, j] * s$H * (1 - s$M * par$spread_red)
      1 - prod(1 - .ina_bp_clip01(pe))
    },
    function(s, j) .ina_bp_set(s, "H", j, 1L))

  # Programmed stop or memoryless information retention.
  for (i in seq_len(n)) {
    nxt <- .ina_bp_dist_new()
    for (z in .ina_bp_dist_items(dist)) {
      s <- z$state; w <- z$prob
      if (s$Info[i] == 0L) {
        .ina_bp_dist_add(nxt, s, w); next
      }
      if (!is.na(par$Kinfo[i])) {
        stop_now <- s$Age[i] < 0L || s$Age[i] >= par$Kinfo[i]
        if (stop_now) {
          s$Info[i] <- 0L; s$Age[i] <- -1L
        }
        .ina_bp_dist_add(nxt, s, w)
      } else {
        pr <- par$info_ret[i]
        if (pr > 0) .ina_bp_dist_add(nxt, s, w * pr)
        if (pr < 1) {
          q <- s; q$Info[i] <- 0L; q$Age[i] <- -1L
          .ina_bp_dist_add(nxt, q, w * (1 - pr))
        }
      }
    }
    dist <- nxt
  }

  # SEAM transfer uses the pre-stop D snapshot.
  dist <- .ina_bp_branch_nodes(dist, seq_len(n),
    function(s, j) {
      if (s$Info[j] == 1L) return(0)
      1 - prod(1 - .ina_bp_clip01(par$seam[, j] * s$D))
    },
    function(s, j) .ina_bp_set(s, "Info", j, 1L))

  # External host and information arrivals.
  dist <- .ina_bp_branch_nodes(dist, seq_len(n),
    function(s, i) if (s$H[i] == 0L) par$ext_host[i] else 0,
    function(s, i) .ina_bp_set(s, "H", i, 1L))
  dist <- .ina_bp_branch_nodes(dist, seq_len(n),
    function(s, i) if (s$Info[i] == 0L) par$ext_info[i] else 0,
    function(s, i) .ina_bp_set(s, "Info", i, 1L))

  # Management and pre-stop detection snapshot are no longer required.
  dist <- .ina_bp_apply_deterministic(dist, function(s) { s$M <- NULL; s$D <- NULL; s })

  # Binary pathogen step: synchronous transmission -> clearance -> resident
  # introduction -> pathogen-associated host extinction. Freeze the infectious
  # source vector before any target is updated so infection cannot traverse two
  # contact edges within the same timestep.
  dist <- .ina_bp_apply_deterministic(dist, function(s) { s$Psource <- s$P; s })
  dist <- .ina_bp_branch_nodes(dist, seq_len(n),
    function(s, j) {
      if (s$H[j] == 0L || s$P[j] == 1L) return(0)
      q <- .ina_bp_clip01(par$contact[, j] * par$trans * s$Psource)
      1 - prod(1 - q)
    },
    function(s, j) .ina_bp_set(s, "P", j, 1L))
  dist <- .ina_bp_apply_deterministic(dist, function(s) { s$Psource <- NULL; s })

  dist <- .ina_bp_branch_nodes(dist, seq_len(n),
    function(s, i) if (s$P[i] == 1L) par$clear[i] else 0,
    function(s, i) .ina_bp_set(s, "P", i, 0L))

  dist <- .ina_bp_branch_nodes(dist, seq_len(n),
    function(s, i) if (s$H[i] == 1L && s$P[i] == 0L) par$intro[i] else 0,
    function(s, i) .ina_bp_set(s, "P", i, 1L))

  dist <- .ina_bp_branch_nodes(dist, seq_len(n),
    function(s, i) if (s$P[i] == 1L) par$path_ext[i] else 0,
    function(s, i) {
      s <- .ina_bp_set(s, "H", i, 0L)
      .ina_bp_set(s, "P", i, 0L)
    })

  # Escape is a first-passage event observed at the complete pathogen-state
  # timestep boundary, consistent with PathogenPresentResults.
  if (length(OutsideNodes)) dist <- .ina_bp_apply_deterministic(dist, function(s) {
    if (any(s$P[OutsideNodes] == 1L)) s$Escaped <- 1L
    s
  })

  # Expected detection events before they are merged through HaveInfo.
  path_det_expected <- numeric(n)
  for (z in .ina_bp_dist_items(dist))
    path_det_expected <- path_det_expected + z$prob * z$state$P * par$path_det

  if (path_acq) dist <- .ina_bp_branch_nodes(dist, seq_len(n),
    function(s, i) if (s$P[i] == 1L) par$path_det[i] else 0,
    function(s, i) {
      s <- .ina_bp_set(s, "Info", i, 1L)
      if (UseInfoPersistence) s <- .ina_bp_set(s, "Age", i, 0L)
      s
    })

  host_det_expected <- numeric(n)
  for (z in .ina_bp_dist_items(dist))
    host_det_expected <- host_det_expected + z$prob * z$state$H * par$host_det

  if (host_acq) dist <- .ina_bp_branch_nodes(dist, seq_len(n),
    function(s, i) if (s$H[i] == 1L) par$host_det[i] else 0,
    function(s, i) {
      s <- .ina_bp_set(s, "Info", i, 1L)
      if (UseInfoPersistence) s <- .ina_bp_set(s, "Age", i, 0L)
      s
    })

  list(distribution = dist,
       PathogenDetectionExpected = path_det_expected,
       HostDetectionExpected = host_det_expected)
}

.ina_bp_summary <- function(dist) {
  it <- .ina_bp_dist_items(dist)
  if (!length(it)) stop("Empty Binary pathogen state distribution")
  n <- length(it[[1L]]$state$H)
  host <- path <- info <- known_host <- numeric(n)
  anypath <- escape <- 0
  for (z in it) {
    s <- z$state; w <- z$prob
    host <- host + w * s$H
    path <- path + w * s$P
    info <- info + w * s$Info
    known_host <- known_host + w * s$H * s$Info
    anypath <- anypath + w * as.numeric(any(s$P == 1L))
    escape <- escape + w * s$Escaped
  }
  list(Mass = sum(vapply(it, function(z) z$prob, numeric(1))),
       HostPresence = host, PathogenPresence = path, Information = info,
       KnownHostPresence = known_host, AnyPathogen = anypath, Escape = escape,
       ExpectedHostNodeCount = sum(host), ExpectedPathogenNodeCount = sum(path),
       ExpectedInformedNodeCount = sum(info), StateCount = length(it))
}

.ina_bp_intrinsic_operator <- function(Pathogen, timestep, n, Ntimesteps,
                                       HostAvailability = 1) {
  C <- .ina_bp_contact(Pathogen, timestep, n, Ntimesteps)
  tr <- .ina_bp_clip01(.ina_bp_resolve_node(Pathogen$TransmissionProb, timestep, n, Ntimesteps, "TransmissionProb"))
  clear <- .ina_bp_clip01(.ina_bp_resolve_node(Pathogen$ClearanceProb, timestep, n, Ntimesteps, "ClearanceProb"))
  pext <- .ina_bp_clip01(.ina_bp_resolve_node(Pathogen$PathogenHostExtinctionProb, timestep, n, Ntimesteps, "PathogenHostExtinctionProb"))
  h <- .ina_bp_clip01(if (length(HostAvailability) == 1L) rep(HostAvailability, n) else HostAvailability)
  if (length(h) != n) stop("HostAvailability must be scalar or length nodes")
  B <- matrix(0, n, n)
  for (i in seq_len(n)) {
    B[i, i] <- (1 - clear[i]) * (1 - pext[i]) * h[i]
    for (j in seq_len(n)) if (j != i) {
      q <- .ina_bp_clip01(C[i, j] * tr[i])[1L]
      B[j, i] <- q * h[j] * (1 - clear[j]) * (1 - pext[j])
    }
  }
  B
}

.ina_bp_spectral_radius <- function(B) {
  if (!length(B)) return(0)
  max(Mod(eigen(B, only.values = TRUE)$values))
}

.ina_bp_branching_extinction <- function(B, horizon = 50L, tolerance = 1e-12) {
  B <- as.matrix(B); n <- nrow(B)
  q <- rep(0, n); hist <- matrix(NA_real_, n, horizon)
  for (tt in seq_len(horizon)) {
    old <- q; nq <- numeric(n)
    for (i in seq_len(n)) nq[i] <- prod((1 - B[, i]) + B[, i] * q)
    q <- .ina_bp_clip01(nq); hist[, tt] <- q
    if (max(abs(q - old)) < tolerance && tt < horizon) {
      hist[, (tt + 1L):horizon] <- q; break
    }
  }
  list(ProbabilityBySource = q, History = hist)
}

.ina_bp_branching_no_escape <- function(B, OutsideNodes, horizon = 10L) {
  B <- as.matrix(B); n <- nrow(B)
  OutsideNodes <- sort(unique(as.integer(OutsideNodes)))
  if (any(OutsideNodes < 1L | OutsideNodes > n)) stop("OutsideNodes contains invalid node indices")
  inside <- setdiff(seq_len(n), OutsideNodes)
  g <- rep(1, n); hist <- matrix(0, n, horizon)
  if (length(OutsideNodes)) g[OutsideNodes] <- 0
  for (tt in seq_len(horizon)) {
    ng <- g
    for (i in inside) {
      outside_no <- if (length(OutsideNodes)) prod(1 - B[OutsideNodes, i]) else 1
      inside_no <- if (length(inside)) prod((1 - B[inside, i]) + B[inside, i] * g[inside]) else 1
      ng[i] <- outside_no * inside_no
    }
    if (length(OutsideNodes)) ng[OutsideNodes] <- 0
    g <- .ina_bp_clip01(ng); hist[, tt] <- 1 - g
  }
  list(NoEscapeBySource = g, EscapeBySource = 1 - g, HistoryEscapeBySource = hist)
}

.ina_bp_exact_joint_one_step_operator <- function(InitialState, InitialInfo,
                                                   Pathogen, InformationAcquisition,
                                                   UseInfoPersistence,
                                                   step_par, n, OutsideNodes = integer(0)) {
  H <- as.numeric(InitialState); I <- as.numeric(InitialInfo)
  if (length(I) == 1L) I <- rep(I, n)
  if (length(H) != n || length(I) != n || any(!(H %in% c(0,1))) || any(!(I %in% c(0,1)))) return(NULL)
  active <- which(H == 1)
  if (!length(active)) return(matrix(0, n, n))
  G <- matrix(0, n, n)
  par0 <- step_par; par0$intro[] <- 0
  for (src in active) {
    Age <- rep(-1L, n); Age[I == 1 & H == 1] <- 0L
    s <- .ina_bp_state(H, as.integer(seq_len(n) == src), I, Age, 0L)
    d <- .ina_bp_dist_new(); .ina_bp_dist_add(d, s, 1)
    z <- .ina_bp_step_exact(d, par0, InformationAcquisition, UseInfoPersistence, OutsideNodes)
    G[, src] <- .ina_bp_summary(z$distribution)$PathogenPresence
  }
  G
}

INApestBinaryPathogenAnalytical <- function(
    Ntimesteps = 10,
    InitialState,
    InitialInfo = 0,
    InitialPathogen = NULL,
    Pathogen,
    InformationAcquisition = NULL,
    ApplyInitialDetection = TRUE,
    DetectionProb = 0,
    ManageProb = 0,
    EradicationProb = 0,
    SpreadReduction = 0,
    InfoRetentionProb = 1,
    InfoPersistenceSteps = NA,
    EnvEstabProb = 1,
    Survival = 1,
    SDDprob,
    LDDprob = 0,
    SEAM = 0,
    InvasionRisk = NA,
    ExternalInfoProb = NA,
    OngoingExternalInvasion = FALSE,
    OngoingExternalInfo = FALSE,
    OutsideNodes = integer(0),
    Exact = c("auto", "always", "never"),
    ExactMaxNodes = 3L,
    ExactMaxStates = 250000L,
    BranchingHorizon = Ntimesteps,
    ReturnStateDistribution = FALSE,
    ReturnOperators = FALSE) {

  if (!inherits(Pathogen, "INApestPathogen") || !identical(Pathogen$Model, "Binary"))
    stop("Pathogen must be INApestPathogen(Model = 'Binary')")
  if (length(Ntimesteps) != 1L || Ntimesteps < 1 || Ntimesteps != floor(Ntimesteps))
    stop("Ntimesteps must be a positive integer")
  Ntimesteps <- as.integer(Ntimesteps)
  Exact <- match.arg(Exact)

  SDD0 <- if (length(dim(SDDprob)) == 3L) SDDprob[, , 1L] else as.matrix(SDDprob)
  n <- nrow(SDD0)
  if (ncol(SDD0) != n) stop("SDDprob must be square")
  Pathogen$binary_validate(n, Ntimesteps)
  InformationAcquisition <- .ina_bp_information_acquisition(InformationAcquisition, Pathogen)
  host_acq <- InformationAcquisition %in% c("host", "both")
  path_acq <- InformationAcquisition %in% c("pathogen", "both")

  OutsideNodes <- sort(unique(as.integer(OutsideNodes)))
  if (length(OutsideNodes) && any(OutsideNodes < 1L | OutsideNodes > n))
    stop("OutsideNodes contains invalid node indices")

  UseInfoPersistence <- .ina_bp_use_persistence(InfoPersistenceSteps, n, Ntimesteps)
  diagnostics <- character(0)
  if (UseInfoPersistence) diagnostics <- c(diagnostics,
    paste0("Programmed stopping uses an explicit local-evidence age state. In InformationAcquisition='",
           InformationAcquisition, "' mode, only the selected biological evidence source(s) refresh that clock."))
  if (path_acq && !isTRUE(Pathogen$DetectionTriggersInfo) && !is.null(InformationAcquisition))
    diagnostics <- c(diagnostics,
      "Explicit InformationAcquisition supersedes the legacy Pathogen$DetectionTriggersInfo switch for this Binary analytical call.")
  diagnostics <- c(diagnostics,
    "Exact comparisons condition on supplied mean probabilities; set Binary INApest DetectionSD, ManageSD, EradicationSD and SpreadReductionSD to zero for stochastic validation.")

  # Intrinsic pathogen-engine operator: hosts are assumed available at the
  # pathogen update. External pathogen introductions are excluded.
  intrinsic_ops <- lapply(seq_len(Ntimesteps), function(tt)
    .ina_bp_intrinsic_operator(Pathogen, tt, n, Ntimesteps, HostAvailability = 1))
  intrinsic_cycle <- Reduce(`%*%`, rev(intrinsic_ops))
  intrinsic_cycle_lambda <- .ina_bp_spectral_radius(intrinsic_cycle)
  intrinsic_step_lambda <- intrinsic_cycle_lambda^(1 / Ntimesteps)

  use_exact <- switch(Exact,
    always = TRUE,
    never = FALSE,
    auto = n <= as.integer(ExactMaxNodes))
  if (Exact == "always" && n > as.integer(ExactMaxNodes)) diagnostics <- c(diagnostics,
    "Exact='always' overrides ExactMaxNodes; state-space growth can be very large.")

  trajectory <- NULL; exact_dist <- NULL; exact_status <- "not run"
  joint_operator <- NULL
  if (use_exact) {
    d <- .ina_bp_initial_distribution(InitialState, InitialInfo, InitialPathogen,
                                      Pathogen, DetectionProb, ApplyInitialDetection,
                                      InformationAcquisition, UseInfoPersistence,
                                      n, Ntimesteps, OutsideNodes)
    s0 <- .ina_bp_summary(d)
    host_mat <- path_mat <- info_mat <- known_host_mat <- matrix(NA_real_, n, Ntimesteps + 1L)
    host_mat[, 1L] <- s0$HostPresence; path_mat[, 1L] <- s0$PathogenPresence; info_mat[, 1L] <- s0$Information
    known_host_mat[, 1L] <- s0$KnownHostPresence
    any_path <- escape <- expected_host <- expected_path <- expected_info <- numeric(Ntimesteps + 1L)
    state_count <- integer(Ntimesteps + 1L)
    any_path[1L] <- s0$AnyPathogen; escape[1L] <- s0$Escape
    expected_host[1L] <- s0$ExpectedHostNodeCount; expected_path[1L] <- s0$ExpectedPathogenNodeCount; expected_info[1L] <- s0$ExpectedInformedNodeCount
    state_count[1L] <- s0$StateCount
    path_det <- host_det <- matrix(0, n, Ntimesteps)

    first_step_par <- NULL
    for (tt in seq_len(Ntimesteps)) {
      par <- .ina_bp_step_parameters(tt, n, Ntimesteps, Pathogen,
        DetectionProb, ManageProb, EradicationProb, SpreadReduction,
        InfoRetentionProb, InfoPersistenceSteps, EnvEstabProb, Survival,
        SDDprob, LDDprob, SEAM, OngoingExternalInvasion, InvasionRisk,
        OngoingExternalInfo, ExternalInfoProb)
      if (tt == 1L) first_step_par <- par
      z <- .ina_bp_step_exact(d, par, InformationAcquisition,
                              UseInfoPersistence, OutsideNodes)
      d <- z$distribution
      sm <- .ina_bp_summary(d)
      if (sm$StateCount > ExactMaxStates) {
        exact_status <- paste0("stopped after timestep ", tt, ": reachable state count ",
                               sm$StateCount, " exceeded ExactMaxStates=", ExactMaxStates)
        diagnostics <- c(diagnostics, exact_status)
        break
      }
      host_mat[, tt + 1L] <- sm$HostPresence
      path_mat[, tt + 1L] <- sm$PathogenPresence
      info_mat[, tt + 1L] <- sm$Information
      known_host_mat[, tt + 1L] <- sm$KnownHostPresence
      any_path[tt + 1L] <- sm$AnyPathogen; escape[tt + 1L] <- sm$Escape
      expected_host[tt + 1L] <- sm$ExpectedHostNodeCount
      expected_path[tt + 1L] <- sm$ExpectedPathogenNodeCount
      expected_info[tt + 1L] <- sm$ExpectedInformedNodeCount
      state_count[tt + 1L] <- sm$StateCount
      path_det[, tt] <- z$PathogenDetectionExpected
      host_det[, tt] <- z$HostDetectionExpected
      exact_status <- "complete"
    }

    completed <- max(which(!is.na(colSums(host_mat)))) - 1L
    keep <- seq_len(completed + 1L)
    trajectory <- list(
      Timestep = 0:completed,
      HostPresenceProb = host_mat[, keep, drop = FALSE],
      PathogenPresenceProb = path_mat[, keep, drop = FALSE],
      InformationProb = info_mat[, keep, drop = FALSE],
      KnownHostPresenceProb = known_host_mat[, keep, drop = FALSE],
      AnyPathogenProb = any_path[keep],
      PathogenAbsenceProb = 1 - any_path[keep],
      EscapeProb = escape[keep],
      ExpectedHostNodeCount = expected_host[keep],
      ExpectedPathogenNodeCount = expected_path[keep],
      ExpectedInformedNodeCount = expected_info[keep],
      ReachableStateCount = state_count[keep],
      ExpectedPathogenDetections = if (completed) path_det[, seq_len(completed), drop = FALSE] else path_det[, FALSE, drop = FALSE],
      ExpectedHostDetections = if (completed) host_det[, seq_len(completed), drop = FALSE] else host_det[, FALSE, drop = FALSE])

    if (!is.null(first_step_par)) joint_operator <-
      .ina_bp_exact_joint_one_step_operator(InitialState, InitialInfo, Pathogen,
        InformationAcquisition, UseInfoPersistence, first_step_par, n, OutsideNodes)
    if (ReturnStateDistribution && identical(exact_status, "complete")) exact_dist <- d
  } else diagnostics <- c(diagnostics,
    paste0("Exact finite-state propagation not run (n=", n, ", ExactMaxNodes=", ExactMaxNodes,
           "). Rare-pathogen operator/branching results are returned instead."))

  # Scalable rare-lineage branch uses the first intrinsic pathogen operator for
  # static interpretation and composes the time-varying operators for expected
  # lineage abundance. It deliberately does not hide host saturation/collisions.
  x0 <- if (is.null(InitialPathogen))
    .ina_bp_clip01(.ina_bp_resolve_node(Pathogen$InitialPresent, 1L, n, Ntimesteps, "InitialPresent")) else
    .ina_bp_clip01(.ina_bp_resolve_node(InitialPathogen, 1L, n, Ntimesteps, "InitialPathogen"))
  h0 <- .ina_bp_clip01(.ina_bp_resolve_node(InitialState, 1L, n, Ntimesteps, "InitialState"))
  x <- x0 * h0
  rare_traj <- matrix(NA_real_, n, Ntimesteps + 1L); rare_traj[, 1L] <- x
  for (tt in seq_len(Ntimesteps)) {
    x <- as.numeric(intrinsic_ops[[tt]] %*% x)
    rare_traj[, tt + 1L] <- x
  }

  B0 <- intrinsic_ops[[1L]]
  bext <- .ina_bp_branching_extinction(B0, as.integer(max(1L, BranchingHorizon)))
  extinction_branch <- if (all(x0 %in% c(0,1))) {
    prod(ifelse(x0 > 0, bext$ProbabilityBySource, 1))
  } else NA_real_

  escape_branch <- NULL
  if (length(OutsideNodes)) {
    bno <- .ina_bp_branching_no_escape(B0, OutsideNodes, as.integer(max(1L, BranchingHorizon)))
    inside <- setdiff(seq_len(n), OutsideNodes)
    no0 <- if (all(x0 %in% c(0,1))) {
      if (any(x0[OutsideNodes] > 0)) 0 else
        prod(ifelse(x0 > 0 & seq_len(n) %in% inside, bno$NoEscapeBySource, 1))
    } else NA_real_
    escape_branch <- list(
      Method = "rare-pathogen Bernoulli multitype branching no-escape approximation",
      ProbabilityByHorizon = if (is.na(no0)) NA_real_ else 1 - no0,
      EscapeBySource = bno$EscapeBySource,
      HistoryEscapeBySource = bno$HistoryEscapeBySource)
  }

  intro_any <- FALSE
  for (tt in seq_len(Ntimesteps)) {
    ip <- .ina_bp_resolve_node(Pathogen$IntroductionProb, tt, n, Ntimesteps, "IntroductionProb")
    if (any(ip > 0)) { intro_any <- TRUE; break }
  }
  exact_extinction <- if (!is.null(trajectory) && !intro_any)
    tail(trajectory$PathogenAbsenceProb, 1L) else NA_real_
  if (intro_any) diagnostics <- c(diagnostics,
    "IntroductionProb is non-zero: exact pathogen absence at the horizon is reported, but it is not labelled lineage extinction because reintroduction is possible.")
  if (length(OutsideNodes)) {
    outside_intro <- FALSE
    for (tt in seq_len(Ntimesteps)) {
      ip <- .ina_bp_resolve_node(Pathogen$IntroductionProb, tt, n, Ntimesteps, "IntroductionProb")
      if (any(ip[OutsideNodes] > 0)) { outside_intro <- TRUE; break }
    }
    if (outside_intro) diagnostics <- c(diagnostics,
      "OutsideNodes have non-zero IntroductionProb: exact Escape includes first outside pathogen presence from exogenous resident infection as well as cross-boundary transmission. Set outside IntroductionProb=0 when the estimand is strictly containment leakage.")
    diagnostics <- c(diagnostics,
      "Exact Escape is first outside pathogen presence at the completed pathogen-state timestep boundary; an outside infection that is created and clears within the same pathogen step is not observable in PathogenPresentResults and is not counted.")
  }

  exact_escape <- if (!is.null(trajectory) && length(OutsideNodes)) tail(trajectory$EscapeProb, 1L) else NA_real_
  escape_increment <- if (!is.null(trajectory) && length(OutsideNodes)) diff(trajectory$EscapeProb) else numeric(0)
  conditional_escape_time <- if (length(escape_increment) && sum(escape_increment) > 0)
    sum(seq_along(escape_increment) * escape_increment) / sum(escape_increment) else NA_real_

  joint_lambda <- if (!is.null(joint_operator)) .ina_bp_spectral_radius(joint_operator) else NA_real_
  result <- list(
    Model = "INApest",
    Variant = "binary_pathogen",
    InformationAcquisition = InformationAcquisition,
    InformationState = "single shared HaveInfo state",
    Exact = list(
      Status = exact_status,
      Method = if (use_exact) "event-order exact finite-state propagation" else "not run",
      Trajectory = trajectory,
      PathogenFreeByHorizon = if (is.null(trajectory)) NA_real_ else tail(trajectory$PathogenAbsenceProb, 1L),
      ExtinctionByHorizon = exact_extinction,
      EscapeByHorizon = exact_escape,
      ConditionalMeanTimeToFirstEscape = conditional_escape_time,
      StateDistribution = exact_dist),
    Growth = list(
      PathogenEngineIntrinsicPerTimestepMultiplier = intrinsic_step_lambda,
      PathogenEngineCycleMultiplier = intrinsic_cycle_lambda,
      JointHostPathogenOneStepMultiplier = joint_lambda,
      IntrinsicOperators = if (ReturnOperators) intrinsic_ops else NULL,
      JointOneStepOperator = if (ReturnOperators) joint_operator else NULL,
      Interpretation = c(
        "PathogenEngineIntrinsic excludes exogenous IntroductionProb and conditions on hosts being available at the pathogen update.",
        "JointHostPathogenOneStepMultiplier is returned only when InitialState and InitialInfo are deterministic; it includes first-step host survival/management before the pathogen update.")),
    RarePathogen = list(
      ExpectedLineageTrajectory = rare_traj,
      BranchingExtinctionByHorizon = extinction_branch,
      Escape = escape_branch,
      Scope = "rare pathogen / independent lineages; finite-node collisions and saturation are not represented"),
    HeadlineEstimands = list(
      GrowthMultiplier = intrinsic_step_lambda,
      PathogenPresenceByHorizon = if (is.null(trajectory)) NA_real_ else tail(trajectory$AnyPathogenProb, 1L),
      PathogenAbsenceByHorizon = if (is.null(trajectory)) NA_real_ else tail(trajectory$PathogenAbsenceProb, 1L),
      ExtinctionByHorizon = exact_extinction,
      EscapeByHorizon = exact_escape),
    Diagnostics = unique(diagnostics),
    ApproximationScope = c(
      "Exact finite-state results are exact for the represented Binary event process and supplied deterministic probabilities; state-space size grows rapidly with nodes and persistence ages.",
      "Rare-pathogen growth/branching results are scalable screening approximations and become less exact as multiple lineages collide in the same finite node set.",
      "Parameter SD integration is not hidden inside the analytical solution; compare to stochastic runs with SD=0 or treat SD uncertainty separately."))
  class(result) <- c("INApestBinaryPathogenAnalyticalResult", "list")
  result
}

print.INApestBinaryPathogenAnalyticalResult <- function(x, ...) {
  cat("INApest Binary pathogen analytical result\n")
  cat("  Information acquisition:", x$InformationAcquisition, "\n")
  cat("  Intrinsic pathogen multiplier:",
      format(round(x$Growth$PathogenEngineIntrinsicPerTimestepMultiplier, 5), nsmall = 5), "\n")
  cat("  Exact finite-state status:", x$Exact$Status, "\n")
  if (is.finite(x$HeadlineEstimands$PathogenPresenceByHorizon))
    cat("  Pathogen present at horizon:", round(x$HeadlineEstimands$PathogenPresenceByHorizon, 5), "\n")
  if (is.finite(x$HeadlineEstimands$EscapeByHorizon))
    cat("  Pathogen escape by horizon:", round(x$HeadlineEstimands$EscapeByHorizon, 5), "\n")
  if (length(x$Diagnostics)) cat("  Diagnostics:", length(x$Diagnostics), "(see $Diagnostics)\n")
  invisible(x)
}

###############################################################################
### Final dispatcher extension: preserve every pre-existing analytical call.
###############################################################################
INApestAnalytical_pre_binary_pathogen <- INApestAnalytical

INApestAnalytical <- function(...) {
  args <- list(...)
  Model <- if (!is.null(args$Model)) as.character(args$Model)[1L] else "INApest"
  Pathogen <- args$Pathogen
  if (!identical(Model, "INApest") || is.null(Pathogen))
    return(do.call(INApestAnalytical_pre_binary_pathogen, args))
  if (!inherits(Pathogen, "INApestPathogen") || !identical(Pathogen$Model, "Binary"))
    stop("For Model='INApest', Pathogen must be NULL or INApestPathogen(Model='Binary')")

  keep <- c("Ntimesteps", "InitialState", "InitialInfo", "InitialPathogen", "Pathogen",
            "InformationAcquisition", "ApplyInitialDetection", "DetectionProb",
            "ManageProb", "EradicationProb", "SpreadReduction", "InfoRetentionProb",
            "InfoPersistenceSteps", "EnvEstabProb", "Survival", "SDDprob", "LDDprob",
            "SEAM", "InvasionRisk", "ExternalInfoProb", "OngoingExternalInvasion",
            "OngoingExternalInfo", "OutsideNodes", "Exact", "ExactMaxNodes",
            "ExactMaxStates", "BranchingHorizon", "ReturnStateDistribution",
            "ReturnOperators")
  unknown <- setdiff(names(args), c("Model", keep))
  if (length(unknown)) {
    # Legacy host-only arguments that have no Binary pathogen meaning are not
    # silently interpreted. Common generic arguments can still be ignored only
    # if they are exact defaults in the old interface.
    harmless <- c("LDDrate", "K", "PropaguleProduction", "PropaguleEstablishment",
                  "Transition", "Nstages", "SeedbankK", "MortalityProb",
                  "DispersalDensityFactor", "ExportProb", "ExportSDDprob",
                  "ExportLDDprob", "OutsideEstablishmentProb", "AssumeResidualExport",
                  "ExtinctionGenerations", "InformationMode")
    bad <- setdiff(unknown, harmless)
    if (length(bad)) stop("Unsupported Binary pathogen analytical argument(s): ", paste(bad, collapse = ", "))
  }
  call_args <- args[intersect(names(args), keep)]
  do.call(INApestBinaryPathogenAnalytical, call_args)
}

###############################################################################
### Meta pathogen analytical extension, 27 August 2026
###############################################################################
###############################################################################
### INApest Meta pathogen analytical extension: SIS, fixed host abundance
### Development milestone: 27 August 2026
###
### Purpose
###   Exact finite-state solution for one-node INApestMeta SIS when total host
###   abundance is fixed, plus the rare-pathogen linear growth operator for a
###   fixed-abundance node network.
###
### This module mirrors the current INApestPathogen SIS event contract:
###   1. transmission is calculated from the start-of-pathogen-step S/I state;
###   2. exogenous IntroductionProb acts after transmission and converts
###      remaining susceptible residents to I without changing N;
###   3. recovery acts on the pre-existing I only; newly infected hosts do not
###      recover in the same disease step;
###   4. PathogenMortalityProb must be zero for this fixed-N exact branch.
###############################################################################

.ina_meta_sis_clip01 <- function(x) pmin(1, pmax(0, x))

.ina_meta_sis_step_value <- function(x, timestep, Ntimesteps, name,
                                     integer = FALSE, positive = FALSE) {
  if (is.function(x)) {
    fm <- names(formals(x))
    a <- list(timestep = timestep, Ntimesteps = Ntimesteps)
    if (!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a), fm)]
    x <- do.call(x, a)
  }
  if (!is.null(dim(x))) stop(name, " must be scalar, length Ntimesteps, or resolver function for the one-node exact SIS solution")
  x <- as.numeric(x)
  if (length(x) == 1L) z <- x
  else if (length(x) == Ntimesteps) z <- x[timestep]
  else stop(name, " must be scalar or length Ntimesteps")
  if (length(z) != 1L || !is.finite(z)) stop(name, " must resolve to one finite value")
  if (integer && (z < 0 || z != floor(z))) stop(name, " must resolve to a non-negative whole number")
  if (positive && z <= 0) stop(name, " must resolve to a positive value")
  z
}

.ina_meta_sis_one_node_operator <- function(
    N,
    Beta,
    RecoveryProb,
    IntroductionProb = 0,
    IntroductionNumber = 1,
    Transmission = c("frequency", "density"),
    DensityScale = 1,
    ContactWeight = 1,
    timestep = 1L,
    Ntimesteps = 1L) {

  Transmission <- match.arg(Transmission)
  N <- as.integer(N)
  if (length(N) != 1L || is.na(N) || N < 0L) stop("N must be one non-negative whole number")

  beta <- .ina_meta_sis_step_value(Beta, timestep, Ntimesteps, "Beta")
  rec <- .ina_meta_sis_step_value(RecoveryProb, timestep, Ntimesteps, "RecoveryProb")
  intro_p <- .ina_meta_sis_step_value(IntroductionProb, timestep, Ntimesteps, "IntroductionProb")
  intro_n <- as.integer(.ina_meta_sis_step_value(IntroductionNumber, timestep, Ntimesteps, "IntroductionNumber", integer = TRUE))
  density_scale <- .ina_meta_sis_step_value(DensityScale, timestep, Ntimesteps, "DensityScale", positive = TRUE)
  cw <- .ina_meta_sis_step_value(ContactWeight, timestep, Ntimesteps, "ContactWeight")

  if (beta < 0) stop("Beta must be non-negative")
  if (rec < 0 || rec > 1) stop("RecoveryProb must be in [0,1]")
  if (intro_p < 0 || intro_p > 1) stop("IntroductionProb must be in [0,1]")
  if (cw < 0) stop("ContactWeight must be non-negative")

  states <- 0:N
  T <- matrix(0, nrow = N + 1L, ncol = N + 1L,
              dimnames = list(paste0("I", states), paste0("I", states)))

  for (ii in states) {
    S <- N - ii
    if (Transmission == "frequency") {
      # With one node, a positive scalar contact weight cancels from I/N.
      foi <- if (N > 0L && cw > 0) beta * ii / N else 0
    } else {
      foi <- beta * ii * cw / density_scale
    }
    p_inf <- .ina_meta_sis_clip01(-expm1(-max(0, foi)))

    px <- dbinom(0:S, size = S, prob = p_inf)
    py <- dbinom(0:ii, size = ii, prob = 1 - rec)

    for (x in 0:S) {
      max_intro <- min(intro_n, S - x)
      intro_values <- if (intro_p <= 0 || max_intro <= 0L) 0L else c(0L, max_intro)
      intro_probs <- if (intro_p <= 0 || max_intro <= 0L) 1 else c(1 - intro_p, intro_p)
      for (y in 0:ii) {
        base_p <- px[x + 1L] * py[y + 1L]
        if (base_p == 0) next
        for (kk in seq_along(intro_values)) {
          jj <- y + x + intro_values[kk]
          T[ii + 1L, jj + 1L] <- T[ii + 1L, jj + 1L] + base_p * intro_probs[kk]
        }
      }
    }
  }

  rs <- rowSums(T)
  if (max(abs(rs - 1)) > 1e-12) stop("Internal SIS operator construction error: transition rows do not sum to 1")
  attr(T, "N") <- N
  attr(T, "timestep") <- timestep
  attr(T, "Transmission") <- Transmission
  T
}

INApestMetaSISExactOneNode <- function(
    Ntimesteps,
    HostPopulation,
    InitialInfected,
    Beta,
    RecoveryProb,
    IntroductionProb = 0,
    IntroductionNumber = 1,
    Transmission = c("frequency", "density"),
    DensityScale = 1,
    ContactWeight = 1,
    DetectionProb = 0,
    ReturnOperators = FALSE) {

  Transmission <- match.arg(Transmission)
  Ntimesteps <- as.integer(Ntimesteps)
  N <- as.integer(HostPopulation)
  I0 <- as.integer(InitialInfected)
  if (length(Ntimesteps) != 1L || is.na(Ntimesteps) || Ntimesteps < 0L) stop("Ntimesteps must be a non-negative whole number")
  if (length(N) != 1L || is.na(N) || N < 0L) stop("HostPopulation must be a non-negative whole number")
  if (length(I0) != 1L || is.na(I0) || I0 < 0L || I0 > N) stop("InitialInfected must be a whole number between 0 and HostPopulation")

  dist <- numeric(N + 1L); dist[I0 + 1L] <- 1
  traj <- matrix(0, nrow = Ntimesteps + 1L, ncol = N + 1L,
                 dimnames = list(timestep = 0:Ntimesteps, I = 0:N))
  traj[1L, ] <- dist
  ops <- vector("list", Ntimesteps)

  expected_I <- numeric(Ntimesteps + 1L)
  prevalence <- numeric(Ntimesteps + 1L)
  extinction <- numeric(Ntimesteps + 1L)
  pdetect <- numeric(Ntimesteps + 1L)
  expected_I[1L] <- I0
  prevalence[1L] <- if (N > 0L) I0 / N else 0
  extinction[1L] <- as.numeric(I0 == 0L)

  detection_at <- function(tt, d) {
    pd <- .ina_meta_sis_step_value(DetectionProb, max(1L, tt), max(1L, Ntimesteps), "DetectionProb")
    if (pd < 0 || pd > 1) stop("DetectionProb must resolve to [0,1]")
    sum(d * (1 - (1 - pd)^(0:N)))
  }
  pdetect[1L] <- detection_at(1L, dist)

  if (Ntimesteps > 0L) {
    for (tt in seq_len(Ntimesteps)) {
      T <- .ina_meta_sis_one_node_operator(
        N = N, Beta = Beta, RecoveryProb = RecoveryProb,
        IntroductionProb = IntroductionProb, IntroductionNumber = IntroductionNumber,
        Transmission = Transmission, DensityScale = DensityScale,
        ContactWeight = ContactWeight, timestep = tt, Ntimesteps = Ntimesteps)
      ops[[tt]] <- T
      dist <- as.numeric(dist %*% T)
      traj[tt + 1L, ] <- dist
      expected_I[tt + 1L] <- sum((0:N) * dist)
      prevalence[tt + 1L] <- if (N > 0L) expected_I[tt + 1L] / N else 0
      extinction[tt + 1L] <- dist[1L]
      pdetect[tt + 1L] <- detection_at(tt, dist)
    }
  }

  # Intrinsic rare-pathogen multiplier excludes external introduction.
  b1 <- .ina_meta_sis_step_value(Beta, 1L, max(1L, Ntimesteps), "Beta")
  r1 <- .ina_meta_sis_step_value(RecoveryProb, 1L, max(1L, Ntimesteps), "RecoveryProb")
  c1 <- .ina_meta_sis_step_value(ContactWeight, 1L, max(1L, Ntimesteps), "ContactWeight")
  ds1 <- .ina_meta_sis_step_value(DensityScale, 1L, max(1L, Ntimesteps), "DensityScale", positive = TRUE)
  lambda <- if (Transmission == "frequency") {
    (1 - r1) + if (N > 0L && c1 > 0) b1 else 0
  } else {
    (1 - r1) + N * b1 * c1 / ds1
  }

  out <- list(
    Model = "INApestMeta",
    PathogenModel = "SIS",
    HostAssumption = "fixed abundance",
    Exact = TRUE,
    HostPopulation = N,
    InitialInfected = I0,
    Transmission = Transmission,
    Growth = list(
      IntrinsicRarePathogenMultiplier = lambda,
      Classification = if (lambda > 1 + 1e-12) "growing" else if (lambda < 1 - 1e-12) "declining" else "threshold"
    ),
    Trajectory = data.frame(
      timestep = 0:Ntimesteps,
      ExpectedInfected = expected_I,
      ExpectedSusceptible = N - expected_I,
      ExpectedPrevalence = prevalence,
      ExtinctionProbability = extinction,
      PathogenPresenceProbability = 1 - extinction,
      DetectionProbability = pdetect,
      stringsAsFactors = FALSE
    ),
    StateDistribution = traj,
    Diagnostics = c(
      "Exact finite-state one-node SIS solution for fixed integer host abundance.",
      "Transmission and recovery use the start-of-pathogen-step state; newly infected hosts do not recover in the same timestep.",
      "IntroductionProb, when non-zero, infects susceptible residents and does not change host abundance.",
      "PathogenMortalityProb is outside this fixed-host branch and must be handled by the host-turnover extension."
    )
  )
  if (ReturnOperators) out$Operators <- ops
  class(out) <- c("INApestMetaSISExactOneNode", "list")
  out
}

INApestMetaSISGrowthOperatorFixedN <- function(
    HostPopulation,
    Beta,
    RecoveryProb,
    ContactMatrix = NULL,
    Transmission = c("frequency", "density"),
    DensityScale = 1) {

  Transmission <- match.arg(Transmission)
  N <- as.numeric(HostPopulation)
  n <- length(N)
  if (!n || any(!is.finite(N)) || any(N < 0)) stop("HostPopulation must contain finite non-negative values")

  recycle_n <- function(x, name) {
    x <- as.numeric(x)
    if (length(x) == 1L) x <- rep(x, n)
    if (length(x) != n || any(!is.finite(x))) stop(name, " must be scalar or length nodes")
    x
  }
  beta <- recycle_n(Beta, "Beta")
  rec <- recycle_n(RecoveryProb, "RecoveryProb")
  ds <- recycle_n(DensityScale, "DensityScale")
  if (any(beta < 0)) stop("Beta must be non-negative")
  if (any(rec < 0 | rec > 1)) stop("RecoveryProb must be in [0,1]")
  if (any(ds <= 0)) stop("DensityScale must be positive")

  C <- if (is.null(ContactMatrix)) diag(n) else as.matrix(ContactMatrix)
  if (!all(dim(C) == c(n, n)) || any(!is.finite(C)) || any(C < 0))
    stop("ContactMatrix must be a finite non-negative nodes x nodes source-by-target matrix")

  A <- matrix(0, n, n)
  diag(A) <- 1 - rec

  if (Transmission == "frequency") {
    denom <- as.numeric(crossprod(N, C))
    for (j in seq_len(n)) {
      if (N[j] <= 0 || denom[j] <= 0 || beta[j] == 0) next
      A[j, ] <- A[j, ] + N[j] * beta[j] * C[, j] / denom[j]
    }
  } else {
    for (j in seq_len(n)) {
      if (N[j] <= 0 || beta[j] == 0) next
      A[j, ] <- A[j, ] + N[j] * beta[j] * C[, j] / ds[j]
    }
  }

  ev <- eigen(A, only.values = TRUE)$values
  lambda <- max(Mod(ev))
  list(
    Operator = A,
    IntrinsicRarePathogenMultiplier = lambda,
    Classification = if (lambda > 1 + 1e-12) "growing" else if (lambda < 1 - 1e-12) "declining" else "threshold",
    Orientation = "rows = recipient nodes; columns = infectious source nodes",
    Diagnostics = "Linearisation is about the pathogen-free fixed-host state; external pathogen introduction is excluded."
  )
}

print.INApestMetaSISExactOneNode <- function(x, ...) {
  cat("INApestMeta SIS exact one-node analytical result\n")
  cat("  Host population:", x$HostPopulation, "\n")
  cat("  Initial infected:", x$InitialInfected, "\n")
  cat("  Transmission:", x$Transmission, "\n")
  cat("  Rare-pathogen multiplier:", format(round(x$Growth$IntrinsicRarePathogenMultiplier, 6), nsmall = 6), "\n")
  cat("  End expected infected:", round(tail(x$Trajectory$ExpectedInfected, 1), 6), "\n")
  cat("  End pathogen presence probability:", round(tail(x$Trajectory$PathogenPresenceProbability, 1), 6), "\n")
  invisible(x)
}

###############################################################################
### Exact one-node SIS with host turnover (joint N, I state)
###############################################################################

.ina_meta_sis_joint_states <- function(K) {
  K <- as.integer(K)
  do.call(rbind, lapply(0:K, function(N) data.frame(N=N, I=0:N)))
}

.ina_meta_sis_external_probs <- function(x) {
  if (is.null(x)) return(c(S=1, I=0))
  if (!is.numeric(x) || is.null(names(x)) || any(!names(x) %in% c("S","I")))
    stop("ExternalPathogenStateProb for SIS must be a named numeric vector using S and/or I")
  out <- c(S=0,I=0); out[names(x)] <- as.numeric(x)
  if (any(!is.finite(out)) || any(out < 0) || abs(sum(out)-1) > 1e-10)
    stop("ExternalPathogenStateProb must contain non-negative probabilities summing to 1")
  out
}

.ina_meta_sis_turnover_operator <- function(
    K,
    Survival = 1,
    RecruitToCapacity = FALSE,
    ExternalHostInvasionProb = 0,
    ExternalHostNumber = 1,
    ExternalPathogenStateProb = c(S=1,I=0),
    Beta = 0,
    RecoveryProb = 0,
    PathogenMortalityProb = 0,
    IntroductionProb = 0,
    IntroductionNumber = 1,
    Transmission = c("frequency","density"),
    DensityScale = 1,
    ContactWeight = 1,
    timestep = 1L,
    Ntimesteps = 1L) {

  Transmission <- match.arg(Transmission)
  K <- as.integer(K)
  if (length(K)!=1L || is.na(K) || K < 0L) stop("K must be one non-negative whole number")
  surv <- .ina_meta_sis_step_value(Survival,timestep,Ntimesteps,"Survival")
  ext_p <- .ina_meta_sis_step_value(ExternalHostInvasionProb,timestep,Ntimesteps,"ExternalHostInvasionProb")
  ext_n <- as.integer(.ina_meta_sis_step_value(ExternalHostNumber,timestep,Ntimesteps,"ExternalHostNumber",integer=TRUE))
  beta <- .ina_meta_sis_step_value(Beta,timestep,Ntimesteps,"Beta")
  rec <- .ina_meta_sis_step_value(RecoveryProb,timestep,Ntimesteps,"RecoveryProb")
  pmort <- .ina_meta_sis_step_value(PathogenMortalityProb,timestep,Ntimesteps,"PathogenMortalityProb")
  intro_p <- .ina_meta_sis_step_value(IntroductionProb,timestep,Ntimesteps,"IntroductionProb")
  intro_n <- as.integer(.ina_meta_sis_step_value(IntroductionNumber,timestep,Ntimesteps,"IntroductionNumber",integer=TRUE))
  ds <- .ina_meta_sis_step_value(DensityScale,timestep,Ntimesteps,"DensityScale",positive=TRUE)
  cw <- .ina_meta_sis_step_value(ContactWeight,timestep,Ntimesteps,"ContactWeight")
  if (surv<0 || surv>1) stop("Survival must be in [0,1]")
  if (ext_p<0 || ext_p>1) stop("ExternalHostInvasionProb must be in [0,1]")
  if (beta<0 || cw<0) stop("Beta and ContactWeight must be non-negative")
  if (rec<0 || pmort<0 || rec+pmort>1+1e-12) stop("RecoveryProb and PathogenMortalityProb must be non-negative and sum to <= 1")
  if (intro_p<0 || intro_p>1) stop("IntroductionProb must be in [0,1]")
  ep <- .ina_meta_sis_external_probs(ExternalPathogenStateProb)

  states <- .ina_meta_sis_joint_states(K)
  key <- paste(states$N,states$I,sep=":")
  idx <- setNames(seq_len(nrow(states)),key)
  T <- matrix(0,nrow(states),nrow(states),dimnames=list(key,key))

  add_pathogen_step <- function(row, Npre, Ipre, host_prob) {
    if (host_prob == 0) return(invisible(NULL))
    Spre <- Npre - Ipre
    if (Transmission=="frequency") foi <- if(Npre>0 && cw>0) beta*Ipre/Npre else 0
    else foi <- beta*Ipre*cw/ds
    pinf <- .ina_meta_sis_clip01(-expm1(-max(0,foi)))
    px <- dbinom(0:Spre,Spre,pinf)
    # Pre-existing infectious hosts have stay/recover/death categorical outcomes.
    for (x in 0:Spre) {
      max_intro <- min(intro_n,Spre-x)
      ivals <- if(intro_p<=0 || max_intro<=0L) 0L else c(0L,max_intro)
      iprobs <- if(intro_p<=0 || max_intro<=0L) 1 else c(1-intro_p,intro_p)
      for (stay in 0:Ipre) {
        for (death in 0:(Ipre-stay)) {
          recover <- Ipre-stay-death
          pfate <- dmultinom(c(stay,recover,death),prob=c(1-rec-pmort,rec,pmort))
          if (pfate==0) next
          for (kk in seq_along(ivals)) {
            Nend <- Npre-death
            Iend <- stay+x+ivals[kk]
            T[row,idx[[paste(Nend,Iend,sep=":")]]] <<-
              T[row,idx[[paste(Nend,Iend,sep=":")]]] + host_prob*px[x+1L]*pfate*iprobs[kk]
          }
        }
      }
    }
    invisible(NULL)
  }

  for (row in seq_len(nrow(states))) {
    N0 <- states$N[row]; I0 <- states$I[row]; S0 <- N0-I0

    # The simulator's Binomial(total N, Survival) followed by hypergeometric
    # reconciliation is distributionally equivalent to independent survival of
    # S and I hosts at the same Survival probability.
    if (isTRUE(RecruitToCapacity)) {
      # Susceptible losses are replaced before the pathogen step, so only the
      # number of surviving infectious hosts matters to the pre-disease state.
      for (ih in 0:I0) {
        phost <- dbinom(ih,I0,surv)
        Nlocal <- K; Ilocal <- ih
        accepted <- min(ext_n,max(0L,K-Nlocal))
        # When filled to capacity no external hosts can be accepted.
        if (accepted==0L || ext_p<=0) add_pathogen_step(row,Nlocal,Ilocal,phost)
        else stop("Internal turnover logic error")
      }
    } else {
      for (ih in 0:I0) for (sh in 0:S0) {
        phost <- dbinom(ih,I0,surv)*dbinom(sh,S0,surv)
        if (phost==0) next
        Nlocal <- ih+sh; Ilocal <- ih
        accepted <- min(ext_n,max(0L,K-Nlocal))
        if (ext_p<=0 || accepted<=0L) {
          add_pathogen_step(row,Nlocal,Ilocal,phost)
        } else {
          # No external event.
          add_pathogen_step(row,Nlocal,Ilocal,phost*(1-ext_p))
          # External event; only actually accepted hosts are assigned state.
          for (ii in 0:accepted) {
            pimm <- dbinom(ii,accepted,ep["I"])
            add_pathogen_step(row,Nlocal+accepted,Ilocal+ii,phost*ext_p*pimm)
          }
        }
      }
    }
  }

  rs <- rowSums(T)
  if(max(abs(rs-1))>1e-11) stop("Internal turnover SIS operator error: rows do not sum to 1; max error=",max(abs(rs-1)))
  attr(T,"States") <- states
  attr(T,"K") <- K
  T
}

INApestMetaSISExactTurnoverOneNode <- function(
    Ntimesteps,
    K,
    InitialHostPopulation,
    InitialInfected,
    Survival = 1,
    RecruitToCapacity = FALSE,
    ExternalHostInvasionProb = 0,
    ExternalHostNumber = 1,
    ExternalPathogenStateProb = c(S=1,I=0),
    Beta = 0,
    RecoveryProb = 0,
    PathogenMortalityProb = 0,
    IntroductionProb = 0,
    IntroductionNumber = 1,
    Transmission = c("frequency","density"),
    DensityScale = 1,
    ContactWeight = 1,
    DetectionProb = 0,
    ReturnOperators = FALSE) {

  Transmission <- match.arg(Transmission)
  K <- as.integer(K); Ntimesteps <- as.integer(Ntimesteps)
  N0 <- as.integer(InitialHostPopulation); I0 <- as.integer(InitialInfected)
  if (K<0 || N0<0 || N0>K || I0<0 || I0>N0) stop("Require 0 <= InitialInfected <= InitialHostPopulation <= K")
  states <- .ina_meta_sis_joint_states(K); keys <- paste(states$N,states$I,sep=":")
  dist <- numeric(nrow(states)); dist[match(paste(N0,I0,sep=":"),keys)] <- 1
  trajdist <- matrix(0,Ntimesteps+1L,nrow(states),dimnames=list(timestep=0:Ntimesteps,state=keys));trajdist[1,]<-dist
  ops <- vector("list",Ntimesteps)

  summarise <- function(d,tt) {
    pd <- .ina_meta_sis_step_value(DetectionProb,max(1L,tt),max(1L,Ntimesteps),"DetectionProb")
    if(pd<0 || pd>1) stop("DetectionProb must be in [0,1]")
    c(ExpectedHost=sum(d*states$N),ExpectedInfected=sum(d*states$I),
      PathogenPresenceProbability=sum(d[states$I>0]),HostPresenceProbability=sum(d[states$N>0]),
      DetectionProbability=sum(d*(1-(1-pd)^states$I)))
  }
  sm <- matrix(0,Ntimesteps+1L,5);sm[1,]<-summarise(dist,1L)
  if(Ntimesteps>0) for(tt in seq_len(Ntimesteps)) {
    T <- .ina_meta_sis_turnover_operator(K,Survival,RecruitToCapacity,ExternalHostInvasionProb,
      ExternalHostNumber,ExternalPathogenStateProb,Beta,RecoveryProb,PathogenMortalityProb,
      IntroductionProb,IntroductionNumber,Transmission,DensityScale,ContactWeight,tt,Ntimesteps)
    ops[[tt]]<-T;dist<-as.numeric(dist%*%T);trajdist[tt+1L,]<-dist;sm[tt+1L,]<-summarise(dist,tt)
  }
  out <- list(Model="INApestMeta",PathogenModel="SIS",HostAssumption="exact one-node turnover",
    Exact=TRUE,K=K,RecruitToCapacity=isTRUE(RecruitToCapacity),
    Trajectory=data.frame(timestep=0:Ntimesteps,ExpectedHost=sm[,1],ExpectedInfected=sm[,2],
      ExpectedPrevalence=ifelse(sm[,1]>0,sm[,2]/sm[,1],0),PathogenPresenceProbability=sm[,3],
      HostPresenceProbability=sm[,4],DetectionProbability=sm[,5]),
    StateTable=states,StateDistribution=trajdist,
    Diagnostics=c("Exact finite-state joint (N,I) solution for one-node SIS host turnover.",
      "Natural host mortality thins susceptible and infectious hosts before local recruitment.",
      if(isTRUE(RecruitToCapacity)) "Local recruitment restores host abundance to K with susceptible hosts before pathogen transmission." else "No local recruitment is applied in this exact branch.",
      "External host immigration occurs before pathogen transmission; accepted immigrants may be S or I.",
      "IntroductionProb remains resident S-to-I conversion and never changes N."))
  if(ReturnOperators) out$Operators<-ops
  class(out)<-c("INApestMetaSISExactTurnoverOneNode","list")
  out
}

###############################################################################
### Exact small-network SIS with fixed host abundance
###############################################################################

.ina_meta_sis_network_states <- function(N) {
  N <- as.integer(N)
  g <- expand.grid(lapply(N, function(z) 0:z), KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  out <- as.matrix(g); storage.mode(out) <- "integer"
  colnames(out) <- paste0("node",seq_along(N))
  out
}

.ina_meta_sis_network_operator_fixedN <- function(
    HostPopulation,
    Beta,
    RecoveryProb,
    ContactMatrix = NULL,
    IntroductionProb = 0,
    IntroductionNumber = 1,
    Transmission = c("frequency","density"),
    DensityScale = 1) {
  Transmission <- match.arg(Transmission)
  N <- as.integer(HostPopulation); n <- length(N)
  if(any(is.na(N)) || any(N<0)) stop("HostPopulation must be non-negative whole numbers")
  recyc <- function(x,name,integer=FALSE){x<-as.numeric(x);if(length(x)==1)x<-rep(x,n);if(length(x)!=n||any(!is.finite(x)))stop(name," must be scalar or length nodes");if(integer&&any(x<0|x!=floor(x)))stop(name," must be non-negative whole numbers");x}
  beta<-recyc(Beta,"Beta");rec<-recyc(RecoveryProb,"RecoveryProb");ip<-recyc(IntroductionProb,"IntroductionProb");inum<-as.integer(recyc(IntroductionNumber,"IntroductionNumber",TRUE));ds<-recyc(DensityScale,"DensityScale")
  if(any(beta<0)||any(rec<0|rec>1)||any(ip<0|ip>1)||any(ds<=0))stop("Invalid SIS parameters")
  C<-if(is.null(ContactMatrix))diag(n) else as.matrix(ContactMatrix)
  if(!all(dim(C)==c(n,n))||any(!is.finite(C))||any(C<0))stop("ContactMatrix must be finite non-negative nodes x nodes")
  states<-.ina_meta_sis_network_states(N); keys<-apply(states,1,paste,collapse=":"); idx<-setNames(seq_len(nrow(states)),keys)
  T<-matrix(0,nrow(states),nrow(states),dimnames=list(keys,keys))
  denom<-as.numeric(crossprod(N,C))
  for(row in seq_len(nrow(states))){
    I<-states[row,];S<-N-I;pressure<-as.numeric(crossprod(I,C))
    foi<-if(Transmission=="frequency") beta*ifelse(denom>0,pressure/denom,0) else beta*pressure/ds
    pinf<-.ina_meta_sis_clip01(-expm1(-pmax(0,foi)))
    node_probs<-vector("list",n)
    for(j in seq_len(n)){
      pnext<-numeric(N[j]+1L)
      px<-dbinom(0:S[j],S[j],pinf[j]);py<-dbinom(0:I[j],I[j],1-rec[j])
      for(x in 0:S[j]){
        mx<-min(inum[j],S[j]-x); iv<-if(ip[j]<=0||mx<=0)0L else c(0L,mx); pp<-if(ip[j]<=0||mx<=0)1 else c(1-ip[j],ip[j])
        for(y in 0:I[j])for(k in seq_along(iv))pnext[y+x+iv[k]+1L]<-pnext[y+x+iv[k]+1L]+px[x+1L]*py[y+1L]*pp[k]
      }
      node_probs[[j]]<-pnext
    }
    # Conditional on the start state, target-node infection/recovery draws are independent.
    dest<-states
    probs<-rep(1,nrow(dest))
    for(j in seq_len(n)) probs<-probs*node_probs[[j]][dest[,j]+1L]
    T[row,]<-probs
  }
  if(max(abs(rowSums(T)-1))>1e-11)stop("Network SIS operator rows do not sum to 1")
  attr(T,"States")<-states;attr(T,"HostPopulation")<-N;T
}

INApestMetaSISExactNetworkFixedN <- function(
    Ntimesteps,
    HostPopulation,
    InitialInfected,
    Beta,
    RecoveryProb,
    ContactMatrix = NULL,
    IntroductionProb = 0,
    IntroductionNumber = 1,
    Transmission = c("frequency","density"),
    DensityScale = 1,
    OutsideNodes = integer(0),
    ReturnOperator = FALSE) {
  Transmission<-match.arg(Transmission);N<-as.integer(HostPopulation);I0<-as.integer(InitialInfected);n<-length(N)
  if(length(I0)!=n||any(I0<0|I0>N))stop("InitialInfected must match HostPopulation and satisfy 0 <= I <= N")
  OutsideNodes<-as.integer(OutsideNodes);if(any(!OutsideNodes%in%seq_len(n)))stop("OutsideNodes contains invalid node indices")
  T<-.ina_meta_sis_network_operator_fixedN(N,Beta,RecoveryProb,ContactMatrix,IntroductionProb,IntroductionNumber,Transmission,DensityScale)
  states<-attr(T,"States");keys<-rownames(T);dist<-numeric(nrow(states));dist[match(paste(I0,collapse=":"),keys)]<-1
  full<-matrix(0,Ntimesteps+1L,nrow(states));full[1,]<-dist
  expI<-matrix(0,Ntimesteps+1L,n);expI[1,]<-as.numeric(dist%*%states)
  ext<-numeric(Ntimesteps+1L);ext[1]<-sum(dist[rowSums(states)==0])
  escape<-numeric(Ntimesteps+1L);safe<-dist
  escaped_states<-if(length(OutsideNodes))apply(states[,OutsideNodes,drop=FALSE],1,function(z)any(z>0)) else rep(FALSE,nrow(states))
  if(length(OutsideNodes)&&any(I0[OutsideNodes]>0)){escape[1]<-1;safe[]<-0}else safe[escaped_states]<-0
  if(Ntimesteps>0)for(tt in seq_len(Ntimesteps)){
    dist<-as.numeric(dist%*%T);full[tt+1,]<-dist;expI[tt+1,]<-as.numeric(dist%*%states);ext[tt+1]<-sum(dist[rowSums(states)==0])
    if(length(OutsideNodes)){
      nxt<-as.numeric(safe%*%T);newesc<-sum(nxt[escaped_states]);escape[tt+1]<-escape[tt]+newesc;nxt[escaped_states]<-0;safe<-nxt
    }
  }
  growth<-INApestMetaSISGrowthOperatorFixedN(N,Beta,RecoveryProb,ContactMatrix,Transmission,DensityScale)
  out<-list(Model="INApestMeta",PathogenModel="SIS",HostAssumption="fixed abundance network",Exact=TRUE,
    Growth=growth,ExpectedInfectedByNode=expI,ExpectedTotalInfected=rowSums(expI),
    ExtinctionProbability=ext,PathogenPresenceProbability=1-ext,
    Escape=if(length(OutsideNodes))list(OutsideNodes=OutsideNodes,ProbabilityByTimestep=escape,ProbabilityByHorizon=tail(escape,1)) else NULL,
    StateTable=states,StateDistribution=full,
    Diagnostics=c("Exact finite-state small-network SIS solution with fixed host counts.","Pathogen transmission uses the start-of-step infectious state, so newly infected nodes/hosts do not transmit again until the next timestep.","Escape is first passage: once any OutsideNodes contains I>0, that realisation is counted as escaped even if infection later clears."))
  if(ReturnOperator)out$Operator<-T;class(out)<-c("INApestMetaSISExactNetworkFixedN","list");out
}

###############################################################################
### Exact one-node SIS + management/information state
###############################################################################

.ina_meta_sis_info_states <- function(Khost, InfoPersistenceSteps=NA) {
  b <- .ina_meta_sis_joint_states(Khost)
  rows <- list(); k <- 0L
  finite <- !is.na(InfoPersistenceSteps)
  for(r in seq_len(nrow(b))) {
    k<-k+1L; rows[[k]]<-data.frame(N=b$N[r],I=b$I[r],H=0L,A=-2L)
    if(finite) {
      k<-k+1L; rows[[k]]<-data.frame(N=b$N[r],I=b$I[r],H=1L,A=-1L) # informed, no local-evidence clock
      ages <- 0:max(0L,as.integer(InfoPersistenceSteps)-1L)
      for(a in unique(ages)){k<-k+1L;rows[[k]]<-data.frame(N=b$N[r],I=b$I[r],H=1L,A=as.integer(a))}
    } else {
      k<-k+1L; rows[[k]]<-data.frame(N=b$N[r],I=b$I[r],H=1L,A=-2L)
    }
  }
  do.call(rbind,rows)
}

.ina_meta_class_mortality_outcomes <- function(count, survival, management_mortality) {
  count<-as.integer(count)
  ps<-survival*(1-management_mortality); pm<-survival*management_mortality; pn<-1-survival
  out<-list();k<-0L
  for(sv in 0:count)for(md in 0:(count-sv)){
    nd<-count-sv-md;p<-dmultinom(c(sv,md,nd),prob=c(ps,pm,pn))
    if(p>0){k<-k+1L;out[[k]]<-c(survive=sv,manage_death=md,prob=p)}
  }
  do.call(rbind,out)
}

.ina_meta_sis_disease_outcomes <- function(Npre,Ipre,beta,rec,pmort,intro_p,intro_n,Transmission,DensityScale,ContactWeight) {
  Spre<-Npre-Ipre
  if(Transmission=="frequency")foi<-if(Npre>0&&ContactWeight>0)beta*Ipre/Npre else 0 else foi<-beta*Ipre*ContactWeight/DensityScale
  pinf<-.ina_meta_sis_clip01(-expm1(-max(0,foi)));px<-dbinom(0:Spre,Spre,pinf)
  acc<-new.env(hash=TRUE,parent=emptyenv())
  add<-function(N,I,p){key<-paste(N,I,sep=":");old<-if(exists(key,acc,inherits=FALSE))get(key,acc) else 0;assign(key,old+p,acc)}
  for(x in 0:Spre){mx<-min(intro_n,Spre-x);iv<-if(intro_p<=0||mx<=0)0L else c(0L,mx);ipp<-if(intro_p<=0||mx<=0)1 else c(1-intro_p,intro_p)
    for(stay in 0:Ipre)for(death in 0:(Ipre-stay)){recover<-Ipre-stay-death;pf<-dmultinom(c(stay,recover,death),prob=c(1-rec-pmort,rec,pmort));if(pf==0)next
      for(kk in seq_along(iv))add(Npre-death,stay+x+iv[kk],px[x+1L]*pf*ipp[kk])}}
  ks<-ls(acc);do.call(rbind,lapply(ks,function(k){z<-strsplit(k,":",fixed=TRUE)[[1]];c(N=as.integer(z[1]),I=as.integer(z[2]),prob=get(k,acc))}))
}

.ina_meta_sis_info_operator <- function(
    Khost,
    Survival=1,
    ManageProb=0,
    MortalityProb=0,
    HostDetectionProb=0,
    PathogenDetectionProb=0,
    InformationAcquisition=c("host","pathogen","both"),
    InfoPersistenceSteps=NA,
    InfoRetentionProb=1,
    Beta=0,
    RecoveryProb=0,
    PathogenMortalityProb=0,
    IntroductionProb=0,
    IntroductionNumber=1,
    Transmission=c("frequency","density"),
    DensityScale=1,
    ContactWeight=1) {
  InformationAcquisition<-match.arg(InformationAcquisition);Transmission<-match.arg(Transmission)
  host_acq<-InformationAcquisition%in%c("host","both");path_acq<-InformationAcquisition%in%c("pathogen","both")
  vals<-c(Survival,ManageProb,MortalityProb,HostDetectionProb,PathogenDetectionProb,InfoRetentionProb,RecoveryProb,PathogenMortalityProb,IntroductionProb)
  if(any(!is.finite(vals))||any(vals<0)||any(vals>1)||RecoveryProb+PathogenMortalityProb>1)stop("Probability inputs must be valid")
  if(Beta<0||DensityScale<=0||ContactWeight<0)stop("Invalid transmission inputs")
  if(!is.na(InfoPersistenceSteps)&&(InfoPersistenceSteps<0||InfoPersistenceSteps!=floor(InfoPersistenceSteps)))stop("InfoPersistenceSteps must be NA or a non-negative whole number")
  st<-.ina_meta_sis_info_states(Khost,InfoPersistenceSteps);key<-paste(st$N,st$I,st$H,st$A,sep=":");idx<-setNames(seq_len(nrow(st)),key);T<-matrix(0,nrow(st),nrow(st),dimnames=list(key,key))
  finite<-!is.na(InfoPersistenceSteps);Kinfo<-if(finite)as.integer(InfoPersistenceSteps) else NA_integer_
  add_state<-function(row,N,I,H,A,p){ky<-paste(N,I,H,A,sep=":");T[row,idx[[ky]]]<<-T[row,idx[[ky]]]+p}
  for(row in seq_len(nrow(st))){N0<-st$N[row];I0<-st$I[row];S0<-N0-I0;H0<-st$H[row];A0<-st$A[row]
    mprob<-ManageProb*H0
    for(managing in 0:1){pmng<-if(managing==1)mprob else 1-mprob;if(pmng==0)next;mm<-MortalityProb*managing
      so<-.ina_meta_class_mortality_outcomes(S0,Survival,mm);io<-.ina_meta_class_mortality_outcomes(I0,Survival,mm)
      for(si in seq_len(nrow(so)))for(ii in seq_len(nrow(io))){ph<-pmng*so[si,"prob"]*io[ii,"prob"];if(ph==0)next
        Ssurv<-as.integer(so[si,"survive"]);Isurv<-as.integer(io[ii,"survive"]);Npre<-Ssurv+Isurv
        host_evidence_now<-host_acq && (so[si,"manage_death"]+io[ii,"manage_death"]>0)
        # Information loss is evaluated after this timestep's management/mortality.
        if(finite){
          if(H0==0){Hmid<-0L;Amid<--2L
          } else if(host_evidence_now){age_stop<-0L;if(age_stop>=Kinfo){Hmid<-0L;Amid<--2L}else{Hmid<-1L;Amid<-0L}
          } else if(A0<0){Hmid<-0L;Amid<--2L
          } else {age_stop<-A0+1L;if(age_stop>=Kinfo){Hmid<-0L;Amid<--2L}else{Hmid<-1L;Amid<-age_stop}}
          info_branches<-matrix(c(Hmid,Amid,1),nrow=1,dimnames=list(NULL,c("H","A","p")))
        } else {
          if(H0==1L&&InfoRetentionProb<1)info_branches<-rbind(c(1,-2,InfoRetentionProb),c(0,-2,1-InfoRetentionProb)) else info_branches<-matrix(c(H0,-2,1),nrow=1)
          colnames(info_branches)<-c("H","A","p")
        }
        douts<-.ina_meta_sis_disease_outcomes(Npre,Isurv,Beta,RecoveryProb,PathogenMortalityProb,IntroductionProb,as.integer(IntroductionNumber),Transmission,DensityScale,ContactWeight)
        for(ib in seq_len(nrow(info_branches)))for(dd in seq_len(nrow(douts))){pbase<-ph*info_branches[ib,"p"]*douts[dd,"prob"];if(pbase==0)next
          Ne<-as.integer(douts[dd,"N"]);Ie<-as.integer(douts[dd,"I"]);Hm<-as.integer(info_branches[ib,"H"]);Am<-as.integer(info_branches[ib,"A"])
          qh<-if(host_acq)1-(1-HostDetectionProb)^Ne else 0;qp<-if(path_acq)1-(1-PathogenDetectionProb)^Ie else 0;qev<-1-(1-qh)*(1-qp)
          if(qev>0)add_state(row,Ne,Ie,1L,if(finite)0L else -2L,pbase*qev)
          if(qev<1)add_state(row,Ne,Ie,Hm,Am,pbase*(1-qev))
        }
      }
    }
  }
  if(max(abs(rowSums(T)-1))>1e-10)stop("Information SIS operator rows do not sum to 1; max error=",max(abs(rowSums(T)-1)))
  attr(T,"States")<-st;T
}

INApestMetaSISExactInformationOneNode <- function(
    Ntimesteps,K,InitialHostPopulation,InitialInfected,InitialInfo=0,
    Survival=1,ManageProb=0,MortalityProb=0,HostDetectionProb=0,PathogenDetectionProb=0,
    InformationAcquisition=c("host","pathogen","both"),InfoPersistenceSteps=NA,InfoRetentionProb=1,
    Beta=0,RecoveryProb=0,PathogenMortalityProb=0,IntroductionProb=0,IntroductionNumber=1,
    Transmission=c("frequency","density"),DensityScale=1,ContactWeight=1,ReturnOperator=FALSE) {
  InformationAcquisition<-match.arg(InformationAcquisition);Transmission<-match.arg(Transmission);K<-as.integer(K);N0<-as.integer(InitialHostPopulation);I0<-as.integer(InitialInfected)
  if(N0<0||N0>K||I0<0||I0>N0||!InitialInfo%in%c(0,1))stop("Invalid initial state")
  T<-.ina_meta_sis_info_operator(K,Survival,ManageProb,MortalityProb,HostDetectionProb,PathogenDetectionProb,InformationAcquisition,InfoPersistenceSteps,InfoRetentionProb,Beta,RecoveryProb,PathogenMortalityProb,IntroductionProb,IntroductionNumber,Transmission,DensityScale,ContactWeight)
  st<-attr(T,"States");keys<-rownames(T);finite<-!is.na(InfoPersistenceSteps);host_acq<-InformationAcquisition%in%c("host","both");path_acq<-InformationAcquisition%in%c("pathogen","both")
  qh<-if(host_acq)1-(1-HostDetectionProb)^N0 else 0;qp<-if(path_acq)1-(1-PathogenDetectionProb)^I0 else 0;qev<-1-(1-qh)*(1-qp)
  dist<-numeric(nrow(st));
  # Direct initial local evidence establishes age zero; user-supplied information
  # without local evidence is deliberately unanchored under programmed stopping.
  if(qev>0)dist[match(paste(N0,I0,1,if(finite)0 else -2,sep=":"),keys)]<-qev
  if(qev<1){H<-as.integer(InitialInfo);A<-if(H==1&&finite)-1L else -2L;dist[match(paste(N0,I0,H,A,sep=":"),keys)]<-dist[match(paste(N0,I0,H,A,sep=":"),keys)]+(1-qev)}
  traj<-matrix(0,Ntimesteps+1L,nrow(st));traj[1,]<-dist;EN<-EI<-PH<-numeric(Ntimesteps+1L);EN[1]<-N0;EI[1]<-I0;PH[1]<-sum(dist*st$H);PM<-numeric(Ntimesteps)
  if(Ntimesteps>0)for(tt in seq_len(Ntimesteps)){PM[tt]<-ManageProb*sum(dist*st$H);dist<-as.numeric(dist%*%T);traj[tt+1,]<-dist;EN[tt+1]<-sum(dist*st$N);EI[tt+1]<-sum(dist*st$I);PH[tt+1]<-sum(dist*st$H)}
  out<-list(Model="INApestMeta",PathogenModel="SIS",Exact=TRUE,InformationAcquisition=InformationAcquisition,
    Trajectory=data.frame(timestep=0:Ntimesteps,ExpectedHost=EN,ExpectedInfected=EI,InformationProbability=PH),ManagingProbability=PM,StateTable=st,StateDistribution=traj,
    Diagnostics=c("Exact one-node joint host/pathogen/information state for small K.","HaveInfo remains the single management-control state; InformationAcquisition determines which direct local evidence creates or refreshes it.",if(finite)"Programmed stopping is represented by explicit time-since-local-evidence states." else "Information loss uses the memoryless InfoRetentionProb pathway."))
  if(ReturnOperator)out$Operator<-T;class(out)<-c("INApestMetaSISExactInformationOneNode","list");out
}

###############################################################################
### Exact small-network SIS + management/information for small host capacities
###############################################################################

.ina_meta_sis_disease_outcomes_pinf <- function(Npre,Ipre,pinf,rec,pmort,intro_p,intro_n) {
  Spre<-Npre-Ipre;px<-dbinom(0:Spre,Spre,.ina_meta_sis_clip01(pinf));acc<-new.env(hash=TRUE,parent=emptyenv())
  add<-function(N,I,p){key<-paste(N,I,sep=":");old<-if(exists(key,acc,inherits=FALSE))get(key,acc) else 0;assign(key,old+p,acc)}
  for(x in 0:Spre){mx<-min(intro_n,Spre-x);iv<-if(intro_p<=0||mx<=0)0L else c(0L,mx);ipp<-if(intro_p<=0||mx<=0)1 else c(1-intro_p,intro_p)
    for(stay in 0:Ipre)for(death in 0:(Ipre-stay)){recover<-Ipre-stay-death;pf<-dmultinom(c(stay,recover,death),prob=c(1-rec-pmort,rec,pmort));if(pf==0)next
      for(kk in seq_along(iv))add(Npre-death,stay+x+iv[kk],px[x+1L]*pf*ipp[kk])}}
  ks<-ls(acc);if(!length(ks))return(matrix(c(Npre,0,1),nrow=1,dimnames=list(NULL,c("N","I","prob"))))
  z<-do.call(rbind,lapply(ks,function(k){q<-strsplit(k,":",fixed=TRUE)[[1]];c(N=as.integer(q[1]),I=as.integer(q[2]),prob=get(k,acc))}));z
}

.ina_meta_sis_info_prebranches <- function(state_row, Survival, ManageProb, MortalityProb,
                                            host_acq, finite, Kinfo, InfoRetentionProb) {
  N0<-state_row[["N"]];I0<-state_row[["I"]];S0<-N0-I0;H0<-state_row[["H"]];A0<-state_row[["A"]]
  out<-list();k<-0L;mprob<-ManageProb*H0
  for(managing in 0:1){pmng<-if(managing==1)mprob else 1-mprob;if(pmng==0)next;mm<-MortalityProb*managing
    so<-.ina_meta_class_mortality_outcomes(S0,Survival,mm);io<-.ina_meta_class_mortality_outcomes(I0,Survival,mm)
    for(si in seq_len(nrow(so)))for(ii in seq_len(nrow(io))){ph<-pmng*so[si,"prob"]*io[ii,"prob"];if(ph==0)next
      Npre<-as.integer(so[si,"survive"]+io[ii,"survive"]);Ipre<-as.integer(io[ii,"survive"]);evidence<-host_acq&&(so[si,"manage_death"]+io[ii,"manage_death"]>0)
      ib<-list()
      if(finite){
        if(H0==0){ib[[1]]<-c(H=0,A=-2,p=1)
        }else if(evidence){if(0>=Kinfo)ib[[1]]<-c(H=0,A=-2,p=1) else ib[[1]]<-c(H=1,A=0,p=1)
        }else if(A0<0){ib[[1]]<-c(H=0,A=-2,p=1)
        }else{ag<-A0+1;if(ag>=Kinfo)ib[[1]]<-c(H=0,A=-2,p=1) else ib[[1]]<-c(H=1,A=ag,p=1)}
      }else{
        if(H0==1&&InfoRetentionProb<1){ib[[1]]<-c(H=1,A=-2,p=InfoRetentionProb);ib[[2]]<-c(H=0,A=-2,p=1-InfoRetentionProb)}else ib[[1]]<-c(H=H0,A=-2,p=1)
      }
      for(q in ib){k<-k+1L;out[[k]]<-c(N=Npre,I=Ipre,H=as.integer(q["H"]),A=as.integer(q["A"]),prob=unname(ph*q["p"]))}
    }}
  do.call(rbind,out)
}

INApestMetaSISExactNetworkInformation <- function(
    Ntimesteps,K,InitialHostPopulation,InitialInfected,InitialInfo=0,
    Survival=1,ManageProb=0,MortalityProb=0,HostDetectionProb=0,PathogenDetectionProb=0,
    InformationAcquisition=c("host","pathogen","both"),InfoPersistenceSteps=NA,InfoRetentionProb=1,
    Beta=0,RecoveryProb=0,PathogenMortalityProb=0,IntroductionProb=0,IntroductionNumber=1,
    ContactMatrix=NULL,Transmission=c("frequency","density"),DensityScale=1,OutsideNodes=integer(0),ReturnOperator=FALSE,
    MaxStates=2000L) {
  InformationAcquisition<-match.arg(InformationAcquisition);Transmission<-match.arg(Transmission)
  K<-as.integer(K);n<-length(K);N0<-as.integer(InitialHostPopulation);I0<-as.integer(InitialInfected);Hinit<-as.integer(InitialInfo)
  if(length(N0)!=n||length(I0)!=n||length(Hinit)==1)Hinit<-rep(Hinit,n)
  if(length(N0)!=n||length(I0)!=n||length(Hinit)!=n||any(N0<0|N0>K)|any(I0<0|I0>N0)|any(!Hinit%in%c(0,1)))stop("Invalid initial network state")
  recyc<-function(x,name){x<-as.numeric(x);if(length(x)==1)x<-rep(x,n);if(length(x)!=n||any(!is.finite(x)))stop(name," must be scalar or length nodes");x}
  surv<-recyc(Survival,"Survival");mp<-recyc(ManageProb,"ManageProb");mort<-recyc(MortalityProb,"MortalityProb");hd<-recyc(HostDetectionProb,"HostDetectionProb");pd<-recyc(PathogenDetectionProb,"PathogenDetectionProb");beta<-recyc(Beta,"Beta");rec<-recyc(RecoveryProb,"RecoveryProb");pmort<-recyc(PathogenMortalityProb,"PathogenMortalityProb");ip<-recyc(IntroductionProb,"IntroductionProb");inum<-as.integer(recyc(IntroductionNumber,"IntroductionNumber"));ds<-recyc(DensityScale,"DensityScale")
  if(any(surv<0|surv>1|mp<0|mp>1|mort<0|mort>1|hd<0|hd>1|pd<0|pd>1|rec<0|pmort<0|rec+pmort>1|ip<0|ip>1)||any(beta<0)|any(ds<=0))stop("Invalid parameters")
  C<-if(is.null(ContactMatrix))diag(n) else as.matrix(ContactMatrix);if(!all(dim(C)==c(n,n))||any(!is.finite(C))||any(C<0))stop("Invalid ContactMatrix")
  finite<-!is.na(InfoPersistenceSteps);if(finite){Kinfo<-as.integer(InfoPersistenceSteps);if(Kinfo<0)stop("InfoPersistenceSteps must be non-negative")}else Kinfo<-NA_integer_
  if(length(InfoRetentionProb)!=1||!is.finite(InfoRetentionProb)||InfoRetentionProb<0||InfoRetentionProb>1)stop("InfoRetentionProb must be one probability")
  host_acq<-InformationAcquisition%in%c("host","both");path_acq<-InformationAcquisition%in%c("pathogen","both")
  local_states<-lapply(K,function(k).ina_meta_sis_info_states(k,InfoPersistenceSteps));dims<-vapply(local_states,nrow,integer(1));nst<-prod(dims);if(nst>MaxStates)stop("Exact network information state space has ",nst," states; increase MaxStates deliberately or use an approximation")
  grid<-as.matrix(expand.grid(lapply(dims,seq_len),KEEP.OUT.ATTRS=FALSE));storage.mode(grid)<-"integer";gkey<-apply(grid,1,paste,collapse=":");gidx<-setNames(seq_len(nrow(grid)),gkey)
  decode<-function(gr){lapply(seq_len(n),function(j)local_states[[j]][gr[j],,drop=FALSE])}
  T<-matrix(0,nst,nst,dimnames=list(gkey,gkey))
  # cache local pre-management branches by node/local-state index
  precache<-lapply(seq_len(n),function(j)lapply(seq_len(dims[j]),function(s).ina_meta_sis_info_prebranches(local_states[[j]][s,],surv[j],mp[j],mort[j],host_acq,finite,Kinfo,InfoRetentionProb)))
  for(row in seq_len(nst)){
    prelists<-lapply(seq_len(n),function(j)precache[[j]][[grid[row,j]]]);pg<-as.matrix(expand.grid(lapply(prelists,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(pg)<-"integer"
    for(pb in seq_len(nrow(pg))){prs<-lapply(seq_len(n),function(j)prelists[[j]][pg[pb,j],]);p_pre<-prod(vapply(prs,function(z)z["prob"],numeric(1)));if(p_pre==0)next
      Nv<-vapply(prs,function(z)as.numeric(z["N"]),numeric(1));Iv<-vapply(prs,function(z)as.numeric(z["I"]),numeric(1));den<-as.numeric(crossprod(Nv,C));pressure<-as.numeric(crossprod(Iv,C));foi<-if(Transmission=="frequency")beta*ifelse(den>0,pressure/den,0)else beta*pressure/ds;pinf<-.ina_meta_sis_clip01(-expm1(-pmax(0,foi)))
      dlists<-lapply(seq_len(n),function(j).ina_meta_sis_disease_outcomes_pinf(as.integer(Nv[j]),as.integer(Iv[j]),pinf[j],rec[j],pmort[j],ip[j],inum[j]));dg<-as.matrix(expand.grid(lapply(dlists,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(dg)<-"integer"
      for(db in seq_len(nrow(dg))){dos<-lapply(seq_len(n),function(j)dlists[[j]][dg[db,j],]);p_dis<-prod(vapply(dos,function(z)z["prob"],numeric(1)));if(p_dis==0)next
        nextlists<-vector("list",n)
        for(j in seq_len(n)){Ne<-as.integer(dos[[j]]["N"]);Ie<-as.integer(dos[[j]]["I"]);Hm<-as.integer(prs[[j]]["H"]);Am<-as.integer(prs[[j]]["A"]);qh<-if(host_acq)1-(1-hd[j])^Ne else 0;qp<-if(path_acq)1-(1-pd[j])^Ie else 0;qev<-1-(1-qh)*(1-qp)
          a<-list();kk<-0L;if(qev>0){kk<-kk+1L;a[[kk]]<-c(N=Ne,I=Ie,H=1,A=if(finite)0 else -2,prob=qev)};if(qev<1){kk<-kk+1L;a[[kk]]<-c(N=Ne,I=Ie,H=Hm,A=Am,prob=1-qev)};nextlists[[j]]<-do.call(rbind,a)}
        ng<-as.matrix(expand.grid(lapply(nextlists,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(ng)<-"integer"
        for(nb in seq_len(nrow(ng))){locidx<-integer(n);pnext<-p_pre*p_dis;for(j in seq_len(n)){zz<-nextlists[[j]][ng[nb,j],];pnext<-pnext*zz["prob"];matchrow<-which(local_states[[j]]$N==zz["N"]&local_states[[j]]$I==zz["I"]&local_states[[j]]$H==zz["H"]&local_states[[j]]$A==zz["A"]);locidx[j]<-matchrow[1]};T[row,gidx[[paste(locidx,collapse=":")]]]<-T[row,gidx[[paste(locidx,collapse=":")]]]+pnext}
      }
    }
  }
  if(max(abs(rowSums(T)-1))>1e-9)stop("Exact network information operator row error ",max(abs(rowSums(T)-1)))
  # initial distribution after direct local detections
  initlists<-vector("list",n);for(j in seq_len(n)){qh<-if(host_acq)1-(1-hd[j])^N0[j] else 0;qp<-if(path_acq)1-(1-pd[j])^I0[j] else 0;qev<-1-(1-qh)*(1-qp);a<-list();kk<-0L;if(qev>0){kk<-kk+1L;a[[kk]]<-c(N=N0[j],I=I0[j],H=1,A=if(finite)0 else -2,prob=qev)};if(qev<1){kk<-kk+1L;a[[kk]]<-c(N=N0[j],I=I0[j],H=Hinit[j],A=if(Hinit[j]==1&&finite)-1 else -2,prob=1-qev)};initlists[[j]]<-do.call(rbind,a)}
  ig<-as.matrix(expand.grid(lapply(initlists,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));dist<-numeric(nst);for(b in seq_len(nrow(ig))){loc<-integer(n);p0<-1;for(j in seq_len(n)){z<-initlists[[j]][ig[b,j],];p0<-p0*z["prob"];loc[j]<-which(local_states[[j]]$N==z["N"]&local_states[[j]]$I==z["I"]&local_states[[j]]$H==z["H"]&local_states[[j]]$A==z["A"])[1]};dist[gidx[[paste(loc,collapse=":")]]]<-dist[gidx[[paste(loc,collapse=":")]]]+p0}
  # decoded matrices for summaries
  Nmat<-Imat<-Hmat<-matrix(0,nst,n);for(r in seq_len(nst))for(j in seq_len(n)){z<-local_states[[j]][grid[r,j],];Nmat[r,j]<-z$N;Imat[r,j]<-z$I;Hmat[r,j]<-z$H}
  EN<-EI<-EH<-matrix(0,Ntimesteps+1,n);EN[1,]<-as.numeric(dist%*%Nmat);EI[1,]<-as.numeric(dist%*%Imat);EH[1,]<-as.numeric(dist%*%Hmat);PM<-matrix(0,Ntimesteps,n)
  OutsideNodes<-as.integer(OutsideNodes);escaped<-if(length(OutsideNodes))apply(Imat[,OutsideNodes,drop=FALSE],1,function(z)any(z>0))else rep(FALSE,nst);escape<-numeric(Ntimesteps+1);safe<-dist;if(length(OutsideNodes)&&any(I0[OutsideNodes]>0)){escape[1]<-1;safe[]<-0}else safe[escaped]<-0
  for(tt in seq_len(Ntimesteps)){PM[tt,]<-as.numeric(dist%*%sweep(Hmat,2,mp,'*'));dist<-as.numeric(dist%*%T);EN[tt+1,]<-as.numeric(dist%*%Nmat);EI[tt+1,]<-as.numeric(dist%*%Imat);EH[tt+1,]<-as.numeric(dist%*%Hmat);if(length(OutsideNodes)){nx<-as.numeric(safe%*%T);escape[tt+1]<-escape[tt]+sum(nx[escaped]);nx[escaped]<-0;safe<-nx}}
  out<-list(Model="INApestMeta",PathogenModel="SIS",Exact=TRUE,InformationAcquisition=InformationAcquisition,ExpectedHostByNode=EN,ExpectedInfectedByNode=EI,InformationProbabilityByNode=EH,ManagingProbabilityByNode=PM,Escape=if(length(OutsideNodes))list(OutsideNodes=OutsideNodes,ProbabilityByTimestep=escape,ProbabilityByHorizon=tail(escape,1))else NULL,StateCount=nst,Diagnostics=c("Exact small-network joint N/I/information operator; intended as validation truth for small capacities.","Escape is first passage into any OutsideNodes."));if(ReturnOperator)out$Operator<-T;class(out)<-c("INApestMetaSISExactNetworkInformation","list");out
}

###############################################################################
### Exact fixed-host SIR / SEIR one-node compartment models
###############################################################################

.ina_meta_compartment_states <- function(Model,N) {
  Model<-match.arg(Model,c("SIR","SEIR"));N<-as.integer(N);out<-list();k<-0L
  if(Model=="SIR"){
    for(I in 0:N)for(R in 0:(N-I)){k<-k+1L;out[[k]]<-c(S=N-I-R,I=I,R=R)}
  }else{
    for(E in 0:N)for(I in 0:(N-E))for(R in 0:(N-E-I)){k<-k+1L;out[[k]]<-c(S=N-E-I-R,E=E,I=I,R=R)}
  }
  z<-as.data.frame(do.call(rbind,out));for(nm in names(z))z[[nm]]<-as.integer(z[[nm]]);z
}

.ina_meta_compartment_operator_one_node <- function(Model,N,Beta,RecoveryProb,ProgressionProb=1,ImmunityLossProb=0,IntroductionProb=0,IntroductionNumber=1,Transmission=c("frequency","density"),DensityScale=1,ContactWeight=1,timestep=1L,Ntimesteps=1L) {
  Model<-match.arg(Model,c("SIR","SEIR"));Transmission<-match.arg(Transmission);N<-as.integer(N)
  beta<-.ina_meta_sis_step_value(Beta,timestep,Ntimesteps,"Beta");rec<-.ina_meta_sis_step_value(RecoveryProb,timestep,Ntimesteps,"RecoveryProb");prog<-.ina_meta_sis_step_value(ProgressionProb,timestep,Ntimesteps,"ProgressionProb");wan<-.ina_meta_sis_step_value(ImmunityLossProb,timestep,Ntimesteps,"ImmunityLossProb");ip<-.ina_meta_sis_step_value(IntroductionProb,timestep,Ntimesteps,"IntroductionProb");inum<-as.integer(.ina_meta_sis_step_value(IntroductionNumber,timestep,Ntimesteps,"IntroductionNumber",integer=TRUE));ds<-.ina_meta_sis_step_value(DensityScale,timestep,Ntimesteps,"DensityScale",positive=TRUE);cw<-.ina_meta_sis_step_value(ContactWeight,timestep,Ntimesteps,"ContactWeight")
  if(beta<0||cw<0||rec<0||rec>1||prog<0||prog>1||wan<0||wan>1||ip<0||ip>1)stop("Invalid SIR/SEIR parameters")
  st<-.ina_meta_compartment_states(Model,N);key<-apply(st,1,paste,collapse=":");idx<-setNames(seq_len(nrow(st)),key);T<-matrix(0,nrow(st),nrow(st),dimnames=list(key,key))
  for(row in seq_len(nrow(st))){S0<-st$S[row];I0<-st$I[row];R0<-st$R[row];E0<-if(Model=="SEIR")st$E[row] else 0L
    foi<-if(Transmission=="frequency")if(N>0&&cw>0)beta*I0/N else 0 else beta*I0*cw/ds;pinf<-.ina_meta_sis_clip01(-expm1(-max(0,foi)));px<-dbinom(0:S0,S0,pinf);pstay<-dbinom(0:I0,I0,1-rec);plose<-dbinom(0:R0,R0,wan);pprog<-if(Model=="SEIR")dbinom(0:E0,E0,prog) else 1
    for(x in 0:S0){mx<-min(inum,S0-x);iv<-if(ip<=0||mx<=0)0L else c(0L,mx);ipp<-if(ip<=0||mx<=0)1 else c(1-ip,ip)
      for(stay in 0:I0){recover<-I0-stay
        for(lose in 0:R0){
          if(Model=="SIR"){
            for(kk in seq_along(iv)){S1<-S0-x-iv[kk]+lose;I1<-stay+x+iv[kk];R1<-R0-lose+recover;ky<-paste(S1,I1,R1,sep=":");T[row,idx[[ky]]]<-T[row,idx[[ky]]]+px[x+1L]*pstay[stay+1L]*plose[lose+1L]*ipp[kk]}
          }else{
            for(pg in 0:E0)for(kk in seq_along(iv)){S1<-S0-x-iv[kk]+lose;E1<-E0-pg+x+iv[kk];I1<-stay+pg;R1<-R0-lose+recover;ky<-paste(S1,E1,I1,R1,sep=":");T[row,idx[[ky]]]<-T[row,idx[[ky]]]+px[x+1L]*pstay[stay+1L]*plose[lose+1L]*pprog[pg+1L]*ipp[kk]}
          }
        }
      }
    }
  }
  if(max(abs(rowSums(T)-1))>1e-11)stop(Model," operator rows do not sum to 1")
  attr(T,"States")<-st;T
}

INApestMetaCompartmentExactOneNode <- function(Model=c("SIR","SEIR"),Ntimesteps,HostPopulation,InitialInfected=0,InitialExposed=0,InitialRecovered=0,Beta,RecoveryProb,ProgressionProb=1,ImmunityLossProb=0,IntroductionProb=0,IntroductionNumber=1,Transmission=c("frequency","density"),DensityScale=1,ContactWeight=1,DetectionProb=0,ReturnOperators=FALSE) {
  Model<-match.arg(Model);Transmission<-match.arg(Transmission);N<-as.integer(HostPopulation);I0<-as.integer(InitialInfected);E0<-as.integer(InitialExposed);R0<-as.integer(InitialRecovered)
  if(Model=="SIR"&&E0!=0)stop("InitialExposed must be zero for SIR");if(any(c(I0,E0,R0)<0)||I0+E0+R0>N)stop("Initial compartment counts exceed HostPopulation")
  st<-.ina_meta_compartment_states(Model,N);key<-apply(st,1,paste,collapse=":");initkey<-if(Model=="SIR")paste(N-I0-R0,I0,R0,sep=":") else paste(N-E0-I0-R0,E0,I0,R0,sep=":");dist<-numeric(nrow(st));dist[match(initkey,key)]<-1
  traj<-matrix(0,Ntimesteps+1,nrow(st));traj[1,]<-dist;ops<-vector("list",Ntimesteps)
  cols<-names(st);means<-matrix(0,Ntimesteps+1,length(cols),dimnames=list(NULL,cols));means[1,]<-as.numeric(dist%*%as.matrix(st));active<-detect<-numeric(Ntimesteps+1)
  activefun<-function(d)if(Model=="SEIR")sum(d[(st$E+st$I)>0])else sum(d[st$I>0]);active[1]<-activefun(dist)
  detfun<-function(d,tt){pd<-.ina_meta_sis_step_value(DetectionProb,max(1L,tt),max(1L,Ntimesteps),"DetectionProb");sum(d*(1-(1-pd)^st$I))};detect[1]<-detfun(dist,1L)
  if(Ntimesteps>0)for(tt in seq_len(Ntimesteps)){T<-.ina_meta_compartment_operator_one_node(Model,N,Beta,RecoveryProb,ProgressionProb,ImmunityLossProb,IntroductionProb,IntroductionNumber,Transmission,DensityScale,ContactWeight,tt,Ntimesteps);ops[[tt]]<-T;dist<-as.numeric(dist%*%T);traj[tt+1,]<-dist;means[tt+1,]<-as.numeric(dist%*%as.matrix(st));active[tt+1]<-activefun(dist);detect[tt+1]<-detfun(dist,tt)}
  b<-.ina_meta_sis_step_value(Beta,1,max(1L,Ntimesteps),"Beta");r<-.ina_meta_sis_step_value(RecoveryProb,1,max(1L,Ntimesteps),"RecoveryProb");beff<-if(Transmission=="frequency")if(N>0&&ContactWeight>0)b else 0 else b*N*ContactWeight/.ina_meta_sis_step_value(DensityScale,1,max(1L,Ntimesteps),"DensityScale",positive=TRUE)
  if(Model=="SIR"){G<-matrix((1-r)+beff,1,1,dimnames=list("I","I"));lambda<-G[1,1]}else{p<-.ina_meta_sis_step_value(ProgressionProb,1,max(1L,Ntimesteps),"ProgressionProb");G<-matrix(c(1-p,p,beff,1-r),2,2,byrow=FALSE,dimnames=list(c("E","I"),c("E","I")));lambda<-max(Mod(eigen(G,only.values=TRUE)$values))}
  tr<-data.frame(timestep=0:Ntimesteps,means,ActivePathogenProbability=active,PathogenFreedomProbability=1-active,DetectionProbability=detect,check.names=FALSE)
  if(Model=="SIR"&&all(as.numeric(ImmunityLossProb)==0))tr$ExpectedEverInfected<-N-tr$S
  out<-list(Model="INApestMeta",PathogenModel=Model,HostAssumption="fixed abundance",Exact=TRUE,Growth=list(ActiveInfectionLinearOperator=G,IntrinsicRarePathogenMultiplier=lambda),Trajectory=tr,StateTable=st,StateDistribution=traj,Diagnostics=c(paste0("Exact fixed-host one-node ",Model," finite-state solution."),"New infections do not progress/recover in the same pathogen update; each biological state transition is at most one compartment step per timestep.",if(Model=="SEIR")"ProgressionProb=1 still represents a one-timestep latent state for newly infected hosts; it does not algebraically collapse SEIR to SIR in this discrete engine." else "ImmunityLossProb acts on hosts already in R at the start of the pathogen step; newly recovered hosts cannot lose immunity in that same step."));if(ReturnOperators)out$Operators<-ops;class(out)<-c("INApestMetaCompartmentExactOneNode","list");out
}
###############################################################################
### Exact one-node SIR / SEIR with host turnover (event-contract validation)
###############################################################################

.ina_meta_compartment_states_uptoK <- function(Model, K) {
  Model <- match.arg(Model, c("SIR", "SEIR"))
  K <- as.integer(K)
  out <- list(); k <- 0L
  if (Model == "SIR") {
    for (S in 0:K) for (I in 0:(K-S)) for (R in 0:(K-S-I)) {
      k <- k + 1L; out[[k]] <- c(S=S, I=I, R=R)
    }
  } else {
    for (S in 0:K) for (E in 0:(K-S)) for (I in 0:(K-S-E)) for (R in 0:(K-S-E-I)) {
      k <- k + 1L; out[[k]] <- c(S=S, E=E, I=I, R=R)
    }
  }
  z <- as.data.frame(do.call(rbind, out))
  for (nm in names(z)) z[[nm]] <- as.integer(z[[nm]])
  z$N <- rowSums(z)
  z
}

.ina_meta_multinom_allocations <- function(size, probs, names_out) {
  size <- as.integer(size)
  probs <- as.numeric(probs)
  names(probs) <- names_out
  if (size <= 0L) {
    z <- as.data.frame(as.list(setNames(rep(0L, length(names_out)), names_out)))
    z$prob <- 1
    return(z)
  }
  if (length(probs) != length(names_out) || any(probs < 0) || abs(sum(probs)-1) > 1e-10)
    stop("Invalid multinomial allocation probabilities")
  out <- list(); k <- 0L
  recurse <- function(pos, remaining, counts) {
    if (pos == length(names_out)) {
      counts[pos] <- remaining
      lp <- lgamma(size+1) - sum(lgamma(counts+1)) + sum(ifelse(counts>0, counts*log(probs), 0))
      if (any(counts > 0 & probs == 0)) p <- 0 else p <- exp(lp)
      k <<- k + 1L
      out[[k]] <<- c(setNames(counts, names_out), prob=p)
      return(invisible(NULL))
    }
    for (v in 0:remaining) {
      counts[pos] <- v
      recurse(pos+1L, remaining-v, counts)
    }
  }
  recurse(1L, size, integer(length(names_out)))
  z <- as.data.frame(do.call(rbind, out))
  for (nm in names_out) z[[nm]] <- as.integer(z[[nm]])
  z$prob <- as.numeric(z$prob)
  z[z$prob > 0, , drop=FALSE]
}

.ina_meta_compartment_host_outcomes <- function(state, Model, K, Survival,
                                                RecruitToCapacity=FALSE,
                                                ExternalInvasionProb=0,
                                                ExternalInvasionNumber=1,
                                                ExternalPathogenStateProb=NULL) {
  Model <- match.arg(Model, c("SIR", "SEIR"))
  states <- if (Model == "SIR") c("S","I","R") else c("S","E","I","R")
  x <- as.integer(state[states]); names(x) <- states
  surv <- as.numeric(Survival)
  if (!is.finite(surv) || surv < 0 || surv > 1) stop("Survival must be in [0,1]")
  ep <- as.numeric(ExternalInvasionProb)
  if (!is.finite(ep) || ep < 0 || ep > 1) stop("ExternalInvasionProb must be in [0,1]")
  en <- as.integer(ExternalInvasionNumber)
  if (en < 0) stop("ExternalInvasionNumber must be non-negative")
  if (is.null(ExternalPathogenStateProb)) {
    extp <- setNames(rep(0, length(states)), states); extp["S"] <- 1
  } else {
    extp <- setNames(rep(0, length(states)), states)
    if (is.null(names(ExternalPathogenStateProb)) || any(!names(ExternalPathogenStateProb) %in% states))
      stop("ExternalPathogenStateProb names do not match model states")
    extp[names(ExternalPathogenStateProb)] <- as.numeric(ExternalPathogenStateProb)
    if (any(extp < 0) || abs(sum(extp)-1) > 1e-10) stop("ExternalPathogenStateProb must sum to 1")
  }
  # Independent class survival is distributionally equivalent to the parent
  # total Binomial host-survival draw followed by unbiased hypergeometric
  # pathogen-state thinning.
  class_lists <- lapply(states, function(nm) data.frame(v=0:x[nm], p=dbinom(0:x[nm], x[nm], surv)))
  g <- as.matrix(expand.grid(lapply(class_lists, function(z) seq_len(nrow(z))), KEEP.OUT.ATTRS=FALSE))
  out <- list(); kk <- 0L
  for (rr in seq_len(nrow(g))) {
    y <- setNames(integer(length(states)), states); p0 <- 1
    for (j in seq_along(states)) {
      z <- class_lists[[j]][g[rr,j],]
      y[j] <- as.integer(z$v); p0 <- p0 * z$p
    }
    if (RecruitToCapacity) y["S"] <- y["S"] + max(0L, as.integer(K-sum(y)))
    ext_events <- if (ep <= 0 || en <= 0) data.frame(event=0L,prob=1) else data.frame(event=c(0L,1L),prob=c(1-ep,ep))
    for (ee in seq_len(nrow(ext_events))) {
      accept <- if (ext_events$event[ee] == 1L) min(en, max(0L, K-sum(y))) else 0L
      alloc <- .ina_meta_multinom_allocations(accept, extp, states)
      for (aa in seq_len(nrow(alloc))) {
        yy <- y
        for (nm in states) yy[nm] <- yy[nm] + as.integer(alloc[aa,nm])
        kk <- kk + 1L
        out[[kk]] <- c(yy, prob=unname(p0 * ext_events$prob[ee] * alloc$prob[aa]))
      }
    }
  }
  z <- as.data.frame(do.call(rbind, out))
  for (nm in states) z[[nm]] <- as.integer(z[[nm]])
  z$prob <- as.numeric(z$prob)
  # aggregate duplicate compartment outcomes
  key <- apply(z[,states,drop=FALSE],1,paste,collapse=":")
  ps <- tapply(z$prob,key,sum)
  first <- match(names(ps),key)
  ans <- z[first,c(states),drop=FALSE]; ans$prob <- as.numeric(ps)
  rownames(ans) <- NULL; ans
}

.ina_meta_compartment_disease_outcomes <- function(state, Model, Beta, RecoveryProb,
                                                   ProgressionProb=1, ImmunityLossProb=0,
                                                   PathogenMortalityProb=0,
                                                   IntroductionProb=0,
                                                   IntroductionNumber=1,
                                                   Transmission=c("frequency","density"),
                                                   DensityScale=1,
                                                   ContactWeight=1) {
  Model <- match.arg(Model,c("SIR","SEIR")); Transmission <- match.arg(Transmission)
  S0 <- as.integer(state["S"]); I0 <- as.integer(state["I"])
  R0 <- as.integer(state["R"]); E0 <- if (Model=="SEIR") as.integer(state["E"]) else 0L
  N0 <- S0+E0+I0+R0
  beta <- as.numeric(Beta); rec <- as.numeric(RecoveryProb); prog <- as.numeric(ProgressionProb)
  wan <- as.numeric(ImmunityLossProb); mort <- as.numeric(PathogenMortalityProb)
  ip <- as.numeric(IntroductionProb); inum <- as.integer(IntroductionNumber)
  if (rec < 0 || mort < 0 || rec+mort > 1+1e-12) stop("RecoveryProb + PathogenMortalityProb must be <= 1")
  foi <- if (Transmission=="frequency") {
    if (N0 > 0 && ContactWeight > 0) beta*I0/N0 else 0
  } else beta*I0*ContactWeight/as.numeric(DensityScale)
  pinf <- .ina_meta_sis_clip01(-expm1(-max(0,foi)))
  out <- list(); kk <- 0L
  for (x in 0:S0) {
    px <- dbinom(x,S0,pinf); remS <- S0-x; mx <- min(inum,remS)
    intro_vals <- if (ip > 0 && mx > 0) c(0L,mx) else 0L
    intro_probs <- if (ip > 0 && mx > 0) c(1-ip,ip) else 1
    for (ii in seq_along(intro_vals)) {
      intro <- intro_vals[ii]; pi <- intro_probs[ii]
      for (prog_n in if (Model=="SEIR") 0:E0 else 0L) {
        pp <- if (Model=="SEIR") dbinom(prog_n,E0,prog) else 1
        for (lose in 0:R0) {
          pl <- dbinom(lose,R0,wan)
          # Multinomial outcomes of pre-existing infectious hosts: stay/recover/die.
          for (stay in 0:I0) for (recover in 0:(I0-stay)) {
            die <- I0-stay-recover
            pr <- dmultinom(c(stay,recover,die),prob=c(1-rec-mort,rec,mort))
            if (pr == 0) next
            if (Model=="SIR") {
              S1 <- S0-x-intro+lose
              I1 <- stay+x+intro
              R1 <- R0-lose+recover
              kk <- kk+1L; out[[kk]] <- c(S=S1,I=I1,R=R1,prob=unname(px*pi*pp*pl*pr))
            } else {
              S1 <- S0-x-intro+lose
              E1 <- E0-prog_n+x+intro
              I1 <- stay+prog_n
              R1 <- R0-lose+recover
              kk <- kk+1L; out[[kk]] <- c(S=S1,E=E1,I=I1,R=R1,prob=unname(px*pi*pp*pl*pr))
            }
          }
        }
      }
    }
  }
  states <- if (Model=="SIR") c("S","I","R") else c("S","E","I","R")
  z <- as.data.frame(do.call(rbind,out)); for(nm in states) z[[nm]] <- as.integer(z[[nm]])
  z$prob <- as.numeric(z$prob); z <- z[z$prob>0,,drop=FALSE]
  key <- apply(z[,states,drop=FALSE],1,paste,collapse=":"); ps <- tapply(z$prob,key,sum); first <- match(names(ps),key)
  ans <- z[first,states,drop=FALSE]; ans$prob <- as.numeric(ps); rownames(ans)<-NULL; ans
}

.ina_meta_compartment_turnover_operator <- function(Model, K, Survival=1,
                                                    RecruitToCapacity=FALSE,
                                                    ExternalInvasionProb=0,
                                                    ExternalInvasionNumber=1,
                                                    ExternalPathogenStateProb=NULL,
                                                    Beta=0, RecoveryProb=0,
                                                    ProgressionProb=1,
                                                    ImmunityLossProb=0,
                                                    PathogenMortalityProb=0,
                                                    IntroductionProb=0,
                                                    IntroductionNumber=1,
                                                    Transmission=c("frequency","density"),
                                                    DensityScale=1,
                                                    ContactWeight=1,
                                                    timestep=1L,Ntimesteps=1L) {
  Model <- match.arg(Model,c("SIR","SEIR")); Transmission<-match.arg(Transmission)
  st <- .ina_meta_compartment_states_uptoK(Model,K)
  states <- if(Model=="SIR")c("S","I","R")else c("S","E","I","R")
  key <- apply(st[,states,drop=FALSE],1,paste,collapse=":"); idx <- setNames(seq_len(nrow(st)),key)
  T <- matrix(0,nrow(st),nrow(st),dimnames=list(key,key))
  vals <- list(
    Survival=.ina_meta_sis_step_value(Survival,timestep,Ntimesteps,"Survival"),
    ExternalInvasionProb=.ina_meta_sis_step_value(ExternalInvasionProb,timestep,Ntimesteps,"ExternalInvasionProb"),
    Beta=.ina_meta_sis_step_value(Beta,timestep,Ntimesteps,"Beta"),
    RecoveryProb=.ina_meta_sis_step_value(RecoveryProb,timestep,Ntimesteps,"RecoveryProb"),
    ProgressionProb=.ina_meta_sis_step_value(ProgressionProb,timestep,Ntimesteps,"ProgressionProb"),
    ImmunityLossProb=.ina_meta_sis_step_value(ImmunityLossProb,timestep,Ntimesteps,"ImmunityLossProb"),
    PathogenMortalityProb=.ina_meta_sis_step_value(PathogenMortalityProb,timestep,Ntimesteps,"PathogenMortalityProb"),
    IntroductionProb=.ina_meta_sis_step_value(IntroductionProb,timestep,Ntimesteps,"IntroductionProb"),
    DensityScale=.ina_meta_sis_step_value(DensityScale,timestep,Ntimesteps,"DensityScale",positive=TRUE)
  )
  for (row in seq_len(nrow(st))) {
    hs <- .ina_meta_compartment_host_outcomes(st[row,],Model,K,vals$Survival,RecruitToCapacity,
                                             vals$ExternalInvasionProb,ExternalInvasionNumber,
                                             ExternalPathogenStateProb)
    for (hh in seq_len(nrow(hs))) {
      ds <- .ina_meta_compartment_disease_outcomes(hs[hh,],Model,vals$Beta,vals$RecoveryProb,
                                                   vals$ProgressionProb,vals$ImmunityLossProb,
                                                   vals$PathogenMortalityProb,vals$IntroductionProb,
                                                   IntroductionNumber,Transmission,vals$DensityScale,
                                                   ContactWeight)
      for (dd in seq_len(nrow(ds))) {
        ky <- paste(as.integer(ds[dd,states]),collapse=":")
        T[row,idx[[ky]]] <- T[row,idx[[ky]]] + hs$prob[hh]*ds$prob[dd]
      }
    }
  }
  if(max(abs(rowSums(T)-1))>1e-10)stop(Model," turnover operator row error ",max(abs(rowSums(T)-1)))
  attr(T,"States")<-st;T
}

INApestMetaCompartmentExactTurnoverOneNode <- function(Model=c("SIR","SEIR"),Ntimesteps,K,
                                                       InitialHostPopulation=K,
                                                       InitialInfected=0,InitialExposed=0,InitialRecovered=0,
                                                       Survival=1,RecruitToCapacity=FALSE,
                                                       ExternalInvasionProb=0,ExternalInvasionNumber=1,
                                                       ExternalPathogenStateProb=NULL,
                                                       Beta=0,RecoveryProb=0,ProgressionProb=1,
                                                       ImmunityLossProb=0,PathogenMortalityProb=0,
                                                       IntroductionProb=0,IntroductionNumber=1,
                                                       Transmission=c("frequency","density"),DensityScale=1,
                                                       ContactWeight=1,ReturnOperators=FALSE) {
  Model<-match.arg(Model);Transmission<-match.arg(Transmission);K<-as.integer(K);N0<-as.integer(InitialHostPopulation)
  I0<-as.integer(InitialInfected);E0<-as.integer(InitialExposed);R0<-as.integer(InitialRecovered)
  if(Model=="SIR"&&E0!=0)stop("InitialExposed must be zero for SIR")
  if(N0<0||N0>K||I0+E0+R0>N0)stop("Invalid initial host/pathogen counts")
  st<-.ina_meta_compartment_states_uptoK(Model,K); states<-if(Model=="SIR")c("S","I","R")else c("S","E","I","R")
  init <- if(Model=="SIR")c(S=N0-I0-R0,I=I0,R=R0)else c(S=N0-E0-I0-R0,E=E0,I=I0,R=R0)
  key<-apply(st[,states,drop=FALSE],1,paste,collapse=":");dist<-numeric(nrow(st));dist[match(paste(init,collapse=":"),key)]<-1
  ops<-vector("list",Ntimesteps);traj<-matrix(0,Ntimesteps+1,nrow(st));traj[1,]<-dist
  means<-matrix(0,Ntimesteps+1,length(states)+1,dimnames=list(NULL,c(states,"N")))
  means[1,]<-as.numeric(dist%*%as.matrix(st[,c(states,"N")]))
  active<-numeric(Ntimesteps+1);active[1]<-if(Model=="SEIR")sum(dist[(st$E+st$I)>0])else sum(dist[st$I>0])
  for(tt in seq_len(Ntimesteps)){
    T<-.ina_meta_compartment_turnover_operator(Model,K,Survival,RecruitToCapacity,ExternalInvasionProb,
      ExternalInvasionNumber,ExternalPathogenStateProb,Beta,RecoveryProb,ProgressionProb,ImmunityLossProb,
      PathogenMortalityProb,IntroductionProb,IntroductionNumber,Transmission,DensityScale,ContactWeight,tt,Ntimesteps)
    ops[[tt]]<-T;dist<-as.numeric(dist%*%T);traj[tt+1,]<-dist
    means[tt+1,]<-as.numeric(dist%*%as.matrix(st[,c(states,"N")]))
    active[tt+1]<-if(Model=="SEIR")sum(dist[(st$E+st$I)>0])else sum(dist[st$I>0])
  }
  tr<-data.frame(timestep=0:Ntimesteps,means,ActivePathogenProbability=active,PathogenFreedomProbability=1-active,check.names=FALSE)
  out<-list(Model="INApestMeta",PathogenModel=Model,HostAssumption="explicit host turnover",Exact=TRUE,
            Trajectory=tr,StateTable=st,StateDistribution=traj,
            Diagnostics=c(paste0("Exact one-node ",Model," joint host-turnover/pathogen finite-state solution."),
              "Natural host mortality thins S/E/I/R without replacement; optional replacement recruitment enters S.",
              "External host immigrants can be assigned explicit pathogen states; resident IntroductionProb leaves N unchanged."))
  if(ReturnOperators)out$Operators<-ops;class(out)<-c("INApestMetaCompartmentExactTurnoverOneNode","list");out
}

###############################################################################
### Rare-pathogen spatial growth operators for SIS / SIR / SEIR Meta
###############################################################################

INApestMetaPathogenGrowthOperator <- function(Model=c("SIS","SIR","SEIR"),
                                              HostPopulation,
                                              Beta,
                                              RecoveryProb,
                                              ProgressionProb=1,
                                              PathogenMortalityProb=0,
                                              ContactMatrix=NULL,
                                              Transmission=c("frequency","density"),
                                              DensityScale=1) {
  Model<-match.arg(Model);Transmission<-match.arg(Transmission)
  N<-as.numeric(HostPopulation);n<-length(N)
  if(!n||any(!is.finite(N))||any(N<0))stop("HostPopulation must contain finite non-negative values")
  recyc<-function(x,name){x<-as.numeric(x);if(length(x)==1L)x<-rep(x,n);if(length(x)!=n||any(!is.finite(x)))stop(name," must be scalar or length nodes");x}
  beta<-recyc(Beta,"Beta");rec<-recyc(RecoveryProb,"RecoveryProb");mort<-recyc(PathogenMortalityProb,"PathogenMortalityProb");prog<-recyc(ProgressionProb,"ProgressionProb");ds<-recyc(DensityScale,"DensityScale")
  if(any(beta<0)||any(rec<0|rec>1)||any(mort<0|mort>1)||any(rec+mort>1+1e-12)||any(prog<0|prog>1)||any(ds<=0))stop("Invalid pathogen growth parameters")
  C<-if(is.null(ContactMatrix))diag(n)else as.matrix(ContactMatrix)
  if(!all(dim(C)==c(n,n))||any(!is.finite(C))||any(C<0))stop("ContactMatrix must be a finite non-negative source x target nodes matrix")
  B<-matrix(0,n,n,dimnames=list(paste0("node",seq_len(n)),paste0("node",seq_len(n))))
  if(Transmission=="frequency"){
    denom<-as.numeric(crossprod(N,C))
    for(j in seq_len(n))if(N[j]>0&&denom[j]>0&&beta[j]>0)B[j,]<-N[j]*beta[j]*C[,j]/denom[j]
  }else{
    for(j in seq_len(n))if(N[j]>0&&beta[j]>0)B[j,]<-N[j]*beta[j]*C[,j]/ds[j]
  }
  stayI<-diag(1-rec-mort,n,n)
  if(Model %in% c("SIS","SIR")){
    G<-stayI+B
    dimnames(G)<-list(paste0("I_node",seq_len(n)),paste0("I_node",seq_len(n)))
    note<-if(Model=="SIR")"At the pathogen-free state SIR and SIS have the same infectious-state linearisation; recovered hosts are absent to first order." else "SIS infectious-state linearisation."
  }else{
    PE<-diag(prog,n,n);stayE<-diag(1-prog,n,n)
    G<-rbind(cbind(stayE,B),cbind(PE,stayI))
    dimnames(G)<-list(c(paste0("E_node",seq_len(n)),paste0("I_node",seq_len(n))),c(paste0("E_node",seq_len(n)),paste0("I_node",seq_len(n))))
    note<-"SEIR linearisation preserves the discrete latent-state delay: new infections enter E, and only E present at the start of a pathogen step can progress to I."
  }
  ev<-eigen(G,only.values=TRUE)$values;lambda<-max(Mod(ev))
  list(Model=Model,Operator=G,TransmissionBlock=B,IntrinsicRarePathogenMultiplier=lambda,
       Classification=if(lambda>1+1e-12)"growing"else if(lambda<1-1e-12)"declining"else "threshold",
       Orientation="rows = recipient active states; columns = source active states",
       Diagnostics=c("Linearisation is about the pathogen-free fixed-host state; external introductions and susceptible depletion are excluded.",note))
}

###############################################################################
### Exact small-network fixed-host SIR / SEIR + first-passage escape
###############################################################################

.ina_meta_compartment_local_outcomes_pinf <- function(state,Model,pinf,RecoveryProb,
                                                      ProgressionProb=1,ImmunityLossProb=0,
                                                      IntroductionProb=0,IntroductionNumber=1) {
  Model<-match.arg(Model,c("SIR","SEIR"));S0<-as.integer(state["S"]);I0<-as.integer(state["I"]);R0<-as.integer(state["R"]);E0<-if(Model=="SEIR")as.integer(state["E"])else 0L
  rec<-as.numeric(RecoveryProb);prog<-as.numeric(ProgressionProb);wan<-as.numeric(ImmunityLossProb);ip<-as.numeric(IntroductionProb);inum<-as.integer(IntroductionNumber)
  out<-list();kk<-0L
  for(x in 0:S0){px<-dbinom(x,S0,pinf);mx<-min(inum,S0-x);iv<-if(ip>0&&mx>0)c(0L,mx)else 0L;ipp<-if(ip>0&&mx>0)c(1-ip,ip)else 1
    for(stay in 0:I0){pr<-dbinom(stay,I0,1-rec);recover<-I0-stay
      for(lose in 0:R0){pl<-dbinom(lose,R0,wan)
        if(Model=="SIR"){
          for(q in seq_along(iv)){kk<-kk+1L;out[[kk]]<-c(S=S0-x-iv[q]+lose,I=stay+x+iv[q],R=R0-lose+recover,prob=unname(px*pr*pl*ipp[q]))}
        }else{
          for(pg in 0:E0){pp<-dbinom(pg,E0,prog);for(q in seq_along(iv)){kk<-kk+1L;out[[kk]]<-c(S=S0-x-iv[q]+lose,E=E0-pg+x+iv[q],I=stay+pg,R=R0-lose+recover,prob=unname(px*pr*pl*pp*ipp[q]))}}
        }
      }
    }
  }
  states<-if(Model=="SIR")c("S","I","R")else c("S","E","I","R");z<-as.data.frame(do.call(rbind,out));for(nm in states)z[[nm]]<-as.integer(z[[nm]]);z$prob<-as.numeric(z$prob);z<-z[z$prob>0,,drop=FALSE]
  key<-apply(z[,states,drop=FALSE],1,paste,collapse=":");ps<-tapply(z$prob,key,sum);first<-match(names(ps),key);ans<-z[first,states,drop=FALSE];ans$prob<-as.numeric(ps);rownames(ans)<-NULL;ans
}

INApestMetaCompartmentExactNetworkFixedN <- function(Model=c("SIR","SEIR"),Ntimesteps,
                                                     HostPopulation,
                                                     InitialInfected,
                                                     InitialExposed=0,
                                                     InitialRecovered=0,
                                                     Beta,
                                                     RecoveryProb,
                                                     ProgressionProb=1,
                                                     ImmunityLossProb=0,
                                                     IntroductionProb=0,
                                                     IntroductionNumber=1,
                                                     ContactMatrix=NULL,
                                                     Transmission=c("frequency","density"),
                                                     DensityScale=1,
                                                     OutsideNodes=integer(0),
                                                     MaxStates=50000L,
                                                     ReturnOperator=FALSE) {
  Model<-match.arg(Model);Transmission<-match.arg(Transmission);N<-as.integer(HostPopulation);n<-length(N)
  recycle_n<-function(x,name){x<-as.numeric(x);if(length(x)==1L)x<-rep(x,n);if(length(x)!=n)stop(name," must be scalar or length nodes");x}
  I0<-as.integer(recycle_n(InitialInfected,"InitialInfected"));E0<-as.integer(recycle_n(InitialExposed,"InitialExposed"));R0<-as.integer(recycle_n(InitialRecovered,"InitialRecovered"));if(Model=="SIR"&&any(E0!=0))stop("InitialExposed must be zero for SIR")
  if(any(I0+E0+R0>N))stop("Initial pathogen counts exceed HostPopulation")
  beta<-recycle_n(Beta,"Beta");rec<-recycle_n(RecoveryProb,"RecoveryProb");prog<-recycle_n(ProgressionProb,"ProgressionProb");wan<-recycle_n(ImmunityLossProb,"ImmunityLossProb");ip<-recycle_n(IntroductionProb,"IntroductionProb");ds<-recycle_n(DensityScale,"DensityScale")
  inum<-as.integer(recycle_n(IntroductionNumber,"IntroductionNumber"));C<-if(is.null(ContactMatrix))diag(n)else as.matrix(ContactMatrix);if(!all(dim(C)==c(n,n)))stop("ContactMatrix dimensions invalid")
  local<-lapply(N,function(nn).ina_meta_compartment_states(Model,nn));dims<-vapply(local,nrow,integer(1));nst<-prod(dims);if(nst>MaxStates)stop("Exact network state space has ",nst," states; increase MaxStates only for deliberate small-network validation")
  grid<-as.matrix(expand.grid(lapply(dims,seq_len),KEEP.OUT.ATTRS=FALSE));storage.mode(grid)<-"integer";gidx<-setNames(seq_len(nrow(grid)),apply(grid,1,paste,collapse=":"));states<-if(Model=="SIR")c("S","I","R")else c("S","E","I","R")
  Imat<-Emat<-matrix(0,nst,n);if(Model=="SIR")Emat[]<-0
  for(rr in seq_len(nst))for(j in seq_len(n)){z<-local[[j]][grid[rr,j],];Imat[rr,j]<-z$I;if(Model=="SEIR")Emat[rr,j]<-z$E}
  makeT<-function(tt){T<-matrix(0,nst,nst);for(row in seq_len(nst)){Iv<-Imat[row,];den<-as.numeric(crossprod(N,C));pressure<-as.numeric(crossprod(Iv,C));foi<-if(Transmission=="frequency")beta*ifelse(den>0,pressure/den,0)else beta*pressure/ds;pinf<-.ina_meta_sis_clip01(-expm1(-pmax(0,foi)))
      lists<-vector("list",n);for(j in seq_len(n)){z<-local[[j]][grid[row,j],];lists[[j]]<-.ina_meta_compartment_local_outcomes_pinf(z,Model,pinf[j],rec[j],prog[j],wan[j],ip[j],inum[j])}
      ng<-as.matrix(expand.grid(lapply(lists,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(ng)<-"integer"
      for(bb in seq_len(nrow(ng))){loc<-integer(n);pp<-1;for(j in seq_len(n)){z<-lists[[j]][ng[bb,j],];pp<-pp*z$prob;ky<-paste(as.integer(z[states]),collapse=":");loc[j]<-match(ky,apply(local[[j]][,states,drop=FALSE],1,paste,collapse=":"))};T[row,gidx[[paste(loc,collapse=":")]]]<-T[row,gidx[[paste(loc,collapse=":")]]]+pp}
    };if(max(abs(rowSums(T)-1))>1e-10)stop("Network compartment operator row error");T}
  # Time-varying pathogen parameters can be represented by extending the resolver; this exact
  # validator currently takes node vectors fixed over the validation horizon.
  T<-makeT(1L)
  initloc<-integer(n);for(j in seq_len(n)){z<-if(Model=="SIR")c(S=N[j]-I0[j]-R0[j],I=I0[j],R=R0[j])else c(S=N[j]-E0[j]-I0[j]-R0[j],E=E0[j],I=I0[j],R=R0[j]);initloc[j]<-match(paste(z,collapse=":"),apply(local[[j]][,states,drop=FALSE],1,paste,collapse=":"))}
  dist<-numeric(nst);dist[gidx[[paste(initloc,collapse=":")]]]<-1
  EI<-matrix(0,Ntimesteps+1,n);EE<-matrix(0,Ntimesteps+1,n);EI[1,]<-as.numeric(dist%*%Imat);EE[1,]<-as.numeric(dist%*%Emat);active<-if(Model=="SEIR")rowSums(Emat+Imat)>0 else rowSums(Imat)>0;presence<-numeric(Ntimesteps+1);presence[1]<-sum(dist[active])
  OutsideNodes<-as.integer(OutsideNodes);escstate<-if(length(OutsideNodes))if(Model=="SEIR")apply((Emat+Imat)[,OutsideNodes,drop=FALSE],1,function(z)any(z>0))else apply(Imat[,OutsideNodes,drop=FALSE],1,function(z)any(z>0)) else rep(FALSE,nst);escape<-numeric(Ntimesteps+1);safe<-dist;if(length(OutsideNodes)&&any((I0+E0)[OutsideNodes]>0)){escape[1]<-1;safe[]<-0}else safe[escstate]<-0
  full<-matrix(0,Ntimesteps+1,nst);full[1,]<-dist
  for(tt in seq_len(Ntimesteps)){dist<-as.numeric(dist%*%T);full[tt+1,]<-dist;EI[tt+1,]<-as.numeric(dist%*%Imat);EE[tt+1,]<-as.numeric(dist%*%Emat);presence[tt+1]<-sum(dist[active]);if(length(OutsideNodes)){nx<-as.numeric(safe%*%T);escape[tt+1]<-escape[tt]+sum(nx[escstate]);nx[escstate]<-0;safe<-nx}}
  growth<-INApestMetaPathogenGrowthOperator(Model,N,beta,rec,prog,0,C,Transmission,ds)
  out<-list(Model="INApestMeta",PathogenModel=Model,HostAssumption="fixed abundance network",Exact=TRUE,Growth=growth,ExpectedInfectedByNode=EI,ExpectedExposedByNode=if(Model=="SEIR")EE else NULL,ActivePathogenProbability=presence,PathogenFreedomProbability=1-presence,Escape=if(length(OutsideNodes))list(OutsideNodes=OutsideNodes,ProbabilityByTimestep=escape,ProbabilityByHorizon=tail(escape,1))else NULL,StateCount=nst,StateDistribution=full,Diagnostics=c(paste0("Exact small-network fixed-host ",Model," finite-state solution."),"Transmission is synchronous and uses infectious hosts present at the start of the disease step.",if(Model=="SEIR")"Escape counts first E or I presence outside containment."else"Escape counts first I presence outside containment."));if(ReturnOperator)out$Operator<-T;class(out)<-c("INApestMetaCompartmentExactNetworkFixedN","list");out
}

###############################################################################
### User-facing Meta pathogen analytical router
###############################################################################

.ina_meta_pathogen_initial_vector <- function(x,n,name) {
  if (is.function(x)) stop(name," resolver functions are not supported for initial counts in the exact Meta pathogen router")
  x<-as.numeric(x);if(length(x)==1L)x<-rep(x,n);if(length(x)!=n||any(!is.finite(x))||any(x<0)||any(x!=floor(x)))stop(name," must be a non-negative whole-number scalar or length nodes");as.integer(x)
}

INApestMetaPathogenAnalytical <- function(
    Ntimesteps,
    HostPopulation,
    Pathogen,
    HostAssumption=c("fixed","event_contract_turnover"),
    K=HostPopulation,
    Survival=1,
    RecruitToCapacity=FALSE,
    ExternalHostInvasionProb=0,
    ExternalHostNumber=1,
    ExternalPathogenStateProb=NULL,
    HostDetectionProb=0,
    InitialInfo=0,
    ManageProb=0,
    MortalityProb=0,
    InformationAcquisition=NULL,
    InfoPersistenceSteps=NA,
    InfoRetentionProb=1,
    OutsideNodes=integer(0),
    Exact=TRUE,
    ExactMaxStates=50000L,
    ReturnOperators=FALSE) {
  if(!inherits(Pathogen,"INApestPathogen"))stop("Pathogen must be returned by INApestPathogen()")
  Model<-as.character(Pathogen$Model)[1L];if(!Model%in%c("SIS","SIR","SEIR"))stop("Meta pathogen analytical router supports SIS, SIR and SEIR")
  HostAssumption<-match.arg(HostAssumption);N<-as.integer(HostPopulation);n<-length(N);if(!n||any(N<0)||any(!is.finite(N)))stop("HostPopulation must contain non-negative whole numbers")
  K<-as.integer(K);if(length(K)==1L)K<-rep(K,n);if(length(K)!=n||any(K<N))stop("K must be scalar or length nodes and >= HostPopulation")
  I0<-.ina_meta_pathogen_initial_vector(Pathogen$InitialInfected,n,"InitialInfected")
  E0<-if(Model=="SEIR").ina_meta_pathogen_initial_vector(Pathogen$InitialExposed,n,"InitialExposed") else rep(0L,n)
  R0<-if(Model%in%c("SIR","SEIR")).ina_meta_pathogen_initial_vector(Pathogen$InitialRecovered,n,"InitialRecovered") else rep(0L,n)
  if(any(I0+E0+R0>N))stop("Initial pathogen-state counts exceed HostPopulation")
  acq<-InformationAcquisition
  if(is.null(acq))acq<-if(isTRUE(Pathogen$DetectionTriggersInfo))"both"else"host"
  acq<-match.arg(as.character(acq),c("host","pathogen","both"))
  info_requested<-any(as.numeric(ManageProb)!=0,na.rm=TRUE)||any(as.numeric(InitialInfo)!=0,na.rm=TRUE)||any(as.numeric(HostDetectionProb)!=0,na.rm=TRUE)||any(as.numeric(Pathogen$DetectionProb)!=0,na.rm=TRUE)||any(!is.na(InfoPersistenceSteps))||any(as.numeric(InfoRetentionProb)!=1,na.rm=TRUE)
  if(!isTRUE(Exact)){
    g<-INApestMetaPathogenGrowthOperator(Model,N,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$PathogenMortalityProb,Pathogen$ContactMatrix,Pathogen$Transmission,Pathogen$DensityScale)
    return(structure(list(Model="INApestMeta",PathogenModel=Model,Exact=FALSE,Growth=g,Diagnostics=c("Scalable rare-pathogen result only; Exact=FALSE does not propagate finite-prevalence state probabilities.","Use stochastic INApestMeta as the final arbiter once susceptible depletion, host turnover or shared information becomes important.")),class=c("INApestMetaPathogenAnalytical","list")))
  }
  if(HostAssumption=="event_contract_turnover"){
    if(n!=1L)stop("Exact event_contract_turnover is currently a one-node validation solution")
    if(info_requested)stop("Exact SIR/SEIR/SIS turnover and information are validated separately; combined turnover+information router is not exposed as a generic exact method")
    if(Model=="SIS"){
      return(INApestMetaSISExactTurnoverOneNode(Ntimesteps,K[1],N[1],I0[1],Survival,RecruitToCapacity,ExternalHostInvasionProb,ExternalHostNumber,
        if(is.null(ExternalPathogenStateProb))c(S=1,I=0)else ExternalPathogenStateProb,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$PathogenMortalityProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,Pathogen$Transmission,Pathogen$DensityScale,1,Pathogen$DetectionProb,ReturnOperators))
    }
    return(INApestMetaCompartmentExactTurnoverOneNode(Model,Ntimesteps,K[1],N[1],I0[1],E0[1],R0[1],Survival,RecruitToCapacity,ExternalHostInvasionProb,ExternalHostNumber,ExternalPathogenStateProb,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$PathogenMortalityProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,Pathogen$Transmission,Pathogen$DensityScale,1,ReturnOperators))
  }
  # Fixed-host branches.
  if(Model=="SIS" && info_requested){
    if(any(as.numeric(Pathogen$PathogenMortalityProb)!=0))stop("Fixed-host SIS information solution requires PathogenMortalityProb=0 unless host loss is explicitly modelled")
    if(n==1L){
      return(INApestMetaSISExactInformationOneNode(Ntimesteps,K[1],N[1],I0[1],as.integer(InitialInfo)[1],Survival,ManageProb,MortalityProb,HostDetectionProb,Pathogen$DetectionProb,acq,InfoPersistenceSteps,InfoRetentionProb,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$PathogenMortalityProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,Pathogen$Transmission,Pathogen$DensityScale,1,ReturnOperators))
    }
    if(any(K!=N))stop("Exact fixed-host network information currently requires K == HostPopulation")
    return(INApestMetaSISExactNetworkInformation(Ntimesteps,K,N,I0,InitialInfo,Survival,ManageProb,MortalityProb,HostDetectionProb,Pathogen$DetectionProb,acq,InfoPersistenceSteps,InfoRetentionProb,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$PathogenMortalityProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,Pathogen$ContactMatrix,Pathogen$Transmission,Pathogen$DensityScale,OutsideNodes,ReturnOperators,ExactMaxStates))
  }
  if(info_requested && Model%in%c("SIR","SEIR")){
    if(n!=1L)stop("Exact SIR/SEIR management-information expansion is currently exposed for one-node validation systems; use stochastic Meta for larger managed networks")
    return(INApestMetaCompartmentExactInformationOneNode(Model,Ntimesteps,K[1],N[1],I0[1],E0[1],R0[1],as.integer(InitialInfo)[1],Survival,ManageProb,MortalityProb,HostDetectionProb,Pathogen$DetectionProb,acq,InfoPersistenceSteps,InfoRetentionProb,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$PathogenMortalityProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,Pathogen$Transmission,Pathogen$DensityScale,1,ReturnOperators))
  }
  if(any(as.numeric(Pathogen$PathogenMortalityProb)!=0))stop("Fixed-host exact branch requires PathogenMortalityProb=0; use HostAssumption='event_contract_turnover' for disease mortality")
  if(n==1L){
    if(Model=="SIS")return(INApestMetaSISExactOneNode(Ntimesteps,N[1],I0[1],Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,Pathogen$Transmission,Pathogen$DensityScale,1,Pathogen$DetectionProb,ReturnOperators))
    return(INApestMetaCompartmentExactOneNode(Model,Ntimesteps,N[1],I0[1],E0[1],R0[1],Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,Pathogen$Transmission,Pathogen$DensityScale,1,Pathogen$DetectionProb,ReturnOperators))
  }
  if(Model=="SIS")return(INApestMetaSISExactNetworkFixedN(Ntimesteps,N,I0,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ContactMatrix,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,Pathogen$Transmission,Pathogen$DensityScale,OutsideNodes,ReturnOperators))
  INApestMetaCompartmentExactNetworkFixedN(Model,Ntimesteps,N,I0,E0,R0,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,Pathogen$ContactMatrix,Pathogen$Transmission,Pathogen$DensityScale,OutsideNodes,ExactMaxStates,ReturnOperators)
}

print.INApestMetaPathogenAnalytical <- function(x,...) {
  cat("INApestMeta pathogen analytical result\n")
  cat("  Pathogen model:",x$PathogenModel,"\n")
  cat("  Exact:",isTRUE(x$Exact),"\n")
  if(!is.null(x$Growth$IntrinsicRarePathogenMultiplier))cat("  Rare-pathogen multiplier:",round(x$Growth$IntrinsicRarePathogenMultiplier,6),"\n")
  invisible(x)
}

###############################################################################
### Exact one-node SIR / SEIR + management/information/programmed stop
###############################################################################

.ina_meta_compartment_info_states <- function(Model,Khost,InfoPersistenceSteps=NA){
  Model<-match.arg(Model,c("SIR","SEIR"));b<-.ina_meta_compartment_states_uptoK(Model,Khost);states<-if(Model=="SIR")c("S","I","R")else c("S","E","I","R");rows<-list();k<-0L;finite<-!is.na(InfoPersistenceSteps)
  for(r in seq_len(nrow(b))){base<-as.list(b[r,states,drop=FALSE]);k<-k+1L;rows[[k]]<-c(base,list(H=0L,A=-2L));if(finite){k<-k+1L;rows[[k]]<-c(base,list(H=1L,A=-1L));ages<-0:max(0L,as.integer(InfoPersistenceSteps)-1L);for(a in unique(ages)){k<-k+1L;rows[[k]]<-c(base,list(H=1L,A=as.integer(a)))}}else{k<-k+1L;rows[[k]]<-c(base,list(H=1L,A=-2L))}}
  z<-do.call(rbind,lapply(rows,as.data.frame));rownames(z)<-NULL;for(nm in c(states,"H","A"))z[[nm]]<-as.integer(z[[nm]]);z$N<-rowSums(z[,states,drop=FALSE]);z
}

.ina_meta_compartment_info_prebranches <- function(state_row,Model,Survival,ManageProb,MortalityProb,host_acq,finite,Kinfo,InfoRetentionProb){
  states<-if(Model=="SIR")c("S","I","R")else c("S","E","I","R");H0<-as.integer(state_row[["H"]]);A0<-as.integer(state_row[["A"]]);mprob<-ManageProb*H0;out<-list();kk<-0L
  for(managing in 0:1){pmng<-if(managing==1)mprob else 1-mprob;if(pmng==0)next;mm<-MortalityProb*managing;lists<-lapply(states,function(nm).ina_meta_class_mortality_outcomes(as.integer(state_row[[nm]]),Survival,mm));g<-as.matrix(expand.grid(lapply(lists,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(g)<-"integer"
    for(rr in seq_len(nrow(g))){surv<-setNames(integer(length(states)),states);md<-0L;p<-pmng;for(j in seq_along(states)){z<-lists[[j]][g[rr,j],];surv[j]<-as.integer(z["survive"]);md<-md+as.integer(z["manage_death"]);p<-p*unname(z["prob"])};if(p==0)next;evidence<-host_acq&&md>0
      ib<-list();if(finite){if(H0==0){ib[[1]]<-c(H=0,A=-2,p=1)}else if(evidence){if(Kinfo<=0)ib[[1]]<-c(H=0,A=-2,p=1)else ib[[1]]<-c(H=1,A=0,p=1)}else if(A0<0){ib[[1]]<-c(H=0,A=-2,p=1)}else{ag<-A0+1L;if(ag>=Kinfo)ib[[1]]<-c(H=0,A=-2,p=1)else ib[[1]]<-c(H=1,A=ag,p=1)}}else{if(H0==1&&InfoRetentionProb<1){ib[[1]]<-c(H=1,A=-2,p=InfoRetentionProb);ib[[2]]<-c(H=0,A=-2,p=1-InfoRetentionProb)}else ib[[1]]<-c(H=H0,A=-2,p=1)}
      for(q in ib){kk<-kk+1L;out[[kk]]<-c(surv,H=as.integer(q["H"]),A=as.integer(q["A"]),prob=unname(p*q["p"]))}
    }}
  z<-as.data.frame(do.call(rbind,out));for(nm in c(states,"H","A"))z[[nm]]<-as.integer(z[[nm]]);z$prob<-as.numeric(z$prob);z
}

.ina_meta_compartment_info_operator <- function(Model,Khost,Survival=1,ManageProb=0,MortalityProb=0,HostDetectionProb=0,PathogenDetectionProb=0,InformationAcquisition=c("host","pathogen","both"),InfoPersistenceSteps=NA,InfoRetentionProb=1,Beta=0,RecoveryProb=0,ProgressionProb=1,ImmunityLossProb=0,PathogenMortalityProb=0,IntroductionProb=0,IntroductionNumber=1,Transmission=c("frequency","density"),DensityScale=1,ContactWeight=1){
  Model<-match.arg(Model,c("SIR","SEIR"));InformationAcquisition<-match.arg(InformationAcquisition);Transmission<-match.arg(Transmission);finite<-!is.na(InfoPersistenceSteps);Kinfo<-if(finite)as.integer(InfoPersistenceSteps)else NA_integer_;if(finite&&Kinfo<0)stop("InfoPersistenceSteps must be non-negative");st<-.ina_meta_compartment_info_states(Model,Khost,InfoPersistenceSteps);states<-if(Model=="SIR")c("S","I","R")else c("S","E","I","R");keys<-apply(st[,c(states,"H","A"),drop=FALSE],1,paste,collapse=":");idx<-setNames(seq_len(nrow(st)),keys);T<-matrix(0,nrow(st),nrow(st),dimnames=list(keys,keys));host_acq<-InformationAcquisition%in%c("host","both");path_acq<-InformationAcquisition%in%c("pathogen","both")
  for(row in seq_len(nrow(st))){pre<-.ina_meta_compartment_info_prebranches(st[row,],Model,Survival,ManageProb,MortalityProb,host_acq,finite,Kinfo,InfoRetentionProb);for(pb in seq_len(nrow(pre))){ds<-.ina_meta_compartment_disease_outcomes(pre[pb,],Model,Beta,RecoveryProb,ProgressionProb,ImmunityLossProb,PathogenMortalityProb,IntroductionProb,IntroductionNumber,Transmission,DensityScale,ContactWeight);for(dd in seq_len(nrow(ds))){N1<-sum(as.integer(ds[dd,states]));I1<-as.integer(ds[dd,"I"]);qh<-if(host_acq)1-(1-HostDetectionProb)^N1 else 0;qp<-if(path_acq)1-(1-PathogenDetectionProb)^I1 else 0;qev<-1-(1-qh)*(1-qp);branches<-list();bb<-0L;if(qev>0){bb<-bb+1L;branches[[bb]]<-c(H=1,A=if(finite)0 else -2,p=qev)};if(qev<1){bb<-bb+1L;branches[[bb]]<-c(H=pre$H[pb],A=pre$A[pb],p=1-qev)};for(q in branches){ky<-paste(c(as.integer(ds[dd,states]),as.integer(q["H"]),as.integer(q["A"])),collapse=":");T[row,idx[[ky]]]<-T[row,idx[[ky]]]+pre$prob[pb]*ds$prob[dd]*unname(q["p"])}}}}
  if(max(abs(rowSums(T)-1))>1e-10)stop(Model," information operator row error ",max(abs(rowSums(T)-1)));attr(T,"States")<-st;T
}

INApestMetaCompartmentExactInformationOneNode <- function(Model=c("SIR","SEIR"),Ntimesteps,K,InitialHostPopulation,InitialInfected=0,InitialExposed=0,InitialRecovered=0,InitialInfo=0,Survival=1,ManageProb=0,MortalityProb=0,HostDetectionProb=0,PathogenDetectionProb=0,InformationAcquisition=c("host","pathogen","both"),InfoPersistenceSteps=NA,InfoRetentionProb=1,Beta=0,RecoveryProb=0,ProgressionProb=1,ImmunityLossProb=0,PathogenMortalityProb=0,IntroductionProb=0,IntroductionNumber=1,Transmission=c("frequency","density"),DensityScale=1,ContactWeight=1,ReturnOperator=FALSE){
  Model<-match.arg(Model);InformationAcquisition<-match.arg(InformationAcquisition);Transmission<-match.arg(Transmission);N0<-as.integer(InitialHostPopulation);I0<-as.integer(InitialInfected);E0<-as.integer(InitialExposed);R0<-as.integer(InitialRecovered);if(Model=="SIR"&&E0!=0)stop("InitialExposed must be zero for SIR");if(I0+E0+R0>N0||N0>K)stop("Invalid initial counts");states<-if(Model=="SIR")c("S","I","R")else c("S","E","I","R");base<-if(Model=="SIR")c(S=N0-I0-R0,I=I0,R=R0)else c(S=N0-E0-I0-R0,E=E0,I=I0,R=R0);finite<-!is.na(InfoPersistenceSteps);host_acq<-InformationAcquisition%in%c("host","both");path_acq<-InformationAcquisition%in%c("pathogen","both");T<-.ina_meta_compartment_info_operator(Model,K,Survival,ManageProb,MortalityProb,HostDetectionProb,PathogenDetectionProb,InformationAcquisition,InfoPersistenceSteps,InfoRetentionProb,Beta,RecoveryProb,ProgressionProb,ImmunityLossProb,PathogenMortalityProb,IntroductionProb,IntroductionNumber,Transmission,DensityScale,ContactWeight);st<-attr(T,"States");keys<-rownames(T);qh<-if(host_acq)1-(1-HostDetectionProb)^N0 else 0;qp<-if(path_acq)1-(1-PathogenDetectionProb)^I0 else 0;qev<-1-(1-qh)*(1-qp);dist<-numeric(nrow(st));if(qev>0){ky<-paste(c(base,1,if(finite)0 else -2),collapse=":");dist[idx<-match(ky,keys)]<-qev};if(qev<1){H<-as.integer(InitialInfo);A<-if(H==1&&finite)-1 else -2;ky<-paste(c(base,H,A),collapse=":");ii<-match(ky,keys);dist[ii]<-dist[ii]+1-qev}
  means<-matrix(0,Ntimesteps+1,length(states)+1,dimnames=list(NULL,c(states,"N")));means[1,]<-as.numeric(dist%*%as.matrix(st[,c(states,"N")])) ;PH<-numeric(Ntimesteps+1);PH[1]<-sum(dist*st$H);PM<-numeric(Ntimesteps);active<-numeric(Ntimesteps+1);active[1]<-if(Model=="SEIR")sum(dist[(st$E+st$I)>0])else sum(dist[st$I>0]);full<-matrix(0,Ntimesteps+1,nrow(st));full[1,]<-dist
  for(tt in seq_len(Ntimesteps)){PM[tt]<-ManageProb*sum(dist*st$H);dist<-as.numeric(dist%*%T);full[tt+1,]<-dist;means[tt+1,]<-as.numeric(dist%*%as.matrix(st[,c(states,"N")]));PH[tt+1]<-sum(dist*st$H);active[tt+1]<-if(Model=="SEIR")sum(dist[(st$E+st$I)>0])else sum(dist[st$I>0])}
  out<-list(Model="INApestMeta",PathogenModel=Model,Exact=TRUE,InformationAcquisition=InformationAcquisition,Trajectory=data.frame(timestep=0:Ntimesteps,means,InformationProbability=PH,ActivePathogenProbability=active,check.names=FALSE),ManagingProbability=PM,StateTable=st,StateDistribution=full,Diagnostics=c(paste0("Exact one-node ",Model," host/pathogen/information state for small K."),"HaveInfo is the single management-control state; acquisition selects host, pathogen or both as direct local evidence.",if(finite)"Programmed stopping uses explicit time-since-local-evidence states."else"Information retention is memoryless."));if(ReturnOperator)out$Operator<-T;class(out)<-c("INApestMetaCompartmentExactInformationOneNode","list");out
}

###############################################################################
### Final Meta-pathogen dispatcher extension.
### Preserve Binary pathogen and all pre-existing host-only analytical calls.
###############################################################################
INApestAnalytical_pre_meta_pathogen <- INApestAnalytical
INApestAnalytical <- function(...) {
  args <- list(...)
  Model <- if (!is.null(args$Model)) as.character(args$Model)[1L] else "INApest"
  Pathogen <- args$Pathogen
  if (!identical(Model, "INApestMeta") || is.null(Pathogen))
    return(do.call(INApestAnalytical_pre_meta_pathogen, args))
  if (is.null(args$HostPopulation)) {
    if (is.null(args$InitialState)) stop("Meta pathogen analytical calls require HostPopulation or InitialState")
    args$HostPopulation <- args$InitialState
  }
  keep <- c("Ntimesteps","HostPopulation","Pathogen","HostAssumption","K","Survival",
            "RecruitToCapacity","ExternalHostInvasionProb","ExternalHostNumber",
            "ExternalPathogenStateProb","HostDetectionProb","InitialInfo","ManageProb",
            "MortalityProb","InformationAcquisition","InfoPersistenceSteps","InfoRetentionProb",
            "OutsideNodes","Exact","ExactMaxStates","ReturnOperators")
  call_args <- args[intersect(names(args), keep)]
  do.call(INApestMetaPathogenAnalytical, call_args)
}
###############################################################################
### INApest multiple-land-use pathogen analytical methods
###
### The MLU pathogen engine treats each node x land-use cell as a pathogen unit.
### R's column-major matrix flattening is preserved exactly: unit index
###   u = node + (land_use - 1) * n_nodes.
### This module therefore reuses the validated Meta exact network operators on
### the flattened pathogen-unit graph, then reshapes results back to node x LU.
###############################################################################

.ina_mlu_check_matrix <- function(x, n_nodes = NULL, n_landuses = NULL, name = deparse(substitute(x)), integer = FALSE) {
  if (!is.matrix(x)) stop(name, " must be a nodes x land-uses matrix")
  if (!is.null(n_nodes) && nrow(x) != n_nodes) stop(name, " has the wrong number of nodes")
  if (!is.null(n_landuses) && ncol(x) != n_landuses) stop(name, " has the wrong number of land uses")
  if (any(!is.finite(x))) stop(name, " must contain finite values")
  if (integer && any(x < 0 | x != floor(x))) stop(name, " must contain non-negative whole numbers")
  x
}

.ina_mlu_unit_index <- function(n_nodes, n_landuses) {
  z <- expand.grid(node = seq_len(n_nodes), land_use = seq_len(n_landuses), KEEP.OUT.ATTRS = FALSE)
  z$unit <- z$node + (z$land_use - 1L) * n_nodes
  z <- z[order(z$unit), , drop = FALSE]
  rownames(z) <- NULL
  z
}

.ina_mlu_resolve_unit <- function(x, n_nodes, n_landuses, name, integer = FALSE) {
  n_units <- n_nodes * n_landuses
  if (is.matrix(x)) {
    if (!all(dim(x) == c(n_nodes, n_landuses))) stop(name, " matrix must be nodes x land uses")
    z <- as.numeric(x)
  } else {
    z <- as.numeric(x)
    if (length(z) == 1L) z <- rep(z, n_units)
    else if (length(z) == n_units) z <- z
    else if (length(z) == n_landuses && n_landuses != n_nodes)
      z <- as.numeric(matrix(rep(z, each = n_nodes), nrow = n_nodes))
    else if (length(z) == n_nodes && n_nodes != n_landuses)
      z <- as.numeric(matrix(rep(z, n_landuses), nrow = n_nodes))
    else stop(name, " must be scalar, nodes x land uses, or an unambiguous node/LU/unit vector")
  }
  if (any(!is.finite(z))) stop(name, " must be finite")
  if (integer && any(z < 0 | z != floor(z))) stop(name, " must contain non-negative whole numbers")
  z
}

INApestMLUPathogenContactMatrix <- function(NodeContact, LandUseMixing = NULL) {
  NodeContact <- as.matrix(NodeContact)
  if (nrow(NodeContact) != ncol(NodeContact) || any(!is.finite(NodeContact)) || any(NodeContact < 0))
    stop("NodeContact must be a finite non-negative square source x target matrix")
  if (is.null(LandUseMixing)) return(NodeContact)
  LandUseMixing <- as.matrix(LandUseMixing)
  if (nrow(LandUseMixing) != ncol(LandUseMixing) || any(!is.finite(LandUseMixing)) || any(LandUseMixing < 0))
    stop("LandUseMixing must be a finite non-negative square source x target matrix")
  kronecker(LandUseMixing, NodeContact)
}

.ina_mlu_contact <- function(n_nodes, n_landuses, ContactMatrix = NULL,
                             NodeContact = NULL, LandUseMixing = NULL) {
  n_units <- n_nodes * n_landuses
  if (!is.null(ContactMatrix)) {
    C <- as.matrix(ContactMatrix)
    if (!all(dim(C) == c(n_units, n_units)) || any(!is.finite(C)) || any(C < 0))
      stop("ContactMatrix must be a finite non-negative pathogen-units x pathogen-units matrix")
    return(C)
  }
  if (is.null(NodeContact)) NodeContact <- diag(n_nodes)
  if (is.null(LandUseMixing)) LandUseMixing <- diag(n_landuses)
  C <- INApestMLUPathogenContactMatrix(NodeContact, LandUseMixing)
  if (!all(dim(C) == c(n_units, n_units))) stop("Combined contact matrix dimensions do not match nodes x land uses")
  C
}

.ina_mlu_reshape_time_unit <- function(x, n_nodes, n_landuses) {
  # x is timesteps x units. Return nodes x LU x timesteps.
  tt <- nrow(x)
  out <- array(0, dim = c(n_nodes, n_landuses, tt),
               dimnames = list(node = seq_len(n_nodes), land_use = seq_len(n_landuses), timestep = 0:(tt - 1L)))
  for (t in seq_len(tt)) out[, , t] <- matrix(x[t, ], nrow = n_nodes, ncol = n_landuses)
  out
}

INApestMLUPathogenGrowthOperator <- function(
    Model = c("SIS", "SIR", "SEIR"),
    HostPopulation,
    Beta,
    RecoveryProb,
    ProgressionProb = 1,
    PathogenMortalityProb = 0,
    ContactMatrix = NULL,
    NodeContact = NULL,
    LandUseMixing = NULL,
    Transmission = c("frequency", "density"),
    DensityScale = 1) {
  Model <- match.arg(Model); Transmission <- match.arg(Transmission)
  HostPopulation <- .ina_mlu_check_matrix(HostPopulation, name = "HostPopulation", integer = FALSE)
  n_nodes <- nrow(HostPopulation); n_landuses <- ncol(HostPopulation); n_units <- n_nodes * n_landuses
  N <- as.numeric(HostPopulation)
  beta <- .ina_mlu_resolve_unit(Beta, n_nodes, n_landuses, "Beta")
  rec <- .ina_mlu_resolve_unit(RecoveryProb, n_nodes, n_landuses, "RecoveryProb")
  prog <- .ina_mlu_resolve_unit(ProgressionProb, n_nodes, n_landuses, "ProgressionProb")
  mort <- .ina_mlu_resolve_unit(PathogenMortalityProb, n_nodes, n_landuses, "PathogenMortalityProb")
  ds <- .ina_mlu_resolve_unit(DensityScale, n_nodes, n_landuses, "DensityScale")
  C <- .ina_mlu_contact(n_nodes, n_landuses, ContactMatrix, NodeContact, LandUseMixing)
  g <- INApestMetaPathogenGrowthOperator(Model, N, beta, rec, prog, mort, C, Transmission, ds)
  map <- .ina_mlu_unit_index(n_nodes, n_landuses)
  active_states <- if (Model == "SEIR") 2L else 1L
  labels <- if (Model == "SEIR") {
    c(paste0("E_n", map$node, "_lu", map$land_use), paste0("I_n", map$node, "_lu", map$land_use))
  } else paste0("I_n", map$node, "_lu", map$land_use)
  dimnames(g$Operator) <- list(labels, labels)
  if (!is.null(g$TransmissionBlock))
    dimnames(g$TransmissionBlock) <- list(paste0("recipient_n", map$node, "_lu", map$land_use), paste0("source_n", map$node, "_lu", map$land_use))
  g$Model <- "INApestMetaMultipleLandUse"
  g$HostArchitecture <- "node x land-use pathogen units"
  g$UnitMap <- map
  g$HostPopulation <- HostPopulation
  g$Diagnostics <- c(g$Diagnostics,
                     "MLU flattening matches the stochastic engine: all nodes in land use 1, then all nodes in land use 2, and so on.",
                     "Land-use mixing and node contact can be supplied separately; the combined source x target contact matrix is kronecker(LandUseMixing, NodeContact).")
  class(g) <- c("INApestMLUPathogenGrowthOperator", "list")
  g
}

INApestMLUPathogenExactFixedN <- function(
    Model = c("SIS", "SIR", "SEIR"),
    Ntimesteps,
    HostPopulation,
    InitialInfected,
    InitialExposed = 0,
    InitialRecovered = 0,
    Beta,
    RecoveryProb,
    ProgressionProb = 1,
    ImmunityLossProb = 0,
    IntroductionProb = 0,
    IntroductionNumber = 1,
    ContactMatrix = NULL,
    NodeContact = NULL,
    LandUseMixing = NULL,
    Transmission = c("frequency", "density"),
    DensityScale = 1,
    OutsideNodes = integer(0),
    OutsideUnits = integer(0),
    MaxStates = 50000L,
    ReturnOperator = FALSE) {
  Model <- match.arg(Model); Transmission <- match.arg(Transmission)
  Nmat <- .ina_mlu_check_matrix(HostPopulation, name = "HostPopulation", integer = TRUE)
  n_nodes <- nrow(Nmat); n_landuses <- ncol(Nmat); n_units <- n_nodes * n_landuses
  N <- as.integer(as.numeric(Nmat))
  I0 <- as.integer(.ina_mlu_resolve_unit(InitialInfected, n_nodes, n_landuses, "InitialInfected", TRUE))
  E0 <- as.integer(.ina_mlu_resolve_unit(InitialExposed, n_nodes, n_landuses, "InitialExposed", TRUE))
  R0 <- as.integer(.ina_mlu_resolve_unit(InitialRecovered, n_nodes, n_landuses, "InitialRecovered", TRUE))
  if (any(I0 + E0 + R0 > N)) stop("Initial pathogen-state counts exceed HostPopulation")
  beta <- .ina_mlu_resolve_unit(Beta, n_nodes, n_landuses, "Beta")
  rec <- .ina_mlu_resolve_unit(RecoveryProb, n_nodes, n_landuses, "RecoveryProb")
  prog <- .ina_mlu_resolve_unit(ProgressionProb, n_nodes, n_landuses, "ProgressionProb")
  wan <- .ina_mlu_resolve_unit(ImmunityLossProb, n_nodes, n_landuses, "ImmunityLossProb")
  ip <- .ina_mlu_resolve_unit(IntroductionProb, n_nodes, n_landuses, "IntroductionProb")
  inum <- as.integer(.ina_mlu_resolve_unit(IntroductionNumber, n_nodes, n_landuses, "IntroductionNumber", TRUE))
  ds <- .ina_mlu_resolve_unit(DensityScale, n_nodes, n_landuses, "DensityScale")
  C <- .ina_mlu_contact(n_nodes, n_landuses, ContactMatrix, NodeContact, LandUseMixing)
  outside <- unique(as.integer(OutsideUnits))
  if (length(OutsideNodes)) {
    if (any(!OutsideNodes %in% seq_len(n_nodes))) stop("OutsideNodes contains invalid node indices")
    map <- .ina_mlu_unit_index(n_nodes, n_landuses)
    outside <- unique(c(outside, map$unit[map$node %in% OutsideNodes]))
  }
  if (any(!outside %in% seq_len(n_units))) stop("OutsideUnits contains invalid pathogen-unit indices")
  if (Model == "SIS") {
    ans <- INApestMetaSISExactNetworkFixedN(Ntimesteps, N, I0, beta, rec, C, ip, inum,
                                            Transmission, ds, OutsideNodes = outside,
                                            ReturnOperator = ReturnOperator)
    EI <- ans$ExpectedInfectedByNode
    EE <- NULL
  } else {
    ans <- INApestMetaCompartmentExactNetworkFixedN(Model, Ntimesteps, N, I0, E0, R0,
                                                     beta, rec, prog, wan, ip, inum, C,
                                                     Transmission, ds, OutsideNodes = outside,
                                                     MaxStates = MaxStates,
                                                     ReturnOperator = ReturnOperator)
    EI <- ans$ExpectedInfectedByNode
    EE <- ans$ExpectedExposedByNode
  }
  map <- .ina_mlu_unit_index(n_nodes, n_landuses)
  out <- ans
  out$Model <- "INApestMetaMultipleLandUse"
  out$HostArchitecture <- "fixed abundance by node x land use"
  out$PathogenModel <- Model
  out$UnitMap <- map
  out$HostPopulation <- Nmat
  out$ExpectedInfectedByUnit <- EI
  out$ExpectedInfected <- .ina_mlu_reshape_time_unit(EI, n_nodes, n_landuses)
  if (!is.null(EE)) {
    out$ExpectedExposedByUnit <- EE
    out$ExpectedExposed <- .ina_mlu_reshape_time_unit(EE, n_nodes, n_landuses)
  }
  if (!is.null(out$Escape)) {
    out$Escape$OutsideUnits <- outside
    out$Escape$OutsideNodes <- as.integer(OutsideNodes)
  }
  out$Growth <- INApestMLUPathogenGrowthOperator(Model, Nmat, beta, rec, prog,
                                                  PathogenMortalityProb = if (Model == "SIS") 0 else 0,
                                                  ContactMatrix = C, Transmission = Transmission,
                                                  DensityScale = ds)
  out$Diagnostics <- c(out$Diagnostics,
                       "Exact MLU solution is obtained on the flattened node x land-use pathogen-unit graph and then reshaped back to nodes x land uses.",
                       "This fixed-host exact branch is intended as a small-system validator; host turnover and management are added in separate MLU analytical branches.")
  class(out) <- c("INApestMLUPathogenExactFixedN", "list")
  out
}

print.INApestMLUPathogenExactFixedN <- function(x, ...) {
  cat("INApest MLU pathogen exact fixed-host solution\n")
  cat("Pathogen model:", x$PathogenModel, "\n")
  cat("Nodes x land uses:", nrow(x$HostPopulation), "x", ncol(x$HostPopulation), "\n")
  cat("Rare-pathogen multiplier:", format(x$Growth$IntrinsicRarePathogenMultiplier, digits = 7), "\n")
  if (!is.null(x$Escape)) cat("Escape probability by horizon:", format(x$Escape$ProbabilityByHorizon, digits = 7), "\n")
  invisible(x)
}

###############################################################################
### Exact one-node MLU SIS with explicit gross host turnover
###############################################################################

.ina_mlu_sis_local_ni_states <- function(K) {
  rows <- list(); k <- 0L
  for (N in 0:as.integer(K)) for (I in 0:N) { k <- k + 1L; rows[[k]] <- c(N=N,I=I) }
  z <- as.data.frame(do.call(rbind, rows)); z$N <- as.integer(z$N); z$I <- as.integer(z$I); z
}

.ina_mlu_sis_turnover_state_table <- function(K) {
  K <- as.integer(K); L <- length(K)
  local <- lapply(K, .ina_mlu_sis_local_ni_states)
  dims <- vapply(local, nrow, integer(1))
  grid <- as.matrix(expand.grid(lapply(dims, seq_len), KEEP.OUT.ATTRS=FALSE)); storage.mode(grid) <- "integer"
  out <- matrix(0L, nrow(grid), 2L*L)
  cn <- character(2L*L)
  for (lu in seq_len(L)) {
    out[,2L*lu-1L] <- local[[lu]]$N[grid[,lu]]
    out[,2L*lu] <- local[[lu]]$I[grid[,lu]]
    cn[2L*lu-1L] <- paste0("N_lu",lu); cn[2L*lu] <- paste0("I_lu",lu)
  }
  colnames(out) <- cn
  attr(out,"LocalStates") <- local
  out
}

.ina_mlu_sis_turnover_operator_one_node <- function(
    K,
    Survival = 1,
    RecruitToCapacity = FALSE,
    Beta = 0,
    RecoveryProb = 0,
    PathogenMortalityProb = 0,
    IntroductionProb = 0,
    IntroductionNumber = 1,
    LandUseMixing = NULL,
    Transmission = c("frequency","density"),
    DensityScale = 1) {
  Transmission <- match.arg(Transmission); K <- as.integer(K); L <- length(K)
  if(any(K<0))stop("K must be non-negative whole numbers")
  res <- function(x,name,int=FALSE){z<-.ina_mlu_resolve_unit(x,1L,L,name,int);z}
  surv<-res(Survival,"Survival");beta<-res(Beta,"Beta");rec<-res(RecoveryProb,"RecoveryProb");pm<-res(PathogenMortalityProb,"PathogenMortalityProb");ip<-res(IntroductionProb,"IntroductionProb");inum<-as.integer(res(IntroductionNumber,"IntroductionNumber",TRUE));ds<-res(DensityScale,"DensityScale")
  if(any(surv<0|surv>1)||any(beta<0)||any(rec<0|rec>1)||any(pm<0|pm>1)||any(rec+pm>1+1e-12)||any(ip<0|ip>1)||any(ds<=0))stop("Invalid turnover/pathogen parameters")
  C <- .ina_mlu_contact(1L,L,NodeContact=matrix(1,1,1),LandUseMixing=if(is.null(LandUseMixing))diag(L)else LandUseMixing)
  st <- .ina_mlu_sis_turnover_state_table(K); keys <- apply(st,1,paste,collapse=":"); idx <- setNames(seq_len(nrow(st)),keys); T <- matrix(0,nrow(st),nrow(st),dimnames=list(keys,keys))
  for(row in seq_len(nrow(st))) {
    N0 <- as.integer(st[row,seq(1,2L*L,by=2)]); I0 <- as.integer(st[row,seq(2,2L*L,by=2)]); S0 <- N0-I0
    host_lists <- vector("list",L)
    for(lu in seq_len(L)) {
      z <- list(); kk <- 0L
      for(ss in 0:S0[lu]) for(ii in 0:I0[lu]) {
        p <- dbinom(ss,S0[lu],surv[lu])*dbinom(ii,I0[lu],surv[lu]); if(p==0)next
        Npre <- ss+ii
        if(isTRUE(RecruitToCapacity)) Npre <- K[lu]
        kk<-kk+1L;z[[kk]]<-c(N=Npre,I=ii,prob=p)
      }
      host_lists[[lu]] <- as.data.frame(do.call(rbind,z))
    }
    hg <- as.matrix(expand.grid(lapply(host_lists,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(hg)<-"integer"
    for(hb in seq_len(nrow(hg))) {
      Npre<-Ipre<-numeric(L);ph<-1
      for(lu in seq_len(L)){z<-host_lists[[lu]][hg[hb,lu],];Npre[lu]<-z$N;Ipre[lu]<-z$I;ph<-ph*z$prob}
      if(ph==0)next
      pressure<-as.numeric(crossprod(Ipre,C));den<-as.numeric(crossprod(Npre,C));foi<-if(Transmission=="frequency")beta*ifelse(den>0,pressure/den,0)else beta*pressure/ds;pinf<-.ina_meta_sis_clip01(-expm1(-pmax(0,foi)))
      dlist<-vector("list",L)
      for(lu in seq_len(L)) dlist[[lu]] <- as.data.frame(.ina_meta_sis_disease_outcomes_pinf(as.integer(Npre[lu]),as.integer(Ipre[lu]),pinf[lu],rec[lu],pm[lu],ip[lu],inum[lu]))
      dg<-as.matrix(expand.grid(lapply(dlist,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(dg)<-"integer"
      for(db in seq_len(nrow(dg))) {
        dest<-integer(2L*L);pd<-ph
        for(lu in seq_len(L)){z<-dlist[[lu]][dg[db,lu],];dest[2L*lu-1L]<-as.integer(z$N);dest[2L*lu]<-as.integer(z$I);pd<-pd*z$prob}
        ky<-paste(dest,collapse=":");T[row,idx[[ky]]]<-T[row,idx[[ky]]]+pd
      }
    }
  }
  er<-max(abs(rowSums(T)-1));if(er>1e-10)stop("MLU turnover operator rows do not sum to 1; max error=",er)
  attr(T,"States")<-st;T
}

INApestMLUSISExactTurnoverOneNode <- function(
    Ntimesteps,
    K,
    InitialHostPopulation,
    InitialInfected,
    Survival = 1,
    RecruitToCapacity = FALSE,
    Beta = 0,
    RecoveryProb = 0,
    PathogenMortalityProb = 0,
    IntroductionProb = 0,
    IntroductionNumber = 1,
    LandUseMixing = NULL,
    Transmission = c("frequency","density"),
    DensityScale = 1,
    ReturnOperator = FALSE) {
  Transmission<-match.arg(Transmission);K<-as.integer(K);L<-length(K)
  N0<-as.integer(.ina_mlu_resolve_unit(InitialHostPopulation,1L,L,"InitialHostPopulation",TRUE));I0<-as.integer(.ina_mlu_resolve_unit(InitialInfected,1L,L,"InitialInfected",TRUE));if(any(N0>K)||any(I0>N0))stop("Invalid initial MLU turnover state")
  T<-.ina_mlu_sis_turnover_operator_one_node(K,Survival,RecruitToCapacity,Beta,RecoveryProb,PathogenMortalityProb,IntroductionProb,IntroductionNumber,LandUseMixing,Transmission,DensityScale)
  st<-attr(T,"States");keys<-rownames(T);init<-as.vector(rbind(N0,I0));dist<-numeric(nrow(st));dist[match(paste(init,collapse=":"),keys)]<-1
  EN<-EI<-matrix(0,Ntimesteps+1L,L);EN[1,]<-N0;EI[1,]<-I0;presence<-numeric(Ntimesteps+1L);presence[1]<-as.numeric(any(I0>0));full<-matrix(0,Ntimesteps+1L,nrow(st));full[1,]<-dist
  Ncols<-seq(1,2L*L,by=2);Icols<-seq(2,2L*L,by=2)
  if(Ntimesteps>0)for(tt in seq_len(Ntimesteps)){dist<-as.numeric(dist%*%T);full[tt+1,]<-dist;EN[tt+1,]<-as.numeric(dist%*%st[,Ncols,drop=FALSE]);EI[tt+1,]<-as.numeric(dist%*%st[,Icols,drop=FALSE]);presence[tt+1]<-sum(dist[rowSums(st[,Icols,drop=FALSE])>0])}
  out<-list(Model="INApestMetaMultipleLandUse",PathogenModel="SIS",HostAssumption=if(RecruitToCapacity)"gross survival followed by susceptible replacement to K in each land use"else"gross survival with no recruitment",Exact=TRUE,ExpectedHostByLandUse=EN,ExpectedInfectedByLandUse=EI,PathogenPresenceProbability=presence,PathogenFreedomProbability=1-presence,StateTable=st,StateDistribution=full,Diagnostics=c("Exact one-node MLU SIS turnover validator.","Host survival is applied to actual S/I hosts before any replacement recruitment; replacement hosts enter S.","Pathogen transmission occurs only after host turnover and uses the post-turnover land-use abundances."));if(ReturnOperator)out$Operator<-T;class(out)<-c("INApestMLUSISExactTurnoverOneNode","list");out
}

###############################################################################
### Shared node-level information / management across multiple land uses
###
### The stochastic MLU engine stores host/pathogen states by node x land use,
### but HaveInfo and LastKnownPresence are node-level states. These exact
### validators preserve that architecture rather than treating land uses as
### independent information units.
###############################################################################

.ina_mlu_model_states <- function(Model) {
  Model <- match.arg(Model, c("SIS","SIR","SEIR"))
  if (Model == "SIS") c("S","I") else if (Model == "SIR") c("S","I","R") else c("S","E","I","R")
}

.ina_mlu_local_compartment_states <- function(Model, K) {
  Model <- match.arg(Model, c("SIS","SIR","SEIR")); K <- as.integer(K)
  if (K < 0L) stop("K must be non-negative")
  if (Model == "SIS") {
    out <- list(); kk <- 0L
    for (N in 0:K) for (I in 0:N) {
      kk <- kk + 1L; out[[kk]] <- c(S=N-I, I=I)
    }
    z <- as.data.frame(do.call(rbind,out)); z$S <- as.integer(z$S); z$I <- as.integer(z$I); z$N <- z$S+z$I
    return(z)
  }
  .ina_meta_compartment_states_uptoK(Model, K)
}

.ina_mlu_info_age_rows <- function(InfoPersistenceSteps=NA) {
  if (is.na(InfoPersistenceSteps)) return(data.frame(H=c(0L,1L), A=c(-2L,-2L)))
  Kinfo <- as.integer(InfoPersistenceSteps)
  if (Kinfo < 0L || Kinfo != InfoPersistenceSteps) stop("InfoPersistenceSteps must be NA or a non-negative whole number")
  ages <- unique(c(-1L, 0:max(0L,Kinfo-1L)))
  rbind(data.frame(H=0L,A=-2L), data.frame(H=1L,A=ages))
}

.ina_mlu_one_node_info_states <- function(Model, K, InfoPersistenceSteps=NA, MaxStates=50000L) {
  Model <- match.arg(Model,c("SIS","SIR","SEIR")); K <- as.integer(K); L <- length(K)
  if (!L || any(K < 0L)) stop("K must contain non-negative land-use capacities")
  local <- lapply(K, function(k) .ina_mlu_local_compartment_states(Model,k))
  dims <- vapply(local,nrow,integer(1)); bio_n <- prod(dims)
  ages <- .ina_mlu_info_age_rows(InfoPersistenceSteps)
  nst <- bio_n*nrow(ages)
  if (nst > MaxStates) stop("Exact one-node MLU information state space has ",nst," states; increase MaxStates deliberately or use stochastic MLU")
  g <- as.matrix(expand.grid(lapply(dims,seq_len),KEEP.OUT.ATTRS=FALSE)); storage.mode(g)<-"integer"
  states <- .ina_mlu_model_states(Model)
  bio <- matrix(0L,nrow(g),length(states)*L)
  cn <- character(ncol(bio))
  for (lu in seq_len(L)) for (j in seq_along(states)) {
    col <- (lu-1L)*length(states)+j
    bio[,col] <- local[[lu]][[states[j]]][g[,lu]]
    cn[col] <- paste0(states[j],"_lu",lu)
  }
  colnames(bio) <- cn
  out <- bio[rep(seq_len(nrow(bio)), each=nrow(ages)),,drop=FALSE]
  aa <- ages[rep(seq_len(nrow(ages)), times=nrow(bio)),,drop=FALSE]
  out <- cbind(out,H=as.integer(aa$H),A=as.integer(aa$A))
  storage.mode(out)<-"integer"
  attr(out,"ModelStates")<-states; attr(out,"K")<-K; out
}

.ina_mlu_compartment_disease_outcomes_pinf <- function(state, Model, p_inf,
                                                        RecoveryProb=0, ProgressionProb=1,
                                                        ImmunityLossProb=0, PathogenMortalityProb=0,
                                                        IntroductionProb=0, IntroductionNumber=1) {
  Model <- match.arg(Model,c("SIR","SEIR")); states <- .ina_mlu_model_states(Model)
  S0 <- as.integer(state[["S"]]); I0 <- as.integer(state[["I"]]); R0 <- as.integer(state[["R"]]); E0 <- if(Model=="SEIR")as.integer(state[["E"]])else 0L
  rec <- as.numeric(RecoveryProb); prog <- as.numeric(ProgressionProb); wan <- as.numeric(ImmunityLossProb); mort <- as.numeric(PathogenMortalityProb)
  ip <- as.numeric(IntroductionProb); inum <- as.integer(IntroductionNumber); p_inf <- .ina_meta_sis_clip01(p_inf)
  if(any(!is.finite(c(rec,prog,wan,mort,ip)))||any(c(rec,prog,wan,mort,ip)<0)||any(c(rec,prog,wan,mort,ip)>1)||rec+mort>1+1e-12) stop("Invalid disease probabilities")
  out<-list();kk<-0L
  for(x in 0:S0){px<-dbinom(x,S0,p_inf);remS<-S0-x;mx<-min(inum,remS);ivals<-if(ip>0&&mx>0)c(0L,mx)else 0L;ips<-if(ip>0&&mx>0)c(1-ip,ip)else 1
    for(iv in seq_along(ivals)){intro<-ivals[iv];pi<-ips[iv]
      for(pg in if(Model=="SEIR")0:E0 else 0L){pp<-if(Model=="SEIR")dbinom(pg,E0,prog)else 1
        for(lose in 0:R0){pl<-dbinom(lose,R0,wan)
          for(stay in 0:I0)for(recover in 0:(I0-stay)){die<-I0-stay-recover;pr<-dmultinom(c(stay,recover,die),prob=c(1-rec-mort,rec,mort));if(pr==0)next
            if(Model=="SIR") z<-c(S=S0-x-intro+lose,I=stay+x+intro,R=R0-lose+recover)
            else z<-c(S=S0-x-intro+lose,E=E0-pg+x+intro,I=stay+pg,R=R0-lose+recover)
            kk<-kk+1L;out[[kk]]<-c(z,prob=unname(px*pi*pp*pl*pr))
  }}}}}
  z<-as.data.frame(do.call(rbind,out));for(nm in states)z[[nm]]<-as.integer(z[[nm]]);z$prob<-as.numeric(z$prob);z<-z[z$prob>0,,drop=FALSE]
  key<-apply(z[,states,drop=FALSE],1,paste,collapse=":");ps<-tapply(z$prob,key,sum);first<-match(names(ps),key);ans<-z[first,states,drop=FALSE];ans$prob<-as.numeric(ps);rownames(ans)<-NULL;ans
}

.ina_mlu_local_management_mortality_branches <- function(counts, states, Survival, ManageProb, MortalityProb, H) {
  out<-list();kk<-0L;mprob<-ManageProb*H
  for(managing in 0:1){pm<-if(managing==1)mprob else 1-mprob;if(pm<=0)next;mm<-MortalityProb*managing
    ls<-lapply(states,function(nm).ina_meta_class_mortality_outcomes(as.integer(counts[[nm]]),Survival,mm))
    g<-as.matrix(expand.grid(lapply(ls,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(g)<-"integer"
    for(rr in seq_len(nrow(g))){surv<-setNames(integer(length(states)),states);md<-0L;p<-pm
      for(j in seq_along(states)){z<-ls[[j]][g[rr,j],];surv[j]<-as.integer(z["survive"]);md<-md+as.integer(z["manage_death"]);p<-p*unname(z["prob"])}
      if(p>0){kk<-kk+1L;out[[kk]]<-c(surv,manage_deaths=md,managing=managing,prob=p)}
    }
  }
  z<-as.data.frame(do.call(rbind,out));for(nm in states)z[[nm]]<-as.integer(z[[nm]]);z$manage_deaths<-as.integer(z$manage_deaths);z$managing<-as.integer(z$managing);z$prob<-as.numeric(z$prob);z
}

.ina_mlu_one_node_info_prebranches <- function(state_row, Model, L, Survival, ManageProb, MortalityProb,
                                                host_acq, InfoPersistenceSteps, InfoRetentionProb) {
  states <- .ina_mlu_model_states(Model); getv<-function(nm)as.integer(if(is.matrix(state_row)||is.data.frame(state_row))state_row[1,nm]else state_row[[nm]]); H0<-getv("H");A0<-getv("A");finite<-!is.na(InfoPersistenceSteps);Kinfo<-if(finite)as.integer(InfoPersistenceSteps)else NA_integer_
  local<-vector("list",L)
  for(lu in seq_len(L)){
    cc<-setNames(vapply(states,function(nm)getv(paste0(nm,"_lu",lu)),integer(1)),states)
    local[[lu]]<-.ina_mlu_local_management_mortality_branches(cc,states,Survival[lu],ManageProb[lu],MortalityProb[lu],H0)
  }
  g<-as.matrix(expand.grid(lapply(local,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(g)<-"integer";out<-list();kk<-0L
  for(rr in seq_len(nrow(g))){vals<-integer(length(states)*L);names(vals)<-as.vector(vapply(seq_len(L),function(lu)paste0(states,"_lu",lu),character(length(states))));p<-1;md<-0L;mg<-integer(L)
    for(lu in seq_len(L)){z<-local[[lu]][g[rr,lu],];for(nm in states)vals[paste0(nm,"_lu",lu)]<-as.integer(z[[nm]]);md<-md+as.integer(z$manage_deaths);mg[lu]<-as.integer(z$managing);p<-p*z$prob}
    if(p==0)next;evidence<-host_acq&&md>0
    info<-list()
    if(finite){
      if(H0==0)info[[1]]<-c(H=0,A=-2,p=1)
      else if(evidence){if(Kinfo<=0)info[[1]]<-c(H=0,A=-2,p=1)else info[[1]]<-c(H=1,A=0,p=1)}
      else if(A0<0)info[[1]]<-c(H=0,A=-2,p=1)
      else {ag<-A0+1L;if(ag>=Kinfo)info[[1]]<-c(H=0,A=-2,p=1)else info[[1]]<-c(H=1,A=ag,p=1)}
    } else {
      if(H0==1L&&InfoRetentionProb<1){info[[1]]<-c(H=1,A=-2,p=InfoRetentionProb);info[[2]]<-c(H=0,A=-2,p=1-InfoRetentionProb)} else info[[1]]<-c(H=H0,A=-2,p=1)
    }
    for(q in info){kk<-kk+1L;out[[kk]]<-c(vals,H=as.integer(q["H"]),A=as.integer(q["A"]),manage_pattern=sum(mg*2^(seq_len(L)-1L)),prob=unname(p*q["p"]))}
  }
  z<-as.data.frame(do.call(rbind,out));for(nm in c(names(vals),"H","A","manage_pattern"))z[[nm]]<-as.integer(z[[nm]]);z$prob<-as.numeric(z$prob);z
}

.ina_mlu_info_operator_one_node <- function(Model,K,Survival=1,ManageProb=0,MortalityProb=0,
                                             HostDetectionProb=0,PathogenDetectionProb=0,
                                             InformationAcquisition=c("host","pathogen","both"),
                                             InfoPersistenceSteps=NA,InfoRetentionProb=1,
                                             Beta=0,RecoveryProb=0,ProgressionProb=1,ImmunityLossProb=0,
                                             PathogenMortalityProb=0,IntroductionProb=0,IntroductionNumber=1,
                                             LandUseMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,
                                             MaxStates=50000L) {
  Model<-match.arg(Model,c("SIS","SIR","SEIR"));InformationAcquisition<-match.arg(InformationAcquisition);Transmission<-match.arg(Transmission);K<-as.integer(K);L<-length(K);states<-.ina_mlu_model_states(Model)
  res<-function(x,nm,int=FALSE).ina_mlu_resolve_unit(x,1L,L,nm,int)
  surv<-res(Survival,"Survival");mng<-res(ManageProb,"ManageProb");mort<-res(MortalityProb,"MortalityProb");hd<-res(HostDetectionProb,"HostDetectionProb");pd<-res(PathogenDetectionProb,"PathogenDetectionProb");beta<-res(Beta,"Beta");rec<-res(RecoveryProb,"RecoveryProb");prog<-res(ProgressionProb,"ProgressionProb");wan<-res(ImmunityLossProb,"ImmunityLossProb");pmort<-res(PathogenMortalityProb,"PathogenMortalityProb");ip<-res(IntroductionProb,"IntroductionProb");inum<-as.integer(res(IntroductionNumber,"IntroductionNumber",TRUE));ds<-res(DensityScale,"DensityScale")
  probs<-c(surv,mng,mort,hd,pd,rec,prog,wan,pmort,ip,InfoRetentionProb);if(any(!is.finite(probs))||any(probs<0)||any(probs>1)||any(rec+pmort>1+1e-12)||any(beta<0)||any(ds<=0))stop("Invalid MLU information/pathogen parameters")
  C<-.ina_mlu_contact(1L,L,NodeContact=matrix(1,1,1),LandUseMixing=if(is.null(LandUseMixing))diag(L)else LandUseMixing)
  st<-.ina_mlu_one_node_info_states(Model,K,InfoPersistenceSteps,MaxStates);bio_cols<-setdiff(colnames(st),c("H","A"));keys<-apply(st[,c(bio_cols,"H","A"),drop=FALSE],1,paste,collapse=":");idx<-setNames(seq_len(nrow(st)),keys);T<-matrix(0,nrow(st),nrow(st),dimnames=list(keys,keys));host_acq<-InformationAcquisition%in%c("host","both");path_acq<-InformationAcquisition%in%c("pathogen","both");finite<-!is.na(InfoPersistenceSteps)
  for(row in seq_len(nrow(st))){pre<-.ina_mlu_one_node_info_prebranches(st[row,],Model,L,surv,mng,mort,host_acq,InfoPersistenceSteps,InfoRetentionProb)
    for(pb in seq_len(nrow(pre))){Ipre<-vapply(seq_len(L),function(lu)as.numeric(pre[pb,paste0("I_lu",lu)]),numeric(1));Npre<-vapply(seq_len(L),function(lu)sum(as.numeric(pre[pb,paste0(states,"_lu",lu)])),numeric(1));pressure<-as.numeric(crossprod(Ipre,C));den<-as.numeric(crossprod(Npre,C));foi<-if(Transmission=="frequency")beta*ifelse(den>0,pressure/den,0)else beta*pressure/ds;pinf<-.ina_meta_sis_clip01(-expm1(-pmax(0,foi)))
      dl<-vector("list",L)
      for(lu in seq_len(L)){
        if(Model=="SIS")dl[[lu]]<-as.data.frame(.ina_meta_sis_disease_outcomes_pinf(as.integer(Npre[lu]),as.integer(Ipre[lu]),pinf[lu],rec[lu],pmort[lu],ip[lu],inum[lu]))
        else {cc<-setNames(vapply(states,function(nm)as.integer(pre[pb,paste0(nm,"_lu",lu)]),integer(1)),states);dl[[lu]]<-.ina_mlu_compartment_disease_outcomes_pinf(cc,Model,pinf[lu],rec[lu],prog[lu],wan[lu],pmort[lu],ip[lu],inum[lu])}
      }
      dg<-as.matrix(expand.grid(lapply(dl,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(dg)<-"integer"
      for(dd in seq_len(nrow(dg))){dest<-integer(length(bio_cols));names(dest)<-bio_cols;pbase<-pre$prob[pb]
        N1<-I1<-numeric(L)
        for(lu in seq_len(L)){z<-dl[[lu]][dg[dd,lu],];pbase<-pbase*z$prob
          if(Model=="SIS"){dest[paste0("S_lu",lu)]<-as.integer(z$N-z$I);dest[paste0("I_lu",lu)]<-as.integer(z$I)}else for(nm in states)dest[paste0(nm,"_lu",lu)]<-as.integer(z[[nm]])
          N1[lu]<-sum(dest[paste0(states,"_lu",lu)]);I1[lu]<-dest[paste0("I_lu",lu)]
        }
        if(pbase==0)next;qh<-if(host_acq)1-prod((1-hd)^N1)else 0;qp<-if(path_acq)1-prod((1-pd)^I1)else 0;qev<-1-(1-qh)*(1-qp);Hm<-as.integer(pre$H[pb]);Am<-as.integer(pre$A[pb])
        if(qev>0){ky<-paste(c(dest,1L,if(finite)0L else -2L),collapse=":");T[row,idx[[ky]]]<-T[row,idx[[ky]]]+pbase*qev}
        if(qev<1){ky<-paste(c(dest,Hm,Am),collapse=":");T[row,idx[[ky]]]<-T[row,idx[[ky]]]+pbase*(1-qev)}
      }
    }
  }
  er<-max(abs(rowSums(T)-1));if(er>1e-10)stop("MLU one-node information operator row error ",er);attr(T,"States")<-st;attr(T,"LandUseCount")<-L;T
}

INApestMLUPathogenExactInformationOneNode <- function(Model=c("SIS","SIR","SEIR"),Ntimesteps,K,
    InitialHostPopulation,InitialInfected=0,InitialExposed=0,InitialRecovered=0,InitialInfo=0,
    Survival=1,ManageProb=0,MortalityProb=0,HostDetectionProb=0,PathogenDetectionProb=0,
    InformationAcquisition=c("host","pathogen","both"),InfoPersistenceSteps=NA,InfoRetentionProb=1,
    Beta=0,RecoveryProb=0,ProgressionProb=1,ImmunityLossProb=0,PathogenMortalityProb=0,
    IntroductionProb=0,IntroductionNumber=1,LandUseMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,
    MaxStates=50000L,ReturnOperator=FALSE) {
  Model<-match.arg(Model);InformationAcquisition<-match.arg(InformationAcquisition);Transmission<-match.arg(Transmission);K<-as.integer(K);L<-length(K);states<-.ina_mlu_model_states(Model)
  resi<-function(x,nm).ina_mlu_resolve_unit(x,1L,L,nm,TRUE);N0<-as.integer(resi(InitialHostPopulation,"InitialHostPopulation"));I0<-as.integer(resi(InitialInfected,"InitialInfected"));E0<-if(Model=="SEIR")as.integer(resi(InitialExposed,"InitialExposed"))else integer(L);R0<-if(Model%in%c("SIR","SEIR"))as.integer(resi(InitialRecovered,"InitialRecovered"))else integer(L)
  if(any(N0>K)||any(I0+E0+R0>N0)||!InitialInfo%in%c(0,1))stop("Invalid initial MLU information state")
  T<-.ina_mlu_info_operator_one_node(Model,K,Survival,ManageProb,MortalityProb,HostDetectionProb,PathogenDetectionProb,InformationAcquisition,InfoPersistenceSteps,InfoRetentionProb,Beta,RecoveryProb,ProgressionProb,ImmunityLossProb,PathogenMortalityProb,IntroductionProb,IntroductionNumber,LandUseMixing,Transmission,DensityScale,MaxStates)
  st<-attr(T,"States");bio_cols<-setdiff(colnames(st),c("H","A"));keys<-rownames(T);finite<-!is.na(InfoPersistenceSteps);hd<-.ina_mlu_resolve_unit(HostDetectionProb,1L,L,"HostDetectionProb");pd<-.ina_mlu_resolve_unit(PathogenDetectionProb,1L,L,"PathogenDetectionProb");host_acq<-InformationAcquisition%in%c("host","both");path_acq<-InformationAcquisition%in%c("pathogen","both")
  init<-integer(length(bio_cols));names(init)<-bio_cols
  for(lu in seq_len(L)){init[paste0("S_lu",lu)]<-N0[lu]-I0[lu]-E0[lu]-R0[lu];init[paste0("I_lu",lu)]<-I0[lu];if(Model=="SEIR")init[paste0("E_lu",lu)]<-E0[lu];if(Model%in%c("SIR","SEIR"))init[paste0("R_lu",lu)]<-R0[lu]}
  qh<-if(host_acq)1-prod((1-hd)^N0)else 0;qp<-if(path_acq)1-prod((1-pd)^I0)else 0;qev<-1-(1-qh)*(1-qp);dist<-numeric(nrow(st));
  if(qev>0){ky<-paste(c(init,1L,if(finite)0L else -2L),collapse=":");dist[match(ky,keys)]<-qev}
  if(qev<1){H<-as.integer(InitialInfo);A<-if(H==1L&&finite)-1L else -2L;ky<-paste(c(init,H,A),collapse=":");ii<-match(ky,keys);dist[ii]<-dist[ii]+1-qev}
  full<-matrix(0,Ntimesteps+1L,nrow(st));full[1,]<-dist;PH<-numeric(Ntimesteps+1L);PH[1]<-sum(dist*st[,"H"]);mng<-.ina_mlu_resolve_unit(ManageProb,1L,L,"ManageProb");PM<-matrix(0,Ntimesteps,L);means<-lapply(states,function(nm)matrix(0,Ntimesteps+1L,L));names(means)<-states
  for(nm in states)for(lu in seq_len(L))means[[nm]][1,lu]<-sum(dist*st[,paste0(nm,"_lu",lu)])
  for(tt in seq_len(Ntimesteps)){PM[tt,]<-mng*sum(dist*st[,"H"]);dist<-as.numeric(dist%*%T);full[tt+1,]<-dist;PH[tt+1]<-sum(dist*st[,"H"]);for(nm in states)for(lu in seq_len(L))means[[nm]][tt+1,lu]<-sum(dist*st[,paste0(nm,"_lu",lu)])}
  EN<-matrix(0,Ntimesteps+1L,L);for(nm in states)EN<-EN+means[[nm]]
  out<-list(Model="INApestMetaMultipleLandUse",PathogenModel=Model,Exact=TRUE,InformationAcquisition=InformationAcquisition,ExpectedHostByLandUse=EN,ExpectedStatesByLandUse=means,ExpectedInfectedByLandUse=means$I,InformationProbability=PH,ManagingProbabilityByLandUse=PM,StateTable=st,StateDistribution=full,Diagnostics=c("Exact one-node MLU joint pathogen/information validator with a single node-level HaveInfo state shared across land uses.","Management adoption is drawn independently by land use conditional on shared HaveInfo; actual management deaths are host evidence only when host acquisition is enabled.",if(finite)"Programmed stop is represented by explicit node-level time-since-evidence states."else"Information loss follows the memoryless retention pathway."))
  if(Model=="SEIR")out$ExpectedExposedByLandUse<-means$E;if(ReturnOperator)out$Operator<-T;class(out)<-c("INApestMLUPathogenExactInformationOneNode","list");out
}

###############################################################################
### Exact one-node MLU SIR / SEIR gross-turnover extension
###############################################################################

.ina_mlu_compartment_turnover_operator_one_node <- function(Model,K,Survival=1,RecruitToCapacity=FALSE,
    Beta=0,RecoveryProb=0,ProgressionProb=1,ImmunityLossProb=0,PathogenMortalityProb=0,
    IntroductionProb=0,IntroductionNumber=1,LandUseMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,MaxStates=50000L) {
  Model<-match.arg(Model,c("SIR","SEIR"));Transmission<-match.arg(Transmission);K<-as.integer(K);L<-length(K);states<-.ina_mlu_model_states(Model)
  local<-lapply(K,function(k).ina_meta_compartment_states_uptoK(Model,k));dims<-vapply(local,nrow,integer(1));nst<-prod(dims);if(nst>MaxStates)stop("Exact MLU turnover state space has ",nst," states")
  g<-as.matrix(expand.grid(lapply(dims,seq_len),KEEP.OUT.ATTRS=FALSE));storage.mode(g)<-"integer";st<-matrix(0L,nst,length(states)*L);cn<-character(ncol(st));for(lu in seq_len(L))for(j in seq_along(states)){cc<-(lu-1L)*length(states)+j;st[,cc]<-local[[lu]][[states[j]]][g[,lu]];cn[cc]<-paste0(states[j],"_lu",lu)};colnames(st)<-cn;keys<-apply(st,1,paste,collapse=":");idx<-setNames(seq_len(nst),keys);T<-matrix(0,nst,nst,dimnames=list(keys,keys))
  res<-function(x,nm,int=FALSE).ina_mlu_resolve_unit(x,1L,L,nm,int);surv<-res(Survival,"Survival");beta<-res(Beta,"Beta");rec<-res(RecoveryProb,"RecoveryProb");prog<-res(ProgressionProb,"ProgressionProb");wan<-res(ImmunityLossProb,"ImmunityLossProb");pmort<-res(PathogenMortalityProb,"PathogenMortalityProb");ip<-res(IntroductionProb,"IntroductionProb");inum<-as.integer(res(IntroductionNumber,"IntroductionNumber",TRUE));ds<-res(DensityScale,"DensityScale");C<-.ina_mlu_contact(1L,L,NodeContact=matrix(1,1,1),LandUseMixing=if(is.null(LandUseMixing))diag(L)else LandUseMixing)
  for(row in seq_len(nst)){hl<-vector("list",L);for(lu in seq_len(L)){cc<-setNames(vapply(states,function(nm)as.integer(st[row,paste0(nm,"_lu",lu)]),integer(1)),states);hl[[lu]]<-.ina_meta_compartment_host_outcomes(cc,Model,K[lu],surv[lu],RecruitToCapacity,0,0,NULL)};hg<-as.matrix(expand.grid(lapply(hl,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(hg)<-"integer"
    for(hb in seq_len(nrow(hg))){pre<-vector("list",L);ph<-1;Ipre<-Npre<-numeric(L);for(lu in seq_len(L)){z<-hl[[lu]][hg[hb,lu],];pre[[lu]]<-z;ph<-ph*z$prob;Ipre[lu]<-z$I;Npre[lu]<-sum(as.numeric(z[states]))};pressure<-as.numeric(crossprod(Ipre,C));den<-as.numeric(crossprod(Npre,C));foi<-if(Transmission=="frequency")beta*ifelse(den>0,pressure/den,0)else beta*pressure/ds;pinf<-.ina_meta_sis_clip01(-expm1(-pmax(0,foi)));dl<-vector("list",L);for(lu in seq_len(L))dl[[lu]]<-.ina_mlu_compartment_disease_outcomes_pinf(pre[[lu]],Model,pinf[lu],rec[lu],prog[lu],wan[lu],pmort[lu],ip[lu],inum[lu]);dg<-as.matrix(expand.grid(lapply(dl,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(dg)<-"integer"
      for(dd in seq_len(nrow(dg))){dest<-integer(length(cn));names(dest)<-cn;p<-ph;for(lu in seq_len(L)){z<-dl[[lu]][dg[dd,lu],];p<-p*z$prob;for(nm in states)dest[paste0(nm,"_lu",lu)]<-as.integer(z[[nm]])};ky<-paste(dest,collapse=":");T[row,idx[[ky]]]<-T[row,idx[[ky]]]+p}
    }
  }
  er<-max(abs(rowSums(T)-1));if(er>1e-10)stop(Model," MLU turnover operator row error ",er);attr(T,"States")<-st;T
}

INApestMLUCompartmentExactTurnoverOneNode <- function(Model=c("SIR","SEIR"),Ntimesteps,K,InitialHostPopulation,InitialInfected=0,InitialExposed=0,InitialRecovered=0,Survival=1,RecruitToCapacity=FALSE,Beta=0,RecoveryProb=0,ProgressionProb=1,ImmunityLossProb=0,PathogenMortalityProb=0,IntroductionProb=0,IntroductionNumber=1,LandUseMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,MaxStates=50000L,ReturnOperator=FALSE){
  Model<-match.arg(Model);Transmission<-match.arg(Transmission);K<-as.integer(K);L<-length(K);states<-.ina_mlu_model_states(Model);resi<-function(x,nm).ina_mlu_resolve_unit(x,1L,L,nm,TRUE);N0<-as.integer(resi(InitialHostPopulation,"InitialHostPopulation"));I0<-as.integer(resi(InitialInfected,"InitialInfected"));E0<-if(Model=="SEIR")as.integer(resi(InitialExposed,"InitialExposed"))else integer(L);R0<-as.integer(resi(InitialRecovered,"InitialRecovered"));if(any(N0>K)||any(I0+E0+R0>N0))stop("Invalid initial counts")
  T<-.ina_mlu_compartment_turnover_operator_one_node(Model,K,Survival,RecruitToCapacity,Beta,RecoveryProb,ProgressionProb,ImmunityLossProb,PathogenMortalityProb,IntroductionProb,IntroductionNumber,LandUseMixing,Transmission,DensityScale,MaxStates);st<-attr(T,"States");keys<-rownames(T);init<-integer(length(states)*L);names(init)<-colnames(st);for(lu in seq_len(L)){init[paste0("S_lu",lu)]<-N0[lu]-I0[lu]-E0[lu]-R0[lu];init[paste0("I_lu",lu)]<-I0[lu];if(Model=="SEIR")init[paste0("E_lu",lu)]<-E0[lu];init[paste0("R_lu",lu)]<-R0[lu]};dist<-numeric(nrow(st));dist[match(paste(init,collapse=":"),keys)]<-1;full<-matrix(0,Ntimesteps+1L,nrow(st));full[1,]<-dist;means<-lapply(states,function(nm)matrix(0,Ntimesteps+1L,L));names(means)<-states;for(nm in states)for(lu in seq_len(L))means[[nm]][1,lu]<-init[paste0(nm,"_lu",lu)];active<-numeric(Ntimesteps+1L);active[1]<-as.numeric(any(I0+E0>0))
  for(tt in seq_len(Ntimesteps)){dist<-as.numeric(dist%*%T);full[tt+1,]<-dist;for(nm in states)for(lu in seq_len(L))means[[nm]][tt+1,lu]<-sum(dist*st[,paste0(nm,"_lu",lu)]);active[tt+1]<-if(Model=="SEIR")sum(dist[rowSums(st[,grep("^(E|I)_lu",colnames(st)),drop=FALSE])>0])else sum(dist[rowSums(st[,grep("^I_lu",colnames(st)),drop=FALSE])>0])};EN<-matrix(0,Ntimesteps+1L,L);for(nm in states)EN<-EN+means[[nm]];out<-list(Model="INApestMetaMultipleLandUse",PathogenModel=Model,HostAssumption="explicit one-node MLU turnover",Exact=TRUE,ExpectedHostByLandUse=EN,ExpectedStatesByLandUse=means,ExpectedInfectedByLandUse=means$I,ActivePathogenProbability=active,PathogenFreedomProbability=1-active,StateTable=st,StateDistribution=full,Diagnostics=c("Exact one-node MLU compartment turnover validator.","Gross host survival thins pathogen states before optional susceptible replacement recruitment.","Pathogen transmission then occurs across land uses using the supplied land-use mixing matrix."));if(Model=="SEIR")out$ExpectedExposedByLandUse<-means$E;if(ReturnOperator)out$Operator<-T;class(out)<-c("INApestMLUCompartmentExactTurnoverOneNode","list");out
}

###############################################################################
### MLU pathogen analytical router and unified dispatcher extension
###############################################################################

.ina_mlu_initial_matrix <- function(x,n_nodes,n_landuses,name){
  if(is.matrix(x)){if(!all(dim(x)==c(n_nodes,n_landuses)))stop(name," must be nodes x land uses");z<-x}else{v<-as.numeric(x);if(length(v)==1L)z<-matrix(v,n_nodes,n_landuses)else if(length(v)==n_nodes*n_landuses)z<-matrix(v,n_nodes,n_landuses)else if(length(v)==n_landuses)z<-matrix(rep(v,each=n_nodes),n_nodes,n_landuses)else stop(name," has incompatible dimensions")};if(any(!is.finite(z))||any(z<0)||any(z!=floor(z)))stop(name," must contain non-negative whole numbers");storage.mode(z)<-"integer";z
}

INApestMLUPathogenAnalytical <- function(Ntimesteps,HostPopulation,Pathogen,
    HostAssumption=c("fixed","event_contract_turnover"),K=HostPopulation,Survival=1,RecruitToCapacity=FALSE,
    HostDetectionProb=0,InitialInfo=0,ManageProb=0,MortalityProb=0,InformationAcquisition=NULL,
    InfoPersistenceSteps=NA,InfoRetentionProb=1,NodeContact=NULL,LandUseMixing=NULL,
    OutsideNodes=integer(0),OutsideUnits=integer(0),Exact=TRUE,ExactMaxStates=50000L,ReturnOperators=FALSE){
  if(!inherits(Pathogen,"INApestPathogen"))stop("Pathogen must be returned by INApestPathogen()");Model<-as.character(Pathogen$Model)[1L];if(!Model%in%c("SIS","SIR","SEIR"))stop("MLU analytical router supports SIS, SIR and SEIR");HostAssumption<-match.arg(HostAssumption);N<-.ina_mlu_check_matrix(HostPopulation,name="HostPopulation",integer=TRUE);n_nodes<-nrow(N);L<-ncol(N);Kmat<-if(is.matrix(K))K else matrix(.ina_mlu_resolve_unit(K,n_nodes,L,"K",TRUE),n_nodes,L);if(!all(dim(Kmat)==dim(N))||any(Kmat<N))stop("K must match HostPopulation dimensions and be >= HostPopulation")
  I0<-.ina_mlu_initial_matrix(Pathogen$InitialInfected,n_nodes,L,"InitialInfected");E0<-if(Model=="SEIR").ina_mlu_initial_matrix(Pathogen$InitialExposed,n_nodes,L,"InitialExposed")else matrix(0L,n_nodes,L);R0<-if(Model%in%c("SIR","SEIR")).ina_mlu_initial_matrix(Pathogen$InitialRecovered,n_nodes,L,"InitialRecovered")else matrix(0L,n_nodes,L);if(any(I0+E0+R0>N))stop("Initial pathogen states exceed HostPopulation")
  acq<-InformationAcquisition;if(is.null(acq))acq<-if(isTRUE(Pathogen$DetectionTriggersInfo))"both"else"host";acq<-match.arg(as.character(acq),c("host","pathogen","both"));info_requested<-any(as.numeric(ManageProb)!=0,na.rm=TRUE)||any(as.numeric(InitialInfo)!=0,na.rm=TRUE)||any(as.numeric(HostDetectionProb)!=0,na.rm=TRUE)||any(as.numeric(Pathogen$DetectionProb)!=0,na.rm=TRUE)||any(!is.na(InfoPersistenceSteps))||any(as.numeric(InfoRetentionProb)!=1,na.rm=TRUE)
  C<-.ina_mlu_contact(n_nodes,L,ContactMatrix=Pathogen$ContactMatrix,NodeContact=NodeContact,LandUseMixing=LandUseMixing)
  if(!isTRUE(Exact)){g<-INApestMLUPathogenGrowthOperator(Model,N,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$PathogenMortalityProb,ContactMatrix=C,Transmission=Pathogen$Transmission,DensityScale=Pathogen$DensityScale);return(structure(list(Model="INApestMetaMultipleLandUse",PathogenModel=Model,Exact=FALSE,Growth=g,Diagnostics=c("Scalable MLU result is the rare-pathogen growth operator on the flattened node x land-use graph.","Use stochastic MLU for finite-prevalence managed landscapes once exact state expansion is impractical.")),class=c("INApestMLUPathogenAnalytical","list")))}
  if(HostAssumption=="event_contract_turnover"){
    if(n_nodes!=1L)stop("Exact MLU event-contract turnover is currently a one-node validation solution")
    if(info_requested)stop("MLU turnover and shared-information validators are exact but exposed separately; combined turnover+information is not routed generically")
    if(Model=="SIS")return(INApestMLUSISExactTurnoverOneNode(Ntimesteps,as.integer(Kmat[1,]),as.integer(N[1,]),as.integer(I0[1,]),Survival,RecruitToCapacity,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$PathogenMortalityProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,if(is.null(LandUseMixing))diag(L)else LandUseMixing,Pathogen$Transmission,Pathogen$DensityScale,ReturnOperators))
    return(INApestMLUCompartmentExactTurnoverOneNode(Model,Ntimesteps,as.integer(Kmat[1,]),as.integer(N[1,]),as.integer(I0[1,]),as.integer(E0[1,]),as.integer(R0[1,]),Survival,RecruitToCapacity,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$PathogenMortalityProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,if(is.null(LandUseMixing))diag(L)else LandUseMixing,Pathogen$Transmission,Pathogen$DensityScale,ExactMaxStates,ReturnOperators))
  }
  if(info_requested){
    if(n_nodes!=1L)stop("Exact shared-information MLU expansion is currently exposed for one-node validation systems; use stochastic MLU for larger managed networks")
    return(INApestMLUPathogenExactInformationOneNode(Model,Ntimesteps,as.integer(Kmat[1,]),as.integer(N[1,]),as.integer(I0[1,]),as.integer(E0[1,]),as.integer(R0[1,]),as.integer(InitialInfo)[1],Survival,ManageProb,MortalityProb,HostDetectionProb,Pathogen$DetectionProb,acq,InfoPersistenceSteps,InfoRetentionProb,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$PathogenMortalityProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,if(is.null(LandUseMixing))diag(L)else LandUseMixing,Pathogen$Transmission,Pathogen$DensityScale,ExactMaxStates,ReturnOperators))
  }
  if(any(as.numeric(Pathogen$PathogenMortalityProb)!=0))stop("Fixed-host exact MLU branch requires PathogenMortalityProb=0; use HostAssumption='event_contract_turnover'")
  INApestMLUPathogenExactFixedN(Model,Ntimesteps,N,I0,E0,R0,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,ContactMatrix=C,Transmission=Pathogen$Transmission,DensityScale=Pathogen$DensityScale,OutsideNodes=OutsideNodes,OutsideUnits=OutsideUnits,MaxStates=ExactMaxStates,ReturnOperator=ReturnOperators)
}

print.INApestMLUPathogenAnalytical <- function(x,...){cat("INApest MLU pathogen analytical result\n");cat("  Pathogen model:",x$PathogenModel,"\n");cat("  Exact:",isTRUE(x$Exact),"\n");if(!is.null(x$Growth$IntrinsicRarePathogenMultiplier))cat("  Rare-pathogen multiplier:",round(x$Growth$IntrinsicRarePathogenMultiplier,6),"\n");invisible(x)}

INApestAnalytical_pre_mlu_pathogen <- INApestAnalytical
INApestAnalytical <- function(...) {
  args<-list(...);Model<-if(!is.null(args$Model))as.character(args$Model)[1L]else"INApest";Pathogen<-args$Pathogen
  if(!identical(Model,"INApestMetaMultipleLandUse")||is.null(Pathogen))return(do.call(INApestAnalytical_pre_mlu_pathogen,args))
  if(is.null(args$HostPopulation)){if(is.null(args$InitialState))stop("MLU pathogen analytical calls require HostPopulation or InitialState");args$HostPopulation<-args$InitialState}
  keep<-c("Ntimesteps","HostPopulation","Pathogen","HostAssumption","K","Survival","RecruitToCapacity","HostDetectionProb","InitialInfo","ManageProb","MortalityProb","InformationAcquisition","InfoPersistenceSteps","InfoRetentionProb","NodeContact","LandUseMixing","OutsideNodes","OutsideUnits","Exact","ExactMaxStates","ReturnOperators")
  do.call(INApestMLUPathogenAnalytical,args[intersect(names(args),keep)])
}

###############################################################################
### Exact tiny-network MLU SIS + shared node information + managed escape
###############################################################################

.ina_mlu_node_vector <- function(x,n_nodes,name,integer=FALSE){z<-as.numeric(x);if(length(z)==1L)z<-rep(z,n_nodes);if(length(z)!=n_nodes||any(!is.finite(z)))stop(name," must be scalar or length nodes");if(integer&&any(z!=floor(z)))stop(name," must contain whole numbers");z}

INApestMLUSISExactNetworkInformation <- function(Ntimesteps,K,InitialHostPopulation,InitialInfected,InitialInfo=0,
    Survival=1,ManageProb=0,MortalityProb=0,HostDetectionProb=0,PathogenDetectionProb=0,
    InformationAcquisition=c("host","pathogen","both"),InfoPersistenceSteps=NA,InfoRetentionProb=1,
    Beta=0,RecoveryProb=0,PathogenMortalityProb=0,IntroductionProb=0,IntroductionNumber=1,
    ContactMatrix=NULL,NodeContact=NULL,LandUseMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,
    OutsideNodes=integer(0),MaxStates=5000L,ReturnOperator=FALSE){
  InformationAcquisition<-match.arg(InformationAcquisition);Transmission<-match.arg(Transmission)
  K<-.ina_mlu_check_matrix(K,name="K",integer=TRUE);N0<-.ina_mlu_check_matrix(InitialHostPopulation,nrow(K),ncol(K),"InitialHostPopulation",TRUE);I0<-.ina_mlu_check_matrix(InitialInfected,nrow(K),ncol(K),"InitialInfected",TRUE);if(any(N0>K)||any(I0>N0))stop("Invalid initial network state")
  n<-nrow(K);L<-ncol(K);n_units<-n*L;H0<-as.integer(.ina_mlu_node_vector(InitialInfo,n,"InitialInfo",TRUE));if(any(!H0%in%c(0L,1L)))stop("InitialInfo must be 0/1")
  ipers<-as.numeric(InfoPersistenceSteps);if(length(ipers)==1L)ipers<-rep(ipers,n);if(length(ipers)!=n||any(!is.na(ipers)&(ipers<0|ipers!=floor(ipers))))stop("InfoPersistenceSteps must be scalar or length nodes")
  iret<-.ina_mlu_node_vector(InfoRetentionProb,n,"InfoRetentionProb");if(any(iret<0|iret>1))stop("InfoRetentionProb must be in [0,1]")
  to_mat<-function(x,nm,int=FALSE){matrix(.ina_mlu_resolve_unit(x,n,L,nm,int),n,L)}
  surv<-to_mat(Survival,"Survival");mng<-to_mat(ManageProb,"ManageProb");mort<-to_mat(MortalityProb,"MortalityProb");hd<-to_mat(HostDetectionProb,"HostDetectionProb");pd<-to_mat(PathogenDetectionProb,"PathogenDetectionProb");beta<-as.numeric(.ina_mlu_resolve_unit(Beta,n,L,"Beta"));rec<-as.numeric(.ina_mlu_resolve_unit(RecoveryProb,n,L,"RecoveryProb"));pmort<-as.numeric(.ina_mlu_resolve_unit(PathogenMortalityProb,n,L,"PathogenMortalityProb"));intro<-as.numeric(.ina_mlu_resolve_unit(IntroductionProb,n,L,"IntroductionProb"));inum<-as.integer(.ina_mlu_resolve_unit(IntroductionNumber,n,L,"IntroductionNumber",TRUE));ds<-as.numeric(.ina_mlu_resolve_unit(DensityScale,n,L,"DensityScale"))
  if(any(c(surv,mng,mort,hd,pd,rec,pmort,intro)<0)||any(c(surv,mng,mort,hd,pd,rec,pmort,intro)>1)||any(rec+pmort>1+1e-12)||any(beta<0)||any(ds<=0))stop("Invalid network parameters")
  C<-.ina_mlu_contact(n,L,ContactMatrix,NodeContact,LandUseMixing);host_acq<-InformationAcquisition%in%c("host","both");path_acq<-InformationAcquisition%in%c("pathogen","both")
  local<-lapply(seq_len(n),function(j).ina_mlu_one_node_info_states("SIS",as.integer(K[j,]),ipers[j],MaxStates));dims<-vapply(local,nrow,integer(1));nst<-prod(dims);if(nst>MaxStates)stop("Exact MLU network information state space has ",nst," states; reduce capacities/nodes or increase MaxStates deliberately")
  grids<-as.matrix(expand.grid(lapply(dims,seq_len),KEEP.OUT.ATTRS=FALSE));storage.mode(grids)<-"integer";gkey<-apply(grids,1,paste,collapse=":");gidx<-setNames(seq_len(nst),gkey);local_keys<-lapply(local,function(z)apply(z,1,paste,collapse=":"));local_idx<-lapply(local,function(z)setNames(seq_len(nrow(z)),apply(z,1,paste,collapse=":")))
  T<-matrix(0,nst,nst)
  for(row in seq_len(nst)){
    preL<-vector("list",n)
    for(j in seq_len(n))preL[[j]]<-.ina_mlu_one_node_info_prebranches(local[[j]][grids[row,j],,drop=FALSE],"SIS",L,surv[j,],mng[j,],mort[j,],host_acq,ipers[j],iret[j])
    pg<-as.matrix(expand.grid(lapply(preL,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(pg)<-"integer"
    for(pb in seq_len(nrow(pg))){Npre<-Ipre<-matrix(0,n,L);Hmid<-Amid<-integer(n);pp<-1
      for(j in seq_len(n)){z<-preL[[j]][pg[pb,j],];for(lu in seq_len(L)){Ipre[j,lu]<-z[[paste0("I_lu",lu)]];Npre[j,lu]<-Ipre[j,lu]+z[[paste0("S_lu",lu)]]};Hmid[j]<-z$H;Amid[j]<-z$A;pp<-pp*z$prob};if(pp==0)next
      Iv<-as.numeric(Ipre);Nv<-as.numeric(Npre);pressure<-as.numeric(crossprod(Iv,C));den<-as.numeric(crossprod(Nv,C));foi<-if(Transmission=="frequency")beta*ifelse(den>0,pressure/den,0)else beta*pressure/ds;pinf<-.ina_meta_sis_clip01(-expm1(-pmax(0,foi)))
      dl<-vector("list",n_units);for(u in seq_len(n_units))dl[[u]]<-as.data.frame(.ina_meta_sis_disease_outcomes_pinf(as.integer(Nv[u]),as.integer(Iv[u]),pinf[u],rec[u],pmort[u],intro[u],inum[u]));dg<-as.matrix(expand.grid(lapply(dl,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(dg)<-"integer"
      for(dd in seq_len(nrow(dg))){N1<-I1<-numeric(n_units);pdis<-pp;for(u in seq_len(n_units)){z<-dl[[u]][dg[dd,u],];N1[u]<-z$N;I1[u]<-z$I;pdis<-pdis*z$prob};if(pdis==0)next;N1m<-matrix(N1,n,L);I1m<-matrix(I1,n,L)
        db<-vector("list",n);for(j in seq_len(n)){qh<-if(host_acq)1-prod((1-hd[j,])^N1m[j,])else 0;qp<-if(path_acq)1-prod((1-pd[j,])^I1m[j,])else 0;qev<-1-(1-qh)*(1-qp);if(qev>0&&qev<1)db[[j]]<-rbind(c(H=1,A=if(is.na(ipers[j]))-2 else 0,p=qev),c(H=Hmid[j],A=Amid[j],p=1-qev))else if(qev>=1)db[[j]]<-matrix(c(1,if(is.na(ipers[j]))-2 else 0,1),1,dimnames=list(NULL,c("H","A","p")))else db[[j]]<-matrix(c(Hmid[j],Amid[j],1),1,dimnames=list(NULL,c("H","A","p")))}
        bg<-as.matrix(expand.grid(lapply(db,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(bg)<-"integer"
        for(bb in seq_len(nrow(bg))){li<-integer(n);pfin<-pdis;for(j in seq_len(n)){q<-db[[j]][bg[bb,j],];bio<-integer(2L*L);for(lu in seq_len(L)){u<-j+(lu-1L)*n;bio[2L*lu-1L]<-as.integer(N1[u]-I1[u]);bio[2L*lu]<-as.integer(I1[u])};ky<-paste(c(bio,as.integer(q["H"]),as.integer(q["A"])),collapse=":");li[j]<-local_idx[[j]][[ky]];pfin<-pfin*unname(q["p"])};gi<-gidx[[paste(li,collapse=":")]];T[row,gi]<-T[row,gi]+pfin}
      }
    }
  }
  er<-max(abs(rowSums(T)-1));if(er>1e-10)stop("MLU network information operator row error ",er)
  # Initial distribution includes direct local evidence before timestep 1.
  local_init<-vector("list",n)
  for(j in seq_len(n)){qh<-if(host_acq)1-prod((1-hd[j,])^N0[j,])else 0;qp<-if(path_acq)1-prod((1-pd[j,])^I0[j,])else 0;qev<-1-(1-qh)*(1-qp);bi<-integer(2L*L);for(lu in seq_len(L)){bi[2L*lu-1L]<-N0[j,lu]-I0[j,lu];bi[2L*lu]<-I0[j,lu]};lst<-list();kk<-0L;if(qev>0){kk<-kk+1L;lst[[kk]]<-c(idx=local_idx[[j]][[paste(c(bi,1,if(is.na(ipers[j]))-2 else 0),collapse=":")]],p=qev)};if(qev<1){A<-if(H0[j]==1&&!is.na(ipers[j]))-1 else -2;kk<-kk+1L;lst[[kk]]<-c(idx=local_idx[[j]][[paste(c(bi,H0[j],A),collapse=":")]],p=1-qev)};local_init[[j]]<-do.call(rbind,lst)}
  ig<-as.matrix(expand.grid(lapply(local_init,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(ig)<-"integer";dist<-numeric(nst);for(rr in seq_len(nrow(ig))){li<-integer(n);p<-1;for(j in seq_len(n)){z<-local_init[[j]][ig[rr,j],];li[j]<-as.integer(z["idx"]);p<-p*z["p"]};dist[gidx[[paste(li,collapse=":")]]]<-dist[gidx[[paste(li,collapse=":")]]]+p}
  EI<-array(0,dim=c(n,L,Ntimesteps+1L));EN<-array(0,dim=c(n,L,Ntimesteps+1L));PH<-matrix(0,Ntimesteps+1L,n);PM<-array(0,dim=c(n,L,Ntimesteps));full<-matrix(0,Ntimesteps+1L,nst);full[1,]<-dist
  summarize<-function(d,tt){for(j in seq_len(n)){h<-numeric(nst);for(r in seq_len(nst)){lr<-local[[j]][grids[r,j],];h[r]<-lr["H"];for(lu in seq_len(L)){EI[j,lu,tt]<<-EI[j,lu,tt]+d[r]*lr[paste0("I_lu",lu)];EN[j,lu,tt]<<-EN[j,lu,tt]+d[r]*(lr[paste0("S_lu",lu)]+lr[paste0("I_lu",lu)])}};PH[tt,j]<<-sum(d*h)}}
  summarize(dist,1L);for(tt in seq_len(Ntimesteps)){for(j in seq_len(n))for(lu in seq_len(L))PM[j,lu,tt]<-mng[j,lu]*PH[tt,j];dist<-as.numeric(dist%*%T);full[tt+1,]<-dist;summarize(dist,tt+1L)}
  outside<-unique(as.integer(OutsideNodes));if(any(!outside%in%seq_len(n)))stop("OutsideNodes contains invalid node indices");escape<-NULL
  if(length(outside)){
    safe<-logical(nst);for(r in seq_len(nst)){ok<-TRUE;for(j in outside){lr<-local[[j]][grids[r,j],];if(sum(lr[grep("^I_lu",names(lr))])>0){ok<-FALSE;break}};safe[r]<-ok};sd<-full[1,safe];Q<-T[safe,safe,drop=FALSE];esc<-numeric(Ntimesteps+1L);esc[1]<-1-sum(sd);if(Ntimesteps>0)for(tt in seq_len(Ntimesteps)){sd<-as.numeric(sd%*%Q);esc[tt+1]<-1-sum(sd)};escape<-list(OutsideNodes=outside,ProbabilityByTimestep=esc,ProbabilityByHorizon=tail(esc,1))
  }
  out<-list(Model="INApestMetaMultipleLandUse",PathogenModel="SIS",Exact=TRUE,InformationAcquisition=InformationAcquisition,ExpectedHost=EN,ExpectedInfected=EI,InformationProbabilityByNode=PH,ManagingProbability=PM,Escape=escape,StateCount=nst,StateDistribution=full,Diagnostics=c("Exact tiny-network MLU SIS operator with land-use host/pathogen states and node-level shared information.","First-passage escape is computed with a safe-state sub-operator, so later clearance outside containment does not erase an earlier escape.","This branch is intentionally state-space limited and serves as validation truth; use stochastic MLU for realistic managed landscapes."));if(ReturnOperator)out$Operator<-T;class(out)<-c("INApestMLUSISExactNetworkInformation","list");out
}

# Replace only the MLU router with the final network-information-aware version.
INApestMLUPathogenAnalytical_v1 <- INApestMLUPathogenAnalytical
INApestMLUPathogenAnalytical <- function(Ntimesteps,HostPopulation,Pathogen,
    HostAssumption=c("fixed","event_contract_turnover"),K=HostPopulation,Survival=1,RecruitToCapacity=FALSE,
    HostDetectionProb=0,InitialInfo=0,ManageProb=0,MortalityProb=0,InformationAcquisition=NULL,
    InfoPersistenceSteps=NA,InfoRetentionProb=1,NodeContact=NULL,LandUseMixing=NULL,
    OutsideNodes=integer(0),OutsideUnits=integer(0),Exact=TRUE,ExactMaxStates=50000L,ReturnOperators=FALSE){
  N<-.ina_mlu_check_matrix(HostPopulation,name="HostPopulation",integer=TRUE);n<-nrow(N);L<-ncol(N);Model<-as.character(Pathogen$Model)[1L];acq<-InformationAcquisition;if(is.null(acq))acq<-if(isTRUE(Pathogen$DetectionTriggersInfo))"both"else"host";info_requested<-any(as.numeric(ManageProb)!=0,na.rm=TRUE)||any(as.numeric(InitialInfo)!=0,na.rm=TRUE)||any(as.numeric(HostDetectionProb)!=0,na.rm=TRUE)||any(as.numeric(Pathogen$DetectionProb)!=0,na.rm=TRUE)||any(!is.na(InfoPersistenceSteps))||any(as.numeric(InfoRetentionProb)!=1,na.rm=TRUE)
  if(isTRUE(Exact)&&match.arg(HostAssumption)=="fixed"&&info_requested&&Model=="SIS"&&n>1L){Kmat<-if(is.matrix(K))K else matrix(.ina_mlu_resolve_unit(K,n,L,"K",TRUE),n,L);I0<-.ina_mlu_initial_matrix(Pathogen$InitialInfected,n,L,"InitialInfected");C<-.ina_mlu_contact(n,L,ContactMatrix=Pathogen$ContactMatrix,NodeContact=NodeContact,LandUseMixing=LandUseMixing);return(INApestMLUSISExactNetworkInformation(Ntimesteps,Kmat,N,I0,InitialInfo,Survival,ManageProb,MortalityProb,HostDetectionProb,Pathogen$DetectionProb,acq,InfoPersistenceSteps,InfoRetentionProb,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$PathogenMortalityProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,ContactMatrix=C,Transmission=Pathogen$Transmission,DensityScale=Pathogen$DensityScale,OutsideNodes=OutsideNodes,MaxStates=ExactMaxStates,ReturnOperator=ReturnOperators))}
  INApestMLUPathogenAnalytical_v1(Ntimesteps,HostPopulation,Pathogen,HostAssumption,K,Survival,RecruitToCapacity,HostDetectionProb,InitialInfo,ManageProb,MortalityProb,InformationAcquisition,InfoPersistenceSteps,InfoRetentionProb,NodeContact,LandUseMixing,OutsideNodes,OutsideUnits,Exact,ExactMaxStates,ReturnOperators)
}


###############################################################################
### INApestMetaTransitionMatrix pathogen analytical methods
### Core mathematics: demographic carrier transport -> pathogen state process
###############################################################################

.ina_tm_ns_matrix <- function(x,n,S,name,prob=FALSE){
  if(is.matrix(x)){
    if(!all(dim(x)==c(n,S))) stop(name," matrix must be nodes x stages")
    z<-x
  } else {
    x<-as.numeric(x)
    if(length(x)==1L) z<-matrix(x,n,S)
    else if(length(x)==S && S!=n) z<-matrix(rep(x,each=n),n,S)
    else if(length(x)==n && n!=S) z<-matrix(rep(x,S),n,S)
    else if(length(x)==n*S) z<-matrix(x,n,S)
    else stop(name," must be scalar, stage vector, node vector, or nodes x stages")
  }
  if(any(!is.finite(z))) stop(name," must be finite")
  if(prob && any(z<0|z>1)) stop(name," must lie in [0,1]")
  z
}

.ina_tm_transition_movement <- function(x,S,n,name){
  if(S<=1L) return(list())
  if(is.null(x)) return(rep(list(NULL),S-1L))
  z<-if(is.list(x))x else rep(list(x),S-1L)
  if(length(z)!=S-1L) stop(name," must be a matrix or list of length Nstages-1")
  for(k in seq_along(z)) if(!is.null(z[[k]])){
    z[[k]]<-as.matrix(z[[k]])
    if(!all(dim(z[[k]])==c(n,n))) stop(name," matrices must be nodes x nodes")
    if(any(!is.finite(z[[k]]))||any(z[[k]]<0)||any(rowSums(z[[k]])>1+1e-12))
      stop(name," rows must contain finite non-negative probabilities summing to <=1")
  }
  z
}

INApestTransitionPathogenCarrierOperator <- function(
    Transition,Nstages,n_nodes,
    ManageProb=0,MortalityProb=0,
    TransitionSDDprob=NULL,TransitionLDDprob=NULL,TransitionLDDrate=0){
  S<-as.integer(Nstages);n<-as.integer(n_nodes)
  if(S<1L||n<1L)stop("Nstages and n_nodes must be positive")
  A<-.transition_list(Transition,n,S)
  a<-.ina_recycle(ManageProb,n,"ManageProb")
  M<-.ina_tm_ns_matrix(MortalityProb,n,S,"MortalityProb",TRUE)
  if(any(a<0|a>1))stop("ManageProb must lie in [0,1]")
  Ps<-.ina_tm_transition_movement(TransitionSDDprob,S,n,"TransitionSDDprob")
  Pl<-.ina_tm_transition_movement(TransitionLDDprob,S,n,"TransitionLDDprob")
  rr<-rep_len(as.numeric(TransitionLDDrate),max(1L,S-1L))
  if(any(!is.finite(rr))||any(rr<0|rr>1))stop("TransitionLDDrate must lie in [0,1]")
  U<-n*S;idx<-function(i,s)(i-1L)*S+s
  H<-matrix(0,U,U)
  # Management occurs before demographic progression. Integrating over the
  # Bernoulli adoption event gives the one-carrier survival factor below.
  q<-1-M*sweep(matrix(1,n,S),1,a,`*`)
  for(i in seq_len(n))for(k in seq_len(S)){
    src<-idx(i,k);Ai<-A[[i]];qq<-q[i,k]
    if(k<S){
      H[idx(i,k),src]<-H[idx(i,k),src]+qq*Ai[k,k]
      tr<-Ai[k+1L,k]
      if(tr>0){
        if(is.null(Ps[[k]])&&is.null(Pl[[k]])) P<-diag(n)
        else if(!is.null(Ps[[k]])&&!is.null(Pl[[k]])) P<-(1-rr[k])*Ps[[k]]+rr[k]*Pl[[k]]
        else if(!is.null(Ps[[k]])) P<-Ps[[k]] else P<-Pl[[k]]
        for(j in seq_len(n))if(P[i,j]>0)
          H[idx(j,k+1L),src]<-H[idx(j,k+1L),src]+qq*tr*P[i,j]
      }
    } else H[idx(i,S),src]<-H[idx(i,S),src]+qq*Ai[S,S]
  }
  map<-expand.grid(stage=seq_len(S),node=seq_len(n),KEEP.OUT.ATTRS=FALSE)
  # expand.grid stage-fast order already matches (node-1)*S+stage.
  labels<-paste0("n",map$node,"_s",map$stage)
  dimnames(H)<-list(target=labels,source=labels)
  list(Operator=H,UnitMap=map,ColumnSurvival=rowSums(t(H)),
       Diagnostics=c(
         "Carrier operator tracks only existing hosts: demographic fecundity is excluded because all offspring enter pathogen state S.",
         "Management mortality is applied before demographic stage survival/progression, matching INApestMetaTransitionMatrix.",
         "Stage-transition movement is included exactly in the uncrowded limit; residual row probability is export/loss.",
         "Capacity blocking and BlockedTransitionMortality are intentionally excluded from this rare-carrier operator and require a crowded-host correction."))
}

INApestTransitionPathogenDiseaseLinearOperator <- function(
    Model=c("SIS","SIR","SEIR"),HostPopulationDiseaseStep,
    Beta,RecoveryProb,ProgressionProb=1,PathogenMortalityProb=0,
    ContactMatrix=NULL,StageMixing=NULL,Transmission=c("frequency","density"),DensityScale=1){
  Model<-match.arg(Model);Transmission<-match.arg(Transmission)
  Nmat<-as.matrix(HostPopulationDiseaseStep);n<-nrow(Nmat);S<-ncol(Nmat);U<-n*S
  if(any(!is.finite(Nmat))||any(Nmat<0))stop("HostPopulationDiseaseStep must be finite non-negative nodes x stages")
  # Stochastic engine is node-major, stage-fast.
  N<-as.numeric(t(Nmat))
  beta<-as.numeric(t(.ina_tm_ns_matrix(Beta,n,S,"Beta")))
  rec<-as.numeric(t(.ina_tm_ns_matrix(RecoveryProb,n,S,"RecoveryProb",TRUE)))
  prog<-as.numeric(t(.ina_tm_ns_matrix(ProgressionProb,n,S,"ProgressionProb",TRUE)))
  mort<-as.numeric(t(.ina_tm_ns_matrix(PathogenMortalityProb,n,S,"PathogenMortalityProb",TRUE)))
  ds<-as.numeric(t(.ina_tm_ns_matrix(DensityScale,n,S,"DensityScale")))
  if(any(beta<0)||any(ds<=0)||any(rec+mort>1+1e-12))stop("Invalid pathogen parameters")
  Cn<-if(is.null(ContactMatrix))diag(n)else as.matrix(ContactMatrix)
  if(!all(dim(Cn)==c(n,n))||any(!is.finite(Cn))||any(Cn<0))stop("ContactMatrix must be finite non-negative nodes x nodes")
  Sm<-if(is.null(StageMixing))matrix(1,S,S)else as.matrix(StageMixing)
  if(!all(dim(Sm)==c(S,S))||any(!is.finite(Sm))||any(Sm<0))stop("StageMixing must be finite non-negative stages x stages")
  C<-kronecker(Cn,Sm) # source x target, node-major stage-fast
  B<-matrix(0,U,U)
  if(Transmission=="frequency"){
    den<-as.numeric(crossprod(N,C))
    for(w in seq_len(U))if(N[w]>0&&den[w]>0&&beta[w]>0)B[w,]<-N[w]*beta[w]*C[,w]/den[w]
  }else for(w in seq_len(U))if(N[w]>0&&beta[w]>0)B[w,]<-N[w]*beta[w]*C[,w]/ds[w]
  stayI<-diag(1-rec-mort,U)
  if(Model%in%c("SIS","SIR")){
    G<-stayI+B
    note<-if(Model=="SIR")"SIR has the same active-I disease-step linearisation as SIS at the pathogen-free state; R is inactive to first order."else"SIS disease-step linearisation."
  } else {
    P<-diag(prog,U);stayE<-diag(1-prog,U)
    G<-rbind(cbind(stayE,B),cbind(P,stayI))
    note<-"SEIR preserves the discrete latent delay: infections generated this step enter E and cannot progress until a later pathogen step."
  }
  list(Operator=G,TransmissionBlock=B,CombinedContact=C,HostPopulationDiseaseStep=Nmat,
       Model=Model,Transmission=Transmission,Diagnostics=c(
         "Rows are recipient active pathogen states; columns are source active pathogen states.",
         "Frequency-transmission denominators use the post-demographic disease-step host distribution, matching the stochastic engine.",note))
}

INApestTransitionPathogenGrowthOperator <- function(
    Model=c("SIS","SIR","SEIR"),Transition,Nstages,HostPopulationDiseaseStep,
    Beta,RecoveryProb,ProgressionProb=1,PathogenMortalityProb=0,
    ContactMatrix=NULL,StageMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,
    ManageProb=0,MortalityProb=0,
    TransitionSDDprob=NULL,TransitionLDDprob=NULL,TransitionLDDrate=0){
  Model<-match.arg(Model);Transmission<-match.arg(Transmission)
  Nmat<-as.matrix(HostPopulationDiseaseStep);n<-nrow(Nmat);S<-ncol(Nmat)
  if(S!=Nstages)stop("HostPopulationDiseaseStep columns must equal Nstages")
  car<-INApestTransitionPathogenCarrierOperator(Transition,S,n,ManageProb,MortalityProb,
                                                 TransitionSDDprob,TransitionLDDprob,TransitionLDDrate)
  dis<-INApestTransitionPathogenDiseaseLinearOperator(Model,Nmat,Beta,RecoveryProb,ProgressionProb,
                                                       PathogenMortalityProb,ContactMatrix,StageMixing,
                                                       Transmission,DensityScale)
  H<-car$Operator;U<-n*S
  if(Model%in%c("SIS","SIR"))G<-dis$Operator%*%H
  else {
    H2<-rbind(cbind(H,matrix(0,U,U)),cbind(matrix(0,U,U),H))
    G<-dis$Operator%*%H2
  }
  ev<-eigen(G,only.values=TRUE)$values;lambda<-max(Mod(ev))
  list(Model="INApestMetaTransitionMatrix",PathogenModel=Model,Operator=G,
       CarrierOperator=H,DiseaseOperator=dis$Operator,TransmissionBlock=dis$TransmissionBlock,
       IntrinsicRarePathogenMultiplier=lambda,
       Classification=if(lambda>1+1e-12)"growing"else if(lambda<1-1e-12)"declining"else"threshold",
       HostPopulationDiseaseStep=Nmat,UnitMap=car$UnitMap,
       Diagnostics=c(car$Diagnostics,dis$Diagnostics,
         "Full-timestep active-pathogen operator is disease-step operator composed after demographic carrier transport.",
         "This is a rare-pathogen/uncrowded carrier result. Demographic births affect the host background and transmission denominators but do not directly create infected offspring."))
}

INApestTransitionPathogenPeriodicGrowth <- function(Operators){
  if(!is.list(Operators)||!length(Operators))stop("Operators must be a non-empty list of same-sized square matrices")
  mats<-lapply(Operators,function(x)if(is.list(x)&&!is.null(x$Operator))x$Operator else as.matrix(x))
  d<-dim(mats[[1L]]);if(d[1]!=d[2]||any(vapply(mats,function(M)!all(dim(M)==d),logical(1))))stop("All operators must have the same square dimension")
  P<-diag(d[1]);for(G in mats)P<-G%*%P
  rho<-max(Mod(eigen(P,only.values=TRUE)$values));g<-rho^(1/length(mats))
  list(CycleOperator=P,CycleMultiplier=rho,PerTimestepMultiplier=g,
       Classification=if(g>1+1e-12)"growing"else if(g<1-1e-12)"declining"else"threshold",
       Diagnostics="For a periodic environment the correct threshold is the spectral radius of the ordered product, not the mean of timestep-specific eigenvalues.")
}

.ina_tm_single_carrier_pinf <- function(v,N,beta,C,Transmission,ds){
  U<-length(N);out<-numeric(U)
  if(Transmission=="frequency"){
    den<-as.numeric(crossprod(N,C))
    for(w in seq_len(U))if(den[w]>0&&beta[w]>0&&C[v,w]>0)out[w]<--expm1(-beta[w]*C[v,w]/den[w])
  } else for(w in seq_len(U))if(beta[w]>0&&C[v,w]>0)out[w]<--expm1(-beta[w]*C[v,w]/ds[w])
  pmin(1,pmax(0,out))
}

INApestTransitionPathogenBranching <- function(
    Model=c("SIS","SIR","SEIR"),timesteps,Transition,Nstages,HostPopulationDiseaseStep,
    Beta,RecoveryProb,ProgressionProb=1,PathogenMortalityProb=0,
    ContactMatrix=NULL,StageMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,
    ManageProb=0,MortalityProb=0,TransitionSDDprob=NULL,TransitionLDDprob=NULL,TransitionLDDrate=0,
    OutsideNodes=integer(0),OutsideUnits=integer(0)){
  Model<-match.arg(Model);Transmission<-match.arg(Transmission);T<-as.integer(timesteps)
  Nmat<-as.matrix(HostPopulationDiseaseStep);n<-nrow(Nmat);S<-ncol(Nmat);U<-n*S;N<-as.numeric(t(Nmat))
  car<-INApestTransitionPathogenCarrierOperator(Transition,S,n,ManageProb,MortalityProb,
                                                 TransitionSDDprob,TransitionLDDprob,TransitionLDDrate)
  H<-car$Operator
  beta<-as.numeric(t(.ina_tm_ns_matrix(Beta,n,S,"Beta")));rec<-as.numeric(t(.ina_tm_ns_matrix(RecoveryProb,n,S,"RecoveryProb",TRUE)))
  prog<-as.numeric(t(.ina_tm_ns_matrix(ProgressionProb,n,S,"ProgressionProb",TRUE)));mort<-as.numeric(t(.ina_tm_ns_matrix(PathogenMortalityProb,n,S,"PathogenMortalityProb",TRUE)))
  ds<-as.numeric(t(.ina_tm_ns_matrix(DensityScale,n,S,"DensityScale")))
  Cn<-if(is.null(ContactMatrix))diag(n)else as.matrix(ContactMatrix);Sm<-if(is.null(StageMixing))matrix(1,S,S)else as.matrix(StageMixing);C<-kronecker(Cn,Sm)
  outside<-unique(as.integer(OutsideUnits));if(length(OutsideNodes))for(i in OutsideNodes)outside<-c(outside,(i-1L)*S+seq_len(S));outside<-unique(outside)
  if(any(!outside%in%seq_len(U)))stop("Outside node/unit index invalid")
  integerN<-all(abs(N-round(N))<1e-10)
  transmission_pgf<-function(pinf,Ssus,z){
    if(integerN)prod((1-pinf+pinf*z)^as.integer(round(Ssus)))
    else exp(sum(Ssus*pinf*(z-1)))
  }
  # q = active-lineage extinction-by-horizon; h = no first-passage outside active state by horizon.
  if(Model%in%c("SIS","SIR")){q<-rep(0,U);h<-rep(1,U);if(length(outside))h[outside]<-0;qh<-hh<-matrix(NA,U,T)
    stepfun<-function(z,noescape=FALSE){zn<-numeric(U);for(u in seq_len(U)){
      val<-max(0,1-sum(H[,u]));for(v in which(H[,u]>0)){
        pinf<-.ina_tm_single_carrier_pinf(v,N,beta,C,Transmission,ds);Ssus<-pmax(0,N-as.numeric(seq_len(U)==v));trans<-transmission_pgf(pinf,Ssus,z);stay<-1-rec[v]-mort[v];carrier<-(1-stay)+stay*z[v];val<-val+H[v,u]*carrier*trans}
      zn[u]<-val};pmin(1,pmax(0,zn))}
    for(tt in seq_len(T)){q<-stepfun(q);h<-stepfun(h,TRUE);if(length(outside))h[outside]<-0;qh[,tt]<-q;hh[,tt]<-h}
    meanop<-matrix(0,U,U);for(u in seq_len(U))for(v in which(H[,u]>0)){
      pinf<-.ina_tm_single_carrier_pinf(v,N,beta,C,Transmission,ds);Ssus<-pmax(0,N-as.numeric(seq_len(U)==v));meanop[,u]<-meanop[,u]+H[v,u]*(Ssus*pinf);meanop[v,u]<-meanop[v,u]+H[v,u]*(1-rec[v]-mort[v])}
  } else {
    q<-rep(0,2*U);h<-rep(1,2*U);if(length(outside))h[c(outside,U+outside)]<-0;qh<-hh<-matrix(NA,2*U,T)
    stepfun<-function(z){zn<-numeric(2*U)
      # E source: demographic carrier then E stay/progress; no transmission.
      for(u in seq_len(U)){val<-max(0,1-sum(H[,u]));for(v in which(H[,u]>0))val<-val+H[v,u]*((1-prog[v])*z[v]+prog[v]*z[U+v]);zn[u]<-val}
      # I source: demographic carrier, transmission creates E, carrier may stay I.
      for(u in seq_len(U)){val<-max(0,1-sum(H[,u]));for(v in which(H[,u]>0)){
        pinf<-.ina_tm_single_carrier_pinf(v,N,beta,C,Transmission,ds);Ssus<-pmax(0,N-as.numeric(seq_len(U)==v));trans<-transmission_pgf(pinf,Ssus,z[seq_len(U)]);stay<-1-rec[v]-mort[v];val<-val+H[v,u]*((1-stay)+stay*z[U+v])*trans};zn[U+u]<-val}
      pmin(1,pmax(0,zn))}
    for(tt in seq_len(T)){q<-stepfun(q);h<-stepfun(h);if(length(outside))h[c(outside,U+outside)]<-0;qh[,tt]<-q;hh[,tt]<-h}
    meanop<-matrix(0,2*U,2*U)
    for(u in seq_len(U))for(v in which(H[,u]>0)){
      meanop[v,u]<-meanop[v,u]+H[v,u]*(1-prog[v]);meanop[U+v,u]<-meanop[U+v,u]+H[v,u]*prog[v]
      pinf<-.ina_tm_single_carrier_pinf(v,N,beta,C,Transmission,ds);Ssus<-pmax(0,N-as.numeric(seq_len(U)==v));meanop[seq_len(U),U+u]<-meanop[seq_len(U),U+u]+H[v,u]*(Ssus*pinf);meanop[U+v,U+u]<-meanop[U+v,U+u]+H[v,u]*(1-rec[v]-mort[v])
    }
  }
  rho<-max(Mod(eigen(meanop,only.values=TRUE)$values))
  inside<-setdiff(seq_len(U),outside)
  if(Model=="SEIR") inside_types<-c(inside,U+inside) else inside_types<-inside
  esc<-1-h
  list(Model="INApestMetaTransitionMatrix",PathogenModel=Model,Exact=FALSE,
       OneCarrierMeanOperator=meanop,OneCarrierBranchingMultiplier=rho,
       ExtinctionByHorizon=q,EscapeByHorizon=esc,
       InsideStartingTypes=inside_types,
       MaxEscapeByHorizonFromInside=if(length(inside_types))max(esc[inside_types])else NA_real_,
       ExtinctionHistory=qh,EscapeHistory=1-hh,OutsideUnits=outside,UnitMap=car$UnitMap,
       Diagnostics=c(
         "Multitype branching approximation is exact for descendants of a single lineage until different pathogen lineages compete for the same susceptible hosts.",
         "The one-carrier offspring distribution uses the stochastic engine's exact exponential infection probability, not only its infinitesimal Jacobian.",
         if(integerN) "Integer host backgrounds use the exact finite-count binomial offspring PGF." else "Fractional analytical host backgrounds use a mean-matched Poisson offspring PGF; this avoids treating a fractional host count as a binomial exponent.",
         "Escape is defined at timestep boundaries by an active pathogen state outside containment; transient carrier movement followed by recovery/death within the same timestep is not counted unless it leaves an active descendant outside."))
}

INApestTransitionHostMeanOperator <- function(
    Transition,Nstages,SDDprob,LDDprob=0,LDDrate=0,
    EnvEstabProb=1,PropaguleEstablishment=1,
    ManageProb=0,MortalityProb=0,SpreadReduction=0,FecundityReduction=0,
    DispersalDensityFactor=0,K=1,SeedbankK=1,
    TransitionSDDprob=NULL,TransitionLDDprob=NULL,TransitionLDDrate=0){
  SDD<-as.matrix(SDDprob);n<-nrow(SDD);S<-as.integer(Nstages)
  z<-transition_components(Transition,S,SDD,LDDprob,LDDrate,EnvEstabProb,
                           PropaguleEstablishment,ManageProb,MortalityProb,
                           SpreadReduction,DispersalDensityFactor,K,SeedbankK,
                           FecundityReduction)
  car<-INApestTransitionPathogenCarrierOperator(Transition,S,n,ManageProb,MortalityProb,
                                                 TransitionSDDprob,TransitionLDDprob,TransitionLDDrate)
  H<-car$Operator;R<-matrix(0,n*S,n*S);idx<-function(i,s)(i-1L)*S+s
  F<-if(!is.null(z$fecundity_reduction))z$fecundity_reduction else matrix(0,n,S)
  for(i in seq_len(n)){
    Ai<-z$A[[i]];a<-z$adoption[i];g<-z$spread_reduction[i]
    for(k in 2:S){
      fec<-Ai[1,k];if(fec<=0)next
      m<-z$mortality[i,k];fr<-F[i,k]
      for(j in seq_len(n)){
        nat<-(1-z$LDDrate)*z$sdd_enabled[i]*z$SDD[i,j]
        hum0<-z$LDDrate*z$LDD[i,j]
        hum1<-z$LDDrate*(1-g)*z$LDD[i,j]
        w<-fec*z$recruit_success[j]*((1-a)*(nat+hum0)+a*(1-m)*(1-fr)*(nat+hum1))
        R[idx(j,1L),idx(i,k)]<-R[idx(j,1L),idx(i,k)]+w
      }
    }
  }
  G<-H+R
  dimnames(R)<-dimnames(H);dimnames(G)<-dimnames(H)
  list(Operator=G,CarrierOperator=H,RecruitmentOperator=R,UnitMap=car$UnitMap,
       Diagnostics=c(car$Diagnostics,
         "Host mean operator is the uncrowded first-moment operator: existing hosts are carried by H and reproduction contributes susceptible stage-1 recruits through R.",
         "Management mortality, fecundity reduction and LDD spread reduction are integrated over the node-level management adoption event in the reproductive term.",
         "Because recruitment and crowding are nonlinear in finite populations, this operator is intended for disease-free mean-background construction in the low-density/large-host regime."))
}

INApestTransitionPathogenGrowthFromHostStart <- function(
    Model=c("SIS","SIR","SEIR"),Transition,Nstages,InitialHostPopulation,
    SDDprob,LDDprob=0,LDDrate=0,EnvEstabProb=1,PropaguleEstablishment=1,
    ManageProb=0,MortalityProb=0,SpreadReduction=0,FecundityReduction=0,
    DispersalDensityFactor=0,K=1,SeedbankK=1,
    TransitionSDDprob=NULL,TransitionLDDprob=NULL,TransitionLDDrate=0,
    Beta,RecoveryProb,ProgressionProb=1,PathogenMortalityProb=0,
    ContactMatrix=NULL,StageMixing=NULL,Transmission=c("frequency","density"),DensityScale=1){
  Model<-match.arg(Model);Transmission<-match.arg(Transmission)
  X<-as.matrix(InitialHostPopulation);n<-nrow(X);S<-ncol(X);if(S!=Nstages)stop("InitialHostPopulation columns must equal Nstages")
  hm<-INApestTransitionHostMeanOperator(Transition,S,SDDprob,LDDprob,LDDrate,EnvEstabProb,
    PropaguleEstablishment,ManageProb,MortalityProb,SpreadReduction,FecundityReduction,
    DispersalDensityFactor,K,SeedbankK,TransitionSDDprob,TransitionLDDprob,TransitionLDDrate)
  x0<-as.numeric(t(X));x1<-as.numeric(hm$Operator%*%x0);N1<-matrix(x1,nrow=n,ncol=S,byrow=TRUE)
  pg<-INApestTransitionPathogenGrowthOperator(Model,Transition,S,N1,Beta,RecoveryProb,ProgressionProb,
    PathogenMortalityProb,ContactMatrix,StageMixing,Transmission,DensityScale,ManageProb,MortalityProb,
    TransitionSDDprob,TransitionLDDprob,TransitionLDDrate)
  pg$InitialHostPopulation<-X;pg$DerivedHostPopulationDiseaseStep<-N1;pg$HostMeanOperator<-hm$Operator
  pg$Diagnostics<-c(pg$Diagnostics,
    "HostPopulationDiseaseStep was derived from the same-timestep uncrowded disease-free host mean operator. For strongly crowded or highly stochastic host populations, provide/estimate the disease-step host background directly instead.")
  pg
}

INApestTransitionPathogenExactDiseaseStepFixedN <- function(
    Model=c("SIS","SIR","SEIR"),Ntimesteps,HostPopulationDiseaseStep,
    InitialInfected,InitialExposed=0,InitialRecovered=0,
    Beta,RecoveryProb,ProgressionProb=1,ImmunityLossProb=0,
    IntroductionProb=0,IntroductionNumber=1,
    ContactMatrix=NULL,StageMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,
    OutsideNodes=integer(0),OutsideUnits=integer(0),MaxStates=50000L,ReturnOperator=FALSE){
  Model<-match.arg(Model);Transmission<-match.arg(Transmission)
  Nmat<-as.matrix(HostPopulationDiseaseStep);n<-nrow(Nmat);S<-ncol(Nmat);U<-n*S
  if(any(Nmat<0|Nmat!=floor(Nmat)))stop("HostPopulationDiseaseStep must contain non-negative whole numbers for exact finite-state analysis")
  N<-as.integer(as.numeric(t(Nmat)))
  rv<-function(x,name,integer=FALSE){z<-.ina_tm_ns_matrix(x,n,S,name);z<-as.numeric(t(z));if(integer)z<-as.integer(z);z}
  I0<-rv(InitialInfected,"InitialInfected",TRUE);E0<-rv(InitialExposed,"InitialExposed",TRUE);R0<-rv(InitialRecovered,"InitialRecovered",TRUE)
  if(any(I0+E0+R0>N))stop("Initial pathogen-state counts exceed host counts")
  beta<-rv(Beta,"Beta");rec<-rv(RecoveryProb,"RecoveryProb");prog<-rv(ProgressionProb,"ProgressionProb");wan<-rv(ImmunityLossProb,"ImmunityLossProb")
  ip<-rv(IntroductionProb,"IntroductionProb");inum<-rv(IntroductionNumber,"IntroductionNumber",TRUE);ds<-rv(DensityScale,"DensityScale")
  Cn<-if(is.null(ContactMatrix))diag(n)else as.matrix(ContactMatrix);Sm<-if(is.null(StageMixing))matrix(1,S,S)else as.matrix(StageMixing);C<-kronecker(Cn,Sm)
  outside<-unique(as.integer(OutsideUnits));if(length(OutsideNodes))for(i in OutsideNodes)outside<-c(outside,(i-1L)*S+seq_len(S));outside<-unique(outside)
  if(Model=="SIS")ans<-INApestMetaSISExactNetworkFixedN(Ntimesteps,N,I0,beta,rec,C,ip,inum,Transmission,ds,OutsideNodes=outside,ReturnOperator=ReturnOperator)
  else ans<-INApestMetaCompartmentExactNetworkFixedN(Model,Ntimesteps,N,I0,E0,R0,beta,rec,prog,wan,ip,inum,C,Transmission,ds,OutsideNodes=outside,MaxStates=MaxStates,ReturnOperator=ReturnOperator)
  ans$Model<-"INApestMetaTransitionMatrix";ans$HostArchitecture<-"conditional fixed node x demographic-stage abundance";ans$PathogenModel<-Model
  ans$HostPopulationDiseaseStep<-Nmat;ans$OutsideUnits<-outside
  ans$Diagnostics<-c(ans$Diagnostics,
    "This exact finite-state branch conditions on a fixed node x stage host background and validates the pathogen step itself.",
    "It deliberately does not apply demographic stage transitions between pathogen steps; full Transition-Matrix invasion growth/escape uses the carrier-composed operators.")
  class(ans)<-c("INApestTransitionPathogenExactDiseaseStepFixedN","list");ans
}

INApestTransitionPathogenNextGeneration <- function(GrowthObject){
  x<-GrowthObject
  if(!is.list(x)||is.null(x$CarrierOperator)||is.null(x$TransmissionBlock)||is.null(x$PathogenModel))stop("Supply an INApestTransitionPathogenGrowthOperator result")
  H<-x$CarrierOperator;B<-x$TransmissionBlock;U<-nrow(H);model<-x$PathogenModel
  # Recover disease-only persistence/progression from the already composed
  # disease operator. This avoids re-resolving parameters and preserves exact
  # orientation.
  D<-x$DiseaseOperator
  if(model%in%c("SIS","SIR")){
    # D = D_I + B
    DI<-D-B;T<-DI%*%H;F<-B%*%H
  } else {
    # disease block rows/cols = E,I. New-infection contribution is only I -> E.
    Z<-matrix(0,U,U);Fdis<-rbind(cbind(Z,B),cbind(Z,Z));Tdis<-D-Fdis
    H2<-rbind(cbind(H,Z),cbind(Z,H));T<-Tdis%*%H2;F<-Fdis%*%H2
  }
  rhoT<-max(Mod(eigen(T,only.values=TRUE)$values))
  if(rhoT>=1-1e-12){
    K<-matrix(NA_real_,nrow(T),ncol(T));R0<-Inf
    note<-"The no-new-infection active-state process is non-transient (spectral radius >= 1), so the expected lifetime next-generation sum diverges."
  } else {
    K<-F%*%solve(diag(nrow(T))-T);R0<-max(Mod(eigen(K,only.values=TRUE)$values))
    note<-"K = F (I-T)^(-1) sums new active infections produced over the full remaining demographic/pathogen lifetime of an initial active-state cohort."
  }
  list(Model="INApestMetaTransitionMatrix",PathogenModel=model,TransitionWithoutNewInfection=T,
       NewInfectionOperator=F,NextGenerationOperator=K,R0=R0,
       PerTimestepLambda=x$IntrinsicRarePathogenMultiplier,
       ThresholdAgreement=if(is.finite(R0))sign(R0-1)==sign(x$IntrinsicRarePathogenMultiplier-1)else x$IntrinsicRarePathogenMultiplier>=1,
       Diagnostics=c(note,
         "R0 and lambda answer different questions: R0 is lifetime secondary active infection production, whereas lambda is the asymptotic per-timestep multiplier.",
         "Under the usual non-negative transient-state conditions their invasion thresholds agree even though their numerical values differ."))
}

INApestTransitionPathogenAnalytical <- function(
    Ntimesteps=10,Transition,Nstages,Pathogen,
    InformationMode=c("none","all_informed"),
    HostPopulationDiseaseStep=NULL,InitialHostPopulation=NULL,
    SDDprob=NULL,LDDprob=0,LDDrate=0,EnvEstabProb=1,PropaguleEstablishment=1,
    K=1,SeedbankK=1,ManageProb=0,MortalityProb=0,SpreadReduction=0,FecundityReduction=0,
    DispersalDensityFactor=0,TransitionSDDprob=NULL,TransitionLDDprob=NULL,TransitionLDDrate=0,
    StageMixing=NULL,InitialInfected=0,InitialExposed=0,InitialRecovered=0,
    OutsideNodes=integer(0),OutsideUnits=integer(0),Exact=TRUE,ExactMaxStates=50000L,
    ReturnOperators=FALSE){
  if(!inherits(Pathogen,"INApestPathogen"))stop("Pathogen must be created by INApestPathogen()")
  InformationMode<-match.arg(InformationMode)
  model<-as.character(Pathogen$Model)[1L];if(!model%in%c("SIS","SIR","SEIR"))stop("Transition-Matrix analytical pathogen methods currently support SIS, SIR and SEIR")
  ManageProbEff<-if(InformationMode=="all_informed")ManageProb else 0
  if(is.null(HostPopulationDiseaseStep)){
    if(is.null(InitialHostPopulation)||is.null(SDDprob))stop("Supply HostPopulationDiseaseStep, or InitialHostPopulation plus SDDprob to derive the disease-step host mean background")
    bg<-INApestTransitionPathogenGrowthFromHostStart(model,Transition,Nstages,InitialHostPopulation,SDDprob,LDDprob,LDDrate,
      EnvEstabProb,PropaguleEstablishment,ManageProbEff,MortalityProb,SpreadReduction,FecundityReduction,
      DispersalDensityFactor,K,SeedbankK,TransitionSDDprob,TransitionLDDprob,TransitionLDDrate,
      Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$PathogenMortalityProb,
      Pathogen$ContactMatrix,StageMixing,Pathogen$Transmission,Pathogen$DensityScale)
    HostPopulationDiseaseStep<-bg$DerivedHostPopulationDiseaseStep;growth<-bg
  } else {
    HostPopulationDiseaseStep<-as.matrix(HostPopulationDiseaseStep)
    growth<-INApestTransitionPathogenGrowthOperator(model,Transition,Nstages,HostPopulationDiseaseStep,
      Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$PathogenMortalityProb,
      Pathogen$ContactMatrix,StageMixing,Pathogen$Transmission,Pathogen$DensityScale,
      ManageProbEff,MortalityProb,TransitionSDDprob,TransitionLDDprob,TransitionLDDrate)
  }
  ng<-INApestTransitionPathogenNextGeneration(growth)
  br<-INApestTransitionPathogenBranching(model,Ntimesteps,Transition,Nstages,HostPopulationDiseaseStep,
    Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$PathogenMortalityProb,
    Pathogen$ContactMatrix,StageMixing,Pathogen$Transmission,Pathogen$DensityScale,
    ManageProbEff,MortalityProb,TransitionSDDprob,TransitionLDDprob,TransitionLDDrate,OutsideNodes,OutsideUnits)
  exact<-NULL
  exact_reason<-NULL
  if(isTRUE(Exact)){
    if(any(as.numeric(Pathogen$PathogenMortalityProb)!=0,na.rm=TRUE)) exact_reason<-"Conditional fixed-host exact disease-step analysis requires PathogenMortalityProb = 0 because pathogen deaths change host abundance."
    else if(any(HostPopulationDiseaseStep!=floor(HostPopulationDiseaseStep))) exact_reason<-"Conditional exact disease-step analysis requires whole-number host counts."
    else {
      exact<-tryCatch(INApestTransitionPathogenExactDiseaseStepFixedN(model,Ntimesteps,HostPopulationDiseaseStep,
        InitialInfected,InitialExposed,InitialRecovered,Pathogen$Beta,Pathogen$RecoveryProb,
        Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,
        Pathogen$ContactMatrix,StageMixing,Pathogen$Transmission,Pathogen$DensityScale,OutsideNodes,OutsideUnits,
        ExactMaxStates,ReturnOperators),error=function(e){exact_reason<<-conditionMessage(e);NULL})
    }
  }
  out<-list(Model="INApestMetaTransitionMatrix",PathogenModel=model,InformationMode=InformationMode,
    HostPopulationDiseaseStep=HostPopulationDiseaseStep,Growth=growth,NextGeneration=ng,Branching=br,
    ExactConditionalDiseaseStep=exact,ExactConditionalUnavailableReason=exact_reason,
    Diagnostics=c(
      "Transition-Matrix pathogen analysis separates demographic carrier transport from the pathogen state process.",
      if(InformationMode=="none") "Management-dependent vital rates are evaluated on an uninformed background (no management adoption)." else "Management-dependent vital rates are evaluated on an all-informed background; ManageProb remains the conditional adoption probability.",
      "Endogenous host/pathogen detection, information transfer and programmed-stop feedback are shared node-level states and are not closed by an independent-lineage branching process; use none/all-informed analytical envelopes and stochastic simulation for that feedback.",
      "The scalable headline invasion metrics are the full-timestep rare-pathogen lambda, lifetime next-generation R0, and multitype branching extinction/escape.",
      "The finite-state exact branch conditions on fixed node x stage host abundance; it validates disease mechanics but does not replace the carrier-composed growth/branching solution for evolving demographic stages."))
  class(out)<-c("INApestTransitionPathogenAnalytical","list");out
}

print.INApestTransitionPathogenAnalytical <- function(x,...){
  cat("INApest Transition-Matrix pathogen analytical result\n")
  cat("  Pathogen model:",x$PathogenModel,"\n")
  cat("  Per-timestep lambda:",format(x$Growth$IntrinsicRarePathogenMultiplier,digits=7),"\n")
  cat("  Lifetime R0:",format(x$NextGeneration$R0,digits=7),"\n")
  if(length(x$Branching$OutsideUnits))cat("  Branching escape by horizon (max inside starting type):",format(x$Branching$MaxEscapeByHorizonFromInside,digits=7),"\n")
  if(!is.null(x$ExactConditionalDiseaseStep))cat("  Conditional finite-state disease-step solution: available\n")
  else if(!is.null(x$ExactConditionalUnavailableReason))cat("  Conditional finite-state disease-step solution:",x$ExactConditionalUnavailableReason,"\n")
  invisible(x)
}


###############################################################################
### Unified dispatcher extension: Transition-Matrix pathogen analysis
###############################################################################
INApestAnalytical_pre_tm_pathogen <- INApestAnalytical
INApestAnalytical <- function(...) {
  args<-list(...)
  Model<-if(!is.null(args$Model))as.character(args$Model)[1L]else"INApest"
  Pathogen<-args$Pathogen
  if(!identical(Model,"INApestMetaTransitionMatrix")||is.null(Pathogen))
    return(do.call(INApestAnalytical_pre_tm_pathogen,args))
  if(is.null(args$Transition))stop("Transition is required for INApestMetaTransitionMatrix pathogen analysis")
  if(is.null(args$Nstages)){
    A0<-if(is.list(args$Transition))args$Transition[[1L]]else args$Transition
    args$Nstages<-nrow(as.matrix(A0))
  }
  # Legacy host analytical calls commonly use InitialState. For pathogen TM
  # analysis this is interpreted as the demographic host matrix at the start
  # of the timestep unless HostPopulationDiseaseStep is supplied explicitly.
  if(is.null(args$InitialHostPopulation)&&!is.null(args$InitialState))args$InitialHostPopulation<-args$InitialState
  keep<-c("Ntimesteps","Transition","Nstages","Pathogen","InformationMode","HostPopulationDiseaseStep","InitialHostPopulation",
          "SDDprob","LDDprob","LDDrate","EnvEstabProb","PropaguleEstablishment","K","SeedbankK",
          "ManageProb","MortalityProb","SpreadReduction","FecundityReduction","DispersalDensityFactor",
          "TransitionSDDprob","TransitionLDDprob","TransitionLDDrate","StageMixing","InitialInfected",
          "InitialExposed","InitialRecovered","OutsideNodes","OutsideUnits","Exact","ExactMaxStates","ReturnOperators")
  do.call(INApestTransitionPathogenAnalytical,args[intersect(names(args),keep)])
}


###############################################################################
### INApest point-pathogen analytical methods
###############################################################################

.ina_point_resolve <- function(x,points,timestep=1L,perm=1L,name="parameter"){
  n<-nrow(points);if(is.function(x)){
    fm<-names(formals(x));a<-list(points=points,timestep=timestep,perm=perm)
    if(!is.null(fm)&&!"..."%in%fm)a<-a[intersect(names(a),fm)]
    z<-do.call(x,a)
  } else z<-if(length(x)==1L)rep(as.numeric(x),n)else as.numeric(x)
  if(length(z)==1L)z<-rep(z,n)
  if(length(z)!=n||any(!is.finite(z)))stop(name," must resolve to one finite value per point")
  z
}

INApestPointPathogenEdgeMatrix <- function(
    Points,Pathogen,ContactRadius=Inf,ContactKernel=NULL,ContactProb=1,
    timestep=1L,perm=1L,IncludeDiagonal=FALSE){
  if(!is.data.frame(Points)||!all(c("x","y")%in%names(Points)))stop("Points must contain x and y")
  if(inherits(Pathogen,"INApestPointPathogenInteraction"))Pathogen<-Pathogen$Pathogen
  if(!inherits(Pathogen,"INApestPathogen"))stop("Pathogen must be created by INApestPathogen() or supplied through INApestPointPathogenInteraction()")
  n<-nrow(Points);Q<-matrix(0,n,n)
  if(n<1L)return(Q)
  if(n<2L&&!isTRUE(IncludeDiagonal))return(Q)
  beta<-.ina_point_resolve(Pathogen$Beta,Points,timestep,perm,"Beta")
  if(any(beta<0))stop("Beta must be non-negative")
  ds<-.ina_point_resolve(Pathogen$DensityScale,Points,timestep,perm,"DensityScale")
  if(any(ds<=0))stop("DensityScale must be positive")
  pairs<-expand.grid(i=seq_len(n),j=seq_len(n),KEEP.OUT.ATTRS=FALSE)
  if(!isTRUE(IncludeDiagonal))pairs<-pairs[pairs$i!=pairs$j,,drop=FALSE]
  dx<-Points$x[pairs$i]-Points$x[pairs$j];dy<-Points$y[pairs$i]-Points$y[pairs$j];d<-sqrt(dx^2+dy^2)
  keep<-d<=ContactRadius;pairs<-pairs[keep,,drop=FALSE];d<-d[keep]
  if(!nrow(pairs))return(Q)
  if(is.function(ContactProb)){
    fm<-names(formals(ContactProb));a<-list(distance=d,source=Points[pairs$i,,drop=FALSE],target=Points[pairs$j,,drop=FALSE],timestep=timestep,perm=perm)
    if(!is.null(fm)&&!"..."%in%fm)a<-a[intersect(names(a),fm)];cp<-do.call(ContactProb,a)
  } else cp<-ContactProb
  if(length(cp)==1L)cp<-rep(cp,nrow(pairs));if(length(cp)!=nrow(pairs)||any(!is.finite(cp)))stop("ContactProb must resolve per candidate pair")
  if(!is.null(ContactKernel)){
    if(!is.function(ContactKernel))stop("ContactKernel must be NULL or a function")
    fm<-names(formals(ContactKernel));a<-list(distance=d,source=Points[pairs$i,,drop=FALSE],target=Points[pairs$j,,drop=FALSE],timestep=timestep,perm=perm)
    if(!is.null(fm)&&!"..."%in%fm)a<-a[intersect(names(a),fm)];km<-do.call(ContactKernel,a)
    if(length(km)==1L)km<-rep(km,nrow(pairs));if(length(km)!=nrow(pairs)||any(!is.finite(km))||any(km<0))stop("ContactKernel must return finite non-negative multipliers")
    cp<-cp*km
  }
  cp<-pmin(1,pmax(0,cp))
  if(Pathogen$Transmission=="frequency")edge<-cp*pmin(1,pmax(0,beta[pairs$j]))
  else edge<-cp*(-expm1(-beta[pairs$j]/ds[pairs$j]))
  Q[cbind(pairs$j,pairs$i)]<-pmin(1,pmax(0,edge)) # row target, col source
  dimnames(Q)<-list(target=if("id"%in%names(Points))Points$id else seq_len(n),source=if("id"%in%names(Points))Points$id else seq_len(n))
  Q
}

INApestPointPathogenGrowthOperator <- function(
    Model=c("SIS","SIR","SEIR"),Points,Pathogen,
    ContactRadius=Inf,ContactKernel=NULL,ContactProb=1,timestep=1L,perm=1L){
  Model<-match.arg(Model);if(Model!=Pathogen$Model)warning("Model differs from Pathogen$Model; using requested analytical compartment structure")
  n<-nrow(Points);Q<-INApestPointPathogenEdgeMatrix(Points,Pathogen,ContactRadius,ContactKernel,ContactProb,timestep,perm)
  rec<-.ina_point_resolve(Pathogen$RecoveryProb,Points,timestep,perm,"RecoveryProb");mort<-.ina_point_resolve(Pathogen$PathogenMortalityProb,Points,timestep,perm,"PathogenMortalityProb")
  prog<-.ina_point_resolve(Pathogen$ProgressionProb,Points,timestep,perm,"ProgressionProb")
  if(any(rec<0|rec>1)||any(mort<0|mort>1)||any(rec+mort>1+1e-12)||any(prog<0|prog>1))stop("Invalid pathogen probabilities")
  stayI<-diag(1-rec-mort,n)
  if(Model%in%c("SIS","SIR"))G<-stayI+Q else G<-rbind(cbind(diag(1-prog,n),Q),cbind(diag(prog,n),stayI))
  lam<-max(Mod(eigen(G,only.values=TRUE)$values))
  list(Model="INApestMetaPoint",PathogenModel=Model,Operator=G,TransmissionEdgeMatrix=Q,
       IntrinsicRarePathogenMultiplier=lam,Classification=if(lam>1+1e-12)"growing"else if(lam<1-1e-12)"declining"else"threshold",
       Diagnostics=c(
         "The realised-contact layer is integrated out exactly: each source-target pair contributes an independent one-source transmission probability.",
         "Newly infected points do not recover/progress in the same timestep, matching the point pathogen engine.",
         if(Model=="SEIR")"SEIR preserves the one-timestep latent-state delay."else"SIS/SIR active-I linearisation."))
}

INApestPointPathogenNextGeneration <- function(GrowthObject){
  x<-GrowthObject;G<-x$Operator;Q<-x$TransmissionEdgeMatrix;n<-nrow(Q);model<-x$PathogenModel;Z<-matrix(0,n,n)
  if(model%in%c("SIS","SIR")){T<-G-Q;F<-Q}else{F<-rbind(cbind(Z,Q),cbind(Z,Z));T<-G-F}
  rhoT<-max(Mod(eigen(T,only.values=TRUE)$values))
  if(rhoT>=1-1e-12){K<-matrix(NA,nrow(T),ncol(T));R0<-Inf}else{K<-F%*%solve(diag(nrow(T))-T);R0<-max(Mod(eigen(K,only.values=TRUE)$values))}
  list(TransitionWithoutNewInfection=T,NewInfectionOperator=F,NextGenerationOperator=K,R0=R0,
       PerTimestepLambda=x$IntrinsicRarePathogenMultiplier,
       ThresholdAgreement=if(is.finite(R0))sign(R0-1)==sign(x$IntrinsicRarePathogenMultiplier-1)else x$IntrinsicRarePathogenMultiplier>=1,
       Diagnostics="Point-model R0 is the lifetime secondary-infection spectral radius after integrating over repeated infectious persistence at the same point.")
}

.ina_point_state_table <- function(Model,n){
  vals<-switch(Model,SIS=c("S","I"),SIR=c("S","I","R"),SEIR=c("S","E","I","R"));m<-length(vals)
  g<-expand.grid(rep(list(seq_len(m)),n),KEEP.OUT.ATTRS=FALSE);z<-matrix(vals[as.matrix(g)],nrow=nrow(g),ncol=n);colnames(z)<-paste0("p",seq_len(n));z
}

.ina_point_state_key <- function(z)paste(z,collapse="|")

.ina_point_exact_operator <- function(Model,Points,Pathogen,Q,MaxStates=100000L,timestep=1L,perm=1L){
  n<-nrow(Points);st<-.ina_point_state_table(Model,n);if(nrow(st)>MaxStates)stop("Exact point pathogen state space exceeds MaxStates")
  keys<-apply(st,1,.ina_point_state_key);idx<-setNames(seq_len(nrow(st)),keys);T<-matrix(0,nrow(st),nrow(st),dimnames=list(keys,keys))
  rec<-.ina_point_resolve(Pathogen$RecoveryProb,Points,timestep,perm,"RecoveryProb");mort<-.ina_point_resolve(Pathogen$PathogenMortalityProb,Points,timestep,perm,"PathogenMortalityProb")
  if(any(mort!=0))stop("Exact fixed-point operator requires PathogenMortalityProb = 0 because pathogen death removes host points")
  prog<-.ina_point_resolve(Pathogen$ProgressionProb,Points,timestep,perm,"ProgressionProb");wan<-.ina_point_resolve(Pathogen$ImmunityLossProb,Points,timestep,perm,"ImmunityLossProb");ip<-.ina_point_resolve(Pathogen$IntroductionProb,Points,timestep,perm,"IntroductionProb")
  if(any(rec<0|rec>1)||any(prog<0|prog>1)||any(wan<0|wan>1)||any(ip<0|ip>1))stop("Invalid pathogen probabilities")
  # Enumerate independent per-point end-state distributions conditional on the start state.
  for(r in seq_len(nrow(st))){s0<-st[r,];inf<-which(s0=="I");outs<-vector("list",n)
    for(j in seq_len(n)){
      if(s0[j]=="S"){
        pinf<-if(length(inf))1-prod(1-Q[j,inf])else 0
        # Introduction applies only if not newly infected.
        if(Model=="SEIR"){
          pE<-pinf+(1-pinf)*ip[j];outs[[j]]<-c(S=1-pE,E=pE)
        } else {pI<-pinf+(1-pinf)*ip[j];outs[[j]]<-c(S=1-pI,I=pI)}
      } else if(s0[j]=="I"){
        if(Model=="SIS")outs[[j]]<-c(S=rec[j],I=1-rec[j])
        else outs[[j]]<-c(I=1-rec[j],R=rec[j])
      } else if(s0[j]=="E")outs[[j]]<-c(E=1-prog[j],I=prog[j])
      else if(s0[j]=="R")outs[[j]]<-c(S=wan[j],R=1-wan[j])
    }
    # Cartesian product of point outcomes.
    namesv<-lapply(outs,names);grid<-expand.grid(lapply(namesv,seq_along),KEEP.OUT.ATTRS=FALSE)
    for(k in seq_len(nrow(grid))){dest<-character(n);pr<-1
      for(j in seq_len(n)){ii<-grid[k,j];dest[j]<-namesv[[j]][ii];pr<-pr*outs[[j]][ii]}
      if(pr>0)T[r,idx[[.ina_point_state_key(dest)]]]<-T[r,idx[[.ina_point_state_key(dest)]]]+pr
    }
  }
  if(max(abs(rowSums(T)-1))>1e-11)stop("Internal exact point operator error")
  attr(T,"States")<-st;T
}

.ina_point_hitting_probability <- function(T,target,avoid=rep(FALSE,nrow(T)),tol=1e-13,maxiter=100000L){
  target<-as.logical(target);avoid<-as.logical(avoid);if(length(target)!=nrow(T)||length(avoid)!=nrow(T))stop("hitting-set length mismatch")
  unknown<-!(target|avoid);u<-as.numeric(target)
  if(any(unknown)){
    Q<-T[unknown,unknown,drop=FALSE];r<-rowSums(T[unknown,target,drop=FALSE]);z<-rep(0,sum(unknown))
    for(k in seq_len(maxiter)){zn<-as.numeric(r+Q%*%z);if(max(abs(zn-z))<tol){z<-zn;break};z<-zn}
    u[unknown]<-pmin(1,pmax(0,z))
  }
  u
}

.ina_point_static_absorption <- function(T,States,model,outside,initial_index){
  active_state<-apply(States,1,function(z)any(z=="I"|(model=="SEIR"&z=="E")))
  extinct<-!active_state
  escaped<-if(any(outside))apply(States[,outside,drop=FALSE],1,function(z)any(z=="I"|(model=="SEIR"&z=="E")))else rep(FALSE,nrow(States))
  pe<-if(any(escaped)).ina_point_hitting_probability(T,escaped) else rep(0,nrow(T))
  # With no exogenous reintroduction, extinction is an absorbing pathogen event;
  # treating escape as the competing absorbing event gives the exact ordering.
  pext_first<-.ina_point_hitting_probability(T,extinct,escaped)
  pesc_first<-if(any(escaped)).ina_point_hitting_probability(T,escaped,extinct) else rep(0,nrow(T))
  transient<-!(extinct|escaped);etime<-NA_real_
  if(transient[initial_index]){
    Q<-T[transient,transient,drop=FALSE];rho<-if(nrow(Q))max(Mod(eigen(Q,only.values=TRUE)$values))else 0
    if(rho<1-1e-12){tt<-solve(diag(nrow(Q))-Q,rep(1,nrow(Q)));etime<-tt[which(which(transient)==initial_index)]}
  } else etime<-0
  list(EventualEscapeProbability=pe[initial_index],EscapeBeforeExtinctionProbability=pesc_first[initial_index],
       ExtinctionBeforeEscapeProbability=pext_first[initial_index],ExpectedTimestepsToEscapeOrExtinction=etime,
       PerStateEventualEscape=pe)
}

INApestPointPathogenExactFixedGeometry <- function(
    Ntimesteps,Points,Pathogen,ContactRadius=Inf,ContactKernel=NULL,ContactProb=1,
    InitialState=NULL,Outside=NULL,MaxStates=100000L,ReturnOperator=FALSE){
  cr_missing<-missing(ContactRadius);ck_missing<-missing(ContactKernel);cp_missing<-missing(ContactProb)
  interaction<-inherits(Pathogen,"INApestPointPathogenInteraction")
  if(interaction&&cr_missing)ContactRadius<-.ina_point_interaction_contact(Pathogen,"ContactRadius",Inf)
  if(interaction&&ck_missing)ContactKernel<-.ina_point_interaction_contact(Pathogen,"ContactKernel",NULL)
  if(interaction&&cp_missing)ContactProb<-.ina_point_interaction_contact(Pathogen,"ContactProb",1)
  Pathogen<-.ina_point_pathogen_spec(Pathogen)
  model<-Pathogen$Model;if(!model%in%c("SIS","SIR","SEIR"))stop("Exact point pathogen analysis supports SIS/SIR/SEIR")
  Ntimesteps<-as.integer(Ntimesteps);if(length(Ntimesteps)!=1L||is.na(Ntimesteps)||Ntimesteps<1L)stop("Ntimesteps must be a positive integer")
  n<-nrow(Points);st<-.ina_point_state_table(model,n);if(nrow(st)>MaxStates)stop("Exact point pathogen state space exceeds MaxStates")
  keys<-apply(st,1,.ina_point_state_key)
  if(is.null(InitialState)){
    z<-rep("S",n);ni<-as.integer(Pathogen$InitialInfected);ne<-if(model=="SEIR")as.integer(Pathogen$InitialExposed)else 0L;nr<-if(model!="SIS")as.integer(Pathogen$InitialRecovered)else 0L
    if(ni+ne+nr>n)stop("Initial pathogen counts exceed points")
    pos<-seq_len(n);if(ne){z[pos[seq_len(ne)]]<-"E";pos<-pos[-seq_len(ne)]};if(ni){z[pos[seq_len(ni)]]<-"I";pos<-pos[-seq_len(ni)]};if(nr)z[pos[seq_len(nr)]]<-"R";InitialState<-z
  }
  InitialState<-as.character(InitialState);if(length(InitialState)!=n)stop("InitialState must have one pathogen state per point")
  dist<-numeric(nrow(st));ii<-match(.ina_point_state_key(InitialState),keys);if(is.na(ii))stop("InitialState contains invalid state");dist[ii]<-1
  full<-matrix(0,Ntimesteps+1,nrow(st),dimnames=list(timestep=0:Ntimesteps,state=keys));full[1,]<-dist
  active<-function(row)row=="I"|(row=="E"&model=="SEIR")
  EI<-matrix(0,Ntimesteps+1,n);EI[1,]<-InitialState=="I";EE<-if(model=="SEIR")matrix(0,Ntimesteps+1,n)else NULL;if(model=="SEIR")EE[1,]<-InitialState=="E"
  ext<-numeric(Ntimesteps+1);ext[1]<-as.numeric(!any(active(InitialState)))
  outside<-if(is.null(Outside))rep(FALSE,n)else as.logical(Outside);if(length(outside)!=n)stop("Outside must be logical length nrow(Points)")
  escstates<-if(any(outside))apply(st[,outside,drop=FALSE],1,function(z)any(z=="I"|(model=="SEIR"&z=="E")))else rep(FALSE,nrow(st))
  escape<-numeric(Ntimesteps+1);safe<-dist
  if(any(escstates&dist>0)){escape[1]<-sum(dist[escstates]);safe[escstates]<-0}else safe[escstates]<-0
  Tlist<-vector("list",Ntimesteps);growths<-vector("list",Ntimesteps);ngs<-vector("list",Ntimesteps)
  for(tt in seq_len(Ntimesteps)){
    Qtt<-INApestPointPathogenEdgeMatrix(Points,Pathogen,ContactRadius,ContactKernel,ContactProb,tt,1L)
    Ttt<-.ina_point_exact_operator(model,Points,Pathogen,Qtt,MaxStates,tt,1L)
    if(!identical(rownames(Ttt),keys))stop("Internal point exact state ordering changed across timesteps")
    Tlist[[tt]]<-Ttt
    growths[[tt]]<-INApestPointPathogenGrowthOperator(model,Points,Pathogen,ContactRadius,ContactKernel,ContactProb,tt,1L)
    ngs[[tt]]<-INApestPointPathogenNextGeneration(growths[[tt]])
    dist<-as.numeric(dist%*%Ttt);full[tt+1,]<-dist
    EI[tt+1,]<-as.numeric(dist%*%(st=="I"));if(model=="SEIR")EE[tt+1,]<-as.numeric(dist%*%(st=="E"))
    ext[tt+1]<-sum(dist[apply(st,1,function(z)!any(z=="I"|(model=="SEIR"&z=="E")))])
    if(any(outside)){
      nxt<-as.numeric(safe%*%Ttt);newesc<-sum(nxt[escstates]);escape[tt+1]<-escape[tt]+newesc;nxt[escstates]<-0;safe<-nxt
    }
  }
  ord<-INApestTransitionPathogenPeriodicGrowth(lapply(growths,`[[`,"Operator"))
  static_inputs<-!any(vapply(list(Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$IntroductionProb,Pathogen$DensityScale,ContactProb,ContactKernel),is.function,logical(1)))
  static_operator<-static_inputs&&all(vapply(Tlist[-1L],function(M)isTRUE(all.equal(M,Tlist[[1L]],tolerance=0)),logical(1)))
  if(length(Tlist)==1L)static_operator<-static_inputs
  intro_zero<-all(.ina_point_resolve(Pathogen$IntroductionProb,Points,1L,1L,"IntroductionProb")==0)
  absorption<-if(static_operator&&intro_zero).ina_point_static_absorption(Tlist[[1L]],st,model,outside,ii)else NULL
  out<-list(Model="INApestMetaPoint",PathogenModel=model,Exact=TRUE,Points=Points,
    StepGrowth=growths,StepNextGeneration=ngs,OrderedGrowth=ord,StaticAbsorption=absorption,
    ExpectedInfected=EI,ExpectedExposed=EE,ExtinctionProbability=ext,
    Escape=if(any(outside))list(Outside=outside,ProbabilityByTimestep=escape,ProbabilityByHorizon=tail(escape,1))else NULL,
    StateTable=st,StateDistribution=full,
    Diagnostics=c(
      "Exact finite-state pathogen solution conditional on a fixed point set and geometry.",
      "The exact operator integrates over realised contacts analytically and preserves the engine's one-state-transition-per-timestep timing.",
      "Timestep-specific pathogen/contact parameters are handled by an ordered sequence of exact transition operators rather than by averaging them.",
      if(!is.null(absorption)) "For static no-reintroduction systems, eventual escape, extinction-before-escape and mean absorption time are solved from the exact finite-state chain." else "Eventual static absorption metrics are omitted when parameters are time-varying or pathogen reintroduction is active.",
      "Pathogen mortality is excluded here because it removes host points and therefore changes geometry."))
  if(ReturnOperator)out$Operators<-Tlist
  class(out)<-c("INApestPointPathogenExactFixedGeometry","list");out
}

INApestPointTransitionPathogenGrowthOperator <- function(
    Model=c("SIS","SIR","SEIR"),Representatives,HostTypePopulationDiseaseStep,
    ParentCarrierOperator,Pathogen,ContactRadius=Inf,ContactKernel=NULL,ContactProb=1,
    timestep=1L,perm=1L){
  Model<-match.arg(Model);R<-as.data.frame(Representatives);nt<-nrow(R);N<-as.numeric(HostTypePopulationDiseaseStep);H<-as.matrix(ParentCarrierOperator)
  if(length(N)!=nt||any(!is.finite(N))||any(N<0))stop("HostTypePopulationDiseaseStep must be finite non-negative length nrow(Representatives)")
  if(!all(dim(H)==c(nt,nt))||any(H<0)||any(!is.finite(H))||any(colSums(H)>1+1e-10))stop("ParentCarrierOperator must be a substochastic target x source matrix")
  Q<-INApestPointPathogenEdgeMatrix(R,Pathogen,ContactRadius,ContactKernel,ContactProb,timestep,perm,IncludeDiagonal=TRUE)
  # One-carrier offspring mean.  A carrier in source type i sees N_j
  # susceptible hosts in every other target type j, but only N_i - 1
  # susceptible hosts in its own type because the carrier itself occupies one host.
  SusceptibleTargets<-matrix(N,nrow=nt,ncol=nt)
  diag(SusceptibleTargets)<-pmax(0,N-1)
  B<-Q*SusceptibleTargets
  rec<-.ina_point_resolve(Pathogen$RecoveryProb,R,timestep,perm,"RecoveryProb");mort<-.ina_point_resolve(Pathogen$PathogenMortalityProb,R,timestep,perm,"PathogenMortalityProb");prog<-.ina_point_resolve(Pathogen$ProgressionProb,R,timestep,perm,"ProgressionProb")
  stayI<-diag(1-rec-mort,nt)
  if(Model%in%c("SIS","SIR")){D<-stayI+B;G<-D%*%H}
  else {D<-rbind(cbind(diag(1-prog,nt),B),cbind(diag(prog,nt),stayI));Z<-matrix(0,nt,nt);H2<-rbind(cbind(H,Z),cbind(Z,H));G<-D%*%H2}
  lam<-max(Mod(eigen(G,only.values=TRUE)$values))
  out<-list(Model="INApestPointTransitionMatrix",PathogenModel=Model,Operator=G,ParentCarrierOperator=H,
            DiseaseOperator=D,TransmissionEdgeMatrix=Q,TransmissionBlock=B,SusceptibleTargetsByCarrier=SusceptibleTargets,HostTypePopulationDiseaseStep=N,
            Representatives=R,Pathogen=Pathogen,
            ResolvedDiseaseParameters=list(RecoveryProb=rec,PathogenMortalityProb=mort,ProgressionProb=prog),
            IntrinsicRarePathogenMultiplier=lam,
            Classification=if(lam>1+1e-12)"growing"else if(lam<1-1e-12)"declining"else"threshold",
            Diagnostics=c(
              "Continuous point-stage movement is represented through the host analytical Parent carrier matrix; reproductive recruits are excluded from pathogen carriage because offspring are susceptible.",
              "Transmission is applied after stage movement/reproduction using the post-host-step type abundance and the point contact geometry.",
              "Within-type transmission is retained for contracted types containing multiple hosts; the infectious carrier itself is removed from its own susceptible target count.",
              "When Representatives are analysis-grid cell centres, within-cell contact is evaluated at zero representative distance and is therefore a grid-resolution approximation to the continuous-space contact process."))
  class(out)<-c("INApestPointTransitionPathogenGrowthOperator","list");out
}

INApestPointTransitionPathogenNextGeneration <- function(GrowthObject){
  x<-GrowthObject;H<-x$ParentCarrierOperator;B<-x$TransmissionBlock;D<-x$DiseaseOperator;nt<-nrow(H);Z<-matrix(0,nt,nt)
  if(x$PathogenModel%in%c("SIS","SIR")){T<-(D-B)%*%H;F<-B%*%H}
  else {Fdis<-rbind(cbind(Z,B),cbind(Z,Z));Tdis<-D-Fdis;H2<-rbind(cbind(H,Z),cbind(Z,H));T<-Tdis%*%H2;F<-Fdis%*%H2}
  rt<-max(Mod(eigen(T,only.values=TRUE)$values));if(rt>=1-1e-12){K<-matrix(NA,nrow(T),ncol(T));R0<-Inf}else{K<-F%*%solve(diag(nrow(T))-T);R0<-max(Mod(eigen(K,only.values=TRUE)$values))}
  list(R0=R0,NextGenerationOperator=K,TransitionWithoutNewInfection=T,NewInfectionOperator=F,
       PerTimestepLambda=x$IntrinsicRarePathogenMultiplier,
       ThresholdAgreement=if(is.finite(R0))sign(R0-1)==sign(x$IntrinsicRarePathogenMultiplier-1)else x$IntrinsicRarePathogenMultiplier>=1)
}

INApestPointTransitionPathogenBranching <- function(GrowthObject,timesteps,OutsideTypes=integer(0)){
  x<-GrowthObject;model<-x$PathogenModel;H<-x$ParentCarrierOperator;Q<-x$TransmissionEdgeMatrix;N<-x$HostTypePopulationDiseaseStep;nt<-nrow(H);Tn<-as.integer(timesteps)
  rec<-x$ResolvedDiseaseParameters$RecoveryProb;mort<-x$ResolvedDiseaseParameters$PathogenMortalityProb;prog<-x$ResolvedDiseaseParameters$ProgressionProb
  outside<-unique(as.integer(OutsideTypes));if(any(!outside%in%seq_len(nt)))stop("OutsideTypes must index analytical host types")
  integerN<-all(abs(N-round(N))<1e-10)
  trans_factor<-function(v,zE){
    Ssus<-pmax(0,N-as.numeric(seq_len(nt)==v));q<-Q[,v]
    if(integerN)prod((1-q+q*zE)^as.integer(round(Ssus))) else exp(sum(Ssus*q*(zE-1)))
  }
  if(model%in%c("SIS","SIR")){
    q<-rep(0,nt);h<-rep(1,nt);if(length(outside))h[outside]<-0;qh<-hh<-matrix(NA,nt,Tn)
    stepfun<-function(z){zn<-numeric(nt);for(u in seq_len(nt)){
      val<-max(0,1-sum(H[,u]));for(v in which(H[,u]>0)){stay<-1-rec[v]-mort[v];val<-val+H[v,u]*((1-stay)+stay*z[v])*trans_factor(v,z)};zn[u]<-val};pmin(1,pmax(0,zn))}
    for(tt in seq_len(Tn)){q<-stepfun(q);h<-stepfun(h);if(length(outside))h[outside]<-0;qh[,tt]<-q;hh[,tt]<-h}
  }else{
    q<-rep(0,2*nt);h<-rep(1,2*nt);if(length(outside))h[c(outside,nt+outside)]<-0;qh<-hh<-matrix(NA,2*nt,Tn)
    stepfun<-function(z){zn<-numeric(2*nt)
      for(u in seq_len(nt)){val<-max(0,1-sum(H[,u]));for(v in which(H[,u]>0))val<-val+H[v,u]*((1-prog[v])*z[v]+prog[v]*z[nt+v]);zn[u]<-val}
      for(u in seq_len(nt)){val<-max(0,1-sum(H[,u]));for(v in which(H[,u]>0)){stay<-1-rec[v]-mort[v];val<-val+H[v,u]*((1-stay)+stay*z[nt+v])*trans_factor(v,z[seq_len(nt)])};zn[nt+u]<-val};pmin(1,pmax(0,zn))}
    for(tt in seq_len(Tn)){q<-stepfun(q);h<-stepfun(h);if(length(outside))h[c(outside,nt+outside)]<-0;qh[,tt]<-q;hh[,tt]<-h}
  }
  inside<-setdiff(seq_len(nt),outside);itypes<-if(model=="SEIR")c(inside,nt+inside)else inside
  esc<-1-h
  list(Model="INApestPointTransitionMatrix",PathogenModel=model,ExtinctionByHorizon=q,EscapeByHorizon=esc,
       ExtinctionHistory=qh,EscapeHistory=1-hh,OutsideTypes=outside,InsideStartingTypes=itypes,
       MaxEscapeByHorizonFromInside=if(length(itypes))max(esc[itypes])else NA_real_,
       Approximation=if(integerN)"finite-count one-lineage branching on contracted point types"else"mean-matched Poisson branching on contracted point types",
       Diagnostics=c("Branching recursion composes stage/space carrier movement with point-contact transmission after movement.",
         "When analytical type abundances are non-integer expectations, transmission offspring use a mean-matched Poisson PGF; integer type abundances use the finite-count binomial PGF.",
         "For containment, include explicit outside analysis-grid cells in OutsideTypes; movement that leaves the analytical grid entirely is treated as loss rather than silently counted as pathogen escape."))
}

INApestPointTransitionPathogenAnalytical <- function(
    Ntimesteps=10,Nstages,Transition,InitialPoints,Pathogen,
    InformationMode=c("none","all_informed"),
    SDDkernel,LDDkernel=NULL,LDDrate=0,PropaguleEstablishment=1,EnvEstabProb=1,
    TransitionKernels=NULL,TransitionHabitatSearch=FALSE,ApplyHabitatToTransitions=FALSE,
    TransitionEstablishment=1,BlockedTransitionMortality=0,
    HabitatSuitability=NULL,HabitatSearchRadius=0,HabitatSearchCandidates=128,LocalK=Inf,KRadius=0,
    DetectionProb=0,DetectionSpatial=NULL,ManageProb=0,ManageSpatial=NULL,
    MortalityProb=0,MortalitySpatial=NULL,FecundityReduction=0,FecundityReductionSpatial=NULL,
    SpreadReduction=0,SpreadReductionSpatial=NULL,SpreadReductionAppliesTo=c("LDD","all"),
    InfoRetentionProb=1,InfoRadius=0,InfoTransferProb=0,InfoKernel=NULL,
    PointAnalysisGrid=NULL,MaxAnalysisCells=400L,KernelSamples=1000L,PointSeed=1L,
    OutsideEstablishmentProb=1,
    ContactRadius=Inf,ContactKernel=NULL,ContactProb=1,
    OutsideCells=integer(0),ReturnOperators=FALSE){
  cr_missing<-missing(ContactRadius);ck_missing<-missing(ContactKernel);cp_missing<-missing(ContactProb)
  interaction<-inherits(Pathogen,"INApestPointPathogenInteraction")
  Pspec<-.ina_point_pathogen_spec(Pathogen)
  if(interaction&&cr_missing)ContactRadius<-.ina_point_interaction_contact(Pathogen,"ContactRadius",Inf)
  if(interaction&&ck_missing)ContactKernel<-.ina_point_interaction_contact(Pathogen,"ContactKernel",NULL)
  if(interaction&&cp_missing)ContactProb<-.ina_point_interaction_contact(Pathogen,"ContactProb",1)
  Pathogen<-Pspec
  InformationMode<-match.arg(InformationMode);SpreadReductionAppliesTo<-match.arg(SpreadReductionAppliesTo)
  Ntimesteps<-as.integer(Ntimesteps);Nstages<-as.integer(Nstages)
  if(length(Ntimesteps)!=1L||is.na(Ntimesteps)||Ntimesteps<1L)stop("Ntimesteps must be a positive integer")
  if(length(Nstages)!=1L||is.na(Nstages)||Nstages<1L)stop("Nstages must be a positive integer")
  model<-Pathogen$Model;if(!model%in%c("SIS","SIR","SEIR"))stop("Point-transition pathogen analytical methods support SIS/SIR/SEIR")
  .ina_pt_require_point_helpers(TRUE)
  analysis_grid<-.ina_pt_analysis_grid(PointAnalysisGrid,HabitatSuitability,DetectionSpatial,ManageSpatial,MortalitySpatial,FecundityReductionSpatial,SpreadReductionSpatial)
  if(is.null(analysis_grid)&&(!is.null(TransitionKernels)||is.function(ContactKernel)||is.finite(ContactRadius)))
    stop("Spatial transition/contact analysis requires PointAnalysisGrid so movement and contact geometry can be represented")
  if(!is.null(analysis_grid)&&is.null(PointAnalysisGrid)&&analysis_grid$nrow*analysis_grid$ncol>MaxAnalysisCells)stop("Inherited analysis grid exceeds MaxAnalysisCells; supply a coarser PointAnalysisGrid")
  reps<-.ina_pt_transition_representatives(InitialPoints,analysis_grid,Nstages)
  x<-.ina_pt_transition_initial(InitialPoints,analysis_grid,Nstages)
  mode<-if(InformationMode=="none")"none"else"all_informed"
  steps<-vector("list",Ntimesteps);growths<-vector("list",Ntimesteps);ngs<-vector("list",Ntimesteps);branches<-vector("list",Ntimesteps)
  hosttraj<-matrix(0,length(x),Ntimesteps+1L);hosttraj[,1L]<-x
  outside_types<-integer(0)
  if(length(OutsideCells)){
    if(is.null(analysis_grid))stop("OutsideCells requires PointAnalysisGrid")
    nc<-analysis_grid$nrow*analysis_grid$ncol;if(any(!OutsideCells%in%seq_len(nc)))stop("OutsideCells contains invalid analysis cells")
    outside_types<-unlist(lapply(seq_len(Nstages),function(s)(s-1L)*nc+OutsideCells),use.names=FALSE)
  }
  for(tt in seq_len(Ntimesteps)){
    st<-.ina_pt_transition_operator_step(tt,Ntimesteps,reps,analysis_grid,Transition,Nstages,
      SDDkernel,LDDkernel,LDDrate,PropaguleEstablishment,EnvEstabProb,
      TransitionKernels,TransitionHabitatSearch,ApplyHabitatToTransitions,TransitionEstablishment,
      BlockedTransitionMortality,HabitatSuitability,HabitatSearchRadius,HabitatSearchCandidates,
      DetectionProb,DetectionSpatial,ManageProb,ManageSpatial,MortalityProb,MortalitySpatial,
      FecundityReduction,FecundityReductionSpatial,SpreadReduction,SpreadReductionSpatial,
      SpreadReductionAppliesTo,InfoRetentionProb,InfoRadius,InfoTransferProb,InfoKernel,
      KernelSamples,PointSeed,mode,OutsideEstablishmentProb)
    Ghost<-if(mode=="none")st$G0 else st$GH;H<-if(mode=="none")st$Parent0 else st$ParentH
    x1<-as.numeric(Ghost%*%x);hosttraj[,tt+1L]<-x1
    pg<-INApestPointTransitionPathogenGrowthOperator(model,reps,x1,H,Pathogen,ContactRadius,ContactKernel,ContactProb,tt,1L)
    growths[[tt]]<-pg;ngs[[tt]]<-INApestPointTransitionPathogenNextGeneration(pg)
    branches[[tt]]<-INApestPointTransitionPathogenBranching(pg,1L,outside_types)
    steps[[tt]]<-st;x<-x1
  }
  ops<-lapply(growths,`[[`,"Operator");cyc<-INApestTransitionPathogenPeriodicGrowth(ops)
  brseq<-INApestPointPathogenBranchingSequence(growths,outside_types)
  result<-list(Model="INApestPointTransitionMatrix",PathogenModel=model,InformationMode=InformationMode,
    Representatives=reps,HostTypeTrajectory=hosttraj,StepGrowth=growths,StepNextGeneration=ngs,
    CycleGrowth=cyc,Branching=brseq,StepBranching=branches,OutsideTypes=outside_types,
    Diagnostics=c(
      "Host stage/space movement and reproduction are contracted with the existing INApest point-transition analytical kernel machinery; only the persisting/progressing Parent operator carries pathogen state.",
      "Reproductive offspring contribute to the susceptible host background before pathogen contact but never inherit pathogen state.",
      "Pathogen contact between contracted stage x grid types is evaluated at representative cell-centre geometry; refine PointAnalysisGrid where contact kernels vary sharply with distance.",
      "Dynamic detection-driven information is intentionally not collapsed into an average management probability here; use InformationMode='none' or 'all_informed' for analytical screening and the stochastic point model for endogenous information feedback.",
      if(!is.infinite(LocalK)||KRadius>0) "Finite LocalK/KRadius is a nonlinear local-crowding process and is excluded from the rare-pathogen carrier/background operator; use stochastic simulation once crowding is material." else "Local crowding is inactive in this analytical screening configuration."))
  if(ReturnOperators)result$HostStepData<-steps
  class(result)<-c("INApestPointTransitionPathogenAnalytical","list");result
}

print.INApestPointTransitionPathogenAnalytical <- function(x,...){
  cat("INApest point-transition pathogen analytical result\n")
  cat("  Pathogen model:",x$PathogenModel,"\n")
  cat("  Equivalent per-timestep multiplier:",format(x$CycleGrowth$PerTimestepMultiplier,digits=7),"\n")
  invisible(x)
}

###############################################################################
### Time-inhomogeneous branching and high-level INApestMetaPoint pathogen path
###############################################################################

.ina_point_pathogen_spec <- function(Pathogen){
  if(inherits(Pathogen,"INApestPointPathogenInteraction")){
    if(is.null(Pathogen$Pathogen)||!inherits(Pathogen$Pathogen,"INApestPathogen"))
      stop("Point pathogen interaction does not contain a valid INApestPathogen specification")
    return(Pathogen$Pathogen)
  }
  if(inherits(Pathogen,"INApestPathogen"))return(Pathogen)
  stop("Pathogen must be an INApestPathogen or INApestPointPathogenInteraction")
}

.ina_point_interaction_contact <- function(Pathogen,name,default){
  if(!inherits(Pathogen,"INApestPointPathogenInteraction"))return(default)
  e<-environment(Pathogen$Contact)
  if(!is.null(e)&&exists(name,envir=e,inherits=TRUE))get(name,envir=e,inherits=TRUE) else default
}

.ina_point_branch_step <- function(GrowthObject,z){
  x<-GrowthObject;model<-x$PathogenModel;H<-x$ParentCarrierOperator;Q<-x$TransmissionEdgeMatrix
  N<-x$HostTypePopulationDiseaseStep;nt<-nrow(H)
  rec<-x$ResolvedDiseaseParameters$RecoveryProb;mort<-x$ResolvedDiseaseParameters$PathogenMortalityProb
  prog<-x$ResolvedDiseaseParameters$ProgressionProb
  integerN<-all(abs(N-round(N))<1e-10)
  trans_factor<-function(v,zE){
    Ssus<-pmax(0,N-as.numeric(seq_len(nt)==v));q<-Q[,v]
    if(integerN)prod((1-q+q*zE)^as.integer(round(Ssus))) else exp(sum(Ssus*q*(zE-1)))
  }
  if(model%in%c("SIS","SIR")){
    if(length(z)!=nt)stop("Branching state length mismatch")
    zn<-numeric(nt)
    for(u in seq_len(nt)){
      val<-max(0,1-sum(H[,u]))
      for(v in which(H[,u]>0)){
        stay<-1-rec[v]-mort[v]
        val<-val+H[v,u]*((1-stay)+stay*z[v])*trans_factor(v,z)
      }
      zn[u]<-val
    }
  } else {
    if(length(z)!=2L*nt)stop("SEIR branching state length mismatch")
    zn<-numeric(2L*nt)
    for(u in seq_len(nt)){
      val<-max(0,1-sum(H[,u]))
      for(v in which(H[,u]>0))val<-val+H[v,u]*((1-prog[v])*z[v]+prog[v]*z[nt+v])
      zn[u]<-val
    }
    for(u in seq_len(nt)){
      val<-max(0,1-sum(H[,u]))
      for(v in which(H[,u]>0)){
        stay<-1-rec[v]-mort[v]
        val<-val+H[v,u]*((1-stay)+stay*z[nt+v])*trans_factor(v,z[seq_len(nt)])
      }
      zn[nt+u]<-val
    }
  }
  pmin(1,pmax(0,zn))
}

INApestPointPathogenBranchingSequence <- function(GrowthObjects,OutsideTypes=integer(0)){
  if(!is.list(GrowthObjects)||!length(GrowthObjects))stop("GrowthObjects must be a non-empty list")
  model<-GrowthObjects[[1L]]$PathogenModel;nt<-nrow(GrowthObjects[[1L]]$ParentCarrierOperator)
  if(any(vapply(GrowthObjects,function(x)x$PathogenModel!=model||nrow(x$ParentCarrierOperator)!=nt,logical(1))))
    stop("All GrowthObjects must use the same pathogen model and analytical type set")
  d<-if(model=="SEIR")2L*nt else nt
  outside<-unique(as.integer(OutsideTypes));if(any(!outside%in%seq_len(nt)))stop("OutsideTypes must index analytical host types")
  outside_active<-if(model=="SEIR")c(outside,nt+outside)else outside
  Tn<-length(GrowthObjects);qh<-hh<-matrix(NA_real_,d,Tn)
  for(hor in seq_len(Tn)){
    q<-rep(0,d)
    for(tt in rev(seq_len(hor)))q<-.ina_point_branch_step(GrowthObjects[[tt]],q)
    qh[,hor]<-q
    z<-rep(1,d);if(length(outside_active))z[outside_active]<-0
    for(tt in rev(seq_len(hor))){z<-.ina_point_branch_step(GrowthObjects[[tt]],z);if(length(outside_active))z[outside_active]<-0}
    hh[,hor]<-z
  }
  inside<-setdiff(seq_len(nt),outside);inside_active<-if(model=="SEIR")c(inside,nt+inside)else inside
  list(Model=GrowthObjects[[1L]]$Model,PathogenModel=model,
       ExtinctionHistory=qh,EscapeHistory=1-hh,
       ExtinctionByHorizon=qh[,Tn],EscapeByHorizon=1-hh[,Tn],
       OutsideTypes=outside,InsideStartingTypes=inside_active,
       MaxEscapeByHorizonFromInside=if(length(inside_active))max((1-hh[,Tn])[inside_active])else NA_real_,
       Diagnostics=c(
         "Time-varying branching PGFs are composed in chronological order: F1(F2(...FT(z)...)), rather than by averaging timestep-specific extinction or escape probabilities.",
         "Integer contracted host counts use finite-count binomial transmission offspring; fractional analytical host backgrounds use a mean-matched Poisson PGF.",
         "Escape is first passage to an active pathogen state in an outside analytical type at a timestep boundary."))
}

INApestMetaPointPathogenAnalytical <- function(
    Ntimesteps=10,InitialPoints,Pathogen,
    InformationMode=c("none","all_informed"),
    Survival=1,PropaguleProduction,PropaguleEstablishment=1,EnvEstabProb=1,
    SDDkernel,LDDkernel=NULL,LDDrate=0,
    HabitatSuitability=NULL,HabitatSearchRadius=0,HabitatSearchCandidates=128,LocalK=Inf,KRadius=0,
    DetectionProb=0,DetectionSpatial=NULL,ManageProb=0,ManageSpatial=NULL,
    MortalityProb=0,MortalitySpatial=NULL,FecundityReduction=0,FecundityReductionSpatial=NULL,
    SpreadReduction=0,SpreadReductionSpatial=NULL,SpreadReductionAppliesTo=c("LDD","all"),
    InfoRetentionProb=1,InfoRadius=0,InfoTransferProb=0,InfoKernel=NULL,
    PointAnalysisGrid=NULL,MaxAnalysisCells=400L,KernelSamples=1000L,PointSeed=1L,
    OutsideEstablishmentProb=1,
    ContactRadius=Inf,ContactKernel=NULL,ContactProb=1,
    OutsideCells=integer(0),ReturnOperators=FALSE){
  cr_missing<-missing(ContactRadius);ck_missing<-missing(ContactKernel);cp_missing<-missing(ContactProb)
  interaction<-inherits(Pathogen,"INApestPointPathogenInteraction")
  Pspec<-.ina_point_pathogen_spec(Pathogen)
  if(interaction&&cr_missing)ContactRadius<-.ina_point_interaction_contact(Pathogen,"ContactRadius",Inf)
  if(interaction&&ck_missing)ContactKernel<-.ina_point_interaction_contact(Pathogen,"ContactKernel",NULL)
  if(interaction&&cp_missing)ContactProb<-.ina_point_interaction_contact(Pathogen,"ContactProb",1)
  Pathogen<-Pspec
  InformationMode<-match.arg(InformationMode);SpreadReductionAppliesTo<-match.arg(SpreadReductionAppliesTo)
  Ntimesteps<-as.integer(Ntimesteps);if(length(Ntimesteps)!=1L||is.na(Ntimesteps)||Ntimesteps<1L)stop("Ntimesteps must be a positive integer")
  .ina_pt_require_point_helpers(FALSE)
  analysis_grid<-.ina_pt_analysis_grid(PointAnalysisGrid,HabitatSuitability,DetectionSpatial,ManageSpatial,MortalitySpatial,FecundityReductionSpatial,SpreadReductionSpatial)
  if(is.null(analysis_grid)&&(is.function(ContactKernel)||is.finite(ContactRadius)))
    stop("Finite-radius or geometry-dependent pathogen contact requires PointAnalysisGrid for the scalable moving-point analytical contraction")
  if(!is.null(analysis_grid)&&is.null(PointAnalysisGrid)&&analysis_grid$nrow*analysis_grid$ncol>MaxAnalysisCells)
    stop("Inherited analysis grid exceeds MaxAnalysisCells; supply a coarser PointAnalysisGrid")
  reps<-.ina_pt_metapoint_representatives(InitialPoints,analysis_grid)
  x<-.ina_pt_metapoint_initial(InitialPoints,analysis_grid)
  mode<-if(InformationMode=="none")"none"else"all_informed"
  steps<-growths<-ngs<-vector("list",Ntimesteps)
  hosttraj<-matrix(0,length(x),Ntimesteps+1L);hosttraj[,1L]<-x
  outside_types<-integer(0)
  if(length(OutsideCells)){
    if(is.null(analysis_grid))stop("OutsideCells requires PointAnalysisGrid")
    nc<-analysis_grid$nrow*analysis_grid$ncol;if(any(!OutsideCells%in%seq_len(nc)))stop("OutsideCells contains invalid analysis cells")
    outside_types<-unique(as.integer(OutsideCells))
  }
  for(tt in seq_len(Ntimesteps)){
    st<-.ina_pt_metapoint_operator_step(tt,Ntimesteps,reps,analysis_grid,
      Survival,PropaguleProduction,PropaguleEstablishment,EnvEstabProb,
      SDDkernel,LDDkernel,LDDrate,HabitatSuitability,HabitatSearchRadius,HabitatSearchCandidates,
      DetectionProb,DetectionSpatial,ManageProb,ManageSpatial,MortalityProb,MortalitySpatial,
      FecundityReduction,FecundityReductionSpatial,SpreadReduction,SpreadReductionSpatial,
      SpreadReductionAppliesTo,InfoRetentionProb,InfoRadius,InfoTransferProb,InfoKernel,
      KernelSamples,PointSeed,mode,OutsideEstablishmentProb)
    Ghost<-if(mode=="none")st$G0 else st$GH;H<-if(mode=="none")st$Parent0 else st$ParentH
    x1<-as.numeric(Ghost%*%x);hosttraj[,tt+1L]<-x1
    pg<-INApestPointTransitionPathogenGrowthOperator(Pathogen$Model,reps,x1,H,Pathogen,ContactRadius,ContactKernel,ContactProb,tt,1L)
    pg$Model<-"INApestMetaPoint";growths[[tt]]<-pg;ngs[[tt]]<-INApestPointTransitionPathogenNextGeneration(pg)
    steps[[tt]]<-st;x<-x1
  }
  cyc<-INApestTransitionPathogenPeriodicGrowth(lapply(growths,`[[`,"Operator"))
  br<-INApestPointPathogenBranchingSequence(growths,outside_types)
  out<-list(Model="INApestMetaPoint",PathogenModel=Pathogen$Model,InformationMode=InformationMode,
    Representatives=reps,HostTypeTrajectory=hosttraj,StepGrowth=growths,StepNextGeneration=ngs,
    OrderedGrowth=cyc,Branching=br,OutsideTypes=outside_types,
    Diagnostics=c(
      "Existing infected points are carried only by the host Parent survival/management operator; Poisson reproductive offspring enter the susceptible host background and never inherit pathogen state.",
      "Pathogen contact is applied after survival, management and recruitment, matching INApestMetaPoint event order.",
      "The scalable moving-point solution supports none/all-informed management backgrounds; endogenous detection-information feedback remains a correlated point process and is intentionally left to stochastic simulation.",
      if(!is.infinite(LocalK)||KRadius>0) "Finite LocalK/KRadius is a nonlinear local-crowding process and is excluded from the rare-pathogen background operator; use stochastic simulation once crowding is material." else "Local crowding is inactive in this analytical screening configuration.",
      if(is.null(analysis_grid)) "Homogeneous contact contraction: all hosts occupy one analytical type." else "Spatial contact and host dispersal are contracted to PointAnalysisGrid cell-centre types; refine the grid where contact kernels or habitat vary sharply."))
  if(ReturnOperators)out$HostStepData<-steps
  class(out)<-c("INApestMetaPointPathogenAnalytical","list");out
}

print.INApestMetaPointPathogenAnalytical <- function(x,...){
  cat("INApest MetaPoint pathogen analytical result\n")
  cat("  Pathogen model:",x$PathogenModel,"\n")
  cat("  Equivalent per-timestep multiplier:",format(x$OrderedGrowth$PerTimestepMultiplier,digits=7),"\n")
  if(length(x$OutsideTypes))cat("  Branching escape by horizon (max inside starting type):",format(x$Branching$MaxEscapeByHorizonFromInside,digits=7),"\n")
  invisible(x)
}
###############################################################################
### Unified dispatcher extension: point pathogen analytical methods
###############################################################################
INApestAnalytical_pre_point_pathogen <- INApestAnalytical
INApestAnalytical <- function(...) {
  args<-list(...)
  Model<-if(!is.null(args$Model))as.character(args$Model)[1L]else"INApest"
  Pathogen<-args$Pathogen
  if(is.null(Pathogen)||!Model%in%c("INApestMetaPoint","INApestPointTransitionMatrix"))
    return(do.call(INApestAnalytical_pre_point_pathogen,args))
  if(Model=="INApestMetaPoint"){
    keep<-c("Ntimesteps","InitialPoints","Pathogen","InformationMode","Survival","PropaguleProduction",
      "PropaguleEstablishment","EnvEstabProb","SDDkernel","LDDkernel","LDDrate","HabitatSuitability",
      "HabitatSearchRadius","HabitatSearchCandidates","LocalK","KRadius","DetectionProb","DetectionSpatial","ManageProb","ManageSpatial",
      "MortalityProb","MortalitySpatial","FecundityReduction","FecundityReductionSpatial","SpreadReduction",
      "SpreadReductionSpatial","SpreadReductionAppliesTo","InfoRetentionProb","InfoRadius","InfoTransferProb","InfoKernel",
      "PointAnalysisGrid","MaxAnalysisCells","KernelSamples","PointSeed","OutsideEstablishmentProb",
      "ContactRadius","ContactKernel","ContactProb","OutsideCells","ReturnOperators")
    return(do.call(INApestMetaPointPathogenAnalytical,args[intersect(names(args),keep)]))
  }
  if(is.null(args$Transition))stop("Transition is required for INApestPointTransitionMatrix pathogen analysis")
  if(is.null(args$Nstages)){
    A0<-if(is.array(args$Transition)&&length(dim(args$Transition))==3L)args$Transition[,,1L] else if(is.matrix(args$Transition))args$Transition else NULL
    if(is.null(A0))stop("Nstages is required when Transition is supplied as a function")
    args$Nstages<-nrow(A0)
  }
  keep<-c("Ntimesteps","Nstages","Transition","InitialPoints","Pathogen","InformationMode","SDDkernel","LDDkernel","LDDrate",
    "PropaguleEstablishment","EnvEstabProb","TransitionKernels","TransitionHabitatSearch","ApplyHabitatToTransitions",
    "TransitionEstablishment","BlockedTransitionMortality","HabitatSuitability","HabitatSearchRadius","HabitatSearchCandidates","LocalK","KRadius",
    "DetectionProb","DetectionSpatial","ManageProb","ManageSpatial","MortalityProb","MortalitySpatial","FecundityReduction",
    "FecundityReductionSpatial","SpreadReduction","SpreadReductionSpatial","SpreadReductionAppliesTo","InfoRetentionProb",
    "InfoRadius","InfoTransferProb","InfoKernel","PointAnalysisGrid","MaxAnalysisCells","KernelSamples","PointSeed",
    "OutsideEstablishmentProb","ContactRadius","ContactKernel","ContactProb","OutsideCells","ReturnOperators")
  do.call(INApestPointTransitionPathogenAnalytical,args[intersect(names(args),keep)])
}
###############################################################################
###############################################################################

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
