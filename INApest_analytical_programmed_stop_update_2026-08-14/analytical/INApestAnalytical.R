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
