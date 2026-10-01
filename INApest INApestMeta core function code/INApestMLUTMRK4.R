###############################################################################
### INApest MLU x Transition Matrix RK4 + linked stochastic bridge
### Repository source-hygiene release v0.5.4 -- 2026-10-01
###
### Scientific behavior: frozen MLUTM RK v0.5.3.
### Source-hygiene changes only:
###   * compatibility RK and compartment helpers are MLUTM-namespaced;
###   * deterministic v0.5 reference implementation is private;
###   * v0.5.3 linked stochastic implementation is the sole public local step;
###   * no assignment can overwrite canonical shared INApestRK4.R,
###     INApestStochasticCompartmentBridge.R or INApestContinuousLocalDynamics.R.
###############################################################################

###############################################################################
### INApest generic RK4 kernel and stochastic state bridge v0.5
### Date: 2026-10-01
###
### Descends from the validated September 2026 ContinuousLocalDynamics contract:
###   * classical RK4 numerical kernel;
###   * default RKMaxStep = 0.025;
###   * invalid negative continuous states fail rather than being clipped;
###   * conversion back to INApest counts is stochastic, never deterministic
###     rounding when a continuous endpoint is fractional.
###
### This source is deliberately architecture-neutral. The MLUTM adapter owns
### state packing and the biological meaning of the joint H/P/B state.
###############################################################################

.INApestMLUTMRK4Step <- function(state, derivative, t, h, ..., NegativeTolerance = 1e-10) {
  y <- as.numeric(state)
  if (!length(y) || any(!is.finite(y))) stop("RK state must be finite", call. = FALSE)
  if (!is.function(derivative)) stop("derivative must be a function", call. = FALSE)
  eval_d <- function(tt, yy) {
    if (any(yy < -NegativeTolerance)) stop("RK intermediate state became negative; reduce RKMaxStep or revise the rate function", call. = FALSE)
    z <- as.numeric(derivative(tt, yy, ...))
    if (length(z) != length(y) || any(!is.finite(z))) stop("RK derivative must return one finite value per state element", call. = FALSE)
    z
  }
  k1 <- eval_d(t, y)
  k2 <- eval_d(t + h/2, y + h*k1/2)
  k3 <- eval_d(t + h/2, y + h*k2/2)
  k4 <- eval_d(t + h, y + h*k3)
  out <- y + h*(k1 + 2*k2 + 2*k3 + k4)/6
  if (any(out < -NegativeTolerance)) stop("RK endpoint became negative; reduce RKMaxStep or revise the rate function", call. = FALSE)
  out[abs(out) <= NegativeTolerance] <- 0
  out
}

.INApestMLUTMRK4Integrate <- function(state, derivative, t0 = 0, dt = 1,
                                RKMaxStep = 0.025, ..., NegativeTolerance = 1e-10) {
  dt <- as.numeric(dt); RKMaxStep <- as.numeric(RKMaxStep); t0 <- as.numeric(t0)
  if (length(dt) != 1L || !is.finite(dt) || dt < 0) stop("dt must be a finite non-negative scalar", call. = FALSE)
  if (length(RKMaxStep) != 1L || !is.finite(RKMaxStep) || RKMaxStep <= 0) stop("RKMaxStep must be finite and > 0", call. = FALSE)
  if (dt == 0) return(as.numeric(state))
  nstep <- max(1L, ceiling(dt/RKMaxStep))
  h <- dt/nstep
  y <- as.numeric(state); tt <- t0
  for (k in seq_len(nstep)) {
    y <- .INApestMLUTMRK4Step(y, derivative, tt, h, ..., NegativeTolerance = NegativeTolerance)
    tt <- tt + h
  }
  y
}

# Independent unbiased stochastic rounding. Exact integer inputs are returned
# without an RNG draw, which is important for zero-rate/feature-off regressions.
.INApestMLUTMStochasticRound <- function(x, tolerance = 1e-10) {
  z <- as.numeric(x)
  if (any(!is.finite(z)) || any(z < -tolerance)) stop("State bridge requires finite non-negative values", call. = FALSE)
  z[z < 0 & z >= -tolerance] <- 0
  lo <- floor(z)
  frac <- z - lo
  out <- as.integer(lo)
  idx <- which(frac > tolerance)
  if (length(idx)) out[idx] <- out[idx] + stats::rbinom(length(idx), 1L, frac[idx])
  out
}

# Retained public name from the validated stochastic-bridge work. For a cohort
# exposed to named competing hazards, allocate each individual exactly once to
# one hazard or survival. Hazard values are continuous-time rates.
.INApestMLUTMStochasticCompetingHazards <- function(n, hazards, dt = 1) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 0L) stop("n must be a non-negative integer", call. = FALSE)
  h <- as.numeric(hazards)
  if (any(!is.finite(h)) || any(h < 0)) stop("hazards must be finite and non-negative", call. = FALSE)
  if (!length(h)) return(c(survive = n))
  total <- sum(h)
  nm <- names(hazards); if (is.null(nm)) nm <- paste0("event", seq_along(h))
  if (n == 0L || total <= 0) return(setNames(c(integer(length(h)), n), c(nm, "survive")))
  p_event <- 1 - exp(-total*dt)
  p <- c(p_event * h/total, 1-p_event)
  z <- as.integer(stats::rmultinom(1L, n, p)[,1L])
  setNames(z, c(nm, "survive"))
}

.INApestMLUTMCombineDerivatives <- function(...) {
  fs <- list(...)
  if (length(fs) == 1L && is.list(fs[[1L]]) && !is.function(fs[[1L]])) fs <- fs[[1L]]
  if (!length(fs) || !all(vapply(fs, is.function, logical(1)))) stop("All derivatives must be functions", call. = FALSE)
  function(t, state, ...) {
    vals <- lapply(fs, function(f) as.numeric(f(t, state, ...)))
    n <- unique(vapply(vals, length, integer(1)))
    if (length(n) != 1L) stop("Combined derivatives returned different state lengths", call. = FALSE)
    Reduce(`+`, vals)
  }
}

# User-facing local-biology specification. The MLUTM adapter supplies the
# standard pathogen and biocontrol attack terms; these functions add optional
# host and agent local rates without owning stage progression or movement.


###############################################################################
# INApestStochasticCompartmentBridge.R
# Generic stochastic RK4 bridge for linked whole-count compartment transfers
# v0.4 -- 30 September 2026
#
# A CompartmentRateFunction returns, for a node x compartment integer state:
#   GainRate          node x compartment gross Poisson entry rates;
#   TransitionHazards node x source x destination per-capita hazards; and
#   ExitHazards       node x source x exit-cause per-capita hazards.
#
# Internal transfers are linked events: one realised S -> I event removes one S
# and creates one I. They are therefore not represented as independent loss and
# gain draws. Within each stochastic substep, RK4 follows the deterministic
# environment and the corresponding time-inhomogeneous compartment generator.
# Starting cohorts are assigned one multinomial endpoint fate; Poisson gains are
# thinned into endpoint compartments or exit causes. This preserves integer
# states and the exact realised identity
#
#   total_end = total_start + realised_gains - realised_exits
#
# while retaining the RK4 continuous-time mean path used to calculate fates.
###############################################################################


# Compatibility helper used by the validated compartment bridge. The MLUTM
# v0.5 kernel predates this helper but has the same supported-argument intent.

# Auxiliary RK4 step for probability/moment states. Unlike the biological
# v0.5 kernel, auxiliary moment vectors legitimately contain signed
# derivatives and must not be clipped using NegativeTolerance.
.inatmlu_scb_rk4_allow <- function(State, Time, Step, RateFunction) {
  y <- as.numeric(State)
  f <- function(tt, yy) {
    z <- as.numeric(RateFunction(tt, yy))
    if (length(z) != length(y) || any(!is.finite(z)))
      stop("Auxiliary RK4 derivative must be finite and match state length")
    z
  }
  k1 <- f(Time, y)
  k2 <- f(Time + Step/2, y + Step*k1/2)
  k3 <- f(Time + Step/2, y + Step*k2/2)
  k4 <- f(Time + Step, y + Step*k3)
  y + Step*(k1 + 2*k2 + 2*k3 + k4)/6
}
if (!exists(".inatmlu_rk4_call_supported", mode = "function")) {
  .inatmlu_rk4_call_supported <- function(fun, args) {
    fm <- names(formals(fun))
    if (!is.null(fm) && !("..." %in% fm)) args <- args[intersect(names(args), fm)]
    do.call(fun, args)
  }
}

.inatmlu_scb_assert_state <- function(State, label = "State") {
  if (is.data.frame(State)) State <- as.matrix(State)
  if (!is.matrix(State)) stop(label, " must be a node x compartment matrix")
  if (!is.numeric(State)) stop(label, " must be numeric")
  if (!nrow(State) || !ncol(State)) stop(label, " must have positive dimensions")
  if (is.null(colnames(State)) || any(!nzchar(colnames(State))) || anyDuplicated(colnames(State)))
    stop(label, " must have unique non-empty compartment column names")
  if (any(!is.finite(State)) || any(State < 0) || any(State != floor(State)))
    stop(label, " must contain finite non-negative whole counts")
  storage.mode(State) <- "integer"
  State
}

.inatmlu_scb_matrix <- function(x, State, label) {
  n <- nrow(State); k <- ncol(State)
  if (is.null(dim(x))) {
    z <- as.numeric(x)
    if (length(z) == 1L) x <- matrix(z, n, k)
    else if (length(z) == n) x <- matrix(rep(z, k), n, k)
    else if (length(z) == n*k) x <- matrix(z, n, k)
    else stop(label, " must be scalar, length nodes, or node x compartment")
  } else x <- as.matrix(x)
  if (!all(dim(x) == c(n,k))) stop(label, " must have node x compartment dimensions")
  if (any(!is.finite(x))) stop(label, " must be finite")
  dimnames(x) <- dimnames(State)
  x
}

.inatmlu_scb_transition_array <- function(x, State) {
  n <- nrow(State); k <- ncol(State); comps <- colnames(State)
  if (is.null(x)) {
    out <- array(0, dim = c(n,k,k), dimnames = list(rownames(State), comps, comps))
    return(out)
  }
  if (length(dim(x)) != 3L || !all(dim(x) == c(n,k,k)))
    stop("TransitionHazards must have dimensions nodes x source compartment x destination compartment")
  out <- array(as.numeric(x), dim = c(n,k,k),
               dimnames = list(rownames(State), comps, comps))
  if (any(!is.finite(out)) || any(out < 0))
    stop("TransitionHazards must be finite and non-negative")
  for (j in seq_len(k)) out[,j,j] <- 0
  out
}

.inatmlu_scb_exit_array <- function(x, State) {
  n <- nrow(State); k <- ncol(State); comps <- colnames(State)
  if (is.null(x))
    return(array(numeric(), dim = c(n,k,0L), dimnames = list(rownames(State), comps, character())))
  d <- dim(x)
  if (length(d) == 2L && all(d == c(n,k))) {
    nm <- if (!is.null(attr(x,"cause"))) as.character(attr(x,"cause"))[1L] else "exit"
    x <- array(as.numeric(x), dim = c(n,k,1L),
               dimnames = list(rownames(State), comps, nm))
  }
  if (length(dim(x)) != 3L || dim(x)[1] != n || dim(x)[2] != k)
    stop("ExitHazards must have dimensions nodes x source compartment x exit cause")
  if (dim(x)[3] == 0L)
    return(array(numeric(), dim = c(n,k,0L), dimnames = list(rownames(State), comps, character())))
  causes <- dimnames(x)[[3]]
  if (is.null(causes) || any(!nzchar(causes)) || anyDuplicated(causes))
    stop("ExitHazards must have unique non-empty exit-cause names")
  out <- array(as.numeric(x), dim = dim(x),
               dimnames = list(rownames(State), comps, causes))
  if (any(!is.finite(out)) || any(out < 0))
    stop("ExitHazards must be finite and non-negative")
  out
}

.inatmlu_scb_rates <- function(fun, t, State, Parameters = NULL, Context = NULL,
                               expected_causes = NULL, extra = list()) {
  ans <- .inatmlu_rk4_call_supported(
    fun,
    c(list(t=t,time=t,state=State,State=State,pars=Parameters,Parameters=Parameters,
           context=Context,Context=Context), extra)
  )
  if (!is.list(ans))
    stop("CompartmentRateFunction must return a list")
  gain <- ans$GainRate
  if (is.null(gain)) gain <- State*0
  gain <- .inatmlu_scb_matrix(gain, State, "GainRate")
  if (any(gain < 0)) stop("GainRate must be non-negative")
  trans <- .inatmlu_scb_transition_array(ans$TransitionHazards, State)
  exits <- .inatmlu_scb_exit_array(ans$ExitHazards, State)
  causes <- dimnames(exits)[[3]]
  if (!is.null(expected_causes) && !identical(causes, expected_causes))
    stop("ExitHazards causes changed during an RK step")
  list(GainRate=gain, TransitionHazards=trans, ExitHazards=exits, Causes=causes)
}

.inatmlu_scb_generator <- function(Tmat, Emat) {
  # Preserve the one-compartment case: array extraction in base R drops a
  # 1 x 1 transition slice to a scalar unless it is normalised here.
  Tmat <- as.matrix(Tmat)
  k <- nrow(Tmat)
  if (is.null(dim(Emat))) Emat <- matrix(Emat, nrow = k)
  else Emat <- as.matrix(Emat)
  Q <- Tmat
  diag(Q) <- 0
  exit_total <- if (ncol(Emat)) rowSums(Emat) else rep(0,k)
  diag(Q) <- -(rowSums(Q) + exit_total)
  Q
}

# Deterministic mean derivative implied by a compartment rate function.
INApestMLUTMCompartmentMeanDerivative <- function(t, State, CompartmentRateFunction,
                                              Parameters=NULL, Context=NULL, ...) {
  State <- as.matrix(State)
  fl <- .inatmlu_scb_rates(CompartmentRateFunction,t,State,Parameters,Context,
                           expected_causes=NULL,extra=list(...))
  out <- fl$GainRate*0
  for (i in seq_len(nrow(State))) {
    E <- if (length(fl$Causes)) fl$ExitHazards[i,,,drop=FALSE][1,,] else matrix(numeric(),ncol(State),0L)
    if (length(fl$Causes)==1L) E <- matrix(E,ncol(State),1L,dimnames=list(colnames(State),fl$Causes))
    Q <- .inatmlu_scb_generator(fl$TransitionHazards[i,,], E)
    out[i,] <- fl$GainRate[i,] + as.numeric(State[i,] %*% Q)
  }
  out
}

# RK4 moments over one stochastic substep.
INApestMLUTMStochasticCompartmentMoments <- function(State, Time, Step,
                                                 CompartmentRateFunction,
                                                 Parameters=NULL, Context=NULL, ...) {
  State <- .inatmlu_scb_assert_state(State)
  if (!is.function(CompartmentRateFunction)) stop("CompartmentRateFunction must be a function")
  Time <- as.numeric(Time)[1L]; Step <- as.numeric(Step)[1L]
  if (!is.finite(Time)) stop("Time must be finite")
  if (!is.finite(Step) || Step <= 0) stop("Step must be finite and > 0")
  n <- nrow(State); k <- ncol(State); extra <- list(...)
  f0 <- .inatmlu_scb_rates(CompartmentRateFunction,Time,State,Parameters,Context,
                           expected_causes=NULL,extra=extra)
  causes <- f0$Causes; cN <- length(causes)

  nx <- n*k; np <- n*k*k; nd <- n*k*cN
  i_x <- seq_len(nx)
  off <- nx
  i_P <- off + seq_len(np); off <- off + np
  i_D <- if (nd) off + seq_len(nd) else integer(); off <- off + nd
  i_B <- off + seq_len(np); off <- off + np
  i_GD <- if (nd) off + seq_len(nd) else integer(); off <- off + nd
  i_GT <- off + seq_len(nx)

  P0 <- array(0,c(n,k,k)); for (j in seq_len(k)) P0[,j,j] <- 1
  D0 <- array(0,c(n,k,cN)); B0 <- array(0,c(n,k,k)); GD0 <- array(0,c(n,k,cN)); GT0 <- matrix(0,n,k)
  z0 <- c(as.numeric(State),as.numeric(P0),as.numeric(D0),as.numeric(B0),as.numeric(GD0),as.numeric(GT0))

  deriv <- function(t,state,...) {
    z <- as.numeric(state)
    X <- matrix(z[i_x],n,k,dimnames=dimnames(State))
    P <- array(z[i_P],c(n,k,k))
    D <- if (nd) array(z[i_D],c(n,k,cN)) else array(numeric(),c(n,k,0L))
    B <- array(z[i_B],c(n,k,k))
    GD <- if (nd) array(z[i_GD],c(n,k,cN)) else array(numeric(),c(n,k,0L))
    fl <- .inatmlu_scb_rates(CompartmentRateFunction,t,X,Parameters,Context,
                             expected_causes=causes,extra=extra)
    dX <- matrix(0,n,k); dP <- array(0,c(n,k,k)); dD <- array(0,c(n,k,cN))
    dB <- array(0,c(n,k,k)); dGD <- array(0,c(n,k,cN)); dGT <- fl$GainRate
    for (ii in seq_len(n)) {
      E <- if (cN) matrix(fl$ExitHazards[ii,,],k,cN) else matrix(numeric(),k,0L)
      Q <- .inatmlu_scb_generator(fl$TransitionHazards[ii,,],E)
      dX[ii,] <- fl$GainRate[ii,] + as.numeric(X[ii,] %*% Q)
      dP[ii,,] <- P[ii,,] %*% Q
      dB[ii,,] <- diag(fl$GainRate[ii,],nrow=k,ncol=k) + B[ii,,] %*% Q
      if (cN) {
        dD[ii,,] <- P[ii,,] %*% E
        dGD[ii,,] <- B[ii,,] %*% E
      }
    }
    c(as.numeric(dX),as.numeric(dP),as.numeric(dD),as.numeric(dB),as.numeric(dGD),as.numeric(dGT))
  }

  # The v0.5 MLUTM RK kernel uses (state, derivative, t, h). Auxiliary
  # probability/moment states may transiently be negative during RK stages, so
  # domain validity is checked below from the completed moment identities.
  zend <- .inatmlu_scb_rk4_allow(z0, Time, Step, deriv)
  X <- matrix(zend[i_x],n,k,dimnames=dimnames(State))
  P <- array(zend[i_P],c(n,k,k),dimnames=list(rownames(State),colnames(State),colnames(State)))
  D <- if (nd) array(zend[i_D],c(n,k,cN),dimnames=list(rownames(State),colnames(State),causes)) else array(numeric(),c(n,k,0L),dimnames=list(rownames(State),colnames(State),character()))
  B <- array(zend[i_B],c(n,k,k),dimnames=list(rownames(State),colnames(State),colnames(State)))
  GD <- if (nd) array(zend[i_GD],c(n,k,cN),dimnames=list(rownames(State),colnames(State),causes)) else array(numeric(),c(n,k,0L),dimnames=list(rownames(State),colnames(State),character()))
  GT <- matrix(zend[i_GT],n,k,dimnames=dimnames(State))

  tol <- 5e-9
  vals <- c(X,P,D,B,GD,GT)
  if (any(!is.finite(vals)) || any(X < -tol) || any(P < -tol) || any(D < -tol) ||
      any(B < -tol) || any(GD < -tol) || any(GT < -tol))
    stop("RK compartment moment step left the valid biological domain; reduce RKMaxStep or revise rates")
  P[P<0 & P>=-tol] <- 0; D[D<0 & D>=-tol] <- 0
  B[B<0 & B>=-tol] <- 0; GD[GD<0 & GD>=-tol] <- 0; GT[GT<0 & GT>=-tol] <- 0

  start_closure <- 0; gain_closure <- 0; state_closure <- 0
  for (ii in seq_len(n)) {
    Pi <- matrix(P[ii,,],k,k)
    Bi <- matrix(B[ii,,],k,k)
    start_exit <- if(cN) rowSums(matrix(D[ii,,],k,cN)) else rep(0,k)
    gain_exit <- if(cN) rowSums(matrix(GD[ii,,],k,cN)) else rep(0,k)
    start_closure <- max(start_closure, max(abs(rowSums(Pi) + start_exit - 1)))
    gain_closure <- max(gain_closure, max(abs(rowSums(Bi) + gain_exit - GT[ii,])))
    expected <- as.numeric(State[ii,] %*% Pi) + colSums(Bi)
    state_closure <- max(state_closure, max(abs(X[ii,]-expected)))
  }
  scale <- max(1,max(abs(c(X,GT,State))))
  if (start_closure > 5e-8 || gain_closure > 5e-8*scale || state_closure > 5e-8*scale)
    stop("RK compartment moment identities failed; reduce RKMaxStep")

  list(ExpectedEnd=X, StartLiveProbabilities=P, StartExitProbabilities=D,
       GainLiveMeans=B, GainExitMeans=GD, TotalGainMeans=GT, Causes=causes,
       Closure=c(start=start_closure,gains=gain_closure,state=state_closure))
}

INApestMLUTMStochasticCompartmentStep <- function(State, Time, Step,
                                              CompartmentRateFunction,
                                              Parameters=NULL, Context=NULL, ...) {
  State <- .inatmlu_scb_assert_state(State)
  mom <- INApestMLUTMStochasticCompartmentMoments(State,Time,Step,CompartmentRateFunction,
                                              Parameters,Context,...)
  n <- nrow(State); k <- ncol(State); cN <- length(mom$Causes)
  final <- matrix(0L,n,k,dimnames=dimnames(State))
  gains <- matrix(0L,n,k,dimnames=dimnames(State))
  exits <- matrix(0L,n,cN,dimnames=list(rownames(State),mom$Causes))

  for (ii in seq_len(n)) for (a in seq_len(k)) {
    probs <- c(mom$StartLiveProbabilities[ii,a,], if(cN) mom$StartExitProbabilities[ii,a,] else numeric())
    if (any(probs < -1e-10) || sum(probs) > 1 + 1e-8) stop("Invalid start-cohort endpoint probabilities")
    probs[probs<0] <- 0
    if (abs(sum(probs)-1) > 1e-8) probs <- probs/sum(probs)
    z <- as.integer(stats::rmultinom(1L,size=State[ii,a],prob=probs))
    final[ii,] <- final[ii,] + z[seq_len(k)]
    if (cN) exits[ii,] <- exits[ii,] + z[k+seq_len(cN)]

    live_mean <- mom$GainLiveMeans[ii,a,]
    dead_mean <- if(cN) mom$GainExitMeans[ii,a,] else numeric()
    gz_live <- stats::rpois(k,pmax(0,live_mean))
    gz_dead <- if(cN) stats::rpois(cN,pmax(0,dead_mean)) else integer()
    final[ii,] <- final[ii,] + gz_live
    gains[ii,a] <- sum(gz_live) + sum(gz_dead)
    if (cN) exits[ii,] <- exits[ii,] + gz_dead
  }

  if (any(final < 0) || any(final != floor(final))) stop("Internal compartment bridge produced invalid counts")
  lhs <- rowSums(final)
  rhs <- rowSums(State) + rowSums(gains) - if(cN) rowSums(exits) else 0
  if (!identical(as.integer(lhs),as.integer(rhs)))
    stop("Internal compartment bridge host bookkeeping failure")
  storage.mode(final) <- "integer"; storage.mode(gains) <- "integer"; storage.mode(exits) <- "integer"
  list(State=final,Gains=gains,Exits=exits,Moments=mom)
}

INApestMLUTMStochasticCompartmentIntegrate <- function(State, CompartmentRateFunction,
                                                   Duration, MaxStep=0.025,
                                                   StartTime=0, Parameters=NULL,
                                                   Context=NULL, ...) {
  State <- .inatmlu_scb_assert_state(State)
  Duration <- as.numeric(Duration)[1L]; MaxStep <- as.numeric(MaxStep)[1L]; StartTime <- as.numeric(StartTime)[1L]
  if (!is.finite(Duration) || Duration < 0) stop("Duration must be finite and >= 0")
  if (!is.finite(MaxStep) || MaxStep <= 0) stop("MaxStep must be finite and > 0")
  if (!is.finite(StartTime)) stop("StartTime must be finite")
  if (Duration == 0) return(list(State=State,Gains=State*0,Exits=matrix(numeric(),nrow(State),0L),NSubsteps=0L,InternalStep=0))
  nsub <- max(1L,as.integer(ceiling(Duration/MaxStep-1e-14)))
  h <- Duration/nsub; x <- State; gains <- State*0; exits <- NULL; causes <- NULL
  max_closure <- c(start=0,gains=0,state=0); tt <- StartTime
  for (jj in seq_len(nsub)) {
    z <- INApestMLUTMStochasticCompartmentStep(x,tt,h,CompartmentRateFunction,Parameters,Context,...)
    if (is.null(causes)) {
      causes <- colnames(z$Exits)
      exits <- matrix(0L,nrow(State),length(causes),dimnames=list(rownames(State),causes))
    } else if (!identical(causes,colnames(z$Exits))) stop("Exit causes changed between compartment substeps")
    x <- z$State; gains <- gains + z$Gains
    if (length(causes)) exits <- exits + z$Exits
    max_closure <- pmax(max_closure,z$Moments$Closure)
    tt <- StartTime + jj*h
  }
  if (is.null(exits)) exits <- matrix(numeric(),nrow(State),0L)
  list(State=x,Gains=gains,Exits=exits,NSubsteps=nsub,InternalStep=h,
       MaxStep=MaxStep,StartTime=StartTime,EndTime=StartTime+Duration,
       MaxClosureError=max_closure)
}


###############################################################################
### INApest MLU x Transition Matrix continuous-local-dynamics RK4 adapter v0.5
### Date: 2026-10-01
###
### Frozen biological parent: MLUTM H+P+B v0.4.2 (native 17/17 PASS).
###
### State contract
### --------------
### H[node, land-use, host-stage]
### P[node, land-use, host-stage, pathogen-state]
### Q[node, agent-stage] for each biocontrol agent
###
### Event ownership in active RK mode
### ---------------------------------
### 1. Response layer: information-gated management and management mortality.
### 2. RK4: continuous local biology within existing host demographic stages.
###    - optional user host rate;
###    - standard continuous pathogen transmission/progression/recovery/mortality;
###    - standard biocontrol attack/recruitment plus optional user agent rate.
### 3. Discrete MLUTM boundary: host stage progression/stasis, transition-associated
###    movement, fecundity, reproductive dispersal and recruitment. Pathogen state
###    is carried with hosts and new host recruits enter S.
### 4. Biocontrol spatial movement.
### 5. Response information persistence/SEAM and surveillance.
###
### RK never owns pest stage progression, transition movement, fecundity or
### recruitment. That boundary remains the frozen MLUTM parent contract.
###
### ContinuousBiology = NULL delegates directly to frozen v0.4.2, preserving
### exact results and RNG state.
###############################################################################

.inatmlu_rk_require <- function() {
  needed <- c(".INApestMLUTMRK4Integrate","INApestMLUTMContinuousBiology",
              "INApestMLUTMPathogenBiocontrolResponseStep",
              "INApestMLUTMPathogenHostStep","INApestMLUTMHostStep",
              "INApestMLUTMPathogenReconcile","INApestMLUTMPathogenInitial")
  miss <- needed[!vapply(needed, exists, logical(1), mode="function")]
  if (length(miss)) stop("Required RK/MLUTM function(s) missing: ",paste(miss,collapse=", "),call.=FALSE)
  invisible(TRUE)
}

.inatmlu_rk_call <- function(f, args) {
  if (is.null(f)) return(NULL)
  fm <- names(formals(f))
  if (!is.null(fm) && !("..." %in% fm)) args <- args[intersect(names(args),fm)]
  do.call(f,args)
}

.inatmlu_rk_zero_pathogen <- function(Pathogen) {
  if (is.null(Pathogen)) return(NULL)
  out <- Pathogen
  for (nm in c("Beta","RecoveryProb","ProgressionProb","PathogenMortalityProb",
               "ImmunityLossProb","IntroductionProb","IntroductionNumber")) out[[nm]] <- 0
  class(out) <- class(Pathogen)
  out
}

.inatmlu_rk_prob_to_hazard <- function(p, dt) {
  p <- as.numeric(p)
  if (any(!is.finite(p)) || any(p < 0 | p > 1)) stop("Probability-to-hazard conversion requires values in [0,1]",call.=FALSE)
  -log(pmax(.Machine$double.eps,1-p))/dt
}

.inatmlu_rk_competing_exit_hazards <- function(a,b,dt) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (any(!is.finite(a)) || any(!is.finite(b)) || any(a<0|b<0|a+b>1+1e-12))
    stop("RecoveryProb + PathogenMortalityProb must resolve to <= 1 in RK mode",call.=FALSE)
  s <- a+b
  total <- -log(pmax(.Machine$double.eps,1-pmin(1,s)))/dt
  ha <- hb <- numeric(length(s)); nz <- s>0
  ha[nz] <- total[nz]*a[nz]/s[nz]; hb[nz] <- total[nz]*b[nz]/s[nz]
  list(a=ha,b=hb)
}

.inatmlu_rk_pack <- function(H,P,Q) {
  vals <- numeric(); meta <- list(); pos <- 0L
  add <- function(name,x) {
    n <- length(x); idx <- if(n) seq.int(pos+1L,pos+n) else integer()
    meta[[name]] <<- list(idx=idx,dim=dim(x),dimnames=dimnames(x)); pos <<- pos+n
    vals <<- c(vals,as.numeric(x))
  }
  if (is.null(P)) add("H",H) else add("P",P)
  if (!is.null(Q)) for (nm in names(Q)) add(paste0("Q::",nm),Q[[nm]])
  list(state=vals,meta=meta)
}

.inatmlu_rk_unpack <- function(y,meta,H_template,P_template,Q_template) {
  take <- function(key) {
    m <- meta[[key]]; z <- y[m$idx]
    if (is.null(m$dim)) return(z)
    array(z,dim=m$dim,dimnames=m$dimnames)
  }
  if (is.null(P_template)) { H <- take("H"); P <- NULL } else {
    P <- take("P"); H <- array(apply(P,c(1,2,3),sum),dim=dim(P)[1:3],dimnames=dimnames(P)[1:3])
  }
  Q <- NULL
  if (!is.null(Q_template)) { Q <- setNames(vector("list",length(Q_template)),names(Q_template)); for(nm in names(Q)) Q[[nm]]<-take(paste0("Q::",nm)) }
  list(H=H,P=P,Q=Q)
}

.inatmlu_rk_bridge_v052 <- function(H,P,Q,tol=1e-10) {
  if (is.null(P)) {
    h <- array(.INApestMLUTMStochasticRound(H,tol),dim=dim(H),dimnames=dimnames(H)); p <- NULL
  } else {
    p <- array(.INApestMLUTMStochasticRound(P,tol),dim=dim(P),dimnames=dimnames(P))
    h <- array(apply(p,c(1,2,3),sum),dim=dim(P)[1:3],dimnames=dimnames(P)[1:3])
  }
  q <- NULL
  if (!is.null(Q)) { q<-Q; for(nm in names(q)) q[[nm]]<-matrix(.INApestMLUTMStochasticRound(q[[nm]],tol),nrow=nrow(q[[nm]]),ncol=ncol(q[[nm]]),dimnames=dimnames(q[[nm]])) }
  list(H=h,P=p,Q=q)
}

.inatmlu_rk_validate_bc <- function(Biocontrol,Prepared,timestep) {
  if (is.null(Biocontrol)) return(invisible(TRUE))
  # Local agent demography must be expressed in BiocontrolRateFunction in RK mode.
  # The existing discrete Transition is therefore required to be identity so a
  # second local-demography step cannot be applied silently.
  for (a in Prepared$Biocontrol$Agents) {
    P <- .inabc_resolve_transition(a$Transition,timestep,Prepared$Context,a)
    if (!isTRUE(all.equal(P,diag(nrow(P)),tolerance=1e-12,check.attributes=FALSE)))
      stop("Active RK mode requires identity Biocontrol agent Transition; encode local agent demography in BiocontrolRateFunction to avoid double-stepping",call.=FALSE)
  }
  invisible(TRUE)
}

.inatmlu_rk_add_releases <- function(Q,Prepared,timestep) {
  if (is.null(Q)) return(NULL)
  out <- Q
  for (nm in names(out)) {
    a <- Prepared$Biocontrol$Agents[[nm]]
    rel <- .inabc_resolve_state_matrix(a$Release,timestep,Prepared$Context,a,"Release")
    out[[nm]] <- out[[nm]] + rel
  }
  out
}

.inatmlu_rk_move_q <- function(Q,Prepared,timestep) {
  if (is.null(Q)) return(NULL)
  out <- Q
  for (nm in names(out)) {
    a <- Prepared$Biocontrol$Agents[[nm]]
    M <- .inabc_resolve_movement(a$Movement,timestep,Prepared$Context,a)
    out[[nm]] <- .inabc_move_agent(out[[nm]],M,a$MovementStages)
  }
  out
}

.inatmlu_rk_pathogen_derivative <- function(H,P,Pathogen,timestep,Ntimesteps,StageMixing,PathogenLandUseMixing,dt) {
  if (is.null(Pathogen)) return(list(dP=NULL, incidence=NULL, deaths=NULL))
  d <- dim(H); nn<-d[1L]; L<-d[2L]; S<-d[3L]; nc<-nn*L; states<-Pathogen$States
  if (is.function(Pathogen$IntroductionProb))
    stop("Active RK v0.5 requires static Pathogen IntroductionProb=0; resolver introductions are reserved for a later gate",call.=FALSE)
  if (any(as.numeric(Pathogen$IntroductionProb)!=0))
    stop("Active RK v0.5 requires Pathogen IntroductionProb=0; discrete introductions are reserved for a later gate",call.=FALSE)
  P2 <- .inatmlu_p_prepare_combined(Pathogen,nn,L,S,PathogenLandUseMixing)
  pf <- .inatmlu_p_flatten_tm(P)
  U <- nc*S; M<-matrix(0,U,length(states),dimnames=list(NULL,states))
  for(i in seq_len(nc)) for(s in seq_len(S)) M[(i-1L)*S+s,] <- pf[i,s,]
  getm <- function(nm,prob=FALSE) {
    z <- .iptm_resolve(P2[[nm]],as.integer(timestep),nc,S,as.integer(Ntimesteps),nm,prob)
    as.numeric(t(z))
  }
  beta<-getm("Beta"); rec<-getm("RecoveryProb",TRUE); mort<-getm("PathogenMortalityProb",TRUE)
  prog<-getm("ProgressionProb",TRUE); wan<-getm("ImmunityLossProb",TRUE); ds<-getm("DensityScale")
  ex <- .inatmlu_rk_competing_exit_hazards(rec,mort,dt); hrec<-ex$a; hmort<-ex$b
  hprog <- .inatmlu_rk_prob_to_hazard(prog,dt); hwan <- .inatmlu_rk_prob_to_hazard(wan,dt)
  if (is.null(StageMixing)) StageMixing <- matrix(1,S,S)
  if (!is.matrix(StageMixing) || !identical(dim(StageMixing),c(S,S)) || any(!is.finite(StageMixing)) || any(StageMixing<0)) stop("StageMixing must be a finite non-negative host-stage x host-stage matrix",call.=FALSE)
  Cnode <- if(is.null(P2$ContactMatrix)) diag(nc) else P2$ContactMatrix
  C <- kronecker(Cnode,StageMixing)
  S0<-M[,"S"]; I0<-M[,"I"]; E0<-if("E"%in%states)M[,"E"] else rep(0,U); R0<-if("R"%in%states)M[,"R"] else rep(0,U)
  live<-rowSums(M); infp<-as.numeric(crossprod(I0,C))
  if(Pathogen$Transmission=="frequency") { den<-as.numeric(crossprod(live,C)); force<-beta*ifelse(den>0,infp/den,0) } else force<-beta*infp/ds
  force <- pmax(0,force); inf_rate<-force*S0
  dM<-matrix(0,U,length(states),dimnames=list(NULL,states))
  dM[,"S"] <- dM[,"S"] - inf_rate
  if("E"%in%states) { dM[,"E"]<-dM[,"E"]+inf_rate-hprog*E0; dM[,"I"]<-dM[,"I"]+hprog*E0 } else dM[,"I"]<-dM[,"I"]+inf_rate
  dM[,"I"] <- dM[,"I"] - (hrec+hmort)*I0
  if("R"%in%states) { dM[,"R"]<-dM[,"R"]+hrec*I0-hwan*R0; dM[,"S"]<-dM[,"S"]+hwan*R0 } else dM[,"S"]<-dM[,"S"]+hrec*I0
  dpf<-array(0,dim=dim(pf),dimnames=dimnames(pf)); for(i in seq_len(nc)) for(s in seq_len(S)) dpf[i,s,]<-dM[(i-1L)*S+s,]
  dP <- .inatmlu_p_unflatten_tm(dpf,nn,L,S,states,P)
  list(dP=dP,incidence=inf_rate,deaths=hmort*I0)
}

.inatmlu_rk_host_rate_to_p <- function(rate,H,P) {
  if (is.null(P)) return(NULL)
  r <- as.numeric(rate); h<-as.numeric(H); dP<-array(0,dim=dim(P),dimnames=dimnames(P)); states<-dimnames(P)[[4L]]
  pos <- pmax(r,0); neg <- pmax(-r,0); hazard <- ifelse(h>0,neg/h,0)
  dP[,,,"S"] <- dP[,,,"S"] + array(pos,dim=dim(H))
  for(q in seq_along(states)) dP[,,,q] <- dP[,,,q] - P[,,,q]*array(hazard,dim=dim(H))
  dP
}

.inatmlu_rk_bc_terms <- function(H,P,Q,Biocontrol,Prepared,timestep,Ntimesteps) {
  if (is.null(Biocontrol)) return(list(dH=array(0,dim(H)),dP=if(is.null(P))NULL else array(0,dim(P)),dQ=NULL))
  d<-dim(H); nn<-d[1L]; L<-d[2L]; S<-d[3L]
  dH<-array(0,dim(H)); dP<-if(is.null(P))NULL else array(0,dim(P)); dQ<-lapply(Q,function(z)matrix(0,nrow(z),ncol(z),dimnames=dimnames(z))); names(dQ)<-names(Q)
  # Simultaneous attack hazards: host loss is driven by summed hazards; recruit
  # production for each agent is its own hazard contribution.
  for(nm in names(Q)) {
    a <- Prepared$OriginalBiocontrol$Agents[[nm]]
    ts <- .inatmlu_bc_base_target_stage(a,Prepared$stage_names,S)
    aa <- Prepared$Biocontrol$Agents[[nm]]
    rate <- .inabc_resolve_node(aa$AttackRate,timestep,Prepared$Context,"AttackRate",aa)
    active <- Q[[nm]][,aa$AttackStage]
    haz_node <- pmax(0,rate*active)
    recruit_node <- numeric(nn)
    for(l in seq_len(L)) for(s in ts) {
      hz <- haz_node; loss <- hz*H[,l,s]
      dH[,l,s] <- dH[,l,s]-loss
      if(!is.null(P)) for(q in seq_len(dim(P)[4L])) dP[,l,s,q] <- dP[,l,s,q]-hz*P[,l,s,q]
      recruit_node <- recruit_node + loss*aa$OffspringPerAttack
    }
    dQ[[nm]][,aa$RecruitStage] <- dQ[[nm]][,aa$RecruitStage]+recruit_node
  }
  list(dH=dH,dP=dP,dQ=dQ)
}

# Integrate local continuous biology only; no pest stage progression or spatial
# boundary occurs here. Returns both continuous endpoint and realised integers.
.INApestMLUTMContinuousLocalStep_v052 <- function(
    n, ContinuousBiology, Pathogen=NULL, PathogenState=NULL,
    StageMixing=NULL, PathogenLandUseMixing=NULL,
    Biocontrol=NULL, BiocontrolState=NULL, BiocontrolPrepared=NULL,
    timestep=1L, Ntimesteps=max(1L,as.integer(timestep)), Realise=TRUE) {
  .inatmlu_rk_require()
  if(!inherits(ContinuousBiology,"INApestContinuousBiology")) stop("ContinuousBiology must be created by INApestMLUTMContinuousBiology()",call.=FALSE)
  dims<-.inatmlu_validate_array_state(n); nn<-dims[1L]; L<-dims[2L]; S<-dims[3L]
  if(!is.null(Pathogen)) {
    if(is.null(PathogenState)) PathogenState<-INApestMLUTMPathogenInitial(n,Pathogen,NULL,Ntimesteps)
    PathogenState<-.inatmlu_p_validate_state(PathogenState,n,Pathogen)
  } else PathogenState<-NULL
  Prepared<-BiocontrolPrepared; Q<-BiocontrolState
  if(!is.null(Biocontrol)) {
    if(is.null(Prepared)) Prepared<-INApestMLUTMPrepareBiocontrol(Biocontrol,n,Ntimesteps)
    .inatmlu_rk_validate_bc(Biocontrol,Prepared,timestep)
    if(is.null(Q)) Q<-INApestMLUTMBiocontrolInitial(Biocontrol,n,Ntimesteps,Prepared)
    Q<-.inatmlu_rk_add_releases(Q,Prepared,timestep)
  } else Q<-NULL
  packed<-.inatmlu_rk_pack(n,PathogenState,Q)
  context<-list(n_nodes=nn,n_landuses=L,n_stages=S,timestep=as.integer(timestep),Ntimesteps=as.integer(Ntimesteps),Pathogen=Pathogen,Biocontrol=Biocontrol,PreparedBiocontrol=Prepared)
  deriv<-function(tt,y) {
    st<-.inatmlu_rk_unpack(y,packed$meta,n,PathogenState,Q); H<-st$H; P<-st$P; QQ<-st$Q
    host_rate<-array(0,dim(H))
    if(!is.null(ContinuousBiology$HostRateFunction)) {
      z<-.inatmlu_rk_call(ContinuousBiology$HostRateFunction,list(t=tt,H=H,P=P,Q=QQ,context=context))
      if(is.null(dim(z)) || !identical(dim(z),dim(H)) || any(!is.finite(z))) stop("HostRateFunction must return a finite array with the same dimensions as H",call.=FALSE)
      host_rate<-z
    }
    if(is.null(P)) dH<-host_rate else dH<-NULL
    dP<-if(is.null(P))NULL else .inatmlu_rk_host_rate_to_p(host_rate,H,P)
    if(!is.null(P)) {
      pp<-.inatmlu_rk_pathogen_derivative(H,P,Pathogen,timestep,Ntimesteps,StageMixing,PathogenLandUseMixing,ContinuousBiology$TimestepLength)
      dP<-dP+pp$dP
    }
    bct<-.inatmlu_rk_bc_terms(H,P,QQ,Biocontrol,Prepared,timestep,Ntimesteps)
    if(is.null(P)) dH<-dH+bct$dH else dP<-dP+bct$dP
    dQ<-bct$dQ
    if(!is.null(QQ) && !is.null(ContinuousBiology$BiocontrolRateFunction)) {
      extra<-.inatmlu_rk_call(ContinuousBiology$BiocontrolRateFunction,list(t=tt,H=H,P=P,Q=QQ,context=context))
      if(!is.list(extra) || !identical(names(extra),names(QQ))) stop("BiocontrolRateFunction must return a named list matching Q",call.=FALSE)
      for(nm in names(QQ)) { if(!identical(dim(extra[[nm]]),dim(QQ[[nm]]))||any(!is.finite(extra[[nm]]))) stop("BiocontrolRateFunction returned invalid dimensions/values for agent ",nm,call.=FALSE); dQ[[nm]]<-dQ[[nm]]+extra[[nm]] }
    }
    pp2<-.inatmlu_rk_pack(if(is.null(P))dH else H,dP,dQ)
    pp2$state
  }
  y1<-.INApestMLUTMRK4Integrate(packed$state,deriv,t0=0,dt=ContinuousBiology$TimestepLength,RKMaxStep=ContinuousBiology$RKMaxStep,NegativeTolerance=ContinuousBiology$NegativeTolerance)
  cont<-.inatmlu_rk_unpack(y1,packed$meta,n,PathogenState,Q)
  if(any(cont$H < -ContinuousBiology$NegativeTolerance) || (!is.null(cont$P)&&any(cont$P < -ContinuousBiology$NegativeTolerance)) || (!is.null(cont$Q)&&any(unlist(cont$Q)< -ContinuousBiology$NegativeTolerance))) stop("Continuous endpoint contains negative abundance",call.=FALSE)
  realised <- if (isTRUE(Realise)) .inatmlu_rk_bridge(cont$H,cont$P,cont$Q,ContinuousBiology$NegativeTolerance) else cont
  list(N=realised$H,PathogenState=realised$P,BiocontrolState=realised$Q,
       Continuous=cont,BiocontrolPrepared=Prepared)
}

# Apply the frozen discrete host boundary while suppressing a second pathogen
# biological update. The zero-process pathogen clone still carries compartments
# through stage progression/movement and assigns new host recruits to S.
INApestMLUTMRKDiscreteHostBoundary <- function(n,biology_args,Pathogen=NULL,PathogenState=NULL,timestep=1L,Ntimesteps=max(1L,as.integer(timestep)),StageMixing=NULL,PathogenLandUseMixing=NULL,managing=0) {
  if(!is.list(biology_args)) stop("biology_args must be a named list",call.=FALSE)
  if(length(biology_args) && (is.null(names(biology_args)) || any(!nzchar(names(biology_args))))) stop("biology_args must be a named list",call.=FALSE)
  if(is.null(Pathogen)) return(do.call(INApestMLUTMHostStep,c(list(n0=n,managing=managing),biology_args)))
  P0<-.inatmlu_rk_zero_pathogen(Pathogen)
  do.call(INApestMLUTMPathogenHostStep,c(list(n0=n,PathogenState=PathogenState,Pathogen=P0,timestep=timestep,Ntimesteps=Ntimesteps,StageMixing=StageMixing,PathogenLandUseMixing=PathogenLandUseMixing,managing=managing),biology_args))
}

INApestMLUTMRKResponseStep <- function(
    n,have_info,timestep=1L,last_known_presence=NULL,
    manage_prob=0,mortality_prob=0,detection_prob=0,info_triggered_detection_prob=0,
    info_retention_prob=1,info_persistence_steps=NA,seam=0,biology_args=list(),
    Pathogen=NULL,PathogenState=NULL,StageMixing=NULL,PathogenLandUseMixing=NULL,
    Biocontrol=NULL,BiocontrolState=NULL,Ntimesteps=max(1L,as.integer(timestep)),BiocontrolPrepared=NULL,
    ContinuousBiology=NULL) {
  if(is.null(ContinuousBiology)) return(INApestMLUTMPathogenBiocontrolResponseStep(
    n=n,have_info=have_info,timestep=timestep,last_known_presence=last_known_presence,
    manage_prob=manage_prob,mortality_prob=mortality_prob,detection_prob=detection_prob,
    info_triggered_detection_prob=info_triggered_detection_prob,info_retention_prob=info_retention_prob,
    info_persistence_steps=info_persistence_steps,seam=seam,biology_args=biology_args,
    Pathogen=Pathogen,PathogenState=PathogenState,StageMixing=StageMixing,PathogenLandUseMixing=PathogenLandUseMixing,
    Biocontrol=Biocontrol,BiocontrolState=BiocontrolState,Ntimesteps=Ntimesteps,BiocontrolPrepared=BiocontrolPrepared))
  .inatmlu_rk_require(); dims<-.inatmlu_validate_array_state(n); nn<-dims[1L];L<-dims[2L];S<-dims[3L]
  have_info<-.inatmlu_binary_info(have_info,nn,"have_info"); if(is.null(last_known_presence))last_known_presence<-rep(NA_real_,nn); last_known_presence<-as.numeric(last_known_presence)
  if(length(last_known_presence)!=nn || any(!is.na(last_known_presence)&!is.finite(last_known_presence))) stop("last_known_presence must be NULL or length nodes",call.=FALSE)
  retention<-.inatmlu_node_prob(info_retention_prob,nn,"info_retention_prob"); persistence<-.inatmlu_persistence(info_persistence_steps,nn); use_persistence<-any(!is.na(persistence))
  if(!is.list(biology_args)||(length(biology_args)&&(is.null(names(biology_args))||any(!nzchar(names(biology_args)))))) stop("biology_args must be a named list",call.=FALSE)
  if(anyDuplicated(names(biology_args))) stop("biology_args names must be unique",call.=FALSE)
  reserved<-intersect(names(biology_args),c("n0","managing","PathogenState","Pathogen","timestep","Ntimesteps","StageMixing","PathogenLandUseMixing"));if(length(reserved))stop("biology_args may not override RK/response core arguments: ",paste(reserved,collapse=", "),call.=FALSE)
  manage_surface<-.inatmlu_response_surface(manage_prob,nn,L,"manage_prob"); mortality_cube<-.inatmlu_response_cube(mortality_prob,nn,L,S,"mortality_prob")
  detection_cube<-.inatmlu_response_cube(detection_prob,nn,L,S,"detection_prob"); info_detection_cube<-.inatmlu_response_cube(info_triggered_detection_prob,nn,L,S,"info_triggered_detection_prob")
  node_invaded_start<-as.integer(vapply(seq_len(nn),function(i)sum(n[i,,])>0,logical(1))); known_occupied_start<-node_invaded_start*have_info
  manage_p<-as.numeric(manage_surface)*rep(have_info,times=L); managing<-matrix(rbinom(nn*L,1L,manage_p),nn,L); managing_cube<-array(rep(as.numeric(managing),times=S),c(nn,L,S))
  n_after_management<-array(rbinom(length(n),size=as.integer(n),prob=as.numeric(1-mortality_cube*managing_cube)),dim=dim(n),dimnames=dimnames(n)); management_deaths<-n-n_after_management
  if(use_persistence){killed<-vapply(seq_len(nn),function(i)sum(management_deaths[i,,])>0,logical(1));last_known_presence[killed]<-timestep}
  if(!is.null(Pathogen)){if(is.null(PathogenState))PathogenState<-INApestMLUTMPathogenInitial(n,Pathogen,NULL,Ntimesteps);p_after_management<-INApestMLUTMPathogenReconcile(PathogenState,n_after_management,Pathogen)}else p_after_management<-NULL
  local<-INApestMLUTMContinuousLocalStep(n_after_management,ContinuousBiology,Pathogen,p_after_management,StageMixing,PathogenLandUseMixing,Biocontrol,BiocontrolState,BiocontrolPrepared,timestep,Ntimesteps,TRUE)
  boundary<-INApestMLUTMRKDiscreteHostBoundary(local$N,biology_args,Pathogen,local$PathogenState,timestep,Ntimesteps,StageMixing,PathogenLandUseMixing,managing)
  if(is.null(Pathogen)){n_after_biology<-boundary;p_after_biology<-NULL}else{n_after_biology<-boundary$N;p_after_biology<-boundary$PathogenState}
  q_after_biology<-.inatmlu_rk_move_q(local$BiocontrolState,local$BiocontrolPrepared,timestep)
  programmed<-which(have_info==1L&!is.na(persistence));if(length(programmed)){elapsed<-timestep-last_known_presence;stop_nodes<-programmed[is.na(last_known_presence[programmed])|elapsed[programmed]>=persistence[programmed]];if(length(stop_nodes))have_info[stop_nodes]<-0L}
  decay<-which(have_info==1L&is.na(persistence)&retention<1);if(length(decay))have_info[decay]<-rbinom(length(decay),1L,retention[decay])
  seam_transferred<-integer(nn);if(!(length(seam)==1L&&identical(as.numeric(seam),0))){if(!is.matrix(seam)||!identical(dim(seam),c(nn,nn)))stop("seam must be 0 or a nodes x nodes matrix",call.=FALSE);if(any(!is.finite(seam))||any(seam<0|seam>1))stop("seam entries must be probabilities in [0,1]",call.=FALSE);Sx<-seam;diag(Sx)<-0;draws<-matrix(rbinom(nn*nn,1L,as.numeric(Sx*known_occupied_start)),nn,nn);seam_transferred<-as.integer(colSums(draws)>0);have_info[have_info==0L]<-seam_transferred[have_info==0L]}
  pathogen_detected<-integer(nn);p_path_det<-rep(0,nn)
  if(!is.null(Pathogen)){p_path_det<-INApestMLUTMPathogenDetectionProbability(p_after_biology,Pathogen$DetectionProb);pathogen_detected<-rbinom(nn,1L,p_path_det);if(isTRUE(Pathogen$DetectionTriggersInfo)){if(use_persistence)last_known_presence[pathogen_detected==1L]<-timestep;have_info[have_info==0L]<-pathogen_detected[have_info==0L]}}
  info_before_surveillance<-as.integer(have_info!=0L);p_background<-INApestMLUTMHostDetectionProbability(n_after_biology,detection_cube);background_detected<-rbinom(nn,1L,p_background);p_info<-INApestMLUTMHostDetectionProbability(n_after_biology,info_detection_cube);info_triggered_detected<-if(any(info_detection_cube>0))rbinom(nn,1L,p_info*info_before_surveillance)else integer(nn);host_detection_evidence<-pmax(background_detected,info_triggered_detected);if(use_persistence)last_known_presence[host_detection_evidence==1L]<-timestep;have_info[have_info==0L]<-host_detection_evidence[have_info==0L]
  invaded_by_landuse<-matrix(0L,nn,L);for(i in seq_len(nn))for(l in seq_len(L))invaded_by_landuse[i,l]<-as.integer(sum(n_after_biology[i,l,])>0);node_invaded_end<-as.integer(rowSums(invaded_by_landuse)>0)
  list(N=n_after_biology,NAfterManagement=n_after_management,NAfterContinuousLocal=local$N,NAfterBiology=n_after_biology,ManagementDeaths=management_deaths,Managing=managing,
       PathogenState=p_after_biology,PathogenStateAfterContinuousLocal=local$PathogenState,PathogenDeaths=if(is.null(Pathogen))NULL else array(NA_real_,dim(n_after_biology)),NewInfections=if(is.null(Pathogen))NULL else array(NA_real_,dim(n_after_biology)),
       PathogenDetectionProbability=p_path_det,PathogenDetected=pathogen_detected,HaveInfo=have_info,LastKnownPresence=last_known_presence,InformationStateBeforeSurveillance=info_before_surveillance,
       BackgroundDetectionProbability=p_background,InfoTriggeredDetectionProbability=p_info,BackgroundDetected=background_detected,InfoTriggeredDetected=info_triggered_detected,HostDetectionEvidence=host_detection_evidence,SEAMTransferred=seam_transferred,
       InvadedByLandUse=invaded_by_landuse,Invaded=node_invaded_end,KnownPresentByLandUse=invaded_by_landuse*have_info,
       BiocontrolState=q_after_biology,BiocontrolImpact=NULL,BiocontrolPrepared=local$BiocontrolPrepared,ContinuousLocalEndpoint=local$Continuous,RKMode=TRUE)
}


###############################################################################
### INApest MLUTM RK linked stochastic bridge v0.5.3
### Date: 2026-10-01
###
### Narrow strengthening of frozen MLUTM RK v0.5.2:
###   * Realise=FALSE retains the v0.5.2 deterministic RK architecture path.
###   * Realise=TRUE no longer independently rounds fractional RK endpoints.
###   * linked whole-count transfers/exits use the validated stochastic
###     compartment bridge used by the main INApest RK architecture chain.
###   * standard pathogen transfers and competing pathogen/biocontrol deaths
###     are realised as linked events; each realised biocontrol attack removes
###     exactly one host and creates the configured integer agent recruits.
###   * optional custom host stochastic biology must be supplied explicitly as
###     gross gains + per-capita loss hazards (HostFluxFunction). A non-zero net
###     HostRateFunction is intentionally not guessed into a stochastic process.
###   * optional continuous biocontrol stage demography must be supplied as
###     explicit transition/exit rates (BiocontrolRates).
###
### Pest demographic stage progression, transition movement, fecundity,
### dispersal/recruitment, biocontrol spatial movement, response and surveillance
### remain owned by the unchanged v0.5.2 / v0.4.2 architecture boundary.
###############################################################################

# Explicit continuous agent-rate contract, identical in meaning to the validated
# Meta/TM RK production chain. TransitionRates orientation follows the existing
# discrete agent Transition matrix: rows=destination, columns=source.
INApestMLUTMContinuousBiocontrolAgentRates <- function(TransitionRates=NULL, ExitRates=0) {
  structure(list(TransitionRates=TransitionRates, ExitRates=ExitRates),
            class=c("INApestMLUTMContinuousBiocontrolAgentRates","list"))
}
INApestMLUTMContinuousBiocontrolRates <- function(AgentRates) {
  if (inherits(AgentRates,"INApestMLUTMContinuousBiocontrolAgentRates"))
    stop("AgentRates must be a named list, one entry per biocontrol agent", call.=FALSE)
  if (!is.list(AgentRates) || !length(AgentRates) || is.null(names(AgentRates)) ||
      any(!nzchar(names(AgentRates))) || anyDuplicated(names(AgentRates)))
    stop("AgentRates must be a uniquely named non-empty list", call.=FALSE)
  if (!all(vapply(AgentRates,inherits,logical(1),"INApestMLUTMContinuousBiocontrolAgentRates")))
    stop("Every AgentRates entry must come from INApestMLUTMContinuousBiocontrolAgentRates()", call.=FALSE)
  structure(AgentRates,class=c("INApestMLUTMContinuousBiocontrolRates","list"))
}

.inatmlu_rk53_call <- function(f,args) {
  fm <- names(formals(f)); if(!is.null(fm) && !("..."%in%fm)) args<-args[intersect(names(args),fm)]
  do.call(f,args)
}

.inatmlu_rk53_expand_H <- function(x,H,label) {
  if (is.null(dim(x))) {
    z<-as.numeric(x)
    if(length(z)==1L) return(array(z,dim=dim(H),dimnames=dimnames(H)))
    if(length(z)==length(H)) return(array(z,dim=dim(H),dimnames=dimnames(H)))
    stop(label," must be scalar or shaped like H",call.=FALSE)
  }
  if(!identical(dim(x),dim(H))) stop(label," must have the same dimensions as H",call.=FALSE)
  out<-array(as.numeric(x),dim=dim(H),dimnames=dimnames(H))
  if(any(!is.finite(out))) stop(label," must be finite",call.=FALSE)
  out
}

.inatmlu_rk53_host_flux <- function(fun,t,H,P,Q,context) {
  if(is.null(fun)) return(list(GainRate=array(0,dim(H),dimnames=dimnames(H)),Hazards=list()))
  ans<-.inatmlu_rk53_call(fun,list(t=t,time=t,H=H,P=P,Q=Q,context=context))
  if(!is.list(ans)) stop("HostFluxFunction must return a list with GainRate and Hazards",call.=FALSE)
  gain<-ans$GainRate; if(is.null(gain)) gain<-ans$Gains
  if(is.null(gain)) stop("HostFluxFunction must return GainRate",call.=FALSE)
  gain<-.inatmlu_rk53_expand_H(gain,H,"Host GainRate")
  if(any(gain<0)) stop("Host GainRate must be non-negative",call.=FALSE)
  hz<-ans$Hazards
  if(is.null(hz)) hz<-list()
  if(!is.list(hz) || is.data.frame(hz)) stop("Host Hazards must be a named list of scalar or H-shaped per-capita hazards",call.=FALSE)
  if(length(hz)) {
    if(is.null(names(hz))||any(!nzchar(names(hz)))||anyDuplicated(names(hz)))
      stop("Host Hazards must have unique non-empty names",call.=FALSE)
    hz<-lapply(seq_along(hz),function(j){z<-.inatmlu_rk53_expand_H(hz[[j]],H,paste0("Host Hazards$",names(hz)[j]));if(any(z<0))stop("Host hazards must be non-negative",call.=FALSE);z}) |>
      setNames(names(ans$Hazards))
  }
  list(GainRate=gain,Hazards=hz)
}

# Re-define the public constructor additively. Positional compatibility with
# v0.5.2 is retained because new arguments are appended.
INApestMLUTMContinuousBiology <- function(HostRateFunction=NULL,
                                     BiocontrolRateFunction=NULL,
                                     RKMaxStep=0.025,
                                     TimestepLength=1,
                                     NegativeTolerance=1e-10,
                                     HostFluxFunction=NULL,
                                     BiocontrolRates=NULL) {
  if(!is.null(HostRateFunction)&&!is.function(HostRateFunction))stop("HostRateFunction must be NULL or a function",call.=FALSE)
  if(!is.null(BiocontrolRateFunction)&&!is.function(BiocontrolRateFunction))stop("BiocontrolRateFunction must be NULL or a function",call.=FALSE)
  if(!is.null(HostFluxFunction)&&!is.function(HostFluxFunction))stop("HostFluxFunction must be NULL or a function",call.=FALSE)
  if(!is.null(BiocontrolRates)&&!inherits(BiocontrolRates,"INApestMLUTMContinuousBiocontrolRates"))stop("BiocontrolRates must be NULL or created by INApestMLUTMContinuousBiocontrolRates()",call.=FALSE)
  if(!is.null(HostRateFunction)&&!is.null(HostFluxFunction))stop("Supply HostRateFunction or HostFluxFunction, not both",call.=FALSE)
  if(!is.null(BiocontrolRateFunction)&&!is.null(BiocontrolRates))stop("Supply BiocontrolRateFunction or BiocontrolRates, not both",call.=FALSE)
  if(!is.finite(RKMaxStep)||RKMaxStep<=0)stop("RKMaxStep must be finite and > 0",call.=FALSE)
  if(!is.finite(TimestepLength)||TimestepLength<=0)stop("TimestepLength must be finite and > 0",call.=FALSE)

  # Auto-generate deterministic mean functions from explicit stochastic
  # decompositions so Realise=FALSE and Realise=TRUE describe the same biology.
  hrate<-HostRateFunction
  if(!is.null(HostFluxFunction)) {
    ff<-HostFluxFunction
    hrate<-function(t,H,P=NULL,Q=NULL,context=NULL,...) {
      fl<-.inatmlu_rk53_host_flux(ff,t,H,P,Q,context)
      total<-array(0,dim(H),dimnames=dimnames(H))
      if(length(fl$Hazards)) for(z in fl$Hazards) total<-total+z
      fl$GainRate-H*total
    }
  }
  brate<-BiocontrolRateFunction
  if(!is.null(BiocontrolRates)) {
    specs<-BiocontrolRates
    brate<-function(t,H,P=NULL,Q,context,...) {
      if(is.null(context$PreparedBiocontrol)) stop("BiocontrolRates require PreparedBiocontrol in context",call.=FALSE)
      agents<-context$PreparedBiocontrol$OriginalBiocontrol$Agents
      if(!identical(names(specs),names(agents)))stop("BiocontrolRates names must match biocontrol agents",call.=FALSE)
      out<-Q
      for(nm in names(Q)) {
        rf<-.inatmlu_rk53_agent_rate_function(specs[[nm]],agents[[nm]],context$timestep,context)
        out[[nm]]<-INApestMLUTMCompartmentMeanDerivative(t,as.matrix(Q[[nm]]),rf)
      }
      out
    }
  }
  structure(list(HostRateFunction=hrate,BiocontrolRateFunction=brate,
                 RKMaxStep=as.numeric(RKMaxStep),TimestepLength=as.numeric(TimestepLength),
                 NegativeTolerance=as.numeric(NegativeTolerance),
                 HostFluxFunction=HostFluxFunction,BiocontrolRates=BiocontrolRates,
                 StochasticBridge="linked-compartment-v0.5.3"),
            class=c("INApestMLUTMContinuousBiology","INApestContinuousBiology"))
}

.inatmlu_rk53_flat_H <- function(H) as.numeric(t(.inatmlu_flatten_state(H)))
.inatmlu_rk53_unflat_H <- function(x,nn,L,S,template=NULL) {
  m<-matrix(as.numeric(x),nrow=nn*L,ncol=S,byrow=TRUE)
  out<-.inatmlu_unflatten_state(m,nn,L,S)
  if(!is.null(template)) dimnames(out)<-dimnames(template)
  out
}
.inatmlu_rk53_flat_P <- function(P) {
  d<-dim(P);nn<-d[1L];L<-d[2L];S<-d[3L];K<-d[4L];pf<-.inatmlu_p_flatten_tm(P)
  M<-matrix(0,nn*L*S,K,dimnames=list(NULL,dimnames(P)[[4L]]))
  for(i in seq_len(nn*L))for(s in seq_len(S))M[(i-1L)*S+s,]<-pf[i,s,]
  M
}
.inatmlu_rk53_unflat_P <- function(M,nn,L,S,states,template=NULL) {
  pf<-array(0,c(nn*L,S,length(states)),dimnames=list(NULL,NULL,states))
  for(i in seq_len(nn*L))for(s in seq_len(S))pf[i,s,]<-M[(i-1L)*S+s,]
  out<-.inatmlu_p_unflatten_tm(pf,nn,L,S,states,template)
  out
}
.inatmlu_rk53_rowmap <- function(nn,L,S) {
  u<-nn*L*S; idx<-seq_len(u); pseudo<-(idx-1L)%/%S+1L
  data.frame(node=(pseudo-1L)%/%L+1L,landuse=(pseudo-1L)%%L+1L,stage=(idx-1L)%%S+1L)
}

.inatmlu_rk53_biocontrol_call <- function(x,t,timestep,context,agent,State=NULL) {
  if(!is.function(x))return(x)
  .inatmlu_rk53_call(x,list(t=t,time=t,timestep=timestep,context=context,agent=agent,state=State,State=State))
}
.inatmlu_rk53_agent_transition_rates <- function(x,t,timestep,context,agent,State) {
  s<-length(agent$Stages);if(is.null(x))return(matrix(0,s,s,dimnames=list(agent$Stages,agent$Stages)))
  z<-.inatmlu_rk53_biocontrol_call(x,t,timestep,context,agent,State)
  if(length(dim(z))!=2L||!identical(dim(z),c(s,s)))stop("Continuous agent TransitionRates must be stages x stages",call.=FALSE)
  out<-as.matrix(z);if(any(!is.finite(out))||any(out<0))stop("Continuous agent TransitionRates must be finite and non-negative",call.=FALSE);diag(out)<-0;dimnames(out)<-list(agent$Stages,agent$Stages);out
}
.inatmlu_rk53_agent_exit_rates <- function(x,t,timestep,context,agent,State) {
  n<-context$n_nodes;s<-length(agent$Stages);z<-.inatmlu_rk53_biocontrol_call(x,t,timestep,context,agent,State);d<-dim(z)
  if(is.null(d)){v<-as.numeric(z);if(length(v)==1L)out<-matrix(v,n,s)else if(length(v)==s)out<-matrix(rep(v,each=n),n,s)else if(length(v)==n&&n!=s)out<-matrix(rep(v,s),n,s)else stop("Continuous agent ExitRates must be scalar, length stages, nodes x stages, or resolver",call.=FALSE)}
  else if(length(d)==2L&&identical(d,c(n,s)))out<-as.matrix(z)else stop("Continuous agent ExitRates must be scalar, length stages, nodes x stages, or resolver",call.=FALSE)
  if(any(!is.finite(out))||any(out<0))stop("Continuous agent ExitRates must be finite and non-negative",call.=FALSE);dimnames(out)<-list(rownames(State),agent$Stages);out
}
.inatmlu_rk53_agent_rate_function <- function(spec,agent,timestep,context) {
  force(spec);force(agent);force(timestep);force(context)
  function(t,State,...) {
    State<-as.matrix(State);n<-nrow(State);s<-ncol(State)
    tr<-.inatmlu_rk53_agent_transition_rates(spec$TransitionRates,t,timestep,context,agent,State)
    T<-array(0,c(n,s,s),dimnames=list(rownames(State),agent$Stages,agent$Stages))
    for(src in seq_len(s))for(dst in seq_len(s))if(src!=dst)T[,src,dst]<-tr[dst,src]
    ex<-.inatmlu_rk53_agent_exit_rates(spec$ExitRates,t,timestep,context,agent,State)
    E<-array(ex,c(n,s,1L),dimnames=list(rownames(State),agent$Stages,"agent_mortality"))
    list(GainRate=matrix(0,n,s,dimnames=dimnames(State)),TransitionHazards=T,ExitHazards=E)
  }
}
.inatmlu_rk53_zero_agent_rate_function <- function(agent) {
  force(agent);function(t,State,...){State<-as.matrix(State);n<-nrow(State);s<-ncol(State);list(GainRate=matrix(0,n,s,dimnames=dimnames(State)),TransitionHazards=array(0,c(n,s,s),dimnames=list(rownames(State),agent$Stages,agent$Stages)),ExitHazards=array(0,c(n,s,0L),dimnames=list(rownames(State),agent$Stages,character())))}
}
.inatmlu_rk53_attacker_exposure <- function(State,Time,Step,RateFunction,AttackStage) {
  State<-as.matrix(State);n<-nrow(State);s<-ncol(State);z0<-c(as.numeric(State),rep(0,n));AttackStage<-as.integer(AttackStage)
  der<-function(t,state,...){X<-matrix(state[seq_len(n*s)],n,s,dimnames=dimnames(State));dX<-INApestMLUTMCompartmentMeanDerivative(t,X,RateFunction);c(as.numeric(dX),as.numeric(X[,AttackStage]))}
  zend<-.inatmlu_scb_rk4_allow(z0,Time,Step,der);exposure<-as.numeric(zend[n*s+seq_len(n)]);if(any(!is.finite(exposure))||any(exposure< -1e-8))stop("Continuous biocontrol attacker exposure left valid domain",call.=FALSE);pmax(0,exposure)
}

.inatmlu_rk53_pathogen_rates <- function(State,H,Pathogen,timestep,Ntimesteps,StageMixing,PathogenLandUseMixing,dt,nn,L,S) {
  comps<-colnames(State);U<-nrow(State);k<-ncol(State)
  T<-array(0,c(U,k,k),dimnames=list(NULL,comps,comps));pmort<-rep(0,U)
  if(is.null(Pathogen))return(list(TransitionHazards=T,PathogenMortality=pmort))
  if(is.function(Pathogen$IntroductionProb)||any(as.numeric(Pathogen$IntroductionProb)!=0))stop("Active RK v0.5.3 requires Pathogen IntroductionProb=0",call.=FALSE)
  P2<-.inatmlu_p_prepare_combined(Pathogen,nn,L,S,PathogenLandUseMixing);nc<-nn*L
  getm<-function(nm,prob=FALSE){z<-.iptm_resolve(P2[[nm]],as.integer(timestep),nc,S,as.integer(Ntimesteps),nm,prob);as.numeric(t(z))}
  beta<-getm("Beta");rec<-getm("RecoveryProb",TRUE);mort<-getm("PathogenMortalityProb",TRUE);prog<-getm("ProgressionProb",TRUE);wan<-getm("ImmunityLossProb",TRUE);ds<-getm("DensityScale")
  ex<-.inatmlu_rk_competing_exit_hazards(rec,mort,dt);hrec<-ex$a;pmort<-ex$b;hprog<-.inatmlu_rk_prob_to_hazard(prog,dt);hwan<-.inatmlu_rk_prob_to_hazard(wan,dt)
  if(is.null(StageMixing))StageMixing<-matrix(1,S,S);if(!is.matrix(StageMixing)||!identical(dim(StageMixing),c(S,S))||any(!is.finite(StageMixing))||any(StageMixing<0))stop("StageMixing must be finite non-negative stage x stage",call.=FALSE)
  Cnode<-if(is.null(P2$ContactMatrix))diag(nc)else P2$ContactMatrix;C<-kronecker(Cnode,StageMixing);I<-State[,"I"];live<-rowSums(State);infp<-as.numeric(crossprod(I,C));if(Pathogen$Transmission=="frequency"){den<-as.numeric(crossprod(live,C));force<-beta*ifelse(den>0,infp/den,0)}else force<-beta*infp/ds;force<-pmax(0,force)
  if("E"%in%comps){T[,"S","E"]<-force;T[,"E","I"]<-hprog}else{T[,"S","I"]<-force}
  if(Pathogen$Model=="SIS"){T[,"I","S"]<-hrec;if(any(hwan>0))stop("ImmunityLossProb must be zero for SIS",call.=FALSE)}else{T[,"I","R"]<-hrec;T[,"R","S"]<-hwan}
  list(TransitionHazards=T,PathogenMortality=pmort)
}

.inatmlu_rk53_host_rate_function <- function(ContinuousBiology,Pathogen,StageMixing,PathogenLandUseMixing,Biocontrol,Prepared,attack_h_cell,timestep,Ntimesteps,nn,L,S,Q,context) {
  comps<-if(is.null(Pathogen))"H"else Pathogen$States;U<-nn*L*S;k<-length(comps);map<-.inatmlu_rk53_rowmap(nn,L,S)
  force(ContinuousBiology);force(Pathogen);force(StageMixing);force(PathogenLandUseMixing);force(Biocontrol);force(Prepared);force(attack_h_cell);force(timestep);force(Ntimesteps);force(nn);force(L);force(S);force(Q);force(context)
  function(t,State,...) {
    State<-as.matrix(State);Hvec<-rowSums(State);H<-.inatmlu_rk53_unflat_H(Hvec,nn,L,S)
    P<-if(is.null(Pathogen))NULL else .inatmlu_rk53_unflat_P(State,nn,L,S,Pathogen$States)
    fl<-.inatmlu_rk53_host_flux(ContinuousBiology$HostFluxFunction,t,H,P,Q,context)
    gvec<-.inatmlu_rk53_flat_H(fl$GainRate);gain<-matrix(0,U,k,dimnames=list(NULL,comps));gain[,if("S"%in%comps)"S"else 1L]<-gvec
    trans<-array(0,c(U,k,k),dimnames=list(NULL,comps,comps));pmort<-rep(0,U)
    if(!is.null(Pathogen)){pr<-.inatmlu_rk53_pathogen_rates(State,H,Pathogen,timestep,Ntimesteps,StageMixing,PathogenLandUseMixing,ContinuousBiology$TimestepLength,nn,L,S);trans<-pr$TransitionHazards;pmort<-pr$PathogenMortality}
    host_names<-names(fl$Hazards);bio_names<-if(is.null(Biocontrol))character()else names(Biocontrol$Agents);causes<-c(if(length(host_names))paste0("host:",host_names)else character(),if(!is.null(Pathogen))"pathogen"else character(),if(length(bio_names))paste0("bio:",bio_names)else character())
    exits<-array(0,c(U,k,length(causes)),dimnames=list(NULL,comps,causes))
    if(length(host_names))for(j in seq_along(host_names)){hv<-.inatmlu_rk53_flat_H(fl$Hazards[[j]]);exits[,,paste0("host:",host_names[j])]<-matrix(rep(hv,k),U,k)}
    if(!is.null(Pathogen))exits[,"I","pathogen"]<-pmort
    if(length(bio_names))for(j in seq_along(bio_names)){hv<-attack_h_cell[,j];exits[,,paste0("bio:",bio_names[j])]<-matrix(rep(hv,k),U,k)}
    list(GainRate=gain,TransitionHazards=trans,ExitHazards=exits)
  }
}

.inatmlu_rk53_check_legacy_rates <- function(ContinuousBiology,t,H,P,Q,context) {
  if(!is.null(ContinuousBiology$HostRateFunction)&&is.null(ContinuousBiology$HostFluxFunction)){
    z<-.inatmlu_rk_call(ContinuousBiology$HostRateFunction,list(t=t,H=H,P=P,Q=Q,context=context));if(any(abs(z)>1e-12))stop("Stochastic Realise=TRUE cannot infer births/deaths from a non-zero HostRateFunction. Supply HostFluxFunction with explicit GainRate and Hazards.",call.=FALSE)
  }
  if(!is.null(ContinuousBiology$BiocontrolRateFunction)&&is.null(ContinuousBiology$BiocontrolRates)){
    z<-.inatmlu_rk_call(ContinuousBiology$BiocontrolRateFunction,list(t=t,H=H,P=P,Q=Q,context=context));if(!is.list(z)||any(abs(unlist(z))>1e-12))stop("Stochastic Realise=TRUE cannot infer agent transfers/losses from a non-zero BiocontrolRateFunction. Supply BiocontrolRates with explicit TransitionRates and ExitRates.",call.=FALSE)
  }
  invisible(TRUE)
}

# The legacy independent endpoint-rounding helper is deliberately disabled in
# the strengthened production path. It remains defined in v0.5 for historical
# reproducibility, but calling it after this hotfix is an explicit error.
.inatmlu_rk_bridge_disabled <- function(...) stop("Independent endpoint stochastic rounding is disabled in MLUTM RK v0.5.3; use INApestMLUTMContinuousLocalStep() linked bridge",call.=FALSE)

INApestMLUTMContinuousLocalStep <- function(
    n,ContinuousBiology,Pathogen=NULL,PathogenState=NULL,
    StageMixing=NULL,PathogenLandUseMixing=NULL,
    Biocontrol=NULL,BiocontrolState=NULL,BiocontrolPrepared=NULL,
    timestep=1L,Ntimesteps=max(1L,as.integer(timestep)),Realise=TRUE) {
  .inatmlu_rk_require()
  if(!inherits(ContinuousBiology,"INApestMLUTMContinuousBiology"))stop("ContinuousBiology must be created by INApestMLUTMContinuousBiology()",call.=FALSE)
  dims<-.inatmlu_validate_array_state(n);nn<-dims[1L];L<-dims[2L];S<-dims[3L]

  # Preserve the exact v0.5.2 deterministic/reference trajectory. Explicit flux
  # constructors auto-generate their deterministic mean functions above.
  det<-.INApestMLUTMContinuousLocalStep_v052(n,ContinuousBiology,Pathogen,PathogenState,StageMixing,PathogenLandUseMixing,Biocontrol,BiocontrolState,BiocontrolPrepared,timestep,Ntimesteps,FALSE)
  if(!isTRUE(Realise))return(det)

  if(!is.null(Pathogen)){if(is.null(PathogenState))PathogenState<-INApestMLUTMPathogenInitial(n,Pathogen,NULL,Ntimesteps);PathogenState<-.inatmlu_p_validate_state(PathogenState,n,Pathogen)}else PathogenState<-NULL
  Prepared<-BiocontrolPrepared;Q<-BiocontrolState
  if(!is.null(Biocontrol)){
    if(is.null(Prepared))Prepared<-INApestMLUTMPrepareBiocontrol(Biocontrol,n,Ntimesteps);.inatmlu_rk_validate_bc(Biocontrol,Prepared,timestep)
    if(is.null(Q))Q<-INApestMLUTMBiocontrolInitial(Biocontrol,n,Ntimesteps,Prepared);Q<-.inatmlu_rk_add_releases(Q,Prepared,timestep)
  }else Q<-NULL
  context<-list(n_nodes=nn,n_landuses=L,n_stages=S,timestep=as.integer(timestep),Ntimesteps=as.integer(Ntimesteps),Pathogen=Pathogen,Biocontrol=Biocontrol,PreparedBiocontrol=Prepared)
  .inatmlu_rk53_check_legacy_rates(ContinuousBiology,0,n,PathogenState,Q,context)

  # Structural zero path preserves exact state and RNG identity.
  no_custom<-is.null(ContinuousBiology$HostFluxFunction)&&is.null(ContinuousBiology$HostRateFunction)&&is.null(ContinuousBiology$BiocontrolRates)&&is.null(ContinuousBiology$BiocontrolRateFunction)
  if(no_custom&&is.null(Pathogen)&&is.null(Biocontrol))return(list(N=n,PathogenState=NULL,BiocontrolState=NULL,Continuous=det$Continuous,BiocontrolPrepared=NULL,StochasticBridge=list(Scheme="linked-compartment-v0.5.3",NSubsteps=0L,Attacks=NULL)))

  Hstate<-if(is.null(Pathogen))matrix(as.integer(.inatmlu_rk53_flat_H(n)),ncol=1L,dimnames=list(NULL,"H"))else{M<-.inatmlu_rk53_flat_P(PathogenState);storage.mode(M)<-"integer";M}
  specs<-ContinuousBiology$BiocontrolRates
  if(!is.null(specs)&&is.null(Biocontrol))stop("BiocontrolRates supplied but Biocontrol is NULL",call.=FALSE)
  if(!is.null(specs)&&!identical(names(specs),names(Biocontrol$Agents)))stop("BiocontrolRates names must exactly match biocontrol agent names",call.=FALSE)
  attacks_node<-if(is.null(Biocontrol))NULL else setNames(lapply(names(Biocontrol$Agents),function(x)integer(nn)),names(Biocontrol$Agents))
  host_exits<-NULL;host_gains<-matrix(0L,nrow(Hstate),ncol(Hstate),dimnames=dimnames(Hstate))
  nsub<-max(1L,as.integer(ceiling(ContinuousBiology$TimestepLength/ContinuousBiology$RKMaxStep-1e-14)));h<-ContinuousBiology$TimestepLength/nsub;tt<-(as.integer(timestep)-1L)*ContinuousBiology$TimestepLength;map<-.inatmlu_rk53_rowmap(nn,L,S)

  for(ss in seq_len(nsub)){
    curH<-.inatmlu_rk53_unflat_H(rowSums(Hstate),nn,L,S,n);curP<-if(is.null(Pathogen))NULL else .inatmlu_rk53_unflat_P(Hstate,nn,L,S,Pathogen$States,PathogenState)
    .inatmlu_rk53_check_legacy_rates(ContinuousBiology,tt,curH,curP,Q,context)
    attack_h<-if(is.null(Biocontrol))matrix(numeric(),nrow(Hstate),0L)else matrix(0,nrow(Hstate),length(Biocontrol$Agents),dimnames=list(NULL,names(Biocontrol$Agents)))
    nextQ<-Q
    if(!is.null(Biocontrol))for(nm in names(Biocontrol$Agents)){
      a<-Prepared$OriginalBiocontrol$Agents[[nm]];spec<-if(is.null(specs))NULL else specs[[nm]];rf<-if(is.null(spec)).inatmlu_rk53_zero_agent_rate_function(a)else .inatmlu_rk53_agent_rate_function(spec,a,timestep,context)
      exposure<-.inatmlu_rk53_attacker_exposure(Q[[nm]],tt,h,rf,a$AttackStage);ar<-.inabc_resolve_node(a$AttackRate,timestep,Prepared$Context,"AttackRate",a);if(any(!is.finite(ar))||any(ar<0))stop("AttackRate must resolve finite non-negative values",call.=FALSE);haz_node<-ar*(exposure/h);targets<-.inatmlu_bc_base_target_stage(a,Prepared$stage_names,S)
      for(u in seq_len(nrow(Hstate)))if(map$stage[u]%in%targets)attack_h[u,nm]<-haz_node[map$node[u]]
      if(is.null(spec))nextQ[[nm]]<-Q[[nm]]else nextQ[[nm]]<-INApestMLUTMStochasticCompartmentStep(Q[[nm]],tt,h,rf)$State
    }
    hfun<-.inatmlu_rk53_host_rate_function(ContinuousBiology,Pathogen,StageMixing,PathogenLandUseMixing,Biocontrol,Prepared,attack_h,timestep,Ntimesteps,nn,L,S,Q,context)
    hz<-INApestMLUTMStochasticCompartmentStep(Hstate,tt,h,hfun);Hstate<-hz$State;host_gains<-host_gains+hz$Gains
    if(is.null(host_exits)){host_exits<-hz$Exits}else{allc<-union(colnames(host_exits),colnames(hz$Exits));expand<-function(x){z<-matrix(0,nrow(x),length(allc),dimnames=list(NULL,allc));if(ncol(x))z[,colnames(x)]<-x;z};host_exits<-expand(host_exits)+expand(hz$Exits)}
    if(!is.null(Biocontrol))for(nm in names(Biocontrol$Agents)){
      cause<-paste0("bio:",nm);killed_cell<-if(cause%in%colnames(hz$Exits))as.integer(hz$Exits[,cause])else integer(nrow(Hstate));killed_node<-integer(nn);for(u in seq_along(killed_cell))killed_node[map$node[u]]<-killed_node[map$node[u]]+killed_cell[u];attacks_node[[nm]]<-attacks_node[[nm]]+killed_node;rec<-killed_node*as.integer(Biocontrol$Agents[[nm]]$OffspringPerAttack);nextQ[[nm]][,Biocontrol$Agents[[nm]]$RecruitStage]<-nextQ[[nm]][,Biocontrol$Agents[[nm]]$RecruitStage]+rec;storage.mode(nextQ[[nm]])<-"integer"
    }
    Q<-nextQ;tt<-tt+h
  }
  Hout<-.inatmlu_rk53_unflat_H(rowSums(Hstate),nn,L,S,n);storage.mode(Hout)<-"integer";Pout<-if(is.null(Pathogen))NULL else .inatmlu_rk53_unflat_P(Hstate,nn,L,S,Pathogen$States,PathogenState);if(!is.null(Pout))storage.mode(Pout)<-"integer"
  if(!is.null(Pout)&&!all(apply(Pout,c(1,2,3),sum)==Hout))stop("Linked bridge internal P-to-H accounting failure",call.=FALSE)
  impact<-NULL;if(!is.null(Biocontrol)){recnode<-lapply(names(Biocontrol$Agents),function(nm)as.integer(attacks_node[[nm]]*Biocontrol$Agents[[nm]]$OffspringPerAttack));names(recnode)<-names(Biocontrol$Agents);impact<-list(AttacksByAgentNode=attacks_node,RecruitsByAgentNode=recnode)}
  list(N=Hout,PathogenState=Pout,BiocontrolState=Q,Continuous=det$Continuous,BiocontrolPrepared=Prepared,BiocontrolImpact=impact,
       StochasticBridge=list(Scheme="linked-compartment-v0.5.3",NSubsteps=nsub,InternalStep=h,HostGains=host_gains,HostExits=host_exits,Attacks=attacks_node))
}
