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

if (!exists("INApestRK4Step", mode = "function") ||
    !exists(".inapest_rk4_call_supported", mode = "function"))
  stop("Source INApestRK4.R before INApestStochasticCompartmentBridge.R")

.inapest_scb_assert_state <- function(State, label = "State") {
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

.inapest_scb_matrix <- function(x, State, label) {
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

.inapest_scb_transition_array <- function(x, State) {
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

.inapest_scb_exit_array <- function(x, State) {
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

.inapest_scb_rates <- function(fun, t, State, Parameters = NULL, Context = NULL,
                               expected_causes = NULL, extra = list()) {
  ans <- .inapest_rk4_call_supported(
    fun,
    c(list(t=t,time=t,state=State,State=State,pars=Parameters,Parameters=Parameters,
           context=Context,Context=Context), extra)
  )
  if (!is.list(ans))
    stop("CompartmentRateFunction must return a list")
  gain <- ans$GainRate
  if (is.null(gain)) gain <- State*0
  gain <- .inapest_scb_matrix(gain, State, "GainRate")
  if (any(gain < 0)) stop("GainRate must be non-negative")
  trans <- .inapest_scb_transition_array(ans$TransitionHazards, State)
  exits <- .inapest_scb_exit_array(ans$ExitHazards, State)
  causes <- dimnames(exits)[[3]]
  if (!is.null(expected_causes) && !identical(causes, expected_causes))
    stop("ExitHazards causes changed during an RK step")
  list(GainRate=gain, TransitionHazards=trans, ExitHazards=exits, Causes=causes)
}

.inapest_scb_generator <- function(Tmat, Emat) {
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
INApestCompartmentMeanDerivative <- function(t, State, CompartmentRateFunction,
                                              Parameters=NULL, Context=NULL, ...) {
  State <- as.matrix(State)
  fl <- .inapest_scb_rates(CompartmentRateFunction,t,State,Parameters,Context,
                           expected_causes=NULL,extra=list(...))
  out <- fl$GainRate*0
  for (i in seq_len(nrow(State))) {
    E <- if (length(fl$Causes)) fl$ExitHazards[i,,,drop=FALSE][1,,] else matrix(numeric(),ncol(State),0L)
    if (length(fl$Causes)==1L) E <- matrix(E,ncol(State),1L,dimnames=list(colnames(State),fl$Causes))
    Q <- .inapest_scb_generator(fl$TransitionHazards[i,,], E)
    out[i,] <- fl$GainRate[i,] + as.numeric(State[i,] %*% Q)
  }
  out
}

# RK4 moments over one stochastic substep.
INApestStochasticCompartmentMoments <- function(State, Time, Step,
                                                 CompartmentRateFunction,
                                                 Parameters=NULL, Context=NULL, ...) {
  State <- .inapest_scb_assert_state(State)
  if (!is.function(CompartmentRateFunction)) stop("CompartmentRateFunction must be a function")
  Time <- as.numeric(Time)[1L]; Step <- as.numeric(Step)[1L]
  if (!is.finite(Time)) stop("Time must be finite")
  if (!is.finite(Step) || Step <= 0) stop("Step must be finite and > 0")
  n <- nrow(State); k <- ncol(State); extra <- list(...)
  f0 <- .inapest_scb_rates(CompartmentRateFunction,Time,State,Parameters,Context,
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
    fl <- .inapest_scb_rates(CompartmentRateFunction,t,X,Parameters,Context,
                             expected_causes=causes,extra=extra)
    dX <- matrix(0,n,k); dP <- array(0,c(n,k,k)); dD <- array(0,c(n,k,cN))
    dB <- array(0,c(n,k,k)); dGD <- array(0,c(n,k,cN)); dGT <- fl$GainRate
    for (ii in seq_len(n)) {
      E <- if (cN) matrix(fl$ExitHazards[ii,,],k,cN) else matrix(numeric(),k,0L)
      Q <- .inapest_scb_generator(fl$TransitionHazards[ii,,],E)
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

  zend <- INApestRK4Step(z0,Time=Time,Step=Step,RateFunction=deriv,NonNegative="allow")
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

INApestStochasticCompartmentStep <- function(State, Time, Step,
                                              CompartmentRateFunction,
                                              Parameters=NULL, Context=NULL, ...) {
  State <- .inapest_scb_assert_state(State)
  mom <- INApestStochasticCompartmentMoments(State,Time,Step,CompartmentRateFunction,
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

INApestStochasticCompartmentIntegrate <- function(State, CompartmentRateFunction,
                                                   Duration, MaxStep=0.025,
                                                   StartTime=0, Parameters=NULL,
                                                   Context=NULL, ...) {
  State <- .inapest_scb_assert_state(State)
  Duration <- as.numeric(Duration)[1L]; MaxStep <- as.numeric(MaxStep)[1L]; StartTime <- as.numeric(StartTime)[1L]
  if (!is.finite(Duration) || Duration < 0) stop("Duration must be finite and >= 0")
  if (!is.finite(MaxStep) || MaxStep <= 0) stop("MaxStep must be finite and > 0")
  if (!is.finite(StartTime)) stop("StartTime must be finite")
  if (Duration == 0) return(list(State=State,Gains=State*0,Exits=matrix(numeric(),nrow(State),0L),NSubsteps=0L,InternalStep=0))
  nsub <- max(1L,as.integer(ceiling(Duration/MaxStep-1e-14)))
  h <- Duration/nsub; x <- State; gains <- State*0; exits <- NULL; causes <- NULL
  max_closure <- c(start=0,gains=0,state=0); tt <- StartTime
  for (jj in seq_len(nsub)) {
    z <- INApestStochasticCompartmentStep(x,tt,h,CompartmentRateFunction,Parameters,Context,...)
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
