###############################################################################
# INApestStochasticFluxBridge.R
# Generic stochastic integer-flux bridge for continuous INApest LocalDynamics
# v0.2 -- 29 September 2026
#
# Purpose
# -------
# Convert continuous-time local biological dynamics into stochastic whole-count
# state changes without rounding a fractional RK state.
#
# Scientific contract
# -------------------
# A FluxFunction must explicitly return:
#   GainRate : gross entry/birth/immigration rate, counts per unit time; and
#   Hazards  : one or more per-capita bounded-loss hazards, per unit time.
#
# The implied deterministic equation for state element i is
#
#   dx_i/dt = GainRate_i - x_i * sum_j Hazards_ij .
#
# This decomposition is required because a net derivative alone cannot identify
# the stochastic process: the same dN/dt can arise from many different birth /
# death combinations with different variances and extinction probabilities.
#
# Within each stochastic internal step, classical RK4 follows the deterministic
# environment and simultaneously integrates:
#   * fate probabilities for the cohort present at the start of the step; and
#   * Poisson-thinning means for gains occurring during the step.
# Initial individuals are realised with one multinomial draw across competing
# loss causes + survival. New gains are realised as independent Poisson-thinned
# fate categories. Therefore no individual can experience two bounded losses,
# all returned states remain whole counts, and there is no deterministic
# rounding of an RK endpoint.
#
# The approximation becomes a continuous-time stochastic process as the
# stochastic internal step is refined. RKMaxStep controls BOTH numerical and
# stochastic event-time resolution.
###############################################################################

if (!exists("INApestRK4Step", mode = "function"))
  stop("Source INApestRK4.R before INApestStochasticFluxBridge.R")

.inapest_sfb_restore <- function(x, template) {
  out <- as.numeric(x)
  dim(out) <- dim(template)
  dimnames(out) <- dimnames(template)
  if (is.null(dim(template))) names(out) <- names(template)
  out
}

.inapest_sfb_assert_whole_state <- function(x, label = "State") {
  if (!is.numeric(x) || any(!is.finite(x)))
    stop(label, " must contain finite numeric values")
  if (any(x < 0)) stop(label, " must contain non-negative whole counts")
  tol <- 1e-9 * pmax(1, abs(x))
  if (any(abs(x - round(x)) > tol))
    stop(label, " must contain whole counts; fractional RK states may not be rounded into production")
  .inapest_sfb_restore(round(x), x)
}

.inapest_sfb_expand <- function(x, state, label) {
  z <- as.numeric(x)
  if (!length(z) || any(!is.finite(z))) stop(label, " must be finite")
  if (length(z) == 1L) z <- rep(z, length(state))
  if (length(z) != length(state))
    stop(label, " must be scalar or have one value per state element")
  .inapest_sfb_restore(z, state)
}

.inapest_sfb_hazard_matrix <- function(x, state) {
  n <- length(state)
  if (is.null(x)) return(matrix(numeric(), nrow = n, ncol = 0L))

  if (is.list(x) && !is.data.frame(x)) {
    if (!length(x)) return(matrix(numeric(), nrow = n, ncol = 0L))
    nm <- names(x)
    if (is.null(nm) || any(!nzchar(nm))) nm <- paste0("loss", seq_along(x))
    out <- do.call(cbind, lapply(seq_along(x), function(i)
      as.numeric(.inapest_sfb_expand(x[[i]], state, paste0("Hazards$", nm[i])))))
    colnames(out) <- nm
    return(out)
  }

  if (is.matrix(x) || is.data.frame(x)) {
    out <- as.matrix(x)
    storage.mode(out) <- "double"
    if (nrow(out) == 1L && n > 1L) out <- out[rep(1L, n), , drop = FALSE]
    if (nrow(out) != n) stop("Hazards must have one row per state element")
    if (is.null(colnames(out))) colnames(out) <- paste0("loss", seq_len(ncol(out)))
    if (any(!nzchar(colnames(out))) || anyDuplicated(colnames(out)))
      stop("Hazard columns must have unique non-empty names")
    return(out)
  }

  out <- as.numeric(.inapest_sfb_expand(x, state, "Hazards"))
  matrix(out, nrow = n, ncol = 1L,
         dimnames = list(NULL, if (!is.null(names(x)) && length(names(x)) == 1L && nzchar(names(x))) names(x) else "loss"))
}

.inapest_sfb_flux <- function(fun, t, state, Parameters = NULL, Context = NULL,
                              expected_causes = NULL, extra = list()) {
  ans <- .inapest_rk4_call_supported(
    fun,
    c(list(t = t, time = t, state = state, State = state,
           pars = Parameters, Parameters = Parameters,
           context = Context, Context = Context), extra)
  )
  if (!is.list(ans))
    stop("FluxFunction must return a list with GainRate and Hazards")
  gain <- ans$GainRate
  if (is.null(gain)) gain <- ans$Gains
  if (is.null(gain)) stop("FluxFunction must return GainRate")
  gain <- .inapest_sfb_expand(gain, state, "GainRate")
  hazards <- .inapest_sfb_hazard_matrix(ans$Hazards, state)
  if (any(gain < 0)) stop("GainRate must be non-negative")
  if (any(!is.finite(hazards)) || any(hazards < 0))
    stop("Hazards must be finite and non-negative")
  if (!is.null(expected_causes) && !identical(colnames(hazards), expected_causes))
    stop("FluxFunction hazard names/number changed during an RK step")
  list(GainRate = gain, Hazards = hazards)
}

# Draw one whole-count cohort across mutually exclusive categories with supplied
# probabilities. The final implicit category is survival. Sequential conditional
# binomials are exactly equivalent to a multinomial draw and vectorise in base R.
INApestStochasticBoundedProbabilities <- function(Count, Probabilities) {
  Count <- .inapest_sfb_assert_whole_state(Count, "Count")
  CountVec <- as.numeric(Count)
  P <- as.matrix(Probabilities)
  if (nrow(P) == 1L && length(CountVec) > 1L)
    P <- P[rep(1L, length(CountVec)), , drop = FALSE]
  if (nrow(P) != length(CountVec))
    stop("Probabilities must have one row per Count value")
  if (any(!is.finite(P)) || any(P < -1e-12))
    stop("Probabilities must be finite and non-negative")
  P[P < 0 & P >= -1e-12] <- 0
  rs <- rowSums(P)
  if (any(rs > 1 + 1e-10))
    stop("Bounded event probabilities exceed one; reduce RKMaxStep or revise the flux model")
  if (any(rs > 1)) P <- P / pmax(1, rs)

  remaining_n <- CountVec
  remaining_p <- rep(1, length(CountVec))
  Draws <- matrix(0, nrow(P), ncol(P), dimnames = dimnames(P))
  for (j in seq_len(ncol(P))) {
    q <- ifelse(remaining_p > 1e-15, P[, j] / remaining_p, 0)
    q <- pmin(1, pmax(0, q))
    z <- stats::rbinom(length(CountVec), size = remaining_n, prob = q)
    Draws[, j] <- z
    remaining_n <- remaining_n - z
    remaining_p <- pmax(0, remaining_p - P[, j])
  }
  list(Events = Draws, Survivors = remaining_n,
       Probabilities = P, SurvivalProbability = pmax(0, 1 - rowSums(P)))
}

# RK4 moment calculation for one stochastic substep. The deterministic state x,
# start-cohort survival s, start-cohort cause probabilities c, surviving-gain
# mean y, new-gain cause means d, and total gain mean G are integrated together.
# For a fixed deterministic environment path, this is the exact multinomial /
# Poisson thinning decomposition of the endpoint population.
INApestStochasticFluxMoments <- function(State, Time, Step, FluxFunction,
                                         Parameters = NULL, Context = NULL, ...) {
  State <- .inapest_sfb_assert_whole_state(State)
  if (!is.function(FluxFunction)) stop("FluxFunction must be a function")
  Time <- as.numeric(Time)[1L]; Step <- as.numeric(Step)[1L]
  if (!is.finite(Time)) stop("Time must be finite")
  if (!is.finite(Step) || Step <= 0) stop("Step must be finite and > 0")
  extra <- list(...)
  n <- length(State)
  f0 <- .inapest_sfb_flux(FluxFunction, Time, State, Parameters, Context,
                          expected_causes = NULL, extra = extra)
  causes <- colnames(f0$Hazards)
  m <- length(causes)

  ix_x <- seq_len(n)
  ix_s <- n + seq_len(n)
  ix_y <- 2L*n + seq_len(n)
  ix_G <- 3L*n + seq_len(n)
  offset <- 4L*n
  ix_c <- if (m) offset + seq_len(n*m) else integer()
  offset <- offset + n*m
  ix_d <- if (m) offset + seq_len(n*m) else integer()

  z0 <- c(as.numeric(State), rep(1, n), rep(0, n), rep(0, n),
          if (m) rep(0, n*m) else numeric(),
          if (m) rep(0, n*m) else numeric())

  deriv <- function(t, state, ...) {
    x <- state[ix_x]
    s <- state[ix_s]
    y <- state[ix_y]
    xt <- .inapest_sfb_restore(x, State)
    fl <- .inapest_sfb_flux(FluxFunction, t, xt, Parameters, Context,
                            expected_causes = causes, extra = extra)
    g <- as.numeric(fl$GainRate)
    Hmat <- fl$Hazards
    H <- if (m) rowSums(Hmat) else rep(0, n)
    dx <- g - H * x
    ds <- -H * s
    dy <- g - H * y
    dG <- g
    dc <- if (m) as.numeric(Hmat * s) else numeric()
    dd <- if (m) as.numeric(Hmat * y) else numeric()
    c(dx, ds, dy, dG, dc, dd)
  }

  zend <- INApestRK4Step(z0, Time = Time, Step = Step,
                         RateFunction = deriv, NonNegative = "allow")
  x <- zend[ix_x]
  s <- zend[ix_s]
  y <- zend[ix_y]
  G <- zend[ix_G]
  C <- if (m) matrix(zend[ix_c], nrow = n, ncol = m,
                     dimnames = list(NULL, causes)) else matrix(numeric(), n, 0L)
  D <- if (m) matrix(zend[ix_d], nrow = n, ncol = m,
                     dimnames = list(NULL, causes)) else matrix(numeric(), n, 0L)

  tol <- 2e-9
  bad <- any(!is.finite(c(x,s,y,G,C,D))) || any(x < -tol) ||
    any(s < -tol | s > 1 + tol) || any(y < -tol) || any(G < -tol) ||
    any(C < -tol) || any(D < -tol)
  if (bad)
    stop("RK stochastic moment step left the valid biological domain; reduce RKMaxStep or revise the flux model")
  s[s < 0] <- 0; s[s > 1] <- 1
  y[y < 0] <- 0; G[G < 0] <- 0
  if (m) { C[C < 0] <- 0; D[D < 0] <- 0 }

  initial_closure <- if (m) max(abs(s + rowSums(C) - 1)) else max(abs(s - 1))
  gain_closure <- if (m) max(abs(G - y - rowSums(D))) else max(abs(G - y))
  decomp <- as.numeric(State) * s + y
  state_closure <- max(abs(x - decomp))
  scale <- max(1, max(abs(c(x,G,as.numeric(State)))))
  if (initial_closure > 2e-8 || gain_closure > 2e-8*scale || state_closure > 2e-8*scale)
    stop("RK stochastic moment identities failed; reduce RKMaxStep")

  list(
    ExpectedEnd = .inapest_sfb_restore(x, State),
    SurvivalProbability = .inapest_sfb_restore(s, State),
    InitialLossProbabilities = C,
    SurvivingGainMean = .inapest_sfb_restore(y, State),
    NewGainLossMeans = D,
    TotalGainMean = .inapest_sfb_restore(G, State),
    Causes = causes,
    Closure = c(initial = initial_closure, gains = gain_closure, state = state_closure)
  )
}

# One stochastic substep. Initial individuals use a bounded multinomial fate
# draw. Gains generated during the interval are a Poisson process; Poisson
# thinning makes endpoint survivors and each loss-cause category independent
# Poisson counts with means supplied by the RK moment calculation.
INApestStochasticFluxStep <- function(State, Time, Step, FluxFunction,
                                      Parameters = NULL, Context = NULL, ...) {
  State <- .inapest_sfb_assert_whole_state(State)
  mom <- INApestStochasticFluxMoments(State, Time, Step, FluxFunction,
                                      Parameters = Parameters, Context = Context, ...)
  n <- length(State); m <- length(mom$Causes)
  init <- INApestStochasticBoundedProbabilities(State, mom$InitialLossProbabilities)
  new_survive <- stats::rpois(n, lambda = as.numeric(mom$SurvivingGainMean))
  new_loss <- if (m) {
    matrix(stats::rpois(n*m, lambda = as.numeric(mom$NewGainLossMeans)),
           nrow = n, ncol = m, dimnames = list(NULL, mom$Causes))
  } else matrix(numeric(), n, 0L)

  initial_loss <- init$Events
  total_loss <- if (m) initial_loss + new_loss else matrix(numeric(), n, 0L)
  total_gain <- new_survive + if (m) rowSums(new_loss) else 0
  final <- init$Survivors + new_survive

  # Exact realised bookkeeping identity for every state element.
  if (any(abs(final - (as.numeric(State) + total_gain - rowSums(total_loss))) > 0))
    stop("Internal stochastic-flux bookkeeping failure")
  final <- .inapest_sfb_restore(final, State)

  list(
    State = final,
    Gains = .inapest_sfb_restore(total_gain, State),
    Losses = total_loss,
    InitialLosses = initial_loss,
    NewGainLosses = new_loss,
    SurvivingGains = .inapest_sfb_restore(new_survive, State),
    Moments = mom
  )
}

# Integrate one arbitrary parent interval using stochastic internal steps no
# larger than MaxStep. Equal subdivision matches the continuous v0.1 kernel and
# makes results independent of parent boundaries when the same internal grid is
# induced over the same elapsed time.
INApestStochasticFluxIntegrate <- function(State, FluxFunction, Duration,
                                           MaxStep = 0.025, StartTime = 0,
                                           Parameters = NULL, Context = NULL, ...) {
  State <- .inapest_sfb_assert_whole_state(State)
  Duration <- as.numeric(Duration)[1L]
  MaxStep <- as.numeric(MaxStep)[1L]
  StartTime <- as.numeric(StartTime)[1L]
  if (!is.finite(Duration) || Duration < 0) stop("Duration must be finite and >= 0")
  if (!is.finite(MaxStep) || MaxStep <= 0) stop("MaxStep must be finite and > 0")
  if (!is.finite(StartTime)) stop("StartTime must be finite")
  if (Duration == 0) return(list(State = State, Gains = State*0,
                                 Losses = matrix(numeric(), length(State), 0L),
                                 NSubsteps = 0L, InternalStep = 0))

  nsub <- max(1L, as.integer(ceiling(Duration / MaxStep - 1e-14)))
  h <- Duration / nsub
  x <- State
  gains <- rep(0, length(State))
  losses <- NULL
  max_closure <- c(initial=0,gains=0,state=0)
  causes <- NULL
  tt <- StartTime
  for (ii in seq_len(nsub)) {
    z <- INApestStochasticFluxStep(x, Time = tt, Step = h,
                                   FluxFunction = FluxFunction,
                                   Parameters = Parameters, Context = Context, ...)
    if (is.null(causes)) {
      causes <- colnames(z$Losses)
      losses <- matrix(0, length(State), length(causes),
                       dimnames = list(NULL, causes))
    } else if (!identical(causes, colnames(z$Losses))) {
      stop("FluxFunction loss causes changed between stochastic substeps")
    }
    x <- z$State
    gains <- gains + as.numeric(z$Gains)
    if (length(causes)) losses <- losses + z$Losses
    max_closure <- pmax(max_closure, z$Moments$Closure)
    tt <- StartTime + ii*h
  }
  if (is.null(losses)) losses <- matrix(numeric(), length(State), 0L)
  list(
    State = x,
    Gains = .inapest_sfb_restore(gains, State),
    Losses = losses,
    NSubsteps = nsub,
    InternalStep = h,
    MaxStep = MaxStep,
    StartTime = StartTime,
    EndTime = StartTime + Duration,
    MaxClosureError = max_closure
  )
}

# LocalDynamics-compatible host-only stochastic adaptor. It returns the ordinary
# whole-count host state and attaches diagnostics; production engines may ignore
# the attribute without changing the numeric state contract.
INApestStochasticLocalDynamics <- function(
    HostFluxFunction,
    TimestepLength = 1,
    RKMaxStep = 0.025,
    Parameters = NULL,
    StartTime = 0) {
  if (!is.function(HostFluxFunction)) stop("HostFluxFunction must be a function")
  TimestepLength <- as.numeric(TimestepLength)[1L]
  RKMaxStep <- as.numeric(RKMaxStep)[1L]
  StartTime <- as.numeric(StartTime)[1L]
  if (!is.finite(TimestepLength) || TimestepLength <= 0)
    stop("TimestepLength must be finite and > 0")
  if (!is.finite(RKMaxStep) || RKMaxStep <= 0)
    stop("RKMaxStep must be finite and > 0")
  if (!is.finite(StartTime)) stop("StartTime must be finite")

  flux_fun <- HostFluxFunction; pars <- Parameters
  dt <- TimestepLength; hmax <- RKMaxStep; t0 <- StartTime
  f <- function(n0, timestep = NULL, Ntimesteps = NULL, ...) {
    step_index <- if (is.null(timestep)) 1L else as.integer(timestep)[1L]
    if (is.na(step_index) || step_index < 1L) stop("timestep must be >= 1 when supplied")
    this_start <- t0 + (step_index - 1L)*dt
    runtime <- list(...)
    reserved <- c("t","time","state","State","pars","Parameters","context","Context")
    collision <- intersect(names(runtime), reserved)
    if (length(collision))
      stop("LocalDynamics runtime arguments may not override reserved stochastic-RK argument(s): ",
           paste(collision, collapse=", "))
    bridge <- function(t, state, pars = NULL, ...) {
      .inapest_rk4_call_supported(
        flux_fun,
        c(list(t=t,time=t,state=state,State=state,pars=pars,Parameters=pars,
               timestep=step_index,Ntimesteps=Ntimesteps), runtime)
      )
    }
    z <- INApestStochasticFluxIntegrate(
      State=n0, FluxFunction=bridge, Duration=dt, MaxStep=hmax,
      StartTime=this_start, Parameters=pars)
    out <- z$State
    attr(out, "INApestStochasticFlux") <- z[setdiff(names(z), "State")]
    out
  }
  attr(f, "INApestStochasticLocalDynamics") <- list(
    Version="0.2", StateMode="integer", TimestepLength=TimestepLength,
    RKMaxStep=RKMaxStep, StochasticScheme="RK4 cohort-thinning flux bridge")
  class(f) <- c("INApestStochasticLocalDynamics","function")
  f
}

print.INApestStochasticLocalDynamics <- function(x, ...) {
  z <- attr(x,"INApestStochasticLocalDynamics")
  cat("INApest stochastic RK4 LocalDynamics adaptor\n")
  cat("  version:",z$Version,"\n")
  cat("  state mode:",z$StateMode,"\n")
  cat("  parent timestep:",z$TimestepLength,"\n")
  cat("  maximum stochastic/RK step:",z$RKMaxStep,"\n")
  invisible(x)
}

# Standard explicit stochastic host flux decompositions ------------------------
INApestHostLinearBirthDeathFlux <- function(BirthRate, DeathRate) {
  force(BirthRate); force(DeathRate)
  function(t, state, ...) {
    b <- .inapest_sfb_expand(BirthRate, state, "BirthRate")
    d <- .inapest_sfb_expand(DeathRate, state, "DeathRate")
    if (any(b < 0) || any(d < 0)) stop("BirthRate and DeathRate must be non-negative")
    list(GainRate = b*state, Hazards = list(death=d))
  }
}

INApestHostImmigrationDeathFlux <- function(ImmigrationRate, DeathRate) {
  force(ImmigrationRate); force(DeathRate)
  function(t, state, ...) {
    g <- .inapest_sfb_expand(ImmigrationRate, state, "ImmigrationRate")
    d <- .inapest_sfb_expand(DeathRate, state, "DeathRate")
    if (any(g < 0) || any(d < 0)) stop("ImmigrationRate and DeathRate must be non-negative")
    list(GainRate = g, Hazards = list(death=d))
  }
}

# Stochastic interpretation of the standard logistic ODE:
#   dN/dt = r*N*(1-N/K) - mu*N
# as gross births r*N, density-dependent mortality hazard r*N/K, and optional
# natural mortality hazard mu. This decomposition is explicit and therefore has
# a defined stochastic variance/extinction process.
INApestHostLogisticBirthDeathFlux <- function(GrowthRate, CarryingCapacity,
                                              MortalityRate = 0) {
  force(GrowthRate); force(CarryingCapacity); force(MortalityRate)
  function(t, state, ...) {
    r <- .inapest_sfb_expand(GrowthRate, state, "GrowthRate")
    K <- .inapest_sfb_expand(CarryingCapacity, state, "CarryingCapacity")
    mu <- .inapest_sfb_expand(MortalityRate, state, "MortalityRate")
    if (any(r < 0) || any(K <= 0) || any(mu < 0))
      stop("GrowthRate must be >=0, CarryingCapacity >0, and MortalityRate >=0")
    list(GainRate = r*state,
         Hazards = list(density = r*state/K, natural = mu))
  }
}
