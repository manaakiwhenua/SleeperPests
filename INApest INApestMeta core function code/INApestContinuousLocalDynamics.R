###############################################################################
# INApestContinuousLocalDynamics.R
# General continuous host LocalDynamics adaptor and standard host rate models
# v0.1 -- 29 September 2026
###############################################################################

if (!exists("INApestRK4Integrate", mode = "function"))
  stop("Source INApestRK4.R before INApestContinuousLocalDynamics.R")

.inapest_host_expand <- function(x, state, label) {
  if (is.function(x)) return(x)
  z <- as.numeric(x)
  if (!length(z) || any(!is.finite(z))) stop(label, " must be finite")
  if (length(z) == 1L) {
    out <- rep(z, length(state))
    dim(out) <- dim(state); dimnames(out) <- dimnames(state)
    if (is.null(dim(state))) names(out) <- names(state)
    return(out)
  }
  if (length(z) != length(state))
    stop(label, " must be scalar or have one value per host state")
  out <- z
  dim(out) <- dim(state); dimnames(out) <- dimnames(state)
  if (is.null(dim(state))) names(out) <- names(state)
  out
}

# Standard host rate constructors ------------------------------------------------
INApestHostExponentialRates <- function(Rate) {
  force(Rate)
  function(t, state, ...) {
    r <- .inapest_host_expand(Rate, state, "Rate")
    r * state
  }
}

INApestHostLogisticRates <- function(GrowthRate, CarryingCapacity,
                                     MortalityRate = 0) {
  force(GrowthRate); force(CarryingCapacity); force(MortalityRate)
  function(t, state, ...) {
    r <- .inapest_host_expand(GrowthRate, state, "GrowthRate")
    K <- .inapest_host_expand(CarryingCapacity, state, "CarryingCapacity")
    mu <- .inapest_host_expand(MortalityRate, state, "MortalityRate")
    if (any(K <= 0)) stop("CarryingCapacity must be > 0")
    r * state * (1 - state / K) - mu * state
  }
}

INApestHostAlleeRates <- function(GrowthRate, CarryingCapacity,
                                  AlleeThreshold, MortalityRate = 0) {
  force(GrowthRate); force(CarryingCapacity); force(AlleeThreshold); force(MortalityRate)
  function(t, state, ...) {
    r <- .inapest_host_expand(GrowthRate, state, "GrowthRate")
    K <- .inapest_host_expand(CarryingCapacity, state, "CarryingCapacity")
    A <- .inapest_host_expand(AlleeThreshold, state, "AlleeThreshold")
    mu <- .inapest_host_expand(MortalityRate, state, "MortalityRate")
    if (any(K <= 0)) stop("CarryingCapacity must be > 0")
    if (any(A <= 0 | A >= K)) stop("AlleeThreshold must be > 0 and < CarryingCapacity")
    r * state * (1 - state / K) * (state / A - 1) - mu * state
  }
}

# Generic adaptor ----------------------------------------------------------------
# Returns an ordinary function suitable for the existing LocalDynamics contract
# when continuous host state is scientifically appropriate. It accepts n0 and
# ... so the parent engine can pass its normal LocalDynamics arguments without
# changing the public call. Only arguments explicitly declared by the supplied
# HostRateFunction are forwarded to that rate function.
#
# IMPORTANT: current stochastic integer INApest engines must not consume the
# fractional result by rounding. Production stochastic use requires a validated
# integer stochastic-flux adaptor. This constructor advertises its state mode
# through attributes so callers/tests can enforce that boundary.
INApestContinuousLocalDynamics <- function(
    HostRateFunction,
    TimestepLength = 1,
    RKMaxStep = 0.025,
    Parameters = NULL,
    NonNegative = c("error", "allow"),
    StartTime = 0) {
  NonNegative <- match.arg(NonNegative)
  if (!is.function(HostRateFunction)) stop("HostRateFunction must be a function")
  TimestepLength <- as.numeric(TimestepLength)[1L]
  RKMaxStep <- as.numeric(RKMaxStep)[1L]
  StartTime <- as.numeric(StartTime)[1L]
  if (!is.finite(TimestepLength) || TimestepLength <= 0)
    stop("TimestepLength must be finite and > 0")
  if (!is.finite(RKMaxStep) || RKMaxStep <= 0)
    stop("RKMaxStep must be finite and > 0")
  if (!is.finite(StartTime)) stop("StartTime must be finite")

  rate_fun <- HostRateFunction
  pars <- Parameters
  dt <- TimestepLength
  hmax <- RKMaxStep
  nonneg <- NonNegative
  t0 <- StartTime

  f <- function(n0, timestep = NULL, Ntimesteps = NULL, ...) {
    step_index <- if (is.null(timestep)) 1L else as.integer(timestep)[1L]
    if (is.na(step_index) || step_index < 1L) stop("timestep must be >= 1 when supplied")
    this_start <- t0 + (step_index - 1L) * dt
    runtime <- list(...)
    reserved <- c("t", "time", "state", "State", "pars", "Parameters")
    collision <- intersect(names(runtime), reserved)
    if (length(collision))
      stop("LocalDynamics runtime arguments may not override reserved RK argument(s): ",
           paste(collision, collapse = ", "))
    bridge <- function(t, state, pars = NULL, ...) {
      .inapest_rk4_call_supported(
        rate_fun,
        c(list(t = t, time = t, state = state, State = state,
               pars = pars, Parameters = pars,
               timestep = step_index, Ntimesteps = Ntimesteps), runtime)
      )
    }
    INApestRK4Integrate(
      State = n0,
      RateFunction = bridge,
      Duration = dt,
      MaxStep = hmax,
      StartTime = this_start,
      Parameters = pars,
      NonNegative = nonneg
    )
  }
  attr(f, "INApestContinuousLocalDynamics") <- list(
    Version = "0.1",
    StateMode = "continuous",
    TimestepLength = TimestepLength,
    RKMaxStep = RKMaxStep,
    NonNegative = NonNegative
  )
  class(f) <- c("INApestContinuousLocalDynamics", "function")
  f
}

print.INApestContinuousLocalDynamics <- function(x, ...) {
  z <- attr(x, "INApestContinuousLocalDynamics")
  cat("INApest continuous LocalDynamics adaptor\n")
  cat("  version:", z$Version, "\n")
  cat("  state mode:", z$StateMode, "\n")
  cat("  parent timestep:", z$TimestepLength, "\n")
  cat("  maximum RK4 step:", z$RKMaxStep, "\n")
  invisible(x)
}
