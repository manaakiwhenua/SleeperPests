###############################################################################
# INApestRK4.R
# Generic continuous-time RK4 numerical kernel for INApest local biology
# v0.1 -- 29 September 2026
#
# Scope
# -----
# This source contains numerical integration only. It does not know about
# hosts, pathogens, biocontrol, movement, surveillance or management.
#
# The production INApest engines are stochastic integer-state simulators.
# Therefore this continuous kernel MUST NOT be inserted into an integer parent
# by silently rounding its output. A separately validated stochastic flux
# bridge is required for that use. The existing biocontrol stochastic
# fixed-delay/RK4 implementation is the parent example for that later bridge.
###############################################################################

.inapest_rk4_call_supported <- function(fun, args) {
  fm <- names(formals(fun))
  if (!is.null(fm) && !("..." %in% fm)) args <- args[intersect(names(args), fm)]
  do.call(fun, args)
}

.inapest_rk4_restore <- function(x, template) {
  out <- as.numeric(x)
  dim(out) <- dim(template)
  dimnames(out) <- dimnames(template)
  if (is.null(dim(template))) names(out) <- names(template)
  out
}

.inapest_rk4_validate_derivative <- function(dx, state) {
  if (!is.numeric(dx)) stop("RateFunction must return a numeric derivative")
  if (length(dx) != length(state))
    stop("RateFunction derivative must have the same length as state")
  if (any(!is.finite(dx))) stop("RateFunction returned non-finite derivative(s)")
  .inapest_rk4_restore(dx, state)
}

# One classical fourth-order Runge-Kutta step.
INApestRK4Step <- function(State, Time, Step, RateFunction,
                          Parameters = NULL, Context = NULL,
                          NonNegative = c("error", "allow"), ...) {
  NonNegative <- match.arg(NonNegative)
  if (!is.numeric(State) || any(!is.finite(State)))
    stop("State must contain finite numeric values")
  Step <- as.numeric(Step)[1L]
  Time <- as.numeric(Time)[1L]
  if (!is.finite(Step) || Step <= 0) stop("Step must be finite and > 0")
  if (!is.finite(Time)) stop("Time must be finite")
  if (!is.function(RateFunction)) stop("RateFunction must be a function")
  extra <- list(...)

  deriv <- function(tt, xx) {
    xx <- .inapest_rk4_restore(xx, State)
    ans <- .inapest_rk4_call_supported(
      RateFunction,
      c(list(t = tt, time = tt, state = xx, State = xx,
             pars = Parameters, Parameters = Parameters,
             context = Context, Context = Context), extra)
    )
    .inapest_rk4_validate_derivative(ans, State)
  }

  y <- as.numeric(State)
  k1 <- as.numeric(deriv(Time, y))
  k2 <- as.numeric(deriv(Time + Step / 2, y + Step * k1 / 2))
  k3 <- as.numeric(deriv(Time + Step / 2, y + Step * k2 / 2))
  k4 <- as.numeric(deriv(Time + Step,     y + Step * k3))
  out <- .inapest_rk4_restore(
    y + Step * (k1 + 2 * k2 + 2 * k3 + k4) / 6,
    State
  )

  if (NonNegative == "error" && any(out < -1e-12))
    stop("RK4 step produced a negative biological state; reduce RKMaxStep or revise the rate model")
  # Numerical round-off near zero is normal. Values are not generally clipped;
  # only tiny negative floating-point noise is mapped to exact zero.
  if (NonNegative == "error") out[out < 0 & out >= -1e-12] <- 0
  out
}

# Integrate over one arbitrary interval using as many equal RK4 substeps as are
# required to keep the internal step <= MaxStep. The final substep is therefore
# never larger than MaxStep and the interval endpoint is reached exactly.
INApestRK4Integrate <- function(State, RateFunction, Duration,
                                MaxStep = 0.025, StartTime = 0,
                                Parameters = NULL, Context = NULL,
                                NonNegative = c("error", "allow"), ...) {
  NonNegative <- match.arg(NonNegative)
  Duration <- as.numeric(Duration)[1L]
  MaxStep <- as.numeric(MaxStep)[1L]
  StartTime <- as.numeric(StartTime)[1L]
  if (!is.finite(Duration) || Duration < 0) stop("Duration must be finite and >= 0")
  if (!is.finite(MaxStep) || MaxStep <= 0) stop("MaxStep must be finite and > 0")
  if (!is.finite(StartTime)) stop("StartTime must be finite")
  if (Duration == 0) return(State)

  nsub <- max(1L, as.integer(ceiling(Duration / MaxStep - 1e-14)))
  h <- Duration / nsub
  x <- State
  tt <- StartTime
  for (ii in seq_len(nsub)) {
    x <- INApestRK4Step(
      State = x, Time = tt, Step = h, RateFunction = RateFunction,
      Parameters = Parameters, Context = Context,
      NonNegative = NonNegative, ...
    )
    tt <- StartTime + ii * h
  }
  attr(x, "INApestRK4") <- list(
    StartTime = StartTime,
    EndTime = StartTime + Duration,
    Duration = Duration,
    NSubsteps = nsub,
    InternalStep = h,
    MaxStep = MaxStep
  )
  x
}

# Compatibility primitive extracted from the accepted biocontrol parent.
# It deliberately preserves the parent's arithmetic exactly so regression can
# establish that genericisation did not change this validated RK calculation.
INApestRK4DecayIntegral <- function(x, rate, dt) {
  x <- as.numeric(x)
  rate <- as.numeric(rate)[1L]
  dt <- as.numeric(dt)[1L]
  if (!is.finite(rate) || rate < 0 || !is.finite(dt) || dt <= 0)
    stop("rate must be >=0 and dt must be >0")

  k1x <- -rate * x; k1i <- x
  x2 <- x + dt * k1x / 2
  k2x <- -rate * x2; k2i <- x2
  x3 <- x + dt * k2x / 2
  k3x <- -rate * x3; k3i <- x3
  x4 <- x + dt * k3x
  k4x <- -rate * x4; k4i <- x4

  list(
    End = x + dt * (k1x + 2*k2x + 2*k3x + k4x) / 6,
    Integral = dt * (k1i + 2*k2i + 2*k3i + k4i) / 6)
}
