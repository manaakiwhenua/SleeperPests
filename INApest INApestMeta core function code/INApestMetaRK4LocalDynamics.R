###############################################################################
# INApestMetaRK4LocalDynamics.R
# Production Meta adaptor for stochastic continuous-time host LocalDynamics
# v0.3 -- 30 September 2026
#
# This adaptor composes the validated generic stochastic RK4 flux bridge with
# the established INApestMeta dispersal/establishment operator. It does not
# replace the default LocalDynamics path.
#
# Event order within the existing Meta biological boundary is:
#   1. parent engine survival / management mortality (unchanged);
#   2. stochastic RK4 local host demographic fluxes;
#   3. established Meta propagule production, dispersal and establishment;
#   4. optional pathogen reconciliation / biocontrol / surveillance downstream
#      in the parent engine (unchanged).
#
# HostFluxFunction therefore describes changes in RESIDENT host abundance.
# PropaguleProduction remains the established dispersal-output process. Users
# are responsible for choosing a biological parameterisation in which those two
# processes have the intended interpretation.
###############################################################################

if (!exists("INApestStochasticLocalDynamics", mode = "function"))
  stop("Source INApestRK4.R and INApestStochasticFluxBridge.R before INApestMetaRK4LocalDynamics.R")

INApestMetaRK4LocalDynamics <- function(
    HostFluxFunction,
    TimestepLength = 1,
    RKMaxStep = 0.025,
    Parameters = NULL,
    StartTime = 0,
    DispersalDynamics = NULL,
    RunWhenEmpty = FALSE,
    CapacityTolerance = 1e-9) {

  if (!is.function(HostFluxFunction)) stop("HostFluxFunction must be a function")
  if (is.null(DispersalDynamics)) {
    if (!exists("local.dynamics", mode = "function"))
      stop("Source the current INApestMeta.r before constructing the Meta RK4 LocalDynamics adaptor")
    DispersalDynamics <- get("local.dynamics", mode = "function")
  }
  if (!is.function(DispersalDynamics)) stop("DispersalDynamics must be a function")
  if (!is.logical(RunWhenEmpty) || length(RunWhenEmpty) != 1L || is.na(RunWhenEmpty))
    stop("RunWhenEmpty must be TRUE or FALSE")
  CapacityTolerance <- as.numeric(CapacityTolerance)[1L]
  if (!is.finite(CapacityTolerance) || CapacityTolerance < 0)
    stop("CapacityTolerance must be finite and >= 0")

  stochastic_local <- INApestStochasticLocalDynamics(
    HostFluxFunction = HostFluxFunction,
    TimestepLength = TimestepLength,
    RKMaxStep = RKMaxStep,
    Parameters = Parameters,
    StartTime = StartTime
  )
  dispersal_fun <- DispersalDynamics
  cap_tol <- CapacityTolerance

  # Match the established Meta LocalDynamics core contract. timestep and
  # Ntimesteps are optional read-only additions passed only when explicitly
  # requested by the custom LocalDynamics function.
  f <- function(
      sddprob,
      nodepropaguleproduction,
      nodeenvestabprob,
      n,
      lddprob,
      lddrate,
      k_is_0,
      nodeK,
      nodepropaguleestablishment,
      nodespreadreduction,
      nodefecundityreduction = 0,
      managing,
      maxinteger,
      timestep = NULL,
      Ntimesteps = NULL) {

    n_before_rk <- n
    n_after_rk <- stochastic_local(
      n0 = n,
      timestep = timestep,
      Ntimesteps = Ntimesteps,
      nodeK = nodeK,
      k_is_0 = k_is_0,
      nodeenvestabprob = nodeenvestabprob,
      nodepropaguleproduction = nodepropaguleproduction,
      nodepropaguleestablishment = nodepropaguleestablishment,
      nodespreadreduction = nodespreadreduction,
      nodefecundityreduction = nodefecundityreduction,
      managing = managing,
      sddprob = sddprob,
      lddprob = lddprob,
      lddrate = lddrate
    )

    # Meta uses nodeK as a hard establishment capacity. The generic stochastic
    # bridge intentionally does not impose a capacity because that is a biology
    # choice. Do not silently clip an RK realisation that is incompatible with
    # the parent Meta capacity contract.
    Kvec <- as.numeric(nodeK)
    Nvec <- as.numeric(n_after_rk)
    if (length(Kvec) == 1L) Kvec <- rep(Kvec, length(Nvec))
    if (length(Kvec) != length(Nvec) || any(!is.finite(Kvec)) || any(Kvec < 0))
      stop("nodeK must resolve to finite non-negative capacity for every Meta node")
    over <- which(Nvec > Kvec + cap_tol * pmax(1, abs(Kvec)))
    if (length(over))
      stop("RK host state exceeded Meta nodeK before dispersal at node(s): ",
           paste(over, collapse = ", "),
           ". Use a capacity-compatible HostFluxFunction or revise the model; the adaptor does not silently clip stochastic host counts.")

    # Preserve the established Meta dispersal/establishment implementation.
    out <- dispersal_fun(
      sddprob = sddprob,
      nodepropaguleproduction = nodepropaguleproduction,
      nodeenvestabprob = nodeenvestabprob,
      n = n_after_rk,
      lddprob = lddprob,
      lddrate = lddrate,
      k_is_0 = k_is_0,
      nodeK = nodeK,
      nodepropaguleestablishment = nodepropaguleestablishment,
      nodespreadreduction = nodespreadreduction,
      nodefecundityreduction = nodefecundityreduction,
      managing = managing,
      maxinteger = maxinteger
    )

    # Diagnostics are available for direct LocalDynamics tests. The parent Meta
    # result arrays intentionally retain their established numeric contract.
    attr(out, "INApestMetaRK4") <- list(
      StateBeforeRK = n_before_rk,
      StateAfterRK = n_after_rk,
      Flux = attr(n_after_rk, "INApestStochasticFlux"),
      Timestep = timestep,
      Ntimesteps = Ntimesteps
    )
    out
  }

  meta <- list(
    Version = "0.3",
    StateMode = "integer",
    TimestepLength = as.numeric(TimestepLength)[1L],
    RKMaxStep = as.numeric(RKMaxStep)[1L],
    RunWhenEmpty = RunWhenEmpty,
    CapacityPolicy = "error-no-clipping",
    Ordering = c("parent survival/management mortality",
                 "RK stochastic local host dynamics",
                 "Meta propagule production/dispersal/establishment")
  )
  attr(f, "INApestMetaRK4LocalDynamics") <- meta
  attr(f, "INApestRunWhenEmpty") <- RunWhenEmpty
  class(f) <- c("INApestMetaRK4LocalDynamics", "function")
  f
}

print.INApestMetaRK4LocalDynamics <- function(x, ...) {
  z <- attr(x, "INApestMetaRK4LocalDynamics")
  cat("INApest Meta stochastic RK4 LocalDynamics adaptor\n")
  cat("  version:", z$Version, "\n")
  cat("  state mode:", z$StateMode, "\n")
  cat("  parent timestep:", z$TimestepLength, "\n")
  cat("  maximum stochastic/RK step:", z$RKMaxStep, "\n")
  cat("  run when globally empty:", z$RunWhenEmpty, "\n")
  invisible(x)
}
