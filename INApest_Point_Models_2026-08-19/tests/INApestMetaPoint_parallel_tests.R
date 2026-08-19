###############################################################################
### Parallel wrapper regression tests
###
### In WebR, true PSOCK multiprocessing is skipped because the browser/WASM
### runtime cannot spawn ordinary R worker processes. In native R the same
### script additionally checks 2-worker PSOCK equality against the 1-core
### wrapper result for a fixed parallel RNG stream schedule.
###############################################################################

source("/mnt-data/INApestMetaPoint.R")
source("/mnt-data/INApestPointTransitionMatrix.R")
source("/mnt-data/INApestMetaPointParallel.R")
source("/mnt-data/INApestPointTransitionMatrixParallel.R")

point_args <- list(
  Ntimesteps = 3,
  InitialPoints = data.frame(x = 0, y = 0),
  Survival = 1,
  PropaguleProduction = 0.5,
  PropaguleEstablishment = 1,
  EnvEstabProb = 1,
  SDDkernel = INApestPointKernelFixed(10),
  DetectionProb = 0, DetectionSD = 0,
  ManageProb = 0, ManageSD = 0,
  MortalityProb = 0, MortalitySD = 0,
  FecundityReduction = 0, FecundityReductionSD = 0,
  SpreadReduction = 0, SpreadReductionSD = 0,
  DoProgress = FALSE
)

p1 <- do.call(INApestMetaPointParallel, c(list(Nperm = 4, Cores = 1, Seed = 123), point_args))
p2 <- do.call(INApestMetaPointParallel, c(list(Nperm = 4, Cores = 1, Seed = 123), point_args))
stopifnot(identical(p1$Summary, p2$Summary))
stopifnot(identical(p1$EventLog, p2$EventLog))
stopifnot(identical(sort(unique(p1$Summary$perm)), 1:4))

A <- matrix(c(0.1, 1.2, 0.8, 0.05), 2, 2, byrow = TRUE)
tm_args <- list(
  Ntimesteps = 3, Nstages = 2, Transition = A,
  InitialPoints = data.frame(x = 0, y = 0, stage = 2),
  SDDkernel = INApestPointKernelFixed(20),
  TransitionKernels = list(INApestPointKernelFixed(2)),
  DetectionProb = c(0, 0), DetectionSD = 0,
  ManageProb = 0, ManageSD = 0,
  MortalityProb = c(0, 0), MortalitySD = 0,
  FecundityReduction = 0, FecundityReductionSD = 0,
  SpreadReduction = 0, SpreadReductionSD = 0,
  DoProgress = FALSE
)

t1 <- do.call(INApestPointTransitionMatrixParallel, c(list(Nperm = 4, Cores = 1, Seed = 321), tm_args))
t2 <- do.call(INApestPointTransitionMatrixParallel, c(list(Nperm = 4, Cores = 1, Seed = 321), tm_args))
stopifnot(identical(t1$Summary, t2$Summary))
stopifnot(identical(t1$EventLog, t2$EventLog))
stopifnot(identical(sort(unique(t1$Summary$perm)), 1:4))

# Native-R cross-backend reproducibility check.
if (!grepl("wasm", R.version$platform, ignore.case = TRUE) && parallel::detectCores() >= 2) {
  pp <- do.call(INApestMetaPointParallel, c(list(Nperm = 4, Cores = 2, Backend = "psock", Seed = 123), point_args))
  stopifnot(identical(p1$Summary, pp$Summary))
  stopifnot(identical(p1$EventLog, pp$EventLog))
  stopifnot(identical(p1$FinalPoints, pp$FinalPoints))

  tp <- do.call(INApestPointTransitionMatrixParallel, c(list(Nperm = 4, Cores = 2, Backend = "psock", Seed = 321), tm_args))
  stopifnot(identical(t1$Summary, tp$Summary))
  stopifnot(identical(t1$EventLog, tp$EventLog))
  stopifnot(identical(t1$FinalPoints, tp$FinalPoints))
}

cat("Parallel wrapper regression tests passed for available backend(s).\n")
