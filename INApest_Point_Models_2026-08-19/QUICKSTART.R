###############################################################################
### INApest point-model quick start
### Run from the extracted package root.
###############################################################################

source(file.path("models", "INApestMetaPoint.R"))
source(file.path("models", "INApestPointTransitionMatrix.R"))

### 1. Generic INApestMeta-style point model ---------------------------------
meta_fit <- INApestMetaPoint(
  Nperm = 1,
  Ntimesteps = 10,
  InitialPoints = data.frame(x = 0, y = 0),
  Survival = 0.9,
  PropaguleProduction = 1.0,
  PropaguleEstablishment = 0.8,
  EnvEstabProb = 1,
  SDDkernel = INApestPointKernelExponential(mean_distance = 100),
  DetectionProb = 0,
  ManageProb = 0,
  MortalityProb = 0,
  FecundityReduction = 0,
  SpreadReduction = 0,
  DoProgress = FALSE,
  Seed = 1
)

### 2. Two-stage transition-matrix point model -------------------------------
A <- matrix(c(
  0.25, 1.30,
  0.60, 0.82
), nrow = 2, byrow = TRUE)

tm_fit <- INApestPointTransitionMatrix(
  Nperm = 1,
  Ntimesteps = 10,
  Nstages = 2,
  Transition = A,
  InitialPoints = data.frame(x = 0, y = 0, stage = 2L),
  SDDkernel = list(
    INApestPointKernelExponential(100),
    INApestPointKernelExponential(250)
  ),
  TransitionKernels = list(
    INApestPointKernelExponential(75)
  ),
  DetectionProb = c(0, 0),
  ManageProb = 0,
  MortalityProb = c(0, 0),
  FecundityReduction = 0,
  SpreadReduction = 0,
  DoProgress = FALSE,
  Seed = 1
)

print(tail(meta_fit$Summary))
print(tail(tm_fit$Summary))

### 3. Run the native-R animated demonstration -------------------------------
# source(file.path("demos", "INApestPointTransitionMatrix_unmanaged_gif_nativeR.R"))
