cat('=== Base-R habitat grid tests ===\n')

source('/mnt-data/INApestMetaPoint.R')

# -----------------------------------------------------------------------------
# 1. Static grid construction, orientation, boundaries, outside extent, NA
# -----------------------------------------------------------------------------
g <- INApestHabitatGrid(
  values = matrix(c(
    0.10, 0.20,
    0.30, 0.90
  ), nrow = 2, byrow = TRUE),
  xmin = 0, xmax = 2,
  ymin = 0, ymax = 2
)

v <- .ipp_habitat_value(
  x = c(0.5, 1.5, 0.5, 1.5, 2.0, -0.1),
  y = c(1.5, 1.5, 0.5, 0.5, 0.0, 1.0),
  HabitatSuitability = g,
  timestep = 1
)
stopifnot(isTRUE(all.equal(v, c(0.10, 0.20, 0.30, 0.90, 0.90, 0))))

g_na <- INApestHabitatGrid(
  values = matrix(c(NA, 0.2, 0.3, 0.9), nrow = 2, byrow = TRUE),
  xmin = 0, xmax = 2, ymin = 0, ymax = 2
)
stopifnot(.ipp_habitat_value(0.5, 1.5, g_na, 1) == 0)

# Explicit reverse row orientation should also work.
g_south <- INApestHabitatGrid(
  values = matrix(c(
    0.30, 0.90,
    0.10, 0.20
  ), nrow = 2, byrow = TRUE),
  xmin = 0, xmax = 2, ymin = 0, ymax = 2,
  row_origin = 'ymin'
)
stopifnot(.ipp_habitat_value(1.5, 0.5, g_south, 1) == 0.90)
stopifnot(.ipp_habitat_value(0.5, 1.5, g_south, 1) == 0.10)

# Invalid suitability should be rejected rather than silently clipped.
err <- try(INApestHabitatGrid(matrix(c(0, 1.1), nrow = 1), 0, 2, 0, 1), silent = TRUE)
stopifnot(inherits(err, 'try-error'))

# -----------------------------------------------------------------------------
# 2. Time-varying grid layers
# -----------------------------------------------------------------------------
a <- array(NA_real_, dim = c(2, 2, 2))
a[, , 1] <- 0.2
a[, , 2] <- 0.8
g_time <- INApestHabitatGrid(a, 0, 2, 0, 2)
stopifnot(.ipp_habitat_value(0.5, 0.5, g_time, 1) == 0.2)
stopifnot(.ipp_habitat_value(0.5, 0.5, g_time, 2) == 0.8)

# -----------------------------------------------------------------------------
# 3. Exact base-R habitat nudge: highest better habitat, nearest tie, radius cap
# -----------------------------------------------------------------------------
g_search <- INApestHabitatGrid(
  values = matrix(c(
    0.1, 0.1,
    0.2, 0.9
  ), nrow = 2, byrow = TRUE),
  xmin = 0, xmax = 2, ymin = 0, ymax = 2
)

s1 <- .ipp_habitat_search(
  x = 0.5, y = 0.5,
  HabitatSuitability = g_search,
  HabitatSearchRadius = 1.01,
  timestep = 1
)
stopifnot(s1$habitat_nudged)
stopifnot(abs(s1$x - 1.5) < 1e-12)
stopifnot(abs(s1$y - 0.5) < 1e-12)
stopifnot(abs(s1$habitat - 0.9) < 1e-12)
stopifnot(sqrt((s1$x - 0.5)^2 + (s1$y - 0.5)^2) <= 1.01)

s2 <- .ipp_habitat_search(
  x = 0.5, y = 0.5,
  HabitatSuitability = g_search,
  HabitatSearchRadius = 0.9,
  timestep = 1
)
stopifnot(!s2$habitat_nudged)
stopifnot(s2$x == 0.5, s2$y == 0.5, s2$habitat == 0.2)

# -----------------------------------------------------------------------------
# 4. INApestPoint: habitat multiplier in establishment
# Expected offspring: 2000 parents * Poisson mean 2 * 0.8 * 0.5 * 0.25 = 400
# -----------------------------------------------------------------------------
g_quarter <- INApestHabitatGrid(matrix(0.25, nrow = 1, ncol = 1), 0, 1, 0, 1)
initial <- data.frame(x = rep(0.5, 2000), y = rep(0.5, 2000))

fit <- INApestMetaPoint(
  Nperm = 1,
  Ntimesteps = 1,
  InitialPoints = initial,
  Survival = 1,
  PropaguleProduction = 2,
  PropaguleEstablishment = 0.8,
  EnvEstabProb = 0.5,
  SDDkernel = INApestPointKernelFixed(0),
  HabitatSuitability = g_quarter,
  HabitatSearchRadius = 0,
  DetectionProb = 0, DetectionSD = 0,
  ManageProb = 0, ManageSD = 0,
  MortalityProb = 0, MortalitySD = 0,
  FecundityReduction = 0, FecundityReductionSD = 0,
  SpreadReduction = 0, SpreadReductionSD = 0,
  Seed = 123,
  DoProgress = FALSE
)

est <- fit$Summary$n_established[1]
expected <- 400
stopifnot(abs(est - expected) <= 5 * sqrt(expected))

# -----------------------------------------------------------------------------
# 5. INApestPoint: model-level habitat nudge
# -----------------------------------------------------------------------------
fit_nudge <- INApestMetaPoint(
  Nperm = 1, Ntimesteps = 1,
  InitialPoints = data.frame(x = rep(0.5, 50), y = rep(0.5, 50)),
  Survival = 1,
  PropaguleProduction = 2,
  PropaguleEstablishment = 1,
  EnvEstabProb = 1,
  SDDkernel = INApestPointKernelFixed(0),
  HabitatSuitability = g_search,
  HabitatSearchRadius = 1.01,
  DetectionProb = 0, DetectionSD = 0,
  ManageProb = 0, ManageSD = 0,
  MortalityProb = 0, MortalitySD = 0,
  FecundityReduction = 0, FecundityReductionSD = 0,
  SpreadReduction = 0, SpreadReductionSD = 0,
  Seed = 456,
  DoProgress = FALSE
)

ev <- fit_nudge$EventLog[fit_nudge$EventLog$event == 'establishment', , drop = FALSE]
stopifnot(nrow(ev) > 0)
stopifnot(all(ev$habitat_nudged))
stopifnot(all(abs(ev$x - 1.5) < 1e-12))
stopifnot(all(abs(ev$y - 0.5) < 1e-12))
stopifnot(all(abs(ev$habitat - 0.9) < 1e-12))

# -----------------------------------------------------------------------------
# 6. Early rejection of insufficient time layers
# -----------------------------------------------------------------------------
err <- try(INApestMetaPoint(
  Nperm = 1, Ntimesteps = 3,
  InitialPoints = data.frame(x = 0.5, y = 0.5),
  Survival = 1, PropaguleProduction = 0, PropaguleEstablishment = 1,
  EnvEstabProb = 1, SDDkernel = INApestPointKernelFixed(0),
  HabitatSuitability = g_time,
  DoProgress = FALSE
), silent = TRUE)
stopifnot(inherits(err, 'try-error'))

cat('INApestPoint base-R habitat tests passed. Established=', est, '\n', sep='')

# -----------------------------------------------------------------------------
# 7. Transition-matrix engine with same native grid
# -----------------------------------------------------------------------------
source('/mnt-data/INApestPointTransitionMatrix.R')

# Reproductive establishment test from stage 2.
# Stage 2 survives and produces 2 stage-1 propagules on average.
A_repro <- matrix(c(
  1, 2,
  0, 1
), nrow = 2, byrow = TRUE)

fit_tm <- INApestPointTransitionMatrix(
  Nperm = 1, Ntimesteps = 1, Nstages = 2,
  Transition = A_repro,
  InitialPoints = data.frame(x = rep(0.5, 2000), y = rep(0.5, 2000), stage = 2L),
  SDDkernel = INApestPointKernelFixed(0),
  PropaguleEstablishment = 0.8,
  EnvEstabProb = 0.5,
  HabitatSuitability = g_quarter,
  HabitatSearchRadius = 0,
  DetectionProb = 0, DetectionSD = 0,
  ManageProb = 0, ManageSD = 0,
  MortalityProb = 0, MortalitySD = 0,
  FecundityReduction = 0, FecundityReductionSD = 0,
  SpreadReduction = 0, SpreadReductionSD = 0,
  Seed = 123,
  DoProgress = FALSE
)

est_tm <- fit_tm$Summary$n_established[1]
stopifnot(abs(est_tm - expected) <= 5 * sqrt(expected))

# Transition movement itself can use the same habitat-search grid.
A_prog <- matrix(c(
  0, 0,
  1, 1
), nrow = 2, byrow = TRUE)

fit_tm_nudge <- INApestPointTransitionMatrix(
  Nperm = 1, Ntimesteps = 1, Nstages = 2,
  Transition = A_prog,
  InitialPoints = data.frame(x = rep(0.5, 50), y = rep(0.5, 50), stage = 1L),
  SDDkernel = INApestPointKernelFixed(0),
  TransitionKernels = list(INApestPointKernelFixed(0)),
  TransitionHabitatSearch = TRUE,
  ApplyHabitatToTransitions = FALSE,
  HabitatSuitability = g_search,
  HabitatSearchRadius = 1.01,
  DetectionProb = 0, DetectionSD = 0,
  ManageProb = 0, ManageSD = 0,
  MortalityProb = 0, MortalitySD = 0,
  FecundityReduction = 0, FecundityReductionSD = 0,
  SpreadReduction = 0, SpreadReductionSD = 0,
  Seed = 789,
  DoProgress = FALSE
)

ev_tm <- fit_tm_nudge$EventLog[fit_tm_nudge$EventLog$event == 'stage_transition', , drop = FALSE]
stopifnot(nrow(ev_tm) == 50)
stopifnot(all(ev_tm$habitat_nudged))
stopifnot(all(abs(ev_tm$x - 1.5) < 1e-12))
stopifnot(all(abs(ev_tm$y - 0.5) < 1e-12))

# Habitat can also block a stage transition when explicitly requested.
g_zero <- INApestHabitatGrid(matrix(0, 1, 1), 0, 1, 0, 1)
fit_tm_block <- INApestPointTransitionMatrix(
  Nperm = 1, Ntimesteps = 1, Nstages = 2,
  Transition = A_prog,
  InitialPoints = data.frame(x = rep(0.5, 20), y = rep(0.5, 20), stage = 1L),
  SDDkernel = INApestPointKernelFixed(0),
  TransitionKernels = list(NULL),
  TransitionHabitatSearch = FALSE,
  ApplyHabitatToTransitions = TRUE,
  TransitionEstablishment = 1,
  BlockedTransitionMortality = 1,
  EnvEstabProb = 1,
  HabitatSuitability = g_zero,
  HabitatSearchRadius = 0,
  DetectionProb = 0, DetectionSD = 0,
  ManageProb = 0, ManageSD = 0,
  MortalityProb = 0, MortalitySD = 0,
  FecundityReduction = 0, FecundityReductionSD = 0,
  SpreadReduction = 0, SpreadReductionSD = 0,
  Seed = 790,
  DoProgress = FALSE
)
stopifnot(fit_tm_block$Summary$n_stage_transitions[1] == 0)
stopifnot(fit_tm_block$Summary$n_end[1] == 0)

cat('INApestPointTransitionMatrix base-R habitat tests passed. Established=', est_tm, '\n', sep='')
cat('ALL BASE-R HABITAT GRID TESTS PASSED\n')
