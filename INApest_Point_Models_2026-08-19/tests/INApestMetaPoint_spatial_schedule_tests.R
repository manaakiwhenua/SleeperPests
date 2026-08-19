###############################################################################
### Regression tests for human-friendly INApestSpatialSchedule
###############################################################################

source('/mnt-data/INApestMetaPoint.R')

k0 <- INApestPointKernelFixed(0)

g0 <- INApestSpatialGrid(matrix(0, 2, 2), xmin=0, xmax=2, ymin=0, ymax=2)
g1 <- INApestSpatialGrid(matrix(1, 2, 2), xmin=0, xmax=2, ymin=0, ymax=2)
g5 <- INApestSpatialGrid(matrix(0.5, 2, 2), xmin=0, xmax=2, ymin=0, ymax=2)

# 1. Recommended change-point syntax: named grids + start timesteps.
s <- INApestSpatialSchedule(
  early = g0,
  response = g5,
  intensive = g1,
  from = c(1, 3, 6)
)
stopifnot(inherits(s, 'INApestSpatialSchedule'))
stopifnot(identical(s$from, c(1L, 3L, 6L)))
stopifnot(identical(s$to, c(2, 5, Inf)))
stopifnot(identical(names(s$grids), c('early','response','intensive')))

v <- vapply(1:8, function(tt) .ipp_spatial_value(1, 1, s, tt), numeric(1))
stopifnot(identical(v, c(0, 0, 0.5, 0.5, 0.5, 1, 1, 1)))

# 2. Explicit ranges remain available when desired.
s_explicit <- INApestSpatialSchedule(
  early = g0,
  later = g1,
  from = c(1, 4),
  to = c(3, 7)
)
stopifnot(.ipp_spatial_value(1,1,s_explicit,3) == 0)
stopifnot(.ipp_spatial_value(1,1,s_explicit,4) == 1)
.ipp_validate_spatial_surface(s_explicit, 7, 'test')

# 3. A list is accepted for programmatic setup.
s_list <- INApestSpatialSchedule(
  list(off = g0, on = g1),
  from = c(1, 2)
)
stopifnot(.ipp_spatial_value(1,1,s_list,1) == 0)
stopifnot(.ipp_spatial_value(1,1,s_list,2) == 1)

# 4. Gaps fail early and explain the missing timestep.
gap <- INApestSpatialSchedule(
  first = g0,
  second = g1,
  from = c(1, 4),
  to = c(2, 6)
)
err <- try(.ipp_validate_spatial_surface(gap, 6, 'DetectionSpatial'), silent=TRUE)
stopifnot(inherits(err, 'try-error'), grepl('timestep 3', as.character(err), fixed=TRUE))

# 5. Overlap, non-increasing change points, incompatible geometry and
# multi-layer scheduled grids are rejected.
err <- try(INApestSpatialSchedule(g0, g1, from=c(1,3), to=c(3,5)), silent=TRUE)
stopifnot(inherits(err, 'try-error'))
err <- try(INApestSpatialSchedule(g0, g1, from=c(3,1)), silent=TRUE)
stopifnot(inherits(err, 'try-error'))
g_bad <- INApestSpatialGrid(matrix(1, 3, 2), xmin=0, xmax=2, ymin=0, ymax=2)
err <- try(INApestSpatialSchedule(g0, g_bad, from=c(1,2)), silent=TRUE)
stopifnot(inherits(err, 'try-error'))
a3 <- array(1, dim=c(2,2,2))
g_time <- INApestSpatialGrid(a3, xmin=0, xmax=2, ymin=0, ymax=2)
err <- try(INApestSpatialSchedule(g_time, from=1), silent=TRUE)
stopifnot(inherits(err, 'try-error'))

# 6. Point engine uses scheduled management surfaces transparently.
# No detection in timestep 1; complete detection from timestep 2 onward.
det_sched <- INApestSpatialSchedule(off = g0, on = g1, from = c(1,2))
init <- data.frame(x=rep(1, 40), y=rep(1,40))
fit <- INApestMetaPoint(
  Nperm=1, Ntimesteps=2, InitialPoints=init,
  Survival=1, PropaguleProduction=0, PropaguleEstablishment=1,
  SDDkernel=k0,
  DetectionProb=1, DetectionSD=0, DetectionSpatial=det_sched,
  ManageProb=0, ManageSD=0,
  MortalityProb=0, MortalitySD=0,
  FecundityReduction=0, FecundityReductionSD=0,
  SpreadReduction=0, SpreadReductionSD=0,
  DoProgress=FALSE, Seed=91
)
stopifnot(fit$Summary$n_detected[1] == 0)
stopifnot(fit$Summary$n_detected[2] == 40)

# 7. Habitat accepts the same schedule architecture.
hab_sched <- INApestSpatialSchedule(unsuitable=g0, suitable=g1, from=c(1,2))
stopifnot(.ipp_habitat_value(1,1,hab_sched,1) == 0)
stopifnot(.ipp_habitat_value(1,1,hab_sched,2) == 1)

cat('INApestPoint spatial schedule tests passed.\n')

###############################################################################
### Same schedule interface in the transition-matrix engine
###############################################################################
source('/mnt-data/INApestPointTransitionMatrix.R')

g0 <- INApestSpatialGrid(matrix(0, 2, 2), xmin=0, xmax=2, ymin=0, ymax=2)
g1 <- INApestSpatialGrid(matrix(1, 2, 2), xmin=0, xmax=2, ymin=0, ymax=2)
det_sched <- INApestSpatialSchedule(off=g0, on=g1, from=c(1,2))

fit_tm <- INApestPointTransitionMatrix(
  Nperm=1, Ntimesteps=2, Nstages=2, Transition=diag(2),
  InitialPoints=data.frame(x=rep(1,40), y=rep(1,40), stage=1L),
  SDDkernel=INApestPointKernelFixed(0),
  DetectionProb=1, DetectionSD=0, DetectionSpatial=det_sched,
  ManageProb=0, ManageSD=0,
  MortalityProb=0, MortalitySD=0,
  FecundityReduction=0, FecundityReductionSD=0,
  SpreadReduction=0, SpreadReductionSD=0,
  DoProgress=FALSE, Seed=92
)
stopifnot(fit_tm$Summary$n_detected[1] == 0)
stopifnot(fit_tm$Summary$n_detected[2] == 40)

cat('INApestPointTransitionMatrix spatial schedule tests passed.\n')
cat('ALL SPATIAL SCHEDULE TESTS PASSED\n')
