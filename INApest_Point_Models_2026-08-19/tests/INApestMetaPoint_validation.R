###############################################################################
### INApestMetaPoint_validation.R
### Neutral demographic + dispersal regression test for INApestMetaPoint
###
### Validation design
### -----------------
### * one biological entity per point
### * 3,000 independent founder lineages, all starting at (0, 0)
### * no mortality, management, detection, information, habitat effects or LDD
### * Poisson propagule production with lambda = 0.5 per point per timestep
### * all propagules establish
### * exponential radial SDD kernel with mean distance = 100 spatial units
###
### Because there are no interactions among points in this neutral case, the
### 3,000 founder lineages are independent Monte Carlo replicates even though
### they are simulated in a single INApestMetaPoint permutation.
###############################################################################

source("/mnt-data/INApestMetaPoint.R")

N0 <- 3000L
Ntimesteps <- 7L
lambda <- 0.5
mean_distance <- 100

initial <- data.frame(
  x = rep(0, N0),
  y = rep(0, N0)
)

fit <- INApestMetaPoint(
  ModelName = "PointValidationBenchmark",
  Nperm = 1,
  Ntimesteps = Ntimesteps,
  InitialPoints = initial,

  Survival = 1,
  PropaguleProduction = lambda,
  PropaguleEstablishment = 1,
  EnvEstabProb = 1,

  SDDkernel = INApestPointKernelExponential(mean_distance),
  LDDkernel = NULL,
  LDDrate = 0,

  HabitatSuitability = NULL,
  HabitatSearchRadius = 0,

  DetectionProb = 0,
  DetectionSD = 0,
  ManageProb = 0,
  ManageSD = 0,
  MortalityProb = 0,
  MortalitySD = 0,
  SpreadReduction = 0,
  SpreadReductionSD = 0,

  InfoRadius = 0,
  OngoingExternalInvasion = FALSE,

  DoProgress = TRUE,
  Seed = 20260818
)

###############################################################################
### Recover independent founder lineages from parent_id
###############################################################################

fp <- fit$FinalPoints[order(fit$FinalPoints$id), ]
root <- integer(max(fp$id))

for (i in seq_len(nrow(fp))) {
  id <- fp$id[i]
  pid <- fp$parent_id[i]
  root[id] <- if (is.na(pid)) id else root[pid]
}

###############################################################################
### 1. Population growth
###
### With survival = 1 and Poisson(lambda) propagules per existing point:
###
###   N[t+1] = N[t] + Poisson(lambda * N[t])
###
### This is a Galton-Watson process with offspring distribution
### 1 + Poisson(lambda), hence
###
###   E[N_t]   = (1 + lambda)^t
###
###   Var[N_t] = lambda * m^(t-1) * (m^t - 1) / (m - 1), m = 1 + lambda
###############################################################################

m <- 1 + lambda
growth_rows <- vector("list", Ntimesteps)

for (tt in seq_len(Ntimesteps)) {

  ids <- fit$PointHistory$id[
    fit$PointHistory$timestep == tt
  ]

  lineage_n <- tabulate(
    root[ids],
    nbins = N0
  )

  q <- quantile(
    lineage_n,
    c(0.05, 0.5, 0.95),
    names = FALSE
  )

  mean_theory <- m^tt

  var_theory <- lambda *
    m^(tt - 1) *
    (m^tt - 1) /
    (m - 1)

  growth_rows[[tt]] <- data.frame(
    timestep = tt,
    mean_sim = mean(lineage_n),
    mean_theory = mean_theory,
    mean_rel_error = (mean(lineage_n) - mean_theory) / mean_theory,
    var_sim = var(lineage_n),
    var_theory = var_theory,
    var_rel_error = (var(lineage_n) - var_theory) / var_theory,
    q05 = q[1],
    median = q[2],
    q95 = q[3]
  )
}

growth_validation <- do.call(
  rbind,
  growth_rows
)

###############################################################################
### 2. Dispersal-kernel validation
###############################################################################

ev <- subset(
  fit$EventLog,
  event == "establishment"
)

parent_x <- fp$x[
  match(ev$parent_id, fp$id)
]

parent_y <- fp$y[
  match(ev$parent_id, fp$id)
]

step_distance <- sqrt(
  (ev$x - parent_x)^2 +
  (ev$y - parent_y)^2
)

emp_q <- quantile(
  step_distance,
  c(0.5, 0.9, 0.95, 0.99),
  names = FALSE
)

theory_q <- qexp(
  c(0.5, 0.9, 0.95, 0.99),
  rate = 1 / mean_distance
)

dispersal_validation <- data.frame(
  statistic = c("mean", "q50", "q90", "q95", "q99"),
  simulated = c(mean(step_distance), emp_q),
  theoretical = c(mean_distance, theory_q)
)

dispersal_validation$relative_error <-
  (dispersal_validation$simulated - dispersal_validation$theoretical) /
  dispersal_validation$theoretical

###############################################################################
### 3. Emergent final spatial footprint
###############################################################################

fp$root <- root[fp$id]

fp$radius_from_origin <- sqrt(
  fp$x^2 +
  fp$y^2
)

max_radius <- tapply(
  fp$radius_from_origin,
  factor(fp$root, levels = seq_len(N0)),
  max
)

lineage_n_final <- tabulate(
  fp$root,
  nbins = N0
)

final_realisation_distribution <- data.frame(
  root_id = seq_len(N0),
  n_final = lineage_n_final,
  max_radius = as.numeric(max_radius)
)

q_n <- quantile(
  lineage_n_final,
  c(0.05, 0.5, 0.95),
  names = FALSE
)

q_r <- quantile(
  max_radius,
  c(0.05, 0.5, 0.95),
  names = FALSE
)

final_summary <- data.frame(
  metric = c(
    "final abundance mean",
    "final abundance q05",
    "final abundance median",
    "final abundance q95",
    "maximum radius mean",
    "maximum radius q05",
    "maximum radius median",
    "maximum radius q95"
  ),
  value = c(
    mean(lineage_n_final),
    q_n[1], q_n[2], q_n[3],
    mean(max_radius),
    q_r[1], q_r[2], q_r[3]
  )
)

###############################################################################
### Regression-test tolerances
###############################################################################

stopifnot(
  max(abs(growth_validation$mean_rel_error)) < 0.03,
  max(abs(growth_validation$var_rel_error)) < 0.10,
  max(abs(dispersal_validation$relative_error)) < 0.05
)

###############################################################################
### Output
###############################################################################

write.csv(
  growth_validation,
  "growth_validation.csv",
  row.names = FALSE
)

write.csv(
  dispersal_validation,
  "dispersal_validation.csv",
  row.names = FALSE
)

write.csv(
  final_realisation_distribution,
  "final_realisation_distribution.csv",
  row.names = FALSE
)

write.csv(
  final_summary,
  "final_summary.csv",
  row.names = FALSE
)

cat("\nPopulation-growth validation\n")
print(growth_validation, row.names = FALSE)

cat("\nDispersal-kernel validation\n")
print(dispersal_validation, row.names = FALSE)

cat("\nFinal spatial footprint\n")
print(final_summary, row.names = FALSE)

cat("\nAll neutral validation checks passed.\n")
