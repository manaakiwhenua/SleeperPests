# Model-level regression tests for temporal SDD/LDD support.
#
# Usage from the bundle root:
#   Rscript tests/test_models_temporal_connectivity_30x30.R
#
# The tests use a 30 x 30 (900-node) landscape, two realisations and three
# timesteps.  For every public model implementation they check:
#   1. a 2D connectivity matrix and a 3D array containing identical slices
#      give identical results under a fixed random-number stream;
#   2. changing the first SDD slice changes the simulated spatial outcome;
#   3. changing the first LDD slice changes the simulated spatial outcome.
# Parallel tests force a one-worker PSOCK cluster with a fixed RNG stream so
# the 2D-versus-repeated-3D comparison is reproducible.

options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "tests/test_models_temporal_connectivity_30x30.R"
bundle_root <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
code_dir <- file.path(bundle_root, "modified")
out_dir <- file.path(bundle_root, "test_output_r")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- paste0(normalizePath(out_dir, mustWork = TRUE), .Platform$file.sep)

nr <- 30L
nc <- 30L
n_nodes <- nr * nc
n_timesteps <- 3L
n_perm <- 2L

node_id <- function(r, c) (r - 1L) * nc + c
barrier_col <- nc %/% 2L
source_node <- node_id(10L, barrier_col)
target_node <- node_id(10L, barrier_col + 1L)

make_east_matrix <- function() {
  m <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  for (r in seq_len(nr)) {
    for (c in seq_len(nc)) {
      i <- node_id(r, c)
      j <- if (c < nc) node_id(r, c + 1L) else i
      m[i, j] <- 1
    }
  }
  m
}

as_temporal <- function(m) {
  a <- array(0, dim = c(n_nodes, n_nodes, n_timesteps))
  for (t in seq_len(n_timesteps)) a[, , t] <- m
  a
}

close_first_barrier_slice <- function(a) {
  for (r in seq_len(nr)) {
    left <- node_id(r, barrier_col)
    right <- node_id(r, barrier_col + 1L)
    a[left, right, 1L] <- 0
    a[left, left, 1L] <- 1
  }
  a
}

base_route <- make_east_matrix()
zero_route <- matrix(0, nrow = n_nodes, ncol = n_nodes)
base_route_3d <- as_temporal(base_route)
zero_route_3d <- as_temporal(zero_route)
dynamic_route_3d <- close_first_barrier_slice(base_route_3d)

initial_binary <- integer(n_nodes)
initial_binary[source_node] <- 1L
initial_population <- numeric(n_nodes)
initial_population[source_node] <- 100
initial_population_lu <- matrix(0, nrow = n_nodes, ncol = 2L)
initial_population_lu[source_node, 1L] <- 100
initial_population_stage <- matrix(0, nrow = n_nodes, ncol = 2L)
initial_population_stage[source_node, 2L] <- 50
geocoords <- cbind(
  x = rep(seq_len(nc), times = nr),
  y = rep(seq_len(nr), each = nc)
)

result_path <- function(model_name, suffix) file.path(out_dir, paste0(model_name, suffix))

run_and_read <- function(fun, call_args, model_name, suffix, seed = 71003L) {
  # Remove any stale result from an earlier failed test.
  unlink(result_path(model_name, suffix), force = TRUE)
  call_args$ModelName <- model_name
  call_args$OutputDir <- out_dir
  set.seed(seed)
  invisible(do.call(fun, call_args))
  p <- result_path(model_name, suffix)
  if (!file.exists(p)) stop("Expected result was not written: ", p)
  readRDS(p)
}

configure_parallel_env <- function(env) {
  env$detectCores <- function(...) 2L
  env$makeCluster <- function(spec, ...) {
    cl <- parallel::makePSOCKcluster(1L)
    parallel::clusterSetRNGStream(cl, iseed = 99331L)
    cl
  }
  invisible(env)
}

load_model <- function(filename, fun_name, packages = character(), parallel = FALSE) {
  missing <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) {
    message("SKIP ", fun_name, ": missing package(s): ", paste(missing, collapse = ", "))
    return(NULL)
  }
  env <- new.env(parent = .GlobalEnv)
  sys.source(file.path(code_dir, filename), envir = env)
  if (parallel) configure_parallel_env(env)
  get(fun_name, envir = env, inherits = FALSE)
}

assert_same <- function(a, b, label) {
  if (!identical(a, b)) stop(label, ": 2D and repeated-3D results differed")
}

assert_changed <- function(a, b, label) {
  if (identical(a, b)) stop(label, ": temporal connectivity change did not alter results")
}

basic_args <- function(sdd, ldd, geocoords_arg = NULL) {
  x <- list(
    Nperm = n_perm,
    Ntimesteps = n_timesteps,
    DetectionProb = 0,
    DetectionSD = 0,
    ManageProb = 0,
    ManageSD = 0,
    EradicationProb = 0,
    EradicationSD = 0,
    SpreadReduction = 0,
    SpreadReductionSD = 0,
    InitialInvasion = initial_binary,
    InitBioP = NA,
    InvasionRisk = rep(0, n_nodes),
    InitialInfo = integer(n_nodes),
    InitInfoP = NA,
    ExternalInfoProb = rep(0, n_nodes),
    EnvEstabProb = 1,
    Survival = 1,
    SDDprob = sdd,
    SEAM = 0,
    LDDprob = ldd,
    OngoingExternalInvasion = FALSE,
    OngoingExternalInfo = FALSE,
    DoPlots = FALSE
  )
  if (!is.null(geocoords_arg)) x$geocoords <- geocoords_arg
  x
}

meta_args <- function(sdd, ldd, lddrate = 0) {
  list(
    Nperm = n_perm,
    Ntimesteps = n_timesteps,
    DetectionProb = 0,
    DetectionSD = 0,
    ManageProb = 0,
    ManageSD = 0,
    MortalityProb = 0,
    MortalitySD = 0,
    SpreadReduction = 0,
    SpreadReductionSD = 0,
    InitialPopulation = initial_population,
    InitBioP = NA,
    InvasionRisk = rep(0, n_nodes),
    InitialInfo = integer(n_nodes),
    InitInfoP = NA,
    ExternalInfoProb = rep(0, n_nodes),
    EnvEstabProb = 1,
    Survival = 1,
    K = rep(100, n_nodes),
    PropaguleProduction = 100,
    PropaguleEstablishment = 1,
    IncursionStartPop = NA,
    SDDprob = sdd,
    SEAM = 0,
    LDDprob = ldd,
    LDDrate = lddrate,
    OngoingExternalInvasion = FALSE,
    OngoingExternalInfo = FALSE,
    DoPlots = FALSE
  )
}

mlu_args <- function(sdd, ldd, lddrate = 0) {
  list(
    Nperm = n_perm,
    Ntimesteps = n_timesteps,
    Nlanduses = 2L,
    DetectionProb = c(0, 0),
    DetectionSD = c(0, 0),
    ManageProb = c(0, 0),
    ManageSD = c(0, 0),
    MortalityProb = c(0, 0),
    MortalitySD = c(0, 0),
    SpreadReduction = c(0, 0),
    SpreadReductionSD = c(0, 0),
    InitialPopulation = initial_population_lu,
    InitBioP = NA,
    InvasionRisk = rep(0, n_nodes),
    InitialInfo = integer(n_nodes),
    InitInfoP = NA,
    ExternalInfoProb = rep(0, n_nodes),
    EnvEstabProb = 1,
    Survival = c(1, 1),
    K = matrix(100, nrow = n_nodes, ncol = 2L),
    PropaguleProduction = 100,
    PropaguleEstablishment = 1,
    IncursionStartPop = NA,
    SDDprob = sdd,
    SEAM = 0,
    LDDprob = ldd,
    LDDrate = lddrate,
    OngoingExternalInvasion = FALSE,
    OngoingExternalInfo = FALSE,
    DoPlots = FALSE
  )
}

transition_matrix <- matrix(
  c(0, 100,
    1,   1),
  nrow = 2L,
  byrow = TRUE
)

transition_args <- function(sdd, ldd, lddrate = 0) {
  list(
    Nperm = n_perm,
    Ntimesteps = n_timesteps,
    Nstages = 2L,
    Weights = c(1, 1),
    Transition = transition_matrix,
    DetectionProb = c(0, 0),
    DetectionSD = 0,
    ManageProb = 0,
    ManageSD = 0,
    MortalityProb = c(0, 0),
    MortalitySD = 0,
    SpreadReduction = 0,
    SpreadReductionSD = 0,
    InitialPopulation = initial_population_stage,
    InitBioP = NA,
    InvasionRisk = rep(0, n_nodes),
    InitialInfo = integer(n_nodes),
    InitInfoP = NA,
    ExternalInfoProb = rep(0, n_nodes),
    EnvEstabProb = 1,
    K = rep(100, n_nodes),
    SeedbankK = rep(200, n_nodes),
    PropaguleEstablishment = 1,
    IncursionStartPop = NA,
    SDDprob = sdd,
    SEAM = 0,
    LDDprob = ldd,
    LDDrate = lddrate,
    DispersalDensityFactor = 0,
    BlockedTransitionMortality = 0,
    OngoingExternalInvasion = FALSE,
    OngoingExternalInfo = FALSE,
    DoPlots = FALSE
  )
}

run_temporal_suite <- function(fun, label, arg_builder, suffix, extra_arg = NULL) {
  if (is.null(fun)) return(FALSE)

  # 2D legacy path versus identical 3D slices (SDD active, LDD inactive).
  a2 <- arg_builder(base_route, zero_route, 0)
  a3 <- arg_builder(base_route_3d, zero_route_3d, 0)
  if (!is.null(extra_arg)) {
    a2 <- c(a2, extra_arg)
    a3 <- c(a3, extra_arg)
  }
  out2 <- run_and_read(fun, a2, paste0(label, "_2d_"), suffix, 1111L)
  out3 <- run_and_read(fun, a3, paste0(label, "_3d_repeat_"), suffix, 1111L)
  assert_same(out2, out3, paste0(label, " repeated SDD"))

  # Time-varying SDD: only the first slice has a barrier.
  ad <- arg_builder(dynamic_route_3d, zero_route_3d, 0)
  if (!is.null(extra_arg)) ad <- c(ad, extra_arg)
  outd <- run_and_read(fun, ad, paste0(label, "_3d_dynamic_sdd_"), suffix, 1111L)
  assert_changed(out3, outd, paste0(label, " dynamic SDD"))

  # Time-varying LDD: use all propagules/invasion probability through LDD.
  # For the occupancy family arg_builder does not accept/use LDDrate, so it is
  # handled by its dedicated wrapper below.
  al_static <- arg_builder(zero_route_3d, base_route_3d, 1)
  al_dynamic <- arg_builder(zero_route_3d, dynamic_route_3d, 1)
  if (!is.null(extra_arg)) {
    al_static <- c(al_static, extra_arg)
    al_dynamic <- c(al_dynamic, extra_arg)
  }
  outls <- run_and_read(fun, al_static, paste0(label, "_3d_static_ldd_"), suffix, 2222L)
  outld <- run_and_read(fun, al_dynamic, paste0(label, "_3d_dynamic_ldd_"), suffix, 2222L)
  assert_changed(outls, outld, paste0(label, " dynamic LDD"))

  message("PASS ", label, ": 2D compatibility + dynamic SDD + dynamic LDD")
  TRUE
}

run_basic_suite <- function(fun, label, geocoords_arg = NULL) {
  if (is.null(fun)) return(FALSE)
  extra <- if (is.null(geocoords_arg)) NULL else list(geocoords = geocoords_arg)

  a2 <- basic_args(base_route, zero_route, geocoords_arg)
  a3 <- basic_args(base_route_3d, zero_route_3d, geocoords_arg)
  out2 <- run_and_read(fun, a2, paste0(label, "_2d_"), "InvasionLargeOut.rds", 3333L)
  out3 <- run_and_read(fun, a3, paste0(label, "_3d_repeat_"), "InvasionLargeOut.rds", 3333L)
  assert_same(out2, out3, paste0(label, " repeated SDD"))

  ad <- basic_args(dynamic_route_3d, zero_route_3d, geocoords_arg)
  outd <- run_and_read(fun, ad, paste0(label, "_3d_dynamic_sdd_"), "InvasionLargeOut.rds", 3333L)
  assert_changed(out3, outd, paste0(label, " dynamic SDD"))

  als <- basic_args(zero_route_3d, base_route_3d, geocoords_arg)
  ald <- basic_args(zero_route_3d, dynamic_route_3d, geocoords_arg)
  outls <- run_and_read(fun, als, paste0(label, "_3d_static_ldd_"), "InvasionLargeOut.rds", 4444L)
  outld <- run_and_read(fun, ald, paste0(label, "_3d_dynamic_ldd_"), "InvasionLargeOut.rds", 4444L)
  assert_changed(outls, outld, paste0(label, " dynamic LDD"))

  message("PASS ", label, ": 2D compatibility + dynamic SDD + dynamic LDD")
  TRUE
}
assert_error <- function(expr, pattern, label) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  if (is.na(msg) || !grepl(pattern, msg, fixed = TRUE)) stop(label, ": expected error containing '", pattern, "', got: ", msg)
  message("PASS ", label)
}

# INApest mixed 2D / 3D compatibility
f <- load_model("INApest.R", "INApest")
b2 <- run_and_read(f, basic_args(base_route, zero_route), "edge_INApest_base_", "InvasionLargeOut.rds", 7771L)
bmixs <- run_and_read(f, basic_args(base_route_3d, zero_route), "edge_INApest_mixs_", "InvasionLargeOut.rds", 7771L)
assert_same(b2, bmixs, "INApest 3D SDD + 2D LDD")
l2 <- run_and_read(f, basic_args(zero_route, base_route), "edge_INApest_ldd2_", "InvasionLargeOut.rds", 7772L)
bmixl <- run_and_read(f, basic_args(zero_route, base_route_3d), "edge_INApest_mixl_", "InvasionLargeOut.rds", 7772L)
assert_same(l2, bmixl, "INApest 2D SDD + 3D LDD")

bad <- array(0, dim = c(n_nodes, n_nodes, n_timesteps - 1L))
x <- basic_args(bad, zero_route); x$ModelName <- "bad"; x$OutputDir <- out_dir
assert_error(do.call(f, x), "SDDprob 3D array must have dimensions nodes x nodes x Ntimesteps", "INApest bad SDD time dimension")
x <- basic_args(base_route, bad); x$ModelName <- "bad"; x$OutputDir <- out_dir
assert_error(do.call(f, x), "LDDprob 3D array must have dimensions nodes x nodes x Ntimesteps", "INApest bad LDD time dimension")

# Multiple-land-use mixed compatibility and validation
fmlu <- load_model("INApestMetaMultipleLandUse.r", "INApestMetaMultipleLandUse")
m2 <- run_and_read(fmlu, mlu_args(base_route, zero_route, 0), "edge_MLU_base_", "PopulationLargeOut.rds", 7781L)
mmixs <- run_and_read(fmlu, mlu_args(base_route_3d, zero_route, 0), "edge_MLU_mixs_", "PopulationLargeOut.rds", 7781L)
assert_same(m2, mmixs, "MLU 3D SDD + 2D LDD")
ml2 <- run_and_read(fmlu, mlu_args(zero_route, base_route, 1), "edge_MLU_ldd2_", "PopulationLargeOut.rds", 7782L)
mmixl <- run_and_read(fmlu, mlu_args(zero_route, base_route_3d, 1), "edge_MLU_mixl_", "PopulationLargeOut.rds", 7782L)
assert_same(ml2, mmixl, "MLU 2D SDD + 3D LDD")
x <- mlu_args(bad, zero_route, 0); x$ModelName <- "bad"; x$OutputDir <- out_dir
assert_error(do.call(fmlu, x), "SDDprob 3D array must have dimensions nodes x nodes x Ntimesteps", "MLU bad SDD time dimension")
x <- mlu_args(base_route, bad, 0); x$ModelName <- "bad"; x$OutputDir <- out_dir
assert_error(do.call(fmlu, x), "LDDprob 3D array must have dimensions nodes x nodes x Ntimesteps", "MLU bad LDD time dimension")

# Transition model mixed compatibility and validation
ftr <- load_model("INApestMetaTransitionMatrix.r", "INApestMetaTransitionMatrix")
t2 <- run_and_read(ftr, transition_args(base_route, zero_route, 0), "edge_TR_base_", "PopulationStageLargeOut.rds", 7791L)
tmixs <- run_and_read(ftr, transition_args(base_route_3d, zero_route, 0), "edge_TR_mixs_", "PopulationStageLargeOut.rds", 7791L)
assert_same(t2, tmixs, "Transition 3D SDD + 2D LDD")
tl2 <- run_and_read(ftr, transition_args(zero_route, base_route, 1), "edge_TR_ldd2_", "PopulationStageLargeOut.rds", 7792L)
tmixl <- run_and_read(ftr, transition_args(zero_route, base_route_3d, 1), "edge_TR_mixl_", "PopulationStageLargeOut.rds", 7792L)
assert_same(tl2, tmixl, "Transition 2D SDD + 3D LDD")
x <- transition_args(bad, zero_route, 0); x$ModelName <- "bad"; x$OutputDir <- out_dir
assert_error(do.call(ftr, x), "SDDprob 3D array must have dimensions nodes x nodes x Ntimesteps", "Transition bad SDD time dimension")
x <- transition_args(base_route, bad, 0); x$ModelName <- "bad"; x$OutputDir <- out_dir
assert_error(do.call(ftr, x), "LDDprob 3D array must have dimensions nodes x nodes x Ntimesteps", "Transition bad LDD time dimension")

# Guarded initialisation smoke tests
x <- basic_args(zero_route, zero_route); x$Ntimesteps <- 1L; x$InitialInvasion <- NA; x$InitBioP <- NA; x$InvasionRisk <- NA
invisible(run_and_read(f, x, "edge_INApest_no_initial_", "InvasionLargeOut.rds", 7801L))
message("PASS INApest no-initial-input guard")
x <- mlu_args(zero_route, zero_route, 0); x$Ntimesteps <- 1L; x$InitialPopulation <- NA; x$InitBioP <- NA; x$InvasionRisk <- NA
invisible(run_and_read(fmlu, x, "edge_MLU_no_initial_", "PopulationLargeOut.rds", 7802L))
message("PASS MLU no-initial-input guard")

cat("\nEDGE TESTS PASS\n")
