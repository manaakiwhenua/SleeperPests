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
# Test-only sequential shims for package/PSOCK infrastructure unavailable in WebR.
# Model source files are not changed.
abind_test <- function(..., along) {
  xs <- list(...)
  d <- dim(xs[[1]])
  if (is.null(d) || along != length(d) + 1L) stop("test abind shim only supports binding on a new final dimension")
  out <- array(NA, dim = c(d, length(xs)))
  out[] <- unlist(xs, use.names = FALSE)
  out
}
foreach_test <- function(..., .combine, .packages = character()) {
  dots <- as.list(substitute(list(...)))[-1L]
  list(iter_expr = dots[[1L]], combine = .combine, env = parent.frame())
}
dopar_test <- function(obj, expr) {
  block <- substitute(expr)
  vals <- eval(obj$iter_expr, obj$env)
  ans <- lapply(vals, function(v) eval(block, envir = new.env(parent = obj$env)))
  comb <- get(obj$combine, envir = obj$env, inherits = TRUE)
  do.call(comb, ans)
}

ina_test <- function(..., initinfo, initbio, bpam, probestabvec, probadoptvec) {
  incoming <- as.numeric(as.vector(initbio) %*% bpam)
  estab <- as.logical((as.vector(initbio) > 0 & as.vector(probestabvec) > 0) | incoming > 0)
  step <- list(vect1cL = list(NULL, as.vector(initinfo)), estabvecL = list(estab))
  list(multdetails = list(list(multout = list(step))))
}

load_parallel_shim <- function(filename, fun_name, ina = FALSE, transition = FALSE) {
  env <- new.env(parent = .GlobalEnv)
  env$abind <- abind_test
  env$foreach <- foreach_test
  env$`%dopar%` <- dopar_test
  env$detectCores <- function(...) 2L
  env$makeCluster <- function(...) structure(list(), class = "testCluster")
  env$registerDoParallel <- function(...) invisible(NULL)
  env$stopCluster <- function(...) invisible(NULL)
  env$clusterExport <- function(...) invisible(NULL)
  env$clusterEvalQ <- function(...) invisible(NULL)
  env$parLapply <- function(cl, X, fun, ...) lapply(X, fun, ...)
  if (ina) env$INAscene <- ina_test

  lines <- readLines(file.path(code_dir, filename), warn = FALSE)
  lines <- lines[!grepl("^[[:space:]]*library\\((abind|doParallel|INA)\\)[[:space:]]*$", lines)]
  if (transition) lines <- gsub("parallel::stopCluster(cl)", "stopCluster(cl)", lines, fixed = TRUE)
  eval(parse(text = lines, keep.source = FALSE), envir = env)
  get(fun_name, envir = env, inherits = FALSE)
}

results <- logical()
unlink(file.path(bundle_root, "shim_progress.log"))
mark <- function(x) cat(x, "\n", file=file.path(bundle_root, "shim_progress.log"), append=TRUE)
mark("START INApestParallel")
f <- load_parallel_shim("INApestParallel.R", "INApestParallel")
results["INApestParallel"] <- run_basic_suite(f, "shim_INApestParallel")

mark("DONE INApestParallel; START INAScene")
f <- load_parallel_shim("INApestParallelInAScene.R", "INApestParallelINAscene", ina = TRUE)
results["INApestParallelINAscene"] <- run_basic_suite(f, "shim_INApestParallelINAscene", geocoords)

mark("DONE INAScene; START MetaParallel")
f <- load_parallel_shim("INApestMetaParallel.r", "INApestMetaParallel")
results["INApestMetaParallel"] <- run_temporal_suite(f, "shim_INApestMetaParallel", meta_args, "PopulationLargeOut.rds")

mark("DONE MetaParallel; START MetaParallelMLU")
f <- load_parallel_shim("INApestMetaParallelMultipleLandUse.r", "INApestMetaParallelMultipleLandUse")
results["INApestMetaParallelMultipleLandUse"] <- run_temporal_suite(f, "shim_INApestMetaParallelMultipleLandUse", mlu_args, "PopulationLargeOut.rds")

mark("DONE MetaParallelMLU; START TransitionParallel")
f <- load_parallel_shim("INApestMetaTransitionMatrixParallel.r", "INApestMetaTransitionMatrixParallel", transition = TRUE)
results["INApestMetaTransitionMatrixParallel"] <- run_temporal_suite(f, "shim_INApestMetaTransitionMatrixParallel", transition_args, "PopulationStageLargeOut.rds")

mark("DONE TransitionParallel")
writeLines(sprintf("%s %s", names(results), ifelse(results, "PASS", "FAIL")), file.path(bundle_root, "shim_results.txt"))
cat("\nPARALLEL-BODY SHIM SUMMARY\n")
for (nm in names(results)) cat(sprintf("%-42s %s\n", nm, if (isTRUE(results[[nm]])) "PASS" else "FAIL"))
