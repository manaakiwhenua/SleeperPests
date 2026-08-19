# Model-level tests for node information decay / retention.
#
# The tests use a 30 x 30 (900-node) landscape, two realisations and three
# timesteps. Every public model implementation is checked for:
#   1. exact backward compatibility when InfoRetentionProb is omitted (= 1);
#   2. stochastic information decay with scalar InfoRetentionProb = 0.5;
#   3. node-specific retention using a vector of length nodes;
#   4. node-by-timestep retention using a nodes x timesteps matrix;
#   5. information refresh after retention through existing detection code;
#   6. rejection of an invalid retention vector length.
# Parallel implementations are exercised using the same test-only sequential
# shims used for the temporal-connectivity regression because WebR lacks the
# external parallel/INA package stack. Production source is not rewritten.

options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
bundle_root <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
baseline_dir <- file.path(bundle_root, "baseline")
modified_dir <- file.path(bundle_root, "modified")
out_dir <- file.path(bundle_root, "test_output_info_decay")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- paste0(normalizePath(out_dir, mustWork = TRUE), .Platform$file.sep)

nr <- 30L
nc <- 30L
n_nodes <- nr * nc
n_timesteps <- 3L
n_perm <- 2L
n_landuses <- 2L
source_node <- 1L

zero_route <- matrix(0, nrow = n_nodes, ncol = n_nodes)
initial_binary <- integer(n_nodes)
initial_binary[source_node] <- 1L
initial_population <- numeric(n_nodes)
initial_population[source_node] <- 20
initial_population_lu <- matrix(0, nrow = n_nodes, ncol = n_landuses)
initial_population_lu[source_node, 1L] <- 20
initial_population_stage <- matrix(0, nrow = n_nodes, ncol = 2L)
initial_population_stage[source_node, 2L] <- 20
all_info <- rep(1, n_nodes)
zero_external <- rep(0, n_nodes)
geocoords <- cbind(
  x = rep(seq_len(nc), times = nr),
  y = rep(seq_len(nr), each = nc)
)

# Use a stable two-stage matrix with no propagule production in the first row
# and persistence of stage 2. The information tests do not depend on spread.
transition_matrix <- matrix(
  c(0, 0,
    0, 1),
  nrow = 2L,
  byrow = TRUE
)

# ----------------------- test infrastructure -----------------------
result_path <- function(model_name, suffix) file.path(out_dir, paste0(model_name, suffix))

run_and_read <- function(fun, call_args, model_name, suffixes, seed = 20260813L) {
  for (s in suffixes) unlink(result_path(model_name, s), force = TRUE)
  call_args$ModelName <- model_name
  call_args$OutputDir <- out_dir
  set.seed(seed)
  invisible(do.call(fun, call_args))
  out <- lapply(suffixes, function(s) {
    p <- result_path(model_name, s)
    if (!file.exists(p)) stop("Expected result was not written: ", p)
    readRDS(p)
  })
  names(out) <- suffixes
  out
}

assert_identical_outputs <- function(a, b, label) {
  if (!identical(names(a), names(b))) stop(label, ": output names differ")
  for (nm in names(a)) {
    if (!identical(a[[nm]], b[[nm]])) stop(label, ": output differed for ", nm)
  }
  TRUE
}

management_profile <- function(x, mlu = FALSE) {
  if (!mlu) {
    if (length(dim(x)) != 3L) stop("Unexpected management result dimensions: ", paste(dim(x), collapse = " x "))
    return(vapply(seq_len(n_timesteps), function(t) sum(x[, t, , drop = FALSE]) / n_perm, numeric(1)))
  }
  if (length(dim(x)) != 4L) stop("Unexpected MLU management result dimensions: ", paste(dim(x), collapse = " x "))
  vapply(seq_len(n_timesteps), function(t) sum(x[, , t, , drop = FALSE]) / n_perm, numeric(1))
}

assert_profile <- function(actual, expected, label, tol = 1e-12) {
  if (length(actual) != length(expected) || any(abs(actual - expected) > tol)) {
    stop(label, ": expected ", paste(expected, collapse = ", "),
         " but got ", paste(round(actual, 3), collapse = ", "))
  }
  TRUE
}

expect_error <- function(expr, pattern, label) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  if (is.na(msg)) stop(label, ": expected an error but none was raised")
  if (!grepl(pattern, msg, fixed = TRUE)) stop(label, ": wrong error: ", msg)
  TRUE
}

# Test-only sequential shims for package/PSOCK infrastructure unavailable in WebR.
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

# Minimal INAscene test double. It preserves initial information and infestation
# and adds any deterministic BPAM arrivals. SEAM is zero in the decay tests.
ina_test <- function(..., initinfo, initbio, bpam, probestabvec, probadoptvec) {
  incoming <- as.numeric(as.vector(initbio) %*% bpam)
  estab <- as.logical((as.vector(initbio) > 0 & as.vector(probestabvec) > 0) | incoming > 0)
  step <- list(vect1cL = list(NULL, as.vector(initinfo)), estabvecL = list(estab))
  list(multdetails = list(list(multout = list(step))))
}

load_model <- function(code_dir, filename, fun_name, parallel = FALSE, ina = FALSE, transition_parallel = FALSE) {
  env <- new.env(parent = .GlobalEnv)
  if (parallel) {
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
    if (transition_parallel) lines <- gsub("parallel::stopCluster(cl)", "stopCluster(cl)", lines, fixed = TRUE)
    eval(parse(text = lines, keep.source = FALSE), envir = env)
  } else {
    sys.source(file.path(code_dir, filename), envir = env)
  }
  get(fun_name, envir = env, inherits = FALSE)
}

# ----------------------- argument builders -----------------------
basic_args <- function(detection = 0) {
  list(
    Nperm = n_perm,
    Ntimesteps = n_timesteps,
    DetectionProb = detection,
    DetectionSD = 0,
    ManageProb = 1,
    ManageSD = 0,
    EradicationProb = 0,
    EradicationSD = 0,
    SpreadReduction = 0,
    SpreadReductionSD = 0,
    InitialInvasion = initial_binary,
    InitBioP = NA,
    InvasionRisk = zero_external,
    InitialInfo = all_info,
    InitInfoP = NA,
    ExternalInfoProb = zero_external,
    EnvEstabProb = 1,
    Survival = 1,
    SDDprob = zero_route,
    SEAM = 0,
    LDDprob = zero_route,
    OngoingExternalInvasion = FALSE,
    OngoingExternalInfo = FALSE,
    DoPlots = FALSE
  )
}

meta_args <- function(detection = 0) {
  list(
    Nperm = n_perm,
    Ntimesteps = n_timesteps,
    DetectionProb = detection,
    DetectionSD = 0,
    ManageProb = 1,
    ManageSD = 0,
    MortalityProb = 0,
    MortalitySD = 0,
    SpreadReduction = 0,
    SpreadReductionSD = 0,
    InitialPopulation = initial_population,
    InitBioP = NA,
    InvasionRisk = zero_external,
    InitialInfo = all_info,
    InitInfoP = NA,
    ExternalInfoProb = zero_external,
    EnvEstabProb = 1,
    Survival = 1,
    K = rep(100, n_nodes),
    PropaguleProduction = 0,
    PropaguleEstablishment = 1,
    IncursionStartPop = NA,
    SDDprob = zero_route,
    SEAM = 0,
    LDDprob = zero_route,
    LDDrate = 0,
    OngoingExternalInvasion = FALSE,
    OngoingExternalInfo = FALSE,
    DoPlots = FALSE
  )
}

mlu_args <- function(detection = c(0, 0)) {
  list(
    Nperm = n_perm,
    Ntimesteps = n_timesteps,
    Nlanduses = n_landuses,
    DetectionProb = detection,
    DetectionSD = c(0, 0),
    ManageProb = c(1, 1),
    ManageSD = c(0, 0),
    MortalityProb = c(0, 0),
    MortalitySD = c(0, 0),
    SpreadReduction = c(0, 0),
    SpreadReductionSD = c(0, 0),
    InitialPopulation = initial_population_lu,
    InitBioP = NA,
    InvasionRisk = zero_external,
    InitialInfo = all_info,
    InitInfoP = NA,
    ExternalInfoProb = zero_external,
    EnvEstabProb = 1,
    Survival = c(1, 1),
    K = matrix(100, nrow = n_nodes, ncol = n_landuses),
    PropaguleProduction = 0,
    PropaguleEstablishment = 1,
    IncursionStartPop = NA,
    SDDprob = zero_route,
    SEAM = 0,
    LDDprob = zero_route,
    LDDrate = 0,
    OngoingExternalInvasion = FALSE,
    OngoingExternalInfo = FALSE,
    DoPlots = FALSE
  )
}

transition_args <- function(detection = c(0, 0)) {
  list(
    Nperm = n_perm,
    Ntimesteps = n_timesteps,
    Nstages = 2L,
    Weights = c(1, 1),
    Transition = transition_matrix,
    DetectionProb = detection,
    DetectionSD = 0,
    ManageProb = 1,
    ManageSD = 0,
    MortalityProb = c(0, 0),
    MortalitySD = 0,
    SpreadReduction = 0,
    SpreadReductionSD = 0,
    InitialPopulation = initial_population_stage,
    InitBioP = NA,
    InvasionRisk = zero_external,
    InitialInfo = all_info,
    InitInfoP = NA,
    ExternalInfoProb = zero_external,
    EnvEstabProb = 1,
    K = rep(100, n_nodes),
    SeedbankK = rep(100, n_nodes),
    PropaguleEstablishment = 1,
    IncursionStartPop = NA,
    SDDprob = zero_route,
    SEAM = 0,
    LDDprob = zero_route,
    LDDrate = 0,
    DispersalDensityFactor = 0,
    BlockedTransitionMortality = 0,
    OngoingExternalInvasion = FALSE,
    OngoingExternalInfo = FALSE,
    DoPlots = FALSE
  )
}

# ----------------------- suite -----------------------
model_specs <- list(
  list(label = "INApest", file = "INApest.R", fun = "INApest", builder = basic_args,
       suffixes = c("InfoLargeOut.rds", "InvasionLargeOut.rds", "DetectedLargeOut.rds"), manage_suffix = "InfoLargeOut.rds", mlu = FALSE),
  list(label = "INApestParallel", file = "INApestParallel.R", fun = "INApestParallel", builder = basic_args,
       suffixes = c("ManagingLargeOut.rds", "InvasionLargeOut.rds", "DetectedLargeOut.rds"), manage_suffix = "ManagingLargeOut.rds", mlu = FALSE, parallel = TRUE),
  list(label = "INApestParallelINAscene", file = "INApestParallelInAScene.R", fun = "INApestParallelINAscene", builder = basic_args,
       suffixes = c("ManagingLargeOut.rds", "InvasionLargeOut.rds", "DetectedLargeOut.rds"), manage_suffix = "ManagingLargeOut.rds", mlu = FALSE, parallel = TRUE, ina = TRUE),
  list(label = "INApestMetaParallel", file = "INApestMetaParallel.r", fun = "INApestMetaParallel", builder = meta_args,
       suffixes = c("InfoLargeOut.rds", "PopulationLargeOut.rds", "InvasionLargeOut.rds", "DetectedLargeOut.rds"), manage_suffix = "InfoLargeOut.rds", mlu = FALSE, parallel = TRUE),
  list(label = "INApestMetaMultipleLandUse", file = "INApestMetaMultipleLandUse.r", fun = "INApestMetaMultipleLandUse", builder = mlu_args,
       suffixes = c("InfoLargeOut.rds", "PopulationLargeOut.rds", "InvasionLargeOut.rds", "DetectedLargeOut.rds"), manage_suffix = "InfoLargeOut.rds", mlu = TRUE),
  list(label = "INApestMetaParallelMultipleLandUse", file = "INApestMetaParallelMultipleLandUse.r", fun = "INApestMetaParallelMultipleLandUse", builder = mlu_args,
       suffixes = c("InfoLargeOut.rds", "PopulationLargeOut.rds", "InvasionLargeOut.rds", "DetectedLargeOut.rds"), manage_suffix = "InfoLargeOut.rds", mlu = TRUE, parallel = TRUE),
  list(label = "INApestMetaTransitionMatrix", file = "INApestMetaTransitionMatrix.r", fun = "INApestMetaTransitionMatrix", builder = transition_args,
       suffixes = c("InfoLargeOut.rds", "PopulationLargeOut.rds", "PopulationStageLargeOut.rds", "InvasionLargeOut.rds", "DetectedLargeOut.rds"), manage_suffix = "InfoLargeOut.rds", mlu = FALSE),
  list(label = "INApestMetaTransitionMatrixParallel", file = "INApestMetaTransitionMatrixParallel.r", fun = "INApestMetaTransitionMatrixParallel", builder = transition_args,
       suffixes = c("InfoLargeOut.rds", "PopulationLargeOut.rds", "PopulationStageLargeOut.rds", "InvasionLargeOut.rds", "DetectedLargeOut.rds"), manage_suffix = "InfoLargeOut.rds", mlu = FALSE, parallel = TRUE, transition_parallel = TRUE)
)

summary_rows <- list()
progress_file <- file.path(bundle_root, "info_decay_progress.txt")
unlink(progress_file)
mark_progress <- function(x) cat(x, "\n", file = progress_file, append = TRUE)

for (spec in model_specs) {
  label <- spec$label
  cat("\n===", label, "===\n")
  mark_progress(paste("START", label))
  is_parallel <- isTRUE(spec$parallel)
  baseline_fun <- load_model(baseline_dir, spec$file, spec$fun, parallel = is_parallel,
                             ina = isTRUE(spec$ina), transition_parallel = isTRUE(spec$transition_parallel))
  modified_fun <- load_model(modified_dir, spec$file, spec$fun, parallel = is_parallel,
                             ina = isTRUE(spec$ina), transition_parallel = isTRUE(spec$transition_parallel))

  # Legacy/default compatibility: omission of InfoRetentionProb must be exact.
  base_args <- spec$builder()
  if (isTRUE(spec$ina)) base_args$geocoords <- geocoords
  old_out <- run_and_read(baseline_fun, base_args, paste0(label, "_baseline_"), spec$suffixes, 10101L)
  new_out <- run_and_read(modified_fun, base_args, paste0(label, "_retention1_"), spec$suffixes, 10101L)
  assert_identical_outputs(old_out, new_out, paste0(label, " backward compatibility"))
  p_no <- management_profile(new_out[[spec$manage_suffix]], spec$mlu)
  mark_progress(paste("COMPAT", label))

  # Scalar stochastic retention.
  scalar_args <- spec$builder()
  if (isTRUE(spec$ina)) scalar_args$geocoords <- geocoords
  scalar_args$InfoRetentionProb <- 0.5
  scalar_out <- run_and_read(modified_fun, scalar_args, paste0(label, "_scalar05_"), spec$suffixes, 20202L)
  p_scalar <- management_profile(scalar_out[[spec$manage_suffix]], spec$mlu)
  mark_progress(paste("SCALAR", label))
  max_units <- n_nodes * if (spec$mlu) n_landuses else 1L
  if (!(abs(p_scalar[1] - max_units) < 1e-12 && p_scalar[2] < p_scalar[1] && p_scalar[3] < p_scalar[2] && p_scalar[3] > 0)) {
    stop(label, ": scalar 0.5 retention did not produce the expected declining management profile: ", paste(p_scalar, collapse = ", "))
  }

  # Node-specific retention: first half retains permanently, second half loses after timestep 1.
  node_retention <- c(rep(1, n_nodes / 2L), rep(0, n_nodes / 2L))
  node_args <- spec$builder()
  if (isTRUE(spec$ina)) node_args$geocoords <- geocoords
  node_args$InfoRetentionProb <- node_retention
  node_out <- run_and_read(modified_fun, node_args, paste0(label, "_nodevector_"), spec$suffixes, 30303L)
  p_node <- management_profile(node_out[[spec$manage_suffix]], spec$mlu)
  mark_progress(paste("NODE", label))
  assert_profile(p_node, c(max_units, max_units / 2, max_units / 2), paste0(label, " node retention"))

  # Node x timestep matrix: retain after t1, lose after t2, so management stops at t3.
  time_retention <- matrix(1, nrow = n_nodes, ncol = n_timesteps)
  time_retention[, 2L] <- 0
  time_args <- spec$builder()
  if (isTRUE(spec$ina)) time_args$geocoords <- geocoords
  time_args$InfoRetentionProb <- time_retention
  time_out <- run_and_read(modified_fun, time_args, paste0(label, "_nodetime_"), spec$suffixes, 40404L)
  p_time <- management_profile(time_out[[spec$manage_suffix]], spec$mlu)
  mark_progress(paste("TIME", label))
  assert_profile(p_time, c(max_units, max_units, 0), paste0(label, " node x timestep retention"))

  # Refresh test: all information is dropped after each step, but the single
  # persistent infestation is detected with probability 1 and is therefore
  # informed again before the next timestep.
  if (grepl("MultipleLandUse", label, fixed = TRUE)) refresh_args <- spec$builder(c(1, 1))
  else if (grepl("TransitionMatrix", label, fixed = TRUE)) refresh_args <- spec$builder(c(1, 1))
  else refresh_args <- spec$builder(1)
  if (isTRUE(spec$ina)) refresh_args$geocoords <- geocoords
  refresh_args$InfoRetentionProb <- 0
  refresh_out <- run_and_read(modified_fun, refresh_args, paste0(label, "_refresh_"), spec$suffixes, 50505L)
  p_refresh <- management_profile(refresh_out[[spec$manage_suffix]], spec$mlu)
  mark_progress(paste("REFRESH", label))
  source_units <- if (spec$mlu) n_landuses else 1L
  assert_profile(p_refresh, c(max_units, source_units, source_units), paste0(label, " detection refresh"))

  # Invalid vector length should be rejected before simulation.
  bad_args <- spec$builder()
  if (isTRUE(spec$ina)) bad_args$geocoords <- geocoords
  bad_args$InfoRetentionProb <- c(0.5, 0.5)
  bad_args$ModelName <- paste0(label, "_invalid_")
  bad_args$OutputDir <- out_dir
  expect_error(do.call(modified_fun, bad_args),
               "InfoRetentionProb must be a single value, vector of length nodes, or matrix nodes x Ntimesteps",
               paste0(label, " invalid retention length"))

  mark_progress(paste("PASS", label))
  cat("PASS", label, "\n")
  cat("  no decay:      ", paste(round(p_no, 1), collapse = " -> "), "\n")
  cat("  scalar 0.5:    ", paste(round(p_scalar, 1), collapse = " -> "), "\n")
  cat("  node vector:   ", paste(round(p_node, 1), collapse = " -> "), "\n")
  cat("  node x time:   ", paste(round(p_time, 1), collapse = " -> "), "\n")
  cat("  detect refresh:", paste(round(p_refresh, 1), collapse = " -> "), "\n")

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    Function = label,
    NoDecay_t1 = p_no[1], NoDecay_t2 = p_no[2], NoDecay_t3 = p_no[3],
    Scalar05_t1 = p_scalar[1], Scalar05_t2 = p_scalar[2], Scalar05_t3 = p_scalar[3],
    NodeVector_t1 = p_node[1], NodeVector_t2 = p_node[2], NodeVector_t3 = p_node[3],
    NodeTime_t1 = p_time[1], NodeTime_t2 = p_time[2], NodeTime_t3 = p_time[3],
    Refresh_t1 = p_refresh[1], Refresh_t2 = p_refresh[2], Refresh_t3 = p_refresh[3],
    BackwardCompatible = TRUE,
    InvalidLengthRejected = TRUE,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(bundle_root, "info_decay_results.csv"), row.names = FALSE)
writeLines(capture.output(print(summary_df, row.names = FALSE)), file.path(bundle_root, "info_decay_results.txt"))
cat("\nALL INFORMATION-DECAY MODEL TESTS PASSED\n")
