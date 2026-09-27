###############################################################################
### INApestVertebrateNodeParallel -- parallel vertebrate-node wrapper
###
### Runs independent INApestVertebrateNode permutations across serial, PSOCK or
### fork workers and combines the standard node x stage x timestep histories.
### The wrapper does not define new vertebrate biology: every worker calls the
### same INApestVertebrateNode() engine, including optional Birth, HomeRange,
### Control, Interaction and Biocontrol modules.
###
### Random-number streams are assigned by global permutation before execution,
### so Cores = 1 and multi-worker PSOCK/fork runs use the same per-permutation
### streams for a fixed Seed.
###############################################################################

# Create one reproducible L'Ecuyer-CMRG random-number stream per permutation.
.ivnp_make_streams <- function(Nperm, Seed = NULL) {
  old_kind <- RNGkind()
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (old_seed_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  if (!is.null(Seed)) set.seed(Seed) else set.seed(sample.int(.Machine$integer.max, 1L))
  streams <- vector("list", Nperm)
  streams[[1L]] <- .Random.seed
  if (Nperm > 1L) {
    for (i in 2:Nperm) streams[[i]] <- parallel::nextRNGStream(streams[[i - 1L]])
  }
  streams
}

# Bind one worker-local result field along its final permutation dimension.
# Every serial worker is called with Nperm = 1, so its final dimension is one.
.ivnp_bind_permutations <- function(worker_results, field) {
  parts <- lapply(worker_results, `[[`, field)
  if (!length(parts) || is.null(parts[[1L]])) return(NULL)

  first <- parts[[1L]]
  d <- dim(first)
  if (is.null(d) || length(d) < 2L || tail(d, 1L) != 1L)
    stop("Worker field '", field, "' must have a final permutation dimension of length 1.")

  dbase <- head(d, -1L)
  dn <- dimnames(first)
  dnbase <- if (is.null(dn)) NULL else head(dn, -1L)
  ans_dn <- if (is.null(dnbase)) NULL else c(dnbase, list(NULL))
  out <- array(vector(typeof(first), prod(c(dbase, length(parts)))),
               dim = c(dbase, length(parts)), dimnames = ans_dn)

  for (pp in seq_along(parts)) {
    z <- parts[[pp]]
    dz <- dim(z)
    if (is.null(dz) || !identical(head(dz, -1L), dbase) || tail(dz, 1L) != 1L)
      stop("Worker field '", field, "' has inconsistent dimensions across permutations.")
    z <- array(z, dim = dbase, dimnames = dnbase)
    idx <- c(rep(list(TRUE), length(dbase)), list(pp))
    out <- do.call(`[<-`, c(list(out), idx, list(value = z)))
  }
  out
}

INApestVertebrateNodeParallel <- function(
  ModelName = "INApestVertebrateNodeParallel", # Model and output name
  Nperm,                                        # Number of stochastic simulation runs
  ...,                                          # Arguments forwarded to INApestVertebrateNode()
  Cores = max(1L, parallel::detectCores(logical = TRUE) - 1L), # Worker processes
  Backend = c("psock", "fork"),               # Parallel backend
  Seed = NULL,                                  # Reproducible stream seed
  Export = NULL,                                # Additional worker objects by name
  OutputDir = NA,                               # Directory for optional combined RDS
  SaveResults = FALSE,                          # Save combined wrapper result
  DoProgress = TRUE                             # Print wrapper progress
) {

  # ---------------------------------------------------------------------------
  # Validate the wrapper call and prepare worker-local serial arguments.
  # ---------------------------------------------------------------------------
  if (!exists("INApestVertebrateNode", mode = "function", inherits = TRUE))
    stop("Source INApestVertebrateNode.R first.")
  if (!requireNamespace("parallel", quietly = TRUE))
    stop("The 'parallel' package is required.")
  if (!is.numeric(Nperm) || length(Nperm) != 1L || !is.finite(Nperm) ||
      Nperm < 1L || Nperm != floor(Nperm))
    stop("Nperm must be a positive integer.")
  Nperm <- as.integer(Nperm)
  if (!is.numeric(Cores) || length(Cores) != 1L || !is.finite(Cores) ||
      Cores < 1L || Cores != floor(Cores))
    stop("Cores must be a positive integer.")
  Cores <- min(as.integer(Cores), Nperm)
  Backend <- match.arg(Backend)
  if (.Platform$OS.type == "windows" && Backend == "fork")
    stop("Backend = 'fork' is not available on Windows; use 'psock'.")

  args <- list(...)
  # The wrapper owns run-level controls; workers are always one-permutation,
  # non-plotting and non-saving serial engine calls.
  for (nm in c("Nperm", "ModelName", "OutputDir", "SaveResults", "DoPlots", "DoProgress"))
    args[[nm]] <- NULL

  streams <- .ivnp_make_streams(Nperm, Seed)
  serial_fun <- get("INApestVertebrateNode", mode = "function", inherits = TRUE)

  worker <- function(i, model_fun, base_args, rng_stream, model_name) {
    assign(".Random.seed", rng_stream[[i]], envir = .GlobalEnv)
    call_args <- c(
      base_args,
      list(
        ModelName = paste0(model_name, "_perm", i),
        Nperm = 1L,
        OutputDir = tempdir(),
        SaveResults = FALSE,
        DoPlots = FALSE,
        DoProgress = FALSE
      )
    )
    do.call(model_fun, call_args)
  }

  # ---------------------------------------------------------------------------
  # Run independent vertebrate-node permutations using the selected backend.
  # ---------------------------------------------------------------------------
  if (DoProgress)
    message("Running ", Nperm, " vertebrate-node permutation",
            if (Nperm == 1L) "" else "s", " with ", Cores, " worker",
            if (Cores == 1L) "" else "s", ".")
  t0 <- proc.time()[[3L]]

  if (Cores == 1L) {
    xs <- lapply(seq_len(Nperm), worker, model_fun = serial_fun,
                 base_args = args, rng_stream = streams, model_name = ModelName)
    backend_used <- "serial"
  } else if (Backend == "fork") {
    xs <- parallel::mclapply(
      seq_len(Nperm), worker, model_fun = serial_fun,
      base_args = args, rng_stream = streams, model_name = ModelName,
      mc.cores = Cores, mc.set.seed = FALSE
    )
    backend_used <- "fork"
  } else {
    cl <- parallel::makeCluster(Cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export only the engine/helper families used by this architecture. User
    # module dependencies that live outside these families can be named in Export.
    env <- environment(serial_fun)
    nms <- ls(env, all.names = TRUE)
    engine_symbols <- nms[grepl(
      "^(\\.iv_|INApestVertebrateNode$|\\.inabc_|INApestBiocontrol)", nms
    )]
    if (length(engine_symbols))
      parallel::clusterExport(cl, engine_symbols, envir = env)
    if (length(Export))
      parallel::clusterExport(cl, Export, envir = parent.frame())

    xs <- parallel::parLapply(
      cl, seq_len(Nperm), worker, model_fun = serial_fun,
      base_args = args, rng_stream = streams, model_name = ModelName
    )
    parallel::stopCluster(cl)
    on.exit(NULL, add = FALSE)
    backend_used <- "psock"
  }
  elapsed <- proc.time()[[3L]] - t0

  # ---------------------------------------------------------------------------
  # Reconstruct the serial vertebrate-node result contract across permutations.
  # ---------------------------------------------------------------------------
  fields <- c(
    "PopulationResults", "PopulationStageResults", "InvasionResults",
    "DetectedResults", "ManagingResults", "BackgroundDetectedResults",
    "BackgroundSurveillanceDetectedResults", "InfoTriggeredDetectedResults",
    "InformationStateBeforeSurveillanceResults", "HaveInfoResults",
    "BackgroundDetectionProbabilityResults", "InfoTriggeredDetectionProbabilityResults",
    "RoutineControlObservationAbundanceResults", "RoutineControlDetectionProbabilityResults",
    "ControlDetectionResults", "ControlDeathResults", "ControlCostResults"
  )
  combined <- setNames(lapply(fields, function(nm) .ivnp_bind_permutations(xs, nm)), fields)

  InvasionProb <- apply(combined$InvasionResults, c(1L, 2L), mean)
  BiocontrolHistory <- if (exists("INApestBiocontrolBindWorkerHistories", mode = "function", inherits = TRUE))
    INApestBiocontrolBindWorkerHistories(xs) else NULL

  out <- c(list(
    ModelName = ModelName,
    # Legacy vertebrate aliases
    Population = combined$PopulationResults,
    PopulationStage = combined$PopulationStageResults,
    Invasion = combined$InvasionResults,
    Detection = combined$DetectedResults,
    Managing = combined$ManagingResults,
    ControlDeaths = combined$ControlDeathResults,
    ControlDetections = combined$ControlDetectionResults,
    ControlCost = combined$ControlCostResults,
    InvasionProbability = InvasionProb
  ), combined, list(
    BiocontrolHistory = BiocontrolHistory,
    ParallelMeta = list(
      Nperm = Nperm, Cores = Cores, Backend = backend_used,
      Seed = Seed, elapsed_seconds = elapsed
    )
  ))
  class(out) <- c("INApestVertebrateNodeParallel", "INApestVertebrateNode", "list")

  if (SaveResults) {
    if (is.na(OutputDir)) OutputDir <- ""
    if (nzchar(OutputDir) && !dir.exists(OutputDir)) dir.create(OutputDir, recursive = TRUE)
    saveRDS(out, file.path(OutputDir, paste0(ModelName, "_VertebrateNodeParallelResults.rds")))
  }
  if (DoProgress) message("Parallel vertebrate-node simulation complete in ", round(elapsed, 2), " s.")
  invisible(out)
}
