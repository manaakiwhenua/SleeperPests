###############################################################################
### INApestVertebratePointParallel -- parallel vertebrate point wrapper
###
### INApestVertebratePoint() is the vertebrate-specialist engine implemented in
### INApestPointTransitionMatrix.R. This file provides its explicit public
### parallel facade without duplicating point, stage-transition or vertebrate
### biology. Independent permutations are run through the same serial engine and
### the point histories, control observations, contacts, interactions and summary
### tables are recombined with global permutation identifiers.
###############################################################################

# Create reproducible L'Ecuyer-CMRG streams for point workers.
.ivpp_make_streams <- function(Nperm, Seed = NULL) {
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

# Row-bind a named data-frame field across worker-local results.
.ivpp_rbind <- function(xs, name) {
  parts <- lapply(xs, function(x) x[[name]])
  parts <- parts[vapply(parts, is.data.frame, logical(1))]
  if (!length(parts)) return(data.frame())
  nonempty <- parts[vapply(parts, ncol, integer(1)) > 0L]
  if (!length(nonempty)) return(data.frame())
  all_names <- unique(unlist(lapply(nonempty, names), use.names = FALSE))
  nonempty <- lapply(nonempty, function(x) {
    miss <- setdiff(all_names, names(x))
    for (nm in miss) x[[nm]] <- rep(NA, nrow(x))
    x[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, nonempty)
  rownames(out) <- NULL
  out
}

# Translate worker-local permutation id 1 back to its global permutation id.
.ivpp_relabel <- function(x, perm) {
  for (nm in c(
    "PointHistory", "EventLog", "FinalPoints", "InfoSites", "ContactHistory",
    "InteractionEvents", "ControlObservationHistory", "Summary",
    "BiocontrolHistory", "BiocontrolPointEvents"
  )) {
    if (is.data.frame(x[[nm]]) && "perm" %in% names(x[[nm]]))
      x[[nm]]$perm <- rep(perm, nrow(x[[nm]]))
  }
  x
}

INApestVertebratePointParallel <- function(
  ModelName = "INApestVertebratePointParallel", # Model and output name
  Nperm,                                         # Number of stochastic simulation runs
  Vertebrate = NULL,                             # Birth/HomeRange/Control/Interaction modules
  ...,                                           # Arguments forwarded to INApestVertebratePoint()
  Cores = max(1L, parallel::detectCores(logical = TRUE) - 1L), # Worker processes
  Backend = c("psock", "fork"),                # Parallel backend
  Seed = NULL,                                   # Reproducible stream seed
  Export = NULL,                                 # Additional worker objects by name
  OutputDir = NA,                                # Directory for optional combined RDS
  SaveResults = FALSE,                           # Save combined wrapper result
  DoProgress = TRUE,                             # Print wrapper progress
  InitialPointGenerator = NULL                   # Optional function(global perm) -> points
) {

  # ---------------------------------------------------------------------------
  # Validate the wrapper call and prepare one-permutation serial arguments.
  # ---------------------------------------------------------------------------
  if (!exists("INApestVertebratePoint", mode = "function", inherits = TRUE))
    stop("Source INApestPointTransitionMatrix.R first.")
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
  args$Vertebrate <- Vertebrate
  args$InitialPointGenerator <- InitialPointGenerator
  for (nm in c("Nperm", "ModelName", "Seed", "OutputDir", "SaveResults", "DoProgress"))
    args[[nm]] <- NULL

  streams <- .ivpp_make_streams(Nperm, Seed)
  serial_fun <- get("INApestVertebratePoint", mode = "function", inherits = TRUE)

  worker <- function(i, model_fun, base_args, rng_stream, model_name) {
    assign(".Random.seed", rng_stream[[i]], envir = .GlobalEnv)

    # InitialPointGenerator is evaluated by each one-permutation worker. Remap
    # its local perm = 1 argument to the original global permutation number.
    if (is.function(base_args$InitialPointGenerator)) {
      original_generator <- base_args$InitialPointGenerator
      global_perm <- i
      base_args$InitialPointGenerator <- local({
        f <- original_generator
        g <- global_perm
        function(perm) f(perm = g)
      })
    }

    call_args <- c(
      base_args,
      list(
        ModelName = paste0(model_name, "_perm", i),
        Nperm = 1L,
        SaveResults = FALSE,
        DoProgress = FALSE,
        Seed = NULL
      )
    )
    .ivpp_relabel(do.call(model_fun, call_args), i)
  }

  # ---------------------------------------------------------------------------
  # Execute independent point permutations using the selected backend.
  # ---------------------------------------------------------------------------
  if (DoProgress)
    message("Running ", Nperm, " vertebrate-point permutation",
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

    env <- environment(serial_fun)
    nms <- ls(env, all.names = TRUE)
    engine_symbols <- nms[grepl(
      "^(\\.ipp_|\\.ipptm_|\\.iv_|\\.ibp_|\\.INApestPointTransitionMatrix_engine$|INApestPoint|INApestSpatial|INApestHabitat|INApestVertebratePoint$|INApestPointBiocontrol|INApestBiocontrolPoint)",
      nms
    )]
    if (length(engine_symbols))
      parallel::clusterExport(cl, engine_symbols, envir = env)
    if (length(Export))
      parallel::clusterExport(cl, Export, envir = parent.frame())
    parallel::clusterExport(cl, ".ivpp_relabel", envir = environment())

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
  # Restore the serial vertebrate-point result contract with global perm ids.
  # ---------------------------------------------------------------------------
  out <- list(
    ModelName = ModelName,
    PointHistory = .ivpp_rbind(xs, "PointHistory"),
    EventLog = .ivpp_rbind(xs, "EventLog"),
    FinalPoints = .ivpp_rbind(xs, "FinalPoints"),
    InfoSites = .ivpp_rbind(xs, "InfoSites"),
    ContactHistory = .ivpp_rbind(xs, "ContactHistory"),
    InteractionEvents = .ivpp_rbind(xs, "InteractionEvents"),
    ControlObservationHistory = .ivpp_rbind(xs, "ControlObservationHistory"),
    Summary = .ivpp_rbind(xs, "Summary"),
    BiocontrolHistory = .ivpp_rbind(xs, "BiocontrolHistory"),
    BiocontrolPointEvents = .ivpp_rbind(xs, "BiocontrolPointEvents"),
    ParallelMeta = list(
      Nperm = Nperm, Cores = Cores, Backend = backend_used,
      Seed = Seed, elapsed_seconds = elapsed
    )
  )

  if (nrow(out$PointHistory))
    out$PointHistory <- out$PointHistory[order(out$PointHistory$perm, out$PointHistory$timestep, out$PointHistory$id), , drop = FALSE]
  if (nrow(out$FinalPoints))
    out$FinalPoints <- out$FinalPoints[order(out$FinalPoints$perm, out$FinalPoints$id), , drop = FALSE]
  if (nrow(out$Summary))
    out$Summary <- out$Summary[order(out$Summary$perm, out$Summary$timestep), , drop = FALSE]
  if (nrow(out$ControlObservationHistory))
    out$ControlObservationHistory <- out$ControlObservationHistory[
      order(out$ControlObservationHistory$perm, out$ControlObservationHistory$timestep,
            out$ControlObservationHistory$id), , drop = FALSE
    ]

  if (nrow(out$BiocontrolHistory))
    out$BiocontrolHistory <- out$BiocontrolHistory[
      order(out$BiocontrolHistory$perm, out$BiocontrolHistory$timestep,
            out$BiocontrolHistory$agent, out$BiocontrolHistory$node, out$BiocontrolHistory$stage),
      , drop = FALSE
    ]
  if (nrow(out$BiocontrolPointEvents) && all(c("perm", "timestep") %in% names(out$BiocontrolPointEvents)))
    out$BiocontrolPointEvents <- out$BiocontrolPointEvents[
      order(out$BiocontrolPointEvents$perm, out$BiocontrolPointEvents$timestep), , drop = FALSE
    ]

  class(out) <- c("INApestVertebratePointParallel", "INApestVertebratePoint",
                  "INApestPointTransitionMatrix", "list")

  if (SaveResults) {
    if (is.na(OutputDir)) OutputDir <- ""
    if (nzchar(OutputDir) && !dir.exists(OutputDir)) dir.create(OutputDir, recursive = TRUE)
    saveRDS(out, file.path(OutputDir, paste0(ModelName, "_VertebratePointParallelResults.rds")))
  }
  if (DoProgress) message("Parallel vertebrate-point simulation complete in ", round(elapsed, 2), " s.")
  out
}
