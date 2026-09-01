###############################################################################
### INApestPointTransitionMatrixParallel
### Definitive thin parallel wrapper around INApestPointTransitionMatrix().
### All population/pathogen biology lives in the serial engine.
###############################################################################

.ipptmp_make_streams <- function(Nperm, Seed = NULL) {
  old_kind <- RNGkind()
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (old_seed_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  if (!is.null(Seed)) set.seed(Seed) else set.seed(sample.int(.Machine$integer.max, 1L))
  streams <- vector("list", Nperm)
  streams[[1L]] <- .Random.seed
  if (Nperm > 1L) for (i in 2:Nperm) streams[[i]] <- parallel::nextRNGStream(streams[[i - 1L]])
  streams
}

.ipptmp_rbind <- function(xs, name) {
  parts <- lapply(xs, function(x) x[[name]])
  parts <- parts[vapply(parts, is.data.frame, logical(1))]
  if (!length(parts)) return(data.frame())
  nonempty <- parts[vapply(parts, ncol, integer(1)) > 0L]
  if (!length(nonempty)) return(data.frame())
  all_names <- unique(unlist(lapply(nonempty, names), use.names = FALSE))
  nonempty <- lapply(nonempty, function(x) {
    miss <- setdiff(all_names, names(x)); for (nm in miss) x[[nm]] <- rep(NA, nrow(x))
    x[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, nonempty); rownames(out) <- NULL; out
}

.ipptmp_relabel <- function(x, perm) {
  for (nm in c("PointHistory", "EventLog", "FinalPoints", "InfoSites", "ContactHistory", "InteractionEvents", "PathogenEvents", "Summary")) {
    if (is.data.frame(x[[nm]]) && "perm" %in% names(x[[nm]])) x[[nm]]$perm <- rep(perm, nrow(x[[nm]]))
  }
  x
}

INApestPointTransitionMatrixParallel <- function(
  ModelName = "INApestPointTransitionMatrixParallel",
  Nperm,
  Pathogen = NULL,
  ...,
  Cores = max(1L, parallel::detectCores(logical = TRUE) - 1L),
  Backend = c("psock", "fork"),
  Seed = NULL,
  Export = NULL,
  OutputDir = NA,
  SaveResults = FALSE,
  DoProgress = TRUE
) {
  if (!exists("INApestPointTransitionMatrix", mode = "function", inherits = TRUE))
    stop("Source INApestPointTransitionMatrix.R first.")
  if (!requireNamespace("parallel", quietly = TRUE)) stop("The 'parallel' package is required.")
  if (!is.numeric(Nperm) || length(Nperm) != 1L || Nperm < 1L || Nperm != floor(Nperm)) stop("Nperm must be a positive integer.")
  Nperm <- as.integer(Nperm)
  if (!is.numeric(Cores) || length(Cores) != 1L || Cores < 1L || Cores != floor(Cores)) stop("Cores must be a positive integer.")
  Cores <- min(as.integer(Cores), Nperm)
  Backend <- match.arg(Backend)
  if (.Platform$OS.type == "windows" && Backend == "fork") stop("Backend = 'fork' is not available on Windows; use 'psock'.")

  args <- list(...)
  args$Pathogen <- Pathogen
  args$Nperm <- NULL; args$ModelName <- NULL; args$Seed <- NULL
  args$OutputDir <- NULL; args$SaveResults <- NULL; args$DoProgress <- NULL
  streams <- .ipptmp_make_streams(Nperm, Seed)
  serial_fun <- get("INApestPointTransitionMatrix", mode = "function", inherits = TRUE)

  worker <- function(i, model_fun, base_args, rng_stream, model_name) {
    assign(".Random.seed", rng_stream[[i]], envir = .GlobalEnv)
    call_args <- c(base_args, list(ModelName = paste0(model_name, "_perm", i), Nperm = 1L,
                                   SaveResults = FALSE, DoProgress = FALSE, Seed = NULL))
    .ipptmp_relabel(do.call(model_fun, call_args), i)
  }

  if (DoProgress) message("Running ", Nperm, " permutations with ", Cores, " worker", if (Cores == 1L) "" else "s", ".")
  t0 <- proc.time()[[3L]]
  if (Cores == 1L) {
    xs <- lapply(seq_len(Nperm), worker, model_fun = serial_fun, base_args = args, rng_stream = streams, model_name = ModelName)
    backend_used <- "serial"
  } else if (Backend == "fork") {
    xs <- parallel::mclapply(seq_len(Nperm), worker, model_fun = serial_fun, base_args = args,
                             rng_stream = streams, model_name = ModelName,
                             mc.cores = Cores, mc.set.seed = FALSE)
    backend_used <- "fork"
  } else {
    cl <- parallel::makeCluster(Cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    env <- environment(serial_fun)
    nms <- ls(env, all.names = TRUE)
    engine_symbols <- nms[grepl("^(\\.ipp_|\\.ipptm_|\\.iv_|\\.INApestPointTransitionMatrix_engine$|INApestPoint|INApestSpatial|INApestHabitat|INApestVertebratePoint$)", nms)]
    if (length(engine_symbols)) parallel::clusterExport(cl, engine_symbols, envir = env)
    if (length(Export)) parallel::clusterExport(cl, Export, envir = parent.frame())
    parallel::clusterExport(cl, c(".ipptmp_relabel"), envir = environment())
    xs <- parallel::parLapply(cl, seq_len(Nperm), worker, model_fun = serial_fun, base_args = args,
                              rng_stream = streams, model_name = ModelName)
    parallel::stopCluster(cl); on.exit(NULL, add = FALSE)
    backend_used <- "psock"
  }
  elapsed <- proc.time()[[3L]] - t0

  out <- list(
    ModelName = ModelName,
    PointHistory = .ipptmp_rbind(xs, "PointHistory"),
    EventLog = .ipptmp_rbind(xs, "EventLog"),
    FinalPoints = .ipptmp_rbind(xs, "FinalPoints"),
    InfoSites = .ipptmp_rbind(xs, "InfoSites"),
    ContactHistory = .ipptmp_rbind(xs, "ContactHistory"),
    InteractionEvents = .ipptmp_rbind(xs, "InteractionEvents"),
    PathogenEvents = .ipptmp_rbind(xs, "PathogenEvents"),
    Summary = .ipptmp_rbind(xs, "Summary"),
    ParallelMeta = list(Nperm = Nperm, Cores = Cores, Backend = backend_used, Seed = Seed, elapsed_seconds = elapsed)
  )
  if (nrow(out$PointHistory)) out$PointHistory <- out$PointHistory[order(out$PointHistory$perm, out$PointHistory$timestep, out$PointHistory$id), , drop = FALSE]
  if (nrow(out$FinalPoints)) out$FinalPoints <- out$FinalPoints[order(out$FinalPoints$perm, out$FinalPoints$id), , drop = FALSE]
  if (nrow(out$Summary)) out$Summary <- out$Summary[order(out$Summary$perm, out$Summary$timestep), , drop = FALSE]
  class(out) <- c("INApestPointTransitionMatrixParallel", "list")

  if (SaveResults) {
    if (is.na(OutputDir)) OutputDir <- ""
    if (nzchar(OutputDir) && !dir.exists(OutputDir)) dir.create(OutputDir, recursive = TRUE)
    saveRDS(out, file.path(OutputDir, paste0(ModelName, "_PointTransitionMatrixParallelResults.rds")))
  }
  if (DoProgress) message("Parallel simulation complete in ", round(elapsed, 2), " s.")
  out
}
