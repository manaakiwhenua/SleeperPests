###############################################################################
### INApestParallel -- parallel wrapper around definitive INApest()
###############################################################################

.inapest_parallel_streams <- function(Nperm, Seed = NULL) {
  old_kind <- RNGkind(); old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({do.call(RNGkind, as.list(old_kind)); if (old_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed", envir=.GlobalEnv)}, add=TRUE)
  RNGkind("L'Ecuyer-CMRG"); if (!is.null(Seed)) set.seed(Seed) else set.seed(sample.int(.Machine$integer.max,1L))
  z <- vector("list",Nperm); z[[1L]] <- .Random.seed
  if (Nperm > 1L) for (i in 2:Nperm) z[[i]] <- parallel::nextRNGStream(z[[i-1L]])
  z
}

INApestParallel <- function(
  ModelName,
  Nperm,
  ...,
  Cores = max(1L, parallel::detectCores(logical = TRUE) - 1L),
  Backend = c("psock", "fork"),
  Seed = NULL,
  Export = NULL,
  OutputDir = NA,
  SaveResults = TRUE,
  DoPlots = FALSE
) {
  if (!exists("INApest", mode="function", inherits=TRUE)) stop("Source INApest.R first")
  if (!requireNamespace("parallel", quietly=TRUE)) stop("The 'parallel' package is required")
  if (Nperm < 1L || Nperm != floor(Nperm)) stop("Nperm must be a positive integer")
  Nperm <- as.integer(Nperm); Cores <- min(as.integer(Cores), Nperm); Backend <- match.arg(Backend)
  if (.Platform$OS.type == "windows" && Backend == "fork") stop("Use Backend='psock' on Windows")
  args <- list(...); args$Nperm <- NULL; args$ModelName <- NULL; args$Seed <- NULL; args$SaveResults <- NULL; args$DoPlots <- NULL; args$OutputDir <- NULL
  streams <- .inapest_parallel_streams(Nperm, Seed)
  fun <- get("INApest", mode="function", inherits=TRUE)
  worker <- function(i, fun, args, streams, model) {
    assign(".Random.seed", streams[[i]], envir=.GlobalEnv)
    do.call(fun, c(args, list(ModelName=paste0(model,"_perm",i), Nperm=1L, SaveResults=FALSE, DoPlots=FALSE, Seed=NULL)))
  }
  if (Cores == 1L) {
    xs <- lapply(seq_len(Nperm), worker, fun=fun, args=args, streams=streams, model=ModelName); backend_used <- "serial"
  } else if (Backend == "fork") {
    xs <- parallel::mclapply(seq_len(Nperm), worker, fun=fun, args=args, streams=streams, model=ModelName, mc.cores=Cores, mc.set.seed=FALSE); backend_used <- "fork"
  } else {
    cl <- parallel::makeCluster(Cores, type="PSOCK"); on.exit(parallel::stopCluster(cl), add=TRUE)
    env <- environment(fun); symbols <- ls(env, all.names=TRUE); symbols <- symbols[grepl("^(\\.inapest_|INApest$)", symbols)]
    if (length(symbols)) parallel::clusterExport(cl, symbols, envir=env)
    if (length(Export)) parallel::clusterExport(cl, Export, envir=parent.frame())
    xs <- parallel::parLapply(cl, seq_len(Nperm), worker, fun=fun, args=args, streams=streams, model=ModelName)
    parallel::stopCluster(cl); on.exit(NULL, add=FALSE); backend_used <- "psock"
  }
  bind3 <- function(name) {
    aa <- lapply(xs, `[[`, name)
    # each serial worker returns nodes x timesteps x 1
    n <- dim(aa[[1]])[1]; tt <- dim(aa[[1]])[2]
    out <- array(0, dim=c(n,tt,Nperm))
    for (i in seq_len(Nperm)) out[,,i] <- aa[[i]][,,1]
    out
  }
  InvasionResults <- bind3("InvasionResults")
  ManagingResults <- bind3("ManagingResults")
  DetectedResults <- bind3("DetectedResults")
  PathogenPresentResults <- bind3("PathogenPresentResults")
  PathogenHostExtinctionResults <- bind3("PathogenHostExtinctionResults")
  PathogenDetectedResults <- bind3("PathogenDetectedResults")
  InvasionProb <- apply(InvasionResults,c(1,2),mean)
  PathogenPresenceProb <- apply(PathogenPresentResults,c(1,2),mean)
  out <- list(ModelName=ModelName, InvasionResults=InvasionResults, ManagingResults=ManagingResults,
              DetectedResults=DetectedResults, InvasionProb=InvasionProb,
              PathogenPresentResults=PathogenPresentResults,
              PathogenHostExtinctionResults=PathogenHostExtinctionResults,
              PathogenDetectedResults=PathogenDetectedResults,
              PathogenPresenceProb=PathogenPresenceProb,
              ParallelMeta=list(Nperm=Nperm,Cores=Cores,Backend=backend_used,Seed=Seed))
  class(out) <- c("INApestParallel","INApest","list")
  if (SaveResults) {
    if (is.na(OutputDir)) OutputDir <- ""; if(nzchar(OutputDir)&&!dir.exists(OutputDir))dir.create(OutputDir,recursive=TRUE)
    stem <- file.path(OutputDir,ModelName)
    saveRDS(ManagingResults,paste0(stem,"InfoLargeOut.rds")); saveRDS(InvasionResults,paste0(stem,"InvasionLargeOut.rds")); saveRDS(DetectedResults,paste0(stem,"DetectedLargeOut.rds")); saveRDS(InvasionProb,paste0(stem,"InvasionProb.rds"))
    if (any(PathogenPresentResults != 0) || any(PathogenHostExtinctionResults != 0)) {saveRDS(PathogenPresentResults,paste0(stem,"PathogenPresentLargeOut.rds"));saveRDS(PathogenHostExtinctionResults,paste0(stem,"PathogenHostExtinctionLargeOut.rds"));saveRDS(PathogenDetectedResults,paste0(stem,"PathogenDetectedLargeOut.rds"));saveRDS(PathogenPresenceProb,paste0(stem,"PathogenPresenceProb.rds"))}
  }
  invisible(out)
}
