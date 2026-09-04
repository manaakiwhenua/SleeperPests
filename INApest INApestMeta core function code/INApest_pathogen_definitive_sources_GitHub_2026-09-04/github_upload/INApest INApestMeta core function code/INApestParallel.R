###############################################################################
### INApestParallel -- parallel binary-node simulation wrapper
###
### Runs independent INApest permutations across serial, PSOCK or fork workers.
### It does not define a new biological model: every worker uses the same Binary
### occupancy, information, response and optional pathogen process as INApest().
###
### Reproducible random-number streams are assigned per permutation, then the
### node x timestep histories are combined into the standard INApest result.
###############################################################################

# Create one reproducible random-number stream per permutation.
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
  ModelName,                                    # Model and output name
  Nperm,                                        # Number of stochastic simulation runs
  ...,                                          # Additional arguments forwarded to the underlying engine
  Cores = max(1L, parallel::detectCores(logical = TRUE) - 1L), # Number of worker processes
  Backend = c("psock", "fork"),                 # Parallel backend: PSOCK or fork
  Seed = NULL,                                  # Random seed for reproducible simulations
  Export = NULL,                                # Additional objects exported to workers
  OutputDir = NA,                               # Directory for saved outputs
  SaveResults = TRUE,                           # Save standard simulation outputs to disk
  DoPlots = FALSE                               # Legacy plotting option; plotting is post-processing
) {

  # ---------------------------------------------------------------------------
  # Set up and validate the parallel run.
  # ---------------------------------------------------------------------------
  if (!exists("INApest", mode="function", inherits=TRUE)) stop("Source INApest.R first")
  if (!requireNamespace("parallel", quietly=TRUE)) stop("The 'parallel' package is required")
  if (Nperm < 1L || Nperm != floor(Nperm)) stop("Nperm must be a positive integer")
  Nperm <- as.integer(Nperm); Cores <- min(as.integer(Cores), Nperm); Backend <- match.arg(Backend)
  if (.Platform$OS.type == "windows" && Backend == "fork") stop("Use Backend='psock' on Windows")
  args <- list(...); args$Nperm <- NULL; args$ModelName <- NULL; args$Seed <- NULL; args$SaveResults <- NULL; args$DoPlots <- NULL; args$OutputDir <- NULL
  streams <- .inapest_parallel_streams(Nperm, Seed)
  fun <- get("INApest", mode="function", inherits=TRUE)
  worker <- function(i, model_fun, args, streams, model) {
    assign(".Random.seed", streams[[i]], envir=.GlobalEnv)
    do.call(model_fun, c(args, list(ModelName=paste0(model,"_perm",i), Nperm=1L, SaveResults=FALSE, DoPlots=FALSE, Seed=NULL)))
  }

  # ---------------------------------------------------------------------------
  # Run independent permutations using the selected backend.
  # ---------------------------------------------------------------------------
  if (Cores == 1L) {
    xs <- lapply(seq_len(Nperm), worker, model_fun=fun, args=args, streams=streams, model=ModelName); backend_used <- "serial"
  } else if (Backend == "fork") {
    xs <- parallel::mclapply(seq_len(Nperm), worker, model_fun=fun, args=args, streams=streams, model=ModelName, mc.cores=Cores, mc.set.seed=FALSE); backend_used <- "fork"
  } else {
    cl <- parallel::makeCluster(Cores, type="PSOCK"); on.exit(parallel::stopCluster(cl), add=TRUE)
    env <- environment(fun); symbols <- ls(env, all.names=TRUE); symbols <- symbols[grepl("^(\\.inapest_|INApest$)", symbols)]
    if (length(symbols)) parallel::clusterExport(cl, symbols, envir=env)
    if (length(Export)) parallel::clusterExport(cl, Export, envir=parent.frame())
    xs <- parallel::parLapply(cl, seq_len(Nperm), worker, model_fun=fun, args=args, streams=streams, model=ModelName)
    parallel::stopCluster(cl); on.exit(NULL, add=FALSE); backend_used <- "psock"
  }

  # ---------------------------------------------------------------------------
  # Combine worker histories into node x timestep x permutation arrays.
  # ---------------------------------------------------------------------------
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
  BackgroundDetectedResults <- bind3("BackgroundDetectedResults")
  InfoTriggeredDetectedResults <- bind3("InfoTriggeredDetectedResults")
  InformationStateBeforeSurveillanceResults <- bind3("InformationStateBeforeSurveillanceResults")
  HaveInfoResults <- bind3("HaveInfoResults")
  PathogenPresentResults <- bind3("PathogenPresentResults")
  PathogenHostExtinctionResults <- bind3("PathogenHostExtinctionResults")
  PathogenDetectedResults <- bind3("PathogenDetectedResults")
  BackgroundDetectionProbabilityResults <- bind3("BackgroundDetectionProbabilityResults")
  InfoTriggeredDetectionProbabilityResults <- bind3("InfoTriggeredDetectionProbabilityResults")
  InvasionProb <- apply(InvasionResults,c(1,2),mean)
  PathogenPresenceProb <- apply(PathogenPresentResults,c(1,2),mean)
  out <- list(ModelName=ModelName, InvasionResults=InvasionResults, ManagingResults=ManagingResults,
              DetectedResults=DetectedResults,
              BackgroundDetectedResults=BackgroundDetectedResults,
              InfoTriggeredDetectedResults=InfoTriggeredDetectedResults,
              InformationStateBeforeSurveillanceResults=InformationStateBeforeSurveillanceResults,
              HaveInfoResults=HaveInfoResults, InvasionProb=InvasionProb,
              PathogenPresentResults=PathogenPresentResults,
              PathogenHostExtinctionResults=PathogenHostExtinctionResults,
              PathogenDetectedResults=PathogenDetectedResults,
              BackgroundDetectionProbabilityResults=BackgroundDetectionProbabilityResults,
              InfoTriggeredDetectionProbabilityResults=InfoTriggeredDetectionProbabilityResults,
              PathogenPresenceProb=PathogenPresenceProb,
              ParallelMeta=list(Nperm=Nperm,Cores=Cores,Backend=backend_used,Seed=Seed))
  class(out) <- c("INApestParallel","INApest","list")

  # ---------------------------------------------------------------------------
  # Save the same standard histories produced by the serial engine.
  # ---------------------------------------------------------------------------
  if (SaveResults) {
    if (is.na(OutputDir)) OutputDir <- ""; if(nzchar(OutputDir)&&!dir.exists(OutputDir))dir.create(OutputDir,recursive=TRUE)
    stem <- file.path(OutputDir,ModelName)
    saveRDS(ManagingResults,paste0(stem,"InfoLargeOut.rds")); saveRDS(InvasionResults,paste0(stem,"InvasionLargeOut.rds")); saveRDS(DetectedResults,paste0(stem,"DetectedLargeOut.rds")); saveRDS(BackgroundDetectedResults,paste0(stem,"BackgroundDetectedLargeOut.rds")); saveRDS(InfoTriggeredDetectedResults,paste0(stem,"InfoTriggeredDetectedLargeOut.rds")); saveRDS(InformationStateBeforeSurveillanceResults,paste0(stem,"InformationStateBeforeSurveillanceLargeOut.rds")); saveRDS(HaveInfoResults,paste0(stem,"HaveInfoLargeOut.rds")); saveRDS(InvasionProb,paste0(stem,"InvasionProb.rds"))
    saveRDS(BackgroundDetectionProbabilityResults,paste0(stem,"BackgroundDetectionProbabilityLargeOut.rds")); saveRDS(InfoTriggeredDetectionProbabilityResults,paste0(stem,"InfoTriggeredDetectionProbabilityLargeOut.rds"))
    if (any(PathogenPresentResults != 0) || any(PathogenHostExtinctionResults != 0)) {saveRDS(PathogenPresentResults,paste0(stem,"PathogenPresentLargeOut.rds"));saveRDS(PathogenHostExtinctionResults,paste0(stem,"PathogenHostExtinctionLargeOut.rds"));saveRDS(PathogenDetectedResults,paste0(stem,"PathogenDetectedLargeOut.rds"));saveRDS(PathogenPresenceProb,paste0(stem,"PathogenPresenceProb.rds"))}
  }
  invisible(out)
}
