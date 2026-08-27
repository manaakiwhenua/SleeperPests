###############################################################################
### INApest -- binary node model with optional binary pathogen occupancy
### Definitive pathogen-enabled refactor, 2026-08-25.
###
### Host/pest state remains binary (`Invaded`). When Pathogen is supplied it
### must be INApestPathogen(Model = "Binary"); pathogen occupancy is tracked in
### a separate binary state constrained by PathogenPresent <= Invaded.
###############################################################################

.inapest_clip01 <- function(x) pmin(1, pmax(0, x))

.inapest_resolve <- function(x, timestep, n, Ntimesteps, name, allow_na = FALSE) {
  if (is.function(x)) {
    z <- x(timestep = timestep, Ntimesteps = Ntimesteps)
  } else if (is.matrix(x)) {
    if (!identical(dim(x), c(n, Ntimesteps))) stop(name, " matrix must have dimensions nodes x Ntimesteps")
    z <- x[, timestep]
  } else if (length(x) == 1L) {
    z <- rep(x, n)
  } else if (length(x) == n) {
    z <- x
  } else {
    stop(name, " must be scalar, length nodes, nodes x Ntimesteps matrix, or resolver function")
  }
  z <- as.numeric(z)
  if (length(z) == 1L) z <- rep(z, n)
  if (length(z) != n) stop(name, " resolved to the wrong length")
  if (!allow_na && any(is.na(z))) stop(name, " may not contain NA")
  z
}

.inapest_prob_draw <- function(mu, sd, n) {
  mu <- rep_len(as.numeric(mu), n)
  sd <- rep_len(as.numeric(sd), n)
  .inapest_clip01(stats::rnorm(n, mean = mu, sd = pmax(0, sd)))
}

.inapest_connectivity <- function(x, timestep, n, Ntimesteps, name) {
  if (is.matrix(x)) {
    if (!identical(dim(x), c(n, n))) stop(name, " matrix must be nodes x nodes")
    return(x)
  }
  if (is.array(x) && length(dim(x)) == 3L) {
    if (!identical(dim(x), c(n, n, Ntimesteps))) stop(name, " 3D array must be nodes x nodes x Ntimesteps")
    return(x[, , timestep])
  }
  if (length(x) == 1L && is.numeric(x)) return(matrix(as.numeric(x), n, n))
  stop(name, " must be scalar, nodes x nodes matrix, or nodes x nodes x Ntimesteps array")
}

INApest <- function(
  ModelName,
  Nperm,
  Ntimesteps,
  DetectionProb,
  DetectionSD = NULL,
  ManageProb,
  ManageSD = NULL,
  EradicationProb,
  EradicationSD = NULL,
  SpreadReduction,
  SpreadReductionSD = NULL,
  InitialInvasion = NA,
  InitBioP = NA,
  InvasionRisk = NA,
  InitialInfo = NA,
  InitInfoP = NA,
  ExternalInfoProb = NA,
  InfoRetentionProb = 1,
  InfoPersistenceSteps = NA,
  EnvEstabProb = 1,
  Survival = 1,
  SDDprob,
  SEAM = 0,
  LDDprob = 0,
  OngoingExternalInvasion = FALSE,
  OngoingExternalInfo = FALSE,
  OutputDir = NA,
  DoPlots = TRUE,
  Pathogen = NULL,
  SaveResults = TRUE,
  Seed = NULL,
  InformationAcquisition = NULL
) {
  if (!is.null(Seed)) set.seed(Seed)
  n <- if (is.matrix(SDDprob)) nrow(SDDprob) else if (is.array(SDDprob)) dim(SDDprob)[1] else stop("SDDprob must be a matrix or 3D array")
  if (Nperm < 1L || Nperm != floor(Nperm)) stop("Nperm must be a positive integer")
  if (Ntimesteps < 1L || Ntimesteps != floor(Ntimesteps)) stop("Ntimesteps must be a positive integer")
  if (is.array(SDDprob) && length(dim(SDDprob)) == 3L && !identical(dim(SDDprob), c(n, n, Ntimesteps))) stop("SDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")
  if (is.matrix(SDDprob) && ncol(SDDprob) != n) stop("SDDprob matrix must be square")
  if (!(length(LDDprob) == 1L || is.matrix(LDDprob) || (is.array(LDDprob) && length(dim(LDDprob)) == 3L))) stop("Invalid LDDprob")

  if (is.null(DetectionSD)) DetectionSD <- mean(DetectionProb, na.rm = TRUE) / 10
  if (is.null(ManageSD)) ManageSD <- mean(ManageProb, na.rm = TRUE) / 10
  if (is.null(EradicationSD)) EradicationSD <- mean(EradicationProb, na.rm = TRUE) / 10
  if (is.null(SpreadReductionSD)) SpreadReductionSD <- (1 - mean(SpreadReduction, na.rm = TRUE)) / 10

  UseInfoPersistence <- any(!is.na(InfoPersistenceSteps))
  if (UseInfoPersistence) {
    tmp <- .inapest_resolve(InfoPersistenceSteps, 1L, n, Ntimesteps, "InfoPersistenceSteps", allow_na = TRUE)
    if (any(!is.na(tmp) & (!is.finite(tmp) | tmp < 0 | tmp != floor(tmp)))) stop("InfoPersistenceSteps must be non-negative whole numbers or NA")
  }
  if (!is.null(Pathogen)) {
    if (!inherits(Pathogen, "INApestPathogen") || !identical(Pathogen$Model, "Binary")) stop("Pathogen must be NULL or INApestPathogen(Model = 'Binary')")
    Pathogen$binary_validate(n, Ntimesteps)
  }

  # InformationAcquisition controls which local biological evidence can create
  # or refresh the single HaveInfo state. NULL preserves the pre-extension
  # contract: host evidence always informs, and pathogen evidence also informs
  # only when Pathogen$DetectionTriggersInfo is TRUE.
  if (is.null(InformationAcquisition)) {
    InformationAcquisition <- if (!is.null(Pathogen) && isTRUE(Pathogen$DetectionTriggersInfo)) "both" else "host"
  } else {
    InformationAcquisition <- match.arg(as.character(InformationAcquisition)[1L], c("host", "pathogen", "both"))
  }
  if (InformationAcquisition %in% c("pathogen", "both") && is.null(Pathogen))
    stop("InformationAcquisition = '", InformationAcquisition, "' requires a Binary Pathogen specification")
  HostTriggersInfo <- InformationAcquisition %in% c("host", "both")
  PathogenTriggersInfo <- InformationAcquisition %in% c("pathogen", "both")

  InvasionResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  DetectedResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  ManagingResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  PathogenPresentResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  PathogenHostExtinctionResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  PathogenDetectedResults <- array(0L, dim = c(n, Ntimesteps, Nperm))

  for (perm in seq_len(Nperm)) {
    InitBio <- integer(n)
    if (length(InitialInvasion) == n && all(!is.na(InitialInvasion))) {
      InitBio <- as.integer(InitialInvasion != 0)
    } else {
      risk <- NULL
      if (is.matrix(InvasionRisk) && nrow(InvasionRisk) == n) risk <- InvasionRisk[, 1L]
      else if (!is.matrix(InvasionRisk) && length(InvasionRisk) == n && !all(is.na(InvasionRisk))) risk <- InvasionRisk
      if (!is.na(InitBioP)[1]) {
        k <- min(n, ceiling(n * InitBioP))
        if (k > 0L) InitBio[sample(seq_len(n), size = k, prob = risk)] <- 1L
      } else if (!is.null(risk)) {
        InitBio <- stats::rbinom(n, 1L, .inapest_clip01(risk))
      }
    }

    InitInfo <- integer(n)
    if (length(InitialInfo) == n && !all(is.na(InitialInfo))) {
      InitInfo <- as.integer(InitialInfo != 0)
    } else if (!is.na(InitInfoP)[1] || (length(ExternalInfoProb) > 1L && !all(is.na(ExternalInfoProb)))) {
      risk_info <- NULL
      if (is.matrix(ExternalInfoProb) && nrow(ExternalInfoProb) == n) risk_info <- ExternalInfoProb[, 1L]
      else if (!is.matrix(ExternalInfoProb) && length(ExternalInfoProb) == n && !all(is.na(ExternalInfoProb))) risk_info <- ExternalInfoProb
      if (!is.na(InitInfoP)[1]) {
        k <- min(n, ceiling(n * InitInfoP))
        if (k > 0L) InitInfo[sample(seq_len(n), size = k, prob = risk_info)] <- 1L
      } else if (!is.null(risk_info)) InitInfo <- stats::rbinom(n, 1L, .inapest_clip01(risk_info))
    }

    mu_det <- .inapest_resolve(DetectionProb, 1L, n, Ntimesteps, "DetectionProb")
    sd_det <- .inapest_resolve(DetectionSD, 1L, n, Ntimesteps, "DetectionSD")
    NodeDetectionProb <- .inapest_prob_draw(mu_det, sd_det, n)
    InitDetection <- stats::rbinom(n, 1L, InitBio * NodeDetectionProb)
    if (HostTriggersInfo) InitInfo[InitInfo == 0L] <- InitDetection[InitInfo == 0L]

    Invaded <- InitBio
    HaveInfo <- InitInfo
    PathogenPresent <- if (is.null(Pathogen)) integer(n) else Pathogen$binary_initial(Invaded, n, Ntimesteps)
    InitPathogenDetection <- integer(n)
    if (!is.null(Pathogen) && PathogenTriggersInfo) {
      InitPathogenDetection <- Pathogen$binary_detect(PathogenPresent, 1L, Ntimesteps)
      HaveInfo[HaveInfo == 0L] <- InitPathogenDetection[HaveInfo == 0L]
    }
    LastKnownPresence <- rep(NA_real_, n)
    if (UseInfoPersistence) {
      # Explicit InitialInfo is accepted as current local evidence when the host
      # is present. Host detections count only in host/both mode; pathogen
      # detections count only in pathogen/both mode.
      explicit_initial <- as.integer(InitInfo != 0L)
      initial_evidence <- explicit_initial
      if (HostTriggersInfo) initial_evidence <- pmax(initial_evidence, InitDetection)
      if (PathogenTriggersInfo) initial_evidence <- pmax(initial_evidence, InitPathogenDetection)
      LastKnownPresence[initial_evidence == 1L & Invaded == 1L] <- 0
    }

    for (timestep in seq_len(Ntimesteps)) {
      SDD <- .inapest_connectivity(SDDprob, timestep, n, Ntimesteps, "SDDprob")
      LDD <- .inapest_connectivity(LDDprob, timestep, n, Ntimesteps, "LDDprob")
      DispProb <- 1 - (1 - SDD) * (1 - LDD)
      env <- .inapest_resolve(EnvEstabProb, timestep, n, Ntimesteps, "EnvEstabProb")
      BPAM <- sweep(DispProb, 2L, .inapest_clip01(env), `*`)

      NodeSurvival <- .inapest_clip01(.inapest_resolve(Survival, timestep, n, Ntimesteps, "Survival"))
      NodeDetectionProb <- .inapest_prob_draw(.inapest_resolve(DetectionProb, timestep, n, Ntimesteps, "DetectionProb"), .inapest_resolve(DetectionSD, timestep, n, Ntimesteps, "DetectionSD"), n)
      NodeManageProb <- .inapest_prob_draw(.inapest_resolve(ManageProb, timestep, n, Ntimesteps, "ManageProb"), .inapest_resolve(ManageSD, timestep, n, Ntimesteps, "ManageSD"), n)
      NodeSpreadReduction <- .inapest_prob_draw(.inapest_resolve(SpreadReduction, timestep, n, Ntimesteps, "SpreadReduction"), .inapest_resolve(SpreadReductionSD, timestep, n, Ntimesteps, "SpreadReductionSD"), n)
      NodeEradicationProb <- .inapest_prob_draw(.inapest_resolve(EradicationProb, timestep, n, Ntimesteps, "EradicationProb"), .inapest_resolve(EradicationSD, timestep, n, Ntimesteps, "EradicationSD"), n)

      # Binary host-mode semantics treat an informed extant host as continuing
      # local host evidence. In pathogen-only mode, host presence does not
      # refresh the programmed-stop clock; only pathogen detection does.
      if (UseInfoPersistence && HostTriggersInfo)
        LastKnownPresence[HaveInfo == 1L & Invaded == 1L] <- timestep
      Managing <- stats::rbinom(n, 1L, NodeManageProb * HaveInfo)
      Invaded <- Invaded * stats::rbinom(n, 1L, NodeSurvival * (1 - NodeEradicationProb * Managing))
      PathogenPresent[Invaded == 0L] <- 0L
      Detected <- Invaded * HaveInfo

      # Host dispersal: source rows -> target columns.
      pspread <- BPAM * Invaded * (1 - Managing * NodeSpreadReduction)
      pspread <- .inapest_clip01(pspread)
      RandBPAM <- matrix(stats::rbinom(n * n, 1L, as.vector(pspread)), n, n)
      NewInvasion <- as.integer(colSums(RandBPAM) > 0)
      Invaded[Invaded == 0L] <- NewInvasion[Invaded == 0L]

      NodeInfoPersistenceSteps <- .inapest_resolve(InfoPersistenceSteps, timestep, n, Ntimesteps, "InfoPersistenceSteps", allow_na = TRUE)
      programmed <- which(HaveInfo == 1L & !is.na(NodeInfoPersistenceSteps))
      if (length(programmed)) {
        age <- timestep - LastKnownPresence
        stop_nodes <- programmed[is.na(LastKnownPresence[programmed]) | age[programmed] >= NodeInfoPersistenceSteps[programmed]]
        if (length(stop_nodes)) HaveInfo[stop_nodes] <- 0L
      }
      NodeInfoRetentionProb <- .inapest_clip01(.inapest_resolve(InfoRetentionProb, timestep, n, Ntimesteps, "InfoRetentionProb"))
      decay <- which(HaveInfo == 1L & is.na(NodeInfoPersistenceSteps) & NodeInfoRetentionProb < 1)
      if (length(decay)) HaveInfo[decay] <- stats::rbinom(length(decay), 1L, NodeInfoRetentionProb[decay])

      if (is.matrix(SEAM)) {
        if (!identical(dim(SEAM), c(n, n))) stop("SEAM must be nodes x nodes")
        pinfo <- .inapest_clip01(SEAM * Detected)
        RandSEAM <- matrix(stats::rbinom(n * n, 1L, as.vector(pinfo)), n, n)
        InfoTransferred <- as.integer(colSums(RandSEAM) > 0)
        HaveInfo[HaveInfo == 0L] <- InfoTransferred[HaveInfo == 0L]
      }

      if (OngoingExternalInvasion) {
        risk <- .inapest_resolve(InvasionRisk, timestep, n, Ntimesteps, "InvasionRisk")
        ext <- stats::rbinom(n, 1L, .inapest_clip01(risk))
        Invaded[Invaded == 0L] <- ext[Invaded == 0L]
      }
      if (OngoingExternalInfo) {
        ep <- .inapest_resolve(ExternalInfoProb, timestep, n, Ntimesteps, "ExternalInfoProb")
        extinfo <- stats::rbinom(n, 1L, .inapest_clip01(ep))
        HaveInfo[HaveInfo == 0L] <- extinfo[HaveInfo == 0L]
      }

      PathogenHostExtinction <- integer(n)
      if (!is.null(Pathogen)) {
        ps <- Pathogen$binary_step(PathogenPresent, Invaded, timestep, Ntimesteps)
        PathogenPresent <- ps$PathogenPresent
        Invaded <- ps$Invaded
        PathogenHostExtinction <- ps$HostExtinction
      } else PathogenPresent[] <- 0L
      PathogenPresent[Invaded == 0L] <- 0L
      PathogenDetected <- integer(n)
      if (!is.null(Pathogen)) {
        PathogenDetected <- Pathogen$binary_detect(PathogenPresent, timestep, Ntimesteps)
        if (PathogenTriggersInfo) {
          if (UseInfoPersistence) LastKnownPresence[PathogenDetected == 1L] <- timestep
          HaveInfo[HaveInfo == 0L] <- PathogenDetected[HaveInfo == 0L]
        }
      }

      ManagingResults[, timestep, perm] <- Managing
      InvasionResults[, timestep, perm] <- Invaded
      PathogenPresentResults[, timestep, perm] <- PathogenPresent
      PathogenHostExtinctionResults[, timestep, perm] <- PathogenHostExtinction
      PathogenDetectedResults[, timestep, perm] <- PathogenDetected

      NewHostDetection <- stats::rbinom(n, 1L, Invaded * NodeDetectionProb)
      if (HostTriggersInfo) {
        if (UseInfoPersistence) LastKnownPresence[NewHostDetection == 1L] <- timestep
        HaveInfo[HaveInfo == 0L] <- NewHostDetection[HaveInfo == 0L]
      }
      DetectedResults[, timestep, perm] <- HaveInfo * Invaded
    }
  }

  InvasionProb <- apply(InvasionResults, c(1L, 2L), mean)
  PathogenPresenceProb <- apply(PathogenPresentResults, c(1L, 2L), mean)

  out <- list(
    ModelName = ModelName,
    InvasionResults = InvasionResults,
    ManagingResults = ManagingResults,
    DetectedResults = DetectedResults,
    InvasionProb = InvasionProb,
    PathogenPresentResults = PathogenPresentResults,
    PathogenHostExtinctionResults = PathogenHostExtinctionResults,
    PathogenDetectedResults = PathogenDetectedResults,
    PathogenPresenceProb = PathogenPresenceProb
  )
  class(out) <- c("INApest", "list")

  if (SaveResults) {
    if (is.na(OutputDir)) OutputDir <- ""
    if (nzchar(OutputDir) && !dir.exists(OutputDir)) dir.create(OutputDir, recursive = TRUE)
    stem <- file.path(OutputDir, ModelName)
    saveRDS(ManagingResults, paste0(stem, "InfoLargeOut.rds"))
    saveRDS(InvasionResults, paste0(stem, "InvasionLargeOut.rds"))
    saveRDS(DetectedResults, paste0(stem, "DetectedLargeOut.rds"))
    saveRDS(InvasionProb, paste0(stem, "InvasionProb.rds"))
    if (!is.null(Pathogen)) {
      saveRDS(PathogenPresentResults, paste0(stem, "PathogenPresentLargeOut.rds"))
      saveRDS(PathogenHostExtinctionResults, paste0(stem, "PathogenHostExtinctionLargeOut.rds"))
      saveRDS(PathogenDetectedResults, paste0(stem, "PathogenDetectedLargeOut.rds"))
      saveRDS(PathogenPresenceProb, paste0(stem, "PathogenPresenceProb.rds"))
    }
  }

  # Legacy plotting is intentionally left to post-hoc plotting code. DoPlots is
  # retained in the signature so existing calls remain valid.
  invisible(out)
}
