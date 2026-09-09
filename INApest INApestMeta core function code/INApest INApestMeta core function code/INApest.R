###############################################################################
### INApest -- binary node invasion simulation engine
###
### INApest represents each spatial node as invaded or not invaded. It is used
### when presence/absence is sufficient rather than abundance, land-use or
### life-stage structure.
###
### Each stochastic simulation links three parts of the response system:
###   Biology     - invasion survives, establishes and spreads among nodes.
###   Information - surveillance and information sharing determine what is known.
###   Response    - informed nodes may be managed to reduce spread or eradicate.
###
### These states are updated repeatedly through space and time. An optional
### binary pathogen process can be coupled to host occupancy. Repeated runs show
### the range and frequency of possible invasion and response outcomes.
###############################################################################

###############################################################################
### INApest -- binary node model with optional binary pathogen occupancy
### Definitive pathogen-enabled refactor, 2026-08-25.
###
### Host/pest state remains binary (`Invaded`). When Pathogen is supplied it
### must be INApestPathogen(Model = "Binary"); pathogen occupancy is tracked in
### a separate binary state constrained by PathogenPresent <= Invaded.
###############################################################################

# Keep probability values within the valid range from 0 to 1.
# Used after draws or calculations that can fall outside this range.
.inapest_clip01 <- function(x) pmin(1, pmax(0, x))

# Resolve a flexible model input to one value per node for the current timestep.
# Inputs can be constant, node-specific, time-varying, or supplied by a function.
.inapest_resolve <- function(x, timestep, n, Ntimesteps, name, allow_na = FALSE) {
  if (is.function(x)) {
    z <- x(timestep = timestep, Ntimesteps = Ntimesteps)
  } else if (is.matrix(x)) {
    if (!all(dim(x) == c(n, Ntimesteps))) stop(name, " matrix must have dimensions nodes x Ntimesteps")
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

# Draw node-specific probabilities around the supplied mean and variation.
# Used to represent stochastic variation in surveillance and response processes.
.inapest_prob_draw <- function(mu, sd, n) {
  mu <- rep_len(as.numeric(mu), n)
  sd <- rep_len(as.numeric(sd), n)
  .inapest_clip01(stats::rnorm(n, mean = mu, sd = pmax(0, sd)))
}

# Return the spatial connectivity used in the current timestep.
# Allows short- and long-distance dispersal pathways to be fixed or time-varying.
.inapest_connectivity <- function(x, timestep, n, Ntimesteps, name) {
  if (is.matrix(x)) {
    if (!all(dim(x) == c(n, n))) stop(name, " matrix must be nodes x nodes")
    return(x)
  }
  if (is.array(x) && length(dim(x)) == 3L) {
    if (!all(dim(x) == c(n, n, Ntimesteps))) stop(name, " 3D array must be nodes x nodes x Ntimesteps")
    return(x[, , timestep])
  }
  if (length(x) == 1L && is.numeric(x)) return(matrix(as.numeric(x), n, n))
  stop(name, " must be scalar, nodes x nodes matrix, or nodes x nodes x Ntimesteps array")
}

INApest <- function(
  ModelName,                     # Model and output name
  Nperm,                         # Number of stochastic simulation runs
  Ntimesteps,                    # Timesteps in each simulation
  DetectionProb,                 # Background-surveillance detection probability
  DetectionSD = NULL,            # Variation in background detection probability
  ManageProb,                    # Management probability when information is available
  ManageSD = NULL,               # Variation in management probability
  EradicationProb,               # Probability management eradicates an invaded node
  EradicationSD = NULL,          # Variation in eradication probability
  SpreadReduction,               # Reduction in spread from a managed node
  SpreadReductionSD = NULL,      # Variation in spread reduction
  InitialInvasion = NA,          # Starting invaded/not-invaded state by node
  InitBioP = NA,                 # Proportion of nodes initially invaded
  InvasionRisk = NA,             # Node-specific risk or weighting for invasion
  InitialInfo = NA,              # Starting information state by node
  InitInfoP = NA,                # Proportion of nodes initially with information
  ExternalInfoProb = NA,        # Probability or weighting for external information
  InfoRetentionProb = 1,        # Probability information is retained between timesteps
  InfoPersistenceSteps = NA,   # Timesteps information persists after local evidence
  EnvEstabProb = 1,             # Establishment probability at receiving nodes
  Survival = 1,                 # Probability an invaded node survives a timestep
  SDDprob,                       # Short-distance source-to-target dispersal
  SEAM = 0,                      # Information transfer between source and target nodes
  LDDprob = 0,                   # Long-distance source-to-target dispersal
  OngoingExternalInvasion = FALSE, # Allow new external invasions after initialisation
  OngoingExternalInfo = FALSE,     # Allow new external information after initialisation
  OutputDir = NA,                # Directory for saved outputs
  DoPlots = TRUE,                # Legacy plotting option; plotting is post-processing
  Pathogen = NULL,               # Optional binary pathogen process
  SaveResults = TRUE,           # Save simulation outputs to disk
  Seed = NULL,                  # Random seed for reproducible simulations
  InformationAcquisition = NULL, # Local evidence used: host, pathogen or both
  InfoTriggeredDetectionProb = 0, # Detection probability where information exists
  InfoTriggeredDetectionSD = NULL # Variation in information-triggered detection
) {
  # ---------------------------------------------------------------------------
  # Set up and validate the simulation
  # ---------------------------------------------------------------------------

  # Set a reproducible random sequence when requested.
  if (!is.null(Seed)) set.seed(Seed)

  # Infer the number of spatial nodes from short-distance connectivity.
  n <- if (is.matrix(SDDprob)) nrow(SDDprob) else if (is.array(SDDprob)) dim(SDDprob)[1] else stop("SDDprob must be a matrix or 3D array")
  # Check run dimensions and spatial connectivity inputs.
  if (Nperm < 1L || Nperm != floor(Nperm)) stop("Nperm must be a positive integer")
  if (Ntimesteps < 1L || Ntimesteps != floor(Ntimesteps)) stop("Ntimesteps must be a positive integer")
  if (is.array(SDDprob) && length(dim(SDDprob)) == 3L && !all(dim(SDDprob) == c(n, n, Ntimesteps))) stop("SDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")
  if (is.matrix(SDDprob) && ncol(SDDprob) != n) stop("SDDprob matrix must be square")
  if (!(length(LDDprob) == 1L || is.matrix(LDDprob) || (is.array(LDDprob) && length(dim(LDDprob)) == 3L))) stop("Invalid LDDprob")

  # Calculate a common mean for equivalent parameter forms.
  # Used only to supply default stochastic variation when no SD is provided.
  # Preserve the historical default-SD rule for resolver functions by taking
  # the mean over the same nodes x timesteps values represented by an equivalent
  # matrix parameterisation. This keeps scalar/vector/matrix/function forms
  # semantically equivalent when the corresponding SD argument is omitted.
  .inapest_default_mean <- function(x, name) {
    if (is.function(x)) {
      vals <- unlist(lapply(seq_len(Ntimesteps), function(tt)
        .inapest_resolve(x, tt, n, Ntimesteps, name)), use.names = FALSE)
      return(mean(vals, na.rm = TRUE))
    }
    mean(x, na.rm = TRUE)
  }

  # Supply default stochastic variation and identify whether targeted
  # information-triggered surveillance is active.
  if (is.null(DetectionSD)) DetectionSD <- .inapest_default_mean(DetectionProb, "DetectionProb") / 10
  if (is.null(InfoTriggeredDetectionSD)) InfoTriggeredDetectionSD <- .inapest_default_mean(InfoTriggeredDetectionProb, "InfoTriggeredDetectionProb") / 10
  UseInfoTriggeredSurveillance <- is.function(InfoTriggeredDetectionProb) || is.function(InfoTriggeredDetectionSD) ||
    any(as.numeric(InfoTriggeredDetectionProb) != 0, na.rm = TRUE) ||
    any(as.numeric(InfoTriggeredDetectionSD) != 0, na.rm = TRUE)
  if (is.null(ManageSD)) ManageSD <- .inapest_default_mean(ManageProb, "ManageProb") / 10
  if (is.null(EradicationSD)) EradicationSD <- .inapest_default_mean(EradicationProb, "EradicationProb") / 10
  if (is.null(SpreadReductionSD)) SpreadReductionSD <- (1 - .inapest_default_mean(SpreadReduction, "SpreadReduction")) / 10

  # Validate programmed information persistence across all timesteps.
  UseInfoPersistence <- is.function(InfoPersistenceSteps) || any(!is.na(InfoPersistenceSteps))
  if (UseInfoPersistence) {
    for (tt in seq_len(Ntimesteps)) {
      tmp <- .inapest_resolve(InfoPersistenceSteps, tt, n, Ntimesteps, "InfoPersistenceSteps", allow_na = TRUE)
      if (any(!is.na(tmp) & (!is.finite(tmp) | tmp < 0 | tmp != floor(tmp))))
        stop("InfoPersistenceSteps must be non-negative whole numbers or NA")
    }
  }
  # Validate the optional binary pathogen process against this simulation.
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

  # ---------------------------------------------------------------------------
  # Allocate arrays used to record simulation histories
  # ---------------------------------------------------------------------------

  # Store biological, information, surveillance, management and pathogen states
  # for every node, timestep and stochastic run.
  InvasionResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  DetectedResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  BackgroundDetectedResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  InfoTriggeredDetectedResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  BackgroundDetectionProbabilityResults <- array(0, dim = c(n, Ntimesteps, Nperm))
  InfoTriggeredDetectionProbabilityResults <- array(0, dim = c(n, Ntimesteps, Nperm))
  InformationStateBeforeSurveillanceResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  HaveInfoResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  ManagingResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  PathogenPresentResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  PathogenHostExtinctionResults <- array(0L, dim = c(n, Ntimesteps, Nperm))
  PathogenDetectedResults <- array(0L, dim = c(n, Ntimesteps, Nperm))

  # ---------------------------------------------------------------------------
  # Run independent stochastic simulations
  # ---------------------------------------------------------------------------

  for (perm in seq_len(Nperm)) {
    # Initialise the hidden biological state. Explicit starting occupancy takes
    # priority; otherwise initial invasion is sampled from the supplied settings.
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

    # Initialise what the response system knows, separately from true occupancy.
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

    # Apply initial background surveillance; host detections create information
    # only when host evidence is an accepted information pathway.
    mu_det <- .inapest_resolve(DetectionProb, 1L, n, Ntimesteps, "DetectionProb")
    sd_det <- .inapest_resolve(DetectionSD, 1L, n, Ntimesteps, "DetectionSD")
    NodeDetectionProb <- .inapest_prob_draw(mu_det, sd_det, n)
    InitDetection <- stats::rbinom(n, 1L, InitBio * NodeDetectionProb)
    if (HostTriggersInfo) InitInfo[InitInfo == 0L] <- InitDetection[InitInfo == 0L]

    # Set the working biological, information and optional pathogen states.
    Invaded <- InitBio
    HaveInfo <- InitInfo
    PathogenPresent <- if (is.null(Pathogen)) integer(n) else Pathogen$binary_initial(Invaded, n, Ntimesteps)
    InitPathogenDetection <- integer(n)
    if (!is.null(Pathogen) && PathogenTriggersInfo) {
      InitPathogenDetection <- Pathogen$binary_detect(PathogenPresent, 1L, Ntimesteps)
      HaveInfo[HaveInfo == 0L] <- InitPathogenDetection[HaveInfo == 0L]
    }
    # Initialise the evidence clock used for programmed stopping of information.
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

    # -------------------------------------------------------------------------
    # Advance the linked biology-information-response cycle through time
    # -------------------------------------------------------------------------

    for (timestep in seq_len(Ntimesteps)) {
      # Resolve current dispersal and environmental establishment to obtain the
      # source-to-target probability of successful spread.
      SDD <- .inapest_connectivity(SDDprob, timestep, n, Ntimesteps, "SDDprob")
      LDD <- .inapest_connectivity(LDDprob, timestep, n, Ntimesteps, "LDDprob")
      DispProb <- 1 - (1 - SDD) * (1 - LDD)
      env <- .inapest_resolve(EnvEstabProb, timestep, n, Ntimesteps, "EnvEstabProb")
      BPAM <- sweep(DispProb, 2L, .inapest_clip01(env), `*`)

      # Draw the node-level biological, surveillance and management probabilities
      # used in this timestep.
      NodeSurvival <- .inapest_clip01(.inapest_resolve(Survival, timestep, n, Ntimesteps, "Survival"))
      NodeDetectionProb <- .inapest_prob_draw(.inapest_resolve(DetectionProb, timestep, n, Ntimesteps, "DetectionProb"), .inapest_resolve(DetectionSD, timestep, n, Ntimesteps, "DetectionSD"), n)
      NodeManageProb <- .inapest_prob_draw(.inapest_resolve(ManageProb, timestep, n, Ntimesteps, "ManageProb"), .inapest_resolve(ManageSD, timestep, n, Ntimesteps, "ManageSD"), n)
      NodeSpreadReduction <- .inapest_prob_draw(.inapest_resolve(SpreadReduction, timestep, n, Ntimesteps, "SpreadReduction"), .inapest_resolve(SpreadReductionSD, timestep, n, Ntimesteps, "SpreadReductionSD"), n)
      NodeEradicationProb <- .inapest_prob_draw(.inapest_resolve(EradicationProb, timestep, n, Ntimesteps, "EradicationProb"), .inapest_resolve(EradicationSD, timestep, n, Ntimesteps, "EradicationSD"), n)

      # Use current information to activate management, then apply natural
      # survival and management-driven eradication before host spread.
      # Binary host-mode semantics treat an informed extant host as continuing
      # local host evidence. In pathogen-only mode, host presence does not
      # refresh the programmed-stop clock; only pathogen detection does.
      if (UseInfoPersistence && HostTriggersInfo)
        LastKnownPresence[HaveInfo == 1L & Invaded == 1L] <- timestep
      Managing <- stats::rbinom(n, 1L, NodeManageProb * HaveInfo)
      Invaded <- Invaded * stats::rbinom(n, 1L, NodeSurvival * (1 - NodeEradicationProb * Managing))
      PathogenPresent[Invaded == 0L] <- 0L
      Detected <- Invaded * HaveInfo

      # Spread surviving invasions through combined SDD and LDD pathways.
      # Management reduces outgoing spread and establishment depends on target habitat.
      # Host dispersal: source rows -> target columns.
      pspread <- BPAM * Invaded * (1 - Managing * NodeSpreadReduction)
      pspread <- .inapest_clip01(pspread)
      RandBPAM <- matrix(stats::rbinom(n * n, 1L, as.vector(pspread)), n, n)
      NewInvasion <- as.integer(colSums(RandBPAM) > 0)
      Invaded[Invaded == 0L] <- NewInvasion[Invaded == 0L]

      # Update information memory. Programmed stopping takes precedence where
      # specified; other informed nodes can lose information stochastically.
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

      # Transfer information from detected nodes through the socioeconomic
      # adjacency matrix (SEAM).
      if (is.matrix(SEAM)) {
        if (!all(dim(SEAM) == c(n, n))) stop("SEAM must be nodes x nodes")
        pinfo <- .inapest_clip01(SEAM * Detected)
        RandSEAM <- matrix(stats::rbinom(n * n, 1L, as.vector(pinfo)), n, n)
        InfoTransferred <- as.integer(colSums(RandSEAM) > 0)
        HaveInfo[HaveInfo == 0L] <- InfoTransferred[HaveInfo == 0L]
      }

      # Add new invasion or information arriving from outside the modelled system.
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

      # Advance the optional binary pathogen process after host invasion dynamics.
      # Pathogen detection can also create or refresh response information.
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

      # Record the response and biological states reached before host surveillance.
      ManagingResults[, timestep, perm] <- Managing
      InvasionResults[, timestep, perm] <- Invaded
      PathogenPresentResults[, timestep, perm] <- PathogenPresent
      PathogenHostExtinctionResults[, timestep, perm] <- PathogenHostExtinction
      PathogenDetectedResults[, timestep, perm] <- PathogenDetected

      # Run the host-surveillance observation round after biological and response
      # processes. Detections update information for later response.
      # Host-surveillance observation round. Information available before either
      # host-surveillance stream is frozen so a background detection cannot
      # activate targeted surveillance retrospectively in the same round.
      InfoBeforeSurveillance <- as.integer(HaveInfo != 0L)
      InformationStateBeforeSurveillanceResults[, timestep, perm] <- InfoBeforeSurveillance

      BackgroundDetectionProbabilityResults[, timestep, perm] <- NodeDetectionProb
      BackgroundDetection <- stats::rbinom(n, 1L, Invaded * NodeDetectionProb)
      InfoTriggeredDetection <- integer(n)
      if (UseInfoTriggeredSurveillance) {
        NodeInfoTriggeredDetectionProb <- .inapest_prob_draw(
          .inapest_resolve(InfoTriggeredDetectionProb, timestep, n, Ntimesteps, "InfoTriggeredDetectionProb"),
          .inapest_resolve(InfoTriggeredDetectionSD, timestep, n, Ntimesteps, "InfoTriggeredDetectionSD"), n
        )
        InfoTriggeredDetectionProbabilityResults[, timestep, perm] <- NodeInfoTriggeredDetectionProb
        InfoTriggeredDetection <- stats::rbinom(
          n, 1L, Invaded * NodeInfoTriggeredDetectionProb * InfoBeforeSurveillance
        )
      }
      BackgroundDetectedResults[, timestep, perm] <- BackgroundDetection
      InfoTriggeredDetectedResults[, timestep, perm] <- InfoTriggeredDetection

      # Combine the two host-surveillance pathways as local detection evidence.
      HostDetectionEvidence <- pmax(BackgroundDetection, InfoTriggeredDetection)
      if (HostTriggersInfo) {
        if (UseInfoPersistence) LastKnownPresence[HostDetectionEvidence == 1L] <- timestep
        HaveInfo[HaveInfo == 0L] <- HostDetectionEvidence[HaveInfo == 0L]
      }
      HaveInfoResults[, timestep, perm] <- HaveInfo
      DetectedResults[, timestep, perm] <- HaveInfo * Invaded
    }
  }

  # ---------------------------------------------------------------------------
  # Summarise and return simulation outputs
  # ---------------------------------------------------------------------------

  # Average binary occupancy across stochastic runs for each node and timestep.
  InvasionProb <- apply(InvasionResults, c(1L, 2L), mean)
  PathogenPresenceProb <- apply(PathogenPresentResults, c(1L, 2L), mean)

  # Return both realised histories and across-run probability summaries.
  out <- list(
    ModelName = ModelName,
    InvasionResults = InvasionResults,
    ManagingResults = ManagingResults,
    DetectedResults = DetectedResults,
    BackgroundDetectedResults = BackgroundDetectedResults,
    InfoTriggeredDetectedResults = InfoTriggeredDetectedResults,
    BackgroundDetectionProbabilityResults = BackgroundDetectionProbabilityResults,
    InfoTriggeredDetectionProbabilityResults = InfoTriggeredDetectionProbabilityResults,
    InformationStateBeforeSurveillanceResults = InformationStateBeforeSurveillanceResults,
    HaveInfoResults = HaveInfoResults,
    InvasionProb = InvasionProb,
    PathogenPresentResults = PathogenPresentResults,
    PathogenHostExtinctionResults = PathogenHostExtinctionResults,
    PathogenDetectedResults = PathogenDetectedResults,
    PathogenPresenceProb = PathogenPresenceProb
  )
  class(out) <- c("INApest", "list")

  # Save the standard simulation histories and summaries when requested.
  if (SaveResults) {
    if (is.na(OutputDir)) OutputDir <- ""
    if (nzchar(OutputDir) && !dir.exists(OutputDir)) dir.create(OutputDir, recursive = TRUE)
    stem <- file.path(OutputDir, ModelName)
    saveRDS(ManagingResults, paste0(stem, "InfoLargeOut.rds"))
    saveRDS(InvasionResults, paste0(stem, "InvasionLargeOut.rds"))
    saveRDS(DetectedResults, paste0(stem, "DetectedLargeOut.rds"))
    saveRDS(BackgroundDetectedResults, paste0(stem, "BackgroundDetectedLargeOut.rds"))
    saveRDS(InfoTriggeredDetectedResults, paste0(stem, "InfoTriggeredDetectedLargeOut.rds"))
    saveRDS(BackgroundDetectionProbabilityResults, paste0(stem, "BackgroundDetectionProbabilityLargeOut.rds"))
    saveRDS(InfoTriggeredDetectionProbabilityResults, paste0(stem, "InfoTriggeredDetectionProbabilityLargeOut.rds"))
    saveRDS(InformationStateBeforeSurveillanceResults, paste0(stem, "InformationStateBeforeSurveillanceLargeOut.rds"))
    saveRDS(HaveInfoResults, paste0(stem, "HaveInfoLargeOut.rds"))
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
  # Return the complete result object without printing it automatically.
  invisible(out)
}
