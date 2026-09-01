###############################################################################
### INApest proof-of-absence wrappers
###
### Common simulation-based Bayesian proof-of-absence inference for the
### maintained INApest model families.  Biological dynamics remain in the
### underlying INApest functions; this file adapts their outputs to one common
### evidence/state contract and conditions stochastic trajectories on observed
### surveillance histories.
###
### Version 2 redesigns the observation contract around two actual surveillance
### pathways: (1) Background surveillance and (2) InfoTriggered surveillance.
### HaveInfo / transferred / external information is response state, not evidence
### by itself.  Legacy Surveillance/Information names remain accepted as aliases.
### Observation likelihoods may be supplied by a separate observation model, so
### biological engines do not need to manufacture PoA-specific q arrays.
###############################################################################


###############################################################################
### Internal validation and array helpers
###############################################################################

.INApestPoAStop <- function(...)
{
  stop(..., call. = FALSE)
}


.INApestPoAClip01 <- function(x, Name = "value")
{
  if(length(x) != 1L || is.na(x) || !is.numeric(x) || !is.finite(x) || x < 0 || x > 1)
    .INApestPoAStop(Name, " must be a single finite numeric value between 0 and 1")
  as.numeric(x)
}


### Canonical evidence streams.  Old names are aliases only; "Information" is
### NOT interpreted as the HaveInfo state.
.INApestPoANormaliseEvidenceSources <- function(EvidenceSources)
{
  if(is.null(EvidenceSources) || !length(EvidenceSources))
    return(character(0))
  x <- as.character(EvidenceSources)
  map <- c(
    Background = "Background",
    InfoTriggered = "InfoTriggered",
    Surveillance = "Background",      # deprecated alias
    Information = "InfoTriggered"      # deprecated alias
  )
  bad <- setdiff(x, names(map))
  if(length(bad))
    .INApestPoAStop(
      "EvidenceSources must use Background and/or InfoTriggered. ",
      "Deprecated aliases Surveillance and Information are also accepted. Unknown: ",
      paste(bad, collapse=", ")
    )
  if(any(x %in% c("Surveillance","Information")))
    warning(
      "EvidenceSources='Surveillance'/'Information' is deprecated. Use 'Background'/'InfoTriggered'. ",
      "HaveInfo itself is not observational evidence.", call.=FALSE
    )
  unique(unname(map[x]))
}


.INApestPoAFirstNonNull <- function(...)
{
  z <- list(...)
  for(x in z) if(!is.null(x)) return(x)
  NULL
}


### Resolve per-particle no-detection probabilities from a separate observation
### model.  This mirrors the pathogen PoF architecture: the biological adapter
### supplies hidden state; the observation model supplies P(observation | state).
### ObservationModel may be:
###   * NULL: use q stored in the adapter (legacy/external contract);
###   * a function(Adapter, ModelResults, Timestep, Source), returning Nperm q's;
###   * a list with Background / InfoTriggered entries.  Each entry may be a
###     function with the same signature, a scalar q, an Nperm vector, or an
###     Nperm x Ntimesteps matrix of q values.
.INApestPoAResolveNoDetection <- function(Adapter, ModelResults, ObservationModel,
                                          Timestep, Source)
{
  if(!(Source %in% c("Background","InfoTriggered")))
    .INApestPoAStop("Unknown observation source: ", Source)

  legacy_q <- if(Source == "Background")
    Adapter$BackgroundNoDetectionProbability else Adapter$InfoTriggeredNoDetectionProbability

  spec <- NULL
  if(is.null(ObservationModel))
    spec <- legacy_q
  else if(is.function(ObservationModel))
    spec <- ObservationModel(Adapter=Adapter, ModelResults=ModelResults,
                             Timestep=Timestep, Source=Source)
  else if(is.list(ObservationModel))
  {
    spec <- ObservationModel[[Source]]
    if(is.null(spec))
    {
      old <- if(Source == "Background") "Surveillance" else "Information"
      spec <- ObservationModel[[old]]
      if(!is.null(spec)) warning("ObservationModel$", old, " is deprecated; use ObservationModel$", Source, ".", call.=FALSE)
    }
    if(is.function(spec))
      spec <- spec(Adapter=Adapter, ModelResults=ModelResults,
                   Timestep=Timestep, Source=Source)
  }
  else
    .INApestPoAStop("ObservationModel must be NULL, a function, or a list")

  if(is.null(spec)) return(NULL)
  if(is.matrix(spec) || (is.array(spec) && length(dim(spec)) == 2L))
  {
    if(!identical(dim(spec), c(Adapter$Nperm, Adapter$Ntimesteps)))
      .INApestPoAStop(Source, " observation q matrix must be Nperm x Ntimesteps")
    q <- spec[, Timestep]
  }
  else
  {
    q <- as.numeric(spec)
    if(length(q) == 1L) q <- rep(q, Adapter$Nperm)
    if(length(q) != Adapter$Nperm)
      .INApestPoAStop(Source, " observation model must return one q per particle")
  }
  if(any(is.na(q) | !is.finite(q) | q < 0 | q > 1))
    .INApestPoAStop(Source, " no-detection probabilities must be finite values in [0,1]")
  q
}


.INApestPoAValidateResults <- function(ModelResults, Required, AdapterName)
{
  if(!is.list(ModelResults))
    .INApestPoAStop(AdapterName, ": ModelResults must be a list returned by an INApest model")

  Missing <- setdiff(Required, names(ModelResults))
  if(length(Missing))
    .INApestPoAStop(
      AdapterName, ": ModelResults is missing required field(s): ",
      paste(Missing, collapse = ", ")
    )
}


### Aggregate an arbitrary numeric/logical array over every dimension except
### TimeDim and PermDim.  Output is always permutations x timesteps.
.INApestPoAAggregateArray <- function(x, TimeDim, PermDim, Name = deparse(substitute(x)))
{
  if(is.null(x)) return(NULL)
  d <- dim(x)
  if(is.null(d))
    .INApestPoAStop(Name, " must be an array")
  if(TimeDim == PermDim || TimeDim < 1L || PermDim < 1L ||
     TimeDim > length(d) || PermDim > length(d))
    .INApestPoAStop("Invalid TimeDim/PermDim for ", Name)

  Keep <- c(PermDim, TimeDim)
  SumDims <- setdiff(seq_along(d), Keep)

  if(length(SumDims) == 0L)
  {
    y <- aperm(x, match(c(PermDim, TimeDim), seq_along(d)))
    return(matrix(as.numeric(y), nrow = d[PermDim], ncol = d[TimeDim]))
  }

  y <- apply(x, Keep, sum, na.rm = TRUE)

  ### apply() can drop dimensions when either kept dimension is length one.
  matrix(as.numeric(y), nrow = d[PermDim], ncol = d[TimeDim])
}


### Convert a unit x timestep x permutation style array (or an array with
### additional unit dimensions) to Unit x Timestep x Permutation.
.INApestPoAUnitArray <- function(x, TimeDim, PermDim, Name = deparse(substitute(x)))
{
  d <- dim(x)
  if(is.null(d)) .INApestPoAStop(Name, " must be an array")

  UnitDims <- setdiff(seq_along(d), c(TimeDim, PermDim))
  NewOrder <- c(UnitDims, TimeDim, PermDim)
  y <- aperm(x, NewOrder)
  Nunit <- if(length(UnitDims)) prod(d[UnitDims]) else 1L
  array(y, dim = c(Nunit, d[TimeDim], d[PermDim]))
}


.INApestPoAInitialAggregate <- function(x, Nperm)
{
  if(is.null(x)) return(rep(0, Nperm))
  if(is.null(dim(x)))
  {
    if(length(x) == Nperm) return(as.numeric(x))
    if(Nperm == 1L) return(sum(x, na.rm = TRUE))
    .INApestPoAStop("Cannot identify permutation dimension for initial detection output")
  }

  d <- dim(x)
  PermCandidates <- which(d == Nperm)
  if(!length(PermCandidates))
  {
    if(Nperm == 1L) return(sum(x, na.rm = TRUE))
    .INApestPoAStop("Cannot identify permutation dimension for initial detection output")
  }

  ### INApest initial detection outputs place permutation last.
  PermDim <- tail(PermCandidates, 1L)
  Keep <- PermDim
  y <- apply(x, Keep, sum, na.rm = TRUE)
  as.numeric(y)
}


.INApestPoAAggregateDetectionEvent <- function(x, TimeDim, PermDim, Nperm, Ntimesteps, Name)
{
  if(is.null(x)) return(NULL)
  d <- dim(x)
  if(is.null(d)) .INApestPoAStop(Name, " must be an array")

  ### Canonical new event contract is unit/node x timestep x permutation, even
  ### when the biological state has extra land-use or stage dimensions.
  if(length(d) == 3L && d[2L] == Ntimesteps && d[3L] == Nperm)
    return(.INApestPoAAggregateArray(x, 2L, 3L, Name))

  ### Backwards-compatible source-specific arrays may still share the dimensions
  ### of the biological state (for example node x land-use x time x permutation).
  if(TimeDim <= length(d) && PermDim <= length(d) &&
     d[TimeDim] == Ntimesteps && d[PermDim] == Nperm)
    return(.INApestPoAAggregateArray(x, TimeDim, PermDim, Name))

  .INApestPoAStop(
    Name, " has unsupported dimensions for the observation-event contract. ",
    "Expected unit/node x timestep x permutation or an architecture-compatible legacy array."
  )
}


.INApestPoAGetDetectionArray <- function(ModelResults, Source,
                                           TimeDim, PermDim, Nperm, Ntimesteps)
{
  if(Source == "Background")
  {
    explicit_names <- c("BackgroundDetectedResults", "SurveillanceDetectedResults")
    for(nm in explicit_names)
      if(!is.null(ModelResults[[nm]]))
        return(list(
          Values=.INApestPoAAggregateDetectionEvent(ModelResults[[nm]], TimeDim, PermDim, Nperm, Ntimesteps, nm),
          Contract="explicit", Field=nm
        ))

    ### Legacy DetectedResults is usually HaveInfo * Invaded: a persistent state,
    ### not a new observation event.  It is kept only for one-round backwards
    ### compatibility and is never accepted for sequential compatible-history
    ### propagation.
    if(!is.null(ModelResults$DetectedResults))
    {
      warning(
        "Explicit BackgroundDetectedResults is unavailable; using legacy DetectedResults for backwards compatibility. ",
        "This is a state, not a reliable new-detection event, and cannot support sequential propagation.",
        call.=FALSE
      )
      return(list(
        Values=.INApestPoAAggregateDetectionEvent(ModelResults$DetectedResults, TimeDim, PermDim, Nperm, Ntimesteps, "DetectedResults"),
        Contract="legacy_state", Field="DetectedResults"
      ))
    }
  }
  else if(Source == "InfoTriggered")
  {
    explicit_names <- c("InfoTriggeredDetectedResults", "InformationDetectedResults")
    for(nm in explicit_names)
      if(!is.null(ModelResults[[nm]]))
      {
        if(nm == "InformationDetectedResults")
          warning("InformationDetectedResults is a deprecated name; use InfoTriggeredDetectedResults.", call.=FALSE)
        return(list(
          Values=.INApestPoAAggregateDetectionEvent(ModelResults[[nm]], TimeDim, PermDim, Nperm, Ntimesteps, nm),
          Contract="explicit", Field=nm
        ))
      }
  }
  else .INApestPoAStop("Unknown detection source: ", Source)

  list(Values=matrix(0, nrow=Nperm, ncol=Ntimesteps), Contract="none", Field=NA_character_)
}


### Standardise hidden state exposed to the observation model.
### ObservationAbundance is observation-unit x timestep x particle.  It can be
### binary occupancy, abundance, node x land-use abundance, or node x stage
### abundance flattened to observation units.
.INApestPoAInformationStateForObservation <- function(ModelResults, ObservationUnitLabels,
                                                       Ntimesteps, Nperm)
{
  x <- ModelResults$InformationStateBeforeSurveillanceResults
  if(is.null(x)) x <- ModelResults$HaveInfoBeforeSurveillanceResults
  if(is.null(x)) return(NULL)
  d <- dim(x)
  if(length(d) != 3L || d[2L] != Ntimesteps || d[3L] != Nperm)
    .INApestPoAStop(
      "InformationStateBeforeSurveillanceResults must be unit/node x timestep x permutation"
    )
  Nobs <- nrow(ObservationUnitLabels)
  if(d[1L] == Nobs) return(array(as.logical(x), dim=d))
  if("Node" %in% names(ObservationUnitLabels))
  {
    node <- as.integer(ObservationUnitLabels$Node)
    if(all(!is.na(node)) && min(node) >= 1L && max(node) <= d[1L])
      return(array(as.logical(x[node,,,drop=FALSE]), dim=c(Nobs,Ntimesteps,Nperm)))
  }
  .INApestPoAStop(
    "Cannot map InformationStateBeforeSurveillanceResults to observation units. ",
    "Supply one row per observation unit or ObservationUnitLabels$Node."
  )
}


.INApestPoAResolveDetectionProb <- function(Spec, Adapter, Timestep, Source)
{
  if(is.function(Spec))
  {
    fm <- names(formals(Spec))
    a <- list(Adapter=Adapter, Timestep=Timestep, Source=Source)
    if(!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a),fm)]
    Spec <- do.call(Spec,a)
  }
  Nobs <- dim(Adapter$ObservationAbundance)[1L]
  if(is.matrix(Spec) || (is.array(Spec) && length(dim(Spec)) == 2L))
  {
    if(!identical(dim(Spec), c(Nobs, Adapter$Ntimesteps)))
      .INApestPoAStop(Source, " DetectionProb matrix must be observation-unit x timestep")
    p <- Spec[,Timestep]
  }
  else
  {
    p <- as.numeric(Spec)
    if(length(p) == 1L) p <- rep(p,Nobs)
    if(length(p) != Nobs)
      .INApestPoAStop(Source, " DetectionProb must be scalar, one value per observation unit, an observation-unit x timestep matrix, or a resolver function")
  }
  if(any(is.na(p) | !is.finite(p) | p < 0 | p > 1))
    .INApestPoAStop(Source, " DetectionProb must resolve to values in [0,1]")
  p
}


### Built-in independent-detection observation model.
### DetectionUnit="individual" gives q = product_u (1-p_u)^N_u.
### DetectionUnit="occupied_unit" gives one detection opportunity per occupied
### unit.  InfoTriggered surveillance is gated by information that existed BEFORE
### the current surveillance round, so a new background detection cannot activate
### targeted surveillance retrospectively within the same round.
INApestPoAIndependentDetectionModel <- function(
  BackgroundDetectionProb,
  InfoTriggeredDetectionProb = 0,
  DetectionUnit = c("individual", "occupied_unit"),
  InfoTriggeredRequiresInformation = TRUE
)
{
  DetectionUnit <- match.arg(DetectionUnit)
  force(BackgroundDetectionProb)
  force(InfoTriggeredDetectionProb)
  force(InfoTriggeredRequiresInformation)

  function(Adapter, ModelResults, Timestep, Source)
  {
    if(is.null(Adapter$ObservationAbundance))
      .INApestPoAStop(
        "The selected adapter does not expose ObservationAbundance. ",
        "Supply a custom ObservationModel function for this architecture."
      )
    A <- Adapter$ObservationAbundance[,Timestep,,drop=FALSE]
    A <- matrix(as.numeric(A), nrow=dim(Adapter$ObservationAbundance)[1L], ncol=Adapter$Nperm)
    if(any(!is.finite(A) | A < 0)) .INApestPoAStop("ObservationAbundance contains invalid values")

    Spec <- if(Source == "Background") BackgroundDetectionProb else InfoTriggeredDetectionProb
    p <- .INApestPoAResolveDetectionProb(Spec, Adapter, Timestep, Source)
    P <- matrix(p, nrow=length(p), ncol=Adapter$Nperm)

    if(Source == "InfoTriggered" && isTRUE(InfoTriggeredRequiresInformation))
    {
      H <- Adapter$InformationStateBeforeSurveillance
      if(is.null(H))
        .INApestPoAStop(
          "InfoTriggered surveillance requires InformationStateBeforeSurveillanceResults ",
          "(or a custom ObservationModel that defines its own gating)."
        )
      gate <- H[,Timestep,,drop=FALSE]
      gate <- matrix(as.numeric(gate > 0), nrow=length(p), ncol=Adapter$Nperm)
      P <- P * gate
    }

    if(DetectionUnit == "individual")
      q <- apply((1-P)^A, 2L, prod)
    else
      q <- apply(ifelse(A > 0, 1-P, 1), 2L, prod)
    as.numeric(q)
  }
}


### Map realised per-unit detection probabilities recorded by the biological
### engine to the same observation-unit x timestep x particle layout as
### ObservationAbundance.  The engine records probabilities, not PoA likelihoods;
### the observation layer still computes P(no detection | hidden state).
.INApestPoAProbabilityForObservation <- function(x, TimeDim, PermDim,
                                                  ObservationAbundance,
                                                  Name)
{
  if(is.null(x)) return(NULL)
  z <- .INApestPoAUnitArray(x, TimeDim, PermDim, Name)
  if(is.null(ObservationAbundance) || !identical(dim(z), dim(ObservationAbundance)))
    .INApestPoAStop(
      Name, " does not map to the ObservationAbundance dimensions. ",
      "Recorded detection probabilities must have one value per observation unit, timestep and particle."
    )
  if(any(is.na(z) | !is.finite(z) | z < 0 | z > 1))
    .INApestPoAStop(Name, " must contain finite probabilities in [0,1]")
  z
}


.INApestPoARecordedDetectionAvailable <- function(Adapter, Source)
{
  if(!(Source %in% c("Background","InfoTriggered"))) return(FALSE)
  P <- if(Source == "Background") Adapter$BackgroundDetectionProbability else Adapter$InfoTriggeredDetectionProbability
  if(!is.null(P))
  {
    if(Source == "InfoTriggered" && is.null(Adapter$InformationStateBeforeSurveillance)) return(FALSE)
    return(TRUE)
  }

  ### Point engines retain the realised per-point detection probability in
  ### PointHistory because the number and attributes of points change through
  ### time.  If no observation opportunities exist anywhere, q is identically 1
  ### and the probability column is not required.
  H <- Adapter$PointHistory
  if(is.data.frame(H))
  {
    if(!is.null(Adapter$ObservationAbundance) && all(Adapter$ObservationAbundance == 0)) return(TRUE)
    need <- if(Source == "Background")
      c("perm","timestep","background_detection_prob") else
      c("perm","timestep","info_triggered_detection_prob","have_info_before_surveillance")
    return(all(need %in% names(H)))
  }
  FALSE
}


### Observation model using the exact detection probabilities realised by each
### INApest simulation.  This preserves DetectionSD uncertainty and any spatial,
### stage, land-use or point heterogeneity while keeping likelihood calculation
### outside the biological engine.
INApestPoARecordedDetectionModel <- function()
{
  function(Adapter, ModelResults, Timestep, Source)
  {
    if(!(Source %in% c("Background","InfoTriggered")))
      .INApestPoAStop("Unknown observation source: ", Source)

    Pstored <- if(Source == "Background") Adapter$BackgroundDetectionProbability else Adapter$InfoTriggeredDetectionProbability
    if(!is.null(Pstored))
    {
      A <- Adapter$ObservationAbundance[,Timestep,,drop=FALSE]
      A <- matrix(as.numeric(A), nrow=dim(Adapter$ObservationAbundance)[1L], ncol=Adapter$Nperm)
      P <- Pstored[,Timestep,,drop=FALSE]
      P <- matrix(as.numeric(P), nrow=dim(Pstored)[1L], ncol=Adapter$Nperm)
      if(Source == "InfoTriggered")
      {
        H <- Adapter$InformationStateBeforeSurveillance
        if(is.null(H))
          .INApestPoAStop("Recorded InfoTriggered probabilities require pre-surveillance information state")
        Gate <- H[,Timestep,,drop=FALSE]
        Gate <- matrix(as.numeric(Gate > 0), nrow=dim(H)[1L], ncol=Adapter$Nperm)
        P <- P * Gate
      }
      return(as.numeric(apply((1-P)^A, 2L, prod)))
    }

    ### Dynamic point models: evaluate each particle directly from the point
    ### snapshot at the surveillance timestep.  One point is one detection
    ### opportunity, so q is the product of its realised no-detection terms.
    H <- Adapter$PointHistory
    if(!is.data.frame(H))
      .INApestPoAStop("Recorded detection probabilities are unavailable for this adapter")
    out <- rep(1, Adapter$Nperm)
    if(!nrow(H)) return(out)
    if(!all(c("perm","timestep") %in% names(H)))
      .INApestPoAStop("PointHistory must contain perm and timestep for recorded observation likelihoods")

    pfield <- if(Source == "Background") "background_detection_prob" else "info_triggered_detection_prob"
    if(!(pfield %in% names(H)))
    {
      if(!is.null(Adapter$ObservationAbundance) && all(Adapter$ObservationAbundance[,Timestep,,drop=FALSE] == 0)) return(out)
      .INApestPoAStop("PointHistory is missing recorded field '", pfield, "'")
    }
    if(Source == "InfoTriggered" && !("have_info_before_surveillance" %in% names(H)))
      .INApestPoAStop("PointHistory is missing have_info_before_surveillance for InfoTriggered likelihoods")

    for(pp in seq_len(Adapter$Nperm))
    {
      z <- H[as.integer(H$perm) == pp & as.integer(H$timestep) == Timestep, , drop=FALSE]
      if(!nrow(z)) { out[pp] <- 1; next }
      pv <- as.numeric(z[[pfield]])
      if(any(is.na(pv) | !is.finite(pv) | pv < 0 | pv > 1))
        .INApestPoAStop("PointHistory ", pfield, " must contain finite probabilities in [0,1]")
      if(Source == "InfoTriggered")
        pv <- pv * as.numeric(as.logical(z$have_info_before_surveillance))
      out[pp] <- prod(1-pv)
    }
    out
  }
}


.INApestPoAAdapterFinish <- function(AdapterName, ModelResults, ResidualPopulation,
                                      Absence, SurveillanceDetections,
                                      InformationDetections, LegacyDetected = NULL,
                                      InitialSurveillanceDetections = NULL,
                                      InitialInformationDetections = NULL,
                                      UnitPresence = NULL, UnitLabels = NULL,
                                      StageResidualPopulation = NULL,
                                      PointHistory = NULL,
                                      ObservationAbundance = NULL,
                                      ObservationUnitLabels = NULL,
                                      InformationStateBeforeSurveillance = NULL,
                                      BackgroundDetectionProbability = NULL,
                                      InfoTriggeredDetectionProbability = NULL,
                                      SurveillanceNoDetectionProbability = NULL,
                                      InformationNoDetectionProbability = NULL,
                                      BackgroundEventContract = "none",
                                      InfoTriggeredEventContract = "none",
                                      BackgroundEventField = NA_character_,
                                      InfoTriggeredEventField = NA_character_)
{
  ResidualPopulation <- as.matrix(ResidualPopulation)
  Absence <- as.matrix(Absence)
  SurveillanceDetections <- as.matrix(SurveillanceDetections)
  InformationDetections <- as.matrix(InformationDetections)

  Nperm <- nrow(Absence)
  Ntimesteps <- ncol(Absence)
  ExpectedDim <- c(Nperm, Ntimesteps)

  for(Name in c("ResidualPopulation", "SurveillanceDetections", "InformationDetections"))
  {
    z <- get(Name)
    if(!identical(dim(z), ExpectedDim))
      .INApestPoAStop(AdapterName, ": ", Name, " does not match absence dimensions")
  }

  if(any(!is.finite(ResidualPopulation)) || any(ResidualPopulation < 0))
    .INApestPoAStop(AdapterName, ": residual population contains invalid values")

  if(!is.null(ObservationAbundance))
  {
    dObs <- dim(ObservationAbundance)
    if(length(dObs) != 3L || dObs[2L] != Ntimesteps || dObs[3L] != Nperm)
      .INApestPoAStop(AdapterName, ": ObservationAbundance must be observation-unit x timestep x permutation")
    if(any(!is.finite(ObservationAbundance)) || any(ObservationAbundance < 0))
      .INApestPoAStop(AdapterName, ": ObservationAbundance contains invalid values")
    if(is.null(ObservationUnitLabels) || !is.data.frame(ObservationUnitLabels) || nrow(ObservationUnitLabels) != dObs[1L])
      .INApestPoAStop(AdapterName, ": ObservationUnitLabels must have one row per observation unit")
  }
  if(!is.null(InformationStateBeforeSurveillance))
  {
    if(is.null(ObservationAbundance) || !identical(dim(InformationStateBeforeSurveillance), dim(ObservationAbundance)))
      .INApestPoAStop(AdapterName, ": InformationStateBeforeSurveillance must match ObservationAbundance dimensions")
    InformationStateBeforeSurveillance <- array(as.logical(InformationStateBeforeSurveillance), dim=dim(ObservationAbundance))
  }

  for(Name in c("BackgroundDetectionProbability", "InfoTriggeredDetectionProbability"))
  {
    z <- get(Name)
    if(!is.null(z))
    {
      if(is.null(ObservationAbundance) || !identical(dim(z), dim(ObservationAbundance)))
        .INApestPoAStop(AdapterName, ": ", Name, " must match ObservationAbundance dimensions")
      if(any(is.na(z) | !is.finite(z) | z < 0 | z > 1))
        .INApestPoAStop(AdapterName, ": ", Name, " must contain finite probabilities in [0,1]")
    }
  }

  for(Name in c("SurveillanceNoDetectionProbability", "InformationNoDetectionProbability"))
  {
    z <- get(Name)
    if(!is.null(z))
    {
      z <- as.matrix(z)
      if(!identical(dim(z), ExpectedDim))
        .INApestPoAStop(AdapterName, ": ", Name, " does not match absence dimensions")
      if(any(!is.na(z) & (!is.finite(z) | z < 0 | z > 1)))
        .INApestPoAStop(AdapterName, ": ", Name, " must contain probabilities in [0,1] or NA")
      assign(Name, z)
    }
  }

  list(
    AdapterName = AdapterName,
    ModelName = if(!is.null(ModelResults$ModelName)) ModelResults$ModelName else NA_character_,
    Nperm = Nperm,
    Ntimesteps = Ntimesteps,
    Absence = matrix(as.logical(Absence), nrow = Nperm, ncol = Ntimesteps),
    ResidualPopulation = ResidualPopulation,
    BackgroundDetections = SurveillanceDetections,
    InfoTriggeredDetections = InformationDetections,
    ### Deprecated aliases retained for frozen analyses.
    SurveillanceDetections = SurveillanceDetections,
    InformationDetections = InformationDetections,
    LegacyDetected = LegacyDetected,
    InitialSurveillanceDetections = if(is.null(InitialSurveillanceDetections)) rep(0, Nperm) else as.numeric(InitialSurveillanceDetections),
    InitialInformationDetections = if(is.null(InitialInformationDetections)) rep(0, Nperm) else as.numeric(InitialInformationDetections),
    UnitPresence = UnitPresence,
    UnitLabels = UnitLabels,
    StageResidualPopulation = StageResidualPopulation,
    PointHistory = PointHistory,
    ObservationAbundance = ObservationAbundance,
    ObservationUnitLabels = ObservationUnitLabels,
    InformationStateBeforeSurveillance = InformationStateBeforeSurveillance,
    BackgroundDetectionProbability = BackgroundDetectionProbability,
    InfoTriggeredDetectionProbability = InfoTriggeredDetectionProbability,
    BackgroundNoDetectionProbability = SurveillanceNoDetectionProbability,
    InfoTriggeredNoDetectionProbability = InformationNoDetectionProbability,
    ### Deprecated aliases retained for frozen analyses.
    SurveillanceNoDetectionProbability = SurveillanceNoDetectionProbability,
    InformationNoDetectionProbability = InformationNoDetectionProbability,
    BackgroundEventContract = BackgroundEventContract,
    InfoTriggeredEventContract = InfoTriggeredEventContract,
    BackgroundEventField = BackgroundEventField,
    InfoTriggeredEventField = InfoTriggeredEventField
  )
}


###############################################################################
### INApest output adapters
###############################################################################

.INApestPoAAdapterBinaryNode <- function(ModelResults)
{
  AdapterName <- ".INApestPoAAdapterBinaryNode"
  .INApestPoAValidateResults(ModelResults, c("InvasionResults"), AdapterName)

  d <- dim(ModelResults$InvasionResults)
  if(length(d) != 3L)
    .INApestPoAStop(AdapterName, ": InvasionResults must be node x timestep x permutation")

  Ntimesteps <- d[2L]
  Nperm <- d[3L]
  Residual <- .INApestPoAAggregateArray(ModelResults$InvasionResults, 2L, 3L, "InvasionResults")
  BackgroundEvent <- .INApestPoAGetDetectionArray(ModelResults, "Background", 2L, 3L, Nperm, Ntimesteps)
  InfoTriggeredEvent <- .INApestPoAGetDetectionArray(ModelResults, "InfoTriggered", 2L, 3L, Nperm, Ntimesteps)
  Legacy <- if(!is.null(ModelResults$DetectedResults)) .INApestPoAAggregateArray(ModelResults$DetectedResults, 2L, 3L, "DetectedResults") else NULL

  .INApestPoAAdapterFinish(
    AdapterName, ModelResults, Residual, Residual <= 0,
    BackgroundEvent$Values, InfoTriggeredEvent$Values, Legacy,
    .INApestPoAInitialAggregate(ModelResults$InitialSurveillanceDetected, Nperm),
    .INApestPoAInitialAggregate(ModelResults$InitialInformationDetected, Nperm),
    UnitPresence = .INApestPoAUnitArray(ModelResults$InvasionResults > 0, 2L, 3L),
    UnitLabels = data.frame(Unit = seq_len(d[1L]), Node = seq_len(d[1L])),
    ObservationAbundance = .INApestPoAUnitArray(ModelResults$InvasionResults, 2L, 3L),
    ObservationUnitLabels = data.frame(Unit=seq_len(d[1L]), Node=seq_len(d[1L])),
    InformationStateBeforeSurveillance = .INApestPoAInformationStateForObservation(
      ModelResults, data.frame(Unit=seq_len(d[1L]), Node=seq_len(d[1L])), Ntimesteps, Nperm),
    BackgroundDetectionProbability = .INApestPoAProbabilityForObservation(
      ModelResults$BackgroundDetectionProbabilityResults, 2L, 3L,
      .INApestPoAUnitArray(ModelResults$InvasionResults, 2L, 3L, "ObservationAbundance"), "BackgroundDetectionProbabilityResults"),
    InfoTriggeredDetectionProbability = .INApestPoAProbabilityForObservation(
      ModelResults$InfoTriggeredDetectionProbabilityResults, 2L, 3L,
      .INApestPoAUnitArray(ModelResults$InvasionResults, 2L, 3L, "ObservationAbundance"), "InfoTriggeredDetectionProbabilityResults"),
    SurveillanceNoDetectionProbability = .INApestPoAFirstNonNull(ModelResults$BackgroundNoDetectionProbability, ModelResults$SurveillanceNoDetectionProbability),
    InformationNoDetectionProbability = .INApestPoAFirstNonNull(ModelResults$InfoTriggeredNoDetectionProbability, ModelResults$InformationNoDetectionProbability),
    BackgroundEventContract = BackgroundEvent$Contract,
    InfoTriggeredEventContract = InfoTriggeredEvent$Contract,
    BackgroundEventField = BackgroundEvent$Field,
    InfoTriggeredEventField = InfoTriggeredEvent$Field
  )
}


.INApestPoAAdapterPopulationNode <- function(ModelResults)
{
  AdapterName <- ".INApestPoAAdapterPopulationNode"
  .INApestPoAValidateResults(ModelResults, c("PopulationResults"), AdapterName)

  d <- dim(ModelResults$PopulationResults)
  if(length(d) != 3L)
    .INApestPoAStop(AdapterName, ": PopulationResults must be node x timestep x permutation")

  Ntimesteps <- d[2L]
  Nperm <- d[3L]
  Residual <- .INApestPoAAggregateArray(ModelResults$PopulationResults, 2L, 3L, "PopulationResults")
  BackgroundEvent <- .INApestPoAGetDetectionArray(ModelResults, "Background", 2L, 3L, Nperm, Ntimesteps)
  InfoTriggeredEvent <- .INApestPoAGetDetectionArray(ModelResults, "InfoTriggered", 2L, 3L, Nperm, Ntimesteps)
  Legacy <- if(!is.null(ModelResults$DetectedResults)) .INApestPoAAggregateArray(ModelResults$DetectedResults, 2L, 3L, "DetectedResults") else NULL

  .INApestPoAAdapterFinish(
    AdapterName, ModelResults, Residual, Residual <= 0,
    BackgroundEvent$Values, InfoTriggeredEvent$Values, Legacy,
    .INApestPoAInitialAggregate(ModelResults$InitialSurveillanceDetected, Nperm),
    .INApestPoAInitialAggregate(ModelResults$InitialInformationDetected, Nperm),
    UnitPresence = .INApestPoAUnitArray(ModelResults$PopulationResults > 0, 2L, 3L),
    UnitLabels = data.frame(Unit = seq_len(d[1L]), Node = seq_len(d[1L])),
    ObservationAbundance = .INApestPoAUnitArray(ModelResults$PopulationResults, 2L, 3L),
    ObservationUnitLabels = data.frame(Unit=seq_len(d[1L]), Node=seq_len(d[1L])),
    InformationStateBeforeSurveillance = .INApestPoAInformationStateForObservation(
      ModelResults, data.frame(Unit=seq_len(d[1L]), Node=seq_len(d[1L])), Ntimesteps, Nperm),
    BackgroundDetectionProbability = .INApestPoAProbabilityForObservation(
      ModelResults$BackgroundDetectionProbabilityResults, 2L, 3L,
      .INApestPoAUnitArray(ModelResults$PopulationResults, 2L, 3L, "ObservationAbundance"), "BackgroundDetectionProbabilityResults"),
    InfoTriggeredDetectionProbability = .INApestPoAProbabilityForObservation(
      ModelResults$InfoTriggeredDetectionProbabilityResults, 2L, 3L,
      .INApestPoAUnitArray(ModelResults$PopulationResults, 2L, 3L, "ObservationAbundance"), "InfoTriggeredDetectionProbabilityResults"),
    SurveillanceNoDetectionProbability = .INApestPoAFirstNonNull(ModelResults$BackgroundNoDetectionProbability, ModelResults$SurveillanceNoDetectionProbability),
    InformationNoDetectionProbability = .INApestPoAFirstNonNull(ModelResults$InfoTriggeredNoDetectionProbability, ModelResults$InformationNoDetectionProbability),
    BackgroundEventContract = BackgroundEvent$Contract,
    InfoTriggeredEventContract = InfoTriggeredEvent$Contract,
    BackgroundEventField = BackgroundEvent$Field,
    InfoTriggeredEventField = InfoTriggeredEvent$Field
  )
}


.INApestPoAAdapterMultipleLandUse <- function(ModelResults)
{
  AdapterName <- ".INApestPoAAdapterMultipleLandUse"
  .INApestPoAValidateResults(ModelResults, c("PopulationResults"), AdapterName)

  d <- dim(ModelResults$PopulationResults)
  if(length(d) != 4L)
    .INApestPoAStop(AdapterName, ": PopulationResults must be node x land-use x timestep x permutation")

  Ntimesteps <- d[3L]
  Nperm <- d[4L]
  Residual <- .INApestPoAAggregateArray(ModelResults$PopulationResults, 3L, 4L, "PopulationResults")
  BackgroundEvent <- .INApestPoAGetDetectionArray(ModelResults, "Background", 3L, 4L, Nperm, Ntimesteps)
  InfoTriggeredEvent <- .INApestPoAGetDetectionArray(ModelResults, "InfoTriggered", 3L, 4L, Nperm, Ntimesteps)
  Legacy <- if(!is.null(ModelResults$DetectedResults))
    .INApestPoAAggregateDetectionEvent(ModelResults$DetectedResults, 3L, 4L, Nperm, Ntimesteps, "DetectedResults") else NULL

  Labels <- expand.grid(
    Node = seq_len(d[1L]),
    LandUse = seq_len(d[2L]),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  Labels$Unit <- seq_len(nrow(Labels))
  Labels <- Labels[, c("Unit", "Node", "LandUse")]

  .INApestPoAAdapterFinish(
    AdapterName, ModelResults, Residual, Residual <= 0,
    BackgroundEvent$Values, InfoTriggeredEvent$Values, Legacy,
    .INApestPoAInitialAggregate(ModelResults$InitialSurveillanceDetected, Nperm),
    .INApestPoAInitialAggregate(ModelResults$InitialInformationDetected, Nperm),
    UnitPresence = .INApestPoAUnitArray(ModelResults$PopulationResults > 0, 3L, 4L),
    UnitLabels = Labels,
    ObservationAbundance = .INApestPoAUnitArray(ModelResults$PopulationResults, 3L, 4L),
    ObservationUnitLabels = Labels,
    InformationStateBeforeSurveillance = .INApestPoAInformationStateForObservation(
      ModelResults, Labels, Ntimesteps, Nperm),
    BackgroundDetectionProbability = .INApestPoAProbabilityForObservation(
      ModelResults$BackgroundDetectionProbabilityResults, 3L, 4L,
      .INApestPoAUnitArray(ModelResults$PopulationResults, 3L, 4L, "ObservationAbundance"), "BackgroundDetectionProbabilityResults"),
    InfoTriggeredDetectionProbability = .INApestPoAProbabilityForObservation(
      ModelResults$InfoTriggeredDetectionProbabilityResults, 3L, 4L,
      .INApestPoAUnitArray(ModelResults$PopulationResults, 3L, 4L, "ObservationAbundance"), "InfoTriggeredDetectionProbabilityResults"),
    SurveillanceNoDetectionProbability = .INApestPoAFirstNonNull(ModelResults$BackgroundNoDetectionProbability, ModelResults$SurveillanceNoDetectionProbability),
    InformationNoDetectionProbability = .INApestPoAFirstNonNull(ModelResults$InfoTriggeredNoDetectionProbability, ModelResults$InformationNoDetectionProbability),
    BackgroundEventContract = BackgroundEvent$Contract,
    InfoTriggeredEventContract = InfoTriggeredEvent$Contract,
    BackgroundEventField = BackgroundEvent$Field,
    InfoTriggeredEventField = InfoTriggeredEvent$Field
  )
}


.INApestPoAAdapterTransitionMatrix <- function(ModelResults)
{
  AdapterName <- ".INApestPoAAdapterTransitionMatrix"
  .INApestPoAValidateResults(ModelResults, c("PopulationStageResults"), AdapterName)

  d <- dim(ModelResults$PopulationStageResults)
  if(length(d) != 4L)
    .INApestPoAStop(AdapterName, ": PopulationStageResults must be node x stage x timestep x permutation")

  Nnode <- d[1L]
  Nstage <- d[2L]
  Ntimesteps <- d[3L]
  Nperm <- d[4L]

  Residual <- .INApestPoAAggregateArray(ModelResults$PopulationStageResults, 3L, 4L, "PopulationStageResults")

  ### Stage totals as permutation x timestep x stage.
  StageTotals <- apply(ModelResults$PopulationStageResults, c(2L, 3L, 4L), sum, na.rm = TRUE)
  if(length(dim(StageTotals)) < 3L)
    StageTotals <- array(StageTotals, dim = c(Nstage, Ntimesteps, Nperm))
  StageResidual <- aperm(StageTotals, c(3L, 2L, 1L))

  BackgroundEvent <- .INApestPoAGetDetectionArray(ModelResults, "Background", 2L, 3L, Nperm, Ntimesteps)
  InfoTriggeredEvent <- .INApestPoAGetDetectionArray(ModelResults, "InfoTriggered", 2L, 3L, Nperm, Ntimesteps)
  Legacy <- if(!is.null(ModelResults$DetectedResults)) .INApestPoAAggregateArray(ModelResults$DetectedResults, 2L, 3L, "DetectedResults") else NULL

  UnitState <- if(!is.null(ModelResults$PopulationResults)) ModelResults$PopulationResults else
    apply(ModelResults$PopulationStageResults, c(1L, 3L, 4L), sum, na.rm = TRUE)
  if(length(dim(UnitState)) < 3L)
    UnitState <- array(UnitState, dim = c(Nnode, Ntimesteps, Nperm))

  ObservationLabels <- expand.grid(
    Node=seq_len(Nnode), Stage=seq_len(Nstage), KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  ObservationLabels$Unit <- seq_len(nrow(ObservationLabels))
  ObservationLabels <- ObservationLabels[,c("Unit","Node","Stage")]

  .INApestPoAAdapterFinish(
    AdapterName, ModelResults, Residual, Residual <= 0,
    BackgroundEvent$Values, InfoTriggeredEvent$Values, Legacy,
    .INApestPoAInitialAggregate(ModelResults$InitialSurveillanceDetected, Nperm),
    .INApestPoAInitialAggregate(ModelResults$InitialInformationDetected, Nperm),
    UnitPresence = .INApestPoAUnitArray(UnitState > 0, 2L, 3L),
    UnitLabels = data.frame(Unit = seq_len(Nnode), Node = seq_len(Nnode)),
    StageResidualPopulation = StageResidual,
    ObservationAbundance = .INApestPoAUnitArray(ModelResults$PopulationStageResults, 3L, 4L),
    ObservationUnitLabels = ObservationLabels,
    InformationStateBeforeSurveillance = .INApestPoAInformationStateForObservation(
      ModelResults, ObservationLabels, Ntimesteps, Nperm),
    BackgroundDetectionProbability = .INApestPoAProbabilityForObservation(
      ModelResults$BackgroundDetectionProbabilityResults, 3L, 4L,
      .INApestPoAUnitArray(ModelResults$PopulationStageResults, 3L, 4L, "ObservationAbundance"), "BackgroundDetectionProbabilityResults"),
    InfoTriggeredDetectionProbability = .INApestPoAProbabilityForObservation(
      ModelResults$InfoTriggeredDetectionProbabilityResults, 3L, 4L,
      .INApestPoAUnitArray(ModelResults$PopulationStageResults, 3L, 4L, "ObservationAbundance"), "InfoTriggeredDetectionProbabilityResults"),
    SurveillanceNoDetectionProbability = .INApestPoAFirstNonNull(ModelResults$BackgroundNoDetectionProbability, ModelResults$SurveillanceNoDetectionProbability),
    InformationNoDetectionProbability = .INApestPoAFirstNonNull(ModelResults$InfoTriggeredNoDetectionProbability, ModelResults$InformationNoDetectionProbability),
    BackgroundEventContract = BackgroundEvent$Contract,
    InfoTriggeredEventContract = InfoTriggeredEvent$Contract,
    BackgroundEventField = BackgroundEvent$Field,
    InfoTriggeredEventField = InfoTriggeredEvent$Field
  )
}


.INApestPoAPointMatrix <- function(Summary, ValueColumn, Nperm, Ntimesteps)
{
  out <- matrix(0, nrow = Nperm, ncol = Ntimesteps)
  if(!nrow(Summary) || !(ValueColumn %in% names(Summary))) return(out)

  p <- as.integer(Summary$perm)
  t <- as.integer(Summary$timestep)
  ok <- !is.na(p) & !is.na(t) & p >= 1L & p <= Nperm & t >= 1L & t <= Ntimesteps
  if(any(ok)) out[cbind(p[ok], t[ok])] <- as.numeric(Summary[[ValueColumn]][ok])
  out
}


.INApestPoAPointObservationStrata <- function(Summary, Nperm, Ntimesteps)
{
  total <- .INApestPoAPointMatrix(Summary, "n_end", Nperm, Ntimesteps)
  if(!("n_info_before_surveillance" %in% names(Summary)))
    return(list(
      Abundance=array(t(total), dim=c(1L,Ntimesteps,Nperm)),
      Information=NULL,
      Labels=data.frame(Unit=1L, InfoStratum="all", stringsAsFactors=FALSE),
      Contract="system_total_without_pre_surveillance_information"
    ))

  informed <- .INApestPoAPointMatrix(Summary, "n_info_before_surveillance", Nperm, Ntimesteps)
  if(any(informed < 0 | informed > total))
    .INApestPoAStop("Point Summary$n_info_before_surveillance must be between zero and n_end")
  uninformed <- total - informed

  abundance <- array(0, dim=c(2L,Ntimesteps,Nperm))
  abundance[1L,,] <- t(uninformed)
  abundance[2L,,] <- t(informed)
  info <- array(FALSE, dim=dim(abundance))
  info[2L,,] <- TRUE

  list(
    Abundance=abundance,
    Information=info,
    Labels=data.frame(Unit=1:2, InfoStratum=c("uninformed","informed"), stringsAsFactors=FALSE),
    Contract="pre_surveillance_information_strata"
  )
}


.INApestPoAAdapterPoint <- function(ModelResults)
{
  AdapterName <- ".INApestPoAAdapterPoint"
  .INApestPoAValidateResults(ModelResults, c("Summary"), AdapterName)

  Summary <- ModelResults$Summary
  if(!is.data.frame(Summary) || !all(c("perm", "timestep", "n_end") %in% names(Summary)))
    .INApestPoAStop(AdapterName, ": Summary must contain perm, timestep and n_end")

  Nperm <- if(nrow(Summary)) max(as.integer(Summary$perm), na.rm = TRUE) else 0L
  Ntimesteps <- if(nrow(Summary)) max(as.integer(Summary$timestep), na.rm = TRUE) else 0L
  if(Nperm < 1L || Ntimesteps < 1L)
    .INApestPoAStop(AdapterName, ": Summary contains no simulation timesteps")

  Residual <- .INApestPoAPointMatrix(Summary, "n_end", Nperm, Ntimesteps)
  BackgroundField <- if("n_new_background_detections" %in% names(Summary)) "n_new_background_detections" else
    if("n_new_surveillance_detections" %in% names(Summary)) "n_new_surveillance_detections" else
    if("n_new_detections" %in% names(Summary)) "n_new_detections" else NA_character_
  Background <- if(!is.na(BackgroundField)) .INApestPoAPointMatrix(Summary, BackgroundField, Nperm, Ntimesteps) else matrix(0,Nperm,Ntimesteps)
  InfoTriggeredField <- if("n_new_info_triggered_detections" %in% names(Summary)) "n_new_info_triggered_detections" else
    if("n_new_information_detections" %in% names(Summary)) "n_new_information_detections" else NA_character_
  InfoTriggered <- if(!is.na(InfoTriggeredField)) .INApestPoAPointMatrix(Summary, InfoTriggeredField, Nperm, Ntimesteps) else matrix(0,Nperm,Ntimesteps)
  Legacy <- .INApestPoAPointMatrix(Summary, "n_new_detections", Nperm, Ntimesteps)

  InitialSurv <- rep(0, Nperm)
  InitialInfo <- rep(0, Nperm)
  EventLog <- ModelResults$EventLog
  if(is.data.frame(EventLog) && nrow(EventLog) && all(c("perm", "timestep", "event") %in% names(EventLog)))
  {
    z <- EventLog[EventLog$timestep == 0 & EventLog$event %in% c("background_detection","detection"), , drop = FALSE]
    if(nrow(z)) InitialSurv <- tabulate(as.integer(z$perm), nbins = Nperm)
    z <- EventLog[EventLog$timestep == 0 & EventLog$event %in% c("info_triggered_detection","information_detection"), , drop = FALSE]
    if(nrow(z)) InitialInfo <- tabulate(as.integer(z$perm), nbins = Nperm)
  }

  PointObservation <- .INApestPoAPointObservationStrata(Summary, Nperm, Ntimesteps)

  .INApestPoAAdapterFinish(
    AdapterName, ModelResults, Residual, Residual <= 0,
    Background, InfoTriggered, Legacy,
    InitialSurv, InitialInfo,
    PointHistory = ModelResults$PointHistory,
    ObservationAbundance = PointObservation$Abundance,
    ObservationUnitLabels = PointObservation$Labels,
    InformationStateBeforeSurveillance = PointObservation$Information,
    SurveillanceNoDetectionProbability = .INApestPoAFirstNonNull(ModelResults$BackgroundNoDetectionProbability, ModelResults$SurveillanceNoDetectionProbability),
    InformationNoDetectionProbability = .INApestPoAFirstNonNull(ModelResults$InfoTriggeredNoDetectionProbability, ModelResults$InformationNoDetectionProbability),
    BackgroundEventContract = if(is.na(BackgroundField)) "none" else "explicit",
    InfoTriggeredEventContract = if(is.na(InfoTriggeredField)) "none" else "explicit",
    BackgroundEventField = BackgroundField,
    InfoTriggeredEventField = InfoTriggeredField
  )
}


.INApestPoAAdapterPointTransitionMatrix <- function(ModelResults)
{
  AdapterName <- ".INApestPoAAdapterPointTransitionMatrix"
  Base <- .INApestPoAAdapterPoint(ModelResults)
  Summary <- ModelResults$Summary

  StageCols <- grep("^stage[0-9]+_end$", names(Summary), value = TRUE)
  StageResidual <- NULL
  if(length(StageCols))
  {
    StageNumber <- as.integer(sub("^stage([0-9]+)_end$", "\\1", StageCols))
    StageCols <- StageCols[order(StageNumber)]
    Nstage <- length(StageCols)
    StageResidual <- array(0, dim = c(Base$Nperm, Base$Ntimesteps, Nstage))
    for(s in seq_len(Nstage))
      StageResidual[, , s] <- .INApestPoAPointMatrix(Summary, StageCols[s], Base$Nperm, Base$Ntimesteps)
  }

  Base$AdapterName <- AdapterName
  Base$StageResidualPopulation <- StageResidual
  Base
}


###############################################################################
### Observation-history and Bayesian conditioning helpers
###############################################################################

.INApestPoANormaliseObservationHistory <- function(ObservationHistory, PoAStartTimestep, Ntimesteps)
{
  if(is.null(ObservationHistory))
    return(data.frame(
      Round=integer(0), Timestep=integer(0),
      BackgroundDetections=numeric(0), InfoTriggeredDetections=numeric(0)
    ))

  if(!is.data.frame(ObservationHistory))
    .INApestPoAStop("ObservationHistory must be a data.frame")
  if(!("Timestep" %in% names(ObservationHistory)))
    .INApestPoAStop("ObservationHistory must contain a Timestep column")

  x <- ObservationHistory
  if(!("BackgroundDetections" %in% names(x)))
  {
    if("SurveillanceDetections" %in% names(x))
    {
      x$BackgroundDetections <- x$SurveillanceDetections
      warning("ObservationHistory$SurveillanceDetections is deprecated; use BackgroundDetections.", call.=FALSE)
    } else x$BackgroundDetections <- NA_real_
  }
  if(!("InfoTriggeredDetections" %in% names(x)))
  {
    if("InformationDetections" %in% names(x))
    {
      x$InfoTriggeredDetections <- x$InformationDetections
      warning(
        "ObservationHistory$InformationDetections is deprecated; use InfoTriggeredDetections. ",
        "The HaveInfo state itself is not a detection.", call.=FALSE
      )
    } else x$InfoTriggeredDetections <- NA_real_
  }

  x$Timestep <- as.integer(x$Timestep)
  if(any(is.na(x$Timestep)) || any(x$Timestep < PoAStartTimestep) || any(x$Timestep > Ntimesteps))
    .INApestPoAStop(
      "ObservationHistory timesteps must be between PoAStartTimestep (",
      PoAStartTimestep, ") and Ntimesteps (", Ntimesteps, ")"
    )
  if(anyDuplicated(x$Timestep))
    .INApestPoAStop("ObservationHistory must contain at most one row per timestep")

  for(Name in c("BackgroundDetections", "InfoTriggeredDetections"))
    if(any(!is.na(x[[Name]]) & (!is.finite(x[[Name]]) | x[[Name]] < 0)))
      .INApestPoAStop(Name, " must contain non-negative finite counts or NA")

  x <- x[order(x$Timestep), , drop=FALSE]
  x$Round <- seq_len(nrow(x))
  out <- x[,c("Round","Timestep","BackgroundDetections","InfoTriggeredDetections"),drop=FALSE]
  ### Backwards-compatible read-only aliases in returned results.
  out$SurveillanceDetections <- out$BackgroundDetections
  out$InformationDetections <- out$InfoTriggeredDetections
  out
}

.INApestPoAReweightPrior <- function(Weights, Absent, PriorPoA)
{
  PriorPoA <- .INApestPoAClip01(PriorPoA, "PriorPoA")
  Present <- !Absent

  if(PriorPoA > 0 && !any(Absent))
    .INApestPoAStop("PriorPoA > 0 cannot be imposed because no simulated particle is absent at PoAStartTimestep")
  if(PriorPoA < 1 && !any(Present))
    .INApestPoAStop("PriorPoA < 1 cannot be imposed because no simulated particle is present at PoAStartTimestep")

  Weights[] <- 0
  if(any(Absent) && PriorPoA > 0)
    Weights[Absent] <- PriorPoA / sum(Absent)
  if(any(Present) && PriorPoA < 1)
    Weights[Present] <- (1 - PriorPoA) / sum(Present)
  Weights / sum(Weights)
}


.INApestPoAWeightedMean <- function(x, w)
{
  if(!length(x) || !length(w) || sum(w) <= 0) return(NA_real_)
  sum(as.numeric(x) * w) / sum(w)
}


.INApestPoAConditionalDetection <- function(DetectionCounts, Present, Weights)
{
  Den <- sum(Weights[Present])
  if(Den <= 0) return(NA_real_)
  sum(Weights[Present] * (DetectionCounts[Present] > 0)) / Den
}


.INApestPoAConditionalSSeLikelihood <- function(NoDetectionProbability, Present, Weights)
{
  if(is.null(NoDetectionProbability)) return(NA_real_)
  Den <- sum(Weights[Present])
  if(Den <= 0) return(NA_real_)
  q <- NoDetectionProbability[Present]
  w <- Weights[Present]
  ok <- !is.na(q) & is.finite(q)
  if(!any(ok) || sum(w[ok]) <= 0) return(NA_real_)
  1 - sum(w[ok] * q[ok]) / sum(w[ok])
}


.INApestPoAObservationLikelihood <- function(NoDetectionProbability, Observed, Match, SourceName)
{
  if(is.na(Observed)) return(rep(1, length(NoDetectionProbability)))
  q <- as.numeric(NoDetectionProbability)
  if(any(is.na(q) | !is.finite(q) | q < 0 | q > 1))
    .INApestPoAStop(
      "Likelihood mode requires finite ", SourceName,
      " no-detection probabilities in [0,1] for every particle at an observed timestep"
    )
  if(Observed == 0) return(q)
  if(Match == "binary" && Observed > 0) return(1 - q)
  .INApestPoAStop(
    "Likelihood mode currently supports zero detections, or positive detections only with ObservationMatch='binary'. ",
    "Exact positive detection counts require a count-likelihood output from the underlying model."
  )
}


.INApestPoARescaleMatchedForPropagation <- function(PriorWeights, Match, Absent,
                                                     PosteriorWeights, Timestep)
{
  Raw <- PriorWeights * as.numeric(Match)
  TargetAbsent <- sum(PosteriorWeights[Absent])
  TargetPresent <- sum(PosteriorWeights[!Absent])
  Out <- rep(0, length(PriorWeights))

  for(group in list(list(idx = Absent, target = TargetAbsent, name = "absent"),
                    list(idx = !Absent, target = TargetPresent, name = "present")))
  {
    idx <- group$idx & Match
    if(group$target > 0)
    {
      den <- sum(Raw[idx])
      if(den <= 0)
        .INApestPoAStop(
          "Likelihood update retained positive posterior mass for the ", group$name,
          " class at timestep ", Timestep,
          " but no realised matching particles are available to propagate that class. ",
          "Increase Nperm or shorten the conditioning horizon."
        )
      Out[idx] <- Raw[idx] / den * group$target
    }
  }
  if(sum(Out) <= 0) .INApestPoAStop("No matched particles available for likelihood propagation")
  Out / sum(Out)
}


.INApestPoAMatchObservation <- function(Simulated, Observed, Match)
{
  if(is.na(Observed)) return(rep(TRUE, length(Simulated)))
  if(Match == "exact") return(Simulated == Observed)
  (Simulated > 0) == (Observed > 0)
}


.INApestPoAFreedomVector <- function(Adapter, Timestep, FreedomFunction)
{
  if(is.null(FreedomFunction)) return(NULL)
  z <- FreedomFunction(Adapter, Timestep)
  if(length(z) != Adapter$Nperm || any(is.na(z)))
    .INApestPoAStop("FreedomFunction must return one non-missing logical value per permutation")
  as.logical(z)
}


.INApestPoAResidualStageSummary <- function(Adapter, Timestep, Weights)
{
  x <- Adapter$StageResidualPopulation
  if(is.null(x)) return(NULL)

  Nstage <- dim(x)[3L]
  out <- data.frame(
    Stage = seq_len(Nstage),
    ExpectedResidualPopulation = NA_real_,
    ProbabilityStagePresent = NA_real_,
    stringsAsFactors = FALSE
  )
  for(s in seq_len(Nstage))
  {
    v <- x[, Timestep, s]
    out$ExpectedResidualPopulation[s] <- .INApestPoAWeightedMean(v, Weights)
    out$ProbabilityStagePresent[s] <- .INApestPoAWeightedMean(v > 0, Weights)
  }
  out
}


.INApestPoAResidualSpatialRisk <- function(Adapter, Timestep, Weights)
{
  if(!is.null(Adapter$UnitPresence))
  {
    UnitPresence <- Adapter$UnitPresence[, Timestep, , drop = FALSE]
    Nunit <- dim(UnitPresence)[1L]
    Risk <- numeric(Nunit)
    for(u in seq_len(Nunit))
      Risk[u] <- .INApestPoAWeightedMean(UnitPresence[u, 1L, ], Weights)

    out <- if(!is.null(Adapter$UnitLabels)) Adapter$UnitLabels else data.frame(Unit = seq_len(Nunit))
    out$PosteriorPresenceProbability <- Risk
    return(out)
  }

  if(is.data.frame(Adapter$PointHistory) && nrow(Adapter$PointHistory) &&
     all(c("perm", "timestep", "x", "y") %in% names(Adapter$PointHistory)))
  {
    out <- Adapter$PointHistory[Adapter$PointHistory$timestep == Timestep, , drop = FALSE]
    if(!nrow(out))
    {
      out$ParticleWeight <- numeric(0)
      return(out)
    }
    out$ParticleWeight <- Weights[as.integer(out$perm)]
    return(out)
  }

  NULL
}


###############################################################################
### Common proof-of-absence inference core
###############################################################################

INApestPoACore <- function(
  ModelResults,
  ModelAdapter,
  PoAStartTimestep = 1,
  ObservationHistory = NULL,
  PriorPoA = NA,
  EvidenceSources = c("Background", "InfoTriggered"),
  ObservationModel = NULL,
  ObservationMatch = c("exact", "binary"),
  PoAMethod = c("simulation", "likelihood"),
  TargetPoA = c(0.95, 0.99),
  AbsenceFunction = NULL,
  FreedomFunction = NULL,
  ReturnParticles = FALSE,
  KeepModelResults = FALSE,
  ResultClass = "INApestPoA"
)
{
  if(!is.function(ModelAdapter))
    .INApestPoAStop("ModelAdapter must be an adapter function")

  Adapter <- ModelAdapter(ModelResults)
  Nperm <- Adapter$Nperm
  Ntimesteps <- Adapter$Ntimesteps

  PoAStartTimestep <- as.integer(PoAStartTimestep)
  if(length(PoAStartTimestep) != 1L || is.na(PoAStartTimestep) ||
     PoAStartTimestep < 1L || PoAStartTimestep > Ntimesteps)
    .INApestPoAStop("PoAStartTimestep must be between 1 and ", Ntimesteps)

  EvidenceSources <- .INApestPoANormaliseEvidenceSources(EvidenceSources)
  ObservationMatch <- match.arg(ObservationMatch)
  PoAMethod <- match.arg(PoAMethod)

  UserObservationModelSupplied <- !is.null(ObservationModel)
  ObservationModelMode <- if(UserObservationModelSupplied) "user_supplied" else "none"
  if(PoAMethod == "likelihood" && is.null(ObservationModel) && length(EvidenceSources))
  {
    Recorded <- vapply(EvidenceSources, function(src)
      .INApestPoARecordedDetectionAvailable(Adapter, src), logical(1))
    Stored <- vapply(EvidenceSources, function(src) {
      q <- if(src == "Background") Adapter$BackgroundNoDetectionProbability else Adapter$InfoTriggeredNoDetectionProbability
      !is.null(q)
    }, logical(1))
    Missing <- !(Recorded | Stored)
    if(any(Missing))
      .INApestPoAStop(
        "Likelihood mode lacks an observation model for: ",
        paste(EvidenceSources[Missing], collapse=", "), ". ",
        "Use patched engine results with recorded detection probabilities, supply ObservationModel, or supply legacy no-detection probabilities."
      )

    if(all(Recorded))
    {
      ObservationModel <- INApestPoARecordedDetectionModel()
      ObservationModelMode <- "recorded_engine_detection_probabilities"
    }
    else if(all(Stored) && !any(Recorded))
      ObservationModelMode <- "stored_no_detection_probabilities"
    else
    {
      RecordedModel <- INApestPoARecordedDetectionModel()
      RecordedBySource <- setNames(Recorded, EvidenceSources)
      ObservationModel <- function(Adapter, ModelResults, Timestep, Source) {
        if(isTRUE(RecordedBySource[[Source]]))
          return(RecordedModel(Adapter, ModelResults, Timestep, Source))
        q <- if(Source == "Background") Adapter$BackgroundNoDetectionProbability else Adapter$InfoTriggeredNoDetectionProbability
        q[,Timestep]
      }
      ObservationModelMode <- "hybrid_recorded_and_stored"
    }
  }
  if(PoAMethod == "simulation") ObservationModelMode <- "realised_detection_events"

  if(!is.logical(ReturnParticles) || length(ReturnParticles) != 1L || is.na(ReturnParticles))
    .INApestPoAStop("ReturnParticles must be TRUE or FALSE")
  if(!is.logical(KeepModelResults) || length(KeepModelResults) != 1L || is.na(KeepModelResults))
    .INApestPoAStop("KeepModelResults must be TRUE or FALSE")

  TargetPoA <- sort(unique(as.numeric(TargetPoA)))
  if(length(TargetPoA) && any(!is.finite(TargetPoA) | TargetPoA < 0 | TargetPoA > 1))
    .INApestPoAStop("TargetPoA values must be between 0 and 1")

  if(!is.null(AbsenceFunction))
  {
    if(!is.function(AbsenceFunction)) .INApestPoAStop("AbsenceFunction must be NULL or a function")
    Absence <- matrix(FALSE, nrow = Nperm, ncol = Ntimesteps)
    for(t in seq_len(Ntimesteps))
    {
      z <- AbsenceFunction(Adapter, t)
      if(length(z) != Nperm || any(is.na(z)))
        .INApestPoAStop("AbsenceFunction must return one non-missing logical value per permutation")
      Absence[, t] <- as.logical(z)
    }
  }
  else
    Absence <- Adapter$Absence

  Obs <- .INApestPoANormaliseObservationHistory(ObservationHistory, PoAStartTimestep, Ntimesteps)

  ### Current-round likelihood weighting can use an observation model alone.
  ### Sequential propagation additionally needs realised source-specific events
  ### because detections can change future information and management.
  if(nrow(Obs) > 1L)
  {
    if("Background" %in% EvidenceSources && Adapter$BackgroundEventContract != "explicit")
      .INApestPoAStop("Sequential PoA requires explicit BackgroundDetectedResults/new background detection events; legacy DetectedResults is not sufficient")
    if("InfoTriggered" %in% EvidenceSources && Adapter$InfoTriggeredEventContract != "explicit")
      .INApestPoAStop("Sequential PoA with InfoTriggered evidence requires explicit InfoTriggeredDetectedResults/new info-triggered detection events")
  }

  Weights <- rep(1 / Nperm, Nperm)
  EmpiricalPriorPoA <- .INApestPoAWeightedMean(Absence[, PoAStartTimestep], Weights)
  if(length(PriorPoA) == 1L && !is.na(PriorPoA))
    Weights <- .INApestPoAReweightPrior(Weights, Absence[, PoAStartTimestep], PriorPoA)

  AppliedPriorPoA <- .INApestPoAWeightedMean(Absence[, PoAStartTimestep], Weights)

  Rows <- vector("list", nrow(Obs) + 1L)
  Rows[[1L]] <- data.frame(
    Round = 0L,
    Timestep = PoAStartTimestep,
    PriorPoA = AppliedPriorPoA,
    BackgroundSSe = NA_real_,
    InfoTriggeredSSe = NA_real_,
    CombinedSSe = NA_real_,
    SurveillanceSSe = NA_real_,
    InformationSSe = NA_real_,
    ObservationProbability = NA_real_,
    PosteriorPoA = AppliedPriorPoA,
    PriorExpectedResidualPopulation = .INApestPoAWeightedMean(Adapter$ResidualPopulation[, PoAStartTimestep], Weights),
    PosteriorExpectedResidualPopulation = .INApestPoAWeightedMean(Adapter$ResidualPopulation[, PoAStartTimestep], Weights),
    PosteriorFunctionalFreedom = if(is.null(FreedomFunction)) NA_real_ else .INApestPoAWeightedMean(.INApestPoAFreedomVector(Adapter, PoAStartTimestep, FreedomFunction), Weights),
    EffectiveParticles = 1 / sum(Weights^2),
    MatchingParticles = Nperm,
    stringsAsFactors = FALSE
  )

  WeightHistory <- if(ReturnParticles) matrix(NA_real_, nrow = Nperm, ncol = nrow(Obs) + 1L) else NULL
  if(ReturnParticles) WeightHistory[, 1L] <- Weights

  DetectionRows <- vector("list", nrow(Obs))

  if(nrow(Obs))
  {
    for(r in seq_len(nrow(Obs)))
    {
      t <- Obs$Timestep[r]
      PriorWeights <- Weights
      Present <- !Absence[, t]

      BackgroundCount <- Adapter$BackgroundDetections[, t]
      InfoTriggeredCount <- Adapter$InfoTriggeredDetections[, t]
      CombinedDetected <- rep(FALSE, Nperm)
      if("Background" %in% EvidenceSources) CombinedDetected <- CombinedDetected | BackgroundCount > 0
      if("InfoTriggered" %in% EvidenceSources) CombinedDetected <- CombinedDetected | InfoTriggeredCount > 0

      BackgroundQ <- .INApestPoAResolveNoDetection(Adapter, ModelResults, ObservationModel, t, "Background")
      InfoTriggeredQ <- .INApestPoAResolveNoDetection(Adapter, ModelResults, ObservationModel, t, "InfoTriggered")
      if(PoAMethod == "likelihood")
      {
        if("Background" %in% EvidenceSources && is.null(BackgroundQ))
          .INApestPoAStop("No Background observation likelihood is available at timestep ", t)
        if("InfoTriggered" %in% EvidenceSources && is.null(InfoTriggeredQ))
          .INApestPoAStop("No InfoTriggered observation likelihood is available at timestep ", t)
      }

      if(PoAMethod == "likelihood")
      {
        BackgroundSSe <- .INApestPoAConditionalSSeLikelihood(BackgroundQ, Present, PriorWeights)
        InfoTriggeredSSe <- .INApestPoAConditionalSSeLikelihood(InfoTriggeredQ, Present, PriorWeights)
        CombinedQ <- rep(1, Nperm)
        if("Background" %in% EvidenceSources) CombinedQ <- CombinedQ * BackgroundQ
        if("InfoTriggered" %in% EvidenceSources) CombinedQ <- CombinedQ * InfoTriggeredQ
        CombinedSSe <- .INApestPoAConditionalSSeLikelihood(CombinedQ, Present, PriorWeights)
      }
      else
      {
        BackgroundSSe <- .INApestPoAConditionalDetection(BackgroundCount, Present, PriorWeights)
        InfoTriggeredSSe <- .INApestPoAConditionalDetection(InfoTriggeredCount, Present, PriorWeights)
        CombinedSSe <- .INApestPoAConditionalDetection(as.numeric(CombinedDetected), Present, PriorWeights)
      }

      PriorPoAThisRound <- .INApestPoAWeightedMean(Absence[, t], PriorWeights)
      PriorExpectedResidual <- .INApestPoAWeightedMean(Adapter$ResidualPopulation[, t], PriorWeights)

      Match <- rep(TRUE, Nperm)
      if("Background" %in% EvidenceSources)
        Match <- Match & .INApestPoAMatchObservation(BackgroundCount, Obs$BackgroundDetections[r], ObservationMatch)
      if("InfoTriggered" %in% EvidenceSources)
        Match <- Match & .INApestPoAMatchObservation(InfoTriggeredCount, Obs$InfoTriggeredDetections[r], ObservationMatch)

      if(PoAMethod == "likelihood")
      {
        Likelihood <- rep(1, Nperm)
        if("Background" %in% EvidenceSources)
          Likelihood <- Likelihood * .INApestPoAObservationLikelihood(
            BackgroundQ, Obs$BackgroundDetections[r], ObservationMatch, "background surveillance"
          )
        if("InfoTriggered" %in% EvidenceSources)
          Likelihood <- Likelihood * .INApestPoAObservationLikelihood(
            InfoTriggeredQ, Obs$InfoTriggeredDetections[r], ObservationMatch, "information-triggered surveillance"
          )
        ObservationProbability <- sum(PriorWeights * Likelihood)
        if(!is.finite(ObservationProbability) || ObservationProbability <= 0)
          .INApestPoAStop("Observed evidence has zero likelihood at timestep ", t)
        PosteriorWeights <- PriorWeights * Likelihood
        PosteriorWeights <- PosteriorWeights / sum(PosteriorWeights)
      }
      else
      {
        ObservationProbability <- sum(PriorWeights[Match])
        if(!is.finite(ObservationProbability) || ObservationProbability <= 0)
          .INApestPoAStop(
            "No simulated particles match ObservationHistory at timestep ", t,
            ". Increase Nperm, use ObservationMatch='binary' where appropriate, or review the model/observation assumptions."
          )
        PosteriorWeights <- PriorWeights * as.numeric(Match)
        PosteriorWeights <- PosteriorWeights / sum(PosteriorWeights)
      }

      Functional <- .INApestPoAFreedomVector(Adapter, t, FreedomFunction)
      PosteriorPoA <- .INApestPoAWeightedMean(Absence[, t], PosteriorWeights)
      PosteriorExpectedResidual <- .INApestPoAWeightedMean(Adapter$ResidualPopulation[, t], PosteriorWeights)

      Rows[[r + 1L]] <- data.frame(
        Round = Obs$Round[r],
        Timestep = t,
        PriorPoA = PriorPoAThisRound,
        BackgroundSSe = BackgroundSSe,
        InfoTriggeredSSe = InfoTriggeredSSe,
        CombinedSSe = CombinedSSe,
        ### Deprecated output aliases.
        SurveillanceSSe = BackgroundSSe,
        InformationSSe = InfoTriggeredSSe,
        ObservationProbability = ObservationProbability,
        PosteriorPoA = PosteriorPoA,
        PriorExpectedResidualPopulation = PriorExpectedResidual,
        PosteriorExpectedResidualPopulation = PosteriorExpectedResidual,
        PosteriorFunctionalFreedom = if(is.null(Functional)) NA_real_ else .INApestPoAWeightedMean(Functional, PosteriorWeights),
        EffectiveParticles = 1 / sum(PosteriorWeights^2),
        MatchingParticles = sum(Match),
        stringsAsFactors = FALSE
      )

      DetectionRows[[r]] <- data.frame(
        Round = Obs$Round[r],
        Timestep = t,
        ObservedBackgroundDetections = Obs$BackgroundDetections[r],
        ObservedInfoTriggeredDetections = Obs$InfoTriggeredDetections[r],
        ExpectedBackgroundDetectionsGivenPrior = .INApestPoAWeightedMean(BackgroundCount, PriorWeights),
        ExpectedInfoTriggeredDetectionsGivenPrior = .INApestPoAWeightedMean(InfoTriggeredCount, PriorWeights),
        ProbabilityAnyBackgroundDetectionGivenPresence = BackgroundSSe,
        ProbabilityAnyInfoTriggeredDetectionGivenPresence = InfoTriggeredSSe,
        ProbabilityAnySelectedDetectionGivenPresence = CombinedSSe,
        stringsAsFactors = FALSE
      )

      if(ReturnParticles) WeightHistory[, r + 1L] <- PosteriorWeights

      ### Direct likelihood weighting gives a lower-variance current PoA, but
      ### detections alter future INApest management. Therefore future state
      ### propagation is restricted to trajectories whose realised evidence
      ### matched the observations, with absent/present class masses rescaled
      ### to the likelihood-updated posterior. This is a Rao-Blackwellized
      ### particle filter rather than an inconsistent post-hoc reweighting of
      ### trajectories with the wrong management history.
      if(r < nrow(Obs))
      {
        if(PoAMethod == "likelihood")
          Weights <- .INApestPoARescaleMatchedForPropagation(
            PriorWeights, Match, Absence[, t], PosteriorWeights, t
          )
        else
          Weights <- PosteriorWeights
      }
      else
        Weights <- PosteriorWeights
    }
  }

  PoASummary <- do.call(rbind, Rows)
  DetectionSummary <- if(length(DetectionRows)) do.call(rbind, DetectionRows) else data.frame()

  if(nrow(DetectionSummary))
  {
    DetectionSummary$ObservedSurveillanceDetections <- DetectionSummary$ObservedBackgroundDetections
    DetectionSummary$ObservedInformationDetections <- DetectionSummary$ObservedInfoTriggeredDetections
    DetectionSummary$ExpectedSurveillanceDetectionsGivenPrior <- DetectionSummary$ExpectedBackgroundDetectionsGivenPrior
    DetectionSummary$ExpectedInformationDetectionsGivenPrior <- DetectionSummary$ExpectedInfoTriggeredDetectionsGivenPrior
    DetectionSummary$ProbabilityAnySurveillanceDetectionGivenPresence <- DetectionSummary$ProbabilityAnyBackgroundDetectionGivenPresence
    DetectionSummary$ProbabilityAnyInformationDetectionGivenPresence <- DetectionSummary$ProbabilityAnyInfoTriggeredDetectionGivenPresence
  }

  StopRows <- lapply(TargetPoA, function(Target)
  {
    Hit <- which(PoASummary$PosteriorPoA >= Target)
    data.frame(
      TargetPoA = Target,
      Reached = length(Hit) > 0L,
      FirstRound = if(length(Hit)) PoASummary$Round[Hit[1L]] else NA_integer_,
      FirstTimestep = if(length(Hit)) PoASummary$Timestep[Hit[1L]] else NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  StoppingSummary <- if(length(StopRows)) do.call(rbind, StopRows) else data.frame()

  FinalTimestep <- if(nrow(Obs)) tail(Obs$Timestep, 1L) else PoAStartTimestep
  ResidualStageSummary <- .INApestPoAResidualStageSummary(Adapter, FinalTimestep, Weights)
  ResidualSpatialRisk <- .INApestPoAResidualSpatialRisk(Adapter, FinalTimestep, Weights)

  SurveillanceSystemSensitivity <- if(nrow(PoASummary) > 1L)
    PoASummary[-1L, c("Round", "Timestep", "BackgroundSSe", "InfoTriggeredSSe", "CombinedSSe"), drop = FALSE]
  else
    data.frame()

  out <- list(
    ModelName = Adapter$ModelName,
    AdapterName = Adapter$AdapterName,
    PoAStartTimestep = PoAStartTimestep,
    EmpiricalPriorPoA = EmpiricalPriorPoA,
    AppliedPriorPoA = AppliedPriorPoA,
    PoASummary = PoASummary,
    SurveillanceSystemSensitivity = SurveillanceSystemSensitivity,
    DetectionSummary = DetectionSummary,
    ResidualStageSummary = ResidualStageSummary,
    ResidualSpatialRisk = ResidualSpatialRisk,
    StoppingSummary = StoppingSummary,
    ObservationHistory = Obs,
    FinalParticleWeights = Weights,
    ParticleWeightHistory = WeightHistory,
    Settings = list(
      EvidenceSources = EvidenceSources,
      ObservationModelSupplied = UserObservationModelSupplied,
      ObservationModelMode = ObservationModelMode,
      BackgroundEventContract = Adapter$BackgroundEventContract,
      InfoTriggeredEventContract = Adapter$InfoTriggeredEventContract,
      ObservationMatch = ObservationMatch,
      PoAMethod = PoAMethod,
      TargetPoA = TargetPoA,
      PriorPoA = PriorPoA,
      ReturnParticles = ReturnParticles
    )
  )

  if(KeepModelResults) out$ModelResults <- ModelResults

  class(out) <- unique(c(ResultClass, "INApestPoA", "list"))
  out
}


###############################################################################
### Common wrapper runner
###############################################################################

.INApestPoARun <- function(
  ModelFunction,
  ModelAdapter,
  ModelArgs = list(),
  ModelResults = NULL,
  PoAStartTimestep = 1,
  ObservationHistory = NULL,
  PriorPoA = NA,
  EvidenceSources = c("Background", "InfoTriggered"),
  ObservationModel = NULL,
  ObservationMatch = c("exact", "binary"),
  PoAMethod = c("simulation", "likelihood"),
  TargetPoA = c(0.95, 0.99),
  AbsenceFunction = NULL,
  FreedomFunction = NULL,
  ReturnParticles = FALSE,
  KeepModelResults = FALSE,
  ResultClass = "INApestPoA"
)
{
  if(is.null(ModelResults))
  {
    if(!is.list(ModelArgs)) .INApestPoAStop("ModelArgs must be a named list")
    if(!exists(ModelFunction, mode = "function", inherits = TRUE))
      .INApestPoAStop("Underlying model function '", ModelFunction, "' is not loaded. Source its INApest file first.")

    Fun <- get(ModelFunction, mode = "function", inherits = TRUE)
    FormalNames <- names(formals(Fun))

    ### Request in-memory node output and suppress optional file/plot/progress
    ### side effects unless the caller explicitly asks for them.
    if("ReturnResults" %in% FormalNames) ModelArgs$ReturnResults <- TRUE
    if("SaveResults" %in% FormalNames && is.null(ModelArgs$SaveResults)) ModelArgs$SaveResults <- FALSE
    if("DoPlots" %in% FormalNames && is.null(ModelArgs$DoPlots)) ModelArgs$DoPlots <- FALSE
    if("DoProgress" %in% FormalNames && is.null(ModelArgs$DoProgress)) ModelArgs$DoProgress <- FALSE

    ModelResults <- do.call(Fun, ModelArgs)
  }

  INApestPoACore(
    ModelResults = ModelResults,
    ModelAdapter = ModelAdapter,
    PoAStartTimestep = PoAStartTimestep,
    ObservationHistory = ObservationHistory,
    PriorPoA = PriorPoA,
    EvidenceSources = EvidenceSources,
    ObservationModel = ObservationModel,
    ObservationMatch = ObservationMatch,
    PoAMethod = PoAMethod,
    TargetPoA = TargetPoA,
    AbsenceFunction = AbsenceFunction,
    FreedomFunction = FreedomFunction,
    ReturnParticles = ReturnParticles,
    KeepModelResults = KeepModelResults,
    ResultClass = ResultClass
  )
}


###############################################################################
### Thin public wrappers
###############################################################################

INApestPoA <- function(ModelArgs = list(), ModelResults = NULL, ...)
  .INApestPoARun("INApest", .INApestPoAAdapterBinaryNode, ModelArgs, ModelResults,
                 ..., ResultClass = "INApestPoAResults")

INApestParallelPoA <- function(ModelArgs = list(), ModelResults = NULL, ...)
  .INApestPoARun("INApestParallel", .INApestPoAAdapterBinaryNode, ModelArgs, ModelResults,
                 ..., ResultClass = "INApestParallelPoA")

INApestMetaPoA <- function(ModelArgs = list(), ModelResults = NULL, ...)
  .INApestPoARun("INApestMeta", .INApestPoAAdapterPopulationNode, ModelArgs, ModelResults,
                 ..., ResultClass = "INApestMetaPoA")

INApestMetaParallelPoA <- function(ModelArgs = list(), ModelResults = NULL, ...)
  .INApestPoARun("INApestMetaParallel", .INApestPoAAdapterPopulationNode, ModelArgs, ModelResults,
                 ..., ResultClass = "INApestMetaParallelPoA")

INApestMetaMultipleLandUsePoA <- function(ModelArgs = list(), ModelResults = NULL, ...)
  .INApestPoARun("INApestMetaMultipleLandUse", .INApestPoAAdapterMultipleLandUse, ModelArgs, ModelResults,
                 ..., ResultClass = "INApestMetaMultipleLandUsePoA")
INApestMetaParallelMultipleLandUsePoA <- function(ModelArgs = list(), ModelResults = NULL, ...)
  .INApestPoARun("INApestMetaParallelMultipleLandUse", .INApestPoAAdapterMultipleLandUse, ModelArgs, ModelResults,
                 ..., ResultClass = "INApestMetaParallelMultipleLandUsePoA")

INApestMetaTransitionMatrixPoA <- function(ModelArgs = list(), ModelResults = NULL, ...)
  .INApestPoARun("INApestMetaTransitionMatrix", .INApestPoAAdapterTransitionMatrix, ModelArgs, ModelResults,
                 ..., ResultClass = "INApestMetaTransitionMatrixPoA")

INApestMetaTransitionMatrixParallelPoA <- function(ModelArgs = list(), ModelResults = NULL, ...)
  .INApestPoARun("INApestMetaTransitionMatrixParallel", .INApestPoAAdapterTransitionMatrix, ModelArgs, ModelResults,
                 ..., ResultClass = "INApestMetaTransitionMatrixParallelPoA")

INApestMetaPointPoA <- function(ModelArgs = list(), ModelResults = NULL, ...)
  .INApestPoARun("INApestMetaPoint", .INApestPoAAdapterPoint, ModelArgs, ModelResults,
                 ..., ResultClass = "INApestMetaPointPoA")

INApestMetaPointParallelPoA <- function(ModelArgs = list(), ModelResults = NULL, ...)
  .INApestPoARun("INApestMetaPointParallel", .INApestPoAAdapterPoint, ModelArgs, ModelResults,
                 ..., ResultClass = "INApestMetaPointParallelPoA")

INApestPointTransitionMatrixPoA <- function(ModelArgs = list(), ModelResults = NULL, ...)
  .INApestPoARun("INApestPointTransitionMatrix", .INApestPoAAdapterPointTransitionMatrix, ModelArgs, ModelResults,
                 ..., ResultClass = "INApestPointTransitionMatrixPoA")

INApestPointTransitionMatrixParallelPoA <- function(ModelArgs = list(), ModelResults = NULL, ...)
  .INApestPoARun("INApestPointTransitionMatrixParallel", .INApestPoAAdapterPointTransitionMatrix, ModelArgs, ModelResults,
                 ..., ResultClass = "INApestPointTransitionMatrixParallelPoA")


###############################################################################
### Print method
###############################################################################

print.INApestPoA <- function(x, ...)
{
  cat("INApest proof-of-absence result\n")
  if(!is.na(x$ModelName)) cat("Model:", x$ModelName, "\n")
  cat("Adapter:", x$AdapterName, "\n")
  cat("Empirical prior PoA:", format(x$EmpiricalPriorPoA, digits = 4), "\n")
  cat("Applied prior PoA:", format(x$AppliedPriorPoA, digits = 4), "\n")
  Last <- x$PoASummary[nrow(x$PoASummary), , drop = FALSE]
  cat("Latest posterior PoA:", format(Last$PosteriorPoA, digits = 4),
      "at timestep", Last$Timestep, "\n")
  invisible(x)
}
