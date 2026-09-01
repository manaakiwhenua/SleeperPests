###############################################################################
### INApest vertebrate proof-of-absence adapters and wrappers
###
### Additive companion to the frozen INApestPoA.R core.  This file does not
### modify the validated Bayesian PoA machinery.  It adapts the two vertebrate
### biological engines and incorporates pre-deployed vertebrate-control
### detections as routine/background observation opportunities.
###############################################################################

.INApestPoAVertebrateRequireCore <- function()
{
  needed <- c(
    ".INApestPoARun", ".INApestPoAAdapterTransitionMatrix",
    ".INApestPoAAdapterPointTransitionMatrix", "INApestPoARecordedDetectionModel"
  )
  missing <- needed[!vapply(needed, exists, logical(1), mode="function", inherits=TRUE)]
  if(length(missing))
    stop("Source INApestPoA.R before INApestPoAVertebrate.R. Missing: ",
         paste(missing, collapse=", "), call.=FALSE)
  invisible(TRUE)
}


### Node-engine routine-control detection is drawn independently of the lethal
### control draw.  Conditional on the recorded pre-control abundance, the
### probability of no routine-control detection is therefore simply
### product((1-p)^N) across node x stage observation units.
.INApestPoAVertebrateNodeControlQ <- function(ModelResults, Timestep, Nperm)
{
  A <- ModelResults$RoutineControlObservationAbundanceResults
  P <- ModelResults$RoutineControlDetectionProbabilityResults
  if(is.null(A) || is.null(P)) return(rep(1, Nperm))
  if(length(dim(A)) != 4L || !identical(dim(A), dim(P)))
    stop("Vertebrate node routine-control abundance and detection probability must be matching node x stage x timestep x permutation arrays", call.=FALSE)
  if(Timestep < 1L || Timestep > dim(A)[3L] || Nperm != dim(A)[4L])
    stop("Invalid vertebrate node routine-control likelihood dimensions", call.=FALSE)
  if(any(!is.finite(A)) || any(A < 0) ||
     any(!is.finite(P)) || any(P < 0 | P > 1))
    stop("Invalid vertebrate node routine-control observation values", call.=FALSE)

  out <- rep(1, Nperm)
  for(pp in seq_len(Nperm))
    out[pp] <- prod((1 - P[,,Timestep,pp]) ^ A[,,Timestep,pp])
  out
}


### The point engine deliberately uses one shared U(0,1) encounter draw for
### routine-control detection and killing.  The particle already contains the
### realised control-death outcome, so the observation likelihood must condition
### on that outcome rather than multiply independent detection and death terms.
###
### For detection probability d and kill probability k:
###   dead:     P(no detect | dead)     = max(k-d,0) / k
###   survived: P(no detect | survived) = [1-max(k,d)] / (1-k)
###
### Impossible recorded outcomes (dead with k=0; survived with k=1) are rejected.
.INApestPoAVertebratePointControlQ <- function(ModelResults, Timestep, Nperm)
{
  H <- ModelResults$ControlObservationHistory
  out <- rep(1, Nperm)
  if(is.null(H) || !is.data.frame(H) || !nrow(H)) return(out)

  needed <- c("perm","timestep","control_detection_prob","control_kill_prob","control_dead")
  missing <- setdiff(needed, names(H))
  if(length(missing))
    stop("ControlObservationHistory is missing: ", paste(missing, collapse=", "), call.=FALSE)

  for(pp in seq_len(Nperm))
  {
    z <- H[as.integer(H$perm) == pp & as.integer(H$timestep) == Timestep, , drop=FALSE]
    if(!nrow(z)) next
    d <- as.numeric(z$control_detection_prob)
    k <- as.numeric(z$control_kill_prob)
    dead <- as.logical(z$control_dead)
    if(any(is.na(d) | !is.finite(d) | d < 0 | d > 1) ||
       any(is.na(k) | !is.finite(k) | k < 0 | k > 1) || any(is.na(dead)))
      stop("Invalid values in ControlObservationHistory", call.=FALSE)

    qi <- numeric(nrow(z))
    id <- which(dead)
    if(length(id))
    {
      if(any(k[id] <= 0))
        stop("ControlObservationHistory contains a death with zero kill probability", call.=FALSE)
      qi[id] <- pmax(k[id] - d[id], 0) / k[id]
    }
    isurv <- which(!dead)
    if(length(isurv))
    {
      if(any(k[isurv] >= 1))
        stop("ControlObservationHistory contains a survivor with unit kill probability", call.=FALSE)
      qi[isurv] <- (1 - pmax(k[isurv], d[isurv])) / (1 - k[isurv])
    }
    if(any(!is.finite(qi) | qi < -1e-12 | qi > 1 + 1e-12))
      stop("Invalid conditional no-detection probability for vertebrate point control", call.=FALSE)
    qi <- pmin(1,pmax(0,qi))
    out[pp] <- prod(qi)
  }
  out
}


### Default vertebrate observation model.
###
### Background = ordinary INApest background surveillance × pre-deployed
### vertebrate-control/device observation.  InfoTriggered is delegated unchanged
### to the already validated recorded-probability PoA observation model.
INApestPoAVertebrateObservationModel <- function()
{
  .INApestPoAVertebrateRequireCore()
  ordinary <- INApestPoARecordedDetectionModel()

  function(Adapter, ModelResults, Timestep, Source)
  {
    q <- ordinary(Adapter=Adapter, ModelResults=ModelResults,
                  Timestep=Timestep, Source=Source)
    if(!identical(Source, "Background")) return(q)

    if(identical(Adapter$AdapterName, ".INApestPoAAdapterVertebrateNode"))
      q <- q * .INApestPoAVertebrateNodeControlQ(ModelResults, Timestep, Adapter$Nperm)
    else if(identical(Adapter$AdapterName, ".INApestPoAAdapterVertebratePoint"))
      q <- q * .INApestPoAVertebratePointControlQ(ModelResults, Timestep, Adapter$Nperm)
    q
  }
}


### Vertebrate node retains the transition-matrix biological state contract.
### The patched engine supplies the same canonical state/surveillance arrays as
### INApestMetaTransitionMatrix plus the pre-control observation phase.
.INApestPoAAdapterVertebrateNode <- function(ModelResults)
{
  .INApestPoAVertebrateRequireCore()
  Base <- .INApestPoAAdapterTransitionMatrix(ModelResults)
  Base$AdapterName <- ".INApestPoAAdapterVertebrateNode"
  Base
}


### Vertebrate point retains the PointTransitionMatrix state contract.  Add the
### routine-control detection count to the Background event stream.  The
### ordinary point adapter continues to provide point/stage residual summaries.
.INApestPoAAdapterVertebratePoint <- function(ModelResults)
{
  .INApestPoAVertebrateRequireCore()
  Base <- .INApestPoAAdapterPointTransitionMatrix(ModelResults)
  S <- ModelResults$Summary
  if(!is.data.frame(S) || !all(c("perm","timestep") %in% names(S)))
    stop("Vertebrate point Summary must contain perm and timestep", call.=FALSE)

  control <- if("n_control_detections" %in% names(S))
    .INApestPoAPointMatrix(S, "n_control_detections", Base$Nperm, Base$Ntimesteps) else
    matrix(0, Base$Nperm, Base$Ntimesteps)
  ordinary <- if("n_new_background_detections" %in% names(S))
    .INApestPoAPointMatrix(S, "n_new_background_detections", Base$Nperm, Base$Ntimesteps) else
    Base$BackgroundDetections

  Base$BackgroundDetections <- control + ordinary
  Base$SurveillanceDetections <- Base$BackgroundDetections
  Base$BackgroundEventContract <- "explicit"
  Base$BackgroundEventField <- "n_control_detections + n_new_background_detections"
  Base$AdapterName <- ".INApestPoAAdapterVertebratePoint"
  Base
}


###############################################################################
### Thin public wrappers
###############################################################################

INApestVertebrateNodePoA <- function(ModelArgs=list(), ModelResults=NULL,
                                      ObservationModel=NULL, ...)
{
  .INApestPoAVertebrateRequireCore()
  if(is.null(ObservationModel))
    ObservationModel <- INApestPoAVertebrateObservationModel()
  .INApestPoARun(
    "INApestVertebrateNode", .INApestPoAAdapterVertebrateNode,
    ModelArgs, ModelResults, ..., ObservationModel=ObservationModel,
    ResultClass="INApestVertebrateNodePoA"
  )
}


INApestVertebratePointPoA <- function(ModelArgs=list(), ModelResults=NULL,
                                       ObservationModel=NULL, ...)
{
  .INApestPoAVertebrateRequireCore()
  if(is.null(ObservationModel))
    ObservationModel <- INApestPoAVertebrateObservationModel()
  .INApestPoARun(
    "INApestVertebratePoint", .INApestPoAAdapterVertebratePoint,
    ModelArgs, ModelResults, ..., ObservationModel=ObservationModel,
    ResultClass="INApestVertebratePointPoA"
  )
}
