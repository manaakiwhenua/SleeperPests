###############################################################################
### INApest pathogen proof-of-freedom (PoF) companion
###
### PoF is defined from latent biological pathogen state, not from detections.
### Surveillance results update the probability of that latent state by
### likelihood weighting across stochastic trajectories (particles).
###############################################################################

.ipof_clip01 <- function(x) pmin(1, pmax(0, x))

.ipof_norm_weights <- function(w) {
  if (!is.numeric(w) || any(!is.finite(w)) || any(w < 0))
    stop("Particle weights must be finite and non-negative.")
  z <- sum(w)
  if (!is.finite(z) || z <= 0)
    stop("Observation history has zero likelihood under all supplied particles.")
  w / z
}

.ipof_pathogen_spec <- function(Pathogen) {
  if (inherits(Pathogen, "INApestPointPathogenInteraction")) return(Pathogen$Pathogen)
  if (inherits(Pathogen, "INApestPathogen")) return(Pathogen)
  stop("Pathogen must be an INApestPathogen or INApestPointPathogenInteraction object.")
}


### Resolve either an in-memory model result or the traditional saved-output
### contract used by node Meta-family functions. Character input may be either
### a complete .rds result object or the filename stem preceding standard
### INApest output suffixes. A list with OutputDir + ModelName is also accepted.
.ipof_model_output <- function(ModelOutput) {
  if (is.list(ModelOutput) && !is.null(ModelOutput$OutputDir) &&
      !is.null(ModelOutput$ModelName) &&
      is.null(ModelOutput$PathogenStateResults) &&
      is.null(ModelOutput$PathogenStageResults) &&
      is.null(ModelOutput$PathogenPresentResults) &&
      is.null(ModelOutput$PointHistory)) {
    od <- ModelOutput$OutputDir
    if (length(od) != 1L || is.na(od)) od <- ""
    ModelOutput <- file.path(od, ModelOutput$ModelName)
  }

  if (is.list(ModelOutput)) return(ModelOutput)

  if (!is.character(ModelOutput) || length(ModelOutput) != 1L || is.na(ModelOutput))
    stop("ModelOutput must be a result object, an .rds result file, a standard output filename stem, or list(OutputDir=..., ModelName=...).")

  if (file.exists(ModelOutput) && grepl("\\.rds$", ModelOutput, ignore.case = TRUE)) {
    z <- readRDS(ModelOutput)
    if (!is.list(z))
      stop("A directly supplied .rds ModelOutput must contain a saved result object. For individual pathogen arrays, supply their common filename stem instead.")
    return(z)
  }

  stem <- ModelOutput
  point_candidates <- c(
    paste0(stem, "_PointResults.rds"),
    paste0(stem, "_MetaPointParallelResults.rds"),
    paste0(stem, "_VertebratePointResults.rds"),
    paste0(stem, "_PointTransitionMatrixParallelResults.rds")
  )
  hit <- point_candidates[file.exists(point_candidates)]
  if (length(hit)) return(readRDS(hit[1L]))

  out <- list()
  candidates <- c(
    PathogenPresentResults = "PathogenPresentLargeOut.rds",
    PathogenDetectedResults = "PathogenDetectedLargeOut.rds",
    PathogenStateResults = "PathogenStateLargeOut.rds",
    PathogenStageResults = "PathogenStageLargeOut.rds",
    PathogenDeathResults = "PathogenDeathLargeOut.rds",
    NewInfectionResults = "NewInfectionLargeOut.rds",
    OriginalPathogenPresentResults = "OriginalPathogenPresentLargeOut.rds",
    ReintroducedPathogenPresentResults = "ReintroducedPathogenPresentLargeOut.rds",
    OriginalPathogenStateResults = "OriginalPathogenStateLargeOut.rds",
    ReintroducedPathogenStateResults = "ReintroducedPathogenStateLargeOut.rds",
    OriginalPathogenStageResults = "OriginalPathogenStageLargeOut.rds",
    ReintroducedPathogenStageResults = "ReintroducedPathogenStageLargeOut.rds"
  )
  for (nm in names(candidates)) {
    f <- paste0(stem, candidates[[nm]])
    if (file.exists(f)) out[[nm]] <- readRDS(f)
  }
  if (!length(out))
    stop("No standard pathogen outputs found for ModelOutput stem: ", stem)
  out$ModelNameStem <- stem
  out
}

.ipof_active_states <- function(Pathogen) {
  Pathogen <- .ipof_pathogen_spec(Pathogen)
  switch(Pathogen$Model,
         Binary = "Present",
         SIS = "I",
         SIR = "I",
         SEIR = c("E", "I"),
         stop("Unsupported pathogen model: ", Pathogen$Model))
}

### Resolve pathogen detection probability using the same parameter schedule as
### the simulator. Returns one probability per spatial surveillance unit.
.ipof_detection_prob <- function(Pathogen, timestep, n_nodes, Ntimesteps,
                                 n_landuses = NULL) {
  if (Pathogen$Model == "Binary") {
    # binary_detect performs a random draw, so resolve its parameter directly.
    x <- Pathogen$DetectionProb
    if (is.function(x)) {
      fm <- names(formals(x)); a <- list(timestep=timestep,n_nodes=n_nodes,Ntimesteps=Ntimesteps)
      if (!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a),fm)]
      x <- do.call(x,a)
    }
    d <- dim(x)
    if (!is.null(d)) {
      if (length(d)==2L && all(d==c(n_nodes,Ntimesteps))) return(.ipof_clip01(as.numeric(x[,timestep])))
      stop("Binary DetectionProb must be scalar, length nodes, nodes x Ntimesteps, or resolver function.")
    }
    x <- as.numeric(x)
    if (length(x)==1L) return(rep(.ipof_clip01(x),n_nodes))
    if (length(x)==n_nodes) return(.ipof_clip01(x))
    if (length(x)==Ntimesteps && Ntimesteps != n_nodes) return(rep(.ipof_clip01(x[timestep]),n_nodes))
    stop("Binary DetectionProb has unsupported or ambiguous shape.")
  }
  context <- list(n_nodes=n_nodes,Ntimesteps=Ntimesteps)
  if (!is.null(n_landuses)) context$n_landuses <- n_landuses
  .ipof_clip01(Pathogen$Engine$Resolve(Pathogen$DetectionProb,timestep,context,"DetectionProb"))
}

### Convert supported INApest outputs into a common latent representation.
### Returned active_count is [surveillance unit x timestep x particle].
### A surveillance unit is a node for node models and a point for point models.
INApestPathogenFreedomState <- function(ModelOutput, Pathogen) {
  ModelOutput <- .ipof_model_output(ModelOutput)
  PathogenInput <- Pathogen
  Pathogen <- .ipof_pathogen_spec(Pathogen)
  active_states <- .ipof_active_states(Pathogen)

  # Binary INApest / INApestParallel.
  if (is.list(ModelOutput) && !is.null(ModelOutput$PathogenPresentResults)) {
    x <- ModelOutput$PathogenPresentResults
    if (length(dim(x)) != 3L) stop("PathogenPresentResults must be nodes x timesteps x particles.")
    active <- array(as.numeric(x > 0), dim=dim(x))
    freedom <- apply(active, c(2,3), sum) == 0
    return(list(Type="binary", ActiveCount=active, Freedom=freedom,
                Nnodes=dim(x)[1], Ntimesteps=dim(x)[2], Nparticles=dim(x)[3]))
  }

  # Aggregate Meta: nodes x states x time x particles.
  if (is.list(ModelOutput) && !is.null(ModelOutput$PathogenStateResults)) {
    x <- ModelOutput$PathogenStateResults
    d <- dim(x); st <- dimnames(x)[[2]]
    if (length(d)==4L && !is.null(st)) {
      active_idx <- match(active_states,st); if(anyNA(active_idx)) stop("Required active pathogen state missing from PathogenStateResults.")
      active <- apply(x[,active_idx,,,drop=FALSE],c(1,3,4),sum)
      if (length(dim(active))==2L) active <- array(active,dim=c(d[1],d[3],d[4]))
      freedom <- apply(active,c(2,3),sum)==0
      return(list(Type="meta",ActiveCount=active,Freedom=freedom,Nnodes=d[1],Ntimesteps=d[3],Nparticles=d[4]))
    }
    # MLU: nodes x landuse x states x time x particles; aggregate land uses to node.
    st <- dimnames(x)[[3]]
    if (length(d)==5L && !is.null(st)) {
      active_idx <- match(active_states,st); if(anyNA(active_idx)) stop("Required active pathogen state missing from MLU PathogenStateResults.")
      active <- apply(x[,,active_idx,,,drop=FALSE],c(1,4,5),sum)
      if(length(dim(active))==2L) active <- array(active,dim=c(d[1],d[4],d[5]))
      freedom <- apply(active,c(2,3),sum)==0
      return(list(Type="mlu",ActiveCount=active,Freedom=freedom,Nnodes=d[1],Nlanduses=d[2],Ntimesteps=d[4],Nparticles=d[5]))
    }
  }

  # Transition matrix: nodes x stages x states x time x particles.
  if (is.list(ModelOutput) && !is.null(ModelOutput$PathogenStageResults)) {
    x <- ModelOutput$PathogenStageResults; d <- dim(x); st <- dimnames(x)[[3]]
    if (length(d)!=5L || is.null(st)) stop("PathogenStageResults must be nodes x stages x states x timesteps x particles with state dimnames.")
    active_idx <- match(active_states,st); if(anyNA(active_idx)) stop("Required active pathogen state missing from PathogenStageResults.")
    active <- apply(x[,,active_idx,,,drop=FALSE],c(1,4,5),sum)
    if(length(dim(active))==2L) active <- array(active,dim=c(d[1],d[4],d[5]))
    freedom <- apply(active,c(2,3),sum)==0
    return(list(Type="transition",ActiveCount=active,Freedom=freedom,Nnodes=d[1],Nstages=d[2],Ntimesteps=d[4],Nparticles=d[5]))
  }

  # Point models: persistent state reconstructed from PointHistory snapshots.
  if (is.list(ModelOutput) && is.data.frame(ModelOutput$PointHistory)) {
    h <- ModelOutput$PointHistory
    state_field <- if (inherits(PathogenInput,"INApestPointPathogenInteraction")) PathogenInput$PathogenStateField else "pathogen_state"
    if (!(state_field %in% names(h))) stop("PointHistory does not contain pathogen state field '",state_field,"'.")
    if (!all(c("perm","timestep") %in% names(h))) stop("PointHistory requires perm and timestep columns.")
    if (is.data.frame(ModelOutput$Summary) && nrow(ModelOutput$Summary) && all(c("perm","timestep") %in% names(ModelOutput$Summary))) {
      np <- max(ModelOutput$Summary$perm); nt <- max(ModelOutput$Summary$timestep)
    } else {
      np <- if(nrow(h)) max(h$perm) else 0L; nt <- if(nrow(h)) max(h$timestep) else 0L
    }
    total_active <- matrix(0,nt,np)
    for (pp in seq_len(np)) for (tt in seq_len(nt)) {
      z <- h[h$perm==pp & h$timestep==tt,,drop=FALSE]
      total_active[tt,pp] <- sum(z[[state_field]] %in% active_states)
    }
    freedom <- total_active==0
    return(list(Type="point",TotalActive=total_active,Freedom=freedom,Ntimesteps=nt,Nparticles=np,StateField=state_field))
  }
  stop("Unsupported ModelOutput. Supply a pathogen-capable INApest result object.")
}

### Build per-particle likelihood for one observation round.
### Observation can be:
###   scalar 0/1: system-wide no-detection/detection;
###   node vector 0/1/NA: node-specific results; NA = node not observed.
### For abundance models DetectionProb is interpreted per infectious host. Exposed
### SEIR hosts count against freedom but are not detectable under the default
### pathogen surveillance model.
INApestPathogenObservationLikelihood <- function(ModelOutput, Pathogen, timestep,
                                                  Observation = 0) {
  ModelOutput <- .ipof_model_output(ModelOutput)
  PathogenInput <- Pathogen
  Pathogen <- .ipof_pathogen_spec(Pathogen)
  fs <- INApestPathogenFreedomState(ModelOutput,PathogenInput)
  if (timestep < 1L || timestep > fs$Ntimesteps) stop("timestep outside model output.")

  if (fs$Type == "point") {
    # Point likelihood is evaluated directly from each point so spatial/function
    # detection schedules can use point attributes in future extensions. Current
    # common Pathogen DetectionProb is treated as per-point probability.
    h <- ModelOutput$PointHistory
    state_field <- fs$StateField
    active_detectable <- "I"
    if (Pathogen$Model=="Binary") active_detectable <- "Present"
    np <- fs$Nparticles; L0 <- numeric(np)
    for(pp in seq_len(np)) {
      z <- h[h$perm==pp & h$timestep==timestep,,drop=FALSE]
      inf <- z[[state_field]] %in% active_detectable
      if(!nrow(z) || !any(inf)) { L0[pp] <- 1; next }
      if (inherits(PathogenInput, "INApestPointPathogenInteraction") &&
          is.function(PathogenInput$ResolvePoint)) {
        pd <- PathogenInput$ResolvePoint(
          Pathogen$DetectionProb, z, timestep, pp, "DetectionProb"
        )
      } else {
        # Backwards-compatible fallback for interaction objects created before
        # ResolvePoint was exposed. Semantics match the point helper resolver.
        pd <- Pathogen$DetectionProb
        if(is.function(pd)) {
          fm <- names(formals(pd)); a <- list(points=z,timestep=timestep,perm=pp)
          if(!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a),fm)]
          pd <- do.call(pd,a)
        } else if(length(pd)==1L) pd <- rep(as.numeric(pd),nrow(z))
        else { pd <- as.numeric(pd); if(length(pd)!=nrow(z)) stop("Point DetectionProb must be scalar, function(points,timestep,perm), or one value per point snapshot.") }
        if(length(pd)==1L) pd <- rep(pd,nrow(z))
        if(length(pd)!=nrow(z) || any(!is.finite(pd))) stop("Point DetectionProb must resolve to one finite value per point.")
      }
      if(any(pd<0|pd>1)) stop("Point DetectionProb must resolve to [0,1] per point.")
      L0[pp] <- prod(1-.ipof_clip01(pd[inf]))
    }
    if(length(Observation)!=1L || is.na(Observation) || !(Observation %in% c(0,1))) stop("Point v1 Observation must be scalar 0/1.")
    return(if(Observation==0) L0 else 1-L0)
  }

  # Construct node-level probability of no pathogen detection for each particle.
  # Surveillance is applied to infectious I only; E counts against freedom but
  # is not detected by the default pathogen observation model.
  if (fs$Type == "mlu") {
    x <- ModelOutput$PathogenStateResults
    I <- x[, , "I", timestep, , drop=FALSE]
    p_unit <- .ipof_detection_prob(Pathogen,timestep,fs$Nnodes,fs$Ntimesteps,n_landuses=fs$Nlanduses)
    pmat <- matrix(p_unit,nrow=fs$Nnodes,ncol=fs$Nlanduses)
    no_node <- matrix(1,nrow=fs$Nnodes,ncol=fs$Nparticles)
    for(pp in seq_len(fs$Nparticles))
      no_node[,pp] <- apply((1-pmat)^matrix(I[,,,,pp],nrow=fs$Nnodes,ncol=fs$Nlanduses),1,prod)
  } else if (fs$Type == "transition") {
    x <- ModelOutput$PathogenStageResults
    I_node <- apply(x[, , "I", timestep, ,drop=FALSE],c(1,5),sum)
    I_node <- matrix(I_node,nrow=fs$Nnodes,ncol=fs$Nparticles)
    p <- .ipof_detection_prob(Pathogen,timestep,fs$Nnodes,fs$Ntimesteps)
    no_node <- (1-p)^I_node
  } else {
    active <- fs$ActiveCount[,timestep,,drop=FALSE]
    active <- matrix(active,nrow=fs$Nnodes,ncol=fs$Nparticles)
    if(Pathogen$Model=="SEIR") {
      x <- ModelOutput$PathogenStateResults
      if(!is.null(x)) active <- matrix(x[,"I",timestep,],nrow=fs$Nnodes,ncol=fs$Nparticles)
    }
    p <- .ipof_detection_prob(Pathogen,timestep,fs$Nnodes,fs$Ntimesteps)
    no_node <- (1-p)^active
  }

  if(length(Observation)==1L) {
    if(is.na(Observation) || !(Observation %in% c(0,1))) stop("Scalar Observation must be 0 (no detection) or 1 (one or more detections).")
    l0 <- apply(no_node,2,prod)
    return(if(Observation==0) l0 else 1-l0)
  }
  if(length(Observation)!=fs$Nnodes) stop("Node Observation must have one value per node.")
  if(any(!is.na(Observation) & !(Observation %in% c(0,1)))) stop("Node Observation values must be 0, 1 or NA.")
  out <- rep(1,fs$Nparticles)
  for(i in seq_len(fs$Nnodes)) if(!is.na(Observation[i])) {
    out <- out * if(Observation[i]==0) no_node[i,] else (1-no_node[i,])
  }
  out
}

.ipof_weighted <- function(x,w) { if(!length(x)||!length(w)||sum(w)<=0)return(NA_real_); sum(as.numeric(x)*w)/sum(w) }

INApestPathogenPoFCore <- function(Freedom, ObservationLikelihood,
                                   PriorWeights = NULL) {
  Freedom <- as.logical(Freedom)
  n <- length(Freedom)
  if(length(ObservationLikelihood)!=n) stop("ObservationLikelihood and Freedom lengths differ.")
  if(is.null(PriorWeights)) PriorWeights <- rep(1/n,n)
  w0 <- .ipof_norm_weights(PriorWeights)
  prior <- sum(w0 * Freedom)
  w1 <- .ipof_norm_weights(w0 * ObservationLikelihood)
  posterior <- sum(w1 * Freedom)
  list(PriorPoF=prior,PosteriorPoF=posterior,PosteriorWeights=w1,
       EffectiveParticles=1/sum(w1^2),ObservationEvidence=sum(w0*ObservationLikelihood))
}

### Main post-processing wrapper. ObservationHistory may be a named/list sequence
### indexed by timestep. Missing timesteps contribute no observation likelihood.
### Weights accumulate through time. This is valid when observation results do
### not alter subsequent simulated dynamics. For pathogen-triggered information
### and response, use ReplayRequired=TRUE as a diagnostic until dynamic replay is
### implemented in the wrapper.
INApestPathogenPoF <- function(ModelOutput, Pathogen,
                               ObservationHistory = NULL,
                               PriorWeights = NULL,
                               ReplayRequired = NULL) {
  ModelOutput <- .ipof_model_output(ModelOutput)
  PathogenInput <- Pathogen
  PathogenSpec <- .ipof_pathogen_spec(Pathogen)
  fs <- INApestPathogenFreedomState(ModelOutput,PathogenInput)
  np <- fs$Nparticles; nt <- fs$Ntimesteps
  w <- if(is.null(PriorWeights)) rep(1/np,np) else .ipof_norm_weights(PriorWeights)
  if(is.null(ReplayRequired)) ReplayRequired <- isTRUE(PathogenSpec$DetectionTriggersInfo)
  ans <- data.frame(timestep=seq_len(nt),PriorPoF=NA_real_,PosteriorPoF=NA_real_,
                    ObservationEvidence=NA_real_,EffectiveParticles=NA_real_)
  for(tt in seq_len(nt)) {
    f <- fs$Freedom[tt,]
    ans$PriorPoF[tt] <- sum(w*f)
    obs <- NULL
    if(!is.null(ObservationHistory)) {
      if(is.list(ObservationHistory) && length(ObservationHistory)>=tt) obs <- ObservationHistory[[tt]]
      else if(is.numeric(ObservationHistory) && length(ObservationHistory)>=tt) obs <- ObservationHistory[tt]
    }
    if(is.null(obs) || (length(obs)==1L && is.na(obs))) L <- rep(1,np)
    else L <- INApestPathogenObservationLikelihood(ModelOutput,PathogenInput,tt,obs)
    ev <- sum(w*L)
    w <- .ipof_norm_weights(w*L)
    ans$PosteriorPoF[tt] <- sum(w*f); ans$ObservationEvidence[tt] <- ev; ans$EffectiveParticles[tt] <- 1/sum(w^2)
  }
  structure(list(Summary=ans,PosteriorWeights=w,Freedom=fs$Freedom,
                 ReplayRequired=ReplayRequired,
                 Interpretation=if(ReplayRequired)
                   "Observation likelihood is valid for the current state, but future management can depend on pathogen detection/information. Dynamic history replay is required for fully conditioned future trajectories."
                 else "Post-hoc sequential likelihood weighting is valid when observations do not alter future simulated dynamics."),
            class=c("INApestPathogenPoF","list"))
}

### Exact one-round binary Bayesian benchmark used by the regression suite.
INApestPathogenPoFExactBinary <- function(StateProbability, PathogenPresent,
                                          DetectionProb, Observation = 0) {
  StateProbability <- .ipof_norm_weights(StateProbability)
  P <- as.matrix(PathogenPresent)
  if(nrow(P)!=length(StateProbability)) stop("One row of PathogenPresent required per exact state.")
  d <- rep_len(DetectionProb,ncol(P)); if(any(d<0|d>1))stop("DetectionProb must be in [0,1].")
  no_det <- apply(sweep(P, 2, 1-d, function(present, q) q^present), 1, prod)
  L <- if(Observation==0) no_det else 1-no_det
  free <- rowSums(P)==0
  INApestPathogenPoFCore(free,L,StateProbability)
}


###############################################################################
### Binary pathogen origin and sequential compatible-history proof of freedom
###############################################################################

### Four mutually exclusive system-level origin classes:
### 0 free; 1 original only; 2 reintroduced only; 3 mixed original+reintroduced.
INApestPathogenOriginState <- function(ModelOutput, Pathogen=NULL) {
  ModelOutput <- .ipof_model_output(ModelOutput)
  if(!is.null(ModelOutput$PathogenPresentResults)) {
    P<-ModelOutput$PathogenPresentResults; O<-ModelOutput$OriginalPathogenPresentResults; R<-ModelOutput$ReintroducedPathogenPresentResults
    if(is.null(O)||is.null(R)) stop("Origin-aware binary PoF requires original and reintroduced pathogen outputs. Run with TrackOrigin=TRUE.")
    if(!identical(dim(P),dim(O))||!identical(dim(P),dim(R))||length(dim(P))!=3L) stop("Binary pathogen origin arrays must have identical nodes x timesteps x particles dimensions.")
    P<-P>0;O<-O>0;R<-R>0;if(any((O|R)!=P))stop("Origin arrays violate union(origin) == PathogenPresent.")
    op<-apply(O,c(2,3),any);rp<-apply(R,c(2,3),any);cls<-matrix(as.integer(op)+2L*as.integer(rp),nrow=dim(P)[2],ncol=dim(P)[3])
    return(list(Type="binary",Class=cls,OriginalPresent=op,ReintroducedPresent=rp,Freedom=!(op|rp),Ntimesteps=dim(P)[2],Nparticles=dim(P)[3]))
  }
  X<-ModelOutput$PathogenStateResults; O<-ModelOutput$OriginalPathogenStateResults; R<-ModelOutput$ReintroducedPathogenStateResults
  if(!is.null(X) && length(dim(X))==4L) {
    if(is.null(Pathogen)) stop("Pathogen is required to interpret compartment origin arrays.")
    Pathogen<-.ipof_pathogen_spec(Pathogen); if(Pathogen$Model=="Binary")stop("Binary origin state expected binary outputs.")
    if(is.null(O)||is.null(R))stop("Origin-aware compartment PoF requires OriginalPathogenStateResults and ReintroducedPathogenStateResults. Run Meta with TrackOrigin=TRUE.")
    if(!identical(dim(X),dim(O))||!identical(dim(X),dim(R))||length(dim(X))!=4L)stop("Compartment origin arrays must have identical nodes x states x timesteps x particles dimensions.")
    st<-dimnames(X)[[2]];if(is.null(st))stop("PathogenStateResults requires pathogen-state dimnames.")
    active<-.ipof_active_states(Pathogen);ai<-match(active,st);if(anyNA(ai))stop("Required active pathogen state missing from origin arrays.")
    if(any(O[,ai,,,drop=FALSE]+R[,ai,,,drop=FALSE]!=X[,ai,,,drop=FALSE]))stop("Compartment origin invariant failed: original + reintroduced counts must equal active pathogen-state counts.")
    inactive<-setdiff(seq_along(st),ai);if(length(inactive)&&(any(O[,inactive,,,drop=FALSE]!=0)||any(R[,inactive,,,drop=FALSE]!=0)))stop("Origin counts must be zero outside active E/I states.")
    op<-apply(O[,ai,,,drop=FALSE],c(3,4),sum)>0;rp<-apply(R[,ai,,,drop=FALSE],c(3,4),sum)>0
    if(is.null(dim(op)))op<-matrix(op,nrow=dim(X)[3],ncol=dim(X)[4]);if(is.null(dim(rp)))rp<-matrix(rp,nrow=dim(X)[3],ncol=dim(X)[4])
    cls<-matrix(as.integer(op)+2L*as.integer(rp),nrow=dim(X)[3],ncol=dim(X)[4])
    return(list(Type="compartment",Class=cls,OriginalPresent=op,ReintroducedPresent=rp,Freedom=!(op|rp),Ntimesteps=dim(X)[3],Nparticles=dim(X)[4]))
  }
  # Multiple-land-use compartment arrays: nodes x land use x states x time x particles.
  if(!is.null(ModelOutput$PathogenStateResults) && length(dim(ModelOutput$PathogenStateResults))==5L) {
    X<-ModelOutput$PathogenStateResults; O<-ModelOutput$OriginalPathogenStateResults; R<-ModelOutput$ReintroducedPathogenStateResults
    if(is.null(Pathogen)) stop("Pathogen is required to interpret MLU origin arrays.")
    Pathogen<-.ipof_pathogen_spec(Pathogen); if(is.null(O)||is.null(R)) stop("Origin-aware MLU PoF requires OriginalPathogenStateResults and ReintroducedPathogenStateResults. Run MLU with TrackOrigin=TRUE.")
    if(!identical(dim(X),dim(O))||!identical(dim(X),dim(R))) stop("MLU origin arrays must have identical dimensions.")
    st<-dimnames(X)[[3]]; active<-.ipof_active_states(Pathogen); ai<-match(active,st); if(anyNA(ai)) stop("Required active pathogen state missing from MLU origin arrays.")
    if(any(O[,,ai,,,drop=FALSE]+R[,,ai,,,drop=FALSE]!=X[,,ai,,,drop=FALSE])) stop("MLU origin invariant failed.")
    inactive<-setdiff(seq_along(st),ai); if(length(inactive)&&(any(O[,,inactive,,,drop=FALSE]!=0)||any(R[,,inactive,,,drop=FALSE]!=0))) stop("MLU origin counts must be zero outside active E/I states.")
    op<-apply(O[,,ai,,,drop=FALSE],c(4,5),sum)>0; rp<-apply(R[,,ai,,,drop=FALSE],c(4,5),sum)>0
    if(is.null(dim(op))) op<-matrix(op,nrow=dim(X)[4],ncol=dim(X)[5]); if(is.null(dim(rp))) rp<-matrix(rp,nrow=dim(X)[4],ncol=dim(X)[5])
    cls<-matrix(as.integer(op)+2L*as.integer(rp),nrow=dim(X)[4],ncol=dim(X)[5])
    return(list(Type="mlu",Class=cls,OriginalPresent=op,ReintroducedPresent=rp,Freedom=!(op|rp),Ntimesteps=dim(X)[4],Nparticles=dim(X)[5]))
  }
  # Demographic transition arrays: nodes x stages x states x time x particles.
  Xs<-ModelOutput$PathogenStageResults; Os<-ModelOutput$OriginalPathogenStageResults; Rs<-ModelOutput$ReintroducedPathogenStageResults
  if(!is.null(Xs)) {
    if(is.null(Pathogen)) stop("Pathogen is required to interpret transition-matrix origin arrays.")
    Pathogen<-.ipof_pathogen_spec(Pathogen); if(is.null(Os)||is.null(Rs)) stop("Origin-aware transition PoF requires OriginalPathogenStageResults and ReintroducedPathogenStageResults. Run transition matrix with TrackOrigin=TRUE.")
    if(!identical(dim(Xs),dim(Os))||!identical(dim(Xs),dim(Rs))||length(dim(Xs))!=5L) stop("Transition origin arrays must have identical nodes x stages x states x timesteps x particles dimensions.")
    st<-dimnames(Xs)[[3]]; active<-.ipof_active_states(Pathogen); ai<-match(active,st); if(anyNA(ai)) stop("Required active pathogen state missing from transition origin arrays.")
    if(any(Os[,,ai,,,drop=FALSE]+Rs[,,ai,,,drop=FALSE]!=Xs[,,ai,,,drop=FALSE])) stop("Transition origin invariant failed.")
    inactive<-setdiff(seq_along(st),ai); if(length(inactive)&&(any(Os[,,inactive,,,drop=FALSE]!=0)||any(Rs[,,inactive,,,drop=FALSE]!=0))) stop("Transition origin counts must be zero outside active E/I states.")
    op<-apply(Os[,,ai,,,drop=FALSE],c(4,5),sum)>0; rp<-apply(Rs[,,ai,,,drop=FALSE],c(4,5),sum)>0
    if(is.null(dim(op))) op<-matrix(op,nrow=dim(Xs)[4],ncol=dim(Xs)[5]); if(is.null(dim(rp))) rp<-matrix(rp,nrow=dim(Xs)[4],ncol=dim(Xs)[5])
    cls<-matrix(as.integer(op)+2L*as.integer(rp),nrow=dim(Xs)[4],ncol=dim(Xs)[5])
    return(list(Type="transition",Class=cls,OriginalPresent=op,ReintroducedPresent=rp,Freedom=!(op|rp),Ntimesteps=dim(Xs)[4],Nparticles=dim(Xs)[5]))
  }
  # Point and point-transition models: persistent lineage/origin fields live
  # directly on each point snapshot. This retains richer lineage roots while the
  # sequential PoF calculation collapses them to the common four origin classes.
  if(is.list(ModelOutput) && is.data.frame(ModelOutput$PointHistory)) {
    h<-ModelOutput$PointHistory
    if(is.null(Pathogen)) stop("Pathogen interaction is required to interpret point origin fields.")
    pi<-Pathogen; ps<-.ipof_pathogen_spec(Pathogen)
    origin_field<-if(inherits(pi,"INApestPointPathogenInteraction")&&!is.null(pi$PathogenOriginField))pi$PathogenOriginField else "pathogen_origin"
    lineage_field<-if(inherits(pi,"INApestPointPathogenInteraction")&&!is.null(pi$PathogenLineageField))pi$PathogenLineageField else "pathogen_lineage"
    state_field<-if(inherits(pi,"INApestPointPathogenInteraction"))pi$PathogenStateField else "pathogen_state"
    if(!all(c("perm","timestep",state_field,origin_field,lineage_field)%in%names(h))) stop("Origin-aware point PoF requires pathogen state, origin and lineage fields in PointHistory. Run with TrackOrigin=TRUE.")
    if(is.data.frame(ModelOutput$Summary)&&nrow(ModelOutput$Summary)&&all(c("perm","timestep")%in%names(ModelOutput$Summary))){np<-max(ModelOutput$Summary$perm);nt<-max(ModelOutput$Summary$timestep)}else{np<-if(nrow(h))max(h$perm)else 0L;nt<-if(nrow(h))max(h$timestep)else 0L}
    active<-.ipof_active_states(ps); op<-matrix(FALSE,nt,np);rp<-matrix(FALSE,nt,np);oc<-matrix(0L,nt,np);rc<-matrix(0L,nt,np)
    for(pp in seq_len(np))for(tt in seq_len(nt)){
      z<-h[h$perm==pp & h$timestep==tt,,drop=FALSE];if(!nrow(z))next
      aa<-z[[state_field]]%in%active
      if(any(aa & (is.na(z[[origin_field]])|is.na(z[[lineage_field]])|!nzchar(as.character(z[[lineage_field]]))))) stop("Active point pathogen states must carry origin and lineage identifiers.")
      bad<-aa & !(z[[origin_field]]%in%c("original","reintroduced"));if(any(bad,na.rm=TRUE))stop("Point pathogen origin must be original or reintroduced for active states.")
      zo<-aa & z[[origin_field]]=="original";zr<-aa & z[[origin_field]]=="reintroduced"
      op[tt,pp]<-any(zo,na.rm=TRUE);rp[tt,pp]<-any(zr,na.rm=TRUE)
      oc[tt,pp]<-length(unique(z[[lineage_field]][zo]));rc[tt,pp]<-length(unique(z[[lineage_field]][zr]))
    }
    cls<-matrix(as.integer(op)+2L*as.integer(rp),nrow=nt,ncol=np)
    return(list(Type="point",Class=cls,OriginalPresent=op,ReintroducedPresent=rp,Freedom=!(op|rp),OriginalLineageCount=oc,ReintroducedLineageCount=rc,Ntimesteps=nt,Nparticles=np,StateField=state_field,OriginField=origin_field,LineageField=lineage_field))
  }
  stop("No supported origin-aware pathogen outputs found.")
}

.ipof_match_binary_detection <- function(SimulatedCounts, Observed, Match=c("binary","exact")) {
  Match <- match.arg(Match)
  if (length(Observed)!=1L || is.na(Observed) || !is.finite(Observed) || Observed < 0)
    stop("Observed pathogen detections must be one non-negative finite count per observed timestep.")
  if (Match=="binary") return((SimulatedCounts > 0) == (Observed > 0))
  SimulatedCounts == Observed
}

### Preserve the likelihood posterior mass of each biological origin class while
### restricting future propagation to realised trajectories with observation
### histories compatible with the real surveillance record.
.ipof_rescale_origin_matched <- function(PriorWeights, Match, OriginClass,
                                         PosteriorWeights, timestep) {
  raw <- PriorWeights * as.numeric(Match)
  out <- rep(0,length(PriorWeights))
  labs <- c("free","original-only","reintroduced-only","mixed")
  for (cc in 0:3) {
    target <- sum(PosteriorWeights[OriginClass==cc])
    idx <- Match & OriginClass==cc
    if (target > 0) {
      den <- sum(raw[idx])
      if (den <= 0)
        stop("Likelihood update retained positive posterior mass for the ",labs[cc+1L],
             " class at timestep ",timestep,
             " but no realised matching particles are available to propagate that class. Increase Nperm or shorten the conditioning horizon.")
      out[idx] <- raw[idx]/den*target
    }
  }
  .ipof_norm_weights(out)
}

### Sequential binary PoF with pathogen-triggered information/management.
### Current-round posterior uses the lower-variance likelihood calculation.
### Between observed rounds only realised detection histories matching the
### observation are propagated, with four origin-class masses rescaled to the
### likelihood posterior. This is the pathogen analogue of the validated PoA
### compatible-history filter, extended so persistence/reintroduction cannot be
### distorted within the broad "present" class.
INApestPathogenPoFSequentialBinary <- function(ModelOutput, Pathogen,
                                                ObservationHistory,
                                                PriorWeights=NULL,
                                                ObservationMatch=c("binary","exact"),
                                                ReturnParticles=FALSE) {
  ModelOutput <- .ipof_model_output(ModelOutput)
  PathogenSpec <- .ipof_pathogen_spec(Pathogen)
  if (!identical(PathogenSpec$Model,"Binary")) stop("Sequential binary PoF requires Pathogen Model='Binary'.")
  origin <- INApestPathogenOriginState(ModelOutput)
  D <- ModelOutput$PathogenDetectedResults
  if (is.null(D) || length(dim(D))!=3L) stop("PathogenDetectedResults nodes x timesteps x particles are required for compatible-history propagation.")
  if (dim(D)[2]!=origin$Ntimesteps || dim(D)[3]!=origin$Nparticles) stop("PathogenDetectedResults dimensions do not match origin arrays.")
  ObservationMatch <- match.arg(ObservationMatch)
  if (!is.data.frame(ObservationHistory) || !all(c("Timestep","PathogenDetections") %in% names(ObservationHistory)))
    stop("ObservationHistory must be a data.frame with Timestep and PathogenDetections columns.")
  obs <- ObservationHistory[order(ObservationHistory$Timestep),,drop=FALSE]
  obs$Timestep <- as.integer(obs$Timestep)
  if (!nrow(obs) || anyDuplicated(obs$Timestep) || any(obs$Timestep<1L | obs$Timestep>origin$Ntimesteps))
    stop("ObservationHistory must contain unique valid observed timesteps.")
  np <- origin$Nparticles
  w <- if(is.null(PriorWeights)) rep(1/np,np) else .ipof_norm_weights(PriorWeights)
  rows <- vector("list",nrow(obs)); wh <- if(ReturnParticles) matrix(NA_real_,np,nrow(obs)) else NULL
  for (rr in seq_len(nrow(obs))) {
    tt <- obs$Timestep[rr]; prior <- w; cls <- origin$Class[tt,]
    L <- INApestPathogenObservationLikelihood(ModelOutput,Pathogen,tt,
                                               as.integer(obs$PathogenDetections[rr] > 0))
    ev <- sum(prior*L)
    if (!is.finite(ev) || ev <= 0) stop("Observed pathogen evidence has zero likelihood at timestep ",tt)
    post <- .ipof_norm_weights(prior*L)
    sim_count <- apply(D[,tt,,drop=FALSE],3,sum)
    match <- .ipof_match_binary_detection(sim_count,obs$PathogenDetections[rr],ObservationMatch)
    mass <- vapply(0:3,function(cc)sum(post[cls==cc]),numeric(1))
    rows[[rr]] <- data.frame(
      Round=rr,Timestep=tt,
      ObservationProbability=ev,
      PriorPoF=sum(prior[cls==0L]), PosteriorPoF=mass[1],
      PosteriorOriginalPersistence=mass[2]+mass[4],
      PosteriorOriginalEliminated=mass[1]+mass[3],
      PosteriorReintroductionPresent=mass[3]+mass[4],
      PosteriorReintroductionOnly=mass[3],
      PosteriorMixed=mass[4],
      ReinvasionGap=(mass[1]+mass[3])-mass[1],
      EffectiveParticles=1/sum(post^2), MatchingParticles=sum(match),
      stringsAsFactors=FALSE)
    if(ReturnParticles) wh[,rr] <- post
    if (rr < nrow(obs))
      w <- .ipof_rescale_origin_matched(prior,match,cls,post,tt)
    else w <- post
  }
  out <- list(Summary=do.call(rbind,rows), FinalParticleWeights=w,
              ObservationHistory=obs, OriginState=origin,
              ParticleWeightHistory=wh,
              Interpretation=paste(
                "Posterior PoF is current freedom (no original or later pathogen lineage).",
                "Original persistence and reintroduction are retained separately.",
                "Between rounds, only realised pathogen-detection histories compatible with the observations are propagated,",
                "with free/original-only/reintroduced-only/mixed masses rescaled to the likelihood posterior."))
  class(out) <- c("INApestPathogenPoFSequential","INApestPathogenPoF","list")
  out
}

INApestPathogenPoFSequentialCompartment <- function(ModelOutput, Pathogen, ObservationHistory,
                                                     PriorWeights=NULL, ObservationMatch=c("binary","exact"),
                                                     ReturnParticles=FALSE) {
  ModelOutput<-.ipof_model_output(ModelOutput); PathogenSpec<-.ipof_pathogen_spec(Pathogen)
  if(!(PathogenSpec$Model %in% c("SIS","SIR","SEIR")))stop("Sequential compartment PoF requires SIS, SIR or SEIR.")
  origin<-INApestPathogenOriginState(ModelOutput,PathogenSpec);D<-ModelOutput$PathogenDetectedResults
  if(is.null(D)||length(dim(D))!=3L)stop("PathogenDetectedResults nodes x timesteps x particles are required.")
  if(dim(D)[2]!=origin$Ntimesteps||dim(D)[3]!=origin$Nparticles)stop("PathogenDetectedResults dimensions do not match origin arrays.")
  ObservationMatch<-match.arg(ObservationMatch)
  if(!is.data.frame(ObservationHistory)||!all(c("Timestep","PathogenDetections")%in%names(ObservationHistory)))stop("ObservationHistory must contain Timestep and PathogenDetections.")
  obs<-ObservationHistory[order(ObservationHistory$Timestep),,drop=FALSE];obs$Timestep<-as.integer(obs$Timestep)
  if(!nrow(obs)||anyDuplicated(obs$Timestep)||any(obs$Timestep<1L|obs$Timestep>origin$Ntimesteps))stop("ObservationHistory must contain unique valid timesteps.")
  np<-origin$Nparticles;w<-if(is.null(PriorWeights))rep(1/np,np)else .ipof_norm_weights(PriorWeights);rows<-vector("list",nrow(obs));wh<-if(ReturnParticles)matrix(NA_real_,np,nrow(obs))else NULL
  for(rr in seq_len(nrow(obs))) {
    tt<-obs$Timestep[rr];prior<-w;cls<-origin$Class[tt,]
    L<-INApestPathogenObservationLikelihood(ModelOutput,PathogenSpec,tt,as.integer(obs$PathogenDetections[rr]>0));ev<-sum(prior*L);if(!is.finite(ev)||ev<=0)stop("Observed pathogen evidence has zero likelihood at timestep ",tt)
    post<-.ipof_norm_weights(prior*L);sim_count<-apply(D[,tt,,drop=FALSE],3,sum);match<-.ipof_match_binary_detection(sim_count,obs$PathogenDetections[rr],ObservationMatch);mass<-vapply(0:3,function(cc)sum(post[cls==cc]),numeric(1))
    rows[[rr]]<-data.frame(Round=rr,Timestep=tt,PathogenModel=PathogenSpec$Model,ObservationProbability=ev,PriorPoF=sum(prior[cls==0L]),PosteriorPoF=mass[1],PosteriorOriginalPersistence=mass[2]+mass[4],PosteriorOriginalEliminated=mass[1]+mass[3],PosteriorReintroductionPresent=mass[3]+mass[4],PosteriorReintroductionOnly=mass[3],PosteriorMixed=mass[4],ReinvasionGap=mass[3],EffectiveParticles=1/sum(post^2),MatchingParticles=sum(match),stringsAsFactors=FALSE)
    if(ReturnParticles)wh[,rr]<-post
    if(rr<nrow(obs))w<-.ipof_rescale_origin_matched(prior,match,cls,post,tt)else w<-post
  }
  out<-list(Summary=do.call(rbind,rows),FinalParticleWeights=w,ObservationHistory=obs,OriginState=origin,ParticleWeightHistory=wh,
            Interpretation=paste("Current freedom requires no active original or reintroduced pathogen lineage.",if(PathogenSpec$Model=="SEIR")"Exposed E hosts count as persistence even though default pathogen detection is based on I."else"Active I hosts determine persistence.","Only compatible realised detection histories are propagated between rounds; four origin-class masses are rescaled to the likelihood posterior."))
  class(out)<-c("INApestPathogenPoFSequentialCompartment","INApestPathogenPoFSequential","INApestPathogenPoF","list");out
}



###############################################################################
### Point / point-transition sequential compatible-history proof of freedom
###############################################################################

.ipof_point_detection_counts <- function(ModelOutput, Ntimesteps, Nparticles) {
  out<-matrix(0L,nrow=Ntimesteps,ncol=Nparticles)
  pe<-ModelOutput$PathogenEvents
  if(!is.data.frame(pe)||!nrow(pe)||!("pathogen_detected"%in%names(pe))) return(out)
  if(!all(c("perm","timestep")%in%names(pe))) stop("PathogenEvents must contain perm and timestep.")
  for(i in seq_len(nrow(pe))){pp<-as.integer(pe$perm[i]);tt<-as.integer(pe$timestep[i]);if(!is.na(pp)&&!is.na(tt)&&pp>=1L&&pp<=Nparticles&&tt>=1L&&tt<=Ntimesteps&&isTRUE(pe$pathogen_detected[i]))out[tt,pp]<-out[tt,pp]+1L}
  out
}

INApestPathogenPoFSequentialPoint <- function(ModelOutput, Pathogen, ObservationHistory,
                                               PriorWeights=NULL, ObservationMatch=c("binary","exact"),
                                               ReturnParticles=FALSE) {
  ModelOutput<-.ipof_model_output(ModelOutput)
  if(!inherits(Pathogen,"INApestPointPathogenInteraction")) stop("Sequential point PoF requires an INApestPointPathogenInteraction object.")
  ps<-.ipof_pathogen_spec(Pathogen);if(!isTRUE(ps$TrackOrigin))stop("Sequential point persistence/reintroduction PoF requires Pathogen TrackOrigin=TRUE.")
  origin<-INApestPathogenOriginState(ModelOutput,Pathogen);D<-.ipof_point_detection_counts(ModelOutput,origin$Ntimesteps,origin$Nparticles)
  ObservationMatch<-match.arg(ObservationMatch)
  if(!is.data.frame(ObservationHistory)||!all(c("Timestep","PathogenDetections")%in%names(ObservationHistory)))stop("ObservationHistory must contain Timestep and PathogenDetections.")
  obs<-ObservationHistory[order(ObservationHistory$Timestep),,drop=FALSE];obs$Timestep<-as.integer(obs$Timestep)
  if(!nrow(obs)||anyDuplicated(obs$Timestep)||any(obs$Timestep<1L|obs$Timestep>origin$Ntimesteps))stop("ObservationHistory must contain unique valid timesteps.")
  np<-origin$Nparticles;w<-if(is.null(PriorWeights))rep(1/np,np)else .ipof_norm_weights(PriorWeights);rows<-vector("list",nrow(obs));wh<-if(ReturnParticles)matrix(NA_real_,np,nrow(obs))else NULL
  for(rr in seq_len(nrow(obs))){
    tt<-obs$Timestep[rr];prior<-w;cls<-origin$Class[tt,]
    L<-INApestPathogenObservationLikelihood(ModelOutput,Pathogen,tt,as.integer(obs$PathogenDetections[rr]>0));ev<-sum(prior*L);if(!is.finite(ev)||ev<=0)stop("Observed pathogen evidence has zero likelihood at timestep ",tt)
    post<-.ipof_norm_weights(prior*L);match<-.ipof_match_binary_detection(D[tt,],obs$PathogenDetections[rr],ObservationMatch);mass<-vapply(0:3,function(cc)sum(post[cls==cc]),numeric(1))
    rows[[rr]]<-data.frame(Round=rr,Timestep=tt,PathogenModel=ps$Model,ObservationProbability=ev,PriorPoF=sum(prior[cls==0L]),PosteriorPoF=mass[1],PosteriorOriginalPersistence=mass[2]+mass[4],PosteriorOriginalEliminated=mass[1]+mass[3],PosteriorReintroductionPresent=mass[3]+mass[4],PosteriorReintroductionOnly=mass[3],PosteriorMixed=mass[4],ReinvasionGap=mass[3],PosteriorExpectedOriginalLineages=.ipof_weighted(origin$OriginalLineageCount[tt,],post),PosteriorExpectedReintroducedLineages=.ipof_weighted(origin$ReintroducedLineageCount[tt,],post),EffectiveParticles=1/sum(post^2),MatchingParticles=sum(match),stringsAsFactors=FALSE)
    if(ReturnParticles)wh[,rr]<-post
    if(rr<nrow(obs))w<-.ipof_rescale_origin_matched(prior,match,cls,post,tt)else w<-post
  }
  out<-list(Summary=do.call(rbind,rows),FinalParticleWeights=w,ObservationHistory=obs,OriginState=origin,ParticleWeightHistory=wh,
            Interpretation=paste("Point pathogen lineage roots are retained explicitly in PointHistory.","PoF collapses those roots to free/original-only/reintroduced-only/mixed for compatible-history rescaling while reporting expected numbers of surviving original and reintroduced lineages."))
  class(out)<-c("INApestPathogenPoFSequentialPoint","INApestPathogenPoFSequential","INApestPathogenPoF","list");out
}

###############################################################################
### Multi-target proof-of-freedom extension: pathogen + biocontrol
###
### This block is additive. The frozen pathogen PoF implementation above is
### retained unchanged. New functions allow the same particle weights to be
### updated from pathogen and biocontrol surveillance in the same timestep.
### Marginal pathogen PoF, marginal biocontrol PoF and their joint intersection
### are calculated from the same posterior particle distribution.
###############################################################################

.ipof_multi_model_output <- function(ModelOutput) {
  if (is.list(ModelOutput) && !is.null(ModelOutput$OutputDir) &&
      !is.null(ModelOutput$ModelName) &&
      is.null(ModelOutput$BiocontrolHistory) &&
      is.null(ModelOutput$PathogenStateResults) &&
      is.null(ModelOutput$PathogenStageResults) &&
      is.null(ModelOutput$PathogenPresentResults) &&
      is.null(ModelOutput$PointHistory)) {
    od <- ModelOutput$OutputDir
    if (length(od) != 1L || is.na(od)) od <- ""
    ModelOutput <- file.path(od, ModelOutput$ModelName)
  }

  if (is.list(ModelOutput)) return(ModelOutput)

  if (!is.character(ModelOutput) || length(ModelOutput) != 1L || is.na(ModelOutput))
    stop("ModelOutput must be a result object, an .rds result file, a standard output filename stem, or list(OutputDir=..., ModelName=...).")

  if (file.exists(ModelOutput) && grepl("\\.rds$", ModelOutput, ignore.case = TRUE)) {
    z <- readRDS(ModelOutput)
    if (!is.list(z)) stop("A directly supplied .rds ModelOutput must contain a saved result object.")
    return(z)
  }

  stem <- ModelOutput
  point_candidates <- c(
    paste0(stem, "_PointResults.rds"),
    paste0(stem, "_MetaPointParallelResults.rds"),
    paste0(stem, "_VertebratePointResults.rds"),
    paste0(stem, "_VertebratePointParallelResults.rds"),
    paste0(stem, "_PointTransitionMatrixParallelResults.rds")
  )
  hit <- point_candidates[file.exists(point_candidates)]
  if (length(hit)) return(readRDS(hit[1L]))

  out <- tryCatch(.ipof_model_output(stem), error = function(e) list())
  bc_file <- paste0(stem, "BiocontrolHistory.rds")
  if (file.exists(bc_file)) out$BiocontrolHistory <- readRDS(bc_file)
  if (!length(out)) stop("No standard pathogen or biocontrol outputs found for ModelOutput stem: ", stem)
  out$ModelNameStem <- stem
  out
}

.ipof_biocontrol_agents <- function(ModelOutput, BiocontrolAgents = NULL) {
  h <- ModelOutput$BiocontrolHistory
  if (is.null(h)) stop("ModelOutput does not contain BiocontrolHistory.")
  available <- if (is.list(h) && is.list(h$State)) names(h$State) else
    if (is.data.frame(h) && "agent" %in% names(h)) unique(as.character(h$agent)) else character()
  available <- available[nzchar(available)]
  if (!length(available)) stop("Could not identify biocontrol agent names from BiocontrolHistory.")
  if (is.null(BiocontrolAgents)) return(available)
  BiocontrolAgents <- as.character(BiocontrolAgents)
  bad <- setdiff(BiocontrolAgents, available)
  if (length(bad)) stop("Requested biocontrol agent(s) not present in history: ", paste(bad, collapse = ", "))
  unique(BiocontrolAgents)
}

.ipof_biocontrol_dims_from_point <- function(ModelOutput, h) {
  if (is.data.frame(ModelOutput$Summary) && nrow(ModelOutput$Summary) &&
      all(c("perm", "timestep") %in% names(ModelOutput$Summary))) {
    return(c(Ntimesteps = max(ModelOutput$Summary$timestep), Nparticles = max(ModelOutput$Summary$perm)))
  }
  if (nrow(h) && all(c("perm", "timestep") %in% names(h)))
    return(c(Ntimesteps = max(h$timestep), Nparticles = max(h$perm)))
  stop("Point BiocontrolHistory requires perm/timestep information, or ModelOutput$Summary with those fields.")
}

### Convert biocontrol history to the common latent freedom representation.
### Freedom means zero abundance across every stage of every selected agent.
### This biological definition is deliberately independent of which stages are
### observable in surveillance.
INApestBiocontrolFreedomState <- function(ModelOutput, BiocontrolAgents = NULL) {
  ModelOutput <- .ipof_multi_model_output(ModelOutput)
  h <- ModelOutput$BiocontrolHistory
  if (is.null(h)) stop("ModelOutput does not contain BiocontrolHistory.")
  agents <- .ipof_biocontrol_agents(ModelOutput, BiocontrolAgents)

  # Node-based Binary, Meta, MLU, Transition Matrix and Vertebrate Node.
  if (is.list(h) && is.list(h$State)) {
    total <- NULL
    base_dim <- NULL
    for (nm in agents) {
      x <- h$State[[nm]]
      d <- dim(x)
      if (length(d) != 4L)
        stop("Node biocontrol state for agent '", nm, "' must be nodes x agent-stage x timesteps x particles.")
      if (is.null(base_dim)) base_dim <- d[c(1,3,4)]
      if (!identical(as.integer(d[c(1,3,4)]), as.integer(base_dim)))
        stop("Selected biocontrol agents have incompatible node/time/particle dimensions.")
      a <- apply(x, c(3,4), sum)
      if (is.null(dim(a))) a <- matrix(a, nrow = d[3], ncol = d[4])
      if (is.null(total)) total <- a else total <- total + a
    }
    freedom <- total == 0
    return(list(
      Type = "node", Freedom = freedom, TotalAbundance = total,
      Nnodes = base_dim[1], Ntimesteps = base_dim[2], Nparticles = base_dim[3],
      Agents = agents
    ))
  }

  # Point-host engines store node-supported agent state in long form.
  if (is.data.frame(h)) {
    need <- c("perm", "timestep", "agent", "node", "stage", "abundance")
    if (!all(need %in% names(h)))
      stop("Point BiocontrolHistory must contain: ", paste(need, collapse = ", "))
    dd <- .ipof_biocontrol_dims_from_point(ModelOutput, h)
    nt <- as.integer(dd[["Ntimesteps"]]); np <- as.integer(dd[["Nparticles"]])
    total <- matrix(0, nrow = nt, ncol = np)
    z <- h[h$agent %in% agents, , drop = FALSE]
    if (nrow(z)) {
      for (i in seq_len(nrow(z))) {
        tt <- as.integer(z$timestep[i]); pp <- as.integer(z$perm[i])
        if (tt >= 1L && tt <= nt && pp >= 1L && pp <= np)
          total[tt, pp] <- total[tt, pp] + as.numeric(z$abundance[i])
      }
    }
    return(list(
      Type = "point", Freedom = total == 0, TotalAbundance = total,
      Ntimesteps = nt, Nparticles = np, Agents = agents
    ))
  }

  stop("Unsupported BiocontrolHistory representation.")
}

### Observation model for direct surveillance of a biocontrol agent.
### DetectionProb is per individual per surveillance round. It may be a scalar,
### node vector, timestep vector, nodes x timesteps matrix, named list by agent,
### or a resolver function(agent, stage, timestep, n_nodes, Ntimesteps).
### DetectStages controls observability only; all stages still count against PoF.
INApestBiocontrolSurveillance <- function(DetectionProb,
                                          DetectStages = NULL,
                                          DetectAgents = NULL) {
  if (is.null(DetectionProb)) stop("DetectionProb is required.")
  structure(list(
    DetectionProb = DetectionProb,
    DetectStages = DetectStages,
    DetectAgents = DetectAgents
  ), class = "INApestBiocontrolSurveillance")
}

.ipof_bio_stages <- function(stage_names, DetectStages, agent) {
  if (is.null(DetectStages)) return(stage_names)
  z <- DetectStages
  if (is.list(z)) {
    if (!is.null(names(z)) && agent %in% names(z)) z <- z[[agent]]
    else if (length(z) == 1L) z <- z[[1L]]
    else stop("DetectStages list must be named by agent, or contain one shared element.")
  }
  if (is.numeric(z)) {
    z <- as.integer(z)
    if (anyNA(z) || any(z < 1L | z > length(stage_names)))
      stop("Numeric DetectStages outside available stage indices for agent '", agent, "'.")
    return(stage_names[z])
  }
  z <- as.character(z)
  bad <- setdiff(z, stage_names)
  if (length(bad)) stop("DetectStages absent from agent '", agent, "': ", paste(bad, collapse = ", "))
  unique(z)
}

.ipof_bio_detection_prob <- function(x, agent, stage, timestep, n_nodes, Ntimesteps) {
  if (is.list(x) && !is.data.frame(x)) {
    if (!is.null(names(x)) && agent %in% names(x))
      return(.ipof_bio_detection_prob(x[[agent]], agent, stage, timestep, n_nodes, Ntimesteps))
    if (length(x) == 1L)
      return(.ipof_bio_detection_prob(x[[1L]], agent, stage, timestep, n_nodes, Ntimesteps))
    stop("DetectionProb list must be named by biocontrol agent, or contain one shared element.")
  }
  if (is.function(x)) {
    a <- list(agent = agent, stage = stage, timestep = timestep,
              n_nodes = n_nodes, Ntimesteps = Ntimesteps)
    fm <- names(formals(x))
    if (!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a), fm)]
    x <- do.call(x, a)
  }
  d <- dim(x)
  if (!is.null(d)) {
    if (length(d) == 2L && identical(as.integer(d), c(as.integer(n_nodes), as.integer(Ntimesteps))))
      return(.ipof_clip01(as.numeric(x[, timestep])))
    stop("Biocontrol DetectionProb array/matrix must be nodes x Ntimesteps; use a resolver function for stage-specific schedules.")
  }
  x <- as.numeric(x)
  if (!length(x) || any(!is.finite(x))) stop("Biocontrol DetectionProb must resolve to finite numeric values.")
  if (length(x) == 1L) return(rep(.ipof_clip01(x), n_nodes))
  if (length(x) == n_nodes && length(x) == Ntimesteps)
    stop("Ambiguous Biocontrol DetectionProb vector: n_nodes equals Ntimesteps. Use a matrix or resolver function.")
  if (length(x) == n_nodes) return(.ipof_clip01(x))
  if (length(x) == Ntimesteps) return(rep(.ipof_clip01(x[timestep]), n_nodes))
  stop("Biocontrol DetectionProb must be scalar, length nodes, length Ntimesteps, nodes x Ntimesteps, named list, or resolver function.")
}

### Likelihood for direct biocontrol surveillance at one timestep.
### Observation can be scalar 0/1 (system-wide) or node vector 0/1/NA.
INApestBiocontrolObservationLikelihood <- function(ModelOutput,
                                                    Surveillance,
                                                    timestep,
                                                    Observation = 0,
                                                    BiocontrolAgents = NULL) {
  ModelOutput <- .ipof_multi_model_output(ModelOutput)
  if (!inherits(Surveillance, "INApestBiocontrolSurveillance"))
    stop("Surveillance must be created by INApestBiocontrolSurveillance().")
  fs <- INApestBiocontrolFreedomState(ModelOutput, BiocontrolAgents)
  if (timestep < 1L || timestep > fs$Ntimesteps) stop("timestep outside model output.")

  agents <- fs$Agents
  if (!is.null(Surveillance$DetectAgents)) {
    detect_agents <- intersect(agents, as.character(Surveillance$DetectAgents))
    if (!length(detect_agents)) stop("No selected PoF agents are included in Surveillance$DetectAgents.")
  } else detect_agents <- agents

  h <- ModelOutput$BiocontrolHistory
  np <- fs$Nparticles

  if (is.list(h) && is.list(h$State)) {
    n_nodes <- fs$Nnodes
    no_node <- matrix(1, nrow = n_nodes, ncol = np)
    for (nm in detect_agents) {
      x <- h$State[[nm]]
      stage_names <- dimnames(x)[[2]]
      if (is.null(stage_names)) stage_names <- as.character(seq_len(dim(x)[2]))
      stages <- .ipof_bio_stages(stage_names, Surveillance$DetectStages, nm)
      for (st in stages) {
        jj <- match(st, stage_names)
        count <- matrix(x[, jj, timestep, , drop = FALSE], nrow = n_nodes, ncol = np)
        p <- .ipof_bio_detection_prob(Surveillance$DetectionProb, nm, st, timestep, n_nodes, fs$Ntimesteps)
        no_node <- no_node * sweep(count, 1, 1 - p, function(n, q) q^n)
      }
    }
  } else if (is.data.frame(h)) {
    z0 <- h[h$timestep == timestep & h$agent %in% detect_agents, , drop = FALSE]
    n_nodes <- if (nrow(h)) max(as.integer(h$node), na.rm = TRUE) else 0L
    if (!is.finite(n_nodes) || n_nodes < 1L) stop("Point BiocontrolHistory has no valid support-node indices.")
    no_node <- matrix(1, nrow = n_nodes, ncol = np)
    for (nm in detect_agents) {
      za <- z0[z0$agent == nm, , drop = FALSE]
      available_stages <- unique(as.character(h$stage[h$agent == nm]))
      stages <- .ipof_bio_stages(available_stages, Surveillance$DetectStages, nm)
      for (st in stages) {
        p <- .ipof_bio_detection_prob(Surveillance$DetectionProb, nm, st, timestep, n_nodes, fs$Ntimesteps)
        zz <- za[za$stage == st, , drop = FALSE]
        if (nrow(zz)) {
          for (i in seq_len(nrow(zz))) {
            node <- as.integer(zz$node[i]); pp <- as.integer(zz$perm[i])
            if (node >= 1L && node <= n_nodes && pp >= 1L && pp <= np)
              no_node[node, pp] <- no_node[node, pp] * (1 - p[node])^as.numeric(zz$abundance[i])
          }
        }
      }
    }
  } else stop("Unsupported BiocontrolHistory representation.")

  if (length(Observation) == 1L) {
    if (is.na(Observation) || !(Observation %in% c(0,1)))
      stop("Scalar Observation must be 0 (no detection) or 1 (one or more detections).")
    l0 <- apply(no_node, 2, prod)
    return(if (Observation == 0) l0 else 1 - l0)
  }
  if (length(Observation) != nrow(no_node)) stop("Node Observation must have one value per biocontrol support node.")
  if (any(!is.na(Observation) & !(Observation %in% c(0,1)))) stop("Node Observation values must be 0, 1 or NA.")
  out <- rep(1, np)
  for (i in seq_len(nrow(no_node))) if (!is.na(Observation[i]))
    out <- out * if (Observation[i] == 0) no_node[i, ] else (1 - no_node[i, ])
  out
}

.ipof_observation_at <- function(ObservationHistory, target, timestep, ntargets) {
  if (is.null(ObservationHistory)) return(NULL)
  h <- ObservationHistory
  if (is.list(h) && !is.null(names(h))) {
    if (target %in% names(h)) h <- h[[target]] else return(NULL)
  } else if (ntargets > 1L)
    stop("With multiple PoF targets, ObservationHistory must be a named list with 'pathogen' and/or 'biocontrol'.")
  if (is.null(h)) return(NULL)
  if (is.list(h)) {
    if (length(h) >= timestep) return(h[[timestep]])
    return(NULL)
  }
  if ((is.numeric(h) || is.logical(h)) && length(h) >= timestep) return(h[timestep])
  NULL
}

.ipof_metric <- function(w, f) sum(w * as.numeric(f))

### General proof-of-freedom wrapper.
###
### Targets may contain "pathogen", "biocontrol", or both. When both are used,
### the joint PoF is the posterior probability that the pathogen is absent AND
### the selected biocontrol agent(s) are absent in the same particle. It is not
### calculated as the product of the two marginal PoFs.
###
### When both observation streams are supplied for one timestep, their
### likelihoods are multiplied conditional on the simulated latent state and the
### particle weights are normalised once. This makes the update order-independent
### under the stated conditional-independence observation assumption.
INApestPoF <- function(ModelOutput,
                       Pathogen = NULL,
                       Biocontrol = NULL,
                       Targets = NULL,
                       ObservationHistory = NULL,
                       PriorWeights = NULL,
                       BiocontrolSurveillance = NULL,
                       BiocontrolAgents = NULL,
                       ReplayRequired = NULL) {
  ModelOutput <- .ipof_multi_model_output(ModelOutput)
  if (is.null(Targets)) {
    Targets <- c(if (!is.null(Pathogen)) "pathogen", if (!is.null(Biocontrol) || !is.null(ModelOutput$BiocontrolHistory)) "biocontrol")
  }
  Targets <- unique(match.arg(Targets, c("pathogen", "biocontrol"), several.ok = TRUE))
  if (!length(Targets)) stop("At least one target is required.")
  if ("pathogen" %in% Targets && is.null(Pathogen)) stop("Pathogen is required when target 'pathogen' is selected.")
  if ("biocontrol" %in% Targets && is.null(ModelOutput$BiocontrolHistory)) stop("BiocontrolHistory is required when target 'biocontrol' is selected.")

  states <- list()
  if ("pathogen" %in% Targets) states$pathogen <- INApestPathogenFreedomState(ModelOutput, Pathogen)
  if ("biocontrol" %in% Targets) states$biocontrol <- INApestBiocontrolFreedomState(ModelOutput, BiocontrolAgents)

  nt <- states[[1L]]$Ntimesteps; np <- states[[1L]]$Nparticles
  for (nm in names(states)) {
    if (states[[nm]]$Ntimesteps != nt || states[[nm]]$Nparticles != np)
      stop("All selected PoF targets must have the same timestep and particle dimensions.")
  }

  w <- if (is.null(PriorWeights)) rep(1 / np, np) else .ipof_norm_weights(PriorWeights)
  if (length(w) != np) stop("PriorWeights must contain one value per particle.")

  if (is.null(ReplayRequired)) {
    ReplayRequired <- FALSE
    if ("pathogen" %in% Targets) {
      ps <- .ipof_pathogen_spec(Pathogen)
      ReplayRequired <- isTRUE(ps$DetectionTriggersInfo)
    }
  }

  cols <- c("timestep", "PriorPoF_Pathogen", "PosteriorPoF_Pathogen",
            "PriorPoF_Biocontrol", "PosteriorPoF_Biocontrol",
            "PriorPoF_Joint", "PosteriorPoF_Joint",
            "ObservationEvidence", "EffectiveParticles")
  ans <- as.data.frame(matrix(NA_real_, nrow = nt, ncol = length(cols)))
  names(ans) <- cols; ans$timestep <- seq_len(nt)

  weights_by_time <- vector("list", nt)
  likelihood_by_time <- vector("list", nt)

  for (tt in seq_len(nt)) {
    fp <- if ("pathogen" %in% Targets) states$pathogen$Freedom[tt, ] else rep(TRUE, np)
    fb <- if ("biocontrol" %in% Targets) states$biocontrol$Freedom[tt, ] else rep(TRUE, np)
    fj <- fp & fb

    if ("pathogen" %in% Targets) ans$PriorPoF_Pathogen[tt] <- .ipof_metric(w, fp)
    if ("biocontrol" %in% Targets) ans$PriorPoF_Biocontrol[tt] <- .ipof_metric(w, fb)
    ans$PriorPoF_Joint[tt] <- .ipof_metric(w, fj)

    L <- rep(1, np)
    components <- list()
    if ("pathogen" %in% Targets) {
      obs <- .ipof_observation_at(ObservationHistory, "pathogen", tt, length(Targets))
      lp <- if (is.null(obs) || (length(obs) == 1L && is.na(obs))) rep(1, np) else
        INApestPathogenObservationLikelihood(ModelOutput, Pathogen, tt, obs)
      components$pathogen <- lp; L <- L * lp
    }
    if ("biocontrol" %in% Targets) {
      obs <- .ipof_observation_at(ObservationHistory, "biocontrol", tt, length(Targets))
      lb <- if (is.null(obs) || (length(obs) == 1L && is.na(obs))) rep(1, np) else {
        if (is.null(BiocontrolSurveillance)) stop("BiocontrolSurveillance is required when biocontrol observations are supplied.")
        INApestBiocontrolObservationLikelihood(ModelOutput, BiocontrolSurveillance, tt, obs, BiocontrolAgents)
      }
      components$biocontrol <- lb; L <- L * lb
    }

    ev <- sum(w * L)
    w <- .ipof_norm_weights(w * L)
    if ("pathogen" %in% Targets) ans$PosteriorPoF_Pathogen[tt] <- .ipof_metric(w, fp)
    if ("biocontrol" %in% Targets) ans$PosteriorPoF_Biocontrol[tt] <- .ipof_metric(w, fb)
    ans$PosteriorPoF_Joint[tt] <- .ipof_metric(w, fj)
    ans$ObservationEvidence[tt] <- ev
    ans$EffectiveParticles[tt] <- 1 / sum(w^2)
    weights_by_time[[tt]] <- w
    likelihood_by_time[[tt]] <- components
  }

  structure(list(
    Summary = ans,
    PosteriorWeights = w,
    PosteriorWeightsByTime = weights_by_time,
    ObservationLikelihoodByTime = likelihood_by_time,
    Freedom = lapply(states, `[[`, "Freedom"),
    TargetState = states,
    Targets = Targets,
    BiocontrolAgents = if ("biocontrol" %in% Targets) states$biocontrol$Agents else NULL,
    ReplayRequired = ReplayRequired,
    ObservationAssumption = if (length(Targets) > 1L)
      "Pathogen and biocontrol observation streams are conditionally independent given the simulated latent state; same-timestep likelihoods are multiplied before one weight normalisation."
      else "Single-target observation likelihood weighting.",
    Interpretation = if (ReplayRequired)
      "Current-state likelihoods are valid, but future trajectories require compatible-history replay when observations trigger later information or management."
      else "Post-hoc sequential likelihood weighting is valid when supplied observations do not alter later simulated dynamics."
  ), class = c("INApestPoF", "list"))
}

### Convenience wrapper for biocontrol-only PoF.
INApestBiocontrolPoF <- function(ModelOutput,
                                 Biocontrol = NULL,
                                 ObservationHistory = NULL,
                                 PriorWeights = NULL,
                                 Surveillance = NULL,
                                 BiocontrolAgents = NULL,
                                 ReplayRequired = FALSE) {
  INApestPoF(
    ModelOutput = ModelOutput, Biocontrol = Biocontrol,
    Targets = "biocontrol", ObservationHistory = ObservationHistory,
    PriorWeights = PriorWeights, BiocontrolSurveillance = Surveillance,
    BiocontrolAgents = BiocontrolAgents, ReplayRequired = ReplayRequired
  )
}
