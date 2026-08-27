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
    NewInfectionResults = "NewInfectionLargeOut.rds"
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
    np <- if(nrow(h)) max(h$perm) else 0L; nt <- if(nrow(h)) max(h$timestep) else 0L
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
