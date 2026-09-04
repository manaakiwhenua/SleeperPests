###############################################################################
### INApestMetaParallelMultipleLandUse -- parallel node x land-use engine
###
### Runs the Multiple Land Use architecture across independent stochastic
### permutations using a common PSOCK/serial parallel contract. Each worker keeps
### the same node x land-use demography, dispersal, information, management and
### optional pathogen processes as the serial engine.
###
### Parallelisation changes execution only; biological semantics are unchanged.
###############################################################################

# Default node x land-use local growth, dispersal and recruitment process.
local.dynamicsLU = function(sddprob = SDDprob, nodepropaguleproduction = NodePropaguleProduction,nodeenvestabprob = NodeEnvEstabProb,n=N0,
lddprob = LDDprob, lddrate = LDDrate,k_is_0 = K_is_0, nodeK = NodeK,nodepropaguleestablishment = NodePropaguleEstablishment,
nodespreadreduction = NodeSpreadReduction,nodefecundityreduction = NodeFecundityReduction,managing = Managing)
{
  # Keep the default LocalDynamics self-contained so it can be serialized to
  # PSOCK workers without relying on objects in the master's global environment.
  AllocateLandUseRecruits <- function(X)
    {
    Vect = vector(length=0)
    for(i in 2:length(X))
      Vect <- c(Vect,rep((i-1),times = X[i]))
    Sample <- sample(Vect,size = X[1],replace = F)
    Out <- vector(length = length(X)-1)
    for(j in 1:length(Out))
      Out[j] <- length(Sample[Sample==j])
    return(Out)
    }
FecundityContrib <- n * (1-nodefecundityreduction*managing)
Propagules <- rpois(nrow(sddprob), nodepropaguleproduction * rowSums(FecundityContrib))# propagules are produced after management effects on fecundity
      # integer per-propagule thinning/allocation.
      # Connectivity entries are applied once; 1-rowSum means no represented target.
      n_nodes <- nrow(sddprob)
      # Allocate realised integer propagules among destinations and outside loss.
      AllocateDispersers <- function(counts, probmat, name) {
        incoming <- numeric(n_nodes)
        if(!is.matrix(probmat)) return(incoming)
        rs <- rowSums(probmat)
        if(any(!is.finite(probmat)) || any(probmat < 0) || any(rs > 1 + 1e-12))
          stop(name, " must contain finite non-negative probabilities with every source-row sum <= 1")
        active <- which(counts > 0 & rs > 0)
        if(!length(active)) return(incoming)
        for(ii in active) {
          probs <- c(probmat[ii,],max(0,1-rs[ii]))
          z <- rmultinom(1,size=counts[ii],prob=probs)
          incoming <- incoming + z[seq_len(n_nodes),1]
        }
        incoming
      }
      if(lddrate <= 0) LDDcount <- integer(n_nodes)
      else if(lddrate >= 1) LDDcount <- Propagules
      else LDDcount <- rbinom(n_nodes,size=Propagules,prob=lddrate)
      SDDcount <- Propagules-LDDcount
      Pin <- AllocateDispersers(SDDcount,sddprob,"SDDprob")
      Qin <- numeric(n_nodes)
      if(is.matrix(lddprob)) {
        contrib_sum <- rowSums(FecundityContrib)
        Pn <- FecundityContrib/contrib_sum
        Pn[!is.finite(Pn)] <- 0
        keep <- rowSums(Pn*(1-nodespreadreduction*managing))
        keep[contrib_sum <= 0] <- 0
        keep <- pmin(1,pmax(0,keep))
        LDDkept <- rbinom(n_nodes,size=LDDcount,prob=keep)
        Qin <- AllocateDispersers(LDDkept,lddprob,"LDDprob")
      }
     
    # propagule success depends on availability of uninfested host plants
    Recruits <- rbinom(nrow(sddprob), rowSums(nodeK)-rowSums(n), 1 - exp(-nodepropaguleestablishment*nodeenvestabprob*(Pin+Qin)))  
    InVector = cbind(Recruits,nodeK-n)
    LUrecruits <- t(apply(InVector,1,FUN = AllocateLandUseRecruits))
    # Possibly due to rounding error when K is large N0 can occasionally = K+1
    Nout <- ifelse(n + LUrecruits>nodeK,nodeK,n + LUrecruits)
return(Nout)
}


# -----------------------------------------------------------------------------
# Optional custom LocalDynamics arguments
# -----------------------------------------------------------------------------
# Ordinary LocalDynamicsArgs entries are passed unchanged. Wrap a vector,
# matrix, array, or list in INApestLocalDynamicsTimeArg() when its final
# dimension/list position is indexed by simulation timestep. This avoids
# guessing whether an arbitrary custom matrix is static or time-varying.
INApestLocalDynamicsTimeArg <- function(x) {
  structure(list(values = x), class = "INApestLocalDynamicsTimeArg")
}

# Resolve custom LocalDynamics arguments for the current timestep.
.resolve_INApest_LocalDynamicsArgs <- function(LocalDynamicsArgs, timestep, Ntimesteps) {
  if (is.null(LocalDynamicsArgs)) LocalDynamicsArgs <- list()
  if (!is.list(LocalDynamicsArgs))
    stop("LocalDynamicsArgs must be a named list")
  if (!length(LocalDynamicsArgs)) return(list())
  if (is.null(names(LocalDynamicsArgs)) || any(!nzchar(names(LocalDynamicsArgs))))
    stop("Every LocalDynamicsArgs entry must have a non-empty name")
  if (anyDuplicated(names(LocalDynamicsArgs)))
    stop("LocalDynamicsArgs names must be unique")

  resolve_one <- function(x, name) {
    if (inherits(x, "INApestLocalDynamicsTimeArg")) {
      values <- x$values
      if (is.list(values)) {
        if (length(values) != Ntimesteps)
          stop(name, " wrapped list must have length Ntimesteps")
        return(values[[timestep]])
      }
      d <- dim(values)
      if (is.null(d)) {
        if (length(values) != Ntimesteps)
          stop(name, " wrapped vector must have length Ntimesteps")
        return(values[[timestep]])
      }
      if (tail(d, 1L) != Ntimesteps)
        stop(name, " wrapped array/matrix must have Ntimesteps in its final dimension")
      index <- lapply(d, seq_len)
      index[[length(d)]] <- timestep
      return(do.call(`[`, c(list(values), index, list(drop = TRUE))))
    }

    # Resolver functions are evaluated by the parent model. The returned
    # current value, not timestep itself, is passed to LocalDynamics.
    if (is.function(x)) {
      fm <- names(formals(x))
      call_args <- list(timestep = timestep, Ntimesteps = Ntimesteps)
      if (!is.null(fm) && !("..." %in% fm))
        call_args <- call_args[intersect(names(call_args), fm)]
      return(do.call(x, call_args))
    }

    x
  }

  out <- Map(resolve_one, LocalDynamicsArgs, names(LocalDynamicsArgs))
  names(out) <- names(LocalDynamicsArgs)
  out
}

###############################################################################
### Generic pathogen-state component for INApest abundance models
###############################################################################

INApestPathogen <- function(
    Model = c("SIS", "SIR", "SEIR", "Binary"),  # Pathogen state model: Binary, SIS, SIR or SEIR
    Beta = 0,                                   # Transmission-rate parameter
    RecoveryProb = 0,                           # Probability an infectious host recovers
    ProgressionProb = 1,                        # Probability an exposed host becomes infectious
    PathogenMortalityProb = 0,                  # Probability infection kills an infectious host
    ImmunityLossProb = 0,                       # Probability a recovered host becomes susceptible again
    InitialInfected = 0,                        # Starting infectious host count/specification
    InitialExposed = 0,                         # Starting exposed host count/specification
    InitialRecovered = 0,                       # Starting recovered host count/specification
    IntroductionProb = 0,                       # Probability an existing susceptible host is externally infected
    IntroductionNumber = 1,                     # Hosts infected when an external introduction occurs
    ContactMatrix = NULL,                       # Source-to-target pathogen contact matrix
    Transmission = c("frequency", "density"),   # Frequency- or density-dependent transmission
    DensityScale = 1,                           # Scale for density-dependent transmission
    InitialPresent = 0,                         # Starting binary pathogen occupancy
    ClearanceProb = 0,                          # Probability pathogen clears without host loss
    TransmissionProb = NULL,                    # Binary source-to-target transmission probability
    PathogenHostExtinctionProb = 0,             # Probability pathogen removes infected host occupancy
    DetectionProb = 0,                          # Background-surveillance detection probability
    DetectionTriggersInfo = FALSE) {            # Whether pathogen detection updates INApest information

  Model <- match.arg(Model)
  Transmission <- match.arg(Transmission)
  if (!is.logical(DetectionTriggersInfo) || length(DetectionTriggersInfo) != 1L || is.na(DetectionTriggersInfo))
    stop("DetectionTriggersInfo must be TRUE or FALSE")

  if (Model == "Binary") {
    if (is.null(TransmissionProb)) TransmissionProb <- Beta
    clip01 <- function(x) pmin(1, pmax(0, x))
    resolve_binary <- function(x, timestep, n_nodes, Ntimesteps, name) {
      if (is.function(x)) {
        fm <- names(formals(x)); a <- list(timestep=timestep,n_nodes=n_nodes,Ntimesteps=Ntimesteps)
        if (!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a),fm)]
        x <- do.call(x,a)
      }
      d <- dim(x)
      if (!is.null(d)) {
        if (length(d)==2L && all(d==c(n_nodes,Ntimesteps))) return(as.numeric(x[,timestep]))
        stop(name," must be scalar, length nodes, nodes x Ntimesteps, or resolver function")
      }
      x <- as.numeric(x)
      if (length(x)==1L) return(rep(x,n_nodes))
      if (length(x)==n_nodes) return(x)
      if (length(x)==Ntimesteps && Ntimesteps != n_nodes) return(rep(x[timestep],n_nodes))
      stop(name," has unsupported or ambiguous shape")
    }
    contact_binary <- function(timestep,n_nodes,Ntimesteps) {
      x <- ContactMatrix
      if (is.null(x)) return(diag(n_nodes))
      if (is.function(x)) {
        fm <- names(formals(x)); a <- list(timestep=timestep,n_nodes=n_nodes,Ntimesteps=Ntimesteps)
        if (!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a),fm)]
        x <- do.call(x,a)
      }
      d <- dim(x)
      if (length(d)==2L && all(d==c(n_nodes,n_nodes))) out <- x
      else if (length(d)==3L && all(d==c(n_nodes,n_nodes,Ntimesteps))) out <- x[,,timestep,drop=TRUE]
      else stop("Binary ContactMatrix must be nodes x nodes or nodes x nodes x Ntimesteps")
      if (any(!is.finite(out)) || any(out < 0)) stop("ContactMatrix entries must be finite and non-negative")
      as.matrix(out)
    }
    validate_binary <- function(n_nodes,Ntimesteps) {
      for (tt in seq_len(Ntimesteps)) {
        for (nm in c("RecoveryProb","PathogenMortalityProb","ImmunityLossProb")) {
          z <- resolve_binary(switch(nm, RecoveryProb=RecoveryProb, PathogenMortalityProb=PathogenMortalityProb, ImmunityLossProb=ImmunityLossProb),tt,n_nodes,Ntimesteps,nm)
          if (any(z != 0)) stop(nm," is not used by Model = 'Binary'; use ClearanceProb or PathogenHostExtinctionProb as appropriate")
        }
        for (nm in c("TransmissionProb","ClearanceProb","IntroductionProb","PathogenHostExtinctionProb","DetectionProb")) {
          x <- switch(nm, TransmissionProb=TransmissionProb, ClearanceProb=ClearanceProb, IntroductionProb=IntroductionProb, PathogenHostExtinctionProb=PathogenHostExtinctionProb, DetectionProb=DetectionProb)
          z <- resolve_binary(x,tt,n_nodes,Ntimesteps,nm)
          if (any(!is.finite(z)) || any(z < 0 | z > 1)) stop(nm," must resolve to [0,1]")
        }
        contact_binary(tt,n_nodes,Ntimesteps)
      }
      for (nm in c("InitialInfected","InitialExposed","InitialRecovered")) {
        z <- resolve_binary(switch(nm, InitialInfected=InitialInfected, InitialExposed=InitialExposed, InitialRecovered=InitialRecovered),1L,n_nodes,Ntimesteps,nm)
        if (any(z != 0)) stop(nm," is not used by Model = 'Binary'; use InitialPresent")
      }
      invisible(TRUE)
    }
    initial_binary <- function(Invaded,n_nodes=length(Invaded),Ntimesteps=1L) {
      validate_binary(n_nodes,Ntimesteps)
      z <- resolve_binary(InitialPresent,1L,n_nodes,Ntimesteps,"InitialPresent")
      if (any(!is.finite(z)) || any(z < 0 | z > 1)) stop("InitialPresent must resolve to [0,1]")
      if (all(z %in% c(0,1))) out <- as.integer(z) else out <- rbinom(n_nodes,1,z)
      as.integer(out * as.integer(Invaded > 0))
    }
    detect_binary <- function(PathogenPresent, timestep, Ntimesteps) {
      n_nodes <- length(PathogenPresent)
      p <- resolve_binary(DetectionProb, timestep, n_nodes, Ntimesteps, "DetectionProb")
      if (any(!is.finite(p)) || any(p < 0 | p > 1)) stop("DetectionProb must resolve to [0,1]")
      as.integer(stats::rbinom(n_nodes, 1L, clip01(p) * as.integer(PathogenPresent > 0)))
    }
    step_binary <- function(PathogenPresent,Invaded,timestep,Ntimesteps) {
      n_nodes <- length(Invaded); validate_binary(n_nodes,Ntimesteps)
      P <- as.integer(PathogenPresent > 0) * as.integer(Invaded > 0)
      tr <- resolve_binary(TransmissionProb,timestep,n_nodes,Ntimesteps,"TransmissionProb")
      clear <- resolve_binary(ClearanceProb,timestep,n_nodes,Ntimesteps,"ClearanceProb")
      intro <- resolve_binary(IntroductionProb,timestep,n_nodes,Ntimesteps,"IntroductionProb")
      ext <- resolve_binary(PathogenHostExtinctionProb,timestep,n_nodes,Ntimesteps,"PathogenHostExtinctionProb")
      C <- contact_binary(timestep,n_nodes,Ntimesteps)
      q <- sweep(C,1,tr*P,"*")
      q[q>1] <- 1
      ptrans <- 1-apply(1-q,2,prod)
      susceptible <- which(Invaded==1 & P==0)
      if (length(susceptible)) P[susceptible] <- rbinom(length(susceptible),1,clip01(ptrans[susceptible]))
      infected <- which(P==1)
      if (length(infected)) P[infected] <- P[infected] * rbinom(length(infected),1,1-clear[infected])
      candidates <- which(Invaded==1 & P==0)
      if (length(candidates)) P[candidates] <- pmax(P[candidates],rbinom(length(candidates),1,intro[candidates]))
      host_extinct <- integer(n_nodes); infected <- which(P==1)
      if (length(infected)) host_extinct[infected] <- rbinom(length(infected),1,ext[infected])
      Invaded[host_extinct==1] <- 0L; P[Invaded==0] <- 0L
      list(PathogenPresent=as.integer(P),Invaded=as.integer(Invaded),HostExtinction=host_extinct,NewPathogen=as.integer(P & !PathogenPresent))
    }
    out <- list(Model="Binary", States=c("Absent","Present"), TransmissionProb=TransmissionProb, ClearanceProb=ClearanceProb,
                InitialPresent=InitialPresent, IntroductionProb=IntroductionProb, ContactMatrix=ContactMatrix,
                PathogenHostExtinctionProb=PathogenHostExtinctionProb, DetectionProb=DetectionProb, DetectionTriggersInfo=DetectionTriggersInfo,
                binary_initial=initial_binary, binary_step=step_binary, binary_detect=detect_binary, binary_validate=validate_binary)
    class(out) <- c("INApestPathogen","list")
    return(out)
  }

  States <- switch(Model,
                   SIS = c("S", "I"),
                   SIR = c("S", "I", "R"),
                   SEIR = c("S", "E", "I", "R"))

  Spec <- list(
    Model = Model, Beta = Beta, RecoveryProb = RecoveryProb,
    ProgressionProb = ProgressionProb,
    PathogenMortalityProb = PathogenMortalityProb,
    ImmunityLossProb = ImmunityLossProb,
    InitialInfected = InitialInfected,
    InitialExposed = InitialExposed,
    InitialRecovered = InitialRecovered,
    IntroductionProb = IntroductionProb,
    IntroductionNumber = IntroductionNumber,
    ContactMatrix = ContactMatrix,
    Transmission = Transmission,
    DensityScale = DensityScale,
    DetectionProb = DetectionProb,
    DetectionTriggersInfo = DetectionTriggersInfo
  )

  clip01 <- function(x) pmin(1, pmax(0, x))

  context_dims <- function(context) {
    if (is.null(context$n_nodes) || is.null(context$Ntimesteps))
      stop("Pathogen context requires n_nodes and Ntimesteps")
    n_nodes <- as.integer(context$n_nodes)
    nts <- as.integer(context$Ntimesteps)
    n_lu <- if (is.null(context$n_landuses)) NULL else as.integer(context$n_landuses)
    list(n_nodes = n_nodes, Ntimesteps = nts, n_landuses = n_lu,
         n_units = if (is.null(n_lu)) n_nodes else n_nodes * n_lu)
  }

  resolve <- function(x, timestep, context, name) {
    dd <- context_dims(context)
    if (is.function(x)) {
      fm <- names(formals(x))
      args <- list(timestep = timestep, n_nodes = dd$n_nodes,
                   n_landuses = dd$n_landuses, Ntimesteps = dd$Ntimesteps)
      if (!is.null(fm) && !("..." %in% fm)) args <- args[intersect(names(args), fm)]
      x <- do.call(x, args)
    }

    d <- dim(x)
    if (is.null(dd$n_landuses)) {
      if (!is.null(d)) {
        if (length(d) == 2L && all(d == c(dd$n_nodes, dd$Ntimesteps)))
          return(as.numeric(x[, timestep]))
        stop(name, " matrix must have dimensions nodes x Ntimesteps")
      }
      x <- as.numeric(x)
      if (length(x) == 1L) return(rep(x, dd$n_units))
      if (length(x) == dd$n_nodes) return(x)
      if (length(x) == dd$Ntimesteps && dd$Ntimesteps != dd$n_nodes)
        return(rep(x[timestep], dd$n_units))
      stop(name, " must be scalar, length nodes, length Ntimesteps (when unambiguous), nodes x Ntimesteps, or a resolver function")
    }

    # Multiple-land-use context. Matrix/array forms are preferred when lengths
    # are ambiguous (e.g. number of nodes equals number of land uses).
    if (!is.null(d)) {
      if (length(d) == 2L && all(d == c(dd$n_nodes, dd$n_landuses)))
        return(as.numeric(x))
      if (length(d) == 3L && all(d == c(dd$n_nodes, dd$n_landuses, dd$Ntimesteps)))
        return(as.numeric(x[, , timestep]))
      stop(name, " for MLU must be nodes x land uses or nodes x land uses x Ntimesteps")
    }
    x <- as.numeric(x)
    if (length(x) == 1L) return(rep(x, dd$n_units))
    if (length(x) == dd$n_units) return(x)
    if (length(x) == dd$n_landuses && dd$n_landuses != dd$n_nodes)
      return(as.numeric(matrix(rep(x, each = dd$n_nodes), nrow = dd$n_nodes)))
    if (length(x) == dd$n_nodes && dd$n_nodes != dd$n_landuses)
      return(as.numeric(matrix(rep(x, dd$n_landuses), nrow = dd$n_nodes)))
    if (length(x) == dd$Ntimesteps && dd$Ntimesteps != dd$n_nodes && dd$Ntimesteps != dd$n_landuses)
      return(rep(x[timestep], dd$n_units))
    stop(name, " has an ambiguous or unsupported MLU shape; use a nodes x land uses matrix or nodes x land uses x Ntimesteps array")
  }

  contact_matrix <- function(timestep, context) {
    dd <- context_dims(context)
    x <- Spec$ContactMatrix
    if (is.null(x)) return(diag(dd$n_units))
    if (is.function(x)) {
      fm <- names(formals(x))
      args <- list(timestep = timestep, n_nodes = dd$n_nodes,
                   n_landuses = dd$n_landuses, Ntimesteps = dd$Ntimesteps)
      if (!is.null(fm) && !("..." %in% fm)) args <- args[intersect(names(args), fm)]
      x <- do.call(x, args)
    }
    d <- dim(x)
    if (is.null(d)) stop("ContactMatrix must be a square matrix, a 3-D source x target x timestep array, NULL, or a resolver function")
    if (length(d) == 2L) {
      if (!all(d == c(dd$n_units, dd$n_units)))
        stop("ContactMatrix must have dimensions pathogen units x pathogen units")
      out <- x
    } else if (length(d) == 3L) {
      if (!all(d == c(dd$n_units, dd$n_units, dd$Ntimesteps)))
        stop("Time-varying ContactMatrix must have dimensions pathogen units x pathogen units x Ntimesteps")
      out <- x[, , timestep, drop = TRUE]
    } else stop("ContactMatrix must be a matrix or 3-D array")
    if (any(!is.finite(out)) || any(out < 0))
      stop("ContactMatrix entries must be finite and non-negative")
    as.matrix(out)
  }

  validate_context <- function(context) {
    dd <- context_dims(context)
    prob_names <- c("RecoveryProb", "PathogenMortalityProb", "ImmunityLossProb", "IntroductionProb", "DetectionProb")
    if (Model == "SEIR") prob_names <- c(prob_names, "ProgressionProb")
    for (tt in seq_len(dd$Ntimesteps)) {
      contact_matrix(tt, context)
      for (nm in prob_names) {
        z <- resolve(Spec[[nm]], tt, context, nm)
        if (any(!is.finite(z)) || any(z < 0 | z > 1))
          stop(nm, " must resolve to values in [0,1] at every timestep")
      }
      rec <- resolve(Spec$RecoveryProb, tt, context, "RecoveryProb")
      mort <- resolve(Spec$PathogenMortalityProb, tt, context, "PathogenMortalityProb")
      if (any(rec + mort > 1 + 1e-12))
        stop("RecoveryProb + PathogenMortalityProb must not exceed 1 at any timestep")
      b <- resolve(Spec$Beta, tt, context, "Beta")
      if (any(!is.finite(b)) || any(b < 0))
        stop("Beta must resolve to finite non-negative values at every timestep")
      ds <- resolve(Spec$DensityScale, tt, context, "DensityScale")
      if (any(!is.finite(ds)) || any(ds <= 0))
        stop("DensityScale must resolve to positive values at every timestep")
      intro_n <- resolve(Spec$IntroductionNumber, tt, context, "IntroductionNumber")
      if (any(!is.finite(intro_n)) || any(intro_n < 0) || any(intro_n != floor(intro_n)))
        stop("IntroductionNumber must resolve to non-negative whole numbers at every timestep")
    }

    # Fail rather than silently ignore state-specific initial conditions or
    # processes that the selected compartment model cannot represent.
    if (!("E" %in% States)) {
      z <- resolve(Spec$InitialExposed, 1L, context, "InitialExposed")
      if (any(z != 0)) stop("InitialExposed must be zero unless Model = 'SEIR'")
    }
    if (!("R" %in% States)) {
      z <- resolve(Spec$InitialRecovered, 1L, context, "InitialRecovered")
      if (any(z != 0)) stop("InitialRecovered must be zero unless the model contains R")
      for (tt in seq_len(dd$Ntimesteps)) {
        z <- resolve(Spec$ImmunityLossProb, tt, context, "ImmunityLossProb")
        if (any(z != 0)) stop("ImmunityLossProb must be zero unless the model contains R")
      }
    }
    invisible(dd)
  }

  initial <- function(N, context) {
    validate_context(context)
    Nvec <- as.integer(as.numeric(N))
    resolve_count <- function(x, name) {
      z <- resolve(x, 1L, context, name)
      if (any(!is.finite(z)) || any(z < 0) || any(z != floor(z)))
        stop(name, " must resolve to non-negative whole-number host counts")
      as.integer(z)
    }
    I <- resolve_count(Spec$InitialInfected, "InitialInfected")
    E <- if ("E" %in% States) resolve_count(Spec$InitialExposed, "InitialExposed") else integer(length(Nvec))
    R <- if ("R" %in% States) resolve_count(Spec$InitialRecovered, "InitialRecovered") else integer(length(Nvec))
    if (any(I + E + R > Nvec)) stop("Initial pathogen-state counts exceed total host abundance")
    out <- matrix(0L, nrow = length(Nvec), ncol = length(States), dimnames = list(NULL, States))
    out[, "S"] <- Nvec - I - E - R
    out[, "I"] <- I
    if ("E" %in% States) out[, "E"] <- E
    if ("R" %in% States) out[, "R"] <- R
    out
  }

  thin_row <- function(counts, target) {
    counts <- as.integer(counts); target <- as.integer(target); total <- sum(counts)
    if (target >= total) return(counts)
    if (target <= 0L) return(integer(length(counts)))
    if (length(counts) == 1L) return(target)
    out <- integer(length(counts)); remain_draw <- target; remain_total <- total
    for (j in seq_len(length(counts) - 1L)) {
      out[j] <- rhyper(1L, counts[j], remain_total - counts[j], remain_draw)
      remain_draw <- remain_draw - out[j]
      remain_total <- remain_total - counts[j]
    }
    out[length(counts)] <- remain_draw
    out
  }

  reconcile <- function(State, N, context = NULL) {
    Nvec <- as.integer(as.numeric(N)); old <- rowSums(State)
    if (length(Nvec) != nrow(State)) stop("Host abundance and PathogenState lengths differ")
    for (i in seq_len(nrow(State))) {
      if (Nvec[i] < old[i]) State[i, ] <- thin_row(State[i, ], Nvec[i])
      else if (Nvec[i] > old[i]) State[i, "S"] <- State[i, "S"] + (Nvec[i] - old[i])
    }
    storage.mode(State) <- "integer"
    State
  }

  step <- function(State, N, timestep, context) {
    validate_context(context)
    n_units <- nrow(State); Nvec <- as.integer(as.numeric(N))
    beta <- pmax(0, resolve(Spec$Beta, timestep, context, "Beta"))
    rec <- resolve(Spec$RecoveryProb, timestep, context, "RecoveryProb")
    mort <- resolve(Spec$PathogenMortalityProb, timestep, context, "PathogenMortalityProb")
    prog <- resolve(Spec$ProgressionProb, timestep, context, "ProgressionProb")
    waning <- resolve(Spec$ImmunityLossProb, timestep, context, "ImmunityLossProb")
    intro_p <- resolve(Spec$IntroductionProb, timestep, context, "IntroductionProb")
    intro_n <- as.integer(resolve(Spec$IntroductionNumber, timestep, context, "IntroductionNumber"))
    density_scale <- resolve(Spec$DensityScale, timestep, context, "DensityScale")
    if (any(rec + mort > 1 + 1e-12)) stop("RecoveryProb + PathogenMortalityProb must not exceed 1")

    S0 <- as.integer(State[, "S"]); I0 <- as.integer(State[, "I"])
    E0 <- if ("E" %in% States) as.integer(State[, "E"]) else integer(n_units)
    R0 <- if ("R" %in% States) as.integer(State[, "R"]) else integer(n_units)
    live <- S0 + E0 + I0 + R0
    contact <- contact_matrix(timestep, context) # rows = infectious source, columns = recipient target
    infectious_pressure <- as.numeric(crossprod(I0, contact))
    if (Transmission == "frequency") {
      contact_population <- as.numeric(crossprod(live, contact))
      prevalence_pressure <- ifelse(contact_population > 0, infectious_pressure / contact_population, 0)
      foi <- beta * prevalence_pressure
    } else {
      foi <- beta * infectious_pressure / density_scale
    }
    p_inf <- clip01(-expm1(-pmax(0, foi)))
    new_inf <- rbinom(n_units, S0, p_inf)

    introduced <- integer(n_units)
    event <- rbinom(n_units, 1L, intro_p)
    idx <- which(event > 0 & S0 - new_inf > 0)
    if (length(idx)) introduced[idx] <- pmin(intro_n[idx], S0[idx] - new_inf[idx])

    progress <- if ("E" %in% States) rbinom(n_units, E0, prog) else integer(n_units)
    stay_i <- recover <- deaths <- integer(n_units)
    for (i in which(I0 > 0)) {
      z <- as.numeric(rmultinom(1L, I0[i], c(1 - rec[i] - mort[i], rec[i], mort[i])))
      stay_i[i] <- z[1]; recover[i] <- z[2]; deaths[i] <- z[3]
    }
    lose <- if ("R" %in% States) rbinom(n_units, R0, waning) else integer(n_units)

    State[, "S"] <- S0 - new_inf - introduced + lose
    if ("E" %in% States) {
      State[, "E"] <- E0 - progress + new_inf + introduced
      State[, "I"] <- stay_i + progress
    } else State[, "I"] <- stay_i + new_inf + introduced
    if ("R" %in% States) State[, "R"] <- R0 - lose + recover
    else State[, "S"] <- State[, "S"] + recover

    Nvec <- Nvec - deaths
    if (any(Nvec < 0) || any(rowSums(State) != Nvec)) stop("Internal pathogen-state conservation error")
    list(State = State, N = Nvec, Deaths = deaths,
         Introduced = introduced, NewInfections = new_inf)
  }

  Engine <- list(States = States, Validate = validate_context, Resolve = resolve,
                 ContactMatrix = contact_matrix, Initial = initial,
                 Reconcile = reconcile, Step = step)
  structure(c(Spec, list(States = States, Engine = Engine)), class = "INApestPathogen")
}

INApestPathogenContactMatrix <- function(NodeContact, LandUseMixing = NULL) {
  NodeContact <- as.matrix(NodeContact)
  if (nrow(NodeContact) != ncol(NodeContact) || any(!is.finite(NodeContact)) || any(NodeContact < 0))
    stop("NodeContact must be a finite non-negative square source x target matrix")
  if (is.null(LandUseMixing)) return(NodeContact)
  LandUseMixing <- as.matrix(LandUseMixing)
  if (nrow(LandUseMixing) != ncol(LandUseMixing) || any(!is.finite(LandUseMixing)) || any(LandUseMixing < 0))
    stop("LandUseMixing must be a finite non-negative square source x target matrix")
  # MLU state is flattened from nodes x land uses with node varying fastest.
  # Kronecker order therefore gives source/target blocks by land-use pair.
  kronecker(LandUseMixing, NodeContact)
}

INApestPathogenOutputs <- function(PathogenStateResults) {
  d <- dim(PathogenStateResults)
  if (length(d) == 4L) {
    st <- dimnames(PathogenStateResults)[[2]]
    get_state <- function(nm) if (nm %in% st) PathogenStateResults[, nm, , , drop = FALSE] else NULL
    I <- get_state("I")
    return(list(S = get_state("S"), E = get_state("E"), I = I, R = get_state("R"),
                PathogenPresence = I > 0,
                InfectiousNodeCount = apply(I > 0, c(3, 4), sum),
                TotalInfectious = apply(I, c(3, 4), sum)))
  }
  if (length(d) == 5L) {
    st <- dimnames(PathogenStateResults)[[3]]
    get_state <- function(nm) if (nm %in% st) PathogenStateResults[, , nm, , , drop = FALSE] else NULL
    I <- get_state("I")
    return(list(S = get_state("S"), E = get_state("E"), I = I, R = get_state("R"),
                PathogenPresence = I > 0,
                InfectiousNodeCount = apply(I > 0, c(4, 5), sum),
                TotalInfectious = apply(I, c(4, 5), sum)))
  }
  stop("Expected Meta or MLU PathogenStateResults array")
}


INApestMetaParallelMultipleLandUse = function(
ModelName,                                      # Model and output name
Nperm,                                          # Number of stochastic simulation runs
Ntimesteps,                                     # Timesteps in each simulation
LocalDynamics = local.dynamicsLU,               # Local population growth, movement and management function
LocalDynamicsArgs = list(),                     # Named extra arguments passed to LocalDynamics
Pathogen = NULL,                                # Optional pathogen process specification
Nlanduses,                                      # Number of land-use classes
DetectionProb,                                  # Background-surveillance detection probability
DetectionSD = NULL,                             # Variation in background detection probability
ManageProb,                                     # Management probability when information is available
ManageSD = NULL,                                # Variation in management probability
MortalityProb,                                  # Management-driven host mortality probability
MortalitySD = NULL,                             # Variation in management mortality probability
FecundityReduction = 0,                         # Proportional reduction in reproduction under management
SpreadReduction,                                # Proportional reduction in dispersal under management
SpreadReductionSD = NULL,                       # Variation in spread reduction
InitialPopulation = NA,                         # Starting host abundance by node x land use
InitBioP = NA,                                  # Proportion of nodes initially invaded
InvasionRisk = NA,                              # External invasion probability or node weighting
InitialInfo = NA,                               # Starting information state
InitInfoP = 0,                                  # Proportion of nodes initially with information
ExternalInfoProb = NA,                          # Information arriving from outside the modelled system
InfoRetentionProb = 1,                          # Probability existing information persists one timestep
InfoPersistenceSteps = NA,                      # Programmed timesteps information persists after evidence
EnvEstabProb = 1,                               # Environmental establishment probability
Survival = 1,                                   # Host survival probability between timesteps
K,                                              # Host carrying capacity by node
PropaguleProduction,                            # Propagules produced per host
PropaguleEstablishment,                         # Establishment probability for arriving propagules
IncursionStartPop=NA,                           # Host abundance assigned to new external incursions
SDDprob,                                        # Short-distance source-to-target dispersal
SEAM = 0,                                       # Information-transfer adjacency between nodes
LDDprob = NA,                                   # Long-distance source-to-target dispersal
LDDrate = 0,                                    # Fraction of propagules entering long-distance dispersal
OngoingExternalInvasion = F,                    # Allow new host incursions after initialisation
OngoingExternalInfo = F,                        # Allow new external information after initialisation
OutputDir = NA,                                 # Directory for saved outputs
DoPlots = TRUE,                                 # Legacy plotting option; plotting is post-processing
Cores = NULL,                                   # Number of worker processes
Seed = NULL,                                    # Random seed for reproducible simulations
ExternalPathogenStateProb = NULL,               # Pathogen-state distribution for external host arrivals
InformationAcquisition = NULL,                  # Local evidence used for information: host, pathogen or both
InfoTriggeredDetectionProb = 0,                 # Detection probability where information already exists
InfoTriggeredDetectionSD = NULL,                # Variation in information-triggered detection
ReturnResults = FALSE                           # Return the in-memory result object
)
{
if(!is.numeric(Nperm) || length(Nperm)!=1L || !is.finite(Nperm) || Nperm < 1 || Nperm != floor(Nperm)) stop("Nperm must be a positive integer")
if(!is.numeric(Ntimesteps) || length(Ntimesteps)!=1L || !is.finite(Ntimesteps) || Ntimesteps < 1 || Ntimesteps != floor(Ntimesteps)) stop("Ntimesteps must be a positive integer")
if(!is.numeric(Nlanduses) || length(Nlanduses)!=1L || !is.finite(Nlanduses) || Nlanduses < 1 || Nlanduses != floor(Nlanduses)) stop("Nlanduses must be a positive integer")

# ---------------------------------------------------------------------------
# Set up and validate the node x land-use simulation.
# ---------------------------------------------------------------------------
if(!is.function(LocalDynamics))
  stop("LocalDynamics must be a function")
if (is.null(LocalDynamicsArgs)) LocalDynamicsArgs <- list()
if (!is.list(LocalDynamicsArgs))
  stop("LocalDynamicsArgs must be a named list")
if (length(LocalDynamicsArgs) &&
    (is.null(names(LocalDynamicsArgs)) || any(!nzchar(names(LocalDynamicsArgs)))))
  stop("Every LocalDynamicsArgs entry must have a non-empty name")
if (anyDuplicated(names(LocalDynamicsArgs)))
  stop("LocalDynamicsArgs names must be unique")
# Force the argument before any parallel worker closure is created. This keeps
# the selected default or user-supplied function as an explicit model input.
force(LocalDynamics)
force(LocalDynamicsArgs)
UserLocalDynamicsArgs <- LocalDynamicsArgs
# collect user lexical bindings needed by custom LocalDynamics on PSOCK workers.
.CollectLocalDynamicsPSOCKBindings <- function(fun) {
  collected <- list(); seen <- character(0)
  find_user_binding <- function(nm, env) {
    ee <- env
    while(!identical(ee, emptyenv())) {
      enm <- environmentName(ee)
      if(isNamespace(ee) || identical(ee, baseenv()) || startsWith(enm, "package:")) return(NULL)
      if(exists(nm, envir=ee, inherits=FALSE)) return(list(found=TRUE,value=get(nm, envir=ee, inherits=FALSE)))
      ee <- parent.env(ee)
    }
    list(found=FALSE,value=NULL)
  }
  collect_fun <- function(f) {
    if(!is.function(f)) return(invisible(NULL))
    fg <- codetools::findGlobals(f, merge=FALSE)
    refs <- unique(c(fg$variables, fg$functions))
    for(nm in refs) {
      if(nm %in% seen) next
      hit <- find_user_binding(nm, environment(f))
      if(!isTRUE(hit$found)) next
      val <- hit$value
      seen <<- c(seen, nm); collected[[nm]] <<- val
      if(is.function(val)) collect_fun(val)
    }
    invisible(NULL)
  }
  collect_fun(fun); collected
}
LocalDynamicsPSOCKBindings <- .CollectLocalDynamicsPSOCKBindings(LocalDynamics)
force(LocalDynamicsPSOCKBindings)
LocalDynamicsHasPreMortalityN <- "pre_mortality_n" %in% names(formals(LocalDynamics))
PathogenOriginal <- Pathogen
if(!is.null(PathogenOriginal) && !inherits(PathogenOriginal, "INApestPathogen")) stop("Pathogen must be NULL or an object returned by INApestPathogen()")
UsePathogen <- !is.null(PathogenOriginal)
if(UsePathogen) { PathogenEngine <- PathogenOriginal$Engine; PathogenContext <- list(n_nodes = nrow(SDDprob), n_landuses = Nlanduses, Ntimesteps = Ntimesteps); PathogenEngine$Validate(PathogenContext); force(PathogenEngine); force(PathogenContext) }

# Choose which direct local biological evidence can create/refresh HaveInfo.
# NULL preserves the legacy pathway: host detection informs, while pathogen
# detection informs only when Pathogen$DetectionTriggersInfo is TRUE.
ExplicitInformationAcquisition <- !is.null(InformationAcquisition)
if(ExplicitInformationAcquisition) {
  if(length(InformationAcquisition) != 1L || is.na(InformationAcquisition))
    stop("InformationAcquisition must be one of: host, pathogen, both")
  InformationAcquisition <- match.arg(as.character(InformationAcquisition), c("host", "pathogen", "both"))
  if(InformationAcquisition == "pathogen" && !UsePathogen)
    stop("InformationAcquisition = 'pathogen' requires Pathogen to be supplied")
  HostInformationAcquisition <- InformationAcquisition %in% c("host", "both")
  PathogenInformationAcquisition <- InformationAcquisition %in% c("pathogen", "both")
} else {
  HostInformationAcquisition <- TRUE
  PathogenInformationAcquisition <- UsePathogen && isTRUE(PathogenOriginal$DetectionTriggersInfo)
  InformationAcquisition <- if(PathogenInformationAcquisition) "both" else "host"
}
force(HostInformationAcquisition); force(PathogenInformationAcquisition); force(InformationAcquisition)

# Optional pathogen-state composition of accepted external host immigrants.
# The existing host invasion and carrying-capacity calculations remain unchanged.
ExternalPathogenStateProbResolved <- NULL
if(!is.null(ExternalPathogenStateProb)) {
  if(!UsePathogen)
    stop("ExternalPathogenStateProb requires Pathogen to be supplied")
  if(!is.numeric(ExternalPathogenStateProb) || is.null(names(ExternalPathogenStateProb)) ||
     any(!nzchar(names(ExternalPathogenStateProb))) || anyDuplicated(names(ExternalPathogenStateProb)))
    stop("ExternalPathogenStateProb must be a uniquely named numeric vector of pathogen-state probabilities")
  if(any(!names(ExternalPathogenStateProb) %in% PathogenEngine$States))
    stop("ExternalPathogenStateProb names must be pathogen states: ", paste(PathogenEngine$States, collapse = ", "))
  if(any(!is.finite(ExternalPathogenStateProb)) || any(ExternalPathogenStateProb < 0))
    stop("ExternalPathogenStateProb values must be finite and non-negative")
  ExternalPathogenStateProbResolved <- setNames(rep(0, length(PathogenEngine$States)), PathogenEngine$States)
  ExternalPathogenStateProbResolved[names(ExternalPathogenStateProb)] <- ExternalPathogenStateProb
  if(abs(sum(ExternalPathogenStateProbResolved) - 1) > 1e-10)
    stop("ExternalPathogenStateProb must sum to 1")
}
force(ExternalPathogenStateProbResolved)
  # Allow SDD and LDD connectivity to vary through time
  # static matrices obey the same source x target geometry as time-varying arrays.
  if(is.matrix(SDDprob) && nrow(SDDprob) != ncol(SDDprob)) stop("SDDprob matrix must be square")
  if(!(is.matrix(SDDprob) || length(dim(SDDprob)) == 3)) stop("SDDprob must be a square matrix or nodes x nodes x Ntimesteps array")
  if(length(dim(SDDprob)) == 3 && (dim(SDDprob)[1] != dim(SDDprob)[2] || dim(SDDprob)[3] != Ntimesteps)) stop("SDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")
  if(is.matrix(LDDprob) && (nrow(LDDprob) != nrow(SDDprob) || ncol(LDDprob) != nrow(SDDprob))) stop("LDDprob matrix must have dimensions nodes x nodes")
  if(length(dim(LDDprob)) == 3 && (dim(LDDprob)[1] != nrow(SDDprob) || dim(LDDprob)[2] != nrow(SDDprob) || dim(LDDprob)[3] != Ntimesteps)) stop("LDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")


# Allow management-induced fecundity reduction to vary by land use, node and time
FecundityReductionDims <- dim(FecundityReduction)
if(is.null(FecundityReductionDims)) {
  if(!(length(FecundityReduction) == 1 || length(FecundityReduction) == Nlanduses))
    stop("FecundityReduction must be scalar, length Nlanduses, nodes x land uses, or nodes x land uses x Ntimesteps")
} else if(length(FecundityReductionDims) == 2) {
  if(!all(FecundityReductionDims == c(nrow(SDDprob), Nlanduses)))
    stop("FecundityReduction matrix must have dimensions nodes x land uses")
} else if(length(FecundityReductionDims) == 3) {
  if(!all(FecundityReductionDims == c(nrow(SDDprob), Nlanduses, Ntimesteps)))
    stop("FecundityReduction array must have dimensions nodes x land uses x Ntimesteps")
} else {
  stop("FecundityReduction must be scalar, length Nlanduses, nodes x land uses, or nodes x land uses x Ntimesteps")
}
if(any(!is.finite(FecundityReduction)) || any(FecundityReduction < 0) || any(FecundityReduction > 1))
  stop("FecundityReduction values must be between 0 and 1")

# Resolve fecundity reduction to node x land-use values.
ResolveFecundityReductionLU <- function(timestep) {
  if(is.null(dim(FecundityReduction))) {
    if(length(FecundityReduction) == 1)
      return(matrix(FecundityReduction, nrow = nrow(SDDprob), ncol = Nlanduses))
    return(matrix(rep(FecundityReduction, each = nrow(SDDprob)), nrow = nrow(SDDprob), ncol = Nlanduses))
  }
  if(length(dim(FecundityReduction)) == 2) return(FecundityReduction)
  FecundityReduction[,,timestep]
}

# Validate targeted-surveillance inputs for node x land-use models.
.ValidateInfoTriggeredDetectionLU <- function(x, name) {
  d <- dim(x)
  if(is.null(d)) { if(!(length(x) == 1 || length(x) == Nlanduses)) stop(name, " must be scalar, length Nlanduses, nodes x land uses, or nodes x land uses x Ntimesteps")
  } else if(length(d) == 2) { if(!all(d == c(nrow(SDDprob),Nlanduses))) stop(name, " matrix must have dimensions nodes x land uses")
  } else if(length(d) == 3) { if(!all(d == c(nrow(SDDprob),Nlanduses,Ntimesteps))) stop(name, " array must have dimensions nodes x land uses x Ntimesteps")
  } else stop(name, " has unsupported dimensions")
  if(any(!is.finite(x)) || any(x < 0)) stop(name, " must contain finite non-negative values")
}
# Expand targeted detection to node x land-use values.
.ExpandInfoTriggeredDetectionLU <- function(x, timestep) {
  d <- dim(x); if(is.null(d)) { if(length(x)==1) return(matrix(x,nrow=nrow(SDDprob),ncol=Nlanduses)); return(matrix(rep(x,each=nrow(SDDprob)),nrow=nrow(SDDprob),ncol=Nlanduses)) }; if(length(d)==2) return(x); x[,,timestep]
}
.ValidateInfoTriggeredDetectionLU(InfoTriggeredDetectionProb, "InfoTriggeredDetectionProb")
if(any(InfoTriggeredDetectionProb > 1)) stop("InfoTriggeredDetectionProb values must be between 0 and 1")
if(is.null(InfoTriggeredDetectionSD)) InfoTriggeredDetectionSD <- InfoTriggeredDetectionProb/10
.ValidateInfoTriggeredDetectionLU(InfoTriggeredDetectionSD, "InfoTriggeredDetectionSD")
UseInfoTriggeredSurveillance <- any(InfoTriggeredDetectionProb != 0) || any(InfoTriggeredDetectionSD != 0)
InfoTriggeredDetectionTimeVarying <- length(dim(InfoTriggeredDetectionProb)) == 3 || length(dim(InfoTriggeredDetectionSD)) == 3
# Draw realised targeted-detection probabilities by node x land use.
.DrawInfoTriggeredDetectionLU <- function(timestep) { mu<-.ExpandInfoTriggeredDetectionLU(InfoTriggeredDetectionProb,timestep); sd<-.ExpandInfoTriggeredDetectionLU(InfoTriggeredDetectionSD,timestep); z<-matrix(rnorm(length(mu),mean=as.numeric(mu),sd=as.numeric(sd)),nrow=nrow(SDDprob),ncol=Nlanduses); z[z<0]<-0;z[z>1]<-1;z }

# Allow information retention to vary by node and through time
if(is.matrix(InfoRetentionProb) == T && (nrow(InfoRetentionProb) != nrow(SDDprob) || ncol(InfoRetentionProb) != Ntimesteps))
  stop("InfoRetentionProb matrix must have dimensions nodes x Ntimesteps")
if(is.matrix(InfoRetentionProb) == F && !(length(InfoRetentionProb) == 1 || length(InfoRetentionProb) == nrow(SDDprob)))
  stop("InfoRetentionProb must be a single value, vector of length nodes, or matrix nodes x Ntimesteps")
if(any(is.na(InfoRetentionProb)) || any(InfoRetentionProb < 0) || any(InfoRetentionProb > 1))
  stop("InfoRetentionProb values must be between 0 and 1")

# Allow programmed information persistence after last known local presence to vary by node and through time
if(is.matrix(InfoPersistenceSteps) == T && (nrow(InfoPersistenceSteps) != nrow(SDDprob) || ncol(InfoPersistenceSteps) != Ntimesteps))
  stop("InfoPersistenceSteps matrix must have dimensions nodes x Ntimesteps")
if(is.matrix(InfoPersistenceSteps) == F && !(length(InfoPersistenceSteps) == 1 || length(InfoPersistenceSteps) == nrow(SDDprob)))
  stop("InfoPersistenceSteps must be a single value, vector of length nodes, or matrix nodes x Ntimesteps")
if(any(!is.na(InfoPersistenceSteps) & (!is.finite(InfoPersistenceSteps) | InfoPersistenceSteps < 0 | InfoPersistenceSteps != floor(InfoPersistenceSteps))))
  stop("InfoPersistenceSteps values must be non-negative whole numbers or NA")
UseInfoPersistence = any(!is.na(InfoPersistenceSteps))
if(UseInfoPersistence == T && any(InfoRetentionProb < 1))
  warning("Both InfoPersistenceSteps and InfoRetentionProb specify information loss. Programmed stopping takes priority where InfoPersistenceSteps is not NA; InfoRetentionProb is only used where InfoPersistenceSteps is NA.",call. = F)

# pre-evaluate some variables for efficiency
if(length(dim(K)) <3)
{
K_is_0 <- rowSums(K)<=0
inv_K <- 1 / sum(colSums(K))
NodeK = K
Pk  = K/rowSums(K)
Pk[is.na(Pk)] = 0
}

if(length(dim(K)) ==3)
    {
    # preserve nodes x land-uses shape when slicing time-varying K.
    NodeK <- matrix(K[,,1],nrow=nrow(SDDprob),ncol=Nlanduses)
    K_is_0 <- rowSums(NodeK)<=0
    Pk <- NodeK/rowSums(NodeK)
    Pk[is.na(Pk)] <- 0
    }
  
  
if(length(dim(PropaguleProduction)) < 3)
  NodePropaguleProduction = PropaguleProduction

if(length(dim(PropaguleEstablishment)) <3)
  NodePropaguleEstablishment = PropaguleEstablishment

if(length(dim(EnvEstabProb)) <3)
  NodeEnvEstabProb <- EnvEstabProb

if(length(dim(Survival)) < 3)
  NodeSurvival <- Survival


# Validate the socioeconomic information-transfer network when supplied.
if(is.matrix(SEAM) == T)
     {
     diag(SEAM) = 0
     RandSEAM <- matrix(NA,nrow = nrow(SDDprob),ncol=nrow(SDDprob))
     }

# Assign standard deviation value to management in no value provided
if(is.null(ManageSD) == T)
	ManageSD = ManageProb/10
if(is.null(SpreadReductionSD) == T)
	SpreadReductionSD = (1-SpreadReduction)/10
if(is.null(DetectionSD) == T)
	DetectionSD = DetectionProb/10
if(is.null(MortalitySD) == T)
    MortalitySD = MortalityProb/10

  # explicit land-use surface contract.
  .ValidateLUSurface <- function(x,name,is_sd=FALSE) {
    d <- dim(x); ok <- FALSE
    if(is.null(d)) ok <- length(x) %in% c(1L,Nlanduses)
    if(!is.null(d) && length(d)==2L) ok <- all(d == c(nrow(SDDprob),Nlanduses))
    if(!is.null(d) && length(d)==3L) ok <- all(d == c(nrow(SDDprob),Nlanduses,Ntimesteps))
    if(!ok) stop(name," must be scalar, length Nlanduses, nodes x land uses, or nodes x land uses x Ntimesteps")
    if(any(!is.finite(x))) stop(name," must contain finite values")
    if(is_sd) { if(any(x < 0)) stop(name," must contain non-negative values") }
    else if(any(x < 0 | x > 1)) stop(name," values must be between 0 and 1")
    invisible(TRUE)
  }
  # Resolve a land-use surface for the current timestep.
  .ResolveLUSurface <- function(x,timestep,name,is_sd=FALSE) {
    .ValidateLUSurface(x,name,is_sd); d <- dim(x)
    if(is.null(d)) {
      if(length(x)==1L) return(matrix(x,nrow=nrow(SDDprob),ncol=Nlanduses))
      return(matrix(rep(x,each=nrow(SDDprob)),nrow=nrow(SDDprob),ncol=Nlanduses))
    }
    if(length(d)==2L) return(x)
    x[,,timestep,drop=TRUE]
  }
  # Draw stochastic values from a resolved land-use surface.
  .DrawLUSurface <- function(mu,sd,timestep,name) {
    mm <- .ResolveLUSurface(mu,timestep,name,FALSE)
    ss <- .ResolveLUSurface(sd,timestep,paste0(name,"SD"),TRUE)
    z <- matrix(rnorm(length(mm),mean=as.numeric(mm),sd=as.numeric(ss)),nrow=nrow(SDDprob),ncol=Nlanduses)
    z[z<0] <- 0; z[z>1] <- 1; z
  }
  for(nm in c("DetectionProb","ManageProb","MortalityProb","SpreadReduction")) .ValidateLUSurface(get(nm),nm,FALSE)
  for(nm in c("DetectionSD","ManageSD","MortalitySD","SpreadReductionSD")) .ValidateLUSurface(get(nm),nm,TRUE)
  # Validate scalar, node or node x timestep inputs.
  .ValidateNodeTime <- function(x,name,unit_interval=FALSE,nonnegative=TRUE) {
    d <- dim(x); ok <- FALSE
    if(is.null(d)) ok <- length(x) %in% c(1L,nrow(SDDprob))
    if(!is.null(d) && length(d)==2L) ok <- all(d == c(nrow(SDDprob),Ntimesteps))
    if(!ok) stop(name," must be scalar, length nodes, or nodes x Ntimesteps matrix")
    if(any(!is.finite(x))) stop(name," must contain finite values")
    if(nonnegative && any(x < 0)) stop(name," must contain non-negative values")
    if(unit_interval && any(x > 1)) stop(name," values must be between 0 and 1")
  }
  .ValidateNodeTime(Survival,"Survival",TRUE,TRUE)
  .ValidateNodeTime(EnvEstabProb,"EnvEstabProb",TRUE,TRUE)
  .ValidateNodeTime(PropaguleProduction,"PropaguleProduction",FALSE,TRUE)
  .ValidateNodeTime(PropaguleEstablishment,"PropaguleEstablishment",TRUE,TRUE)
  kd <- dim(K)
  if(is.null(kd) || !(length(kd) %in% c(2L,3L)) ||
     (length(kd)==2L && !all(kd==c(nrow(SDDprob),Nlanduses))) ||
     (length(kd)==3L && !all(kd==c(nrow(SDDprob),Nlanduses,Ntimesteps))))
    stop("K must have dimensions nodes x land uses or nodes x land uses x Ntimesteps")
  if(any(!is.finite(K)) || any(K < 0) || any(K != floor(K))) stop("K must contain non-negative whole-number carrying capacities")
  if(!(length(InitialPopulation)==1L && is.na(InitialPopulation))) {
    if(!is.matrix(InitialPopulation) || !all(dim(InitialPopulation)==c(nrow(SDDprob),Nlanduses))) stop("InitialPopulation must be NA or a nodes x land uses matrix")
    if(any(!is.finite(InitialPopulation)) || any(InitialPopulation < 0) || any(InitialPopulation != floor(InitialPopulation))) stop("InitialPopulation must contain non-negative whole-number abundances")
  }
  HasExplicitInitialInfo <- !(length(InitialInfo)==1L && is.na(InitialInfo))
  if(HasExplicitInitialInfo) {
    if(!is.numeric(InitialInfo) || length(InitialInfo)!=nrow(SDDprob) || any(!is.finite(InitialInfo)) || any(!InitialInfo %in% c(0,1))) stop("InitialInfo must be scalar NA or a binary vector of length nodes")
  }
  if(!is.numeric(LDDrate) || length(LDDrate)!=1L || !is.finite(LDDrate) || LDDrate < 0 || LDDrate > 1) stop("LDDrate must be one finite probability between 0 and 1")
  if(!(length(IncursionStartPop)==1L && (is.na(IncursionStartPop) || (is.finite(IncursionStartPop) && IncursionStartPop >= 0 && IncursionStartPop == floor(IncursionStartPop))))) stop("IncursionStartPop must be NA or one non-negative whole number")


# Normalise scalar uncertainty inputs across land-use classes.
if(is.null(dim(DetectionSD)) && length(DetectionSD) == 1L) DetectionSD <- rep(DetectionSD, Nlanduses)
if(is.null(dim(ManageSD)) && length(ManageSD) == 1L) ManageSD <- rep(ManageSD, Nlanduses)
if(is.null(dim(MortalitySD)) && length(MortalitySD) == 1L) MortalitySD <- rep(MortalitySD, Nlanduses)
if(is.null(dim(SpreadReductionSD)) && length(SpreadReductionSD) == 1L) SpreadReductionSD <- rep(SpreadReductionSD, Nlanduses)


###########################################################
### Start of simulation
###########################################################
    
detected_cores <- parallel::detectCores()
if (is.na(detected_cores)) detected_cores <- 2L
if (is.null(Cores)) {
  n_cores <- max(1L, min(Nperm, detected_cores - 1L))
} else {
  if(!is.numeric(Cores) || length(Cores)!=1L || !is.finite(Cores) || Cores < 1 || Cores != floor(Cores)) stop("Cores must be a positive integer or NULL")
  n_cores <- max(1L, min(Nperm, as.integer(Cores)))
}
if(!is.null(Seed)) {
  if(!is.numeric(Seed) || length(Seed)!=1L || !is.finite(Seed)) stop("Seed must be one finite number or NULL")
  set.seed(as.integer(Seed))
}

# Run one stochastic realisation. Function arguments and local helpers are
# captured in this closure, avoiding fragile manual worker export lists.
# Capture the LocalDynamicsArgs resolver in this call environment so Windows
# PSOCK workers do not depend on a helper that exists only in the master session.
LocalDynamicsArgsResolver <- .resolve_INApest_LocalDynamicsArgs
force(LocalDynamicsArgsResolver)

# ---------------------------------------------------------------------------
# Run one complete stochastic land-use history on a worker.
# ---------------------------------------------------------------------------
PermutationWorker <- function(i_perm)
  {
  # Set initial dispersal connectivity
  NodeSDDprob = SDDprob
  if(length(dim(SDDprob)) == 3)
    NodeSDDprob = SDDprob[,,1]
  NodeLDDprob = LDDprob
  if(length(dim(LDDprob)) == 3)
    NodeLDDprob = LDDprob[,,1]
    # Allocate an exact number of recruits among land-use classes.
    # Sampling is explicit so the requested total is preserved exactly.
    #
    # Allocate an exact count among land-use classes without replacement.
    SampleVector <- function(X)
      {
      Vect = vector(length=0)
      for(i  in 2:length(X))
        Vect <- c(Vect,rep((i-1),times = X[i]))
      Sample <- sample(Vect,size = X[1],replace = F)
      Out <- vector(length = length(X)-1)
      for(j in 1:length(Out))
        Out[j] <- length(Sample[Sample==j])
      return(Out)  
      }
    
    
  InvasionResultsLoop <- array(dim = c(nrow(SDDprob),Nlanduses,Ntimesteps))
  PopulationResultsLoop <- InvasionResultsLoop
  ManagingResultsLoop <- array(dim = c(nrow(SDDprob),Nlanduses,Ntimesteps))
  DetectedResultsLoop <- array(dim = c(nrow(SDDprob),Nlanduses,Ntimesteps))
  BackgroundDetectedResultsLoop <- matrix(0L,nrow=nrow(SDDprob),ncol=Ntimesteps)
  InfoTriggeredDetectedResultsLoop <- matrix(0L,nrow=nrow(SDDprob),ncol=Ntimesteps)
  BackgroundDetectionProbabilityResultsLoop <- array(0,dim=c(nrow(SDDprob),Nlanduses,Ntimesteps))
  InfoTriggeredDetectionProbabilityResultsLoop <- array(0,dim=c(nrow(SDDprob),Nlanduses,Ntimesteps))
  InformationStateBeforeSurveillanceResultsLoop <- matrix(0L,nrow=nrow(SDDprob),ncol=Ntimesteps)
  HaveInfoResultsLoop <- matrix(0L,nrow=nrow(SDDprob),ncol=Ntimesteps)
  if(UsePathogen) {
    PathogenStates <- PathogenEngine$States
    PathogenStateResultsLoop <- array(0L, dim = c(nrow(SDDprob), Nlanduses, length(PathogenStates), Ntimesteps), dimnames = list(NULL, NULL, PathogenStates, NULL))
    PathogenDetectedResultsLoop <- matrix(0L,nrow=nrow(SDDprob),ncol=Ntimesteps)
  }
# Initialise host abundance from explicit starting values or sampled invasion settings.
InitBio = matrix(ncol = Nlanduses, nrow = nrow(SDDprob))
InitBio[,] = 0
InintInfested = rep(0,times = nrow(SDDprob))

if(is.matrix(InitialPopulation) == T && nrow(InitialPopulation) == nrow(SDDprob) && ncol(InitialPopulation) == Nlanduses)
  InitBio = InitialPopulation

if(is.matrix(InitialPopulation) == F || nrow(InitialPopulation) != nrow(SDDprob) || ncol(InitialPopulation) != Nlanduses)
{
risk = NULL
if(is.matrix(InvasionRisk) == T && nrow(InvasionRisk) == nrow(SDDprob))
  risk = InvasionRisk[,1]
if(is.matrix(InvasionRisk) == F && length(InvasionRisk) == nrow(SDDprob))
  risk = InvasionRisk

if(is.na(InitBioP) == F)
  Infested = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitBioP),prob = risk)
if(is.na(InitBioP) == T && is.null(risk) == F)
  {
  Infested = rbinom(1:nrow(SDDprob),size = 1,prob = risk)
  Infested = which(Infested == 1)
  }
if(is.na(InitBioP) == T && is.null(risk) == T)
  Infested = integer(0)

if(is.na(IncursionStartPop) == T)
  InintInfested[Infested] = 1
if(is.na(IncursionStartPop) == F)
  InintInfested[Infested] = IncursionStartPop
# Find alternative to for loop
InintInfested = pmin(InintInfested,rowSums(NodeK))
    InVector = cbind(InintInfested,NodeK)
InitBio <- t(apply(InVector,1,FUN = SampleVector))
}

# Cap starting host abundance at local carrying capacity.
for(i in 1:Nlanduses)
  InitBio[,i] = apply(cbind(NodeK[,i], InitBio[,i]),MARGIN = 1,FUN = min)

# Set the working host abundance and initialise pathogen state when present.
N <- InitBio
if(UsePathogen) { PathogenStateFlat <- PathogenEngine$Initial(N, PathogenContext) }
if(sum(N) == 0 && OngoingExternalInvasion == F)
  warning("No initial populations and no future external invasions")

# Initialise response information independently of true host abundance.
# Without an initial information specification, all nodes start uninformed.
InitInfo = rep(0,times = nrow(SDDprob))
if(HasExplicitInitialInfo || (is.na(InitInfoP) == F && InitInfoP>0) || is.na(sum(ExternalInfoProb)) == F )
{
if(!HasExplicitInitialInfo)
  {
  if(length(ExternalInfoProb) == nrow(SDDprob))
    {
    if(is.na(InitInfoP) == F)
      Info = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitInfoP),prob = ExternalInfoProb)
    if(is.na(InitInfoP) == T)
      {
      Info = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb) 
      Info = which(Info == 1)
      } 
    }
  if(length(ExternalInfoProb) != nrow(SDDprob))
    {
    if(is.matrix(ExternalInfoProb) == F)
      Info = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitInfoP))
    if(is.matrix(ExternalInfoProb) == T)
      {
      # first timestep supplies weights when InitInfoP fixes initial informed share.
      if(is.na(InitInfoP) == F) Info = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitInfoP),prob = ExternalInfoProb[,1])
      if(is.na(InitInfoP) == T) { Info = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb[,1]); Info = which(Info == 1) }
      }
    }
  InitInfo[Info] = 1
  
  }
if(HasExplicitInitialInfo)
  InitInfo = InitialInfo  
}

# Draw node-level background-surveillance detection probabilities.
# unified first-timestep detection surface.
  NodeDetectionProb <- .DrawLUSurface(DetectionProb,DetectionSD,1L,"DetectionProb")
  
  
if(UseInfoTriggeredSurveillance && !InfoTriggeredDetectionTimeVarying)
  NodeInfoTriggeredDetectionProb <- .DrawInfoTriggeredDetectionLU(1L)

# Draw node-level management-adoption probabilities where information is available.
# unified first-timestep management surface.
  NodeManageProb <- .DrawLUSurface(ManageProb,ManageSD,1L,"ManageProb")
  
  
# Draw the management effect on outgoing spread.
# unified first-timestep spread-reduction surface.
  NodeSpreadReduction <- .DrawLUSurface(SpreadReduction,SpreadReductionSD,1L,"SpreadReduction")
  
  
# Draw management-driven host mortality probabilities.
# unified first-timestep mortality surface.
  NodeMortalityProb <- .DrawLUSurface(MortalityProb,MortalitySD,1L,"MortalityProb")
  
  
# Record current host/pest presence from starting abundance.
Invaded = ifelse(InitBio>0,1,0) 


# Apply initial surveillance to the starting host population.
# Accepted detections can add response information.

LUdetectionProb = 1-(1-NodeDetectionProb)^(InitBio)
InitDetection = rbinom(1:nrow(SDDprob),size = 1,prob = 1-apply(1-LUdetectionProb,1,prod))
if(HostInformationAcquisition)
  InitInfo[InitInfo == 0] = InitDetection[InitInfo == 0]
# Set the working information state for the simulation.
HaveInfo = InitInfo
InitPathogenDetection <- integer(nrow(SDDprob))
if(UsePathogen && PathogenInformationAcquisition)
  {
  PathogenDetectionP <- PathogenEngine$Resolve(PathogenOriginal$DetectionProb, 1L, PathogenContext, "DetectionProb")
  p_not <- (1-PathogenDetectionP)^PathogenStateFlat[,"I"]
  p_node <- 1-apply(matrix(p_not,nrow=nrow(SDDprob),ncol=Nlanduses),1,prod)
  InitPathogenDetection <- rbinom(nrow(SDDprob),1,p_node)
  HaveInfo[HaveInfo == 0] <- InitPathogenDetection[HaveInfo == 0]
  }

# Track the most recent timestep with accepted direct local evidence.
LastKnownPresence = rep(NA,nrow(SDDprob))
if(UseInfoPersistence == T)
  {
  InitialKnownPresence = integer(0)
  if(HostInformationAcquisition)
    InitialKnownPresence = union(InitialKnownPresence, which(InitDetection == 1))
  if(PathogenInformationAcquisition && ExplicitInformationAcquisition)
    InitialKnownPresence = union(InitialKnownPresence, which(InitPathogenDetection == 1))
  # For InformationAcquisition=NULL, preserve historical initial pathogen-clock
  # semantics exactly; later pathogen detections can still refresh the clock.
  if(length(InitialKnownPresence) > 0)
    LastKnownPresence[InitialKnownPresence] = 0
  }

 
# Declare matrices outside loop
Managing = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
N0 = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
Propagules = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
Recruits = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))  
  # run simulation

# ---------------------------------------------------------------------------
# Advance land-use biology, information and response through time.
# ---------------------------------------------------------------------------
for (timestep in 1:Ntimesteps) 
  { 
 
  # Resolve short- and long-distance connectivity for the current timestep.
  if(length(dim(SDDprob)) == 3)
    NodeSDDprob = SDDprob[,,timestep]
  if(length(dim(LDDprob)) == 3)
    NodeLDDprob = LDDprob[,,timestep]
  NodeFecundityReduction <- ResolveFecundityReductionLU(timestep)

  # Resolve time-varying environmental establishment inputs.
  if(is.matrix(EnvEstabProb) == T)
    NodeEnvEstabProb <- EnvEstabProb[,timestep]
   
  if(is.matrix(Survival) == T)
    NodeSurvival <- Survival[,timestep]
    
    
  # If carrying capacity provided as 3D array assign values for relevant timestep
    if(length(dim(K)) == 3)
      {
      # preserve nodes x land-uses shape when slicing time-varying K.
      NodeK <- matrix(K[,,timestep],nrow=nrow(SDDprob),ncol=Nlanduses)
      K_is_0 <- rowSums(NodeK)<=0
      inv_K <- 1 / sum(NodeK)
      }
    
     # Resolve propagule production for the current timestep.
  if(is.matrix(PropaguleProduction) == TRUE)
    NodePropaguleProduction = PropaguleProduction[,timestep] 
  
  if(is.matrix(PropaguleEstablishment) == TRUE)
    NodePropaguleEstablishment = PropaguleEstablishment[,timestep]
      
  # redraw detection when mean or SD is time-varying.
    if(timestep > 1L && (length(dim(DetectionProb))==3 || length(dim(DetectionSD))==3)) NodeDetectionProb <- .DrawLUSurface(DetectionProb,DetectionSD,timestep,"DetectionProb")
    
     if(UseInfoTriggeredSurveillance && InfoTriggeredDetectionTimeVarying)
    NodeInfoTriggeredDetectionProb <- .DrawInfoTriggeredDetectionLU(timestep)

  # Draw node-level management-adoption probabilities where information is available.
  # redraw management when mean or SD is time-varying.
    if(timestep > 1L && (length(dim(ManageProb))==3 || length(dim(ManageSD))==3)) NodeManageProb <- .DrawLUSurface(ManageProb,ManageSD,timestep,"ManageProb")
    
     # redraw spread reduction when mean or SD is time-varying.
    if(timestep > 1L && (length(dim(SpreadReduction))==3 || length(dim(SpreadReductionSD))==3)) NodeSpreadReduction <- .DrawLUSurface(SpreadReduction,SpreadReductionSD,timestep,"SpreadReduction")
    
     # Randomly assign timestep management-mortality probability when management applied
  # redraw mortality when mean or SD is time-varying.
    if(timestep > 1L && (length(dim(MortalityProb))==3 || length(dim(MortalitySD))==3)) NodeMortalityProb <- .DrawLUSurface(MortalityProb,MortalitySD,timestep,"MortalityProb")
    
     # Use current information to activate management at each node.
  Managing[] = rbinom(Nlanduses*nrow(SDDprob),size = 1,prob = NodeManageProb*HaveInfo)
  # Identify occupied nodes currently known to the response system.
  Detected = Invaded*HaveInfo
  
  
  
  # Adjust starting population for natural and managed mortality.
  # The authoritative host draw is unchanged. With Pathogen active, realised
  # deaths are conditionally separated by cause before LocalDynamics so gross
  # turnover cannot be hidden by later recruitment.
  NBeforeMortality <- N
  ManagementMortality = NodeMortalityProb*Managing
  N0[] = rbinom(Nlanduses*nrow(SDDprob),size = c(N),prob = NodeSurvival*(1-ManagementMortality))

  ConditionalManagementMortality = matrix(0,nrow = nrow(SDDprob),ncol = Nlanduses)
  if(UsePathogen || UseInfoPersistence == T) {
    TotalMortalityProb = 1-NodeSurvival*(1-ManagementMortality)
    ManagementMortalityCells = which((N-N0) > 0 & ManagementMortality > 0 & TotalMortalityProb > 0)
    if(length(ManagementMortalityCells) > 0)
      ConditionalManagementMortality[ManagementMortalityCells] =
        (NodeSurvival*ManagementMortality)[ManagementMortalityCells]/TotalMortalityProb[ManagementMortalityCells]
  }

  ManagementDeaths <- NaturalDeaths <- NULL
  if(UsePathogen) {
    TotalDeathsFlat <- as.integer(c(N-N0))
    ConditionalManagementMortalityFlat <- as.numeric(ConditionalManagementMortality)
    ManagementDeathsFlat <- integer(length(TotalDeathsFlat))
    CellsWithPossibleManagementDeaths <- which(TotalDeathsFlat > 0 & ConditionalManagementMortalityFlat > 0)
    if(length(CellsWithPossibleManagementDeaths) > 0)
      ManagementDeathsFlat[CellsWithPossibleManagementDeaths] <- rbinom(
        n = length(CellsWithPossibleManagementDeaths),
        size = TotalDeathsFlat[CellsWithPossibleManagementDeaths],
        prob = ConditionalManagementMortalityFlat[CellsWithPossibleManagementDeaths]
      )
    NaturalDeathsFlat <- TotalDeathsFlat-ManagementDeathsFlat
    NAfterNaturalMortality <- matrix(
      as.integer(c(N))-NaturalDeathsFlat,
      nrow = nrow(SDDprob), ncol = Nlanduses
    )
    PathogenStateFlat <- PathogenEngine$Reconcile(PathogenStateFlat, NAfterNaturalMortality, PathogenContext)
    PathogenStateFlat <- PathogenEngine$Reconcile(PathogenStateFlat, N0, PathogenContext)
    ManagementDeaths <- matrix(ManagementDeathsFlat,nrow = nrow(SDDprob),ncol = Nlanduses)
    NaturalDeaths <- matrix(NaturalDeathsFlat,nrow = nrow(SDDprob),ncol = Nlanduses)
  }

  # Track known local presence from actual management mortality. Host-only
  # runs retain the previous conditional at-least-one-kill calculation.
  if(UseInfoPersistence == T)
    {
    if(UsePathogen) {
      if(HostInformationAcquisition) {
        KnownPresence = which(rowSums(ManagementDeaths) > 0)
        if(length(KnownPresence) > 0)
          LastKnownPresence[KnownPresence] = timestep
      }
    } else {
      ManagementKillProb = 1-apply((1-ConditionalManagementMortality)^(N-N0),1,prod)
      ManagementKillProb[ManagementKillProb < 0] = 0
      ManagementKillProb[ManagementKillProb > 1] = 1
      CertainManagementKillNodes = which(ManagementKillProb >= 1)
      if(length(CertainManagementKillNodes) > 0)
        LastKnownPresence[CertainManagementKillNodes] = timestep
      PotentialManagementKillNodes = which(ManagementKillProb > 0 & ManagementKillProb < 1)
      if(length(PotentialManagementKillNodes) > 0)
        {
        ManagementKilled = rbinom(n = length(PotentialManagementKillNodes),size = 1,prob = ManagementKillProb[PotentialManagementKillNodes])
        KnownPresence = PotentialManagementKillNodes[ManagementKilled == 1]
        if(length(KnownPresence) > 0)
          LastKnownPresence[KnownPresence] = timestep
        }
      }
    }
  if(sum(N0)<=0 )
    N = N0
  Pin <-0
  Qin <- 0  
    # natural dispersal
  if(sum(N0)>0 || (LocalDynamicsHasPreMortalityN && sum(NBeforeMortality)>0)) 
  {
  CoreLocalDynamicsArgs <- list(
    sddprob = NodeSDDprob,
    nodepropaguleproduction = NodePropaguleProduction,
    nodeenvestabprob = NodeEnvEstabProb,
    n = N0,
    lddprob = NodeLDDprob,
    lddrate = LDDrate,
    k_is_0 = K_is_0,
    nodeK = NodeK,
    nodepropaguleestablishment = NodePropaguleEstablishment,
    nodespreadreduction = NodeSpreadReduction,
    managing = Managing
  )
  LocalDynamicsFormals <- names(formals(LocalDynamics))
  # Optional read-only start-of-step host abundance for custom dynamics that
  # explicitly request it. It is not passed merely because a function has ...,
  # so existing custom/default LocalDynamics calls remain unchanged.
  if(LocalDynamicsHasPreMortalityN)
    CoreLocalDynamicsArgs$pre_mortality_n <- NBeforeMortality
  LocalDynamicsAcceptsFecundityReduction <-
    "nodefecundityreduction" %in% LocalDynamicsFormals || "..." %in% LocalDynamicsFormals
  if (LocalDynamicsAcceptsFecundityReduction)
    CoreLocalDynamicsArgs$nodefecundityreduction <- NodeFecundityReduction
  else if (any(NodeFecundityReduction * Managing > 0))
    stop("Custom LocalDynamics must accept a 'nodefecundityreduction' argument (or ...) when FecundityReduction is active")
  ResolvedLocalDynamicsArgs <- LocalDynamicsArgsResolver(
      UserLocalDynamicsArgs, timestep = timestep, Ntimesteps = Ntimesteps
    )
    if (length(ResolvedLocalDynamicsArgs)) {
      duplicate_args <- intersect(names(ResolvedLocalDynamicsArgs), names(CoreLocalDynamicsArgs))
      if (length(duplicate_args))
        stop("LocalDynamicsArgs may not override INApest core LocalDynamics argument(s): ",
             paste(duplicate_args, collapse = ", "))
      LocalDynamicsFormals <- names(formals(LocalDynamics))
      unknown_args <- setdiff(names(ResolvedLocalDynamicsArgs), LocalDynamicsFormals)
      if (!("..." %in% LocalDynamicsFormals) && length(unknown_args))
        stop("Custom LocalDynamics does not accept LocalDynamicsArgs argument(s): ",
             paste(unknown_args, collapse = ", "))
      CoreLocalDynamicsArgs <- c(CoreLocalDynamicsArgs, ResolvedLocalDynamicsArgs)
    }
    N <- do.call(LocalDynamics, CoreLocalDynamicsArgs)
    # Legacy LocalDynamics returns only total host abundance by node x land use.
    # Reconcile at this event boundary so positive recruitment enters S before
    # external host immigration and pathogen transmission.
    if(UsePathogen)
      PathogenStateFlat <- PathogenEngine$Reconcile(PathogenStateFlat, N, PathogenContext)
  } 
 # Apply programmed stopping after last known local presence
NodeInfoPersistenceSteps = InfoPersistenceSteps
if(is.matrix(InfoPersistenceSteps) == T)
  NodeInfoPersistenceSteps = InfoPersistenceSteps[,timestep]
if(length(NodeInfoPersistenceSteps) == 1)
  NodeInfoPersistenceSteps = rep(NodeInfoPersistenceSteps,nrow(SDDprob))
ProgrammedInfoNodes = which(HaveInfo == 1 & !is.na(NodeInfoPersistenceSteps))
if(length(ProgrammedInfoNodes) > 0)
  {
  TimeSinceKnownPresence = timestep-LastKnownPresence
  InfoStopNodes = ProgrammedInfoNodes[is.na(LastKnownPresence[ProgrammedInfoNodes]) | TimeSinceKnownPresence[ProgrammedInfoNodes] >= NodeInfoPersistenceSteps[ProgrammedInfoNodes]]
  if(length(InfoStopNodes) > 0)
    HaveInfo[InfoStopNodes] = 0
  }

# Allow information to decay after management and spread where no programmed stop is supplied
NodeInfoRetentionProb = InfoRetentionProb
if(is.matrix(InfoRetentionProb) == T)
  NodeInfoRetentionProb = InfoRetentionProb[,timestep]
if(length(NodeInfoRetentionProb) == 1)
  NodeInfoRetentionProb = rep(NodeInfoRetentionProb,nrow(SDDprob))
InfoDecayNodes = which(HaveInfo == 1 & is.na(NodeInfoPersistenceSteps) & NodeInfoRetentionProb < 1)
if(length(InfoDecayNodes) > 0)
  HaveInfo[InfoDecayNodes] = rbinom(n = length(InfoDecayNodes),size = 1,prob = NodeInfoRetentionProb[InfoDecayNodes])

# Update info vector for any info spread (if SEAM supplied)
 # Only zero values updated here so information can refresh nodes that lost information
 if(is.matrix(SEAM) == T)
  {
  RandSEAM[] <- rbinom(n=nrow(SDDprob)^2, size=1, prob = SEAM*apply(Detected,1,max))
  InfoTransferred = ifelse(colSums(RandSEAM)>0,1,0)
  HaveInfo[HaveInfo == 0] = InfoTransferred[HaveInfo == 0]
  }
 
 
 # Add invasion resulting from colonisation from external sources.
 # The host invasion draw/distribution and capacity behaviour are left unchanged;
 # pathogen state follows only the host increase actually accepted by that logic.
 NBeforeExternalInvasion <- N
 if(OngoingExternalInvasion == T)
  {
  if(is.matrix(InvasionRisk) == F) ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk)
  if(is.matrix(InvasionRisk) == T) ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk[,timestep])
  # IncursionStartPop is TOTAL node-level abundance allocated across land uses by available capacity.
  ExternalInvasionAdd <- matrix(0,nrow=nrow(SDDprob),ncol=Nlanduses)
  ExternalStartN <- if(is.na(IncursionStartPop)) 1L else as.integer(IncursionStartPop)
  for(ii in which(ExternalInvasion == 1L & rowSums(N) == 0)) {
    AcceptedStartN <- min(ExternalStartN,as.integer(sum(NodeK[ii,])))
    if(AcceptedStartN > 0L) ExternalInvasionAdd[ii,] <- SampleVector(c(AcceptedStartN,NodeK[ii,]))
  }
  N = N + ExternalInvasionAdd
  N[N > NodeK] = NodeK[N > NodeK]
  }

   if(UsePathogen)
    {
    if(is.null(ExternalPathogenStateProbResolved)) {
      PathogenStateFlat <- PathogenEngine$Reconcile(PathogenStateFlat, N, PathogenContext)
    } else {
      ExternalBeforeFlat <- as.integer(c(NBeforeExternalInvasion))
      ExternalAfterFlat <- as.integer(c(N))
      ExternalBaseN <- pmin(ExternalBeforeFlat, ExternalAfterFlat)
      PathogenStateFlat <- PathogenEngine$Reconcile(PathogenStateFlat, ExternalBaseN, PathogenContext)
      ExternalAccepted <- pmax(0L, ExternalAfterFlat-ExternalBeforeFlat)
      for(ii in which(ExternalAccepted > 0L)) {
        ExternalByState <- as.integer(rmultinom(1L, size = ExternalAccepted[ii], prob = ExternalPathogenStateProbResolved))
        PathogenStateFlat[ii, PathogenEngine$States] <-
          PathogenStateFlat[ii, PathogenEngine$States] + ExternalByState
      }
      storage.mode(PathogenStateFlat) <- "integer"
      if(any(rowSums(PathogenStateFlat) != ExternalAfterFlat))
        stop("External host pathogen-state assignment violated S/E/I/R = N")
    }

    # Pathogen transmission/progression/recovery occurs once after host events.
    PathogenStep <- PathogenEngine$Step(PathogenStateFlat, N, timestep, PathogenContext)
    PathogenStateFlat <- PathogenStep$State
    N[] <- PathogenStep$N
    for(ss in seq_along(PathogenStates)) PathogenStateResultsLoop[,,ss,timestep] <- matrix(PathogenStateFlat[,ss], nrow = nrow(SDDprob), ncol = Nlanduses)
    PathogenDetectionP <- PathogenEngine$Resolve(PathogenOriginal$DetectionProb, timestep, PathogenContext, "DetectionProb")
    p_not <- (1-PathogenDetectionP)^PathogenStateFlat[,"I"]
    p_node <- 1-apply(matrix(p_not,nrow=nrow(SDDprob),ncol=Nlanduses),1,prod)
    PathogenDetectedNow <- rbinom(nrow(SDDprob),1,p_node)
    PathogenDetectedResultsLoop[,timestep] <- PathogenDetectedNow
    if(PathogenInformationAcquisition)
      {
      if(UseInfoPersistence == T) LastKnownPresence[PathogenDetectedNow == 1] <- timestep
      HaveInfo[HaveInfo == 0] <- PathogenDetectedNow[HaveInfo == 0]
      }
    }
 
  # Add nodes with information resulting from external sources
  if(OngoingExternalInfo == T)
    {
    if(is.matrix(ExternalInfoProb) == F)
      ExternalInfo = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb)
    if(is.matrix(ExternalInfoProb) == T)
      ExternalInfo = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb[,timestep])
    HaveInfo[HaveInfo == 0] = ExternalInfo[HaveInfo==0]
     }
 # Update host/pest presence from current abundance.
 Invaded = ifelse(N>0,1,0)
 head(Invaded)
 # Record nodes adopting management
 ManagingResultsLoop[,,timestep] = Managing
  
 # Record infested nodes
 InvasionResultsLoop[,,timestep] = Invaded

 # Record populations
 PopulationResultsLoop[,,timestep] = N

 # Two node-level surveillance streams over the land-use populations.
 InfoBeforeSurveillance = as.integer(HaveInfo != 0)
 InformationStateBeforeSurveillanceResultsLoop[,timestep] = InfoBeforeSurveillance
 BackgroundDetectionProbabilityResultsLoop[,,timestep] = NodeDetectionProb
 LUdetectionProb = 1-(1-NodeDetectionProb)^(N)
 BackgroundDetection = rbinom(1:nrow(SDDprob),size = 1,prob = 1-apply(1-LUdetectionProb,1,prod))
 InfoTriggeredDetection = integer(nrow(SDDprob))
 if(UseInfoTriggeredSurveillance)
   {
   InfoTriggeredDetectionProbabilityResultsLoop[,,timestep] = NodeInfoTriggeredDetectionProb
   LUinfoDetectionProb = 1-(1-NodeInfoTriggeredDetectionProb)^(N)
   InfoTriggeredDetection = rbinom(1:nrow(SDDprob),size = 1,prob = (1-apply(1-LUinfoDetectionProb,1,prod)) * InfoBeforeSurveillance)
   }
 BackgroundDetectedResultsLoop[,timestep] = BackgroundDetection
 InfoTriggeredDetectedResultsLoop[,timestep] = InfoTriggeredDetection
 HostDetectionEvidence = pmax(BackgroundDetection,InfoTriggeredDetection)
 if(HostInformationAcquisition && UseInfoPersistence == T)
   {
   KnownPresence = which(HostDetectionEvidence == 1)
   if(length(KnownPresence) > 0) LastKnownPresence[KnownPresence] = timestep
   }
 if(HostInformationAcquisition) HaveInfo[HaveInfo==0] = HostDetectionEvidence[HaveInfo==0]
 HaveInfoResultsLoop[,timestep] = HaveInfo
 # Legacy detection status remains persistent known-present state by land use.
 DetectedResultsLoop[,,timestep] = HaveInfo*Invaded 
 }
 return(list(Invasion = InvasionResultsLoop, Population = PopulationResultsLoop, Managing = ManagingResultsLoop, Detected = DetectedResultsLoop,
   BackgroundDetected = BackgroundDetectedResultsLoop, InfoTriggeredDetected = InfoTriggeredDetectedResultsLoop,
   BackgroundDetectionProbability = BackgroundDetectionProbabilityResultsLoop, InfoTriggeredDetectionProbability = InfoTriggeredDetectionProbabilityResultsLoop,
   InformationStateBeforeSurveillance = InformationStateBeforeSurveillanceResultsLoop, HaveInfo = HaveInfoResultsLoop,
   PathogenState = if(UsePathogen) PathogenStateResultsLoop else NULL, PathogenDetected = if(UsePathogen) PathogenDetectedResultsLoop else NULL))
}

# Use a common PSOCK/parLapply architecture across parallel INApest variants.
# Static scheduling and L'Ecuyer-CMRG worker streams support reproducible
# parallel simulations when the caller fixes the R seed.
if(n_cores == 1L)
  {
  PermutationResults <- lapply(seq_len(Nperm), PermutationWorker)
  } else {
  cluster <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(if(inherits(cluster, "cluster")) parallel::stopCluster(cluster), add = TRUE)
  if (is.null(Seed)) parallel::clusterSetRNGStream(cluster) else parallel::clusterSetRNGStream(cluster, iseed=as.integer(Seed))
  # explicitly materialise custom LocalDynamics lexical bindings on PSOCK workers.
  if(length(LocalDynamicsPSOCKBindings))
    parallel::clusterCall(cluster, function(bindings) { list2env(bindings, envir=.GlobalEnv); NULL }, LocalDynamicsPSOCKBindings)
  PermutationResults <- parallel::parLapply(cluster, seq_len(Nperm), PermutationWorker)
  parallel::stopCluster(cluster)
  cluster <- NULL
  }
InvasionResults <- simplify2array(lapply(PermutationResults, `[[`, "Invasion"), higher = TRUE)
PopulationResults <- simplify2array(lapply(PermutationResults, `[[`, "Population"), higher = TRUE)
ManagingResults <- simplify2array(lapply(PermutationResults, `[[`, "Managing"), higher = TRUE)
DetectedResults <- simplify2array(lapply(PermutationResults, `[[`, "Detected"), higher = TRUE)
BackgroundDetectedResults <- simplify2array(lapply(PermutationResults, `[[`, "BackgroundDetected"), higher = TRUE)
InfoTriggeredDetectedResults <- simplify2array(lapply(PermutationResults, `[[`, "InfoTriggeredDetected"), higher = TRUE)
InformationStateBeforeSurveillanceResults <- simplify2array(lapply(PermutationResults, `[[`, "InformationStateBeforeSurveillance"), higher = TRUE)
HaveInfoResults <- simplify2array(lapply(PermutationResults, `[[`, "HaveInfo"), higher = TRUE)
BackgroundDetectionProbabilityResults <- simplify2array(lapply(PermutationResults, `[[`, "BackgroundDetectionProbability"), higher = TRUE)
InfoTriggeredDetectionProbabilityResults <- simplify2array(lapply(PermutationResults, `[[`, "InfoTriggeredDetectionProbability"), higher = TRUE)
dim(InvasionResults) <- c(nrow(SDDprob),Nlanduses,Ntimesteps,Nperm)
dim(PopulationResults) <- c(nrow(SDDprob),Nlanduses,Ntimesteps,Nperm)
dim(ManagingResults) <- c(nrow(SDDprob),Nlanduses,Ntimesteps,Nperm)
dim(DetectedResults) <- c(nrow(SDDprob),Nlanduses,Ntimesteps,Nperm)
dim(BackgroundDetectedResults) <- c(nrow(SDDprob),Ntimesteps,Nperm)
dim(InfoTriggeredDetectedResults) <- c(nrow(SDDprob),Ntimesteps,Nperm)
dim(InformationStateBeforeSurveillanceResults) <- c(nrow(SDDprob),Ntimesteps,Nperm)
dim(HaveInfoResults) <- c(nrow(SDDprob),Ntimesteps,Nperm)
dim(BackgroundDetectionProbabilityResults) <- c(nrow(SDDprob),Nlanduses,Ntimesteps,Nperm)
dim(InfoTriggeredDetectionProbabilityResults) <- c(nrow(SDDprob),Nlanduses,Ntimesteps,Nperm)
if(UsePathogen) {
  PathogenStates <- PathogenEngine$States
  PathogenStateResults <- array(0L, dim = c(nrow(SDDprob),Nlanduses,length(PathogenStates),Ntimesteps,Nperm), dimnames = list(NULL,NULL,PathogenStates,NULL,NULL))
  PathogenDetectedResults <- array(0L,dim=c(nrow(SDDprob),Ntimesteps,Nperm))
  for(pp in seq_len(Nperm)) { PathogenStateResults[,,,,pp] <- PermutationResults[[pp]]$PathogenState; PathogenDetectedResults[,,pp] <- PermutationResults[[pp]]$PathogenDetected }
}
###########################################################
### End of Simulation
###########################################################

###########################################################
### Save results for post-hoc analyses
###########################################################
### ModelName used to generate filenames
# Use standard format for ease of reading results to produce heat maps
# Support post-processing comparisons among management scenarios.
if(is.na(OutputDir) == T)
	OutputDir = ""
FileNameStem = paste0(OutputDir,ModelName)

# These are 3D arrays with dimensions (Nodes,Timesteps,Realisations)
saveRDS(ManagingResults, paste0(FileNameStem,"InfoLargeOut.rds"))
saveRDS(PopulationResults, paste0(FileNameStem,"PopulationLargeOut.rds"))
saveRDS(InvasionResults, paste0(FileNameStem,"InvasionLargeOut.rds"))
saveRDS(DetectedResults, paste0(FileNameStem,"DetectedLargeOut.rds"))
saveRDS(BackgroundDetectedResults, paste0(FileNameStem,"BackgroundDetectedLargeOut.rds"))
saveRDS(InfoTriggeredDetectedResults, paste0(FileNameStem,"InfoTriggeredDetectedLargeOut.rds"))
saveRDS(InformationStateBeforeSurveillanceResults, paste0(FileNameStem,"InformationStateBeforeSurveillanceLargeOut.rds"))
saveRDS(HaveInfoResults, paste0(FileNameStem,"HaveInfoLargeOut.rds"))
saveRDS(BackgroundDetectionProbabilityResults, paste0(FileNameStem,"BackgroundDetectionProbabilityLargeOut.rds"))
saveRDS(InfoTriggeredDetectionProbabilityResults, paste0(FileNameStem,"InfoTriggeredDetectionProbabilityLargeOut.rds"))
if(UsePathogen) {
  saveRDS(PathogenStateResults, paste0(FileNameStem,"PathogenStateLargeOut.rds"))
  saveRDS(PathogenDetectedResults, paste0(FileNameStem,"PathogenDetectedLargeOut.rds"))
}

##########################################################
### Store node-level invasion probabilities for each timestep.
### and estimation of invasion threat to other regions
##########################################################

InvasionProb = matrix(ncol = Ntimesteps, nrow = nrow(SDDprob))
for(timestep in 1:Ntimesteps)
{
TimestepData = InvasionResults[,,timestep,,drop=FALSE]
dim(TimestepData)
NodeInvaded <- apply(TimestepData, c(1,4), max)
if(is.null(dim(NodeInvaded))) NodeInvaded = matrix(NodeInvaded,nrow=nrow(SDDprob),ncol=Nperm)
InvasionProb[,timestep] = rowSums(NodeInvaded)/Nperm
}
saveRDS(InvasionProb, paste0(FileNameStem,"InvasionProb.rds"))
if(DoPlots == T)
{
###########################################################
### Produce summary figs when processing completed
###########################################################

Title = ModelName


# Change in total population with time
# Plots of raw values for each realisation and summaries (median and 95% CI) provided

PopulationSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(PopulationSummary) = c("Realisation",   "Timestep",  "NodesInfested")

# if(is.matrix(K) == TRUE)
# inv_K <- 1 / colSums(K)
inv_K
for(perm in 1:Nperm)
{
PopulationData = PopulationResults[,,,perm]
dim(PopulationData)
NodesInfested = apply(PopulationData,3,sum)*inv_K
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesInfested)
PopulationSummary = rbind(PopulationSummary,Results)
}


Filename = paste0(FileNameStem,"PopulationRaw.png")
png(Filename)
plot(PopulationSummary$Timestep,PopulationSummary$NodesInfested,ylim = c(0,1),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Total population (proportion of K)", main = Title)

for(perm in 1:Nperm)
{
Sub = PopulationSummary[PopulationSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesInfested,col  = perm)
}
dev.off()

Quantiles = as.data.frame(aggregate(PopulationSummary$NodesInfested, by = list(PopulationSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"PopulationSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,1), xlab = "Time since incursion detected (timesteps)",
ylab = "Total population (proportion of K)", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()

# Change in number of nodes infested with time
# Plots of raw values for each realisation and summaries (median and 95% CI) provided

InvasionSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(InvasionSummary) = c("Realisation",   "Timestep",  "NodesInfested")

for(perm in 1:Nperm)
{
InvasionData = InvasionResults[,,,perm]
NodesInfested = colSums(apply(InvasionData,c(1,3),max))

Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesInfested)
InvasionSummary = rbind(InvasionSummary,Results)
}


Filename = paste0(FileNameStem,"InvasionRaw.png")
png(Filename)
plot(InvasionSummary$Timestep,InvasionSummary$NodesInfested,ylim = c(0,max(InvasionSummary$NodesInfested)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Number of nodes infested", main = Title)

for(perm in 1:Nperm)
{
Sub = InvasionSummary[InvasionSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesInfested,col  = perm)
}
dev.off()

Quantiles = as.data.frame(aggregate(InvasionSummary$NodesInfested, by = list(InvasionSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"InvasionSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Number of nodes infested", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


# Summarise the number of nodes under management through time.
# Plots of raw values for each realisation and summaries (median and 95% CI) provided

ManagingSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(ManagingSummary) = c("Realisation",   "Timestep",  "NodesManaging")

for(perm in 1:Nperm)
{
ManagingData = ManagingResults[,,,perm]
dim(ManagingData)
NodesManaging = apply(ManagingData,3,sum)/Nlanduses
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesManaging)
ManagingSummary = rbind(ManagingSummary,Results)
}



 
Filename = paste0(FileNameStem,"ManagingRaw.png")
png(Filename)
plot(ManagingSummary$Timestep,ManagingSummary$NodesManaging,ylim = c(0,max(ManagingSummary$NodesManaging)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes under management", main = Title)

for(perm in 1:Nperm)
{
Sub = ManagingSummary[ManagingSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesManaging,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(ManagingSummary$NodesManaging, by = list(ManagingSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"ManagingSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes under management", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


# Summarise occupied nodes known to the response system through time.
# Plots of raw values for each realisation and summaries (median and 95% CI) provided

DetectedSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(DetectedSummary) = c("Realisation",   "Timestep",  "NodesDetected")

for(perm in 1:Nperm)
{
DetectedData = DetectedResults[,,,perm]
dim(DetectedData)
NodesDetected = colSums(apply(DetectedData,c(1,3),max))
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesDetected)
DetectedSummary = rbind(DetectedSummary,Results)
}


Filename = paste0(FileNameStem,"DetectedRaw.png")
png(Filename)
plot(DetectedSummary$Timestep,DetectedSummary$NodesDetected,ylim = c(0,max(DetectedSummary$NodesDetected)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes pest detected", main = Title)

for(perm in 1:Nperm)
{
Sub = DetectedSummary[DetectedSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesDetected,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(DetectedSummary$NodesDetected, by = list(DetectedSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"DetectedSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes pest detected", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


# Summarise the proportion of occupied nodes known through time.
# Plots of raw values for each realisation and summaries (median and 95% CI) provided
 
DetectedProportionSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(DetectedProportionSummary) = c("Realisation",   "Timestep",  "DetectedProportion")

for(perm in 1:Nperm)
{
InvasionData = InvasionResults[,,,perm]
DetectedData = DetectedResults[,,,perm]
NodesDetected = colSums(apply(DetectedData,c(1,3),max))
NodesInvaded = colSums(apply(InvasionData,c(1,3),max))
DetectedProportion = NodesDetected/NodesInvaded
DetectedProportion[is.na(DetectedProportion)==T] = 1
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,DetectedProportion)
DetectedProportionSummary = rbind(DetectedProportionSummary,Results)
}


Filename = paste0(FileNameStem,"DetectedProportionRaw.png")
png(Filename)
plot(DetectedProportionSummary$Timestep,DetectedProportionSummary$DetectedProportion,ylim = c(0,max(DetectedProportionSummary$DetectedProportion)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Proportion of infested nodes detected", main = Title)

for(perm in 1:Nperm)
{
Sub = DetectedProportionSummary[DetectedProportionSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$DetectedProportion,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(DetectedProportionSummary$DetectedProportion, by = list(DetectedProportionSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"DetectedProportionSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Proportion infested nodes detected", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()
}
if(ReturnResults)
  {
  ResultObject <- list(ModelName=ModelName, PopulationResults=PopulationResults, InvasionResults=InvasionResults,
    ManagingResults=ManagingResults, DetectedResults=DetectedResults,
    BackgroundDetectedResults=BackgroundDetectedResults, InfoTriggeredDetectedResults=InfoTriggeredDetectedResults,
    BackgroundDetectionProbabilityResults=BackgroundDetectionProbabilityResults, InfoTriggeredDetectionProbabilityResults=InfoTriggeredDetectionProbabilityResults,
    InformationStateBeforeSurveillanceResults=InformationStateBeforeSurveillanceResults, HaveInfoResults=HaveInfoResults,
    InvasionProb=InvasionProb)
  if(UsePathogen) { ResultObject$PathogenStateResults <- PathogenStateResults; ResultObject$PathogenDetectedResults <- PathogenDetectedResults }
  class(ResultObject) <- c("INApestMetaParallelMultipleLandUse","list")
  return(invisible(ResultObject))
  }
invisible(NULL)
}

################################################################
################################################################
### End of function
################################################################
################################################################
